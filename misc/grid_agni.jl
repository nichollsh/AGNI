using AGNI
using .Iterators
using LoggingExtras
using Printf

const ROOT_DIR::String = abspath(dirname(abspath(@__FILE__)), "../")

# =============================================================================
# User config
# -------------------------------

# Base parameters
cfg_base = "res/config/greygas.toml"
@info "Using base config: $cfg_base"
cfg::Dict = AGNI.open_config(joinpath(ROOT_DIR,cfg_base))

# Define grid
grid::Dict = Dict((
    "vmr_H2S" => range(start=0.05,  stop=0.95,   length=6),
    "p_surf"  => 10.0 .^ range(start=0.0, stop=4.0, length=6),
))

# Variables to record
output_keys = ["p_boa", "tmp_surf", "transspec_r", "transspec_rho", "transspec_μ", "transspec_tmp", "Kzz_kbreak"]

# Other options
save_netcdfs = false
save_plots   = false

# =============================================================================
# Parse keys and flatten grid
# -------------------------------

# Output folder
output_dir = joinpath(ROOT_DIR, cfg["files"]["output_dir"])
@info "Output folder: $output_dir"
rm(output_dir,force=true,recursive=true)
mkdir(output_dir)

# Logging
AGNI.setup_logging(joinpath(output_dir, "runner.log"), cfg["execution"]["verbosity"])

# Parse parameters
incl_convect     = cfg["execution"]["convection"]
incl_conduct     = cfg["execution"]["conduction"]
incl_sens        = cfg["execution"]["sensible_heat"]
incl_latent      = cfg["execution"]["latent_heat"]
sol_type         = cfg["execution"]["solution_type"]
conv_atol        = cfg["execution"]["converge_atol"]
conv_rtol        = cfg["execution"]["converge_rtol"]
perturb_all      = cfg["execution"]["perturb_all"]
plt_tmp          = cfg["plots"]["temperature"]
plt_ani          = cfg["plots"]["animate"]
p_top            = cfg["composition"]["p_top"]
real_gas         = cfg["execution"]["real_gas"]
chem_type        = cfg["composition"]["chemistry"]
condensates      = cfg["composition"]["condensates"]
turb_coeff       = cfg["planet"]["turb_coeff"]
wind_speed       = cfg["planet"]["wind_speed"]
flux_int         = cfg["planet"]["flux_int"]
surface_mat      = cfg["planet"]["surface_material"]
p_surf           = cfg["composition"]["p_surf"]
mf_dict          = cfg["composition"]["vmr_dict"]
star_Teff        = cfg["planet"]["star_Teff"]

radius  = cfg["planet"]["radius"]
mass    = cfg["planet"]["mass"]
gravity = phys.grav_accel(mass, radius)

# Get the keys from the grid dictionary
input_keys = collect(keys(grid))

# Create a vector of dictionaries for all parameter combinations
grid_flat = [
    Dict(zip(input_keys, values_combination))
    for values_combination in Iterators.product(values(grid)...)
]
gridsize = length(grid_flat)

# Write combinations to file
open(joinpath(output_dir,"gridpoints.csv"), "w") do hdl
    # Header
    head = "index," * join(input_keys, ",") * "\n"
    write(hdl,head)

    # Each row
    for (i,p) in enumerate(grid_flat)
        row = @sprintf("%07d,",i) * join([@sprintf("%.6e",v) for v in values(p)], ",") * "\n"
        write(hdl,row)
    end
end

# Create dictionary to record results in
result::Array{Dict, 1} = [Dict{String,Float64}() for _ in 1:gridsize]

@info "Generated grid of $(length(input_keys)) dimensions, with $gridsize points"

# =============================================================================
# Setup atmosphere object
# -------------------------------

# AGNI struct, create
atmos = atmosphere.Atmos_t()

# AGNI struct, setup
atmosphere.setup!(atmos, ROOT_DIR, output_dir,
                                String(cfg["files" ]["input_sf"]),
                                Float64(cfg["planet"]["instellation"]),
                                Float64(cfg["planet"]["s0_fact"]),
                                Float64(cfg["planet"]["albedo_b"]),
                                Float64(cfg["planet"]["zenith_angle"]),
                                Float64(cfg["planet"]["tmp_surf"]),
                                gravity, radius,
                                Int(cfg["execution"]["num_levels"]),
                                p_surf,
                                p_top,
                                mf_dict, "";

                                condensates=condensates,
                                flag_gcontinuum   = cfg["execution"]["continua"],
                                flag_rayleigh     = cfg["execution"]["rayleigh"],
                                flag_cloud        = cfg["execution"]["cloud"],
                                overlap_method    = cfg["execution"]["overlap_method"],
                                real_gas          = real_gas,
                                thermo_functions  = cfg["execution"]["thermo_funct"],
                                use_all_gases     = true,
                                C_d=turb_coeff, U=wind_speed,
                                flux_int=flux_int,
                                surface_material=surface_mat,
                                mlt_criterion=only(cfg["execution"]["convection_crit"][1]),
                        )

# AGNI struct, allocate
atmosphere.allocate!(atmos, ""; stellar_Teff=star_Teff)

# Set PT
setpt.request!(atmos, cfg["execution"]["initial_state"])

# =============================================================================
# Iterate over parameters
# -------------------------------

for (i,p) in enumerate(grid_flat)
    @info @sprintf("Grid point %d / %-d = %.1f%%",i,gridsize,100*i/gridsize)

    succ = true

    # Set all VMRs to zero
    for gas in atmos.gas_names
        atmos.gas_vmr[gas][:] .= 0.0
        atmos.gas_ovmr[gas][:] .= 0.0
    end

    # Update parameters
    for (key,val) in p
        @info "    set $key = $val"
        if key == "p_surf"
            atmos.p_boa = val * 1e5 # convert bar to Pa
            atmosphere.generate_pgrid!(atmos)

        elseif key == "mass"
            atmos.grav_surf = phys.grav_accel(val, atmos.rp)

        elseif key == "radius"
            atmos.rp = val

        elseif key == "instellation"
            atmos.instellation = val

        elseif startswith(key, "vmr_")
            gas = split(key,"_")[2]
            atmos.gas_vmr[gas][:]  .= val

        else
            @error "Unhandled parameter: $key"
            succ = false
        end
    end
    if !succ
        break
    end

    # Ensure VMRs sum to unity
    tot_vmr::Float64 = 0.0
    for i in 1:atmos.nlev_c
        # get total
        tot_vmr = 0.0
        for g in atmos.gas_names
            tot_vmr += atmos.gas_vmr[g][i]
        end

        # normalise to 1 if greater than 1, otherwise fill with H2
        if tot_vmr > 1
            for g in atmos.gas_names
                atmos.gas_vmr[g][i] /= tot_vmr
            end
        else
            atmos.gas_vmr["H2"][i] += 1-tot_vmr
        end
    end
    # set original vmr arrays
    for g in atmos.gas_names
        atmos.gas_ovmr[g][:] .= atmos.gas_vmr[g][:]
    end

    # Solve for RCE
    succ = solver.solve_energy!(atmos, sol_type=sol_type,
                                            conduct=incl_conduct, chem_type=chem_type,
                                            convect=incl_convect, latent=incl_latent,
                                            sens_heat=incl_sens,
                                            max_steps=Int(cfg["execution"]["max_steps"]),
                                            max_runtime=Float64(cfg["execution"]["max_runtime"]),
                                            conv_atol=conv_atol,
                                            conv_rtol=conv_rtol,
                                            method=1,
                                            rainout=Bool(cfg["execution"]["rainout"]),
                                            dx_max=Float64(cfg["execution"]["dx_max"]),
                                            ls_method=Int(cfg["execution"]["linesearch"]),
                                            easy_start=Bool(cfg["execution"]["easy_start"]),
                                            modplot=0,
                                            save_frames=false,
                                            perturb_all=perturb_all
                                            )

    # Write NetCDF file
    if save_netcdfs
        save.write_ncdf(atmos, joinpath(atmos.OUT_DIR,@sprintf("%07d.nc",i)))
    end

    # Make plot
    if save_plots
        plotting.plot_pt(atmos, joinpath(atmos.OUT_DIR,@sprintf("%07d_pt.png",i)))
    end

    # Record keys
    for k in output_keys
        result[i][k] = Float64(getfield(atmos, Symbol(k)))
    end

    # Iterate
    @info "  "
end

atmosphere.deallocate!(atmos)

# =============================================================================
# Write result to file, and exit
# -------------------------------

# Write results
@info "Writing results to file..."
open(joinpath(output_dir,"result.csv"), "w") do hdl
    # Header
    head = "index," * join(output_keys, ",") * "\n"
    write(hdl,head)

    # Each row
    for i in 1:gridsize
        row = @sprintf("%07d",i)
        for k in output_keys
            row *= @sprintf(",%.6e",result[i][k])
        end
        write(hdl,row*"\n")
    end
end

# Done
@info "Done!"
