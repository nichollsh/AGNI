using AGNI
using .Iterators
using LoggingExtras
using Printf
using NCDatasets

const ROOT_DIR::String = abspath(dirname(abspath(@__FILE__)), "../")
const R_earth::Float64 = 6.371e6
const M_earth::Float64 = 5.972e24

# =============================================================================
# User config
# -------------------------------

# Base parameters
cfg_base = "res/config/structure_grid.toml"
@info "Using base config: $cfg_base"
cfg::Dict = AGNI.open_config(joinpath(ROOT_DIR,cfg_base))

# Define grid
grid::Dict = Dict((
    "mass_tot"      =>       range(start=1.00,  stop=10.00, length=4),  # M_earth
    "frac_atm"      =>       range(start=0.01,  stop=0.10,  length=4),
    "frac_core"     =>       range(start=0.10,  stop=0.80,  length=4),
    "metal_C"       => 10 .^ range(start=-3,    stop=2,     length=2),
    # "metal_S"       => 10 .^ range(start=-3,    stop=2,     length=2),
    # "metal_O"       => 10 .^ range(start=-3,    stop=2,     length=2),
    "instellation"  => 10 .^ range(start=-0.5,  stop=3.5,   length=3)
))

# Variables to record
output_keys = ["succ", "p_surf", "t_surf", "r_surf", "μ_surf", "r_phot", "μ_phot", "t_phot",  "Kzz_max"]

# Grid management options
save_netcdfs = false        # NetCDF file for each case
save_plots   = false        # plots for each case
save_ncdf_tp = true         # a NetCDF containing all T(p) solutions

# Runtime options
transspec_p   = 2e3    # Pa
fc_floor      = 80.0   # K

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
metallicities    = cfg["composition"]["metallicities"]
turb_coeff       = cfg["planet"]["turb_coeff"]
wind_speed       = cfg["planet"]["wind_speed"]
flux_int         = cfg["planet"]["flux_int"]
surface_mat      = cfg["planet"]["surface_material"]
p_surf           = cfg["composition"]["p_surf"]
mf_dict          = cfg["composition"]["vmr_dict"]
star_Teff        = cfg["planet"]["star_Teff"]
stellar_spectrum = cfg["files"]["input_star"]
nlev_c           = cfg["execution"]["num_levels"]

# Intial values for interior structure
radius   = cfg["planet"]["radius"]
mass     = cfg["planet"]["mass"] # interior mass
gravity  = phys.grav_accel(mass, radius)
mass_tot  = mass * 1.1
frac_core = 0.325
frac_atm  = 0.01

# Get the keys from the grid dictionary
input_keys = collect(keys(grid))

# Create a vector of dictionaries for all parameter combinations
grid_flat = [
    Dict(zip(input_keys, values_combination))
    for values_combination in Iterators.product(values(grid)...)
]
gridsize = length(grid_flat)
numfails = 0

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

# Create output variables to record results in
result_table::Array{Dict, 1}  = [Dict{String,Float64}() for _ in 1:gridsize]
result_profs::Array{Dict, 1}  = [Dict{String,Array}()    for _ in 1:gridsize] # array of dicts (p, t, r)

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
                                nlev_c,
                                p_surf,
                                p_top,
                                mf_dict, "";

                                condensates=condensates,
                                metallicities=metallicities,
                                flag_gcontinuum   = cfg["execution"]["continua"],
                                flag_rayleigh     = cfg["execution"]["rayleigh"],
                                flag_cloud        = cfg["execution"]["cloud"],
                                overlap_method    = cfg["execution"]["overlap_method"],
                                real_gas          = real_gas,
                                thermo_functions  = cfg["execution"]["thermo_funct"],
                                use_all_gases     = true,
                                C_d=turb_coeff, U=wind_speed,
                                fastchem_floor = fc_floor,
                                Kzz_floor = 0.0,
                                flux_int=flux_int,
                                surface_material=surface_mat,
                                mlt_criterion=only(cfg["execution"]["convection_crit"][1]),
                        )

# AGNI struct, allocate
atmosphere.allocate!(atmos, stellar_spectrum; stellar_Teff=star_Teff)
atmos.transspec_p = transspec_p

# Set PT
setpt.request!(atmos, cfg["execution"]["initial_state"])

# =============================================================================
# Iterate over parameters
# -------------------------------

"""
Calc interior radius as a function of: interior mass, interior core mass fraction.
    Earth units.
"""
function calc_Rint(m, c)
    # fit coefficients
    m0 = 1.10190624
    m1 = 0.27968094
    c0 = -0.20955757
    c1 = 2.07892259
    e0 = -0.08667749
    e1 = 0.66494823
    o1 = -0.0593633

    # evaluate fit
    return  m0*m^m1 + c0*c^c1 + e0*(m*c)^e1  + o1
end

"""
Update interior structure
"""
function update_structure!(atmos, mass_tot, frac_atm, frac_core)
    if frac_atm + frac_core > 1
        @warn "Core fraction ($frac_core) and atm fraction ($frac_atm) sum to >1"
    end

    # atmosphere mass
    mass_atm = mass_tot * frac_atm

    # interior mass from remainder
    atmos.interior_mass = mass_tot - mass_atm

    # get interior radius from fit to Zalmoxis
    frac_core_int = frac_core / (1-frac_atm)
    atmos.rp = calc_Rint(atmos.interior_mass/M_earth, frac_core_int) * R_earth

    # surface gravity
    atmos.grav_surf = phys.grav_accel(atmos.interior_mass, atmos.rp)

    # surface pressure
    atmos.p_boa = mass_atm * atmos.grav_surf / (4 * pi * atmos.rp^2)
    atmosphere.generate_pgrid!(atmos)
end

for (i,p) in enumerate(grid_flat)
    @info @sprintf("Grid point %d / %-d = %.1f%%",i,gridsize,100*i/gridsize)

    succ = true

    # Set all VMRs to zero
    # for gas in atmos.gas_names
    #     atmos.gas_vmr[gas][:] .= 0.0
    #     atmos.gas_ovmr[gas][:] .= 0.0
    # end

    # Update parameters
    for (k,val) in p
        @info "    set $k = $val"
        if k == "p_surf"
            atmos.p_boa = val * 1e5 # convert bar to Pa
            atmosphere.generate_pgrid!(atmos)

        elseif k == "mass"
            atmos.interior_mass = val * M_earth
            atmos.grav_surf = phys.grav_accel(val, atmos.rp)

        elseif k == "radius"
            atmos.rp = val
            atmos.grav_surf = phys.grav_accel(atmos.interior_mass, atmos.rp)

        elseif k == "frac_atm"
            global frac_atm = val
            update_structure!(atmos, mass_tot, frac_atm, frac_core)

        elseif k == "frac_core"
            global frac_core = val
            update_structure!(atmos, mass_tot, frac_atm, frac_core)

        elseif k == "mass_tot"
            global mass_tot = val * M_earth
            update_structure!(atmos, mass_tot, frac_atm, frac_core)

        elseif k == "instellation"
            atmos.instellation = val * 1361.0 # W/m^2

        elseif startswith(k, "vmr_")
            gas = split(k,"_")[2]
            atmos.gas_vmr[gas][:]  .= val

        elseif startswith(k, "metal_")
            gas = split(k,"_")[2]
            atmos.metal_orig[gas] = val

            # remove FC input file to force update
            rm(atmos.fastchem_elem, force=true)

        else
            @error "Unhandled parameter: $k"
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
        # if atmos.gas_vmr[g][end] > 0
        #     println("$g: $(atmos.gas_vmr[g][end])")
        # end
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

    # Record keys (all in SI)
    for k in output_keys
        field = Symbol(k)
        if hasfield(atmosphere.Atmos_t, field)
            result_table[i][k] = Float64(getfield(atmos, Symbol(k)))
        elseif k == "succ"
            if succ
                result_table[i][k] = 1.0 # success
            else
                global numfails += 1
                result_table[i][k] = -1.0 # failure
            end
        elseif k == "p_surf"
            result_table[i][k] = atmos.p_boa
        elseif k == "t_surf"
            result_table[i][k] = atmos.tmp_surf
        elseif k == "r_surf"
            result_table[i][k] = atmos.rp
        elseif k == "μ_surf"
            result_table[i][k] = atmos.layer_μ[end]
        elseif k == "r_phot"
            result_table[i][k] = atmos.transspec_r
        elseif k == "μ_phot"
            result_table[i][k] = atmos.transspec_μ
        elseif k == "t_phot"
            result_table[i][k] = atmos.transspec_tmp
        elseif k == "Kzz_max"
            result_table[i][k] = maximum(atmos.Kzz)
        else
            @error "Unhandled variable: $k"
            result_table[i][k]  = 0.0
        end
    end

    # Record profile (also in SI)
    result_profs[i] = Dict("p"=>deepcopy(atmos.p),
                            "t"=>deepcopy(atmos.tmp),
                            "r"=>deepcopy(atmos.r)
                        )

    # Iterate
    @info "  "
end

atmosphere.deallocate!(atmos)

# =============================================================================
# Write result_table to CSV, result_profs to NetCDF, and exit
# -------------------------------

if numfails >0
    @warn "Number of failed grid points: $numfails"
end

# Write results
@info "Writing result_table to CSV..."
open(joinpath(output_dir,"result_table.csv"), "w") do hdl
    # Header
    head = "index," * join(output_keys, ",") * "\n"
    write(hdl,head)

    # Each row
    for i in 1:gridsize
        row = @sprintf("%07d",i)
        for k in output_keys
            row *= @sprintf(",%.6e",result_table[i][k])
        end
        write(hdl,row*"\n")
    end
end

# Write NetCDF of profiles
@info "Writing result_profs to NetCDF.."
ds = Dataset(joinpath(output_dir,"result_profs.nc"),"c")
ds.attrib["description"]        = "AGNI grid profiles (TPR)"
ds.attrib["hostname"]           = gethostname()
ds.attrib["username"]           = ENV["USER"]
ds.attrib["AGNI_version"]       = atmos.AGNI_VERSION
ds.attrib["SOCRATES_version"]   = atmos.SOCRATES_VERSION

defDim(ds, "nlev_c",   atmos.nlev_c)
defDim(ds, "gridsize", gridsize)

var_p = defVar(ds, "p", Float64, ("nlev_c","gridsize",) ) # saved in python dimension order
var_t = defVar(ds, "t", Float64, ("nlev_c","gridsize",) )
var_r = defVar(ds, "r", Float64, ("nlev_c","gridsize",) )

for i in 1:gridsize
    for j in 1:nlev_c
        var_p[j,i] = result_profs[i]["p"][j]
        var_t[j,i] = result_profs[i]["t"][j]
        var_r[j,i] = result_profs[i]["r"][j]
    end
    # @info @sprintf("%d : Ri=%.3f",i,result_profs[i]["r"][end]/R_earth)
end

close(ds)

# Done
@info "Done!"
