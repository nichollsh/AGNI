#!/usr/bin/env -S julia --color=yes --startup-file=no

# To be run from within the AGNI root directory
# Run as: julia --project=. misc/utilities/grid_agni.jl

using AGNI
using .Iterators
using LoggingExtras
using Printf
using DataStructures
using NCDatasets
using Dates

const ROOT_DIR::String = abspath(dirname(abspath(@__FILE__)), "../../")
const R_earth::Float64 = 6.371e6
const M_earth::Float64 = 5.972e24

# =============================================================================
# User config
# -------------------------------

# Base parameters
cfg_base = "res/config/structure_grid.toml"
@info "Using base config: $cfg_base"

# Mass array with custom spacing
mass_arr::Array{Float64, 1} = 10.0 .^ vcat( range(start=log10(0.5),  stop=log10(4.5),   length=7),
                                            range(start=log10(5.0),  stop=log10(10.0),  length=9)
                                          )


# Define grid
#    parameters will be varied in the same order as these keys
#    enter the least-important parameters first
grid::OrderedDict = OrderedDict{String,Array{Float64,1}}((

    # metallicities here are by MASS fraction relative to hydrogen (converted to mole below)
    "metal_C"       => 10 .^ range(start=-1.0,  stop=3.0,   step=2.0),
    # "metal_S"       => 10 .^ range(start=-1.0,  stop=3.0,     step=2.0),
    # "metal_O"       => 10 .^ range(start=-1.0,  stop=3.0,     step=2.0),

    "frac_core"     =>       range(start=0.2,   stop=0.7,   step=0.1),
    "frac_atm"      =>       range(start=0.00,  stop=0.15,  step=0.03),

    # "Teff"          =>       range(start=2500,  stop=6000,  step=700.0),
    "instellation"  => 10 .^ range(start=log10(1.0),  stop=log10(2500.0),  length=5), # S_earth

    "mass_tot"      =>       mass_arr,  # M_earth
))

# Variables to record
output_keys =  ["succ", "flux_loss",
                "p_surf", "t_surf", "r_surf", "μ_surf",
                "t_phot", "r_phot", "μ_phot", "g_phot",
                "Kzz_max", "conv_ptop", "conv_pbot",]

# Grid management options
save_netcdfs = false        # NetCDF file for each case
save_plots   = false        # plots for each case
save_ncdf_tp = true         # a single NetCDF containing all T(p) solutions

# Runtime options
AGNI.solver.ls_increase= 1.02
modwrite::Int          = 1           # frequency to write CSV file
modplot::Int           = 2            # Plot during runtime (debug)
frac_min::Float64      = 0.001        # 0.001 -> 1170 bar for Earth
frac_max::Float64      = 1.0
transspec_p::Float64   = 2e3    # Pa
fc_floor::Float64      = 500.0   # K

# =============================================================================
# Parse keys and flatten grid
# -------------------------------

# Parse config file
cfg::Dict = AGNI.open_config(joinpath(ROOT_DIR,cfg_base))

# Output folder
output_dir = joinpath(ROOT_DIR, cfg["files"]["output_dir"])
@info "Output folder: $output_dir"
rm(output_dir,force=true,recursive=true)
mkdir(output_dir)

# Backup config to output dir
cp(joinpath(ROOT_DIR,cfg_base), joinpath(output_dir,"base_config.toml"))

# Results path
result_table_path::String = joinpath(output_dir,"result_table.csv")
result_emits_path::String = joinpath(output_dir,"result_emits.csv")
result_profs_path::String = joinpath(output_dir,"result_profs.nc")

# Logging
AGNI.setup_logging(joinpath(output_dir, "manager.log"), cfg["execution"]["verbosity"])

# Parse parameters
incl_convect     = cfg["physics"]["convection"]
incl_conduct     = cfg["physics"]["conduction"]
incl_sens        = cfg["physics"]["sensible_heat"]
incl_latent      = cfg["physics"]["latent_heat"]
sol_type         = cfg["execution"]["solution_type"]
conv_atol        = cfg["execution"]["converge_atol"]
conv_rtol        = cfg["execution"]["converge_rtol"]
perturb_all      = cfg["execution"]["perturb_all"]
plt_tmp          = cfg["plots"]["temperature"]
plt_ani          = cfg["plots"]["animate"]
p_top            = cfg["composition"]["p_top"]
chem             = cfg["physics"]["chemistry"]
roughness        = cfg["planet"]["roughness"]
windspeed        = cfg["planet"]["wind_speed"]
condensates      = cfg["composition"]["condensates"]
metallicities    = cfg["composition"]["metallicities"]
roughness        = cfg["planet"]["roughness"]
wind_speed       = cfg["planet"]["wind_speed"]
flux_int         = cfg["planet"]["flux_int"]
surface_mat      = cfg["planet"]["surface_material"]
p_surf           = cfg["composition"]["p_surf"]
star_Teff        = cfg["planet"]["star_Teff"]
stellar_spectrum = cfg["files"]["input_star"]
nlev_c           = cfg["execution"]["num_levels"]
grey_lw          = cfg["physics"]["grey_lw"]
grey_sw          = cfg["physics"]["grey_sw"]

# Get number of spectral bands
bands::Int = 1
if cfg["files"]["input_sf"] != "greygas"
    bands = parse(Int,split(cfg["files"]["input_sf"],"/")[end-1])
end
@info "Spectral bands: $bands"

# Intial values for interior structure
radius   = cfg["planet"]["radius"]
mass     = cfg["planet"]["mass"] # interior mass
gravity  = phys.grav_accel(mass, radius)
mass_tot  = mass * 1.1
frac_core = 0.325
frac_atm  = 0.01
mf_dict = Dict("H2"=>0.6, "H2O"=>0.1, "CO2"=>0.1, "N2"=>0.1, "H2S"=>0.1)

# Get the keys from the grid dictionary
input_keys = collect(keys(grid))

# Warn about size of grid
gz_est = prod([length(_v) for _v in values(grid)])
if gz_est > 1e6
    @info "Grid has $(gz_est/1e6)M points"
elseif gz_est > 1e3
    @info "Grid has $(gz_est/1e3)k points"
else
    @info "Grid has $gz_est points"
end
rt_est = gz_est * 5.0 # seconds
if rt_est > 60*60
    rt_est /= 60*60 # hrs
    if rt_est > 24
        rt_est /= 24 # days
        @info "Single-core runtime will be approximately $(rt_est) days"
    else
        @info "Single-core runtime will be approximately $(rt_est) hours"
    end
else
    @info "Single-core runtime will be approximately $(rt_est) seconds"
end


# Tidy grid
@info "Grid axes:"
#    round total mass to 2dp
if "mass_tot" in keys(grid)
    grid["mass_tot"] = round.(grid["mass_tot"]; digits=2)
end
#    limit atmosphere mass fraction
if "frac_atm" in keys(grid)
    grid["frac_atm"] = clamp.(grid["frac_atm"], frac_min, frac_max)
end
#    limit core mass fraction
if "frac_core" in keys(grid)
    grid["frac_core"] = clamp.(grid["frac_core"], frac_min, frac_max)
end
#    round instellation to 1dp
#    sort in descending order
if "instellation" in keys(grid)
    grid["instellation"] = reverse(sort(round.(grid["instellation"]; digits=1)))
end

# Print  gridpoints for user
for k in keys(grid)
    grid[k] = collect(Float64, grid[k])
    @info "  $k : $(grid[k])"
end

# Create a vector of dictionaries for all parameter combinations
@info "Flattening grid"
grid_flat = [
    OrderedDict(zip(input_keys, values_combination))
    for values_combination in Iterators.product(values(grid)...)
]
gridsize = length(grid_flat)
numfails = 0

# Write combinations to file
@info "Writing flattened grid to file"
gpfile = joinpath(output_dir,"gridpoints.csv")
@info "    $gpfile"
open(gpfile, "w") do hdl
    # Header
    head = "index," * join(input_keys, ",") * "\n"
    write(hdl,head)

    # Each row
    for (i,p) in enumerate(grid_flat)
        row = @sprintf("%08d,",i) * join([@sprintf("%.6e",v) for v in values(p)], ",") * "\n"
        write(hdl,row)
    end
end

# Create output variables to record results in
result_table::Array{Dict,    1}  = [Dict{String,Float64}(k => Float64(-1.0) for k in output_keys)   for _ in 1:gridsize]
result_profs::Array{Dict,    1}  = [Dict{String,Array}()                    for _ in 1:gridsize] # array of dicts (p, t, r)
result_emits::Array{Float64, 2}  = zeros(Float64, (gridsize,bands)) # array of fluxes

wlarr::Array{Float64,1} = zeros(Float64, bands)

@info "Generated grid of $(length(input_keys)) dimensions, with $gridsize points"
sleep(3)

# Write results table to disk
function write_table(res_tab::Array{Dict,1}, fpath::String, nrows::Int)

    global output_keys

    @info "Writing results table '$fpath' (nrows=$nrows)"

    # Remove old file
    if isfile(fpath)
        rm(fpath)
    end

    # Write file
    open(fpath, "w") do hdl
        # Header
        head = "index," * join(output_keys, ",") * "\n"
        write(hdl,head)

        # Each row
        for i in 1:nrows
            row = @sprintf("%08d",i)
            for k in output_keys
                row *= @sprintf(",%.6e",res_tab[i][k])
            end
            write(hdl,row*"\n")
        end
    end
end

# Write emission fluxes to disk
function write_emits(emi_arr::Array, fpath::String, nrows::Int)

    @info "Writing fluxes array '$fpath' (nrows=$nrows)"

    global wlarr

    # Remove old file
    if isfile(fpath)
        rm(fpath)
    end

    # Write file
    open(fpath, "w") do hdl
        # Header of wavelength values
        head = "index"
        for b in 1:bands
            head *= @sprintf(",%.6e",wlarr[b])
        end
        head *= "\n"
        write(hdl,head)

        # Each row
        for i in 1:nrows
            row = @sprintf("%08d",i)
            for b in 1:bands
                row *= @sprintf(",%.6e",emi_arr[i,b])
            end
            write(hdl,row*"\n")
        end
    end
end

# Write P-T-R profiles to NetCDF file
function write_profs(res_pro::Array, fpath::String)

    if isfile(fpath)
        rm(fpath)
    end

    ds = Dataset(fpath,"c")

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
            var_p[j,i] = res_pro[i]["p"][j]
            var_t[j,i] = res_pro[i]["t"][j]
            var_r[j,i] = res_pro[i]["r"][j]
        end
    end

    close(ds)
end


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
                                κ_grey_lw=grey_lw,
                                κ_grey_sw=grey_sw,
                                metallicities=metallicities,
                                flag_gcontinuum   = cfg["physics"]["continua"],
                                flag_rayleigh     = cfg["physics"]["rayleigh"],
                                flag_cloud        = cfg["physics"]["cloud"],
                                overlap_method    = cfg["physics"]["overlap_method"],
                                real_gas          = cfg["physics"]["real_gas"],
                                thermo_functions  = cfg["physics"]["thermo_funct"],
                                use_all_gases     = true,
                                surf_roughness=roughness, surf_windspeed=windspeed,
                                fastchem_floor = fc_floor,
                                Kzz_floor = 0.0,
                                flux_int=flux_int,
                                surface_material=surface_mat,
                                mlt_criterion=only(cfg["physics"]["convection_crit"][1]),
                        )

# AGNI struct, allocate
atmosphere.allocate!(atmos, stellar_spectrum; stellar_Teff=star_Teff)
atmos.transspec_p = transspec_p

# Check ok
if !atmos.is_alloc
    @error "Could not allocate atmosphere struct"
    exit(1)
end

wlarr[:] .= atmos.bands_cen[:]

# Set PT
setpt.request!(atmos, cfg["execution"]["initial_state"])

# =============================================================================
# Iterate over parameters
# -------------------------------

"""
Calc interior radius as a function of: interior mass, interior core mass fraction.
All in Earth units, derived from Zalmoxis.
"""
function calc_Rint(m, c)
    # fit coefficients
    m0 =  1.2034502662
    m1 =  0.2638026977
    c0 = -0.2115696893
    c1 =  1.992280927
    e0 = -0.1028476164
    e1 =  0.5909898648
    o1 = -0.1505066123

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

# Get start time [seconds]
time_start::Float64 = time()
@info "Start time: $(now())"

# Run the grid of models in series
succ = true
succ_last = true
for (i,p) in enumerate(grid_flat)
    @info @sprintf("Grid point %d / %-d (%2.1f%%)",i,gridsize,100*i/gridsize)

    global succ_last
    global succ
    succ_last = succ

    # Set all VMRs to zero
    # for gas in atmos.gas_names
    #     atmos.gas_vmr[gas][:] .= 0.0
    #     atmos.gas_ovmr[gas][:] .= 0.0
    # end

    # Update parameters
    for (k,val) in p
        @info "    set $k = $val"
        if k == "p_surf"
            atmos.p_oboa = val * 1e5 # convert bar to Pa
            atmos.p_boa = atmos.p_oboa
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

            # metallicity key is by mass frac, but atmosphere stores value by mol frac
            # convert these via scaling with factor mu_H/mu_gas
            gas = String(split(k,"_")[2])
            atmos.metal_orig[gas] = val * phys._get_mmw("H") / phys._get_mmw(gas)

            # remove FC input file to force update
            rm(atmos.fastchem_elem, force=true)

        else
            @error "Unhandled input parameter: $k"
            exit(1)
        end
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

    @info @sprintf("    using p_surf = %.2e bar",atmos.pl[end]/1e5)

    # Set temperature array based on interpolation from last solution
    max_steps = Int(cfg["execution"]["max_steps"])
    easy_start = Bool(cfg["execution"]["easy_start"])
    if i > 1
        if succ_last
            # last iter was successful
            setpt.fromarrays!(atmos, result_profs[i-1]["p"], result_profs[i-1]["t"])
        else
            # last iter failed -> restore initial guess for T(p)
            setpt.request!(atmos, cfg["execution"]["initial_state"])
            easy_start = true
        end
    end

    # Solve for RCE
    succ = solver.solve_energy!(atmos, sol_type=sol_type,
                                            conduct=incl_conduct, chem=chem,
                                            convect=incl_convect, latent=incl_latent,
                                            sens_heat=incl_sens,
                                            max_steps=max_steps,
                                            max_runtime=Float64(cfg["execution"]["max_runtime"]),
                                            conv_atol=conv_atol,
                                            conv_rtol=conv_rtol,
                                            method=1,
                                            rainout=Bool(cfg["physics"]["rainout"]),
                                            oceans=Bool(cfg["physics"]["oceans"]),
                                            dx_max=Float64(cfg["execution"]["dx_max"]),
                                            ls_method=Int(cfg["execution"]["linesearch"]),
                                            easy_start=easy_start,
                                            modplot=modplot,
                                            save_frames=false,
                                            radiative_Kzz=false,
                                            perturb_all=perturb_all,

                                            )

    # Report radius
    @info @sprintf("    found r_phot = %.3f R_earth",atmos.transspec_r/R_earth)

    # Write NetCDF file for this case
    if save_netcdfs
        save.write_ncdf(atmos, joinpath(atmos.OUT_DIR,@sprintf("%08d.nc",i)))
    end

    # Make plot for this case
    if save_plots
        plotting.plot_pt(atmos, joinpath(atmos.OUT_DIR,@sprintf("%08d_pt.png",i)))
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
        elseif k == "flux_loss"
            result_table[i][k] = maximum(abs.(atmos.flux_tot)) - minimum(abs.(atmos.flux_tot))

        elseif k == "p_surf"
            result_table[i][k] = atmos.p_boa
        elseif k == "t_surf"
            result_table[i][k] = atmos.tmp_surf
        elseif k == "r_surf"
            result_table[i][k] = atmos.rp
        elseif k == "μ_surf"
            result_table[i][k] = atmos.layer_μ[end]
        elseif k == "g_surf"
            result_table[i][k] = atmos.grav_surf

        elseif k == "t_phot"
            result_table[i][k] = atmos.transspec_tmp
        elseif k == "r_phot"
            result_table[i][k] = atmos.transspec_r
        elseif k == "μ_phot"
            result_table[i][k] = atmos.transspec_μ
        elseif k == "g_phot"
            result_table[i][k] = atmos.transspec_grav

        elseif k == "Kzz_max"
            result_table[i][k] = maximum(atmos.Kzz)
        elseif k == "conv_ptop"
            result_table[i][k] = atmosphere.estimate_convective_zone(atmos)[1]
        elseif k == "conv_pbot"
            result_table[i][k] = atmosphere.estimate_convective_zone(atmos)[2]
        else
            @error "Unhandled output variable: $k"
            exit(1)
        end
    end

    # Record fluxes
    result_emits[i,:] .= atmos.band_u_lw[1, :] .+ atmos.band_u_sw[1, :]

    # Record profile (also in SI)
    result_profs[i] = Dict("p"=>deepcopy(atmos.p),
                            "t"=>deepcopy(atmos.tmp),
                            "r"=>deepcopy(atmos.r)
                        )


    # Update results file on the go?
    if mod(i,modwrite) == 0
        write_table(result_table, result_table_path, i)
        write_emits(result_emits, result_emits_path, i)
    end

    # Iterate
    @info "  "
end

@info "===================================="
@info " "

# Get end time
time_end::Float64 = time()
@info "Finish time: $(now())"

# Tidy up
atmosphere.deallocate!(atmos)

# =============================================================================
# Write result_table to CSV, result_profs to NetCDF, and exit
# -------------------------------

if numfails >0
    @warn "Number of failed grid points: $numfails"
else
    @info "No failures recorded"
end

# Print model statistics
duration = (time_end - time_start)
@info "Average iteration duration: $(duration/gridsize) seconds"

duration /= 60 # minutes
if duration > 60
    @info "Total runtime: $(duration/60) hours"
else
    @info "Total runtime: $duration minutes"
end

# Write results table
@info "Writing final results table to CSV..."
write_table(result_table, result_table_path, gridsize)

@info "Writing final fluxes array to CSV..."
write_emits(result_emits, result_emits_path, gridsize)

# Write NetCDF of profiles
@info "Writing final profiles to NetCDF.."
write_profs(result_profs, result_profs_path)

# Done
@info "Done!"
