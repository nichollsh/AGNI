#!/usr/bin/env -S julia --color=yes --startup-file=no

# To be run from within the AGNI root directory
# First command line argument must be the worker ID

# e.g. to run worker 1 of 3
#   julia --project=. misc/grid/worker.jl 1 3

using .Iterators
using LoggingExtras
using Printf
using DataStructures
using NCDatasets
using Dates
using AGNI

const ROOT_DIR::String = abspath(dirname(abspath(@__FILE__)), "../../")
const R_earth::Float64 = 6.371e6    # m
const M_earth::Float64 = 5.972e24   # kg
const DEFAULT_FILL::Float64 = 0.0   # fill value for arrays
const SGL_RUNTIME::Float64 = 10.0   # estimated runtime for a single gridpoint [seconds]

# =============================================================================
#                        ALL CONFIGURATION HERE
# -----------------------------------------------------------------------------

# Base parameters
const cfg_base = "res/config/struct_grid.toml"

# Mass array with custom spacing
mass_arr::Array{Float64, 1} = 10.0 .^ vcat( range(start=log10(1.0),  stop=log10(4.5),   length=8)[1:end-1],
                                            range(start=log10(4.5),  stop=log10(10.0),  length=8)
                                          )
mass_arr = reverse(sort(mass_arr))

# Define grid
#    parameters will be varied in the same order as these keys
const grid::OrderedDict = OrderedDict{String,Array{Float64,1}}((

    "frac_atm"      =>  10.0 .^ range(start=-3.0,  stop=-1.0,  length=8),

    # "frac_core"     =>  range(start=0.2,   stop=0.7,   step=0.1),
    "frac_core"     =>  Float64[0.325, 0.5, 0.7],

    "logCO"         =>  range(start=-3.0,  stop=0.0,   step=1.0),  # C/O mass ratio

    "mass_tot"      =>  mass_arr,  # M_earth

    "logZ"          =>  range(start=1.5,  stop=-1.0,   step=-0.5),  # total metallicity

    # "flux_int"      => Float64[0.0, 1.0],   # internal heat flux
    "instellation"  =>  Float64[1000.0, 300.0, 100.0, 10.0, 1.0 ], # S_earth

    "Teff"          =>  range(start=2500,  stop=5500,  step=600.0),
))

# Variables to record
const output_keys =  ["succ", "flux_loss", "r_bound",
                        "p_surf",
                        "t_surf", "r_surf", "μ_surf", "g_surf",
                        "t_phot", "r_phot", "μ_phot", "g_phot",
                        "vmr_H2", "vmr_H2O", "vmr_CO2", "vmr_CO", "vmr_O2",
                        "Kzz_max", "conv_ptop", "conv_pbot",]

# Grid management options
const save_netcdfs           = false        # NetCDF file for each case
const save_plots             = false        # plots for each case
const modwrite::Int          = 40           # Write CSV file every `modwrite` gridpoints
const modplot::Int           = 0            # Plot every `modplot` solver steps (debug)
const use_tmpdir             = true
const frac_min::Float64      = 1e-7         # 0.001 -> 1170 bar for Earth
const frac_max::Float64      = 0.999
const transspec_p::Float64   = 2e3          # Pa
const fc_floor::Float64      = 800.0       # K
const fc_wellmixed::Bool     = false      # calculate abundances as well-mixed ?
const mlt_asymptotic::Bool   = true

atmosphere.HYDROGRAV_selfg  = true
atmosphere.HYDROGRAV_constg = false

# solver.ls_increase = 1.02
solver.easy_incr  = 1/solver.easy_ini

# energy.CONVECT_MIN_PRESSURE = 1e-3 * 1e5    # 1 mbar -> Pa


# =============================================================================
#                      WORKER EXECUTION BELOW
# -----------------------------------------------------------------------------


# Define work requirements
num_work::Int = 0
id_work::Int = 0
if length(ARGS) == 2
    id_work = parse(Int, ARGS[1])
    num_work = parse(Int, ARGS[2])
else
    println(stderr, "Got invalid command line arguments: $ARGS")
    println(stderr, "Expected: ID_WORK NUM_WORK")
    exit(1)
end

# Root output folder
if (id_work < 1) || (id_work > num_work)
    println(stderr, "Invalid worker ID=$id_work")
    exit(1)
end

# Parse config file
println("Using base config: $cfg_base")
cfg::Dict = AGNI.open_config(joinpath(ROOT_DIR,cfg_base))
output_dir = joinpath(ROOT_DIR, cfg["files"]["output_dir"])

# Clean output folder
if id_work==1
    println("Creating output folder: $output_dir")
    rm(output_dir, force=true, recursive=true)
    mkdir(output_dir)
end
if !isdir(output_dir)
    println(stderr, "Could not find output folder '$output_dir'")
    exit(1)
end
output_dir = realpath(output_dir)

# Auxiliary folders
save_netcdfs && mkpath(joinpath(output_dir,"nc"))
save_plots && mkpath(joinpath(output_dir,"pl"))

# Results path and base config
cp(joinpath(ROOT_DIR,cfg_base), joinpath(output_dir,"base.toml"), force=true)

# Output dir for this particular worker
OUT_DIR = joinpath(output_dir,"wk_$id_work")
rm(OUT_DIR, recursive=true, force=true)
mkdir(OUT_DIR)

# IO dir for this worker
if haskey(ENV, "TMPDIR") && use_tmpdir
    tmpdir = realpath(ENV["TMPDIR"])
    IO_DIR = joinpath(tmpdir,"AGNIgrid","wk_$id_work")
    rm(IO_DIR, recursive=true, force=true)
    mkpath(IO_DIR)
else
    IO_DIR = OUT_DIR
end


# Setup logging ASAP
AGNI.setup_logging(
    joinpath(OUT_DIR, "wk_$(id_work).log"),
    cfg["execution"]["verbosity"]
)

@info "This process is operating worker ID=$id_work (of $num_work total)"
@info "    HOSTNAME=$(gethostname())"
@info "    OUT_DIR=$OUT_DIR"
@info "    IO_DIR=$IO_DIR"

# Output files
result_table_path::String = joinpath(OUT_DIR,"result_table.csv")
result_emits_path::String = joinpath(OUT_DIR,"result_emits.csv")
result_profs_path::String = joinpath(OUT_DIR,"result_profs.nc")

# Get number of spectral bands
bands::Int = 1
if cfg["files"]["input_sf"] != "greygas"
    bands = parse(Int,split(cfg["files"]["input_sf"],"/")[end-1])
end
@info "Spectral bands for RT: $bands"

# Shared structure variables, assigned temporary values
radius    = cfg["planet"]["radius"]
mass      = cfg["planet"]["mass"] # interior mass
gravity   = phys.grav_accel(mass, radius)
mass_tot  = mass * 1.1
frac_core = 0.325
frac_atm  = 0.01
stellar_Teff = cfg["planet"]["star_Teff"]
input_star   = cfg["files"]["input_star"]
logZ       = 2.0
logCO      = 0.0
easy_start = Bool(cfg["execution"]["easy_start"])

# Get the keys from the grid dictionary
input_keys = collect(keys(grid))

# Tidy grid
@info "Grid axes:"
#    limit atmosphere mass fraction
if "frac_atm" in keys(grid)
    grid["frac_atm"] = clamp.(grid["frac_atm"], frac_min, frac_max)
    grid["frac_atm"] = round.(grid["frac_atm"]; sigdigits=3)
end
#    limit core mass fraction
if "frac_core" in keys(grid)
    grid["frac_core"] = clamp.(grid["frac_core"], frac_min, frac_max)
    grid["frac_core"] = round.(grid["frac_core"]; sigdigits=3)
end
# round other quantities to 3dp...
for k in String["mass_tot", "Teff", "instellation"]
    if k in keys(grid)
        grid[k] = round.(grid[k]; digits=3)
    end
end

# Print gridpoints for user
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

# Warn about total grid size
gz_est = prod([length(_v) for _v in values(grid)])
if gz_est > 1e6
    @info "Grid size is $(gz_est/1e6)M points"
elseif gz_est > 1e4
    @info "Grid size is $(gz_est/1e3)k points"
else
    @info "Grid size is $gz_est points"
end

# Assign workers to grid points
#    Last worker will take the "remainder" points if chunks do not fit wholly
chunksize = round(Int,gridsize/(num_work),RoundDown)
@info "Chunk size is $chunksize points"
grid_flat[1]["worker"] = 1
for i in 2:gridsize
    if (mod(i-1,chunksize) == 0) && (i>1)
        grid_flat[i]["worker"] = min(1+grid_flat[i-1]["worker"], num_work)
    else
        grid_flat[i]["worker"] = grid_flat[i-1]["worker"]
    end
end

# Estimate worker runtime
rt_est = chunksize * SGL_RUNTIME # seconds
if rt_est > 60*60
    rt_est /= 60*60 # hrs
    if rt_est > 24
        rt_est /= 24 # days
        @info @sprintf("Worker runtime will be approx %.2f days",rt_est)
    else
        @info @sprintf("Worker runtime will be approx %.2f hours",rt_est)
    end
else
    @info @sprintf("Worker runtime will be approx %.2f seconds",rt_est)
end


# Write combinations to file
gpfile = joinpath(output_dir,"gridpoints.csv")
if (id_work == 1) || !isfile(gpfile)
    @info "Writing flattened grid to file"
    @info "    $gpfile"
    open(gpfile, "w") do hdl
        # Header
        head = "index,worker," * join(input_keys, ",") * "\n"
        write(hdl,head)

        # Each row
        for (i,p) in enumerate(grid_flat)

            # Index
            row = @sprintf("%08d,",i)

            # Worker
            row *= @sprintf("%08d,",p["worker"])

            # Other keys
            row *= join([@sprintf("%.5e",p[k]) for k in input_keys], ",") * "\n"

            # Write out
            write(hdl,row)
        end
    end
end

# Create output variables to record results in
result_table::Array{Dict,    1}  = [Dict{String,Float64}(k => Float64(DEFAULT_FILL) for k in output_keys) for _ in 1:gridsize]
result_profs::Array{Dict,    1}  = [Dict{String,Array}()  for _ in 1:gridsize] # array of dicts (p, t, r)
result_emits::Array{Float64, 2}  = fill(DEFAULT_FILL, (gridsize,bands)) # array of fluxes
wlarr::Array{Float64,1}          = fill(DEFAULT_FILL, bands)

@info "Generated grid of $(length(input_keys)) dimensions, with $gridsize points"
sleep(3)

# Write results table to disk
function write_table()

    global output_keys
    global result_table_path
    global result_table

    @info "Writing results table '$result_table_path'"

    # Remove old file
    if isfile(result_table_path)
        rm(result_table_path)
    end

    # Write file
    open(result_table_path, "w") do hdl
        # Header
        head = "index," * join(output_keys, ",") * "\n"
        write(hdl,head)

        # Each row
        for i in 1:gridsize

            # this worker?
            if grid_flat[i]["worker"] != id_work
                continue
            end

            row = @sprintf("%08d",i)
            for k in output_keys
                row *= @sprintf(",%.5e",result_table[i][k])
            end
            write(hdl,row*"\n")
        end
    end
end

# Write emission fluxes to disk
function write_emits()

    global wlarr
    global result_emits
    global result_emits_path

    @info "Writing fluxes array '$result_emits_path'"

    # Remove old file
    if isfile(result_emits_path)
        rm(result_emits_path)
    end

    # Write file
    open(result_emits_path, "w") do hdl
        # Header of wavelength values
        head = "0.0"
        for b in 1:bands
            head *= @sprintf(",%.5e",wlarr[b]*1e9) # convert to nm
        end
        head *= "\n"
        write(hdl,head)

        # Each row
        for i in 1:gridsize

            # this worker?
            if grid_flat[i]["worker"] != id_work
                continue
            end

            row = @sprintf("%08d",i)
            for b in 1:bands
                row *= @sprintf(",%.5e",result_emits[i,b])
            end
            write(hdl,row*"\n")
        end
    end
end

# Write P-T-R profiles to NetCDF file
function write_profs(nlev::Int)

    global result_profs
    global result_profs_path

    @info "Writing profiles netcdf '$result_profs_path'"

    if isfile(result_profs_path)
        rm(result_profs_path)
    end

    ds = Dataset(result_profs_path,"c")

    ds.attrib["description"]        = "AGNI grid worker profiles (t-p-r)"
    ds.attrib["hostname"]           = gethostname()
    ds.attrib["username"]           = ENV["USER"]
    ds.attrib["AGNI_version"]       = atmos.AGNI_VERSION
    ds.attrib["SOCRATES_version"]   = atmos.SOCRATES_VERSION
    ds.attrib["id_work"]            = id_work

    defDim(ds, "nlev_c",   nlev)
    defDim(ds, "gridsize", gridsize)

    var_p = defVar(ds, "p", Float64, ("nlev_c","gridsize",) ) # saved in python dimension order
    var_t = defVar(ds, "t", Float64, ("nlev_c","gridsize",) )
    var_r = defVar(ds, "r", Float64, ("nlev_c","gridsize",) )

    for i in 1:gridsize

        # this worker?
        if !haskey(result_profs[i],"p")
            continue
        end

        for j in 1:nlev
            var_p[j,i] = result_profs[i]["p"][j]
            var_t[j,i] = result_profs[i]["t"][j]
            var_r[j,i] = result_profs[i]["r"][j]
        end
    end

    close(ds)
end

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
function update_structure!(atmos, mtot, fatm, fcor)
    if fatm + fcor > 1
        @warn "Core fraction ($fcor) and atm fraction ($fatm) sum to >1"
    end

    # atmosphere mass [SI units]
    mass_atm = mtot * fatm

    # interior mass from remainder
    atmos.interior_mass = mtot - mass_atm

    # get interior radius from fit to Zalmoxis
    fcor_i = fcor / (1-fatm)
    atmos.rp = calc_Rint(atmos.interior_mass/M_earth, fcor_i) * R_earth

    # surface gravity
    atmos.grav_surf = phys.grav_accel(atmos.interior_mass, atmos.rp)

    # surface pressure
    atmos.p_oboa = mass_atm * atmos.grav_surf / (4 * pi * atmos.rp^2)
    atmos.p_boa = atmos.p_oboa
    atmosphere.generate_pgrid!(atmos)
end


"""
Update metallicity (by */H mol) from Z and C/O mass ratios
"""
function update_metallicity!(atmos, _logZ, _logCO)

    # Convert from log-scale to linear
    Z  = 10^_logZ
    CO = 10^_logCO

    # Need to work out what the X/H mass ratios are...

    # C/H by mass
    C_to_H = (CO * Z) / (1+CO)

    # O/H by mass
    O_to_H = Z - C_to_H

    # metallicity is by mass frac, but atmosphere stores value by mol frac
    atmos.metal_orig["C"] = C_to_H * phys._get_mmw("H") / phys._get_mmw("C")
    atmos.metal_orig["O"] = O_to_H * phys._get_mmw("H") / phys._get_mmw("O")
    atmos.metal_orig["H"] = 1.0
    atmos.metal_orig["S"] = 0.0

    # remove FC input file to force update
    rm(atmos.fastchem_elem, force=true)
end

"""
Initialise atmosphere object

Arguments:
 - `OUT_DIR::String`    worker output folder
 - `IO_DIR::String`    worker temp folder
"""
function init_atmos(OUT_DIR::String, IO_DIR::String)

    global radius
    global mass_tot
    global gravity
    global stellar_Teff
    global input_star

    # temp values
    mf_dict = Dict("H2"=>0.6, "H2O"=>0.1, "CO2"=>0.1, "N2"=>0.1, "H2S"=>0.1)
    gravity  = phys.grav_accel(mass_tot, radius)

    # Instantiate
    atm = atmosphere.Atmos_t()

    # AGNI struct, setup
    atmosphere.setup!(atm, ROOT_DIR, OUT_DIR,
                                    String(cfg["files" ]["input_sf"]),
                                    Float64(cfg["planet"]["instellation"]),
                                    Float64(cfg["planet"]["s0_fact"]),
                                    Float64(cfg["planet"]["albedo_b"]),
                                    Float64(cfg["planet"]["zenith_angle"]),
                                    Float64(cfg["planet"]["tmp_surf"]),
                                    gravity, radius,
                                    cfg["execution"]["num_levels"],
                                    cfg["composition"]["p_surf"],
                                    cfg["composition"]["p_top"],
                                    mf_dict, "";

                                    IO_DIR=IO_DIR,
                                    condensates=cfg["composition"]["condensates"],
                                    κ_grey_lw=cfg["physics"]["grey_lw"],
                                    κ_grey_sw=cfg["physics"]["grey_sw"],
                                    metallicities     = cfg["composition"]["metallicities"],
                                    flag_gcontinuum   = cfg["physics"]["continua"],
                                    flag_rayleigh     = cfg["physics"]["rayleigh"],
                                    flag_cloud        = cfg["physics"]["cloud"],
                                    overlap_method    = cfg["physics"]["overlap_method"],
                                    real_gas          = cfg["physics"]["real_gas"],
                                    thermo_functions  = cfg["physics"]["thermo_funct"],
                                    use_all_gases     = true,
                                    surf_roughness=cfg["planet"]["roughness"],
                                    surf_windspeed=cfg["planet"]["wind_speed"],
                                    mlt_asymptotic=mlt_asymptotic,
                                    fastchem_floor = fc_floor,
                                    fastchem_wellmixed = fc_wellmixed,
                                    Kzz_floor = 0.0,
                                    flux_int=cfg["planet"]["flux_int"],
                                    surface_material=cfg["planet"]["surface_material"],
                                    mlt_criterion=only(cfg["physics"]["convection_crit"][1]),
                            )

    # AGNI struct, allocate
    atmosphere.allocate!(atm, input_star; stellar_Teff=stellar_Teff)
    atm.transspec_p = transspec_p

    # Check ok
    if !atm.is_alloc
        @error "Could not allocate atmosphere struct"
        exit(1)
    end

    # Set PT
    atm.tmp_surf = Float64(cfg["planet"]["tmp_surf"])
    setpt.request!(atm, cfg["execution"]["initial_state"])

    return atm

end

# Initialise new atmosphere
atmos = init_atmos(OUT_DIR, IO_DIR)
wlarr[:] .= atmos.bands_cen[:]

# Get start time
@info "Operating worker $id_work of $num_work"
time_start::Float64 = time()
@info "Start time: $(now())"

# Write initial files
write_table()
write_emits()
write_profs(atmos.nlev_c)

# Run the grid of models
succ = false
succ_last = false
i_counter = 0
for (i,p) in enumerate(grid_flat)

    # i = gridpoint index
    # p = parameters at this gridpoint

    # globals, defined outside the loop
    global i_counter
    global succ
    global succ_last
    global atmos
    global mass_tot
    global radius
    global frac_core
    global frac_atm
    global gravity
    global result_emits
    global result_profs
    global result_table
    global output_dir
    global gridsize
    global save_netcdfs
    global save_plots
    global wlarr
    global stellar_Teff
    global logZ
    global logCO
    global easy_start

    # Check that this worker is assigned to this grid point, by index
    if grid_flat[i]["worker"] != id_work
        continue
    end
    i_counter += 1
    @info @sprintf("Grid point %d of %d total (%d of %d in chunk)",
                                i,gridsize,   i_counter,chunksize)

    succ_last = succ
    updated_k = false
    easy_start = false

    # Update parameters
    for (idx,(k,val)) in enumerate(p)

        # idx = index of parameter
        # k   = name of parameter
        # val = value of parameter

        # inform user
        if k == "worker"
            continue
        end

        # updating the stellar spectrum requires creating a whole new atmos object
        #   so do this first...
        if ("Teff" in keys(p)) && (idx == 1)
            if (i > 1) && (grid_flat[i-1]["Teff"] != grid_flat[i]["Teff"])
                @info "Updating Teff parameter..."
                stellar_Teff = Float64(grid_flat[i]["Teff"])

                # make new atmosphere
                atmosphere.deallocate!(atmos)
                atmos = init_atmos(OUT_DIR, IO_DIR)
                wlarr[:] .= atmos.bands_cen[:]

                easy_start = true
            end
        end

        # set other parameter...
        updated_k = (i_counter ==1) || (i==1) || Bool(grid_flat[i-1][k] != grid_flat[i][k])
        if updated_k
            @info "    updated $k = $val"
        else
            @info "    use     $k = $val"
        end
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
            frac_atm = val
            update_structure!(atmos, mass_tot, frac_atm, frac_core)

        elseif k == "frac_core"
            frac_core = val
            update_structure!(atmos, mass_tot, frac_atm, frac_core)

        elseif k == "mass_tot"
            mass_tot = val * M_earth # convert to SI
            update_structure!(atmos, mass_tot, frac_atm, frac_core)

        elseif k == "flux_int"
            atmos.flux_int = val

            # set T(p)
            if updated_k
                atmos.tmp_surf = Float64(cfg["planet"]["tmp_surf"])
                setpt.request!(atmos, cfg["execution"]["initial_state"])
                easy_start = true
            end

        elseif k == "instellation"
            atmos.instellation = val * 1361.0 # W/m^2

            # set T(p)
            if updated_k
                atmos.tmp_surf = Float64(cfg["planet"]["tmp_surf"])
                setpt.request!(atmos, cfg["execution"]["initial_state"])
                easy_start = true
            end

        elseif startswith(k, "vmr_")
            gas = split(k,"_")[2]
            atmos.gas_vmr[gas][:]  .= val

        elseif k == "logZ"
            logZ = val
            update_metallicity!(atmos, logZ, logCO)

        elseif k == "logCO"
            logCO = val
            update_metallicity!(atmos, logZ, logCO)

        elseif k == "Teff"
            # already handled


        else
            @error "Unhandled input parameter: $k"
            exit(1)
        end
    end

    @info @sprintf("    using p_surf = %.2e bar",atmos.p_boa/1e5)

    # Set temperature array based on interpolation from last solution
    max_steps   = Int(cfg["execution"]["max_steps"])
    max_runtime = Float64(cfg["execution"]["max_runtime"])
    if succ_last && (i>1) && haskey(result_profs[i-1],"p") && (i_counter != 1)
        # last iter was successful
        setpt.fromarrays!(atmos, result_profs[i-1]["p"], result_profs[i-1]["t"]; extrap=true)
        # easy_start = false
    else
        # last iter failed
        atmos.tmp_surf = Float64(cfg["planet"]["tmp_surf"])
        setpt.request!(atmos, cfg["execution"]["initial_state"])
        easy_start = true
    end

    # Allow more steps for first solution
    if i_counter == 1
        max_steps   *= 2
        max_runtime *= 2
    end

    # Solve for RCE
    succ = solver.solve_energy!(atmos, sol_type=cfg["execution"]["solution_type"],
                                            conduct=cfg["physics"]["conduction"],
                                            convect=cfg["physics"]["convection"],
                                            latent=cfg["physics"]["latent_heat"],
                                            sens_heat= cfg["physics"]["sensible_heat"],
                                            chem=cfg["physics"]["chemistry"],
                                            max_steps=max_steps,
                                            max_runtime=max_runtime,
                                            conv_atol= cfg["execution"]["converge_atol"],
                                            conv_rtol=cfg["execution"]["converge_rtol"],
                                            method=1,
                                            rainout=Bool(cfg["physics"]["rainout"]),
                                            oceans=Bool(cfg["physics"]["oceans"]),
                                            dx_max=Float64(cfg["execution"]["dx_max"]),
                                            dx_min=Float64(cfg["execution"]["dx_min"]),
                                            ls_method=Int(cfg["execution"]["linesearch"]),
                                            easy_start=easy_start,
                                            modplot=modplot,
                                            save_frames=false,
                                            radiative_Kzz=false,
                                            perturb_all=cfg["execution"]["perturb_all"],
                                            )
    # Report radius
    @info @sprintf("    found r_phot = %.3f R_earth",atmos.transspec_r/R_earth)

    # Write NetCDF file for this case
    if save_netcdfs
        @info "    write netcdf file"
        save.write_ncdf(atmos, joinpath(output_dir,"nc",@sprintf("%08d.nc",i)))
    end

    # Make plot for this case
    if save_plots
        plotting.combined(plotting.plot_pt(atmos,""), plotting.plot_fluxes(atmos, ""),
                            plotting.plot_vmr(atmos,""), plotting.plot_radius(atmos, ""),
                            "Index = $i     Success = $succ",
                            joinpath(output_dir,"pl",@sprintf("%08d_pl.png",i)))
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
                result_table[i][k] = -1.0 # failure
            end
        elseif k == "flux_loss"
            result_table[i][k] = maximum(abs.(atmos.flux_tot)) - minimum(abs.(atmos.flux_tot))

        elseif k == "r_bound"
            if all(atmos.layer_isbound)
                # fully bound by gravity
                result_table[i][k] = -1
            else
                # unbound at some layer
                result_table[i][k] = minimum(atmos.r[.!atmos.layer_isbound])
            end

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

        elseif startswith(k, "vmr_")
            gas = split(k,"_")[2]
            result_table[i][k] = atmos.gas_vmr[gas][end]

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
    end # end keys  loop

    # Record fluxes
    result_emits[i,:] .= atmos.band_u_lw[1, :] .+ atmos.band_u_sw[1, :]

    # Record profile (also in SI)
    result_profs[i] = Dict( "p"=>deepcopy(atmos.p),
                            "t"=>deepcopy(atmos.tmp),
                            "r"=>deepcopy(atmos.r)
                          )


    # Update results file on the go?
    if mod(i,modwrite) == 0
        write_table()
        write_emits()
        write_profs(atmos.nlev_c)
    end

    # Iterate
    @info "  "
end

@info "Worker (ID=$id_work) completed all allocated points"
@info "------------------------------"
@info " "

@info "Writing final data..."
write_table()
write_emits()
write_profs(atmos.nlev_c)

@info "Deallocating atmosphere"
atmosphere.deallocate!(atmos)

@info "===================================="
@info " "

# Get end time
time_end::Float64 = time()
@info "Finish time: $(now())"

# Print model statistics
duration = (time_end - time_start) / 60 # minutes
if duration > 60
    @info @sprintf("Total runtime: %.2f hours", duration/60)
else
    @info @sprintf("Total runtime: %.2f minutes", duration)
end

# Done
@info "Done!"
