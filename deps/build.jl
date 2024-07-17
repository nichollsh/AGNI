# Build AGNI
ROOT_DIR=abspath(joinpath(dirname(abspath(PROGRAM_FILE)),".."))

# Find socrates
RAD_DIR = abspath(ENV["RAD_DIR"])

# Generate wrappers
wrap = joinpath(RAD_DIR, "julia/src/generate_wrappers.jl")
include(wrap)

# Build libSOCRATES
cd(joinpath(RAD_DIR,"julia/lib/")) do 
    run(`make`)
end 

# Download basic data 
run(`bash $(joinpath(ROOT_DIR,get_data.sh)) basic`)

