# Build AGNI

# Find socrates
RAD_DIR = abspath(ENV["RAD_DIR"])

# Generate wrappers
wrap = joinpath(RAD_DIR, "julia/src/generate_wrappers.jl")
include(wrap)

# Build libSOCRATES
cd(joinpath(RAD_DIR,"julia/lib/")) do 
    run(`make`)
end 


