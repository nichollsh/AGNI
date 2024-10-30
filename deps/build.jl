# Build AGNI
ROOT_DIR=abspath(joinpath(dirname(abspath(PROGRAM_FILE)),".."))
println("ROOT_DIR = $ROOT_DIR")

# Find socrates
RAD_DIR = abspath(ENV["RAD_DIR"])
println("RAD_DIR = $RAD_DIR")

# Get socrates version
SOC_VER = readchomp(joinpath(RAD_DIR,"version"))
println("SOC_VER = $SOC_VER")

# Generate wrappers
println("Generate wrappers")
wrap = joinpath(RAD_DIR, "julia/src/generate_wrappers.jl")
include(wrap)

# Build libSOCRATES
println("Build libSOCRATES")
cd(joinpath(RAD_DIR,"julia/lib/")) do
    run(`make`)
end

# Download basic data
println("Get data")
get_data = joinpath(ROOT_DIR,"src/get_data.sh")
run(`bash $get_data basic`)

println("Build completed")
