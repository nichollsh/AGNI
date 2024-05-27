#!/usr/bin/env -S julia --color=yes --startup-file=no

# Don't show plot windows
ENV["GKSwstype"] = "100"

# Activate environment
import Pkg
Pkg.activate(dirname(abspath(@__FILE__)))

# Include AGNI
include("src/AGNI.jl")
import .AGNI 

# Run
if AGNI.main()
    exit(0)
end 
exit(1)
