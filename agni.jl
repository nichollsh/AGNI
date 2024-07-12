#!/usr/bin/env -S julia -O2 -t auto --color=yes --startup-file=no

# Don't show plot windows
ENV["GKSwstype"] = "100"

# Check RAD_DIR 
if !("RAD_DIR" in keys(ENV))
    error("Cannot find SOCRATES! Have you set RAD_DIR?")
end 

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
