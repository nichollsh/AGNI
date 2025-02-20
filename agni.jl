#!/usr/bin/env -S julia -O2 --color=yes --startup-file=no

# Don't show plot windows
ENV["GKSwstype"] = "100"
AGNI_DIR = dirname(abspath(@__FILE__))

# Check RAD_DIR
if !("RAD_DIR" in keys(ENV))
    error("Cannot find SOCRATES! Have you set RAD_DIR?")
end

# Check SOCRATES.jl
SOCjl = joinpath(abspath(ENV["RAD_DIR"]),"julia","src","SOCRATES.jl")
if !isfile(SOCjl)
    error("Cannot find SOCRATES library! Tried: '$SOCjl'")
end

# Activate environment
import Pkg
Pkg.activate(AGNI_DIR)

# Include AGNI
import AGNI

# Run
if AGNI.main()
    exit(0)
end
exit(1)
