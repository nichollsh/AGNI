#!/usr/bin/env -S julia --color=yes --startup-file=no

# Get AGNI root directory
ENV["GKSwstype"] = "100"

include("src/AGNI.jl")
import .AGNI 

if AGNI.main()
    exit(0)
end 
exit(1)
