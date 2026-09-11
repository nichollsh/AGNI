#!/usr/bin/env -S julia -O2 --color=yes --startup-file=no

# Configure how plotting windows are displayed
ENV["WAYLAND_DISPLAY"] = ""
ENV["QT_PLUGIN_PATH"] = ""
ENV["QML2_IMPORT_PATH"] = ""
ENV["QT_QPA_PLATFORM"] = "xcb;offscreen"
if get(ENV,"SHOW_PLOTS","0") == "1"
    ENV["GKSwstype"] = "411"
else
    ENV["GKSwstype"] = "100"
end

# Path to AGNI directory
const AGNI_DIR::String = dirname(abspath(@__FILE__))

# Check RAD_DIR
if !("RAD_DIR" in keys(ENV))
    error("Cannot find SOCRATES! Have you set RAD_DIR?")
end

# Check SOCRATES.jl
const SOCjl::String = joinpath(abspath(ENV["RAD_DIR"]),"julia","src","SOCRATES.jl")
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
