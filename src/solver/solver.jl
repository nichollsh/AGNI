# This file is part of AGNI. License is GPL-3.0: https://www.gnu.org/licenses

"""
**Wrapper module around the various solver modules.**
"""
module solver

    using Printf
    using LoggingExtras
    using Statistics: median, mean
    using LinearAlgebra

    import ..atmosphere
    import ..diagnostics
    import ..setpt
    import ..energy
    import ..consts
    import ..phys
    import ..plotting
    import ..chemistry
    import ..ocean
    import ..multicol

    include("golden.jl"); using .golden
    include("energy.jl"); import .solve_energy: solve_energy!
    include("prescribed.jl"); import .solve_prescribed: solve_prescribed!
    include("transparent.jl"); import .solve_transparent: solve_transparent!
    include("globe.jl"); import .solve_globe: solve_globe!

end
