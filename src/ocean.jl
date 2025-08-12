# Contains module handling surface oceans

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

"""
**This module handles ocean formation**

Assumes a binary area of ocean basin (OB) and contintenal shelf (CS). Rained-out liquid is
deposited into these two reservoirs, filling the oceans first. The functions are designed
to handle multiple condensing liquids. The exposed liquid is the one with the lowest
density, assuming that they do not mix miscibly.

At the moment this assumes a top-hat continent shape. But these could be expanded to allow
for a shelf area distribution which changes with height (i.e. sloped shelf edges).
"""
module ocean

    import ..phys

    """
    **Determine the layering structure of surface liquids**

    Output is an array of length=4 tuples, one tuple per liquid. 1st tuple index is the
    location of the liquid (i=1 is bottom of ocean). 2nd index is the liquid name.
    3rd index is its height [m] in the ocean bains. 4th index is height above oceans.

    Arguments:
    - `sigs::Dict{String, Float64}`    dictionary of rainout densities for liquids [kg m-2]
    - `fOB::Float64`                   ocean basin area, as fraction of planet surface
    - `hCS::Float64`                   continental shelf height [m]
    - `Rpl::Float64`                   planet radius at bottom of ocean [m]

    Returns:
    - `output::Array`                  array containing tuples as described above.
    """
    function dist_surf_liq(sigs::Dict{String, Float64},
                            fOB::Float64, hCS::Float64, Rpl::Float64)::Array{Tuple,1}

        # Areas [m2]
        aPL::Float64 = 4 * pi * Rpl^2    # planet area
        aOB::Float64 = aPL * fOB         # ocean basin

        # Total capacity of oceans [m3]
        v_cap_OB::Float64 = aOB * hCS

        # Get liquid densities
        liqs::Array{String, 1} = collect(keys(sigs))
        rhos::Array{Float64,1} = Float64[ phys.liquid_rho(l) for l in liqs]

        # Sort liquids by decreasing density
        mask::Array{Int, 1} = reverse(sortperm(rhos))

        # Loop through liquids, filling oceans first
        v_tot_OB::Float64 = 0.0 # current amount of all things in oceans
        v_fill_OB::Array{Float64,1} = zero(rhos)  # volume of each liquid in oceans
        v_fill_CS::Array{Float64,1} = zero(rhos)  # volume of each liquid on continents
        output::Array = []
        for (j,i) in enumerate(mask)
            # total volume of rain from this liquid
            v_liq = sigs[liqs[i]] * aPL / rhos[i]

            # volume currently contained within oceans
            v_tot_OB = sum(v_fill_OB)

            # oceans full already?
            if v_tot_OB >= v_cap_OB
                # add all liquid to surface
                v_fill_CS[i] = v_liq

            # oceans not full yet
            else
                # oceans *will* be filled by this liquid
                if v_liq > v_cap_OB - v_tot_OB
                    v_fill_OB[i] = v_cap_OB - v_tot_OB    # first into oceans
                    v_fill_CS[i] = v_liq - v_fill_OB[i]   # remainder onto surface

                # oceans can contain all of this liquid, without filling them
                else
                    v_fill_OB[i] = v_liq
                end
            end

            # output array with consolidated info
            #   surface liquid is divided by total area, since it covers oceans+continents
            #   do not include empty layers
            if v_fill_CS[i] + v_fill_OB[i] > 1.0
                push!(output, (j, liqs[i], v_fill_OB[i]/aOB, v_fill_CS[i]/aPL))
            end
        end

        return output
    end # end ocean_layers

    """
    **Get component of topmost ocean layer.**

    This is the liquified component exposed to the atmosphere.
    """
    function get_topliq(layers::Array{Tuple,1})::String
        if length(layers) == 0
            return "none"
        else
            return layers[end][2]
        end
    end

    """
    **Get ocean depth [m] at deepest point of ocean**
    """
    function get_maxdepth(layers::Array{Tuple,1})::Float64
        if length(layers) == 0
            return 0.0
        else
            # sum over all layers
            return sum([la[3]+la[4] for la in layers])
        end
    end

    """
    **Get area-fraction of the planet that is covered by oceans.**

    The area depends on how full the ocean basins are.
    """
    function get_areacov(layers::Array{Tuple,1},
                            fOB::Float64)::Float64

        # basins empty => desert planet
        if get_maxdepth(layers) < 1e-2
            return 0.0

        # basins not empty
        else
            # liquid on shelf => aqua planet
            if sum([la[4] for la in layers]) > 1e-2
                return 1.0

            # otherwise, must be contintental planet
            else
                return fOB
            end
        end
    end
end
