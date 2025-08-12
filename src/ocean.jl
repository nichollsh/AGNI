# Contains module handling surface oceans

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module ocean

    import ..phys

    # This module handles ocean formation given a binary distribution of ocean basin (OB)
    #    and contintenal shelf (CS). Rained-out liquid is deposited into these two
    #    reservoirs, filling the oceans first. The functions are designed to handle
    #    multiple condensing liquids. The exposed liquid is the one with the lowest
    #    density, assuming that they do not mix miscibly.

    """
    **Determine the layering structure of surface liquids**

    Rained-out liquid is deposited into two reservoirs. Ocean bains are filled first.
    Once the OBs are full, remaining liquid spills onto the surface and covers the entire
    planet area. Can handle multiple liquids, assuming that they don't mix.

    Output is an array of length=4 tuples, one tuple per liquid. First tuple index is the
    location of the liquid (denser )

    Arguments:
    - `sigs::Dict{String, Float64}`    dictionary of rainout densities for liquids [kg m-2]
    - `fOB::Float64`                   ocean basin area, as fraction of planet surface
    - `hCS::Float64`                   continental shelf height [m]
    - `Rpl::Float64`                   planet radius at bottom of ocean [m]

    Returns:
    - `output::Array`                   Array containing tuples as described above.
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
        rhos::Array{Float64,1} = Float64[ liquid_rho(l) for l in liqs]

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
            push!(output, (j, liqs[i], v_fill_OB[i]/aOB, v_fill_CS[i]/aPL))
        end

        return output
    end # end ocean_layers

end
