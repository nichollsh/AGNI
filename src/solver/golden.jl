# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0

"""
**Module for performing a golden section search**
"""
module golden

    using Printf
    using LoggingExtras

    """
    **Golden section search algorithm**

    Minimises a function `f` between the bounds `a` and `b`.
    The function `f` must return a positive scalar.

    Arguments:
    - `f`           Function to be minimised
    - `a`           Initial bracket, lower value
    - `b`           Initial bracket, upper valuie
    - `dxtol`       Convergence: exit when bracket is smaller than this size
    - `atol`        Convergence: absolute tolerance on minimum
    - `max_steps`   Maximum number of iterations
    - `warnings`    Print warning if search does not converge before `max_steps`` are taken.

    Returns:
    - `sol`         Best solution found by search.
    - `succ`        Did search converge?
    """
    function gs_search(f::Function,a::Float64,b::Float64,
                            dxtol::Float64,atol::Float64,max_steps::Int64;
                            warnings::Bool=false)::Tuple{Float64,Bool}
        c::Float64 = (-1+sqrt(5))/2

        x1::Float64 = c*a + (1-c)*b
        x2::Float64 = (1-c)*a + c*b

        fx1::Float64 = f(x1)
        fx2::Float64 = f(x2)

        midp::Float64 = 0.5*(a+b)

        for i in 1:max_steps
            if fx1 < fx2
                b = x2
                x2 = x1
                fx2 = fx1
                x1 = c*a + (1-c)*b
                fx1 = f(x1)
            else
                a = x1
                x1 = x2
                fx1 = fx2
                x2 = (1-c)*a + c*b
                fx2 = f(x2)
            end


            if fx1 < atol
                # @debug "GS search found minimum. $(i+2) function evaluations, best = $x1"
                return (x1, true)
            end

            if fx2 < atol
                # @debug "GS search found minimum. $(i+2) function evaluations, best = $x2"
                return (x2, true)
            end

            midp = 0.5*(a+b)

            if abs(b-a) < dxtol
                # @debug "GS search reached minimum bracket size. $(i+2) function evaluations, best = $best"
                return (midp, false)
            end
        end

        if warnings
            @warn "GS search did not converge"
        end
        return (midp, false)
    end
    export gs_search

end
