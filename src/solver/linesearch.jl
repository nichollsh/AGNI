# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0

"""
**Module for step-length control (globalisation) of Newton-type solvers.**

These are pure functions which act on user-provided closures. They do not depend on the
atmosphere model, so that they can be tested on analytic problems.
"""
module linesearch

    using LinearAlgebra

    """
    **Backtracking linesearch satisfying the Armijo (sufficient decrease) condition.**

    Finds a step scale `α` along a search direction such that the merit function along
    that direction, `φ(α)`, satisfies the Armijo condition
    ```
    φ(α) <= φ_ref + c1 * α * dφ0
    ```
    where `dφ0 < 0` is the directional derivative at `α=0`. Setting `φ_ref = φ0` gives the
    usual monotone condition. Setting `φ_ref` to the largest of the recent merit values
    gives the non-monotone condition of Grippo, Lampariello & Lucidi (1986).

    Starts from `α0` (normally the full step), then backtracks using a safeguarded
    quadratic (first backtrack) or cubic (subsequent backtracks) interpolation of `φ`,
    restricting each new `α` to lie within `[0.1, 0.5]` times the previous value. Trial
    points returning a non-finite merit are halved. See Nocedal & Wright (2006, Sec. 3.5).

    For a Newton step `d = -J⁻¹ r` and merit `φ(α) = ‖r(x + α d)‖ₚ`, the directional
    derivative is `dφ0 = -φ0`, for any p-norm.

    Arguments:
    - `φ::Function`         merit function along the search direction, `φ(α)::Float64`
    - `φ0::Float64`         merit at `α=0`
    - `dφ0::Float64`        directional derivative of the merit at `α=0`

    Optional arguments:
    - `φ_ref::Float64`      reference merit in the sufficient decrease condition
    - `α0::Float64`         initial step scale
    - `c1::Float64`         sufficient decrease parameter, 0 < c1 < 1
    - `α_min::Float64`      smallest step scale to try before giving up
    - `max_iter::Int64`     maximum number of merit evaluations

    Returns:
    - `α::Float64`          accepted step scale, or the best one found if unsuccessful
    - `φα::Float64`         merit at `α`
    - `nfev::Int64`         number of merit evaluations
    - `ok::Bool`            whether the sufficient decrease condition was satisfied
    """
    function armijo(φ::Function, φ0::Float64, dφ0::Float64;
                        φ_ref::Float64=φ0, α0::Float64=1.0, c1::Float64=1e-4,
                        α_min::Float64=1e-4, max_iter::Int64=30
                        )::Tuple{Float64,Float64,Int64,Bool}

        # Validate inputs
        if !(0.0 < c1 < 1.0)
            throw(ArgumentError("Armijo parameter c1 must be in (0,1), got $c1"))
        end
        if !(α0 > 0.0) || !(α_min > 0.0)
            throw(ArgumentError("Step scales must be positive, got α0=$α0, α_min=$α_min"))
        end

        # Not a descent direction: no step can be guaranteed to decrease the merit
        if !(dφ0 < 0.0) || !isfinite(φ0)
            return (0.0, φ0, 0, false)
        end

        # Best point found so far
        α_best::Float64 = 0.0
        φ_best::Float64 = φ0

        # Previous trial, for cubic interpolation
        α_prv::Float64 = 0.0
        φ_prv::Float64 = φ0

        α::Float64  = α0
        φα::Float64 = Inf
        nfev::Int64 = 0
        while nfev < max_iter

            # Evaluate merit at this trial
            φα = φ(α)
            nfev += 1
            if isfinite(φα) && (φα < φ_best)
                α_best = α
                φ_best = φα
            end

            # Sufficient decrease?
            if isfinite(φα) && (φα <= φ_ref + c1 * α * dφ0)
                return (α, φα, nfev, true)
            end

            # Choose next trial
            if !isfinite(φα)
                # Merit evaluation failed, so halve the step
                α_new = 0.5 * α

            elseif α_prv == 0.0
                # Minimiser of quadratic through φ0, dφ0, and φ(α)
                α_new = -dφ0 * α^2 / (2.0 * (φα - φ0 - dφ0 * α))

            else
                # Minimiser of cubic through φ0, dφ0, φ(α), and φ(α_prv)
                t1 = φα    - φ0 - dφ0 * α
                t2 = φ_prv - φ0 - dφ0 * α_prv
                den = α - α_prv
                a = ( t1 / α^2 - t2 / α_prv^2) / den
                b = (-α_prv * t1 / α^2 + α * t2 / α_prv^2) / den
                if abs(a) < eps(Float64)
                    # Cubic is degenerate, so use the quadratic
                    α_new = -dφ0 / (2.0 * b)
                else
                    disc = b^2 - 3.0 * a * dφ0
                    α_new = (disc < 0.0) ? 0.5 * α : (-b + sqrt(disc)) / (3.0 * a)
                end
            end

            # Safeguard the new step scale
            if !isfinite(α_new)
                α_new = 0.5 * α
            end
            α_new = clamp(α_new, 0.1 * α, 0.5 * α)

            # Record this trial, and move to the next
            if isfinite(φα)
                α_prv = α
                φ_prv = φα
            end
            α = α_new

            # Step scale too small?
            if α < α_min
                break
            end
        end

        # Unsuccessful: return best point found
        return (α_best, φ_best, nfev, false)
    end # end armijo

    """
    **Reference merit for a non-monotone linesearch.**

    Returns the largest of the last `M` merit values, which is the reference value used
    by the non-monotone Armijo condition of Grippo, Lampariello & Lucidi (1986). With
    `M=1` this reduces to the most recent merit value, and the usual monotone condition.

    Arguments:
    - `history::Array{Float64,1}`   merit values at previous iterates, oldest first
    - `M::Int64`                    memory length, M >= 1

    Returns:
    - `Float64`                     reference merit value
    """
    function nonmonotone_ref(history::Array{Float64,1}, M::Int64)::Float64
        if M < 1
            throw(ArgumentError("Non-monotone memory length must be >= 1, got $M"))
        end
        if isempty(history)
            throw(ArgumentError("Merit history is empty"))
        end
        return maximum(@view history[max(1, end-M+1):end])
    end # end nonmonotone_ref

    """
    **Broyden's ("good") rank-one update of an approximate Jacobian.**

    Updates `B` in-place so that it satisfies the secant equation `B s = y`, while making
    the smallest change to `B` in the Frobenius norm (Broyden 1965; Nocedal & Wright 2006,
    Sec. 11.1).

    Arguments:
    - `B::Array{Float64,2}`     approximate Jacobian, updated in-place
    - `s::Array{Float64,1}`     step taken in the solution vector
    - `y::Array{Float64,1}`     corresponding change in the residuals

    Returns:
    - `Bool`                    whether the update was applied (not if `s` is zero or non-finite)
    """
    function broyden_update!(B::Array{Float64,2}, s::Array{Float64,1},
                                y::Array{Float64,1})::Bool
        sts = dot(s, s)
        if !(sts > 0.0) || !isfinite(sts) || !all(isfinite, y)
            return false
        end
        u = (y .- B * s) ./ sts
        B .+= u * s'
        return true
    end # end broyden_update!

    export armijo, nonmonotone_ref, broyden_update!

end
