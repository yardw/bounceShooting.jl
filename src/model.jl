using DifferentialEquations, StaticArrays

# eom for Lorenz system
function eom_Lorentz(u, p, t)
    x, y, z = u
    dx = 10 * (y - x)
    dy = x * (28 - z) - y
    dz = x * y - (8/3) * z
    return SA[dx, dy, dz]
end

# eom for standard bounce problem
function eom_bounce(fields, params, r)
    ϕ, δϕ = fields # ϕ is the scalar field, ϕ' is its derivative over r, where r is the Euclidean time assuming O(4) symmetry of the solution
    ϕ0, λ, ϵ = params # ϕ0 is the vev, λ is the quartic coupling
    # the scalar potential is V(ϕ) = λ((ϕ+ϕ0)^2 - ϕ0^2)^2 + ϵ ϕ
    dϕ = δϕ # the derivative of ϕ with respect to r
    dδϕ = 4λ*(ϕ+ϕ0)*((ϕ+ϕ0)^2 - ϕ0^2) + ϵ - 3δϕ/r # the bounce equation
    # the factor of 3/r comes from the radial coordinate in spherical coordinates
    return δfields = SA[dϕ, dδϕ] 
end
# callback when the solution reaches the boundary condition: 1) ϕ > 0, 2) ϕ' < 0, 3) ϕ < -2ϕ0
@inline function is_falling_out(fields, r, integrator) #terminate the integration
    ϕ, δϕ = fields
    vev, λ, ϵ = integrator.p
    return  ϕ > 0.2vev && δϕ > 0
end
function is_falling_out(sol, p)
    ϕ, δϕ = sol.u[end]
    vev, λ, ϵ = p
    return ϕ > 0.2vev && δϕ > 0
end
@inline function is_overshot(fields, r, integrator)
    ϕ, δϕ = fields
    vev, λ, ϵ = integrator.p
    return ϕ < -2.2vev && δϕ < 0
end
function is_overshot(sol, p)
    ϕ, δϕ = sol.u[end]
    vev, λ, ϵ = p
    return ϕ < -2.2vev && δϕ < 0
end
@inline function is_undershot(fields, r, integrator)
    ϕ, δϕ = fields
    vev, λ, ϵ = integrator.p
    return -1.8vev < ϕ < -0.2vev && δϕ > 0
end
function is_undershot(sol, p)
    ϕ, δϕ = sol.u[end]
    vev, λ, ϵ = p
    return -1.8vev < ϕ < -0.2vev && δϕ > 0
end
# Callback functions for the bounce problem
cbs_bounce = CallbackSet(
    DiscreteCallback(is_falling_out, terminate!),
    DiscreteCallback(is_overshot, terminate!),
    DiscreteCallback(is_undershot, terminate!)
)

"""
Function to parse parameters for the bounce problem
# structure of params:
    - vev: the vacuum expectation value (vev) of the scalar field
    - λ: the quartic coupling constant
    - ϵ: the linear term in the scalar potential
    - ϕ0: the initial value of the scalar field at r0
    - dϕ0: the initial derivative of the scalar field at r0
    - r0: the initial value of the radial coordinate
    - r1: the final value of the radial coordinate
"""
function param_parser_bounce(params)
    p = SA[params.vev, params.λ, params.ϵ]  # Convert parameters to a StaticArray
    u0 = SA[params.ϕ0, params.dϕ0]  # Initial conditions u0 = [ϕ, ϕ'], where ϕ is the scalar field and ϕ' is its derivative over r
    tspan = (params.r0, params.r1)  # Time span for the simulation
    return u0, tspan, p
end

function shoot(inputs, ϕ0_scan_range, n=50)
    # check ϕ0_scan_range is a vector of length 2 and of type Number but not NaN
    if length(ϕ0_scan_range) != 2 || !all(isa.(ϕ0_scan_range, Number)) || any(isnan.(ϕ0_scan_range))
        error("please check ϕ0_scan_range = [lowerbound, upperbound]")
        return 
    end
    u0, tspan, params = param_parser_bounce(inputs)
    prob = ODEProblem(eom_bounce, u0, tspan, params)
    ϕ0_sections = sort(ϕ0_scan_range)
    ϕ0mid = sum(ϕ0_sections) / length(ϕ0_sections)
    for i in 1:n
        sol = solve(remake(prob, u0=SA[ϕ0mid, 0.0]), Tsit5(), callback = cbs_bounce, save_on=false)

        if is_falling_out(sol, params) || is_overshot(sol, params)
            ϕ0_sections[2] = ϕ0mid
        elseif is_undershot(sol, params)
            ϕ0_sections[1] = ϕ0mid
        else
            @info sol.retcode
            break
        end
        ϕ0mid = sum(ϕ0_sections) / length(ϕ0_sections)
    end
    return ϕ0mid, ϕ0_sections[2] - ϕ0mid
end