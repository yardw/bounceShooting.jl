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
    ϕ0, λ, ϵ = integrator.p
    return  ϕ > 0.2ϕ && δϕ > 0
end
function is_falling_out(sol, p)
    ϕ, δϕ = sol.u[end]
    ϕ0, λ, ϵ = p
    return ϕ > 0.2ϕ && δϕ > 0
end
@inline function is_overshot(fields, r, integrator)
    ϕ, δϕ = fields
    ϕ0, λ, ϵ = integrator.p
    return ϕ < -2.2ϕ0 && δϕ < 0
end
function is_overshot(sol, p)
    ϕ, δϕ = sol.u[end]
    ϕ0, λ, ϵ = p
    return ϕ < -2.2ϕ0 && δϕ < 0
end
@inline function is_undershot(fields, r, integrator)
    ϕ, δϕ = fields
    ϕ0, λ, ϵ = integrator.p
    return -1.8ϕ0 < ϕ < -0.2ϕ0 && δϕ > 0
end
function is_undershot(sol, p)
    ϕ, δϕ = sol.u[end]
    ϕ0, λ, ϵ = p
    return -1.8ϕ0 < ϕ < -0.2ϕ0 && δϕ > 0
end
# Callback functions for the bounce problem
cbs_bounce = CallbackSet(
    DiscreteCallback(is_falling_out, terminate!),
    DiscreteCallback(is_overshot, terminate!),
    DiscreteCallback(is_undershot, terminate!)
)

function param_parser_bounce(params)
    p = SA[params.ϕ0, params.λ, params.ϵ]  # Convert parameters to a StaticArray
    u0 = SA[params.ϕ0, params.dϕ0]  # Initial conditions u0 = [ϕ, ϕ'], where ϕ is the scalar field and ϕ' is its derivative over r
    tspan = (params.r0, params.r1)  # Time span for the simulation
    return u0, tspan, p
end