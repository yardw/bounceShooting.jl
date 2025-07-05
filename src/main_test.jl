#include model.jl file
include("model.jl")
using BenchmarkTools

# Test the Lorentz model
u0 = SA[1.0, 0.0, 0.0]  # Initial conditions
tspan = (0.0, 50.0)  # Time span for the simulation
param = []  # Parameters (empty for this model)
prob = ODEProblem(eom_Lorentz, u0, tspan, param)
@btime sol = solve(prob, Tsit5())  # Solve the ODE problem using Tsit5 method


# Visualize the solution
using Plots
plot(sol, vars=(1, 2), xlabel="x", ylabel="y", title="Lorenz Attractor: x vs y")


# test the bounce problem
params = SA[1.0, 0.1, 1e-1]  # Parameters: ϕ0 = 1.0, λ = 0.1, ϵ = 1e-3 where ϕ0 is the vacuum expectation value (vev), λ is the quartic coupling and ϵ is a small perturbation
u0 = SA[-0.16243456471667703, 0.0]  # Initial conditions u0 = [ϕ, ϕ'], where ϕ is the scalar field and ϕ' is its derivative over r
tspan = (1e-4, 150.0)  # Time span for the simulation. Note that the initial time is set to a small value to avoid singularity at r=0

prob = ODEProblem(eom_bounce, u0, tspan, params)
sol = solve(prob, Tsit5())  # Solve the ODE problem using Tsit5 method
sol = solve(prob, Tsit5(), save_on=false)  # Solve the ODE problem using Tsit5 method
sol = solve(prob, Tsit5(), callback = cbs_bounce)  # Solve the ODE problem using Tsit5 method
let params = SA[1.0, 0.1, 1e-1], ϕ0_scan_range = [-1., 0.], tspan = (1e-4, 200)
    prob = ODEProblem(eom_bounce, SA[ϕ0_scan_range[1], 0.0], tspan, params)
    ϕ0_sections = [ϕ0_scan_range[1], ϕ0_scan_range[2]]
    ϕ0mid = (ϕ0_sections[1] + ϕ0_sections[2]) / 2
    # ϕ0s = [(ϕ0mid,NaN,NaN)]
    for i in 1:100

        sol = solve(remake(prob, u0=SA[ϕ0mid, 0.0]), Tsit5(), callback = cbs_bounce, save_on=false)

        if is_falling_out(sol, params) || is_overshot(sol, params)
            # push!(ϕ0s, (ϕ0mid,1, sol.t[end]))
            # @info "overshot at ϕ0 = $ϕ0mid"
            ϕ0_sections[2] = ϕ0mid
        elseif is_undershot(sol, params)
            # push!(ϕ0s, (ϕ0mid,-1, sol.t[end]))
            # @info "undershot at ϕ0 = $ϕ0mid"
            ϕ0_sections[1] = ϕ0mid
        else
            @info sol.retcode
            return ϕ0_sections
        end
        ϕ0mid = (ϕ0_sections[1] + ϕ0_sections[2]) / 2
    end
    return ϕ0_sections
end

# sol
# if sol.retcode != :Success
#     error("ODE solver failed to find a solution: $(sol.retcode)")
# end

# Visualize the solution
using Plots
plot(sol, idxs=(1, 2), xlabel="ϕ", ylabel="ϕ'", title="Bounce Solution: ϕ vs ϕ'")
plot(sol, vars=(0, 1), xlabel="r", ylabel="ϕ", title="Bounce Solution: ϕ vs r")