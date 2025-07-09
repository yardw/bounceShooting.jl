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
# let inputs = (r0 = 1e-4, r1 = 200.0, ϕ0 = 1.0, dϕ0 = 0.0, λ = 0.1, ϵ = 1e-1), ϕ0_scan_range = [0.0, -1.0]
#     # check ϕ0_scan_range is a vector of length 2 and of type Number but not NaN
#     if length(ϕ0_scan_range) != 2 || !all(isa.(ϕ0_scan_range, Number)) || any(isnan.(ϕ0_scan_range))
#         error("please check ϕ0_scan_range = [lowerbound, upperbound]")
#         return 
#     end
#     u0, tspan, params = param_parser_bounce(inputs)
#     prob = ODEProblem(eom_bounce, u0, tspan, params)
#     ϕ0_sections = sort(ϕ0_scan_range)
#     ϕ0mid = sum(ϕ0_sections) / length(ϕ0_sections)
#     for i in 1:50

#         sol = solve(remake(prob, u0=SA[ϕ0mid, 0.0]), Tsit5(), callback = cbs_bounce, save_on=false)

#         if is_falling_out(sol, params) || is_overshot(sol, params)
#             ϕ0_sections[2] = ϕ0mid
#         elseif is_undershot(sol, params)
#             ϕ0_sections[1] = ϕ0mid
#         else
#             @info sol.retcode
#             break
#         end
#         ϕ0mid = sum(ϕ0_sections) / length(ϕ0_sections)
#     end
#     return ϕ0mid, ϕ0_sections .- ϕ0mid
# end
using Plots
let inputs = (r0 = 1e-4, r1 = 200.0, ϕ0 = 1.0, dϕ0 = 0.0, λ = 0.1, ϵ = 1e-1, vev=1.0), ϕ0_scan_range = [-0.1, -1.0]
    mean, var = shoot(inputs, ϕ0_scan_range)
    @info "Mean ϕ0: $mean, Variance: $var"
    # Update inputs with the mean ϕ0
    refined_inputs = (inputs..., ϕ0 = mean+var)
    plot(solve(ODEProblem(eom_bounce, param_parser_bounce(refined_inputs)...)), idxs=(0, 1), xlabel="r", ylabel="ϕ", title="Bounce Solution: ϕ vs r", ylims=(-1.5,0.5), label = "overshot"
    )
    refined_inputs = (inputs..., ϕ0 = mean-var)
    plot!(solve(ODEProblem(eom_bounce, param_parser_bounce(refined_inputs)...)), idxs=(0, 1), label="undershot")
end

# test whether is_falling_out is working
let inputs = (r0 = 1e-4, r1 = 200.0, ϕ0 = 1.0, dϕ0 = 0.0, λ = 0.1, ϵ = 1e-1, vev=1.0), ϕ0_scan_range = [-0.1, -1.0]
    mean, var = shoot(inputs, ϕ0_scan_range)
    @info "Mean ϕ0: $mean, Variance: $var"
    # Update inputs with the mean ϕ0
    refined_inputs = (inputs..., ϕ0 = mean-var)
    sol = solve(ODEProblem(eom_bounce, param_parser_bounce(refined_inputs)...), Tsit5(), callback = DiscreteCallback(is_falling_out, terminate!), save_on=false)
    @info "Final state: $(sol.u[end]), retcode: $(sol.retcode)"
    plot(sol, idxs=(0, 1), xlabel="r", ylabel="ϕ", title="Bounce Solution: ϕ vs r", ylims=(-20.5,1.5))
end


# Visualize the solution
using Plots
plot(sol, idxs=(1, 2), xlabel="ϕ", ylabel="ϕ'", title="Bounce Solution: ϕ vs ϕ'")
plot(sol, idxs=(0, 1), xlabel="r", ylabel="ϕ", title="Bounce Solution: ϕ vs r")

