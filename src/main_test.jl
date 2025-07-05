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
# params = SA[1.0, 0.1, 1e-3]  # Parameters: ϕ0 = 1.0, λ = 0.1, ϵ = 1e-3 where ϕ0 is the vacuum expectation value (vev), λ is the quartic coupling and ϵ is a small perturbation
# u0 = SA[1.0-1e-3, 0.0]  # Initial conditions u0 = [ϕ, ϕ'], where ϕ is the scalar field and ϕ' is its derivative over r
# tspan = (1e-3, 100.0)  # Time span for the simulation. Note that the initial time is set to a small value to avoid singularity at r=0

# prob = ODEProblem(eom_bounce, u0, tspan, params)
# sol = solve(prob, Tsit5(), callback = cbs_bounce);  # Solve the ODE problem using Tsit5 method
# sol.retcode

@btime let params = SA[1.0, 0.1, 1e-3], ϕ0_scan_range = [-1e-2, -1e-3], tspan = (1e-4, 100.0)
    # Scan for the bounce solution by varying ϕ0 in the range ϕ0_scan_range
    prob = ODEProblem(eom_bounce, SA[ϕ0_scan_range[1], 0.0], tspan, params) # initialize the ODE problem with the first value of ϕ0_scan_range
    get_t_end(prob, ϕ0) = solve(remake(prob, u0=SA[ϕ0, 0.0]), Tsit5(), callback = cbs_bounce).t[end] # function to get the end time of the solution
    function gradient_descent(xmin, xmax, f, f_goal; α=0.1, max_iter=10) 
        # Gradient descent to find the minimum of f in the range [xmin, xmax]
        x0 = xmax
        x1 = (xmin + xmax) / 2 
        f0 = f(x0)  
        f1 = f(x1) 
        xs = [(x0, f0), (x1, f1)]  # Store the initial points
        iter = 0
        status = "Unfinished"
        while iter < max_iter
            df = (f1 - f0)  # Numerical derivative
            x_new = x1 - α * df/f1 * (x1-x0)  # Update x using the gradient
            if x_new < xmin || x_new > xmax  # Ensure x stays within bounds
                @error "x_new out of bounds: $x_new"
                break
            end
            if f1>f_goal
                status = "Goal Reached"
                break
            end
            x0 = x1  
            x1 = x_new  
            f0 = f1  
            f1 = f(x1)  
            push!(xs, (x1,f1))  # Store the current x
            iter += 1
        end
        return (x=x1, f=f1, iter=iter, status=status, records = xs)  # Return the optimal x and its function value
    end
    xmin, xmax = sort(ϕ0_scan_range)
    return gradient_descent(xmin, xmax, x->get_t_end(prob, x), 20)
end

# sol
# if sol.retcode != :Success
#     error("ODE solver failed to find a solution: $(sol.retcode)")
# end

# Visualize the solution
using Plots
plot(sol, idxs=(1, 2), xlabel="ϕ", ylabel="ϕ'", title="Bounce Solution: ϕ vs ϕ'")
plot(sol, vars=(0, 1), xlabel="r", ylabel="ϕ", title="Bounce Solution: ϕ vs r")