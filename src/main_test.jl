#include model.jl file
include("model.jl")
using StaticArrays
using Benchmarktools
# Test the model
u0 = SA[1.0, 0.0, 0.0]  # Initial conditions
tspan = (0.0, 50.0)  # Time span for the simulation
param = []  # Parameters (empty for this model)
prob = ODEProblem(eom, u0, tspan, param)
@btime sol = solve(prob, Tsit5())  # Solve the ODE problem using Tsit5 method

# Visualize the solution
using Plots
plot(sol, vars=(1, 2), xlabel="x", ylabel="y", title="Lorenz Attractor: x vs y")