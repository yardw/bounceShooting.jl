using DifferentialEquations, StaticArrays

function eom(u, p, t)
    x, y, z = u
    dx = 10 * (y - x)
    dy = x * (28 - z) - y
    dz = x * y - (8/3) * z
    return SA[dx, dy, dz]
end