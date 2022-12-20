
module NBodyProblem
using LinearAlgebra

include("SolarSystemDynamics.jl")
# Write your package code here.

export twobody!

# Definition of Ordinary Differential Equation
function twobody!(dxdt, x, p, t)
    # dx/dt
    r = norm(x[1:3])
    dxdt[1:3] = x[4:6]
    dxdt[4:6] = -ssd.GM["EARTH"] * x[1:3] / (r^3)
end

end
