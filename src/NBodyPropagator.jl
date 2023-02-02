module NBodyPropagator
using LinearAlgebra

include("SolarSystemDynamics.jl")
# Write your package code here.

export twobody!

# Define the Solar System Dynamics Constant
ssd = SolarSystemDynamics()

# Definition of Ordinary Differential Equation
function twobody!(dxdt, x, p, t)
    # dx/dt
    r = norm(x[1:3])
    dxdt[1:3] = x[4:6]
    dxdt[4:6] = -ssd.GM["EARTH"] * x[1:3] / (r^3)
end


"""
NBodyProblem{K}(x0, tspan, list_bodies; kwargs...)

Constructor of the Struct `NBodyProblem{K}`.

# Arguments
- `x0`: Initial state vector.
- `tspan`: Time span.
- `list_bodies`: List of bodies.
- `kwargs`: Keyword arguments.

# Keyword Arguments
- `id_center`: ID of the center body.
- `ref_frame`: Reference frame.
- `lsf`: Length scale factor.
- `tsf`: Time scale factor.
- `msf`: Mass scale factor.
- `need_stm`: Whether the STM is needed.
- `need_transitional_state`: Whether the transitional state is needed.
"""
struct NBodyProblem{K}
    # Fields
    x0::Array{Float64,1}
    tspan::Tuple{Float64,Float64}
    list_bodies::Array{Int,1}
    kwargs::K

    # Constructor of the Struct
    function NBodyProblem(x0, tspan, list_bodies; kwargs...)

        # Merge keyword argments
        defaults = (; id_center=0, #= SOLAR SYSTEM BARYCENTER =#
            ref_frame="ECLIPJ2000",
            lsf=1.0e6, #= km =#
            tsf=1.0e6, #= s =#
            msf=1.0, #= kg =#
            need_stm=false,
            need_transitional_state=false
        )
        kwargs_ = merge(defaults, kwargs)

        # Substitute the value to the fields
        new{typeof(kwargs_)}(x0, tspan, list_bodies, kwargs_)
    end
end

end
