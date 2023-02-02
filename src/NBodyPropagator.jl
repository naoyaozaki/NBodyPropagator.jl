module NBodyPropagator
using LinearAlgebra
using DifferentialEquations

include("SolarSystemDynamics.jl")
# Write your package code here.

export NBodyProblem
export propagate

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
"""
struct NBodyProblem{K}
    # Fields
    x0::Array{Float64,1}
    tspan::Tuple{Float64,Float64}
    list_bodies::Array{Int,1}
    ssd::SolarSystemDynamics
    kwargs::K

    # Constructor of the Struct
    function NBodyProblem(x0, tspan, list_bodies; kwargs...)
        # Define the Solar System Dynamics Constant
        ssd_ = SolarSystemDynamics()

        # Merge keyword argments
        defaults = (;
            id_center=0, #= SOLAR SYSTEM BARYCENTER =#
            ref_frame="ECLIPJ2000",
            lsf=1.0e6, #= km =#
            tsf=1.0e6, #= s =#
            msf=1.0, #= kg =#
            need_stm=false
        )
        kwargs_ = merge(defaults, kwargs)

        # Substitute the value to the fields
        new{typeof(kwargs_)}(x0, tspan, list_bodies, ssd_, kwargs_)
    end
end

"""
   eom_nbp!(dxdt, x, p, t)

Computes the equation of motion for the N-body problem.

# Arguments
- `dxdt`: The derivative of the state vector.
- `x`: The state vector.
- `p`: The parameter sets.
- `t`: The time.
"""
function eom_nbp!(dxdt, x, p, t)
    # Initialize dx/dt
    dxdt[1:3] = x[4:6]
    dxdt[4:6] = zeros(3)

    # Calculate Gravitational Acceleration
    for id in p.list_bodies
        # Get planetary ephemeris
        x_body, _ = SPICE.spkez(id, t * p.tsf, p.ref_frame, "NONE", p.id_center)
        r_body = x_body[1:3] / p.lsf
        v_body = x_body[4:6] / (p.lsf / p.tsf)
        gm_scaled = p.ssd.GM[p.ssd.NAME[id]] / (p.lsf^3 / p.tsf^2)

        # Acceleration due to the body id
        r_rel = x[1:3] - r_body
        dxdt[4:6] -= gm_scaled * r_rel / norm(r_rel)^3

        #     // For STM Pre-computation
        #     if (x.size() > dim_state_) // if (std::get<bool>(options_["need_stm"]))
        #     {
        #         dfdr_id_temp = gm_id_scaled * (-RowMatrixXd::Identity(3, 3) / std::pow(r_rel.norm(), 3) + 3.0 * r_rel * r_rel.transpose() / std::pow(r_rel.norm(), 5));
        #         dfdr += dfdr_id_temp;
        #         dfdt0 -= dfdr_id_temp * v_body;
        #     }

        # Acceleration of the central body relative to an inertial frame
        if p.id_center != 0 && p.id_center != id
            dxdt[4:6] -= gm_scaled * r_body / norm(r_body)^3 # TODO: This acceleration should be calculated via ephemeris

            #         // For STM Pre-computation
            #         if (x.size() > dim_state_) // if (std::get<bool>(options_["need_stm"]))
            #         {
            #             dfdr_id_temp = gm_id_scaled * (-RowMatrixXd::Identity(3, 3) / std::pow(r_body.norm(), 3) + 3.0 * r_body * r_body.transpose() / std::pow(r_body.norm(), 5));
            #             dfdt0 += dfdr_id_temp * v_body;
            #         }
        end
    end
end


"""
   propagate(prob, kwargs...)

Computes the equation of motion for the N-body problem.

# Arguments
- `nbp`: NBodyProblem definition
- `kwargs`: Keyword arguments.
"""
function propagate(nbp::NBodyProblem; kwargs...)
    # Merge keyword argments
    defaults = (;
        reltol=1.0e-10, # Relative Tolerance
        abstol=1.0e-10 # Absolute Tolerance
    )
    kwargs_ = merge(defaults, kwargs)

    # Scaling
    x0_ = [nbp.x0[1:3] / nbp.kwargs.lsf; nbp.x0[4:6] / (nbp.kwargs.lsf / nbp.kwargs.tsf)]
    tspan_ = (nbp.tspan[1] / nbp.kwargs.tsf, nbp.tspan[2] / nbp.kwargs.tsf)

    # Definition of ODEProblem
    prob = ODEProblem(eom_nbp!, x0_, tspan_, merge((; list_bodies=nbp.list_bodies, ssd=nbp.ssd), nbp.kwargs))
    sol_ = solve(prob, Vern7(); kwargs_...)

    # Unscaling
    sol = vcat(sol_[1:3, :] .* nbp.kwargs.lsf, sol_[4:6, :] .* (nbp.kwargs.lsf / nbp.kwargs.tsf))

    # Returns
    return sol
end



# Definition of Ordinary Differential Equation
function twobody!(dxdt, x, p, t)
    # dx/dt
    r = norm(x[1:3])
    dxdt[1:3] = x[4:6]
    dxdt[4:6] = -ssd.GM["EARTH"] * x[1:3] / (r^3)
end

end
