module NBodyPropagator
using LinearAlgebra
using DifferentialEquations

include("SolarSystemDynamics.jl")
# Write your package code here.

export NBodyProblem
export propagate

"""
Initialization function of the NBodyPropagator module, which is called immediately after the module is loaded (e.g., by import, using, or require) at runtime for the first time 
    (i.e., __init__ is only called once, and only after all statements in the module have been executed).
"""
function __init__()
    # Furnish SPICE Kernels
    init_spice_kernels()
end

"""
NBodyProblem{K}(x0, tspan, list_bodies; kwargs...)

Constructor of the Struct `NBodyProblem{K}`.

# Arguments
- `x0`: Initial state vector.
- `tspan`: Time span.
- `list_bodies`: List of bodies.
- `ssd`: Solar system dynamics parameter sets.
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
    if p.need_stm
        dfdr = zeros(3, 3)
        dfdt0 = zeros(3)
    end

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

        # For STM pre-computation
        if p.need_stm
            dfdr_id = gm_scaled * (-Matrix{Float64}(I, 3, 3) / norm(r_rel)^3 + 3.0 * kron(r_rel, r_rel') / norm(r_rel)^5)
            dfdr += dfdr_id
            dfdt0 -= dfdr_id * v_body
        end

        # Acceleration of the central body relative to an inertial frame
        if p.id_center != 0 && p.id_center != id
            dxdt[4:6] -= gm_scaled * r_body / norm(r_body)^3 # TODO: This acceleration should be calculated via ephemeris

            # For STM pre-computation
            if p.need_stm
                dfdr_id = gm_scaled * (-Matrix{Float64}(I, 3, 3) / norm(r_body)^3 + 3.0 * kron(r_body, r_body') / norm(r_body)^5)
                dfdt0 += dfdr_id * v_body
            end
        end
    end

    # For STM computation
    if p.need_stm
        ## Calculate dxdx0
        # Preprocess
        dxdx0 = reshape(x[7:42], 6, 6)
        # Calculate df/dx and d(dxd0)/dt
        dfdx = [
            zeros(3, 3) Matrix{Float64}(I, 3, 3)
            dfdr zeros(3, 3)
        ]
        d_dxdx0_dt = dfdx * dxdx0
        # Postprocess
        dxdt[7:42] = vec(d_dxdx0_dt)

        ## Calcualte dxdt0
        dxdt0 = x[43:48]
        dxdt[43:48] = dfdx * dxdt0
        dxdt[46:48] += dfdt0
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

    # Preprocess (Scaling, initial STM)
    x0_ = [nbp.x0[1:3] / nbp.kwargs.lsf; nbp.x0[4:6] / (nbp.kwargs.lsf / nbp.kwargs.tsf)]
    tspan_ = (nbp.tspan[1] / nbp.kwargs.tsf, nbp.tspan[2] / nbp.kwargs.tsf)
    if nbp.kwargs.need_stm
        append!(x0_, [vec(Matrix{Float64}(I, 6, 6)); zeros(6, 1)])
    end

    # Definition of ODEProblem
    prob = ODEProblem(eom_nbp!, x0_, tspan_, merge((; list_bodies=nbp.list_bodies, ssd=nbp.ssd), nbp.kwargs))
    sol_ = solve(prob, Vern7(); kwargs_...)

    # Postprocess (Unscaling, reshape STMs)
    state_all = vcat(sol_[1:3, :] .* nbp.kwargs.lsf, sol_[4:6, :] .* (nbp.kwargs.lsf / nbp.kwargs.tsf))
    if nbp.kwargs.need_stm
        stm_all = reshape(sol_[7:42, :], 6, 6, :)
        dxdt0_all = sol_[43:48, :]
        # Unscaling
        # stm_all[1:3,1:3,:] *= nbp.kwargs.lsf/nbp.kwargs.lsf #= dr/dr =#
        stm_all[1:3, 4:6, :] *= nbp.kwargs.tsf #= dr/dv =#
        stm_all[4:6, 1:3, :] /= nbp.kwargs.tsf #= dv/dr =#
        # stm_all[4:6,4:6,:] *= (nbp.kwargs.lsf/nbp.kwargs.tsf) / (nbp.kwargs.lsf/nbp.kwargs.tsf) #= dv/dv =#
        dxdt0_all[1:3, :] *= nbp.kwargs.lsf / nbp.kwargs.tsf #= dr/dt0 =#
        dxdt0_all[4:6, :] *= nbp.kwargs.lsf / nbp.kwargs.tsf^2 #= dv/dt0 =#
    end

    # Returns
    if nbp.kwargs.need_stm
        return state_all, stm_all, dxdt0_all
    else
        return state_all
    end
end

end
