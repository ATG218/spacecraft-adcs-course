# Coupled translational (two-body + J2 + optional env force) + gyrostat attitude
# with optional gravity-gradient torque and LEO drag (force + torque) from
# EnvironmentalPerturbations. SRP omitted by default (LEO + drag choice).
#
# Load after:
#   include("environmental_perturbations.jl")
#   using .EnvironmentalPerturbations

using LinearAlgebra
using StaticArrays
using DifferentialEquations

Base.@kwdef struct OrbitGyrostatEnvParams
    I::SMatrix{3,3,Float64}
    H_r::Float64
    e::SVector{3,Float64}
    env::EarthEnvironment
    aero::SpacecraftAeroProperties
    include_gg::Bool
    include_drag::Bool
end

function orbit_gyrostat_env_ode!(dx, x, p::OrbitGyrostatEnvParams, t)
    r = SVector(x[1], x[2], x[3])
    v = SVector(x[4], x[5], x[6])
    q = SVector(x[7], x[8], x[9], x[10])
    ω = SVector(x[11], x[12], x[13])

    mu = 398600.4418
    RE = 6378.1363
    J2 = 1.08262668e-3
    rnorm = norm(r)
    a_tb = -mu / rnorm^3 * r
    z2 = r[3]^2 / rnorm^2
    a_J2 = (3 * J2 * mu * RE^2 / (2 * rnorm^5)) .* SVector(
        r[1] * (5 * z2 - 1), r[2] * (5 * z2 - 1), r[3] * (5 * z2 - 3))

    w = environmental_wrench_from_km_state(
        p.I, q, r, v;
        sun_eci_m=nothing,
        env=p.env,
        aero=p.aero,
        srp=nothing,
        include_gravity_gradient=p.include_gg,
        include_drag=p.include_drag,
        include_srp=false,
    )
    F_i = SVector{3,Float64}(w.force_eci_N)
    τ_b = SVector{3,Float64}(w.torque_body_Nm)
    m = p.aero.mass_kg
    a_env_km_s2 = (F_i ./ m) ./ 1000.0  # N -> m/s^2 -> km/s^2

    Ω = @SMatrix [
        0.0   -ω[1]  -ω[2]  -ω[3]
        ω[1]   0.0    ω[3]  -ω[2]
        ω[2]  -ω[3]   0.0    ω[1]
        ω[3]   ω[2]  -ω[1]   0.0
    ]
    dq = 0.5 .* (Ω * q)
    h_r = p.H_r .* p.e
    ω_dot = p.I \ (τ_b - cross(ω, p.I * ω + h_r))

    dx[1] = v[1]
    dx[2] = v[2]
    dx[3] = v[3]
    dv = a_tb + a_J2 + a_env_km_s2
    dx[4] = dv[1]
    dx[5] = dv[2]
    dx[6] = dv[3]
    dx[7] = dq[1]
    dx[8] = dq[2]
    dx[9] = dq[3]
    dx[10] = dq[4]
    dx[11] = ω_dot[1]
    dx[12] = ω_dot[2]
    dx[13] = ω_dot[3]
    return nothing
end

function integrate_orbit_gyrostat_env!(
    x0::AbstractVector{Float64},
    p::OrbitGyrostatEnvParams;
    tspan::Tuple{Float64,Float64},
    reltol::Float64=1e-9,
    abstol::Float64=1e-12,
    saveat::Union{Nothing,Real}=30.0,
)
    prob = ODEProblem(orbit_gyrostat_env_ode!, Vector(x0), tspan, p)
    sol = isnothing(saveat) ? solve(prob, Vern9(); reltol=reltol, abstol=abstol) :
          solve(
              prob,
              Vern9();
              reltol=reltol,
              abstol=abstol,
              save_everystep=false,
              save_start=true,
              save_end=true,
              saveat=saveat,
          )
    for u in sol.u
        u[7:10] ./= norm(u[7:10])
    end
    return sol
end
