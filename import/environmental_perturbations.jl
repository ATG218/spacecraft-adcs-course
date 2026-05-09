# Environmental perturbation models aligned with Lecture 13 (drag, SRP, gravity gradient torque)
# and Curtis-style orbital atmosphere conventions (rho0, h_ref, scale height).
# Drag: exponential rho(h); v_rel = v - omega_E x r (rotating atmosphere).
# SRP: pressure flux at 1 AU ~ 4.5e-6 N/m² [solar flux ~1368 W/m²]; scales with reflectivity Cr and AU distance.
# Gravity gradient: tau = (3 mu / r^3) * r_hat_b x (I * r_hat_b) with r_hat in body frame via quaternion DCM.

module EnvironmentalPerturbations

using LinearAlgebra
using StaticArrays

export EarthEnvironment,
       SpacecraftAeroProperties,
       SpacecraftSRPProperties,
       quat_to_dcm_body_to_inertial,
       atmospheric_density_exponential,
       relative_atmosphere_velocity,
       gravity_gradient_torque,
       atmospheric_drag_force,
       atmospheric_drag_torque,
       solar_radiation_pressure_force,
       solar_radiation_pressure_torque,
       environmental_wrench,
       environmental_wrench_from_km_state

const MU_EARTH = 3.986004418e14
const R_EARTH = 6_378_136.3
const OMEGA_EARTH = @SVector [0.0, 0.0, 7.2921159e-5]
const SOLAR_PRESSURE_1AU = 4.56e-6

Base.@kwdef struct EarthEnvironment
    mu::Float64 = MU_EARTH
    radius::Float64 = R_EARTH
    omega_planet::SVector{3,Float64} = OMEGA_EARTH
    rho_ref::Float64 = 3.614e-13
    h_ref::Float64 = 700e3
    H_scale::Float64 = 88.667e3
    solar_pressure_1au::Float64 = SOLAR_PRESSURE_1AU
end

Base.@kwdef struct SpacecraftAeroProperties
    Cd::Float64 = 2.2
    area_m2::Float64 = 0.03
    mass_kg::Float64 = 12.0
    area_normal_body::SVector{3,Float64} = @SVector [1.0, 0.0, 0.0]
    center_of_pressure_body_m::SVector{3,Float64} = @SVector [0.0, 0.0, 0.0]
end

Base.@kwdef struct SpacecraftSRPProperties
    Cr::Float64 = 1.3
    area_m2::Float64 = 0.03
    area_normal_body::SVector{3,Float64} = @SVector [1.0, 0.0, 0.0]
    center_of_pressure_body_m::SVector{3,Float64} = @SVector [0.0, 0.0, 0.0]
    sun_distance_au::Float64 = 1.0
    eclipse::Bool = false
end

function quat_to_dcm_body_to_inertial(q::SVector{4,<:Real})
    qn = SVector{4,Float64}(q ./ norm(q))
    q0, q1, q2, q3 = qn
    return @SMatrix [
        1 - 2 * (q2^2 + q3^2)   2 * (q1*q2 - q0*q3)       2 * (q1*q3 + q0*q2)
        2 * (q1*q2 + q0*q3)     1 - 2 * (q1^2 + q3^2)     2 * (q2*q3 - q0*q1)
        2 * (q1*q3 - q0*q2)     2 * (q2*q3 + q0*q1)       1 - 2 * (q1^2 + q2^2)
    ]
end

function atmospheric_density_exponential(altitude_m::Real, env::EarthEnvironment=EarthEnvironment())
    return env.rho_ref * exp(-(altitude_m - env.h_ref) / env.H_scale)
end

function relative_atmosphere_velocity(r_eci_m::SVector{3,<:Real},
                                      v_eci_mps::SVector{3,<:Real},
    env::EarthEnvironment=EarthEnvironment())
    r = SVector{3,Float64}(r_eci_m)
    v = SVector{3,Float64}(v_eci_mps)
    return v - cross(env.omega_planet, r)
end

function gravity_gradient_torque(I_body::SMatrix{3,3,<:Real},
                                 q_b_to_i::SVector{4,<:Real},
                                 r_eci_m::SVector{3,<:Real};
                                 mu::Real=MU_EARTH)
    Q_bi = quat_to_dcm_body_to_inertial(q_b_to_i)
    r_hat_i = normalize(SVector{3,Float64}(r_eci_m))
    r_hat_b = Q_bi' * r_hat_i
    I_b = SMatrix{3,3,Float64}(I_body)
    return (3.0 * Float64(mu) / norm(r_eci_m)^3) * cross(r_hat_b, I_b * r_hat_b)
end

function atmospheric_drag_force(r_eci_m::SVector{3,<:Real},
                                v_eci_mps::SVector{3,<:Real},
                                q_b_to_i::SVector{4,<:Real},
                                aero::SpacecraftAeroProperties=SpacecraftAeroProperties();
                                env::EarthEnvironment=EarthEnvironment(),
                                use_projected_area::Bool=true)
    r = SVector{3,Float64}(r_eci_m)
    v_rel_i = relative_atmosphere_velocity(r, v_eci_mps, env)
    speed = norm(v_rel_i)
    speed == 0.0 && return @SVector [0.0, 0.0, 0.0]

    altitude = norm(r) - env.radius
    rho = atmospheric_density_exponential(altitude, env)

    area = aero.area_m2
    if use_projected_area
        Q_bi = quat_to_dcm_body_to_inertial(q_b_to_i)
        n_hat_i = Q_bi * normalize(aero.area_normal_body)
        flow_hat_i = -v_rel_i / speed
        area *= max(0.0, dot(n_hat_i, flow_hat_i))
    end

    return -0.5 * rho * aero.Cd * area * speed * v_rel_i
end

function atmospheric_drag_torque(r_eci_m::SVector{3,<:Real},
                                 v_eci_mps::SVector{3,<:Real},
                                 q_b_to_i::SVector{4,<:Real},
                                 aero::SpacecraftAeroProperties=SpacecraftAeroProperties();
                                 env::EarthEnvironment=EarthEnvironment(),
                                 use_projected_area::Bool=true)
    F_i = atmospheric_drag_force(r_eci_m, v_eci_mps, q_b_to_i, aero;
                                 env=env, use_projected_area=use_projected_area)
    Q_bi = quat_to_dcm_body_to_inertial(q_b_to_i)
    F_b = Q_bi' * F_i
    return cross(aero.center_of_pressure_body_m, F_b)
end

function solar_radiation_pressure_force(q_b_to_i::SVector{4,<:Real},
                                        sun_eci_m::SVector{3,<:Real},
                                        srp::SpacecraftSRPProperties=SpacecraftSRPProperties();
                                        env::EarthEnvironment=EarthEnvironment(),
                                        use_projected_area::Bool=true)
    srp.eclipse && return @SVector [0.0, 0.0, 0.0]

    s_i = SVector{3,Float64}(sun_eci_m)
    s_norm = norm(s_i)
    s_norm == 0.0 && return @SVector [0.0, 0.0, 0.0]

    sun_hat_i = s_i / s_norm
    pressure = env.solar_pressure_1au * srp.Cr / srp.sun_distance_au^2

    area = srp.area_m2
    if use_projected_area
        Q_bi = quat_to_dcm_body_to_inertial(q_b_to_i)
        n_hat_i = Q_bi * normalize(srp.area_normal_body)
        area *= max(0.0, dot(n_hat_i, sun_hat_i))
    end

    # sun_eci_m is interpreted as the inertial direction from the spacecraft to the Sun.
    # Radiation pressure therefore pushes away from the Sun, opposite this direction.
    return -pressure * area * sun_hat_i
end

function solar_radiation_pressure_torque(q_b_to_i::SVector{4,<:Real},
                                         sun_eci_m::SVector{3,<:Real},
                                         srp::SpacecraftSRPProperties=SpacecraftSRPProperties();
                                         env::EarthEnvironment=EarthEnvironment(),
                                         use_projected_area::Bool=true)
    F_i = solar_radiation_pressure_force(q_b_to_i, sun_eci_m, srp;
                                         env=env, use_projected_area=use_projected_area)
    Q_bi = quat_to_dcm_body_to_inertial(q_b_to_i)
    F_b = Q_bi' * F_i
    return cross(srp.center_of_pressure_body_m, F_b)
end

function environmental_wrench(I_body::SMatrix{3,3,<:Real},
                              q_b_to_i::SVector{4,<:Real},
                              r_eci_m::SVector{3,<:Real},
                              v_eci_mps::SVector{3,<:Real};
                              sun_eci_m::Union{Nothing,SVector{3,Float64}}=nothing,
                              env::EarthEnvironment=EarthEnvironment(),
                              aero::Union{Nothing,SpacecraftAeroProperties}=nothing,
                              srp::Union{Nothing,SpacecraftSRPProperties}=nothing,
                              include_gravity_gradient::Bool=true,
                              include_drag::Bool=true,
                              include_srp::Bool=true)
    F_i = @SVector [0.0, 0.0, 0.0]
    tau_b = @SVector [0.0, 0.0, 0.0]

    if include_gravity_gradient
        tau_b += gravity_gradient_torque(I_body, q_b_to_i, r_eci_m; mu=env.mu)
    end
    if include_drag && aero !== nothing
        F_i += atmospheric_drag_force(r_eci_m, v_eci_mps, q_b_to_i, aero; env=env)
        tau_b += atmospheric_drag_torque(r_eci_m, v_eci_mps, q_b_to_i, aero; env=env)
    end
    if include_srp && srp !== nothing && sun_eci_m !== nothing
        F_i += solar_radiation_pressure_force(q_b_to_i, sun_eci_m, srp; env=env)
        tau_b += solar_radiation_pressure_torque(q_b_to_i, sun_eci_m, srp; env=env)
    end

    return (; force_eci_N=F_i, torque_body_Nm=tau_b)
end

function environmental_wrench_from_km_state(I_body::SMatrix{3,3,<:Real},
                                            q_b_to_i::SVector{4,<:Real},
                                            r_eci_km::SVector{3,<:Real},
                                            v_eci_kmps::SVector{3,<:Real};
                                            kwargs...)
    r_eci_m = 1000.0 .* SVector{3,Float64}(r_eci_km)
    v_eci_mps = 1000.0 .* SVector{3,Float64}(v_eci_kmps)
    return environmental_wrench(I_body, q_b_to_i, r_eci_m, v_eci_mps; kwargs...)
end

end
