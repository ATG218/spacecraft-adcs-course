using LinearAlgebra
using StaticArrays
using DifferentialEquations
using Rotations
using Random
using PythonPlot

include("safe_mode_gyrostat.jl")

# Quaternion rotating unit vector a onto unit vector b
function rotation_quat_between(a::SVector{3,Float64}, b::SVector{3,Float64})
    ua = normalize(a); ub = normalize(b)
    c = dot(ua, ub)
    if c >= 0.999999
        return SVector{4,Float64}(1.0, 0.0, 0.0, 0.0)
    elseif c <= -0.999999
        v = abs(ua[1]) < 0.9 ? SVector(1.0,0.0,0.0) : abs(ua[2]) < 0.9 ? SVector(0.0,1.0,0.0) : SVector(0.0,0.0,1.0)
        axis = normalize(cross(ua, v))
        return SVector{4,Float64}(0.0, axis[1], axis[2], axis[3])
    else
        axis = normalize(cross(ua, ub))
        angle = acos(clamp(c, -1.0, 1.0))
        qv = axis * sin(angle/2)
        return SVector{4,Float64}(cos(angle/2), qv[1], qv[2], qv[3])
    end
end

# Classical orbital elements -> ECI position and velocity (km, km/s)
function coe_to_eci(a_km, e, inc, raan, argp, nu)
    mu = 398600.4418
    p = a_km * (1 - e^2)
    rn = p / (1 + e*cos(nu))
    r_pqw = rn .* SVector{3}(cos(nu), sin(nu), 0.0)
    v_pqw = sqrt(mu/p) .* SVector{3}(-sin(nu), e+cos(nu), 0.0)

    cr, sr = cos(raan), sin(raan)
    ci, si = cos(inc),  sin(inc)
    cp, sp = cos(argp), sin(argp)
    R = @SMatrix [
        cr*cp - sr*ci*sp   -cr*sp - sr*ci*cp    sr*si
        sr*cp + cr*ci*sp   -sr*sp + cr*ci*cp   -cr*si
        si*sp               si*cp               ci
    ]
    return R*r_pqw, R*v_pqw
end

# SSO inclination at given altitude (km)
function compute_sso_inclination(altitude_km)
    mu = 398600.4418; RE = 6378.1363; J2 = 1.08262668e-3
    omega_dot = 0.985647 * (π/180) / 86400
    a = RE + altitude_km
    n = sqrt(mu / a^3)
    return acos(clamp(omega_dot / (-1.5 * J2 * n * (RE/a)^2), -1.0, 1.0))
end

function main()

    I_body = @SMatrix [
        0.249  0.0    0.0
        0.0    0.187  0.02
        0.0    0.02   0.148
    ]
    e_panel   = @SVector [0.0, 1.0, 0.0]
    omega_des = rpm2rad(10.0)
    w_des     = omega_des .* e_panel

    I_pert = perturb_inertia(I_body; sig_eig=0.03, sig_rot=deg2rad(2.0))
    H_r, I_rotor_eff = required_rotor_momentum(I_pert, e_panel; ratio=1.2, omega=omega_des)
    println("H_r = $(round(H_r, digits=4)) N·m·s,  I_rotor = $(round(I_rotor_eff, digits=4)) kg·m²")

    # Attitude-only simulation
    q0_att    = @SVector [1.0, 0.0, 0.0, 0.0]
    delta_w   = 0.01 * omega_des
    xi        = randn(3); xi -= dot(xi, e_panel)*e_panel; xi /= norm(xi)
    w0_att    = w_des + delta_w .* xi

    sol_att = integrate_gyrostat!(I_pert, w0_att, q0_att, H_r, e_panel;
        tspan=(0.0, 800.0), reltol=1e-9, abstol=1e-12)

    t = sol_att.t
    fig, ax = subplots(3, 1, figsize=(8, 6), sharex=true)
    ax[0].plot(t, [u[5] for u in sol_att.u], color="steelblue")
    ax[0].set_ylabel("ω₁ (rad/s)"); ax[0].set_title("Safe-mode gyrostat spin"); ax[0].grid(true, alpha=0.3)
    ax[1].plot(t, [u[6] for u in sol_att.u], color="tomato")
    ax[1].set_ylabel("ω₂ (rad/s)"); ax[1].grid(true, alpha=0.3)
    ax[2].plot(t, [u[7] for u in sol_att.u], color="green")
    ax[2].set_ylabel("ω₃ (rad/s)"); ax[2].set_xlabel("Time (s)"); ax[2].grid(true, alpha=0.3)
    fig.tight_layout()
    savefig("figs/safe_mode_angular_velocity.png"; dpi=150, bbox_inches="tight")

    # Orbit-attitude simulation
    alt = 480.0; RE = 6378.1363
    r0, v0   = coe_to_eci(RE + alt, 0.0, compute_sso_inclination(alt), 0.0, 0.0, 0.0)
    sun_vec  = @SVector [1.0, 0.0, 0.0]
    q0_orbit = rotation_quat_between(e_panel, sun_vec)
    xi2      = randn(3); xi2 -= dot(xi2, e_panel)*e_panel; xi2 /= norm(xi2)
    w0_orbit = w_des + delta_w .* xi2

    x0 = Vector{Float64}(undef, 13)
    x0[1:3] = r0; x0[4:6] = v0; x0[7:10] = q0_orbit; x0[11:13] = w0_orbit

    sol_orb = integrate_orbit_gyrostat!(x0, I_pert, H_r, e_panel;
        tspan=(0.0, 2000.0), drag=false, reltol=1e-9, abstol=1e-12)

    t = sol_orb.t
    theta_err = pointing_error(sol_orb, e_panel, sun_vec)
    fig2, ax2 = subplots(2, 1, figsize=(8, 6), sharex=true)
    ax2[0].plot(t, [u[7]  for u in sol_orb.u], label="q0", color="black")
    ax2[0].plot(t, [u[8]  for u in sol_orb.u], label="q1", color="steelblue")
    ax2[0].plot(t, [u[9]  for u in sol_orb.u], label="q2", color="tomato")
    ax2[0].plot(t, [u[10] for u in sol_orb.u], label="q3", color="green")
    ax2[0].set_ylabel("Quaternion components"); ax2[0].legend(loc="upper right"); ax2[0].grid(true, alpha=0.3)
    ax2[0].set_title("Orbit–gyrostat: quaternion and pointing error")
    ax2[1].plot(t, theta_err, color="purple")
    ax2[1].set_xlabel("Time (s)"); ax2[1].set_ylabel("Pointing error (deg)"); ax2[1].grid(true, alpha=0.3)
    fig2.tight_layout()
    savefig("figs/orbit_attitude_and_pointing.png"; dpi=150, bbox_inches="tight")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end