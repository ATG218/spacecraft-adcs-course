# Starling attitude control + MEKF. Load after include("environmental_perturbations.jl").
# Experiments: attitude_control_starling.ipynb
using LinearAlgebra
using StaticArrays
using Random
using Statistics

using .EnvironmentalPerturbations: EarthEnvironment,
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

const _STARLING_IMPL_DIR = dirname(@__FILE__)
include(joinpath(_STARLING_IMPL_DIR, "mekf.jl"))

const MU_KM = 398600.4418
const RE_KM = 6378.1363

I_body = @SMatrix [
    0.249  0.0    0.0
    0.0    0.187  0.02
    0.0    0.02   0.148
]
const J_dyn = Matrix{Float64}(I_body)
const Ji_dyn = inv(J_dyn)

e_panel = @SVector [0.0, 1.0, 0.0]
alt_km = 480.0
a_km = RE_KM + alt_km

function compute_sso_inclination(altitude_km)
    mu = MU_KM
    re = RE_KM
    j2 = 1.08262668e-3
    omega_dot = 0.985647 * (pi / 180) / 86400
    a = re + altitude_km
    n = sqrt(mu / a^3)
    return acos(clamp(omega_dot / (-1.5 * j2 * n * (re / a)^2), -1.0, 1.0))
end

function coe_to_eci(a_km, e, inc, raan, argp, nu)
    mu = MU_KM
    p = a_km * (1 - e^2)
    rn = p / (1 + e * cos(nu))
    r_pqw = rn .* @SVector [cos(nu), sin(nu), 0.0]
    v_pqw = sqrt(mu / p) .* @SVector [-sin(nu), e + cos(nu), 0.0]
    cr, sr = cos(raan), sin(raan)
    ci, si = cos(inc), sin(inc)
    cp, sp = cos(argp), sin(argp)
    R = @SMatrix [
        cr * cp - sr * ci * sp   -cr * sp - sr * ci * cp    sr * si
        sr * cp + cr * ci * sp   -sr * sp + cr * ci * cp   -cr * si
        si * sp                   si * cp                    ci
    ]
    return R * r_pqw, R * v_pqw
end

inc = compute_sso_inclination(alt_km)
mean_motion_rad_s = sqrt(MU_KM / a_km^3)

env = EarthEnvironment()
aero = SpacecraftAeroProperties(
    Cd=2.2,
    area_m2=0.03,
    mass_kg=12.0,
    area_normal_body=e_panel,
    center_of_pressure_body_m=SVector(0.02, 0.0, 0.0),
)
srp = SpacecraftSRPProperties(
    Cr=1.3,
    area_m2=0.03,
    area_normal_body=e_panel,
    center_of_pressure_body_m=SVector(0.02, 0.0, 0.0),
    sun_distance_au=1.0,
    eclipse=false,
)

τ_rw_max = 6e-3
rho_rw_max = 50e-3
sun_vec_eci = @SVector [1.0e8, 0.0, 0.0]

Hq = [zeros(1, 3); I(3)]

function Lmat(q::AbstractVector{<:Real})
    q0, q1, q2, q3 = q[1], q[2], q[3], q[4]
    return @SMatrix [q0 -q1 -q2 -q3; q1 q0 -q3 q2; q2 q3 q0 -q1; q3 -q2 q1 q0]
end

function Rmat(q::AbstractVector{<:Real})
    q0, q1, q2, q3 = q[1], q[2], q[3], q[4]
    return @SMatrix [q0 -q1 -q2 -q3; q1 q0 q3 -q2; q2 -q3 q0 q1; q3 q2 -q1 q0]
end

function Gmat(q::AbstractVector{<:Real})
    return Lmat(q) * Hq
end

function Qbody_from_quat(q::AbstractVector{<:Real})
    Tq = [1.0 zeros(1, 3); zeros(3, 1) -I(3)]
    return Hq' * (Rmat(q)' * Lmat(q)) * Hq
end

function quat_conj(q)
    return @SVector [q[1], -q[2], -q[3], -q[4]]
end

function quat_normalize(q)
    n = norm(q)
    return SVector{4,Float64}(q ./ n)
end

function quat_err(q::AbstractVector, qd::AbstractVector)
    qe = Lmat(quat_conj(qd)) * SVector{4,Float64}(q)
    qe = quat_normalize(qe)
    if qe[1] < 0
        qe = -qe
    end
    return qe
end

function expq(ϕ::AbstractVector{<:Real})
    θ = norm(ϕ)
    if θ < 1e-14
        return @SVector [1.0, 0.0, 0.0, 0.0]
    end
    return SVector{4,Float64}([cos(θ); ϕ ./ θ * sin(θ)])
end

function logq(q::AbstractVector{<:Real})
    qv = SVector(q[2], q[3], q[4])
    s = norm(qv)
    if s < 1e-14
        return @SVector [0.0, 0.0, 0.0]
    end
    θ = atan(s, q[1])
    return (θ / s) .* qv
end

function principal_angle_deg(qe::AbstractVector)
    return 2 * rad2deg(acos(clamp(qe[1], -1.0, 1.0)))
end

function random_quat_max_angle(theta_max_rad)
    ax = randn(3)
    ax = ax / norm(ax)
    θ = rand() * theta_max_rad
    return quat_normalize(@SVector [cos(θ / 2), ax[1] * sin(θ / 2), ax[2] * sin(θ / 2), ax[3] * sin(θ / 2)])
end

function sat_vec(u, umax)
    return clamp.(u, -umax, umax)
end

function env_torque_body(q_bi::SVector{4,Float64}, r_km::SVector{3,Float64}, v_kmps::SVector{3,Float64})
    r_m = 1000.0 .* r_km
    v_mps = 1000.0 .* v_kmps
    w = environmental_wrench(
        I_body,
        q_bi,
        r_m,
        v_mps;
        sun_eci_m=sun_vec_eci,
        env=env,
        aero=aero,
        srp=srp,
        include_gravity_gradient=true,
        include_drag=true,
        include_srp=false,
    )
    return SVector{3,Float64}(w.torque_body_Nm)
end

function dynamics_rhs(q, ω, rho, τ_dist, u_cmd)
    qn = quat_normalize(q)
    ωv = SVector{3,Float64}(ω)
    rhov = SVector{3,Float64}(rho)
    u = SVector{3,Float64}(u_cmd)
    τd = SVector{3,Float64}(τ_dist)
    qdot = 0.5 .* (Lmat(qn) * SVector{4,Float64}(0.0, ωv[1], ωv[2], ωv[3]))
    # Lecture convention: H = J*ω + rho is total angular momentum, while rho
    # is the internal reaction-wheel momentum; rho_dot = -u for bus torque u.
    ωdot = Ji_dyn * (τd + u - cross(ωv, J_dyn * ωv + rhov))
    rhodot = -Vector(u)
    return Vector(qdot), Vector(ωdot), Vector(rhodot)
end

function rk4_step(q, ω, rho, τ_dist, u_cmd, dt)
    function f(qv, ωv, rhov)
        return dynamics_rhs(qv, ωv, rhov, τ_dist, u_cmd)
    end
    k1q, k1w, k1rho = f(q, ω, rho)
    k2q, k2w, k2rho = f(q .+ 0.5 * dt .* k1q, ω .+ 0.5 * dt .* k1w, rho .+ 0.5 * dt .* k1rho)
    k3q, k3w, k3rho = f(q .+ 0.5 * dt .* k2q, ω .+ 0.5 * dt .* k2w, rho .+ 0.5 * dt .* k2rho)
    k4q, k4w, k4rho = f(q .+ dt .* k3q, ω .+ dt .* k3w, rho .+ dt .* k3rho)
    qn = q .+ (dt / 6) .* (k1q .+ 2 .* k2q .+ 2 .* k3q .+ k4q)
    ωn = ω .+ (dt / 6) .* (k1w .+ 2 .* k2w .+ 2 .* k3w .+ k4w)
    rhon = rho .+ (dt / 6) .* (k1rho .+ 2 .* k2rho .+ 2 .* k3rho .+ k4rho)
    qn = Vector(quat_normalize(qn))
    rhon = sat_vec(rhon, rho_rw_max)
    return qn, Vector(ωn), Vector(rhon)
end

struct PDGains
    kp::Float64
    kd::Float64
end

function torque_pd(q::AbstractVector, ω::AbstractVector, qd::AbstractVector, g::PDGains)
    qe = quat_err(q, qd)
    ϕ = logq(qe)
    τ = -g.kp .* ϕ .- g.kd .* ω
    return τ, ϕ, qe
end

function meas_stack(q, rN_mat::Matrix{Float64})
    m = size(rN_mat, 2)
    y = zeros(3 * m)
    Qk = Qbody_from_quat(q)
    for k in 1:m
        r = rN_mat[:, k]
        # Body measurement of inertial unit vector: b = Q(q)' * r (see mekf-simple.ipynb)
        y[(3*(k-1)+1):(3*k)] .= Qk' * r
    end
    return y
end

Tq = [1.0 zeros(1, 3); zeros(3, 1) -I(3)]

function meas_jacobian(q, rN_mat::Matrix{Float64})
    m = size(rN_mat, 2)
    C = zeros(3 * m, 3)
    for k in 1:m
        r = SVector{3,Float64}(rN_mat[:, k])
        rows = (3*(k-1)+1):(3*k)
        Hr = Hq * r
        C[rows, :] .= Hq' * (Lmat(q)' * Lmat(Hr) + Rmat(q) * Rmat(Hr) * Tq) * Gmat(q)
    end
    return C
end

function mekf_predict(q, P, ω_gyro, dt, V)
    Δq = expq(0.5 * dt .* ω_gyro)
    qpred = quat_normalize(Lmat(q) * Δq)
    A = Matrix(Gmat(qpred)' * Rmat(Δq) * Gmat(q))
    Ppred = A * P * A' + V
    return Vector(qpred), Ppred, A
end

function mekf_update(qpred, Ppred, y_meas, rN_mat, W)
    yhat = meas_stack(qpred, rN_mat)
    z = y_meas - yhat
    C = meas_jacobian(qpred, rN_mat)
    S = C * Ppred * C' + W
    S = S + (1e-9 + 1e-9 * tr(S) / size(S, 1)) * I(size(S, 1))
    K = Ppred * C' / S
    δ = K * z
    ϕ = δ[1:3]
    ns = max(0.0, 1.0 - dot(ϕ, ϕ))
    dq = [sqrt(ns); ϕ[1]; ϕ[2]; ϕ[3]]
    dq = Vector(quat_normalize(dq))
    qn = Vector(quat_normalize(Lmat(qpred) * dq))
    I3 = I(3)
    Pn = (I3 - K * C) * Ppred * (I3 - K * C)' + K * W * K'
    return qn, Pn
end

function inertial_refs(r_km::SVector{3,Float64})
    ri = normalize(-1000.0 .* r_km)
    si = normalize(sun_vec_eci)
    return Matrix(hcat(si, ri))
end

function simulate_closed_loop(; Tfinal, dt, q0, omega0, rho0, qd,
        g::PDGains,
        sigma_gyro,
        sigma_vec,
        V_mekf,
        rng::AbstractRNG,
        control_from_truth::Bool=false,
        control_gyro_lpf_tau::Float64=0.0,
        professor_mekf::Bool=true,
        sigma_bias_walk::Float64=1e-5,
    )
    nst = max(1, round(Int, Tfinal / dt))
    n = nst + 1
    q = Vector{Float64}(quat_normalize(q0))
    ω = Vector{Float64}(omega0)
    rho = Vector{Float64}(rho0)
    β_true = zeros(3)
    ν = 0.0
    ω_lpf = Vector{Float64}(omega0)

    # Initial estimate: small multiplicative perturbation of truth (MEKF will correct)
    dq0 = expq(0.08 .* randn(rng, 3))
    qhat = Vector(quat_normalize(Lmat(q) * dq0))
    P = Matrix(0.25 .* I(3))
    # Course 7-state MEKF (mekf.jl, module MEKF) — quaternion + gyro bias
    qhat_s = SVector{4,Float64}(qhat[1], qhat[2], qhat[3], qhat[4])
    x_prof = SVector{7,Float64}(qhat_s[1], qhat_s[2], qhat_s[3], qhat_s[4], 0.0, 0.0, 0.0)
    P_prof = Matrix{Float64}(0.5 * I(6))
    Qproc_prof = zeros(6, 6)
    Qproc_prof[1:3, 1:3] .= (sigma_gyro^2) * (dt^2) * I(3)
    Qproc_prof[4:6, 4:6] .= (sigma_bias_walk^2) * dt * I(3)
    R3 = Matrix{Float64}((sigma_vec^2) * I(3))
    R_vec_list = Matrix{Float64}[R3, R3]

    r_km, v_kmps = coe_to_eci(a_km, 0.0, inc, 0.0, 0.0, ν)
    r_km = SVector{3,Float64}(r_km)
    v_kmps = SVector{3,Float64}(v_kmps)
    rN_init = inertial_refs(r_km)
    mref = size(rN_init, 2)
    Wmekf = (sigma_vec^2) * I(3 * mref)

    q_hist = zeros(4, n)
    qh_hist = zeros(4, n)
    ω_hist = zeros(3, n)
    beta_hist = zeros(3, n)
    betahat_hist = zeros(3, n)
    rho_hist = zeros(3, n)
    tau_hist = zeros(3, n)
    ang_hist = zeros(n)
    ang_est_hist = zeros(n)

    for k in 1:n
        q_hist[:, k] .= q
        qh_hist[:, k] .= qhat
        ω_hist[:, k] .= ω
        beta_hist[:, k] .= β_true
        rho_hist[:, k] .= rho
        if professor_mekf
            betahat_hist[:, k] .= x_prof[5:7]
        end
        qe_t = quat_err(q, qd)
        ang_hist[k] = principal_angle_deg(qe_t)
        ang_est_hist[k] = principal_angle_deg(quat_err(q, qhat))

        τ_dist = env_torque_body(quat_normalize(SVector{4,Float64}(q)), r_km, v_kmps)
        # Lecture gyro model: measured rate = true rate + bias + white noise.
        # The course MEKF estimates beta; the simple MEKF path is retained only
        # for comparison/debugging and does not remove bias.
        ω_meas = ω .+ β_true .+ sigma_gyro .* randn(rng, 3)
        if control_gyro_lpf_tau > 0.0
            α_lpf = 1.0 - exp(-dt / control_gyro_lpf_tau)
            ω_lpf .= (1.0 - α_lpf) .* ω_lpf .+ α_lpf .* ω_meas
        else
            ω_lpf .= ω_meas
        end
        q_ctrl = control_from_truth ? q : qhat
        if control_from_truth
            ω_ctrl = ω
        elseif professor_mekf
            β̂ = SVector{3,Float64}(x_prof[5], x_prof[6], x_prof[7])
            ω_ctrl = ω_lpf .- Vector(β̂)
        else
            ω_ctrl = ω_lpf
        end

        τ, _, _ = torque_pd(q_ctrl, ω_ctrl, qd, g)
        τ_sat = sat_vec(τ, τ_rw_max)
        tau_hist[:, k] .= τ_sat

        if k < n
            q, ω, rho = rk4_step(q, ω, rho, τ_dist, τ_sat, dt)
            β_true .= β_true .+ sigma_bias_walk * sqrt(dt) .* randn(rng, 3)
            ν = ν + mean_motion_rad_s * dt
            r_km, v_kmps = coe_to_eci(a_km, 0.0, inc, 0.0, 0.0, ν)
            r_km = SVector{3,Float64}(r_km)
            v_kmps = SVector{3,Float64}(v_kmps)
            rN_post = inertial_refs(r_km)
            y_meas = meas_stack(q, rN_post) .+ sigma_vec .* randn(rng, 3 * mref)
            if professor_mekf
                ω_meas_sv = SVector{3,Float64}(ω_meas[1], ω_meas[2], ω_meas[3])
                x_pred, P_pred = MEKF.mekf_predict(x_prof, P_prof, ω_meas_sv, dt, Qproc_prof)
                q_pred_sv = SVector{4,Float64}(x_pred[1], x_pred[2], x_pred[3], x_pred[4])
                vec_refs_k = SVector{3,Float64}[
                    normalize(SVector{3,Float64}(rN_post[:, 1])),
                    normalize(SVector{3,Float64}(rN_post[:, 2])),
                ]
                vec_meas_k = SVector{3,Float64}[
                    normalize(SVector{3,Float64}(y_meas[1:3])),
                    normalize(SVector{3,Float64}(y_meas[4:6])),
                ]
                zv, Cv, Wv = MEKF.build_vector_measurement(q_pred_sv, vec_meas_k, vec_refs_k, R_vec_list)
                # Tiny diagonal load on W for numerical conditioning (optional in course script)
                Wv = Wv + (1e-12 + 1e-12 * tr(Wv) / size(Wv, 1)) * I(size(Wv, 1))
                x_prof, P_prof = MEKF.mekf_update(x_pred, P_pred, zv, Cv, Wv)
                qhat[1], qhat[2], qhat[3], qhat[4] = x_prof[1], x_prof[2], x_prof[3], x_prof[4]
                nq = norm(qhat)
                qhat ./= nq
            else
                qpred, Ppred, _ = mekf_predict(qhat, P, ω_meas, dt, V_mekf)
                qhat, P = mekf_update(qpred, Ppred, y_meas, rN_post, Wmekf)
            end
        end
    end

    t = range(0.0; step=dt, length=n)
    return (; t, q_hist, qh_hist, omega_hist=ω_hist, beta_hist, betahat_hist, rho_hist, tau_hist, ang_hist, ang_est_hist)
end

function rms_steady(ang_deg::Vector{Float64}; skip_frac=0.25)
    n = length(ang_deg)
    i0 = Int(floor(skip_frac * n)) + 1
    v = @view ang_deg[i0:end]
    return sqrt(mean(abs2.(v)))
end

# Aliases for notebook convenience after `include` (same names as report4)
Ji = Ji_dyn
J = J_dyn
