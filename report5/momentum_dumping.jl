# Momentum dumping helpers for report5.
#
# Load this after the report4 modules:
#   include("../report4/environmental_perturbations.jl")
#   using .EnvironmentalPerturbations
#   include("../report4/attitude_control_starling_impl.jl")

using LinearAlgebra
using StaticArrays
using Random
using Statistics

const MD_EARTH_DIPOLE_EQUATOR_T = 3.12e-5

function orbit_period_seconds_md(cfg::StarlingConfig)
    return 2π / cfg.mean_motion_rad_s
end

function environmental_torque_body_md(q_bi::SVector{4,Float64},
                                      r_km::SVector{3,Float64},
                                      v_kmps::SVector{3,Float64},
                                      cfg::StarlingConfig;
                                      include_srp::Bool=false)
    r_m = 1000.0 .* r_km
    v_mps = 1000.0 .* v_kmps
    wrench = environmental_wrench(
        cfg.I_body, q_bi, r_m, v_mps;
        sun_eci_m=cfg.sun_vec_eci,
        env=cfg.env,
        aero=cfg.aero,
        srp=cfg.srp,
        include_gravity_gradient=true,
        include_drag=true,
        include_srp=include_srp,
    )
    return SVector{3,Float64}(wrench.torque_body_Nm)
end

function mean_environment_torque_body(cfg::StarlingConfig;
                                      q_ref::SVector{4,Float64}=SVector{4,Float64}(1.0, 0.0, 0.0, 0.0),
                                      samples::Int=720,
                                      include_srp::Bool=false)
    torque_sum = zeros(3)
    for j in 0:(samples - 1)
        nu = 2π * j / samples
        r_km, v_kmps = coe_to_eci(cfg.a_km, 0.0, cfg.inc_rad, 0.0, 0.0, nu)
        tau = environmental_torque_body_md(
            q_ref,
            SVector{3,Float64}(r_km),
            SVector{3,Float64}(v_kmps),
            cfg;
            include_srp=include_srp,
        )
        torque_sum .+= Vector(tau)
    end
    return SVector{3,Float64}(torque_sum ./ samples)
end

function estimate_stress_scale_for_saturation(cfg::StarlingConfig;
                                              q_ref::SVector{4,Float64}=SVector{4,Float64}(1.0, 0.0, 0.0, 0.0),
                                              Tfinal::Float64,
                                              target_fraction::Float64=1.15,
                                              samples::Int=720,
                                              include_srp::Bool=false)
    mean_tau = mean_environment_torque_body(
        cfg;
        q_ref=q_ref,
        samples=samples,
        include_srp=include_srp,
    )
    nominal_axis_growth = abs.(mean_tau) .* Tfinal
    denom = maximum(nominal_axis_growth)
    scale = denom > 0.0 ? max(1.0, target_fraction * cfg.rho_rw_max / denom) : Inf
    return (;
        scale,
        mean_torque_body_Nm=mean_tau,
        nominal_axis_growth_Nms=nominal_axis_growth,
        target_axis_growth_Nms=target_fraction * cfg.rho_rw_max,
    )
end

function dipole_magnetic_field_eci(r_km::SVector{3,Float64};
                                   equator_field_T::Float64=MD_EARTH_DIPOLE_EQUATOR_T,
                                   dipole_axis_eci::SVector{3,Float64}=SVector{3,Float64}(0.0, 0.0, 1.0),
                                   earth_radius_m::Float64=6_378_136.3)
    r_m = 1000.0 .* r_km
    r_norm = norm(r_m)
    r_norm == 0.0 && return @SVector [0.0, 0.0, 0.0]

    r_hat = r_m / r_norm
    m_hat = normalize(dipole_axis_eci)
    return equator_field_T * (earth_radius_m / r_norm)^3 *
           (3.0 * dot(m_hat, r_hat) * r_hat - m_hat)
end

function magnetic_field_body_md(q_bi::SVector{4,Float64},
                                r_km::SVector{3,Float64};
                                equator_field_T::Float64=MD_EARTH_DIPOLE_EQUATOR_T,
                                dipole_axis_eci::SVector{3,Float64}=SVector{3,Float64}(0.0, 0.0, 1.0))
    b_eci = dipole_magnetic_field_eci(
        r_km;
        equator_field_T=equator_field_T,
        dipole_axis_eci=dipole_axis_eci,
    )
    return quat_to_dcm_body_to_inertial(q_bi)' * b_eci
end

function magnetorquer_unloading_command(B_body::SVector{3,Float64},
                                        rho::SVector{3,Float64};
                                        dump_gain::Float64=2.0e-4,
                                        dipole_max::Float64=0.1,
                                        deadband::Float64=0.0)
    B2 = dot(B_body, B_body)
    if B2 < eps(Float64) || norm(rho) <= deadband
        z = @SVector [0.0, 0.0, 0.0]
        return (; dipole_body_Am2=z, torque_body_Nm=z, desired_torque_body_Nm=z)
    end

    tau_des = -dump_gain .* rho
    B_hat = B_body / sqrt(B2)
    tau_perp = tau_des - dot(tau_des, B_hat) .* B_hat
    m_cmd = cross(B_body, tau_perp) ./ B2
    m_sat = SVector{3,Float64}(clamp.(m_cmd, -dipole_max, dipole_max))
    tau_mag = cross(m_sat, B_body)

    return (;
        dipole_body_Am2=m_sat,
        torque_body_Nm=tau_mag,
        desired_torque_body_Nm=tau_des,
    )
end

function dynamics_rhs_md(q, omega, rho, external_torque, wheel_torque, cfg::StarlingConfig)
    return dynamics_rhs(q, omega, rho, external_torque, wheel_torque, cfg)
end

function momentum_limited_wheel_torque(tau_cmd::AbstractVector{<:Real},
                                       rho::AbstractVector{<:Real},
                                       dt::Real,
                                       cfg::StarlingConfig)
    dt <= 0 && throw(ArgumentError("dt must be positive for wheel momentum limiting"))
    tau = Vector{Float64}(sat_vec(tau_cmd, cfg.τ_rw_max))
    rho_max = cfg.rho_rw_max
    for i in eachindex(tau)
        rho_i = clamp(float(rho[i]), -rho_max, rho_max)
        tau_min = (rho_i - rho_max) / dt
        tau_max = (rho_i + rho_max) / dt
        tau[i] = clamp(tau[i], tau_min, tau_max)
    end
    return tau
end

function rk4_step_md(q, omega, rho, external_torque, wheel_torque, dt, cfg::StarlingConfig)
    function f(qv, omegav, rhov)
        return dynamics_rhs_md(qv, omegav, rhov, external_torque, wheel_torque, cfg)
    end
    k1q, k1w, k1rho = f(q, omega, rho)
    k2q, k2w, k2rho = f(q .+ 0.5 * dt .* k1q, omega .+ 0.5 * dt .* k1w, rho .+ 0.5 * dt .* k1rho)
    k3q, k3w, k3rho = f(q .+ 0.5 * dt .* k2q, omega .+ 0.5 * dt .* k2w, rho .+ 0.5 * dt .* k2rho)
    k4q, k4w, k4rho = f(q .+ dt .* k3q, omega .+ dt .* k3w, rho .+ dt .* k3rho)

    qn = q .+ (dt / 6) .* (k1q .+ 2 .* k2q .+ 2 .* k3q .+ k4q)
    omegan = omega .+ (dt / 6) .* (k1w .+ 2 .* k2w .+ 2 .* k3w .+ k4w)
    rhon = rho .+ (dt / 6) .* (k1rho .+ 2 .* k2rho .+ 2 .* k3rho .+ k4rho)

    qn = Vector(quat_normalize(qn))
    rhon = sat_vec(rhon, cfg.rho_rw_max)
    return qn, Vector(omegan), Vector(rhon)
end

function simulate_momentum_dumping(cfg::StarlingConfig; Tfinal, dt, q0, omega0, rho0, qd,
        g::PDGains,
        sigma_gyro,
        sigma_vec,
        V_mekf,
        rng::AbstractRNG,
        control_from_truth::Bool=false,
        control_gyro_lpf_tau::Float64=0.0,
        professor_mekf::Bool=true,
        sigma_bias_walk::Float64=1e-5,
        disturbance_scale::Float64=1.0,
        include_srp::Bool=false,
        enable_dumping::Bool=false,
        feedforward_external_torque::Bool=false,
        wheel_torque_mode::Symbol=:controlled,
        enforce_wheel_momentum_limit::Bool=true,
        dump_gain::Float64=2.0e-4,
        dipole_max::Float64=0.1,
        dump_deadband::Float64=0.0,
        magnetic_equator_T::Float64=MD_EARTH_DIPOLE_EQUATOR_T,
    )
    nst = max(1, round(Int, Tfinal / dt))
    n = nst + 1
    q = Vector{Float64}(quat_normalize(q0))
    omega = Vector{Float64}(omega0)
    rho = Vector{Float64}(rho0)
    beta_true = zeros(3)
    nu = 0.0
    omega_lpf = Vector{Float64}(omega0)

    dq0 = expq(0.08 .* randn(rng, 3))
    qhat = Vector(quat_normalize(Lmat(q) * dq0))
    P = Matrix(0.25 .* I(3))
    qhat_s = SVector{4,Float64}(qhat[1], qhat[2], qhat[3], qhat[4])
    x_prof = SVector{7,Float64}(qhat_s[1], qhat_s[2], qhat_s[3], qhat_s[4], 0.0, 0.0, 0.0)
    P_prof = zeros(6, 6)
    P_prof[1:3, 1:3] .= 0.25 * I(3)
    initial_bias_var = max((10.0 * sigma_gyro)^2, (sigma_bias_walk * sqrt(dt))^2)
    P_prof[4:6, 4:6] .= initial_bias_var * I(3)
    Qproc_prof = zeros(6, 6)
    Qproc_prof[1:3, 1:3] .= (sigma_gyro^2) * (dt^2) * I(3)
    Qproc_prof[4:6, 4:6] .= (sigma_bias_walk^2) * dt * I(3)
    R3 = Matrix{Float64}((sigma_vec^2) * I(3))
    R_vec_list = Matrix{Float64}[R3, R3]
    omega_meas = omega .+ beta_true .+ sigma_gyro .* randn(rng, 3)

    r_km, v_kmps = coe_to_eci(cfg.a_km, 0.0, cfg.inc_rad, 0.0, 0.0, nu)
    r_km = SVector{3,Float64}(r_km)
    v_kmps = SVector{3,Float64}(v_kmps)
    rN_init = inertial_refs(r_km, cfg)
    mref = size(rN_init, 2)
    Wmekf = (sigma_vec^2) * I(3 * mref)

    q_hist = zeros(4, n)
    qh_hist = zeros(4, n)
    omega_hist = zeros(3, n)
    beta_hist = zeros(3, n)
    betahat_hist = zeros(3, n)
    rho_hist = zeros(3, n)
    tau_hist = zeros(3, n)
    tau_env_hist = zeros(3, n)
    tau_mag_hist = zeros(3, n)
    tau_external_hist = zeros(3, n)
    dipole_hist = zeros(3, n)
    B_body_hist = zeros(3, n)
    ang_hist = zeros(n)
    ang_est_hist = zeros(n)

    for k in 1:n
        q_hist[:, k] .= q
        qh_hist[:, k] .= qhat
        omega_hist[:, k] .= omega
        beta_hist[:, k] .= beta_true
        rho_hist[:, k] .= rho
        if professor_mekf
            betahat_hist[:, k] .= x_prof[5:7]
        end

        qe_t = quat_err(q, qd)
        ang_hist[k] = principal_angle_deg(qe_t)
        ang_est_hist[k] = principal_angle_deg(quat_err(q, qhat))

        q_sv = quat_normalize(SVector{4,Float64}(q[1], q[2], q[3], q[4]))
        tau_env = disturbance_scale .* environmental_torque_body_md(
            q_sv, r_km, v_kmps, cfg;
            include_srp=include_srp,
        )
        B_body = magnetic_field_body_md(q_sv, r_km; equator_field_T=magnetic_equator_T)
        rho_sv = SVector{3,Float64}(rho[1], rho[2], rho[3])
        dump = enable_dumping ? magnetorquer_unloading_command(
            B_body,
            rho_sv;
            dump_gain=dump_gain,
            dipole_max=dipole_max,
            deadband=dump_deadband,
        ) : magnetorquer_unloading_command(
            B_body,
            SVector{3,Float64}(0.0, 0.0, 0.0);
            dump_gain=dump_gain,
            dipole_max=dipole_max,
            deadband=Inf,
        )
        tau_mag = dump.torque_body_Nm
        tau_external = tau_env + tau_mag

        tau_env_hist[:, k] .= tau_env
        tau_mag_hist[:, k] .= tau_mag
        tau_external_hist[:, k] .= tau_external
        dipole_hist[:, k] .= dump.dipole_body_Am2
        B_body_hist[:, k] .= B_body

        if control_gyro_lpf_tau > 0.0
            alpha_lpf = 1.0 - exp(-dt / control_gyro_lpf_tau)
            omega_lpf .= (1.0 - alpha_lpf) .* omega_lpf .+ alpha_lpf .* omega_meas
        else
            omega_lpf .= omega_meas
        end

        q_ctrl = control_from_truth ? q : qhat
        if control_from_truth
            omega_ctrl = omega
        elseif professor_mekf
            beta_hat = SVector{3,Float64}(x_prof[5], x_prof[6], x_prof[7])
            omega_ctrl = omega_lpf .- Vector(beta_hat)
        else
            omega_ctrl = omega_lpf
        end

        tau_pd, _, _ = torque_pd(q_ctrl, omega_ctrl, qd, g)
        tau_cmd = feedforward_external_torque ? tau_pd .- tau_external : tau_pd
        tau_sat = if wheel_torque_mode === :none
            zeros(3)
        elseif enforce_wheel_momentum_limit
            momentum_limited_wheel_torque(tau_cmd, rho, dt, cfg)
        else
            sat_vec(tau_cmd, cfg.τ_rw_max)
        end
        tau_hist[:, k] .= tau_sat

        if k < n
            q, omega, rho = rk4_step_md(q, omega, rho, tau_external, tau_sat, dt, cfg)
            beta_true .= beta_true .+ sigma_bias_walk * sqrt(dt) .* randn(rng, 3)
            nu = nu + cfg.mean_motion_rad_s * dt
            r_km, v_kmps = coe_to_eci(cfg.a_km, 0.0, cfg.inc_rad, 0.0, 0.0, nu)
            r_km = SVector{3,Float64}(r_km)
            v_kmps = SVector{3,Float64}(v_kmps)
            rN_post = inertial_refs(r_km, cfg)
            y_meas = meas_stack(q, rN_post) .+ sigma_vec .* randn(rng, 3 * mref)
            omega_meas_next = omega .+ beta_true .+ sigma_gyro .* randn(rng, 3)
            omega_meas_interval = 0.5 .* (omega_meas .+ omega_meas_next)
            if professor_mekf
                omega_meas_sv = SVector{3,Float64}(
                    omega_meas_interval[1],
                    omega_meas_interval[2],
                    omega_meas_interval[3],
                )
                x_pred, P_pred = MEKF.mekf_predict(x_prof, P_prof, omega_meas_sv, dt, Qproc_prof)
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
                Wv = Wv + (1e-12 + 1e-12 * tr(Wv) / size(Wv, 1)) * I(size(Wv, 1))
                x_prof, P_prof = MEKF.mekf_update(x_pred, P_pred, zv, Cv, Wv)
                qhat[1], qhat[2], qhat[3], qhat[4] = x_prof[1], x_prof[2], x_prof[3], x_prof[4]
                qhat ./= norm(qhat)
            else
                qpred, Ppred, _ = mekf_predict(qhat, P, omega_meas_interval, dt, V_mekf)
                qhat, P = mekf_update(qpred, Ppred, y_meas, rN_post, Wmekf)
            end
            omega_meas = omega_meas_next
        end
    end

    t = range(0.0; step=dt, length=n)
    return (;
        t,
        q_hist,
        qh_hist,
        omega_hist,
        beta_hist,
        betahat_hist,
        rho_hist,
        tau_hist,
        tau_env_hist,
        tau_mag_hist,
        tau_external_hist,
        dipole_hist,
        B_body_hist,
        ang_hist,
        ang_est_hist,
        disturbance_scale,
        enable_dumping,
        feedforward_external_torque,
        wheel_torque_mode,
        enforce_wheel_momentum_limit,
        dump_gain,
        dipole_max,
    )
end

function max_abs_wheel_momentum(sol)
    return maximum(abs.(sol.rho_hist))
end

function saturation_fraction(sol, limit::Float64; threshold::Float64=0.999)
    saturated = [any(abs.(sol.rho_hist[:, k]) .>= threshold * limit) for k in axes(sol.rho_hist, 2)]
    return count(saturated) / length(saturated)
end

function first_saturation_time(sol, limit::Float64; threshold::Float64=0.999)
    for k in axes(sol.rho_hist, 2)
        if any(abs.(sol.rho_hist[:, k]) .>= threshold * limit)
            return sol.t[k]
        end
    end
    return nothing
end

function run_summary(label::AbstractString, sol, cfg::StarlingConfig; skip_frac::Float64=0.2)
    B_norm = [norm(sol.B_body_hist[:, k]) for k in axes(sol.B_body_hist, 2)]
    return (;
        label,
        disturbance_scale=sol.disturbance_scale,
        dumping_enabled=sol.enable_dumping,
        feedforward_external_torque=sol.feedforward_external_torque,
        max_abs_rho_Nms=max_abs_wheel_momentum(sol),
        rho_limit_Nms=cfg.rho_rw_max,
        first_saturation_s=first_saturation_time(sol, cfg.rho_rw_max),
        saturation_fraction=saturation_fraction(sol, cfg.rho_rw_max),
        rms_pointing_deg=rms_steady(collect(sol.ang_hist); skip_frac=skip_frac),
        max_dipole_Am2=maximum(abs.(sol.dipole_hist)),
        dipole_limit_Am2=sol.dipole_max,
        mean_B_uT=1e6 * mean(B_norm),
        min_B_uT=1e6 * minimum(B_norm),
        max_B_uT=1e6 * maximum(B_norm),
    )
end

function print_run_summary(summary)
    sat_time = summary.first_saturation_s === nothing ? "none" :
               string(round(summary.first_saturation_s; digits=1), " s")
    println(summary.label)
    println("  disturbance scale: ", round(summary.disturbance_scale; sigdigits=5))
    println("  dumping enabled: ", summary.dumping_enabled)
    println("  max |rho_i|: ", round(summary.max_abs_rho_Nms; sigdigits=5),
            " N*m*s / ", summary.rho_limit_Nms)
    println("  first saturation: ", sat_time)
    println("  saturated time fraction: ", round(100 * summary.saturation_fraction; digits=2), "%")
    println("  RMS pointing: ", round(summary.rms_pointing_deg; digits=4), " deg")
    println("  max |m_i|: ", round(summary.max_dipole_Am2; sigdigits=5),
            " A*m^2 / ", summary.dipole_limit_Am2)
    println("  B-field range: ",
            round(summary.min_B_uT; digits=2), " to ",
            round(summary.max_B_uT; digits=2), " uT")
end
