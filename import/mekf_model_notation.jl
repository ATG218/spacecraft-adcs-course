# MEKF implementation using notation from hw2_simulation.jl and report2.ipynb
#
# This script provides a complete implementation of the multiplicative
# extended Kalman filter (MEKF) for spacecraft attitude estimation.  It
# follows the quaternion notation and helper functions defined in the
# provided homework files (hw2_simulation.jl, safe_mode_gyrostat.jl)
# and the report notebook report2.ipynb.  Quaternion vectors are
# ordered as [scalar; vector] = [q0; q1; q2; q3], and rotation
# matrices map body vectors to inertial coordinates via Q(q): r_N = Q(q)*r_B.
#
# The filter estimates the attitude quaternion and gyro bias using
# gyroscope, vector sensor (sun and magnetometer) and star tracker
# measurements.  A basic simulation is included at the end to
# illustrate the filter in action.

using LinearAlgebra
using StaticArrays
using Random

# Bring in dynamics helpers (integrate_gyrostat!, rpm2rad, etc.)
include("safe_mode_gyrostat.jl")

##############################
# Quaternion helper functions
##############################

"Skew‑symmetric matrix for a 3‑vector v"
@inline function hat(v::SVector{3,Float64})::SMatrix{3,3,Float64,9}
    return @SMatrix [   0.0    -v[3]   v[2];
                       v[3]     0.0   -v[1];
                      -v[2]    v[1]    0.0]
end

"Inverse of hat: extract the vector from a skew‑symmetric matrix"
@inline function unhat(S::SMatrix{3,3,Float64,9})::SVector{3,Float64}
    return @SVector [0.5*(S[3,2] - S[2,3]),
                     0.5*(S[1,3] - S[3,1]),
                     0.5*(S[2,1] - S[1,2])]
end

# Matrix H picks out the vector part of a quaternion.  See report2.ipynb.
const H = @SMatrix [0.0 0.0 0.0;
                    1.0 0.0 0.0;
                    0.0 1.0 0.0;
                    0.0 0.0 1.0]

"Left quaternion multiplication matrix: L(q)*p == q ⊗ p"
function L(q::SVector{4,Float64})::SMatrix{4,4,Float64,16}
    q0 = q[1]; qv = @SVector [q[2], q[3], q[4]]
    return @SMatrix [ q0       -qv[1]   -qv[2]   -qv[3];
                      qv[1]     q0       qv[3]   -qv[2];
                      qv[2]    -qv[3]    q0       qv[1];
                      qv[3]     qv[2]   -qv[1]    q0     ]
end

"Right quaternion multiplication matrix: R(p)*q == q ⊗ p"
function R(p::SVector{4,Float64})::SMatrix{4,4,Float64,16}
    p0 = p[1]; pv = @SVector [p[2], p[3], p[4]]
    return @SMatrix [ p0       -pv[1]   -pv[2]   -pv[3];
                      pv[1]     p0      -pv[3]    pv[2];
                      pv[2]     pv[3]     p0     -pv[1];
                      pv[3]    -pv[2]     pv[1]    p0   ]
end

"Quaternion to 3×3 rotation matrix mapping body → inertial: r_N = Q(q) * r_B"
function Q(q::SVector{4,Float64})::SMatrix{3,3,Float64,9}
    # Using the definition from report2.ipynb: Q(q) = H' * (R(q)' * L(q)) * H
    return transpose(H) * (transpose(R(q)) * L(q)) * H
end

"Quaternion exponential map: convert a small rotation vector φ into a unit quaternion.
   The input φ should be half the desired rotation vector (see usage below)."
function expq(phi::SVector{3,Float64})::SVector{4,Float64}
    θ = norm(phi)
    if θ < 1e-12
        return @SVector [1.0, 0.0, 0.0, 0.0]
    else
        # sinc(x) in Julia is sin(πx)/(πx), so sinc(θ/π) = sin(θ)/θ
        vv = phi * sinc(θ / π)
        return SVector{4,Float64}(cos(θ), vv[1], vv[2], vv[3])
    end
end

function expq(phi::AbstractVector{<:Real})::SVector{4,Float64}
    return expq(SVector{3,Float64}(phi[1], phi[2], phi[3]))
end

"Quaternion logarithm map: returns the rotation vector for a unit quaternion q.
   For small rotations, this map is well‑behaved."
function logq(q::SVector{4,Float64})::SVector{3,Float64}
    c = q[1]
    qv = @SVector [q[2], q[3], q[4]]
    s = norm(qv)
    θ = atan(s, c)
    # Use sinc(θ/π) = sin(θ)/θ to avoid division by zero at small θ
    return qv / sinc(θ / π)
end

"Quaternion conjugate (inverse for unit quaternions)"
@inline function qconj(q::SVector{4,Float64})::SVector{4,Float64}
    return @SVector [q[1], -q[2], -q[3], -q[4]]
end

##############################
# Sensor noise models
##############################

"Sample Gaussian noise with covariance R_cov (3×3 matrix)"
function sample_noise(R_cov::AbstractMatrix{Float64})::SVector{3,Float64}
    Lchol = cholesky(Symmetric(Matrix(R_cov))).L
    return SVector{3,Float64}(Lchol * randn(3))
end

"Generate a noisy unit vector measurement: v_true + noise, then normalize."
function noisy_vector_meas(v_true::SVector{3,Float64}, R_cov::AbstractMatrix{Float64})::SVector{3,Float64}
    v = v_true + sample_noise(R_cov)
    return v / norm(v)
end

"Generate a noisy star‑tracker quaternion: true quaternion q_true perturbed
   by a small rotation δφ ~ N(0, R_cov).  The noise vector δφ is applied
   in the body frame and converted to a quaternion via expq(0.5*δφ)."
function noisy_star_tracker_meas(q_true::SVector{4,Float64}, R_cov::AbstractMatrix{Float64})::SVector{4,Float64}
    δφ = sample_noise(R_cov)
    δq = expq(0.5 * δφ)
    q_meas = L(q_true) * δq
    return q_meas / norm(q_meas)
end

##############################
# MEKF prediction and update
##############################

"State prediction for the MEKF with gyro bias.  The state x consists of
   the attitude quaternion q (4 elements) and gyro bias β (3 elements).
   The measured angular rate u is used along with the bias estimate to
   propagate the attitude over a step of duration h.  The bias is
   assumed constant over the step."
function state_prediction(x::SVector{7,Float64}, u::SVector{3,Float64}, h::Float64)
    q = SVector{4,Float64}(x[1], x[2], x[3], x[4])
    β = SVector{3,Float64}(x[5], x[6], x[7])
    ωc = u - β
    # The rotation vector over the interval is φ = h*ωc.  Pass half of it to expq.
    Δq = expq(0.5 * h * ωc)
    qn = L(q) * Δq
    qn /= norm(qn)
    return SVector{7,Float64}(qn[1], qn[2], qn[3], qn[4], β[1], β[2], β[3])
end

"Linearised state transition Jacobian for the error state [φ; δβ].  This
   function computes the discrete‑time matrix A such that
   δx_{k+1} ≈ A * δx_k, where δx = [φ; δβ].  The derivation follows
   Lecture 12: δφ̇ = −ωc× δφ − δβ.  We approximate expm for a small step h."
function state_transition_jacobian(x::SVector{7,Float64}, u::SVector{3,Float64}, h::Float64)
    # Continuous‑time linearised dynamics:
    # δφ̇ = −hat(ωc) * δφ − δβ
    # δβ̇ = 0
    β = SVector{3,Float64}(x[5], x[6], x[7])
    ωc = u - β
    F_c = zeros(6,6)
    F_c[1:3,1:3] .= -Matrix(hat(ωc))
    F_c[1:3,4:6] .= -Matrix(I(3))
    # Discrete‑time approximation A ≈ I + F_c*h
    return I(6) + F_c * h
end

"Prediction step: propagate the state and its covariance over one time step.
   Qproc is the 6×6 process noise covariance on the error state."
function mekf_predict(x::SVector{7,Float64}, P::Matrix{Float64}, u::SVector{3,Float64}, h::Float64, Qproc::Matrix{Float64})
    xpred = state_prediction(x, u, h)
    A = state_transition_jacobian(x, u, h)
    Ppred = A * P * transpose(A) + Qproc
    return xpred, Ppred
end

"Form the measurement innovation, Jacobian and covariance for vector sensors.
   Each vector sensor measures a known inertial vector r_N expressed in
   body coordinates.  The predicted measurement is y_pred = Q(q_pred)' * r_N.
   The innovation is y_meas − y_pred.  The Jacobian with respect to φ is
   −hat(y_pred) and zeros with respect to bias."
function build_vector_measurement(q_pred::SVector{4,Float64}, vec_meas::Vector{SVector{3,Float64}}, vec_refs::Vector{SVector{3,Float64}}, R_vec::Vector{Matrix{Float64}})
    n = length(vec_meas)
    z = zeros(3*n)
    C = zeros(3*n, 6)
    W = zeros(3*n, 3*n)
    for i in 1:n
        r_N = vec_refs[i]
        y_pred = transpose(Q(q_pred)) * r_N
        z[(3*(i-1)+1):(3*i)] .= vec_meas[i] - y_pred
        C[(3*(i-1)+1):(3*i), 1:3] .= -Matrix(hat(y_pred))
        C[(3*(i-1)+1):(3*i), 4:6] .= 0.0
        W[(3*(i-1)+1):(3*i), (3*(i-1)+1):(3*i)] .= R_vec[i]
    end
    return z, C, W
end

"Form the measurement innovation, Jacobian and covariance for star trackers.
   Each star tracker returns a quaternion measurement q_meas that equals
   q_true ⊗ δq, where δq ≈ expq(0.5*δφ).  We form the error quaternion
   q_err = conj(q_pred) ⊗ q_meas, extract its rotation vector φ = logq(q_err),
   and treat this as the measurement.  The predicted measurement is 0.
   The Jacobian is identity in φ and zero in bias."
function build_star_measurement(q_pred::SVector{4,Float64}, star_meas::Vector{SVector{4,Float64}}, R_st::Vector{Matrix{Float64}})
    m = length(star_meas)
    z = zeros(3*m)
    C = zeros(3*m, 6)
    W = zeros(3*m, 3*m)
    for i in 1:m
        q_m = star_meas[i]
        # Error quaternion: δq = conj(q_pred) ⊗ q_meas
        q_err = L(qconj(q_pred)) * q_m
        φ = logq(q_err)
        z[(3*(i-1)+1):(3*i)] .= φ
        C[(3*(i-1)+1):(3*i), 1:3] .= I(3)
        C[(3*(i-1)+1):(3*i), 4:6] .= 0.0
        W[(3*(i-1)+1):(3*i), (3*(i-1)+1):(3*i)] .= R_st[i]
    end
    return z, C, W
end

"Update step: incorporate the measurement innovations z with Jacobian C and
   measurement noise covariance W.  This function updates the estimate
   (x_pred, P_pred) to (x_upd, P_upd)."
function mekf_update(xpred::SVector{7,Float64}, Ppred::Matrix{Float64}, z::Vector{Float64}, C::Matrix{Float64}, W::Matrix{Float64})
    S = C * Ppred * transpose(C) + W
    K = Ppred * transpose(C) * inv(S)
    δx = K * z
    φ = δx[1:3]
    δβ = δx[4:6]
    q_pred = SVector{4,Float64}(xpred[1], xpred[2], xpred[3], xpred[4])
    β_pred = SVector{3,Float64}(xpred[5], xpred[6], xpred[7])
    Δq = expq(φ)
    q_upd = L(q_pred) * Δq
    q_upd /= norm(q_upd)
    β_upd = β_pred + δβ
    x_upd = SVector{7,Float64}(q_upd[1], q_upd[2], q_upd[3], q_upd[4], β_upd[1], β_upd[2], β_upd[3])
    I6 = I(6)
    P_upd = (I6 - K * C) * Ppred * transpose(I6 - K * C) + K * W * transpose(K)
    return x_upd, P_upd
end

##############################
# Example simulation
##############################

function run_mekf_model_notation_example()
    # Seed RNG for reproducibility
    Random.seed!(20240409)

    # Define the true inertia and spin parameters (same as hw2)
    I_body = @SMatrix [0.249  0.0    0.0;
                       0.0    0.187  0.02;
                       0.0    0.02   0.148]
    e_panel = @SVector [0.0, 1.0, 0.0]
    omega_des = rpm2rad(10.0)
    w_des = omega_des .* e_panel

    # Perturb inertia and compute required rotor momentum
    I_pert = perturb_inertia(I_body; sig_eig=0.03, sig_rot=deg2rad(2.0))
    H_r, _ = required_rotor_momentum(I_pert, e_panel; ratio=1.2, omega=omega_des)

    # Simulate the true attitude dynamics using integrate_gyrostat!
    q0_true = @SVector [1.0, 0.0, 0.0, 0.0]
    delta_w = 0.01 * omega_des
    xi = randn(3); xi -= dot(xi, e_panel) * e_panel; xi /= norm(xi)
    w0_true = w_des + delta_w * xi
    sol_true = integrate_gyrostat!(I_pert, w0_true, q0_true, H_r, e_panel;
                                  tspan=(0.0, 400.0), reltol=1e-9, abstol=1e-12)

    # Measurement noise parameters from report2.ipynb
    sigma_st_cross = 6.0 * (π/180) / 3600    # 6 arcsec → rad
    sigma_st_bore  = 40.0 * (π/180) / 3600   # 40 arcsec → rad
    R_st_covs = [Diagonal([sigma_st_cross^2, sigma_st_cross^2, sigma_st_bore^2]) for _ in 1:2]
    sigma_sun = 0.3 * (π/180)
    sigma_mag = 1.0 * (π/180)
    R_vec_covs = [sigma_sun^2 * I(3), sigma_mag^2 * I(3)]
    # Inertial reference vectors (sun and magnetic field)
    r_sun = normalize(@SVector [1.0, 0.0, 0.0])
    r_mag = normalize(@SVector [0.0, 0.0, 1.0])

    # Extract time step and number of samples
    n = length(sol_true.t)
    h = sol_true.t[2] - sol_true.t[1]

    # True bias (for simulation) and gyro noise parameters
    β_true = @SVector [0.0, 0.0, 0.0]
    sigma_gyro_noise = 0.001   # rad/s standard deviation
    sigma_bias_walk  = 1e-5     # rad/s^1.5 random walk per √s
    # Process noise covariance on the error state
    Qproc = zeros(6,6)
    Qproc[1:3,1:3] .= sigma_gyro_noise^2 * h^2 * I(3)
    Qproc[4:6,4:6] .= sigma_bias_walk^2 * h * I(3)

    # Preallocate arrays for measurements
    gyro_meas = Vector{SVector{3,Float64}}(undef, n)
    vec_meas = [Vector{SVector{3,Float64}}(undef, n) for _ in 1:2]
    star_meas = [Vector{SVector{4,Float64}}(undef, n) for _ in 1:2]

    # Generate synthetic measurements from the truth solution
    for k in 1:n
        u = sol_true.u[k]
        qk = SVector{4,Float64}(u[1], u[2], u[3], u[4])
        ωk = SVector{3,Float64}(u[5], u[6], u[7])
        # Gyro measurement = true rate + bias + white noise
        gyro_meas[k] = ωk + β_true + sigma_gyro_noise * randn(SVector{3,Float64})
        # Vector sensors measure body directions of inertial vectors
        y_sun = transpose(Q(qk)) * r_sun
        y_mag = transpose(Q(qk)) * r_mag
        vec_meas[1][k] = noisy_vector_meas(y_sun, R_vec_covs[1])
        vec_meas[2][k] = noisy_vector_meas(y_mag, R_vec_covs[2])
        # Star tracker measurements (two sensors)
        star_meas[1][k] = noisy_star_tracker_meas(qk, R_st_covs[1])
        star_meas[2][k] = noisy_star_tracker_meas(qk, R_st_covs[2])
    end

    # Initialise the state estimate and covariance
    x_est = SVector{7,Float64}(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    P_est = 0.5 * I(6)        # fairly large initial uncertainty
    # Storage for estimates
    est_quat = Vector{SVector{4,Float64}}(undef, n)
    est_beta = Vector{SVector{3,Float64}}(undef, n)
    est_quat[1] = SVector{4,Float64}(x_est[1], x_est[2], x_est[3], x_est[4])
    est_beta[1] = SVector{3,Float64}(x_est[5], x_est[6], x_est[7])

    # Run the MEKF over the time history
    for k in 1:(n-1)
        # Prediction step
        x_pred, P_pred = mekf_predict(x_est, P_est, gyro_meas[k], h, Qproc)
        # Build vector measurements at the next time step
        vec_meas_kp1 = [vec_meas[1][k+1], vec_meas[2][k+1]]
        vec_refs = [r_sun, r_mag]
        q_pred_k = SVector{4,Float64}(x_pred[1], x_pred[2], x_pred[3], x_pred[4])
        zv, Cv, Wv = build_vector_measurement(q_pred_k, vec_meas_kp1, vec_refs, R_vec_covs)
        # Build star tracker measurements at the next step
        star_meas_kp1 = [star_meas[1][k+1], star_meas[2][k+1]]
        zst, Cst, Wst = build_star_measurement(q_pred_k, star_meas_kp1, R_st_covs)
        # Concatenate innovations, Jacobians and covariances
        z = vcat(zv, zst)
        C = vcat(Cv, Cst)
        W = zeros(length(z), length(z))
        W[1:length(zv), 1:length(zv)] = Wv
        W[(length(zv)+1):end, (length(zv)+1):end] = Wst
        # Update step
        x_est, P_est = mekf_update(x_pred, P_pred, z, C, W)
        # Store results
        est_quat[k+1] = SVector{4,Float64}(x_est[1], x_est[2], x_est[3], x_est[4])
        est_beta[k+1] = SVector{3,Float64}(x_est[5], x_est[6], x_est[7])
    end

    # Compute final attitude error relative to truth
    uf = sol_true.u[end]
    q_true_final = SVector{4,Float64}(uf[1], uf[2], uf[3], uf[4])
    q_est_final = est_quat[end]
    q_err_final = L(qconj(q_est_final)) * q_true_final
    φ_err_final = logq(q_err_final)
    println("Final attitude error magnitude (deg): ", norm(φ_err_final) * 180/π)
    println("Final bias estimate: ", est_beta[end])
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_mekf_model_notation_example()
end