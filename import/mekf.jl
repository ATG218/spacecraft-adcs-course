# MEKF: 7-state (q + gyro bias β), 6×6 error covariance.
# Self-contained — no external includes needed.
# Used by attitude_control_starling_impl.jl when professor_mekf=true.
module MEKF

using LinearAlgebra
using StaticArrays

# H picks out the vector part of a quaternion: [scalar; vector] → vector
const H = @SMatrix [0.0 0.0 0.0;
                    1.0 0.0 0.0;
                    0.0 1.0 0.0;
                    0.0 0.0 1.0]

@inline function hat(v::SVector{3,Float64})::SMatrix{3,3,Float64,9}
    return @SMatrix [  0.0   -v[3]   v[2];
                       v[3]   0.0   -v[1];
                      -v[2]   v[1]   0.0]
end

function L(q::SVector{4,Float64})::SMatrix{4,4,Float64,16}
    q0 = q[1]; qv = @SVector [q[2], q[3], q[4]]
    return @SMatrix [ q0      -qv[1]  -qv[2]  -qv[3];
                      qv[1]    q0      qv[3]  -qv[2];
                      qv[2]   -qv[3]   q0      qv[1];
                      qv[3]    qv[2]  -qv[1]   q0   ]
end

function R(p::SVector{4,Float64})::SMatrix{4,4,Float64,16}
    p0 = p[1]; pv = @SVector [p[2], p[3], p[4]]
    return @SMatrix [ p0      -pv[1]  -pv[2]  -pv[3];
                      pv[1]    p0     -pv[3]   pv[2];
                      pv[2]    pv[3]   p0     -pv[1];
                      pv[3]   -pv[2]   pv[1]   p0   ]
end

# Rotation matrix body→inertial: r_N = Q(q) * r_B
function Q(q::SVector{4,Float64})::SMatrix{3,3,Float64,9}
    return transpose(H) * (transpose(R(q)) * L(q)) * H
end

function expq(phi::SVector{3,Float64})::SVector{4,Float64}
    θ = norm(phi)
    θ < 1e-12 && return @SVector [1.0, 0.0, 0.0, 0.0]
    vv = phi * sinc(θ / π)
    return SVector{4,Float64}(cos(θ), vv[1], vv[2], vv[3])
end

function expq(phi::AbstractVector{<:Real})::SVector{4,Float64}
    return expq(SVector{3,Float64}(phi[1], phi[2], phi[3]))
end

function state_prediction(x::SVector{7,Float64}, u::SVector{3,Float64}, h::Float64)
    q = SVector{4,Float64}(x[1], x[2], x[3], x[4])
    β = SVector{3,Float64}(x[5], x[6], x[7])
    ωc = u - β
    Δq = expq(0.5 * h * ωc)
    qn = L(q) * Δq
    qn /= norm(qn)
    return SVector{7,Float64}(qn[1], qn[2], qn[3], qn[4], β[1], β[2], β[3])
end

function state_transition_jacobian(x::SVector{7,Float64}, u::SVector{3,Float64}, h::Float64)
    β = SVector{3,Float64}(x[5], x[6], x[7])
    ωc = u - β
    F_c = zeros(6, 6)
    F_c[1:3, 1:3] .= -Matrix(hat(ωc))
    F_c[1:3, 4:6] .= -Matrix(I(3))
    return I(6) + F_c * h
end

function mekf_predict(x::SVector{7,Float64}, P::Matrix{Float64}, u::SVector{3,Float64}, h::Float64, Qproc::Matrix{Float64})
    xpred = state_prediction(x, u, h)
    A = state_transition_jacobian(x, u, h)
    Ppred = A * P * transpose(A) + Qproc
    return xpred, Ppred
end

function build_vector_measurement(q_pred::SVector{4,Float64}, vec_meas::Vector{SVector{3,Float64}}, vec_refs::Vector{SVector{3,Float64}}, R_vec::Vector{Matrix{Float64}})
    n = length(vec_meas)
    z = zeros(3 * n)
    C = zeros(3 * n, 6)
    W = zeros(3 * n, 3 * n)
    for i in 1:n
        r_N = vec_refs[i]
        y_pred = transpose(Q(q_pred)) * r_N
        z[(3*(i-1)+1):(3*i)] .= vec_meas[i] - y_pred
        C[(3*(i-1)+1):(3*i), 1:3] .= -Matrix(hat(y_pred))
        W[(3*(i-1)+1):(3*i), (3*(i-1)+1):(3*i)] .= R_vec[i]
    end
    return z, C, W
end

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

end # module MEKF
