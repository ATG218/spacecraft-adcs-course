using LinearAlgebra
using StaticArrays
using DifferentialEquations
using Rotations
using Random

rpm2rad(rpm) = rpm * 2π / 60
rad2rpm(rad) = rad * 60 / (2π)

# Perturb inertia eigenvalues by sig_eig (fractional) and eigenvectors by
# a random rotation with std sig_rot (rad)
function perturb_inertia(I_mat::SMatrix{3,3,Float64}; sig_eig=0.03, sig_rot=deg2rad(2.0))
    E = eigen(Matrix(I_mat))
    V = SMatrix{3,3,Float64,9}(E.vectors)
    D = E.values

    Dp = D .* (1 .+ sig_eig * randn(3))
    D_tilde = @SMatrix [Dp[1] 0.0 0.0; 0.0 Dp[2] 0.0; 0.0 0.0 Dp[3]]

    v = sig_rot * randn(3)
    theta = norm(v)
    I3 = SMatrix{3,3,Float64}(I(3))
    if theta > 0
        axis = v / theta
        K = @SMatrix [
             0.0      -axis[3]  axis[2]
             axis[3]   0.0     -axis[1]
            -axis[2]   axis[1]  0.0
        ]
        Rdelta = I3 + sin(theta)*K + (1 - cos(theta))*(K*K)
    else
        Rdelta = I3
    end

    V_tilde = Rdelta * V
    return V_tilde * D_tilde * V_tilde'
end

# Required rotor momentum so that J_spin_eff / J_transverse >= ratio
function required_rotor_momentum(I::SMatrix{3,3,Float64}, e::SVector{3,Float64}; ratio=1.2, omega::Float64=1.0)
    en = e / norm(e)
    J_s0 = dot(en, I * en)

    # Pick a basis vector not parallel to en for Gram-Schmidt
    abs_e = abs.(en)
    if abs_e[1] < 0.9
        w0 = SVector{3,Float64}(1.0, 0.0, 0.0)
    elseif abs_e[2] < 0.9
        w0 = SVector{3,Float64}(0.0, 1.0, 0.0)
    else
        w0 = SVector{3,Float64}(0.0, 0.0, 1.0)
    end
    v1 = normalize(w0 - dot(w0, en)*en)
    v2 = normalize(cross(en, v1))
    J_t = max(dot(v1, I*v1), dot(v2, I*v2))

    I_rotor = max(0.0, ratio*J_t - J_s0)
    return I_rotor * omega, I_rotor
end

function gyrostat_attitude_ode!(dx, x, p, t)
    q = SVector(x[1], x[2], x[3], x[4])
    ω = SVector(x[5], x[6], x[7])
    h_r = p.H_r .* p.e

    Ω = @SMatrix [
         0.0   -ω[1]  -ω[2]  -ω[3]
         ω[1]   0.0    ω[3]  -ω[2]
         ω[2]  -ω[3]   0.0    ω[1]
         ω[3]   ω[2]  -ω[1]   0.0
    ]
    dq = 0.5 .* (Ω * q)
    ω_dot = p.I \ (-cross(ω, p.I*ω + h_r))

    dx[1] = dq[1]; dx[2] = dq[2]; dx[3] = dq[3]; dx[4] = dq[4]
    dx[5] = ω_dot[1]; dx[6] = ω_dot[2]; dx[7] = ω_dot[3]
end

function integrate_gyrostat!(I_mat, ω0, q0, H_r, e; tspan=(0.0, 1000.0), reltol=1e-9, abstol=1e-12)
    x0 = [q0/norm(q0); ω0]
    sol = solve(ODEProblem(gyrostat_attitude_ode!, x0, tspan, (; I=I_mat, H_r=H_r, e=e)),
                Vern9(); reltol=reltol, abstol=abstol)
    for u in sol.u
        u[1:4] ./= norm(u[1:4])
    end
    return sol
end

function orbit_gyrostat_ode!(dx, x, p, t)
    r = SVector(x[1], x[2], x[3])
    v = SVector(x[4], x[5], x[6])
    q = SVector(x[7], x[8], x[9], x[10])
    ω = SVector(x[11], x[12], x[13])

    mu = 398600.4418; RE = 6378.1363; J2 = 1.08262668e-3
    rnorm = norm(r)
    a_tb = -mu / rnorm^3 * r
    z2 = r[3]^2 / rnorm^2
    a_J2 = (3*J2*mu*RE^2 / (2*rnorm^5)) .* SVector(
        r[1]*(5*z2 - 1), r[2]*(5*z2 - 1), r[3]*(5*z2 - 3))

    a_drag = SVector(0.0, 0.0, 0.0)
    if p.drag
        r_m = 1000*r; v_m = 1000*v
        h = norm(r_m) - RE*1000
        rho = 3.614e-13 * exp(-(h - 700e3) / 88.667e3)
        v_rel = v_m - cross(SVector(0.0, 0.0, 7.2921159e-5), r_m)
        a_drag = (-0.5 * rho * (2.2*0.03/12.0) * norm(v_rel) / 1000) .* v_rel
    end

    Ω = @SMatrix [
         0.0   -ω[1]  -ω[2]  -ω[3]
         ω[1]   0.0    ω[3]  -ω[2]
         ω[2]  -ω[3]   0.0    ω[1]
         ω[3]   ω[2]  -ω[1]   0.0
    ]
    dq = 0.5 .* (Ω * q)
    h_r = p.H_r .* p.e
    ω_dot = p.I \ (-cross(ω, p.I*ω + h_r))

    dx[1]=v[1];  dx[2]=v[2];  dx[3]=v[3]
    dv = a_tb + a_J2 + a_drag
    dx[4]=dv[1]; dx[5]=dv[2]; dx[6]=dv[3]
    dx[7]=dq[1]; dx[8]=dq[2]; dx[9]=dq[3]; dx[10]=dq[4]
    dx[11]=ω_dot[1]; dx[12]=ω_dot[2]; dx[13]=ω_dot[3]
end

function integrate_orbit_gyrostat!(x0, I, H_r, e; tspan=(0.0, 3600.0), drag=false, reltol=1e-9, abstol=1e-12)
    sol = solve(ODEProblem(orbit_gyrostat_ode!, x0, tspan, (; I=I, H_r=H_r, e=e, drag=drag)),
                Vern9(); reltol=reltol, abstol=abstol)
    for u in sol.u
        u[7:10] ./= norm(u[7:10])
    end
    return sol
end

function quat_to_dcm(q::SVector{4,Float64})
    q = q / norm(q)
    q0, q1, q2, q3 = q
    return @SMatrix [
        1-2*(q2^2+q3^2)   2*(q1*q2-q0*q3)  2*(q1*q3+q0*q2)
        2*(q1*q2+q0*q3)   1-2*(q1^2+q3^2)  2*(q2*q3-q0*q1)
        2*(q1*q3-q0*q2)   2*(q2*q3+q0*q1)  1-2*(q1^2+q2^2)
    ]
end

function pointing_error(sol, e_panel::SVector{3,Float64}, sun_vec::SVector{3,Float64})
    en = e_panel / norm(e_panel)
    s  = sun_vec  / norm(sun_vec)
    return [acos(clamp(dot(quat_to_dcm(SVector(u[7],u[8],u[9],u[10])) * en, s), -1.0, 1.0)) * (180/π)
            for u in sol.u]
end
