# MeshCat visualization for momentum dumping scenarios (report5).
#
# Load after momentum_dumping.jl and report4 stack:
#   include("momentum_dumping.jl")
#   include("momentum_dumping_meshcat.jl")
#
# Requires: MeshCat, GeometryBasics, CoordinateTransformations, Rotations, Colors

using MeshCat
using GeometryBasics
using CoordinateTransformations
using Rotations
using Colors
using LinearAlgebra
using Random
using StaticArrays

# 6U-like rectangular prism: ~30 x 10 x 10 cm, body frame origin at geometric center.
const MESHCAT_BUS_LX = 0.30
const MESHCAT_BUS_LY = 0.10
const MESHCAT_BUS_LZ = 0.10

function build_bus_visualizer(; bus_name=:bus)
    vis = Visualizer()
    lx, ly, lz = MESHCAT_BUS_LX, MESHCAT_BUS_LY, MESHCAT_BUS_LZ
    rect = HyperRectangle(
        Vec3(-lx / 2, -ly / 2, -lz / 2),
        Vec3(lx, ly, lz),
    )
    setobject!(
        vis[bus_name],
        rect,
        MeshPhongMaterial(color=RGBA(0.25, 0.55, 0.85, 0.85)),
    )
    return vis
end

# q = [q0, q1, q2, q3] scalar-first, matches Rotations.QuatRotation(w, x, y, z).
function quat_to_linear_map(q::AbstractVector{<:Real})
    return LinearMap(QuatRotation(q[1], q[2], q[3], q[4]))
end

function animate_quaternion_history!(
    vis,
    t,
    q_hist::AbstractMatrix{<:Real};
    bus_name=:bus,
    subsample::Int=1,
    frame_pause::Float64=0.02,
    playback_slowdown::Float64=1.0,
)
    subsample < 1 && throw(ArgumentError("subsample must be >= 1"))
    playback_slowdown < 0.0 && throw(ArgumentError("playback_slowdown must be >= 0"))
    dt_sleep = frame_pause * playback_slowdown
    n = size(q_hist, 2)
    for k in 1:subsample:n
        settransform!(vis[bus_name], quat_to_linear_map(@view q_hist[:, k]))
        dt_sleep > 0.0 && sleep(dt_sleep)
    end
    return nothing
end

function prepare_three_meshcat_solutions(
    cfg;
    T_demo::Float64,
    dt_ctrl::Float64,
    qd::Vector{Float64},
    stress_scale::Float64,
    g_pd::PDGains,
    rng::AbstractRNG=MersenneTwister(101),
    sigma_gyro::Float64=2e-4,
    sigma_vec::Float64=1e-3,
    professor_mekf::Bool=true,
)
    q0 = copy(qd)
    omega0 = zeros(3)
    rho0 = zeros(3)
    V_mekf = (sigma_gyro^2) * I(3)
    common = (;
        Tfinal=T_demo,
        dt=dt_ctrl,
        q0=q0,
        omega0=omega0,
        rho0=rho0,
        qd=qd,
        g=g_pd,
        sigma_gyro=sigma_gyro,
        sigma_vec=sigma_vec,
        V_mekf=V_mekf,
        control_from_truth=true,
        control_gyro_lpf_tau=2.0,
        sigma_bias_walk=1e-7,
        feedforward_external_torque=true,
        professor_mekf=professor_mekf,
    )

    sol_disturb = simulate_momentum_dumping(
        cfg;
        common...,
        rng=rng,
        disturbance_scale=stress_scale,
        enable_dumping=false,
        wheel_torque_mode=:none,
    )
    sol_wheels = simulate_momentum_dumping(
        cfg;
        common...,
        rng=rng,
        disturbance_scale=stress_scale,
        enable_dumping=false,
        wheel_torque_mode=:controlled,
    )
    sol_dump = simulate_momentum_dumping(
        cfg;
        common...,
        rng=rng,
        disturbance_scale=stress_scale,
        enable_dumping=true,
        wheel_torque_mode=:controlled,
        dump_gain=2.0e-4,
        dipole_max=0.1,
        dump_deadband=2e-3,
    )

    return (; sol_disturb, sol_wheels, sol_dump)
end

# Runs three simulations and animates sequentially in one MeshCat window:
#   1. Disturbance only (wheel_torque_mode=:none)
#   2. Reaction wheels, no dumping
#   3. Wheels + magnetorquers (momentum dumping)
function run_three_meshcat_demos(
    cfg;
    T_demo::Float64,
    dt_ctrl::Float64,
    qd::Vector{Float64},
    stress_scale::Float64,
    g_pd::PDGains,
    rng::AbstractRNG=MersenneTwister(101),
    sigma_gyro::Float64=2e-4,
    sigma_vec::Float64=1e-3,
    professor_mekf::Bool=true,
    subsample::Int=0,
    frame_pause::Float64=0.02,
    playback_slowdown::Float64=1.0,
    pause_between_scenarios::Float64=2.0,
)
    playback_slowdown < 0.0 && throw(ArgumentError("playback_slowdown must be >= 0"))

    sols = prepare_three_meshcat_solutions(
        cfg;
        T_demo=T_demo,
        dt_ctrl=dt_ctrl,
        qd=qd,
        stress_scale=stress_scale,
        g_pd=g_pd,
        rng=rng,
        sigma_gyro=sigma_gyro,
        sigma_vec=sigma_vec,
        professor_mekf=professor_mekf,
    )

    vis = build_bus_visualizer()
    open(vis)

    nframes = size(sols.sol_disturb.q_hist, 2)
    sub = subsample > 0 ? subsample : max(1, div(nframes, 200))

    pause_scenarios = pause_between_scenarios * playback_slowdown

    @info "MeshCat scenario 1/3: disturbance only (no wheel torque)"
    animate_quaternion_history!(
        vis, collect(sols.sol_disturb.t), sols.sol_disturb.q_hist;
        subsample=sub, frame_pause=frame_pause, playback_slowdown=playback_slowdown,
    )
    sleep(pause_scenarios)

    @info "MeshCat scenario 2/3: reaction wheels (saturate under stress, no magnetic dumping)"
    animate_quaternion_history!(
        vis, collect(sols.sol_wheels.t), sols.sol_wheels.q_hist;
        subsample=sub, frame_pause=frame_pause, playback_slowdown=playback_slowdown,
    )
    sleep(pause_scenarios)

    @info "MeshCat scenario 3/3: wheels + magnetorquers (momentum dumping, attitude recovery)"
    animate_quaternion_history!(
        vis, collect(sols.sol_dump.t), sols.sol_dump.q_hist;
        subsample=sub, frame_pause=frame_pause, playback_slowdown=playback_slowdown,
    )

    return (; vis, sol_disturb=sols.sol_disturb, sol_wheels=sols.sol_wheels, sol_dump=sols.sol_dump)
end
