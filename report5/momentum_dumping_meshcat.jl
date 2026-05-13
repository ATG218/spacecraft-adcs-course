# MeshCat visualization for momentum dumping scenarios (report5).
#
# Load after momentum_dumping.jl and report4 stack:
#   include("momentum_dumping.jl")
#   include("momentum_dumping_meshcat.jl")
#
# Requires: MeshCat, GeometryBasics, CoordinateTransformations, Rotations, Colors
#
# Optional (interactive UI only): WebIO, Interact — see `interactive_meshcat_scenarios`.

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

const _MESHCAT_INTERACTIVE_PKG_READY = Ref(false)
const _MESHCAT_ANIM_LOCK = ReentrantLock()

"""
    build_bus_visualizer(; bus_name=:bus)

Create a MeshCat `Visualizer`, attach a translucent rectangular prism at `vis[bus_name]`.
Returns `vis`.
"""
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

"""
Quaternion convention matches report1 / attitude_control: `q = [q0, q1, q2, q3]` scalar-first,
same as `Rotations.QuatRotation(w, x, y, z)`.
"""
function quat_to_linear_map(q::AbstractVector{<:Real})
    return LinearMap(QuatRotation(q[1], q[2], q[3], q[4]))
end

"""
    animate_quaternion_history!(vis, t, q_hist; bus_name=:bus, subsample=1, frame_pause=0.02, playback_slowdown=1.0)

Replay columns of `q_hist` (4 x N) on `vis[bus_name]`. `subsample > 1` skips frames for long runs.

**Playback speed:** `playback_slowdown` multiplies the sleep after each frame (`sleep(frame_pause * playback_slowdown)`).
Use values `> 1.0` to slow down (e.g. `4.0` ≈ 4× slower). `0.0` disables pausing between frames (fastest scrub).
"""
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

"""
    prepare_three_meshcat_solutions(cfg; T_demo, dt_ctrl, qd, stress_scale, g_pd, ...)

Run the three momentum-dumping comparison simulations used for MeshCat (no viewer).

Returns `(; sol_disturb, sol_wheels, sol_dump)`:

1. **Disturbance only** — `wheel_torque_mode=:none`, no dumping.
2. **Reaction wheels** — `wheel_torque_mode=:controlled`, no dumping.
3. **Magnetorquers** — wheels + magnetic dumping enabled.
"""
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

"""
    run_three_meshcat_demos(cfg; T_demo, dt_ctrl, qd, stress_scale, g_pd, ..., playback_slowdown=1.0)

Run three simulations and animate sequentially in one MeshCat window:

1. **Disturbance only** — scaled environmental torque, no reaction-wheel torque (`wheel_torque_mode=:none`).
2. **Wheels saturated** — same stress, wheels with feedforward, no magnetorquers.
3. **Magnetorquers** — same stress, wheels + magnetic momentum dumping.

Opens the browser viewer (`open(vis)`). Run the notebook from `report5/` so includes resolve.

**Playback:** `frame_pause` sets the base delay per displayed frame; `playback_slowdown > 1` stretches
that delay and the pause between scenarios by the same factor (e.g. `playback_slowdown=4` for ~4× slower).
"""
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

function _ensure_interactive_packages!()
    _MESHCAT_INTERACTIVE_PKG_READY[] && return nothing
    try
        Core.eval(Main, :(import Interact, WebIO))
    catch
        error(
            "Interactive MeshCat needs WebIO and Interact. In Julia run:\n" *
            "  import Pkg; Pkg.add([\"WebIO\", \"Interact\"])\n" *
            "Then restart the kernel and include this file again.",
        )
    end
    _MESHCAT_INTERACTIVE_PKG_READY[] = true
    return nothing
end

# Julia 1.12+: if Interact is loaded via `Core.eval` after this file is included, call Interact
# APIs through `invokelatest` so methods from the newer "world" are reachable (otherwise Widgets
# can error with "slider ... overflow" / MethodError world age).
function _ix_interact()
    return Base.invokelatest(getproperty, Main, Symbol("Interact"))
end

function _ix_call(sym::Symbol, args...; kwargs...)
    ix = _ix_interact()
    f = Base.invokelatest(getproperty, ix, sym)
    return Base.invokelatest(f, args...; kwargs...)
end

# Install WebIO into Jupyter (see momentum_dumping_notes.md). Optional if you only use MeshCat in a browser tab.
function install_webio_jupyter!()
    _ensure_interactive_packages!()
    wx = Base.invokelatest(getproperty, Main, Symbol("WebIO"))
    Base.invokelatest(getproperty(wx, :install_jupyter_nbextension))()
    Base.invokelatest(getproperty(wx, :install_notebook_config))()
    @info "WebIO nbextension installed; restart Jupyter / VS Code, then re-open the notebook."
    return nothing
end

function _meshcat_try_animate!(f)
    if !trylock(_MESHCAT_ANIM_LOCK)
        @warn "MeshCat animation already running; wait for it to finish before clicking another scenario."
        return nothing
    end
    try
        f()
    finally
        unlock(_MESHCAT_ANIM_LOCK)
    end
    return nothing
end

# Interactive MeshCat (WebIO + Interact). Docs: momentum_dumping_notes.md
function interactive_meshcat_scenarios(
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
    default_playback_slowdown::Float64=6.0,
    open_external::Bool=true,
    bus_name=:bus,
)
    _ensure_interactive_packages!()

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

    vis = build_bus_visualizer(; bus_name=bus_name)
    if open_external
        open(vis)
    end
    disp = render(vis)

    nframes = size(sols.sol_disturb.q_hist, 2)
    sub_hi = max(1, min(200, nframes))
    sub_default = subsample > 0 ? subsample : max(1, div(nframes, 200))
    sub_default = clamp(sub_default, 1, sub_hi)

    stp = 0.5
    slow_val = clamp(round(Int, default_playback_slowdown / stp) * stp, stp, 40.0)
    slowdown_w = _ix_call(:slider, stp:stp:40.0; label="playback slowdown", value=slow_val)
    subsample_w = _ix_call(:slider, 1:sub_hi; label="subsample (every Nth frame)", value=sub_default)

    btn1 = _ix_call(:button, "1 Disturbance only (no wheel torque)")
    btn2 = _ix_call(:button, "2 Reaction wheels, no magnetic dumping")
    btn3 = _ix_call(:button, "3 Wheels + magnetorquers (dumping)")

    _read_slowdown() = Float64(_ix_call(:observe, slowdown_w)[])
    _read_sub() = Int(round(_ix_call(:observe, subsample_w)[]))

    _ix_call(
        :on,
        (_) -> _meshcat_try_animate!() do
            @info "MeshCat replay: disturbance only"
            animate_quaternion_history!(
                vis, collect(sols.sol_disturb.t), sols.sol_disturb.q_hist;
                bus_name=bus_name,
                subsample=_read_sub(),
                frame_pause=frame_pause,
                playback_slowdown=_read_slowdown(),
            )
        end,
        _ix_call(:observe, btn1),
    )
    _ix_call(
        :on,
        (_) -> _meshcat_try_animate!() do
            @info "MeshCat replay: reaction wheels (no dumping)"
            animate_quaternion_history!(
                vis, collect(sols.sol_wheels.t), sols.sol_wheels.q_hist;
                bus_name=bus_name,
                subsample=_read_sub(),
                frame_pause=frame_pause,
                playback_slowdown=_read_slowdown(),
            )
        end,
        _ix_call(:observe, btn2),
    )
    _ix_call(
        :on,
        (_) -> _meshcat_try_animate!() do
            @info "MeshCat replay: wheels + magnetorquers"
            animate_quaternion_history!(
                vis, collect(sols.sol_dump.t), sols.sol_dump.q_hist;
                bus_name=bus_name,
                subsample=_read_sub(),
                frame_pause=frame_pause,
                playback_slowdown=_read_slowdown(),
            )
        end,
        _ix_call(:observe, btn3),
    )

    controls = _ix_call(
        :vbox,
        _ix_call(:hbox, btn1),
        _ix_call(:hbox, btn2),
        _ix_call(:hbox, btn3),
        slowdown_w,
        subsample_w,
    )
    ui = _ix_call(:vbox, controls, disp)

    return (;
        ui,
        vis,
        sol_disturb=sols.sol_disturb,
        sol_wheels=sols.sol_wheels,
        sol_dump=sols.sol_dump,
    )
end
