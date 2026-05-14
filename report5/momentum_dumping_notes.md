# Momentum Dumping Implementation Notes

## Goal

Create a rigorous report 5 scenario where reaction wheels saturate under a documented disturbance stress test, then show magnetic torque rods unloading wheel momentum while attitude remains controlled.

## Approach

- Keep the nominal report 4 environmental torque model as the source of the disturbance direction and orbit dependence.
- Apply an explicit stress-test scale factor computed from the nominal orbit-averaged torque and the desired demonstration duration.
- Use a simple Earth-centered dipole field for the magnetorquer demonstration.
- Use `m x B` magnetic torque with per-axis dipole saturation at 0.1 A m^2.
- Use external-torque feedforward in the comparison runs so wheel momentum growth reflects the disturbance balance rather than large transient PD control effort.

## Important Tuning Lesson

A first attempt with pure PD control saturated the wheels almost immediately and produced large attitude error. That was not the desired phenomenon: the wheel momentum was dominated by controller transient torques rather than accumulated environmental disturbance torque.

The current helper therefore supports `feedforward_external_torque=true`. With this enabled, the wheel torque cancels the known external disturbance plus magnetorquer torque while the PD term remains available for attitude regulation. This makes the no-dump wheel momentum grow approximately from the stressed environmental torque, and makes the dumping case demonstrate the intended unloading law.

## Current Demonstration Behavior

Using a four-orbit stress window and a stress scale computed by `estimate_stress_scale_for_saturation`, the no-dump run reaches the 50 mN m s per-axis wheel limit near the end of the window. With magnetic dumping enabled, the same stress case keeps wheel momentum below the limit while respecting the 0.1 A m^2 rod command limit.

## Limitations

- The magnetic field model is a dipole approximation, not IGRF.
- The stress-test scale is not a nominal Starling environmental claim.
- The notebook isolates momentum management by using truth-fed control and known external-torque feedforward. A later estimator-in-the-loop demonstration should retune gains and noise assumptions separately.

## MeshCat visualization

- **File:** [momentum_dumping_meshcat.jl](momentum_dumping_meshcat.jl) — `build_bus_visualizer`, `animate_quaternion_history!`, `prepare_three_meshcat_solutions` (three sims, no viewer), `run_three_meshcat_demos` (auto sequential replay), `interactive_meshcat_scenarios` (WebIO + Interact buttons/sliders; optional).
- **Dependencies:** `MeshCat`, `GeometryBasics`, `CoordinateTransformations`, `Colors`, and `Rotations` (already loaded in the report5 notebook from the report4 stack). If `using MeshCat` fails, run once: `import Pkg; Pkg.add(["MeshCat", "GeometryBasics", "CoordinateTransformations", "Colors"])`.
- **Interactive UI:** `interactive_meshcat_scenarios` needs WebIO and Interact (`Pkg.add(["WebIO", "Interact"])`), then restart IJulia. The returned `.ui` should be the last expression in the notebook cell. It calls `prepare_three_meshcat_solutions`, opens MeshCat, and builds buttons plus sliders; each button replays one trajectory from t = 0. `open_external=true` (default) also opens a browser tab. Return value: `(; ui, vis, sol_disturb, sol_wheels, sol_dump)`. **Julia 1.12:** Interact is loaded with `import` after `include`; the implementation uses `Base.invokelatest` so slider/button methods are callable (avoids “slider overflow” / world-age `MethodError`). **Jupyter:** the red “WebIO extension not detected” banner is separate from MeshCat—use `install_webio_jupyter!()` once (defined in `momentum_dumping_meshcat.jl` after include) or skip embedded UI and rely on the browser tab.
- **Run location:** execute the notebook from the `report5/` directory so `include("momentum_dumping_meshcat.jl")` resolves.
- **Quaternion convention:** attitude uses `q = [q0, q1, q2, q3]` scalar-first, matching `Rotations.QuatRotation(q0, q1, q2, q3)` as in [report1/starling_adcs.ipynb](../report1/starling_adcs.ipynb).
- **Ports:** MeshCat starts a local HTTP server; if the browser does not open, read the printed URL. If a port is in use, close other MeshCat sessions or change the default in MeshCat settings.
- **Performance:** long `T_demo` runs produce thousands of frames. Use `subsample` (e.g. 25) and `frame_pause` in `run_three_meshcat_demos` to keep playback length reasonable.
- **Playback speed:** `playback_slowdown` (default `1.0`) multiplies both the per-frame sleep and `pause_between_scenarios`. Use `playback_slowdown=3` or higher to slow the animation without retuning `frame_pause`.
- **Scenario A:** `simulate_momentum_dumping(..., wheel_torque_mode=:none)` — no wheel torque; body rotates under scaled environmental torque only (no magnetorquers in that replay).

## Troubleshooting (MeshCat / WebIO)

- **`interactive_meshcat_scenarios` errors on missing WebIO/Interact:** run `import Pkg; Pkg.add(["WebIO", "Interact"])` and restart the Jupyter kernel so IJulia reloads extensions.
- **Blank iframe but MeshCat URL works:** set `open_external=true` (default) and use the browser tab; some VS Code / Jupyter setups need the [WebIO Jupyter extension](https://juliagizmos.github.io/WebIO.jl/latest/providers/ijulia/) or a plain browser for the embedded iframe.
- **“Animation already running” warning:** wait for the current replay to finish before clicking another scenario; only one blocking `animate_quaternion_history!` loop runs at a time.
- **Second `include("momentum_dumping_meshcat.jl")` fails on `const`:** restart the kernel; top-level `const` in included files cannot be redefined in the same session.
- **`Function slider was about to overflow` or world-age `MethodError` on Julia 1.12:** `interactive_meshcat_scenarios` uses `Base.invokelatest` for all Interact calls when packages are brought in via `Core.eval` during `_ensure_interactive_packages!`. Pull the latest `momentum_dumping_meshcat.jl`. Alternatively, run `using WebIO, Interact` in a cell **before** `include("momentum_dumping_meshcat.jl")` so bindings exist in the same world as your other code.
- **WebIO “nbextension” / “jupyter-lab-provider” not found or widgets blank:** WebIO must register with the **same** Jupyter install you actually launch (conda vs pip vs Julia’s Conda.jl often diverge). After `include("momentum_dumping_meshcat.jl")`, run `install_webio_jupyter!()` once, or manually: `using WebIO, IJulia` then `WebIO.install_jupyter_nbextension()` and `WebIO.install_notebook_config()` for **classic Notebook**; for **JupyterLab**, `WebIO.install_jupyter_labextension()`. In the shell that matches your `jupyter`, run `jupyter labextension list` and confirm `@webio/jupyter-lab-provider` appears. Restart Jupyter / VS Code after installing. Official notes: [WebIO IJulia provider](https://juliagizmos.github.io/WebIO.jl/stable/providers/ijulia/).
