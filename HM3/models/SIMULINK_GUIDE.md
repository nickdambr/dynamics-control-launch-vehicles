# HM3 — Simulink guided track

This folder hosts the block-diagram mirror of the HM3 control design. The
**MATLAB scripts are the source of truth**: gains and filters are designed and
validated in `main_task1/2/3.m`; the Simulink model `hm3_closed_loop.slx` only
*reproduces* them so the closed loop can be inspected and simulated as a block
diagram. `run_simulink_closed_loop.m` overlays the two and should give matching
traces.

> A `.slx` is a binary file and cannot be generated reliably from a text agent,
> so this model is built **interactively** by following the checklist below. The
> companion script `init_simulink_hm3.m` pre-computes every parameter the blocks
> reference, so the wiring is the only manual step.

## 0. One-time setup

```matlab
cd HM3
S = init_simulink_hm3(2);     % task = 1 | 2 | 3 ; pushes variables to base ws
```

`init_simulink_hm3` exports (among others):

| Variable | Meaning |
|----------|---------|
| `A_rigid`, `Bdelta_rigid`, `Bwind_rigid`, `C_meas_rigid`, `C_plot_rigid` | rigid plant (Task 1) |
| `A_full`, `Bdelta_full`, `Bwind_full`, `C_meas_full`, `C_plot_full` | full 6-state plant (Task 2/3) |
| `Kp_th`, `Kd_th`, `Kp_z`, `Kd_z` | controller gains |
| `tvc_num`, `tvc_den` | TVC + 20 ms delay transfer function (Eq. 3) |
| `notch_num`, `notch_den` | bending notch (Eq. 4) |
| `wind_ts` | `alpha_w(t)` timeseries for a *From Workspace* block |
| `Tstop` | suggested simulation stop time |

The measurement vector convention everywhere is
`y = [theta_m, thetadot_m, z_m, zdot_m]` and the plot vector is
`[theta, z, zdot]`.

---

## 1. Task 1 — rigid plant + PD (ideal actuator)

Create `hm3_closed_loop.slx` and add:

1. **`Plant_Rigid`** — *State-Space* block
   - A = `A_rigid`, B = `[Bdelta_rigid Bwind_rigid]`, C = `[C_meas_rigid; C_plot_rigid]`, D = `zeros(7,2)`
   - Inputs (in order): `delta`, `alpha_w`
   - Outputs (in order): `theta_m, thetadot_m, z_m, zdot_m, theta, z, zdot`
2. **`Controller_PD`** — build `u_pd` from the measurements:
   ```
   u_pd = Kp_th*(theta_ref - theta_m) - Kd_th*thetadot_m
          - Kp_z*z_m - Kd_z*zdot_m
   ```
   Use four *Gain* blocks (`Kp_th`, `Kd_th`, `Kp_z`, `Kd_z`) + *Sum* blocks, or a
   single *Gain* (matrix `[Kp_th -Kp_th -Kd_th -Kp_z -Kd_z]`, multiplication
   "Matrix(K*u)") fed by a *Mux* of `[theta_ref; theta_m; thetadot_m; z_m; zdot_m]`.
   Set `theta_ref = 0` (a *Constant*).
3. **Ideal actuator** — wire `u_pd` straight into `delta` (no TVC in Task 1).
4. **`Wind`** — *From Workspace* block, data `wind_ts`, output → `alpha_w`.
5. **Logging** — name the signals `theta_sl`, `z_sl`, `zdot_sl`, `delta_sl`
   (right-click → Properties → *Log signal data*, or *To Workspace* blocks with
   those variable names, format *Timeseries*).
6. **Solver** — variable-step `ode45`, stop time `Tstop`.

Validate:

```matlab
run_simulink_closed_loop(1);   % writes figures/task1_simulink_vs_script.png
```

---

## 2. Task 2 — full model (TVC + delay + notch)

Reuse the same file. The simplest robust approach is a **Variant Subsystem** (or
just a second sheet) so Task 1 and Task 2 coexist:

1. **`Plant_Full`** — *State-Space* with `A_full`, `[Bdelta_full Bwind_full]`,
   `[C_meas_full; C_plot_full]`, `zeros(7,2)`. Same I/O names as Task 1; this
   block now carries the bending state, and the INS outputs include the bending
   contamination (Eq. 2).
2. **`Notch_Hx`** — *Transfer Fcn* block, numerator `notch_num`, denominator
   `notch_den` (Eq. 4). Place it on `u_pd`.
3. **`TVC`** — *Transfer Fcn* block, numerator `tvc_num`, denominator `tvc_den`
   (2nd-order actuator + Pade delay, Eq. 3). Cascade after the notch.
   Chain: `u_pd → Notch_Hx → TVC → delta`.
4. Keep the PD controller and wind source unchanged.
5. Validate:
   ```matlab
   run_simulink_closed_loop(2);   % figures/task2_simulink_vs_script.png
   ```

### Coupling the professor's wind generator (`strong_wind.slx`)

`General/hw3-v3/strong_wind.slx` is a **wind generator only** (output `v_w`
[m/s]); it contains no plant/controller. To drive the closed loop with it
instead of the `From Workspace` block:

1. `load('General/hw3-v3/drywind.mat')` and the `GreensiteLPV` struct
   (`load('General/hw3-v3/GreensiteLPV_DATA.mat')`) into the base workspace —
   the generator's lookup tables need them. Stop time of that model is 140 s.
2. Add a **Model** block referencing `strong_wind.slx`; take its `v_w` output.
3. Convert wind velocity to wind angle of attack with a *Gain* `1/V`
   (`alpha_w = v_w / V`, `V = S.p.V`), then feed `alpha_w` into `Plant_Full`.
4. Optionally save the generated profile to `figures/task2_wind_profile.png`.

---

## 3. Task 3 — corner cases (±30 %)

The plant matrices depend on `mu_alpha` (A6) and `mu_c` (K1). Re-initialise the
workspace per corner and re-run:

```matlab
corners = {1.0 1.0; 0.7 1.0; 1.3 1.0; 1.0 0.7; 1.0 1.3};
for k = 1:size(corners,1)
    init_simulink_hm3(3,'mu_alpha_scale',corners{k,1},'mu_c_scale',corners{k,2});
    sim('hm3_closed_loop','StopTime',num2str(Tstop));
    % collect theta_sl / z_sl for the overlay
end
```

Because `init_simulink_hm3` rebuilds `A_full`/`Bdelta_full` with the scaled
coefficients and pushes them to the base workspace, the *State-Space* block
automatically picks up the corner without any rewiring (provided the block
references the workspace variables by name and the model is configured to
*re-evaluate* parameters between runs — the default for tunable parameters).

---

## Signal-naming contract (must match `run_simulink_closed_loop.m`)

| Signal | Source |
|--------|--------|
| `theta_sl` | `theta` plant output |
| `z_sl`     | `z` plant output |
| `zdot_sl`  | `zdot` plant output |
| `delta_sl` | actual TVC deflection `delta` |

Log them as *Timeseries* (signal logging or *To Workspace*). The overlay script
reads them from `logsout` first, then from the SimulationOutput fields.
