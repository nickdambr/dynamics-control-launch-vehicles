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

`init_simulink_hm3` exports to the base workspace:

| Variable | Meaning |
|----------|---------|
| `A_rigid` | rigid plant A matrix (4×4) |
| `Bdelta_rigid`, `Bwind_rigid` | rigid plant input columns for `delta` and `alpha_w` (4×1 each) |
| `C_meas_rigid` | rigid plant measurement output matrix (4×4), rows = `[theta_m, thetadot_m, z_m, zdot_m]` |
| `C_plot_rigid` | rigid plant plot output matrix (3×4), rows = `[theta, z, zdot]` |
| `A_full` | full 6-state plant A matrix (6×6) |
| `Bdelta_full`, `Bwind_full` | full plant input columns (6×1 each) |
| `C_meas_full` | full plant measurement matrix (4×6), includes INS bending contamination |
| `C_plot_full` | full plant plot matrix (3×6), true (uncontaminated) states |
| `Kp_th`, `Kd_th` | pitch PD gains (tuned on the **rigid plant** for all tasks) |
| `Kp_z`, `Kd_z` | lateral drift gains (fixed at −1×10⁻³) |
| `tvc_num`, `tvc_den` | TVC 2nd-order + 3rd-order Padé delay (row vectors for Transfer Fcn block) |
| `notch_num`, `notch_den` | bending notch — **minimum-phase variant** (`numSign=+1`, Eq. 4) |
| `wind_ts` | `alpha_w(t)` as a `timeseries` for a *From Workspace* block |
| `Tstop` | suggested simulation stop time (end of wind profile) |
| `p` | full parameter struct (use `p.V` for flight velocity, etc.) |

The measurement vector convention everywhere is
`y_meas = [theta_m, thetadot_m, z_m, zdot_m]` and the plot vector is
`[theta, z, zdot]` (true states).

Optional arguments (same name/value pairs as `init_simulink_hm3`):

```matlab
S = init_simulink_hm3(3, 'mu_alpha_scale', 0.7, 'mu_c_scale', 1.3);
S = init_simulink_hm3(2, 'severity', 'severe');   % default
```

---

## 1. Task 1 — rigid plant + PD (ideal actuator)

Create `hm3_closed_loop.slx` (Save As from a blank model) and add:

1. **`Plant_Rigid`** — *State-Space* block
   - A = `A_rigid`
   - B = `[Bdelta_rigid, Bwind_rigid]`
   - C = `[C_meas_rigid; C_plot_rigid]`
   - D = `zeros(7, 2)`
   - Inputs (in order): `delta`, `alpha_w`
   - Outputs (in order): `theta_m, thetadot_m, z_m, zdot_m, theta, z, zdot`

2. **`Controller_PD`** — computes `u_pd` from measurements:
   ```
   u_pd = Kp_th*(theta_ref - theta_m) - Kd_th*thetadot_m
          - Kp_z*z_m - Kd_z*zdot_m
   ```
   Simplest implementation: one *Gain* block with row-vector gain
   `[Kp_th, -Kp_th, -Kd_th, -Kp_z, -Kd_z]` (multiplication "Matrix(K*u)"),
   fed by a *Mux* of `[theta_ref; theta_m; thetadot_m; z_m; zdot_m]`.
   Set `theta_ref = 0` via a *Constant* block.

3. **Ideal actuator** — wire `u_pd` directly to `delta` input of `Plant_Rigid`
   (no TVC block in Task 1).

4. **`Wind`** — *From Workspace* block
   - Data: `wind_ts`
   - Output after final data value: *Extrapolation* (hold last)
   - Output → `alpha_w` input of `Plant_Rigid`

5. **Logging** — right-click each signal → *Properties* → enable *Log signal
   data* and set the name; or use *To Workspace* blocks (format *Timeseries*):

   | Signal | Source | Logged name |
   |--------|--------|-------------|
   | `theta` | plant output 5 | `theta_sl` |
   | `z` | plant output 6 | `z_sl` |
   | `zdot` | plant output 7 | `zdot_sl` |
   | `delta` | wire between controller and plant | `delta_sl` |

6. **Solver** — Simulation → Model Configuration Parameters:
   - Solver: variable-step `ode45`
   - Stop time: `Tstop`

Validate:

```matlab
run_simulink_closed_loop(1);   % writes figures/task1_simulink_vs_script.png
```

---

## 2. Task 2 — full plant + TVC + delay + notch

Reuse the same file. The simplest robust approach is a **Variant Subsystem** (or
a second sheet/canvas area) so Task 1 and Task 2 coexist:

1. **`Plant_Full`** — *State-Space* block
   - A = `A_full`
   - B = `[Bdelta_full, Bwind_full]`
   - C = `[C_meas_full; C_plot_full]`
   - D = `zeros(7, 2)`
   - Same I/O names as `Plant_Rigid`; this block carries the two bending states
     (`eta`, `etadot`) and the INS outputs include bending contamination (Eq. 2).

2. **`Notch_Hx`** — *Transfer Fcn* block (minimum-phase notch, Eq. 4)
   - Numerator: `notch_num`
   - Denominator: `notch_den`
   - Place on the `u_pd` signal, before the TVC.

3. **`TVC`** — *Transfer Fcn* block (2nd-order actuator + Padé delay, Eq. 3)
   - Numerator: `tvc_num`
   - Denominator: `tvc_den`
   - Chain: `u_pd → Notch_Hx → TVC → delta`

4. Keep the PD controller (`Controller_PD`) and wind source unchanged.

5. Validate:
   ```matlab
   run_simulink_closed_loop(2);   % figures/task2_simulink_vs_script.png
   ```

### Coupling the professor's wind generator (`strong_wind.slx`)

`General/hw3-v3/strong_wind.slx` is a wind generator only (output `v_w` [m/s]);
it contains no plant or controller. To drive the closed loop with it instead of
the *From Workspace* block:

1. Load its lookup-table data into the base workspace:
   ```matlab
   load('General/hw3-v3/drywind.mat')
   load('General/hw3-v3/GreensiteLPV_DATA.mat')   % GreensiteLPV struct
   ```
   The model stop time in `strong_wind.slx` is 140 s.

2. Add a *Model* block referencing `strong_wind.slx`; connect its `v_w` output.

3. Convert wind velocity to angle of attack with a *Gain* `1/p.V`
   (`alpha_w = v_w / V`) and feed the result into the `alpha_w` input of
   `Plant_Full`.

4. Optionally save the generated profile to `figures/task2_wind_profile.png`.

---

## 3. Task 3 — corner cases (±30 %)

The plant matrices depend on `p.mu_alpha` (coefficient A6) and `p.K1`. Re-run
`init_simulink_hm3` with corner scales and re-simulate for each corner. The
*State-Space* blocks pick up the updated workspace variables automatically
(Simulink re-evaluates tunable parameters between runs by default).

```matlab
corners = {1.0 1.0; 0.7 1.0; 1.3 1.0; 1.0 0.7; 1.0 1.3};
for k = 1:size(corners,1)
    init_simulink_hm3(3, 'mu_alpha_scale', corners{k,1}, 'mu_c_scale', corners{k,2});
    so = sim('hm3_closed_loop', 'StopTime', num2str(Tstop));
    % collect theta_sl / z_sl from so for the overlay plot
end
```

---

## Signal-naming contract (must match `run_simulink_closed_loop.m`)

| Signal | Source | Logged name |
|--------|--------|-------------|
| `theta` | plant output 5 (`C_plot` row 1) | `theta_sl` |
| `z` | plant output 6 (`C_plot` row 2) | `z_sl` |
| `zdot` | plant output 7 (`C_plot` row 3) | `zdot_sl` |
| `delta` | actual TVC deflection (after TVC block, or `u_pd` in Task 1) | `delta_sl` |

Log them as *Timeseries* (signal logging or *To Workspace*). `run_simulink_closed_loop`
reads them from `logsout` first, then from the `SimulationOutput` struct fields.
