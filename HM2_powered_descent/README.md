# HM2 — Powered Descent and Landing of Reusable Launch Vehicles

Direct-collocation solution to the minimum-fuel powered-descent problem for a
2D point-mass model (no aerodynamics, flat Earth). Ports the SpaceX-style
"divert-and-soft-land" planning problem to a small NLP solvable in MATLAB
without external dependencies.

> **Assignment:** see
> [`Homework 2 - Powered Descent Landing.pdf`](../DCLV_MATERIALE_CORSO_26042026/Homework%202%20-%20Powered%20Descent%20Landing.pdf)
> in the course material folder.

## Problem at a glance

| Quantity              | Value                          |
| --------------------- | ------------------------------ |
| Initial position      | (1000, 3000) m                 |
| Initial velocity      | (300, −200) m/s                |
| Initial mass          | 2000 kg                        |
| Final state           | (0, 0, 0, 0) — pinpoint, soft  |
| Flight time `tf`      | 38 s (fixed)                   |
| Thrust bounds         | 0 ≤ |T| ≤ 70 kN                |
| Glide-slope half-angle| 60°                            |
| `Isp`                 | 225 s                          |

Cost: **minimize fuel** ≡ maximize `m(tf)`.

## Approach

- **Transcription:** trapezoidal direct collocation on `N = 50` evenly spaced
  nodes. Decision vector stacks `[x, y, vx, vy, m, Tx, Ty]` per node.
- **Solver:** `fmincon` with the SQP algorithm. No external optimization
  toolbox required (CasADi / YALMIP not used at this stage — see roadmap).
- **Initial guess:** linear interpolation between the boundary states for the
  state variables; constant hover thrust `(0, m₀·g)` for the controls.
- **Glide-slope constraint** rewritten as the linear pair
  `±x − tan(θmax)·y ≤ 0`, which is convex.
- **Sensitivity sweep:** the script re-solves the problem for `tf ∈ {0.95, 1.00, 1.05} · 38 s`
  and overlays the three solutions on the same plots.

## How to run

From this folder:

```matlab
main_task1
```

Or headless:

```bash
matlab -batch "run('main_task1.m')"
```

Expected runtime: ~30 s per `tf` value on a modern laptop (three runs total).

## Files

| File              | Role                                                 |
| ----------------- | ---------------------------------------------------- |
| `main_task1.m`    | Top-level script: data, sensitivity sweep, plots.    |
|                   | `solve_trapcol` — builds and solves the NLP.         |
|                   | `trap_nonlcon` — defects + thrust bounds + glide-slope. |
|                   | `dyn_rhs` — continuous dynamics (Eq. 2-6 of PDF).    |
|                   | `plot_results` — trajectory, thrust, mass, glide-slope plots. |

## Results

| `tf` [s] | `m_f` [kg] | fuel [kg] |
| -------- | ---------- | --------- |
| 36.10    | 1406.33    | 593.67    |
| 38.00    | 1403.20    | 596.80    |
| 39.90    | 1398.82    | 601.18    |

Fuel consumption grows monotonically with `tf` over the swept window
(longer hover ⇒ more gravity losses). All three solutions respect the
glide-slope corridor and the thrust-magnitude bounds. The two shorter-`tf`
runs terminate at the `MaxIterations` cap of 1000 with first-order
optimality of order `1e-4` in non-dim units; constraint violations are
below `1e-6` in all cases.

| Trajectory (3 sensitivity runs + glide-slope corridor) | Thrust magnitude |
|:-:|:-:|
| ![Trajectory](figures/task1_trajectory.png) | ![Thrust](figures/task1_thrust_magnitude.png) |

| Mass | Glide-slope angle |
|:-:|:-:|
| ![Mass](figures/task1_mass.png) | ![Glide-slope](figures/task1_glide_slope.png) |

## Roadmap / TODO

- [x] **Task 2 (PDF Appendix A) — Nonlinear ZOH + RK4.** Multiple-shooting NLP
      with `x_{k+1} = RK4(x_k, u_k, dt)`. Validated against `ode45` forward
      integration (max node error `1.4e-8` non-dim).
- [x] **Task 2 — LTV-linearised ZOH with SCvx (fmincon/SQP).** Literal
      Appendix-A construction (augmented ODE for `Φ`, `B̂`, `ĉ`), adaptive
      trust-region ratio test, warm-started from the trapezoidal solution.
- [x] **Task 2 — SCvx with conic inner solver (YALMIP + ECOS).** Same outer
      loop; inner sub-problem cast as an SOCP and solved by ECOS. Per-step
      fidelity is an order of magnitude tighter than the fmincon path.
- [ ] **Free-time variant.** Lift the `tf`-fixed assumption and minimise fuel
      over a variable horizon — the right framing for a GFOLD-like
      formulation.
- [ ] **Lossless convexification.** Handle a non-zero lower thrust bound
      `T_min > 0` via the slack-variable change of variable
      (Açikmeşe–Ploen), preserving the SOCP structure of variant (c).
- [ ] **Tighten convergence at short `tf`.** The 36.10 s and 38.00 s
      sensitivity runs currently hit the `MaxIterations = 1000` cap with
      `OptimalityTolerance ≈ 1e-4` non-dim. Either raise the cap further,
      supply analytical gradients, or warm-start from the nominal solution.
