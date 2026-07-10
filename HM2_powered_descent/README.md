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
- **Solver:** `fmincon` with the SQP algorithm for Task 1 and the Task-2
  variants (a)/(b); the conic Task-2 variants (c)/(d) use YALMIP + ECOS and
  are skipped gracefully if those packages are absent.
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

Expected runtime: ~1–1.5 min per `tf` value on a modern laptop (three sweep
runs plus a three-point grid-convergence study).

## Files

| File              | Role                                                 |
| ----------------- | ---------------------------------------------------- |
| `main_task1.m`    | Top-level script: data, sensitivity sweep, grid study, plots. |
|                   | `solve_trapcol` — builds and solves the NLP.         |
|                   | `trap_nonlcon` — defects + thrust bounds + glide-slope. |
|                   | `dyn_rhs` — continuous dynamics (Eq. 2-6 of PDF).    |
|                   | `diagnostics` — switching times, glide-slope margin, KKT activity. |
|                   | `fwd_integrate_pwl` / `node_err` — ode45 replay fidelity metric. |
|                   | `plot_results` — trajectory, thrust, mass, glide-slope plots. |

## Results

| `tf` [s] | `m_f` [kg] | fuel [kg] |
| -------- | ---------- | --------- |
| 36.10    | 1406.33    | 593.67    |
| 38.00    | 1403.20    | 596.80    |
| 39.90    | 1398.82    | 601.18    |

Fuel consumption grows monotonically with `tf` over the swept window
(longer coast ⇒ more gravity losses). The thrust profile is the classic
**max–coast–max** (bang-off-bang): for the nominal run the burns switch at
`t ≈ 14.0 s` and `t ≈ 33.1 s`, and the `tf` variation is absorbed almost
entirely by the coast arc. The KKT multipliers confirm that the upper
thrust bound is the **only active path constraint** (glide-slope multipliers
are numerically zero; minimum corridor margin 1.0–2.1° across the sweep).
The two shorter-`tf` runs terminate at the `MaxIterations` cap of 1000 with
first-order optimality of order `1e-3`–`1e-4` in non-dim units; constraint
violations are below `1e-6` in all cases.

A grid-convergence study at the nominal `tf` (N = 25/50/100) shows the
ode45-replay node error dropping by ×4 per mesh halving — the expected
`O(Δt²)` of the trapezoidal rule — while `m_f` fluctuates by less than
±0.8 kg. Replaying the optimized controls open-loop through ode45 lands
the trapezoidal solution 4.3 m from the pad at 0.11 m/s; the ZOH (RK4)
variant replays to micrometre accuracy, the LTV SCvx variants land within
centimetres to decimetres, and the GFOLD log-mass variant reaches the
integrator floor (see the report for the full tables).

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
- [x] **Task 2 — GFOLD log-mass ZOH (exact LTI) with SCvx.** Change of
      variables `z = ln m`, `u = T/m` (Açıkmese & Blackmore): the dynamics
      become exactly LTI and are discretised by a single matrix exponential.
      Self-starting (no trapezoidal warm start), converges in 3 SCvx
      iterations (~5 s wall time) and replays to the integrator floor
      (`7.3e-12` non-dim node error).
- [ ] **Lossless convexification → single SOCP.** Log-mass change
      of variables + slack `Γ`: the whole OCP becomes one convex SOCP — no
      SCvx loop, no warm start, sub-second ECOS solve with a global-optimality
      certificate. Also covers `T_min > 0` (quadratic lower bound) and the
      free-`tf` variant (1-D search over the SOCP). Full plan:
      [`tickets/open/T006_hm2-lossless-socp.md`](../tickets/open/T006_hm2-lossless-socp.md).
- [ ] **Tighten convergence at short `tf`.** The 36.10 s and 38.00 s
      sensitivity runs hit the `MaxIterations = 1000` cap with first-order
      optimality stalled at `1e-3`–`1e-4` non-dim. Root cause: the coast arc
      sits exactly at the nonsmooth point `T = 0` of `‖T‖` under
      finite-difference gradients. Fixes, in order of elegance: the
      lossless-convexification slack `Γ ≥ ‖T‖` (removes the norm from the
      dynamics), analytical gradients with a smoothed norm, or continuation
      (warm-start each sweep run from the nominal solution).
