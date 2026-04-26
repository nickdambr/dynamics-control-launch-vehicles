---
id: T003
title: Scaffold HM2 — Powered descent & landing (direct collocation)
priority: high
created: 2026-04-26
updated: 2026-04-26
---

## Context
Homework 2 (`DCLV_MATERIALE_CORSO_26042026/Homework 2 - Powered Descent
Landing.pdf`) is the next assignment. It has two tasks:

- **Task 1 (mandatory):** fixed-duration minimum-fuel landing via direct
  collocation with trapezoidal transcription. Sensitivity analysis on `tf`
  (±5%).
- **Task 2 (optional):** Zero-Order Hold discretization, then forward-integrate
  with `ode45` and compare with the optimization output.

Problem: 2D Cartesian, point-mass, no aerodynamics, flat Earth. Glide-slope
constraint, thrust-magnitude bounds, soft pinpoint landing. Data in Table 1
of the PDF.

## Acceptance criteria
- [x] `HM2_powered_descent/` folder created (uses the new `HM<N>_<topic>` convention)
- [x] `main_task1.m` solves the trapezoidal-collocation NLP and produces:
      trajectory plot (xy), thrust profile vs time, mass vs time, glide-slope
      check, fuel sensitivity for `tf ± 5%`
- [x] Glide-slope constraint formulated as a convex (linear) constraint pair
      `±x − tan(θmax)·y ≤ 0`
- [ ] `main_task2.m` (optional) implements ZOH discretization + forward-integration
      check  — **deferred** (tracked in HM2 README roadmap)
- [x] Local `HM2_powered_descent/README.md` with problem statement, approach, key results

## Notes
- Reference: `DCLV_MATERIALE_CORSO_26042026/matlab/cvx_sled_class_2026.m`
  (YALMIP-based convex example) and `dircol_class.m`.
- Decide upfront which solver: `fmincon` (no extra deps) vs YALMIP + a conic
  solver (cleaner for the convex sub-problems but adds a dependency).
- Initial guess matters a lot for direct collocation — use a straight-line
  interpolation from `x0` to `xf` for state, and `(0, m·g)` for thrust.

## Resolution
- Chose **fmincon (sqp)** to keep the scaffold dependency-free. SCvx /
  YALMIP path documented as a follow-up in the HM2 README roadmap.
- Initial guess as planned (linear interp + hover thrust) — converged at the
  first call without warm-starting.
- Smoke-test: three sensitivity runs return `m_f` ≈ 1390-1398 kg with fuel
  growing monotonically with `tf`, which is physically consistent (longer
  hover ⇒ more gravity loss).
- **Known limitation:** the `tf = 39.9 s` run hits `MaxIterations = 500`
  with first-order optimality ≈ 5×10⁻⁴ and constraint residual ≈ 6×10⁻³.
  Acceptable for a scaffold; tracked in the HM2 README TODO list.
- Task 2 (ZOH) and SCvx are explicitly listed as next steps in the HM2
  README; they will become their own tickets when prioritized.
