---
id: T006
title: HM2 variant (d) — lossless convexification, single-SOCP GFOLD-style solve
priority: medium
created: 2026-06-12
updated: 2026-06-12
---

## Context

HM2 Task 2 currently has four transcriptions; the best one (SCvx with
YALMIP+ECOS inner SOCP) still needs the ~90 s trapezoidal warm start plus a
15-iteration outer loop (~29 s). For this exact problem class (flat Earth, no
aerodynamics) the industry answer — GFOLD / Açikmeşe-Ploen — removes the outer
loop entirely: a change of variables plus a lossless slack relaxation makes the
whole OCP **one convex SOCP**, solvable by ECOS in well under a second, with a
global-optimality certificate. This is also the natural vehicle for two items
already on the HM2 roadmap (lossless convexification with `T_min > 0`, and the
free-final-time variant) and matches the course notes (lezione 11:
convessificazione, log-mass, lossless slack; lezione 12: why no SCvx is needed
when the problem is convex).

Baseline numbers to beat (same laptop, batch session, `tf = 38 s`, `N = 50`):

| Variant                  | m_f [kg] | wall [s] | replay pos err [m] |
|--------------------------|----------|----------|--------------------|
| Trapezoidal (fmincon)    | 1403.20  | ~90      | 4.3                |
| ZOH RK4 (fmincon)        | 1403.37  | ~100     | 6e-6               |
| SCvx (fmincon inner)     | 1399.54  | ~82 (*)  | 0.18               |
| SCvx (YALMIP+ECOS inner) | 1400.84  | ~29 (*)  | 0.019              |

(*) plus the trapezoidal warm start.

## Mathematical formulation

Change of variables (all in the existing non-dim scheme, `a_ref = g`):

- thrust acceleration `a = T/m`, slack `sigma >= ||a||`, log-mass `z = ln m`
  (non-dim `m0 = 1` so `z(0) = 0`);
- dynamics become **linear**: `dv/dt = a + g_vec`, `dz/dtau = -Vc * sigma`
  with `Vc = V_ref/c` (same residual parameter as Task 1/2);
- cost: maximize `z_N` (monotone in `m_f`);
- thrust magnitude: `||a|| <= sigma` is a second-order cone (lossless: the
  optimizer pushes `sigma` down onto `||a||` because `sigma` burns mass);
- upper bound `sigma <= Tmax_nd * exp(-z)` is non-convex; impose the tangent
  (affine) bound about a fixed profile `z0(t)`:
  `sigma_k <= Tmax_nd * exp(-z0_k) * (1 - (z_k - z0_k))`.
  Convexity of `exp(-z)` makes the tangent an **under**-estimate, so the
  affine constraint is an inner (safe, slightly conservative) approximation —
  feasibility is never lost. Use the max-burn mass profile
  `z0(tau) = ln(1 - Vc * Tmax_nd * tau)` as reference; optionally one
  re-linearization solve about the obtained `z(t)` to shave the conservatism
  (2 solves total, still no SCvx loop).
- `T_min = 0` makes the lower bound trivial (`sigma >= 0`). For the
  `T_min > 0` showcase the quadratic lower bound of Açikmeşe-Ploen
  (`Tmin_nd * exp(-z0)*(1 - (z-z0) + (z-z0)^2/2) <= sigma`) is convex and
  drops straight in.
- glide-slope, altitude, boundary conditions: unchanged, already linear.

Discretization: ZOH on `(a, sigma)`. In the new variables the dynamics are
LTI, so the discrete matrices are **exact closed form** (double integrator +
gravity; `z_{k+1} = z_k - Vc*sigma_k*dt`) — no STM integration, no ode45 in
the transcription at all.

Note the modeling shift: ZOH on *acceleration*, not on thrust. The physical
thrust within an interval is `T(t) = m(t)*a_k`, decreasing with `m(t)`; since
the affine bound is imposed at the left node (where `m` is largest), the
intra-interval thrust bound is automatically respected.

## Implementation plan

1. New script `HM2_powered_descent/main_task2_socp.m` (keeps `main_task2.m`
   runtime manageable; reuse `nondim`/`dim_sol` helpers by copy, same style).
2. Build the exact ZOH LTI matrices in closed form; assemble the SOCP in
   YALMIP (`X` 4xN pos/vel, `Z` 1xN log-mass, `A` 2x(N-1), `S` 1x(N-1));
   solve with ECOS; measure wall time (parse + solve, single call).
3. Post-checks:
   - **losslessness**: report `max_k (sigma_k - ||a_k||)` — must be ~solver
     tolerance; this is the empirical proof of the relaxation theorem;
   - reconstruct `T = exp(z) .* a`, verify `||T|| <= Tmax` at nodes and
     between nodes (replay);
   - ode45 replay with `T(t) = m(t)*a_k` (consistent with the a-ZOH model):
     touchdown dispersion + node errors, same metric as Task 2;
   - optional second solve re-linearized about the obtained `z(t)`.
4. Figures (`task2d_*` prefix, light theme): trajectory/thrust/mass overlay
   vs the four existing variants; **sigma vs ||a|| overlay** (the signature
   lossless plot) replaces the SCvx convergence trace (single solve — there
   is no iteration history).
5. Report: new subsection in `Task2.tex` ("Variant (d): lossless
   convexification — a single SOCP") with the derivation above, the
   losslessness argument (cite acikmese2007, already in ref.bib; lezione 11),
   extended comparison tables (m_f, fidelity, dispersion, wall time), and the
   punchline: certificate + sub-second solve vs minutes of NLP.
6. README: results row + tick the "Lossless convexification" roadmap box
   (point it at this ticket).
7. Stretch goals (separate acceptance, can slip):
   - `T_min = 20 kN` re-solve with the quadratic lower bound — shows the
     coast arc being replaced by a min-throttle arc (the whole point of
     lossless convexification);
   - free final time: outer 1-D search on `tf` (golden section over the SOCP,
     each evaluation < 1 s) — closes the "free-time variant" roadmap item.

## Acceptance criteria

- [ ] Single ECOS solve (no SCvx loop, no NLP warm start) returns
      `m_f` within ~1 kg of the four existing transcriptions
- [ ] Total wall time (YALMIP parse + ECOS) under ~2 s at `N = 50`
- [ ] Losslessness verified: `max(sigma - ||a||)` at solver-tolerance level,
      reported in script output and in the report
- [ ] ode45 replay touchdown dispersion at the ZOH-RK4 level (sub-mm/cm),
      thrust bound respected between nodes
- [ ] `Task2.tex` subsection + comparison tables updated, PDF recompiles
      with no new overfull boxes / undefined refs
- [ ] HM2 README results + roadmap updated
- [ ] (stretch) `T_min > 0` showcase with quadratic lower bound
- [ ] (stretch) free-`tf` golden-section search

## Notes

- Pitfall: do **not** reuse the Appendix-A STM machinery here — the LTI
  closed form is exact and the whole point is that no integration is needed.
- Pitfall: the replay must drive the dynamics with `T = m(t)*a_k` (a-ZOH),
  not `T_k` constant, or the fidelity check will measure the wrong model.
- The affine upper bound is an inner approximation: expect `m_f` a hair
  below the true optimum on the first solve; the optional re-linearization
  pass should close the gap to < 0.1 kg.
- Dependencies: YALMIP + ECOS (already used by Task 2 variant (c)).
- References: Açikmeşe & Ploen 2007, Blackmore et al. 2010 (both in
  `ref.bib`), course notes lezione 11 (log-mass + slack) and 12.
- Related: HM2 README roadmap items "Lossless convexification" and
  "Free-time variant"; this ticket supersedes both.
