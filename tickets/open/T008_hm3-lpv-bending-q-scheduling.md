---
id: T008
title: HM3 LPV showcase — bending Varying Notch + q-scheduled gains (T007 stretch)
priority: low
created: 2026-06-13
updated: 2026-06-13
---

## Context

The full-ascent LPV showcase from [`T007`](../done/T007_hm3-lpv-full-ascent.md)
landed in [`HM3/LTV_FULL_ASCENT/`](../../HM3/LTV_FULL_ASCENT/) with the **rigid**
plant and **time-scheduled** gains. Two stretch goals were deliberately scoped
out of that round and are tracked here. They are independent and can be done
separately. Like T007, this is a portfolio showcase, **not** part of the HM3
deliverable.

Everything needed already exists: `GreensiteLPV_DATA.mat` carries `omega(t)`
(bending frequency), `sigma_ins(t)`, `phi_ins(t)`, `phi_tvc(t)`, `Tc(t)`,
`Q(t)`; [`build_plant_full`](../../HM3/build_plant_full.m) builds the 6-state
flexible plant; [`build_notch_filter`](../../HM3/build_notch_filter.m) and
[`build_tvc`](../../HM3/build_tvc.m) build the notch and TVC. The LPV scaffolding
([`init_simulink_lpv`](../../HM3/LTV_FULL_ASCENT/init_simulink_lpv.m),
[`build_hm3_full_ascent`](../../HM3/LTV_FULL_ASCENT/build_hm3_full_ascent.m),
[`ode_lpv_ascent`](../../HM3/LTV_FULL_ASCENT/ode_lpv_ascent.m)) is the starting
point — extend it, don't fork it.

## Goal 1 — bending mode with time-varying ω(t) + Varying Notch

Lift the rigid LPV plant to the **flexible** 6-state model over the ascent:
add `eta, etadot` with `omega(t)` from the dataset, the INS bending leakage
(`sigma_ins(t)`, `phi_ins(t)`) and the TVC bending forcing (`-phi_tvc*Tc`,
already in the dataset as `aqk`). The HM3 fixed min-phase notch is centered on
`omega(72) = 18.9 rad/s`; as `omega(t)` sweeps the lookup range
(~16.5–31.8 rad/s) it **detunes**, so the natural tool is a **Varying Notch
Filter** retuned on `omega(t)` (cstextras, same family the professor uses for
the Dryden filter). Compare:
- fixed HM3 notch held over the ascent (shows the detuning),
- `omega(t)`-tracking Varying Notch (recovers gain stabilisation).

The LTV ode45 baseline gains two states; the Simulink build adds the bending
block + the Varying Notch in the actuator path (before the TVC, as in HM3).

## Goal 2 — q(t)-scheduled gains (true LPV scheduling)

Replace the **time**-scheduled gains with gains scheduled on the **measurable**
parameter `q(t)` (dynamic pressure) — true LPV scheduling on a physical
variable rather than the clock. Reuse the existing `design_controller` grid but
key the lookup on `Q` instead of `t`: build `Kp_th(q)`, `Kd_th(q)` and, in both
the ode baseline and Simulink, drive the schedule from the live `Q(t)` signal
(a 1-D lookup `Q -> gains`). Discuss the t↔q mapping ambiguity (Q is
non-monotonic: it rises then falls, so a single q maps to two flight times with
different plants — quantify the resulting schedule error and whether ascending
vs descending q branches need separating).

## Implementation plan

1. `init_simulink_lpv`: add the bending coefficient tables (`omega`,
   `sigma_ins`, `phi_ins`, `aqk`) and a `q`-keyed gain schedule alongside the
   existing `t`-keyed one (option flag to pick).
2. `ode_lpv_ascent`: extend the state to 6 (bending), add the Varying-notch
   actuator dynamics, and allow `q`-scheduled gains.
3. `build_hm3_full_ascent`: add the bending integrators + INS leakage outputs,
   the Varying Notch block, and switch the gain lookup key to `Q` for Goal 2.
4. Figures: fixed-vs-varying notch `|L(omega(t))|` over the ascent; t- vs
   q-scheduled response/margins; refresh the overlay validation (~1e-6).
5. Extend `HM3/LTV_FULL_ASCENT/README.md` with the two results.

## Acceptance criteria

- [ ] Flexible LPV plant (6-state, `omega(t)`) simulates over 0–140 s in both
      the ode45 baseline and `hm3_full_ascent.slx`
- [ ] Varying Notch tracking `omega(t)` vs fixed HM3 notch compared (the fixed
      notch detuning shown explicitly)
- [ ] `q(t)`-scheduled gains implemented and compared against the `t`-scheduled
      schedule; non-monotonic-q ambiguity discussed
- [ ] ode45 vs Simulink overlay still ~1e-6 on theta with the new blocks
- [ ] README section + figures; `HM3/` frozen-time deliverable untouched

## Notes

- Depends on **cstextras** (Varying Notch / Varying Transfer Function blocks) —
  confirm availability before relying on it (Control System Toolbox is present;
  these blocks ship with Simulink Control Design, which is installed).
- The Varying Notch needs `omega(t)` as a signal input — feed it from the same
  Clock-driven lookup used for the plant coefficients.
- Keep `strong_wind.slx` read-only (copy the generator subsystem, never save).
- If the bending-augmented loop with the live generator proves stiff, bound the
  solver `MaxStep` as T007 did (the 0.1 s wind noise already forced this).
