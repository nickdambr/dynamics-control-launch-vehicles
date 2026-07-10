---
id: T002
title: Write public-facing README.md for GitHub
priority: high
created: 2026-04-26
updated: 2026-04-26
---

## Context
The repo will be published on GitHub as a portfolio piece for GNC engineering
roles. We need a `README.md` at the root targeted at:
1. **Recruiters / hiring managers** scrolling fast — clear value prop in 30 s.
2. **GNC engineers / professors** evaluating depth — one click to runnable code
   and to a representative figure per homework.

## Acceptance criteria
- [x] `README.md` at repo root
- [x] Top section: one-paragraph pitch (course, what it demonstrates, link to
      author profile)
- [x] Homework table with: topic, methods (PMP, SCvx, …), one hero figure per
      HM, link to the folder
- [x] "How to run" section: MATLAB version + toolboxes required
- [x] License + author attribution
- [x] Looks good rendered on github.com (markdown sanity-checked locally)

## Notes
- Pick the strongest plot from each completed HM as the hero figure.
- Keep tone professional but not stiff — this is a portfolio, not a thesis.
- Depends on T001 (folder rename) so README links use the final names.

## Resolution
- **Hero figures**:
  - HM0 → `3d_trajectory.png` (visual story of the entire ascent + events).
  - HM1 → `task1a_final_mass_vs_q.png` (clear interior maximum, single image
    that motivates the indirect-optimization approach).
  - HM2 → `task1_trajectory.png` (three sensitivity runs + glide-slope
    corridor in one frame).
- **HM2 figures**: hadn't been exported yet (HM2 README from T003 only
  described what the script would plot). Added the same `EXPORT FIGURES`
  block to `HM2_powered_descent/main_task1.m`, regenerated, then updated
  the HM2 README with a 2×2 figure grid.
- **Badges**: status badges for license, MATLAB version, and per-HM
  completion. Plain shields.io URLs — no signup, no rate-limit risk.
- **Toolbox claims**: confirmed `fmincon` and `fsolve` are both in the
  Optimization Toolbox; nothing else from Robust Control / Global Opt is
  used yet. Updated when HM3 lands.
- **License footer**: explicit note that course PDFs are © Prof. Zavoli
  and intentionally excluded — pre-empts any redistribution complaints.
