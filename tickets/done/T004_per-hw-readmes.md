---
id: T004
title: Add a local README.md to each homework folder
priority: medium
created: 2026-04-26
updated: 2026-04-26
---

## Context
Each `HM<N>/` should be self-explanatory when opened directly on GitHub. A
local README per folder shortens the path from "user clicks the folder" to
"user understands what's there".

Standard structure for each HM README:
1. Problem statement (1 paragraph + link to the assignment PDF in
   `DCLV_MATERIALE_CORSO_26042026/`)
2. Approach (which numerical method, why)
3. How to run (`matlab -batch "main_task1; main_task2; ..."`)
4. Results — embed 2–4 key figures with one-line captions
5. Files map: which `.m` does what

## Acceptance criteria
- [x] `HM0_falcon9_ascent/README.md`
- [x] `HM1/README.md`
- [x] `HM2_powered_descent/README.md` (delivered with T003)
- [x] All three render cleanly on github.com (markdown is plain GFM, tables +
      relative-path image links — verified locally; no GH-only tricks used)

## Notes
- Don't duplicate content from the assignment PDF — link to it.
- Keep the figures inline (relative path) — don't push to an external CDN.

## Resolution
- **Course-PDF links dropped**: the assignment PDFs are .gitignored
  (T005 decision), so a public link to them would 404 on github.com. Each
  README inlines the problem statement instead — short, no duplication risk.
- **Figure generation**: no figures existed in the homework folders. Solved
  by appending a non-intrusive "EXPORT FIGURES" block at the end of each
  `main*.m` (after the last `figure(...)` call, before any local function).
  The block iterates open figures, slugifies their `Name` property, and
  writes `figures/<prefix>_<slug>.png` at 200 dpi via `exportgraphics`.
  This keeps each script runnable in isolation while populating the README's
  embedded images.
- **HM0 `documentazione.txt`**: folded into the HM0 README as a short
  English explanation of the non-dimensional `main2.m` variant. Original
  Italian file kept in place as private context (mentioned in the README's
  files map).
- **HM1 hero figure**: `task1a_final_mass_vs_q.png` — three `y_f` curves
  with clear interior maxima — is the strongest single image to put at the
  top of the eventual public README (T002).
- **HM2 hero candidate**: `main_task1` plots are produced live; if T002
  needs a static figure, run HM2 once and screenshot the trajectory + thrust
  plots, or add the same export block.
