---
id: T005
title: Repo hygiene — .gitignore, LICENSE, and git init check
priority: low
created: 2026-04-26
updated: 2026-04-26
---

## Context
Before publishing on GitHub the repo needs basic hygiene. Current state:
- `git init` already done by the user (branch `master`, zero commits yet).
- No `.gitignore` — MATLAB produces `*.asv` autosaves, `slprj/` Simulink build
  artifacts, `*.slxc` cache files, etc. that should not be tracked.
- No `LICENSE`.
- Default branch is `master`; GitHub default is `main` — consider renaming
  before the first push.

## Acceptance criteria
- [x] `git init` run at repo root
- [x] Default branch decision: renamed to `main` (done by user)
- [x] `.gitignore` covering MATLAB / Simulink / OS junk
- [x] `LICENSE` file — MIT, © Niccolò D'Ambrosio
- [x] Decided to `.gitignore` `DCLV_MATERIALE_CORSO_26042026/` plus the
      `**/LVdynamics_*.pdf`, `**/Homework*.pdf`, `**/Classwork*.pdf` patterns
      to catch copies leaked into HM folders
- [x] Initial commit clean — `6d6e20e chore: initial repo scaffold`

## Notes
- Check size of any committed `.mat` / `.png` / `.svg` before the first push;
  GitHub warns above 50 MB and rejects above 100 MB.
- The `DCLV_MATERIALE_CORSO_26042026/` PDFs are course material — confirm with
  the user whether redistribution is OK before publishing them publicly.

## Resolution
- User had already done `git init` and `git branch -M main` before the chat.
- `.gitignore` covers MATLAB autosaves/MEX, Simulink build cache, IDE configs,
  OS junk, plus the entire course material folder + any duplicated PDFs.
- LICENSE: MIT, copyright "Niccolò D'Ambrosio" 2026.
- Removed scratch file `HM0_falcon9_ascent/Untitled-2.txt` per user direction.
- Kept `main2.m`, `main2_backup.m`, `documentazione.txt` in HM0 for now —
  user will revisit the secondary `main2.m` variant later; documentation
  format will be normalized when T004 (per-HM READMEs) lands.
