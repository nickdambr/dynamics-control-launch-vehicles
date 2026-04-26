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
- [ ] Default branch decision: keep `master` or rename to `main`
- [ ] `.gitignore` covering MATLAB / Simulink / OS junk (Thumbs.db, .DS_Store,
      `*.asv`, `slprj/`, `*.slxc`, `*.mexw64`, `*.mat` if heavy, etc.)
- [ ] `LICENSE` file (decision pending — ask user MIT vs Apache-2.0)
- [ ] Decide whether to commit `DCLV_MATERIALE_CORSO_26042026/` (course PDFs:
      redistribution rights unclear) or `.gitignore` it
- [ ] Initial commit clean (no IDE config files, no large binaries leaked)

## Notes
- Check size of any committed `.mat` / `.png` / `.svg` before the first push;
  GitHub warns above 50 MB and rejects above 100 MB.
- The `DCLV_MATERIALE_CORSO_26042026/` PDFs are course material — confirm with
  the user whether redistribution is OK before publishing them publicly.
