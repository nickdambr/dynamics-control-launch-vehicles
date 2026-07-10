---
id: T001
title: Rename DV_HM0/ to HM0_falcon9_ascent/
priority: medium
created: 2026-04-26
updated: 2026-04-26
---

## Context
The HM0 folder is named `DV_HM0/`, while HM1 follows the cleaner `HM<N>/`
scheme. The agreed convention is `HM<N>_<short_topic>/`. Renaming aligns the
repo with the convention and makes the homework index in the future public
README easier to scan.

## Acceptance criteria
- [x] Folder renamed to `HM0_falcon9_ascent/`
- [x] All relative paths inside the folder's `.m` files still work (figures
      are saved alongside the script — should be fine, but verify)
- [x] Homework index updated (remove the legacy-name caveat)
- [x] Any reference to `DV_HM0` elsewhere in the repo updated (grep first)
- [x] HM0 still runs end-to-end after the rename

## Notes
- Do this BEFORE T002 (public README) so the README links to the final name.
- If git history matters for the move, use `git mv` rather than delete + create.
- **Resolution:** repo not yet under git (T005), so used a plain `mv`. Smoke-test:
  ran `main.m` headless via `matlab -batch` from `HM0_falcon9_ascent/` — all
  computations and event detection (Mach 1, max-Q) reproduce the expected
  values; figures saved correctly into the renamed folder.
- Grep for `DV_HM0` confirmed clean post-edit (only this archived ticket
  retains the legacy name in its title, by design).
