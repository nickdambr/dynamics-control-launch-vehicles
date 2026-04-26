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
- [ ] `HM0_falcon9_ascent/README.md`
- [ ] `HM1/README.md`
- [x] `HM2_powered_descent/README.md` (delivered with T003)
- [ ] All three render cleanly on github.com

## Notes
- Don't duplicate content from the assignment PDF — link to it.
- Keep the figures inline (relative path) — don't push to an external CDN.
