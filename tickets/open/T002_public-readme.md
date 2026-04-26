---
id: T002
title: Write public-facing README.md for GitHub
priority: high
created: 2026-04-26
updated: 2026-04-26
---

## Context
The repo will be published on GitHub as a portfolio piece for GNC engineering
roles. `CLAUDE.md` is for AI/maintainer context, not for human visitors. We
need a `README.md` at the root targeted at:
1. **Recruiters / hiring managers** scrolling fast — clear value prop in 30 s.
2. **GNC engineers / professors** evaluating depth — one click to runnable code
   and to a representative figure per homework.

## Acceptance criteria
- [ ] `README.md` at repo root
- [ ] Top section: one-paragraph pitch (course, what it demonstrates, link to
      author profile)
- [ ] Homework table with: topic, methods (PMP, SCvx, …), one hero figure per
      HM, link to the folder
- [ ] "How to run" section: MATLAB version + toolboxes required
- [ ] License + author attribution
- [ ] Looks good rendered on github.com (preview before declaring done)

## Notes
- Pick the strongest plot from each completed HM as the hero figure.
- Keep tone professional but not stiff — this is a portfolio, not a thesis.
- Depends on T001 (folder rename) so README links use the final names.
