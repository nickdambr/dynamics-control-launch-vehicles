# Tickets

Lightweight ticket system for tracking work on the DCLV repo. Each ticket is a
markdown file with YAML frontmatter; **status is encoded by which folder the
file lives in** — moving the file *is* the state transition.

## Folders

```
tickets/
├── open/          # not started
├── in-progress/   # currently being worked on
└── done/          # completed (kept for history)
```

## Ticket format

Filename: `T<NNN>_short-slug.md` (zero-padded id, e.g. `T001_rename-hm0.md`).

```markdown
---
id: T001
title: One-line summary
priority: high | medium | low
created: YYYY-MM-DD
updated: YYYY-MM-DD
---

## Context
Why this ticket exists. Background, motivation, links to relevant files.

## Acceptance criteria
- [ ] concrete, checkable item
- [ ] another concrete, checkable item

## Notes
Free-form: dependencies, decisions, open questions.
```

## Priority

- **high** — blocking the GitHub publication, or due to course deadlines.
- **medium** — improves the repo but not blocking.
- **low** — nice-to-have polish.

## Workflow

1. Open: create the file in `open/`, fill out context + acceptance criteria.
2. Start: move the file to `in-progress/`, bump `updated:`.
3. Finish: tick all acceptance criteria, move the file to `done/`, bump `updated:`.

Don't delete done tickets — they're the project history.

## Index

Quick view of what's where: `ls tickets/open tickets/in-progress tickets/done`.
