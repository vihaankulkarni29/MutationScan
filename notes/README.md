# Notes

Project reports and developer notes. Everything here is committed and intended to
be read by people reviewing the work.

| Document | What it is |
|---|---|
| [`benchmarks.md`](benchmarks.md) | Measured performance, with the machine described and commands to reproduce each figure |
| [`developer-notes.md`](developer-notes.md) | Dated record of what changed and why, plus open items |

Historical v2 → v3 rebuild documents are in
[`docs/archive/`](../docs/archive/README.md). That code was removed in `d5f0bb0`,
so those documents describe a system that no longer exists.

## What belongs here

- **Benchmark and performance output**, including the machine it was measured on.
  A number without its environment is not a result.
- **Developer notes** for substantive changes: what moved, why, what it affected,
  and how it was verified.
- **Anything a reviewer would otherwise have to ask for twice.**

## What does not belong here

- Anything the user-facing [`README`](../README.md) should say instead. If a
  feature needs explaining to someone running the tool, it goes there.
- Scratch work, half-finished analysis, or numbers not yet checked.
- Agent-internal working memory. That is kept strictly local and outside the
  repository.

## Adding a report

State the machine, the cohort, and the command. Report variance rather than a
single best figure — if a measurement moved between runs, say so and explain
which number you trust. `notes/benchmarks.md` is the worked example: it labels
best-of-N results as such and states the 2x swing observed on identical work,
because a reader who finds those numbers differ needs to know why before they
conclude something is broken.