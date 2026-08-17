# dft-learn (Claude Code)

Follow **[AGENTS.md](AGENTS.md)** as the project spec: package layout (`src/dftlearn/`), `dftrun` CLI, Python 3.12+ / **uv** / **Ruff** / **ty** / **pytest**, NumPy-style public docstrings, and Igor/StoBe context.

## Skills

Project skills are canonical under `.agents/skills/`. This repo links that tree to **`.claude/skills/`**, so Claude Code discovers the same `SKILL.md` files it would under `.claude/skills/<name>/SKILL.md`. Cursor uses `.cursor/skills/` for the same tree.

Load the matching skill **before** implementing. Each skill’s `SKILL.md` plus `references/` is the long-form contract.

| Skill | Use when |
|-------|----------|
| [`general-python`](.claude/skills/general-python/SKILL.md) | uv, ruff, ty, pytest, dataclasses, typing, scientific defaults |
| [`numpy-scientific`](.claude/skills/numpy-scientific/SKILL.md) | NumPy dtypes, views, broadcasting, ufuncs, linalg |
| [`dataframes`](.claude/skills/dataframes/SKILL.md) | pandas / Polars tables, joins, lazy plans, I/O |
| [`numpy-docstrings`](.claude/skills/numpy-docstrings/SKILL.md) | numpydoc on public APIs |
| [`matplotlib-scientific`](.claude/skills/matplotlib-scientific/SKILL.md) | publication Matplotlib figures |
| [`lab-instrumentation`](.claude/skills/lab-instrumentation/SKILL.md) | PyVISA, sockets, HAL, instrument tests |
| [`general`](.claude/skills/general/SKILL.md) | completeness, public contracts, no placeholder code |

## Subagents and rules

Subagents live in `.agents/agents/` and are linked as `.claude/agents/` (`python-reviewer`, `python-types`, `python-refactor`, `dotagent-general-standards-auditor`).

Always-on editor rules live in `.agents/rules/` (linked as `.claude/rules/` and `.cursor/rules/`). The Python rule applies to `**/*.py`: uv, Ruff, ty, NumPy docstrings.

## Tooling

- Environments and tests: `uv sync`, `uv run pytest`, `uv run ruff check`, `uvx ty check`
- Dependencies: `uv add` / `uv remove` only; do not hand-edit version pins in `pyproject.toml`
