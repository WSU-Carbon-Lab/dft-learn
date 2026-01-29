# StoBe Transition-Potential DFT Workflow: Issue Breakdown and Remediations

This document provides a structured analysis of issues in `generate.py` and `schedule.py` relative to the StoBe transition-potential (TP) NEXAFS workflow, and concrete remediations to restore parity with `StoBe_Input_Generator.py` / `molConfig.py` and support the desired output layout (GND/EXC/TP/NEXAFS, optional packaging).

---

## Summary

| Area | Critical issues | Remediation focus |
|------|-----------------|-------------------|
| **generate.py** | No auxiliary/orbital/MCP/augmentation basis blocks; XAS `.xas` path wrong; SPIN omitted | Implement full basis emission (geometry order); fix XAS symlink; add SPIN |
| **schedule.py** | Outputs stay in C1/...; no GND/EXC/TP/NEXAFS layout; run scripts executed with bash (csh shebang); `all` mis-calls `calculate_optimal_workers` | Run-file redirection or post-run organize; use csh or convert to bash; pass `calc_types` to workers |
| **Workflow** | alfaOcc/alfaOccTP correct; LUMO index often updated manually after GND | Keep; optional helper to update from Molden |

---

## 1. Workflow Context

### 1.1 Transition-potential NEXAFS sequence

1. **Ground state (gnd)**
   SCF to obtain reference orbitals and energies.

2. **Excited state (exc)**
   Core hole in 1s, electron in LUMO; no relaxation.
   Uses `ALFA 0 1 <n_virt> 0.0` (e.g. `0 1 14 0.0`). The integer is the **LUMO index**; often updated after inspecting the ground-state Molden (e.g. 1s -> LUMO).

3. **Transition potential (tp)**
   Half core hole: 0.5 in core, 0.5 in same virtual.
   Uses `ALFA 0 1 <n_virt> 0.5` (e.g. `0 1 14 0.5`). Same `<n_virt>` as exc.

4. **XAS (xas)**
   `xrayspec.x` post-processes the `.xas` from TP to produce the NEXAFS spectrum.

### 1.2 Reference implementations

- **Input generation**: `StoBe_Input_Generator.py` + `molConfig.py` (geometry, effective charges, basis sets, alfaOcc/alfaOccTP).
- **Example run files**: `.cursor/skills/stobe-pro/references/examples/` (C1gnd.run, C1exc.run, C1tp.run, C1xas.run).
- **Downstream consumers**: `streamlit_demo/stobe_loader.py` expects outputs under `GND/`, `EXC/`, `TP/`, `NEXAFS/`.

---

## 2. Issues in `generate.py`

### 2.1 Orbital / basis sections missing (critical)

**Observed behavior**

After the electronic-state block (FSYM, ALFA, BETA, FILE, etc.), `makeRunFile` writes only:

```python
# Add basis set sections (simplified - could be expanded)
# Auxiliary basis sets
f.write("END\n")
```

Then it closes the heredoc and appends the StoBe invocation. No auxiliary, orbital, MCP, or augmentation blocks are emitted.

**Required structure (see stobe-pro examples and xray-spectroscopy.md)**

Between the electronic-state `END` and the final `END` that closes the input:

1. **Auxiliary basis**
   One line per atom in geometry order (`A-CARBON (5,2;5,2)`, `A-CARBON(+4) (3,3;3,3)`, etc.).

2. **Orbital basis**
   One line per atom (`O-CARBON iii_iglo`, `O-CARBON(+4) (311/211/1)`, etc.).

3. **MCP basis** (if used)
   One line per target atom, e.g. `P-CARBON(+4) (3,1:8,0)`.

4. **Augmentation**
   `X-FIRST` for the core-hole atom, `X-DUMMY` for all other target atoms.

5. **Final `END`**
   To terminate the basis block.

**Impact**

Generated `.inp` files are invalid: StoBe expects basis sets in geometry order. Calculations fail or are ill-defined.

**Remediation**

- Reuse the **element-order logic** from `StoBe_Input_Generator.makeRunFile` (including `difElements`, `nElem1..6`, `element1a`/`element1b`, and basis strings from molConfig).
- Emit, in order:
  - Auxiliary basis lines (one per atom, `elem1_Abasis_a` vs `elem1_Abasis_b` for the excited atom).
  - Orbital basis lines (same ordering, `elem1_Obasis_a` vs `elem1_Obasis_b`).
  - MCP lines if `elem1_MCPbasis` (and similar) are defined.
  - Augmentation: `X-FIRST` for the excited atom, `X-DUMMY` for the rest.
- Ensure one-to-one correspondence with geometry lines (same order as in the heredoc). Prefer iterating over the parsed geometry explicitly rather than inferring from generic “run file” lines.

---

### 2.2 SPIN keyword omitted

**Observed behavior**

`makeRunFile` writes ORBI, MULLIKEN, VIRT but never `SPIN`. `StoBe_Input_Generator` writes `SPIN FULL` before the FSYM block.

**Remediation**

- Add `SPIN {spin}` (or default `SPIN FULL`) in the same place as in the reference (before FSYM), using `molConfig.spin` when available.

---

### 2.3 Geometry parsing and formatting differences

**Observed behavior**

- `generate.py` uses `line.split()` and `len(parts) >= 4`, and always appends `\n` to geometry lines.
- `StoBe_Input_Generator` has special cases (e.g. last atom of a given element) where the final newline is omitted for certain `difElements` branches.

**Impact**

Mostly cosmetic, but last-line formatting can matter for strict parsers. Inconsistent with the original generator.

**Remediation**

- Align geometry writing with `StoBe_Input_Generator`: same handling of “last atom” per element group and newline rules.
- Use a single, shared geometry iterator so that basis-set emission can use the same order unambiguously.

---

### 2.4 XAS run file path for `.xas` input

**Observed behavior**

`makeXASrun` uses:

```python
f.write(f"ln -s ~/{mname}/{aname}{i}/{aname}{i}.xas fort.1\n")
```

So it expects `~/<mname>/<aname><i>/<aname><i>.xas` (e.g. `~/tcta_flat/C1/C1.xas`).

**Actual layout**

- Run files live in `C1/`, `C2/`, ... (same as `aname`+index).
- TP writes `fort.11` -> `C1.xas` in the **same** directory as the run (e.g. `C1/`).
- XAS runs **in** that directory too.

So `C1.xas` is `./C1.xas` or `${PWD}/C1.xas` from `C1/`, not `~/mname/C1/C1.xas`.

**Impact**

XAS run fails: `fort.1` symlink target does not exist.

**Remediation**

- Use a path relative to the run directory, e.g.
  `ln -s ${PWD}/${aname}${i}.xas fort.1`
  or `ln -s ./${aname}${i}.xas fort.1`, consistent with the stobe-pro example and `StoBe_Input_Generator` (`ln -s ${PWD}/...`).

---

### 2.5 StoBe / xrayspec paths and run script interpreter

**Observed behavior**

- `generate.py` uses `StoBe.x` and `xrayspec.x` from `PATH` or `STOBE_HOME`, and writes direct invocations (e.g. `StoBe.x`, `{XASRUN}`).
- Run files use `#!/bin/csh -f` (csh).

**Remediation**

- Keep using `PATH` / `STOBE_HOME` if that is the desired deployment model; alternatively allow molConfig overrides.
- Ensure `schedule.py` invokes run scripts with the **same shell** as the shebang (`csh`), or that generated scripts are valid for the shell actually used (see Section 3.3).

---

## 3. Issues in `schedule.py`

### 3.1 Output location: no GND/EXC/TP/NEXAFS layout

**Observed behavior**

- `schedule.py` only **runs** `.run` files; it does not create directories or move outputs.
- Generated run files write outputs next to the script (e.g. `C1/C1gnd.out`, `C1/C1gnd.molden`, `C1/C1.xas`).

**Desired layout (e.g. polystyrene “working” case, streamlit_demo)**

- Outputs grouped by calculation type: `GND/`, `EXC/`, `TP/`, `NEXAFS/`.
- Optional `packaged_output.tar.gz` (or similar) for archiving.

**Impact**

- All outputs stay under `C1/`, `C2/`, ...
- `stobe_loader` (and similar tools) expect `GND/`, `EXC/`, `TP/`, `NEXAFS/` and will not find files.

**Remediation (choose one or combine)**

**Option A – Generate run files that write into type-specific dirs**

- Emit run scripts that create `../GND/`, `../EXC/`, `../TP/`, `../NEXAFS/` (or configurable roots) and redirect outputs there, e.g.
  `StoBe.x <C1gnd.inp> &> ../GND/C1gnd.out`,
  `mv Molden.molf ../.../C1gnd.molden`, etc.
- Same idea as in stobe-pro examples (`../gnd/`, `../exc/`, `../tp/`, `../molden/`). NEXAFS can go to `../NEXAFS/` or equivalent.
- `schedule.py` continues to only run scripts; no extra moving of files.

**Option B – Post-run organization in `schedule.py`**

- After each run, move `*gnd.out`, `*gnd.molden`, `*exc.out`, etc. into `GND/`, `EXC/`, `TP/`, and `*xas.out` / `*.out` (NEXAFS) into `NEXAFS/`.
- Require a single “output root” (e.g. project dir) so that `GND/`, `EXC/`, etc. are shared across `C1/`, `C2/`, ...

**Option C – Dedicated “organize” / “package” command**

- New subcommand, e.g. `schedule.py organize <directory>`, that:
  - Scans for `*gnd.out`, `*exc.out`, `*tp.out`, `*xas.out`, NEXAFS `*.out`, etc.
  - Creates `GND/`, `EXC/`, `TP/`, `NEXAFS/` under a configurable root.
  - Moves (or copies) files accordingly.
- Optional: build `packaged_output.tar.gz` from those directories.

---

### 3.2 “Copies in same directory” / execution directory

**Observed behavior**

- `schedule.py` finds `**/*{calc_type}.run` under the given directory and runs each script with `cwd=run_file.parent` (e.g. `C1/`).
- It does **not** copy run files elsewhere; it executes them in place.

**Interpretation**

- “Copies of each run file in the same directory it is executed in” may mean:
  - Outputs (`.out`, `.molden`, `.xas`) are written next to the `.run` in `C1/`, etc., instead of a dedicated output tree (GND/EXC/TP/NEXAFS), or
  - Some other duplicate setup when using both generator and scheduler.

**Remediation**

- Clarify intended layout: either (a) outputs under `C1/`, ... and then organization via Option B or C above, or (b) run files already redirect into `GND/`, `EXC/`, etc. (Option A).
- Avoid overloading “run directory” with “output directory”; make output root explicit (flag or config).

---

### 3.3 Run script interpreter: bash vs csh

**Observed behavior**

```python
subprocess.run(["bash", run_script], cwd=run_dir, ...)
```

Run files have `#!/bin/csh -f` and may use csh-only syntax.

**Impact**

Csh-specific constructs can fail under `bash`.

**Remediation**

- Run scripts with `csh -f` (or `tcsh -f` if that’s the deployment default), e.g.
  `["csh", "-f", run_script]`,
  and ensure `csh` is available on target systems, **or**
- Change generated run files to `#!/bin/bash` and plain POSIX/bash syntax so that `bash` remains valid.

---

### 3.4 Job key overwrite when multiple calc types per directory

**Observed behavior**

In `run_calculations_with_progress`:

```python
for dir_name, run_files in directories.items():
    for run_file in run_files:
        job = DirectoryJob(...)
        self.jobs[dir_name] = job
```

Jobs are keyed only by `dir_name`. Multiple run files per directory (e.g. gnd, exc, tp, xas) overwrite the same key; only one job per directory is retained.

**Current usage**

- `schedule.py all` calls `run_calculations_with_progress` once per calc type (`gnd`, `exc`, `tp`, `xas`), so each call sees at most one run file per directory. The overwrite does not currently change behavior.

**Remediation**

- If you add a mode that runs multiple calc types in a single call, key jobs by `(dir_name, calc_type)` or `(dir_name, run_file)` so all run files are executed.
- Optionally, unify “run all types for each directory” in one pass and use the extended keying.

---

### 3.5 `calculate_optimal_workers` and `--atom`

**Observed behavior**

- With `--atom C1`, patterns like `**/C1gnd.run` are used; workers are chosen from directories that contain those runs.
- Document that `--atom` expects full atom IDs (e.g. `C1`, `C2`), not only element names.

**Remediation**

- Keep current behavior; add a short note in CLI help and in this doc.

---

### 3.6 `all` command: wrong `calculate_optimal_workers` call

**Observed behavior**

```python
workers = calculate_optimal_workers(str(directory), atom or [])
```

`calculate_optimal_workers` signature is `(scan_dir, calc_types, atoms=None)`. The `all` command passes `(directory, atom or [])`, so `calc_types` is effectively `atom or []` (often `[]`) and `atoms` is never passed. With `calc_types == []`, no run files are found, so `num_directories == 0` and `optimal_workers == 0`.

**Impact**

Worker count can be zero when using `schedule.py all` with auto workers, breaking parallel execution.

**Remediation**

- Call `calculate_optimal_workers(str(directory), calc_types, atom)` in the `all` command, with `calc_types = ["gnd", "exc", "tp", "xas"]`.

---

## 4. AlphaOcc / orbital indexing (alfaOcc, alfaOccTP)

### 4.1 Role of `0 1 <n> 0.0` and `0 1 <n> 0.5`

- **Exc**: `ALFA 0 1 <n> 0.0` – core hole, electron in orbital `<n>` (LUMO).
- **TP**: `ALFA 0 1 <n> 0.5` – half in core, half in orbital `<n>`. Same `<n>` as exc.

`<n>` is often determined from the ground-state Molden (1s -> LUMO) and may differ from the default `1`.

### 4.2 Current generator behavior

- `generate.py` forwards `alfaOcc` and `alfaOccTP` from `molConfig` into `makeRunFile` and writes them into the FSYM/ALFA blocks. No bug here.
- The **omission of basis blocks** is what prevents correct calculations; fixing that is prioritary.

### 4.3 Post–ground-state updates

- If users update alfaOcc/alfaOccTP (e.g. `0 1 8 0.0` / `0 1 8 0.5`) after inspecting GND, that is done **manually** in molConfig or by editing generated inputs.
- Optional future improvement: a small helper or `schedule` subcommand that updates alfaOcc/alfaOccTP in existing `.inp` or `.run` from a simple mapping (e.g. LUMO index per atom).

---

## 5. Package / streamlit_demo alignment

- **Folders**: `stobe_loader` expects `GND/`, `EXC/`, `TP/`, `NEXAFS/`. Resolving Section 3.1 (output layout) is required for compatibility.
- **Naming**: Match expected suffixes (`*gnd.out`, `*exc.out`, `*tp.out`, `*nexafs.out` or equivalent) and atom extraction logic.
- **`packaged_output.tar.gz`**: Not implemented in the repo. Can be added as part of an “organize”/“package” step (Section 3.1, Option C).

---

## 6. Remediation priority

| Priority | Item | Script | Effort |
|----------|------|--------|--------|
| P0 | Emit full basis blocks (auxiliary, orbital, MCP, augmentation) in geometry order | `generate.py` | High |
| P0 | Fix XAS `.xas` path (use `${PWD}` or `./`) | `generate.py` | Low |
| P1 | Add `SPIN` keyword | `generate.py` | Low |
| P1 | Align geometry formatting with `StoBe_Input_Generator` | `generate.py` | Medium |
| P1 | Output layout: run files write to GND/EXC/TP/NEXAFS or add organize/package step | `generate.py` and/or `schedule.py` | Medium |
| P2 | Use csh for `.run` execution or switch scripts to bash | `schedule.py` / `generate.py` | Low |
| P2 | Key jobs by `(dir, calc_type)` or `(dir, run_file)` if multiple types per dir | `schedule.py` | Low |
| P3 | Optional “organize”/“package” CLI and `packaged_output.tar.gz` | `schedule.py` | Medium |

---

## 7. Suggested implementation order

1. **Basis and XAS path in `generate.py`**
   Implement full basis-set emission (Section 2.1) and fix `makeXASrun` (Section 2.4). Add `SPIN` (Section 2.2).
   This restores valid StoBe inputs and a working XAS step.

2. **Output layout**
   Decide between run-file redirection (Option A) vs post-run organization (B) or explicit organize/package (C). Implement the chosen strategy so that `GND/`, `EXC/`, `TP/`, `NEXAFS/` exist and are populated.

3. **`schedule.py`**
   Fix `all` command `calculate_optimal_workers` call (pass `calc_types`). Switch run script execution to csh if keeping csh shebangs; otherwise convert scripts to bash. Fix job keying if you add multi-type-per-dir execution.

4. **Geometry and packaging**
   Align geometry handling with `StoBe_Input_Generator`; add organize/package and `packaged_output.tar.gz` if desired.

---

## 8. References

- `scripts/StoBe_Input_Generator.py` – geometry, basis, and run-file structure.
- `scripts/molConfig.py` – molConfig variables and defaults.
- `.cursor/skills/stobe-pro/references/examples/` – C1gnd.run, C1exc.run, C1tp.run, C1xas.run.
- `.cursor/skills/stobe-pro/references/xray-spectroscopy.md` – TP setup, alfaOcc, basis, XAS.
- `.cursor/skills/stobe-pro/references/calculation-workflows.md` – workflow order and layout.
- `streamlit_demo/stobe_loader.py` – expected folder layout and file naming.
