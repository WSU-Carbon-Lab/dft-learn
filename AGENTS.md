# Contributor Quick Start Guide

This guide is intended to quickly get you up to speed on the project and the conventions used.

## High Level Project Design Philosophy

This project is designed as a library of simple statistical learning theory algorithms applied to the analysis of DFT calculations. This libary is deeply integrated into the [STOBE](https://www.fhi.mpg.de/1022673/StoBe) dft calculation platform, and designed to be used as a tool for the analysis of transition dipole moments (TDMs) between core 1s orbitals and excited unoccupied molecular orbitals (UMOs). However, it is designed to be general purpose, and can be used to analyze any DFT calculation output, so long as it contains a set of transition dipole moments, transition energies, orbitals, atomic coordinates, and a calculated near edge X-ray absorption spectrum (NEXAFS).

## Project Structure

The project is an attempt to solidify and expand uppon the work of [@victormurcia](https://github.com/victormurcia) to develop a clustering and refinement tool within the Igor Pro programing language. This project will continue to maintain that initial tool, but will expand into a more general purpose library within the python programming language and the python ecosystem.

### Folder Structure

The project is organized into the following folders:

- **docs/**: This folder contains some simple documentation for the project. For now it is mostly just a placeholder, but will be expanded uppon in the future to contain a properly generated documentation for the project.
- **igor/**: Legacy Igor Pro implementation of the clustering and refinement tool (not the active development target). See **Igor Pro Implementation** below for WaveMetrics folder roles, file map, and a suggested reading order when tracing the pipeline.
- **notebooks/**: This folder contains the Jupyter notebooks for the project. These notebooks are used to explore the data and the algorithms, and to test and validate the code.
- **src/**: Standard layout parent directory; the **Python implementation** is entirely under **`src/dftlearn/`** (import name **`dftlearn`**). See **Python package (`src/dftlearn/`)** below.
- **streamlit_demo/**: Streamlit application that exercises core workflows against StoBe-style data. It is **not** the product surface for library users; see **Streamlit demo (`streamlit_demo/`)** below.

### Python package (`src/dftlearn/`)

**Canonical location.** All first-party Python code for this project ships from **`src/dftlearn/`** (Hatch **`packages = ["src/dftlearn"]`** in `pyproject.toml`). The **`dftrun`** CLI already lives here as **`dftlearn.cli`**. Additional tools (I/O, clustering and analysis, visualization, and shared utilities) are intended to live **in the same package tree** as further **`dftlearn.*`** subpackages or modules, not as separate install roots.

**Current layout.** Besides **`cli/`**, **`python_pipeline/`** is **staging** code ported from the Igor workflow: StoBe loading (`stobeLoader.py`), clustering and overlap (`clusteringProcs.py`, `overlapProcs.py`), step-edge handling (`stepEdgeProcs.py`), plotting helpers (`plotUtils.py`), multi-spectrum fitting (`multiSpecFitProcs.py`), and related notebooks. Expect these modules to **migrate into** dedicated **`dftlearn.*`** subpackages as APIs stabilize.

**CLI (`dftlearn.cli`, console script `dftrun`).** Entry point is `dftlearn.cli:app` (see `pyproject.toml` **`[project.scripts]`**). Top-level commands today:

- **`dftrun build`** -- generate StoBe input files from **`molConfig`**-style configuration and **XYZ** geometry (see `cli/build.py` and `cli/shared.py`).
- **`dftrun run`** -- schedule and run StoBe DFT jobs (logs, progress, packaging of outputs; see `cli/run.py` and `cli/config.py`).
- **`dftrun init`** -- placeholder for future project initialization.

Legacy or superseded generator helpers live under **`cli/_deprecated/`** (excluded from default Ruff/pyrefly paths). New work should not spread into that tree without an explicit cleanup plan.

**Where the library should grow (all under `src/dftlearn/`).** Target shape is a **small, testable core** that the CLI, notebooks, demos, and future apps import by name:

- **`dftlearn.io`** (or equivalent) -- parsers and loaders for StoBe outputs (and, where practical, other DFT stacks) into well-defined tables or arrays; path and schema handling; no Streamlit or Typer inside this layer.
- **`dftlearn.clustering`** / analysis -- filtering, overlap construction, merge and refinement aligned with the Igor reference; optional **scikit-learn**-compatible estimators or pure functions with explicit inputs and outputs.
- **`dftlearn.visualization`** (or equivalent) -- plotting and optional 3D structure views (optional **`viz`** extras in `pyproject.toml` such as **py3dmol** and **seaborn**); keep rendering separate from numerics for headless and batch use.

Until those subpackages exist, **`dftrun`** and **`python_pipeline`** remain the practical entry points; new features should add or extend **`dftlearn.*`** modules rather than new top-level packages outside **`src/dftlearn/`**.

### Streamlit demo (`streamlit_demo/`)

**Role.** A browser UI to load StoBe-related calculation outputs and drive much of the **core** interactive workflow that the Igor panel originally supported (exploration, plotting, and analysis in one place).

**What it is not.** The demo **does not** deliver a **high-quality, reusable clustering library** for consumers of **`dftlearn`**. It couples behavior to Streamlit session state and app layout, does not define a stable, documented public API for clustering and refinement, and is not a substitute for composable Python APIs, automated tests, and packaging contracts that library users expect. Contributors should assume that **library-quality interfaces live in `src/dftlearn/`** (and future subpackages as above), while **`streamlit_demo/`** is a **reference client** and **exploratory shell**.

**Lifecycle.** The project still intends to retire or replace this demo once the Python library and a clearer application or example story exist; until then it remains useful for parity checks against Igor and for manual validation.

### Igor Pro Implementation

**What this section is for.** It maps the **`igor/`** tree for contributors and agents: WaveMetrics procedure roles, the single panel entry point, how **eleven** project modules plus the **three** Chem3D procedure files (**fourteen** `.ipf` files total) are grouped, and a practical order for reading code. Use it when porting behavior to **`src/dftlearn/`**, comparing against **`streamlit_demo/`**, or answering where a step lives in Igor. The implementation is **legacy** (not the primary development target) but remains the most complete reference for the original clustering, filtering, tensor-based NEXAFS modeling, and StoBe loading workflow. Copying files into a live Igor install is described in [`igor/README.md`](igor/README.md).

#### Conventions when reading Igor code

- The **language** treats identifiers case-insensitively, but **file names** for `#include` are resolved on a case-sensitive filesystem: `SomeFile.ipf` and `somefile.ipf` are distinct paths. String-based wave names follow the same pitfall as file names when case differs.
- **State** is mostly visible in the **data browser**: waves, variables, and strings in folders are often global in practice. Many routines **side-effect** those objects instead of returning packed results; follow `GetDataFolder` / `SetDataFolder` and who last wrote each wave.
- **`SetDataFolder`** selects Igor's **logical** current folder for resolving wave names; it is not a path into the host operating system.
- **`Wave`** objects hold array data (`Make/O/N=(...) name`); use `$name` or `$"name"` to resolve dynamically named waves as needed.
- Identifiers in this tree mostly use **camelCase**, unlike typical **snake_case** in the Python package.

#### WaveMetrics: global procedures versus `#include`

Legacy layout matches [Experiments, Files and Folders](https://docs.wavemetrics.com/igorpro/igor-basics/experiments-files-and-folders). **`Igor Pro .../Igor Procedures`** (application folder and **User Files**) is scanned at startup: each procedure file there is a **global procedure file**, compiled for every experiment without an `#include`. **`.../User Procedures`** is the search path for `#include "SomeName"` (including shortcuts and subfolders); those units compile only when a parent file includes them or when you place a shortcut under **Igor Procedures**. Typical setup: put **`clusteringPanel v1.ipf`** (or a shortcut) under **`Igor Pro User Files/Igor Procedures`** and mirror this repo's **`User Procedures`** under **`Igor Pro User Files/User Procedures`**, unless you extend Igor's procedure search path to point at a checkout.

#### `igor/Igor Procedures/` (entry and global-style panel)

**What you will find here:** the one file intended to behave like a **global entry** when installed in WaveMetrics **Igor Procedures**.

```
igor/Igor Procedures/
└── clusteringPanel v1.ipf
```

**`clusteringPanel v1.ipf`** is the **UI shell**: `#include`s WaveMetrics widgets (`WaveSelectorWidget`, `PopupWaveSelector`) plus `Chem3dhooks1_1`, `Chem3dprocs1_3a`, **`FilteringMain`**, and **`StoBeImportProcs`**. It registers **`Macros → Clustering Algorithm`** → **`clusterPanel()`**, ensures **`root:Packages:DFTClustering`**, and builds the **`ClusteringAlgorithm`** window. Regions align with the pipeline stages: StoBe load (atom count, element, molecule name); optional **frame reorientation**; **experiment** NEXAFS paths and anchors; **filtering/clustering** thresholds and symmetry; **broadening** and energy grid; **tilt / angle-resolved** model and fit; **step-edge** and bare-atom controls; **cluster and transition** inspection; **overlap and OS** displays; optional **parameter-space** widgets; **`Run Algorithm`** (**`RunAlgorithmButton`**) gathering control state and calling into **`filterDFT`** and related routines. The same file holds **`clusterPanelSimple()`** (reduced UI), many control callbacks, and auxiliary plot panels (for example BIC and overlap views). **No heavy numerics here** beyond glue; algorithms live under **`User Procedures`**.

#### `igor/User Procedures/` (included libraries)

**What you will find here:** **eleven** DFT-clustering-related **`.ipf`** modules plus the **three** Chem3D procedure files (**fourteen** `.ipf` files total), all brought in through **`#include`** from the panel and from **`FilteringMain.ipf`**. Nothing here runs at Igor startup until included.

**Why WaveMetrics uses this folder:** **`#include "SomeName"`** resolves against **`User Procedures`** trees; that keeps large libraries out of the global compile set and lets the panel pull a defined dependency cone. **`FilteringMain.ipf`** is the main hub that assembles the clustering stack.

**Organization of the subsections below:** fourteen file names are split into **three** clusters by primary responsibility. Overlap is normal (for example **`SymmetryOperations.ipf`** feeds filtering decisions and tensor logic); headings state the main reason the file exists.

**Suggested reading order for one full interactive run:** (1) **`clusteringPanel v1.ipf`** control procedures, especially the run handler; (2) **`StoBeImportProcs.ipf`** when following **Load DFT**; (3) **`FilteringMain.ipf`** **`filterDFT`** and its `#include` list; (4) **`FilteringWrapper.ipf`** into **`ClusteringProcs.ipf`** / **`peakOverlap.ipf`** / **`fittingAmplitudes.ipf`**; (5) **`SymmetryOperations.ipf`** when symmetry-aware paths run; (6) **`alphaModeling.ipf`** for **`simDFT`** and tensor or experiment-facing spectra; (7) **`sequentialThresholding.ipf`** only for OS/OVP **grids**; (8) **`DFT Plotting Utilities.ipf`** and **Chem3D** when validating from the UI.

##### 1. Core clustering, filtering, and tensor-based NEXAFS modeling

**Scope:** filter DFT transitions, merge by overlap, refit amplitudes to the pre-merge spectrum, then tensor-based film and angle-resolved modeling tied to experiment.

```
igor/User Procedures/
├── alphaModeling.ipf
├── FilteringMain.ipf
├── FilteringWrapper.ipf
├── sequentialThresholding.ipf
├── ClusteringProcs.ipf
├── peakOverlap.ipf
└── fittingAmplitudes.ipf
```

- **`alphaModeling.ipf`** is the largest **numerical** module: tensor construction and normalization, film absorption after tilt and azimuthal averaging, building-block **angle-resolved** spectra (**`simDFT`** and preparation or refinement paths), and wiring clustered Gaussians to experiment-oriented outputs. After you understand **`filterDFT`**, trace modeling here.
- **`FilteringMain.ipf`** states the **pipeline contract** in its header (energy window, OS cutoff, overlap threshold, transition symmetry, clustering, amplitude refit, film tensor, step) and implements **`filterDFT`**. It is the central **`#include`** hub for the files in this cluster and several dependencies in clusters 2 and 3.
- **`FilteringWrapper.ipf`** performs **filtering** mechanics: OS scaling, transition culling, parameter and Gaussian waves, optional full overlap matrices, folder layout, and calls to **`clusteringTransitions`** or an alternate percent-difference clusterer.
- **`sequentialThresholding.ipf`** exposes **`seqThresholds`**, repeating the pipeline over **waves** of OS and overlap settings for parameter exploration without editing the panel each time.
- **`ClusteringProcs.ipf`** owns **iterative merging**: sort transitions, build summed Gaussians, **`simpleCluster2`** / **`simpleCluster3`** loops, and per-stage combined parameter waves.
- **`peakOverlap.ipf`** fills **overlap matrices** (Gaussian intersection / erf-based numerics) that drive merge decisions.
- **`fittingAmplitudes.ipf`** **re-fits** merged-peak parameters so the clustered spectrum matches the **original** DFT NEXAFS (Igor fit functions over the merged parameter wave).

##### 2. StoBe ingest, reference absorption, and transition symmetry

**Scope:** read StoBe outputs into waves, elemental and step-edge infrastructure, and symmetry classification on transition dipoles before and during filtering.

```
igor/User Procedures/
├── StoBeImportProcs.ipf
├── makeBareAtomAbsorption.ipf
└── SymmetryOperations.ipf
```

- **`StoBeImportProcs.ipf`** is the primary **StoBe I/O** layer: e.g. **`DFTwrapper1`** loads ground, excited, and transition-potential data, builds **energy corrections** and transition waves consumed downstream.
- **`makeBareAtomAbsorption.ipf`** sets up **package folders**, **element** and bare-atom tables (including **`AtomicWeight.txt`** from a user-side element library path), and UI helpers for **molecule-level** step and reference absorption under modeled spectra.
- **`SymmetryOperations.ipf`** applies **molecular symmetry** to TDMs: rotations, per-transition tensors, **isotropic / uniaxial / biaxial / triaxial** classification, and waves feeding symmetrized oscillator-strength logic in the filter.

##### 3. Diagnostics plotting and molecular structure viewing

**Scope:** optional graphs for validating stages and third-party **3D** viewing; not required to understand the numerical pipeline in isolation.

```
igor/User Procedures/
├── DFT Plotting Utilities.ipf
├── Chem3Dhooks1_1.ipf
├── Chem3Dprocs1_3.ipf
└── Chem3Dprocs1_3a.ipf
```

- **`DFT Plotting Utilities.ipf`** builds **stage-comparison** graphs (for example OS versus energy for original, first-filtered, and merged-cluster parameters).
- **`Chem3Dhooks1_1.ipf`**, **`Chem3Dprocs1_3.ipf`**, and **`Chem3Dprocs1_3a.ipf`** are the **Chem3D** package by **Richard Knochenmuss** (independent Igor module): Gizmo-style **3D** structures, experiment hooks, menus. They support **geometry** next to spectra, not StoBe file parsing. The panel **`#include`s** the **`1_3a`** file; **`1_3`** is the sibling variant from the same lineage.


<!-- DO NOT EDIT THIS BLOCK IT IS MANAGED BY DOTAGENTS -->

# General

## Audience

You are assisting someone who holds a physics PhD and has extensive experience writing software for scientific and engineering work: numerical computing, data analysis, instrumentation, simulation, and research-grade reproducibility expectations. Assume strong mathematical literacy, comfort with linear algebra and statistics, and low tolerance for hand-wavy numerics or silent type coercion in scientific code.

## Operating principles

- Prefer the smallest coherent change set that satisfies the stated specification. Avoid drive-by refactors, unrelated formatting sweeps, and scope expansion.
- Treat the repository’s existing patterns as the default contract. Match naming, module boundaries, error-handling style, and test layout unless the user explicitly requests a migration.
- Default to production-grade output: complete, runnable, and reviewable. Do not ship placeholder text such as ellipses, “the rest of the implementation here”, “TODO: implement”, or “fill in later” inside code or patches unless the user explicitly authorizes a stub.
- Remain non-lazy: if a command fails, diagnose, adjust, and retry with a different approach when reasonable. Do not stop after the first error without analysis.
- Do not use emoji in code, comments, documentation strings, commit messages, or user-facing text unless the user explicitly requests emoji.

## Communication and documentation outside code

- Avoid standalone documentation files or long narrative write-ups unless the user asks for them or the repository already uses them for the same purpose.
- Prefer editing the code and tests that enforce behavior over adding parallel prose that can drift out of date.
- When the user asks for explanation, keep it precise and tied to the change set.

## Public API documentation (language-agnostic)

- Every **public** function, method, or type exported from a library module carries documentation appropriate to the language (for example Python docstrings, Rust `///` on public items, TSDoc/JSDoc on exported symbols).
- Documentation states the **surface**: name, purpose, parameters, return value, and thrown or returned error shapes when that is part of the contract.
- Prefer **prescriptive** voice that states what the symbol **does** and **means** for callers. Prefer “Maximum `foo` grouped by `bar` using a stable sort on `bar`.” over “Returns the max of foo by bar.” or "Computes the maximum `foo` grouped by `bar` using a stable sort on `bar`." or "Returns the maximum `foo` grouped by `bar` using a stable sort on `bar`."
- For each parameter: name, type as used in the project, allowed ranges or invariants when non-obvious, and interaction with other parameters.
- For results: type, semantics, units when relevant, ordering guarantees, and stability promises when they matter for science or reproducibility.
- Describe **what** the function does at the abstraction level of the API, **how** only when algorithmic choices affect correctness, performance contracts, or numerical stability, and **why** that approach is chosen when trade-offs are non-obvious (for example streaming vs materializing, online vs batch statistics).
- **Internal** helpers omit the long-form contract unless complexity warrants a short note. **Private** helpers keep documentation minimal (a phrase or single sentence at most).

## Module and package documentation

- Each library module (or the closest equivalent in the language’s module system) includes a short module-level description of responsibility: what problems it solves, what it explicitly does not handle, and which invariants callers should respect.
- Module docs are prescriptive about intent and boundaries so new contributors and agents do not duplicate concerns across modules.

## Tooling, skills, and continued learning

- Use project rules, agent skills, and MCP documentation tools when they apply to the task. Prefer authoritative library and framework documentation over memory when behavior, defaults, or breaking changes matter.
- When touching unfamiliar APIs, verify signatures, deprecations, and error modes against current docs or source in the dependency before guessing.
- Prefer file-scoped or package-scoped commands when the repository documents them (typecheck, lint, format, test on a single path) to shorten feedback loops.
- State permission-sensitive actions clearly (dependency installs, destructive commands, credential access) and follow the user’s safety expectations for the workspace.

## Task shape (goal, context, constraints, completion)

- Restate the goal in terms of observable outcomes: behavior, tests, or interfaces that change.
- Ground work in the relevant files, modules, and existing tests named by the user or discovered through search.
- Honor explicit constraints (performance, numerics, compatibility, style) before proposing alternatives.
- Stop when the completion criteria are met: tests pass where applicable, edge cases called out by the user are handled, and no placeholder implementation remains.

## Quality bar for agent output

- Do not substitute templates, pseudocode, or abbreviated implementations when the user asked for working code.
- If scope is too large for one pass, propose a staged plan and complete the first stage fully rather than leaving partial files full of omissions.

# Python

## Python

The following applies to **Python** work in this repository: scientific and general-purpose code, with emphasis on clear structure, reproducible tooling, and documentation that matches how the team uses Cursor (skills, subagents, and editor rules).

### Conventions

- Follow **PEP 8** surface style; treat **Ruff** configuration in `pyproject.toml` as the enforced interpretation of those conventions.
- Prefer **Python 3.12+** unless the repo pins an older interpreter.
- Prefer **readability** over micro-optimizations; vectorize numerics when a library primitive exists instead of tight Python loops over large data.
- **NumPy-style docstrings** on **public** APIs (parameters, returns, and examples where they clarify behavior). Keep implementation bodies clear **without** long narrative **inline comments**—use names, small functions, and docstrings instead.
- **Tables and time series**: explicit **index/column** semantics when using **pandas**; **lazy** queries when standardizing on **Polars** for heavy pipelines.
- **Lab / instruments**: separate **resource lifecycle** (open, configure, close) from **command strings** and parsing (e.g. **pyvisa** patterns).

### Tooling

Use **[uv](https://docs.astral.sh/uv/latest/)** for environments, runs, and dependency changes. Pair it with the **[Astral](https://astral.sh/)** stack as configured in this project.

- **Dependencies**: add, upgrade, and remove packages with **`uv add`**, **`uv add … --upgrade`**, **`uv remove`**—do **not** hand-edit version pins in `pyproject.toml`.
- **Environment**: **`uv sync`** after cloning or when the lockfile changes; **`uv run …`** to execute Python, tools, and tests.
- **Dev tools**: keep **`ruff`** and **`ty`** in the development (or project) dependency group; run **`ruff check`** (and project formatting if applicable) plus **`ty check`** on changed code. **`uvx`** remains an option for one-off tool runs.
- **pytest**: install via **`uv add --dev pytest`** (or the project’s dev group); run with **`uv run pytest`**.

If a **`uv`** subcommand differs by version, use **`uv --help`** or the [uv docs](https://docs.astral.sh/uv/latest/).

### Testing

- Prefer **fast, deterministic** unit tests; isolate I/O and timing-sensitive checks when the team uses markers or separate jobs.
- **Regression tests** for fixed bugs; for numerics, assert **shapes**, **dtypes**, and stability expectations when science or reproducibility requires it.

### Cursor: skills

Load these **skills** by **name** when the task matches (each skill’s own `SKILL.md` and references hold the full detail). Installed skills usually live under `.cursor/skills/` (or your editor’s equivalent).

| Skill | Use it for |
|-------|------------|
| **general-python** | Hub: **uv** / **ruff** / **ty** workflow, builtins and collections, functions and classes, **dataclasses**, typing boundaries, **pytest**, scientific defaults, and pointers to the other skills. |
| **numpy-scientific** | **NumPy**: dtypes, views vs copies, broadcasting, ufuncs and reductions, **linalg** / **einsum**, **`Generator`**, I/O, interop with tables and plotting. |
| **dataframes** | **pandas** and **Polars**: when to use which, indexing, joins, lazy execution, I/O, nulls, Arrow interop. |
| **numpy-docstrings** | **Numpydoc**-style docstrings: section order, semantics (what belongs in docstrings vs types vs tests), anti-patterns, **Parameters** / **Returns** / **Examples** / classes / modules. |
| **matplotlib-scientific** | Publication-style **Matplotlib**: OO API, axes and legends, layout, export, journal widths, optional **SciencePlots**. |
| **lab-instrumentation** | **PyVISA** / VISA sessions, **sockets** vs VISA, **hardware abstraction**, **input validation** before I/O, **testing** without hardware, **PDF** extraction for datasheets and manuals. |

### Cursor: subagents

Delegate by **subagent name** when a focused pass is better than inline editing. Subagents usually live under `.cursor/agents/` (or your editor’s equivalent).

| Subagent | Use it for |
|----------|------------|
| **python-reviewer** | Reviewing changes: **uv** hygiene, typing, numerics footguns, tests, docstring quality. |
| **python-types** | Deep **typing** for **ty**: annotations, PEP 695-style generics, exhaustive **`match`**, fixing checker output. |
| **python-refactor** | **Structure**: unclear multi-value returns, composition vs inheritance, oversized functions or classes, deterministic boundaries. |

### Cursor: rules

- A **Python** Cursor **rule** applies to Python sources (typically `**/*.py` when the rule is configured for those globs). It restates **interpreter preference**, **uv** usage, **ruff** / **ty** expectations, numerics and docstring defaults, and points to **general-python**, domain skills such as **lab-instrumentation** when editing drivers or lab I/O, and the subagents above.
- **Rule text is authoritative for “always on” editor hints**; **skills** carry the long-form patterns and examples. When the two differ on a detail, follow **this spec** and **`pyproject.toml`**, then the **rule**, then skill nuance.

### External references

- [numpydoc format](https://numpydoc.readthedocs.io/en/latest/format.html)
- [uv](https://docs.astral.sh/uv/latest/), [Ruff](https://docs.astral.sh/ruff/), [ty](https://docs.astral.sh/ty/)

# Python - Jupyter

This workspace extends Python with **Jupyter notebook** expectations.

### General Guidelines
This project will make strong use of jupyter notebooks to solve a number of problems, primarily with a scientific focus, but may also have a general purpose focus. Notebooks should be written and well documented. But keep in mind that a notebook allows for a lot of flexibility, and as such, the code should be written in a way that is easy to understand and maintain. As a general rule, we will use notebooks for one of the following purposes:

#### Lightweight Exploration Notebooks
Here we will have a few cells that will build the idea or concept, test it with some data, and then present the results in a clean and easy way. Use matplotlib for static plotting, or hvplot/altair/plotly for interactive plotting.

It is important to note that these notebooks are not really designed to be robust, and as such we should not focus too hard on making them so. Write the minimal code to get the job done, and then move on to the next notebook.

#### Prototyping Notebooks
This is where we will prototype a robust coding solution to a problem. We will use this to test and validate each chunk of code in the complex workflow. As such it is important that we use many small cells to test and validate each chunk of code. Eventually, we will want to move this code into a final script/library.

It is important to note that these notebooks, are not really designed to showcase the code, but rather to test and validate each chunk of code in the complex workflow. As such, we should not focus too hard on making them look nice, but rather ensure that we atomize the code into small cells that can be tested and validated individually. Testing and validation might be done using simple displays, or plotting, or other cases. But assertstetements are not necessary, and should be avoided if possible in favor of displaying the results to the user.

#### Demonstration Notebooks
These notebooks are designed to showcase a workflow of production ready code. Ideally, after a library is complete and ready to use, the user will be able to import the library and use the code treating the notebook as a production environment. They will have a minimal ammount of cells, mixing in documentation and examples of how to use the library.

It is important to note that these notebooks are designed to be robust. These should mix in a healthy ammount of markdown documentation and explaination. but not be too heavy handed. Keep in mind that the goal of these notebooks is that a user can copy them into their own notebooks and know how to use the library.

### Use of Cells
- Use markdown cells sparingly to explain the code, or why it works the way it does. But avoid long narrative markdown cells that are overly robust.
- Use code cells to write the code. Ensure that each cell is small and atomic. Each cell should have one responsibility and deterministically produce a result. If you define a function, call it in the same cell. Avoid using global variables, or mutable state if possible.
- Class definitions should idealy be avoided in notebooks, unless we are prototyping a library. If we need a class, it should be defined in a separate cell and have a minimal footprint.
- Variables should be defined in the cells that use them, and should usually be displayed to the user. Avoid using global variables, or mutable state if possible.
- When prototping a library, use the `%autoreload 2` magic to ensure that the code is reloaded when it is changed.
<!-- END OF MANAGED SECTION -->

## Learned User Preferences

- For RDKit 2D structure layouts in this repository, use `AllChem.Compute2DCoords` rather than `Compute2DMCoords`, which is not consistently available across RDKit builds.

## Learned Workspace Facts

- `dftrun postprocess` accepts a StoBe run root directory and writes `packaged_output/xray_spectra_long.csv` (long-form spectra with columns `energy_ev`, `abs`, and `site`) plus `xas_site_summary.png` by aggregating per-site `XrayT*.out` tables together with geometry from an XYZ file resolved like `dftrun build`.
- XYZ parsing in `dftlearn` supports optional standard two-line XYZ headers or headerless StoBe-style coordinate lines; row order defines atom indices for site labels such as `C1` and `C2`, while tokens such as `C01` encode element type without relying on the comment line.
