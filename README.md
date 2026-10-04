# dft-learn

[![PyPI](https://img.shields.io/pypi/v/dft-learn?style=flat-square&logo=pypi&logoColor=white&label=PyPI)](https://pypi.org/project/dft-learn/)
[![Python](https://img.shields.io/badge/python-3.12%2B-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/downloads/)
[![CI](https://img.shields.io/badge/CI-GitHub%20Actions-2088FF?style=flat-square&logo=githubactions&logoColor=white)](https://github.com/WSU-Carbon-Lab/dft-learn/actions)
[![arXiv](https://img.shields.io/badge/arXiv-2509.01734-b31b1b?style=flat-square&logo=arxiv&logoColor=white)](https://arxiv.org/abs/2509.01734)
[![Hugging Face](https://img.shields.io/badge/Hugging%20Face-carbon--lab-FFD21E?style=flat-square&logo=huggingface&logoColor=black)](https://huggingface.co/carbon-lab)

Python tools for StoBe-style core-level spectra: parse transitions, cluster by
Gaussian peak overlap, and build bond-traceable resonant X-ray optical tensors
for angle-resolved NEXAFS (RSoXS, XRR).

[Carbon Lab](https://labs.wsu.edu/carbon/) ·
[Hugging Face](https://huggingface.co/carbon-lab) ·
[X-ray Atlas](https://xrayatlas.wsu.edu/) ·
[arXiv](https://arxiv.org/abs/2509.01734) ·
[PRL](https://doi.org/10.1103/rfgg-ffyz)

## Install

```bash
pip install -U dft-learn          # or: uv add dft-learn
uv tool install dft-learn         # optional CLI: dftrun
```

Requires Python 3.12+.

## Quick start

```python
import numpy as np
from dftlearn.clustering import TransitionSticks, cluster_by_overlap

sticks = TransitionSticks(
    energy_ev=np.array([284.0, 284.15, 295.0]),
    oscillator_strength=np.array([1.0, 0.9, 0.5]),
    sigma_ev=np.full(3, 0.2),
    site=np.array(["C1", "C1", "C2"]),
    os_xx=np.array([0.1, 0.1, 0.0]),
    os_yy=np.zeros(3),
    os_zz=np.array([0.2, 0.2, 0.5]),
)
result = cluster_by_overlap(sticks, overlap_threshold=50.0)
print(result.n_iterations, result.energy_ev.shape[0])
```

StoBe workflow CLI:

```bash
dftrun build --help
dftrun run --help
dftrun postprocess --help
```

## Demos

| Space | Dataset |
|-------|---------|
| [CuPc optical model](https://huggingface.co/spaces/carbon-lab/cupc-optical-model) (reference) | [optical-cupc](https://huggingface.co/carbon-lab/optical-cupc) |
| [ZnPc optical model](https://huggingface.co/spaces/carbon-lab/znpc-optical-model) | [optical-znpc](https://huggingface.co/carbon-lab/optical-znpc) |

CuPc walks through the publication workflow: sites, DFT sticks (isotropic / xx / zz),
clusters, and refinement against experiment.

## Package layout

| Area | Role |
|------|------|
| `dftlearn.io` | StoBe / XYZ parsers |
| `dftlearn.clustering` | Overlap matrices, merge, OS elbow, threshold selection |
| `dftlearn.xas` | Spectrum reconstruction, C3 symmetry helpers |
| `dftlearn.visualization` | Optional figures (`viz` extras) |
| `dftrun` | Build inputs, schedule runs, package spectra |
| [`igor/`](igor/README.md) | Legacy Igor reference (not the install target) |

<details>
<summary><strong>Development</strong></summary>

```bash
git clone https://github.com/WSU-Carbon-Lab/dft-learn.git
cd dft-learn
make install && make verify    # ruff + format check + pytest
```

| Target | Action |
|--------|--------|
| `make test` | pytest |
| `make lint` | ruff check |
| `make type-check` | ty (advisory) |
| `make fix` | ruff check --fix + format |
| `make build` | sdist + wheel |

Conventions: [`AGENTS.md`](AGENTS.md).

**Release (Trusted Publishing).** Configure a GitHub Environment named `pypi`
for the PyPI project [`dft-learn`](https://pypi.org/project/dft-learn/), bump
the version in `pyproject.toml`, then:

```bash
git tag v0.2.0 && git push origin v0.2.0
```

The [Release](https://github.com/WSU-Carbon-Lab/dft-learn/actions/workflows/release.yml)
workflow publishes to PyPI and creates a GitHub Release.

</details>

## Igor → Python port

Status of the [`igor/`](igor/) clustering pipeline in `dftlearn`
([#1](https://github.com/WSU-Carbon-Lab/dft-learn/issues/1) ·
[#2](https://github.com/WSU-Carbon-Lab/dft-learn/issues/2) ·
[#3](https://github.com/WSU-Carbon-Lab/dft-learn/issues/3)).
Audited against `FilteringMain` / `filterDFT` and current library + `python_pipeline/` staging.

Summaries show one box per capability and **done / total**
(done = tested or done-untested; out of scope omitted from the total).

<img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> tested ·
<img src="docs/readme/box-done.svg" alt="■" width="12" height="12"> done, not tested ·
<img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"> in progress ·
<img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"> not started ·
<img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"> out of scope

<details>
<summary><strong>Ingest</strong> · <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"> 3/4</summary>

| | Capability | Completion | Location | Track |
|:-:|---|---|---|---|
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | `XrayT*.out` / TP `*.xas` sticks | tested | `dftlearn.io` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | XYZ geometry & site labels | tested | `dftlearn.io` |  |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Delta-KS / FINAL ENERGY / TP LUMO (`E^c`) | tested | `dftlearn.io` |  |
| <img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"> | Full ground / excited / TP directory load | in progress | `python_pipeline/stobeLoader` | [#2](https://github.com/WSU-Carbon-Lab/dft-learn/issues/2) |

</details>

<details>
<summary><strong>Filter & cluster</strong> · <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"><img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"> 5/8</summary>

| | Capability | Completion | Location | Track |
|:-:|---|---|---|---|
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Energy window + OS% cull | tested | `dftlearn.clustering` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Gaussian peak-overlap matrices | tested | `dftlearn.clustering` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Iterative overlap merge (`simpleCluster3`) | tested | `dftlearn.clustering` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | OS% elbow cutoff | tested | `dftlearn.clustering` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | OVP threshold selection (BIC / GP) | tested | `dftlearn.clustering` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"> | OS × OVP sequential grids (`seqThresholds`) | not started | — |  |
| <img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"> | Amplitude refit to pre-merge DFT NEXAFS | not started | — | [#1](https://github.com/WSU-Carbon-Lab/dft-learn/issues/1) |
| <img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"> | `filterDFT`-style end-to-end orchestration | in progress | `dftrun postprocess` (partial) | [#1](https://github.com/WSU-Carbon-Lab/dft-learn/issues/1) |
| <img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"> | Alternate %-difference clusterer (`pDiff`) | out of scope | — |  |

</details>

<details>
<summary><strong>Symmetry</strong> · <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"> 1/3</summary>

| | Capability | Completion | Location | Track |
|:-:|---|---|---|---|
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | C3 dipole fold / site OS from tensors | tested | `dftlearn.xas` |  |
| <img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"> | General n-fold TDM symmetry (iso / uni / bi / tri) | not started | — | [#1](https://github.com/WSU-Carbon-Lab/dft-learn/issues/1) |
| <img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"> | TDM frame reorientation (Euler on transitions) | not started | — |  |

</details>

<details>
<summary><strong>Tensors & experiment</strong> · <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"> 2/5</summary>

| | Capability | Completion | Location | Track |
|:-:|---|---|---|---|
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | TP XAS reconstruction & broadening schedule | tested | `dftlearn.xas` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Molecular dipole / Cartesian OS tables | tested | `dftlearn.xas` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"> | Film tilt / `simDFT` / angle-resolved model | in progress | `python_pipeline/multiSpecFitProcs` | [#3](https://github.com/WSU-Carbon-Lab/dft-learn/issues/3) |
| <img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"> | Bare-atom / Henke step edge | in progress | `python_pipeline/stepEdgeProcs` | [#3](https://github.com/WSU-Carbon-Lab/dft-learn/issues/3) |
| <img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"> | Multi-spectrum experiment fit | in progress | `python_pipeline/multiSpecFitProcs` | [#3](https://github.com/WSU-Carbon-Lab/dft-learn/issues/3) |

</details>

<details>
<summary><strong>Visualization & UI</strong> · <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-done.svg" alt="■" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"><img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"> 2/3</summary>

| | Capability | Completion | Location | Track |
|:-:|---|---|---|---|
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Cluster & reconstruction report figures | tested | `dftlearn.visualization` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-done.svg" alt="■" width="12" height="12"> | Site / SCF / orbital summary figures | done | `dftlearn.visualization` |  |
| <img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"> | Interactive clustering panel parity | not started | HF Space / demos | [#5](https://github.com/WSU-Carbon-Lab/dft-learn/pull/5) |
| <img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"> | Chem3D molecular viewer | out of scope | — |  |
| <img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"> | Igor panel chrome / wave selectors | out of scope | — |  |

</details>

Overall: <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"><img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-done.svg" alt="■" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"><img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"> **13/23**

## Citation

```bibtex
@article{murcia2025opticaltensors,
  title   = {Quantitative and bond-traceable resonant X-ray optical tensors of organic molecules},
  author  = {Murcia, Victor and Alqahtani, Obaid and Heilman, Harlan and Collins, Brian A.},
  journal = {Phys. Rev. Lett.},
  year    = {2025},
  doi     = {10.1103/rfgg-ffyz},
  eprint  = {2509.01734},
  archivePrefix = {arXiv}
}
```

[`CITATION.cff`](CITATION.cff) · [lab publications](https://labs.wsu.edu/carbon/publications/)

## License

MIT ([`LICENSE`](LICENSE)). Maintainer: [Harlan Heilman](mailto:harlan.heilman@wsu.edu).
