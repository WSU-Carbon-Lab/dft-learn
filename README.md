# dft-learn

[![PyPI](https://img.shields.io/pypi/v/dft-learn?style=flat-square&logo=pypi&logoColor=white&label=PyPI)](https://pypi.org/project/dft-learn/)
[![Python](https://img.shields.io/badge/python-3.12%2B-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/downloads/)
[![CI](https://img.shields.io/badge/CI-GitHub%20Actions-2088FF?style=flat-square&logo=githubactions&logoColor=white)](https://github.com/WSU-Carbon-Lab/dft-learn/actions)
[![arXiv](https://img.shields.io/badge/arXiv-2509.01734-b31b1b?style=flat-square&logo=arxiv&logoColor=white)](https://arxiv.org/abs/2509.01734)
[![Hugging Face](https://img.shields.io/badge/Hugging%20Face-carbon--lab-FFD21E?style=flat-square&logo=huggingface&logoColor=black)](https://huggingface.co/carbon-lab)

**dft-learn** is a Python library for analyzing DFT / StoBe-style core-level spectra:
filter transitions, cluster by peak overlap, and build bond-traceable resonant X-ray
optical tensors for angle-resolved NEXAFS and optical-constant work
(RSoXS, XRR).

Website / org: [huggingface.co/carbon-lab](https://huggingface.co/carbon-lab)·
Lab: [labs.wsu.edu/carbon](https://labs.wsu.edu/carbon/) ·
Atlas: [xrayatlas.wsu.edu](https://xrayatlas.wsu.edu/)

---

## Installation

```bash
pip install -U dft-learn
```

Or with [uv](https://docs.astral.sh/uv/) (recommended):

```bash
uv add dft-learn
# CLI
uv tool install dft-learn
dftrun --help
```

Requires Python 3.12+.

## Quick start

```python
import dftlearn

print(dftlearn.__version__)
```

Build StoBe inputs and schedule runs with the CLI:

```bash
dftrun build --help
dftrun run --help
dftrun postprocess --help
```

## Interactive demos

| Role | Demo | Data |
|------|------|------|
| **Primary** | [CuPc optical model](https://huggingface.co/spaces/carbon-lab/cupc-optical-model) | [optical-cupc](https://huggingface.co/carbon-lab/optical-cupc) |
| Demo result | [ZnPc optical model](https://huggingface.co/spaces/carbon-lab/znpc-optical-model) | [optical-znpc](https://huggingface.co/carbon-lab/optical-znpc) |

The CuPc Space is the reference walkthrough for the publication workflow
([arXiv:2509.01734](https://arxiv.org/abs/2509.01734) /
[PRL](https://doi.org/10.1103/rfgg-ffyz)): molecule sites, DFT sticks
(isotropic / xx / zz), clusters, and refinement next to experiment.

## Library goals

- **Composable Python APIs** under `dftlearn` for I/O, clustering / overlap, and analysis
- **scikit-learn-friendly** estimators and plain functions with explicit inputs / outputs
- **`dftrun`** for StoBe input generation, job scheduling, and spectrum packaging
- Headless-friendly numerics; visualization stays optional (`viz` extras)

Legacy Igor Pro procedures that inspired the clustering pipeline live under
[`igor/`](igor/README.md) for reference only — they are not the install target.

## Development

```bash
git clone https://github.com/WSU-Carbon-Lab/dft-learn.git
cd dft-learn
make install
make verify          # ruff + format check + pytest
```

Useful targets:

```bash
make test
make lint
make type-check      # ty (advisory while the tree is typed incrementally)
make fix             # ruff check --fix + format
make build           # sdist + wheel via uv
```

Contributor conventions: [`AGENTS.md`](AGENTS.md).

### Releasing to PyPI

CI runs on every push / PR. Publishing uses
[Trusted Publishing](https://docs.pypi.org/trusted-publishers/) (OIDC):

1. Configure a GitHub Environment named `pypi` linked to the PyPI project
   [`dft-learn`](https://pypi.org/project/dft-learn/).
2. Bump the version in `pyproject.toml`.
3. Tag and push: `git tag v0.2.0 && git push origin v0.2.0`.

The [Release](https://github.com/WSU-Carbon-Lab/dft-learn/actions/workflows/release.yml)
workflow builds the sdist/wheel, uploads to PyPI, and creates a GitHub Release.

<details>
<summary><strong>Igor → Python conversion</strong></summary>

Port status for the [`igor/`](igor/) clustering pipeline → `dftlearn`.
See also [#1](https://github.com/WSU-Carbon-Lab/dft-learn/issues/1) ·
[#2](https://github.com/WSU-Carbon-Lab/dft-learn/issues/2) ·
[#3](https://github.com/WSU-Carbon-Lab/dft-learn/issues/3).

`✅` tested &nbsp;·&nbsp; `☑️` implemented &nbsp;·&nbsp; `🔄` in progress &nbsp;·&nbsp; `⬜` not started &nbsp;·&nbsp; `➖` out of scope

#### Ingest

| | Capability | Python | Track |
|:-:|---|---|---|
| ✅ | `XrayT*.out` / XAS sticks | `io.xray_out`, `io.stobe_xas_sticks` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| ✅ | XYZ geometry & site labels | `io.xyz_structure` | |
| ☑️ | Aligned reconstruction tables | `xas.spectrum`, `dftrun postprocess` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| 🔄 | Full ground / excited / TP load | `python_pipeline.stobeLoader` | [#2](https://github.com/WSU-Carbon-Lab/dft-learn/issues/2) |

#### Clustering & filtering

| | Capability | Python | Track |
|:-:|---|---|---|
| ✅ | Peak-overlap matrices | `clustering.overlap` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| ✅ | Iterative overlap merge | `clustering` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| ✅ | OS% elbow cutoff | `clustering.os_elbow` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| ✅ | Overlap-threshold selection | `clustering.selection` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| 🔄 | `filterDFT` orchestration | — | [#1](https://github.com/WSU-Carbon-Lab/dft-learn/issues/1) |
| 🔄 | OS × OVP parameter grids | `clustering.selection` | |
| ⬜ | Amplitude refit (pre-merge) | — | [#1](https://github.com/WSU-Carbon-Lab/dft-learn/issues/1) |

#### Symmetry, tensors & experiment

| | Capability | Python | Track |
|:-:|---|---|---|
| ✅ | C3 dipole fold / site OS | `xas.c3_symmetry` | |
| 🔄 | General TDM symmetry classes | — | [#1](https://github.com/WSU-Carbon-Lab/dft-learn/issues/1) |
| 🔄 | Film tensors / `simDFT` / tilt | `python_pipeline.multiSpecFitProcs` | [#3](https://github.com/WSU-Carbon-Lab/dft-learn/issues/3) |
| 🔄 | Bare-atom / Henke step edge | `python_pipeline.stepEdgeProcs` | [#3](https://github.com/WSU-Carbon-Lab/dft-learn/issues/3) |
| 🔄 | Multi-spectrum experiment fit | `python_pipeline.multiSpecFitProcs` | [#3](https://github.com/WSU-Carbon-Lab/dft-learn/issues/3) |
| ☑️ | Cluster / site / orbital figures | `visualization` | |

#### UI

| | Capability | Python | Track |
|:-:|---|---|---|
| ⬜ | Interactive clustering panel | HF Space / demos | [#5](https://github.com/WSU-Carbon-Lab/dft-learn/pull/5) |
| ➖ | Chem3D viewer | — | |

</details>

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

See also [`CITATION.cff`](CITATION.cff) and the lab
[publications list](https://labs.wsu.edu/carbon/publications/).

## License

MIT — see [`LICENSE`](LICENSE). Maintainer: Harlan Heilman
\<harlan.heilman@wsu.edu\>.
