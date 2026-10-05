# dft-learn

[![PyPI](https://img.shields.io/pypi/v/dft-learn?style=flat-square&logo=pypi&logoColor=white&label=PyPI)](https://pypi.org/project/dft-learn/)
[![Python](https://img.shields.io/badge/python-3.12%2B-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/downloads/)
[![PRL](https://img.shields.io/badge/PRL-137%2C%20158001-003087?style=flat-square&logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAHBUlEQVR42nWXW2wdVxWGv7X3zDknvhwTJy6O47ZporhO45SYQNskUmpCXxASCRJq1YJUQQFFqlQQAvEGNAhV9AXxgoSoEFILbR6A8kAloEWK1AuFSBWkuTRuU4XE5OLYjm%2FnNjN78bBnzsw5dkcazZw9M2v969%2F%2F%2Bvc%2BwjeOg%2FAIhuMYdmIAETAiGAUjGCuogLoE4ghEIQwgtGAFjPgxHKgDlwDqn4cBWIMEVq0xOFEUnVbhBygnAlQfAXnJZ1BFEQRQBRUMgqs3wQgjtw1y8K5R7hsdYfy2TYxWqwxUyohAK3HU4oj5RoMry8tcXFzg%2FPwc%2F745y8XFW7QadYmTGEqhEgZjYs1LPUGI8MTT72FkDCMRoiFGFKOIMaKiELU4ML6dbx26jyO7xygHAQBnb9zkP9dvMLO8zGKziUPZEAZsrFQYqfaztb%2BfLX19DJTLLLeaXFhY0DevXuG1y5fkH9euRrXaSoiLLwhfe9ph8JRLllx8coFnjz7E9w7vB2B2pcazJ9%2FixdNnmVm8BXEMLvZTYMTTbkipDwnLJW4fqLJvyzAHtm7l4Miojg0OEjknp25c1T99MI3w1eMOo%2BIrFxBELGgc86vHj%2FD1%2FZMAvHrhQx773cvMzs1BpQyhwViD8d%2BQa0dxgCPThHqQ6qAUsrla1ak77uDIjp3y2Tu3qQcgGQDEWkPSaPDlB%2FfxwuNHATh3%2FSaf%2BtmvqdVqhD1lYhQVCok1TU7h6k8xYIxBBGLUCzSOFFTCvl4NUPUlqAoqJImDconvHH6A7Hjmb29QW1giHOglipNU9fir05x%2B5xnApV1hBHVKookHJp5eUw4EgagVS4CmL2sq%2Fihm151b2DMyBMBKs8Wr730IpRJx7HwgVVLQ%2FroeiGxcUrCpplQg0WzcEGTtBmAQkjhhx%2BaNhNYCMLdaZ2G1ASI%2BX3elhvVBSOE36Tuo%2FynSBhSg2TspEAcla9r0B8YgkCYoBsqSpVPRDeIjGEA6gRgPwAdU5xMsrDbaAAZ7Kgz19kDiEE0Tahpc02lw2o6BSxlw2nVKYZz2vfEfpqBUwRimr81Ta0UAbCiFTN4%2BDFGMEZOjXwOCPJEW3lHNzyKQtEtNHkg9A9ZyZX6JM%2F%2BbbbPw6J57oAFaTyBWiqDXJu2utnDffp6DMe1A6UMrBlbr%2FOX0BwBcbq1yfQSCA0O4sR7oD6DpIMmqLyZdj3Y%2BYlr8d0FHSznBiYIEnPjXWZYnK%2FzixnlWkiZycABxoE0H06twahFiB4HxILLWg1yExQ7IhjJRarsLiqoGTRRKIe%2Fuinj32jvgLEYsrum8swUCk1W0GsBrc5C41JCkU%2Bnr3WuhI8i6wBXUC57eT%2FRi7u3HaogxXiaSFqeq6GqCbO%2BBiX5oafYgnZKCHjpop1OIacsXNCB%2BXssG9vShLYdD%2FaKSmQcgmddHDu7uhYr1wmzPdVHpWXEFIEqHFkwbFUAMVAN%2FJrk%2FZVMrIr5V8S0kPca%2FGxcDs1aE2tn7RVYKGvBBCSTX0jpHBwgRvyVLgIBO4Wm3DjqLyRzRWzEpUhSWEogUKpLHUvXUdx3GgVtO0OL3ReF1A6EbjBScMPP0xRj%2B24SSDyq5%2BtpgDIKGQnJpGZ1vIRavn%2FXstmg8XfPvRSiy1tneuIXWHVIS35Ztn1KMCi6JKKnh57sOMdTTizZjDLmj5kAorB3rGJEXIfka79Qby%2BUm%2FPGmB9FrkUAQK0hgcBuEj%2FX28dzIAZ6a2Mtvjh0FMbg46QTRIbyCQFU6mLGMT%2F2oPW%2BZIAMD11owXffUWvHL%2B1KCvrPAd3t28%2B2JvdRaEfeMDDE2vJk%2FvHUGNeKnLBMpBTPMdiPFMYWAdfYMOIWSgdkIXpmDDRYNxO9kGgnPvP53Pr1xiC9M3k29FfHYAxOsNlt885cvo%2BUAYwXntMPx8vtiF4Bl%2FMEf5rBF2kyor5wgpTVORV62uDjh96fOM7VrGzs%2BPkijFXP%2Fjq0MD1b58z%2FPoWK8g3Zs3bJCBTRtSEVzEfr20LxX08RJep%2BxFzmMDYjqEZ%2F76W95%2B%2F0ZKqWARhRz7DP7eO7YFyFKcFGCweR7hrb9qhY7w6TJ1e%2BO6QQhxc0G7VZyicMGAatLdR76yfO8%2Ff4VKmFAoxXzxKFJTjz1MKjgmjEWk7Zo%2BmfBX%2F2yq4IBptN0UQ5CNdfF%2BotLEjtsGLKyVOPwj5%2Fn5LlLbSYevn83f%2F3%2BV%2BiplElqTb%2FHcKqpPQtOIg9Cpy3jU7PAlxCxbSZEpMM6RboWBv9QnWIDS6ve4oU3zzC5bZjdo0O04oSx4U18fu9OXjl9kVtzi5gwFNXsQ2xa4JOW8akzwHngXmBTQROyZn8hhRUuBaJOMYEliWJefP0024cH%2BeRdW2jGMaODVR7dP8HJ6Rlmri%2BoWJOxOQ08CXLi%2FwziXzwZxpbOAAAAAElFTkSuQmCC)](https://journals.aps.org/prl/abstract/10.1103/rfgg-ffyz)
[![arXiv](https://img.shields.io/badge/arXiv-2509.01734-b31b1b?style=flat-square&logo=arxiv&logoColor=white)](https://arxiv.org/abs/2509.01734)
[![Hugging Face](https://img.shields.io/badge/Hugging%20Face-carbon--lab-FFD21E?style=flat-square&logo=huggingface&logoColor=black)](https://huggingface.co/carbon-lab)

Python tools for StoBe-style core-level spectra: parse transitions, cluster by
Gaussian peak overlap, and build bond-traceable resonant X-ray optical tensors
for angle-resolved NEXAFS (RSoXS, XRR).

[Carbon Lab](https://labs.wsu.edu/carbon/) ·
[Hugging Face](https://huggingface.co/carbon-lab) ·
[X-ray Atlas](https://xrayatlas.wsu.edu/) ·
[PRL](https://journals.aps.org/prl/abstract/10.1103/rfgg-ffyz) ·
[arXiv](https://arxiv.org/abs/2509.01734)

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
<summary><strong>Ingest</strong> · <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> 4/4</summary>

| | Capability | Completion | Location | Track |
|:-:|---|---|---|---|
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | `XrayT*.out` / TP `*.xas` sticks | tested | `dftlearn.io` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | XYZ geometry & site labels | tested | `dftlearn.io` |  |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Delta-KS / FINAL ENERGY / TP LUMO (`E^c`) | tested | `dftlearn.io` |  |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Full ground / excited / TP directory load | tested | `dftlearn.io.stobe_run` (`load_stobe_run`; staging `stobeLoader` remains) | [#2](https://github.com/WSU-Carbon-Lab/dft-learn/issues/2) |

</details>

<details>
<summary><strong>Filter & cluster</strong> · <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"> 5/8</summary>

| | Capability | Completion | Location | Track |
|:-:|---|---|---|---|
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Energy window + OS% cull | tested | `dftlearn.clustering` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Gaussian peak-overlap matrices | tested | `dftlearn.clustering` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Iterative overlap merge (`simpleCluster3`) | tested | `dftlearn.clustering` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | OS% elbow cutoff | tested | `dftlearn.clustering` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | OVP threshold selection (BIC / GP) | tested | `dftlearn.clustering` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"> | `filterDFT`-style end-to-end orchestration | in progress | `dftrun postprocess` (partial) | [#1](https://github.com/WSU-Carbon-Lab/dft-learn/issues/1) |
| <img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"> | OS × OVP sequential grids (`seqThresholds`) | not started | — |  |
| <img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"> | Amplitude refit to pre-merge DFT NEXAFS | not started | — | [#1](https://github.com/WSU-Carbon-Lab/dft-learn/issues/1) |
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
<summary><strong>Tensors & experiment</strong> · <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"> 3/5</summary>

| | Capability | Completion | Location | Track |
|:-:|---|---|---|---|
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | TP XAS reconstruction & broadening schedule | tested | `dftlearn.xas` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Molecular dipole / Cartesian OS tables | tested | `dftlearn.xas` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Bare-atom / Henke step edge | tested | `dftlearn.xas.step_edge` | [#3](https://github.com/WSU-Carbon-Lab/dft-learn/issues/3) |
| <img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"> | Film tilt / `simDFT` / angle-resolved model | in progress | `python_pipeline/multiSpecFitProcs` | [#3](https://github.com/WSU-Carbon-Lab/dft-learn/issues/3) |
| <img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"> | Multi-spectrum experiment fit | in progress | `python_pipeline/multiSpecFitProcs` | [#3](https://github.com/WSU-Carbon-Lab/dft-learn/issues/3) |

</details>

<details>
<summary><strong>Visualization & UI</strong> · <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"><img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"> 2/3</summary>

| | Capability | Completion | Location | Track |
|:-:|---|---|---|---|
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Cluster & reconstruction report figures | tested | `dftlearn.visualization` | [#4](https://github.com/WSU-Carbon-Lab/dft-learn/pull/4) |
| <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"> | Site / SCF / orbital summary figures | tested | `dftlearn.visualization` (`package_stobe_run`) |  |
| <img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"> | Interactive clustering panel parity | not started | HF Space / demos | [#5](https://github.com/WSU-Carbon-Lab/dft-learn/pull/5) |
| <img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"> | Chem3D molecular viewer | out of scope | — |  |
| <img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"> | Igor panel chrome / wave selectors | out of scope | — |  |

</details>

Overall: <img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-tested.svg" alt="■" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"><img src="docs/readme/box-wip.svg" alt="☒" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-todo.svg" alt="☐" width="12" height="12"><img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"><img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"><img src="docs/readme/box-skip.svg" alt="■" width="12" height="12"> **15/23**

## Citation

```bibtex
@article{Murcia2026OpticalTensors,
  title   = {Quantitative and Bond-Traceable Resonant X-Ray Optical Tensors of Organic Molecules},
  author  = {Murcia, Victor and Alqahtani, Obaid and Heilman, Harlan and Collins, Brian A.},
  journal = {Phys. Rev. Lett.},
  volume  = {137},
  issue   = {15},
  pages   = {158001},
  year    = {2026},
  doi     = {10.1103/rfgg-ffyz},
  url     = {https://journals.aps.org/prl/abstract/10.1103/rfgg-ffyz}
}
```

[`CITATION.cff`](CITATION.cff) · [lab publications](https://labs.wsu.edu/carbon/publications/)

## License

MIT ([`LICENSE`](LICENSE)). Maintainer: [Harlan Heilman](mailto:harlan.heilman@wsu.edu).
