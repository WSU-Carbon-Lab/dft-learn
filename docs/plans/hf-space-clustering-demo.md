# Plan: Hugging Face Space demo for DFT clustering

**Status:** planning only (no implementation in this PR)
**Branch:** `plan/hf-space-clustering-demo`
**Date:** 2026-10-03
**Space:** https://huggingface.co/spaces/carbon-lab/dft-learn
**Library repo:** https://github.com/WSU-Carbon-Lab/dft-learn

---

## 1. Goal

Ship a public Hugging Face Space that demos the Carbon Lab DFT clustering / optical-tensor workflow for CuPc (and related molecules), while hosting the Igor project artifacts people need to reproduce or inspect the legacy pipeline.

### Done when (finish line)

1. **Space UI** loads curated clustered-result demos in-browser (spectra, overlays vs experiment, cluster sticks / OS plots, key parameters).
2. **Igor procedures** (canonical `.ipf` tree) are downloadable from the Hub with a clear install note matching `igor/README.md`.
3. **Curated Igor experiment packages** (selected `.pxp` and exported tables) live on the Hub without blowing past free storage / LFS limits.
4. **No paid compute required** for visitors: browse precomputed demos; optional live clustering is out of scope for v1 (or CPU-only on a tiny fixed sample only if free `cpu-basic` is acceptable later).
5. **Docs cross-link** GitHub library, Space, arXiv paper, and asset inventory so newcomers know what is code vs data vs UI.

### Explicit non-goals (v1)

- Interactive full Igor parity (panel, Chem3D, sequential OS/OVP grids).
- Running StoBe DFT inside the Space.
- Uploading every historical `.pxp` (many are multi-GB and currently OneDrive-dataless).
- Replacing `streamlit_demo/` as the long-term product; Space is the public demo surface.

---

## 2. Current state (inventory)

### 2.1 Hugging Face Space `carbon-lab/dft-learn`

| Item | Value |
|------|--------|
| SDK | `static` |
| Hardware | none (static) |
| Runtime | RUNNING |
| Contents | placeholder `index.html`, `style.css`, stub `README.md` frontmatter |
| License (card) | gpl-2.0 |
| Short description | Machine learning - DFT integrations |

Space is a blank static template. No demo data, no Igor assets, no Gradio/Streamlit app yet.

Authenticated Hub user for this work: `hduvallh` (org `carbon-lab`).

### 2.2 GitHub `WSU-Carbon-Lab/dft-learn`

Already present and relevant:

| Path | Role |
|------|------|
| `igor/` | Canonical Igor procedures (panel + User Procedures); ~740 KB, 14 `.ipf` + README |
| `src/dftlearn/clustering/` | Python overlap-merge core (not full Igor film/tensor stack) |
| `src/dftlearn/python_pipeline/` | Staging port of StoBe load / clustering / overlap / step-edge |
| `src/dftlearn/visualization/` | e.g. XAS cluster summary figures |
| `streamlit_demo/` | Legacy exploratory Streamlit shell (`stobe_loader.py`); not a stable library API |
| `docs/stobe/example-run/` | ZnPc StoBe example inputs / packaged outputs |
| Paper cite (streamlit README) | arXiv:2509.01734 |

### 2.3 OneDrive source trees (Victor Murcia)

**A. Previous Group Members / Projects**

```
.../Victor Murcia/Projects/
  DFT Calculations/     # AFRL, Assorted IGOR Files, CuPc_DFT, P3HT, data, backups
  IGOR Files/           # Clustering Files, Old Versions, Blade Coater, Misc, XRR
```

Rough counts (paths visible locally):

- ~232 `.pxp` under `DFT Calculations`
- ~53 `.pxp` + ~162 `.ipf` under `IGOR Files` (includes dated copies of the clustering code)
- `IGOR Files/Old Versions of Clustering Code/` has dated snapshots (2021-11 through 2022-03+) of the same procedure set already mirrored in repo `igor/`

**B. Group Publishing / Victor CuPc Optical Model Letter**

Figure and model `.pxp` files, paper drafts, PRL response folders, exported CSVs:

- `cupc-cui-data.csv` / `cupc-cui-data-edited.csv`
- `cupc-si-data.csv` / `cupc-si-data-edited.csv`
- `verify_and_plot_cupc_data.py`
- Figure `.pxp` (RawDFT vs Cluster, Model vs Exp, chi-squared heatmaps, overlap matrices, GIWAXS)

### 2.4 Critical hosting constraint: OneDrive + size

Many large `.pxp` paths report tens of GB and `dataless` / `blocks=0` (macOS OneDrive placeholders, not downloaded). Example: `CuPc Clustering - Reworked Scaling V1.pxp` (~13 GB logical, zero local blocks).

Small publication figure experiments are often fully local (e.g. `Figure 4 - RawDFTvsClusterDFT.pxp` ~437 KB).

**Implication:** do not plan a blind bulk upload. Curate, hydrate from OneDrive only what we keep, and prefer exported CSV/Parquet/JSON + PNG over giant `.pxp` for the browser demo.

HF free Space / Dataset storage and Git LFS quotas make dumping all param-space `.pxp` unrealistic. Prefer:

- **Dataset repo** (or Space `data/` + LFS) for curated binaries
- **GitHub `igor/`** as source of truth for procedures; Space mirrors or submodules a pinned copy

---

## 3. Recommended architecture

### Split Hub surface (preferred)

| Hub repo | Purpose |
|----------|---------|
| `spaces/carbon-lab/dft-learn` | Demo UI (browse + plot precomputed results) |
| `datasets/carbon-lab/dft-learn-igor` (new) | Igor `.ipf` archive + curated `.pxp` + exported demo tables/figures |
| GitHub `WSU-Carbon-Lab/dft-learn` | Library, CLI, tests, plan docs; optional sync workflow to Hub |

### Space SDK choice for v1

**Recommendation: Gradio on free `cpu-basic` (or keep Static if UI stays pure JS).**

| Option | Pros | Cons |
|--------|------|------|
| **A. Gradio + cpu-basic** | Fast Python plots, load CSV/Parquet from Dataset, reuse `dftlearn` plotting helpers | Needs Space SDK change from `static`; cold start |
| **B. Static + Plotly/Vega in browser** | Zero hardware; free forever | Must pre-bake JSON; harder to reuse Python viz |
| **C. Docker Streamlit** | Closest to `streamlit_demo/` | Heavier image; more ops; still no GPU needed |

**Decision for finish line:** Option A (Gradio) unless we want zero Python on HF, then B. Do **not** attach paid GPU. Live full clustering of visitor uploads is deferred.

Optional later (cute, cheap): Gradio tab that runs `dftlearn.clustering` on a **bundled tiny stick table** only (CPU seconds), not on arbitrary user StoBe trees.

### Data flow (v1)

```
OneDrive curated .pxp / CSV  -->  export scripts (local, GitHub)
                              -->  datasets/carbon-lab/dft-learn-igor
                              -->  Space loads via huggingface_hub / datasets
Visitor  -->  Space (select molecule / figure / OS-OVP case)
         -->  plots + download links (ipf zip, selected pxp, CSV)
```

---

## 4. What to host (asset curation)

### Tier 0 -- always host (small)

- Canonical `.ipf` tree = current GitHub `igor/` (not every Old Versions snapshot)
- Manifest JSON listing every curated asset, source path, license/attribution, paper figure mapping
- CuPc experiment CSVs already exported in the Optical Model Letter folder
- README install instructions for Igor User Procedures / Igor Procedures

### Tier 1 -- curated demo payloads (browser)

Export from publication-critical Igor experiments (hydrate from OneDrive first):

| Demo | Suggested sources | Exports needed |
|------|-------------------|----------------|
| Raw DFT vs clustered DFT | Figure 4 / RawDFT vs ClusterDFT `.pxp` | energy, abs raw, abs clustered |
| Model vs Exp (CuPc-CuI, CuPc-Si) | Figure 6 / Model vs Exp `.pxp`, `cupc-*-data*.csv` | theta series, model, experiment |
| Parameter change / chi-sq heatmap | Figure 5 / 7 `.pxp` | heatmap grids as CSV + PNG |
| OS/OVP refined params story | PPTX notes + Final_Comparison\* `.pxp` | stick tables, OS/OVP metadata |
| ZnPc StoBe package (optional) | `docs/stobe/example-run/znpc-cif/packaged_output/` | already in GitHub; mirror or link |

Target formats: Parquet or CSV + sidecar `meta.json` (OS%, OVP%, FWHM, alpha, molecule, citation).

### Tier 2 -- downloadable Igor experiments (not for in-browser parse)

Selected `.pxp` under a few MB to low tens of MB after hydration. Skip GIWAXS/image-heavy and multi-GB param-space experiments unless someone explicitly funds storage and LFS.

### Tier 3 -- archive / deferred

- Full `GaussianClustering CuPc very Fine Grain ParamSpace V*` series
- Duplicate Old Versions of Clustering Code trees
- Blade coater / XRR unrelated to clustering demo
- 100+ intermediate `CuPc Clustering vN` experiments

Document these in the manifest as "source-only on OneDrive; not mirrored".

---

## 5. Work breakdown (build to finish line)

### Phase 0 -- Decisions and access (blocking)

- [ ] Confirm Hub layout: Space-only vs Space + Dataset
- [ ] Confirm SDK: Gradio vs Static
- [ ] Confirm license on rehosted Victor Murcia materials (repo GPL-2.0 vs paper data; Chem3D third-party attribution)
- [ ] Hydrate OneDrive files needed for Tier 1/2 (mark Files On-Demand available offline)
- [ ] Measure hydrated sizes; set hard cap (e.g. Dataset < 5 GB for v1)

### Phase 1 -- Asset inventory and manifest

- [ ] Script a local inventory of OneDrive trees (path, logical size, dataless?, sha256 when hydrated)
- [ ] Tag each path: `procedure` / `demo-export` / `pxp-download` / `skip`
- [ ] Map paper figures (PRL / arXiv) to assets
- [ ] Write `manifest.json` schema and first filled copy under Dataset (or `docs/plans/assets/`)

### Phase 2 -- Export pipeline (GitHub, offline)

- [ ] Define export schema for clustered sticks and spectra (column names, units eV)
- [ ] Igor export checklist OR Python extractors when waves already CSV-exported
- [ ] Reuse / extend `verify_and_plot_cupc_data.py` patterns for CuPc-CuI / CuPc-Si
- [ ] Produce golden PNG previews for each demo card
- [ ] Unit tests on schema + smoke plots (`uv run pytest`)

### Phase 3 -- Hub data repo

- [ ] Create `datasets/carbon-lab/dft-learn-igor` (or agreed name)
- [ ] Upload Tier 0 procedures + Tier 1 tables/figures
- [ ] Upload Tier 2 selected `.pxp` via LFS as needed
- [ ] Dataset card: citation, install, file map, what is not included

### Phase 4 -- Space demo app

- [ ] Replace static placeholder README frontmatter (title, emoji, sdk, license, short_description)
- [ ] Implement Gradio (or Static) app:
  - Landing: one-line pitch + paper link + GitHub link
  - Demo picker: CuPc raw vs cluster, model vs exp (substrates), heatmap
  - Plot panel: interactive overlay
  - Downloads: procedures zip, CSV, selected `.pxp`
  - About: pipeline diagram (StoBe -> filter -> overlap cluster -> refit -> tensor)
- [ ] Load data from Dataset (not giant files in Space git if using Dataset)
- [ ] Pin dependency versions; prefer installing `dftlearn` from GitHub tag/commit or vendored plotting helpers
- [ ] Smoke test on HF after push

### Phase 5 -- Igor hosting polish

- [ ] Ensure Space/Dataset `igor/` matches GitHub `igor/` (CI sync or documented copy step)
- [ ] Publish install snippet (WaveMetrics User Procedures / Igor Procedures paths)
- [ ] Note Chem3D provenance (Knochenmuss) in Dataset card

### Phase 6 -- Optional cheap "cluster for people"

Only if Phase 4 is solid and free CPU is enough:

- [ ] Bundle one small preloaded stick table (ZnPc or CuPc subset)
- [ ] Gradio controls: OS cutoff, overlap threshold
- [ ] Call `dftlearn.clustering.cluster_by_overlap` / `select_overlap_threshold`
- [ ] Hard reject large uploads; no StoBe `.out` tree processing in v1

### Phase 7 -- Docs and release

- [ ] Link Space from GitHub README homepage
- [ ] Short `docs/` page: how the demo relates to Igor and `dftlearn.clustering`
- [ ] Zenodo / DOI note if packaging data release
- [ ] Checklist against arXiv:2509.01734 figures for demo coverage

---

## 6. Suggested Space IA (information architecture)

1. **Home** -- brand + one sentence + CTA "Browse CuPc demos"
2. **Demos** -- tabs or cards for Tier 1 plots
3. **Downloads** -- procedures, CSVs, curated `.pxp`
4. **Methods** -- high-level OS/OVP clustering explanation (no novel theory dump)
5. **Links** -- GitHub, paper, StoBe, lab

Keep v1 thin: Home + Demos + Downloads is enough.

---

## 7. Alternatives

| Alternative | When to choose |
|-------------|----------------|
| Keep Space **static** and put Python viz only in GitHub notebooks | Want zero HF runtime; demos are static HTML/JSON |
| Port `streamlit_demo` to HF Docker | Fastest reuse of existing UI; higher ops cost |
| Host only on GitHub Pages + release assets | Avoid HF; lose Hub discovery / org Spaces |
| Full interactive clustering + user uploads on ZeroGPU/paid CPU | Explicitly rejected for v1 cost reasons |

---

## 8. Risks and open questions

1. **Which `.pxp` are canonical for the paper figures?** (Figure 4/5/6/7 vs later `3-25-2025` Model vs Exp V3)
2. **Can we redistribute Chem3D `.ipf` on HF under the Dataset license?**
3. **Do we need WSU / Carbon Lab branding assets** beyond current placeholder emoji?
4. **Dataset name:** `dft-learn-igor` vs `dft-clustering-cupc`?
5. **Should GitHub Actions push to the Space**, or is manual `hf upload` enough for v1?
6. **Hydration time / disk:** confirming which multi-GB files are even worth fetching

---

## 9. Proposed implementation order (after this PR)

1. Phase 0 decisions (short meeting / reply on this PR)
2. Phase 1 manifest + hydrate Tier 1 sources
3. Phase 2 exports + tests in GitHub
4. Phase 3 Dataset create + upload
5. Phase 4 Space Gradio demo
6. Phase 5/7 polish; Phase 6 only if desired

---

## 10. Out-of-repo paths (reference)

```
# Procedures + large clustering experiments
/Users/hduva/Library/CloudStorage/OneDrive-SharedLibraries-WashingtonStateUniversity(email.wsu.edu)/Carbon Lab Research Group - Documents/Previous Group Members/Victor Murcia/Projects/IGOR Files
/Users/hduva/Library/CloudStorage/OneDrive-SharedLibraries-WashingtonStateUniversity(email.wsu.edu)/Carbon Lab Research Group - Documents/Previous Group Members/Victor Murcia/Projects/DFT Calculations

# Publication figures, model vs exp, exported CSVs
/Users/hduva/Library/CloudStorage/OneDrive-SharedLibraries-WashingtonStateUniversity(email.wsu.edu)/Carbon Lab Research Group - Documents/Group Publishing/Victor CuPc Optical Model Letter
```

---

## 11. Success metrics

- Space cold load shows a real CuPc demo in under ~30 s on cpu-basic (or instant if static)
- Visitor can download Igor procedures and open them in Igor following README
- At least three Tier 1 demos match paper figures qualitatively
- No paid HF hardware attached
- Manifest documents every hosted file and every intentionally excluded giant `.pxp`
