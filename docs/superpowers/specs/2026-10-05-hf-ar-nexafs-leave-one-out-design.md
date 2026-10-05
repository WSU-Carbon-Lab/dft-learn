# Browser AR-NEXAFS leave-one-out overlay (HF Spaces)

**Date:** 2026-10-05  
**Status:** approved for implementation  
**Surfaces:** `carbon-lab/cupc-optical-model`, `carbon-lab/znpc-optical-model`

## Goal

Let visitors see that refined cluster Gaussians quantitatively reproduce measured angle-resolved NEXAFS, and that leaving one peak at its **initial** parameters (others **final**) produces visible disagreement with experiment.

## UX

### Clusters tab
- New plot: model AR-NEXAFS (all peaks **final**, film α_f, i0, φ) overlaid on measured θ-series (same colors as Overview).
- Optional toggle: show all-**initial** model for context.
- Uses shared energy window / rigidShift where applicable.

### Refinement tab
1. Click heatmap cell (cluster × Δ-parameter) → focus that cluster.
2. Peak detail panel: selected peak **initial** vs **final** Gaussian + DFT member sticks.
3. Leave-one-out AR plot: all peaks final except selected peak initial → model vs experiment (+ residual note).
4. Clear selection → all-final match.

## Forward model (browser)

Port Igor `simDFTfit2` uniaxial film path:

1. Molecular tensor diagonal from cluster `xx`, `yy`, `zz` (abs).
2. Scale by `i0 * amp`.
3. Film average at sample tilt `α` (degrees → rad):
   - `Tin = (Txx*(1+c²) + Tzz*s²)/2`
   - `Tout = Txx*s² + Tzz*c²`
4. Mass absorption at measurement θ, φ: `MA = Tin*sin²θ + Tout*cos²θ` (φ=0 → e_y=0).
5. `I(E) += MA * exp(-((E-pos)/wid)²)` (Igor `gauss` convention).

Step-edge continuum is **out of scope for v1** (resonances-only overlay; caption notes this). IPs remain available for a later step.

## Data

- Peaks / globals: `param_change.json` + `cluster_gaussians.csv`
- Experiment: existing edited long-form CSVs (`energy,mu,theta,...`)
- Sticks: `*_sticks_slim.csv` for DFT members

## Non-goals

- Live least-squares refit
- Porting full `simDFT` step-edge / mask machinery
- Streamlit UI in this pass
