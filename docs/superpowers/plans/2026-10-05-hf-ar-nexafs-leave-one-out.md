# HF AR-NEXAFS leave-one-out Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** In-browser film-tensor AR-NEXAFS vs experiment on CuPc/ZnPc Spaces; leave-one-out on Refinement matrix click.

**Architecture:** Shared `nexafsModel.js` (Igor `simDFTfit2` uniaxial path) + Space `app.js` / `index.html` wiring. Deploy to both HF Spaces.

**Tech Stack:** Vanilla JS, Plotly, existing Space CSV/JSON assets.

---

### Task 1: `nexafsModel.js` forward model

**Files:**
- Create: Space root `nexafsModel.js`

**Steps:**
1. Implement `filmMassAbsorption({xx,zz,amp,i0,alphaDeg,thetaDeg,phiDeg})`
2. Implement `igorGauss(energy, pos, wid)` and `spectrumFromPeaks(energy, peaks, globals, {leaveOneOutCluster})`
3. Export on `window.NexafsModel`

### Task 2: CuPc Space UI wiring

**Files:**
- Modify: `index.html`, `app.js`, include `nexafsModel.js`

**Steps:**
1. Clusters: add `#cluster-ar-plot` + all-final (optional initial) overlay
2. Refinement: matrix click → peak detail + leave-one-out AR plot
3. Smoke-test locally with `python -m http.server`

### Task 3: ZnPc Space parity

**Files:** same pattern under ZnPc clone (θ list 20–90, ZnPc CSV names)

### Task 4: Push Spaces + repo design doc

**Steps:**
1. Commit design/plan in dft-learn
2. `git push` both HF Space repos
