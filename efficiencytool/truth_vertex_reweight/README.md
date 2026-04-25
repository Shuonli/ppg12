# Truth Vertex Reweight (PPG12)

Standalone fit of a truth-level vertex reweight `w(z)` for the PPG12
isolated-photon analysis. Produces one reweight function per crossing-angle
period that, when applied per-vertex to the MC (once for single-interaction
events, twice for double-interaction events), makes the mixed MC reco-vertex
distribution match the data reco-vertex distribution.

This directory is **standalone** — it does not modify the analysis pipeline.
The output is a ROOT file containing a `TH1` that downstream code can read
via `TH1::Interpolate` or equivalent. Integrating the reweight into
`ShowerShapeCheck.C` and `RecoEffCalculator_TTreeReader.C` is a separate,
later step gated on the closure quality produced here.

## Background

The MC was generated with the full RHIC beam profile
(`σ_gen ≈ 62 cm`), which is much wider than the actual 1.5 mrad luminous
region (`σ_tgt ≈ 22 cm`). The current production applies a 1D reco-vertex
reweight `f(z_r) = D(z_r)/M(z_r)` that closes the reco axis exactly but
leaves the underlying truth vertices at a width intermediate between MC and
target. For double-interaction MC, where two truth vertices share the same
generator profile and both feed the reconstructed vertex via MBD
charge-weighted averaging, this mismatch produces a ~+1.4 % cross-section
bias at 1.5 mrad. See `reports/truth_vertex_reco_check.tex` for the full
motivation.

The fix is a **factorized 2D Gaussian reweight** `w(z_h, z_mb) = w(z_h)·w(z_mb)`
with a single function `w(z)`, applied per truth vertex:

- Single-interaction events: `event_weight *= w(z_h)`
- Double-interaction events: `event_weight *= w(z_h) · w(z_mb)`

The function `w(z)` is determined by requiring the **mixed reco-vertex
distribution** (single + double MC at physical pile-up fractions,
plus the MC→reco smearing automatically provided by the slimtree's
`vertexz` branch) to match the **data reco-vertex distribution**.

## Two-stage fit

1. **Stage 1 — Parametric.** Fit a single-parameter Gaussian ratio
   `w(z; σ_tgt) ∝ exp(z²·(1/σ_gen² − 1/σ_tgt²)/2)` by scanning `σ_tgt`
   over a grid and picking the minimum χ² against data reco.

2. **Stage 2 — Residual correction.** After stage 1, measure the
   residual `c(v_r) = D(v_r)/M_stage1(v_r)` at reco level. Treat `c` as
   a truth-level function via the identity map `z ↔ v_r` (justified
   because MBD vertex resolution ≪ bin width). Apply `c(z)` as a
   per-vertex multiplier alongside `w(z)`. The residual is typically
   small after stage 1, so the over-correction incurred for double MC
   (where `v_r ≈ (z_h + z_mb)/2`) is second-order and numerically
   negligible at 1.5 mrad. At 0 mrad it is sub-percent; if the closure
   plot shows visible non-closure there, a half-power tweak is available
   (see `fit_truth_vertex_reweight.py` config flag).

## Factorization prerequisite (gate)

The whole "single function `w(z)`" idea rests on the assumption that in
the generator, `z_h` and `z_mb` are statistically independent and drawn
from the same distribution. This holds by construction (two independent
Pythia events thrown through the same beam profile), but must be verified
before running the fit. `check_factorization.py` reports three gates:

- `|corr(z_h, z_mb)| < 0.02`
- `|std(z_h) − std(z_mb)| < 2σ_stat`
- `|std(z_h in double) − std(z_h in nominal)| < 2σ_stat`

The fit script refuses to run unless these pass.

## Inputs

**Only the double-interaction MC needs to be reprocessed** through the
rebuilt anatreemaker. Single-interaction MC and data slimtrees can be
reused as-is:

| Sample | Reprocess? | Reason |
|---|---|---|
| `jet12_double` | **yes** | needs new `vertexz_truth_mb` branch |
| `photon10_double` | **yes** | needs new `vertexz_truth_mb` branch |
| `jet12` (nominal) | no | only uses `vertexz`, `vertexz_truth`, `runnumber` — all present in existing production |
| `photon10` (nominal) | no | same |
| data | no | data events never populate truth branches; `vertexz` and `runnumber` already present |

The loader handles "missing `vertexz_truth_mb` branch" gracefully (it is
only read if the branch is present in the tree). `build_T_s` never
consumes `z_mb` — single MC only contributes through the `(z_h, z_r)`
matrix. Double MC is the only consumer of the new branch via `build_T_d`.

Suggested minimum stats for the double reprocess: ~500 k events each,
~1 M is comfortable. The existing production has ~10 M events per
sample, so reprocessing a fraction is fine for a first pass.

Data is split by `runnumber` at fit time into 1.5 mrad and 0 mrad periods.

## Outputs

Per period (`output/1p5mrad/` and `output/0mrad/`):

- `fit_result.yaml` — σ_gen, best-fit σ_tgt, χ²/NDF, full iteration log,
  git hash of inputs, command-line arguments
- `reweight.root` — TH1 `h_c_truth` (the stage-2 residual correction on
  a truth-z axis), plus `TNamed` metadata: `sigma_gen`, `sigma_tgt`,
  `period`, `sample_list`, `git_hash`, `timestamp`. Downstream code
  consumes this one file.
- `closure.pdf` — 4-panel closure plots:
  1. Reco vertex: data vs mixed MC, before / after stage 1 / after stage 2
  2. Truth vertex: `z_h` and `z_mb` after full reweight, vs Gaussian(σ_tgt)
  3. `w(z)` parametric curve and `c(z)` residual TH1 overlaid
  4. Final reco residual `D/M_final` — should be flat at 1 within stats
- `fit.log` — textual log

## Python dependencies

- Python 3.8+
- `uproot` (tree I/O)
- `numpy`
- `scipy` (minimization, smoothing)
- `pyyaml`
- `matplotlib` (closure plots)

Already available in the sPHENIX python environment used elsewhere in
this repository.

## How to run

```bash
cd efficiencytool/truth_vertex_reweight

# Edit config.yaml to point at the reprocessed slimtree files.

# Stage 0 — factorization check (gate)
python3 check_factorization.py --config config.yaml

# Full fit (stage 1 + stage 2), per period
python3 fit_truth_vertex_reweight.py --config config.yaml --period 1p5mrad
python3 fit_truth_vertex_reweight.py --config config.yaml --period 0mrad

# Closure plots only (reads existing fit output)
python3 plot_closure.py --config config.yaml --period 1p5mrad
```

## Downstream consumption (future)

Once the closure is validated, a helper function in
`ShowerShapeCheck.C` and `RecoEffCalculator_TTreeReader.C` will look like:

```cpp
// Load once in Init()
TFile *f = TFile::Open("reweight.root");
TH1 *h_c = (TH1*) f->Get("h_c_truth");
float sigma_gen = ((TNamed*) f->Get("sigma_gen"))->GetTitle()->Atof();
float sigma_tgt = ((TNamed*) f->Get("sigma_tgt"))->GetTitle()->Atof();

// Per event
auto TruthVertexWeight = [&](float z) {
    float wG = std::exp(0.5f*z*z*(1.f/(sigma_gen*sigma_gen) - 1.f/(sigma_tgt*sigma_tgt)));
    float cZ = h_c ? h_c->Interpolate(z) : 1.f;
    return wG * cZ;
};
float w_truth = TruthVertexWeight(*vertexz_truth);
if (*vertexz_truth_mb > -9000.f) {
    w_truth *= TruthVertexWeight(*vertexz_truth_mb);
}
event_weight *= w_truth;
```

That hookup is **not** done here. This directory only produces the
reweight file.

## Interpretation of failures

- `check_factorization.py` gate fails → investigate why the two truth
  vertices are not statistically equivalent. Factorized `w(z)` cannot
  be used; the full 2D non-parametric approach would be needed.
- Stage 1 χ²/NDF ≫ 1 → Gaussian ansatz is too restrictive. Likely
  causes: non-Gaussian tails in data, wrong period split, broken
  cross-section weights. Look at `closure.pdf` panel 1.
- Stage 2 max residual > 2 % → identity-map approximation breaking
  down. Enable the half-power-for-double tweak in config, or escalate
  to the non-parametric iterative reweight (not implemented here).
- Parametric and non-parametric `w(z)` disagree visibly in panel 3 →
  the underlying truth distribution is not quite Gaussian. Use the
  non-parametric form as nominal and the parametric as a systematic.
