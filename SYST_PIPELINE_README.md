# BDT Systematic Uncertainty Pipeline

Automated pipeline for computing photon cross-section systematic uncertainties
in the sPHENIX PPG12 analysis. Replaces the old hand-written `syst_*.C` ROOT macros.

---

## Directory Layout

```
ppg12/
├── efficiencytool/
│   ├── make_bdt_variations.py     # Step 1 — generates config_bdt_*.yaml
│   ├── config_bdt_nom.yaml        # nominal base config (hand-maintained)
│   ├── config_bdt_<name>.yaml     # generated variation configs (23 total)
│   ├── oneforall.sh               # runs MergeSim + CalculatePhotonYield for one config
│   ├── oneforall.sub              # HTCondor submit file (queues all config_bdt_*.yaml)
│   └── results/
│       └── Photon_final_bdt_<name>.root   # output of the efficiency tool per variant
└── plotting/
    ├── calc_syst_bdt.py           # Step 3 — computes deviations, writes ROOT + plots
    ├── rootFiles/
    │   ├── syst_bdt_<type>.root   # per syst_type deviation histograms
    │   ├── syst_bdt_<group>.root  # quadrature-summed group uncertainties
    │   └── syst_bdt_total.root    # total systematic (all groups + lumi)
    └── figures/
        ├── syst_bdt_rel_<type>.pdf    # per-type two-pad spectrum + relative deviation
        └── syst_bdt_breakdown.pdf     # all groups + total on one canvas
```

---

## Overview: Three-Step Workflow

```
Step 1                  Step 2                       Step 3
make_bdt_variations.py  oneforall.sub (HTCondor)     calc_syst_bdt.py
      │                       │                            │
      ▼                       ▼                            ▼
config_bdt_*.yaml  →  Photon_final_bdt_*.root  →  syst_bdt_*.root + *.pdf
```

---

## Step 1 — Generate Variation Configs

`make_bdt_variations.py` reads the nominal config and produces one `config_bdt_<name>.yaml`
per variant, with a single parameter overridden from the nominal.

```bash
cd efficiencytool
python make_bdt_variations.py config_bdt_nom.yaml
```

This writes 23 configs (including `config_bdt_nom.yaml` itself as a sanity check).
The nominal base config `config_bdt_nom.yaml` must exist and be hand-maintained before running.

### Variant list

| Config name | Syst type | Role | What changes |
|---|---|---|---|
| nom | — | — | identical to base (sanity check) |
| 0rad | — | — | 0-rad run range only (cross-check) |
| all | — | — | full run range (cross-check) |
| split | — | — | split cluster node (cross-check) |
| tightbdt50 | tight_bdt | down | BDT min threshold lowered (looser) |
| tightbdt70 | tight_bdt | up | BDT min threshold raised (tighter) |
| ntbdtmin02 | nt_bdt | one_sided | non-tight BDT lower boundary = 0.02 |
| noniso04 | noniso | down | noniso sideband window tightened |
| noniso10 | noniso | up | noniso sideband window widened |
| npb03 | npb_cut | down | NPB score cut = 0.3 (looser) |
| npb07 | npb_cut | up | NPB score cut = 0.7 (tighter) |
| purity_pade | purity_fit | one_sided | Padé purity fit function |
| vtxreweight0 | vtx_reweight | one_sided | vertex reweighting disabled |
| etbin_v3E_v3E | bdt_model | max | ET-binned model: v3E / v3E |
| etbin_E_E | bdt_model | max | ET-binned model: E / E |
| etbin_v1E_v3E | bdt_model | max | ET-binned model: v1E / v3E |
| b2bjet | b2bjet | one_sided | back-to-back jet cut applied |
| timingcut_2 | timing | down | timing window ±2 ns (tighter) |
| timingcut_5 | timing | up | timing window ±5 ns (looser) |
| energyscale26 | escale | one_sided | cluster energy scale × 1.026 |
| energyresolution7 | eres | max | energy resolution smearing σ = 7% |
| energyresolution8 | eres | max | energy resolution smearing σ = 8% |
| energyresolution10 | eres | max | energy resolution smearing σ = 9% |

---

## Step 2 — Run the Efficiency Tool on Every Variant

Each `config_bdt_<name>.yaml` must be processed through the efficiency tool
(`MergeSim.C` + `CalculatePhotonYield.C`) to produce `results/Photon_final_bdt_<name>.root`.

### Option A: HTCondor (recommended — runs all 23 variants in parallel)

```bash
cd efficiencytool
mkdir -p logs
condor_submit oneforall.sub
```

The submit file (`oneforall.sub`) queues one job per `config_bdt_*.yaml`. Each job runs
`oneforall.sh <config>` which executes `MergeSim.C` then `CalculatePhotonYield.C`.

Monitor jobs:
```bash
condor_q
# or watch specific logs:
tail -f logs/config_bdt_nom.out
```

### Option B: Run a single variant locally

```bash
cd efficiencytool
bash oneforall.sh config_bdt_nom.yaml
```

### Expected output

After all jobs complete, `results/` should contain one ROOT file per variant:
```
results/Photon_final_bdt_nom.root
results/Photon_final_bdt_tightbdt50.root
results/Photon_final_bdt_tightbdt70.root
... (23 files total)
```

---

## Step 3 — Compute Systematics

`calc_syst_bdt.py` reads the nominal and variation ROOT files, computes per-type
deviations, aggregates into groups by quadrature sum, and writes output ROOT files
and plots.

```bash
cd plotting
python calc_syst_bdt.py \
    --results /sphenix/user/shuhangli/ppg12/efficiencytool/results \
    --outdir  rootFiles \
    --figdir  figures \
    --skip-missing
```

`--skip-missing` silently skips any variation whose result file is absent
(useful during partial runs). Remove it once all 23 result files exist.

### All CLI options

| Option | Default | Description |
|---|---|---|
| `--results DIR` | `/sphenix/user/shuhangli/ppg12/efficiencytool/results` | Directory with `Photon_final_*.root` files |
| `--outdir DIR` | `rootFiles/` | Output directory for ROOT files |
| `--figdir DIR` | `figures/` | Output directory for PDF plots |
| `--nom VAR_TYPE` | `bdt_nom` | var_type string for the nominal file |
| `--histogram NAME` | `h_unfold_sub_result` | Histogram name inside each ROOT file |
| `--skip-missing` | off | Silently skip syst types with missing files |

### What the script does

1. **Per-type deviations** — for each `syst_type` in `SYST_TYPES`, loads the variation
   file(s), computes `(h_var - h_nom)` and `(h_var - h_nom)/h_nom` bin-by-bin, and
   aggregates according to the mode:
   - `two_sided`: separate up/down bands; mirrors if only one direction available
   - `one_sided`: symmetric band = `|delta|`
   - `max`: takes the bin-wise maximum over all "max" variants
   - `placeholder`: skipped with a warning (mbd, nor — not yet implemented)

2. **Group sums** — quadrature-sums the per-type results within each group:
   - `purity` = noniso ⊕ nt_bdt ⊕ purity_fit
   - `eff` = tight_bdt ⊕ npb_cut ⊕ vtx_reweight ⊕ bdt_model ⊕ b2bjet ⊕ timing
   - `escale` = escale
   - `eres` = eres
   - `mbd`, `nor` = placeholder (zero until variants are added)

3. **Total** — quadrature-sums all groups and adds the luminosity uncertainty
   (−6.7% / +9.1%, asymmetric) in quadrature.

4. **Plots** — saves per-type two-pad PDFs and a breakdown summary PDF.

### Output ROOT files

Each ROOT file contains four histograms with the same binning as the nominal spectrum:

| Histogram name | Content |
|---|---|
| `h_dev_low` | Absolute downward deviation (pb/GeV) |
| `h_dev_high` | Absolute upward deviation (pb/GeV) |
| `h_dev_rel_low` | Relative downward deviation (fractional) |
| `h_dev_rel_high` | Relative upward deviation (fractional) |

Files produced:

```
rootFiles/
  syst_bdt_tight_bdt.root       # two-sided efficiency syst (BDT cut)
  syst_bdt_nt_bdt.root          # one-sided purity syst (non-tight BDT boundary)
  syst_bdt_noniso.root          # two-sided purity syst (sideband window)
  syst_bdt_npb_cut.root         # two-sided efficiency syst (NPB score cut)
  syst_bdt_purity_fit.root      # one-sided purity syst (fit function)
  syst_bdt_vtx_reweight.root    # one-sided efficiency syst (vertex reweighting)
  syst_bdt_bdt_model.root       # max efficiency syst (BDT model choice)
  syst_bdt_b2bjet.root          # one-sided efficiency syst (back-to-back jet)
  syst_bdt_timing.root          # two-sided efficiency syst (timing cut)
  syst_bdt_escale.root          # one-sided energy scale syst
  syst_bdt_eres.root            # max energy resolution syst
  syst_bdt_purity.root          # purity group quadrature sum
  syst_bdt_eff.root             # efficiency group quadrature sum
  syst_bdt_escale.root          # (also group file, single member)
  syst_bdt_eres.root            # (also group file, single member)
  syst_bdt_total.root           # total syst = all groups + lumi in quadrature
```

---

## Adding a New Systematic Variant

1. Add a new entry to `VARIANTS` in `make_bdt_variations.py`:
   ```python
   dict(name="myvariant", some_override_key=value,
        syst_type="existing_or_new_type", syst_role="up"),
   ```

2. If `syst_type` is new, add it to `SYST_TYPES` and `SYST_GROUPS` in the same file.

3. Re-run Step 1 to regenerate configs, then re-run Step 2 for the new variant only:
   ```bash
   cd efficiencytool
   python make_bdt_variations.py config_bdt_nom.yaml
   bash oneforall.sh config_bdt_myvariant.yaml
   ```

4. Re-run Step 3:
   ```bash
   cd plotting
   python calc_syst_bdt.py --skip-missing
   ```

### Adding mbd / nor (currently placeholder)

The `mbd` and `nor` systematic types are defined in `SYST_TYPES` with
`"mode": "placeholder"` and are skipped by the pipeline with a warning.
To activate them:

1. Determine the parameter values for MBD efficiency up/down variants and the
   normalization variation.
2. Add the corresponding `dict(name=..., syst_type="mbd", syst_role=...)` entries
   to `VARIANTS` in `make_bdt_variations.py`.
3. Change `"mode": "placeholder"` to `"mode": "two_sided"` (or `"one_sided"`) in
   `SYST_TYPES` for the relevant type.
4. Follow the steps above to regenerate, run, and recalculate.

---

## Verification Checklist

After running all three steps, verify:

- [ ] 23 `config_bdt_*.yaml` files exist in `efficiencytool/`
- [ ] 23 `Photon_final_bdt_*.root` files exist in `efficiencytool/results/`
- [ ] 11 per-type `syst_bdt_*.root` files exist in `plotting/rootFiles/`
- [ ] 4–6 group `syst_bdt_{group}.root` files exist (fewer if mbd/nor still placeholder)
- [ ] `syst_bdt_total.root` exists and `h_dev_rel_high` values are ~5–30% in the 10–26 GeV range
- [ ] `figures/syst_bdt_breakdown.pdf` shows all available groups and total band
- [ ] Per-type PDFs `figures/syst_bdt_rel_*.pdf` look sensible (no runaway bins)
