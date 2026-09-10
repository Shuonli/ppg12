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
│   ├── config_bdt_<name>.yaml     # generated all-range configs (63, one per variant, 175 files in total with the feeders below, see Step 1)
│   ├── config_bdt_<name>_{0rad,1p5mrad}.yaml  # per-period merge-feeders (auto-generated)
│   ├── oneforall_tree_double_dispatch.sh  # tree stage for one config (RecoEffCalculator_TTreeReader, SI+DI mix)
│   ├── oneforall.sh               # yield stage for one config (merge_periods/MergeSim + CalculatePhotonYield)
│   ├── oneforall_tree_double.sub  # HTCondor submit for the tree stage (Phase 1)
│   ├── oneforall_all.sub          # HTCondor submit for the all-range yield stage (Phase 2)
│   └── results/
│       └── Photon_final_bdt_<name>.root   # output of the efficiency tool per variant
└── plotting/
    ├── calc_syst_bdt.py           # Step 3 — computes deviations, writes ROOT + plots
    ├── rootFiles/
    │   ├── syst_bdt_<type>.root   # per syst_type deviation histograms
    │   ├── syst_bdt_<group>.root  # quadrature-summed group uncertainties
    │   ├── syst_bdt_total.root    # total systematic (all groups + flat lumi and MBD-vertex terms)
    │   └── syst_sum.root          # same total as h_sum_* histograms, read by plot_final_selection.C
    └── figures/
        ├── syst_bdt_rel_<type>.pdf        # per-type two-pad spectrum + relative deviation
        ├── syst_bdt_rel_group_<group>.pdf # same per group
        └── syst_bdt_breakdown.pdf         # all groups + total on one canvas
```

---

## Overview: Three-Step Workflow

```
Step 1                  Step 2                       Step 3
make_bdt_variations.py  condor Phase 1 + Phase 2     calc_syst_bdt.py
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

This writes 175 configs from the 63 config-producing entries in `VARIANTS` (64 entries, of
which `iso_generator` is `aggregate_only` and writes no config). Every entry writes its
bare-name config `config_bdt_<name>.yaml` (the all-range analysis target, identical to the
base for `nom`). The 56 entries that do not pin their own `run_min`/`run_max` are also
expanded into two per-period merge-feeders, `config_bdt_<name>_0rad.yaml` and
`config_bdt_<name>_1p5mrad.yaml`, via `PER_PERIOD_OVERRIDES` (runs 47289–51273 with
`lumi` 47.2076 and runs 51274–54000 with `lumi` 17.1642, both with `lumi_target` 64.3718 so a
plain hadd of the two feeders reproduces the all-range MC). The 7 period-pinned entries
(`unfold_iter1/3/4`, `unfold_prior_flat`, `di_frac_fit`, `di_frac_fit_0rad`,
`di_frac_fit_1p5mrad`) write one config each.
The nominal base config `config_bdt_nom.yaml` must exist and be hand-maintained before running.

### Variants in the systematic budget (21 as of 2026-09)

Only entries with a non-null `syst_type` enter `calc_syst_bdt.py`. Role is the `syst_role`
the aggregator buckets on; mode and group come from `SYST_TYPES`.

| Config name | Syst type | Role | Mode | Group | What changes |
|---|---|---|---|---|---|
| mciso_no_shift | iso_resolution | one_sided | one_sided | efficiency | `mc_iso_shift = 0` (additive iso-ET pedestal off, scale stays nominal) |
| iso_generator | iso_generator | one_sided | one_sided | efficiency | HERWIG7 vs PYTHIA8 iso efficiency. No config: `Photon_final_bdt_iso_generator.root` is built post-hoc by `plotting/build_iso_generator_variant.C` |
| noniso04 | noniso | down | two_sided | purity | `reco_noniso_min_shift = 0.1` (non-iso sideband starts 0.1 GeV above the iso cut, nominal 0.8) |
| noniso10 | noniso | up | two_sided | purity | `reco_noniso_min_shift = 1.0` (non-iso sideband starts 1.0 GeV above the iso cut, nominal 0.8) |
| npb03 | npb_cut | down | two_sided | npb | `common.npb_score_cut = 0.3` (looser) |
| npb07 | npb_cut | up | two_sided | npb | `common.npb_score_cut = 0.7` (tighter) |
| purity_pade | purity_fit | one_sided | one_sided | purity | `fit_option = 0` (Erf purity fit `[0]*Erf((x-[1])/[2])` in place of the nominal `fit_option = 1` Padé form, so the variant name is a legacy misnomer) |
| purity_fit_ci_up | purity_fit_ci | up | two_sided | purity | `fittingerror = 1` (upper edge of the 68.3% fit confidence band) |
| purity_fit_ci_down | purity_fit_ci | down | two_sided | purity | `fittingerror = -1` (lower edge of the 68.3% fit confidence band) |
| mc_purity_correction | mc_purity_correction | one_sided | one_sided | purity | `mc_purity_correction = 1`, `mc_purity_corr_fitorder = 1` (pol1-smoothed MC non-closure ratio) |
| no_unfolding_reweighting | reweight | one_sided | one_sided | unfolding | `unfold.reweight = 0` (unweighted response prior) |
| tightup_p05 | photon_id_tight | one_sided | one_sided | purity | `tight.bdt_min_intercept = 0.865625` (+0.05 on the nominal line, slope unchanged) |
| ntdown_m10 | photon_id_nontight | one_sided | one_sided | purity | `non_tight.bdt_min_intercept = 0.6333` (−0.10 on the nominal line, slope unchanged) |
| energyscale148up | escale | down | two_sided | escale | `cluster_escale = 1.0148` on MC cluster ET (moves the cross-section down, so the role follows the cross-section direction) |
| energyscale148down | escale | up | two_sided | escale | `cluster_escale = 0.9852` on MC cluster ET |
| escale_nl | escale_nl | one_sided | one_sided | escale | `cluster_escale_nl_slope = 5e-4` per GeV, `cluster_escale_nl_ref_ET = 10` (0% at 10 GeV, +1% at 30 GeV) |
| eres_smear_none | eres | max | max | eres | `cluster_eres_data_p0/p1/p2 = 0` (no extra response smearing) |
| eres_smear_cE0p08 | eres | max | max | eres | `cluster_eres_data_p0/p1/p2 = 0.13/0.08/0.08` (wider data-resolution fit) |
| unfold_iter3 | unfold_iter | max | max | unfolding | `unfold.resultit = 3` (nominal 2). Result file written by `postprocess_unfold_iter_scan.py` from the nominal output, no re-run |
| unfold_iter4 | unfold_iter | max | max | unfolding | `unfold.resultit = 4`, same post-processing |
| di_frac_fit | di_fraction | one_sided | one_sided | di_fraction | DI blending fraction at the fitted values, carried by the two period-pinned feeder entries `di_frac_fit_0rad` / `di_frac_fit_1p5mrad` (`double_frac_override` 0.290 at 0 mrad, 0.000 at 1.5 mrad) |

The other 43 entries carry `syst_type=None` and are cross-checks that never enter the
quadrature: `nom`, `mciso_noscaleshift`, `mciso_no_scale`, the six `flat_t*_nt*`
flat-threshold partitions, `iso_p49_s37`, `mc_purity_correction_pol2`, `vtxreweight0`,
`bdtmodel_v0`, `bdtmodel_v0E`, `etbin_E_E`, the four `b2bjet_*`, the ±1.1/1.5/2.6%
`energyscale*` pairs, `unfold_iter1`, `unfold_prior_flat`, `di_frac_fit_0rad`,
`di_frac_fit_1p5mrad`, `pjlike`, the nine `ntbdtpair_*` anchor scans and the four
`mask_phisymm_*` tower masks.

---

## Step 2 — Run the Efficiency Tool on Every Variant

Each variant must be processed through the efficiency tool to produce
`results/Photon_final_bdt_<name>.root`: first the tree stage per period
(`RecoEffCalculator_TTreeReader.C` on the 8 SI/DI sample pairs + data, then `MergeSim.C`),
then the all-range merge + `CalculatePhotonYield.C`. `MergeSim`/`oneforall.sub` alone is
not sufficient whenever cuts or per-sample MC content changed.

### Two condor phases (production)

```bash
cd efficiencytool
mkdir -p logs

# Phase 1 — tree stage for every *_0rad.yaml / *_1p5mrad.yaml (and for standalone
# configs without feeders). oneforall_tree_double_dispatch.sh derives DOUBLE_FRAC
# from analysis.run_min (0.224 below run 51274, 0.079 otherwise, or
# analysis.double_frac_override when set), runs oneforall_tree_double.sh and
# MergeSim. Bare-name configs whose two feeders exist on disk are skipped here.
condor_submit oneforall_tree_double.sub   # queue: config_bdt*.yaml glob, 16 GB per job

# Phase 2 — after Phase 1 drains: all-range merge + yield for the bare-name configs
# (oneforall.sh: merge_periods.sh hadd of the two feeder MC outputs, hadd of the two
# per-period data_histo files, CalculatePhotonYield x 2, selection plots).
condor_submit oneforall_all.sub           # queue: all_range_configs.list, 4 GB per job
```

The schedd rejects submit files whose event `log` is not under `/tmp/`. Both files
above still carry `log = logs/$(cfg)...`, so change that line to
`log = /tmp/<name>_$(Cluster).log` before submitting (`output`/`error` may stay under
`logs/`). `oneforall_tree_double_rerun.sub` already follows this pattern and queues
from `syst_configs_rerun.list` (62 entries: nominal, the budget variants with their
feeders, and the `mc_purity_correction_pol2` cross-check triplet); use it when only the
systematic set needs a re-run.

Two budget entries need no condor job: `unfold_iter3/4` come from
`python3 postprocess_unfold_iter_scan.py` on the nominal output, and `iso_generator`
from `root -l -b -q ../plotting/build_iso_generator_variant.C` (the macro uses absolute
input paths). Run both after `Photon_final_bdt_nom.root` is final.

Monitor jobs:
```bash
condor_q
# or watch specific logs:
tail -f logs/config_bdt_nom_0rad.yaml.tree_double.out
```

### Run a single variant locally

```bash
cd efficiencytool
bash oneforall_tree_double_dispatch.sh config_bdt_noniso04_0rad.yaml
bash oneforall_tree_double_dispatch.sh config_bdt_noniso04_1p5mrad.yaml
bash oneforall.sh config_bdt_noniso04.yaml   # merge + yield, needs both feeders done
```

### Expected output

After both phases, `results/` holds one all-range file per variant plus the per-period feeders:
```
results/Photon_final_bdt_nom.root              # nominal (all-range)
results/Photon_final_bdt_nom_0rad.root         # 0 mrad feeder
results/Photon_final_bdt_nom_1p5mrad.root      # 1.5 mrad feeder
results/Photon_final_bdt_noniso04.root
... one Photon_final_bdt_<name>.root per budget variant (21 files including
    iso_generator and unfold_iter3/4), each with a *_mc.root closure companion
    except iso_generator, whose builder writes only the data-side file
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
(useful during partial runs). Remove it once every budget variant has its result file.

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

1. **Per-type deviations** — for each of the 15 `syst_type` keys in `SYST_TYPES`, loads the
   variation file(s) `Photon_final_bdt_<name>.root` (histogram `h_unfold_sub_result`),
   computes `(h_var - h_nom)` and `(h_var - h_nom)/h_nom` bin-by-bin, and
   aggregates according to the mode:
   - `two_sided`: low band = `|delta|` of the `down` variant, high band = `|delta|` of the
     `up` variant; mirrors if only one direction is available
   - `one_sided`: symmetric band = `|delta|`
   - `max`: low = high = bin-wise maximum of `|delta|` over all `max` variants
   - `placeholder`: skipped with a warning (no type uses it at present)

2. **Group sums** — quadrature-sums the per-type results within each `SYST_GROUPS` entry:
   - `escale` = escale ⊕ escale_nl
   - `eres` = eres
   - `purity` = photon_id_tight ⊕ photon_id_nontight ⊕ noniso ⊕ purity_fit ⊕ purity_fit_ci ⊕ mc_purity_correction
   - `efficiency` = iso_resolution ⊕ iso_generator
   - `unfolding` = reweight ⊕ unfold_iter
   - `di_fraction` = di_fraction
   - `npb` = npb_cut

3. **Total** — quadrature-sums the seven groups in `FINAL_SYSTS` (all of the above), then
   adds the two multiplicative-flat terms in `FLAT_SYSTS` in quadrature, post-unfold
   (`add_flat()`):
   - `lumi`: −6.75% / +9.13% (asymmetric, from the MBD inelastic cross section 25.2 +2.3/−1.7 mb)
   - `mbd_vertex_eff`: ±6% (symmetric MBD-vertex-efficiency envelope)

   Per bin and side: `total = sqrt( sum_groups dev^2 + (f_lumi * nom)^2 + (f_mbd * nom)^2 )`.

4. **Plots** — saves per-type two-pad PDFs (`syst_bdt_rel_<type>.pdf`), per-group PDFs
   (`syst_bdt_rel_group_<group>.pdf`) and the breakdown summary PDF.

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
  syst_bdt_reweight.root               # one_sided  (unfolding group)
  syst_bdt_photon_id_tight.root        # one_sided  (purity group)
  syst_bdt_photon_id_nontight.root     # one_sided  (purity group)
  syst_bdt_noniso.root                 # two_sided  (purity group)
  syst_bdt_purity_fit.root             # one_sided  (purity group)
  syst_bdt_purity_fit_ci.root          # two_sided  (purity group)
  syst_bdt_mc_purity_correction.root   # one_sided  (purity group)
  syst_bdt_npb_cut.root                # two_sided  (npb group)
  syst_bdt_iso_resolution.root         # one_sided  (efficiency group)
  syst_bdt_iso_generator.root          # one_sided  (efficiency group)
  syst_bdt_escale_nl.root              # one_sided  (escale group)
  syst_bdt_unfold_iter.root            # max        (unfolding group)
  syst_bdt_escale.root                 # per-type file, then OVERWRITTEN by the escale group
                                       #   sum (escale ⊕ escale_nl), so on disk it is the group
  syst_bdt_eres.root                   # per-type and group file (single member)
  syst_bdt_di_fraction.root            # per-type and group file (single member)
  syst_bdt_purity.root                 # group quadrature sums
  syst_bdt_efficiency.root
  syst_bdt_unfolding.root
  syst_bdt_npb.root
  syst_bdt_total.root                  # total = all groups + lumi + mbd_vertex_eff in quadrature
  syst_sum.root                        # same total as h_sum_low/high, h_sum_rel_low/high
                                       #   (read by plot_final_selection.C)
```

Stale files from retired types (`tight_bdt`, `nt_bdt`, `nor`, `timing`, ...) may still sit in
`rootFiles/` from older runs; the script never deletes them, so go by the list above.

---

## Adding a New Systematic Variant

1. Add a new entry to `VARIANTS` in `make_bdt_variations.py`:
   ```python
   dict(name="myvariant", some_override_key=value,
        syst_type="existing_or_new_type", syst_role="up"),
   ```

2. If `syst_type` is new, add it to `SYST_TYPES` and `SYST_GROUPS` in the same file.

3. Re-run Step 1 to regenerate configs, then re-run Step 2 for the new variant only
   (both feeders through the tree stage, then the bare name for merge + yield):
   ```bash
   cd efficiencytool
   python make_bdt_variations.py config_bdt_nom.yaml
   bash oneforall_tree_double_dispatch.sh config_bdt_myvariant_0rad.yaml
   bash oneforall_tree_double_dispatch.sh config_bdt_myvariant_1p5mrad.yaml
   bash oneforall.sh config_bdt_myvariant.yaml
   ```

4. Re-run Step 3:
   ```bash
   cd plotting
   python calc_syst_bdt.py --skip-missing
   ```

### Flat terms and post-hoc variants

The flat terms are not variants. `LUMI_SYST` and `MBD_VTX_SYST` live in
`make_bdt_variations.py` and are combined through `FLAT_SYSTS` by `add_flat()` in
`calc_syst_bdt.py`. To change one, edit the fraction there and re-run Step 3 only.

Two budget entries bypass Steps 1–2:

- `iso_generator` is declared with `aggregate_only=True`, so Step 1 writes no config.
  Its result file is built by `plotting/build_iso_generator_variant.C`, which scales the
  nominal spectrum by the PYTHIA/HERWIG iso-efficiency ratio bin by bin.
- `unfold_iter3` / `unfold_iter4` are period-pinned single configs whose result files are
  written by `efficiencytool/postprocess_unfold_iter_scan.py` from the nominal
  `h_unfold_sub_leak_{N}` histograms (no condor job).

---

## Verification Checklist

After running all three steps, verify:

- [ ] the 175 generated `config_bdt_*.yaml` files exist in `efficiencytool/` (hand-made extras may sit alongside)
- [ ] `Photon_final_bdt_<name>.root` exists in `efficiencytool/results/` for `nom` and all 21 budget variants
- [ ] 15 per-type and 7 group `syst_bdt_*.root` files were written to `plotting/rootFiles/` (19 distinct names, since `escale`, `eres` and `di_fraction` share a per-type and a group name)
- [ ] `syst_bdt_total.root` and `syst_sum.root` exist and `h_dev_rel_high` / `h_sum_rel_high` are roughly 15–35% in the 10–36 GeV bins
- [ ] `figures/syst_bdt_breakdown.pdf` shows all seven groups and the total band
- [ ] Per-type and per-group PDFs `figures/syst_bdt_rel_*.pdf` look sensible (no runaway bins)
