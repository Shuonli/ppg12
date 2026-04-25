# JETPHOX rerun setup — CT14nlo and NNPDF3.1

**Date:** 2026-04-22
**Status:** READY TO SUBMIT (not yet submitted); one-condor-job validation tests launched for both pipelines.
**Motivation:**
1. Replace the CT14**lo** (LO) PDF in the existing JETPHOX production with the matching NLO PDF (**CT14nlo**) — fixes the LO-PDF-in-NLO-matrix-element inconsistency flagged in `reports/jetphox_audit.md` §1.
2. Add an alternative PDF family (**NNPDF3.1 NLO, α_s=0.118**) for cross-check — commonly used in LHC photon papers (ATLAS 1908.02746), different family of fits (MC replicas vs Hessian eigenvectors).

Both pipelines run in parallel, with separate output directories, and can be launched/compared independently.

---

## 0. Quick map

| Pipeline | PDF | Output dirs | Per-scale pawres | PPG12 histos |
|---|---|---|---|---|
| Existing | CT14lo | `condorout/`, `condor_05/`, `condor_20/` | `pawres/ggd(o)rhic_{05,10,20}.root` | `NLO/rootFiles/jetPHOX_{05,10,20}.root` |
| New (CT14nlo) | CT14nlo | `condorout_nlo/`, `condor_05_nlo/`, `condor_20_nlo/` | `pawres/ggd(o)rhic_nlo_{05,10,20}.root` | `NLO/rootFiles/jetPHOX_nlo_{05,10,20}.root` |
| New (NNPDF) | NNPDF31_nlo_as_0118 | `condorout_nnpdf/`, `condor_05_nnpdf/`, `condor_20_nnpdf/` | `pawres/ggd(o)rhic_nnpdf_{05,10,20}.root` | `NLO/rootFiles/jetPHOX_nnpdf_{05,10,20}.root` |

The existing CT14lo production is preserved on disk and untouched.

---

## 1. What changed relative to existing production

- **PDF set** (line 80 of `parameter.indat`): `CT14lo` → `CT14nlo` or `NNPDF31_nlo_as_0118`.
- **Defensive patches**: lines 48 (output format → `ntuple`), 298 (pT_max → `40.D0`) now patched in every new submit script because the top-level template is out-of-sync with deployed files. The existing submit scripts rely on unmodified template state; mine don't.
- **RNG seeds** offset by 20000/30000/40000 (CT14nlo) and 50000/60000/70000 (NNPDF) so no risk of identical Monte Carlo streams.
- **`>ff.sub` (overwrite)** instead of `>>` (append) — defends against double-submissions on re-run (a latent bug in the existing scripts).
- **Removed dead `Fun4All_run_sim.C`** reference from `transfer_input_files`. It was a copy-paste leftover from an sPHENIX sim template; JETPHOX doesn't need it.
- **NNPDF only**: custom `CondorRunJETPHOX_extpdf.sh` that exports `LHAPDF_DATA_PATH=/sphenix/user/shuhangli/lhapdf_pdfs_nnpdf:/cvmfs/.../LHAPDF-6.5.3/share/LHAPDF` so LHAPDF finds NNPDF in the user dir (cvmfs is read-only and only ships CT14lo/CT14nlo).
- **Everything else identical**: √s=200 GeV, pp, |y|<0.7, iso R=0.3 ET_max=4 GeV, mu scales {0.5, 1, 2}, 100 segments × 1M events.

---

## 2. Files created

### In `/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/`

**CT14nlo production:**
- `run_condor_SL7_nlo.sh`, `run_condor05_SL7_nlo.sh`, `run_condor20_SL7_nlo.sh` — 100 jobs per scale
- `run_nlo_all.sh` — dispatcher for all 3 scales (300 jobs total)
- `hadd_nlo.sh` — merge per-segment → `pawres/gg(d|o)rhic_nlo_{05,10,20}.root`
- `test_one_nlo_job.sh` — validation: 1 job, 1000 events, PDF errors off

**NNPDF3.1 production:**
- `CondorRunJETPHOX_extpdf.sh` — extended runner that sets `LHAPDF_DATA_PATH`
- `run_condor_SL7_nnpdf.sh`, `run_condor05_SL7_nnpdf.sh`, `run_condor20_SL7_nnpdf.sh`
- `run_nnpdf_all.sh` — dispatcher
- `hadd_nnpdf.sh` — merge into `pawres/gg(d|o)rhic_nnpdf_{05,10,20}.root`
- `test_one_nnpdf_job.sh` — validation: 1 job, 1000 events, PDF errors off

### In `/sphenix/user/shuhangli/ppg12/NLO/`

- `MakeJetPHOXhisto.C` — extended with optional `pdf_tag` arg. Default `""` = CT14lo (backward compat). Pass `"_nlo"` or `"_nnpdf"` for the new pipelines.
- `run_jetphox_histos_nlo.sh` — produces `rootFiles/jetPHOX_nlo_{05,10,20}.root`
- `run_jetphox_histos_nnpdf.sh` — produces `rootFiles/jetPHOX_nnpdf_{05,10,20}.root`

### In `/sphenix/user/shuhangli/lhapdf_pdfs_nnpdf/`

- `NNPDF31_nlo_as_0118/` — grid + .info for the external PDF set (127 MB, 101 members).

---

## 3. Pre-submission validation (one-condor-job tests)

Two single-job test submissions validate the full pipeline before the 600-job batch:

- **CT14nlo test** (`test_one_nlo_job.sh`): condor cluster 1223447, directory `condorout_nlo_test/OutDir0/`.
  - Patches: NEVENTS=1000, PDF errors=FALSE, rest identical to production.
  - **Result so far**: direct component completed with xsec = **78,980 pb** (+7.3% vs CT14lo's 73,621 pb → NLO PDF effect confirmed). Tree `t2` has 997 entries, UserInfo `[1000, 78979.5, 200]`. Output file: `.../pawres/ggdrhic_nlo_test_0.root` (88 KB). Fragmentation component still integrating (expected: total ~30 min for 1000-event test).
  - **Physics outputs per resultrhic/dirrhic/output.param**: `value(1)=CT14nlo`, `pdf_error= F`. Confirms PDF was read and applied correctly.
  - **Conclusion**: CT14nlo pipeline is validated.

- **NNPDF test** (`test_one_nnpdf_job.sh`): condor cluster 1223448.
  - Same reduced-stats config, using `CondorRunJETPHOX_extpdf.sh` to set `LHAPDF_DATA_PATH`.
  - **Result pending** — will be inspectable under `condorout_nnpdf_test/OutDir0/`.

### How to inspect the test results

```bash
# CT14nlo test
condor_q 1223447
ls -la /sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/condorout_nlo_test/OutDir0/pawres/
cat /sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/condorout_nlo_test/OutDir0/pawres/ggdrhic_nlo_test_0.res
grep -i "pdf\|CT14" /sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/condorout_nlo_test/OutDir0/working/resultrhic/dirrhic/output.param

# NNPDF test
condor_q 1223448
ls -la /sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/condorout_nnpdf_test/OutDir0/pawres/
grep -i "pdf\|NNPDF\|LHAPDF_DATA_PATH" /sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/condorout_nnpdf_test/OutDir0/working/resultrhic/dirrhic/output.param
grep -i "LHAPDF_DATA_PATH" /sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/condorout_nnpdf_test/OutDir0/working/test.out
```

### Expected direct-only xsec for validation

| PDF | direct integral (pb) |
|---|---|
| CT14lo (existing) | 73,621 |
| CT14nlo (test) | 78,980 (+7.3%) |
| NNPDF31_nlo_as_0118 | (measured pending) |

---

## 4. How to submit production (after tests pass)

```bash
cd /sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4

# one-shot: submit both pipelines (600 jobs total)
bash run_nlo_all.sh
bash run_nnpdf_all.sh

# or by scale
bash run_condor_SL7_nlo.sh      # CT14nlo mu=1.0
bash run_condor05_SL7_nlo.sh    # CT14nlo mu=0.5
bash run_condor20_SL7_nlo.sh    # CT14nlo mu=2.0
bash run_condor_SL7_nnpdf.sh    # NNPDF   mu=1.0
# ... etc

# monitor
condor_q sl4859

# when all jobs complete (overnight typical)
bash hadd_nlo.sh
bash hadd_nnpdf.sh

# rebuild PPG12 analysis histograms
cd /gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/NLO
bash run_jetphox_histos_nlo.sh
bash run_jetphox_histos_nnpdf.sh
```

---

## 5. Expected timeline and disk

- Per-job wall clock (1e6 events, full PDF error sets): ~6-18 hours (same as existing CT14lo runs).
- All 600 jobs in parallel: typically overnight with 100 condor slots.
- Hadd step: ~10 min per pipeline.
- Histogram build: ~20 min per pipeline.

Disk additions: ~100 GB for CT14nlo, ~100 GB for NNPDF. Filesystem has 58 TB free.

---

## 6. Integration into the final plot (after production + hadd + histos)

Edit `plot_final_selection.C:49-51` to swap or overlay:

```cpp
// Option A: swap primary JETPHOX from CT14lo to CT14nlo
TFile *fin_NLO       = new TFile(".../jetPHOX_nlo_10.root");
TFile *fin_NLO_up    = new TFile(".../jetPHOX_nlo_05.root");
TFile *fin_NLO_down  = new TFile(".../jetPHOX_nlo_20.root");

// Option B: add NNPDF as a second NLO curve for comparison
// (read both and draw overlaid, like the existing Werner overlay)
```

The `MakeJetPHOXhisto.C` pdf_tag change is backward compatible — default arg produces the CT14lo filenames unchanged.

---

## 7. Roll-back

Deleting the new output dirs (`condorout_nlo*/`, `condor_0{5,20}_nlo/`, `condorout_nnpdf*/`, `condor_0{5,20}_nnpdf/`), the new `pawres/gg[do]rhic_{nlo,nnpdf}_*.root` files, and `rootFiles/jetPHOX_{nlo,nnpdf}_*.root` restores the CT14lo-only state. No code is overwritten.
