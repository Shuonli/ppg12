---
description: Run the 3-step BDT systematic uncertainty pipeline (generate configs, submit condor, aggregate)
---

Run the BDT systematic uncertainty pipeline. See SYST_PIPELINE_README.md for details.

Arguments: $ARGUMENTS
- No arguments → run all 3 steps interactively (ask before condor submit)
- `--step 1` → generate configs only
- `--step 2` → submit condor jobs only
- `--step 3` → aggregate systematics only

## Step 1 — Generate variation configs

```
cd /sphenix/user/shuhangli/ppg12/efficiencytool
python make_bdt_variations.py config_bdt_nom.yaml
```
Report how many `config_bdt_*.yaml` files were created.

## Step 2 — Submit to HTCondor

Ask the user to confirm before submitting.
```
cd /sphenix/user/shuhangli/ppg12/efficiencytool
condor_submit oneforall.sub
```
Report the cluster ID and number of jobs queued.

## Step 3 — Aggregate systematics

First check that all expected `Photon_final_bdt_*.root` files exist in `results/`. List any missing results and warn the user.
```
cd /sphenix/user/shuhangli/ppg12/plotting
python calc_syst_bdt.py --results ../efficiencytool/results --outdir rootFiles --figdir figures
```
Report output files created in `rootFiles/` and `figures/`.
