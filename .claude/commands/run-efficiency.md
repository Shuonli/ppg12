---
description: Run the efficiency pipeline (MergeSim + CalculatePhotonYield) for a single config
---

Run the efficiency pipeline for a single config in efficiencytool/.

Config: $ARGUMENTS (default: config_bdt_nom.yaml)

## Steps

1. Verify the config file exists at `/sphenix/user/shuhangli/ppg12/efficiencytool/$ARGUMENTS`. If no argument given, use `config_bdt_nom.yaml`.

2. Read the config and extract `var_type` from `output.var_type`.

3. Check if `results/Photon_final_{var_type}.root` already exists — warn the user it will be overwritten.

4. Run from `/sphenix/user/shuhangli/ppg12/efficiencytool/`:
   ```
   source /sphenix/u/shuhang98/setup.sh && bash oneforall.sh $ARGUMENTS
   ```

5. After completion, verify output files exist:
   - `results/Photon_final_{var_type}.root`
   - `results/Photon_final_{var_type}_mc.root`

6. Report file sizes. Flag any output < 1KB as likely failed.
