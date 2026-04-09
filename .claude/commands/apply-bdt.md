---
description: Apply trained BDT scores to simulation and data samples via apply_BDT.C
---

Apply BDT scores to simulation/data samples.

Samples: $ARGUMENTS
- If `all`: use photon5 photon10 photon20 jet10 jet15 jet20 jet30 jet50
- Otherwise: use the space-separated sample names provided (e.g., `photon5 jet20`)

Working directory: `/sphenix/user/shuhangli/ppg12/FunWithxgboost`
Config: `config_nom.yaml`

## Steps

1. Verify `config_nom.yaml` exists and read `input.photon_jet_file_root_dir`.

2. For each sample, verify the input ROOT file exists at `{root_dir}/{sample}/`.

3. For each sample, run:
   ```
   source /sphenix/u/shuhang98/setup.sh && cd /sphenix/user/shuhangli/ppg12/FunWithxgboost && root -l -b -q 'apply_BDT.C("config_nom.yaml", "{sample}")'
   ```

4. Report exit status for each sample (success/failure).

5. For any failures, show the last 20 lines of output to help diagnose.
