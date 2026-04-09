---
description: Inspect a ROOT file's contents (histograms, trees, graphs) using uproot
---

Inspect a ROOT file and summarize its contents.

Arguments: $ARGUMENTS (path to a ROOT file)

## Steps

1. Open the file with `python3` + `uproot` and list all keys with their types:
   ```python
   import uproot
   f = uproot.open('$ARGUMENTS')
   for k in sorted(f.keys()):
       obj = f[k]
       print(f'{k}: {type(obj).__name__}')
   ```

2. For each TH1 or TH2 histogram, print:
   - Number of bins, axis range
   - Integral (sum of bin contents)
   - Non-zero bin count
   - Min and max bin content (excluding under/overflow)
   - Limit to 20 histograms. If more exist, prioritize those matching `h_unfold*`, `h_data*`, `h_purity*`, `gpurity*`.

3. For TTree objects, print:
   - Entry count
   - Branch names (first 30, with total count if more)

4. For TGraph objects, print the number of points.

5. If the file is a `Photon_final_*.root`, additionally extract the cross-section spectrum:
   ```python
   h = f['h_unfold_sub_result']
   vals = h.values(); errs = h.errors(); edges = h.axis().edges()
   print('Cross-section spectrum (d sigma / d pT):')
   for lo, hi, v, e in zip(edges[:-1], edges[1:], vals, errs):
       print(f'  pT [{lo:.0f}, {hi:.0f}] GeV: {v:.4e} +/- {e:.4e} pb/GeV')
   ```

This is a **read-only** inspection. Do not modify any files.
