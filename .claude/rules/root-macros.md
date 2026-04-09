---
globs: "**/*.C"
---

# ROOT Macro Conventions (PPG12)

- Tree name: `slimtree` — load from config `configYaml["input"]["tree"]`, never hardcode
- Cluster node: from config `configYaml["input"]["cluster_node_name"]` (currently `CLUSTERINFO_CEMC`)
- Config loading: `gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so"); YAML::Node configYaml = YAML::LoadFile(configname);`
- Cross-section weights: `{sample}cross` variables (photon5cross, jet50cross, etc.) must be identical across MergeSim.C, RecoEffCalculator.C, and CalculatePhotonYield.C
- Output filenames must include `var_type` from config — prevents silent overwrites between systematic variations
- Parametric cuts are ET-dependent: `threshold = intercept + slope * ET` — never hardcode flat thresholds for BDT or isolation
- ABCD regions: tight_iso (A=signal), tight_noniso (B), nontight_iso (C), nontight_noniso (D). When editing cuts in the tight block, always check the non_tight and common blocks for consistency
- `SaveYamlToRoot()` pattern: preserve config archival in output ROOT files for reproducibility
- Plotting macros: `#include "plotcommon.h"` provides shared pT bins (`NptBins=12`, range 8-35 GeV), frames, legend strings. Use `SetsPhenixStyle()` and `BlairUtils.C` for error bands
- pT bins in `plotcommon.h` (`ptRanges[13]`) must stay in sync with analysis config pT bins
