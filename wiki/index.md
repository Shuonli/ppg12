# PPG12 Codebase Wiki

Isolated photon cross-section in pp at sqrt(s) = 200 GeV, sPHENIX (Run 24).

## Pipeline

- [Overview](pipeline/overview.md) -- End-to-end analysis flow from DSTs to final cross-section
- [Tree Making](pipeline/01-tree-making.md) -- CaloAna24.cc: DST to slimtree conversion
- [BDT Training](pipeline/02-bdt-training.md) -- XGBoost photon-ID and NPB training pipeline
- [BDT Application](pipeline/03-bdt-application.md) -- apply_BDT.C: scoring clusters in slimtrees
- [Efficiency and Yield](pipeline/04-efficiency-yield.md) -- RecoEffCalculator, MergeSim, CalculatePhotonYield, ABCD, unfolding
- [Plotting and Systematics](pipeline/05-plotting-systematics.md) -- Final plots, systematic aggregation, report generators

## Concepts

- [ABCD Method](concepts/abcd-method.md) -- Background subtraction with signal leakage corrections
- [Isolation Cuts](concepts/isolation-cuts.md) -- Parametric isolation, cone sizes, tower vs topo, truth isolation
- [Shower Shape Variables](concepts/shower-shape-variables.md) -- All shower shape variables with definitions
- [Unfolding](concepts/unfolding.md) -- Bayesian unfolding, response matrix, closure tests
- [Systematic Variations](concepts/systematic-variations.md) -- Variation taxonomy, quadrature groups, adding new ones
- [Double-Interaction Efficiency](concepts/double-interaction-efficiency.md) -- Efficiency drop investigation: truth-matching artifact from vertex displacement
- [Eta Edge Migration](concepts/eta-edge-migration.md) -- Fiducial boundary fake rate from eta migration (1.2% signal contamination)
- [Delta-R Truth Matching](concepts/deltaR-truth-matching.md) -- Truth-matching deltaR study: vertex-driven artifact, cut removed
- [Efficiency Delta-R Comparison](concepts/efficiency-deltaR-comparison.md) -- Per-stage efficiency with deltaR removed: single, double, mixed overlay
- [MC Purity Correction](concepts/mc-purity-correction.md) -- ABCD independence violation, R-factor, correction factor investigation

## Reference

- [Config Schema](reference/config-schema.md) -- Every YAML field with type, meaning, and which code reads it
- [Constants Sync](reference/constants-sync.md) -- Constants that must match across files, with current discrepancies
- [Data Flow](reference/data-flow.md) -- File naming conventions, branch contracts between stages
- [MC Samples](reference/mc-samples.md) -- All MC samples, cross-sections, pT ranges, file locations
- [File Inventory](reference/file-inventory.md) -- Every significant file with one-line purpose

## Physics Knowledge Base

- [Physics Wiki Index](physics/index.md) -- Domain knowledge from ~70 research papers covering measurements, theory, techniques, and detector

## Guides

- [Running the Full Pipeline](guides/running-full-pipeline.md) -- Step-by-step commands
- [Adding a Systematic](guides/adding-a-systematic.md) -- Config, code, and aggregation steps
- [Common Issues](guides/common-issues.md) -- Known bugs, gotchas, and workarounds
