# Shower Shape Plotting Update Summary

## Overview
Updated `plot_showershapes_selections.C` to include ALL plots from the original `plot_showershapes.C` PLUS new BDT and isolation selection plots.

## Plot Structure

### Part 0: Standard Data/MC Comparison Plots (from original)
**Location**: Lines 101-296

Creates exactly the same plots as `plot_showershapes.C`:

1. **Projection Plots** (saved as `dis_*.pdf`):
   - Normalized 1D distributions of shower shape variables
   - Overlays data, signal MC, and background MC
   - For all shower shape variables: prob, CNN_prob, weta_cogx, wphi_cogx, et1-et4, e11_to_e33, e17_to_e77, e32_to_e35, bdt
   - For all 4 cuts: cut0 (no nbkg), cut1 (w/ nbkg), cut2 (tight), cut3 (non-tight)
   - For all eta and pT bins

2. **Background Profile Plots** (saved as `pfx_*.pdf`):
   - Mean isolation ET vs shower shape variable
   - Shows background MC only
   - Displays correlation factor
   - Same variables, cuts, and bins as projection plots

### Part 1: BDT Overlay Plots (NEW)
**Location**: Lines 298-403

Creates isolation ET distributions for different BDT score bins:
- Three BDT bins: 0.0-0.3, 0.3-0.7, 0.7-1.0
- Shows how isolation ET distribution changes with BDT score
- Saved as `iso_BDT_overlay_eta*_pt*.pdf`
- Gracefully handles missing histograms (if ShowerShapeCheck.C hasn't been rerun)

### Part 2: Isolation Region Overlay Plots (NEW)
**Location**: Lines 405-494

Compares shower shapes for isolated vs non-isolated regions:
- Tight (isolated, cut2) vs Non-tight (sideband, cut3)
- For all shower shape variables
- Shows background MC only
- Saved as `*_iso_overlay_eta*_pt*.pdf`

### Part 3: Data vs MC Isolation ET Comparison (NEW)
**Location**: Lines 496-577

Compares isolation ET distributions between data and inclusive MC:
- Data vs background MC (jet)
- For tight selection only
- Saved as `isoET_datamc_eta*_pt*.pdf`

## File Outputs

The script now creates:

1. **From Part 0 (Original plots):**
   - `dis_prob_eta*_pt*_cut*.pdf` - Projection plots for prob
   - `dis_CNN_prob_eta*_pt*_cut*.pdf` - Projection plots for CNN_prob
   - `dis_weta_cogx_eta*_pt*_cut*.pdf` - Projection plots for weta_cogx
   - `dis_wphi_cogx_eta*_pt*_cut*.pdf` - Projection plots for wphi_cogx
   - `dis_et1_eta*_pt*_cut*.pdf` - Projection plots for et1
   - `dis_et2_eta*_pt*_cut*.pdf` - Projection plots for et2
   - `dis_et3_eta*_pt*_cut*.pdf` - Projection plots for et3
   - `dis_et4_eta*_pt*_cut*.pdf` - Projection plots for et4
   - `dis_e11_to_e33_eta*_pt*_cut*.pdf` - Projection plots for e11_to_e33
   - `dis_e17_to_e77_eta*_pt*_cut*.pdf` - Projection plots for e17_to_e77
   - `dis_e32_to_e35_eta*_pt*_cut*.pdf` - Projection plots for e32_to_e35
   - `dis_bdt_eta*_pt*_cut*.pdf` - Projection plots for bdt
   - `pfx_*_eta*_pt*_cut*.pdf` - Profile plots (all variables above)

2. **From Part 1 (BDT overlays):**
   - `iso_BDT_overlay_eta*_pt*.pdf`

3. **From Part 2 (Isolation region overlays):**
   - `prob_iso_overlay_eta*_pt*.pdf`
   - `CNN_prob_iso_overlay_eta*_pt*.pdf`
   - `weta_cogx_iso_overlay_eta*_pt*.pdf`
   - `wphi_cogx_iso_overlay_eta*_pt*.pdf`
   - `et1_iso_overlay_eta*_pt*.pdf`
   - `e11_to_e33_iso_overlay_eta*_pt*.pdf`
   - `e32_to_e35_iso_overlay_eta*_pt*.pdf`
   - `bdt_iso_overlay_eta*_pt*.pdf`

4. **From Part 3 (Data/MC isolation):**
   - `isoET_datamc_eta*_pt*.pdf`

**Total plots**: ~500+ PDFs (depends on number of eta and pT bins)

## h2d_bdt Histogram Warning - RESOLVED

The warning "Could not retrieve histograms for h2d_bdt" occurred because:

1. **Root Cause**: File naming mismatch due to inclusive flag
   - The plotting script was loading non-inclusive files: `MC_efficiencyshower_shape_photon20.root` and `MC_efficiencyshower_shape_jet10.root`
   - These files were created on Mar 12, 2025 with older code that didn't create `h2d_bdt` histograms
   - The newer inclusive files (`*_inclusive.root`) created on Dec 29, 2024 DO have these histograms
   - The `hadd` commands in `run_showershape.sh` were combining non-inclusive files

2. **Solution**: Updated both scripts:

   **A. run_showershape.sh**: Updated `hadd` commands to combine `_inclusive` files:
   ```bash
   hadd -f MC_efficiencyshower_shape_signal_inclusive.root MC_efficiencyshower_shape_photon5_inclusive.root MC_efficiencyshower_shape_photon10_inclusive.root MC_efficiencyshower_shape_photon20_inclusive.root

   hadd -f MC_efficiencyshower_shape_jet_inclusive.root MC_efficiencyshower_shape_jet30_inclusive.root MC_efficiencyshower_shape_jet10_inclusive.root MC_efficiencyshower_shape_jet15_inclusive.root MC_efficiencyshower_shape_jet20_inclusive.root
   ```

   **B. plot_showershapes_selections.C**: Updated to use combined `_inclusive` files:
   - Signal MC: `MC_efficiencyshower_shape_signal_inclusive.root`
   - Background MC: `MC_efficiencyshower_shape_jet_inclusive.root`

3. **Note**: After running the updated `run_showershape.sh`, the combined files with `_inclusive` suffix will be created with all required histograms

## Key Features

- ✅ All original plots preserved exactly
- ✅ New BDT overlay plots added
- ✅ New isolation region comparison plots added
- ✅ New data/MC isolation ET comparison added
- ✅ Graceful error handling for missing histograms
- ✅ Consistent plotting style across all plots
- ✅ Clear console output showing progress

## Usage

```bash
cd /gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting
root -l -b -q 'plot_showershapes_selections.C()'
```

Output directory: `../PPG12-analysis-note/Figures/showershapes_selections/`

## Next Steps

To generate the BDT-binned histograms (Part 1 plots), run the updated ShowerShapeCheck.C:

```bash
cd /gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool
root -l -b -q 'ShowerShapeCheck.C("config_nom.yaml", "photon20")'
root -l -b -q 'ShowerShapeCheck.C("config_nom.yaml", "jet10")'
root -l -b -q 'ShowerShapeCheck.C("config_nom.yaml", "data")'
```

This will create the BDT-binned histograms needed for Part 1 plots.
