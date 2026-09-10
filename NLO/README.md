# NLO — JETPHOX theory reference for PPG12

The isolated-photon NLO curves drawn by `plotting/plot_final_selection.C` (lines 51-70) and the
`plotting/paper/plot_paper_final*.C` macros are histograms produced here from JETPHOX ntuples.
This file records how those histograms were made so that the chain can be followed from HEAD.
Everything stated here was read off the build directory and the scripts on 2026-09-10.

## Build

- JETPHOX build: `/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/`, unpacked from
  `/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4.tar.gz` (116315358 B, dated 2024-12-19).
  The version string comes from the tarball and directory name (1.3.1_4). The bundled
  `Readme_jetphox.html` is the generic "JETPHOX version 1.0" document dated 2016-10-20.
- Layout: `src/` (Fortran sources), `basesv5.1/` (BASES/SPRING integrator), `frag/` (fragmentation
  tables), `pdfa/` (PDF tables), `working/parameter.indat` (template card, 2025-06-18),
  `pawres/` (merged ntuples, a symlink to `/sphenix/tg/tg01/commissioning/CaloCalibWG/sli/jetphox/pawres`,
  248 files on 2026-09-10).
- External PDF grids: `/sphenix/user/shuhangli/lhapdf_pdfs_nnpdf/` (CT18NLO has 59 members there).
  The submit wrapper `CondorRunJETPHOX_extpdf.sh` puts that directory on `LHAPDF_DATA_PATH`.
- Each condor segment copies `basesv5.1 frag pdfa src working` into its own `OutDir<i>/`, edits
  `parameter.indat` with `sed` on fixed line numbers, and writes `pawres/gg{d,o}rhic_<name>.root`.
  Card lines touched: 38 (ntuple name), 48 (`ntuple` output), 80 (PDF set), 173/177/181 (scales),
  298 (ptmax), 302 (ptmin), 401 (seed), and 339 (isolation Etmax, iso scan only).

## Physics settings (template card `working/parameter.indat` unless noted)

- Process: direct (`ggdrhic`) plus one-fragmentation (`ggorhic`) prompt photon at NLO.
  The two channels are produced separately and summed in the histogram macros.
- Fragmentation function: BFG set II (card line 156 = `0200`).
- Photon rapidity window: -0.7 < y < 0.7 (lines 290 and 294).
- Isolation: flag 1 (fixed transverse energy), cone radius R = 0.3 (line 327),
  Etmax = 4.0 GeV (line 339). This matches the PPG12 truth isolation definition
  (R = 0.3, iso ET < 4 GeV). The `ct18iso200` production sets line 339 to 200 GeV, which
  removes the isolation requirement and is used for the isolated/inclusive correction.
- Scales: the initial-state factorisation, renormalisation and final-state factorisation
  scales (lines 173, 177, 181) are always set to the same multiple of the photon pT.
  Tag `10` = 1.0 pT (nominal), `05` = 0.5 pT, `20` = 2.0 pT.
- Generated events per segment: 1,000,000 (line 398).
- Chunked photon-pT production (all `*_chunked.root` files): five pT windows, each an independent
  set of condor segments with its own seed range, so that Vegas adapts inside a narrow range.
  A [8,12], B [12,16], C [16,22], D [22,30], E [30,40] GeV.
  CT18 nominal layout (`run_condor_SL7_ct18_chunked.sh`): 1/2/3/3/3 segments per chunk, seed
  bases 60000/61000/62000/63000/64000. The other PDF sets and scales were submitted with
  `run_chunked.sh <pdf_tag> <scale> <seed_base> [seg_csv]` (same five windows, default seg_csv
  1,2,3,3,3, the usage text shows 3,6,9,9,9 as an example) plus `run_chunked_more.sh` top-ups
  (4/6/6/6 extra segments for chunks B-E), so their per-chunk segment counts are not necessarily
  the CT18 layout and are not recorded anywhere.

## PDF sets

LHAPDF6 names as mapped in `run_chunked.sh` (line 80 of the card):

| tag | LHAPDF set | role |
|---|---|---|
| `ct18` | CT18NLO | nominal curve and scale band in `plot_final_selection.C` (`jetPHOX_ct18_{05,10,20}_chunked.root`) |
| `nlo` | CT14nlo | alternative central value (`jetPHOX_nlo_10_chunked.root`) |
| `nnpdf` | NNPDF31_nlo_as_0118 | alternative central value plus scale band |
| `nnpdf4` | NNPDF40_nlo_as_01180 | alternative central value plus scale band |
| `cteq` | cteq66 | alternative central value plus scale band |
| `msht` | MSHT20nlo_as118 | alternative central value plus scale band |
| `ct18nnlo`, `pdf4lhc21`, `nnpdfh`, `nnpdf4h` | CT18NNLO, PDF4LHC21_40, NNPDF31_nlo_as_0118_hessian, NNPDF40_nnlo_as_01180_hessian | PDF-member studies only (`rootFiles/jetPHOX_<tag>_10_pdfmem.root`, `compute_pdf_unc.py`) |

PDF uncertainty recipe (`compute_pdf_unc.py`, `PDF_INFO`): Hessian sets use the PDF4LHC15
asymmetric formula with C = 1.645 for CT/CTEQ (90% CL to 68% CL) and C = 1 for MSHT, the
symmetric-Hessian sets (nnpdfh, nnpdf4h, pdf4lhc21) use C = 1, and the NNPDF replica sets use the
replica mean and standard deviation. Output `rootFiles/jetPHOX_pdfunc.root`, read by
`plotting/paper/plot_paper_final_yj.C`, `_yj_v2.C` and `_yj_v3.C`.

## Scripts in the build directory (not in this repository)

| script | purpose |
|---|---|
| `run_condor_SL7_ct18_chunked.sh` | CT18NLO nominal-scale chunked submit (12 segments) |
| `run_chunked.sh`, `run_chunked_more.sh` | generic chunked submit for the other PDF tags and scales, plus extra segments |
| `run_condor_SL7_ct18iso_chunked.sh` | isolation-scan production (`ct18iso200`) |
| `hadd_ct18_chunked.sh`, `hadd_chunked.sh`, `hadd_ct18iso_chunked.sh` | per-chunk merges into `pawres/gg{d,o}rhic_<tag>_chunk{A..E}_<scale>.root` (never merged across chunks, each chunk keeps its own xsec and nb_evt in the TTree UserInfo) |
| `hadd_ct18.sh`, `hadd_nlo.sh`, `hadd_nnpdf.sh` | merges for the older non-chunked 100-segment productions into `pawres/gg{d,o}rhic_<tag>_{05,10,20}.root` |
| `run_condor{,05,20}_SL7_{ct18,nlo,nnpdf}.sh` | non-chunked 100-segment submits |
| `CT18_README.md` | staging notes for the non-chunked CT18 rerun, including the `plot_final_selection.C` switch from CT14nlo to CT18NLO |

About `hadd_ct18.sh`: `run_jetphox_histos_ct18.sh` (line 4) says to run it first. It is not part of
this repository. It exists only in the build directory (1545 B, 2026-04-30) and belongs to the
non-chunked CT18 production that feeds `MakeJetPHOXhisto.C`, not to the chunked files the final
plot reads. The chunked files come from `hadd_ct18_chunked.sh` (CT18) and `hadd_chunked.sh` (other tags).

## Scripts in this directory

- `MakeJetPHOXhisto.C` (non-chunked) → `rootFiles/jetPHOX<tag>_<scale>.root`, driven by `run_jetphox_histos{,_ct18,_nlo,_nnpdf}.sh`.
- `MakeJetPHOXhisto_chunked.C` → `rootFiles/jetPHOX<tag>_<scale>_chunked.root`. Per chunk file the
  normalisation is `xsec / nb_evt / nseg` from the TTree UserInfo, photons with |y| < 0.7 are filled
  with `pdf_weight[0]`, direct and fragmentation are summed, the result is divided by bin width and
  binned in `pT_bins_truth` of `efficiencytool/config_bdt_nom.yaml`. The `chunks` vector
  (label, nseg) is hard-coded to the CT18 layout 1/2/3/3/3 and must be edited to the number of
  segments actually merged for any other production.
- `MakeJetPHOXhisto_chunked_iso.C` — same for the iso scan (nseg 1/4/6/6/6, tags `_ct18iso200`, `_ct18iso4`).
- `MakeJetPHOXhisto_pdfmem.C`, `compute_pdf_unc.py` — per-member histograms and 68% CL PDF bands.
- `compute_iso_correction.py`, `make_iso_corr_hist.C` → `rootFiles/iso_correction_ct18nlo.{txt,root}` (isolated/inclusive ratio, read by `plotting/paper/plot_paper_final_yj.C` for the PHENIX comparison).
- `compute_eta_correction.py` → `rootFiles/truth_eta_ratio_jetphox.root` (read by `plotting/paper/plot_paper_final_yj.C`).
- `PlotTruthIso.C` — truth isolation plots.

## Histogram files kept in git (`rootFiles/`, force-added over the `*.root` ignore)

| file | date |
|---|---|
| `jetPHOX_ct18_{05,10,20}_chunked.root` | 2026-05-05 |
| `jetPHOX_nlo_{05,10,20}_chunked.root` | 2026-05-05 |
| `jetPHOX_nnpdf_{05,10,20}_chunked.root` | 2026-05-05 |
| `jetPHOX_nnpdf4_{05,10,20}_chunked.root` | 2026-05-10 |
| `jetPHOX_cteq_{05,10,20}_chunked.root` | 2026-05-05 / 05-07 |
| `jetPHOX_msht_{05,10,20}_chunked.root` | 2026-05-05 / 05-07 |
| `jetPHOX_ct18iso200_10_chunked.root` | 2026-06-02 |

Each holds `h_truth_pT` (dσ/dpT summed over direct and fragmentation, |y| < 0.7, pb/GeV per
the JETPHOX xsec units) and `h_truth_eta_pT`.

## State of the inputs on 2026-09-10

- `pawres/` holds the merged per-chunk inputs for `ct18` (05/10/20), `cteq` (05/10/20), `msht`
  (05/10/20), `nlo` (05/10/20), `nnpdf4` (05/10/20), `ct18iso200` (10), `ct18nnlo` (10),
  `nnpdfh` (10), `nnpdf4h` (10) and `pdf4lhc21` (10), ten files per tag and scale.
- There are no `*_nnpdf_chunk*` files left in `pawres/` (`cleanup_disk.sh` deletes
  `pawres/*_nnpdf_*.root`), so the three `jetPHOX_nnpdf_*_chunked.root` (NNPDF3.1) histograms
  cannot be regenerated from what is on disk.
- The per-segment condor output directories of the `ct18`, `cteq`, `msht`, `nlo`, `nnpdf` and
  `nnpdf4` chunked productions are gone (`cleanup_disk.sh` removes leftover `condorout_*`
  directories). Only `condorout_ct18iso200_chunked`, `condorout_ct18nnlo_chunked_10` and
  `condorout_pdf4lhc21_chunked_10` remain.
- No run logs of the chunked histogram macros are kept under `NLO/logs/`, so the per-chunk
  segment counts used for the non-CT18 chunked histograms are not recorded anywhere.
