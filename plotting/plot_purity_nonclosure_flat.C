#include "plotcommon.h"

// --------------------------------------------------------------------------
// plot_purity_nonclosure_flat
//
// Visualises the ABCD purity non-closure on inclusive MC for 4 flat-threshold
// BDT partition variants, with the parametric nominal shown as a reference.
//
// Naming: flat_t{T*100}_nt{nt_min*100}
//   tight = [T, 1.0], non-tight = [nt_min, T]  (no gap at T)
//
// Variants (tight / non-tight):
//   nom               -- parametric: tight = 0.80 - 0.015*ET, non-tight = [0.02, tight]
//   flat_t85_nt50     -- tight = [0.85, 1.0], non-tight = [0.50, 0.85]
//   flat_t90_nt50     -- tight = [0.90, 1.0], non-tight = [0.50, 0.90]  (primary)
//   flat_t90_nt70     -- tight = [0.90, 1.0], non-tight = [0.70, 0.90]
//   flat_t95_nt50     -- tight = [0.95, 1.0], non-tight = [0.50, 0.95]
//
// Key difference vs the parametric nt_bdt scan: here the tight cut differs
// between variants, so region A differs and EACH variant has its OWN truth
// curve.  We therefore plot 4 truth curves (dashed, same color as variant),
// and 4 non-closure curves, each using its own truth.
//
// The parametric nominal is shown on the TOP panel as a light-gray reference
// (open circles) but excluded from the bottom panel to keep it readable.
//
// Inputs (all 12-bin, current CLUSTERINFO_CEMC pipeline):
//   Photon_final_bdt_nom_mc.root             -- reference (parametric)
//   Photon_final_bdt_flat_t85_nt50_mc.root
//   Photon_final_bdt_flat_t90_nt50_mc.root
//   Photon_final_bdt_flat_t90_nt70_mc.root
//   Photon_final_bdt_flat_t95_nt50_mc.root
//
// Each file contains
//   gpurity_leak   (TGraphErrors)      -- ABCD purity with signal-leak correction
//   g_purity_truth (TGraphAsymmErrors) -- truth purity (region A)
//
// Output:
//   figures/purity_nonclosure_flat.pdf    -- 2-panel canvas
//   rootFiles/purity_nonclosure_flat.root -- archived TGraphs
// --------------------------------------------------------------------------

namespace
{
  struct FlatEntry
  {
    std::string tag;       // file-tag, e.g. "flat_t90_nt50" or "nom"
    std::string suffix;    // short suffix used for ROOT object names
    std::string label;     // legend label
    int         color;
    int         marker;
    bool        is_ref;    // true for parametric nominal (reference only)
  };

  TGraphErrors *shiftGraphX(const TGraphErrors *g, double dx, const char *name)
  {
    const int n = g->GetN();
    TGraphErrors *gs = new TGraphErrors(n);
    for (int i = 0; i < n; ++i)
    {
      double x, y;
      g->GetPoint(i, x, y);
      gs->SetPoint(i, x + dx, y);
      gs->SetPointError(i, g->GetErrorX(i), g->GetErrorY(i));
    }
    gs->SetName(name);
    return gs;
  }

  // Non-closure = gABCD - gTruth, bin-by-bin on a shared 12-bin x-grid.
  // Errors added in quadrature; asymmetric truth errors symmetrised.
  TGraphErrors *buildNonClosure(const TGraphErrors *gABCD,
                                const TGraphAsymmErrors *gTruth,
                                const char *name)
  {
    const int n = gABCD->GetN();
    TGraphErrors *gdiff = new TGraphErrors(n);
    for (int i = 0; i < n; ++i)
    {
      double x, yA;
      gABCD->GetPoint(i, x, yA);
      int    jbest = -1;
      double dbest = 1e9;
      for (int j = 0; j < gTruth->GetN(); ++j)
      {
        double xj, yj;
        gTruth->GetPoint(j, xj, yj);
        double d = std::fabs(xj - x);
        if (d < dbest) { dbest = d; jbest = j; }
      }
      double xT, yT;
      gTruth->GetPoint(jbest, xT, yT);
      const double eA  = gABCD->GetErrorY(i);
      const double eTl = gTruth->GetErrorYlow(jbest);
      const double eTh = gTruth->GetErrorYhigh(jbest);
      const double eT  = 0.5 * (eTl + eTh);
      const double err = std::sqrt(eA * eA + eT * eT);
      gdiff->SetPoint(i, x, yA - yT);
      gdiff->SetPointError(i, gABCD->GetErrorX(i), err);
    }
    gdiff->SetName(name);
    return gdiff;
  }
}

void plot_purity_nonclosure_flat(const char *outpdf =
                                 "figures/purity_nonclosure_flat.pdf")
{
  init_plot();

  const std::string basepath =
      "/sphenix/user/shuhangli/ppg12/efficiencytool/results/";
  const std::string outroot =
      "/sphenix/user/shuhangli/ppg12/plotting/rootFiles/"
      "purity_nonclosure_flat.root";

  // Index 0 is the parametric nominal reference (top panel only).
  // Indices 1..4 are the flat variants (top + bottom).
  std::vector<FlatEntry> scan = {
    {"nom",           "nom",   "parametric nom (ref)",
                                 kGray + 2,    24, true },   // open circle
    {"flat_t85_nt50", "t85n50", "flat t=[0.85,1], nt=[0.50,0.85]",
                                 kRed + 1,     20, false},   // full circle
    {"flat_t90_nt50", "t90n50", "flat t=[0.90,1], nt=[0.50,0.90]",
                                 kBlue + 1,    21, false},   // full square
    {"flat_t90_nt70", "t90n70", "flat t=[0.90,1], nt=[0.70,0.90]",
                                 kGreen + 2,   22, false},   // full triangle up
    {"flat_t95_nt50", "t95n50", "flat t=[0.95,1], nt=[0.50,0.95]",
                                 kMagenta + 2, 23, false},   // full triangle down
  };

  // ------------------------------------------------------------------------
  // Load input graphs.
  // ------------------------------------------------------------------------
  std::vector<TGraphErrors *>      gABCD(scan.size(), nullptr);
  std::vector<TGraphAsymmErrors *> gTruth(scan.size(), nullptr);
  std::vector<TGraphErrors *>      gDiff(scan.size(), nullptr);

  for (size_t i = 0; i < scan.size(); ++i)
  {
    const std::string fname =
        basepath + "Photon_final_bdt_" + scan[i].tag + "_mc.root";
    TFile *f = TFile::Open(fname.c_str());
    if (!f || f->IsZombie())
    {
      std::cerr << "ERROR: cannot open " << fname << std::endl;
      return;
    }
    TGraphErrors      *gl = (TGraphErrors *)     f->Get("gpurity_leak");
    TGraphAsymmErrors *gt = (TGraphAsymmErrors *)f->Get("g_purity_truth");
    if (!gl || !gt)
    {
      std::cerr << "ERROR: missing gpurity_leak or g_purity_truth in "
                << fname << std::endl;
      return;
    }
    gl = (TGraphErrors *)     gl->Clone(Form("gpurity_leak_%s",  scan[i].suffix.c_str()));
    gt = (TGraphAsymmErrors *)gt->Clone(Form("g_purity_truth_%s", scan[i].suffix.c_str()));
    gl->SetTitle("");
    gt->SetTitle("");

    gABCD[i]  = gl;
    gTruth[i] = gt;
  }

  // ------------------------------------------------------------------------
  // Build non-closure curves -- each variant uses its OWN truth, since region
  // A depends on the tight cut.  (We still build one for the reference so it
  // is archived, but it will not be drawn in the bottom panel.)
  // ------------------------------------------------------------------------
  for (size_t i = 0; i < scan.size(); ++i)
  {
    gDiff[i] = buildNonClosure(gABCD[i], gTruth[i],
                               Form("g_nonclosure_%s", scan[i].suffix.c_str()));
  }

  // ------------------------------------------------------------------------
  // Print the 12 NC values per variant for user inspection.
  // ------------------------------------------------------------------------
  std::cout << "---------- ABCD non-closure (gpurity_leak - truth) ----------\n";
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "  pT[GeV]";
  for (size_t i = 0; i < scan.size(); ++i)
    std::cout << "   " << std::setw(22) << scan[i].suffix;
  std::cout << "\n";
  const int nbins = gDiff[0]->GetN();
  for (int k = 0; k < nbins; ++k)
  {
    double x0, y0;
    gDiff[0]->GetPoint(k, x0, y0);
    std::cout << "  " << std::setw(6) << x0;
    for (size_t i = 0; i < scan.size(); ++i)
    {
      double xi, yi;
      gDiff[i]->GetPoint(k, xi, yi);
      std::cout << "   " << std::setw(22) << yi;
    }
    std::cout << "\n";
  }
  std::cout << "-------------------------------------------------------------\n";

  // ------------------------------------------------------------------------
  // Slightly x-shift each curve to avoid marker overlap.
  // 5 variants: offsets -0.30, -0.15, 0.00, +0.15, +0.30 GeV.
  // Order follows the scan vector so colors/markers stay matched.
  // ------------------------------------------------------------------------
  const std::array<double, 5> dx_map = { 0.0, -0.30, -0.15, +0.15, +0.30 };
  std::vector<TGraphErrors *> gABCDs(scan.size());
  std::vector<TGraphErrors *> gDiffs(scan.size());
  for (size_t i = 0; i < scan.size(); ++i)
  {
    const double dx = dx_map[i];
    gABCDs[i] = shiftGraphX(gABCD[i], dx,
                            Form("gpurity_leak_%s_sh", scan[i].suffix.c_str()));
    gDiffs[i] = shiftGraphX(gDiff[i], dx,
                            Form("g_nonclosure_%s_sh", scan[i].suffix.c_str()));
  }

  // ------------------------------------------------------------------------
  // Canvas with stacked pads (shared x-axis).
  // ------------------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c_nonclosure_flat", "", 700, 900);
  c1->Divide(1, 2);

  TPad *pad_top = (TPad *)c1->cd(1);
  pad_top->SetPad(0, 0.38, 1, 1);
  pad_top->SetTopMargin(0.10);
  pad_top->SetBottomMargin(0.004);
  pad_top->SetLeftMargin(0.17);
  pad_top->SetRightMargin(0.04);

  TH1F *frame_top = (TH1F *)frame_et_rec->Clone("frame_top_flat");
  frame_top->SetTitle("");
  frame_top->SetYTitle("Purity");
  frame_top->SetXTitle("");
  frame_top->GetYaxis()->SetRangeUser(0.0, 1.20);
  frame_top->GetXaxis()->SetRangeUser(10, 35);
  frame_top->GetYaxis()->SetTitleOffset(1.2);
  frame_top->GetYaxis()->SetTitleSize(0.055);
  frame_top->GetYaxis()->SetLabelSize(0.050);
  frame_top->GetXaxis()->SetLabelSize(0.0);
  frame_top->GetXaxis()->SetTitleSize(0.0);
  frame_top->GetYaxis()->SetNdivisions(508);
  frame_top->GetXaxis()->SetNdivisions(510);
  frame_top->Draw("axis");

  lineone->SetLineColor(kGray + 2);
  lineone->SetLineStyle(2);
  lineone->Draw("L same");

  // Draw truth curves first (dashed lines, same color as variant).
  // Skip the parametric-nominal truth to avoid visual clutter -- the nominal
  // entry exists in the top panel for its ABCD points only.
  std::vector<TGraphAsymmErrors *> gTruthDraw(scan.size(), nullptr);
  for (size_t i = 0; i < scan.size(); ++i)
  {
    if (scan[i].is_ref) continue;
    TGraphAsymmErrors *gt = (TGraphAsymmErrors *)gTruth[i]->Clone(
        Form("g_purity_truth_%s_draw", scan[i].suffix.c_str()));
    gt->SetTitle("");
    gt->SetLineColor(scan[i].color);
    gt->SetLineStyle(2);
    gt->SetLineWidth(3);
    gt->SetMarkerStyle(0);
    gt->SetFillStyle(0);
    gt->Draw("L same");
    gTruthDraw[i] = gt;
  }

  // Draw ABCD points.
  for (size_t i = 0; i < scan.size(); ++i)
  {
    gABCDs[i]->SetMarkerColor(scan[i].color);
    gABCDs[i]->SetLineColor(scan[i].color);
    gABCDs[i]->SetMarkerStyle(scan[i].marker);
    gABCDs[i]->SetMarkerSize(1.1);
    gABCDs[i]->SetLineWidth(2);
    gABCDs[i]->Draw("P same");
  }

  // sPHENIX header block (pad top-margin = 0.10 leaves NDC [0.90,1.0]).
  myText(0.20, 0.955, 1, strleg1.c_str(), 0.050, 0);
  myText(0.20, 0.908, 1,
         "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, jet MC (ABCD self-closure)",
         0.042, 0);
  myText(0.95, 0.955, 1, strleg3.c_str(), 0.046, 1);
  // strleg4 placed inside frame at y=0.855 to avoid frame-edge clipping.
  myText(0.20, 0.810, 1, strleg4.c_str(), 0.038, 0);

  // Legend in the bottom-right of the top panel (below the data bulk).
  TLegend *leg = new TLegend(0.50, 0.11, 0.94, 0.38);
  legStyle(leg, 0.14, 0.032);
  for (size_t i = 0; i < scan.size(); ++i)
    leg->AddEntry(gABCDs[i], scan[i].label.c_str(), "pl");
  // One representative truth entry (dashed line, color-neutral label).
  if (gTruthDraw[1])
    leg->AddEntry(gTruthDraw[1], "truth (per-variant region A)", "l");
  leg->Draw("same");

  // ------------------------------------------------------------------------
  // Bottom pad: non-closure = gpurity_leak - truth (per variant).
  // ------------------------------------------------------------------------
  TPad *pad_bot = (TPad *)c1->cd(2);
  pad_bot->SetPad(0, 0, 1, 0.38);
  pad_bot->SetTopMargin(0.05);
  pad_bot->SetBottomMargin(0.28);
  pad_bot->SetLeftMargin(0.17);
  pad_bot->SetRightMargin(0.04);

  TH1F *frame_bot = (TH1F *)frame_et_rec->Clone("frame_bot_flat");
  frame_bot->SetTitle("");
  frame_bot->SetXTitle("#it{E}_{T}^{#gamma,rec} [GeV]");
  frame_bot->SetYTitle("#varepsilon_{ABCD} #minus #varepsilon_{truth}");
  frame_bot->GetXaxis()->SetRangeUser(10, 35);
  frame_bot->GetYaxis()->SetRangeUser(-0.50, 0.45);
  frame_bot->GetYaxis()->SetNdivisions(505);
  frame_bot->GetXaxis()->SetNdivisions(510);
  const double s = 6.2 / 3.8;
  frame_bot->GetXaxis()->SetTitleSize(0.055 * s);
  frame_bot->GetYaxis()->SetTitleSize(0.055 * s);
  frame_bot->GetXaxis()->SetLabelSize(0.050 * s);
  frame_bot->GetYaxis()->SetLabelSize(0.050 * s);
  frame_bot->GetYaxis()->SetTitleOffset(0.95);
  frame_bot->GetXaxis()->SetTitleOffset(1.05);
  frame_bot->Draw("axis");

  linezero->SetLineColor(kGray + 2);
  linezero->SetLineStyle(2);
  linezero->Draw("L same");

  // Bottom panel: 4 flat variants only (skip parametric reference).
  for (size_t i = 0; i < scan.size(); ++i)
  {
    if (scan[i].is_ref) continue;
    gDiffs[i]->SetMarkerColor(scan[i].color);
    gDiffs[i]->SetLineColor(scan[i].color);
    gDiffs[i]->SetMarkerStyle(scan[i].marker);
    gDiffs[i]->SetMarkerSize(1.1);
    gDiffs[i]->SetLineWidth(2);
    gDiffs[i]->Draw("P same");
  }

  // Annotation clarifying per-variant truth.
  TPaveText *pt = new TPaveText(0.18, 0.82, 0.72, 0.94, "brNDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextFont(42);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.055);
  pt->AddText("truth uses each variant's own region A");
  pt->Draw("same");

  // ------------------------------------------------------------------------
  // Save PDF + archive.
  // ------------------------------------------------------------------------
  c1->SaveAs(outpdf);

  TFile *fout = TFile::Open(outroot.c_str(), "RECREATE");
  for (size_t i = 0; i < scan.size(); ++i)
  {
    gABCD[i]->Write(Form("gpurity_leak_%s",    scan[i].suffix.c_str()));
    gTruth[i]->Write(Form("g_purity_truth_%s", scan[i].suffix.c_str()));
    gDiff[i]->Write(Form("g_nonclosure_%s",    scan[i].suffix.c_str()));
  }
  fout->Close();

  std::cout << "[plot_purity_nonclosure_flat] wrote PDF:  " << outpdf  << std::endl;
  std::cout << "[plot_purity_nonclosure_flat] wrote ROOT: " << outroot << std::endl;
}
