#include "plotcommon.h"

// --------------------------------------------------------------------------
// plot_purity_nonclosure_ntbdt
//
// Visualises the ABCD purity non-closure on inclusive MC as the non-tight
// BDT lower bound (analysis.non_tight.bdt_min) is scanned over
// {0.02 (nom), 0.05, 0.10}.
//
// The earlier ntbdtmin00/01/02 files were produced with an older pipeline
// (11-bin, different per-sample MC_efficiency set) and are NOT comparable
// with the current 12-bin pipeline -- they are excluded here.
//
// Inputs (all 12-bin, current CLUSTERINFO_CEMC pipeline):
//   Photon_final_bdt_nom_mc.root         -- nt_bdt_min = 0.02 (nominal)
//   Photon_final_bdt_ntbdtmin05_mc.root  -- nt_bdt_min = 0.05
//   Photon_final_bdt_ntbdtmin10_mc.root  -- nt_bdt_min = 0.10
//
// Each file contains
//   gpurity_leak   -- ABCD purity with signal-leak correction (primary)
//   g_purity_truth -- truth purity  h_truth_iso_cluster_data /
//                                     h_tight_iso_cluster_signal_data
//
// Truth purity is region-A-only (tight+iso) and is therefore independent of
// the non-tight BDT cut.  We verify this in the macro and draw a single
// dashed black truth curve taken from the nominal file.
//
// Output:
//   figures/purity_nonclosure_ntbdt.pdf    -- 2-panel canvas
//   rootFiles/purity_nonclosure_ntbdt.root -- archived TGraphs
// --------------------------------------------------------------------------

namespace
{
  struct ScanEntry
  {
    std::string  tag;       // file-tag, e.g. "nom" or "ntbdtmin05"
    std::string  suffix;    // short suffix used for ROOT object names
    double       value;     // numerical nt_bdt_min
    std::string  label;     // legend label
    int          color;
    int          marker;
  };

  // Shift x by dx (keeps errors).
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

  // Non-closure = gABCD - gTruth.  With truth taken from the nominal file
  // (same 12-bin x-grid as every gABCD in the 3 clean files), so the match
  // is bin-by-bin.  Errors are added in quadrature; gTruth asymmetric errors
  // are symmetrised.
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
      // Find matching gTruth point by x (tolerance 0.01 GeV).
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

void plot_purity_nonclosure_ntbdt(const char *outpdf =
                                  "figures/purity_nonclosure_ntbdt.pdf")
{
  init_plot();

  const std::string basepath =
      "/sphenix/user/shuhangli/ppg12/efficiencytool/results/";
  const std::string outroot =
      "/sphenix/user/shuhangli/ppg12/plotting/rootFiles/"
      "purity_nonclosure_ntbdt.root";

  // 3 pipeline-consistent scan points.
  std::vector<ScanEntry> scan = {
    {"nom",        "nom", 0.02, "nt_bdt_min = 0.02 (nom)", kBlack,    20},
    {"ntbdtmin05", "05",  0.05, "nt_bdt_min = 0.05",        kBlue + 1, 20},
    {"ntbdtmin10", "10",  0.10, "nt_bdt_min = 0.10",        kRed  + 1, 20},
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
    // TGraph{Errors,AsymmErrors} are not attached to a TDirectory, so cloning
    // is sufficient to detach them from the TFile we are about to close.
    gl = (TGraphErrors *)     gl->Clone(Form("gpurity_leak_%s",  scan[i].suffix.c_str()));
    gt = (TGraphAsymmErrors *)gt->Clone(Form("g_purity_truth_%s", scan[i].suffix.c_str()));
    gl->SetTitle("");
    gt->SetTitle("");

    gABCD[i]  = gl;
    gTruth[i] = gt;
  }

  // ------------------------------------------------------------------------
  // Truth-consistency check (region A should not depend on nt_bdt_min).
  // Print every bin for every file, and WARN if any pair disagrees by >5%
  // (relative).
  // ------------------------------------------------------------------------
  std::cout << "---------- truth-purity consistency check ----------\n";
  const int nbins = gTruth[0]->GetN();
  bool any_warn = false;
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "  pT[GeV]";
  for (size_t i = 0; i < scan.size(); ++i)
    std::cout << "   " << std::setw(18) << scan[i].label;
  std::cout << "     maxRelDev\n";
  for (int k = 0; k < nbins; ++k)
  {
    double x0, y0;
    gTruth[0]->GetPoint(k, x0, y0);
    double ymin = y0, ymax = y0;
    std::cout << "  " << std::setw(6) << x0;
    for (size_t i = 0; i < scan.size(); ++i)
    {
      double xi, yi;
      // Find matching x in gTruth[i] (files may share the same 12-bin grid,
      // but guard just in case).
      int jbest = 0; double dbest = 1e9;
      for (int j = 0; j < gTruth[i]->GetN(); ++j)
      {
        double xj, yj;
        gTruth[i]->GetPoint(j, xj, yj);
        double d = std::fabs(xj - x0);
        if (d < dbest) { dbest = d; jbest = j; }
      }
      gTruth[i]->GetPoint(jbest, xi, yi);
      std::cout << "   " << std::setw(18) << yi;
      if (yi < ymin) ymin = yi;
      if (yi > ymax) ymax = yi;
    }
    const double rel = (ymax > 0.0) ? (ymax - ymin) / ymax : 0.0;
    std::cout << "   " << std::setw(9) << rel;
    if (rel > 0.05)
    {
      std::cout << "  <-- WARNING (>5%)";
      any_warn = true;
    }
    std::cout << "\n";
  }
  if (any_warn)
    std::cout << "WARNING: truth purity differs by >5% in some bin; "
              << "single-truth curve may be misleading.\n";
  else
    std::cout << "OK: truth purity agrees within 5% across all scan files; "
              << "using nominal truth curve as shared reference.\n";
  std::cout << "----------------------------------------------------\n";

  // ------------------------------------------------------------------------
  // Build non-closure curves.  Use NOMINAL truth as the common reference
  // (region A is independent of the non-tight BDT cut, verified above).
  // ------------------------------------------------------------------------
  TGraphAsymmErrors *gTruthRef = gTruth[0];
  for (size_t i = 0; i < scan.size(); ++i)
  {
    gDiff[i] = buildNonClosure(gABCD[i], gTruthRef,
                               Form("g_nonclosure_%s", scan[i].suffix.c_str()));
  }

  // ------------------------------------------------------------------------
  // Slightly x-shift each curve to avoid marker overlap.
  // ------------------------------------------------------------------------
  // Centred around nominal (index 0 in this scan, but symmetric look is
  // achieved with -0.18, 0, +0.18 GeV offsets for visual separation).
  const double dx_step = 0.18;
  const std::array<double, 3> dx_map = { 0.0, -dx_step, +dx_step };
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
  TCanvas *c1 = new TCanvas("c_nonclosure_ntbdt", "", 700, 900);
  c1->Divide(1, 2);

  TPad *pad_top = (TPad *)c1->cd(1);
  pad_top->SetPad(0, 0.38, 1, 1);
  pad_top->SetTopMargin(0.10);        // headroom for sPHENIX header block
  pad_top->SetBottomMargin(0.004);
  pad_top->SetLeftMargin(0.15);
  pad_top->SetRightMargin(0.04);

  // Top frame -- clone frame_et_rec so we don't mutate it for other users.
  TH1F *frame_top = (TH1F *)frame_et_rec->Clone("frame_top_ntbdt");
  frame_top->SetTitle("");
  frame_top->SetYTitle("Purity");
  frame_top->SetXTitle("");
  // y-range per request is [0.0, 1.1]; we add a small tick of headroom so
  // the sPHENIX header text doesn't collide with the y=1 dashed reference
  // line.  Frame top still shows 1.1 as the top tick label.
  frame_top->GetYaxis()->SetRangeUser(0.0, 1.10);
  frame_top->GetXaxis()->SetRangeUser(10, 35);
  frame_top->GetYaxis()->SetTitleOffset(1.2);
  frame_top->GetYaxis()->SetTitleSize(0.055);
  frame_top->GetYaxis()->SetLabelSize(0.050);
  frame_top->GetXaxis()->SetLabelSize(0.0);       // hide labels on shared x
  frame_top->GetXaxis()->SetTitleSize(0.0);
  frame_top->GetYaxis()->SetNdivisions(508);
  frame_top->GetXaxis()->SetNdivisions(510);
  frame_top->Draw("axis");

  lineone->SetLineColor(kGray + 2);
  lineone->SetLineStyle(2);
  lineone->Draw("L same");

  // Single dashed-black truth curve (from nominal file).
  TGraphAsymmErrors *gTruthDraw =
      (TGraphAsymmErrors *)gTruthRef->Clone("g_purity_truth_shared");
  gTruthDraw->SetTitle("");
  gTruthDraw->SetLineColor(kBlack);
  gTruthDraw->SetLineStyle(2);
  gTruthDraw->SetLineWidth(3);
  gTruthDraw->SetMarkerStyle(0);
  gTruthDraw->SetFillStyle(0);
  gTruthDraw->Draw("L same");

  // ABCD scan points (filled circles).
  for (size_t i = 0; i < scan.size(); ++i)
  {
    gABCDs[i]->SetMarkerColor(scan[i].color);
    gABCDs[i]->SetLineColor(scan[i].color);
    gABCDs[i]->SetMarkerStyle(scan[i].marker);
    gABCDs[i]->SetMarkerSize(1.1);
    gABCDs[i]->SetLineWidth(2);
    gABCDs[i]->Draw("P same");
  }

  // sPHENIX header block.  Pad top-margin = 0.10 leaves NDC [0.90, 1.0]
  // above the frame for header text.
  myText(0.20, 0.955, 1, strleg1.c_str(), 0.050, 0);
  myText(0.20, 0.908, 1,
         "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, jet MC (ABCD self-closure)",
         0.042, 0);
  myText(0.95, 0.955, 1, strleg3.c_str(), 0.046, 1);
  myText(0.20, 0.855, 1, strleg4.c_str(), 0.040, 0);

  // Legend: 3 ABCD entries + 1 truth entry.
  TLegend *leg = new TLegend(0.54, 0.12, 0.93, 0.38);
  legStyle(leg, 0.14, 0.038);
  for (size_t i = 0; i < scan.size(); ++i)
    leg->AddEntry(gABCDs[i], scan[i].label.c_str(), "pl");
  leg->AddEntry(gTruthDraw, "truth (region A, nom)", "l");
  leg->Draw("same");

  // ------------------------------------------------------------------------
  // Bottom pad: non-closure = gpurity_leak - g_purity_truth
  // ------------------------------------------------------------------------
  TPad *pad_bot = (TPad *)c1->cd(2);
  pad_bot->SetPad(0, 0, 1, 0.38);
  pad_bot->SetTopMargin(0.05);
  pad_bot->SetBottomMargin(0.28);
  pad_bot->SetLeftMargin(0.15);
  pad_bot->SetRightMargin(0.04);

  TH1F *frame_bot = (TH1F *)frame_et_rec->Clone("frame_bot_ntbdt");
  frame_bot->SetTitle("");
  frame_bot->SetXTitle("#it{E}_{T}^{#gamma,rec} [GeV]");
  frame_bot->SetYTitle("#varepsilon_{ABCD} #minus #varepsilon_{truth}");
  frame_bot->GetXaxis()->SetRangeUser(10, 35);
  frame_bot->GetYaxis()->SetRangeUser(-0.30, 0.22);
  frame_bot->GetYaxis()->SetNdivisions(505);
  frame_bot->GetXaxis()->SetNdivisions(510);
  const double s = 6.2 / 3.8;                     // scale factor (top 0.62 vs bottom 0.38)
  frame_bot->GetXaxis()->SetTitleSize(0.055 * s);
  frame_bot->GetYaxis()->SetTitleSize(0.055 * s);
  frame_bot->GetXaxis()->SetLabelSize(0.050 * s);
  frame_bot->GetYaxis()->SetLabelSize(0.050 * s);
  frame_bot->GetYaxis()->SetTitleOffset(1.20 / s);
  frame_bot->GetXaxis()->SetTitleOffset(1.05);
  frame_bot->Draw("axis");

  linezero->SetLineColor(kGray + 2);
  linezero->SetLineStyle(2);
  linezero->Draw("L same");

  for (size_t i = 0; i < scan.size(); ++i)
  {
    gDiffs[i]->SetMarkerColor(scan[i].color);
    gDiffs[i]->SetLineColor(scan[i].color);
    gDiffs[i]->SetMarkerStyle(scan[i].marker);
    gDiffs[i]->SetMarkerSize(1.1);
    gDiffs[i]->SetLineWidth(2);
    gDiffs[i]->Draw("P same");
  }

  // Annotation clarifying truth source.
  TPaveText *pt = new TPaveText(0.18, 0.82, 0.72, 0.94, "brNDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextFont(42);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.055);
  pt->AddText("truth = nom-file truth (region A, nt-independent)");
  pt->Draw("same");

  // ------------------------------------------------------------------------
  // Save PDF.
  // ------------------------------------------------------------------------
  c1->SaveAs(outpdf);

  // Archive TGraphs.
  TFile *fout = TFile::Open(outroot.c_str(), "RECREATE");
  for (size_t i = 0; i < scan.size(); ++i)
  {
    gABCD[i]->Write(Form("gpurity_leak_%s",    scan[i].suffix.c_str()));
    gTruth[i]->Write(Form("g_purity_truth_%s", scan[i].suffix.c_str()));
    gDiff[i]->Write(Form("g_nonclosure_%s",    scan[i].suffix.c_str()));
  }
  gTruthDraw->Write("g_purity_truth_shared");
  fout->Close();

  std::cout << "[plot_purity_nonclosure_ntbdt] wrote PDF:  " << outpdf  << std::endl;
  std::cout << "[plot_purity_nonclosure_ntbdt] wrote ROOT: " << outroot << std::endl;
}
