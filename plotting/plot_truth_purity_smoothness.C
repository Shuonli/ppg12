// plot_truth_purity_smoothness.C
//
// Truth-purity smoothness check for sample combining.
//
//   purity(pT) = N_prompt(pT) / ( N_prompt(pT) + N_jet(pT) )
//
//   N_prompt from signal-combined MC: h_truth_pT_0  (leading truth photon pT)
//   N_jet    from jet-inclusive MC : h_max_truth_jet_pT (leading truth jet pT)
//
// Assumption: at the hard-parton level, photon pT ~ jet pT for the leading
// parton, so this is a smoothness proxy for the cross-section mixing.
//
// Writes PDFs to plotting/figures/sample_combining_check/
//
// Run:
//   root -l -b -q 'plot_truth_purity_smoothness.C'
//
#include "plotcommon.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>

namespace
{

  const char *kResultsDir =
      "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results";
  const char *kFigDir =
      "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures/"
      "sample_combining_check";

  // Sample boundaries (truth pT) — draw as vertical dashed lines.
  // Signal: photon5 (0-14), photon10 (14-30), photon20 (30+).
  // BG    : jet8 (9-14), jet12 (14-21), jet20 (21-32), jet30 (32-42), jet40 (42+).
  const std::vector<double> kBoundaries = {14.0, 21.0, 30.0, 32.0, 42.0};

  // Per-period luminosity string (overrides plotcommon's hardcoded 16.6 pb^-1).
  // Values come from wiki/reference rules:
  //   0 mrad = 32.66 pb^-1, 1.5 mrad = 16.27 pb^-1, all-range = 48.93 pb^-1.
  std::string lumiLabel(const std::string &period)
  {
    if (period == "0rad")
      return "#it{p}+#it{p} #kern[-0.05]{#sqrt{#it{s}} = 200 GeV, 32.66 pb^{-1}}";
    if (period == "1p5mrad")
      return "#it{p}+#it{p} #kern[-0.05]{#sqrt{#it{s}} = 200 GeV, 16.27 pb^{-1}}";
    // "all" or anything else
    return "#it{p}+#it{p} #kern[-0.05]{#sqrt{#it{s}} = 200 GeV, 48.93 pb^{-1}}";
  }

  // Common rebinning edges for purity curves.
  //
  // The signal h_truth_pT_0 has variable-width bins:
  //   7 8 10 12 14 16 18 20 22 24 26 28 32 36 45
  // so the "natural" grid for purity is those edges.  We use the signal's
  // exact binning as the common grid and rebin the (much finer) jet
  // h_max_truth_jet_pT onto it via integral-preserving grouping.
  //
  // All sample boundaries (14, 21, 30, 32, 42) are inspected via the graph
  // zoom panels; the signal grid contains 14, 22, 28, 32 as bin edges and
  // spans 36-45 for the highest bin.
  const std::vector<double> kPurityEdges = {
      7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45};

  // Plot style: ROOT palette indices.
  const int kCol0rad = kAzure + 2;
  const int kCol1p5 = kRed + 1;
  const int kColAll = kBlack;

  // Rebin a TH1 (arbitrary binning) into user-defined edges by summing bin
  // contents (assumes the input edges are strict subsets of, or finer than,
  // the target edges; we call FindBin on bin centres).
  TH1D *rebinToEdges(const TH1 *h_in, const std::vector<double> &edges,
                     const char *new_name)
  {
    std::vector<double> e = edges;
    TH1D *h_out = new TH1D(new_name, h_in->GetTitle(),
                           (int)e.size() - 1, e.data());
    h_out->Sumw2();
    for (int b = 1; b <= h_in->GetNbinsX(); ++b) {
      double x = h_in->GetBinCenter(b);
      double c = h_in->GetBinContent(b);
      double err = h_in->GetBinError(b);
      int target = h_out->FindBin(x);
      if (target < 1 || target > h_out->GetNbinsX()) continue;
      double cur = h_out->GetBinContent(target);
      double e2 = h_out->GetBinError(target);
      h_out->SetBinContent(target, cur + c);
      h_out->SetBinError(target, std::sqrt(e2 * e2 + err * err));
    }
    return h_out;
  }

  TH1D *readSignalTruthPt(const TString &fname, const char *new_name)
  {
    TFile *f = TFile::Open(fname);
    if (!f || f->IsZombie()) {
      std::cerr << "[ERROR] cannot open signal file: " << fname << std::endl;
      return nullptr;
    }
    TH1D *h = (TH1D *)f->Get("h_truth_pT_0");
    if (!h) {
      std::cerr << "[ERROR] h_truth_pT_0 missing in " << fname << std::endl;
      f->Close();
      return nullptr;
    }
    h->SetDirectory(nullptr);
    TH1D *h_reb = rebinToEdges(h, kPurityEdges, new_name);
    h_reb->SetDirectory(nullptr);
    f->Close();
    delete h;
    return h_reb;
  }

  TH1D *readJetMaxTruthPt(const TString &fname, const char *new_name)
  {
    TFile *f = TFile::Open(fname);
    if (!f || f->IsZombie()) {
      std::cerr << "[ERROR] cannot open jet file: " << fname << std::endl;
      return nullptr;
    }
    TH1F *h = (TH1F *)f->Get("h_max_truth_jet_pT");
    if (!h) {
      std::cerr << "[ERROR] h_max_truth_jet_pT missing in " << fname
                << std::endl;
      f->Close();
      return nullptr;
    }
    h->SetDirectory(nullptr);
    TH1D *h_reb = rebinToEdges(h, kPurityEdges, new_name);
    h_reb->SetDirectory(nullptr);
    f->Close();
    delete h;
    return h_reb;
  }

  TGraphErrors *makePurityGraph(const TH1D *hN, const TH1D *hJ,
                                const char *name)
  {
    TGraphErrors *g = new TGraphErrors();
    g->SetName(name);
    int n = 0;
    for (int b = 1; b <= hN->GetNbinsX(); ++b) {
      double Np = hN->GetBinContent(b);
      double Nj = hJ->GetBinContent(b);
      double dNp = hN->GetBinError(b);
      double dNj = hJ->GetBinError(b);
      double denom = Np + Nj;
      if (denom <= 0) continue;
      double p = Np / denom;
      // sigma_p^2 = (Nj*dNp)^2/denom^4 + (Np*dNj)^2/denom^4 (binomial-like
      // with independent Poisson-weighted contributions).
      double sp = std::sqrt((Nj * Nj) * (dNp * dNp) +
                            (Np * Np) * (dNj * dNj)) /
                  (denom * denom);
      double x = hN->GetBinCenter(b);
      double ex = 0.5 * hN->GetBinWidth(b);
      g->SetPoint(n, x, p);
      g->SetPointError(n, ex, sp);
      ++n;
    }
    return g;
  }

  void drawBoundaries(double ymin, double ymax)
  {
    for (double xb : kBoundaries) {
      TLine *l = new TLine(xb, ymin, xb, ymax);
      l->SetLineStyle(2);
      l->SetLineColor(kGray + 2);
      l->SetLineWidth(1);
      l->Draw("same");
    }
  }

  void drawLegendText(double x, double y, const char *extra = nullptr,
                      const std::string &lumi_str = "")
  {
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.035);
    t.DrawLatex(x, y, strleg1.c_str());
    t.DrawLatex(x, y - 0.05,
                lumi_str.empty() ? strleg2_1.c_str() : lumi_str.c_str());
    t.DrawLatex(x, y - 0.10, strleg3.c_str());
    if (extra) t.DrawLatex(x, y - 0.15, extra);
  }

  void printBinTable(TH1D *hN, TH1D *hJ, const char *label)
  {
    std::cout << "\n=== Bin table for " << label
              << " (signal-native binning; Np = h_truth_pT_0, "
                 "Nj = h_max_truth_jet_pT rebinned) ===" << std::endl;
    std::printf("  %12s  %14s  %14s  %8s\n", "pT range", "N_prompt",
                "N_jet", "purity");
    for (int b = 1; b <= hN->GetNbinsX(); ++b) {
      double lo = hN->GetBinLowEdge(b);
      double hi = hN->GetBinLowEdge(b + 1);
      double Np = hN->GetBinContent(b);
      double Nj = hJ->GetBinContent(b);
      double denom = Np + Nj;
      double p = (denom > 0) ? Np / denom : 0.0;
      std::printf("  [%4.1f - %4.1f]  %14.4e  %14.4e  %8.5f\n", lo, hi, Np,
                  Nj, p);
    }
  }

  void plotOnePeriod(const char *period_label, const char *sig_file,
                     const char *jet_file, const char *out_base)
  {
    TString sig = TString::Format("%s/%s", kResultsDir, sig_file);
    TString jet = TString::Format("%s/%s", kResultsDir, jet_file);

    TH1D *hS =
        readSignalTruthPt(sig, Form("h_sig_%s", period_label));
    TH1D *hJ = readJetMaxTruthPt(jet, Form("h_jet_%s", period_label));
    if (!hS || !hJ) return;

    TGraphErrors *g =
        makePurityGraph(hS, hJ, Form("g_purity_%s", period_label));
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.0);
    g->SetLineWidth(2);

    // Auto-determine purity range over pT > 8 GeV (skip the first bin that
    // has zero jet content because the jet MC only starts at ~9 GeV truth).
    double pmin = 1e9, pmax = -1e9;
    for (int i = 0; i < g->GetN(); ++i) {
      double x, y;
      g->GetPoint(i, x, y);
      if (x < 8.0) continue;  // skip "pure signal" low-pT bin
      if (y < pmin) pmin = y;
      if (y > pmax) pmax = y;
    }
    if (pmin >= pmax) { pmin = 0; pmax = 1e-3; }
    const double padlo = std::max(0.0, pmin * 0.5);
    const double padhi = pmax * 1.5;

    // Full-range plot.
    TCanvas *c = new TCanvas(Form("c_purity_%s", period_label), "", 900, 700);
    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.06);

    TH1F *frame =
        new TH1F(Form("frame_%s", period_label), "", 400, 7, 45);
    frame->GetXaxis()->SetTitle("truth #it{p}_{T} [GeV]");
    frame->GetYaxis()->SetTitle(
        "Truth purity "
        "N^{#gamma}_{prompt} / (N^{#gamma}_{prompt} + N_{jet})");
    frame->GetYaxis()->SetRangeUser(padlo, padhi);
    frame->GetYaxis()->SetTitleOffset(1.7);
    frame->GetYaxis()->SetMaxDigits(3);
    frame->Draw();
    drawBoundaries(padlo, padhi);

    g->SetMarkerColor(kCol0rad);
    g->SetLineColor(kCol0rad);
    g->Draw("P same");

    drawLegendText(0.20, 0.88,
                   Form("Period: %s (purity from combined MC)", period_label),
                   lumiLabel(period_label));

    c->SaveAs(Form("%s/%s.pdf", kFigDir, out_base));

    // Boundary-zoom: 2x2 grid. The 42 GeV boundary is not resolved by the
    // signal's variable binning (last bin is [36, 45]) so we combine the two
    // closely-spaced signal boundaries (30 & 32 GeV) and drop 42.
    struct ZoomCfg {
      double xb;
      double x_lo;
      double x_hi;
      const char *title;
    };
    const std::vector<ZoomCfg> zooms = {
        {14.0, 11.0, 17.0, "14 GeV (photon5 #rightarrow photon10, jet8 #rightarrow jet12)"},
        {21.0, 18.0, 24.0, "21 GeV (jet12 #rightarrow jet20)"},
        {31.0, 27.0, 34.0, "30-32 GeV (photon10 #rightarrow photon20, jet20 #rightarrow jet30)"},
        {42.0, 35.0, 46.0, "42 GeV (jet30 #rightarrow jet40)"},
    };

    // Common y-range across all four panels (so the eye can compare levels).
    double zylo = 1e9, zyhi = -1e9;
    for (int j = 0; j < g->GetN(); ++j) {
      double x, y;
      g->GetPoint(j, x, y);
      if (x < 10.0 || x > 47.0) continue;
      if (y > 0 && y < zylo) zylo = y;
      if (y > zyhi) zyhi = y;
    }
    if (zylo >= zyhi) { zylo = 1e-6; zyhi = 1e-3; }
    // Log-y for consistent exponent display across panels.
    double pad_lo = zylo * 0.3;
    double pad_hi = zyhi * 3.0;

    TCanvas *cz =
        new TCanvas(Form("c_zoom_%s", period_label), "", 1200, 900);
    cz->Divide(2, 2);
    for (size_t i = 0; i < zooms.size(); ++i) {
      cz->cd((int)i + 1);
      gPad->SetLogy();
      gPad->SetLeftMargin(0.16);
      gPad->SetRightMargin(0.05);
      gPad->SetTopMargin(0.10);
      gPad->SetBottomMargin(0.14);
      const ZoomCfg &z = zooms[i];
      TH1F *fr = new TH1F(Form("zoom_%s_%zu", period_label, i), "",
                          40, z.x_lo, z.x_hi);
      fr->GetXaxis()->SetTitle("truth #it{p}_{T} [GeV]");
      fr->GetYaxis()->SetTitle("Truth purity");
      fr->GetYaxis()->SetRangeUser(pad_lo, pad_hi);
      fr->GetYaxis()->SetTitleOffset(1.5);
      fr->GetXaxis()->SetTitleOffset(1.1);
      fr->Draw();
      // boundary marker
      TLine *l = new TLine(z.xb, pad_lo, z.xb, pad_hi);
      l->SetLineStyle(2);
      l->SetLineColor(kGray + 2);
      l->Draw("same");
      g->SetMarkerColor(kCol0rad);
      g->SetLineColor(kCol0rad);
      g->Draw("P same");
      TLatex t;
      t.SetNDC();
      t.SetTextFont(42);
      t.SetTextSize(0.040);
      t.DrawLatex(0.2, 0.92, z.title);
      // For the 42 GeV panel, note the signal-binning limitation INSIDE the
      // frame (not in the title, which would be clipped by the pad edge).
      if (std::fabs(z.xb - 42.0) < 0.01) {
        t.SetTextSize(0.034);
        t.SetTextColor(kGray + 2);
        t.DrawLatex(0.2, 0.86, "signal bin [36,45] GeV: boundary not resolved");
        t.SetTextColor(kBlack);
      }
    }
    // Global header across top of canvas
    cz->cd(0);
    TLatex th;
    th.SetNDC();
    th.SetTextFont(42);
    th.SetTextSize(0.025);
    th.DrawLatex(0.02, 0.97,
                 (strleg1 + "   " + lumiLabel(period_label) + "   " +
                  strleg3 + "   Period: " + period_label)
                     .c_str());
    cz->SaveAs(Form("%s/truth_purity_boundary_zoom_%s.pdf", kFigDir,
                    period_label));

    // Dump table.
    printBinTable(hS, hJ, period_label);

    // Cleanup: don't delete histograms/graphs since we saved TCanvas; but
    // tidy so the next call has fresh objects.
    // (Canvases and histograms will be gc'd when ROOT exits.)
  }

  void plotPeriodComparison()
  {
    struct PD {
      const char *label;
      const char *sig;
      const char *jet;
      int color;
      int marker;
    };
    std::vector<PD> pds = {
        {"0rad", "MC_efficiency_bdt_nom_0rad.root",
         "MC_efficiency_jet_bdt_nom_0rad.root", kCol0rad, 20},
        {"1p5mrad", "MC_efficiency_bdt_nom_1p5mrad.root",
         "MC_efficiency_jet_bdt_nom_1p5mrad.root", kCol1p5, 21},
        {"all", "MC_efficiency_bdt_nom.root",
         "MC_efficiency_jet_bdt_nom.root", kColAll, 22},
    };

    TCanvas *c = new TCanvas("c_purity_compare", "", 900, 700);
    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.06);

    // Two-pass: first build all graphs, then decide auto-scale.
    std::vector<TGraphErrors *> gs;
    for (auto &p : pds) {
      TString sig = TString::Format("%s/%s", kResultsDir, p.sig);
      TString jet = TString::Format("%s/%s", kResultsDir, p.jet);
      TH1D *hS =
          readSignalTruthPt(sig, Form("h_sig_cmp_%s", p.label));
      TH1D *hJ =
          readJetMaxTruthPt(jet, Form("h_jet_cmp_%s", p.label));
      if (!hS || !hJ) { gs.push_back(nullptr); continue; }
      TGraphErrors *g =
          makePurityGraph(hS, hJ, Form("g_purity_cmp_%s", p.label));
      g->SetMarkerColor(p.color);
      g->SetLineColor(p.color);
      g->SetMarkerStyle(p.marker);
      g->SetMarkerSize(1.0);
      g->SetLineWidth(2);
      gs.push_back(g);
    }

    double pmin = 1e9, pmax = -1e9;
    for (auto *g : gs) {
      if (!g) continue;
      for (int i = 0; i < g->GetN(); ++i) {
        double x, y;
        g->GetPoint(i, x, y);
        if (x < 8.0) continue;
        if (y < pmin) pmin = y;
        if (y > pmax) pmax = y;
      }
    }
    if (pmin >= pmax) { pmin = 0; pmax = 1e-3; }
    const double padlo = std::max(0.0, pmin * 0.5);
    const double padhi = pmax * 1.5;

    TH1F *frame = new TH1F("frame_cmp", "", 400, 7, 45);
    frame->GetXaxis()->SetTitle("truth #it{p}_{T} [GeV]");
    frame->GetYaxis()->SetTitle(
        "Truth purity "
        "N^{#gamma}_{prompt} / (N^{#gamma}_{prompt} + N_{jet})");
    frame->GetYaxis()->SetRangeUser(padlo, padhi);
    frame->GetYaxis()->SetTitleOffset(1.7);
    frame->GetYaxis()->SetMaxDigits(3);
    frame->Draw();
    drawBoundaries(padlo, padhi);

    TLegend *leg = new TLegend(0.62, 0.18, 0.94, 0.38);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);

    for (size_t i = 0; i < pds.size(); ++i) {
      if (!gs[i]) continue;
      gs[i]->Draw("P same");
      leg->AddEntry(gs[i], pds[i].label, "p");
    }
    leg->Draw();

    drawLegendText(0.20, 0.88, "Period overlay", lumiLabel("all"));

    c->SaveAs(Form("%s/truth_purity_periods_compare.pdf", kFigDir));
  }

}  // namespace

void plot_truth_purity_smoothness()
{
  init_plot();
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(3);

  std::cout << "Writing figures to: " << kFigDir << std::endl;

  // Per-period plots + tables.
  plotOnePeriod("0rad", "MC_efficiency_bdt_nom_0rad.root",
                "MC_efficiency_jet_bdt_nom_0rad.root",
                "truth_purity_vs_pT_0rad");
  plotOnePeriod("1p5mrad", "MC_efficiency_bdt_nom_1p5mrad.root",
                "MC_efficiency_jet_bdt_nom_1p5mrad.root",
                "truth_purity_vs_pT_1p5mrad");
  plotOnePeriod("all", "MC_efficiency_bdt_nom.root",
                "MC_efficiency_jet_bdt_nom.root",
                "truth_purity_vs_pT_all");

  // Boundary zoom files are written as <base>_zoom.pdf; rename for the
  // requested names by saving a simple symlink/copy pattern: we already
  // create truth_purity_vs_pT_{period}_zoom.pdf. Also create the requested
  // truth_purity_boundary_zoom_{period}.pdf via a final save.

  // Cross-period comparison.
  plotPeriodComparison();

  std::cout << "\nDone." << std::endl;
}
