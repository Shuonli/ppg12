// Validate the signal sample combining for the PPG12 nominal analysis.
//
// For each period {0rad, 1p5mrad, ""} the macro reads the 6 per-sample
// MC_efficiency_photonX_{nom,double}_bdt_nom_{period}.root files and the
// hadd-combined MC_efficiency_bdt_nom_{period}.root, pulls the full-eta
// truth photon pT histogram (h_truth_pT_0), and produces:
//
//   (A) signal_truth_pT_overlay_{period}.pdf   — log-y overlay of the 6
//        per-sample spectra + the combined spectrum
//   (B) signal_ratio_per_sample_{period}.pdf   — individual fractions
//        and the sum-of-per-sample / combined sanity ratio
//   (C) signal_boundary_zoom_{period}.pdf      — zoom at 14 and 30 GeV
//
// All figures follow PPG12 plotcommon.h conventions.

#include "plotcommon.h"

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TColor.h>

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

const char *kResultsDir = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results";
const char *kFigDir     = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures/sample_combining_check";

struct SampleInfo {
  std::string name;     // e.g. "photon5_nom"
  std::string label;    // e.g. "photon5 nom"
  int         color;
  int         style;    // line style
};

const std::vector<SampleInfo> kSamples = {
  {"photon5_nom",    "photon5 nom",    kAzure + 2,   1},
  {"photon5_double", "photon5 DI",     kAzure + 2,   2},
  {"photon10_nom",   "photon10 nom",   kGreen + 2,   1},
  {"photon10_double","photon10 DI",    kGreen + 2,   2},
  {"photon20_nom",   "photon20 nom",   kOrange + 7,  1},
  {"photon20_double","photon20 DI",    kOrange + 7,  2},
};

std::string makePath(const std::string &tag, const std::string &period)
{
  // tag="photon5_nom", period="0rad" -> .../MC_efficiency_photon5_nom_bdt_nom_0rad.root
  // tag="",            period="0rad" -> .../MC_efficiency_bdt_nom_0rad.root (combined)
  std::string name = "MC_efficiency_";
  if (!tag.empty()) name += tag + "_";
  name += "bdt_nom";
  if (!period.empty()) name += "_" + period;
  name += ".root";
  return std::string(kResultsDir) + "/" + name;
}

TH1D *fetchHist(TFile *f, const std::string &hname, const std::string &cloneName)
{
  TH1D *h = (TH1D *)f->Get(hname.c_str());
  if (!h) {
    std::cerr << "[fetchHist] " << hname << " not found in " << f->GetName() << std::endl;
    return nullptr;
  }
  TH1D *c = (TH1D *)h->Clone(cloneName.c_str());
  c->SetDirectory(nullptr);
  return c;
}

void drawHeader(double x = 0.18)
{
  TLatex t;
  t.SetNDC();
  t.SetTextFont(42);
  t.SetTextSize(0.040);
  t.DrawLatex(x, 0.88, strleg1.c_str());
  t.SetTextSize(0.034);
  t.DrawLatex(x, 0.83, strleg2.c_str());
  t.DrawLatex(x, 0.79, strleg3.c_str());
}

// Produce a PDF for a single period. Returns true on success.
// Also writes a compact per-sample integral table to stdout.
bool processPeriod(const std::string &period, const std::string &periodLabel)
{
  std::cout << "\n=====================================================" << std::endl;
  std::cout << "Processing period: " << (period.empty() ? "all-range" : period) << std::endl;
  std::cout << "=====================================================" << std::endl;

  // Open combined file
  std::string combinedPath = makePath("", period);
  TFile *fComb = TFile::Open(combinedPath.c_str(), "READ");
  if (!fComb || fComb->IsZombie()) {
    std::cerr << "[processPeriod] cannot open combined file: " << combinedPath << std::endl;
    return false;
  }

  // Per-sample files + histograms
  std::vector<TFile *>   files(kSamples.size(), nullptr);
  std::vector<TH1D *>    hSamples(kSamples.size(), nullptr);
  bool anyMissing = false;
  for (size_t i = 0; i < kSamples.size(); ++i) {
    std::string path = makePath(kSamples[i].name, period);
    files[i] = TFile::Open(path.c_str(), "READ");
    if (!files[i] || files[i]->IsZombie()) {
      std::cerr << "[processPeriod] cannot open " << path << std::endl;
      anyMissing = true;
      continue;
    }
    std::string clone = "h_" + kSamples[i].name + "_" + (period.empty() ? "all" : period);
    hSamples[i] = fetchHist(files[i], "h_truth_pT_0", clone);
    if (!hSamples[i]) anyMissing = true;
  }
  if (anyMissing) {
    std::cerr << "[processPeriod] missing per-sample inputs for period '"
              << period << "' — skipping." << std::endl;
    for (auto f : files) if (f) f->Close();
    fComb->Close();
    return false;
  }

  TH1D *hComb = fetchHist(fComb, "h_truth_pT_0",
                          std::string("h_combined_") + (period.empty() ? "all" : period));
  if (!hComb) {
    std::cerr << "[processPeriod] missing combined h_truth_pT_0" << std::endl;
    fComb->Close();
    return false;
  }

  const int    nb   = hComb->GetNbinsX();
  const double xmin = hComb->GetXaxis()->GetXmin();
  const double xmax = hComb->GetXaxis()->GetXmax();
  std::cout << "Combined binning: " << nb << " bins, [" << xmin << ", " << xmax << "]" << std::endl;

  // Build bin-by-bin sum of per-sample contents for sanity check
  TH1D *hSum = (TH1D *)hComb->Clone(("h_sum_of_samples_" + (period.empty() ? std::string("all") : period)).c_str());
  hSum->Reset();
  hSum->SetDirectory(nullptr);
  for (int b = 1; b <= nb; ++b) {
    double v  = 0.0;
    double e2 = 0.0;
    for (size_t i = 0; i < kSamples.size(); ++i) {
      v  += hSamples[i]->GetBinContent(b);
      e2 += hSamples[i]->GetBinError(b) * hSamples[i]->GetBinError(b);
    }
    hSum->SetBinContent(b, v);
    hSum->SetBinError(b, std::sqrt(e2));
  }

  // Print per-sample integral table
  std::cout << "\nIntegrals over full range (" << xmin << "-" << xmax << " GeV):" << std::endl;
  printf("  %-18s  %14s  %14s\n", "sample", "integral", "sum");
  double grandSum = 0.0;
  for (size_t i = 0; i < kSamples.size(); ++i) {
    double I = hSamples[i]->Integral(0, nb + 1);
    grandSum += I;
    printf("  %-18s  %14.6g\n", kSamples[i].name.c_str(), I);
  }
  double ICombined = hComb->Integral(0, nb + 1);
  double ISum      = hSum->Integral(0, nb + 1);
  printf("  %-18s  %14.6g\n", "SUM of per-sample", grandSum);
  printf("  %-18s  %14.6g\n", "COMBINED (hadd)",   ICombined);
  printf("  %-18s  %14.6g\n", "ratio sum/comb",    ICombined > 0 ? ISum / ICombined : -1);

  // Per-bin sum/combined ratio diagnostic
  std::cout << "\nBin-by-bin (sum_of_samples / combined):" << std::endl;
  double maxDev = 0.0; int maxDevBin = -1;
  for (int b = 1; b <= nb; ++b) {
    double c = hComb->GetBinContent(b);
    double s = hSum->GetBinContent(b);
    double r = c != 0 ? s / c : 0.0;
    if (std::fabs(r - 1.0) > maxDev) { maxDev = std::fabs(r - 1.0); maxDevBin = b; }
    printf("  bin %2d [%5.1f-%5.1f]  comb=%11.4g  sum=%11.4g  ratio=%.6f\n",
           b, hComb->GetBinLowEdge(b), hComb->GetBinLowEdge(b + 1), c, s, r);
  }
  if (maxDevBin > 0) {
    printf("Max |ratio-1| = %.2e at bin %d [%.1f, %.1f]\n",
           maxDev, maxDevBin, hComb->GetBinLowEdge(maxDevBin),
           hComb->GetBinLowEdge(maxDevBin + 1));
  }

  // =====================================================================
  // (A) Overlay canvas: log-y per-sample + combined
  // =====================================================================
  {
    TCanvas c("c_overlay", "overlay", 800, 700);
    c.SetLogy();
    c.SetRightMargin(0.05);
    c.SetLeftMargin(0.15);
    c.SetTopMargin(0.05);
    c.SetBottomMargin(0.13);

    // Prepare samples for drawing: use as-is (weights already baked in)
    for (size_t i = 0; i < kSamples.size(); ++i) {
      hSamples[i]->SetLineColor(kSamples[i].color);
      hSamples[i]->SetLineStyle(kSamples[i].style);
      hSamples[i]->SetLineWidth(2);
      hSamples[i]->SetMarkerColor(kSamples[i].color);
      hSamples[i]->SetMarkerStyle(20 + (i % 6));
      hSamples[i]->SetMarkerSize(0.6);
    }
    hComb->SetLineColor(kBlack);
    hComb->SetLineWidth(3);
    hComb->SetMarkerStyle(0);
    hComb->SetMarkerSize(0);

    // Y-axis range
    double ymax = hComb->GetMaximum();
    for (auto *h : hSamples) if (h->GetMaximum() > ymax) ymax = h->GetMaximum();
    double ymin = 1.0;
    for (int b = 1; b <= nb; ++b) {
      double c2 = hComb->GetBinContent(b);
      if (c2 > 0 && c2 < ymin) ymin = c2;
    }
    if (ymin <= 0) ymin = 1.0;

    hComb->SetTitle("");
    hComb->GetYaxis()->SetTitle("weighted counts / bin");
    hComb->GetXaxis()->SetTitle("truth #it{p}_{T}^{#gamma} [GeV]");
    hComb->GetYaxis()->SetTitleOffset(1.5);
    hComb->GetXaxis()->SetTitleOffset(1.1);
    hComb->GetYaxis()->SetRangeUser(ymin * 0.3, ymax * 20);

    hComb->Draw("HIST");
    for (size_t i = 0; i < kSamples.size(); ++i) {
      hSamples[i]->Draw("HIST SAME");
    }
    hComb->Draw("HIST SAME"); // put combined on top

    TLegend leg(0.55, 0.55, 0.93, 0.88);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.028);
    leg.AddEntry(hComb, "combined (hadd)", "L");
    for (size_t i = 0; i < kSamples.size(); ++i)
      leg.AddEntry(hSamples[i], kSamples[i].label.c_str(), "L");
    leg.Draw();

    drawHeader();
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.030);
    t.DrawLatex(0.18, 0.74, ("period: " + periodLabel).c_str());

    std::string out = std::string(kFigDir) + "/signal_truth_pT_overlay_" +
                      (period.empty() ? "all" : period) + ".pdf";
    c.SaveAs(out.c_str());
    std::cout << "Wrote " << out << std::endl;
  }

  // =====================================================================
  // (B) Ratio canvas: per-sample/combined + sum/combined
  // =====================================================================
  {
    TCanvas c("c_ratio", "ratio", 800, 900);
    c.Divide(1, 2, 0, 0);

    // --- Top pad: individual per-sample fractions (linear, 0..1) ------
    c.cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.02);

    std::vector<TH1D *> hFrac(kSamples.size(), nullptr);
    for (size_t i = 0; i < kSamples.size(); ++i) {
      hFrac[i] = (TH1D *)hSamples[i]->Clone(Form("hFrac_%s_%s", kSamples[i].name.c_str(),
                                                 period.empty() ? "all" : period.c_str()));
      hFrac[i]->SetDirectory(nullptr);
      hFrac[i]->Divide(hComb);
      hFrac[i]->SetLineColor(kSamples[i].color);
      hFrac[i]->SetLineStyle(kSamples[i].style);
      hFrac[i]->SetLineWidth(2);
      hFrac[i]->SetMarkerColor(kSamples[i].color);
      hFrac[i]->SetMarkerStyle(20 + (i % 6));
      hFrac[i]->SetMarkerSize(0.6);
    }
    hFrac[0]->SetTitle("");
    hFrac[0]->GetYaxis()->SetTitle("per-sample / combined");
    hFrac[0]->GetYaxis()->SetRangeUser(0.0, 1.05);
    hFrac[0]->GetXaxis()->SetLabelSize(0);
    hFrac[0]->GetXaxis()->SetTitleSize(0);
    hFrac[0]->GetYaxis()->SetTitleOffset(1.3);
    hFrac[0]->GetYaxis()->SetTitleSize(0.05);
    hFrac[0]->GetYaxis()->SetLabelSize(0.04);
    hFrac[0]->Draw("HIST");
    for (size_t i = 1; i < kSamples.size(); ++i) hFrac[i]->Draw("HIST SAME");

    TLine one(xmin, 1.0, xmax, 1.0);
    one.SetLineStyle(7);
    one.SetLineColor(kGray + 1);
    one.Draw("SAME");

    // Legend under top-right to avoid overlap with the y=1 reference line.
    TLegend leg(0.60, 0.18, 0.94, 0.50);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.030);
    for (size_t i = 0; i < kSamples.size(); ++i)
      leg.AddEntry(hFrac[i], kSamples[i].label.c_str(), "L");
    leg.Draw();

    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.045);
    t.DrawLatex(0.18, 0.86, strleg1.c_str());
    t.SetTextSize(0.032);
    t.DrawLatex(0.18, 0.81, ("period: " + periodLabel).c_str());
    t.DrawLatex(0.18, 0.76, strleg3.c_str());

    // --- Bottom pad: sum / combined (sanity; should be 1.0) ----------
    c.cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.15);

    TH1D *hSumOverComb = (TH1D *)hSum->Clone(("hSumOverComb_" + (period.empty() ? std::string("all") : period)).c_str());
    hSumOverComb->SetDirectory(nullptr);
    hSumOverComb->Divide(hComb);
    hSumOverComb->SetLineColor(kBlack);
    hSumOverComb->SetMarkerStyle(20);
    hSumOverComb->SetMarkerColor(kBlack);
    hSumOverComb->SetMarkerSize(0.9);
    hSumOverComb->SetLineWidth(2);
    hSumOverComb->SetTitle("");
    hSumOverComb->GetYaxis()->SetTitle("#Sigma(samples) / combined");
    hSumOverComb->GetYaxis()->SetRangeUser(0.8, 1.2);
    hSumOverComb->GetXaxis()->SetTitle("truth #it{p}_{T}^{#gamma} [GeV]");
    hSumOverComb->GetYaxis()->SetTitleOffset(1.3);
    hSumOverComb->GetYaxis()->SetTitleSize(0.05);
    hSumOverComb->GetYaxis()->SetLabelSize(0.04);
    hSumOverComb->GetXaxis()->SetTitleSize(0.05);
    hSumOverComb->GetXaxis()->SetLabelSize(0.04);
    hSumOverComb->Draw("E1");
    TLine one2(xmin, 1.0, xmax, 1.0);
    one2.SetLineStyle(7);
    one2.SetLineColor(kGray + 1);
    one2.Draw("SAME");

    std::string out = std::string(kFigDir) + "/signal_ratio_per_sample_" +
                      (period.empty() ? "all" : period) + ".pdf";
    c.SaveAs(out.c_str());
    std::cout << "Wrote " << out << std::endl;
  }

  // =====================================================================
  // (C) Boundary zoom: 12-16 GeV (photon5->photon10) and 28-32 GeV (photon10->photon20)
  // =====================================================================
  {
    TCanvas c("c_zoom", "zoom", 1200, 600);
    // Small pad gap so each panel keeps its own log-y axis visible.
    c.Divide(2, 1, 0.001, 0.001);

    auto drawZoom = [&](int pad, double xlo, double xhi, const char *title) {
      c.cd(pad);
      // Use independent pads with full margins on both so the log-y axis
      // is visible on each panel (the previous squeezed-margin layout on
      // pad 2 made the axis unreadable).
      gPad->SetLogy(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.05);
      gPad->SetTopMargin(0.08);
      gPad->SetBottomMargin(0.14);

      TH1D *hCombZ = (TH1D *)hComb->Clone(Form("hCombZ%d", pad));
      hCombZ->SetDirectory(nullptr);
      hCombZ->SetTitle("");
      hCombZ->GetXaxis()->SetRangeUser(xlo, xhi);
      hCombZ->GetXaxis()->SetTitle("truth #it{p}_{T}^{#gamma} [GeV]");
      hCombZ->GetYaxis()->SetTitle("weighted counts / bin");
      hCombZ->GetYaxis()->SetTitleOffset(1.3);
      hCombZ->GetYaxis()->SetTitleSize(0.05);
      hCombZ->GetYaxis()->SetLabelSize(0.04);
      hCombZ->GetXaxis()->SetTitleSize(0.05);
      hCombZ->GetXaxis()->SetLabelSize(0.04);
      hCombZ->SetLineColor(kBlack);
      hCombZ->SetLineWidth(3);

      // y-range based on window
      int blo = hCombZ->GetXaxis()->FindBin(xlo + 1e-6);
      int bhi = hCombZ->GetXaxis()->FindBin(xhi - 1e-6);
      double ymax = 0, ymin = 1e30;
      for (int b = blo; b <= bhi; ++b) {
        double v = hCombZ->GetBinContent(b);
        if (v > ymax) ymax = v;
        for (size_t i = 0; i < kSamples.size(); ++i) {
          double vv = hSamples[i]->GetBinContent(b);
          if (vv > ymax) ymax = vv;
          if (vv > 0 && vv < ymin) ymin = vv;
        }
        if (v > 0 && v < ymin) ymin = v;
      }
      if (ymin <= 0) ymin = 1;
      hCombZ->GetYaxis()->SetRangeUser(ymin * 0.3, ymax * 10);

      hCombZ->Draw("HIST");
      for (size_t i = 0; i < kSamples.size(); ++i) hSamples[i]->Draw("HIST SAME");
      hCombZ->Draw("HIST SAME");

      TLatex t;
      t.SetNDC();
      t.SetTextFont(42);
      t.SetTextSize(0.040);
      t.DrawLatex(0.18, 0.93, title);
      t.SetTextSize(0.030);
      t.DrawLatex(0.18, 0.89, ("period: " + periodLabel).c_str());
      t.DrawLatex(0.18, 0.85, strleg3.c_str());

      if (pad == 1) {
        TLegend leg(0.55, 0.55, 0.97, 0.84);
        leg.SetFillStyle(0);
        leg.SetBorderSize(0);
        leg.SetTextFont(42);
        leg.SetTextSize(0.028);
        leg.AddEntry(hCombZ, "combined", "L");
        for (size_t i = 0; i < kSamples.size(); ++i)
          leg.AddEntry(hSamples[i], kSamples[i].label.c_str(), "L");
        leg.Draw();
      }
    };

    drawZoom(1, 12.0, 16.0, "photon5 #rightarrow photon10 boundary");
    drawZoom(2, 28.0, 36.0, "photon10 #rightarrow photon20 boundary");

    std::string out = std::string(kFigDir) + "/signal_boundary_zoom_" +
                      (period.empty() ? "all" : period) + ".pdf";
    c.SaveAs(out.c_str());
    std::cout << "Wrote " << out << std::endl;
  }

  // Cleanup
  for (auto f : files) if (f) f->Close();
  fComb->Close();
  return true;
}

} // namespace

// Cross-period check: does (0rad combined) + (1p5mrad combined) == (all-range combined)?
void checkAllRangeMerge()
{
  std::cout << "\n=====================================================" << std::endl;
  std::cout << "Cross-period merge check (0rad + 1p5mrad vs all-range combined)" << std::endl;
  std::cout << "=====================================================" << std::endl;

  TFile *f0 = TFile::Open((std::string(kResultsDir) + "/MC_efficiency_bdt_nom_0rad.root").c_str(), "READ");
  TFile *f1 = TFile::Open((std::string(kResultsDir) + "/MC_efficiency_bdt_nom_1p5mrad.root").c_str(), "READ");
  TFile *fA = TFile::Open((std::string(kResultsDir) + "/MC_efficiency_bdt_nom.root").c_str(), "READ");
  if (!f0 || f0->IsZombie() || !f1 || f1->IsZombie() || !fA || fA->IsZombie()) {
    std::cerr << "[checkAllRangeMerge] cannot open one of the combined files" << std::endl;
    return;
  }
  TH1D *h0 = (TH1D *)f0->Get("h_truth_pT_0");
  TH1D *h1 = (TH1D *)f1->Get("h_truth_pT_0");
  TH1D *hA = (TH1D *)fA->Get("h_truth_pT_0");
  if (!h0 || !h1 || !hA) { std::cerr << "missing hist" << std::endl; return; }

  double I0 = h0->Integral(), I1 = h1->Integral(), IA = hA->Integral();
  std::cout << "0rad combined integral:       " << I0 << std::endl;
  std::cout << "1p5mrad combined integral:    " << I1 << std::endl;
  std::cout << "SUM 0rad + 1p5mrad:           " << I0 + I1 << std::endl;
  std::cout << "all-range combined integral:  " << IA << std::endl;
  std::cout << "ratio (sum / all-range):      " << (IA != 0 ? (I0 + I1) / IA : 0) << std::endl;

  double maxDev = 0.0; int maxDevBin = -1;
  for (int b = 1; b <= h0->GetNbinsX(); ++b) {
    double s = h0->GetBinContent(b) + h1->GetBinContent(b);
    double a = hA->GetBinContent(b);
    if (a > 0) {
      double r = s / a;
      if (std::fabs(r - 1.0) > maxDev) { maxDev = std::fabs(r - 1.0); maxDevBin = b; }
    }
  }
  if (maxDevBin > 0) {
    printf("Max |sum/all-range - 1| = %.2e at bin %d [%.1f, %.1f]\n",
           maxDev, maxDevBin, h0->GetBinLowEdge(maxDevBin), h0->GetBinLowEdge(maxDevBin + 1));
  }

  f0->Close(); f1->Close(); fA->Close();
}

void plot_sample_combining_signal()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);

  gSystem->mkdir(kFigDir, kTRUE);

  // Period labels show the lumi used in each config
  std::map<std::string, std::string> labels = {
    {"0rad",    "0 mrad (32.66 pb^{-1})"},
    {"1p5mrad", "1.5 mrad (16.27 pb^{-1})"},
    {"",        "all-range (48.93 pb^{-1})"},
  };

  bool ok0 = processPeriod("0rad",    labels["0rad"]);
  bool ok1 = processPeriod("1p5mrad", labels["1p5mrad"]);
  // Skip all-range per-sample plotting: bare per-sample files are stale copies of
  // 1p5mrad (Apr 20, before the Apr 21 re-run). The all-range combined file is
  // built by merge_periods.sh from the two period combined files; instead we do
  // a cross-period merge sanity check.
  checkAllRangeMerge();

  std::cout << "\n=== SUMMARY ===" << std::endl;
  std::cout << "  0rad:      " << (ok0 ? "OK" : "FAILED") << std::endl;
  std::cout << "  1p5mrad:   " << (ok1 ? "OK" : "FAILED") << std::endl;
}
