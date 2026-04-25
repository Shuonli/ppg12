// plot_sample_combining_background.C
//
// Validate that the five background samples (jet8/12/20/30/40, each with
// _nom + _double DI partners = 10 per-sample files) combine correctly into
// the jet-inclusive spectrum `MC_efficiency_jet_bdt_nom_{period}.root`.
//
// Diagnostics produced per period in {0rad, 1p5mrad, all-range (bare)}:
//   A. Per-sample leading truth-jet pT overlay + combined (log-y).
//   B. Ratio sum(10 per-sample) / combined (should be 1.0 everywhere).
//   C. Boundary zoom (13-15, 20-22, 31-33, 41-43 GeV) to expose any step.
//
// Output:
//   plotting/figures/sample_combining_check/
//     background_truth_jet_pT_overlay_{period}.pdf
//     background_ratio_per_sample_{period}.pdf
//     background_boundary_zoom_{period}.pdf

#include "plotcommon.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>

#include <array>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

const std::string kResultsDir =
    "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/";
const std::string kOutDir =
    "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures/"
    "sample_combining_check/";
const std::string kHistName = "h_max_truth_jet_pT";

// Per-jet-sample truth-pT windows (for reference in the plots).
struct JetWindow {
  std::string label;  // e.g. "jet8"
  double lo;          // truth pT lower edge [GeV]
  double hi;          // truth pT upper edge [GeV] (negative = open-ended)
  int color;
};

const std::vector<JetWindow> kJetSamples = {
    {"jet8",  9,  14, kBlue + 1},
    {"jet12", 14, 21, kGreen + 2},
    {"jet20", 21, 32, kOrange + 7},
    {"jet30", 32, 42, kRed + 1},
    {"jet40", 42, -1, kMagenta + 2},
};

struct PeriodCfg {
  std::string key;       // for filenames: 0rad, 1p5mrad, all
  std::string suffix;    // suffix on ROOT filenames: _0rad, _1p5mrad, ""
  std::string legend;    // legend label
};

// ---------------------------------------------------------------------------
// Utility: open file + fetch histogram as TH1D (clone; caller owns).
TH1D *FetchHist(const std::string &fn, const std::string &hname,
                const std::string &tag) {
  TFile *f = TFile::Open(fn.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[WARN] Cannot open " << fn << std::endl;
    return nullptr;
  }
  TH1 *h = dynamic_cast<TH1 *>(f->Get(hname.c_str()));
  if (!h) {
    std::cerr << "[WARN] Missing " << hname << " in " << fn << std::endl;
    f->Close();
    return nullptr;
  }
  auto *clone = dynamic_cast<TH1D *>(
      h->Clone(("_loaded_" + tag).c_str()));
  if (!clone) {
    // Promote TH1F -> TH1D so we have a uniform type.
    clone = new TH1D(("_loaded_" + tag).c_str(), h->GetTitle(),
                      h->GetNbinsX(), h->GetXaxis()->GetXmin(),
                      h->GetXaxis()->GetXmax());
    for (int i = 0; i <= h->GetNbinsX() + 1; ++i) {
      clone->SetBinContent(i, h->GetBinContent(i));
      clone->SetBinError(i, h->GetBinError(i));
    }
  }
  clone->SetDirectory(nullptr);
  f->Close();
  return clone;
}

// Rebin a 1000-bin, 0-100 GeV histogram to 0.5 GeV bins (200 bins).
TH1D *RebinHalfGeV(TH1D *h, const std::string &name) {
  if (!h) return nullptr;
  auto *r = dynamic_cast<TH1D *>(h->Rebin(5, name.c_str()));
  r->SetDirectory(nullptr);
  return r;
}

// Integral within [lo, hi] GeV (full-precision bins; caller uses pre-rebin).
double IntegralRange(TH1D *h, double lo, double hi) {
  if (!h) return 0.0;
  int b1 = h->FindBin(lo + 1e-6);
  int b2 = h->FindBin(hi - 1e-6);
  return h->Integral(b1, b2);
}

void PrintTableHeader(std::ostream &os) {
  os << "\n=================================================================\n";
  os << std::left << std::setw(18) << "sample" << std::setw(18)
     << "integral" << std::setw(14) << "entries" << std::setw(14)
     << "in-window" << "\n";
  os << "-----------------------------------------------------------------\n";
}

void PrintTableRow(std::ostream &os, const std::string &name, double integ,
                    double entries, double in_window) {
  std::ostringstream ss1, ss2, ss3;
  ss1 << std::scientific << std::setprecision(4) << integ;
  ss2 << std::scientific << std::setprecision(4) << entries;
  ss3 << std::scientific << std::setprecision(4) << in_window;
  os << std::left << std::setw(18) << name << std::setw(18) << ss1.str()
     << std::setw(14) << ss2.str() << std::setw(14) << ss3.str() << "\n";
}

// ---------------------------------------------------------------------------
// (A) Per-sample overlay + combined, log-y.
void MakeOverlayPlot(const PeriodCfg &period,
                      const std::vector<TH1D *> &per_sample,
                      TH1D *combined) {
  const double x_lo = 5.0;
  const double x_hi = 60.0;

  TCanvas c("c_overlay", "c_overlay", 900, 700);
  c.SetLogy(true);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.06);
  c.SetBottomMargin(0.13);

  TH1F *frame = new TH1F(
      "frame_overlay", ";truth leading jet #it{p}_{T} [GeV];Weighted counts",
      100, x_lo, x_hi);
  double ymax = combined ? combined->GetMaximum() * 5.0 : 1.0;
  // Scan ALL per-sample histograms (including _double, which at 1.5 mrad are
  // ~0.086 x _nom) to anchor ymin low enough to show the DI tail.
  double ymin = ymax;
  for (auto *h : per_sample) {
    if (!h) continue;
    for (int b = h->FindBin(x_lo); b <= h->FindBin(x_hi); ++b) {
      double v = h->GetBinContent(b);
      if (v > 0 && v < ymin) ymin = v;
    }
  }
  // Pad the floor by another 0.3x to leave visual breathing room.
  ymin *= 0.3;
  if (ymin < ymax * 1e-10) ymin = ymax * 1e-10;
  if (ymin >= ymax) ymin = ymax * 1e-8;
  frame->GetYaxis()->SetRangeUser(ymin, ymax);
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitleOffset(1.4);
  frame->Draw("AXIS");

  TLegend *leg = new TLegend(0.55, 0.55, 0.93, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.028);

  // Per-sample: solid for _nom, dashed for _double, same color per jet pT.
  const int ntypes = 2;
  const std::array<const char *, ntypes> suffixes = {"_nom", "_double"};
  const std::array<int, ntypes> line_styles = {1, 2};

  for (size_t i = 0; i < kJetSamples.size(); ++i) {
    for (int k = 0; k < ntypes; ++k) {
      size_t idx = 2 * i + k;
      if (idx >= per_sample.size() || !per_sample[idx]) continue;
      TH1D *h = per_sample[idx];
      h->SetLineColor(kJetSamples[i].color);
      h->SetMarkerColor(kJetSamples[i].color);
      h->SetLineStyle(line_styles[k]);
      h->SetLineWidth(2);
      h->Draw("HIST SAME");
      std::string lbl = kJetSamples[i].label + suffixes[k];
      leg->AddEntry(h, lbl.c_str(), "l");
    }
  }

  if (combined) {
    combined->SetLineColor(kBlack);
    combined->SetLineWidth(3);
    combined->SetLineStyle(1);
    combined->Draw("HIST SAME");
    leg->AddEntry(combined, "Combined jet", "l");
  }
  leg->Draw();

  // Header text: sPHENIX internal + period + strleg3 eta + strIncMC.
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.033);
  lat.DrawLatex(0.17, 0.88, strleg1.c_str());
  lat.DrawLatex(0.17, 0.84, strleg2.c_str());
  lat.DrawLatex(0.17, 0.80, ("Period: " + period.legend).c_str());
  lat.DrawLatex(0.17, 0.76, strleg3.c_str());
  lat.DrawLatex(0.17, 0.72, strIncMC.c_str());

  std::string out =
      kOutDir + "background_truth_jet_pT_overlay_" + period.key + ".pdf";
  c.SaveAs(out.c_str());
  std::cout << "  wrote " << out << std::endl;

  delete frame;
  delete leg;
}

// ---------------------------------------------------------------------------
// (B) Ratio: sum(10 per-sample) / combined.
void MakeRatioPlot(const PeriodCfg &period,
                    const std::vector<TH1D *> &per_sample, TH1D *combined) {
  if (!combined) return;
  TH1D *sum = dynamic_cast<TH1D *>(per_sample[0]->Clone("h_sum_all"));
  sum->SetDirectory(nullptr);
  sum->Reset();
  for (auto *h : per_sample) {
    if (h) sum->Add(h);
  }

  TH1D *ratio = dynamic_cast<TH1D *>(sum->Clone("h_ratio_sum_combined"));
  ratio->SetDirectory(nullptr);
  ratio->Divide(combined);

  // Report min / max / mean of the ratio in the sensitive range.
  const double x_lo = 5.0;
  const double x_hi = 55.0;
  int b_lo = ratio->FindBin(x_lo + 1e-6);
  int b_hi = ratio->FindBin(x_hi - 1e-6);
  double rmin = 1e30, rmax = -1e30, rsum = 0.0;
  int nfill = 0;
  for (int b = b_lo; b <= b_hi; ++b) {
    double v = ratio->GetBinContent(b);
    if (v == 0.0) continue;  // skip empty combined bins
    if (v < rmin) rmin = v;
    if (v > rmax) rmax = v;
    rsum += v;
    ++nfill;
  }
  double rmean = (nfill > 0) ? rsum / nfill : 0.0;
  std::cout << "  [" << period.key << "] sum/combined over [" << x_lo
            << ", " << x_hi << "] GeV  min=" << rmin << "  max=" << rmax
            << "  mean=" << rmean << "  nfilled=" << nfill << std::endl;

  TCanvas c("c_ratio", "c_ratio", 900, 500);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.08);
  c.SetBottomMargin(0.15);

  TH1F *frame = new TH1F(
      "frame_ratio",
      ";truth leading jet #it{p}_{T} [GeV];#sum_{samples} / Combined",
      100, x_lo, x_hi);
  frame->GetYaxis()->SetRangeUser(0.95, 1.05);
  frame->GetXaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetTitleOffset(1.4);
  frame->Draw("AXIS");

  ratio->SetLineColor(kBlack);
  ratio->SetMarkerColor(kBlack);
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerSize(0.6);
  ratio->SetLineWidth(2);
  ratio->Draw("HIST P SAME");

  TLine *l1 = new TLine(x_lo, 1.0, x_hi, 1.0);
  l1->SetLineColor(kRed);
  l1->SetLineStyle(2);
  l1->Draw();

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.17, 0.88, strleg1.c_str());
  lat.DrawLatex(0.17, 0.83, ("Period: " + period.legend).c_str());
  lat.DrawLatex(0.17, 0.78, strleg3.c_str());
  std::ostringstream stats;
  stats << "min=" << std::setprecision(5) << rmin
        << "  max=" << rmax << "  mean=" << rmean;
  lat.DrawLatex(0.17, 0.73, stats.str().c_str());

  std::string out =
      kOutDir + "background_ratio_per_sample_" + period.key + ".pdf";
  c.SaveAs(out.c_str());
  std::cout << "  wrote " << out << std::endl;

  delete frame;
  delete l1;
  delete sum;
  delete ratio;
}

// ---------------------------------------------------------------------------
// (C) Boundary zoom: 4 pads around sample boundaries.
void MakeBoundaryZoom(const PeriodCfg &period,
                       const std::vector<TH1D *> &per_sample,
                       TH1D *combined) {
  struct Zoom {
    std::string title;
    double x_lo;
    double x_hi;
  };
  const std::vector<Zoom> zooms = {
      {"jet8 #rightarrow jet12 (14 GeV)",  12.0, 16.0},
      {"jet12 #rightarrow jet20 (21 GeV)", 19.0, 23.0},
      {"jet20 #rightarrow jet30 (32 GeV)", 30.0, 34.0},
      {"jet30 #rightarrow jet40 (42 GeV)", 40.0, 44.0},
  };

  TCanvas c("c_bzoom", "c_bzoom", 1200, 900);
  c.Divide(2, 2);

  for (size_t iz = 0; iz < zooms.size(); ++iz) {
    const Zoom &z = zooms[iz];
    c.cd(iz + 1);
    gPad->SetLogy(true);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.14);

    std::string fname = "frame_bz_" + std::to_string(iz);
    TH1F *frame = new TH1F(fname.c_str(),
                            (z.title +
                             ";truth leading jet #it{p}_{T} [GeV];Weighted counts")
                                .c_str(),
                            100, z.x_lo, z.x_hi);

    // Range: scan combined AND all per-sample histograms (_double are
    // ~0.086 x _nom at 1.5 mrad, so anchor ymin on them to stay visible).
    double ymax = 0.0, ymin_nz = 1e30;
    int b1 = combined ? combined->FindBin(z.x_lo + 1e-6) : 1;
    int b2 = combined ? combined->FindBin(z.x_hi - 1e-6) : 1;
    if (combined) {
      for (int b = b1; b <= b2; ++b) {
        double v = combined->GetBinContent(b);
        if (v > ymax) ymax = v;
        if (v > 0 && v < ymin_nz) ymin_nz = v;
      }
    }
    for (auto *h : per_sample) {
      if (!h) continue;
      for (int b = b1; b <= b2; ++b) {
        double v = h->GetBinContent(b);
        if (v > ymax) ymax = v;
        if (v > 0 && v < ymin_nz) ymin_nz = v;
      }
    }
    if (ymax <= 0.0) ymax = 1.0;
    if (ymin_nz > ymax) ymin_nz = ymax * 1e-4;
    frame->GetYaxis()->SetRangeUser(ymin_nz * 0.3, ymax * 5.0);
    frame->GetXaxis()->SetTitleOffset(1.2);
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->Draw("AXIS");

    const std::array<const char *, 2> suffixes = {"_nom", "_double"};
    const std::array<int, 2> line_styles = {1, 2};
    for (size_t i = 0; i < kJetSamples.size(); ++i) {
      for (int k = 0; k < 2; ++k) {
        size_t idx = 2 * i + k;
        if (idx >= per_sample.size() || !per_sample[idx]) continue;
        TH1D *h = per_sample[idx];
        h->SetLineColor(kJetSamples[i].color);
        h->SetLineStyle(line_styles[k]);
        h->SetLineWidth(2);
        h->Draw("HIST SAME");
      }
    }
    if (combined) {
      combined->SetLineColor(kBlack);
      combined->SetLineWidth(3);
      combined->Draw("HIST SAME");
    }

    // Mark the boundary with a vertical line.
    double boundary = (z.x_lo + z.x_hi) / 2.0;
    TLine *vline = new TLine(boundary, ymin_nz * 0.3, boundary, ymax * 5.0);
    vline->SetLineColor(kGray + 2);
    vline->SetLineStyle(3);
    vline->SetLineWidth(2);
    vline->Draw();

    // Draw the legend only in pad 4 (the 42 GeV zoom has the sparsest data,
    // so the legend does not overlap any curves there).
    if (iz == 3) {
      TLegend *leg = new TLegend(0.22, 0.40, 0.95, 0.82);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.032);
      leg->SetNColumns(3);
      for (size_t i = 0; i < kJetSamples.size(); ++i) {
        for (int k = 0; k < 2; ++k) {
          size_t idx = 2 * i + k;
          if (idx >= per_sample.size() || !per_sample[idx]) continue;
          std::string lbl = kJetSamples[i].label + suffixes[k];
          leg->AddEntry(per_sample[idx], lbl.c_str(), "l");
        }
      }
      if (combined) leg->AddEntry(combined, "Combined jet", "l");
      leg->Draw();
    }
  }

  // Global title.
  c.cd(0);
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.022);
  lat.DrawLatex(0.02, 0.97, (strleg1 + "   " + strleg2 + "   Period: " +
                              period.legend + "   " + strIncMC)
                                 .c_str());

  std::string out =
      kOutDir + "background_boundary_zoom_" + period.key + ".pdf";
  c.SaveAs(out.c_str());
  std::cout << "  wrote " << out << std::endl;
}

// ---------------------------------------------------------------------------
void ProcessPeriod(const PeriodCfg &period) {
  std::cout << "\n================================================\n";
  std::cout << "Period: " << period.key << "  (suffix='" << period.suffix
            << "')\n";
  std::cout << "================================================\n";

  std::vector<TH1D *> per_sample_raw;  // 1000-bin originals for integrals
  std::vector<TH1D *> per_sample_reb;  // 0.5-GeV rebin for plotting
  std::vector<std::string> per_sample_names;

  for (const auto &w : kJetSamples) {
    for (const char *suf : {"_nom", "_double"}) {
      std::string sname = w.label + suf;
      std::string fn = kResultsDir + "MC_efficiency_" + sname + "_bdt_nom" +
                       period.suffix + ".root";
      TH1D *h = FetchHist(fn, kHistName, sname + "_" + period.key);
      if (!h) {
        std::cerr << "[ERROR] missing " << fn << " — skipping sample\n";
      }
      per_sample_raw.push_back(h);
      per_sample_reb.push_back(
          h ? RebinHalfGeV(h, (sname + "_reb_" + period.key).c_str())
            : nullptr);
      per_sample_names.push_back(sname);
    }
  }

  std::string combfn =
      kResultsDir + "MC_efficiency_jet_bdt_nom" + period.suffix + ".root";
  TH1D *combined_raw =
      FetchHist(combfn, kHistName, "combined_" + period.key);
  TH1D *combined_reb =
      combined_raw ? RebinHalfGeV(combined_raw, ("combined_reb_" + period.key).c_str())
                   : nullptr;
  if (!combined_raw) {
    std::cerr << "[ERROR] missing combined file " << combfn
              << " — aborting period\n";
    return;
  }

  // -------- Table --------
  PrintTableHeader(std::cout);
  double sum_all = 0.0;
  for (size_t i = 0; i < per_sample_raw.size(); ++i) {
    TH1D *h = per_sample_raw[i];
    if (!h) {
      PrintTableRow(std::cout, per_sample_names[i], 0.0, 0.0, 0.0);
      continue;
    }
    double integ = h->Integral(0, h->GetNbinsX() + 1);
    double entries = h->GetEntries();
    // in-window integral for this sample's nominal truth pT window
    size_t wi = i / 2;
    double lo = kJetSamples[wi].lo;
    double hi = kJetSamples[wi].hi > 0 ? kJetSamples[wi].hi : 100.0;
    double in_w = IntegralRange(h, lo, hi);
    PrintTableRow(std::cout, per_sample_names[i], integ, entries, in_w);
    sum_all += integ;
  }
  double comb_integ = combined_raw->Integral(0, combined_raw->GetNbinsX() + 1);
  double comb_entries = combined_raw->GetEntries();
  std::cout << "-----------------------------------------------------------------\n";
  PrintTableRow(std::cout, "SUM(10)", sum_all, 0.0, 0.0);
  PrintTableRow(std::cout, "COMBINED", comb_integ, comb_entries, 0.0);
  std::cout << "\n  sum / combined = "
            << std::setprecision(8) << (sum_all / comb_integ) << std::endl;

  // -------- Plots --------
  MakeOverlayPlot(period, per_sample_reb, combined_reb);
  MakeRatioPlot(period, per_sample_reb, combined_reb);
  MakeBoundaryZoom(period, per_sample_reb, combined_reb);
}

}  // anonymous namespace

void plot_sample_combining_background() {
  init_plot();
  gStyle->SetOptStat(0);
  gSystem->mkdir(kOutDir.c_str(), kTRUE);

  std::vector<PeriodCfg> periods = {
      {"0rad", "_0rad", "0 mrad, L = 32.66 pb^{-1}"},
      {"1p5mrad", "_1p5mrad", "1.5 mrad, L = 16.27 pb^{-1}"},
      // All-range per-sample combining is not run. The bare
      // MC_efficiency_jet{8..40}_{nom,double}_bdt_nom.root files on disk
      // (Apr 20) are stale copies of the 1.5 mrad content; they are not used
      // by merge_periods.sh, which combines at the *combined* file level.
      // Cross-period additivity is verified separately in the signal macro
      // (checkAllRangeMerge) and at the hadd/TFileMerger level.
  };
  for (const auto &p : periods) ProcessPeriod(p);

  std::cout << "\nDone. Output in " << kOutDir << std::endl;
}
