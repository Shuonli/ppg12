// x_T-scaling compilation -- REPRODUCES Bock HP2018 (arXiv:1901.10950) Jacobian bug
// for ATLAS data. Multiplies each ATLAS point's Y by (pT_bin / pT_fixed), with
// pT_fixed pinned to the per-table value inferred from EPS coordinate decode:
//   ATLAS 13 TeV: pT_fixed = 132.6 GeV
//   ATLAS  7 TeV: pT_fixed = 135.3 GeV
//   ATLAS  8 TeV: pT_fixed = 194.0 GeV  (different data source possibly)
// PHENIX 200 GeV is left alone (Bock plotted it correctly because PHENIX
// publishes invariant E d^3sigma/dp^3 directly, no Jacobian needed).
//
// Output: xtscaling_compilation_bockbug.pdf
//
// Run from the plotting/ directory:
//   cd plotting && root -l -b -q xtscaling/plot_xtscaling_bockbug.C

#include "plotcommon.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

namespace {
const char *XDIR = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/xtscaling";

struct Pt { double x, y, ylo, yhi; };
struct DS {
  std::string id, label;
  int color, marker;
  double msize;
  bool ours;
  double sqrt_s;   // [GeV] for back-computing pT = x_T * sqrt_s / 2
  double pt_fixed; // [GeV], 0 = no bug applied; otherwise Y *= pT_bin / pt_fixed
  std::vector<Pt> pts;
};

std::vector<Pt> readCSV(const std::string &path) {
  std::vector<Pt> v;
  std::ifstream in(path.c_str());
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::stringstream ss(line);
    std::string c; std::vector<double> f;
    while (std::getline(ss, c, ',')) f.push_back(atof(c.c_str()));
    if (f.size() >= 4) v.push_back({f[0], f[1], f[2], f[3]});
  }
  return v;
}
}  // namespace

void plot_xtscaling_bockbug() {
  init_plot();
  gStyle->SetOptStat(0);

  // style table; only ids present on disk are drawn. order = legend order after PPG12.
  // marker styles match Bock HP2018 Fig.1 (1901.10950) as closely as practical.
  std::vector<DS> all = {
    {"ALICE_13TeV",           "ALICE (13 TeV)",   kPink+8,    20, 1.0, false, 13000,    0.0, {}},
    {"ALICE_7TeV",            "ALICE (7 TeV)",    kViolet-6,  21, 1.0, false,  7000,    0.0, {}},
    {"ATLAS_13TeV",           "ATLAS (13 TeV)",   kViolet+1,  21, 1.0, false, 13000,  132.6, {}},  // BOCK BUG
    {"ATLAS_8TeV",            "ATLAS (8 TeV)",    kGreen+2,   27, 1.1, false,  8000,  194.0, {}},  // BOCK BUG
    {"ATLAS_7TeV",            "ATLAS (7 TeV)",    kBlue,      24, 1.0, false,  7000,  135.3, {}},  // BOCK BUG
    {"CMS_7TeV",              "CMS (7 TeV)",      kAzure+1,   25, 1.0, false,  7000,    0.0, {}},
    {"CDF_1800",              "CDF (1.8 TeV)",    kTeal+3,    34, 1.1, false,  1800,    0.0, {}},
    {"D0_1800",               "D0 (1.8 TeV)",     kTeal-1,    33, 1.1, false,  1800,    0.0, {}},
    {"UA1_630",               "UA1 (630 GeV)",    kOrange+1,  29, 1.2, false,   630,    0.0, {}},
    {"UA2_630",               "UA2 (630 GeV)",    kOrange+7,  25, 1.0, false,   630,    0.0, {}},
    {"PHENIX_200GeV",         "PHENIX (200 GeV)", kViolet-1,  20, 1.1, false,   200,    0.0, {}},
    {"PHENIX_510GeV_isolated","PHENIX iso. (510 GeV)", kAzure+2, 22, 1.0, false, 510,   0.0, {}},
    {"R110_63",               "R110 (63 GeV)",    kGray+2,    5,  1.2, false,    63,    0.0, {}},
    {"E706_38p8GeV",          "E706 (38.8 GeV)",  kGray+2,    30, 1.1, false,  38.8,    0.0, {}},
    {"E706_31p6GeV",          "E706 (31.6 GeV)",  kGray+2,    25, 1.0, false,  31.6,    0.0, {}},
    {"NA24_23p8",             "NA24 (23.8 GeV)",  kGray+2,    20, 0.9, false,  23.8,    0.0, {}},
    {"UA6_24p3",              "UA6 (24.3 GeV)",   kGray+2,    28, 1.0, false,  24.3,    0.0, {}},
    {"WA70_23",               "WA70 (23 GeV)",    kGray+2,    21, 0.9, false,    23,    0.0, {}},
    {"E704_19p4",             "E704 (19.4 GeV)",  kBlack,     24, 0.9, false,  19.4,    0.0, {}},
    {"PPG12_200GeV",          "sPHENIX (200 GeV)",kRed,       29, 2.4, true,    200,    0.0, {}},
  };

  std::vector<DS> ds;
  double ymin = 1e300, ymax = -1e300, xmin = 1e300, xmax = -1e300;
  for (auto d : all) {
    d.pts = readCSV(std::string(XDIR) + "/converted/" + d.id + ".csv");
    if (d.pts.empty()) { printf("  (absent) %s\n", d.id.c_str()); continue; }
    // Apply Bock's apparent Jacobian bug if pt_fixed > 0
    if (d.pt_fixed > 0 && d.sqrt_s > 0) {
      for (auto &p : d.pts) {
        double pt = p.x * d.sqrt_s / 2.0;
        double scale = pt / d.pt_fixed;
        p.y *= scale; p.ylo *= scale; p.yhi *= scale;
      }
      printf("  applied BOCK BUG to %-18s pt_fixed=%.1f GeV\n", d.id.c_str(), d.pt_fixed);
    }
    for (auto &p : d.pts) {
      if (p.y > 0) { ymin = std::min(ymin, p.ylo > 0 ? p.ylo : p.y); ymax = std::max(ymax, p.yhi); }
      xmin = std::min(xmin, p.x); xmax = std::max(xmax, p.x);
    }
    ds.push_back(d);
    printf("  loaded %-24s n=%zu\n", d.id.c_str(), d.pts.size());
  }
  if (ds.empty()) { printf("no datasets found\n"); return; }

  // Match Bock HP2018 (arXiv 1901.10950) axis range exactly for side-by-side check
  double xlo = 1e-4, xhi = 1.5;
  double ylo = 1e4,  yhi = 1e29;

  TCanvas *c = new TCanvas("c_xt", "", 820, 840);
  c->SetLogx(); c->SetLogy(); c->SetTicks(1, 1);
  c->SetLeftMargin(0.165); c->SetRightMargin(0.04);
  c->SetTopMargin(0.05); c->SetBottomMargin(0.115);

  TH2F *fr = new TH2F("fr_xt", "", 100, xlo, xhi, 100, ylo, yhi);
  fr->GetXaxis()->SetTitle("#it{x}_{T} = 2#it{p}_{T} / #sqrt{#it{s}}");
  fr->GetYaxis()->SetTitle("(#sqrt{#it{s}}/GeV)^{#it{n}} #it{E} d^{3}#it{#sigma}/d#it{p}^{3}  [pb GeV^{-2} #it{c}^{3}]");
  fr->GetXaxis()->SetTitleOffset(1.15);
  fr->GetYaxis()->SetTitleOffset(1.7);
  // decade-only x labels (avoid 6e-3 colliding with 1e-2 at the low end)
  fr->GetXaxis()->SetNdivisions(510);
  fr->Draw();

  // draw non-PPG12 first, PPG12 last; keep (DS,graph) in draw order for the legend
  std::vector<std::pair<const DS *, TGraphAsymmErrors *>> drawn;
  for (int pass = 0; pass < 2; ++pass) {
    for (auto &d : ds) {
      if (d.ours != (pass == 1)) continue;
      auto *g = new TGraphAsymmErrors();
      for (size_t i = 0; i < d.pts.size(); ++i) {
        const Pt &p = d.pts[i];
        g->SetPoint(i, p.x, p.y);
        g->SetPointError(i, 0, 0, std::max(0.0, p.y - p.ylo), std::max(0.0, p.yhi - p.y));
      }
      g->SetMarkerColor(d.color); g->SetLineColor(d.color);
      g->SetMarkerStyle(d.marker); g->SetMarkerSize(d.msize);
      g->SetLineWidth(d.ours ? 2 : 1);
      g->Draw("P SAME");
      drawn.push_back({&d, g});
    }
  }

  TLegend *leg = new TLegend(0.175, 0.115, 0.75, 0.43);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.026);
  leg->SetNColumns(2); leg->SetColumnSeparation(0.02);
  for (auto &pr : drawn) if (pr.first->ours)
    leg->AddEntry(pr.second, ("#bf{" + pr.first->label + "}").c_str(), "p");
  for (auto &pr : drawn) if (!pr.first->ours)
    leg->AddEntry(pr.second, pr.first->label.c_str(), "p");
  leg->Draw();

  TLatex tl; tl.SetNDC(); tl.SetTextFont(42);
  tl.SetTextSize(0.038); tl.DrawLatex(0.21, 0.945, "Bock-style Jacobian (constant pT) reproduction");
  tl.SetTextSize(0.041); tl.DrawLatex(0.585, 0.885, "#bf{#it{sPHENIX}} Internal");
  tl.SetTextSize(0.032);
  tl.DrawLatex(0.585, 0.840, "#it{p}+#it{p}(#bar{#it{p}}) #rightarrow #it{#gamma} + X");
  tl.DrawLatex(0.585, 0.804, "#it{y} #approx 0,  #it{n} = 4.5");

  gSystem->mkdir("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures", true);
  c->SaveAs("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures/xtscaling_compilation_bockbug.pdf");
  c->SaveAs(Form("%s/xtscaling_compilation_bockbug.pdf", XDIR));
  printf("saved xtscaling_compilation_bockbug.pdf  (x[%.3f,%.3f] y[%.1e,%.1e], %zu datasets)\n",
         xlo, xhi, ylo, yhi, ds.size());
}
