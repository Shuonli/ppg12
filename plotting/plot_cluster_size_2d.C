// plot_cluster_size_2d.C
//
// Companion to plot_cluster_size_cuts.C. Produces two new figure families for
// the cluster_size_short report, using the TH2F `h2_{metric}_{nkey}{suffix}`
// histograms added to ClusterSizeStudy.C:
//
//   (A) 2D heat-maps of cluster-size metric vs cluster E_T
//       per metric x cut x period, 3-panel row (data | photon pool | jet pool)
//       filename: cluster_size_2d_{metric}_{variant}_split_{period}.pdf
//
//   (B) 1D cluster-size distributions binned by cluster E_T
//       per metric x cut x period, 4x3 grid of the 12 reco pT bins from
//       plotcommon.h (ptRanges[] = 8,10,...,36), each panel overlays
//       data / photon pool / jet pool, normalized to unit area.
//       filename: cluster_size_1d_byet_{metric}_{variant}_split_{period}.pdf
//
// Sample pooling and cut-variant suffixes mirror plot_cluster_size_cuts.C.

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "plotcommon.h"  // ptRanges[13], NptBins=12, SetsPhenixStyle()

static const char* kResultsDir =
    "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results";
static const char* kFigDir =
    "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/reports/figures/cluster_size";

namespace {

struct SampleFile {
  std::string filetype;
  std::unique_ptr<TFile> file;
};

std::unique_ptr<SampleFile> LoadSample(const std::string& period,
                                       const std::string& filetype) {
  auto sf = std::make_unique<SampleFile>();
  sf->filetype = filetype;
  const std::string path =
      std::string(kResultsDir) + "/cluster_size_cuts_" + period +
      "_" + filetype + ".root";
  sf->file.reset(TFile::Open(path.c_str(), "READ"));
  if (!sf->file || sf->file->IsZombie()) {
    std::cerr << "Cannot open " << path << std::endl;
    sf->file.reset();
  }
  return sf;
}

// Pool TH2F across a list of samples by TH2::Add, returning heap ownership.
TH2F* PoolH2(const std::vector<const SampleFile*>& samples,
             const std::string& hname,
             const std::string& newname)
{
  TH2F* out = nullptr;
  for (const SampleFile* sf : samples) {
    if (!sf || !sf->file) continue;
    auto* src = dynamic_cast<TH2F*>(sf->file->Get(hname.c_str()));
    if (!src) continue;
    if (!out) {
      out = (TH2F*)src->Clone(newname.c_str());
      out->SetDirectory(nullptr);
    } else {
      out->Add(src, 1.0);
    }
  }
  return out;
}

TLatex* MakeLabel(const std::string& txt, double x, double y,
                  double size = 0.045) {
  auto* l = new TLatex(x, y, txt.c_str());
  l->SetNDC();
  l->SetTextFont(42);
  l->SetTextSize(size);
  return l;
}

void sPhenixInternal(double x, double y, const std::string& period,
                     const std::string& luminosity, double size = 0.045) {
  auto* t1 = MakeLabel("#bf{#it{sPHENIX}} Internal", x, y, size);
  t1->Draw();
  const std::string period_txt =
      "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, " + period +
      ", " + luminosity + " pb^{-1}";
  auto* t2 = MakeLabel(period_txt, x, y - 0.055, size * 0.78);
  t2->Draw();
}

void EnsureDir(const char* path) { gSystem->mkdir(path, true); }

struct PeriodParams {
  std::string period_key;    // "1p5mrad" / "0mrad"
  std::string period_label;  // "1.5 mrad" / "0 mrad"
  std::string lumi;          // "16.86" / "32.66"
};

PeriodParams MakePeriodParams(const std::string& pk) {
  PeriodParams p;
  p.period_key = pk;
  if (pk == "1p5mrad") {
    p.period_label = "1.5 mrad";
    p.lumi         = "16.86";
  } else {
    p.period_label = "0 mrad";
    p.lumi         = "32.66";
  }
  return p;
}

struct MetricSpec {
  std::string key;      // "n_owned" / "width_eta" / "width_phi"
  std::string ytitle;   // axis title for 2D and 1D
  std::string shortname;
  double y_min;         // for 2D display range
  double y_max;
  int rebin_y;          // 1 unless we want to coarsen
};

const std::vector<MetricSpec> kMetrics = {
    {"n_owned",   "n_{owned}",      "n_{owned}",     0.5,  25.5, 1},
    {"width_eta", "#Deltai_{#eta}", "#Deltai_{#eta}", 0.5,  7.5, 1},
    {"width_phi", "#Deltai_{#phi}", "#Deltai_{#phi}", 0.5,  7.5, 1},
};

struct VariantSpec {
  std::string suffix;   // "" / "_common" / "_tight"
  std::string key;      // "baseline" / "common" / "tight"
  std::string label;    // displayed
};

const std::vector<VariantSpec> kVariants = {
    {"",        "baseline", "baseline"},
    {"_common", "common",   "common cut"},
    {"_tight",  "tight",    "tight cut"},
};

// -----------------------------------------------------------------------
// Type A:  3-panel TH2 row (data | photon pool | jet pool)
// -----------------------------------------------------------------------
void Figure_2D_MetricCutPeriod(
    const std::map<std::string, std::unique_ptr<SampleFile>>& all,
    const MetricSpec& m, const VariantSpec& v, const PeriodParams& p,
    const std::string& nkey = "split")
{
  const std::string hname = "h2_" + m.key + "_" + nkey + v.suffix;

  auto get = [&](const std::string& k) -> const SampleFile* {
    auto it = all.find(k);
    return it == all.end() ? nullptr : it->second.get();
  };

  const SampleFile* data = get("data");
  std::vector<const SampleFile*> photon_pool = {
      get("photon5"), get("photon10"), get("photon20")};
  std::vector<const SampleFile*> jet_pool = {
      get("jet5"), get("jet8"), get("jet12"),
      get("jet20"), get("jet30"), get("jet40")};

  TH2F* h_data   = data ? dynamic_cast<TH2F*>(data->file->Get(hname.c_str())) : nullptr;
  if (h_data) h_data->SetDirectory(nullptr);
  TH2F* h_photon = PoolH2(photon_pool, hname,
      Form("ph_pool_%s_%s%s_%s",
           m.key.c_str(), nkey.c_str(), v.suffix.c_str(), p.period_key.c_str()));
  TH2F* h_jet    = PoolH2(jet_pool, hname,
      Form("jt_pool_%s_%s%s_%s",
           m.key.c_str(), nkey.c_str(), v.suffix.c_str(), p.period_key.c_str()));

  auto styleH2 = [&](TH2F* h, const std::string& title_tag) {
    if (!h) return;
    h->SetStats(0);
    h->SetTitle("");
    h->GetXaxis()->SetRangeUser(4, 36);
    h->GetYaxis()->SetRangeUser(m.y_min, m.y_max);
    h->GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
    h->GetYaxis()->SetTitle(m.ytitle.c_str());
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.042);
    h->GetYaxis()->SetLabelSize(0.042);
    h->GetYaxis()->SetTitleOffset(1.2);
    h->GetZaxis()->SetLabelSize(0.035);
  };
  styleH2(h_data,   "data");
  styleH2(h_photon, "photon");
  styleH2(h_jet,    "jet");

  auto* c = new TCanvas(
      Form("c2d_%s_%s_%s_%s", m.key.c_str(), v.key.c_str(),
           nkey.c_str(), p.period_key.c_str()),
      Form("cluster size 2D %s %s %s %s",
           m.key.c_str(), v.key.c_str(), nkey.c_str(), p.period_key.c_str()),
      1800, 600);
  c->Divide(3, 1, 0.003, 0.003);

  auto drawPanel = [&](int idx, TH2F* h, const std::string& tag) {
    c->cd(idx);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.17);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.14);
    if (!h || h->GetEntries() <= 0) {
      auto* l = MakeLabel("no entries", 0.30, 0.50, 0.060);
      l->Draw();
      return;
    }
    h->Draw("COLZ");
    sPhenixInternal(0.17, 0.92, p.period_label, p.lumi, 0.045);
    auto* t = MakeLabel(Form("%s (%s, %s)", tag.c_str(), v.label.c_str(),
                             nkey == "split" ? "split clusters" : "no-split"),
                        0.17, 0.81, 0.040);
    t->Draw();
  };
  drawPanel(1, h_data,   "Data");
  drawPanel(2, h_photon, "Photon MC pool");
  drawPanel(3, h_jet,    "Jet MC pool");

  EnsureDir(kFigDir);
  const std::string out = std::string(kFigDir) + "/cluster_size_2d_" +
                          m.key + "_" + v.key + "_" + nkey + "_" +
                          p.period_key + ".pdf";
  c->SaveAs(out.c_str());
  std::cout << "Wrote " << out << std::endl;
  delete c;
}

// -----------------------------------------------------------------------
// Type B:  12-panel 1D-binned-by-ET grid, overlay data/photon/jet
// -----------------------------------------------------------------------
// Project TH2F onto Y for a given ET range [et_lo, et_hi) and normalize to
// unit area.  Returns nullptr if total integral <= 0.
TH1D* ProjectYNorm(TH2F* h2, double et_lo, double et_hi,
                   const std::string& newname)
{
  if (!h2) return nullptr;
  const int bx_lo = h2->GetXaxis()->FindBin(et_lo + 1e-6);
  const int bx_hi = h2->GetXaxis()->FindBin(et_hi - 1e-6);
  TH1D* proj = h2->ProjectionY(newname.c_str(), bx_lo, bx_hi, "e");
  if (!proj) return nullptr;
  proj->SetDirectory(nullptr);
  const double integ = proj->Integral();
  if (integ > 0) proj->Scale(1.0 / integ);
  else { delete proj; return nullptr; }
  return proj;
}

void Figure_1D_ByET(
    const std::map<std::string, std::unique_ptr<SampleFile>>& all,
    const MetricSpec& m, const VariantSpec& v, const PeriodParams& p,
    const std::string& nkey = "split")
{
  const std::string hname = "h2_" + m.key + "_" + nkey + v.suffix;

  auto get = [&](const std::string& k) -> const SampleFile* {
    auto it = all.find(k);
    return it == all.end() ? nullptr : it->second.get();
  };

  const SampleFile* data = get("data");
  std::vector<const SampleFile*> photon_pool = {
      get("photon5"), get("photon10"), get("photon20")};
  std::vector<const SampleFile*> jet_pool = {
      get("jet5"), get("jet8"), get("jet12"),
      get("jet20"), get("jet30"), get("jet40")};

  TH2F* h_data   = data ? dynamic_cast<TH2F*>(data->file->Get(hname.c_str())) : nullptr;
  if (h_data) h_data->SetDirectory(nullptr);
  TH2F* h_photon = PoolH2(photon_pool, hname,
      Form("ph_byet_%s_%s%s_%s",
           m.key.c_str(), nkey.c_str(), v.suffix.c_str(), p.period_key.c_str()));
  TH2F* h_jet    = PoolH2(jet_pool, hname,
      Form("jt_byet_%s_%s%s_%s",
           m.key.c_str(), nkey.c_str(), v.suffix.c_str(), p.period_key.c_str()));

  // 12-panel canvas: 4 cols x 3 rows (ET bins 0..11 from ptRanges)
  auto* c = new TCanvas(
      Form("c1d_%s_%s_%s_%s", m.key.c_str(), v.key.c_str(),
           nkey.c_str(), p.period_key.c_str()),
      Form("cluster size 1D-byET %s %s %s %s",
           m.key.c_str(), v.key.c_str(), nkey.c_str(), p.period_key.c_str()),
      1800, 1400);
  c->Divide(4, 3, 0.002, 0.002);

  for (int ipt = 0; ipt < NptBins; ++ipt) {
    c->cd(ipt + 1);
    gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.09);
    gPad->SetBottomMargin(0.15);

    const double et_lo = ptRanges[ipt];
    const double et_hi = ptRanges[ipt + 1];

    const std::string tag =
        Form("_%s_%s%s_%s_pt%d", m.key.c_str(), nkey.c_str(),
             v.suffix.c_str(), p.period_key.c_str(), ipt);
    TH1D* p_data   = ProjectYNorm(h_data,   et_lo, et_hi, "pd" + tag);
    TH1D* p_photon = ProjectYNorm(h_photon, et_lo, et_hi, "pp" + tag);
    TH1D* p_jet    = ProjectYNorm(h_jet,    et_lo, et_hi, "pj" + tag);

    auto style = [&](TH1D* h, int color, int marker, int line_style = 1) {
      if (!h) return;
      h->SetStats(0);
      h->SetTitle("");
      h->SetLineColor(color);
      h->SetMarkerColor(color);
      h->SetMarkerStyle(marker);
      h->SetMarkerSize(0.9);
      h->SetLineWidth(2);
      h->SetLineStyle(line_style);
      h->GetXaxis()->SetTitle(m.ytitle.c_str());
      h->GetYaxis()->SetTitle("Normalized");
      h->GetXaxis()->SetTitleSize(0.060);
      h->GetYaxis()->SetTitleSize(0.060);
      h->GetXaxis()->SetLabelSize(0.050);
      h->GetYaxis()->SetLabelSize(0.050);
      h->GetYaxis()->SetTitleOffset(1.15);
      h->GetXaxis()->SetRangeUser(m.y_min, m.y_max);
    };
    style(p_data,   kBlack,  20);
    style(p_photon, kRed+1,  24, 1);
    style(p_jet,    kBlue+1, 25, 1);

    // Frame: whichever exists first among data/photon/jet.
    TH1D* frame = p_data ? p_data : (p_photon ? p_photon : p_jet);
    if (!frame) {
      auto* l = MakeLabel("no entries", 0.30, 0.50, 0.060);
      l->Draw();
      continue;
    }

    double ymax = 0;
    for (TH1D* h : {p_data, p_photon, p_jet}) {
      if (!h) continue;
      if (h->GetMaximum() > ymax) ymax = h->GetMaximum();
    }
    frame->SetMinimum(0.0);
    frame->SetMaximum(ymax * 1.35);
    frame->Draw("E1");
    if (p_photon && p_photon != frame) p_photon->Draw("E1 SAME");
    if (p_jet    && p_jet    != frame) p_jet->Draw("E1 SAME");
    if (p_data   && p_data   != frame) p_data->Draw("E1 SAME");

    auto* etlab = MakeLabel(
        Form("%.0f < E_{T} < %.0f GeV", et_lo, et_hi),
        0.22, 0.85, 0.050);
    etlab->Draw();

    if (ipt == 0) {
      sPhenixInternal(0.22, 0.77, p.period_label, p.lumi, 0.045);
      auto* sel = MakeLabel(
          Form("%s, %s", nkey == "split" ? "split clusters" : "no-split",
               v.label.c_str()),
          0.22, 0.66, 0.042);
      sel->Draw();

      auto* leg = new TLegend(0.55, 0.62, 0.95, 0.90);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.045);
      if (p_data)   leg->AddEntry(p_data,   "Data",            "lp");
      if (p_photon) leg->AddEntry(p_photon, "Photon MC pool",  "lp");
      if (p_jet)    leg->AddEntry(p_jet,    "Jet MC pool",     "lp");
      leg->Draw();
    }
  }

  EnsureDir(kFigDir);
  const std::string out = std::string(kFigDir) + "/cluster_size_1d_byet_" +
                          m.key + "_" + v.key + "_" + nkey + "_" +
                          p.period_key + ".pdf";
  c->SaveAs(out.c_str());
  std::cout << "Wrote " << out << std::endl;
  delete c;
}

}  // namespace

void plot_cluster_size_2d() {
  init_plot();
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  const std::vector<std::string> filetypes = {
      "data",
      "photon5", "photon10", "photon20",
      "jet5", "jet8", "jet12", "jet20", "jet30", "jet40",
  };

  for (const std::string& period : {std::string("1p5mrad"), std::string("0mrad")}) {
    std::cout << "\n=== period " << period << " ===\n";
    PeriodParams p = MakePeriodParams(period);

    std::map<std::string, std::unique_ptr<SampleFile>> all;
    for (const std::string& ft : filetypes) {
      auto sf = LoadSample(period, ft);
      if (sf && sf->file) {
        all[ft] = std::move(sf);
      } else {
        std::cerr << "Missing " << ft << " for period " << period << std::endl;
      }
    }

    for (const auto& m : kMetrics) {
      for (const auto& v : kVariants) {
        Figure_2D_MetricCutPeriod(all, m, v, p, "split");
        Figure_1D_ByET           (all, m, v, p, "split");
      }
    }
  }
}
