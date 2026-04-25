// plot_cluster_size_cuts.C
//
// Per-period cluster-size figures for cluster_size_short.pdf, showing
// data vs pooled MC (photon inclusive, jet inclusive) with:
//   - baseline (no cuts, just |eta_reco|<0.7 fiducial and vertex cut)
//   - common cut
//   - tight cut
// with a quadratic trend line (pol2) overlaid on each curve.
//
// Produces:
//   cluster_size_period_baseline_{period}.pdf     3 rows x 1 col
//   cluster_size_period_cuts_{period}.pdf         3 rows x 2 cols
//   cluster_size_period_di_photon10_{period}.pdf  3 rows x 2 cols, no-cut + common
//   cluster_size_period_di_jet12_{period}.pdf     3 rows x 2 cols, no-cut + common
// for period in {1p5mrad, 0mrad}.
//
// DI mixing uses bin-by-bin  <y>_mix = (1-f)*<y>_single + f*<y>_double
// with f from the period configuration.

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "../efficiencytool/CrossSectionWeights.h"
#include "sPhenixStyle.C"

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
    return sf;
  }
  return sf;
}

// Pool TProfiles across a list of samples (Reset + Add). Returns heap
// ownership; caller must manage.  Empty result is returned as nullptr.
TProfile* PoolProfs(
    const std::vector<const SampleFile*>& samples,
    const std::string& hname,
    const std::string& newname)
{
  TProfile* out = nullptr;
  for (const SampleFile* sf : samples) {
    if (!sf || !sf->file) continue;
    auto* src = dynamic_cast<TProfile*>(sf->file->Get(hname.c_str()));
    if (!src) continue;
    if (!out) {
      out = (TProfile*)src->Clone(newname.c_str());
      out->SetDirectory(nullptr);
      out->Reset();
      out->Add(src, 1.0);
    } else {
      out->Add(src, 1.0);
    }
  }
  return out;
}

// Blend two TProfile means bin-by-bin at fixed double-interaction fraction.
// Result is a TH1D of per-bin (1-f)*<y>_nom + f*<y>_double. When one side
// has no entries in a bin, the other side's value is returned.
TH1D* BlendLinear(const TProfile* p_nom, const TProfile* p_double, double f,
                  const std::string& newname)
{
  if (!p_nom) return nullptr;
  auto* h = new TH1D(newname.c_str(), p_nom->GetTitle(),
                     p_nom->GetNbinsX(),
                     p_nom->GetXaxis()->GetXmin(),
                     p_nom->GetXaxis()->GetXmax());
  h->SetDirectory(nullptr);
  for (int b = 1; b <= p_nom->GetNbinsX(); ++b) {
    const double mn  = p_nom->GetBinContent(b);
    const double en  = p_nom->GetBinError(b);
    const double nn  = p_nom->GetBinEntries(b);
    double md = 0, ed = 0, nd = 0;
    if (p_double) {
      md = p_double->GetBinContent(b);
      ed = p_double->GetBinError(b);
      nd = p_double->GetBinEntries(b);
    }
    if (nn <= 0 && nd <= 0) { h->SetBinContent(b, 0); h->SetBinError(b, 0); continue; }
    double mean, err;
    if (nd <= 0)      { mean = mn; err = en; }
    else if (nn <= 0) { mean = md; err = ed; }
    else {
      mean = (1.0 - f) * mn + f * md;
      err  = std::sqrt((1.0 - f)*(1.0 - f)*en*en + f*f*ed*ed);
    }
    h->SetBinContent(b, mean);
    h->SetBinError(b, err);
  }
  return h;
}

void StyleProfile(TProfile* p, int color, int marker, int line_style = 1) {
  if (!p) return;
  p->SetLineColor(color);
  p->SetMarkerColor(color);
  p->SetMarkerStyle(marker);
  p->SetMarkerSize(0.9);
  p->SetLineWidth(2);
  p->SetLineStyle(line_style);
  p->SetStats(0);
}

void StyleH1(TH1* h, int color, int marker, int line_style = 1) {
  if (!h) return;
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(0.9);
  h->SetLineWidth(2);
  h->SetLineStyle(line_style);
  h->SetStats(0);
}

// Return a TF1 "pol2" fit to the non-empty bins of a TH1 (or TProfile) in
// ET-range [et_lo, et_hi], or nullptr if fewer than 3 good points.
TF1* FitTrend(TH1* h, double et_lo, double et_hi, int color,
              const std::string& name)
{
  if (!h) return nullptr;
  int n_good = 0;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    if (h->GetBinCenter(b) < et_lo || h->GetBinCenter(b) > et_hi) continue;
    if (h->GetBinError(b) <= 0) continue;
    if (!std::isfinite(h->GetBinContent(b))) continue;
    ++n_good;
  }
  if (n_good < 3) return nullptr;
  auto* f = new TF1(name.c_str(), "pol2", et_lo, et_hi);
  f->SetLineColor(color);
  f->SetLineStyle(2);
  f->SetLineWidth(2);
  h->Fit(f, "QN0R");
  return f;
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
                     const std::string& luminosity) {
  auto* t1 = MakeLabel("#bf{#it{sPHENIX}} Internal", x, y);
  t1->Draw();
  const std::string period_txt =
      "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, " + period +
      ", " + luminosity + " pb^{-1}";
  auto* t2 = MakeLabel(period_txt, x, y - 0.055, 0.034);
  t2->Draw();
}

void EnsureDir(const char* path) { gSystem->mkdir(path, true); }

// Adjust Y axis to bracket a set of histograms' bin contents (bin_entries>=1
// for TProfiles).  The `pad_frac_low / pad_frac_high` bracket the min/max.
void SetAutoYRange(TH1* frame,
                   const std::vector<TH1*>& hs,
                   double pad_frac_low = 0.2,
                   double pad_frac_high = 0.35,
                   double et_lo = 4, double et_hi = 36)
{
  double ymin = 1e9, ymax = -1e9;
  for (auto* h : hs) {
    if (!h) continue;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      const double ec = h->GetBinCenter(b);
      if (ec < et_lo || ec > et_hi) continue;
      if (h->GetBinError(b) <= 0) continue;
      const double v = h->GetBinContent(b);
      if (!std::isfinite(v) || v == 0) continue;
      if (v > ymax) ymax = v;
      if (v < ymin) ymin = v;
    }
  }
  if (ymax < ymin) { ymin = 0; ymax = 1; }
  const double span = ymax - ymin;
  frame->SetMinimum(std::max(0.0, ymin - pad_frac_low * span));
  frame->SetMaximum(ymax + pad_frac_high * span);
}

struct PeriodParams {
  std::string period_key;       // "1p5mrad" or "0mrad"
  std::string period_label;     // "1.5 mrad" / "0 mrad"
  std::string lumi;             // "16.86" / "32.66"
  double f_double;              // DI fraction
};

PeriodParams MakePeriodParams(const std::string& pk) {
  PeriodParams p;
  p.period_key = pk;
  if (pk == "1p5mrad") {
    p.period_label = "1.5 mrad";
    p.lumi         = "16.86";
    p.f_double     = 0.079;
  } else {
    p.period_label = "0 mrad";
    p.lumi         = "32.66";
    p.f_double     = 0.224;
  }
  return p;
}

}  // namespace

// -----------------------------------------------------------------------
// Figure: 3 rows (metrics) x 1 col (split container), data vs pooled MC
// (single + DI-blended) with trend lines.  Used both for baseline and
// for each cut variant.
// -----------------------------------------------------------------------
void DrawMetricPanel(
    int pad_index,
    const std::string& mkey,
    const std::string& mtitle,
    const std::string& nkey,                // "split" or "nosplit"
    const std::string& suffix,              // "" or "_common" or "_tight"
    const PeriodParams& p,
    const SampleFile* data,
    const std::vector<const SampleFile*>& photon_pool,
    const std::vector<const SampleFile*>& jet_pool,
    const SampleFile* ph_nom,               // for DI blend
    const SampleFile* ph_dbl,
    const SampleFile* jt_nom,
    const SampleFile* jt_dbl,
    TLegend*& out_leg,
    bool draw_legend_here)
{
  const std::string hname = "prof_" + mkey + "_" + nkey + suffix;

  TProfile* p_data   = data ? (TProfile*)data->file->Get(hname.c_str()) : nullptr;

  TProfile* p_photon_single = PoolProfs(
      photon_pool, hname,
      Form("ph_pool_%s_%s%s", mkey.c_str(), nkey.c_str(), suffix.c_str()));
  TProfile* p_jet_single    = PoolProfs(
      jet_pool, hname,
      Form("jt_pool_%s_%s%s", mkey.c_str(), nkey.c_str(), suffix.c_str()));

  // DI blending at the pool level: use photon10_double (or jet12_double) as
  // the DI "addend" in its ET range; outside that range the single-MC pool
  // is used unchanged (the blend weight becomes 0 there).
  TProfile* ph_nom_p    = ph_nom    ? (TProfile*)ph_nom   ->file->Get(hname.c_str()) : nullptr;
  TProfile* ph_dbl_p    = ph_dbl    ? (TProfile*)ph_dbl   ->file->Get(hname.c_str()) : nullptr;
  TProfile* jt_nom_p    = jt_nom    ? (TProfile*)jt_nom   ->file->Get(hname.c_str()) : nullptr;
  TProfile* jt_dbl_p    = jt_dbl    ? (TProfile*)jt_dbl   ->file->Get(hname.c_str()) : nullptr;

  // Replace single pool mean with blended mean in bins where DI addend has stats.
  TH1D* p_photon_mix = nullptr;
  TH1D* p_jet_mix    = nullptr;
  if (p_photon_single) {
    p_photon_mix = new TH1D(Form("ph_mix_%s_%s%s", mkey.c_str(), nkey.c_str(), suffix.c_str()),
                            "",
                            p_photon_single->GetNbinsX(),
                            p_photon_single->GetXaxis()->GetXmin(),
                            p_photon_single->GetXaxis()->GetXmax());
    p_photon_mix->SetDirectory(nullptr);
    for (int b = 1; b <= p_photon_single->GetNbinsX(); ++b) {
      const double m_s = p_photon_single->GetBinContent(b);
      const double e_s = p_photon_single->GetBinError(b);
      const double n_s = p_photon_single->GetBinEntries(b);
      if (n_s <= 0) { p_photon_mix->SetBinContent(b, 0); p_photon_mix->SetBinError(b, 0); continue; }
      const double m_d = ph_dbl_p ? ph_dbl_p->GetBinContent(b) : 0;
      const double e_d = ph_dbl_p ? ph_dbl_p->GetBinError(b)   : 0;
      const double n_d = ph_dbl_p ? ph_dbl_p->GetBinEntries(b) : 0;
      const double m_n = ph_nom_p ? ph_nom_p->GetBinContent(b) : 0;
      const double n_n = ph_nom_p ? ph_nom_p->GetBinEntries(b) : 0;
      if (n_d > 0 && n_n > 0) {
        // DI-induced shift is (m_d - m_n); apply on top of the single pool.
        // Error carried in quadrature from (single, nom, double) contributions.
        const double e_n = (ph_nom_p ? ph_nom_p->GetBinError(b) : 0.0);
        const double shift = (m_d - m_n);
        p_photon_mix->SetBinContent(b, m_s + p.f_double * shift);
        p_photon_mix->SetBinError(b,
            std::sqrt(e_s*e_s + p.f_double*p.f_double*(e_d*e_d + e_n*e_n)));
      } else {
        p_photon_mix->SetBinContent(b, m_s);
        p_photon_mix->SetBinError(b, e_s);
      }
    }
  }
  if (p_jet_single) {
    p_jet_mix = new TH1D(Form("jt_mix_%s_%s%s", mkey.c_str(), nkey.c_str(), suffix.c_str()),
                         "",
                         p_jet_single->GetNbinsX(),
                         p_jet_single->GetXaxis()->GetXmin(),
                         p_jet_single->GetXaxis()->GetXmax());
    p_jet_mix->SetDirectory(nullptr);
    for (int b = 1; b <= p_jet_single->GetNbinsX(); ++b) {
      const double m_s = p_jet_single->GetBinContent(b);
      const double e_s = p_jet_single->GetBinError(b);
      const double n_s = p_jet_single->GetBinEntries(b);
      if (n_s <= 0) { p_jet_mix->SetBinContent(b, 0); p_jet_mix->SetBinError(b, 0); continue; }
      const double m_d = jt_dbl_p ? jt_dbl_p->GetBinContent(b) : 0;
      const double e_d = jt_dbl_p ? jt_dbl_p->GetBinError(b)   : 0;
      const double n_d = jt_dbl_p ? jt_dbl_p->GetBinEntries(b) : 0;
      const double m_n = jt_nom_p ? jt_nom_p->GetBinContent(b) : 0;
      const double n_n = jt_nom_p ? jt_nom_p->GetBinEntries(b) : 0;
      if (n_d > 0 && n_n > 0) {
        const double e_n = (jt_nom_p ? jt_nom_p->GetBinError(b) : 0.0);
        const double shift = (m_d - m_n);
        p_jet_mix->SetBinContent(b, m_s + p.f_double * shift);
        p_jet_mix->SetBinError(b,
            std::sqrt(e_s*e_s + p.f_double*p.f_double*(e_d*e_d + e_n*e_n)));
      } else {
        p_jet_mix->SetBinContent(b, m_s);
        p_jet_mix->SetBinError(b, e_s);
      }
    }
  }

  StyleProfile(p_data,          kBlack,   20);
  StyleProfile(p_photon_single, kRed+1,   24);
  StyleProfile(p_jet_single,    kBlue+1,  25);
  StyleH1(p_photon_mix,         kRed+1,   26, 2);
  StyleH1(p_jet_mix,             kBlue+1, 32, 2);

  gPad->SetLeftMargin(0.17);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.14);

  // Frame: use whichever of data/photon/jet has most points
  TH1* frame = nullptr;
  if (p_data)              frame = p_data;
  else if (p_photon_single) frame = p_photon_single;
  else if (p_jet_single)   frame = p_jet_single;
  if (!frame) return;

  frame->GetXaxis()->SetRangeUser(4, 36);
  frame->GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
  frame->GetYaxis()->SetTitle(mtitle.c_str());
  frame->GetXaxis()->SetTitleSize(0.055);
  frame->GetYaxis()->SetTitleSize(0.055);
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.3);
  frame->SetTitle("");
  SetAutoYRange(frame, {p_data, p_photon_single, p_jet_single, p_photon_mix, p_jet_mix});

  frame->Draw("E1");
  if (p_photon_single && p_photon_single != frame) p_photon_single->Draw("E1 SAME");
  if (p_jet_single    && p_jet_single    != frame) p_jet_single->Draw("E1 SAME");
  if (p_photon_mix && p_photon_mix->GetEntries() > 0) p_photon_mix->Draw("HIST SAME");
  if (p_jet_mix    && p_jet_mix->GetEntries()    > 0) p_jet_mix->Draw("HIST SAME");
  if (p_data && p_data != frame) p_data->Draw("E1 SAME");

  // Trend lines (pol2) for data, single photon, single jet.
  if (p_data)          FitTrend(p_data,          6, 32, kBlack,  Form("tr_data_%s_%s%s",   mkey.c_str(), nkey.c_str(), suffix.c_str()))->Draw("LSAME");
  if (p_photon_single) FitTrend(p_photon_single, 6, 32, kRed+1,  Form("tr_ph_%s_%s%s",     mkey.c_str(), nkey.c_str(), suffix.c_str()))->Draw("LSAME");
  if (p_jet_single)    FitTrend(p_jet_single,    6, 32, kBlue+1, Form("tr_jt_%s_%s%s",     mkey.c_str(), nkey.c_str(), suffix.c_str()))->Draw("LSAME");

  sPhenixInternal(0.20, 0.88, p.period_label, p.lumi);

  const std::string sel_label =
      suffix.empty() ? "baseline" :
      (suffix == "_common" ? "common cut" : "tight cut");
  auto* tn = MakeLabel(
      Form("%s, %s", nkey == "split" ? "split clusters" : "no-split clusters", sel_label.c_str()),
      0.20, 0.76, 0.038);
  tn->Draw();

  if (draw_legend_here) {
    auto* leg = new TLegend(0.55, 0.62, 0.94, 0.91);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.034);
    if (p_data)          leg->AddEntry(p_data,          "Data",                            "lp");
    if (p_photon_single) leg->AddEntry(p_photon_single, "Photon MC (single)",              "lp");
    if (p_jet_single)    leg->AddEntry(p_jet_single,    "Jet MC (single)",                 "lp");
    if (p_photon_mix)    leg->AddEntry(p_photon_mix,    Form("Photon MC (DI-blend %.1f%%)", p.f_double*100), "l");
    if (p_jet_mix)       leg->AddEntry(p_jet_mix,       Form("Jet MC (DI-blend %.1f%%)",    p.f_double*100), "l");
    leg->Draw();
    out_leg = leg;
  }
}

void Figure_OneCutVariant(
    const std::map<std::string, std::unique_ptr<SampleFile>>& all,
    const PeriodParams& p,
    const std::string& suffix,       // "" | "_common" | "_tight"
    const std::string& variant_label,// "baseline" | "common" | "tight"
    const std::string& nkey)         // "split" | "nosplit"
{
  const std::vector<std::pair<std::string, std::string>> metrics = {
      {"n_owned",   "#LTn_{owned}#GT"},
      {"width_eta", "#LT#Deltai_{#eta}#GT"},
      {"width_phi", "#LT#Deltai_{#phi}#GT"},
  };

  auto get = [&](const std::string& k) -> const SampleFile* {
    auto it = all.find(k);
    return it == all.end() ? nullptr : it->second.get();
  };

  std::vector<const SampleFile*> photon_pool = {
      get("photon5"), get("photon10"), get("photon20")};
  std::vector<const SampleFile*> jet_pool = {
      get("jet5"), get("jet8"), get("jet12"), get("jet20"), get("jet30"), get("jet40")};

  auto* c = new TCanvas(
      Form("c_%s_%s_%s", variant_label.c_str(), nkey.c_str(), p.period_key.c_str()),
      Form("cluster size %s %s %s", variant_label.c_str(), nkey.c_str(), p.period_key.c_str()),
      700, 1500);
  c->Divide(1, 3, 0.003, 0.003);

  TLegend* leg = nullptr;
  int pad = 0;
  for (const auto& m : metrics) {
    c->cd(++pad);
    DrawMetricPanel(pad, m.first, m.second, nkey, suffix, p,
                    get("data"),
                    photon_pool, jet_pool,
                    get("photon10_nom"), get("photon10_double"),
                    get("jet12_nom"),    get("jet12_double"),
                    leg, /*draw_legend_here=*/ (pad == 1));
  }

  EnsureDir(kFigDir);
  std::string out = std::string(kFigDir) + "/cluster_size_period_" +
                    variant_label + "_" + nkey + "_" + p.period_key + ".pdf";
  c->SaveAs(out.c_str());
  std::cout << "Wrote " << out << std::endl;
}

// -----------------------------------------------------------------------
// Overlay figure: for a given period and node, each panel (one per metric)
// overlays:
//   Data: baseline, common-pass, tight-pass
//   Jet MC: baseline, common-pass, tight-pass
//   Photon MC: baseline, tight-pass (as requested, photon MC with the same
//              tight selection gives the direct comparator for data-tight)
// Trend lines (pol2) are drawn on the DATA curves to emphasize the
// progression across selections.  DI-blended curves are added on top where
// the DI samples span the ET range.
// -----------------------------------------------------------------------
void DrawOverlayPanel(
    int pad_index,
    const std::string& mkey,
    const std::string& mtitle,
    const std::string& nkey,                // "split" / "nosplit"
    const PeriodParams& p,
    const SampleFile* data,
    const std::vector<const SampleFile*>& photon_pool,
    const std::vector<const SampleFile*>& jet_pool,
    const SampleFile* ph_dbl,
    const SampleFile* jt_dbl,
    bool draw_legend_here)
{
  gPad->SetLeftMargin(0.17);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.14);

  // ---- Data: three selections
  auto load_data_prof = [&](const std::string& suffix) -> TProfile* {
    if (!data || !data->file) return nullptr;
    const std::string hname = "prof_" + mkey + "_" + nkey + suffix;
    return dynamic_cast<TProfile*>(data->file->Get(hname.c_str()));
  };
  TProfile* d_base   = load_data_prof("");
  TProfile* d_common = load_data_prof("_common");
  TProfile* d_tight  = load_data_prof("_tight");

  // ---- Pooled jet MC: three selections
  auto pool_suffix = [&](const std::vector<const SampleFile*>& pool,
                         const std::string& label,
                         const std::string& suffix) -> TProfile* {
    const std::string hname = "prof_" + mkey + "_" + nkey + suffix;
    return PoolProfs(pool, hname,
        Form("%s_pool_%s_%s%s", label.c_str(), mkey.c_str(), nkey.c_str(), suffix.c_str()));
  };

  TProfile* j_base   = pool_suffix(jet_pool,    "jt", "");
  TProfile* j_common = pool_suffix(jet_pool,    "jt", "_common");
  TProfile* j_tight  = pool_suffix(jet_pool,    "jt", "_tight");
  TProfile* ph_base  = pool_suffix(photon_pool, "ph", "");
  TProfile* ph_tight = pool_suffix(photon_pool, "ph", "_tight");

  // ---- Styles
  // Data: filled black markers increasing in "emphasis" as the selection tightens.
  StyleProfile(d_base,   kBlack,        24);  // open circle
  StyleProfile(d_common, kGray+2,       25);  // open square
  StyleProfile(d_tight,  kBlack,        20);  // filled circle

  // Jet MC: blue, open → filled as selection tightens
  StyleProfile(j_base,   kBlue+1,       26);  // open triangle
  StyleProfile(j_common, kAzure+2,      32);  // open diamond
  StyleProfile(j_tight,  kBlue+1,       22);  // filled triangle

  // Photon MC: red, tight-passing only (the direct comparator for data-tight);
  // baseline is shown as a lighter curve for context.
  StyleProfile(ph_base,  kRed-7,        27);  // open star
  StyleProfile(ph_tight, kRed+1,        33);  // filled diamond

  // ---- Frame
  TH1* frame = d_base ? (TH1*)d_base : (TH1*)j_base;
  if (!frame) return;
  frame->GetXaxis()->SetRangeUser(4, 36);
  frame->GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
  frame->GetYaxis()->SetTitle(mtitle.c_str());
  frame->GetXaxis()->SetTitleSize(0.055);
  frame->GetYaxis()->SetTitleSize(0.055);
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.3);
  frame->SetTitle("");
  SetAutoYRange(frame,
      {d_base, d_common, d_tight, j_base, j_common, j_tight, ph_base, ph_tight},
      0.15, 0.55);

  frame->Draw("E1");
  if (ph_base   && ph_base   != frame) ph_base  ->Draw("E1 SAME");
  if (j_base    && j_base    != frame) j_base   ->Draw("E1 SAME");
  if (j_common  && j_common  != frame) j_common ->Draw("E1 SAME");
  if (j_tight   && j_tight   != frame) j_tight  ->Draw("E1 SAME");
  if (ph_tight  && ph_tight  != frame) ph_tight ->Draw("E1 SAME");
  if (d_common  && d_common  != frame) d_common ->Draw("E1 SAME");
  if (d_tight   && d_tight   != frame) d_tight  ->Draw("E1 SAME");
  if (d_base    && d_base    != frame) d_base   ->Draw("E1 SAME");

  // Trend lines: pol2 on the key curves so the reader sees the progression
  // (data baseline, data tight, photon tight, jet tight).
  auto drawTrend = [&](TProfile* h, int color, const std::string& tag) {
    if (!h) return;
    TF1* f = FitTrend(h, 6, 32, color,
        Form("tr_%s_%s_%s_%s", tag.c_str(), mkey.c_str(), nkey.c_str(), p.period_key.c_str()));
    if (f) f->Draw("LSAME");
  };
  drawTrend(d_base,  kBlack,  "d_base");
  drawTrend(d_tight, kBlack,  "d_tight");
  drawTrend(ph_tight, kRed+1, "ph_tight");
  drawTrend(j_tight, kBlue+1, "j_tight");

  sPhenixInternal(0.20, 0.88, p.period_label, p.lumi);
  auto* tn = MakeLabel(
      Form("%s", nkey == "split" ? "split clusters" : "no-split clusters"),
      0.20, 0.76, 0.038);
  tn->Draw();

  if (draw_legend_here) {
    auto* leg = new TLegend(0.47, 0.53, 0.95, 0.91);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.029);
    leg->SetNColumns(2);
    // Column 1: data selections
    if (d_base)   leg->AddEntry(d_base,   "Data (all)",       "lp");
    // Column 2: jet MC selections
    if (j_base)   leg->AddEntry(j_base,   "Jet MC (all)",     "lp");
    if (d_common) leg->AddEntry(d_common, "Data (common)",    "lp");
    if (j_common) leg->AddEntry(j_common, "Jet MC (common)",  "lp");
    if (d_tight)  leg->AddEntry(d_tight,  "Data (tight)",     "lp");
    if (j_tight)  leg->AddEntry(j_tight,  "Jet MC (tight)",   "lp");
    if (ph_base)  leg->AddEntry(ph_base,  "Photon MC (all)",  "lp");
    if (ph_tight) leg->AddEntry(ph_tight, "Photon MC (tight)","lp");
    leg->Draw();
  }
}

void Figure_Overlay(
    const std::map<std::string, std::unique_ptr<SampleFile>>& all,
    const PeriodParams& p,
    const std::string& nkey)
{
  const std::vector<std::pair<std::string, std::string>> metrics = {
      {"n_owned",   "#LTn_{owned}#GT"},
      {"width_eta", "#LT#Deltai_{#eta}#GT"},
      {"width_phi", "#LT#Deltai_{#phi}#GT"},
  };

  auto get = [&](const std::string& k) -> const SampleFile* {
    auto it = all.find(k);
    return it == all.end() ? nullptr : it->second.get();
  };

  std::vector<const SampleFile*> photon_pool = {
      get("photon5"), get("photon10"), get("photon20")};
  std::vector<const SampleFile*> jet_pool = {
      get("jet5"), get("jet8"), get("jet12"), get("jet20"), get("jet30"), get("jet40")};

  auto* c = new TCanvas(
      Form("c_overlay_%s_%s", nkey.c_str(), p.period_key.c_str()),
      Form("cluster size overlay %s %s", nkey.c_str(), p.period_key.c_str()),
      900, 1500);
  c->Divide(1, 3, 0.003, 0.003);

  int pad = 0;
  for (const auto& m : metrics) {
    c->cd(++pad);
    DrawOverlayPanel(pad, m.first, m.second, nkey, p,
                     get("data"), photon_pool, jet_pool,
                     get("photon10_double"), get("jet12_double"),
                     /*draw_legend_here=*/ (pad == 1));
  }

  EnsureDir(kFigDir);
  std::string out = std::string(kFigDir) + "/cluster_size_period_overlay_" +
                    nkey + "_" + p.period_key + ".pdf";
  c->SaveAs(out.c_str());
  std::cout << "Wrote " << out << std::endl;
}

// 3 rows x 2 cols (common, tight), split container — primary cut figure
void Figure_CutsCombined(
    const std::map<std::string, std::unique_ptr<SampleFile>>& all,
    const PeriodParams& p)
{
  const std::vector<std::pair<std::string, std::string>> metrics = {
      {"n_owned",   "#LTn_{owned}#GT"},
      {"width_eta", "#LT#Deltai_{#eta}#GT"},
      {"width_phi", "#LT#Deltai_{#phi}#GT"},
  };
  const std::vector<std::pair<std::string, std::string>> variants = {
      {"_common", "common"},
      {"_tight",  "tight"},
  };

  auto get = [&](const std::string& k) -> const SampleFile* {
    auto it = all.find(k);
    return it == all.end() ? nullptr : it->second.get();
  };

  std::vector<const SampleFile*> photon_pool = {
      get("photon5"), get("photon10"), get("photon20")};
  std::vector<const SampleFile*> jet_pool = {
      get("jet5"), get("jet8"), get("jet12"), get("jet20"), get("jet30"), get("jet40")};

  auto* c = new TCanvas(
      Form("c_cuts_combined_%s", p.period_key.c_str()),
      Form("cluster size cuts %s", p.period_key.c_str()),
      1400, 1500);
  c->Divide(2, 3, 0.003, 0.003);

  TLegend* leg = nullptr;
  int pad = 0;
  for (const auto& m : metrics) {
    for (const auto& v : variants) {
      ++pad;
      c->cd(pad);
      DrawMetricPanel(pad, m.first, m.second, "split", v.first, p,
                      get("data"),
                      photon_pool, jet_pool,
                      get("photon10_nom"), get("photon10_double"),
                      get("jet12_nom"),    get("jet12_double"),
                      leg, /*draw_legend_here=*/ (pad == 1));
    }
  }

  EnsureDir(kFigDir);
  std::string out = std::string(kFigDir) + "/cluster_size_period_cuts_split_" +
                    p.period_key + ".pdf";
  c->SaveAs(out.c_str());
  std::cout << "Wrote " << out << std::endl;
}

void plot_cluster_size_cuts() {
  SetsPhenixStyle();
  gStyle->SetOptStat(0);

  const std::vector<std::string> filetypes = {
      "data",
      "photon5", "photon10", "photon20",
      "photon10_nom", "photon10_double",
      "jet5", "jet8", "jet12", "jet20", "jet30", "jet40",
      "jet12_nom", "jet12_double",
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

    // Overlay figure: primary — data and MC at baseline / common / tight, plus
    // photon-MC tight-pass curve for the direct signal-MC comparator.
    Figure_Overlay(all, p, "split");
    Figure_Overlay(all, p, "nosplit");

    // Baseline (no cuts), for reference / continuity with the March report.
    Figure_OneCutVariant(all, p, "",        "baseline", "split");
    Figure_OneCutVariant(all, p, "",        "baseline", "nosplit");

    // Combined cuts figure (common + tight side-by-side) — retained for
    // readers who prefer separated panels.
    Figure_CutsCombined(all, p);

    // Secondary: individual cut figures for nosplit container (sim single
    // MC only — DI MC has no cut branches in nosplit).  Provides context
    // for the split-ratio discussion.
    Figure_OneCutVariant(all, p, "_common", "common",   "nosplit");
    Figure_OneCutVariant(all, p, "_tight",  "tight",    "nosplit");
  }
}
