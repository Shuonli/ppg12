#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TSystem.h>

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
  float weight_rel    = 1.0f;
  float ref_xsec_pb   = 1.0f;
  bool  isbackground  = false;
};

std::unique_ptr<SampleFile> LoadSample(const std::string& filetype) {
  auto sf = std::make_unique<SampleFile>();
  sf->filetype = filetype;
  const std::string path =
      std::string(kResultsDir) + "/cluster_size_" + filetype + ".root";
  sf->file.reset(TFile::Open(path.c_str(), "READ"));
  if (!sf->file || sf->file->IsZombie()) {
    std::cerr << "Cannot open " << path << std::endl;
    sf->file.reset();
    return sf;
  }
  auto get_f = [&](const char* name) {
    auto* p = (TParameter<float>*)sf->file->Get(name);
    return p ? p->GetVal() : 1.0f;
  };
  auto get_i = [&](const char* name) {
    auto* p = (TParameter<int>*)sf->file->Get(name);
    return p ? p->GetVal() : 0;
  };
  sf->weight_rel   = get_f("weight_rel");
  sf->ref_xsec_pb  = get_f("ref_xsec_pb");
  sf->isbackground = (get_i("isbackground") == 1);
  return sf;
}

// Clone a TProfile from a sample.
TProfile* CloneProf(const SampleFile& sf,
                    const std::string& hname,
                    const std::string& newname)
{
  if (!sf.file) return nullptr;
  auto* src = dynamic_cast<TProfile*>(sf.file->Get(hname.c_str()));
  if (!src) {
    std::cerr << "Missing " << hname << " in " << sf.filetype << std::endl;
    return nullptr;
  }
  auto* cp = (TProfile*)src->Clone(newname.c_str());
  cp->SetDirectory(nullptr);
  return cp;
}

// Pool a list of profiles via sum_w-weighted averaging (TProfile::Add).
// The per-sample scale multiplies the sample's sum_w, sum_wy contribution.
// Each sample already contains xsec*entries weighting from the analysis macro,
// so passing unit scale gives a physics-weighted pooled mean across samples.
// NOTE: TProfile::Scale scales bin means (wrong for pooling); we Reset + Add
// for the first sample to avoid that pitfall.
TProfile* PoolProfs(
    const std::vector<std::pair<const SampleFile*, double>>& samples,
    const std::string& hname,
    const std::string& newname)
{
  TProfile* out = nullptr;
  for (const auto& [sf, scale] : samples) {
    if (!sf || !sf->file) continue;
    auto* src = dynamic_cast<TProfile*>(sf->file->Get(hname.c_str()));
    if (!src) continue;
    if (!out) {
      out = (TProfile*)src->Clone(newname.c_str());
      out->SetDirectory(nullptr);
      out->Reset();
      out->Add(src, scale);
    } else {
      out->Add(src, scale);
    }
  }
  return out;
}

// Linear combination of TProfile means at fixed physical fractions:
//   mean_out(bin) = (1-f) * mean_nom(bin) + f * mean_double(bin)
//   err_out(bin)  = sqrt((1-f)^2 * err_nom^2 + f^2 * err_double^2)
// Returned as TH1D since the result is a function of bin-wise means, not a
// pooling of underlying fills.
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
    if (nn <= 0 && nd <= 0) {
      h->SetBinContent(b, 0);
      h->SetBinError(b, 0);
      continue;
    }
    // If one side has no stats, fall back to the other.
    double mean, err;
    if (nd <= 0) { mean = mn;  err = en; }
    else if (nn <= 0) { mean = md;  err = ed; }
    else {
      mean = (1.0 - f) * mn + f * md;
      err  = std::sqrt((1.0 - f) * (1.0 - f) * en * en +
                       f * f * ed * ed);
    }
    h->SetBinContent(b, mean);
    h->SetBinError(b, err);
  }
  return h;
}

// Draw helper. Applies style, range, and labels for a cluster-size profile.
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

TLatex* MakeLabel(const std::string& txt, double x, double y,
                  double size = 0.045) {
  auto* l = new TLatex(x, y, txt.c_str());
  l->SetNDC();
  l->SetTextFont(42);
  l->SetTextSize(size);
  return l;
}

void sPhenixInternal(double x = 0.2, double y = 0.86) {
  auto* t1 = MakeLabel("#bf{#it{sPHENIX}} Internal", x, y);
  t1->Draw();
  auto* t2 = MakeLabel(
      "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV", x, y - 0.055, 0.04);
  t2->Draw();
}

void EnsureDir(const char* path) {
  gSystem->mkdir(path, true);
}

}  // namespace

// Figure 1: nominal comparison (data + photon MC inclusive + jet MC inclusive)
// 3 rows (n_owned, width_eta, width_phi) x 2 cols (split, nosplit)
void Figure_Nominal(
    const std::map<std::string, std::unique_ptr<SampleFile>>& all)
{
  const std::vector<std::pair<std::string, std::string>> metrics = {
      {"prof_n_owned",   "#LTn_{owned}#GT"},
      {"prof_width_eta", "#LT#Deltai_{#eta}#GT"},
      {"prof_width_phi", "#LT#Deltai_{#phi}#GT"},
  };
  const std::vector<std::pair<std::string, std::string>> nodes = {
      {"split",   "split clusters"},
      {"nosplit", "no-split clusters"},
  };

  auto get_sf = [&](const std::string& k) -> const SampleFile* {
    auto it = all.find(k);
    return it == all.end() ? nullptr : it->second.get();
  };

  auto TTo_2D_YRange = [&](TProfile* p) {
    if (!p) return;
    double maxv = p->GetMaximum();
    p->SetMinimum(0.0);
    p->SetMaximum(maxv * 1.25);
  };

  // Photon inclusive: photon5 + photon10 + photon20 (already stitched in macro)
  // Both in units where ref = photon20cross.
  // Jet inclusive: jet5 + jet8 + jet12 + jet20 + jet30 + jet40, ref = jet50cross.
  // To put them on a common "cluster-averaged" footing (photon-inclusive, jet-inclusive
  // each as mean cluster size), we do NOT need to re-scale because TProfile mean is
  // invariant under global scale; we Add with weight=1 within each class.

  auto photon_samples = std::vector<std::pair<const SampleFile*, double>>{
      {get_sf("photon5"),  1.0},
      {get_sf("photon10"), 1.0},
      {get_sf("photon20"), 1.0},
  };
  auto jet_samples = std::vector<std::pair<const SampleFile*, double>>{
      {get_sf("jet5"),  1.0},
      {get_sf("jet8"),  1.0},
      {get_sf("jet12"), 1.0},
      {get_sf("jet20"), 1.0},
      {get_sf("jet30"), 1.0},
      {get_sf("jet40"), 1.0},
  };

  auto* c = new TCanvas("c_nominal", "Nominal cluster size", 1400, 1400);
  c->Divide(2, 3, 0.003, 0.003);

  int pad = 0;
  for (const auto& [mkey, ytitle] : metrics) {
    for (const auto& [nkey, ntitle] : nodes) {
      ++pad;
      c->cd(pad);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.04);
      gPad->SetTopMargin(0.08);
      gPad->SetBottomMargin(0.13);

      const std::string hname = mkey + "_" + nkey;

      auto* data = get_sf("data");
      auto* p_data   = data ? CloneProf(*data, hname,
                         Form("data_%s_%s", mkey.c_str(), nkey.c_str())) : nullptr;

      auto* p_photon = PoolProfs(photon_samples, hname,
          Form("photon_%s_%s", mkey.c_str(), nkey.c_str()));
      auto* p_jet    = PoolProfs(jet_samples, hname,
          Form("jet_%s_%s", mkey.c_str(), nkey.c_str()));

      StyleProfile(p_data,   kBlack, 20);
      StyleProfile(p_photon, kRed+1, 24);
      StyleProfile(p_jet,    kBlue+1, 25);

      TProfile* first = p_data ? p_data : (p_photon ? p_photon : p_jet);
      if (!first) continue;

      first->GetXaxis()->SetRangeUser(4, 36);
      first->GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
      first->GetYaxis()->SetTitle(ytitle.c_str());
      first->GetXaxis()->SetTitleSize(0.055);
      first->GetYaxis()->SetTitleSize(0.055);
      first->GetXaxis()->SetLabelSize(0.045);
      first->GetYaxis()->SetLabelSize(0.045);
      first->GetYaxis()->SetTitleOffset(1.25);
      first->SetTitle("");

      // Y range: auto from data if available, else from photon
      double ymax = 0, ymin = 999;
      for (auto* h : {p_data, p_photon, p_jet}) {
        if (!h) continue;
        double hmax = h->GetMaximum();
        double hmin = h->GetMinimum();
        if (hmax > ymax) ymax = hmax;
        if (hmin < ymin) ymin = hmin;
      }
      if (ymax > 0) {
        first->SetMinimum(std::max(0.0, ymin - 0.15 * (ymax - ymin)));
        first->SetMaximum(ymax + 0.35 * (ymax - ymin));
      }

      first->Draw("E1");
      if (p_photon && p_photon != first) p_photon->Draw("E1 SAME");
      if (p_jet    && p_jet    != first) p_jet->Draw("E1 SAME");
      if (p_data   && p_data   != first) p_data->Draw("E1 SAME");

      sPhenixInternal(0.18, 0.88);
      auto* tn = MakeLabel(ntitle, 0.18, 0.72, 0.04);
      tn->Draw();

      if (pad <= 2) {
        auto* leg = new TLegend(0.60, 0.68, 0.93, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.038);
        if (p_data)   leg->AddEntry(p_data,   "Data",       "lp");
        if (p_photon) leg->AddEntry(p_photon, "Photon MC",  "lp");
        if (p_jet)    leg->AddEntry(p_jet,    "Jet MC",     "lp");
        leg->Draw();
      }
    }
  }

  EnsureDir(kFigDir);
  std::string out = std::string(kFigDir) + "/cluster_size_nominal.pdf";
  c->SaveAs(out.c_str());
  std::cout << "Wrote " << out << std::endl;
}

// Figure 2/3: DI comparison for a given base sample ("photon10" or "jet12").
// Compares: single (photon10_nom), double (photon10_double),
// mixed at 7.9%, mixed at 22.4%.
void Figure_DI(
    const std::map<std::string, std::unique_ptr<SampleFile>>& all,
    const std::string& base_label,
    const std::string& nom_key,      // e.g. "photon10_nom"
    const std::string& double_key)   // e.g. "photon10_double"
{
  const std::vector<std::pair<std::string, std::string>> metrics = {
      {"prof_n_owned",   "#LTn_{owned}#GT"},
      {"prof_width_eta", "#LT#Deltai_{#eta}#GT"},
      {"prof_width_phi", "#LT#Deltai_{#phi}#GT"},
  };
  const std::vector<std::pair<std::string, std::string>> nodes = {
      {"split",   "split clusters"},
      {"nosplit", "no-split clusters"},
  };

  auto get_sf = [&](const std::string& k) -> const SampleFile* {
    auto it = all.find(k);
    return it == all.end() ? nullptr : it->second.get();
  };

  const SampleFile* sf_nom    = get_sf(nom_key);
  const SampleFile* sf_double = get_sf(double_key);
  if (!sf_nom || !sf_double) {
    std::cerr << "DI pair missing for " << base_label << std::endl;
    return;
  }

  const double f_7p9  = 0.079;
  const double f_22p4 = 0.224;

  auto* c = new TCanvas(
      Form("c_di_%s", base_label.c_str()),
      Form("DI cluster size — %s", base_label.c_str()),
      1400, 1400);
  c->Divide(2, 3, 0.003, 0.003);

  int pad = 0;
  for (const auto& [mkey, ytitle] : metrics) {
    for (const auto& [nkey, ntitle] : nodes) {
      ++pad;
      c->cd(pad);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.04);
      gPad->SetTopMargin(0.08);
      gPad->SetBottomMargin(0.13);

      const std::string hname = mkey + "_" + nkey;

      auto* p_nom    = CloneProf(*sf_nom,    hname,
          Form("%s_nom_%s_%s",    base_label.c_str(), mkey.c_str(), nkey.c_str()));
      auto* p_double = CloneProf(*sf_double, hname,
          Form("%s_double_%s_%s", base_label.c_str(), mkey.c_str(), nkey.c_str()));

      // Linear mix of means at fixed physical fraction
      TH1D* h_mix_7  = BlendLinear(p_nom, p_double, f_7p9,
          Form("%s_mix_7p9_%s_%s", base_label.c_str(), mkey.c_str(), nkey.c_str()));
      TH1D* h_mix_22 = BlendLinear(p_nom, p_double, f_22p4,
          Form("%s_mix_22p4_%s_%s", base_label.c_str(), mkey.c_str(), nkey.c_str()));

      auto style_h1 = [](TH1D* h, int color, int marker, int line_style = 1) {
        if (!h) return;
        h->SetLineColor(color);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(marker);
        h->SetMarkerSize(0.9);
        h->SetLineWidth(2);
        h->SetLineStyle(line_style);
        h->SetStats(0);
      };

      StyleProfile(p_nom,    kBlue+1,  24);
      StyleProfile(p_double, kRed+1,   25);
      style_h1(h_mix_7,  kGreen+2, 26, 2);
      style_h1(h_mix_22, kMagenta+1, 32, 2);

      TProfile* first = p_nom;
      if (!first) continue;

      first->GetXaxis()->SetRangeUser(4, 36);
      first->GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
      first->GetYaxis()->SetTitle(ytitle.c_str());
      first->GetXaxis()->SetTitleSize(0.055);
      first->GetYaxis()->SetTitleSize(0.055);
      first->GetXaxis()->SetLabelSize(0.045);
      first->GetYaxis()->SetLabelSize(0.045);
      first->GetYaxis()->SetTitleOffset(1.25);
      first->SetTitle("");

      double ymax = 0, ymin = 999;
      auto update_range_p = [&](TProfile* h) {
        if (!h) return;
        for (int b = 1; b <= h->GetNbinsX(); ++b) {
          if (h->GetBinEntries(b) <= 0) continue;
          double v = h->GetBinContent(b);
          if (v > ymax) ymax = v;
          if (v < ymin && v > 0) ymin = v;
        }
      };
      auto update_range_h = [&](TH1D* h) {
        if (!h) return;
        for (int b = 1; b <= h->GetNbinsX(); ++b) {
          double v = h->GetBinContent(b);
          if (v <= 0) continue;
          if (v > ymax) ymax = v;
          if (v < ymin) ymin = v;
        }
      };
      update_range_p(p_nom);
      update_range_p(p_double);
      update_range_h(h_mix_7);
      update_range_h(h_mix_22);
      if (ymax > 0) {
        first->SetMinimum(std::max(0.0, ymin - 0.20 * (ymax - ymin)));
        first->SetMaximum(ymax + 0.45 * (ymax - ymin));
      }

      first->Draw("E1");
      if (p_double) p_double->Draw("E1 SAME");
      if (h_mix_7)  h_mix_7->Draw("E1 SAME");
      if (h_mix_22) h_mix_22->Draw("E1 SAME");
      p_nom->Draw("E1 SAME");

      sPhenixInternal(0.18, 0.88);
      auto* tn = MakeLabel(
          Form("%s: %s", base_label.c_str(), ntitle.c_str()),
          0.18, 0.72, 0.04);
      tn->Draw();

      if (pad <= 2) {
        auto* leg = new TLegend(0.56, 0.63, 0.93, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.034);
        if (p_nom)    leg->AddEntry(p_nom,    "single",              "lp");
        if (p_double) leg->AddEntry(p_double, "double",              "lp");
        if (h_mix_7)  leg->AddEntry(h_mix_7,  "mix 7.9% (1.5 mrad)", "lp");
        if (h_mix_22) leg->AddEntry(h_mix_22, "mix 22.4% (0 mrad)",  "lp");
        leg->Draw();
      }
    }
  }

  EnsureDir(kFigDir);
  std::string out = std::string(kFigDir) + "/cluster_size_di_" + base_label + ".pdf";
  c->SaveAs(out.c_str());
  std::cout << "Wrote " << out << std::endl;
}

// Figure 4: split/nosplit ratio
void Figure_SplitRatio(
    const std::map<std::string, std::unique_ptr<SampleFile>>& all)
{
  const std::vector<std::pair<std::string, std::string>> metrics = {
      {"prof_n_owned",   "n_{owned}"},
      {"prof_width_eta", "#Deltai_{#eta}"},
      {"prof_width_phi", "#Deltai_{#phi}"},
  };

  auto get_sf = [&](const std::string& k) -> const SampleFile* {
    auto it = all.find(k);
    return it == all.end() ? nullptr : it->second.get();
  };

  auto photon_samples = std::vector<std::pair<const SampleFile*, double>>{
      {get_sf("photon5"),  1.0},
      {get_sf("photon10"), 1.0},
      {get_sf("photon20"), 1.0},
  };
  auto jet_samples = std::vector<std::pair<const SampleFile*, double>>{
      {get_sf("jet5"),  1.0},
      {get_sf("jet8"),  1.0},
      {get_sf("jet12"), 1.0},
      {get_sf("jet20"), 1.0},
      {get_sf("jet30"), 1.0},
      {get_sf("jet40"), 1.0},
  };

  // Convert TProfile to TH1D of mean/error-on-mean, skipping empty/low-stat bins.
  auto prof_to_h1 = [](TProfile* p, const std::string& name) -> TH1D* {
    if (!p) return nullptr;
    TH1D* h = new TH1D(name.c_str(), p->GetTitle(),
                       p->GetNbinsX(),
                       p->GetXaxis()->GetXmin(),
                       p->GetXaxis()->GetXmax());
    h->SetDirectory(nullptr);
    for (int b = 1; b <= p->GetNbinsX(); ++b) {
      double n = p->GetBinEntries(b);
      if (n < 10) {
        h->SetBinContent(b, 0);
        h->SetBinError(b, 0);
        continue;
      }
      h->SetBinContent(b, p->GetBinContent(b));
      h->SetBinError(b, p->GetBinError(b));
    }
    return h;
  };

  auto* c = new TCanvas("c_split_ratio", "no-split / split ratio", 1400, 500);
  c->Divide(3, 1, 0.003, 0.003);

  int pad = 0;
  for (const auto& [mkey, ytitle] : metrics) {
    ++pad;
    c->cd(pad);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.14);

    auto build_ratio = [&](const std::string& label,
                           const std::vector<std::pair<const SampleFile*, double>>& samples,
                           int color, int marker) -> TH1D* {
      auto* p_split   = PoolProfs(samples, mkey + "_split",
                         Form("%s_split_%s",   label.c_str(), mkey.c_str()));
      auto* p_nosplit = PoolProfs(samples, mkey + "_nosplit",
                         Form("%s_nosplit_%s", label.c_str(), mkey.c_str()));
      auto* h_split   = prof_to_h1(p_split,
                         Form("h_split_%s_%s",   label.c_str(), mkey.c_str()));
      auto* h_nosplit = prof_to_h1(p_nosplit,
                         Form("h_nosplit_%s_%s", label.c_str(), mkey.c_str()));
      if (!h_split || !h_nosplit) return nullptr;
      auto* ratio = (TH1D*)h_nosplit->Clone(
          Form("ratio_%s_%s", label.c_str(), mkey.c_str()));
      ratio->Divide(h_split);
      ratio->SetLineColor(color);
      ratio->SetMarkerColor(color);
      ratio->SetMarkerStyle(marker);
      ratio->SetMarkerSize(0.9);
      ratio->SetLineWidth(2);
      ratio->SetStats(0);
      return ratio;
    };

    auto build_ratio_single = [&](const std::string& key,
                                  int color, int marker) -> TH1D* {
      const SampleFile* sf = get_sf(key);
      if (!sf || !sf->file) return nullptr;
      std::vector<std::pair<const SampleFile*, double>> v{{sf, 1.0}};
      return build_ratio(key, v, color, marker);
    };

    auto* r_data   = build_ratio_single("data", kBlack, 20);
    auto* r_photon = build_ratio("photon", photon_samples, kRed+1, 24);
    auto* r_jet    = build_ratio("jet",    jet_samples,    kBlue+1, 25);

    // Always use r_jet as the frame so the y-axis is set from its extents
    TH1D* first = r_jet ? r_jet : (r_photon ? r_photon : r_data);
    if (!first) continue;

    first->GetXaxis()->SetRangeUser(4, 36);
    first->GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
    first->GetYaxis()->SetTitle(Form("%s (no-split) / %s (split)",
                                      ytitle.c_str(), ytitle.c_str()));
    first->GetXaxis()->SetTitleSize(0.055);
    first->GetYaxis()->SetTitleSize(0.05);
    first->GetXaxis()->SetLabelSize(0.045);
    first->GetYaxis()->SetLabelSize(0.045);
    first->GetYaxis()->SetTitleOffset(1.35);
    first->SetTitle("");
    first->SetMinimum(0.95);
    first->SetMaximum(1.45);

    first->Draw("E1");
    if (r_photon && r_photon != first) r_photon->Draw("E1 SAME");
    if (r_data   && r_data   != first) r_data->Draw("E1 SAME");
    if (r_jet    && r_jet    != first) r_jet->Draw("E1 SAME");

    auto* line = new TLine(4, 1.0, 36, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray+2);
    line->Draw();

    sPhenixInternal(0.2, 0.88);

    if (pad == 1) {
      auto* leg = new TLegend(0.55, 0.68, 0.96, 0.92);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.04);
      if (r_data)   leg->AddEntry(r_data,   "Data",       "lp");
      if (r_photon) leg->AddEntry(r_photon, "Photon MC",  "lp");
      if (r_jet)    leg->AddEntry(r_jet,    "Jet MC",     "lp");
      leg->Draw();
    }
  }

  EnsureDir(kFigDir);
  std::string out = std::string(kFigDir) + "/cluster_size_split_ratio.pdf";
  c->SaveAs(out.c_str());
  std::cout << "Wrote " << out << std::endl;
}

void plot_cluster_size() {
  SetsPhenixStyle();
  gStyle->SetOptStat(0);

  const std::vector<std::string> filetypes = {
      "data",
      "photon5", "photon10", "photon20",
      "photon10_nom", "photon10_double",
      "jet5", "jet8", "jet12", "jet20", "jet30", "jet40",
      "jet12_nom", "jet12_double",
  };

  std::map<std::string, std::unique_ptr<SampleFile>> all;
  for (const auto& ft : filetypes) {
    auto sf = LoadSample(ft);
    if (sf && sf->file) {
      std::cout << "Loaded " << ft << "  weight=" << sf->weight_rel
                << " isbg=" << sf->isbackground << std::endl;
      all[ft] = std::move(sf);
    } else {
      std::cerr << "Missing " << ft << std::endl;
    }
  }

  Figure_Nominal(all);
  Figure_DI(all, "photon10", "photon10_nom", "photon10_double");
  Figure_DI(all, "jet12",    "jet12_nom",    "jet12_double");
  Figure_SplitRatio(all);
}
