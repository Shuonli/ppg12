#include "plotcommon.h"

namespace
{
  TH1* getHist(TFile* f, const std::string& name)
  {
    if(!f)
    {
      std::cout << "[plot_sideband_compare_left_right] ERROR: null TFile for " << name << std::endl;
      return nullptr;
    }
    TObject* obj = f->Get(name.c_str());
    if(!obj)
    {
      std::cout << "[plot_sideband_compare_left_right] WARNING: missing object '" << name
                << "' in file " << f->GetName() << std::endl;
      return nullptr;
    }
    TH1* h = dynamic_cast<TH1*>(obj);
    if(!h)
    {
      std::cout << "[plot_sideband_compare_left_right] WARNING: object '" << name
                << "' is not a TH1 (class=" << obj->ClassName() << ") in file " << f->GetName() << std::endl;
      return nullptr;
    }
    return h;
  }

  TH1F* makeFrameFrom(const TH1* href, const std::string& frameName, const std::string& ytitle)
  {
    // Clone to inherit binning (works for variable binning too), then reset contents.
    TH1F* fr = (TH1F*)href->Clone(frameName.c_str());
    fr->Reset("ICES");
    fr->SetTitle("");
    fr->GetYaxis()->SetTitle(ytitle.c_str());
    return fr;
  }

  void styleHist(TH1* h, int color, int markerStyle)
  {
    if(!h) return;
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(markerStyle);
    h->SetMarkerSize(1.0);
    h->SetLineWidth(2);
  }

  void ensureSumw2(TH1* h)
  {
    if(!h) return;
    if(h->GetSumw2N() == 0) h->Sumw2();
  }

  void drawCompare1D(
      TH1* h_right,
      TH1* h_left,
      const std::string& outBase,
      const std::string& plotTitle,
      bool logy,
      const std::string& ytitle,
      double ymin,
      double ymax)
  {
    if(!h_right || !h_left) return;

    // Clone so we can style safely.
    TH1* hr = (TH1*)h_right->Clone(Form("%s_right_clone", outBase.c_str()));
    TH1* hl = (TH1*)h_left->Clone(Form("%s_left_clone", outBase.c_str()));
    hr->SetDirectory(nullptr);
    hl->SetDirectory(nullptr);
    ensureSumw2(hr);
    ensureSumw2(hl);

    styleHist(hr, kBlack, 20);
    styleHist(hl, kRed, 24);

    TH1F* frame = makeFrameFrom(hr, Form("frame_%s", outBase.c_str()), ytitle);
    frame->GetXaxis()->SetTitle(hr->GetXaxis()->GetTitle());
    frame->GetXaxis()->SetRangeUser(hr->GetXaxis()->GetXmin(), hr->GetXaxis()->GetXmax());

    double maxv = std::max(hr->GetMaximum(), hl->GetMaximum());
    double minv = 0.0;
    for(int ib=1; ib<=hr->GetNbinsX(); ++ib)
    {
      double v = hr->GetBinContent(ib);
      if(v>0 && (minv==0.0 || v<minv)) minv = v;
    }
    for(int ib=1; ib<=hl->GetNbinsX(); ++ib)
    {
      double v = hl->GetBinContent(ib);
      if(v>0 && (minv==0.0 || v<minv)) minv = v;
    }

    if(ymax <= ymin)
    {
      if(logy)
      {
        frame->GetYaxis()->SetRangeUser(std::max(1e-6, minv*0.5), maxv*5.0);
      }
      else
      {
        frame->GetYaxis()->SetRangeUser(0.0, maxv*1.3);
      }
    }
    else
    {
      frame->GetYaxis()->SetRangeUser(ymin, ymax);
    }

    TCanvas* c = new TCanvas(Form("c_%s", outBase.c_str()), Form("c_%s", outBase.c_str()), 650, 600);
    frame->Draw("axis");
    hr->Draw("same E1X0");
    hl->Draw("same E1X0");

    myText(0.18, 0.90, 1, strleg1.c_str(), 0.04);
    myText(0.18, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.18, 0.80, 1, strleg3.c_str(), 0.04);
    myText(0.18, 0.75, 1, plotTitle.c_str(), 0.04);

    myMarkerLineText(0.58, 0.86, 0, kBlack, 20, kBlack, 1, "#eta > 0", 0.045, true);
    myMarkerLineText(0.58, 0.80, 0, kRed,   24, kRed,   1, "#eta < 0", 0.045, true);

    gPad->SetLogy(logy ? 1 : 0);

    gSystem->mkdir("figures", kTRUE);
    c->SaveAs(Form("figures/%s.pdf", outBase.c_str()));

    delete frame;
    delete hr;
    delete hl;
  }

  void drawRatio1D(
      TH1* h_right,
      TH1* h_left,
      const std::string& outBase,
      const std::string& plotTitle,
      const std::string& ytitle,
      double ymin,
      double ymax)
  {
    if(!h_right || !h_left) return;

    TH1* hr = (TH1*)h_right->Clone(Form("%s_right_clone_ratio", outBase.c_str()));
    TH1* hl = (TH1*)h_left->Clone(Form("%s_left_clone_ratio", outBase.c_str()));
    hr->SetDirectory(nullptr);
    hl->SetDirectory(nullptr);
    ensureSumw2(hr);
    ensureSumw2(hl);

    TH1* hratio = (TH1*)hl->Clone(Form("%s_ratio", outBase.c_str())); // (eta<0)/(eta>0) = left/right
    hratio->SetDirectory(nullptr);
    ensureSumw2(hratio);
    hratio->Divide(hr);
    styleHist(hratio, kBlack, 20);

    TH1F* frame = makeFrameFrom(hratio, Form("frame_ratio_%s", outBase.c_str()), ytitle);
    frame->GetXaxis()->SetTitle(hratio->GetXaxis()->GetTitle());
    frame->GetXaxis()->SetRangeUser(hratio->GetXaxis()->GetXmin(), hratio->GetXaxis()->GetXmax());
    frame->GetYaxis()->SetRangeUser(ymin, ymax);
    frame->SetXTitle("#it{E}_{T}^{#gamma,rec} [GeV]");

    TCanvas* c = new TCanvas(Form("c_ratio_%s", outBase.c_str()), Form("c_ratio_%s", outBase.c_str()), 650, 600);
    frame->Draw("axis");
    const double x1 = frame->GetXaxis()->GetXmin();
    const double x2 = frame->GetXaxis()->GetXmax();
    TLine* l1 = new TLine(x1, 1.0, x2, 1.0);
    l1->SetLineStyle(2);
    l1->SetLineWidth(2);
    l1->SetLineColor(kGray+2);
    l1->Draw("same");
    hratio->Draw("same E1X0");

    myText(0.18, 0.90, 1, strleg1.c_str(), 0.04);
    myText(0.18, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.18, 0.80, 1, strleg3.c_str(), 0.04);
    myText(0.18, 0.75, 1, plotTitle.c_str(), 0.04);
    myText(0.18, 0.70, 1, "Ratio: (#eta < 0)/(#eta > 0)", 0.04);

    gSystem->mkdir("figures", kTRUE);
    c->SaveAs(Form("figures/%s_ratio.pdf", outBase.c_str()));

    delete frame;
    delete hratio;
    delete hr;
    delete hl;
  }

  void drawRatioABCDOverlay(
      const std::vector<TH1*>& ratios,
      const std::vector<std::string>& labels,
      const std::string& outBase,
      const std::string& plotTitle,
      const std::string& ytitle,
      double ymin,
      double ymax)
  {
    if(ratios.empty()) return;
    if(ratios.size() != labels.size()) return;

    TH1* href = ratios.front();
    if(!href) return;

    TH1F* frame = makeFrameFrom(href, Form("frame_ratio_abcd_%s", outBase.c_str()), ytitle);
    frame->GetXaxis()->SetTitle(href->GetXaxis()->GetTitle());
    frame->GetXaxis()->SetRangeUser(href->GetXaxis()->GetXmin(), href->GetXaxis()->GetXmax());
    frame->GetYaxis()->SetRangeUser(ymin, ymax);

    TCanvas* c = new TCanvas(Form("c_ratio_abcd_%s", outBase.c_str()), Form("c_ratio_abcd_%s", outBase.c_str()), 650, 600);
    frame->Draw("axis");
    const double x1 = frame->GetXaxis()->GetXmin();
    const double x2 = frame->GetXaxis()->GetXmax();
    TLine* l1 = new TLine(x1, 1.0, x2, 1.0);
    l1->SetLineStyle(2);
    l1->SetLineWidth(2);
    l1->SetLineColor(kGray+2);
    l1->Draw("same");
    for(auto* h : ratios)
    {
      if(!h) continue;
      h->Draw("same E1X0");
    }

    myText(0.18, 0.90, 1, strleg1.c_str(), 0.04);
    myText(0.18, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.18, 0.80, 1, strleg3.c_str(), 0.04);
    myText(0.18, 0.75, 1, plotTitle.c_str(), 0.04);
    myText(0.18, 0.70, 1, "Ratio: (#eta < 0)/(#eta > 0)", 0.04);

    // Legend-style markers
    double y0 = 0.86;
    double dy = 0.06;
    double x0 = 0.58;
    for(size_t i = 0; i < ratios.size(); ++i)
    {
      if(!ratios[i]) continue;
      myMarkerLineText(x0, y0 - dy * i, 0,
                       ratios[i]->GetMarkerColor(), ratios[i]->GetMarkerStyle(),
                       ratios[i]->GetLineColor(),   ratios[i]->GetLineStyle(),
                       labels[i].c_str(), 0.045, true);
    }

    gSystem->mkdir("figures", kTRUE);
    c->SaveAs(Form("figures/%s.pdf", outBase.c_str()));

    delete frame;
  }
} // namespace

// Compare (eta<0 vs eta>0) sideband region spectra stored in:
//   efficiencytool/results/data_histo_<tag>_left.root   (eta < 0)  [useMC=false]
//   efficiencytool/results/data_histo_<tag>_right.root  (eta > 0)  [useMC=false]
// or for simulation:
//   efficiencytool/results/MC_efficiency_<tag>_left.root   (eta < 0)  [useMC=true]
//   efficiencytool/results/MC_efficiency_<tag>_right.root  (eta > 0)  [useMC=true]
//
// For each region (A/B/C/D), it produces:
// - overlay plot: left vs right
// - ratio plot: (left/right) = (eta<0)/(eta>0)
void plot_sideband_compare_left_right(const std::string tag = "bdt_none", bool useMC = false)
{
  init_plot();

  const std::string prefix = useMC ? "MC_efficiency_jet" : "data_histo";
  const std::string sampleLabel = useMC ? "MC" : "data";
  const std::string fRightName = Form("/sphenix/user/shuhangli/ppg12/efficiencytool/results/%s_%s_right.root", prefix.c_str(), tag.c_str());
  const std::string fLeftName  = Form("/sphenix/user/shuhangli/ppg12/efficiencytool/results/%s_%s_left.root",  prefix.c_str(), tag.c_str());

  TFile* fR = TFile::Open(fRightName.c_str(), "READ");
  TFile* fL = TFile::Open(fLeftName.c_str(), "READ");
  if(!fR || fR->IsZombie() || !fL || fL->IsZombie())
  {
    std::cout << "[plot_sideband_compare_left_right] ERROR: cannot open input files:" << std::endl;
    std::cout << "  right: " << fRightName << std::endl;
    std::cout << "  left : " << fLeftName << std::endl;
    return;
  }

  // Sideband region cluster spectra (as used in plot_sideband_selection.C)
  struct Item
  {
    std::string hname;
    std::string out;
    std::string title;
    bool logy;
    std::string ytitle;
    double ymin;
    double ymax;
  };

  std::vector<Item> items = {
      {"h_tight_iso_cluster_0",          Form("sbs_eta_%s_%s_A_tightIso",      sampleLabel.c_str(), tag.c_str()), "Region A: tight iso",        true,  "Counts", 0, 0},
      {"h_tight_noniso_cluster_0",       Form("sbs_eta_%s_%s_B_tightNonIso",   sampleLabel.c_str(), tag.c_str()), "Region B: tight non-iso",    true,  "Counts", 0, 0},
      {"h_nontight_iso_cluster_0",       Form("sbs_eta_%s_%s_C_nonTightIso",   sampleLabel.c_str(), tag.c_str()), "Region C: non-tight iso",    true,  "Counts", 0, 0},
      {"h_nontight_noniso_cluster_0",    Form("sbs_eta_%s_%s_D_nonTightNonIso",sampleLabel.c_str(), tag.c_str()), "Region D: non-tight non-iso",true,  "Counts", 0, 0},
  };

  // Collect ratios for an ABCD overlay plot.
  std::vector<TH1*> ratios_abcd;
  std::vector<std::string> ratio_labels;

  for(const auto& it : items)
  {
    TH1* hr = getHist(fR, it.hname);
    TH1* hl = getHist(fL, it.hname);
    if(!hr || !hl)
    {
      std::cout << "[plot_sideband_compare_left_right] Skipping " << it.hname << " (missing in one file)" << std::endl;
      continue;
    }

    // overlay
    drawCompare1D(hr, hl, it.out, it.title, it.logy, it.ytitle, it.ymin, it.ymax);
    // ratio: (eta<0)/(eta>0)
    drawRatio1D(hr, hl, it.out, it.title, "Ratio", 0.5, 1.8);

    // For the combined ABCD ratio plot, build a standalone ratio histogram.
    TH1* hr_cl = (TH1*)hr->Clone(Form("%s_right_clone_ratioABCD", it.out.c_str()));
    TH1* hl_cl = (TH1*)hl->Clone(Form("%s_left_clone_ratioABCD",  it.out.c_str()));
    hr_cl->SetDirectory(nullptr);
    hl_cl->SetDirectory(nullptr);
    ensureSumw2(hr_cl);
    ensureSumw2(hl_cl);

    TH1* hratio = (TH1*)hl_cl->Clone(Form("%s_ratioABCD", it.out.c_str())); // (eta<0)/(eta>0)
    hratio->SetDirectory(nullptr);
    ensureSumw2(hratio);
    hratio->Divide(hr_cl);

    // Style by region (A/B/C/D) so all four can be overlaid.
    if(it.out.find("_A_") != std::string::npos)      styleHist(hratio, kBlack,   20);
    else if(it.out.find("_B_") != std::string::npos) styleHist(hratio, kRed,     24);
    else if(it.out.find("_C_") != std::string::npos) styleHist(hratio, kBlue+1,  21);
    else if(it.out.find("_D_") != std::string::npos) styleHist(hratio, kGreen+2, 25);
    else                                             styleHist(hratio, kBlack,   20);

    ratios_abcd.push_back(hratio);
    ratio_labels.push_back(it.title);

    delete hr_cl;
    delete hl_cl;
  }

  // One canvas: ratio for A/B/C/D overlaid.
  // Use a stable output base (independent of the per-region naming).
  const std::string outABCD = Form("sbs_eta_ratio_ABCD_%s_%s", sampleLabel.c_str(), tag.c_str());
  const std::string abcdTitle = useMC ? "MC: Sideband regions A/B/C/D" : "Data: Sideband regions A/B/C/D";
  drawRatioABCDOverlay(ratios_abcd, ratio_labels, outABCD, abcdTitle, "Ratio", 0.5, 1.8);
  for(auto* h : ratios_abcd) delete h;

  fR->Close();
  fL->Close();
}


