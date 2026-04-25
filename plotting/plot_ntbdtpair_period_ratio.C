#include "plotcommon.h"

void plot_ntbdtpair_period_ratio()
{
  init_plot();

  const std::string basepath =
      "/sphenix/user/shuhangli/ppg12/efficiencytool/results/";
  const std::string figDir =
      "/sphenix/user/shuhangli/ppg12/plotting/figures";
  gSystem->mkdir(figDir.c_str(), true);

  std::vector<std::string> suffixes = {
      "t80_70_h80_70_l70_40", "t80_75_h80_75_l70_50", "t80_75_h80_75_l70_40",
      "t80_70_h80_70_l60_40", "t80_75_h80_75_l60_50", "t80_75_h80_75_l60_40",
      "t80_70_h70_70_l60_40", "t80_75_h70_75_l60_50", "t80_75_h70_75_l60_40",
  };

  auto getHist = [&](const std::string &fname) -> TH1F *
  {
    TFile *f = TFile::Open(fname.c_str(), "READ");
    if (!f || f->IsZombie())
    {
      std::cerr << "Warning: cannot open " << fname << std::endl;
      if (f) { f->Close(); delete f; }
      return nullptr;
    }
    TH1F *h = (TH1F *) f->Get("h_unfold_sub_result");
    if (!h)
    {
      std::cerr << "Warning: h_unfold_sub_result missing in " << fname << std::endl;
      return nullptr;
    }
    TH1F *hc = (TH1F *) h->Clone(Form("%s_clone", fname.c_str()));
    hc->SetDirectory(0);
    return hc;
  };

  for (const auto &sfx : suffixes)
  {
    TH1F *h0   = getHist(basepath + "Photon_final_bdt_ntbdtpair_" + sfx + "_0mrad.root");
    TH1F *h1p5 = getHist(basepath + "Photon_final_bdt_ntbdtpair_" + sfx + "_1p5mrad.root");
    if (!h0 || !h1p5) continue;

    // ROOT's default TH1::Divide treats numerator and denominator as statistically
    // independent — standard propagation: (sigma/h)^2 = (sigma_num/num)^2 + (sigma_den/den)^2
    TH1F *hratio = (TH1F *) h0->Clone(Form("h_ratio_%s", sfx.c_str()));
    hratio->SetDirectory(0);
    hratio->Divide(h1p5);

    TCanvas *c = new TCanvas(Form("c_ratio_%s", sfx.c_str()), "", 700, 600);
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.06);

    TH1F *frame = (TH1F *) frame_et_truth->Clone(Form("frame_%s", sfx.c_str()));
    frame->SetTitle("");
    frame->SetYTitle("#sigma(0 mrad) / #sigma(1.5 mrad)");
    frame->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");
    frame->GetYaxis()->SetRangeUser(0.0, 2.5);
    frame->GetXaxis()->SetRangeUser(pTmin, pTmax);
    frame->GetYaxis()->SetNdivisions(506);
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->Draw("axis");

    lineone->SetLineColor(kGray + 2);
    lineone->SetLineStyle(2);
    lineone->Draw("L same");

    hratio->SetMarkerStyle(20);
    hratio->SetMarkerColor(kBlack);
    hratio->SetLineColor(kBlack);
    hratio->SetMarkerSize(1.3);
    hratio->SetLineWidth(2);
    hratio->Draw("PE1 same");

    const float xpos = 0.22, ypos = 0.88, dy = 0.055, fs = 0.038;
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fs, 0);
    myText(xpos, ypos - 1 * dy, 1,
           "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV", fs, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fs, 0);
    myText(xpos, ypos - 3 * dy, 1,
           Form("ntbdtpair %s", sfx.c_str()), fs, 0);

    c->SaveAs(Form("%s/xsec_ratio_0over1p5mrad_bdt_ntbdtpair_%s.pdf",
                   figDir.c_str(), sfx.c_str()));
    delete c;
  }
}
