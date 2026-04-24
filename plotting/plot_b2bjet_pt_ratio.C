#include "plotcommon.h"

// Cross-section ratio plot: b2bjet pt7/pt5 and pt9/pt5
// vs reco photon ET. Uses h_unfold_sub_result from the three Photon_final
// files for the b2bjet pT-threshold scan.

void plot_b2bjet_pt_ratio()
{
    init_plot();

    const char *base = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";

    TFile *f5 = TFile::Open(Form("%s/Photon_final_bdt_b2bjet_pt5.root", base));
    TFile *f7 = TFile::Open(Form("%s/Photon_final_bdt_b2bjet_pt7.root", base));
    TFile *f9 = TFile::Open(Form("%s/Photon_final_bdt_b2bjet_pt9.root", base));

    TH1F *h5 = (TH1F*)f5->Get("h_unfold_sub_result");
    TH1F *h7 = (TH1F*)f7->Get("h_unfold_sub_result");
    TH1F *h9 = (TH1F*)f9->Get("h_unfold_sub_result");

    // deta cancels in ratios; no need to scale
    TH1F *r7  = (TH1F*)h7->Clone("r7");   r7 ->Divide(h5);
    TH1F *r9  = (TH1F*)h9->Clone("r9");   r9 ->Divide(h5);

    TCanvas *c = new TCanvas("c_b2b_pt_ratio", "", 720, 600);
    c->cd();
    c->SetTickx();
    c->SetTicky();

    TH2F *frame = new TH2F("frame", ";#it{E}_{T}^{#gamma} [GeV];d#sigma(b2b #it{p}_{T}^{jet}#geqX) / d#sigma(b2b #it{p}_{T}^{jet}#geq5 GeV)",
                           100, 10, 36, 100, 0, 1.4);
    frame->GetYaxis()->SetTitleOffset(1.2);
    frame->Draw();

    TLine *one = new TLine(10, 1.0, 36, 1.0);
    one->SetLineStyle(2);
    one->SetLineColor(kGray + 1);
    one->Draw("same");

    r7->SetLineColor(kAzure + 2);
    r7->SetMarkerColor(kAzure + 2);
    r7->SetMarkerStyle(20);
    r7->SetMarkerSize(1.2);
    r7->SetLineWidth(2);
    r7->Draw("E1 same");

    r9->SetLineColor(kRed - 4);
    r9->SetMarkerColor(kRed - 4);
    r9->SetMarkerStyle(21);
    r9->SetMarkerSize(1.2);
    r9->SetLineWidth(2);
    r9->Draw("E1 same");

    TLegend *leg = new TLegend(0.22, 0.18, 0.55, 0.34);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(r7, "#it{p}_{T}^{jet} #geq 7 GeV / #geq 5 GeV", "pl");
    leg->AddEntry(r9, "#it{p}_{T}^{jet} #geq 9 GeV / #geq 5 GeV", "pl");
    leg->Draw();

    myText(0.22, 0.87, 1, strleg1.c_str(), 0.05);
    myText(0.22, 0.82, 1, strleg2_1.c_str(), 0.04);
    myText(0.22, 0.78, 1, strleg3.c_str(), 0.04);
    myText(0.22, 0.74, 1, "b2b-jet requirement cross-check", 0.04);

    c->SaveAs("/sphenix/user/shuhangli/ppg12/plotting/figures/b2bjet_pt_ratio.pdf");
    c->SaveAs("/sphenix/user/shuhangli/ppg12/plotting/figures/b2bjet_pt_ratio.png");

    f5->Close(); f7->Close(); f9->Close();
}
