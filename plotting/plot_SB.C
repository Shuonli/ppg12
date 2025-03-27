#include "plotcommon.h"

void plot_SB()
{
    init_plot();

    TFile *fin_sig = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_signal.root", "READ");
    TFile *fin_bg = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_jet.root", "READ");

    const float photon20cross = 1.571e+05 * 0.000673448;

    const float jet30cross = 2.505e-9 * 1e12;

    float jet_scale = jet30cross / photon20cross;

    TH2D* h_sig = (TH2D *)fin_sig->Get("h_ET_isoET_eta0");
    TH2D* h_bg = (TH2D *)fin_bg->Get("h_ET_isoET_eta0");

    int rebinx = 16;

    h_sig->RebinX(rebinx);
    h_bg->RebinX(rebinx);

    TH1D* h_sig_proj = h_sig->ProjectionX("h_sig_proj");
    TH1D* h_bg_proj = h_bg->ProjectionX("h_bg_proj");

    h_sig_proj->Sumw2();
    h_bg_proj->Sumw2();

    //scale background
    h_bg_proj->Scale(jet_scale);

    //calculate s/b
    TH1D* h_sb = (TH1D *)h_sig_proj->Clone("h_sb");

    h_sb->Divide(h_bg_proj);

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    frame_et_rec->SetYTitle("S/B");
    frame_et_rec->GetYaxis()->SetRangeUser(0.0, 1.1);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 30);
    frame_et_rec->Draw("axis");

    h_sb->SetMarkerColor(kBlack);
    h_sb->SetMarkerStyle(20);
    h_sb->SetMarkerSize(1.5);
    h_sb->SetLineColor(kBlack);
    h_sb->Draw("P same");
    //h_bg_proj->Draw();

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5, 0.80, 1, strMC.c_str(), 0.04);

    c1->SaveAs("figures/SB.pdf");










}