#include "plotcommon.h"

void plot_efficiency()
{

    init_plot();

    string savePath = "figures/";
    TFile *fmc = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_nom.root");
    TFile *fdata = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nom.root");

    TEfficiency *eff_reco = (TEfficiency *)fmc->Get("eff_reco_eta_0");
    TEfficiency *eff_iso = (TEfficiency *)fmc->Get("eff_iso_eta_0");
    TEfficiency *eff_id = (TEfficiency *)fmc->Get("eff_id_eta_0");
    TEfficiency *eff_all = (TEfficiency *)fmc->Get("eff_all_eta_0");

    TGraphAsymmErrors *g_mbd_eff = (TGraphAsymmErrors *)fdata->Get("g_mbd_eff");

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    frame_et_truth->SetYTitle("Reco Efficiency");
    frame_et_truth->GetXaxis()->SetRangeUser(10, 35);
    frame_et_truth->Draw("axis");

    eff_reco->SetMarkerColor(kBlack);
    eff_reco->SetMarkerStyle(20);
    eff_reco->SetLineColor(kBlack);
    eff_reco->Draw("same");
    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5, 0.80, 1, strleg3.c_str(), 0.04);
    myText(0.2, 0.9, 1, strSigMC.c_str(), 0.04);

    c1->SaveAs(Form("%s/eff_reco.pdf", savePath.c_str()));

    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
    frame_et_truth->SetYTitle("Isolation Efficiency");
    frame_et_truth->Draw("axis");

    eff_iso->SetMarkerColor(kBlack);
    eff_iso->SetMarkerStyle(20);
    eff_iso->SetLineColor(kBlack);
    eff_iso->Draw("same");
    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5, 0.80, 1, strleg3.c_str(), 0.04);
    myText(0.2, 0.9, 1, strSigMC.c_str(), 0.04);

    c2->SaveAs(Form("%s/eff_iso.pdf", savePath.c_str()));

    TCanvas *c3 = new TCanvas("c3", "c3", 600, 600);
    frame_et_truth->SetYTitle("Identification Efficiency");
    frame_et_truth->Draw("axis");

    eff_id->SetMarkerColor(kBlack);
    eff_id->SetMarkerStyle(20);
    eff_id->SetLineColor(kBlack);
    eff_id->Draw("same");
    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5, 0.80, 1, strleg3.c_str(), 0.04);
    myText(0.2, 0.9, 1, strSigMC.c_str(), 0.04);


    c3->SaveAs(Form("%s/eff_id.pdf", savePath.c_str()));

    TCanvas *c4 = new TCanvas("c4", "c4", 600, 600);
    frame_et_truth->SetYTitle("Total Efficiency");
    frame_et_truth->Draw("axis");

    eff_all->SetMarkerColor(kBlack);
    eff_all->SetMarkerStyle(20);
    eff_all->SetLineColor(kBlack);
    eff_all->Draw("same");
    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5, 0.80, 1, strleg3.c_str(), 0.04);
    myText(0.2, 0.9, 1, strSigMC.c_str(), 0.04);

    c4->SaveAs(Form("%s/eff_total.pdf", savePath.c_str()));

    TCanvas *c5 = new TCanvas("c5", "c5", 600, 600);
    frame_et_truth->SetYTitle("MBD Efficiency");
    frame_et_truth->Draw("axis");

    g_mbd_eff->SetMarkerColor(kBlack);
    g_mbd_eff->SetMarkerStyle(20);
    g_mbd_eff->SetLineColor(kBlack);
    g_mbd_eff->Draw("same p");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5, 0.80, 1, strleg3.c_str(), 0.04);
    myText(0.2, 0.9, 1, strSigMC.c_str(), 0.04);

    c5->SaveAs(Form("%s/eff_mbd.pdf", savePath.c_str()));
}