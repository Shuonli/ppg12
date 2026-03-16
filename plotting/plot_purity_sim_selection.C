#include "plotcommon.h"

void plot_purity_sim_selection(const std::string suffix = "bdt_nom")
{

    init_plot();
    
    string savePath = "figures/";

    std::string dataname = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_" + suffix + "_mc.root";

    TFile *fdata = new TFile(dataname.c_str());

    TGraphErrors *gpurity = (TGraphErrors *)fdata->Get("gpurity");
    TGraphErrors *gpurity_leak = (TGraphErrors *)fdata->Get("gpurity_leak");
    TGraphAsymmErrors *g_purity_truth = (TGraphAsymmErrors *)fdata->Get("g_purity_truth");
    TF1 *f_purity_leak_fit = (TF1 *)fdata->Get("f_purity_leak_fit");
    TGraphErrors *grFineConf_leak = (TGraphErrors *)fdata->Get("grFineConf_leak");

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    frame_et_rec->SetYTitle("Purity");
    frame_et_rec->GetYaxis()->SetRangeUser(0.0, 1.2);
    frame_et_rec->GetXaxis()->SetRangeUser(8, 35);
    frame_et_rec->Draw("axis");

    if (grFineConf_leak)
    {
        grFineConf_leak->SetMarkerColor(kAzure + 2);
        grFineConf_leak->SetLineColor(kAzure + 2);
        grFineConf_leak->SetFillColorAlpha(kAzure + 2, 0.2);
        grFineConf_leak->Draw("e3 same");
    }

    if (f_purity_leak_fit)
    {
        f_purity_leak_fit->SetLineColor(kAzure + 2);
        f_purity_leak_fit->SetLineWidth(2);
        f_purity_leak_fit->Draw("same");
    }

    gpurity->SetMarkerColor(kBlack);
    gpurity->SetMarkerStyle(20);
    gpurity->SetMarkerSize(1.5);
    gpurity->SetLineColor(kBlack);
    gpurity->Draw("P same");

    gpurity_leak->SetMarkerColor(kBlue);
    gpurity_leak->SetMarkerStyle(20);
    gpurity_leak->SetMarkerSize(1.5);
    gpurity_leak->SetLineColor(kBlue);
    gpurity_leak->Draw("P same");

    g_purity_truth->SetMarkerColor(kRed);
    g_purity_truth->SetMarkerStyle(20);
    g_purity_truth->SetMarkerSize(1.5);
    g_purity_truth->SetLineColor(kRed);
    g_purity_truth->Draw("P same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5, 0.80, 1, strMC.c_str(), 0.04);
    myText(0.5, 0.75, 1, suffix.c_str(), 0.04);

    myMarkerLineText(0.30, 0.25, 1, kBlack, 20, kBlack, 1, "w/o signal leakage correction", 0.05, true);
    myMarkerLineText(0.30, 0.20, 1, kBlue, 20, kBlue, 1, "w/ signal leakage correction", 0.05, true);
    myMarkerLineText(0.30, 0.30, 1, kRed, 20, kRed, 1, "truth", 0.05, true);

    c1->SaveAs(Form("%s/purity_sim_%s.pdf", savePath.c_str(), suffix.c_str()));

}
