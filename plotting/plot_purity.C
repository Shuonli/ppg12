#include "plotcommon.h"

void plot_purity()
{
    init_plot();

    string savePath = "figures/";

    TFile *fdata = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nom.root");

    TGraphErrors *gpurity = (TGraphErrors *)fdata->Get("gpurity");
    TGraphErrors *gpurity_leak = (TGraphErrors *)fdata->Get("gpurity_leak");

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    frame_et_rec->SetYTitle("Purity");
    frame_et_rec->GetYaxis()->SetRangeUser(0.0, 1.1);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 30);
    frame_et_rec->Draw("axis");

    gpurity->SetMarkerColor(kBlack);
    gpurity->SetMarkerStyle(20);
    gpurity->SetMarkerSize(1.5);
    gpurity->SetLineColor(kBlack);
    gpurity->Draw("P same");

    gpurity_leak->SetMarkerColor(kBlue);
    gpurity_leak->SetMarkerStyle(24);
    gpurity_leak->SetMarkerSize(1.5);
    gpurity_leak->SetLineColor(kBlue);
    gpurity_leak->Draw("P same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);

    myMarkerLineText(0.26, 0.25, 1, kBlack, 20, kBlack, 1, "w/o signal leakage correction", 0.05, true);
    myMarkerLineText(0.26, 0.20, 1, kBlue, 24, kBlue, 1, "w/ signal leakage correction", 0.05, true);

    c1->SaveAs(Form("%s/purity.pdf", savePath.c_str()));

}