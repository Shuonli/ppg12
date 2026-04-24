#include "plotcommon.h"

void plot_purity()
{
    init_plot();

    string savePath = "figures/";

    bool plotMC_truth = true;

    TFile *fdata = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nomtestv3.root");

    TFile *fMC = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nomtestv3_mc.root");

    TGraphErrors *gpurity = (TGraphErrors *)fdata->Get("gpurity");
    TGraphErrors *gpurity_leak = (TGraphErrors *)fdata->Get("gpurity_leak");
    TGraphAsymmErrors *g_purity_truth = (TGraphAsymmErrors *)fMC->Get("g_purity_truth");
    
    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    frame_et_rec->SetTitle(";#it{E}_{T}^{#gamma} [GeV];Purity");
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

    if(plotMC_truth)
    {
        g_purity_truth->SetMarkerColor(kRed);
        g_purity_truth->SetMarkerStyle(20);
        g_purity_truth->SetMarkerSize(1.5);
        g_purity_truth->SetLineColor(kRed);
        g_purity_truth->Draw("P same");
    }

    float xpos(0.2), xpos2(0.915), ypos(0.885), ypos2(0.19), dy(0.054), dy1(0.065), fontsize(0.042);
    myText(xpos,ypos-0*dy,1,strleg1.c_str(),fontsize,0);
    myText(xpos,ypos-1*dy,1,strleg2.c_str(),fontsize,0);
    myText(xpos2,ypos-0*dy,1,strleg3.c_str(),fontsize,1);
    if(plotMC_truth) 
      myText(xpos2,ypos-1*dy,1,strIncMC.c_str(),fontsize,1);

    int nEntry = plotMC_truth ? 3 : 2;
    TLegend* l1 = new TLegend(xpos, ypos2, xpos2, ypos2+nEntry*dy1);
    legStyle(l1, 0.14, 0.05);
    if(plotMC_truth) 
      l1->AddEntry(g_purity_truth, "truth purity", "pl");
    l1->AddEntry(gpurity_leak, "w/ signal leakage correction", "pl");
    l1->AddEntry(gpurity, "w/o signal leakage correction", "pl");
    l1->Draw("same");

    c1->SaveAs(Form("%s/purity.pdf", savePath.c_str()));

}
