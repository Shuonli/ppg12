#include "plotcommon.h"

void plot_purity_selection(const std::string suffix = "bdt_nom")
{
    init_plot();

    string savePath = "figures/";

    bool plotMC_truth = false;

    std::string dataname = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_"+suffix+".root";

    TFile *fdata = new TFile(dataname.c_str());

    std::string mcname = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_"+suffix+"_mc.root";
    TFile *fMC = new TFile(mcname.c_str());

    TGraphErrors *gpurity = (TGraphErrors *)fdata->Get("gpurity");
    TGraphErrors *gpurity_leak = (TGraphErrors *)fdata->Get("gpurity_leak");
    TGraphErrors *grFineConf = (TGraphErrors *)fdata->Get("grFineConf");
    TGraphErrors *grFineConf_leak = (TGraphErrors *)fdata->Get("grFineConf_leak");
    TF1 *f_purity_fit = (TF1 *)fdata->Get("f_purity_fit");
    TF1 *f_purity_leak_fit = (TF1 *)fdata->Get("f_purity_leak_fit");
    TGraphAsymmErrors *g_purity_truth = (TGraphAsymmErrors *)fMC->Get("g_purity_truth");
    TF1 *f_purity_mc_fit = (TF1 *)fMC->Get("f_purity_leak_fit");
    TGraphErrors *grFineConf_mc = (TGraphErrors *)fMC->Get("grFineConf_leak");

    if (f_purity_leak_fit)
    {
        std::cout << "data fit range: " << f_purity_leak_fit->GetXmin() << " to " << f_purity_leak_fit->GetXmax() << std::endl;
    }
    else
    {
        std::cout << "data fit range: f_purity_leak_fit is null" << std::endl;
    }

    if (f_purity_mc_fit)
    {
        std::cout << "mc fit range:   " << f_purity_mc_fit->GetXmin() << " to " << f_purity_mc_fit->GetXmax() << std::endl;
    }
    else
    {
        std::cout << "mc fit range:   f_purity_mc_fit is null" << std::endl;
    }
    
    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    frame_et_rec->SetTitle(";#it{E}_{T}^{#gamma} [GeV];Purity");
    frame_et_rec->GetYaxis()->SetRangeUser(0.0, 1.1);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 35);
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

    
    grFineConf_leak->SetMarkerColor(kAzure + 2);
    grFineConf_leak->SetLineColor(kAzure + 2);
    grFineConf_leak->SetFillColorAlpha(kAzure + 2, 0.2);

    grFineConf_leak->Draw("e3 same");

    f_purity_leak_fit->SetLineColor(kAzure + 2);

    f_purity_leak_fit->Draw("same");
    if(plotMC_truth)
    {
        grFineConf_mc->SetMarkerColor(kRed);
        grFineConf_mc->SetLineColor(kRed);
        grFineConf_mc->SetFillColorAlpha(kRed, 0.2);
        grFineConf_mc->Draw("e3 same");

        f_purity_mc_fit->SetLineColor(kRed);
        f_purity_mc_fit->Draw("same");

        g_purity_truth->SetMarkerColor(kRed);
        g_purity_truth->SetMarkerStyle(24);
        g_purity_truth->SetMarkerSize(1.5);
        g_purity_truth->SetLineColor(kRed);
        g_purity_truth->Draw("P same");
    }

    float xpos(0.2), xpos2(0.915), ypos(0.88), ypos2(0.19), dy(0.056), dy1(0.075), fontsize(0.048);
    myText(xpos,ypos-0*dy,1,strleg1.c_str(),fontsize,0);
    myText(xpos,ypos-1*dy,1,strleg2.c_str(),fontsize,0);
    myText(xpos2,ypos-0*dy,1,strleg3.c_str(),fontsize,1);
    myText(xpos,ypos-2*dy,1,suffix.c_str(),fontsize,0);
    if(plotMC_truth)
      myText(xpos2,ypos-1*dy,1,strIncMC.c_str(),fontsize,1);

    int nEntry = plotMC_truth ? 5 : 3;
    TLegend* l1 = new TLegend(xpos, ypos2, xpos2, ypos2+nEntry*dy1);
    legStyle(l1, 0.14, fontsize);
    if(plotMC_truth)
    {
      l1->AddEntry(g_purity_truth, "Inc. MC purity", "pl");
      l1->AddEntry(grFineConf_mc, "Inc. MC fit w/ 68\% C.L.", "fl");
    }
    l1->AddEntry(gpurity_leak, "w/ sig. leak. corr.", "pl");
    l1->AddEntry(gpurity, "w/o sig. leak. corr.", "pl");
    l1->AddEntry(grFineConf_leak, "fit w/ 68\% C.L.", "fl");
    l1->Draw("same");

    c1->SaveAs(Form("%s/purity_%s.pdf", savePath.c_str(), suffix.c_str()));

}
