#include "plotcommon.h"

const int col[]={kBlack,kPink+5,kGreen-2,kAzure+7,kRed-4,kBlack,kBlue-3,kYellow+2,kPink-5,kGreen+3,kBlue-3};
const int mkStyle[]={21,20,34,33,25,27,28,24,29,28,22};
const float mkSize[]={1.2,1.2,1.6,1.1,1,1,1,1,1,1,1,1};
void plot_efficiency()
{

    init_plot();

    // string savePath = "test/";
    string savePath = "figures/";
    TFile *fmc = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_nom.root");
    TFile *fdata = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nom.root");

    TEfficiency *eff_reco = (TEfficiency *)fmc->Get("eff_reco_eta_0");
    TEfficiency *eff_iso = (TEfficiency *)fmc->Get("eff_iso_eta_0");
    TEfficiency *eff_id = (TEfficiency *)fmc->Get("eff_id_eta_0");
    TEfficiency *eff_all = (TEfficiency *)fmc->Get("eff_all_eta_0");

    TGraphAsymmErrors *g_mbd_eff = (TGraphAsymmErrors *)fdata->Get("g_mbd_eff");

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    frame_et_truth->SetYTitle("Reconstruction Efficiency");
    frame_et_truth->GetXaxis()->SetRangeUser(10, 26);
    frame_et_truth->GetYaxis()->SetRangeUser(0.45, 1.1);
    frame_et_truth->Draw("axis");

    eff_reco->SetMarkerColor(kBlack);
    eff_reco->SetMarkerStyle(20);
    eff_reco->SetLineColor(kBlack);
    eff_reco->Draw("same");

    float xpos(0.2), xpos2(0.915), ypos(0.885), ypos2(0.22), dy(0.054), dy1(0.08), fontsize(0.046), fontsize1(0.048);
    myText(xpos2,ypos-0*dy,1,strMC.c_str(),fontsize,1);
    myText(xpos2,ypos-1*dy,1,strleg3.c_str(),fontsize,1);
    myText(xpos,ypos-0*dy,1,strleg1.c_str(),fontsize,0);
    myText(xpos,ypos-1*dy,1,strleg2.c_str(),fontsize,0);

    c1->SaveAs(Form("%s/eff_reco.pdf", savePath.c_str()));

    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
    frame_et_truth->SetYTitle("Isolation Efficiency");
    frame_et_truth->Draw("axis");

    eff_iso->SetMarkerColor(kBlack);
    eff_iso->SetMarkerStyle(20);
    eff_iso->SetLineColor(kBlack);
    eff_iso->Draw("same");
    myText(xpos2,ypos-0*dy,1,strMC.c_str(),fontsize,1);
    myText(xpos2,ypos-1*dy,1,strleg3.c_str(),fontsize,1);
    myText(xpos,ypos-0*dy,1,strleg1.c_str(),fontsize,0);
    myText(xpos,ypos-1*dy,1,strleg2.c_str(),fontsize,0);

    c2->SaveAs(Form("%s/eff_iso.pdf", savePath.c_str()));

    TCanvas *c3 = new TCanvas("c3", "c3", 600, 600);
    frame_et_truth->SetYTitle("Identification Efficiency");
    frame_et_truth->Draw("axis");

    eff_id->SetMarkerColor(kBlack);
    eff_id->SetMarkerStyle(20);
    eff_id->SetLineColor(kBlack);
    eff_id->Draw("same");
    myText(xpos2,ypos-0*dy,1,strMC.c_str(),fontsize,1);
    myText(xpos2,ypos-1*dy,1,strleg3.c_str(),fontsize,1);
    myText(xpos,ypos-0*dy,1,strleg1.c_str(),fontsize,0);
    myText(xpos,ypos-1*dy,1,strleg2.c_str(),fontsize,0);

    c3->SaveAs(Form("%s/eff_id.pdf", savePath.c_str()));

    TCanvas *c4 = new TCanvas("c4", "c4", 600, 600);
    frame_et_truth->SetYTitle("Total Efficiency");
    frame_et_truth->Draw("axis");

    eff_all->SetMarkerColor(kBlack);
    eff_all->SetMarkerStyle(20);
    eff_all->SetLineColor(kBlack);
    eff_all->Draw("same");
    myText(xpos2,ypos-0*dy,1,strMC.c_str(),fontsize,1);
    myText(xpos2,ypos-1*dy,1,strleg3.c_str(),fontsize,1);
    myText(xpos,ypos-0*dy,1,strleg1.c_str(),fontsize,0);
    myText(xpos,ypos-1*dy,1,strleg2.c_str(),fontsize,0);

    c4->SaveAs(Form("%s/eff_total.pdf", savePath.c_str()));

    TCanvas *c5 = new TCanvas("c5", "c5", 600, 600);
    frame_et_truth->SetYTitle("MBD Trigger Efficiency");
    frame_et_truth->GetYaxis()->SetRangeUser(0.45, 1.1);
    frame_et_truth->Draw("axis");

    g_mbd_eff->SetMarkerColor(kBlack);
    g_mbd_eff->SetMarkerStyle(20);
    g_mbd_eff->SetLineColor(kBlack);
    g_mbd_eff->Draw("same p");

    myText(xpos2,ypos-0*dy,1,strMC.c_str(),fontsize,1);
    myText(xpos2,ypos-1*dy,1,strleg3.c_str(),fontsize,1);
    myText(xpos,ypos-0*dy,1,strleg1.c_str(),fontsize,0);
    myText(xpos,ypos-1*dy,1,strleg2.c_str(),fontsize,0);

    c5->SaveAs(Form("%s/eff_mbd.pdf", savePath.c_str()));

    TGraphAsymmErrors *g_reco_id_product = (TGraphAsymmErrors *)g_mbd_eff->Clone("g_reco_id_product");

    // loop over bins
    for (int i = 0; i < eff_reco->GetTotalHistogram()->GetNbinsX(); i++)
    {
        double reco = eff_reco->GetEfficiency(i+1);
        double id = eff_id->GetEfficiency(i+1);
        double reco_err = eff_reco->GetEfficiencyErrorLow(i+1);
        double id_err = eff_id->GetEfficiencyErrorLow(i+1);

        double product = 0;
        double product_err = 0;

        // double bin_

        product = reco * id;
        product_err = product * sqrt(pow(reco_err / reco, 2) + pow(id_err / id, 2));
        float x = g_mbd_eff->GetX()[i];
        g_reco_id_product->SetPoint(i, x, product);
        // g_reco_id_ratio->SetPointError(i, 0, 0, ratio_err, ratio_err);
    }

    TCanvas *c6 = new TCanvas("c6", "c6", 600, 600);
    frame_et_truth->SetYTitle("Efficiency");
    frame_et_truth->GetYaxis()->SetRangeUser(0.35, 1.1);
    frame_et_truth->Draw("axis");

    eff_reco->SetMarkerColor(col[0]);
    eff_reco->SetMarkerStyle(mkStyle[0]);
    eff_reco->SetMarkerSize(mkSize[0]);
    eff_reco->SetLineColor(col[0]);
    eff_reco->Draw("same");

    g_reco_id_product->SetMarkerColor(col[1]);
    g_reco_id_product->SetMarkerStyle(mkStyle[1]);
    g_reco_id_product->SetMarkerSize(mkSize[1]);
    g_reco_id_product->SetLineColor(col[1]);
    g_reco_id_product->Draw("same p");

    eff_all->SetMarkerColor(col[2]);
    eff_all->SetMarkerStyle(mkStyle[2]);
    eff_all->SetMarkerSize(mkSize[2]);
    eff_all->SetLineColor(col[2]);
    eff_all->SetLineWidth(2);
    eff_all->Draw("same");

    // float xpos(0.2), xpos2(0.915), ypos(0.885), ypos2(0.225), dy(0.054), dy1(0.065), fontsize(0.042);
    myText(xpos,ypos-0*dy,1,strleg1.c_str(),fontsize1,0);
    myText(xpos,ypos-1*dy,1,strleg2.c_str(),fontsize,0);
    myText(xpos2,ypos-0*dy,1,strMC.c_str(),fontsize,1);
    myText(xpos2,ypos-1*dy,1,strleg3.c_str(),fontsize,1);

    TLegend* l1 = new TLegend(0.45, ypos2, xpos2, ypos2+3*dy1);
    legStyle(l1, 0.20, 0.06);
    l1->AddEntry(eff_reco, "#varepsilon_{reco}", "pl");
    l1->AddEntry(g_reco_id_product, "#varepsilon_{reco}#times#varepsilon_{ID}", "pl");
    l1->AddEntry(eff_all, "#varepsilon_{reco}#times#varepsilon_{ID}#times#varepsilon_{iso}", "pl");
    l1->Draw("same");

    c6->SaveAs(Form("%s/eff_photon.pdf", savePath.c_str()));
}