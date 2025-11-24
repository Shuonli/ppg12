#include "plotcommon.h"

const int col[] = {kBlack, kPink + 5, kGreen - 2, kAzure + 7, kRed - 4, kBlack, kBlue - 3, kYellow + 2, kPink - 5, kGreen + 3, kBlue - 3};
const int mkStyle[] = {21, 20, 34, 33, 25, 27, 28, 24, 29, 28, 22};
const float mkSize[] = {1.2, 1.2, 1.6, 1.1, 1, 1, 1, 1, 1, 1, 1, 1};
void plot_xjg()
{

    init_plot();

    float minpt = 10;
    float maxpt = 35;
    int rebin = 5;

    // string savePath = "test/";
    string savePath = "figures/";
    TFile *fmc = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_nom.root");
    TFile *fmc_inclusive = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_nom.root");
    TFile *fdata = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_nom.root");

    TH2D *h_xjg_pt_signalmc = (TH2D *)fmc->Get("h_tight_iso_xjgamma_signal_0");
    TH2D *h_xjg_pt_inclusive = (TH2D *)fmc_inclusive->Get("h_tight_iso_xjgamma_0");
    TH2D *h_xjg_pt_data = (TH2D *)fdata->Get("h_tight_iso_xjgamma_0");

    TCanvas *c1 = new TCanvas("c1", "c1", 900, 700);
    // set right margin
    c1->SetRightMargin(0.12);
    h_xjg_pt_signalmc->SetXTitle("x_{J#gamma}");
    h_xjg_pt_signalmc->SetYTitle("p_{T}^{#gamma} [GeV]");
    h_xjg_pt_signalmc->Draw("colz");

    float xpos(0.2), xpos2(0.835), ypos(0.885), ypos2(0.65), dy(0.054), dy1(0.08), fontsize(0.046), fontsize1(0.048);
    myText(xpos2, ypos - 0 * dy, 1, strSigMC.c_str(), fontsize, 1);
    myText(xpos2, ypos - 1 * dy, 1, strlegphotonjet.c_str(), fontsize, 1);
    myText(xpos2, ypos - 2 * dy, 1, strdijet.c_str(), fontsize, 1);
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);

    c1->SaveAs(Form("%s/xjg_signalmc.pdf", savePath.c_str()));

    TCanvas *c2 = new TCanvas("c2", "c2", 900, 700);
    // set right margin
    c2->SetRightMargin(0.12);
    h_xjg_pt_inclusive->SetXTitle("x_{J#gamma}");
    h_xjg_pt_inclusive->SetYTitle("p_{T}^{#gamma} [GeV]");
    h_xjg_pt_inclusive->Draw("colz");

    myText(xpos2, ypos - 0 * dy, 1, strIncMC.c_str(), fontsize, 1);
    myText(xpos2, ypos - 1 * dy, 1, strlegphotonjet.c_str(), fontsize, 1);
    myText(xpos2, ypos - 2 * dy, 1, strdijet.c_str(), fontsize, 1);
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);

    c2->SaveAs(Form("%s/xjg_incmc.pdf", savePath.c_str()));

    TCanvas *c3 = new TCanvas("c3", "c3", 900, 700);
    // set right margin
    c3->SetRightMargin(0.12);
    h_xjg_pt_data->SetXTitle("x_{J#gamma}");
    h_xjg_pt_data->SetYTitle("p_{T}^{#gamma} [GeV]");
    h_xjg_pt_data->Draw("colz");
    myText(xpos2, ypos - 0 * dy, 1, "Data", fontsize, 1);
    myText(xpos2, ypos - 1 * dy, 1, strlegphotonjet.c_str(), fontsize, 1);
    myText(xpos2, ypos - 2 * dy, 1, strdijet.c_str(), fontsize, 1);
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
    c3->SaveAs(Form("%s/xjg_data.pdf", savePath.c_str()));

    h_xjg_pt_data->GetYaxis()->SetRangeUser(minpt, maxpt);
    TH1D *h_xjg_data = h_xjg_pt_data->ProjectionX("h_xjg_data", 1, h_xjg_pt_data->GetNbinsY());
    h_xjg_data->Rebin(rebin);
    h_xjg_data->SetLineColor(col[0]);
    h_xjg_data->SetMarkerColor(col[0]);
    h_xjg_data->SetMarkerStyle(mkStyle[0]);
    h_xjg_data->SetMarkerSize(mkSize[0]);

    h_xjg_pt_signalmc->GetYaxis()->SetRangeUser(minpt, maxpt);
    TH1D *h_xjg_signalmc = h_xjg_pt_signalmc->ProjectionX("h_xjg_signalmc", 1, h_xjg_pt_signalmc->GetNbinsY());
    h_xjg_signalmc->Rebin(rebin);
    h_xjg_signalmc->SetLineColor(col[1]);
    h_xjg_signalmc->SetMarkerColor(col[1]);
    h_xjg_signalmc->SetMarkerStyle(mkStyle[1]);
    h_xjg_signalmc->SetMarkerSize(mkSize[1]);

    h_xjg_pt_inclusive->GetYaxis()->SetRangeUser(minpt, maxpt);
    TH1D *h_xjg_inclusive = h_xjg_pt_inclusive->ProjectionX("h_xjg_incmc", 1, h_xjg_pt_inclusive->GetNbinsY());
    h_xjg_inclusive->Rebin(rebin);
    h_xjg_inclusive->SetLineColor(col[2]);
    h_xjg_inclusive->SetMarkerColor(col[2]);
    h_xjg_inclusive->SetMarkerStyle(mkStyle[2]);
    h_xjg_inclusive->SetMarkerSize(mkSize[2]);

    // normalize by integral
    h_xjg_data->Scale(1.0 / h_xjg_data->Integral()/rebin);
    h_xjg_signalmc->Scale(1.0 / h_xjg_signalmc->Integral()/rebin);
    h_xjg_inclusive->Scale(1.0 / h_xjg_inclusive->Integral()/rebin);

    TCanvas *c4 = new TCanvas("c4", "c4", 700, 700);
    //logy
    //c4->SetLogy();
    // c4->SetRightMargin(0.12);
    frame_isoET->GetXaxis()->SetRangeUser(0.6, 1.5);
    frame_isoET->GetYaxis()->SetRangeUser(0.000, 0.035);
    frame_isoET->SetXTitle("x_{J#gamma}");
    frame_isoET->SetYTitle("Arbitrary Units");
    frame_isoET->Draw("axis");

    h_xjg_data->Draw("same pl");
    h_xjg_signalmc->Draw("same pl");
    h_xjg_inclusive->Draw("same pl");

    myText(xpos2+0.08, ypos - 0 * dy, 1, strlegphotonjet.c_str(), fontsize, 1);
    myText(xpos2+0.08, ypos - 1 * dy, 1, strdijet.c_str(), fontsize, 1);
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);

    TLegend* l1 = new TLegend(0.5, ypos2, xpos2, ypos2+3*dy1);
    legStyle(l1, 0.20, 0.04);
    l1->AddEntry(h_xjg_data, "Data", "pl");
    l1->AddEntry(h_xjg_signalmc, strSigMC.c_str(), "pl");
    l1->AddEntry(h_xjg_inclusive, strIncMC.c_str(), "pl");
    l1->Draw("same");

    c4->SaveAs(Form("%s/xjg_data_signalmc_incmc.pdf", savePath.c_str()));



    /*
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
    */
}