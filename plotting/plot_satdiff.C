#include "plotcommon.h"

void plot_satdiff()
{
    const int col[] = {kAzure + 2, kPink + 5, kGreen - 2, kBlack, kRed - 4, kBlack, kBlue - 3, kYellow + 2, kPink - 5, kGreen + 3, kBlue - 3};

    init_plot();

    string savePath = "figures";

    TFile *fdata = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_nom468.root");
    TFile *fdata_new = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_468newcalib.root");

    TH1D *h_tight_iso_cluster = (TH1D *)fdata->Get("h_tight_iso_cluster_0");
    TH1D *h_tight_noniso_cluster = (TH1D *)fdata->Get("h_tight_noniso_cluster_0");
    TH1D *h_nontight_iso_cluster = (TH1D *)fdata->Get("h_nontight_iso_cluster_0");
    TH1D *h_nontight_noniso_cluster = (TH1D *)fdata->Get("h_nontight_noniso_cluster_0");
    TH1D *h_cluster_common_ET = (TH1D *)fdata->Get("h_cluster_common_E");
    h_cluster_common_ET->Rebin(20);

    TH1D *h_tight_iso_cluster_new = (TH1D *)fdata_new->Get("h_tight_iso_cluster_0");
    TH1D *h_tight_noniso_cluster_new = (TH1D *)fdata_new->Get("h_tight_noniso_cluster_0");
    TH1D *h_nontight_iso_cluster_new = (TH1D *)fdata_new->Get("h_nontight_iso_cluster_0");
    TH1D *h_nontight_noniso_cluster_new = (TH1D *)fdata_new->Get("h_nontight_noniso_cluster_0");
    TH1D *h_cluster_common_ET_new = (TH1D *)fdata_new->Get("h_cluster_common_E");
    h_cluster_common_ET_new->Rebin(20);

    TH1D *h_tight_iso_cluster_ratio = (TH1D *)h_tight_iso_cluster_new->Clone("h_tight_iso_cluster_ratio");
    h_tight_iso_cluster_ratio->Divide(h_tight_iso_cluster);
    TH1D *h_tight_noniso_cluster_ratio = (TH1D *)h_tight_noniso_cluster_new->Clone("h_tight_noniso_cluster_ratio");
    h_tight_noniso_cluster_ratio->Divide(h_tight_noniso_cluster);
    TH1D *h_nontight_iso_cluster_ratio = (TH1D *)h_nontight_iso_cluster_new->Clone("h_nontight_iso_cluster_ratio");
    h_nontight_iso_cluster_ratio->Divide(h_nontight_iso_cluster);
    TH1D *h_nontight_noniso_cluster_ratio = (TH1D *)h_nontight_noniso_cluster_new->Clone("h_nontight_noniso_cluster_ratio");
    h_nontight_noniso_cluster_ratio->Divide(h_nontight_noniso_cluster);

    TH1D *h_cluster_common_ET_ratio = (TH1D *)h_cluster_common_ET_new->Clone("h_cluster_common_ET_ratio");
    h_cluster_common_ET_ratio->Divide(h_cluster_common_ET);


    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    frame_et_rec->Draw("axis");
    frame_et_rec->GetXaxis()->SetRangeUser(8, 35);
    frame_et_rec->GetYaxis()->SetRangeUser(0.8, 1.3);
    //frame_et_rec->SetYTitle("ana468/ana462");
    frame_et_rec->SetYTitle("new/old calib");
    h_tight_iso_cluster_ratio->SetLineColor(col[0]);
    h_tight_iso_cluster_ratio->SetLineWidth(4);
    h_tight_iso_cluster_ratio->Draw("same hist");
    h_tight_noniso_cluster_ratio->SetLineColor(col[1]);
    h_tight_noniso_cluster_ratio->SetLineWidth(3);
    h_tight_noniso_cluster_ratio->Draw("same hist");
    h_nontight_iso_cluster_ratio->SetLineColor(col[2]);
    h_nontight_iso_cluster_ratio->SetLineWidth(2);
    h_nontight_iso_cluster_ratio->Draw("same hist");
    h_nontight_noniso_cluster_ratio->SetLineColor(col[3]);
    h_nontight_noniso_cluster_ratio->SetLineWidth(1);
    h_nontight_noniso_cluster_ratio->Draw("same hist");

    lineone->Draw("same");

    float dx = 0.3;
    float dy = 0.3;

    myText(0.5 - dx, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5 - dx, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5 - dx, 0.80, 1, Form("Data   %s", strleg3.c_str()), 0.04);
    myMarkerLineText(0.55 - dx, 0.75 - dy, 0, col[0], 0, col[0], 1, "tight iso", 0.05, true);
    myMarkerLineText(0.55 - dx, 0.70 - dy, 0, col[1], 0, col[1], 1, "tight noniso", 0.05, true);
    myMarkerLineText(0.55 - dx, 0.65 - dy, 0, col[2], 0, col[2], 1, "nontight iso", 0.05, true);
    myMarkerLineText(0.55 - dx, 0.60 - dy, 0, col[3], 0, col[3], 1, "nontight noniso", 0.05, true);

    c1->SaveAs(Form("%s/et_sat.pdf", savePath.c_str())) ;

    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
    frame_et_rec->Draw("axis");
    frame_et_rec->GetXaxis()->SetRangeUser(8, 35);
    frame_et_rec->GetYaxis()->SetRangeUser(0.2, 1.8);

    h_cluster_common_ET_ratio->SetLineColor(kBlack);
    h_cluster_common_ET_ratio->Draw("same hist");

    lineone->Draw("same");

    myText(0.5 - dx, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5 - dx, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5 - dx, 0.80, 1, Form("Data   %s", strleg3.c_str()), 0.04);
    myMarkerLineText(0.55 - dx, 0.75, 0, kBlack, 0, kBlack, 1, "cluster w/npb cut", 0.05, true);


}
