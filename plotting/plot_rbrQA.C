#include "plotcommon.h"

void plot_rbrQA()
{
    init_plot();

    string savePath = "figures/";

    float minrun = 48200;
    float maxrun = 48280;

    TFile *fdata = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/rbrQA.root");
    TFile *fdata_nc = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/rbrQA_nc.root");

    TGraphErrors *gr_common = (TGraphErrors *)fdata->Get("gr_common");
    TGraphErrors *gr_tight_iso = (TGraphErrors *)fdata->Get("gr_tight_iso");
    TGraphErrors *gr_tight_noniso = (TGraphErrors *)fdata->Get("gr_tight_noniso");
    TGraphErrors *gr_nontight_iso = (TGraphErrors *)fdata->Get("gr_nontight_iso");
    //TGraphErrors *gr_nontight_noniso = (TGraphErrors *)fdata->Get("gr_nontight_noniso");

    TGraphErrors *gr_common_nc = (TGraphErrors *)fdata_nc->Get("gr_common")->Clone("gr_common_nc");
    TGraphErrors *gr_tight_iso_nc = (TGraphErrors *)fdata_nc->Get("gr_tight_iso")->Clone("gr_tight_iso_nc");
    TGraphErrors *gr_tight_noniso_nc = (TGraphErrors *)fdata_nc->Get("gr_tight_noniso")->Clone("gr_tight_noniso_nc");
    TGraphErrors *gr_nontight_iso_nc = (TGraphErrors *)fdata_nc->Get("gr_nontight_iso")->Clone("gr_nontight_iso_nc");
    //TGraphErrors *gr_nontight_noniso_nc = (TGraphErrors *)fdata_nc->Get("gr_nontight_noniso")->Clone("gr_nontight_noniso_nc");

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);

    TH1D *frame = new TH1D("frame", ";Run Number;Common Clusters", 1, minrun, maxrun);
    frame->GetYaxis()->SetRangeUser(60000, 140000);
    frame->GetXaxis()->SetNdivisions(505);

    frame->Draw("axis");

    gr_common->SetMarkerColor(kBlack);
    gr_common->SetMarkerStyle(20);
    gr_common->SetLineColor(kBlack);
    gr_common->Draw("P same");

    gr_common_nc->SetMarkerColor(kRed);
    gr_common_nc->SetMarkerStyle(20);
    gr_common_nc->SetLineColor(kRed);
    gr_common_nc->Draw("P same");

    //gr_common->Draw("P same");

    myMarkerLineText(0.25, 0.9, 1, kBlack, 20, kBlack, 1, "w/ pileup correction", 0.05, true);
    myMarkerLineText(0.25, 0.85, 1, kRed, 20, kRed, 1, "w/o pileup correction", 0.05, true);


    TH1D *h_common_rbr = new TH1D("h_common_rbr", "", 50, 2e4, 1.8e5);
    TH1D *h_common_rbr_nc = new TH1D("h_common_rbr_nc", "", 50, 2e4, 1.8e5);
    h_common_rbr->SetXTitle("Common Clusters/lumi");

    //loop over all the points
    for (int i = 0; i < gr_common->GetN(); i++)
    {
        double x, y;
        gr_common->GetPoint(i, x, y);
        double yerr = gr_common->GetErrorY(i);
        h_common_rbr->Fill(y);
    }
    
    for (int i = 0; i < gr_common_nc->GetN(); i++)
    {
        double x, y;
        gr_common_nc->GetPoint(i, x, y);
        double yerr = gr_common_nc->GetErrorY(i);
        h_common_rbr_nc->Fill(y);
    }
    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
    //print histogram mean and rms
    cout << "h_common_rbr mean: " << h_common_rbr->GetMean() << " rms: " << h_common_rbr->GetRMS()/h_common_rbr->GetMean() << endl;
    cout << "h_common_rbr_nc mean: " << h_common_rbr_nc->GetMean() << " rms: " << h_common_rbr_nc->GetRMS()/h_common_rbr_nc->GetMean() << endl;
    //fit each with an gaussian
    float xlow = h_common_rbr->GetMean() - 1. * h_common_rbr->GetRMS();
    float xhigh = h_common_rbr->GetMean() + 1.5 * h_common_rbr->GetRMS();
    float xlow_nc = h_common_rbr_nc->GetMean() - 1. * h_common_rbr_nc->GetRMS();
    float xhigh_nc = h_common_rbr_nc->GetMean() + 1.5 * h_common_rbr_nc->GetRMS();
    TF1 *f1 = new TF1("f1", "gaus", xlow, xhigh);
    h_common_rbr->Fit("f1", "REM", "", xlow, xhigh);
    TF1 *f1_nc = new TF1("f1_nc", "gaus", xlow_nc, xhigh_nc);
    h_common_rbr_nc->Fit("f1_nc", "REM", "", xlow_nc, xhigh_nc);
    //print the fitted mean and rms
    cout << "fitted mean: " << f1->GetParameter(1) << " rms: " << f1->GetParameter(2)/f1->GetParameter(1) << endl;
    cout << "fitted mean_nc: " << f1_nc->GetParameter(1) << " rms_nc: " << f1_nc->GetParameter(2)/f1_nc->GetParameter(1) << endl;
    
    h_common_rbr->Draw();
    h_common_rbr_nc->SetLineColor(kRed);
    h_common_rbr_nc->Draw("same");

}