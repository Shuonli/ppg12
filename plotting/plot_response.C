#include "plotcommon.h"

void plot_response()
{

    init_plot();

    string savePath = "figures";

    string matrix = "response matrix";

    TFile *fdata = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nom.root");
    TFile *fdata_nr = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nr.root");

    TH2D* h_response_leak_reweighted = (TH2D*) fdata->Get("h_response_full_0");

    TH2D* h_response = (TH2D*) fdata_nr->Get("h_response_full_0");


    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    //frame_response->Draw("axis");
    h_response_leak_reweighted->SetXTitle(frame_response->GetXaxis()->GetTitle());
    h_response_leak_reweighted->SetYTitle(frame_response->GetYaxis()->GetTitle());
    h_response_leak_reweighted->Draw("colz");

    //logz
    c1->SetLogz();

     gPad->SetRightMargin(0.15);

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.3, 0.75, 1, Form("%s  reweighted", matrix.c_str()), 0.04);

    c1->SaveAs(Form("%s/response_reweighted.pdf", savePath.c_str()));


    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);

    //frame_response->Draw("axis");

    h_response->SetXTitle(frame_response->GetXaxis()->GetTitle());
    h_response->SetYTitle(frame_response->GetYaxis()->GetTitle());
    h_response->Draw("colz");

    //logz
    c2->SetLogz();

     gPad->SetRightMargin(0.15);

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5, 0.75, 1, matrix.c_str(), 0.04);

    c2->SaveAs(Form("%s/response.pdf", savePath.c_str()));


}
