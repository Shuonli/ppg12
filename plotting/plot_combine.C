#include "plotcommon.h"

void plot_combine()
{
    init_plot();

    TFile *f_photon5 = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_photon5_mbdeffup.root", "READ");
    TFile *f_photon10 = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_photon10_mbdeffup.root", "READ");
    TFile *f_photon20 = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_photon20_mbdeffup.root", "READ");

    int rebinx = 10;

    TH1F *h_max_photon_pT_5 = (TH1F *)f_photon5->Get("h_max_photon_pT");
    h_max_photon_pT_5->Rebin(rebinx);
    TH1F *h_max_photon_pT_10 = (TH1F *)f_photon10->Get("h_max_photon_pT");
    h_max_photon_pT_10->Rebin(rebinx);
    h_max_photon_pT_10->Scale(1.00);
    TH1F *h_max_photon_pT_20 = (TH1F *)f_photon20->Get("h_max_photon_pT");
    h_max_photon_pT_20->Rebin(rebinx);

    TH1F *h_max_photon_pT_sum = (TH1F *)h_max_photon_pT_5->Clone("h_max_photon_pT_sum");
    h_max_photon_pT_sum->Add(h_max_photon_pT_10);
    h_max_photon_pT_sum->Add(h_max_photon_pT_20);

    std::vector<int> colors = {kPink + 5, kGreen - 2, kAzure + 7, kRed - 4, kBlue - 3, kYellow + 2, kPink - 5, kGreen + 3, kBlue - 3, kBlack};

    TH1F *h_20_over_10 = (TH1F *)h_max_photon_pT_20->Clone("h_20_over_10");
    TH1F *h_10_over_05 = (TH1F *)h_max_photon_pT_10->Clone("h_10_over_05");

    h_20_over_10->Divide(h_max_photon_pT_10);
    h_10_over_05->Divide(h_max_photon_pT_5);
    
    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    h_20_over_10->GetXaxis()->SetRangeUser(20, 40);
    h_20_over_10->GetYaxis()->SetRangeUser(0.6, 1.4);
    h_20_over_10->SetYTitle("photon20/photon10");
    h_20_over_10->SetXTitle("Leading #it{E}_{T}^{#gamma} [GeV]");
    h_20_over_10->Draw();
    lineone->Draw("L");
    // h_10_over_05->Draw();
    // h_histo_over_fit->Draw();
    // h_max_photon_pT_sum->Draw();
    // f1->SetLineColor(kRed);
    // f1->Draw("same");

    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
    h_10_over_05->GetXaxis()->SetRangeUser(10, 30);
    h_10_over_05->GetYaxis()->SetRangeUser(0.6, 1.4);
    h_10_over_05->SetYTitle("photon10/photon5");
    h_10_over_05->SetXTitle("Leading #it{E}_{T}^{#gamma} [GeV]");
    h_10_over_05->Draw();
    lineone->Draw("L"); 

    c1->SaveAs("figures/20_over_10.pdf");
    c2->SaveAs("figures/10_over_05.pdf");

    // cut the stuff
    float photon5min = 0.0;
    float photon5max = 14.0;
    

    // fit the sum with Modified Hagedorn Function
    float xlower = 10;
    float xupper = 50;
    //TF1 *f1 = new TF1("f1", "[0]*(1+x/[1])^(-[2])", xlower, xupper);
    //f1->SetParameters(2.09375e+11, 1.54186e+01, 1.32747e+01);

    TF1 *f1 = new TF1("f1", "[0]*pow([1]/x,[2]+[3]*log(x/[1])+ [4]*x)", xlower, xupper);
    f1->SetParameters(2.09375e+9, 1.0, 1.0, 2.0, 0.01);
    //f1->SetParameter(0, 2.09375e+11);  
    //[0]pow([1]/x,[2]+[3]log(x/[1]))

    h_max_photon_pT_sum->Fit("f1", "REMN", "",xlower, xupper);
    h_max_photon_pT_sum->Fit("f1", "REMN", "",xlower, xupper);

    TH1F *h_histo_over_fit = (TH1F *)h_max_photon_pT_sum->Clone("h_histo_over_fit");
    TH1F *h_histo_over_fit_05 = (TH1F *)h_max_photon_pT_sum->Clone("h_histo_over_fit");
    TH1F *h_histo_over_fit_10 = (TH1F *)h_max_photon_pT_sum->Clone("h_histo_over_fit");
    TH1F *h_histo_over_fit_20 = (TH1F *)h_max_photon_pT_sum->Clone("h_histo_over_fit");

    h_histo_over_fit->Divide(f1);


    h_histo_over_fit->GetXaxis()->SetRangeUser(10, 40);
    h_histo_over_fit->SetYTitle("Sim/Fit");
    h_histo_over_fit->SetXTitle("Leading #it{E}_{T}^{#gamma} [GeV]");

    TCanvas *c3 = new TCanvas("c3", "c3", 600, 600);
    h_histo_over_fit->GetYaxis()->SetRangeUser(0.9, 1.1);
    h_histo_over_fit->Draw();
    lineone->Draw("L");
    //h_max_photon_pT_sum->SetXTitle("E_{T}^{#gamma}_{max} [GeV]");
    //h_max_photon_pT_sum->SetYTitle("Counts");
    //h_max_photon_pT_sum->Draw();
    f1->SetLineColor(kRed);
    f1->Draw("same");

    TCanvas *c4 = new TCanvas("can", "", 800, 889);
    c4->Divide(1, 2);

    TPad *pad_1 = (TPad *)c4->cd(1);
    pad_1->SetPad(0, 0.4, 1, 1);
    pad_1->SetTopMargin(0.05);
    pad_1->SetLeftMargin(0.13);
    pad_1->SetBottomMargin(0.002);
    pad_1->SetRightMargin(0.08);
    pad_1->SetLogy();


    frame_et_rec->SetYTitle("counts");
    // frame_et_rec->SetYTitle("d#sigma/d#eta/dE_{T} [pb/GeV]");
    frame_et_rec->GetYaxis()->SetRangeUser(5e3, 1e9);
    // frame_et_rec->GetYaxis()->SetRangeUser(0.2, 4e3);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 40);

    frame_et_rec->GetXaxis()->SetTitleOffset(1.05);
    frame_et_rec->GetYaxis()->SetTitleOffset(1.05);
    frame_et_rec->GetYaxis()->SetTitleSize(0.053);
    frame_et_rec->GetXaxis()->SetLabelSize(0.050);
    frame_et_rec->GetYaxis()->SetLabelSize(0.050);
    frame_et_rec->GetXaxis()->SetLabelOffset(2);
    // frame_et_rec->GetXaxis()->CenterTitle();
    // frame_et_rec->GetYaxis()->CenterTitle();
    frame_et_rec->GetXaxis()->SetNdivisions(505);

    frame_et_rec->Draw("axis");

    h_max_photon_pT_5->SetMarkerStyle(20);
    h_max_photon_pT_5->SetMarkerColor(colors[0]);
    h_max_photon_pT_5->SetLineColor(colors[0]);
    h_max_photon_pT_5->Draw("same");

    h_max_photon_pT_10->SetMarkerStyle(20);
    h_max_photon_pT_10->SetMarkerColor(colors[1]);
    h_max_photon_pT_10->SetLineColor(colors[1]);
    h_max_photon_pT_10->Draw("same");

    h_max_photon_pT_20->SetMarkerStyle(20);
    h_max_photon_pT_20->SetMarkerColor(colors[2]);
    h_max_photon_pT_20->SetLineColor(colors[2]);
    h_max_photon_pT_20->Draw("same p");

    f1->SetLineColor(colors[3]);
    f1->Draw("same");


    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, strMC.c_str(), 0.05);

    myMarkerLineText(0.6, 0.25 + 0.5, 1, colors[0], 20, colors[0], 1, "photon5", 0.05, true);
    myMarkerLineText(0.6, 0.20 + 0.5, 1, colors[1], 20, colors[1], 1, "photon10", 0.05, true);
    myMarkerLineText(0.6, 0.15 + 0.5, 1, colors[2], 20, colors[2], 1, "photon20", 0.05, true);

    TPad *pad_2 = (TPad *)c4->cd(2);
    pad_2->SetPad(0, 0, 1, 0.4);
    pad_2->SetTopMargin(0.023);
    pad_2->SetLeftMargin(0.13);
    pad_2->SetBottomMargin(0.25);
    pad_2->SetRightMargin(0.08);

    frame_et_truth->SetYTitle("Data / Fit");
    frame_et_truth->SetXTitle("Leading #it{E}_{T}^{#gamma} [GeV]");
    frame_et_truth->GetYaxis()->SetNdivisions(506);
    frame_et_truth->GetYaxis()->SetRangeUser(0.9, 1.1);
    frame_et_truth->GetXaxis()->SetRangeUser(10, 40);
    frame_et_truth->GetXaxis()->SetTitleOffset(frame_et_rec->GetXaxis()->GetTitleOffset() * 4 / 6.* 1.4);
    frame_et_truth->GetYaxis()->SetTitleOffset(frame_et_rec->GetYaxis()->GetTitleOffset() * 4 / 6.);
    frame_et_truth->GetYaxis()->SetLabelOffset(frame_et_rec->GetYaxis()->GetLabelOffset() * 4 / 6.);
    frame_et_truth->GetXaxis()->SetLabelSize(frame_et_rec->GetXaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetYaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetTitleSize(frame_et_rec->GetXaxis()->GetTitleSize() * 6 / 4. * 1.2);
    frame_et_truth->GetYaxis()->SetTitleSize(frame_et_rec->GetYaxis()->GetTitleSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetNdivisions(505);
    frame_et_truth->Draw("axis");

    h_histo_over_fit->Draw("same");
    lineone->Draw("L");

    c4->SaveAs("figures/combine.pdf");


}