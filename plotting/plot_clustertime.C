#include "plotcommon.h"

void plot_clustertime()
{
    init_plot();

    TFile *f_in_time = new TFile("/sphenix/user/shuhangli/ppg12/showershapecheck/clustertime.root");

    float clusterminET = 10;
    float clustermaxET = 50;

    float rebinx = 2;

    std::vector<float> eta_bins = {-0.7, -0.2, 0.2, 0.7};

    TH2F *h_cluster_time = (TH2F *)f_in_time->Get("h_cluster_time");
    TH2F *h_cluster_time_withjet = (TH2F *)f_in_time->Get("h_cluster_time_withjet");
    TH2F *h_cluster_time_large_wetacogx = (TH2F *)f_in_time->Get("h_cluster_time_large_wetacogx");

    TH3F *h_cluster_time_eta = (TH3F *)f_in_time->Get("h_cluster_time_eta");
    TH3F *h_cluster_time_withjet_eta = (TH3F *)f_in_time->Get("h_cluster_time_withjet_eta");
    TH3F *h_cluster_time_large_wetacogx_eta = (TH3F *)f_in_time->Get("h_cluster_time_large_wetacogx_eta");

    h_cluster_time_eta->GetYaxis()->SetRangeUser(clusterminET, clustermaxET);
    h_cluster_time_withjet_eta->GetYaxis()->SetRangeUser(clusterminET, clustermaxET);
    h_cluster_time_large_wetacogx_eta->GetYaxis()->SetRangeUser(clusterminET, clustermaxET);

    TH2F *h_cluster_eta_time = (TH2F *)h_cluster_time_eta->Project3D("xz");
    TH2F *h_cluster_eta_time_withjet = (TH2F *)h_cluster_time_withjet_eta->Project3D("xz");
    TH2F *h_cluster_eta_time_large_wetacogx = (TH2F *)h_cluster_time_large_wetacogx_eta->Project3D("xz");

    h_cluster_eta_time->GetYaxis()->SetRangeUser(-1, 1);
    h_cluster_eta_time->SetYTitle("Cluster Time [samples]");
    h_cluster_eta_time->GetXaxis()->SetRangeUser(-0.7, 0.7);
    h_cluster_eta_time->SetXTitle("Cluster #eta");
    h_cluster_eta_time->GetZaxis()->SetRangeUser(10, 300);

    std::vector<TH1F *> h_cluster_eta_time_withjet_list;
    std::vector<TH1F *> h_cluster_eta_time_large_wetacogx_list;

    std::vector<int> colors = {kPink + 8, kSpring - 7, kAzure - 3, kViolet + 3, kOrange + 10, kYellow + 4};

    for (int i = 0; i < (int)eta_bins.size() - 1; i++)
    {
        h_cluster_eta_time_withjet->GetXaxis()->SetRangeUser(eta_bins[i], eta_bins[i + 1]);
        h_cluster_eta_time_withjet_list.push_back((TH1F *)h_cluster_eta_time_withjet->ProjectionY(Form("h_cluster_eta_time_withjet_%d", i)));

        h_cluster_eta_time_large_wetacogx->GetXaxis()->SetRangeUser(eta_bins[i], eta_bins[i + 1]);
        h_cluster_eta_time_large_wetacogx_list.push_back((TH1F *)h_cluster_eta_time_large_wetacogx->ProjectionY(Form("h_cluster_eta_time_large_wetacogx_%d", i)));
    }
    h_cluster_eta_time_withjet->GetXaxis()->SetRangeUser(-0.7, 0.7);
    h_cluster_eta_time_large_wetacogx->GetXaxis()->SetRangeUser(-0.7, 0.7);

    TProfile *h_cluster_eta_time_pfx = h_cluster_eta_time->ProfileX("h_cluster_eta_time_pfx");
    TProfile *h_cluster_eta_time_withjet_pfx = h_cluster_eta_time_withjet->ProfileX("h_cluster_eta_time_withjet_pfx");
    TProfile *h_cluster_eta_time_large_wetacogx_pfx = h_cluster_eta_time_large_wetacogx->ProfileX("h_cluster_eta_time_large_wetacogx_pfx");

    h_cluster_time->RebinX(rebinx);
    h_cluster_time_withjet->RebinX(rebinx);
    h_cluster_time_large_wetacogx->RebinX(rebinx);

    h_cluster_time->GetYaxis()->SetRangeUser(clusterminET, clustermaxET);
    h_cluster_time_withjet->GetYaxis()->SetRangeUser(clusterminET, clustermaxET);
    h_cluster_time_large_wetacogx->GetYaxis()->SetRangeUser(clusterminET, clustermaxET);

    TH1D *h_cluster_time_proj = h_cluster_time->ProjectionX("h_cluster_time_proj");
    TH1D *h_cluster_time_withjet_proj = h_cluster_time_withjet->ProjectionX("h_cluster_time_withjet_proj");
    TH1D *h_cluster_time_large_wetacogx_proj = h_cluster_time_large_wetacogx->ProjectionX("h_cluster_time_large_wetacogx_proj");
    //find fraction of h_cluster_time_proj with time < -0.5
    float fraction = h_cluster_time_proj->Integral(1, h_cluster_time_proj->FindBin(-0.5)) / h_cluster_time_proj->Integral();
    std::cout << "fraction of h_cluster_time_proj with time < -0.5: " << fraction << std::endl;

    // normalization
    h_cluster_time_proj->Scale(1.0 / h_cluster_time_proj->Integral());
    h_cluster_time_withjet_proj->Scale(1.0 / h_cluster_time_withjet_proj->Integral());
    h_cluster_time_large_wetacogx_proj->Scale(1.0 / h_cluster_time_large_wetacogx_proj->Integral());

    h_cluster_time_proj->GetXaxis()->SetRangeUser(-2, 1);
    h_cluster_time_proj->GetYaxis()->SetRangeUser(0, 0.06 * rebinx);
    h_cluster_time_proj->SetXTitle("Cluster Time [samples]");

    TCanvas *c1 = new TCanvas("c1", "", 600, 600);

    h_cluster_time_proj->SetLineColor(kBlack);
    h_cluster_time_proj->SetMarkerColor(kBlack);

    h_cluster_time_withjet_proj->SetLineColor(kRed);
    h_cluster_time_withjet_proj->SetMarkerColor(kRed);

    h_cluster_time_large_wetacogx_proj->SetLineColor(kBlue);
    h_cluster_time_large_wetacogx_proj->SetMarkerColor(kBlue);

    h_cluster_time_proj->Draw("HIST");
    h_cluster_time_withjet_proj->Draw("same HIST");
    h_cluster_time_large_wetacogx_proj->Draw("same HIST");

    myText(0.5 - 0.3, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5 - 0.3, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5 - 0.3, 0.80, 1, Form("%2.0f<E_{T}^{#gamma}<%2.0fGeV, |#eta^{#gamma}|<0.7", clusterminET, clustermaxET), 0.05);

    myMarkerLineText(0.55 - 0.3, 0.7, 1, kBlack, 20, kBlack, 1, "All clusters", 0.05, true);
    myMarkerLineText(0.55 - 0.3, 0.65, 1, kRed, 20, kRed, 1, "With jet#Delta#phi<7#pi/8", 0.05, true);
    myMarkerLineText(0.55 - 0.3, 0.6, 1, kBlue, 20, kBlue, 1, "w_{#eta}^{COGX}>1", 0.05, true);

    TCanvas *c2 = new TCanvas("c2", "", 600, 600);
    // right margin
    c2->SetRightMargin(0.15);

    h_cluster_eta_time->Draw("colz");

    h_cluster_eta_time_withjet_pfx->SetLineColor(kRed);
    h_cluster_eta_time_large_wetacogx_pfx->SetLineColor(kBlue);

    h_cluster_eta_time_withjet_pfx->Draw("hist same");
    h_cluster_eta_time_large_wetacogx_pfx->Draw("hist same");

    myText(0.5 - 0.3, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5 - 0.3, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5 - 0.3, 0.80, 1, Form("%2.0f<E_{T}^{#gamma}<%2.0fGeV", clusterminET, clustermaxET), 0.05);

    myMarkerLineText(0.55 - 0.3, 0.65 + 0.1, 1, kRed, 20, kRed, 1, "With jet#Delta#phi<7#pi/8", 0.05, true);
    myMarkerLineText(0.55 - 0.3, 0.6 + 0.1, 1, kBlue, 20, kBlue, 1, "w_{#eta}^{COGX}>1", 0.05, true);

    TCanvas *c3 = new TCanvas("c3", "", 600, 600);

    for (int i = 0; i < (int)eta_bins.size() - 1; i++)
    {
        h_cluster_eta_time_withjet_list[i]->Rebin(rebinx);
        h_cluster_eta_time_withjet_list[i]->SetXTitle("Cluster Time [samples]");
        h_cluster_eta_time_withjet_list[i]->Scale(1.0 / h_cluster_eta_time_withjet_list[i]->Integral());
        h_cluster_eta_time_withjet_list[i]->GetYaxis()->SetRangeUser(0, 0.07 * rebinx);
        h_cluster_eta_time_withjet_list[i]->GetXaxis()->SetRangeUser(-1.0, 0.5);
        h_cluster_eta_time_withjet_list[i]->SetLineColor(colors[i]);
        h_cluster_eta_time_withjet_list[i]->SetMarkerColor(colors[i]);
        h_cluster_eta_time_withjet_list[i]->SetLineStyle(5);
        h_cluster_eta_time_withjet_list[i]->Draw(i == 0 ? "HIST" : "same HIST");
    }
    for (int i = 0; i < (int)eta_bins.size() - 1; i++)
    {
        h_cluster_eta_time_large_wetacogx_list[i]->Rebin(rebinx);
        h_cluster_eta_time_large_wetacogx_list[i]->Scale(1.0 / h_cluster_eta_time_large_wetacogx_list[i]->Integral());
        h_cluster_eta_time_large_wetacogx_list[i]->GetYaxis()->SetRangeUser(0, 0.06 * rebinx);
        h_cluster_eta_time_large_wetacogx_list[i]->GetXaxis()->SetRangeUser(-2, 2);
        h_cluster_eta_time_large_wetacogx_list[i]->SetLineColor(colors[i]);
        h_cluster_eta_time_large_wetacogx_list[i]->SetMarkerColor(colors[i]);
        h_cluster_eta_time_large_wetacogx_list[i]->SetLineStyle(1);
        h_cluster_eta_time_large_wetacogx_list[i]->Draw("same HIST");

        myText(0.5 - 0.3, 0.9, 1, strleg1.c_str(), 0.04);
        myText(0.5 - 0.3, 0.85, 1, strleg2.c_str(), 0.04);
        myText(0.5 - 0.3, 0.80, 1, Form("%2.0f<E_{T}^{#gamma}<%2.0fGeV, |#eta^{#gamma}|<0.7", clusterminET, clustermaxET), 0.04);

        myMarkerLineText(0.55 - 0.3, 0.7 + 0.05, 0, colors[0], 20, colors[0], 1, Form("%.1f<|#eta^{#gamma}|<%.1f", eta_bins[0], eta_bins[1]), 0.04, true);
        myMarkerLineText(0.55 - 0.3, 0.65 + 0.05, 0, colors[1], 20, colors[1], 1, Form("%.1f<|#eta^{#gamma}|<%.1f", eta_bins[1], eta_bins[2]), 0.04, true);
        myMarkerLineText(0.55 - 0.3, 0.6 + 0.05, 0, colors[2], 20, colors[2], 1, Form("%.1f<|#eta^{#gamma}|<%.1f", eta_bins[2], eta_bins[3]), 0.04, true);

        myMarkerLineText(0.55 + 0.1, 0.75 + 0.1, 0, kBlack, 20, kBlack, 5, Form("With jet#Delta#phi<7#pi/8"), 0.04, true);
        myMarkerLineText(0.55 + 0.1, 0.7 + 0.1, 0, kBlack, 20, kBlack, 1, Form("w_{#eta}^{COGX}>1"), 0.04, true);
    }
}
