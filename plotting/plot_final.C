#include "plotcommon.h"

void plot_final()
{
    init_plot();

    std::string MCstring = "JETPHOX";

    std::string reweightedstring = "";

    std::string datastring = "data";
    std::string bg_MCstring = "Inclusive Sim";

    float datalumi = 16.8468;

    float deta = 1.4;

    TFile *fin_data = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nom.root");

    TFile *fin_syst = new TFile("/sphenix/user/shuhangli/ppg12/plotting/rootFiles/syst_sum.root");
    TFile *fin_NLO = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX.root");
    TFile *fin_mc = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nom_mc.root");

    TH1F *h_data = (TH1F *)fin_data->Get("h_unfold_sub_result");
    h_data->Scale(1.0 / deta);
    TH1F *h_data_cp = (TH1F* )h_data->Clone("h_data_cp");
    TH1F *h_pythia = (TH1F *)fin_data->Get("h_truth_pT_0");
    h_pythia->Scale(1.0 / deta);
    TH1F *h_common_cluster_data = (TH1F *)fin_data->Get("h_common_cluster_0");
    h_common_cluster_data->Scale(1.0 / deta);
    TH1F *h_common_cluster_mc = (TH1F *)fin_mc->Get("h_common_cluster_0");
    h_common_cluster_mc->Scale(1.0 / deta);
    TH1F *h_tight_iso_cluster_data = (TH1F *)fin_data->Get("h_tight_iso_cluster_0");
    h_tight_iso_cluster_data->Scale(1.0 / deta);
    TH1F *h_sub_data = (TH1F *)fin_data->Get("h_data_sub_leak");
    TH1F *h_sub_data_unfold = (TH1F *)fin_data->Get("h_unfold_sub_result_woeff");
    TH1F *h_tight_iso_cluster_mc = (TH1F *)fin_mc->Get("h_tight_iso_cluster_0");

    TH1F *h_NLO = (TH1F *)fin_NLO->Get("h_truth_pT");
    h_NLO->Scale(1.0 / deta);

    TH1F *h_data_NLO = (TH1F *)h_data->Clone("h_data_NLO");
    h_data_NLO->Divide(h_NLO);

    TH1F *h_NLO_data = (TH1F *)h_NLO->Clone("h_NLO_data");
    h_NLO_data->Divide(h_data);

    TH1F *h_pythia_data = (TH1F *)h_pythia->Clone("h_pythia_data");
    h_pythia_data->Divide(h_data);

    TH1F *h_syst_low = (TH1F *)fin_syst->Get("h_sum_low");
    TH1F *h_syst_high = (TH1F *)fin_syst->Get("h_sum_high");

    TH1F *h_syst_rel_low = (TH1F *)fin_syst->Get("h_sum_rel_low");
    TH1F *h_syst_rel_high = (TH1F *)fin_syst->Get("h_sum_rel_high");

    TGraphAsymmErrors *g_syst = new TGraphAsymmErrors(h_data);
    TGraphAsymmErrors *g_syst_rel_NLO = new TGraphAsymmErrors(h_NLO_data);
    TGraphAsymmErrors *g_syst_rel_pythia = new TGraphAsymmErrors(h_pythia_data);
    TGraphAsymmErrors *g_syst_rel = new TGraphAsymmErrors(h_data_NLO);

    for (int i = 0; i < h_data->GetNbinsX(); i++)
    {
        float xlowerror = g_syst->GetErrorXlow(i);
        float xuperror = g_syst->GetErrorXhigh(i);

        g_syst->SetPointError(i, xlowerror, xuperror, h_syst_low->GetBinContent(i + 1), h_syst_high->GetBinContent(i + 1));

        float y = g_syst_rel_NLO->GetY()[i];
        std::cout << "y: " << y << std::endl;
        g_syst_rel_NLO->SetPointError(i, xlowerror, xuperror, y * h_syst_rel_high->GetBinContent(i + 1), y * h_syst_rel_low->GetBinContent(i + 1));

        std::cout << "h_syst_rel_low->GetBinContent(i + 1): " << h_syst_rel_low->GetBinContent(i + 1) << std::endl;

        float y_pythia = g_syst_rel_pythia->GetY()[i];

        std::cout << "y_pythia: " << y_pythia << std::endl;
        std::cout << "xlowerror" << xlowerror << std::endl;
        g_syst_rel_pythia->SetPointError(i, xlowerror, xuperror, y_pythia * h_syst_rel_high->GetBinContent(i + 1), y_pythia * h_syst_rel_low->GetBinContent(i + 1));

        g_syst_rel->SetPoint(i, g_syst_rel_NLO->GetX()[i], 1);
        g_syst_rel->SetPointError(i, xlowerror, xuperror, h_syst_rel_low->GetBinContent(i + 1), h_syst_rel_high->GetBinContent(i + 1));
    }

    TCanvas *c1 = new TCanvas("can", "", 800, 889);
    c1->Divide(1, 2);

    TPad *pad_1 = (TPad *)c1->cd(1);
    pad_1->SetPad(0, 0.4, 1, 1);
    pad_1->SetTopMargin(0.05);
    pad_1->SetLeftMargin(0.13);
    pad_1->SetBottomMargin(0.03);
    pad_1->SetRightMargin(0.08);
    pad_1->SetLogy();

    frame_et_rec->SetYTitle("d#sigma/d#eta/dE_{T} [pb/GeV]");
    frame_et_rec->GetYaxis()->SetRangeUser(2, 4e3);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 30);

    frame_et_rec->GetXaxis()->SetTitleOffset(0.98);
    frame_et_rec->GetYaxis()->SetTitleOffset(1.15);
    frame_et_rec->GetXaxis()->SetLabelSize(0.045);
    frame_et_rec->GetYaxis()->SetLabelSize(0.045);
    frame_et_rec->GetXaxis()->SetLabelOffset(2);
    // frame_et_rec->GetXaxis()->CenterTitle();
    // frame_et_rec->GetYaxis()->CenterTitle();
    frame_et_rec->GetXaxis()->SetNdivisions(505);

    frame_et_rec->Draw("axis");

    // PHENIX data
    //  Number of data points

    const Int_t n = 18;

    // Central pT
    Double_t x[n] = {
        5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75,
        11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0};

    // pT low edges
    Double_t xLow[n] = {
        5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
        10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0};

    // pT high edges
    Double_t xHigh[n] = {
        5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
        12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0};

    // Cross section values, E d^3sigma/dp^3 [pb * GeV^-2 * c^3]
    Double_t y[n] = {
        1140, 613, 348, 231, 136, 93, 67, 48.3, 32.1, 20.4,
        9.8, 2.97, 1.06, 0.34, 0.173, 0.088, 0.042, 0.029};

    // Stat errors (asymmetric: + and -)
    Double_t statUp[n] = {
        30, 19, 13, 9, 6, 4, 3.2, 2.5, 1.9, 1.5,
        0.4, 0.19, 0.1, 0.06, 0.034, 0.02, 0.014, 0.009};

    Double_t statDown[n] = {
        30, 19, 13, 9, 6, 4, 3.2, 2.5, 1.9, 1.5,
        0.4, 0.19, 0.1, 0.06, 0.034, 0.021, 0.015, 0.014};

    // Sys errors (asymmetric: + and -)
    Double_t sysUp[n] = {
        478, 221, 101, 62, 31, 20, 13.4, 9.2, 6.1, 3.7,
        1.7, 0.48, 0.17, 0.05, 0.028, 0.015, 0.007, 0.004};

    Double_t sysDown[n] = {
        478, 221, 101, 62, 31, 20, 13.4, 9.2, 6.1, 3.7,
        1.7, 0.48, 0.17, 0.05, 0.028, 0.015, 0.007, 0.004};

    // Create the TGraphAsymmErrors for STAT
    TGraphAsymmErrors *gStat_PHENIX = new TGraphAsymmErrors(n);
    gStat_PHENIX->SetName("gStat");
    gStat_PHENIX->SetTitle("Cross Section (scaled) with Stat. Errors");

    // Create the TGraphAsymmErrors for SYS
    TGraphAsymmErrors *gSys_PHENIX = new TGraphAsymmErrors(n);
    gSys_PHENIX->SetName("gSys");
    gSys_PHENIX->SetTitle("Cross Section (scaled) with Syst. Errors");

    // Common factor that does NOT include pT, since pT changes per bin
    // Overall factor: 2*pi * 2 * deta, with deta=0.7
    //Double_t factorCommon = 2.0 * TMath::Pi() * 2.0 * 0.7; // ~ 8.7964
    Double_t factorCommon = 2.0 * TMath::Pi(); 

    for (Int_t i = 0; i < n; i++)
    {
        // x errors = half bin widths
        Double_t exl = x[i] - xLow[i];
        Double_t exh = xHigh[i] - x[i];

        // Scale factor for the i-th point = pT * factorCommon
        Double_t scaleFactor = x[i] * factorCommon;

        // Scale the central value and uncertainties
        Double_t yScaled = y[i] * scaleFactor;
        Double_t statUpScaled = statUp[i] * scaleFactor;
        Double_t statDownScaled = statDown[i] * scaleFactor;
        Double_t sysUpScaled = sysUp[i] * scaleFactor;
        Double_t sysDownScaled = sysDown[i] * scaleFactor;

        // Fill the TGraphAsymmErrors (Stat)
        gStat_PHENIX->SetPoint(i, x[i], yScaled);
        gStat_PHENIX->SetPointError(i, exl, exh, statDownScaled, statUpScaled);

        // Fill the TGraphAsymmErrors (Sys)
        gSys_PHENIX->SetPoint(i, x[i], yScaled);
        gSys_PHENIX->SetPointError(i, exl, exh, sysDownScaled, sysUpScaled);
    }

    // getting a different NLO
    ifstream myfile;
    TGraph *tphoton = new TGraph();
    TGraph *tphoton05 = new TGraph();
    TGraph *tphoton02 = new TGraph();
    myfile.open("sphenix_nlo/photons_newphenix_sc1.dat");
    if (!myfile)
    {
        cout << "No file found!" << endl;
        return;
    }

    for (int index = 0; index < 38; index++)
    {
        double pt;
        double yield;
        double fragyield;
        double foo;
        myfile >> pt >> yield >> fragyield >> foo;
        // rescale points
        yield = (yield + fragyield) * factorCommon * pt;
        tphoton->SetPoint(index, pt, yield);
    }
    myfile.close();

    myfile.open("sphenix_nlo/photons_newphenix_sc05.dat");
    if (!myfile)
    {
        cout << "No file found!" << endl;
        return;
    }

    for (int index = 0; index < 38; index++)
    {
        double pt;
        double yield;
        double fragyield;
        double foo;
        myfile >> pt >> yield >> fragyield >> foo;
        // rescale points
        yield = (yield + fragyield) * factorCommon * pt;
        tphoton05->SetPoint(index, pt, yield);
    }
    myfile.close();

    myfile.open("sphenix_nlo/photons_newphenix_sc2.dat");
    if (!myfile)
    {
        cout << "No file found!" << endl;
        return;
    }

    for (int index = 0; index < 38; index++)
    {
        double pt;
        double yield;
        double fragyield;
        double foo;
        myfile >> pt >> yield >> fragyield >> foo;
        // rescale points
        yield = (yield + fragyield) * factorCommon * pt;
        tphoton02->SetPoint(index, pt, yield);
    }
    myfile.close();

    tphoton->SetMarkerStyle(21);
    tphoton->SetMarkerColor(4);
    tphoton->SetLineColor(4);
    tphoton->Draw(" l,same");

    tphoton05->SetMarkerStyle(21);
    tphoton05->SetMarkerColor(4);
    tphoton05->SetLineColor(4);
    tphoton05->SetLineStyle(7);
    tphoton05->Draw("l,same");

    tphoton02->SetMarkerStyle(21);
    tphoton02->SetMarkerColor(4);
    tphoton02->SetLineColor(4);
    tphoton02->SetLineStyle(7);
    tphoton02->Draw("l,same");

    gStat_PHENIX->SetMarkerStyle(28);
    gStat_PHENIX->SetMarkerSize(1.5);
    gStat_PHENIX->SetMarkerColor(kRed + 1);
    gStat_PHENIX->SetLineColor(kRed + 1);
    //gStat_PHENIX->Draw("P");

    gSys_PHENIX->SetMarkerStyle(28);
    gSys_PHENIX->SetMarkerSize(1.5);
    gSys_PHENIX->SetMarkerColor(kRed);
    gSys_PHENIX->SetLineColor(kRed);
    gSys_PHENIX->SetFillColorAlpha(kRed, 0.25);
    //gSys_PHENIX->Draw("2 same");

    g_syst->SetMarkerStyle(20);
    g_syst->SetMarkerColor(kAzure + 2);
    g_syst->SetLineColor(kAzure + 2);
    g_syst->SetFillColorAlpha(kAzure + 2, 0.5);

    g_syst->Draw("2 same");

    h_NLO->SetMarkerStyle(25);
    h_NLO->SetMarkerColor(kPink + 8);
    h_NLO->SetLineColor(kPink + 8);

    h_NLO->Draw("same");

    h_pythia->SetMarkerStyle(27);
    h_pythia->SetMarkerColor(kSpring - 7);
    h_pythia->SetLineColor(kSpring - 7);
    h_pythia->SetMarkerSize(2);

    h_pythia->Draw("same");
    
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerColor(kBlack);
    h_data->SetLineColor(kBlack);

    h_data->Draw("same");

    //--------------------------------------------------lower panel

    // myBoxTextAlpha(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor, Double_t falpha,Int_t lcolor,Int_t lstyle, const char *text)

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, "|#eta^{#gamma}|<0.7, E_{T}^{iso, R=0.3} < 4 GeV", 0.05);

    myBoxTextAlpha(0.55, 0.25 + 0.5, 0.025, kAzure + 2, 0.5, kAzure + 2, 1, "");
    myMarkerLineText(0.55, 0.25 + 0.5, 1, kBlack, 20, kBlack, 1, Form("Data %s", reweightedstring.c_str()), 0.05, true);
    myMarkerLineText(0.55, 0.20 + 0.5, 1, kPink + 8, 25, kPink + 8, 1, MCstring.c_str(), 0.05, true);
    myMarkerLineText(0.55, 0.15 + 0.5, 2, kSpring - 7, 27, kSpring - 7, 1, "Pythia", 0.05, true);
    //myBoxTextAlpha(0.55, 0.10 + 0.5, 0.025, kRed, 0.25, kRed, 1, "");
    //myMarkerLineText(0.55, 0.10 + 0.5, 2, kRed + 1, 28, kRed + 1, 1, "PHENIX", 0.05, true);
    //isolation creteria
    //myText(0.55, 0.10 + 0.5, 1, "Isolation: E_{T}^{iso}<4 GeV ", 0.05);
    //myText(0.55, 0.55, 1, "Phys. Rev. D 86 072008", 0.05);
    myMarkerLineText(0.55, 0.60, 0, kBlue, 21, kBlue, 1, "NLO pQCD by Werner V.", 0.05, true);
    myText(0.50, 0.55, 1, "(no E_{T}^{iso} requirement)", 0.05);

    TPad *pad_2 = (TPad *)c1->cd(2);
    pad_2->SetPad(0, 0, 1, 0.4);
    pad_2->SetTopMargin(0.02);
    pad_2->SetLeftMargin(0.13);
    pad_2->SetBottomMargin(0.25);
    pad_2->SetRightMargin(0.08);

    frame_et_truth->SetYTitle("theory/data");
    frame_et_truth->GetYaxis()->SetNdivisions(506);
    frame_et_truth->GetYaxis()->SetRangeUser(0.2, 2);
    frame_et_truth->GetXaxis()->SetRangeUser(10, 30);
    frame_et_truth->GetYaxis()->SetTitleOffset(frame_et_rec->GetYaxis()->GetTitleOffset() * 4 / 6.);
    frame_et_truth->GetYaxis()->SetLabelOffset(frame_et_rec->GetYaxis()->GetLabelOffset() * 4 / 6.);
    frame_et_truth->GetXaxis()->SetLabelSize(frame_et_rec->GetXaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetYaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetTitleSize(frame_et_rec->GetXaxis()->GetTitleSize() * 6 / 4.);
    frame_et_truth->GetYaxis()->SetTitleSize(frame_et_rec->GetYaxis()->GetTitleSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetNdivisions(505);
    frame_et_truth->Draw("axis");

    g_syst_rel->SetMarkerStyle(20);
    g_syst_rel->SetMarkerColor(kAzure + 2);
    g_syst_rel->SetLineColor(kAzure + 2);
    g_syst_rel->SetFillColorAlpha(kAzure + 2, 0.5);

    lineone->Draw("L");

    g_syst_rel->Draw("2 same");

    h_NLO_data->SetMarkerStyle(25);
    h_NLO_data->SetMarkerColor(kPink + 8);
    h_NLO_data->SetLineColor(kPink + 8);

    h_NLO_data->Draw("same");

    h_pythia_data->SetMarkerStyle(27);
    h_pythia_data->SetMarkerColor(kSpring - 7);
    h_pythia_data->SetLineColor(kSpring - 7);
    h_pythia_data->SetMarkerSize(2);

    h_pythia_data->Draw("same");

    std::string outputname = "figures/final.pdf";

    c1->SaveAs(outputname.c_str());



    pad_1->cd();
    frame_et_rec->GetYaxis()->SetRangeUser(1, 2e4);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 30);
    frame_et_rec->Draw("axis");

    h_common_cluster_data->SetMarkerStyle(20);
    h_common_cluster_data->SetMarkerColor(kBlack);
    h_common_cluster_data->SetLineColor(kBlack);

    h_common_cluster_data->Draw("same");

    h_common_cluster_mc->SetMarkerStyle(24);
    h_common_cluster_mc->SetMarkerColor(kPink + 8);
    h_common_cluster_mc->SetLineColor(kPink + 8);

    h_common_cluster_mc->Draw("same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, "Clusters pass preliminary cuts", 0.05);

    myMarkerLineText(0.25, 0.25, 1, kBlack, 20, kBlack, 1, Form("%s", datastring.c_str()), 0.05, true);
    myMarkerLineText(0.25, 0.20, 1, kPink + 8, 24, kPink + 8, 1, bg_MCstring.c_str(), 0.05, true);

    pad_2->cd();
    frame_et_truth->GetYaxis()->SetRangeUser(0.5, 1.5);
    frame_et_truth->GetXaxis()->SetRangeUser(10, 30);
    frame_et_truth->SetYTitle("data/MC");
    frame_et_truth->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
    frame_et_truth->Draw("axis");

    lineone->Draw("L");

    TH1F *h_common_cluster_data_ratio = (TH1F *)h_common_cluster_data->Clone("h_common_cluster_data_ratio");
    h_common_cluster_data_ratio->Divide(h_common_cluster_mc);

    h_common_cluster_data_ratio->SetMarkerStyle(20);
    h_common_cluster_data_ratio->SetMarkerColor(kBlack);
    h_common_cluster_data_ratio->SetLineColor(kBlack);

    h_common_cluster_data_ratio->Draw("same");

    c1->SaveAs("figures/final_common_cluster.pdf");

    pad_1->cd();
    frame_et_rec->GetYaxis()->SetRangeUser(0.1, 1e3);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 30);
    frame_et_rec->Draw("axis");

    h_tight_iso_cluster_data->SetMarkerStyle(20);
    h_tight_iso_cluster_data->SetMarkerColor(kBlack);
    h_tight_iso_cluster_data->SetLineColor(kBlack);

    h_tight_iso_cluster_data->Draw("same");

    h_tight_iso_cluster_mc->SetMarkerStyle(24);
    h_tight_iso_cluster_mc->SetMarkerColor(kPink + 8);
    h_tight_iso_cluster_mc->SetLineColor(kPink + 8);

    h_tight_iso_cluster_mc->Draw("same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, "Tight iso clusters", 0.05);

    myMarkerLineText(0.25, 0.25, 1, kBlack, 20, kBlack, 1, Form("%s", datastring.c_str()), 0.05, true);
    myMarkerLineText(0.25, 0.20, 1, kPink + 8, 24, kPink + 8, 1, bg_MCstring.c_str(), 0.05, true);

    pad_2->cd();
    frame_et_truth->GetYaxis()->SetRangeUser(0.2, 1.5);
    frame_et_truth->GetXaxis()->SetRangeUser(10, 30);
    frame_et_truth->SetYTitle("data/MC");
    frame_et_truth->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
    frame_et_truth->Draw("axis");

    lineone->Draw("L");

    TH1F *h_tight_iso_cluster_data_ratio = (TH1F *)h_tight_iso_cluster_data->Clone("h_tight_iso_cluster_data_ratio");
    h_tight_iso_cluster_data_ratio->Divide(h_tight_iso_cluster_mc);

    h_tight_iso_cluster_data_ratio->SetMarkerStyle(20);
    h_tight_iso_cluster_data_ratio->SetMarkerColor(kBlack);
    h_tight_iso_cluster_data_ratio->SetLineColor(kBlack);

    h_tight_iso_cluster_data_ratio->Draw("same");

    c1->SaveAs("figures/final_tight_iso_cluster.pdf");

    /*
    TCanvas *c2 = new TCanvas("can", "", 800, 900);
    //log y
    c2->SetLogy();
    TH2F *frame = new TH2F("frame", "", 1, 10, 30, 1, 0.1, 1e5);
    frame->SetXTitle("E_{T}^{#gamma} [GeV]");
    frame->SetYTitle("d#sigma/dE_{T} [pb/GeV]");
    frame->GetYaxis()->SetRangeUser(0.5, 5e4);
    frame->GetXaxis()->SetRangeUser(10, 30);

    frame->Draw("axis");
    */
    pad_1->cd();
    // common, tight iso, sub, sub with unfold and final result
    frame_et_rec->GetYaxis()->SetRangeUser(5, 2e5);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 30);
    frame_et_rec->SetYTitle("dN/dE_{T} [GeV^{-1}]");
    frame_et_rec->Draw("axis");

    std::vector<int> marker_styles = {20, 24, 25, 27, 28};
    std::vector<float> marker_sizes = {1, 1, 1.0, 1.5, 1.5};
    std::vector<int> colors = {kPink + 8, kSpring - 7, kAzure - 3, kOrange + 10, kViolet + 3, kYellow + 4};

    // scale things up by lumi
    h_common_cluster_data->Scale(datalumi);
    h_tight_iso_cluster_data->Scale(datalumi);
    h_sub_data->Scale(datalumi);
    h_sub_data_unfold->Scale(datalumi);
    h_data->Scale(datalumi);

    h_common_cluster_data->SetMarkerStyle(marker_styles[0]);
    h_common_cluster_data->SetMarkerColor(colors[0]);
    h_common_cluster_data->SetLineColor(colors[0]);
    h_common_cluster_data->SetMarkerSize(marker_sizes[0]);

    h_common_cluster_data->Draw("same");

    h_tight_iso_cluster_data->SetMarkerStyle(marker_styles[1]);
    h_tight_iso_cluster_data->SetMarkerColor(colors[1]);
    h_tight_iso_cluster_data->SetLineColor(colors[1]);
    h_tight_iso_cluster_data->SetMarkerSize(marker_sizes[1]);
    h_tight_iso_cluster_data->Draw("same");

    h_sub_data->SetMarkerStyle(marker_styles[2]);
    h_sub_data->SetMarkerColor(colors[2]);
    h_sub_data->SetLineColor(colors[2]);
    h_sub_data->SetMarkerSize(marker_sizes[2]);
    h_sub_data->Draw("same");

    h_sub_data_unfold->SetMarkerStyle(marker_styles[3]);
    h_sub_data_unfold->SetMarkerColor(colors[3]);
    h_sub_data_unfold->SetLineColor(colors[3]);
    h_sub_data_unfold->SetMarkerSize(marker_sizes[3]);
    h_sub_data_unfold->Draw("same");

    h_data->SetMarkerStyle(marker_styles[4]);
    h_data->SetMarkerColor(colors[4]);
    h_data->SetLineColor(colors[4]);
    h_data->SetMarkerSize(marker_sizes[4]);
    h_data->Draw("same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, "|#eta^{#gamma}|<0.7", 0.05);

    myMarkerLineText(0.6, 0.25 + 0.5, 1, colors[0], marker_styles[0], colors[0], 1, "w/ ps cut", 0.05, true);
    myMarkerLineText(0.6, 0.20 + 0.5, 1, colors[1], marker_styles[1], colors[1], 1, "tight iso", 0.05, true);
    myMarkerLineText(0.6, 0.15 + 0.5, 1, colors[2], marker_styles[2], colors[2], 1, "purity corrected", 0.05, true);
    myMarkerLineText(0.6, 0.10 + 0.5, 1.5, colors[3], marker_styles[3], colors[3], 1, "unfolded", 0.05, true);
    myMarkerLineText(0.6, 0.05 + 0.5, 1.5, colors[4], marker_styles[4], colors[4], 1, "efficiency corrected", 0.05, true);

    pad_2->cd();
    frame_et_truth->GetYaxis()->SetRangeUser(0.8, 1.1);
    frame_et_truth->GetXaxis()->SetRangeUser(10, 30);
    frame_et_truth->SetYTitle("after/before unfolding");
    frame_et_truth->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");
    frame_et_truth->Draw("axis");

    lineone->Draw("L");

    TH1F *h_unfold_ratio = (TH1F *)h_sub_data_unfold->Clone("h_unfold_ratio");
    for (int i = 2; i <= h_unfold_ratio->GetNbinsX(); i++)
    {
        float centerX = h_unfold_ratio->GetBinCenter(i);
        int bin = h_sub_data->FindBin(centerX);
        float value = h_sub_data->GetBinContent(bin);
        float error = h_sub_data->GetBinError(bin);

        if (value != 0)
        {
            h_unfold_ratio->SetBinContent(i, h_sub_data_unfold->GetBinContent(i) / value);
            std::cout << "bincenter: " << centerX << " " << h_sub_data_unfold->GetBinContent(i) << " " << value << std::endl;
            h_unfold_ratio->SetBinError(i, 0);
        }
    }

    h_unfold_ratio->SetMarkerStyle(marker_styles[3]);
    h_unfold_ratio->SetMarkerColor(colors[3]);
    h_unfold_ratio->SetLineColor(colors[3]);
    h_unfold_ratio->SetMarkerSize(marker_sizes[3]);
    h_unfold_ratio->Draw("same");

    c1->SaveAs("figures/final_all.pdf");


    //-----------------------------------------------------------------
    // only compare to PHENIX
    init_plot();
    h_data = h_data_cp;
    TCanvas *c2 = new TCanvas("can2", "", 800, 800);
    //log y 
    c2->SetLogy();

    frame_et_rec->SetYTitle("d#sigma/d#eta/dE_{T} [pb/GeV]");
    frame_et_rec->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
    frame_et_rec->GetYaxis()->SetRangeUser(2, 4e3);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 30);


    frame_et_rec->Draw("axis");


    gStat_PHENIX->SetMarkerStyle(28);
    gStat_PHENIX->SetMarkerSize(1.5);
    gStat_PHENIX->SetMarkerColor(kRed + 1);
    gStat_PHENIX->SetLineColor(kRed + 1);
    gStat_PHENIX->Draw("P");

    gSys_PHENIX->SetMarkerStyle(28);
    gSys_PHENIX->SetMarkerSize(1.5);
    gSys_PHENIX->SetMarkerColor(kRed);
    gSys_PHENIX->SetLineColor(kRed);
    gSys_PHENIX->SetFillColorAlpha(kRed, 0.25);
    gSys_PHENIX->Draw("2 same");

    g_syst->SetMarkerStyle(20);
    g_syst->SetMarkerColor(kAzure + 2);
    g_syst->SetLineColor(kAzure + 2);
    g_syst->SetFillColorAlpha(kAzure + 2, 0.5);

    g_syst->Draw("2 same");
    
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerColor(kBlack);
    h_data->SetLineColor(kBlack);

    h_data->Draw("same");

    myText(0.5 - 0.05, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5 - 0.05, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5 - 0.05, 0.80, 1, "|#eta^{#gamma}|<0.7, E_{T}^{iso, R=0.3} < 4 GeV", 0.04);
    myBoxTextAlpha(0.557 - 0.05, 0.25 + 0.5, 0.028, kAzure + 2, 0.5, kAzure + 2, 1, "");
    myMarkerLineText(0.55 - 0.05, 0.25 + 0.5, 1, kBlack, 20, kBlack, 1, Form("Data %s", reweightedstring.c_str()), 0.04, true);
    myText(0.50 - 0.05, 0.60, 1, "(|#eta^{#gamma}|<0.25, no E_{T}^{iso} requirement)", 0.04);
    //myText(0.50 - 0.05, 0.65, 1, "", 0.04);
    myBoxTextAlpha(0.557 - 0.05, 0.20 + 0.5, 0.028, kRed, 0.25, kRed, 1, "");
    myMarkerLineText(0.55 - 0.05, 0.20 + 0.5, 2, kRed + 1, 28, kRed + 1, 1, "PHENIX", 0.04, true);
    myText(0.55 - 0.05, 0.65, 1, "Phys. Rev. D 86 072008", 0.04);

    c2->SaveAs("figures/final_phenix.pdf");

    //-----------------------------------------------------------------
}
