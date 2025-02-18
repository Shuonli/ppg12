#include "plotcommon.h"

void plot_final()
{
    init_plot();

    bool reweighted = true;
    bool usepythia = true;

    std::string MCstring = usepythia ? "Pythia" : "JETPHOX";

    std::string reweightedstring = reweighted ? "w/ reweight" : "w/o reweight";

    std::string datastring = "data";
    std::string bg_MCstring = "Inclusive Sim";

    TFile *fin_data = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nr.root");

    if (reweighted)
        fin_data = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nom.root");

    TFile *fin_syst = new TFile("/sphenix/user/shuhangli/ppg12/plotting/rootFiles/syst_sum.root");
    TFile *fin_NLO = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX.root");
    TFile *fin_mc = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nom_mc.root");

    TH1F *h_data = (TH1F *)fin_data->Get("h_unfold_sub_result");
    TH1F *h_pythia = (TH1F *)fin_data->Get("h_truth_pT_0");
    TH1F *h_common_cluster_data = (TH1F *)fin_data->Get("h_common_cluster_0");
    TH1F *h_common_cluster_mc = (TH1F *)fin_mc->Get("h_common_cluster_0");
    TH1F *h_tight_iso_cluster_data = (TH1F *)fin_data->Get("h_tight_iso_cluster_0");
    TH1F *h_tight_iso_cluster_mc = (TH1F *)fin_mc->Get("h_tight_iso_cluster_0");

    TH1F *h_NLO = (TH1F *)fin_NLO->Get("h_truth_pT");
    if (usepythia)
        h_NLO = h_pythia;

    TH1F *h_data_NLO = (TH1F *)h_data->Clone("h_data_NLO");
    h_data_NLO->Divide(h_NLO);

    TH1F *h_syst_low = (TH1F *)fin_syst->Get("h_sum_low");
    TH1F *h_syst_high = (TH1F *)fin_syst->Get("h_sum_high");

    TH1F *h_syst_rel_low = (TH1F *)fin_syst->Get("h_sum_rel_low");
    TH1F *h_syst_rel_high = (TH1F *)fin_syst->Get("h_sum_rel_high");

    TGraphAsymmErrors *g_syst = new TGraphAsymmErrors(h_data);
    TGraphAsymmErrors *g_syst_rel = new TGraphAsymmErrors(h_data_NLO);

    for (int i = 0; i < h_data->GetNbinsX(); i++)
    {
        float xlowerror = g_syst->GetErrorXlow(i);
        float xuperror = g_syst->GetErrorXhigh(i);

        g_syst->SetPointError(i, xlowerror, xuperror, h_syst_low->GetBinContent(i + 1), h_syst_high->GetBinContent(i + 1));

        float y = g_syst_rel->GetY()[i];
        std::cout << "y: " << y << std::endl;
        g_syst_rel->SetPointError(i, xlowerror, xuperror, y * h_syst_rel_low->GetBinContent(i + 1), y * h_syst_rel_high->GetBinContent(i + 1));

        std::cout << "h_syst_rel_low->GetBinContent(i + 1): " << h_syst_rel_low->GetBinContent(i + 1) << std::endl;
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

    frame_et_rec->SetYTitle("d#sigma/dE_{T} [pb/GeV]");
    frame_et_rec->GetYaxis()->SetRangeUser(0.1, 2e3);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 30);

    frame_et_rec->GetXaxis()->SetTitleOffset(0.98);
    frame_et_rec->GetYaxis()->SetTitleOffset(1.15);
    frame_et_rec->GetXaxis()->SetLabelSize(0.045);
    frame_et_rec->GetYaxis()->SetLabelSize(0.045);
    frame_et_rec->GetXaxis()->SetLabelOffset(2);
    frame_et_rec->GetXaxis()->CenterTitle();
    frame_et_rec->GetYaxis()->CenterTitle();
    frame_et_rec->GetXaxis()->SetNdivisions(505);

    frame_et_rec->Draw("axis");

    g_syst->SetMarkerStyle(20);
    g_syst->SetMarkerColor(kAzure + 2);
    g_syst->SetLineColor(kAzure + 2);
    g_syst->SetFillColorAlpha(kAzure + 2, 0.3);

    g_syst->Draw("5 same");

    h_NLO->SetMarkerStyle(24);
    h_NLO->SetMarkerColor(kPink + 8);
    h_NLO->SetLineColor(kPink + 8);

    h_NLO->Draw("same");

    h_data->SetMarkerStyle(20);
    h_data->SetMarkerColor(kBlack);
    h_data->SetLineColor(kBlack);

    h_data->Draw("same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);

    myMarkerLineText(0.25, 0.25, 1, kBlack, 20, kBlack, 1, Form("data %s", reweightedstring.c_str()), 0.05, true);
    myMarkerLineText(0.25, 0.20, 1, kPink + 8, 24, kPink + 8, 1, MCstring.c_str(), 0.05, true);

    TPad *pad_2 = (TPad *)c1->cd(2);
    pad_2->SetPad(0, 0, 1, 0.4);
    pad_2->SetTopMargin(0.02);
    pad_2->SetLeftMargin(0.13);
    pad_2->SetBottomMargin(0.25);
    pad_2->SetRightMargin(0.08);

    frame_et_truth->SetYTitle("data/theory");
    frame_et_truth->GetYaxis()->SetNdivisions(506);
    frame_et_truth->GetYaxis()->SetRangeUser(0.2, 1.5);
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
    g_syst_rel->SetFillColorAlpha(kAzure + 2, 0.3);

    lineone->Draw("L");

    g_syst_rel->Draw("5 same");

    h_data_NLO->SetMarkerStyle(20);
    h_data_NLO->SetMarkerColor(kBlack);
    h_data_NLO->SetLineColor(kBlack);

    h_data_NLO->Draw("same");

    std::string outputname = "figures/final_" + MCstring + (reweighted ? "_rewe" : "") + ".pdf";

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
    frame_et_truth->SetXTitle("#it{E}_{T}^{iso} [GeV]");
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
    frame_et_truth->SetXTitle("#it{E}_{T}^{iso} [GeV]");
    frame_et_truth->Draw("axis");

    lineone->Draw("L");

    TH1F *h_tight_iso_cluster_data_ratio = (TH1F *)h_tight_iso_cluster_data->Clone("h_tight_iso_cluster_data_ratio");
    h_tight_iso_cluster_data_ratio->Divide(h_tight_iso_cluster_mc);

    h_tight_iso_cluster_data_ratio->SetMarkerStyle(20);
    h_tight_iso_cluster_data_ratio->SetMarkerColor(kBlack);
    h_tight_iso_cluster_data_ratio->SetLineColor(kBlack);

    h_tight_iso_cluster_data_ratio->Draw("same");

    c1->SaveAs("figures/final_tight_iso_cluster.pdf");






}
