#include "plotcommon.h"

void plot_reweight()
{
    init_plot();
    TFile *fin_data = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nom.root");

    TH1F* h_data_sub_leak = (TH1F*)fin_data->Get("h_data_sub_leak_copy");
    TH1F* h_tight_iso_cluster_signal = (TH1F*)fin_data->Get("h_tight_iso_cluster_signal_copy");

    TH1F* h_data_sub_leak_ratio = (TH1F*)h_data_sub_leak->Clone("h_data_sub_leak_ratio");
    h_data_sub_leak_ratio->Divide(h_tight_iso_cluster_signal);


    TCanvas *c1 = new TCanvas("can", "", 800, 889);
    c1->Divide(1, 2);

    TPad *pad_1 = (TPad *)c1->cd(1);
    pad_1->SetPad(0, 0.4, 1, 1);
    pad_1->SetTopMargin(0.05);
    pad_1->SetLeftMargin(0.13);
    pad_1->SetBottomMargin(0.03);
    pad_1->SetRightMargin(0.08);
    pad_1->SetLogy();

    frame_et_rec->SetYTitle("dN/N/dE_{T}");
    frame_et_rec->GetYaxis()->SetRangeUser(1e-6, 1);
    frame_et_rec->GetXaxis()->SetRangeUser(8, 35);

    frame_et_rec->GetXaxis()->SetTitleOffset(0.98);
    frame_et_rec->GetYaxis()->SetTitleOffset(1.15);
    frame_et_rec->GetXaxis()->SetLabelSize(0.045);
    frame_et_rec->GetYaxis()->SetLabelSize(0.045);
    frame_et_rec->GetXaxis()->SetLabelOffset(2);
    //frame_et_rec->GetXaxis()->CenterTitle();
    //frame_et_rec->GetYaxis()->CenterTitle();
    frame_et_rec->GetXaxis()->SetNdivisions(505);

    frame_et_rec->Draw("axis");

    //loop over h_data_sub_leak and scale it by the sum of the bins
    float sum = 0;
    for (int i = 1; i <= h_data_sub_leak->GetNbinsX(); i++)
    {
        sum += h_data_sub_leak->GetBinContent(i);
    }
    h_data_sub_leak->Scale(1.0 / sum);
    h_data_sub_leak->SetMarkerStyle(20);
    h_data_sub_leak->SetMarkerColor(kBlack);
    h_data_sub_leak->SetLineColor(kBlack);
    h_data_sub_leak->Draw("same");
    float sum_signal = 0;
    for (int i = 1; i <= h_tight_iso_cluster_signal->GetNbinsX(); i++)
    {
        sum_signal += h_tight_iso_cluster_signal->GetBinContent(i);
    }
    h_tight_iso_cluster_signal->Scale(1.0 / sum_signal);
    h_tight_iso_cluster_signal->SetMarkerStyle(24);
    h_tight_iso_cluster_signal->SetMarkerColor(kPink + 8);
    h_tight_iso_cluster_signal->SetLineColor(kPink + 8);
    h_tight_iso_cluster_signal->Draw("same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, "|#eta^{#gamma}|<0.7", 0.05);

    myMarkerLineText(0.25, 0.25, 1, kBlack, 20, kBlack, 1, "Data (purity corrected)", 0.05, true);
    myMarkerLineText(0.25, 0.20, 1, kPink + 8, 25, kPink + 8, 1, "Pythia signal reco", 0.05, true);


    TPad *pad_2 = (TPad *)c1->cd(2);
    pad_2->SetPad(0, 0, 1, 0.4);
    pad_2->SetTopMargin(0.02);
    pad_2->SetLeftMargin(0.13);
    pad_2->SetBottomMargin(0.25);
    pad_2->SetRightMargin(0.08);

    frame_et_truth->SetYTitle("Reweighting factor");
    frame_et_truth->GetYaxis()->SetNdivisions(506);
    frame_et_truth->GetYaxis()->SetRangeUser(0.1, 1.2);
    frame_et_truth->GetXaxis()->SetRangeUser(8, 35);
    frame_et_truth->GetYaxis()->SetTitleOffset(frame_et_rec->GetYaxis()->GetTitleOffset() * 4 / 6.);
    frame_et_truth->GetYaxis()->SetLabelOffset(frame_et_rec->GetYaxis()->GetLabelOffset() * 4 / 6.);
    frame_et_truth->GetXaxis()->SetLabelSize(frame_et_rec->GetXaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetYaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetTitleSize(frame_et_rec->GetXaxis()->GetTitleSize() * 6 / 4.);
    frame_et_truth->GetYaxis()->SetTitleSize(frame_et_rec->GetYaxis()->GetTitleSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetNdivisions(505);
    frame_et_truth->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");
    frame_et_truth->Draw("axis");

    // polynomial 3 fit for ratio
    TF1 *fit_ratio = new TF1("fit_ratio", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 8, 35);
    fit_ratio->SetParameter(0, 1);
    fit_ratio->SetParameter(1, 0);
    fit_ratio->SetParameter(2, 0);
    fit_ratio->SetParameter(3, 0);
    h_data_sub_leak_ratio->Fit(fit_ratio, "REM", "", 8, 35);

    //print the fit parameters
    std::cout << "fit parameters: " << fit_ratio->GetParameter(0) << " " << fit_ratio->GetParameter(1) << " " << fit_ratio->GetParameter(2) << " " << fit_ratio->GetParameter(3) << std::endl;



    h_data_sub_leak_ratio->SetMarkerStyle(20);
    h_data_sub_leak_ratio->SetMarkerColor(kBlack);
    h_data_sub_leak_ratio->SetLineColor(kBlack);
    h_data_sub_leak_ratio->Draw("same");

    c1->SaveAs("figures/response_reweight.pdf");

}