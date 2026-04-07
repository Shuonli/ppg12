#include "plotcommon.h"

void plot_reweight(const char *infile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_no_unfolding_reweighting.root")
{
    init_plot();
    TFile *fin = new TFile(infile);

    // Use unfolded data (truth pT, efficiency-corrected) vs MC response-matrix prior (truth pT).
    // These are the correct inputs for reweighting the Bayesian unfolding prior.
    TH1D *h_unfolded = (TH1D *)fin->Get("h_unfold_sub_result");
    TH1D *h_mc_prior = (TH1D *)fin->Get("h_pT_truth_response_0");

    if (!h_unfolded || !h_mc_prior)
    {
        std::cerr << "ERROR: could not retrieve histograms from " << infile << std::endl;
        return;
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

    frame_et_truth->SetYTitle("dN/N/dp_{T}");
    frame_et_truth->GetYaxis()->SetRangeUser(1e-6, 1);
    frame_et_truth->GetXaxis()->SetRangeUser(8, 35);
    frame_et_truth->GetXaxis()->SetTitleOffset(0.98);
    frame_et_truth->GetYaxis()->SetTitleOffset(1.15);
    frame_et_truth->GetXaxis()->SetLabelSize(0.045);
    frame_et_truth->GetYaxis()->SetLabelSize(0.045);
    frame_et_truth->GetXaxis()->SetLabelOffset(2);
    frame_et_truth->GetXaxis()->SetNdivisions(505);
    frame_et_truth->SetXTitle("#it{p}_{T}^{#gamma, truth} [GeV]");
    frame_et_truth->Draw("axis");

    // normalize to unit area for shape comparison
    h_unfolded->Scale(1.0 / h_unfolded->Integral());
    h_unfolded->SetMarkerStyle(20);
    h_unfolded->SetMarkerColor(kBlack);
    h_unfolded->SetLineColor(kBlack);
    h_unfolded->Draw("same");

    h_mc_prior->Scale(1.0 / h_mc_prior->Integral());
    h_mc_prior->SetMarkerStyle(25);
    h_mc_prior->SetMarkerColor(kPink + 8);
    h_mc_prior->SetLineColor(kPink + 8);
    h_mc_prior->Draw("same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, "|#eta^{#gamma}|<0.7", 0.05);

    myMarkerLineText(0.25, 0.25, 1, kBlack, 20, kBlack, 1, "Data (unfolded)", 0.05, true);
    myMarkerLineText(0.25, 0.20, 1, kPink + 8, 25, kPink + 8, 1, "Pythia prior", 0.05, true);

    // --- ratio panel ---
    TPad *pad_2 = (TPad *)c1->cd(2);
    pad_2->SetPad(0, 0, 1, 0.4);
    pad_2->SetTopMargin(0.02);
    pad_2->SetLeftMargin(0.13);
    pad_2->SetBottomMargin(0.25);
    pad_2->SetRightMargin(0.08);

    TH1D *h_ratio = (TH1D *)h_unfolded->Clone("h_ratio");
    h_ratio->Divide(h_mc_prior);

    // [2/2] Pade approximant fit
    TF1 *fit_ratio = new TF1("fit_ratio", "([0] + [1]*x + [3]*x*x) / (1 + [2]*x + [4]*x*x)", 8, 35);
    fit_ratio->SetParameters(1.0, 0.1, 0.1, 0.1, 0.1);
    h_ratio->Fit(fit_ratio, "REM", "", 8, 35);
    h_ratio->Fit(fit_ratio, "REM", "", 8, 35);

    std::cout << "Pade fit parameters (copy into RecoEffCalculator_TTreeReader.C f_reweight):" << std::endl;
    std::cout << "  p0=" << fit_ratio->GetParameter(0)
              << " p1=" << fit_ratio->GetParameter(1)
              << " p2=" << fit_ratio->GetParameter(2)
              << " p3=" << fit_ratio->GetParameter(3)
              << " p4=" << fit_ratio->GetParameter(4) << std::endl;

    TH1F *frame_ratio = new TH1F("frame_ratio", "", 430, 0, 100);
    frame_ratio->SetYTitle("Data / Pythia prior");
    frame_ratio->GetYaxis()->SetNdivisions(506);
    frame_ratio->GetYaxis()->SetRangeUser(0.0, 2.0);
    frame_ratio->GetXaxis()->SetRangeUser(8, 35);
    frame_ratio->GetYaxis()->SetTitleOffset(frame_et_truth->GetYaxis()->GetTitleOffset() * 4 / 6.);
    frame_ratio->GetYaxis()->SetLabelOffset(frame_et_truth->GetYaxis()->GetLabelOffset() * 4 / 6.);
    frame_ratio->GetXaxis()->SetLabelSize(frame_et_truth->GetXaxis()->GetLabelSize() * 6 / 4.);
    frame_ratio->GetYaxis()->SetLabelSize(frame_et_truth->GetYaxis()->GetLabelSize() * 6 / 4.);
    frame_ratio->GetXaxis()->SetTitleSize(frame_et_truth->GetXaxis()->GetTitleSize() * 6 / 4.);
    frame_ratio->GetYaxis()->SetTitleSize(frame_et_truth->GetYaxis()->GetTitleSize() * 6 / 4.);
    frame_ratio->GetXaxis()->SetNdivisions(505);
    frame_ratio->SetXTitle("#it{p}_{T}^{#gamma, truth} [GeV]");
    frame_ratio->Draw("axis");

    TLine *line_one = new TLine(8, 1, 35, 1);
    line_one->SetLineStyle(2);
    line_one->SetLineColor(kGray + 1);
    line_one->Draw();

    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerColor(kBlack);
    h_ratio->SetLineColor(kBlack);
    h_ratio->Draw("same");

    fit_ratio->Draw("same");

    c1->SaveAs("figures/response_reweight.pdf");
}
