#include "plotcommon.h"

// Plot the MC purity correction factor g_mc_purity_fit_ratio
// (ratio of MC truth purity to ABCD-fit purity) vs E_T^gamma for a
// given variant suffix. Written to figures/mc_purity_correction_{tune}.pdf.
//
// Called per-variant by make_selection_plots.sh; referenced by
// make_selection_report.py under the "MC purity correction" plot group.

void plot_mc_purity_correction(const std::string &tune = "bdt_nom")
{
    init_plot();

    // g_mc_purity_fit_ratio and the smooth fit are produced in the MC-closure
    // pass of CalculatePhotonYield. Read them from the mc_purity_correction
    // variant's _mc file: same closure points as the nominal MC pass, and the
    // stored fit is the one actually applied as the correction (order set by
    // analysis.mc_purity_corr_fitorder, pol1 since 2026-09-08).
    TString infile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_mc_purity_correction_mc.root";
    TFile *fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "ERROR: cannot open input file: " << infile << std::endl;
        return;
    }

    TGraphErrors *g_ratio = (TGraphErrors *)fin->Get("g_mc_purity_fit_ratio");
    TF1 *f_corr = (TF1 *)fin->Get("f_mc_purity_corr_fit");
    if (!g_ratio)
    {
        std::cerr << "WARNING: g_mc_purity_fit_ratio not found in " << infile
                  << " -- skipping plot" << std::endl;
        fin->Close();
        return;
    }

    TCanvas *c = new TCanvas("c_mc_purity_corr", "", 700, 600);
    frame_et_rec->SetTitle(";#it{E}_{T}^{#gamma} [GeV];MC purity correction");
    // Reported range only: the 10-12 and 32-36 bins are unfolding under/overflow
    // and are not quoted, so they are not shown.
    frame_et_rec->GetXaxis()->SetRangeUser(12, 32);
    frame_et_rec->GetYaxis()->SetRangeUser(0.6, 1.4);
    frame_et_rec->Draw("axis");

    lineone->SetLineColor(kGray + 2);
    lineone->SetLineStyle(2);
    lineone->Draw("L same");

    g_ratio->SetMarkerStyle(20);
    g_ratio->SetMarkerSize(1.4);
    g_ratio->SetMarkerColor(kBlue + 1);
    g_ratio->SetLineColor(kBlue + 1);
    g_ratio->SetLineWidth(2);
    g_ratio->Draw("P same");

    if (f_corr)
    {
        f_corr->SetLineColor(kRed + 1);
        f_corr->SetLineWidth(3);
        f_corr->SetLineStyle(1);
        f_corr->SetNpx(400);
        f_corr->Draw("L same");
    }

    myText(0.22, 0.89, 1, strleg1.c_str(), 0.042, 0);
    myText(0.22, 0.84, 1, strleg2.c_str(), 0.042, 0);
    myText(0.88, 0.89, 1, strleg3.c_str(), 0.042, 1);
    myText(0.22, 0.78, 1, Form("tune: %s", tune.c_str()), 0.036, 0);

    TLegend *leg = new TLegend(0.22, 0.18, 0.62, 0.30);
    legStyle(leg, 0.15, 0.040);
    leg->AddEntry(g_ratio, "bin-by-bin: P_{truth}^{MC} / P_{ABCD}^{MC}", "pl");
    if (f_corr) leg->AddEntry(f_corr, Form("pol%d fit (applied as correction)", f_corr->GetNpar() - 1), "l");
    leg->Draw("same");

    TString outfile = Form("figures/mc_purity_correction_%s.pdf", tune.c_str());
    c->SaveAs(outfile);
    fin->Close();
}
