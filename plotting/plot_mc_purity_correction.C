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

    // g_mc_purity_fit_ratio is produced only in the MC-closure pass of
    // CalculatePhotonYield (gated by `if (isMC)` at line ~706). The jet MC
    // serves as the inclusive-MC reference for this analysis, so the ratio
    // is read from the _mc.root file.
    TString infile = Form("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_%s_mc.root", tune.c_str());
    TFile *fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "ERROR: cannot open input file: " << infile << std::endl;
        return;
    }

    TGraphErrors *g_ratio = (TGraphErrors *)fin->Get("g_mc_purity_fit_ratio");
    if (!g_ratio)
    {
        std::cerr << "WARNING: g_mc_purity_fit_ratio not found in " << infile
                  << " -- skipping plot" << std::endl;
        fin->Close();
        return;
    }

    TCanvas *c = new TCanvas("c_mc_purity_corr", "", 700, 600);
    frame_et_rec->SetTitle(";#it{E}_{T}^{#gamma} [GeV];MC purity correction");
    frame_et_rec->GetXaxis()->SetRangeUser(10, 35);
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

    myText(0.22, 0.89, 1, strleg1.c_str(), 0.042, 0);
    myText(0.22, 0.84, 1, strleg2.c_str(), 0.042, 0);
    myText(0.88, 0.89, 1, strleg3.c_str(), 0.042, 1);
    myText(0.22, 0.78, 1, Form("tune: %s", tune.c_str()), 0.036, 0);

    TLegend *leg = new TLegend(0.22, 0.18, 0.62, 0.26);
    legStyle(leg, 0.15, 0.045);
    leg->AddEntry(g_ratio, "g_{mc}^{truth} / f_{purity}^{fit}", "pl");
    leg->Draw("same");

    TString outfile = Form("figures/mc_purity_correction_%s.pdf", tune.c_str());
    c->SaveAs(outfile);
    fin->Close();
}
