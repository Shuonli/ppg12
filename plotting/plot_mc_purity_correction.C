#include "plotcommon.h"

void plot_mc_purity_correction(
    const char *infile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom_mc.root",
    const char *graph_name = "g_mc_purity_fit_ratio")
{
    init_plot();

    TFile *fin = new TFile(infile, "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "ERROR: cannot open input file: " << infile << std::endl;
        return;
    }

    TGraphErrors *g_ratio = (TGraphErrors *)fin->Get(graph_name);
    if (!g_ratio)
    {
        std::cerr << "ERROR: graph '" << graph_name << "' not found in " << infile << std::endl;
        fin->Close();
        return;
    }

    TCanvas *c = new TCanvas("c_mc_purity_corr", "", 700, 600);
    frame_et_rec->SetTitle(";#it{E}_{T}^{#gamma} [GeV];MC purity correction");
    frame_et_rec->GetXaxis()->SetRangeUser(8, 35);
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
    myText(0.22, 0.78, 1, "MC truth / purity-fit ratio", 0.038, 0);

    TLegend *leg = new TLegend(0.22, 0.18, 0.62, 0.26);
    legStyle(leg, 0.15, 0.045);
    leg->AddEntry(g_ratio, "g_{mc}^{truth} / f_{purity}^{fit}", "pl");
    leg->Draw("same");

    c->SaveAs("figures/mc_purity_correction_vs_et.pdf");
    fin->Close();
}
