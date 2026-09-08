#include "plotcommon.h"

// =====================================================================
// plot_mc_purity_correction_compare.C
//
// Side-by-side NEW vs OLD overlay of the MC purity correction factor
// g_mc_purity_fit_ratio (bin-by-bin P_truth^MC / P_ABCD^MC) and its
// applied pol2 fit f_mc_purity_corr_fit, for two analysis variants.
//
// Built by reusing the single-variant plot_mc_purity_correction.C idiom
// (same objects, same frame, same sPHENIX labels) so the comparison is a
// drop-in overlay of two of those plots.
//
//   NEW = regenerated waveform-sim trees   (default: bdt_nom    _mc.root)
//   OLD = clean reference trees            (default: bdt_nomold _mc.root)
//
// Both objects live only in the MC-closure pass of CalculatePhotonYield
// (gated by `if (isMC)`), so they are read from the *_mc.root files.
//
// Output: figures/compare/mc_purity_correction_compare.pdf
//
// Usage:
//   root -l -b -q 'plot_mc_purity_correction_compare.C(
//       "/.../Photon_final_bdt_nom_mc.root",
//       "/.../Photon_final_bdt_nomold_mc.root")'
// =====================================================================

void plot_mc_purity_correction_compare(
    const std::string &newfile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom_mc.root",
    const std::string &oldfile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nomold_mc.root")
{
    init_plot();

    const std::string outDir = "/sphenix/user/shuhangli/ppg12/plotting/figures/compare";
    gSystem->mkdir(outDir.c_str(), true);

    TFile *fnew = TFile::Open(newfile.c_str(), "READ");
    if (!fnew || fnew->IsZombie())
    {
        std::cerr << "ERROR: cannot open NEW file: " << newfile << std::endl;
        return;
    }
    TFile *fold = TFile::Open(oldfile.c_str(), "READ");
    if (!fold || fold->IsZombie())
    {
        std::cerr << "ERROR: cannot open OLD file: " << oldfile << std::endl;
        return;
    }

    TGraphErrors *g_new = (TGraphErrors *)fnew->Get("g_mc_purity_fit_ratio");
    TGraphErrors *g_old = (TGraphErrors *)fold->Get("g_mc_purity_fit_ratio");
    TF1 *f_new = (TF1 *)fnew->Get("f_mc_purity_corr_fit");
    TF1 *f_old = (TF1 *)fold->Get("f_mc_purity_corr_fit");

    if (!g_new)
    {
        std::cerr << "ERROR: g_mc_purity_fit_ratio not found in NEW file " << newfile << std::endl;
        return;
    }
    if (!g_old)
    {
        std::cerr << "ERROR: g_mc_purity_fit_ratio not found in OLD file " << oldfile << std::endl;
        return;
    }

    TCanvas *c = new TCanvas("c_mc_purity_corr_cmp", "", 700, 600);
    frame_et_rec->SetTitle(";#it{E}_{T}^{#gamma} [GeV];MC purity correction");
    frame_et_rec->GetXaxis()->SetRangeUser(10, 36);
    frame_et_rec->GetYaxis()->SetRangeUser(0.6, 1.4);
    frame_et_rec->Draw("axis");

    lineone->SetLineColor(kGray + 2);
    lineone->SetLineStyle(2);
    lineone->Draw("L same");

    // NEW: blue filled circles + solid red... keep the two variants on a
    // consistent NEW=blue / OLD=red colour code (markers and fit share the
    // variant colour so the overlay reads cleanly without a colour key).
    g_new->SetMarkerStyle(20);
    g_new->SetMarkerSize(1.4);
    g_new->SetMarkerColor(kBlue + 1);
    g_new->SetLineColor(kBlue + 1);
    g_new->SetLineWidth(2);

    g_old->SetMarkerStyle(24);
    g_old->SetMarkerSize(1.4);
    g_old->SetMarkerColor(kRed + 1);
    g_old->SetLineColor(kRed + 1);
    g_old->SetLineWidth(2);

    if (f_new)
    {
        f_new->SetLineColor(kBlue + 1);
        f_new->SetLineWidth(3);
        f_new->SetLineStyle(1);
        f_new->SetNpx(400);
        f_new->Draw("L same");
    }
    if (f_old)
    {
        f_old->SetLineColor(kRed + 1);
        f_old->SetLineWidth(3);
        f_old->SetLineStyle(2);
        f_old->SetNpx(400);
        f_old->Draw("L same");
    }

    g_new->Draw("P same");
    g_old->Draw("P same");

    myText(0.22, 0.89, 1, strleg1.c_str(), 0.042, 0);
    myText(0.22, 0.84, 1, strleg2.c_str(), 0.042, 0);
    myText(0.88, 0.89, 1, strleg3.c_str(), 0.042, 1);

    TLegend *leg = new TLegend(0.22, 0.16, 0.66, 0.34);
    legStyle(leg, 0.15, 0.038);
    leg->AddEntry(g_new, "NEW (waveform sim): bin-by-bin", "pl");
    if (f_new) leg->AddEntry(f_new, "NEW pol2 fit", "l");
    leg->AddEntry(g_old, "OLD (clean ref): bin-by-bin", "pl");
    if (f_old) leg->AddEntry(f_old, "OLD pol2 fit", "l");
    leg->Draw("same");

    TString outfile = Form("%s/mc_purity_correction_compare.pdf", outDir.c_str());
    c->SaveAs(outfile);

    fnew->Close();
    fold->Close();
}
