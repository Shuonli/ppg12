// plot_paper_efficiency.C -- paper Fig.~\ref{fig:eff_photon}
//
// Lifts the c6 canvas (L199-236) from plotting/plot_efficiency.C: the
// stacked photon-efficiency panel showing eps_reco, eps_reco x eps_ID, and
// eps_reco x eps_ID x eps_iso vs ETg. The single-component panels (c1-c4)
// and the MBD-eff / MBD-fraction panels (c5, c_mbdfrac) are not used in
// the paper and are intentionally NOT carried over.
//
// Inputs:
//   MC_efficiency_<tune>.root  -- TEfficiency objects eff_reco_eta_0,
//                                 eff_id_eta_0, eff_all_eta_0
//   Photon_final_<tune>.root   -- (only used for the x-binning of g_mbd_eff)
//
// Output: PPG12-Paper/figures/eff_photon_<tune>.pdf

#include <TFile.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <cmath>
#include <string>

#include "paper_style.h"

namespace
{
const int   kCol[]    = {kBlack,  kPink + 5, kGreen - 2};
const int   kMkStyle[]= {21,      20,        34};
const float kMkSize[] = {1.2f,    1.2f,      1.6f};
}  // namespace

void plot_paper_efficiency_yj(const std::string &tune = "bdt_nom")
{
    paper_init();

    const std::string resdir = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    TFile *fmc   = TFile::Open(Form("%s/MC_efficiency_%s.root", resdir.c_str(), tune.c_str()));
    TFile *fdata = TFile::Open(Form("%s/Photon_final_%s.root",  resdir.c_str(), tune.c_str()));
    if (!fmc || fmc->IsZombie() || !fdata || fdata->IsZombie())
    {
        std::cerr << "[plot_paper_efficiency] missing MC_efficiency / Photon_final for tune=" << tune << std::endl;
        return;
    }

    TEfficiency *eff_reco = (TEfficiency *)fmc->Get("eff_reco_eta_0");
    TEfficiency *eff_id   = (TEfficiency *)fmc->Get("eff_id_eta_0");
    TEfficiency *eff_all  = (TEfficiency *)fmc->Get("eff_all_eta_0");
    TGraphAsymmErrors *g_mbd_eff = (TGraphAsymmErrors *)fdata->Get("g_mbd_eff");
    if (!eff_reco || !eff_id || !eff_all || !g_mbd_eff)
    {
        std::cerr << "[plot_paper_efficiency] missing efficiency object(s)." << std::endl;
        return;
    }

    // Build eps_reco x eps_ID by hand using the per-bin product, with the
    // fractional uncertainties summed in quadrature.
    TGraphAsymmErrors *g_reco_id_product = (TGraphAsymmErrors *)g_mbd_eff->Clone("g_reco_id_product");
    for (int i = 0; i < eff_reco->GetTotalHistogram()->GetNbinsX(); ++i)
    {
        double reco     = eff_reco->GetEfficiency(i + 1);
        double id       = eff_id  ->GetEfficiency(i + 1);
        double reco_err = eff_reco->GetEfficiencyErrorLow(i + 1);
        double id_err   = eff_id  ->GetEfficiencyErrorLow(i + 1);
        double product  = reco * id;
        // Quadrature of relative errors (kept for future re-enabling, not
        // drawn on the paper figure -- markers only).
        (void) std::sqrt(std::pow(reco_err / std::max(reco, 1e-9), 2) +
                         std::pow(id_err   / std::max(id,   1e-9), 2));
        float x = g_mbd_eff->GetX()[i];
        g_reco_id_product->SetPoint(i, x, product);
    }

    TCanvas *c6 = new TCanvas("c6_paper_eff", "", 600, 600);
    frame_et_truth->SetYTitle("Efficiency");
    // y from 0 so the lower-left legend sits comfortably below the
    // lowest curve (eps_reco*ID*iso ~ 0.38 -> NDC ~0.39, legend top
    // NDC ~0.30 -> safe gap). x covers the full unfolding range
    // (10 < ETg < 36 GeV), wider than the reported analysis range
    // (12-32) so the truth-spectrum overflow bins are visible.
    frame_et_truth->GetYaxis()->SetRangeUser(0.0, 1.15);
    frame_et_truth->GetXaxis()->SetRangeUser(10, 36);
    frame_et_truth->Draw("axis");

    eff_reco->SetMarkerColor(kCol[0]);
    eff_reco->SetMarkerStyle(kMkStyle[0]);
    eff_reco->SetMarkerSize(kMkSize[0]);
    eff_reco->SetLineColor(kCol[0]);
    eff_reco->Draw("same");

    g_reco_id_product->SetMarkerColor(kCol[1]);
    g_reco_id_product->SetMarkerStyle(kMkStyle[1]);
    g_reco_id_product->SetMarkerSize(kMkSize[1]);
    g_reco_id_product->SetLineColor(kCol[1]);
    g_reco_id_product->Draw("same p");

    eff_all->SetMarkerColor(kCol[2]);
    eff_all->SetMarkerStyle(kMkStyle[2]);
    eff_all->SetMarkerSize(kMkSize[2]);
    eff_all->SetLineColor(kCol[2]);
    eff_all->SetLineWidth(2);
    eff_all->Draw("same");

    float xpos(0.2), xpos2(0.915), ypos(0.885),
          dy(0.054), fontsize(0.046), fontsize1(0.048);
    myText(xpos,  ypos - 0 * dy, 1, strleg1.c_str(), fontsize1, 0);
    myText(xpos,  ypos - 1 * dy, 1, strleg2.c_str(), fontsize,  0);
    myText(xpos2, ypos - 0 * dy, 1, strMC.c_str(),   fontsize,  1);
    myText(xpos2, ypos - 1 * dy, 1, strleg3.c_str(), fontsize,  1);

    // Legend at lower-left in the empty band below the lowest curve
    // (eps_reco*ID*iso ~ 0.38 -> NDC ~0.27 with the new y-range), away
    // from data points across the whole 12-32 GeV span.
    TLegend *l1 = new TLegend(0.22, 0.25, 0.62, 0.50);
    legStyle(l1, 0.20, 0.058);
    l1->AddEntry(eff_reco,           "#varepsilon_{reco}",                                                 "pl");
    l1->AddEntry(g_reco_id_product,  "#varepsilon_{reco}#times#varepsilon_{ID}",                            "pl");
    l1->AddEntry(eff_all,            "#varepsilon_{reco}#times#varepsilon_{ID}#times#varepsilon_{iso}",     "pl");
    l1->Draw("same");

    // Canonical paper name (no tune suffix) for main.tex stability.
    c6->SaveAs(Form("%s/eff_photon.pdf", paper_savepath().c_str()));
}
