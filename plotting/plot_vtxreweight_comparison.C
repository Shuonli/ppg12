#include "plotcommon.h"

// =====================================================================
// plot_vtxreweight_comparison.C
//
// Quick 2-panel comparison: nominal vs "no vertex reweighting"
//   Left : g_mbd_eff      (TGraphAsymmErrors) overlay
//   Right: h_unfold_sub_result ratio (no-reweight / nominal)
//
// Inputs (both already exist):
//   efficiencytool/results/Photon_final_bdt_nom.root
//   efficiencytool/results/Photon_final_bdt_vtxreweight0.root
// Output:
//   plotting/figures/vtxreweight_comparison.pdf
// =====================================================================

void plot_vtxreweight_comparison()
{
    init_plot();

    const std::string fnom_path = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root";
    const std::string fnovtx_path = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_vtxreweight0.root";
    const std::string savePath = "/sphenix/user/shuhangli/ppg12/plotting/figures";
    gSystem->mkdir(savePath.c_str(), true);

    TFile *fnom = TFile::Open(fnom_path.c_str());
    TFile *fnovtx = TFile::Open(fnovtx_path.c_str());
    if (!fnom || fnom->IsZombie() || !fnovtx || fnovtx->IsZombie())
    {
        std::cerr << "Cannot open one of the input files" << std::endl;
        return;
    }

    // ---- MBD efficiency graphs ----
    TGraphAsymmErrors *g_mbd_nom = (TGraphAsymmErrors *)fnom->Get("g_mbd_eff");
    TGraphAsymmErrors *g_mbd_novtx = (TGraphAsymmErrors *)fnovtx->Get("g_mbd_eff");
    if (!g_mbd_nom || !g_mbd_novtx)
    {
        std::cerr << "Missing g_mbd_eff in one of the files" << std::endl;
        return;
    }

    g_mbd_nom->SetLineColor(kRed + 1);
    g_mbd_nom->SetMarkerColor(kRed + 1);
    g_mbd_nom->SetMarkerStyle(20);
    g_mbd_nom->SetMarkerSize(1.2);
    g_mbd_nom->SetLineWidth(2);
    g_mbd_nom->SetLineStyle(1);

    g_mbd_novtx->SetLineColor(kBlue + 1);
    g_mbd_novtx->SetMarkerColor(kBlue + 1);
    g_mbd_novtx->SetMarkerStyle(24);
    g_mbd_novtx->SetMarkerSize(1.2);
    g_mbd_novtx->SetLineWidth(2);
    g_mbd_novtx->SetLineStyle(2);

    // ---- Cross-section ratio (no-reweight / nominal) ----
    TH1D *h_xs_nom = (TH1D *)fnom->Get("h_unfold_sub_result");
    TH1D *h_xs_novtx = (TH1D *)fnovtx->Get("h_unfold_sub_result");
    if (!h_xs_nom || !h_xs_novtx)
    {
        std::cerr << "Missing h_unfold_sub_result in one of the files" << std::endl;
        return;
    }
    TH1D *h_ratio = (TH1D *)h_xs_novtx->Clone("h_xs_ratio_novtx_over_nom");
    h_ratio->Divide(h_xs_nom);
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetMarkerColor(kBlack);
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerSize(1.2);
    h_ratio->SetLineWidth(2);

    // ---- Canvas: 1 row x 2 columns ----
    TCanvas *c = new TCanvas("c_vtxreweight_cmp", "", 1300, 600);
    c->Divide(2, 1, 0.005, 0.005);

    const float xpos = 0.22, ypos = 0.885, dy = 0.054, fontsize = 0.042;
    const std::string nomLabel = "Nominal (vtx reweight ON)";
    const std::string novtxLabel = "No vertex reweighting";

    // ================================================================
    // Left panel: MBD efficiency overlay
    // ================================================================
    c->cd(1);
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.16);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.06);

    // Clone shared frame so we can customise Y-axis without polluting globals
    TH1F *frameL = (TH1F *)frame_et_rec->Clone("frame_mbd_eff");
    frameL->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
    frameL->SetYTitle("MBD #times vertex efficiency");
    frameL->GetYaxis()->SetRangeUser(0.0, 1.25);
    frameL->GetXaxis()->SetRangeUser(10, 40);
    frameL->Draw("axis");

    g_mbd_nom->Draw("PE same");
    g_mbd_novtx->Draw("PE same");

    TLegend *legL = new TLegend(0.35, 0.20, 0.93, 0.40);
    legStyle(legL, 0.20, 0.040);
    legL->AddEntry(g_mbd_nom, nomLabel.c_str(), "pl");
    legL->AddEntry(g_mbd_novtx, novtxLabel.c_str(), "pl");
    legL->Draw("same");

    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(), fontsize, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 0);
    myText(xpos, ypos - 3 * dy, 1, strleg4.c_str(), fontsize, 0);

    // ================================================================
    // Right panel: cross-section ratio
    // ================================================================
    c->cd(2);
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.16);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.06);

    TH1F *frameR = (TH1F *)frame_et_rec->Clone("frame_xs_ratio");
    frameR->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
    frameR->SetYTitle("#sigma(no vtx rwt) / #sigma(nominal)");
    frameR->GetYaxis()->SetRangeUser(0.8, 2.4);
    frameR->GetXaxis()->SetRangeUser(10, 40);
    frameR->Draw("axis");

    // Unity reference line
    lineone->Draw("L same");

    h_ratio->Draw("E1 same");

    TLegend *legR = new TLegend(0.35, 0.20, 0.93, 0.35);
    legStyle(legR, 0.20, 0.040);
    legR->AddEntry(h_ratio, "no reweight / nominal", "pl");
    legR->AddEntry(lineone, "unity", "l");
    legR->Draw("same");

    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(), fontsize, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 0);
    myText(xpos, ypos - 3 * dy, 1, strleg4.c_str(), fontsize, 0);

    // ---- Save ----
    const std::string outpdf = savePath + "/vtxreweight_comparison.pdf";
    c->SaveAs(outpdf.c_str());
    std::cout << "\nSaved: " << outpdf << std::endl;

    fnom->Close();
    fnovtx->Close();
}
