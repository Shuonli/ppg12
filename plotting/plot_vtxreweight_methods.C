#include "plotcommon.h"

// =====================================================================
// plot_vtxreweight_methods.C
//
// 3-way comparison: reco-vertex reweight (OLD, nominal) vs
//                   no reweight vs
//                   truth-vertex reweight (NEW, recommended)
//   Left : g_mbd_eff      (TGraphAsymmErrors) overlay
//   Right: h_unfold_sub_result ratios (each variant / OLD)
//
// Inputs (should all exist after pipeline run completes):
//   efficiencytool/results/Photon_final_bdt_nom.root            (OLD reco-vtx reweight)
//   efficiencytool/results/Photon_final_bdt_vtxreweight0.root   (no reweight)
//   efficiencytool/results/Photon_final_bdt_truthvtxreweight.root (truth-vtx reweight)
// Output:
//   plotting/figures/vtxreweight_methods_comparison.pdf
// =====================================================================

namespace {
    TGraphAsymmErrors *getMbd(TFile *f) { return (TGraphAsymmErrors *)f->Get("g_mbd_eff"); }
    TH1D *getXs(TFile *f) { return (TH1D *)f->Get("h_unfold_sub_result"); }
    void styleG(TGraphAsymmErrors *g, int color, int marker, int lineStyle) {
        g->SetLineColor(color); g->SetMarkerColor(color);
        g->SetMarkerStyle(marker); g->SetMarkerSize(1.2);
        g->SetLineWidth(2); g->SetLineStyle(lineStyle);
    }
    void styleH(TH1D *h, int color, int marker, int lineStyle) {
        h->SetLineColor(color); h->SetMarkerColor(color);
        h->SetMarkerStyle(marker); h->SetMarkerSize(1.2);
        h->SetLineWidth(2); h->SetLineStyle(lineStyle);
    }
}

void plot_vtxreweight_methods()
{
    init_plot();

    const std::string base = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/";
    const std::string fold_path   = base + "Photon_final_bdt_nom.root";               // reco-vtx reweight (OLD)
    const std::string fno_path    = base + "Photon_final_bdt_vtxreweight0.root";      // no reweight
    const std::string ftruth_path = base + "Photon_final_bdt_truthvtxreweight.root";  // truth-vtx reweight (NEW)
    const std::string savePath    = "/sphenix/user/shuhangli/ppg12/plotting/figures";
    gSystem->mkdir(savePath.c_str(), true);

    TFile *fold   = TFile::Open(fold_path.c_str());
    TFile *fno    = TFile::Open(fno_path.c_str());
    TFile *ftruth = TFile::Open(ftruth_path.c_str());
    if (!fold   || fold->IsZombie())   { std::cerr << "Cannot open OLD file "   << fold_path   << std::endl; return; }
    if (!fno    || fno->IsZombie())    { std::cerr << "Cannot open NO file "    << fno_path    << std::endl; return; }
    if (!ftruth || ftruth->IsZombie()) { std::cerr << "Cannot open TRUTH file " << ftruth_path << std::endl; return; }

    // ---- MBD efficiency graphs ----
    TGraphAsymmErrors *g_old   = getMbd(fold);
    TGraphAsymmErrors *g_no    = getMbd(fno);
    TGraphAsymmErrors *g_truth = getMbd(ftruth);
    if (!g_old || !g_no || !g_truth) { std::cerr << "Missing g_mbd_eff" << std::endl; return; }

    styleG(g_old,   kRed     + 1, 20, 1);  // solid red circle
    styleG(g_no,    kBlue    + 1, 24, 2);  // dashed blue open circle
    styleG(g_truth, kGreen   + 2, 22, 1);  // solid green triangle

    // ---- Cross-section ratios (each variant / OLD) ----
    TH1D *h_old   = getXs(fold);
    TH1D *h_no    = getXs(fno);
    TH1D *h_truth = getXs(ftruth);
    if (!h_old || !h_no || !h_truth) { std::cerr << "Missing h_unfold_sub_result" << std::endl; return; }

    TH1D *r_no    = (TH1D *)h_no   ->Clone("h_xs_ratio_no_over_old");    r_no   ->Divide(h_old);
    TH1D *r_truth = (TH1D *)h_truth->Clone("h_xs_ratio_truth_over_old"); r_truth->Divide(h_old);
    styleH(r_no,    kBlue  + 1, 24, 2);
    styleH(r_truth, kGreen + 2, 22, 1);

    TCanvas *c = new TCanvas("c_vtxreweight_methods", "", 1300, 600);
    c->Divide(2, 1, 0.005, 0.005);

    const float xpos = 0.22, ypos = 0.885, dy = 0.054, fontsize = 0.042;
    const std::string labOld   = "Reco-vtx reweight (OLD nominal)";
    const std::string labNo    = "No reweight";
    const std::string labTruth = "Truth-vtx reweight (NEW)";

    // --- Left: MBD efficiency overlay
    c->cd(1);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.16);
    gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.06);

    TH1F *frameL = (TH1F *)frame_et_rec->Clone("frame_mbd_eff_methods");
    frameL->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
    frameL->SetYTitle("MBD #times vertex efficiency");
    frameL->GetYaxis()->SetRangeUser(0.0, 1.25);
    frameL->GetXaxis()->SetRangeUser(10, 40);
    frameL->Draw("axis");

    g_old  ->Draw("PE same");
    g_no   ->Draw("PE same");
    g_truth->Draw("PE same");

    TLegend *legL = new TLegend(0.30, 0.18, 0.95, 0.42);
    legStyle(legL, 0.20, 0.038);
    legL->AddEntry(g_old,   labOld.c_str(),   "pl");
    legL->AddEntry(g_no,    labNo.c_str(),    "pl");
    legL->AddEntry(g_truth, labTruth.c_str(), "pl");
    legL->Draw("same");

    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),   fontsize, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(), fontsize, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),   fontsize, 0);
    myText(xpos, ypos - 3 * dy, 1, strleg4.c_str(),   fontsize, 0);

    // --- Right: cross-section ratios to OLD
    c->cd(2);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.16);
    gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.06);

    TH1F *frameR = (TH1F *)frame_et_rec->Clone("frame_xs_ratio_methods");
    frameR->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
    frameR->SetYTitle("#sigma(variant) / #sigma(OLD)");
    frameR->GetYaxis()->SetRangeUser(0.7, 2.4);
    frameR->GetXaxis()->SetRangeUser(10, 40);
    frameR->Draw("axis");

    lineone->Draw("L same");
    r_no   ->Draw("E1 same");
    r_truth->Draw("E1 same");

    TLegend *legR = new TLegend(0.30, 0.18, 0.95, 0.42);
    legStyle(legR, 0.20, 0.038);
    legR->AddEntry(r_no,    "No reweight / OLD",    "pl");
    legR->AddEntry(r_truth, "Truth-vtx rwt / OLD",  "pl");
    legR->AddEntry(lineone, "unity",                "l");
    legR->Draw("same");

    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),   fontsize, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(), fontsize, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),   fontsize, 0);
    myText(xpos, ypos - 3 * dy, 1, strleg4.c_str(),   fontsize, 0);

    const std::string outpdf = savePath + "/vtxreweight_methods_comparison.pdf";
    c->SaveAs(outpdf.c_str());
    std::cout << "\nSaved: " << outpdf << std::endl;

    fold->Close(); fno->Close(); ftruth->Close();
}
