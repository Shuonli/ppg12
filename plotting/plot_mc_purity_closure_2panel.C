#include "plotcommon.h"

// Two-panel MC-purity-closure figure for the analysis note (fig:syst_mc_purity).
//   LEFT  : actual closure correction C_MC(E_T^g) = P_truth^MC / P_ABCD^MC
//           (bin-by-bin points + pol2 fit + unity line), from Photon_final_{tune}_mc.root
//   RIGHT : derived systematic = signed (sigma_var - sigma_nom)/sigma_nom per E_T bin,
//           from h_unfold_sub_result in the mc_purity_correction variant vs nominal.
// Both panels share one SetsPhenixStyle() and a common 12-32 GeV reported range.
// Written to figures/syst_mc_purity_2panel_{tune}.pdf.

void plot_mc_purity_closure_2panel(const std::string &tune = "bdt_nom")
{
    init_plot();
    gStyle->SetOptStat(0);

    const double xlo = 12.0, xhi = 32.0;

    // ---------------- LEFT: C_MC actual values ----------------
    TString fmc = Form("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_%s_mc.root", tune.c_str());
    TFile *finL = TFile::Open(fmc, "READ");
    if (!finL || finL->IsZombie()) { std::cerr << "ERROR: cannot open " << fmc << std::endl; return; }
    TGraphErrors *g_ratio = (TGraphErrors *)finL->Get("g_mc_purity_fit_ratio");
    TF1 *f_corr = (TF1 *)finL->Get("f_mc_purity_corr_fit");
    if (!g_ratio) { std::cerr << "ERROR: g_mc_purity_fit_ratio not found" << std::endl; return; }

    // ---------------- RIGHT: derived systematic ----------------
    TString fnom = Form("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_%s.root", tune.c_str());
    TString fvar = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_mc_purity_correction.root";
    TFile *finN = TFile::Open(fnom, "READ");
    TFile *finV = TFile::Open(fvar, "READ");
    if (!finN || finN->IsZombie() || !finV || finV->IsZombie()) { std::cerr << "ERROR: cannot open nom/var" << std::endl; return; }
    TH1 *hnom = (TH1 *)finN->Get("h_unfold_sub_result");
    TH1 *hvar = (TH1 *)finV->Get("h_unfold_sub_result");
    if (!hnom || !hvar) { std::cerr << "ERROR: h_unfold_sub_result missing" << std::endl; return; }
    TH1D *hdev = (TH1D *)hnom->Clone("h_mc_purity_dev_rel");
    hdev->SetDirectory(0);
    for (int i = 1; i <= hnom->GetNbinsX(); ++i)
    {
        double n = hnom->GetBinContent(i);
        hdev->SetBinContent(i, (n != 0.0) ? (hvar->GetBinContent(i) - n) / n : 0.0);
        hdev->SetBinError(i, 0.0);
    }

    // ---------------- Canvas: two equal, top-aligned pads ----------------
    TCanvas *c = new TCanvas("c_mcpc2", "", 1000, 470);
    TPad *pL = new TPad("pL", "", 0.00, 0.0, 0.50, 1.0);
    TPad *pR = new TPad("pR", "", 0.50, 0.0, 1.00, 1.0);
    for (TPad *p : {pL, pR})
    {
        p->SetLeftMargin(0.17);
        p->SetRightMargin(0.04);
        p->SetBottomMargin(0.14);
        p->SetTopMargin(0.05);
        p->Draw();
    }

    // ----- LEFT panel -----
    pL->cd();
    TH1F *frL = new TH1F("frL", "", 100, 8, 40);
    frL->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");
    frL->SetYTitle("MC purity correction #it{C}_{MC}");
    frL->GetXaxis()->SetRangeUser(xlo, xhi);
    frL->GetYaxis()->SetRangeUser(0.6, 1.45);
    frL->GetYaxis()->SetTitleOffset(1.6);
    frL->Draw("axis");

    TLine *lu = new TLine(xlo, 1.0, xhi, 1.0);
    lu->SetLineColor(kGray + 2); lu->SetLineStyle(2); lu->Draw("L same");

    if (f_corr)
    {
        f_corr->SetLineColor(kRed + 1); f_corr->SetLineWidth(3);
        f_corr->SetNpx(400); f_corr->SetRange(xlo, xhi); f_corr->Draw("L same");
    }
    g_ratio->SetMarkerStyle(20); g_ratio->SetMarkerSize(1.3);
    g_ratio->SetMarkerColor(kBlue + 1); g_ratio->SetLineColor(kBlue + 1); g_ratio->SetLineWidth(2);
    g_ratio->Draw("P same");

    myText(0.22, 0.90, 1, strleg1.c_str(), 0.045, 0);
    myText(0.22, 0.845, 1, strleg2.c_str(), 0.045, 0);
    myText(0.94, 0.90, 1, strleg3.c_str(), 0.045, 1);

    TLegend *legL = new TLegend(0.22, 0.19, 0.70, 0.32);
    legStyle(legL, 0.15, 0.042);
    legL->AddEntry(g_ratio, "bin-by-bin ratio", "p");
    if (f_corr) legL->AddEntry(f_corr, "pol2 fit (applied as correction)", "l");
    legL->Draw("same");

    // ----- RIGHT panel -----
    pR->cd();
    TH1F *frR = new TH1F("frR", "", 100, 8, 40);
    frR->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");
    frR->SetYTitle("Relative difference #Delta#sigma/#sigma");
    frR->GetXaxis()->SetRangeUser(xlo, xhi);
    frR->GetYaxis()->SetRangeUser(-0.10, 0.22);
    frR->GetYaxis()->SetTitleOffset(1.6);
    frR->Draw("axis");

    TLine *lz = new TLine(xlo, 0.0, xhi, 0.0);
    lz->SetLineColor(kGray + 2); lz->SetLineStyle(2); lz->Draw("L same");

    hdev->SetLineColor(kRed + 1); hdev->SetLineWidth(3);
    hdev->GetXaxis()->SetRangeUser(xlo, xhi);
    hdev->Draw("HIST same ][");

    myText(0.94, 0.90, 1, strleg1.c_str(), 0.045, 1);
    myText(0.94, 0.845, 1, strleg2.c_str(), 0.045, 1);

    myText(0.22, 0.315, 1, "one-sided; |#delta| symmetrized in budget", 0.036, 0);
    TLegend *legR = new TLegend(0.22, 0.20, 0.74, 0.28);
    legStyle(legR, 0.15, 0.042);
    legR->AddEntry(hdev, "MC closure correction", "l");
    legR->Draw("same");

    TString outfile = Form("figures/syst_mc_purity_2panel_%s.pdf", tune.c_str());
    c->SaveAs(outfile);
    std::cout << "WROTE " << outfile << std::endl;

    // per-bin dump for verification
    std::cout << "=== per-bin (E_T center | C_MC point | dev_rel %) ===" << std::endl;
    for (int i = 1; i <= hnom->GetNbinsX(); ++i)
    {
        double cx = hnom->GetBinCenter(i);
        if (cx < xlo || cx > xhi) continue;
        std::cout << Form("  %5.1f  |  dev_rel = %+6.2f%%", cx, 100.0 * hdev->GetBinContent(i)) << std::endl;
    }
    finL->Close(); finN->Close(); finV->Close();
}
