// plot_purity_pade_vs_erf.C
//
// Overlay the nominal Padé [1/1] purity fit and the alternative
// error-function fit on the leakage-corrected data purity vs E_T.
// Used as the "actual fit" diagnostic figure for Section 4 of the
// PPG12 analysis note (reviewer comment 25/29).
//
// Output: figures/purity_fit_compare_pade_erf.pdf

#include "plotcommon.h"

void plot_purity_pade_vs_erf()
{
    init_plot();

    const std::string inFile  = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/"
                                "ppg12/efficiencytool/results/"
                                "Photon_final_bdt_nom.root";
    const std::string outFile = "figures/purity_fit_compare_pade_erf.pdf";

    TFile *f = TFile::Open(inFile.c_str(), "READ");
    if (!f || f->IsZombie()) { printf("cannot open %s\n", inFile.c_str()); return; }

    TGraphErrors *g = dynamic_cast<TGraphErrors*>(f->Get("gpurity_leak"));
    if (!g) { printf("missing gpurity_leak\n"); return; }

    // Match nominal fit range in CalculatePhotonYield.C:641-642
    const double xMin = 10.0;
    const double xMax = 36.0;

    // Fit 1: Padé [1/1] (analysis nominal).  Initial parameters must match
    // CalculatePhotonYield.C (constants pade_purity_leak_p{0,1,2}) — the
    // (0.5, 0.5, 0.5) defaults converge to a degenerate horizontal-line
    // local minimum that gives chi^2/ndf ~ 150/8 instead of the true fit.
    TF1 *f_pade = new TF1("f_pade",
                          "([0] + [1]*x) / (1 + [2]*[2]*x)",
                          xMin, xMax);
    f_pade->SetParameters(0.085956, 0.052318, -0.164694);  // pade_purity_leak_p{0,1,2}
    g->Fit(f_pade, "REMN", "", xMin, xMax);
    g->Fit(f_pade, "REMN", "", xMin, xMax);
    const double pade_chi2 = f_pade->GetChisquare();
    const double pade_ndf  = f_pade->GetNDF();

    // Fit 2: Error function (alternative, used as fit-form systematic).
    TF1 *f_erf = new TF1("f_erf",
                         "[0]*TMath::Erf((x - [1])/[2])",
                         xMin, xMax);
    f_erf->SetParameters(1.0, 8.0, 8.0);
    g->Fit(f_erf, "REMN", "", xMin, xMax);
    g->Fit(f_erf, "REMN", "", xMin, xMax);
    const double erf_chi2 = f_erf->GetChisquare();
    const double erf_ndf  = f_erf->GetNDF();

    TCanvas *c = new TCanvas("c_purity_compare", "Purity fit comparison", 600, 600);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.13);

    TH1F *frame = new TH1F("frame_purity_compare", "", 1, 8, 36);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("#it{E}_{T}^{#gamma} [GeV]");
    frame->GetYaxis()->SetTitle("Purity");
    frame->GetYaxis()->SetTitleOffset(1.20);
    frame->GetXaxis()->SetTitleOffset(1.10);
    frame->GetYaxis()->SetRangeUser(0.0, 1.20);
    frame->Draw("axis");

    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlack);
    g->SetMarkerSize(1.0);
    g->SetLineColor(kBlack);
    g->Draw("same p");

    f_pade->SetLineColor(kRed + 1);
    f_pade->SetLineWidth(2);
    f_pade->SetLineStyle(1);
    f_pade->Draw("same");

    f_erf->SetLineColor(kAzure + 2);
    f_erf->SetLineWidth(2);
    f_erf->SetLineStyle(2);
    f_erf->Draw("same");

    const float xpos = 0.20, ypos = 0.85;
    const float dy = 0.050, fontsize = 0.038;
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),   fontsize, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(), fontsize, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),   fontsize, 0);

    TLegend *leg = new TLegend(0.42, 0.16, 0.93, 0.34);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.036);
    leg->AddEntry(g, "Data purity (leakage-corr.)", "lep");
    leg->AddEntry(f_pade,
                  Form("Pad#acute{e}[1/1] nom., #chi^{2}=%.1f/%d",
                       pade_chi2, (int)pade_ndf), "l");
    leg->AddEntry(f_erf,
                  Form("Erf alt., #chi^{2}=%.1f/%d",
                       erf_chi2, (int)erf_ndf), "l");
    leg->Draw();

    c->SaveAs(outFile.c_str());
    printf("Pade chi2/ndf = %.3f/%d\n", pade_chi2, (int)pade_ndf);
    printf("Erf  chi2/ndf = %.3f/%d\n", erf_chi2, (int)erf_ndf);
    printf("saved %s\n", outFile.c_str());
}
