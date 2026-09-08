// plot_trigger_bit30_zoom.C
//
// Zoomed-in version of the photon-4 GeV (bit 30) L1 trigger turn-on,
// using the same visual style as the middle panel of Fig. 6 in the PPG12
// analysis note: a TEfficiency overlay (Clopper-Pearson 68% CI errors)
// fitted with the saturating asymmetric Gumbel form
//   eps(ET) = p0 * exp(-exp(-(ET-mu)/beta))
// and the 1-sigma fit-covariance band drawn as a shaded region.
//
// The numerator/denominator histograms are pulled from the precomputed
// high-statistics file
//   efficiencytool/results/trigger_bit30_turnon.root
// (h_num_percluster / h_den_percluster, 6.6M events) so that the residual
// sub-percent inefficiency in the plateau is actually resolved instead of
// clamping to 1.0 like the live 10-file refit does.
//
// Output: figures/Photon_4_GeV_PlateauZoom.pdf

#include "plotcommon.h"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"

#include <array>
#include <vector>

void plot_trigger_bit30_zoom()
{
    init_plot();

    const std::string inFile  =
        "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/"
        "results/trigger_bit30_turnon.root";
    const std::string outFile = "figures/Photon_4_GeV_PlateauZoom.pdf";

    TFile *f = TFile::Open(inFile.c_str(), "READ");
    if (!f || f->IsZombie()) { printf("cannot open %s\n", inFile.c_str()); return; }

    TH1D *h_num = (TH1D*) f->Get("h_num_percluster");
    TH1D *h_den = (TH1D*) f->Get("h_den_percluster");
    if (!h_num || !h_den) { printf("missing num/den\n"); return; }
    h_num->SetDirectory(nullptr);
    h_den->SetDirectory(nullptr);

    TEfficiency *eff = new TEfficiency(*h_num, *h_den);
    eff->SetStatisticOption(TEfficiency::kFCP);
    eff->SetConfidenceLevel(0.683);
    eff->SetMarkerStyle(20);
    eff->SetMarkerColor(kRed + 1);
    eff->SetMarkerSize(0.85);
    eff->SetLineColor(kRed + 1);

    TCanvas *c = new TCanvas("c_trig_zoom", "Bit-30 plateau zoom", 600, 600);
    c->SetLeftMargin(0.17);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.15);

    // Frame: zoom y-range to expose plateau inefficiency, x-range
    // covering the turn-on plus the analysis range.
    TH2F *frame = new TH2F("frame_trig_zoom", "",
                            100, 7.0, 30.0, 100, 0.985, 1.005);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("#it{E}_{T}^{cluster} [GeV]");
    frame->GetYaxis()->SetTitle("#varepsilon_{L1}(Photon 4 GeV | MBD N&S)");
    frame->GetYaxis()->SetTitleOffset(1.40);
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetYaxis()->SetNdivisions(505);
    frame->Draw("axis");

    eff->Draw("same p");

    // Fit: same Gumbel form as the middle panel of Fig. 6.
    const double fit_lo = 7.0, fit_hi = 30.0;
    TF1 *f_fit = new TF1("fgumbel_bit30",
                         "[0]*TMath::Exp(-TMath::Exp(-(x-[1])/[2]))",
                         fit_lo, fit_hi);
    f_fit->SetParameters(1.0, -2.0, 3.33);
    f_fit->SetParLimits(0, 0.5, 1.0);
    f_fit->SetParLimits(1, -50.0, 8.0);
    f_fit->SetParLimits(2, 0.1, 50.0);
    f_fit->SetLineColor(kRed + 1);
    f_fit->SetLineWidth(2);

    // Build a TH1F mirror of the per-bin efficiencies (with averaged
    // Clopper-Pearson errors) so we can run a proper TGraphErrors-style
    // fit and recover the parameter covariance for the CI band.
    const int nb = h_num->GetNbinsX();
    TH1F *h_eff_for_fit = new TH1F("h_eff_for_fit_bit30", "",
                                   h_num->GetXaxis()->GetNbins(),
                                   h_num->GetXaxis()->GetXbins()->GetSize() > 0
                                       ? h_num->GetXaxis()->GetXbins()->GetArray()
                                       : nullptr);
    if (h_num->GetXaxis()->GetXbins()->GetSize() == 0) {
        // Uniform binning fallback — re-create with the original limits
        delete h_eff_for_fit;
        h_eff_for_fit = new TH1F("h_eff_for_fit_bit30", "", nb,
                                 h_num->GetXaxis()->GetXmin(),
                                 h_num->GetXaxis()->GetXmax());
    }
    h_eff_for_fit->SetDirectory(nullptr);
    for (int i = 1; i <= nb; ++i) {
        const double tot = h_den->GetBinContent(i);
        if (tot <= 0) continue;
        const double eps = eff->GetEfficiency(i);
        const double el  = eff->GetEfficiencyErrorLow(i);
        const double eh  = eff->GetEfficiencyErrorUp(i);
        const double err = 0.5 * (el + eh);
        h_eff_for_fit->SetBinContent(i, eps);
        h_eff_for_fit->SetBinError(i, err > 0 ? err : 1.0 / std::sqrt(tot));
    }
    TFitResultPtr fr = h_eff_for_fit->Fit(f_fit, "R Q N S", "", fit_lo, fit_hi);

    // 1-sigma confidence band from the fit covariance.
    const int n_ci = 400;
    std::vector<double> x_arr(n_ci), ci_arr(n_ci);
    for (int i = 0; i < n_ci; ++i)
        x_arr[i] = fit_lo + i * (fit_hi - fit_lo) / (n_ci - 1);
    fr->GetConfidenceIntervals(n_ci, 1, 1, x_arr.data(), ci_arr.data(),
                                0.683, false);
    TGraphErrors *gci = new TGraphErrors(n_ci);
    for (int i = 0; i < n_ci; ++i) {
        gci->SetPoint(i, x_arr[i], f_fit->Eval(x_arr[i]));
        gci->SetPointError(i, 0.0, ci_arr[i]);
    }
    gci->SetFillColorAlpha(kRed + 1, 0.30);
    gci->SetFillStyle(1001);
    gci->SetLineWidth(0);
    gci->Draw("3 same");
    f_fit->Draw("same");
    eff->Draw("same p");  // re-draw points on top of band

    const float xpos = 0.22, ypos = 0.32;
    const float dy = 0.050, fontsize = 0.038;
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),   fontsize, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(), fontsize, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),   fontsize, 0);

    TLegend *leg = new TLegend(0.42, 0.76, 0.94, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.038);
    leg->AddEntry(eff, "Photon 4 GeV (bit 30, nom.)", "lep");
    leg->AddEntry((TObject*)gci,
                  Form("Gumbel fit, plateau = %.4f", f_fit->GetParameter(0)),
                  "f");
    leg->Draw();

    c->SaveAs(outFile.c_str());
    printf("plateau = %.4f #pm %.4f  chi2/ndf=%.2f/%d\n",
           f_fit->GetParameter(0), f_fit->GetParError(0),
           f_fit->GetChisquare(), f_fit->GetNDF());
    printf("saved %s\n", outFile.c_str());
}
