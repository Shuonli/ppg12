#include "plotcommon.h"
#include "../efficiencytool/CrossSectionWeights.h"
using namespace PPG12;

void plot_SB()
{
    init_plot();

    // Updated 2026-04 to current shower-shape filename convention.
    TFile *fin_sig = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_signal_combined_showershape.root", "READ");
    TFile *fin_bg  = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_jet_inclusive_combined_showershape.root", "READ");

    // Convert the bg histogram from jet50-normalised units (the
    // SampleConfig::weight reference for jet samples) into the same
    // photon20-normalised units as the signal histogram, so S/B is a
    // proper physical ratio.
    float jet_scale = jet50cross / photon20cross;

    TH2D* h_sig = (TH2D *)fin_sig->Get("h_ET_isoET_eta0");
    TH2D* h_bg = (TH2D *)fin_bg->Get("h_ET_isoET_eta0");

    int rebinx = 16;
    h_sig->RebinX(rebinx);
    h_bg->RebinX(rebinx);

    // Project the iso-ET axis up to the parametric reco-iso ceiling
    // applied in the analysis, reco_iso_max(ET) = 0.502 + 0.0433 ET
    // (approximated bin-by-bin). Without this projection cut the
    // S/B ratio is dominated by the loose tail of the jet sample
    // and the figure looks empty on a [0, 1] axis.
    auto project_iso = [](TH2D *h2, const char *name) {
        TH1D *h1 = (TH1D *)h2->ProjectionX(name, 0, 0);  // template
        h1->Reset();
        for (int ix = 1; ix <= h2->GetNbinsX(); ++ix) {
            double et = h2->GetXaxis()->GetBinCenter(ix);
            double iso_max = 0.502095 + 0.0433036 * et;
            int iy_lo = h2->GetYaxis()->FindBin(-1.0);  // include negative
            int iy_hi = h2->GetYaxis()->FindBin(iso_max);
            double sum = 0, sumw2 = 0;
            for (int iy = iy_lo; iy <= iy_hi; ++iy) {
                sum   += h2->GetBinContent(ix, iy);
                sumw2 += h2->GetBinError(ix, iy) * h2->GetBinError(ix, iy);
            }
            h1->SetBinContent(ix, sum);
            h1->SetBinError(ix, std::sqrt(sumw2));
        }
        return h1;
    };
    TH1D *h_sig_proj = project_iso(h_sig, "h_sig_proj");
    TH1D *h_bg_proj  = project_iso(h_bg,  "h_bg_proj");

    h_sig_proj->Sumw2();
    h_bg_proj->Sumw2();

    //scale background
    h_bg_proj->Scale(jet_scale);

    //calculate s/b
    TH1D* h_sb = (TH1D *)h_sig_proj->Clone("h_sb");

    h_sb->Divide(h_bg_proj);

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    frame_et_rec->SetYTitle("S/B");
    double y_max = 1.2 * h_sb->GetMaximum();
    frame_et_rec->GetYaxis()->SetRangeUser(0.0, y_max < 1.0 ? 1.0 : y_max);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 32);
    frame_et_rec->Draw("axis");

    h_sb->SetMarkerColor(kBlack);
    h_sb->SetMarkerStyle(20);
    h_sb->SetMarkerSize(1.5);
    h_sb->SetLineColor(kBlack);
    h_sb->Draw("P same");
    //h_bg_proj->Draw();

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5, 0.80, 1, strMC.c_str(), 0.04);

    c1->SaveAs("figures/SB.pdf");










}