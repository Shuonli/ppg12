// plot_resolution_vs_data.C
//
// Compare PPG12 cluster-response MC resolution against Brennan's data-driven
// resolution fit. Reads the MC TH2 (cluster_Et / truth_pT vs truth_pT) and
// projects in fine 1-GeV truth-pT bins to avoid bin-width broadening, fits
// Gaussian sigma per bin, then overlays:
//   - sigma_MC (this analysis, native = no artificial smearing)
//   - sigma_data (Brennan f_energy_global_minimum)
//   - sqrt(sigma_data^2 - sigma_MC^2)  [additional smearing needed]
//
// Brennan's file: /gpfs/mnt/gpfs02/sphenix/user/bseidlitz/emcPerf/macros/figMaker/output/function_compare_wide.root
//   f_energy_global_minimum     : data-driven best-fit (p0,p1,p2) = (0.18, 0.02, 0.015)
//   f_energy_single_particle_gun: GEANT single-photon gun (0.13, 0.05, 0.02)
//
// Usage: root -l -b -q plot_resolution_vs_data.C

#include "plotcommon.h"

#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <iostream>
#include <vector>

void plot_resolution_vs_data()
{
    init_plot();
    gStyle->SetOptStat(0);

    const char *outdir = "/sphenix/user/shuhangli/ppg12/plotting/figures/cluster_response_tightiso";
    gSystem->mkdir(outdir, true);

    // --- Open MC TH2 (lumi-weighted all-range) ---
    TFile *fmc = TFile::Open("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_signal_combined_bdt_clusterresp_all.root");
    TH2D *h2 = (TH2D *)fmc->Get("h_pT_truth_reco_tightiso_0");
    if (!h2) { std::cerr << "TH2 missing\n"; return; }

    // --- Open external resolution fits (function_compare_wide.root) ---
    TFile *fdata = TFile::Open("/gpfs/mnt/gpfs02/sphenix/user/bseidlitz/emcPerf/macros/figMaker/output/function_compare_wide.root");
    TF1 *fE_gmin     = (TF1 *)fdata->Get("f_energy_global_minimum");
    TF1 *fE_cE0p05   = (TF1 *)fdata->Get("f_energy_example_cE_0p05");
    TF1 *fE_cE0p08   = (TF1 *)fdata->Get("f_energy_example_cE_0p08");
    TF1 *fP_gmin     = (TF1 *)fdata->Get("f_position_global_minimum");
    TF1 *fP_cE0p05   = (TF1 *)fdata->Get("f_position_example_cE_0p05");
    TF1 *fP_cE0p08   = (TF1 *)fdata->Get("f_position_example_cE_0p08");

    // PPG12 PYTHIA prompt-photon MC parametric fit — must match config_bdt_nom.yaml
    // cluster_eres_mc_p{0,1,2}. Form: sigma(pT) = sqrt(p0^2/pT + p1^2/pT^2 + p2^2).
    // This is what the live pipeline uses for the MC side of sigma_extra(pT).
    const double mc_p0 = 0.185, mc_p1 = 0.000, mc_p2 = 0.040;
    TF1 *fE_mc = new TF1("fE_mc",
        "sqrt([0]*[0]/x + [1]*[1]/(x*x) + [2]*[2])", 5, 45);
    fE_mc->SetParameters(mc_p0, mc_p1, mc_p2);
    fE_mc->SetLineColor(kBlack);
    fE_mc->SetLineStyle(1);
    fE_mc->SetLineWidth(3);

    // --- Analysis truth-pT bins: 2 GeV from 8 to 28, then [28,32), [32,36), [36,45) ---
    const int NBINS = 13;
    const double pt_edges_ana[NBINS + 1] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45};
    double bin_centers[NBINS], sigma_mc[NBINS], sigma_mc_err[NBINS];
    int n_valid = 0;
    double centers_v[NBINS], sigmas_v[NBINS], errs_v[NBINS];
    double sigma_tight_v[NBINS], sigma_tight_err_v[NBINS];
    double sigma_dscb_v[NBINS], sigma_dscb_err_v[NBINS];
    double sigma_data_v[NBINS], sigma_extra_v[NBINS], sigma_extra_err_v[NBINS];

    for (int i = 0; i < NBINS; ++i) {
        double pt_lo = pt_edges_ana[i];
        double pt_hi = pt_edges_ana[i + 1];
        int ix_lo = h2->GetXaxis()->FindBin(pt_lo + 1e-6);
        int ix_hi = h2->GetXaxis()->FindBin(pt_hi - 1e-6);
        TH1D *hp = h2->ProjectionY(Form("p_%d", i), ix_lo, ix_hi);
        hp->SetDirectory(nullptr);
        if (hp->Integral() < 200) { delete hp; continue; }
        double mean = hp->GetMean();
        double rms  = hp->GetStdDev();
        if (rms <= 0) { delete hp; continue; }
        double lo = std::max(0.3, mean - 2.0 * rms);
        double hi = std::min(1.7, mean + 2.0 * rms);
        TF1 fg("fg", "gaus", lo, hi);
        fg.SetParameter(0, hp->GetMaximum());
        fg.SetParameter(1, mean);
        fg.SetParameter(2, std::max(0.03, rms));
        int s = hp->Fit(&fg, "RQN0");
        if (s != 0) { delete hp; continue; }

        // Iteratively tighten window to ±1.5*sigma_fit (refit twice) — captures core only
        double mu_t = fg.GetParameter(1);
        double sg_t = fg.GetParameter(2);
        TF1 fg_tight("fg_tight", "gaus", mu_t - 1.5*sg_t, mu_t + 1.5*sg_t);
        fg_tight.SetParameters(fg.GetParameter(0), mu_t, sg_t);
        hp->Fit(&fg_tight, "RQN0");
        mu_t = fg_tight.GetParameter(1); sg_t = fg_tight.GetParameter(2);
        TF1 fg_tight2("fg_tight2", "gaus", mu_t - 1.5*sg_t, mu_t + 1.5*sg_t);
        fg_tight2.SetParameters(fg_tight.GetParameter(0), mu_t, sg_t);
        hp->Fit(&fg_tight2, "RQN0");

        // DSCB fit on [0.3, 1.6]
        double sg_dscb = 0, sg_dscb_err = 0;
        {
            TF1 fdscb("fdscb",
                "[0]*( ((x-[1])/[2] > -[3] && (x-[1])/[2] < [5]) ? exp(-0.5*((x-[1])/[2])*((x-[1])/[2])) : "
                "( (x-[1])/[2] <= -[3] ? pow([4]/abs([3]),[4])*exp(-0.5*[3]*[3])*pow([4]/abs([3])-abs([3]) - (x-[1])/[2], -[4]) : "
                "pow([6]/abs([5]),[6])*exp(-0.5*[5]*[5])*pow([6]/abs([5])-abs([5]) + (x-[1])/[2], -[6]) ))",
                0.3, 1.6);
            fdscb.SetParameters(hp->GetMaximum(), mean, std::max(0.03, rms), 1.0, 5.0, 1.5, 5.0);
            fdscb.SetParLimits(2, 0.01, 0.4);
            fdscb.SetParLimits(3, 0.1, 10); fdscb.SetParLimits(4, 1.05, 100);
            fdscb.SetParLimits(5, 0.1, 10); fdscb.SetParLimits(6, 1.05, 100);
            int sd = hp->Fit(&fdscb, "RQN0");
            if (sd == 0) { sg_dscb = fdscb.GetParameter(2); sg_dscb_err = fdscb.GetParError(2); }
        }
        delete hp;

        double pt_c    = 0.5 * (pt_lo + pt_hi);
        double sg      = fg.GetParameter(2);
        double sg_err  = fg.GetParError(2);
        double sd      = fE_gmin->Eval(pt_c);

        centers_v[n_valid]         = pt_c;
        sigmas_v[n_valid]          = sg;
        errs_v[n_valid]            = sg_err;
        sigma_tight_v[n_valid]     = fg_tight2.GetParameter(2);
        sigma_tight_err_v[n_valid] = fg_tight2.GetParError(2);
        sigma_dscb_v[n_valid]      = sg_dscb;
        sigma_dscb_err_v[n_valid]  = sg_dscb_err;
        sigma_data_v[n_valid]      = sd;

        double sg_use = sigma_tight_v[n_valid];  // use tight-window σ for quad diff
        double diff2 = sd*sd - sg_use*sg_use;
        sigma_extra_v[n_valid]     = (diff2 > 0) ? sqrt(diff2) : 0.0;
        sigma_extra_err_v[n_valid] = (diff2 > 0) ? sg_use * sigma_tight_err_v[n_valid] / sqrt(diff2) : sigma_tight_err_v[n_valid];
        n_valid++;
    }

    // --- Build TGraphs ---
    auto *g_mc = new TGraphErrors(n_valid, centers_v, sigmas_v, nullptr, errs_v);
    g_mc->SetMarkerStyle(20); g_mc->SetMarkerColor(kBlack); g_mc->SetLineColor(kBlack);
    g_mc->SetMarkerSize(1.1);  g_mc->SetLineWidth(2);

    // ----- Plot 1: ENERGY resolution overlay (single panel) -----
    TCanvas c("c", "", 900, 700);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.16);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.04);

    auto style_fn = [](TF1 *f, int col, int ls, int lw) {
        if (!f) return;
        f->SetLineColor(col); f->SetLineStyle(ls); f->SetLineWidth(lw);
        f->SetRange(5, 45);
    };
    style_fn(fE_gmin,   kRed+1,    2, 3);
    style_fn(fE_cE0p05, kGreen+2,  9, 2);
    style_fn(fE_cE0p08, kViolet+2, 6, 2);

    TH1F *frE = new TH1F("frE", "", 1, 5, 45);
    frE->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
    frE->SetYTitle("Energy resolution  #sigma_{E}/E");
    frE->GetYaxis()->SetRangeUser(0.02, 0.14);
    frE->GetXaxis()->SetRangeUser(5, 45);
    frE->GetYaxis()->SetTitleOffset(1.30);
    frE->GetXaxis()->SetTitleOffset(1.15);
    frE->Draw();
    fE_cE0p05->Draw("SAME"); fE_cE0p08->Draw("SAME"); fE_gmin->Draw("SAME");
    fE_mc    ->Draw("SAME");
    g_mc     ->Draw("P SAME");

    TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.038);
    l.DrawLatex(0.18, 0.91, strleg1.c_str());
    l.DrawLatex(0.18, 0.85, strleg2.c_str());
    l.DrawLatex(0.18, 0.79, strleg3.c_str());
    l.DrawLatex(0.18, 0.73, strleg4.c_str());
    l.DrawLatex(0.18, 0.67, "|#it{z}_{reco}| < 60 cm");

    TLegend legE(0.50, 0.66, 0.96, 0.93);
    legE.SetFillStyle(0); legE.SetBorderSize(0); legE.SetTextFont(42); legE.SetTextSize(0.030);
    legE.AddEntry(g_mc,      "PYTHIA prompt #gamma #sigma_{Gauss}", "lp");
    legE.AddEntry(fE_mc,     "MC param fit (0.185, 0, 0.040)", "l");
    legE.AddEntry(fE_gmin,   "f_energy_global_minimum", "l");
    legE.AddEntry(fE_cE0p05, "f_energy_example_cE_0p05", "l");
    legE.AddEntry(fE_cE0p08, "f_energy_example_cE_0p08", "l");
    legE.Draw();

    c.SaveAs(Form("%s/resolution_vs_data.pdf", outdir));

    // ----- Plot 2: quadrature differences (one line per data fit) -----
    //   sqrt(sigma_data^2 - sigma_MC^2) at each analysis bin center.
    //   Where this is imaginary (sigma_MC > sigma_data), the point is set to 0.
    struct DataFit { TF1 *fn; int col; int mks; const char *name; };
    DataFit fits[] = {
        {fE_gmin,   kRed+1,    20, "f_energy_global_minimum"},
        {fE_cE0p05, kGreen+2,  21, "f_energy_example_cE_0p05"},
        {fE_cE0p08, kViolet+2, 22, "f_energy_example_cE_0p08"},
    };
    const int NFITS = sizeof(fits)/sizeof(fits[0]);

    double cx[NBINS];
    for (int i = 0; i < n_valid; i++) cx[i] = centers_v[i];

    std::vector<TGraphErrors *> g_qd_v;
    for (int f = 0; f < NFITS; ++f) {
        double cy[NBINS], cye[NBINS];
        for (int i = 0; i < n_valid; i++) {
            double m  = sigmas_v[i];
            double me = errs_v[i];
            double d  = fits[f].fn->Eval(cx[i]);
            double diff2 = d*d - m*m;
            cy [i] = (diff2 > 0) ? sqrt(diff2) : 0.0;
            cye[i] = (cy[i] > 0) ? m * me / cy[i] : me;
        }
        auto *g = new TGraphErrors(n_valid, cx, cy, nullptr, cye);
        g->SetMarkerStyle(fits[f].mks); g->SetMarkerColor(fits[f].col); g->SetLineColor(fits[f].col);
        g->SetMarkerSize(1.1); g->SetLineWidth(2);
        g_qd_v.push_back(g);
    }

    // Smooth analytical sigma_extra(pT) = sqrt(sigma_data^2 - sigma_MC,param^2) using
    // the MC param fit (matches what the pipeline applies per cluster, modulo binning).
    std::vector<TGraph *> g_qd_curve_v;
    const int NS = 200;
    double pt_lo_s = 5.0, pt_hi_s = 45.0;
    for (int f = 0; f < NFITS; ++f) {
        double xs[NS], ys[NS];
        for (int i = 0; i < NS; ++i) {
            double pt = pt_lo_s + (pt_hi_s - pt_lo_s) * i / (NS - 1);
            double d  = fits[f].fn->Eval(pt);
            double d2 = d * d;
            double m2 = (mc_p0 * mc_p0) / pt + (mc_p1 * mc_p1) / (pt * pt) + mc_p2 * mc_p2;
            xs[i] = pt;
            ys[i] = (d2 > m2) ? std::sqrt(d2 - m2) : 0.0;
        }
        auto *g = new TGraph(NS, xs, ys);
        g->SetLineColor(fits[f].col); g->SetLineStyle(1); g->SetLineWidth(2);
        g_qd_curve_v.push_back(g);
    }

    TCanvas c2("c2", "", 900, 700);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.16);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.04);

    TH1F *frQ = new TH1F("frQ", "", 1, 5, 45);
    frQ->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
    frQ->SetYTitle("#sqrt{#sigma_{data}^{2} #minus #sigma_{MC}^{2}}");
    frQ->GetYaxis()->SetRangeUser(0.0, 0.14);
    frQ->GetXaxis()->SetRangeUser(5, 45);
    frQ->GetYaxis()->SetTitleOffset(1.30);
    frQ->GetXaxis()->SetTitleOffset(1.15);
    frQ->Draw();

    for (auto *g : g_qd_curve_v) g->Draw("L SAME");
    for (auto *g : g_qd_v)       g->Draw("P SAME");

    // dashed reference at 0
    TLine *L0 = new TLine(5, 0, 45, 0);
    L0->SetLineStyle(7); L0->SetLineColor(kGray+2); L0->Draw("SAME");

    TLatex l2; l2.SetNDC(); l2.SetTextFont(42); l2.SetTextSize(0.038);
    l2.DrawLatex(0.18, 0.91, strleg1.c_str());
    l2.DrawLatex(0.18, 0.85, strleg2.c_str());
    l2.DrawLatex(0.18, 0.79, strleg3.c_str());
    l2.DrawLatex(0.18, 0.73, strleg4.c_str());
    l2.DrawLatex(0.18, 0.67, "|#it{z}_{reco}| < 60 cm");

    TLegend legQ(0.35, 0.55, 0.95, 0.92);
    legQ.SetFillStyle(0); legQ.SetBorderSize(0); legQ.SetTextFont(42); legQ.SetTextSize(0.026);
    legQ.SetHeader("Lines: smooth using MC param fit (0.185, 0, 0.040)", "C");
    for (int f = 0; f < NFITS; ++f) {
        legQ.AddEntry(g_qd_v[f],       Form("%s (per-bin)", fits[f].name), "p");
        legQ.AddEntry(g_qd_curve_v[f], Form("%s (param fit)",  fits[f].name), "l");
    }
    legQ.Draw();

    c2.SaveAs(Form("%s/resolution_quad_diff.pdf", outdir));

    // --- Print numbers ---
    FILE *fout = fopen(Form("%s/resolution_vs_data_numbers.txt", outdir), "w");
    fprintf(fout, "# pt_center  sg_2RMS  sg_2RMS_err  sg_tight  sg_tight_err  sg_dscb  sg_dscb_err  sigma_data  sigma_extra_needed (using sg_tight)  sigma_extra_err\n");
    std::cout << "\n# pt_c   sg_2RMS  sg_tight  sg_dscb  sg_data  sg_extra(tight)\n";
    for (int i = 0; i < n_valid; i++) {
        fprintf(fout, "%5.1f  %.4f %.4f  %.4f %.4f  %.4f %.4f  %.4f  %.4f %.4f\n",
                centers_v[i], sigmas_v[i], errs_v[i], sigma_tight_v[i], sigma_tight_err_v[i],
                sigma_dscb_v[i], sigma_dscb_err_v[i],
                sigma_data_v[i], sigma_extra_v[i], sigma_extra_err_v[i]);
        std::cout << Form("%5.1f  %.4f  %.4f  %.4f  %.4f  %.4f",
                          centers_v[i], sigmas_v[i], sigma_tight_v[i], sigma_dscb_v[i],
                          sigma_data_v[i], sigma_extra_v[i]) << "\n";
    }
    fclose(fout);
    std::cout << "[plot_resolution_vs_data] PDF + numbers written to " << outdir << "\n";
}
