// plot_cluster_response.C
//
// Photon cluster ET response and resolution for the tight + isolated
// prompt-photon arm of the analysis, using the double-interaction blended
// MC pipeline at 0 mrad, 1.5 mrad, and the lumi-weighted combination.
//
// Input: per-period hadded TH2 of (truth_pT, cluster_ET/truth_pT)
//        h_pT_truth_reco_tightiso_0 filled by RecoEffCalculator_TTreeReader.C
//        inside the tight && iso arm of the prompt-photon truth-match block.
//
// Outputs (plotting/figures/cluster_response_tightiso/):
//   resp_grid_{0rad,1p5mrad,all}.pdf          — 4x3 grid of 1D projections
//                                                per truth-ET bin with Gauss
//                                                + DSCB fit overlays.
//   resp_mean_vs_pt.pdf                       — Gauss mean vs truth ET, 3 overlays.
//   resp_sigma_vs_pt.pdf                      — Gauss sigma vs truth ET.
//   resp_rms_vs_pt.pdf                        — Histogram RMS vs truth ET (cross-check).
//   resp_summary_4panel.pdf                   — 2x2 summary: mu, sigma, RMS, mean.
//
// Usage:
//   root -l -b -q 'plot_cluster_response.C()'
//
// Companion hadd: hadd_clusterresp.sh produces the input files.

#include "plotcommon.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TSystem.h>
#include <TStyle.h>
#include <iostream>
#include <vector>

namespace {

// =========================================================================
// Double-sided Crystal-Ball PDF — matches the one in compare_energy_response.C
// =========================================================================
static Double_t dscb_fn(Double_t *x, Double_t *par)
{
    Double_t N  = par[0];
    Double_t mu = par[1];
    Double_t sg = par[2];
    Double_t aL = par[3], nL = par[4];
    Double_t aH = par[5], nH = par[6];
    if (sg <= 0) return 0;
    Double_t t = (x[0] - mu) / sg;
    if (t > -aL && t < aH) return N * TMath::Exp(-0.5 * t * t);
    if (t <= -aL) {
        Double_t A = TMath::Power(nL / TMath::Abs(aL), nL) * TMath::Exp(-0.5 * aL * aL);
        Double_t B = nL / TMath::Abs(aL) - TMath::Abs(aL);
        return N * A * TMath::Power(B - t, -nL);
    }
    Double_t A = TMath::Power(nH / TMath::Abs(aH), nH) * TMath::Exp(-0.5 * aH * aH);
    Double_t B = nH / TMath::Abs(aH) - TMath::Abs(aH);
    return N * A * TMath::Power(B + t, -nH);
}

struct PeriodInfo {
    const char *tag;
    const char *label;
    int  color;
    int  marker;
};

static const PeriodInfo kPeriods[] = {
    {"0rad",    "0 mrad (22.4% DI)",        kRed+1,    21},
    {"1p5mrad", "1.5 mrad (7.9% DI)",       kBlue+1,   22},
    {"all",     "Lumi-weighted combined",   kBlack,    20}
};
static const int kNper = sizeof(kPeriods) / sizeof(kPeriods[0]);

// Truth pT bins from config_bdt_nom_{0rad,1p5mrad}.yaml (`pT_bins_truth`):
//   {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45}  → 13 bins.
// Extends below the reco range on the low side (8-10 covers truth pT just
// below the reco ET>8 GeV threshold) and above on the high side for the
// unfolding overflow.
// The truth-pT axis of h_pT_truth_reco_tightiso is 400 bins over [5, 40] GeV;
// the [36, 45) bin is truncated at 40 GeV but still captures the high-ET tail.
static const int   NptBinsTruth = 13;
static const float ptRangesTruth[NptBinsTruth + 1] =
    {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45};

static void draw_sphenix_labels(double x = 0.22, double y0 = 0.88,
                                double dy = 0.052, double sz = 0.040,
                                const char *extra = nullptr)
{
    TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(sz);
    l.DrawLatex(x, y0,        strleg1.c_str());
    l.DrawLatex(x, y0 -   dy, strleg2_1.c_str());
    l.DrawLatex(x, y0 - 2*dy, strleg3.c_str());
    l.DrawLatex(x, y0 - 3*dy, strleg4.c_str());
    if (extra) l.DrawLatex(x, y0 - 4*dy, extra);
}

// Find the truth-pT axis bin range covering [pt_lo, pt_hi)
static std::pair<int,int> pt_bin_range(TH2D *h2, double pt_lo, double pt_hi)
{
    int ix_lo = h2->GetXaxis()->FindBin(pt_lo + 1e-6);
    int ix_hi = h2->GetXaxis()->FindBin(pt_hi - 1e-6);
    // Truth-pT axis only extends to 40 GeV; clip the high edge of the last bin.
    int nbins = h2->GetXaxis()->GetNbins();
    if (ix_hi > nbins) ix_hi = nbins;
    return {ix_lo, ix_hi};
}

// Project [pt_lo, pt_hi) on the y-axis (cluster_Et / truth_pT)
static TH1D *project_response(TH2D *h2, int ibin, const char *suffix)
{
    auto rng = pt_bin_range(h2, ptRangesTruth[ibin], ptRangesTruth[ibin+1]);
    TH1D *hp = h2->ProjectionY(Form("%s_pt%g_%g", suffix, ptRangesTruth[ibin], ptRangesTruth[ibin+1]),
                               rng.first, rng.second);
    hp->SetDirectory(nullptr);
    return hp;
}

// Gaussian fit in a window around the histogram mean.
// Returns success flag; out_mu / out_sg / out_mu_err / out_sg_err filled on success.
static bool fit_gauss(TH1D *hp, double &out_mu, double &out_mu_err,
                      double &out_sg, double &out_sg_err)
{
    if (!hp || hp->Integral() < 100) return false;
    double mean = hp->GetMean();
    double rms  = hp->GetStdDev();
    if (rms <= 0) return false;
    double lo = std::max(0.3, mean - 2.0 * rms);
    double hi = std::min(1.7, mean + 2.0 * rms);
    TF1 fg("fg", "gaus", lo, hi);
    fg.SetParameter(0, hp->GetMaximum());
    fg.SetParameter(1, mean);
    fg.SetParameter(2, std::max(0.03, rms));
    int status = hp->Fit(&fg, "RQN0");
    if (status != 0) return false;
    out_mu     = fg.GetParameter(1);
    out_mu_err = fg.GetParError(1);
    out_sg     = fg.GetParameter(2);
    out_sg_err = fg.GetParError(2);
    return true;
}

// Double-sided Crystal-Ball fit — captures low-side tail from punch-through
static TF1 *fit_dscb(TH1D *hp, const char *name)
{
    if (!hp || hp->Integral() < 100) return nullptr;
    TF1 *fcb = new TF1(name, dscb_fn, 0.3, 1.6, 7);
    double mean = hp->GetMean();
    double rms  = std::max(0.03, (double)hp->GetStdDev());
    fcb->SetParameters(hp->GetMaximum(), mean, rms, 1.0, 5.0, 1.5, 5.0);
    fcb->SetParLimits(2, 0.01, 0.4);
    fcb->SetParLimits(3, 0.1, 10.0);
    fcb->SetParLimits(4, 1.05, 100.0);
    fcb->SetParLimits(5, 0.1, 10.0);
    fcb->SetParLimits(6, 1.05, 100.0);
    fcb->SetLineColor(kGreen+2);
    fcb->SetLineWidth(2);
    int s = hp->Fit(fcb, "RQN0");
    if (s != 0) {
        delete fcb;
        return nullptr;
    }
    return fcb;
}

} // namespace

void plot_cluster_response(const char *inputdir = "/sphenix/user/shuhangli/ppg12/efficiencytool/results",
                            const char *outdir   = "/sphenix/user/shuhangli/ppg12/plotting/figures/cluster_response_tightiso")
{
    init_plot();
    gStyle->SetOptStat(0);
    gSystem->mkdir(outdir, true);

    // --- Load per-period TH2s ---------------------------------------------
    std::vector<TFile *> fin(kNper, nullptr);
    std::vector<TH2D *>  h2(kNper, nullptr);
    for (int ip = 0; ip < kNper; ++ip) {
        TString path = Form("%s/MC_efficiency_signal_combined_bdt_clusterresp_%s.root",
                            inputdir, kPeriods[ip].tag);
        fin[ip] = TFile::Open(path);
        if (!fin[ip] || fin[ip]->IsZombie()) {
            std::cerr << "[plot_cluster_response] FATAL: cannot open " << path << std::endl;
            return;
        }
        h2[ip] = (TH2D *) fin[ip]->Get("h_pT_truth_reco_tightiso_0");
        if (!h2[ip]) {
            std::cerr << "[plot_cluster_response] FATAL: missing h_pT_truth_reco_tightiso_0 in "
                      << path << std::endl;
            return;
        }
        std::cout << "[plot_cluster_response] " << kPeriods[ip].tag
                  << "  integral=" << h2[ip]->Integral() << std::endl;
    }

    // --- Per-period: bin-by-bin Gaussian fits ----------------------------
    // Allocate summary histos: mean and sigma vs truth pT, one per period.
    double pt_edges[NptBinsTruth + 1];
    for (int i = 0; i <= NptBinsTruth; ++i) pt_edges[i] = ptRangesTruth[i];

    std::vector<TH1D *> h_mean(kNper);
    std::vector<TH1D *> h_sigma(kNper);
    std::vector<TH1D *> h_rms(kNper);
    std::vector<TH1D *> h_histmean(kNper);
    for (int ip = 0; ip < kNper; ++ip) {
        h_mean[ip]     = new TH1D(Form("h_mean_%s",  kPeriods[ip].tag),  ";#it{E}_{T}^{#gamma, truth} [GeV];#mu (Gauss) #LT #it{E}_{T}^{cluster}/#it{E}_{T}^{#gamma, truth}#GT", NptBinsTruth, pt_edges);
        h_sigma[ip]    = new TH1D(Form("h_sigma_%s", kPeriods[ip].tag),  ";#it{E}_{T}^{#gamma, truth} [GeV];#sigma (Gauss) of #it{E}_{T}^{cluster}/#it{E}_{T}^{#gamma, truth}", NptBinsTruth, pt_edges);
        h_rms[ip]      = new TH1D(Form("h_rms_%s",   kPeriods[ip].tag),  ";#it{E}_{T}^{#gamma, truth} [GeV];RMS of #it{E}_{T}^{cluster}/#it{E}_{T}^{#gamma, truth}", NptBinsTruth, pt_edges);
        h_histmean[ip] = new TH1D(Form("h_histmean_%s", kPeriods[ip].tag), ";#it{E}_{T}^{#gamma, truth} [GeV];Hist mean #LT #it{E}_{T}^{cluster}/#it{E}_{T}^{#gamma, truth}#GT", NptBinsTruth, pt_edges);
        h_mean[ip]->SetDirectory(nullptr);
        h_sigma[ip]->SetDirectory(nullptr);
        h_rms[ip]->SetDirectory(nullptr);
        h_histmean[ip]->SetDirectory(nullptr);
    }

    // Store per-bin projections (and DSCB fits) so we can draw the QC grid.
    std::vector<std::vector<TH1D *>> proj(kNper, std::vector<TH1D *>(NptBinsTruth, nullptr));
    std::vector<std::vector<TF1  *>> dscb(kNper, std::vector<TF1  *>(NptBinsTruth, nullptr));

    for (int ip = 0; ip < kNper; ++ip) {
        for (int ib = 0; ib < NptBinsTruth; ++ib) {
            TH1D *hp = project_response(h2[ip], ib, Form("p_%s", kPeriods[ip].tag));
            proj[ip][ib] = hp;
            double mu = 0, mu_err = 0, sg = 0, sg_err = 0;
            bool ok = fit_gauss(hp, mu, mu_err, sg, sg_err);
            if (ok) {
                h_mean[ip]->SetBinContent(ib + 1, mu);
                h_mean[ip]->SetBinError(ib + 1, mu_err);
                h_sigma[ip]->SetBinContent(ib + 1, sg);
                h_sigma[ip]->SetBinError(ib + 1, sg_err);
            } else {
                std::cerr << "[plot_cluster_response] Gauss fit failed for "
                          << kPeriods[ip].tag << " bin " << ib << " (pT " << ptRangesTruth[ib] << "-" << ptRangesTruth[ib+1] << ")" << std::endl;
            }
            // Histogram mean / RMS — restrict to [0.3, 1.7] to remove far-tail outliers
            int b_lo = hp->FindBin(0.3 + 1e-6);
            int b_hi = hp->FindBin(1.7 - 1e-6);
            hp->GetXaxis()->SetRange(b_lo, b_hi);
            h_histmean[ip]->SetBinContent(ib + 1, hp->GetMean());
            h_histmean[ip]->SetBinError  (ib + 1, hp->GetMeanError());
            h_rms[ip]     ->SetBinContent(ib + 1, hp->GetStdDev());
            h_rms[ip]     ->SetBinError  (ib + 1, hp->GetStdDevError());
            hp->GetXaxis()->SetRange(0, 0);  // restore

            dscb[ip][ib] = fit_dscb(hp, Form("dscb_%s_%d", kPeriods[ip].tag, ib));
        }
    }

    // -- Style helper
    auto style_summary = [](TH1D *h, const PeriodInfo &p, double marker_size = 1.2) {
        h->SetMarkerColor(p.color);
        h->SetLineColor(p.color);
        h->SetMarkerStyle(p.marker);
        h->SetMarkerSize(marker_size);
        h->SetLineWidth(2);
    };

    // --- Summary overlays --------------------------------------------------
    auto draw_summary = [&](std::vector<TH1D *> &hs, double ylo, double yhi,
                            const char *ytitle, const char *outname,
                            bool add_unity_line)
    {
        TCanvas c("c", "", 800, 650);
        gPad->SetLeftMargin(0.17);
        gPad->SetBottomMargin(0.17);
        TH1F *fr = new TH1F("fr_summary", "", 1, 8, 45);
        fr->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
        fr->SetYTitle(ytitle);
        fr->GetYaxis()->SetRangeUser(ylo, yhi);
        fr->GetXaxis()->SetRangeUser(8, 45);
        fr->GetYaxis()->SetTitleOffset(1.55);
        fr->Draw();
        for (int ip = 0; ip < kNper; ++ip) {
            style_summary(hs[ip], kPeriods[ip]);
            hs[ip]->Draw("E1 SAME");
        }
        if (add_unity_line) {
            TLine *L = new TLine(8, 1.0, 45, 1.0);
            L->SetLineStyle(7); L->SetLineColor(kGray+2);
            L->Draw("SAME");
        }
        draw_sphenix_labels(0.22, 0.88, 0.052, 0.040, "Tight + iso, prompt photon");
        TLegend leg(0.55, 0.20, 0.92, 0.34);
        leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextFont(42); leg.SetTextSize(0.035);
        for (int ip = 0; ip < kNper; ++ip) leg.AddEntry(hs[ip], kPeriods[ip].label, "lp");
        leg.Draw();
        c.SaveAs(Form("%s/%s", outdir, outname));
        delete fr;
    };

    draw_summary(h_mean,     0.88, 1.16,
                 "Gauss fit #mu of #it{E}_{T}^{cluster}/#it{E}_{T}^{#gamma, truth}",
                 "resp_mean_vs_pt.pdf", true);
    draw_summary(h_sigma,    0.02, 0.12,
                 "Gauss fit #sigma of #it{E}_{T}^{cluster}/#it{E}_{T}^{#gamma, truth}",
                 "resp_sigma_vs_pt.pdf", false);
    draw_summary(h_rms,      0.02, 0.16,
                 "Histogram RMS of #it{E}_{T}^{cluster}/#it{E}_{T}^{#gamma, truth}",
                 "resp_rms_vs_pt.pdf", false);
    draw_summary(h_histmean, 0.88, 1.16,
                 "Histogram mean of #it{E}_{T}^{cluster}/#it{E}_{T}^{#gamma, truth}",
                 "resp_histmean_vs_pt.pdf", true);

    // --- 4-panel summary ---------------------------------------------------
    {
        TCanvas c4("c4", "", 1200, 950);
        c4.Divide(2, 2, 0.001, 0.001);
        struct Panel { std::vector<TH1D *> &hs; double ylo, yhi; const char *ytit; bool unity; };
        std::vector<Panel> panels = {
            {h_mean,     0.88, 1.16, "Gauss #mu",       true},
            {h_sigma,    0.02, 0.12, "Gauss #sigma",    false},
            {h_histmean, 0.88, 1.16, "Hist mean",       true},
            {h_rms,      0.02, 0.16, "Hist RMS",        false}
        };
        for (size_t k = 0; k < panels.size(); ++k) {
            c4.cd(k + 1);
            gPad->SetLeftMargin(0.17);
            gPad->SetBottomMargin(0.17);
            TH1F *fr = new TH1F(Form("fr4_%zu", k), "", 1, 8, 45);
            fr->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
            fr->SetYTitle(Form("%s of #it{E}_{T}^{cluster}/#it{E}_{T}^{#gamma, truth}", panels[k].ytit));
            fr->GetYaxis()->SetRangeUser(panels[k].ylo, panels[k].yhi);
            fr->GetXaxis()->SetRangeUser(8, 45);
            fr->Draw();
            for (int ip = 0; ip < kNper; ++ip) {
                style_summary(panels[k].hs[ip], kPeriods[ip]);
                panels[k].hs[ip]->Draw("E1 SAME");
            }
            if (panels[k].unity) {
                TLine *L = new TLine(10, 1.0, 36, 1.0);
                L->SetLineStyle(7); L->SetLineColor(kGray+2);
                L->Draw("SAME");
            }
            if (k == 0) {
                draw_sphenix_labels(0.22, 0.88, 0.060, 0.045, "Tight + iso, prompt #gamma");
            }
            if (k == 1) {
                TLegend *leg = new TLegend(0.55, 0.60, 0.92, 0.85);
                leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextFont(42); leg->SetTextSize(0.04);
                for (int ip = 0; ip < kNper; ++ip) leg->AddEntry(panels[k].hs[ip], kPeriods[ip].label, "lp");
                leg->Draw();
            }
        }
        c4.SaveAs(Form("%s/resp_summary_4panel.pdf", outdir));
    }

    // --- Per-period QC grids: 4x3 (11 panels + 1 legend) ------------------
    auto draw_grid = [&](int ip) {
        // 13 truth bins + 1 legend → 4 columns × 4 rows (15 slots, 14 used).
        TCanvas cg(Form("cg_%s", kPeriods[ip].tag), "", 1600, 1400);
        cg.Divide(4, 4, 0.0008, 0.0008);
        for (int ib = 0; ib < NptBinsTruth; ++ib) {
            cg.cd(ib + 1);
            gPad->SetLeftMargin(0.18);
            gPad->SetBottomMargin(0.20);
            gPad->SetRightMargin(0.04);
            gPad->SetTopMargin(0.06);
            TH1D *hp = proj[ip][ib];
            if (!hp || hp->Integral() <= 0) continue;
            // Frame
            TH1F *fr = new TH1F(Form("frgrid_%s_%d", kPeriods[ip].tag, ib), "", 1, 0.2, 1.8);
            fr->SetXTitle("#it{E}_{T}^{cluster} / #it{E}_{T}^{#gamma, truth}");
            fr->SetYTitle("Counts (per bin)");
            double ymax = 1.3 * hp->GetMaximum();
            fr->GetYaxis()->SetRangeUser(0, ymax);
            fr->GetYaxis()->SetTitleOffset(1.6);
            fr->GetXaxis()->SetTitleOffset(1.2);
            fr->GetXaxis()->SetLabelSize(0.05);
            fr->GetYaxis()->SetLabelSize(0.05);
            fr->GetXaxis()->SetTitleSize(0.055);
            fr->GetYaxis()->SetTitleSize(0.055);
            fr->Draw();
            hp->SetLineColor(kBlack);
            hp->SetMarkerColor(kBlack);
            hp->SetMarkerStyle(20);
            hp->SetMarkerSize(0.7);
            hp->SetLineWidth(1);
            hp->Draw("E1 SAME");

            // Gauss + DSCB overlays
            double mu = h_mean[ip]->GetBinContent(ib + 1);
            double sg = h_sigma[ip]->GetBinContent(ib + 1);
            if (sg > 0) {
                TF1 *fg = new TF1(Form("fg_grid_%s_%d", kPeriods[ip].tag, ib), "gaus", 0.3, 1.6);
                fg->SetParameters(hp->GetMaximum(), mu, sg);
                fg->SetLineColor(kBlue+1);
                fg->SetLineStyle(2);
                fg->SetLineWidth(2);
                fg->Draw("SAME");
            }
            if (dscb[ip][ib]) {
                dscb[ip][ib]->SetRange(0.3, 1.6);
                dscb[ip][ib]->Draw("SAME");
            }

            TLine *L = new TLine(1.0, 0, 1.0, ymax);
            L->SetLineStyle(7); L->SetLineColor(kGray+2);
            L->Draw("SAME");

            // Per-panel label: ET range + Gauss mu / sigma
            TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.055);
            t.DrawLatex(0.22, 0.88, Form("%g < #it{E}_{T}^{truth} < %g GeV", ptRangesTruth[ib], ptRangesTruth[ib+1]));
            t.SetTextSize(0.045);
            t.DrawLatex(0.22, 0.80, Form("#mu = %.3f, #sigma = %.3f", mu, sg));
        }
        // 14th panel: legend (15th and 16th panels stay empty)
        cg.cd(NptBinsTruth + 1);
        TLatex lab; lab.SetNDC(); lab.SetTextFont(42);
        lab.SetTextSize(0.075);
        lab.DrawLatex(0.10, 0.85, strleg1.c_str());
        lab.SetTextSize(0.060);
        lab.DrawLatex(0.10, 0.75, strleg2_1.c_str());
        lab.DrawLatex(0.10, 0.67, strleg3.c_str());
        lab.DrawLatex(0.10, 0.59, strleg4.c_str());
        lab.SetTextSize(0.055);
        lab.DrawLatex(0.10, 0.49, Form("MC: %s", kPeriods[ip].label));
        lab.DrawLatex(0.10, 0.41, "Tight + iso, prompt photon matched");

        // Use a TLegend with TGraph stubs so the line swatches always render.
        TGraph *gG = new TGraph(); gG->SetLineColor(kBlue+1);  gG->SetLineStyle(2); gG->SetLineWidth(3);
        TGraph *gC = new TGraph(); gC->SetLineColor(kGreen+2); gC->SetLineStyle(1); gC->SetLineWidth(3);
        TLegend *legLines = new TLegend(0.08, 0.18, 0.92, 0.34);
        legLines->SetFillStyle(0); legLines->SetBorderSize(0);
        legLines->SetTextFont(42); legLines->SetTextSize(0.052);
        legLines->AddEntry(gG, "Gaussian fit ([#mu-2#sigma, #mu+2#sigma])", "l");
        legLines->AddEntry(gC, "Double-sided Crystal Ball fit", "l");
        legLines->Draw();

        cg.SaveAs(Form("%s/resp_grid_%s.pdf", outdir, kPeriods[ip].tag));
    };

    for (int ip = 0; ip < kNper; ++ip) draw_grid(ip);

    // --- Save summary numbers to a text file for re-derivation ------------
    FILE *fout = fopen(Form("%s/resp_summary_numbers.txt", outdir), "w");
    fprintf(fout, "# period  pT_lo  pT_hi  gauss_mu  gauss_mu_err  gauss_sigma  gauss_sigma_err  hist_mean  hist_RMS  N_entries\n");
    for (int ip = 0; ip < kNper; ++ip) {
        for (int ib = 0; ib < NptBinsTruth; ++ib) {
            TH1D *hp = proj[ip][ib];
            fprintf(fout, "%-10s %5.1f %5.1f  %.4f %.4f  %.4f %.4f  %.4f %.4f  %.0f\n",
                    kPeriods[ip].tag, ptRangesTruth[ib], ptRangesTruth[ib+1],
                    h_mean[ip]->GetBinContent(ib+1), h_mean[ip]->GetBinError(ib+1),
                    h_sigma[ip]->GetBinContent(ib+1), h_sigma[ip]->GetBinError(ib+1),
                    h_histmean[ip]->GetBinContent(ib+1), h_rms[ip]->GetBinContent(ib+1),
                    hp->Integral());
        }
    }
    fclose(fout);

    std::cout << "[plot_cluster_response] Done. Outputs in " << outdir << std::endl;
}
