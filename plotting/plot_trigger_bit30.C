// plot_trigger_bit30.C
//
// Measure and plot the L1 Photon-4-GeV trigger (bit 30,
// `Photon_4_GeV_plus_MBD_NS`) turn-on vs cluster ET, conditioned on
// bit 10 (`MBD_N&S_geq_1`). Produces the 3 PDFs used in PPG12
// analysis-note Section 2 (trigger efficiency):
//
//   figures/h_maxEnergyClus_NewTriggerFilling_doNotScale_Overlay.pdf
//   figures/h_maxEnergyClus_NewTriggerFilling_doNotScale_PhotonTurnOn.pdf
//   figures/Photon_4_GeV_ConstantFit.pdf
//
// Methodology:
//   - 10 of 77 slimtree data files (6.6M events).
//   - Fiducial: |cluster_Eta_CLUSTERINFO_CEMC| < 0.7 and |vertexz| < 200 cm.
//   - Per-cluster fill of ET = cluster_Et_CLUSTERINFO_CEMC.
//   - Denominator: scaledtrigger[10] == 1 (MBD N&S, prescale-selected reference).
//   - Numerator:   scaledtrigger[10] == 1 && livetrigger[30] == 1
//                  (bit 30 fired at L1 in the unbiased bit-10 sample;
//                   independent of bit-30's own prescale wheel).
//   - Binning: 0.5 GeV from 5->20, 1 GeV from 20->40 (50 variable bins).
//
// Run:  root -l -b -q plot_trigger_bit30.C

#include "plotcommon.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <chrono>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace {

// 0.5 GeV binning from 5 -> 20, then 1 GeV from 20 -> 40. 50 bins, 51 edges.
// Bin width is non-uniform; plotting code below scales raw counts by
// 1/binwidth before drawing the left-panel overlay so the bin-width
// transition at 20 GeV does not produce a step-up artifact.
std::vector<double> make_turnon_edges()
{
    std::vector<double> edges;
    edges.reserve(51);
    for (int i = 0; i <= 30; ++i) edges.push_back(5.0 + 0.5 * i);
    for (int i = 1; i <= 20; ++i) edges.push_back(20.0 + 1.0 * i);
    return edges;
}

}  // namespace

void plot_trigger_bit30()
{
    init_plot();

    const std::string savePath = "figures";
    const std::string indir = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/"
                              "anatreemaker/macro_maketree/data/ana521/condorout";
    const int kNfiles = 10;  // per reports/trigger_bit30_turnon.md

    TChain *t = new TChain("slimtree");
    for (int i = 0; i < kNfiles; ++i)
    {
        t->Add(Form("%s/part_%d_with_bdt_split.root", indir.c_str(), i));
    }
    std::cout << "added " << kNfiles << " files; entries = "
              << t->GetEntries() << std::endl;

    // ---- Binning ---------------------------------------------------------
    auto edges_vec = make_turnon_edges();
    const int nbins = static_cast<int>(edges_vec.size()) - 1;
    std::vector<double> edges = edges_vec;  // keep storage alive

    TH1F *h_MBDns = new TH1F("h_MBDns",
                              ";#it{E}_{T}^{cluster} [GeV];Clusters",
                              nbins, edges.data());
    TH1F *h_Photon4GeV = new TH1F("h_Photon4GeV",
                                   ";#it{E}_{T}^{cluster} [GeV];Clusters",
                                   nbins, edges.data());

    // TEfficiency with variable binning.
    TEfficiency *eff_bit30 = new TEfficiency(
        "eff_bit30",
        ";#it{E}_{T}^{cluster} [GeV];#varepsilon(bit 30 | bit 10)",
        nbins, edges.data());
    eff_bit30->SetStatisticOption(TEfficiency::kFCP);  // Clopper-Pearson 68% CI

    // ---- Readers ---------------------------------------------------------
    TTreeReader R(t);
    TTreeReaderValue<int>   ncluster(R, "ncluster_CLUSTERINFO_CEMC");
    TTreeReaderArray<float> cEt (R, "cluster_Et_CLUSTERINFO_CEMC");
    TTreeReaderArray<float> cEta(R, "cluster_Eta_CLUSTERINFO_CEMC");
    TTreeReaderValue<float> vtxz(R, "vertexz");
    TTreeReaderArray<bool>  scaledtrig(R, "scaledtrigger");
    TTreeReaderArray<bool>  livetrig(R, "livetrigger");

    // ---- Event / cluster loop -------------------------------------------
    auto tstart = std::chrono::steady_clock::now();
    long long nev = 0, nev_vz = 0, nev_mbd = 0;
    long long nclus_den = 0, nclus_num = 0;

    while (R.Next())
    {
        ++nev;
        if (std::fabs(*vtxz) > 200.0) continue;
        ++nev_vz;
        if (scaledtrig.GetSize() <= 30)   continue;
        if (livetrig.GetSize()   <= 30)   continue;
        if (scaledtrig[10] == 0)          continue;
        ++nev_mbd;
        // Numerator: bit 30 fired at L1 in the unbiased scaled[10] reference.
        // Use livetrigger[30] (not scaledtrigger[30]) so the result is
        // independent of bit-30's own prescale wheel.
        const bool photon_on = (livetrig[30] != 0);
        for (int i = 0; i < *ncluster; ++i)
        {
            const float et  = cEt[i];
            const float eta = cEta[i];
            if (std::fabs(eta) >= 0.7) continue;
            if (et < edges.front() || et >= edges.back()) continue;
            ++nclus_den;
            h_MBDns->Fill(et);
            if (photon_on)
            {
                ++nclus_num;
                h_Photon4GeV->Fill(et);
            }
            eff_bit30->Fill(photon_on, et);
        }
    }

    auto tend = std::chrono::steady_clock::now();
    const double dt_sec = std::chrono::duration<double>(tend - tstart).count();

    std::cout << "events total      = " << nev     << "\n"
              << "events |vz|<60    = " << nev_vz  << "\n"
              << "events + bit10 on = " << nev_mbd << "\n"
              << "clusters den      = " << nclus_den << "\n"
              << "clusters num      = " << nclus_num << "\n"
              << "integrated eff    = "
              << (nclus_den > 0 ? double(nclus_num) / nclus_den : 0.0)
              << std::endl;
    std::cout << "reader loop time  = " << dt_sec << " s" << std::endl;

    // =====================================================================
    // PDF 1: Overlay (log-y raw counts)
    // =====================================================================
    {
        TCanvas *c1 = new TCanvas("c1_trig_overlay",
                                  "Trigger counts overlay", 600, 600);
        c1->SetLogy();
        c1->SetLeftMargin(0.17);
        c1->SetRightMargin(0.04);
        c1->SetTopMargin(0.06);
        c1->SetBottomMargin(0.15);

        // Density scaling: divide raw counts by bin width so the 20 GeV
        // bin-width transition (0.5 GeV -> 1.0 GeV) does not produce a
        // step-up artifact in the log-y display. Sumw2 is set so the
        // per-bin Poisson errors are propagated correctly through Scale.
        TH1F *h_MBDns_density       = (TH1F *)h_MBDns->Clone("h_MBDns_density");
        TH1F *h_Photon4GeV_density  = (TH1F *)h_Photon4GeV->Clone("h_Photon4GeV_density");
        h_MBDns_density->Sumw2();
        h_Photon4GeV_density->Sumw2();
        h_MBDns_density->Scale(1.0, "width");
        h_Photon4GeV_density->Scale(1.0, "width");

        h_MBDns_density->SetStats(0);
        h_MBDns_density->GetXaxis()->SetRangeUser(5, 15);
        h_MBDns_density->GetYaxis()->SetRangeUser(
            0.5, 5.0 * h_MBDns_density->GetMaximum());
        h_MBDns_density->GetXaxis()->SetTitle("#it{E}_{T}^{cluster} [GeV]");
        h_MBDns_density->GetYaxis()->SetTitle("Clusters / GeV");
        h_MBDns_density->GetYaxis()->SetTitleOffset(1.45);
        h_MBDns_density->GetXaxis()->SetTitleOffset(1.15);

        // MBD N&S in solid black (drawn first as the wider denominator);
        // Photon 4 GeV dashed red drawn on top so the dashes cut through
        // the black line at the ~99.6% plateau where the two histograms
        // are nearly identical.
        h_MBDns_density->SetLineColor(kBlack);
        h_MBDns_density->SetLineStyle(1);
        h_MBDns_density->SetLineWidth(2);
        h_MBDns_density->SetMarkerColor(kBlack);
        h_MBDns_density->SetMarkerStyle(20);
        h_MBDns_density->SetMarkerSize(0.7);

        h_Photon4GeV_density->SetLineColor(kRed + 1);
        h_Photon4GeV_density->SetLineStyle(2);   // dashed
        h_Photon4GeV_density->SetLineWidth(3);
        h_Photon4GeV_density->SetMarkerColor(kRed + 1);
        h_Photon4GeV_density->SetMarkerStyle(24);
        h_Photon4GeV_density->SetMarkerSize(0.7);

        // Draw with E1 to show statistical (Poisson) error bars.
        h_MBDns_density->Draw("E1");
        h_Photon4GeV_density->Draw("E1 same");

        TLegend *l = new TLegend(0.42, 0.66, 0.95, 0.82);
        legStyle(l, 0.20, 0.038);
        l->AddEntry(h_MBDns_density,      "scaled[10]  (MBD N&S)",        "lep");
        l->AddEntry(h_Photon4GeV_density, "scaled[10] & live[30]",        "lep");
        l->Draw("same");

        const float xpos = 0.22, xpos2 = 0.93, ypos = 0.88;
        const float dy = 0.054, fontsize = 0.040;
        myText(xpos,  ypos - 0 * dy, 1, strleg1.c_str(),   fontsize, 0);
        myText(xpos,  ypos - 1 * dy, 1, strleg2_1.c_str(), fontsize, 0);
        myText(xpos2, ypos - 0 * dy, 1, strleg3.c_str(),   fontsize, 1);

        c1->SaveAs(Form("%s/"
                        "h_maxEnergyClus_NewTriggerFilling_doNotScale_"
                        "Overlay.pdf",
                        savePath.c_str()));
    }

    // =====================================================================
    // PDF 2: Turn-on (TEfficiency)
    // =====================================================================
    {
        TCanvas *c2 = new TCanvas("c2_trig_turnon",
                                  "Bit-30 turn-on", 600, 600);
        c2->SetLeftMargin(0.17);
        c2->SetRightMargin(0.04);
        c2->SetTopMargin(0.06);
        c2->SetBottomMargin(0.15);

        // TH2F frame so the x-range is enforced for the TEfficiency overlay
        // (TH1F + SetRangeUser is ignored by TEfficiency::Draw("same")).
        TH2F *frame = new TH2F("frame_trig_turnon", "",
                                100, 5, 15, 100, 0.0, 1.1);
        frame->SetStats(0);
        frame->GetXaxis()->SetTitle("#it{E}_{T}^{cluster} [GeV]");
        frame->GetYaxis()->SetTitle(
            "#varepsilon_{L1}(Photon 4 GeV | MBD N&S)");
        frame->GetYaxis()->SetTitleOffset(1.35);
        frame->GetXaxis()->SetTitleOffset(1.15);
        frame->GetYaxis()->SetNdivisions(510);
        frame->Draw("axis");

        eff_bit30->SetMarkerStyle(20);
        eff_bit30->SetMarkerColor(kBlack);
        eff_bit30->SetMarkerSize(0.9);
        eff_bit30->SetLineColor(kBlack);
        eff_bit30->Draw("same p");

        // Gumbel-CDF turn-on fit:
        //   eps(E_T) = p0 * exp( -exp( -(E_T - p1) / p2 ) )
        // p0 = plateau, p1 = location (mode), p2 = scale (>0).
        // Asymmetric turn-on (slow rise at low ET, fast asymptote) often
        // describes calo trigger turn-ons better than the symmetric erf.
        // Fit on the visible range 5 < E_T^cluster < 15 GeV.
        TF1 *ferf = new TF1("fgumbel_trig_turnon",
            "[0]*TMath::Exp(-TMath::Exp(-(x-[1])/[2]))",
            5.0, 15.0);
        // All three parameters free; initial guess from the previous
        // plateau-locked fit (mu = -2, beta = 3.33).
        ferf->SetParameters(1.0, -2.0, 3.33); // (plateau, mu, beta)
        ferf->SetParLimits(0, 0.5, 1.0);     // plateau capped at 1 (physical upper bound)
        ferf->SetParLimits(1, -50.0, 8.0);
        ferf->SetParLimits(2, 0.1, 50.0);
        ferf->SetLineColor(kRed + 1);
        ferf->SetLineWidth(2);

        // Build a TH1 of efficiencies with symmetrized CP errors, then fit.
        const auto edges_local = make_turnon_edges();
        const int    nb_local  = static_cast<int>(edges_local.size()) - 1;
        TH1F *h_eff_for_fit = new TH1F("h_eff_for_fit_turnon", "",
                                       nb_local, edges_local.data());
        for (int i = 1; i <= nb_local; ++i) {
            const double tot = eff_bit30->GetTotalHistogram()
                                          ->GetBinContent(i);
            if (tot <= 0) continue;
            const double e   = eff_bit30->GetEfficiency(i);
            const double el  = eff_bit30->GetEfficiencyErrorLow(i);
            const double eh  = eff_bit30->GetEfficiencyErrorUp(i);
            const double err = 0.5 * (el + eh);
            h_eff_for_fit->SetBinContent(i, e);
            h_eff_for_fit->SetBinError(i, err > 0 ? err : 1.0/std::sqrt(tot));
        }
        TFitResultPtr fit_result = h_eff_for_fit->Fit(ferf, "R Q N S");

        // 1-sigma (68.3%) CI band from the fit covariance, evaluated on a
        // fine grid over the fit range and propagated via
        // TFitResult::GetConfidenceIntervals.
        const int n_ci = 200;
        std::vector<double> x_ci_arr(n_ci), ci_err_arr(n_ci);
        for (int i = 0; i < n_ci; ++i) {
            x_ci_arr[i] = 5.0 + i * (15.0 - 5.0) / (n_ci - 1);
        }
        fit_result->GetConfidenceIntervals(n_ci, 1, 1,
                                           x_ci_arr.data(), ci_err_arr.data(),
                                           0.683, false);
        TGraphErrors *g_ci = new TGraphErrors(n_ci);
        for (int i = 0; i < n_ci; ++i) {
            g_ci->SetPoint(i, x_ci_arr[i], ferf->Eval(x_ci_arr[i]));
            g_ci->SetPointError(i, 0.0, ci_err_arr[i]);
        }
        g_ci->SetFillColorAlpha(kRed + 1, 0.30);
        g_ci->SetFillStyle(1001);
        g_ci->SetLineWidth(0);
        g_ci->Draw("3 same");
        ferf->Draw("same");

        const double p0 = ferf->GetParameter(0), e0 = ferf->GetParError(0);
        const double p1 = ferf->GetParameter(1), e1 = ferf->GetParError(1);
        const double p2 = ferf->GetParameter(2), e2 = ferf->GetParError(2);
        std::cout << "Gumbel fit: plateau=" << p0 << " ± " << e0
                  << "  mu=" << p1 << " ± " << e1
                  << "  beta=" << p2 << " ± " << e2
                  << "  chi2/ndf=" << ferf->GetChisquare()
                  << "/" << ferf->GetNDF() << std::endl;

        TLegend *l = new TLegend(0.50, 0.36, 0.93, 0.48);
        legStyle(l, 0.20, 0.036);
        l->AddEntry(eff_bit30, "#varepsilon = N(bit 10 & 30) / N(bit 10)", "pl");
        l->AddEntry(ferf,      "fit (1#sigma CI)",                         "lf");
        l->GetListOfPrimitives()->Last(); // place CI marker after fit
        // To show both line and fill in the legend entry, replace with
        // explicit dual-icon construction:
        l->GetListOfPrimitives()->RemoveLast();
        l->AddEntry(g_ci,      "fit (1#sigma CI)",                         "lf");
        l->Draw("same");

        // Fit parameters in middle of panel (data has empty mid-y region
        // because eps jumps quickly from ~0.9 at low ET to ~1.0 plateau).
        myText(0.50, 0.66, 1, Form("plateau = %.3f #pm %.3f",   p0, e0), 0.034, 0);
        myText(0.50, 0.62, 1, Form("#mu = %.2f #pm %.2f GeV",   p1, e1), 0.034, 0);
        myText(0.50, 0.58, 1, Form("#beta = %.2f #pm %.2f GeV", p2, e2), 0.034, 0);

        // sPHENIX labels at bottom-left, fit params at top-right (data
        // saturates near y=1 so top area is also empty above the fit curve).
        const float xpos = 0.22, ypos = 0.30;
        const float dy = 0.054, fontsize = 0.040;
        myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),   fontsize, 0);
        myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(), fontsize, 0);
        myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),   fontsize, 0);

        c2->SaveAs(Form("%s/"
                        "h_maxEnergyClus_NewTriggerFilling_doNotScale_"
                        "PhotonTurnOn.pdf",
                        savePath.c_str()));
    }

    // =====================================================================
    // PDF 3: Two-panel (top = full turn-on, bottom = plateau constant fit)
    // =====================================================================
    double plateau_val = 0.0;
    double plateau_err = 0.0;
    {
        TCanvas *c3 = new TCanvas("c3_trig_plateau",
                                  "Bit-30 plateau fit", 600, 800);
        c3->SetLeftMargin(0.0);
        c3->SetRightMargin(0.0);
        c3->SetTopMargin(0.0);
        c3->SetBottomMargin(0.0);

        TPad *ptop = new TPad("ptop", "", 0.0, 0.52, 1.0, 1.00);
        ptop->SetLeftMargin(0.17);
        ptop->SetRightMargin(0.04);
        ptop->SetTopMargin(0.06);
        ptop->SetBottomMargin(0.16);
        ptop->Draw();

        TPad *pbot = new TPad("pbot", "", 0.0, 0.00, 1.0, 0.52);
        pbot->SetLeftMargin(0.17);
        pbot->SetRightMargin(0.04);
        pbot->SetTopMargin(0.06);
        pbot->SetBottomMargin(0.17);
        pbot->Draw();

        // ---- Top panel: full turn-on ------------------------------------
        ptop->cd();
        TH1F *frame_top = new TH1F("frame_trig_fit_top", "",
                                   1, 7, 40);
        frame_top->SetStats(0);
        frame_top->GetXaxis()->SetRangeUser(10, 40);
        frame_top->GetYaxis()->SetRangeUser(0.0, 1.1);
        frame_top->GetXaxis()->SetTitle("#it{E}_{T}^{cluster} [GeV]");
        frame_top->GetYaxis()->SetTitle(
            "#varepsilon_{L1}");
        frame_top->GetYaxis()->SetTitleOffset(1.25);
        frame_top->GetXaxis()->SetTitleSize(0.055);
        frame_top->GetXaxis()->SetLabelSize(0.050);
        frame_top->GetYaxis()->SetTitleSize(0.055);
        frame_top->GetYaxis()->SetLabelSize(0.050);
        frame_top->Draw("axis");

        eff_bit30->Draw("same p");

        const float xpos = 0.22, xpos2 = 0.93, ypos = 0.88;
        const float dy = 0.07, fontsize = 0.050;
        myText(xpos,  ypos - 0 * dy, 1, strleg1.c_str(),   fontsize, 0);
        myText(xpos,  ypos - 1 * dy, 1, strleg2_1.c_str(), fontsize, 0);
        myText(xpos2, ypos - 0 * dy, 1, strleg3.c_str(),   fontsize, 1);

        // ---- Bottom panel: zoom 10-20 GeV + constant fit ----------------
        pbot->cd();
        TH1F *frame_bot = new TH1F("frame_trig_fit_bot", "",
                                   1, 10, 20);
        frame_bot->SetStats(0);
        frame_bot->GetXaxis()->SetRangeUser(10, 20);
        frame_bot->GetYaxis()->SetRangeUser(0.985, 1.005);
        frame_bot->GetXaxis()->SetTitle("#it{E}_{T}^{cluster} [GeV]");
        frame_bot->GetYaxis()->SetTitle(
            "#varepsilon_{L1}(plateau)");
        frame_bot->GetYaxis()->SetTitleOffset(1.20);
        frame_bot->GetXaxis()->SetTitleSize(0.055);
        frame_bot->GetXaxis()->SetLabelSize(0.050);
        frame_bot->GetYaxis()->SetTitleSize(0.055);
        frame_bot->GetYaxis()->SetLabelSize(0.050);
        frame_bot->GetYaxis()->SetNdivisions(505);
        frame_bot->Draw("axis");

        eff_bit30->Draw("same p");

        // Constant fit on 10 < ET < 20 GeV.  Build a TH1F of passed/total
        // per bin with binomial errors and fit with [0].
        TH1F *h_eff_plot = new TH1F("h_eff_plot", "",
                                    nbins, edges.data());
        for (int i = 1; i <= nbins; ++i)
        {
            const double tot = eff_bit30->GetTotalHistogram()
                                         ->GetBinContent(i);
            if (tot <= 0) continue;
            const double eff = eff_bit30->GetEfficiency(i);
            const double err = 0.5 * (eff_bit30->GetEfficiencyErrorLow(i) +
                                      eff_bit30->GetEfficiencyErrorUp(i));
            h_eff_plot->SetBinContent(i, eff);
            h_eff_plot->SetBinError(i,   err > 0 ? err : 1.0 / sqrt(tot));
        }

        TF1 *fconst = new TF1("fconst", "[0]", 10.0, 20.0);
        fconst->SetParameter(0, 0.996);
        fconst->SetLineColor(kRed + 1);
        fconst->SetLineWidth(2);
        h_eff_plot->Fit(fconst, "R Q N");
        fconst->Draw("same");

        plateau_val = fconst->GetParameter(0);
        plateau_err = fconst->GetParError(0);

        TLegend *l = new TLegend(0.22, 0.20, 0.80, 0.38);
        legStyle(l, 0.20, 0.045);
        l->SetHeader("plateau fit, 10 < #it{E}_{T} < 20 GeV", "L");
        l->AddEntry(fconst,
                    Form("constant = %.4f #pm %.4f",
                         plateau_val, plateau_err),
                    "l");
        l->Draw("same");

        std::cout << "plateau constant fit: "
                  << plateau_val << " +/- " << plateau_err
                  << "  (chi2/ndf = " << fconst->GetChisquare()
                  << " / " << fconst->GetNDF() << ")" << std::endl;

        c3->SaveAs(Form("%s/Photon_4_GeV_ConstantFit.pdf",
                        savePath.c_str()));
    }

    std::cout << "DONE. Plateau = " << plateau_val
              << " +/- " << plateau_err << std::endl;
}
