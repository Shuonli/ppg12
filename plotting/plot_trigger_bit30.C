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
// Methodology follows reports/trigger_bit30_turnon.md:
//   - 10 of 77 slimtree data files (6.6M events).
//   - Fiducial: |cluster_Eta_CLUSTERINFO_CEMC| < 0.7 and |vertexz| < 60 cm.
//   - Per-cluster fill of ET = cluster_Et_CLUSTERINFO_CEMC.
//   - Denominator: livetrigger[10] == 1 (MBD NandS coincidence).
//   - Numerator:   livetrigger[10] == 1 && livetrigger[30] == 1.
//   - Binning: 0.5 GeV from 7->20, 1 GeV from 20->40 (46 variable bins).
//   - Plateau constant fit on 10-20 GeV returns ~0.9958 per the report.
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

// 0.5 GeV binning from 7 -> 20, then 1 GeV from 20 -> 40. 46 bins, 47 edges.
std::vector<double> make_turnon_edges()
{
    std::vector<double> edges;
    edges.reserve(47);
    for (int i = 0; i <= 26; ++i) edges.push_back(7.0 + 0.5 * i);
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
    eff_bit30->SetStatisticOption(TEfficiency::kFCP);  // Clopper-Pearson

    // ---- Readers ---------------------------------------------------------
    TTreeReader R(t);
    TTreeReaderValue<int>   ncluster(R, "ncluster_CLUSTERINFO_CEMC");
    TTreeReaderArray<float> cEt (R, "cluster_Et_CLUSTERINFO_CEMC");
    TTreeReaderArray<float> cEta(R, "cluster_Eta_CLUSTERINFO_CEMC");
    TTreeReaderValue<float> vtxz(R, "vertexz");
    TTreeReaderArray<bool>  livetrig(R, "livetrigger");

    // ---- Event / cluster loop -------------------------------------------
    auto tstart = std::chrono::steady_clock::now();
    long long nev = 0, nev_vz = 0, nev_mbd = 0;
    long long nclus_den = 0, nclus_num = 0;

    while (R.Next())
    {
        ++nev;
        if (std::fabs(*vtxz) > 60.0) continue;
        ++nev_vz;
        if (livetrig.GetSize() <= 30)   continue;
        if (livetrig[10] == 0)          continue;
        ++nev_mbd;
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

        h_MBDns->SetStats(0);
        h_MBDns->GetXaxis()->SetRangeUser(7, 40);
        h_MBDns->GetYaxis()->SetRangeUser(1, 5.0 * h_MBDns->GetMaximum());
        h_MBDns->GetXaxis()->SetTitle("#it{E}_{T}^{cluster} [GeV]");
        h_MBDns->GetYaxis()->SetTitle("Clusters");
        h_MBDns->GetYaxis()->SetTitleOffset(1.45);
        h_MBDns->GetXaxis()->SetTitleOffset(1.15);

        // MBD N&S drawn on top (thick solid black). Photon 4 GeV underneath
        // (thinner dashed red). Since bit-30 eff ~99.6%, the two histograms
        // sit on top of each other above the turn-on; differentiating line
        // style makes both visible at a glance.
        h_MBDns->SetLineColor(kBlack);
        h_MBDns->SetLineStyle(1);   // solid
        h_MBDns->SetLineWidth(3);
        h_MBDns->SetMarkerColor(kBlack);
        h_MBDns->SetMarkerStyle(20);
        h_MBDns->SetMarkerSize(0.7);

        h_Photon4GeV->SetLineColor(kRed + 1);
        h_Photon4GeV->SetLineStyle(2);   // dashed
        h_Photon4GeV->SetLineWidth(2);
        h_Photon4GeV->SetMarkerColor(kRed + 1);
        h_Photon4GeV->SetMarkerStyle(24);
        h_Photon4GeV->SetMarkerSize(0.7);

        // Draw Photon4GeV first (underneath), then MBDns on top so both are
        // distinguishable in the overlap region.
        h_Photon4GeV->Draw("hist");
        h_MBDns->Draw("hist same");

        TLegend *l = new TLegend(0.45, 0.66, 0.93, 0.82);
        legStyle(l, 0.20, 0.038);
        l->AddEntry(h_MBDns,      "MBD N&S #geq 1 (bit 10)",                "l");
        l->AddEntry(h_Photon4GeV, "Photon 4 GeV + MBD NS (bit 30)",         "l");
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

        TH1F *frame = new TH1F("frame_trig_turnon", "",
                               1, 7, 40);
        frame->SetStats(0);
        frame->GetXaxis()->SetRangeUser(7, 40);
        frame->GetYaxis()->SetRangeUser(0.0, 1.1);
        frame->GetXaxis()->SetTitle("#it{E}_{T}^{cluster} [GeV]");
        frame->GetYaxis()->SetTitle(
            "#varepsilon_{L1}(Photon 4 GeV | MBD N&S)");
        frame->GetYaxis()->SetTitleOffset(1.35);
        frame->GetXaxis()->SetTitleOffset(1.15);
        frame->Draw("axis");

        eff_bit30->SetMarkerStyle(20);
        eff_bit30->SetMarkerColor(kBlack);
        eff_bit30->SetMarkerSize(0.9);
        eff_bit30->SetLineColor(kBlack);
        eff_bit30->Draw("same p");

        TLegend *l = new TLegend(0.45, 0.20, 0.93, 0.32);
        legStyle(l, 0.20, 0.038);
        l->AddEntry(eff_bit30,
                    "#varepsilon = N(bit 10 & 30) / N(bit 10)",
                    "pl");
        l->Draw("same");

        const float xpos = 0.22, xpos2 = 0.93, ypos = 0.88;
        const float dy = 0.054, fontsize = 0.040;
        myText(xpos,  ypos - 0 * dy, 1, strleg1.c_str(),   fontsize, 0);
        myText(xpos,  ypos - 1 * dy, 1, strleg2_1.c_str(), fontsize, 0);
        myText(xpos2, ypos - 0 * dy, 1, strleg3.c_str(),   fontsize, 1);

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
        frame_top->GetXaxis()->SetRangeUser(7, 40);
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
