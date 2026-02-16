#include "plotcommon.h"
#include <yaml-cpp/yaml.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>

void plot_double_interaction(const std::string &suffix = "double_interaction_check")
{
    init_plot();
    string savePath = "figures";

    // Load pT bins from config
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node config = YAML::LoadFile("../efficiencytool/config_nom.yaml");
    std::vector<double> pT_bin_edges = config["analysis"]["pT_bins"].as<std::vector<double>>();
    int n_pT_bins = pT_bin_edges.size() - 1;

    // Open files
    TFile *f_data = TFile::Open(Form("../efficiencytool/results/data_histo_%s.root", suffix.c_str()));
    TFile *f_sig  = TFile::Open(Form("../efficiencytool/results/MC_efficiency_%s_signal.root", suffix.c_str()));
    TFile *f_jet  = TFile::Open(Form("../efficiencytool/results/MC_efficiency_%s_jet.root", suffix.c_str()));

    int ieta = 0;

    // ---------------------------------------------------------------
    // Plot 1: ABCD yields for NPB < 0.5 (from Task 1 TH2Ds)
    // ---------------------------------------------------------------
    auto plotABCDYields = [&](TFile *fin, const std::string &label, const std::string &tag)
    {
        if (!fin || fin->IsZombie()) return;

        TH2D *h2_A = (TH2D *)fin->Get(Form("h_tight_iso_cluster_npb_%d", ieta));
        TH2D *h2_B = (TH2D *)fin->Get(Form("h_tight_noniso_cluster_npb_%d", ieta));
        TH2D *h2_C = (TH2D *)fin->Get(Form("h_nontight_iso_cluster_npb_%d", ieta));
        TH2D *h2_D = (TH2D *)fin->Get(Form("h_nontight_noniso_cluster_npb_%d", ieta));
        if (!h2_A || !h2_B || !h2_C || !h2_D) { std::cout << "Missing TH2D for " << tag << std::endl; return; }

        // Project X for NPB < 0.5 (Y-axis bins 1 to FindBin(0.5)-1)
        int npb_cut_bin = h2_A->GetYaxis()->FindBin(0.5) - 1;
        if (npb_cut_bin < 1) npb_cut_bin = 1;

        TH1D *h_A = h2_A->ProjectionX(Form("hA_%s", tag.c_str()), 1, npb_cut_bin);
        TH1D *h_B = h2_B->ProjectionX(Form("hB_%s", tag.c_str()), 1, npb_cut_bin);
        TH1D *h_C = h2_C->ProjectionX(Form("hC_%s", tag.c_str()), 1, npb_cut_bin);
        TH1D *h_D = h2_D->ProjectionX(Form("hD_%s", tag.c_str()), 1, npb_cut_bin);

        h_A->SetLineColor(kBlack);
        h_B->SetLineColor(kRed);
        h_C->SetLineColor(kBlue);
        h_D->SetLineColor(kMagenta);

        TCanvas *c = new TCanvas(Form("c_abcd_%s", tag.c_str()), "", 600, 600);
        frame_et_rec->Draw("axis");
	//frame_et_rec->GetYaxis()->SetRangeUser(0,0.5);
        frame_et_rec->GetXaxis()->SetRangeUser(10, 35);
        float ymax = h_A->GetMaximum() * 5;
        float ymin = 0.5;
        frame_et_rec->GetYaxis()->SetRangeUser(ymin, ymax);

        h_A->Draw("same hist");
        h_B->Draw("same hist");
        h_C->Draw("same hist");
        h_D->Draw("same hist");

        myText(0.50, 0.90, 1, strleg1.c_str(), 0.04);
        myText(0.50, 0.85, 1, strleg2.c_str(), 0.04);
        myText(0.50, 0.80, 1, label.c_str(), 0.04);
        myText(0.18, 0.90, 1, "NPB < 0.5", 0.04);
        myMarkerLineText(0.55, 0.75, 0, kBlack, 0, kBlack, 1, "A: tight iso", 0.05, true);
        myMarkerLineText(0.55, 0.70, 0, kRed, 0, kRed, 1, "B: tight noniso", 0.05, true);
        myMarkerLineText(0.55, 0.65, 0, kBlue, 0, kBlue, 1, "C: nontight iso", 0.05, true);
        myMarkerLineText(0.55, 0.60, 0, kMagenta, 0, kMagenta, 1, "D: nontight noniso", 0.05, true);

        gPad->SetLogy();
        c->SaveAs(Form("%s/abcd_yields_npb_lt05_%s.pdf", savePath.c_str(), tag.c_str()));
        delete c;
    };

    plotABCDYields(f_data, "Data", "data");
    plotABCDYields(f_sig, strSigMC, "signal");
    plotABCDYields(f_jet, strIncMC, "jet");

    // ---------------------------------------------------------------
    // Plot 1b: Fraction of NPB < 0.5 vs ET for each ABCD region
    // ---------------------------------------------------------------
    auto plotNPBFraction = [&](TFile *fin, const std::string &label, const std::string &tag)
    {
        if (!fin || fin->IsZombie()) return;

        TH2D *h2_A = (TH2D *)fin->Get(Form("h_tight_iso_cluster_npb_%d", ieta));
        TH2D *h2_B = (TH2D *)fin->Get(Form("h_tight_noniso_cluster_npb_%d", ieta));
        TH2D *h2_C = (TH2D *)fin->Get(Form("h_nontight_iso_cluster_npb_%d", ieta));
        TH2D *h2_D = (TH2D *)fin->Get(Form("h_nontight_noniso_cluster_npb_%d", ieta));
        if (!h2_A || !h2_B || !h2_C || !h2_D) return;

        int npb_cut_bin = h2_A->GetYaxis()->FindBin(0.5) - 1;
        if (npb_cut_bin < 1) npb_cut_bin = 1;
        int npb_max_bin = h2_A->GetYaxis()->GetNbins();

        // NPB < 0.5 projections
        TH1D *h_A_low = h2_A->ProjectionX(Form("hA_low_%s", tag.c_str()), 1, npb_cut_bin);
        TH1D *h_B_low = h2_B->ProjectionX(Form("hB_low_%s", tag.c_str()), 1, npb_cut_bin);
        TH1D *h_C_low = h2_C->ProjectionX(Form("hC_low_%s", tag.c_str()), 1, npb_cut_bin);
        TH1D *h_D_low = h2_D->ProjectionX(Form("hD_low_%s", tag.c_str()), 1, npb_cut_bin);

        // All NPB projections (total)
        TH1D *h_A_all = h2_A->ProjectionX(Form("hA_all_%s", tag.c_str()), 1, npb_max_bin);
        TH1D *h_B_all = h2_B->ProjectionX(Form("hB_all_%s", tag.c_str()), 1, npb_max_bin);
        TH1D *h_C_all = h2_C->ProjectionX(Form("hC_all_%s", tag.c_str()), 1, npb_max_bin);
        TH1D *h_D_all = h2_D->ProjectionX(Form("hD_all_%s", tag.c_str()), 1, npb_max_bin);

        // Fraction = low / all
        TGraphAsymmErrors *g_frac_A = new TGraphAsymmErrors(h_A_low, h_A_all);
        TGraphAsymmErrors *g_frac_B = new TGraphAsymmErrors(h_B_low, h_B_all);
        TGraphAsymmErrors *g_frac_C = new TGraphAsymmErrors(h_C_low, h_C_all);
        TGraphAsymmErrors *g_frac_D = new TGraphAsymmErrors(h_D_low, h_D_all);

        g_frac_A->SetMarkerStyle(20); g_frac_A->SetMarkerColor(kBlack);   g_frac_A->SetLineColor(kBlack);
        g_frac_B->SetMarkerStyle(21); g_frac_B->SetMarkerColor(kRed);     g_frac_B->SetLineColor(kRed);
        g_frac_C->SetMarkerStyle(22); g_frac_C->SetMarkerColor(kBlue);    g_frac_C->SetLineColor(kBlue);
        g_frac_D->SetMarkerStyle(23); g_frac_D->SetMarkerColor(kMagenta); g_frac_D->SetLineColor(kMagenta);

        TCanvas *c = new TCanvas(Form("c_npbfrac_%s", tag.c_str()), "", 600, 600);
        TH1F *frame_frac = new TH1F(Form("frame_npbfrac_%s", tag.c_str()),
            ";#it{E}_{T}^{#gamma,rec} [GeV];Fraction with NPB < 0.5", 100, 8, 40);
        frame_frac->GetYaxis()->SetRangeUser(0, 0.5);
        frame_frac->GetXaxis()->SetRangeUser(10, 35);
        frame_frac->Draw("axis");

        lineone->Draw("same");

        g_frac_A->Draw("P same");
        g_frac_B->Draw("P same");
        g_frac_C->Draw("P same");
        g_frac_D->Draw("P same");

        myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
        myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
        myText(0.20, 0.80, 1, label.c_str(), 0.04);
        myMarkerLineText(0.55, 0.90, 1.0, kBlack, 20, kBlack, 1, "A: tight iso", 0.04, true);
        myMarkerLineText(0.55, 0.85, 1.0, kRed, 21, kRed, 1, "B: tight noniso", 0.04, true);
        myMarkerLineText(0.55, 0.80, 1.0, kBlue, 22, kBlue, 1, "C: nontight iso", 0.04, true);
        myMarkerLineText(0.55, 0.75, 1.0, kMagenta, 23, kMagenta, 1, "D: nontight noniso", 0.04, true);

        c->SaveAs(Form("%s/npb_lt05_fraction_%s.pdf", savePath.c_str(), tag.c_str()));
        delete c;
    };

    plotNPBFraction(f_data, "Data", "data");
    plotNPBFraction(f_sig, strSigMC, "signal");
    plotNPBFraction(f_jet, strIncMC, "jet");

    // ---------------------------------------------------------------
    // Plot 2: MBD avg sigma vs run number (data only)
    // ---------------------------------------------------------------
    if (f_data && !f_data->IsZombie())
    {
        TProfile *h_run_all    = (TProfile *)f_data->Get("h_mbd_avgsigma_vs_run_all");
        TProfile *h_run_0mrad  = (TProfile *)f_data->Get("h_mbd_avgsigma_vs_run_0mrad");
        TProfile *h_run_1p5mrad = (TProfile *)f_data->Get("h_mbd_avgsigma_vs_run_1p5mrad");

        if (h_run_all && h_run_0mrad && h_run_1p5mrad)
        {
            TCanvas *c2 = new TCanvas("c_mbd_run", "", 600, 600);

            h_run_all->SetMarkerStyle(20);
            h_run_all->SetMarkerSize(0.4);
            h_run_all->SetMarkerColor(kBlack);
            h_run_all->SetLineColor(kBlack);
            h_run_all->GetYaxis()->SetRangeUser(0, 1.0);
            h_run_all->GetXaxis()->SetTitle("Run Number");
            h_run_all->GetYaxis()->SetTitle("MBD Avg #sigma_{t} [ns]");
            h_run_all->Draw("ex0");

            h_run_0mrad->SetMarkerStyle(21);
            h_run_0mrad->SetMarkerSize(0.4);
            h_run_0mrad->SetMarkerColor(kRed);
            h_run_0mrad->SetLineColor(kRed);
            //h_run_0mrad->Draw("same ex0");

            h_run_1p5mrad->SetMarkerStyle(22);
            h_run_1p5mrad->SetMarkerSize(0.4);
            h_run_1p5mrad->SetMarkerColor(kBlue);
            h_run_1p5mrad->SetLineColor(kBlue);
            //h_run_1p5mrad->Draw("same ex0");

            myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
            myMarkerLineText(0.55, 0.90, 0.4, kBlack, 20, kBlack, 1, "All runs", 0.04, true);
            //myMarkerLineText(0.55, 0.85, 0.4, kRed, 21, kRed, 1, "0 mrad", 0.04, true);
            //myMarkerLineText(0.55, 0.80, 0.4, kBlue, 22, kBlue, 1, "1.5 mrad", 0.04, true);

            c2->SaveAs(Form("%s/mbd_avgsigma_vs_run.pdf", savePath.c_str()));
            delete c2;
        }
        else
        {
            std::cout << "Missing MBD TProfiles in data file" << std::endl;
        }
    }

    // ---------------------------------------------------------------
    // Plot 3: NPB score distribution under tight selection (iso + noniso)
    //         Data vs inclusive jet MC, wider pT bins
    // pT_bin_edges: [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35]
    // Wider bins:   [8-12, 12-16, 16-20, 20-26, 26-35]
    //               bins 1-2, 3-4, 5-6, 7-9, 10
    // ---------------------------------------------------------------
    {
        bool have_data_npb = (f_data && !f_data->IsZombie());
        bool have_jet_npb  = (f_jet && !f_jet->IsZombie());

        // Tight iso + noniso TH2Ds
        TH2D *h2_data_tight_iso    = have_data_npb ? (TH2D *)f_data->Get(Form("h_tight_iso_cluster_npb_%d", ieta))    : nullptr;
        TH2D *h2_data_tight_noniso = have_data_npb ? (TH2D *)f_data->Get(Form("h_tight_noniso_cluster_npb_%d", ieta)) : nullptr;
        TH2D *h2_jet_tight_iso     = have_jet_npb  ? (TH2D *)f_jet->Get(Form("h_tight_iso_cluster_npb_%d", ieta))     : nullptr;
        TH2D *h2_jet_tight_noniso  = have_jet_npb  ? (TH2D *)f_jet->Get(Form("h_tight_noniso_cluster_npb_%d", ieta))  : nullptr;

        if (h2_data_tight_iso || h2_jet_tight_iso)
        {
            // Wider pT bin ranges (indices into pT_bin_edges)
            std::vector<int> wide_bin_lo = {1, 3, 5, 7, 10};  // X-axis bin numbers (1-indexed)
            std::vector<int> wide_bin_hi = {2, 4, 6, 9, 10};
            std::vector<double> wide_pt_lo = {8, 12, 16, 20, 26};
            std::vector<double> wide_pt_hi = {12, 16, 20, 26, 35};
            int n_wide = wide_bin_lo.size();

            for (int iw = 0; iw < n_wide; iw++)
            {
                // Project tight iso + noniso for data
                TH1D *h_npb_data = nullptr;
                if (h2_data_tight_iso)
                {
                    h_npb_data = h2_data_tight_iso->ProjectionY(Form("h_npb_data_w%d", iw), wide_bin_lo[iw], wide_bin_hi[iw]);
                    if (h2_data_tight_noniso)
                    {
                        TH1D *htmp = h2_data_tight_noniso->ProjectionY(Form("h_npb_data_noniso_w%d", iw), wide_bin_lo[iw], wide_bin_hi[iw]);
                        h_npb_data->Add(htmp);
                        delete htmp;
                    }
                }

                // Project tight iso + noniso for jet MC
                TH1D *h_npb_jet = nullptr;
                if (h2_jet_tight_iso)
                {
                    h_npb_jet = h2_jet_tight_iso->ProjectionY(Form("h_npb_jet_w%d", iw), wide_bin_lo[iw], wide_bin_hi[iw]);
                    if (h2_jet_tight_noniso)
                    {
                        TH1D *htmp = h2_jet_tight_noniso->ProjectionY(Form("h_npb_jet_noniso_w%d", iw), wide_bin_lo[iw], wide_bin_hi[iw]);
                        h_npb_jet->Add(htmp);
                        delete htmp;
                    }
                }

                // Normalize to unit area
                if (h_npb_data && h_npb_data->Integral() > 0) h_npb_data->Scale(1.0 / h_npb_data->Integral());
                if (h_npb_jet  && h_npb_jet->Integral() > 0)  h_npb_jet->Scale(1.0 / h_npb_jet->Integral());

                float ymax3 = 0;
                if (h_npb_data) ymax3 = h_npb_data->GetMaximum();
                if (h_npb_jet && h_npb_jet->GetMaximum() > ymax3) ymax3 = h_npb_jet->GetMaximum();
                ymax3 *= 1.4;

                TCanvas *c3 = new TCanvas(Form("c_npb_w%d", iw), "", 600, 600);
                TH1F *frame_npb = new TH1F(Form("frame_npb_w%d", iw), ";NPB Score;Normalized counts", 50, 0.7, 1);
                gPad->SetLogy();
                frame_npb->GetYaxis()->SetRangeUser(0.01, ymax3);
                frame_npb->Draw("axis");

                if (h_npb_data)
                {
                    h_npb_data->SetLineColor(kBlack);
                    h_npb_data->SetLineWidth(2);
                    h_npb_data->Draw("same hist");
                }
                if (h_npb_jet)
                {
                    h_npb_jet->SetLineColor(kRed);
                    h_npb_jet->SetLineWidth(2);
                    h_npb_jet->Draw("same hist");
                }

                myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                myText(0.20, 0.80, 1, Form("%0.0f < #it{E}_{T}^{#gamma,rec} < %0.0f GeV", wide_pt_lo[iw], wide_pt_hi[iw]), 0.04);
                myText(0.20, 0.75, 1, "Tight (iso + noniso)", 0.04);
                myMarkerLineText(0.55, 0.90, 0, kBlack, 0, kBlack, 1, "Data", 0.05, true);
                myMarkerLineText(0.55, 0.85, 0, kRed, 0, kRed, 1, "Inclusive MC", 0.05, true);

                c3->SaveAs(Form("%s/npb_score_tight_pt%d.pdf", savePath.c_str(), iw));
                delete c3;
            }
        }
        else
        {
            std::cout << "Missing NPB TH2Ds for Plot 3" << std::endl;
        }
    }

    // ---------------------------------------------------------------
    // Plot 4: Truth purity single vs double interaction (inclusive jet MC)
    // Uses _incsig histograms (signal tagged for all MC including jet)
    // Both numerator and denominator from the inclusive jet MC file
    // Overlays additional smearing levels (5, 10, 15, 20 cm)
    // ---------------------------------------------------------------
    bool have_jet = (f_jet && !f_jet->IsZombie());
    if (have_jet)
    {
        TH1D *h_single_signal = (TH1D *)f_jet->Get(Form("h_tight_iso_cluster_single_incsig_%d", ieta));
        TH1D *h_single_total  = (TH1D *)f_jet->Get(Form("h_tight_iso_cluster_single_%d", ieta));
        TH1D *h_double_signal = (TH1D *)f_jet->Get(Form("h_tight_iso_cluster_double_incsig_%d", ieta));
        TH1D *h_double_total  = (TH1D *)f_jet->Get(Form("h_tight_iso_cluster_double_%d", ieta));

        if (h_single_total && h_single_signal && h_double_total && h_double_signal)
        {
            TGraphAsymmErrors *g_purity_single = new TGraphAsymmErrors(h_single_signal, h_single_total);
            g_purity_single->SetName("g_purity_single_plot");

            TGraphAsymmErrors *g_purity_double = new TGraphAsymmErrors(h_double_signal, h_double_total);
            g_purity_double->SetName("g_purity_double_plot");

            g_purity_single->SetMarkerStyle(20);
            g_purity_single->SetMarkerSize(1.0);
            g_purity_single->SetMarkerColor(kBlack);
            g_purity_single->SetLineColor(kBlack);

            g_purity_double->SetMarkerStyle(21);
            g_purity_double->SetMarkerSize(1.0);
            g_purity_double->SetMarkerColor(kRed);
            g_purity_double->SetLineColor(kRed);

            // Smeared double interaction purity curves
            std::vector<int> smear_sigmas = {5, 10, 15, 20};
            int smear_markers[] = {22, 23, 33, 34};
            int smear_colors[]  = {kBlue, kGreen + 2, kMagenta, kOrange + 1};
            std::vector<TGraphAsymmErrors *> g_purity_smear;

            for (int ism = 0; ism < (int)smear_sigmas.size(); ism++)
            {
                TH1D *h_sm_signal = (TH1D *)f_jet->Get(Form("h_tight_iso_cluster_double_smear%d_incsig_%d", smear_sigmas[ism], ieta));
                TH1D *h_sm_total  = (TH1D *)f_jet->Get(Form("h_tight_iso_cluster_double_smear%d_%d", smear_sigmas[ism], ieta));
                if (h_sm_signal && h_sm_total && h_sm_total->Integral() > 0)
                {
                    TGraphAsymmErrors *g = new TGraphAsymmErrors(h_sm_signal, h_sm_total);
                    g->SetName(Form("g_purity_smear%d", smear_sigmas[ism]));
                    g->SetMarkerStyle(smear_markers[ism]);
                    g->SetMarkerSize(1.0);
                    g->SetMarkerColor(smear_colors[ism]);
                    g->SetLineColor(smear_colors[ism]);
                    g_purity_smear.push_back(g);
                }
                else
                {
                    g_purity_smear.push_back(nullptr);
                }
            }

            TCanvas *c4 = new TCanvas("c_purity", "", 600, 600);
            TH1F *frame_purity = new TH1F("frame_purity", ";#it{E}_{T}^{#gamma,rec} [GeV];Truth Purity", 100, 8, 40);
            frame_purity->GetYaxis()->SetRangeUser(0, 1.2);
            frame_purity->GetXaxis()->SetRangeUser(10, 35);
            frame_purity->Draw("axis");

            lineone->Draw("same");

            g_purity_single->Draw("P same");
            g_purity_double->Draw("P same");
            for (auto *g : g_purity_smear)
            {
                if (g) g->Draw("P same");
            }

            myText(0.20, 0.40, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.35, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.30, 1, strleg3.c_str(), 0.04);
            myText(0.20, 0.25, 1, strIncMC.c_str(), 0.04);

            float leg_y = 0.90;
            myMarkerLineText(0.55, leg_y, 1.0, kBlack, 20, kBlack, 1, "Single interaction", 0.04, true);
            leg_y -= 0.05;
            myMarkerLineText(0.55, leg_y, 1.0, kRed, 21, kRed, 1, "Double (no smear)", 0.04, true);
            for (int ism = 0; ism < (int)smear_sigmas.size(); ism++)
            {
                if (g_purity_smear[ism])
                {
                    leg_y -= 0.05;
                    myMarkerLineText(0.55, leg_y, 1.0, smear_colors[ism], smear_markers[ism],
                        smear_colors[ism], 1, Form("Double +%d cm smear", smear_sigmas[ism]), 0.04, true);
                }
            }

            c4->SaveAs(Form("%s/truth_purity_single_vs_double.pdf", savePath.c_str()));
            delete c4;
        }
        else
        {
            std::cout << "Missing single/double or incsig histograms in jet MC file" << std::endl;
        }
    }

    std::cout << "Done. Plots saved to " << savePath << "/" << std::endl;
}
