#include "plotcommon.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <THStack.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLatex.h>

#include <iostream>
#include <iomanip>

void plot_eta_migration(
    const std::string &single_file = "../efficiencytool/results/eta_migration_signal_nom.root",
    const std::string &double_file = "../efficiencytool/results/eta_migration_double_signal_nom.root",
    const std::string &background_file = "../efficiencytool/results/eta_migration_background_nom.root",
    const std::string &background_double_file = "")
{
    init_plot();
    gSystem->Exec("mkdir -p figures");

    TFile *f_single = TFile::Open(single_file.c_str());
    TFile *f_double = TFile::Open(double_file.c_str());
    TFile *f_bkg    = TFile::Open(background_file.c_str());
    TFile *f_bkg_d  = background_double_file.empty()
                          ? nullptr
                          : TFile::Open(background_double_file.c_str());

    if (!f_single || f_single->IsZombie())
    {
        std::cout << "ERROR: Cannot open single-interaction file: " << single_file << std::endl;
        return;
    }

    bool have_double = (f_double && !f_double->IsZombie());
    if (!have_double)
    {
        std::cout << "WARNING: Cannot open double-interaction file: " << double_file << std::endl;
        std::cout << "         Will skip double-interaction overlays." << std::endl;
    }

    bool have_bkg = (f_bkg && !f_bkg->IsZombie());
    if (!have_bkg)
    {
        std::cout << "WARNING: Cannot open background file: " << background_file << std::endl;
        std::cout << "         Will skip background sideband flow plots." << std::endl;
    }

    bool have_bkg_double = (f_bkg_d && !f_bkg_d->IsZombie());
    if (!have_bkg_double)
    {
        std::cout << "WARNING: Cannot open background double-interaction file: " << background_double_file << std::endl;
        std::cout << "         Will skip background double-int overlays." << std::endl;
    }

    const float eta_fid = 0.7;  // fiducial eta boundary

    // Helper: draw dashed lines at +/-0.7 on both axes of a 2D plot
    auto drawEtaLines2D = [&]()
    {
        TLine *l1 = new TLine(-eta_fid, -1.5, -eta_fid, 1.5);
        TLine *l2 = new TLine( eta_fid, -1.5,  eta_fid, 1.5);
        TLine *l3 = new TLine(-1.5, -eta_fid, 1.5, -eta_fid);
        TLine *l4 = new TLine(-1.5,  eta_fid, 1.5,  eta_fid);
        for (auto *l : {l1, l2, l3, l4})
        {
            l->SetLineStyle(7);
            l->SetLineWidth(2);
            l->SetLineColor(kRed);
            l->Draw("same");
        }
    };

    // Helper: draw vertical lines at +/-0.7 on a 1D plot
    auto drawEtaLines1D = [&](float ymin, float ymax)
    {
        TLine *l1 = new TLine(-eta_fid, ymin, -eta_fid, ymax);
        TLine *l2 = new TLine( eta_fid, ymin,  eta_fid, ymax);
        for (auto *l : {l1, l2})
        {
            l->SetLineStyle(7);
            l->SetLineWidth(2);
            l->SetLineColor(kRed);
            l->Draw("same");
        }
    };

    // =================================================================
    // Plot 1: 2D Migration Matrix (single interaction)
    // =================================================================
    {
        TH2D *h2 = (TH2D *)f_single->Get("h2_eta_truth_vs_reco");
        if (h2)
        {
            TCanvas *c = new TCanvas("c_migration_single", "", 900, 800);
            c->SetRightMargin(0.15);

            h2->SetStats(0);
            h2->SetTitle("");
            h2->GetXaxis()->SetTitle("Reco #eta");
            h2->GetYaxis()->SetTitle("Truth #eta");
            h2->GetXaxis()->SetRangeUser(-1.2, 1.2);
            h2->GetYaxis()->SetRangeUser(-1.2, 1.2);
            h2->Draw("COLZ");

            drawEtaLines2D();

            myText(0.20, 0.92, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.87, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.82, 1, "Single Interaction", 0.035);

            c->SaveAs("figures/eta_migration_matrix.pdf");
            delete c;
        }
        else
        {
            std::cout << "Missing h2_eta_truth_vs_reco in single file" << std::endl;
        }
    }

    // =================================================================
    // Plot 2: 2D Migration Matrix (double interaction)
    // =================================================================
    if (have_double)
    {
        TH2D *h2 = (TH2D *)f_double->Get("h2_eta_truth_vs_reco_double");
        if (!h2) h2 = (TH2D *)f_double->Get("h2_eta_truth_vs_reco");
        if (h2)
        {
            TCanvas *c = new TCanvas("c_migration_double", "", 900, 800);
            c->SetRightMargin(0.15);

            h2->SetStats(0);
            h2->SetTitle("");
            h2->GetXaxis()->SetTitle("Reco #eta");
            h2->GetYaxis()->SetTitle("Truth #eta");
            h2->GetXaxis()->SetRangeUser(-1.2, 1.2);
            h2->GetYaxis()->SetRangeUser(-1.2, 1.2);
            h2->Draw("COLZ");

            drawEtaLines2D();

            myText(0.20, 0.92, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.87, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.82, 1, "Double Interaction", 0.035);

            c->SaveAs("figures/eta_migration_matrix_double.pdf");
            delete c;
        }
        else
        {
            std::cout << "Missing h2_eta_truth_vs_reco_double in double file" << std::endl;
        }
    }

    // =================================================================
    // Plot 3: Inward Migration Fake Rate vs ET (single vs double)
    //   fake rate = truth outside fiducial / all matched reco
    // =================================================================
    {
        TH1D *h_fake_num_s = (TH1D *)f_single->Get("h_reco_ET_truth_outside");
        TH1D *h_fake_den_s = (TH1D *)f_single->Get("h_reco_ET_matched_all");

        TH1D *h_fake_num_d = have_double ? (TH1D *)f_double->Get("h_reco_ET_truth_outside_double") : nullptr;
        TH1D *h_fake_den_d = have_double ? (TH1D *)f_double->Get("h_reco_ET_matched_all_double") : nullptr;
        // Fallback naming without _double suffix
        if (!h_fake_num_d && have_double) h_fake_num_d = (TH1D *)f_double->Get("h_reco_ET_truth_outside");
        if (!h_fake_den_d && have_double) h_fake_den_d = (TH1D *)f_double->Get("h_reco_ET_matched_all");

        if (h_fake_num_s && h_fake_den_s)
        {
            TGraphAsymmErrors *g_fake_s = new TGraphAsymmErrors(h_fake_num_s, h_fake_den_s);
            g_fake_s->SetMarkerStyle(20);
            g_fake_s->SetMarkerSize(1.0);
            g_fake_s->SetMarkerColor(kBlack);
            g_fake_s->SetLineColor(kBlack);

            TGraphAsymmErrors *g_fake_d = nullptr;
            if (h_fake_num_d && h_fake_den_d)
            {
                g_fake_d = new TGraphAsymmErrors(h_fake_num_d, h_fake_den_d);
                g_fake_d->SetMarkerStyle(21);
                g_fake_d->SetMarkerSize(1.0);
                g_fake_d->SetMarkerColor(kRed);
                g_fake_d->SetLineColor(kRed);
            }

            TCanvas *c = new TCanvas("c_fake_rate", "", 800, 600);
            TH1F *frame = new TH1F("frame_fake_rate",
                ";Reco #it{E}_{T} [GeV];Inward Migration Fake Rate", 100, 8, 40);
            float ymax_fake = 0.10;
            if (g_fake_d) { // accommodate double interaction values (~25%)
                for (int ip = 0; ip < g_fake_d->GetN(); ip++) {
                    double xx, yy; g_fake_d->GetPoint(ip, xx, yy);
                    if (xx >= 8 && xx <= 36 && yy > ymax_fake) ymax_fake = yy;
                }
            }
            frame->GetYaxis()->SetRangeUser(0, ymax_fake * 1.3);
            frame->GetXaxis()->SetRangeUser(10, 36);
            frame->Draw("axis");

            g_fake_s->Draw("P same");
            if (g_fake_d) g_fake_d->Draw("P same");

            float leg_y = 0.90;
            myText(0.20, leg_y, 1, strleg1.c_str(), 0.04);
            myText(0.20, leg_y - 0.05, 1, strleg2.c_str(), 0.04);
            myText(0.20, leg_y - 0.10, 1, strSigMC.c_str(), 0.035);

            myMarkerLineText(0.55, 0.90, 1.0, kBlack, 20, kBlack, 1, "Single interaction", 0.04, true);
            if (g_fake_d) myMarkerLineText(0.55, 0.85, 1.0, kRed, 21, kRed, 1, "Double interaction", 0.04, true);

            c->SaveAs("figures/eta_migration_fake_rate.pdf");
            delete c;
        }
        else
        {
            std::cout << "Missing fake rate histograms in single file" << std::endl;
        }
    }

    // =================================================================
    // Plot 4: Outward Migration Loss Rate vs pT (single vs double)
    //   loss rate = reco outside fiducial / all truth inside fiducial
    // =================================================================
    {
        TH1D *h_loss_num_s = (TH1D *)f_single->Get("h_truth_pT_reco_outside");
        TH1D *h_loss_den_s = (TH1D *)f_single->Get("h_truth_pT_inside_all");

        TH1D *h_loss_num_d = have_double ? (TH1D *)f_double->Get("h_truth_pT_reco_outside_double") : nullptr;
        TH1D *h_loss_den_d = have_double ? (TH1D *)f_double->Get("h_truth_pT_inside_all_double") : nullptr;
        if (!h_loss_num_d && have_double) h_loss_num_d = (TH1D *)f_double->Get("h_truth_pT_reco_outside");
        if (!h_loss_den_d && have_double) h_loss_den_d = (TH1D *)f_double->Get("h_truth_pT_inside_all");

        if (h_loss_num_s && h_loss_den_s)
        {
            TGraphAsymmErrors *g_loss_s = new TGraphAsymmErrors(h_loss_num_s, h_loss_den_s);
            g_loss_s->SetMarkerStyle(20);
            g_loss_s->SetMarkerSize(1.0);
            g_loss_s->SetMarkerColor(kBlack);
            g_loss_s->SetLineColor(kBlack);

            TGraphAsymmErrors *g_loss_d = nullptr;
            if (h_loss_num_d && h_loss_den_d)
            {
                g_loss_d = new TGraphAsymmErrors(h_loss_num_d, h_loss_den_d);
                g_loss_d->SetMarkerStyle(21);
                g_loss_d->SetMarkerSize(1.0);
                g_loss_d->SetMarkerColor(kRed);
                g_loss_d->SetLineColor(kRed);
            }

            TCanvas *c = new TCanvas("c_loss_rate", "", 800, 600);
            TH1F *frame = new TH1F("frame_loss_rate",
                ";Truth #it{p}_{T} [GeV];Outward Migration Loss Rate", 100, 8, 40);
            float ymax_loss = 0.05;
            if (g_loss_d) {
                for (int ip = 0; ip < g_loss_d->GetN(); ip++) {
                    double xx, yy; g_loss_d->GetPoint(ip, xx, yy);
                    if (xx >= 8 && xx <= 36 && yy > ymax_loss) ymax_loss = yy;
                }
            }
            frame->GetYaxis()->SetRangeUser(0, ymax_loss * 1.3);
            frame->GetXaxis()->SetRangeUser(10, 36);
            frame->Draw("axis");

            g_loss_s->Draw("P same");
            if (g_loss_d) g_loss_d->Draw("P same");

            myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.80, 1, strSigMC.c_str(), 0.035);

            myMarkerLineText(0.55, 0.90, 1.0, kBlack, 20, kBlack, 1, "Single interaction", 0.04, true);
            if (g_loss_d) myMarkerLineText(0.55, 0.85, 1.0, kRed, 21, kRed, 1, "Double interaction", 0.04, true);

            c->SaveAs("figures/eta_migration_loss_rate.pdf");
            delete c;
        }
        else
        {
            std::cout << "Missing loss rate histograms in single file" << std::endl;
        }
    }

    // =================================================================
    // Plot 5: Response Contamination Fraction vs Reco ET
    //   fraction = contaminated / (contaminated + clean)
    // =================================================================
    {
        TH2D *h2_contam_s = (TH2D *)f_single->Get("h2_response_contaminated");
        TH2D *h2_clean_s  = (TH2D *)f_single->Get("h2_response_clean");

        TH2D *h2_contam_d = have_double ? (TH2D *)f_double->Get("h2_response_contaminated_double") : nullptr;
        TH2D *h2_clean_d  = have_double ? (TH2D *)f_double->Get("h2_response_clean_double") : nullptr;
        if (!h2_contam_d && have_double) h2_contam_d = (TH2D *)f_double->Get("h2_response_contaminated");
        if (!h2_clean_d && have_double)  h2_clean_d  = (TH2D *)f_double->Get("h2_response_clean");

        if (h2_contam_s && h2_clean_s)
        {
            // Project onto reco ET axis (X-axis)
            TH1D *h_contam_s = h2_contam_s->ProjectionX("h_contam_proj_s");
            TH1D *h_clean_s  = h2_clean_s->ProjectionX("h_clean_proj_s");
            TH1D *h_total_s  = (TH1D *)h_contam_s->Clone("h_total_proj_s");
            h_total_s->Add(h_clean_s);

            TGraphAsymmErrors *g_frac_s = new TGraphAsymmErrors(h_contam_s, h_total_s);
            g_frac_s->SetMarkerStyle(20);
            g_frac_s->SetMarkerSize(1.0);
            g_frac_s->SetMarkerColor(kBlack);
            g_frac_s->SetLineColor(kBlack);

            TGraphAsymmErrors *g_frac_d = nullptr;
            if (h2_contam_d && h2_clean_d)
            {
                TH1D *h_contam_d = h2_contam_d->ProjectionX("h_contam_proj_d");
                TH1D *h_clean_d  = h2_clean_d->ProjectionX("h_clean_proj_d");
                TH1D *h_total_d  = (TH1D *)h_contam_d->Clone("h_total_proj_d");
                h_total_d->Add(h_clean_d);

                g_frac_d = new TGraphAsymmErrors(h_contam_d, h_total_d);
                g_frac_d->SetMarkerStyle(21);
                g_frac_d->SetMarkerSize(1.0);
                g_frac_d->SetMarkerColor(kRed);
                g_frac_d->SetLineColor(kRed);
            }

            TCanvas *c = new TCanvas("c_response_contam", "", 800, 600);
            TH1F *frame = new TH1F("frame_resp_contam",
                ";Reco #it{E}_{T} [GeV];Response Contamination Fraction", 100, 8, 40);
            float ymax_resp = 0.10;
            if (g_frac_d) {
                for (int ip = 0; ip < g_frac_d->GetN(); ip++) {
                    double xx, yy; g_frac_d->GetPoint(ip, xx, yy);
                    if (xx >= 8 && xx <= 36 && yy > ymax_resp) ymax_resp = yy;
                }
            }
            frame->GetYaxis()->SetRangeUser(0, ymax_resp * 1.3);
            frame->GetXaxis()->SetRangeUser(10, 36);
            frame->Draw("axis");

            g_frac_s->Draw("P same");
            if (g_frac_d) g_frac_d->Draw("P same");

            myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.80, 1, strSigMC.c_str(), 0.035);

            myMarkerLineText(0.55, 0.90, 1.0, kBlack, 20, kBlack, 1, "Single interaction", 0.04, true);
            if (g_frac_d) myMarkerLineText(0.55, 0.85, 1.0, kRed, 21, kRed, 1, "Double interaction", 0.04, true);

            c->SaveAs("figures/eta_migration_response_contamination.pdf");
            delete c;
        }
        else
        {
            std::cout << "Missing response contamination histograms in single file" << std::endl;
        }
    }

    // =================================================================
    // Plot 5b: Tight+Iso Response Matrix Contamination (inward + outward)
    //   Inward: fake_tight_iso / (good_tight_iso + fake_tight_iso) vs reco ET
    //   Outward: loss / (good + loss) vs truth pT — all clusters, not tight+iso specific
    // =================================================================
    {
        // --- Inward: tight+iso contamination fraction ---
        TH1D *h_good_ti_s = (TH1D *)f_single->Get("h_good_tight_iso");
        TH1D *h_fake_ti_s = (TH1D *)f_single->Get("h_reco_ET_truth_outside_tight_iso");

        TH1D *h_good_ti_d = have_double ? (TH1D *)f_double->Get("h_good_tight_iso_double") : nullptr;
        TH1D *h_fake_ti_d = have_double ? (TH1D *)f_double->Get("h_reco_ET_truth_outside_tight_iso_double") : nullptr;

        // --- Outward: loss rate ---
        TH1D *h_loss_num_s = (TH1D *)f_single->Get("h_truth_pT_reco_outside");
        TH1D *h_loss_den_s = (TH1D *)f_single->Get("h_truth_pT_inside_all");
        TH1D *h_loss_num_d = have_double ? (TH1D *)f_double->Get("h_truth_pT_reco_outside_double") : nullptr;
        TH1D *h_loss_den_d = have_double ? (TH1D *)f_double->Get("h_truth_pT_inside_all_double") : nullptr;

        if (h_good_ti_s && h_fake_ti_s)
        {
            // Inward: denominator = good_ti + fake_ti
            TH1D *rGood_s = rebinToPtBins(h_good_ti_s, "rGood_ti_s");
            TH1D *rFake_s = rebinToPtBins(h_fake_ti_s, "rFake_ti_s");
            TH1D *rDen_s = (TH1D *)rGood_s->Clone("rDen_ti_s");
            rDen_s->Add(rFake_s);

            TH1D *frac_in_s = (TH1D *)rFake_s->Clone("frac_in_s");
            frac_in_s->Divide(rFake_s, rDen_s, 1, 1, "B");

            TH1D *frac_in_d = nullptr;
            if (h_good_ti_d && h_fake_ti_d)
            {
                TH1D *rGood_d = rebinToPtBins(h_good_ti_d, "rGood_ti_d");
                TH1D *rFake_d = rebinToPtBins(h_fake_ti_d, "rFake_ti_d");
                TH1D *rDen_d = (TH1D *)rGood_d->Clone("rDen_ti_d");
                rDen_d->Add(rFake_d);
                frac_in_d = (TH1D *)rFake_d->Clone("frac_in_d");
                frac_in_d->Divide(rFake_d, rDen_d, 1, 1, "B");
            }

            // Outward: rebin loss rate
            TH1D *frac_out_s = nullptr;
            TH1D *frac_out_d = nullptr;
            if (h_loss_num_s && h_loss_den_s)
            {
                TH1D *rLossN_s = rebinToPtBins(h_loss_num_s, "rLossN_s");
                TH1D *rLossD_s = rebinToPtBins(h_loss_den_s, "rLossD_s");
                frac_out_s = (TH1D *)rLossN_s->Clone("frac_out_s");
                frac_out_s->Divide(rLossN_s, rLossD_s, 1, 1, "B");
            }
            if (h_loss_num_d && h_loss_den_d)
            {
                TH1D *rLossN_d = rebinToPtBins(h_loss_num_d, "rLossN_d");
                TH1D *rLossD_d = rebinToPtBins(h_loss_den_d, "rLossD_d");
                frac_out_d = (TH1D *)rLossN_d->Clone("frac_out_d");
                frac_out_d->Divide(rLossN_d, rLossD_d, 1, 1, "B");
            }

            TCanvas *c = new TCanvas("c_tightiso_response", "", 800, 600);
            TH1F *frame = new TH1F("frame_tightiso_response",
                ";#it{E}_{T} (#it{p}_{T}) [GeV];Fraction", 100, 8, 36);
            // Auto y-range
            float ymax_ti = 0.05;
            if (frac_out_d) {
                for (int ib = 1; ib <= frac_out_d->GetNbinsX(); ib++)
                    if (frac_out_d->GetBinContent(ib) > ymax_ti) ymax_ti = frac_out_d->GetBinContent(ib);
            }
            if (frac_out_s) {
                for (int ib = 1; ib <= frac_out_s->GetNbinsX(); ib++)
                    if (frac_out_s->GetBinContent(ib) > ymax_ti) ymax_ti = frac_out_s->GetBinContent(ib);
            }
            frame->GetYaxis()->SetRangeUser(0, ymax_ti * 1.4);
            frame->Draw("axis");

            // Inward: filled circles
            frac_in_s->SetMarkerStyle(20);
            frac_in_s->SetMarkerSize(1.2);
            frac_in_s->SetMarkerColor(kBlack);
            frac_in_s->SetLineColor(kBlack);
            frac_in_s->Draw("E1 same");

            if (frac_in_d) {
                frac_in_d->SetMarkerStyle(20);
                frac_in_d->SetMarkerSize(1.2);
                frac_in_d->SetMarkerColor(kRed);
                frac_in_d->SetLineColor(kRed);
                frac_in_d->Draw("E1 same");
            }

            // Outward: open squares
            if (frac_out_s) {
                frac_out_s->SetMarkerStyle(24);
                frac_out_s->SetMarkerSize(1.2);
                frac_out_s->SetMarkerColor(kBlack);
                frac_out_s->SetLineColor(kBlack);
                frac_out_s->Draw("E1 same");
            }
            if (frac_out_d) {
                frac_out_d->SetMarkerStyle(24);
                frac_out_d->SetMarkerSize(1.2);
                frac_out_d->SetMarkerColor(kRed);
                frac_out_d->SetLineColor(kRed);
                frac_out_d->Draw("E1 same");
            }

            myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.80, 1, strSigMC.c_str(), 0.035);

            TLegend *leg = new TLegend(0.40, 0.62, 0.93, 0.92);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.033);
            leg->AddEntry(frac_in_s, "Inward (tight+iso contam.) - single", "lp");
            if (frac_out_s) leg->AddEntry(frac_out_s, "Outward (loss rate) - single", "lp");
            if (frac_in_d) leg->AddEntry(frac_in_d, "Inward (tight+iso contam.) - double", "lp");
            if (frac_out_d) leg->AddEntry(frac_out_d, "Outward (loss rate) - double", "lp");
            leg->Draw();

            c->SaveAs("figures/eta_migration_tightiso_response.pdf");
            std::cout << "[Plot 5b] Saved figures/eta_migration_tightiso_response.pdf" << std::endl;
            delete c;
        }
        else
        {
            std::cout << "Missing tight+iso histograms for response contamination plot" << std::endl;
        }
    }

    // =================================================================
    // Plot 6: ABCD Distribution of Fakes (single interaction, stacked)
    // =================================================================
    {
        TH1D *h_A = (TH1D *)f_single->Get("h_reco_ET_truth_outside_tight_iso");
        TH1D *h_B = (TH1D *)f_single->Get("h_reco_ET_truth_outside_tight_noniso");
        TH1D *h_C = (TH1D *)f_single->Get("h_reco_ET_truth_outside_nontight_iso");
        TH1D *h_D = (TH1D *)f_single->Get("h_reco_ET_truth_outside_nontight_noniso");

        if (h_A && h_B && h_C && h_D)
        {
            h_A->SetFillColor(kRed);
            h_A->SetLineColor(kRed);
            h_B->SetFillColor(kBlue);
            h_B->SetLineColor(kBlue);
            h_C->SetFillColor(kGreen + 2);
            h_C->SetLineColor(kGreen + 2);
            h_D->SetFillColor(kGray + 1);
            h_D->SetLineColor(kGray + 1);

            THStack *hs = new THStack("hs_abcd_fakes", "");
            hs->Add(h_D);
            hs->Add(h_C);
            hs->Add(h_B);
            hs->Add(h_A);

            TCanvas *c = new TCanvas("c_abcd_fakes", "", 800, 600);
            TH1F *frame = new TH1F("frame_abcd_fakes",
                ";Reco #it{E}_{T} [GeV];Inward Migration Fakes", 100, 8, 40);
            float ymax = hs->GetMaximum("nostack") * 1.5;
            if (ymax <= 0) ymax = 10;
            frame->GetYaxis()->SetRangeUser(0, ymax);
            frame->GetXaxis()->SetRangeUser(10, 36);
            frame->Draw("axis");

            hs->Draw("HIST same");

            myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.80, 1, strSigMC.c_str(), 0.035);
            myText(0.20, 0.75, 1, "Single Interaction", 0.035);

            TLegend *leg = new TLegend(0.55, 0.70, 0.90, 0.92);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.035);
            leg->AddEntry(h_A, "A: tight iso (signal region)", "f");
            leg->AddEntry(h_B, "B: tight noniso", "f");
            leg->AddEntry(h_C, "C: nontight iso", "f");
            leg->AddEntry(h_D, "D: nontight noniso", "f");
            leg->Draw();

            c->SaveAs("figures/eta_migration_abcd_fakes.pdf");
            delete c;
        }
        else
        {
            std::cout << "Missing ABCD fake histograms in single file" << std::endl;
        }
    }

    // =================================================================
    // Plot 6b: ABCD Distribution of Fakes (double interaction, stacked)
    // =================================================================
    if (have_double)
    {
        TH1D *h_A_d = (TH1D *)f_double->Get("h_reco_ET_truth_outside_tight_iso_double");
        TH1D *h_B_d = (TH1D *)f_double->Get("h_reco_ET_truth_outside_tight_noniso_double");
        TH1D *h_C_d = (TH1D *)f_double->Get("h_reco_ET_truth_outside_nontight_iso_double");
        TH1D *h_D_d = (TH1D *)f_double->Get("h_reco_ET_truth_outside_nontight_noniso_double");

        if (h_A_d && h_B_d && h_C_d && h_D_d)
        {
            TH1D *hAd = (TH1D *)h_A_d->Clone("hAd_stack");
            TH1D *hBd = (TH1D *)h_B_d->Clone("hBd_stack");
            TH1D *hCd = (TH1D *)h_C_d->Clone("hCd_stack");
            TH1D *hDd = (TH1D *)h_D_d->Clone("hDd_stack");

            hAd->SetFillColor(kRed);    hAd->SetLineColor(kRed);
            hBd->SetFillColor(kBlue);   hBd->SetLineColor(kBlue);
            hCd->SetFillColor(kGreen + 2); hCd->SetLineColor(kGreen + 2);
            hDd->SetFillColor(kGray + 1);  hDd->SetLineColor(kGray + 1);

            THStack *hs_d = new THStack("hs_abcd_fakes_double", "");
            hs_d->Add(hDd);
            hs_d->Add(hCd);
            hs_d->Add(hBd);
            hs_d->Add(hAd);

            TCanvas *c = new TCanvas("c_abcd_fakes_double", "", 800, 600);
            TH1F *frame = new TH1F("frame_abcd_fakes_d",
                ";Reco #it{E}_{T} [GeV];Inward Migration Fakes (Double Interaction)", 100, 8, 40);
            float ymax = hs_d->GetMaximum("nostack") * 1.5;
            if (ymax <= 0) ymax = 10;
            frame->GetYaxis()->SetRangeUser(0, ymax);
            frame->GetXaxis()->SetRangeUser(10, 36);
            frame->Draw("axis");

            hs_d->Draw("HIST same");

            myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.80, 1, strSigMC.c_str(), 0.035);
            myText(0.20, 0.75, 1, "Double Interaction", 0.035);

            TLegend *leg = new TLegend(0.55, 0.70, 0.90, 0.92);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.035);
            leg->AddEntry(hAd, "A: tight iso (signal region)", "f");
            leg->AddEntry(hBd, "B: tight noniso", "f");
            leg->AddEntry(hCd, "C: nontight iso", "f");
            leg->AddEntry(hDd, "D: nontight noniso", "f");
            leg->Draw();

            c->SaveAs("figures/eta_migration_abcd_fakes_double.pdf");
            std::cout << "[Plot 6b] Saved figures/eta_migration_abcd_fakes_double.pdf" << std::endl;
            delete c;
        }
        else
        {
            std::cout << "Missing ABCD fake histograms in double file" << std::endl;
        }
    }

    // =================================================================
    // Plot 7: Truth Eta of Fakes (single vs double)
    // =================================================================
    {
        TH1D *h_eta_s = (TH1D *)f_single->Get("h_truth_eta_of_fakes");
        TH1D *h_eta_d = have_double ? (TH1D *)f_double->Get("h_truth_eta_of_fakes_double") : nullptr;
        if (!h_eta_d && have_double) h_eta_d = (TH1D *)f_double->Get("h_truth_eta_of_fakes");

        if (h_eta_s)
        {
            h_eta_s->SetLineColor(kBlack);
            h_eta_s->SetLineWidth(2);
            h_eta_s->SetStats(0);
            h_eta_s->SetTitle("");

            if (h_eta_d)
            {
                h_eta_d->SetLineColor(kRed);
                h_eta_d->SetLineWidth(2);
                h_eta_d->SetStats(0);
            }

            TCanvas *c = new TCanvas("c_truth_eta_fakes", "", 800, 600);
            float ymax_eta = h_eta_s->GetMaximum();
            if (h_eta_d && h_eta_d->GetMaximum() > ymax_eta) ymax_eta = h_eta_d->GetMaximum();
            ymax_eta *= 1.5;

            TH1F *frame = new TH1F("frame_truth_eta_fakes",
                ";Truth #eta of Inward Fakes;Counts", 100, -1.5, 1.5);
            frame->GetYaxis()->SetRangeUser(0, ymax_eta);
            frame->Draw("axis");

            h_eta_s->Draw("HIST same");
            if (h_eta_d) h_eta_d->Draw("HIST same");

            drawEtaLines1D(0, ymax_eta);

            myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.80, 1, strSigMC.c_str(), 0.035);

            myMarkerLineText(0.55, 0.90, 0, kBlack, 0, kBlack, 1, "Single interaction", 0.04, true);
            if (h_eta_d) myMarkerLineText(0.55, 0.85, 0, kRed, 0, kRed, 1, "Double interaction", 0.04, true);

            c->SaveAs("figures/eta_migration_truth_eta_fakes.pdf");
            delete c;
        }
        else
        {
            std::cout << "Missing h_truth_eta_of_fakes in single file" << std::endl;
        }
    }

    // =================================================================
    // Plot 8: Delta-eta vs ET (2D with profile overlay)
    // =================================================================
    {
        TH2D *h2_deta = (TH2D *)f_single->Get("h2_deta_vs_ET");
        if (h2_deta)
        {
            TCanvas *c = new TCanvas("c_deta_vs_ET", "", 900, 800);
            c->SetRightMargin(0.15);

            h2_deta->SetStats(0);
            h2_deta->SetTitle("");
            h2_deta->GetXaxis()->SetTitle("Reco #it{E}_{T} [GeV]");
            h2_deta->GetYaxis()->SetTitle("#Delta#eta = #eta_{reco} - #eta_{truth}");
            h2_deta->GetXaxis()->SetRangeUser(10, 36);
            h2_deta->Draw("COLZ");

            // Profile overlay: mean delta-eta vs ET
            TProfile *prof = h2_deta->ProfileX("prof_deta_vs_ET");
            prof->SetMarkerStyle(20);
            prof->SetMarkerSize(0.8);
            prof->SetMarkerColor(kBlack);
            prof->SetLineColor(kBlack);
            prof->SetLineWidth(2);
            prof->Draw("same");

            myText(0.20, 0.92, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.87, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.82, 1, "Single Interaction", 0.035);

            c->SaveAs("figures/eta_migration_deta_vs_ET.pdf");
            delete c;
        }
        else
        {
            std::cout << "Missing h2_deta_vs_ET in single file" << std::endl;
        }
    }

    // =================================================================
    // Plot 9: Signal leakage fraction — single vs double
    // Of all inward-migrating fakes, what fraction leaks into each ABCD region?
    // =================================================================
    {
        // Single interaction ABCD fakes
        TH1D *h_total_s = (TH1D *)f_single->Get("h_reco_ET_truth_outside");
        TH1D *h_A_s = (TH1D *)f_single->Get("h_reco_ET_truth_outside_tight_iso");
        TH1D *h_B_s = (TH1D *)f_single->Get("h_reco_ET_truth_outside_tight_noniso");
        TH1D *h_C_s = (TH1D *)f_single->Get("h_reco_ET_truth_outside_nontight_iso");
        TH1D *h_D_s = (TH1D *)f_single->Get("h_reco_ET_truth_outside_nontight_noniso");

        // Double interaction ABCD fakes
        TH1D *h_total_d = have_double ? (TH1D *)f_double->Get("h_reco_ET_truth_outside_double") : nullptr;
        TH1D *h_A_d = have_double ? (TH1D *)f_double->Get("h_reco_ET_truth_outside_tight_iso_double") : nullptr;
        TH1D *h_B_d = have_double ? (TH1D *)f_double->Get("h_reco_ET_truth_outside_tight_noniso_double") : nullptr;
        TH1D *h_C_d = have_double ? (TH1D *)f_double->Get("h_reco_ET_truth_outside_nontight_iso_double") : nullptr;
        TH1D *h_D_d = have_double ? (TH1D *)f_double->Get("h_reco_ET_truth_outside_nontight_noniso_double") : nullptr;

        if (h_total_s && h_A_s && h_C_s)
        {
            // Denominator: fakes passing any ABCD region (A+B+C+D), not all fakes
            // This avoids the ~98% of double-interaction fakes that fail common cuts
            TH1D *rA_s = rebinToPtBins(h_A_s, "rA_s");
            TH1D *rB_s = rebinToPtBins(h_B_s, "rB_s");
            TH1D *rC_s = rebinToPtBins(h_C_s, "rC_s");
            TH1D *rD_s = rebinToPtBins(h_D_s, "rD_s");
            TH1D *rABCD_s = (TH1D *)rA_s->Clone("rABCD_s");
            rABCD_s->Add(rB_s);
            rABCD_s->Add(rC_s);
            rABCD_s->Add(rD_s);

            // Compute signal leakage fraction: A / (A+B+C+D) per bin
            TH1D *frac_A_s = (TH1D *)rA_s->Clone("frac_A_s");
            frac_A_s->Divide(rA_s, rABCD_s, 1, 1, "B"); // binomial errors
            TH1D *frac_C_s = (TH1D *)rC_s->Clone("frac_C_s");
            frac_C_s->Divide(rC_s, rABCD_s, 1, 1, "B");

            TH1D *frac_A_d = nullptr;
            TH1D *frac_C_d = nullptr;
            if (h_total_d && h_A_d && h_B_d && h_C_d && h_D_d)
            {
                TH1D *rA_d = rebinToPtBins(h_A_d, "rA_d");
                TH1D *rB_d = rebinToPtBins(h_B_d, "rB_d");
                TH1D *rC_d = rebinToPtBins(h_C_d, "rC_d");
                TH1D *rD_d = rebinToPtBins(h_D_d, "rD_d");
                TH1D *rABCD_d = (TH1D *)rA_d->Clone("rABCD_d");
                rABCD_d->Add(rB_d);
                rABCD_d->Add(rC_d);
                rABCD_d->Add(rD_d);

                frac_A_d = (TH1D *)rA_d->Clone("frac_A_d");
                frac_A_d->Divide(rA_d, rABCD_d, 1, 1, "B");
                frac_C_d = (TH1D *)rC_d->Clone("frac_C_d");
                frac_C_d->Divide(rC_d, rABCD_d, 1, 1, "B");
            }

            // --- Plot 9a: Region A (tight+iso) leakage fraction ---
            TCanvas *c9a = new TCanvas("c_signal_leakage", "", 800, 600);
            TH1F *frame9a = new TH1F("frame_signal_leakage",
                ";Reco #it{E}_{T} [GeV];Fraction of ABCD-classified fakes", 100, 8, 36);
            frame9a->GetYaxis()->SetRangeUser(0, 1.0);
            frame9a->Draw("axis");

            // Region A: single = black circle, double = red circle
            frac_A_s->SetMarkerStyle(20);
            frac_A_s->SetMarkerSize(1.2);
            frac_A_s->SetMarkerColor(kBlack);
            frac_A_s->SetLineColor(kBlack);
            frac_A_s->Draw("E1 same");

            // Region C: single = black open square, double = red open square
            frac_C_s->SetMarkerStyle(24);
            frac_C_s->SetMarkerSize(1.2);
            frac_C_s->SetMarkerColor(kBlack);
            frac_C_s->SetLineColor(kBlack);
            frac_C_s->Draw("E1 same");

            if (frac_A_d)
            {
                frac_A_d->SetMarkerStyle(20);
                frac_A_d->SetMarkerSize(1.2);
                frac_A_d->SetMarkerColor(kRed);
                frac_A_d->SetLineColor(kRed);
                frac_A_d->Draw("E1 same");
            }
            if (frac_C_d)
            {
                frac_C_d->SetMarkerStyle(24);
                frac_C_d->SetMarkerSize(1.2);
                frac_C_d->SetMarkerColor(kRed);
                frac_C_d->SetLineColor(kRed);
                frac_C_d->Draw("E1 same");
            }

            myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.80, 1, strSigMC.c_str(), 0.035);

            TLegend *leg9 = new TLegend(0.45, 0.68, 0.92, 0.92);
            leg9->SetBorderSize(0);
            leg9->SetFillStyle(0);
            leg9->SetTextSize(0.033);
            leg9->AddEntry(frac_A_s, "Region A (tight+iso) - single", "lp");
            leg9->AddEntry(frac_C_s, "Region C (nontight+iso) - single", "lp");
            if (frac_A_d)
                leg9->AddEntry(frac_A_d, "Region A (tight+iso) - double", "lp");
            if (frac_C_d)
                leg9->AddEntry(frac_C_d, "Region C (nontight+iso) - double", "lp");
            leg9->Draw();

            c9a->SaveAs("figures/eta_migration_signal_leakage.pdf");
            std::cout << "[Plot 9] Saved figures/eta_migration_signal_leakage.pdf" << std::endl;
            delete c9a;
        }
        else
        {
            std::cout << "Missing ABCD fake histograms for signal leakage plot" << std::endl;
        }
    }

    // Leakage fractions cX = good_signal_in_X / good_signal_in_A
    // (matches the ABCD definition used by the main pipeline).
    {
        TH1D *h_gA_s = (TH1D *)f_single->Get("h_good_tight_iso");
        TH1D *h_gB_s = (TH1D *)f_single->Get("h_good_tight_noniso");
        TH1D *h_gC_s = (TH1D *)f_single->Get("h_good_nontight_iso");
        TH1D *h_gD_s = (TH1D *)f_single->Get("h_good_nontight_noniso");

        TH1D *h_gA_d = have_double ? (TH1D *)f_double->Get("h_good_tight_iso_double") : nullptr;
        TH1D *h_gB_d = have_double ? (TH1D *)f_double->Get("h_good_tight_noniso_double") : nullptr;
        TH1D *h_gC_d = have_double ? (TH1D *)f_double->Get("h_good_nontight_iso_double") : nullptr;
        TH1D *h_gD_d = have_double ? (TH1D *)f_double->Get("h_good_nontight_noniso_double") : nullptr;

        if (h_gA_s && h_gB_s && h_gC_s && h_gD_s)
        {
            // Single interaction leakage fractions
            TH1D *rA_s = rebinToPtBins(h_gA_s, "leak_rA_s");
            TH1D *rB_s = rebinToPtBins(h_gB_s, "leak_rB_s");
            TH1D *rC_s = rebinToPtBins(h_gC_s, "leak_rC_s");
            TH1D *rD_s = rebinToPtBins(h_gD_s, "leak_rD_s");

            TH1D *cB_s = (TH1D *)rB_s->Clone("cB_s");
            cB_s->Divide(rB_s, rA_s, 1, 1, "B");
            TH1D *cC_s = (TH1D *)rC_s->Clone("cC_s");
            cC_s->Divide(rC_s, rA_s, 1, 1, "B");
            TH1D *cD_s = (TH1D *)rD_s->Clone("cD_s");
            cD_s->Divide(rD_s, rA_s, 1, 1, "B");

            // Double interaction leakage fractions
            TH1D *cB_d = nullptr, *cC_d = nullptr, *cD_d = nullptr;
            if (h_gA_d && h_gB_d && h_gC_d && h_gD_d)
            {
                TH1D *rA_d = rebinToPtBins(h_gA_d, "leak_rA_d");
                TH1D *rB_d = rebinToPtBins(h_gB_d, "leak_rB_d");
                TH1D *rC_d = rebinToPtBins(h_gC_d, "leak_rC_d");
                TH1D *rD_d = rebinToPtBins(h_gD_d, "leak_rD_d");

                cB_d = (TH1D *)rB_d->Clone("cB_d");
                cB_d->Divide(rB_d, rA_d, 1, 1, "B");
                cC_d = (TH1D *)rC_d->Clone("cC_d");
                cC_d->Divide(rC_d, rA_d, 1, 1, "B");
                cD_d = (TH1D *)rD_d->Clone("cD_d");
                cD_d->Divide(rD_d, rA_d, 1, 1, "B");
            }

            // --- Upper panel: absolute leakage fractions ---
            TCanvas *c10 = new TCanvas("c_leakage_fracs", "", 800, 800);
            TPad *pad1 = new TPad("pad1_leak", "", 0, 0.35, 1, 1.0);
            pad1->SetBottomMargin(0.02);
            pad1->Draw();
            TPad *pad2 = new TPad("pad2_leak", "", 0, 0.0, 1, 0.35);
            pad2->SetTopMargin(0.02);
            pad2->SetBottomMargin(0.30);
            pad2->Draw();

            pad1->cd();
            float ymax_leak = 0.20;
            if (cB_d)
            {
                for (int ib = 1; ib <= cB_d->GetNbinsX(); ib++)
                    if (cB_d->GetBinContent(ib) > ymax_leak) ymax_leak = cB_d->GetBinContent(ib);
            }
            TH1F *frame10 = new TH1F("frame_leak_fracs",
                ";; c_{X} = N_{X}^{sig} / N_{A}^{sig}", 100, 8, 36);
            frame10->GetYaxis()->SetRangeUser(0, ymax_leak * 1.4);
            frame10->GetYaxis()->SetTitleSize(0.055);
            frame10->GetYaxis()->SetLabelSize(0.05);
            frame10->GetXaxis()->SetLabelSize(0);
            frame10->Draw("axis");

            // Single: black, filled markers
            cB_s->SetMarkerStyle(20); cB_s->SetMarkerSize(1.2);
            cB_s->SetMarkerColor(kBlack); cB_s->SetLineColor(kBlack);
            cB_s->Draw("E1 same");

            cC_s->SetMarkerStyle(21); cC_s->SetMarkerSize(1.2);
            cC_s->SetMarkerColor(kBlack); cC_s->SetLineColor(kBlack);
            cC_s->Draw("E1 same");

            cD_s->SetMarkerStyle(22); cD_s->SetMarkerSize(1.3);
            cD_s->SetMarkerColor(kBlack); cD_s->SetLineColor(kBlack);
            cD_s->Draw("E1 same");

            // Double: red, open markers
            if (cB_d)
            {
                cB_d->SetMarkerStyle(24); cB_d->SetMarkerSize(1.2);
                cB_d->SetMarkerColor(kRed); cB_d->SetLineColor(kRed);
                cB_d->Draw("E1 same");
            }
            if (cC_d)
            {
                cC_d->SetMarkerStyle(25); cC_d->SetMarkerSize(1.2);
                cC_d->SetMarkerColor(kRed); cC_d->SetLineColor(kRed);
                cC_d->Draw("E1 same");
            }
            if (cD_d)
            {
                cD_d->SetMarkerStyle(26); cD_d->SetMarkerSize(1.3);
                cD_d->SetMarkerColor(kRed); cD_d->SetLineColor(kRed);
                cD_d->Draw("E1 same");
            }

            myText(0.20, 0.92, 1, strleg1.c_str(), 0.045);
            myText(0.20, 0.86, 1, strleg2.c_str(), 0.045);
            myText(0.20, 0.80, 1, strSigMC.c_str(), 0.04);

            TLegend *leg10 = new TLegend(0.50, 0.55, 0.93, 0.92);
            leg10->SetBorderSize(0);
            leg10->SetFillStyle(0);
            leg10->SetTextSize(0.037);
            leg10->SetHeader("Single (filled)  Double (open)");
            leg10->AddEntry(cB_s, "c_{B} (tight, non-iso) - single", "lp");
            leg10->AddEntry(cC_s, "c_{C} (non-tight, iso) - single", "lp");
            leg10->AddEntry(cD_s, "c_{D} (non-tight, non-iso) - single", "lp");
            if (cB_d) leg10->AddEntry(cB_d, "c_{B} - double", "lp");
            if (cC_d) leg10->AddEntry(cC_d, "c_{C} - double", "lp");
            if (cD_d) leg10->AddEntry(cD_d, "c_{D} - double", "lp");
            leg10->Draw();

            // --- Lower panel: double / single ratio ---
            pad2->cd();
            TH1F *frame10r = new TH1F("frame_leak_ratio",
                ";Reco #it{E}_{T} [GeV];Double / Single", 100, 8, 36);
            frame10r->GetYaxis()->SetRangeUser(0, 5.0);
            frame10r->GetYaxis()->SetTitleSize(0.09);
            frame10r->GetYaxis()->SetTitleOffset(0.55);
            frame10r->GetYaxis()->SetLabelSize(0.08);
            frame10r->GetXaxis()->SetTitleSize(0.10);
            frame10r->GetXaxis()->SetLabelSize(0.08);
            frame10r->GetYaxis()->SetNdivisions(505);
            frame10r->Draw("axis");

            TLine *unity = new TLine(8, 1, 36, 1);
            unity->SetLineStyle(7);
            unity->SetLineColor(kGray + 1);
            unity->Draw("same");

            auto makeRatio = [](TH1D *num, TH1D *den, const char *name,
                                int marker, double size, int color) -> TH1D * {
                if (!num) return nullptr;
                TH1D *r = (TH1D *)num->Clone(name);
                r->Divide(den);
                r->SetMarkerStyle(marker); r->SetMarkerSize(size);
                r->SetMarkerColor(color);  r->SetLineColor(color);
                r->Draw("E1 same");
                return r;
            };
            TH1D *rat_B = makeRatio(cB_d, cB_s, "rat_B", 20, 1.0, kBlue + 1);
            TH1D *rat_C = makeRatio(cC_d, cC_s, "rat_C", 21, 1.0, kGreen + 2);
            TH1D *rat_D = makeRatio(cD_d, cD_s, "rat_D", 22, 1.1, kMagenta + 1);

            TLegend *leg10r = new TLegend(0.18, 0.65, 0.55, 0.95);
            leg10r->SetBorderSize(0);
            leg10r->SetFillStyle(0);
            leg10r->SetTextSize(0.07);
            if (rat_B) leg10r->AddEntry(rat_B, "c_{B}", "lp");
            if (rat_C) leg10r->AddEntry(rat_C, "c_{C}", "lp");
            if (rat_D) leg10r->AddEntry(rat_D, "c_{D}", "lp");
            leg10r->Draw();

            c10->SaveAs("figures/eta_migration_leakage_fractions.pdf");
            std::cout << "[Plot 10] Saved figures/eta_migration_leakage_fractions.pdf" << std::endl;
            delete c10;
        }
        else
        {
            std::cout << "Missing good signal ABCD histograms for leakage fraction plot" << std::endl;
        }
    }

    // _noNPB histograms are required: the NPB model is trained only on
    // in-fiducial clusters, so its score is a sentinel for outside-fiducial
    // ones — leaving the NPB cut on would zero the entire outflow.
    if (have_bkg)
    {
        struct RegionHists {
            const char *key;      // "A", "B", "C", "D"
            const char *label;    // legend label
            TH1D *native = nullptr;
            TH1D *inflow = nullptr;
            TH1D *outflow = nullptr;
            TH1D *inflow_rate = nullptr;
            TH1D *outflow_rate = nullptr;
            // Double interaction (full GEANT jet12_double)
            TH1D *native_d = nullptr;
            TH1D *inflow_d = nullptr;
            TH1D *outflow_d = nullptr;
            TH1D *inflow_rate_d = nullptr;
            TH1D *outflow_rate_d = nullptr;
        };

        RegionHists regs[4] = {
            {"A", "A: tight + iso"},
            {"B", "B: tight + non-iso"},
            {"C", "C: non-tight + iso"},
            {"D", "D: non-tight + non-iso"}
        };

        auto loadRatesForFile = [&](TFile *f, RegionHists &r,
                                     TH1D *&rb_nat, TH1D *&rb_in, TH1D *&rb_out,
                                     TH1D *&in_rate, TH1D *&out_rate,
                                     const char *tag)
        {
            TH1D *h_nat  = (TH1D *)f->Get(Form("h_bkg_%s_native_noNPB",  r.key));
            TH1D *h_in   = (TH1D *)f->Get(Form("h_bkg_%s_inflow_noNPB",  r.key));
            TH1D *h_out  = (TH1D *)f->Get(Form("h_bkg_%s_outflow_noNPB", r.key));
            if (!h_nat || !h_in || !h_out) return false;
            rb_nat = rebinToPtBins(h_nat, Form("rb_nat_%s_%s", r.key, tag));
            rb_in  = rebinToPtBins(h_in,  Form("rb_in_%s_%s",  r.key, tag));
            rb_out = rebinToPtBins(h_out, Form("rb_out_%s_%s", r.key, tag));
            TH1D *den_in = (TH1D *)rb_nat->Clone(Form("den_in_%s_%s",  r.key, tag));
            den_in->Add(rb_in);
            TH1D *den_out = (TH1D *)rb_nat->Clone(Form("den_out_%s_%s", r.key, tag));
            den_out->Add(rb_out);
            in_rate = (TH1D *)rb_in->Clone(Form("in_rate_%s_%s", r.key, tag));
            in_rate->Divide(rb_in, den_in, 1, 1, "B");
            out_rate = (TH1D *)rb_out->Clone(Form("out_rate_%s_%s", r.key, tag));
            out_rate->Divide(rb_out, den_out, 1, 1, "B");
            return true;
        };

        for (auto &r : regs)
        {
            TH1D *rb_nat = nullptr, *rb_in = nullptr, *rb_out = nullptr;
            if (!loadRatesForFile(f_bkg, r, rb_nat, rb_in, rb_out, r.inflow_rate, r.outflow_rate, "s"))
            {
                std::cout << "  [Plot 11] missing _noNPB histograms for region " << r.key << std::endl;
                continue;
            }
            r.native = rb_nat;
            r.inflow = rb_in;
            r.outflow = rb_out;

            if (have_bkg_double)
            {
                TH1D *rb_nat_d = nullptr, *rb_in_d = nullptr, *rb_out_d = nullptr;
                if (loadRatesForFile(f_bkg_d, r, rb_nat_d, rb_in_d, rb_out_d, r.inflow_rate_d, r.outflow_rate_d, "d"))
                {
                    r.native_d = rb_nat_d;
                    r.inflow_d = rb_in_d;
                    r.outflow_d = rb_out_d;
                }
            }
        }

        // --- 2x2 panel figure: inflow + outflow rates per region ---
        TCanvas *c11 = new TCanvas("c_bkg_sideband_flow", "", 1100, 900);
        c11->Divide(2, 2);
        for (int ir = 0; ir < 4; ir++)
        {
            c11->cd(ir + 1);
            RegionHists &r = regs[ir];
            if (!r.inflow_rate || !r.outflow_rate) continue;

            // Auto y-range across all graphs
            float ymax_panel = 0.10;
            auto scan = [&](TH1D *h) {
                if (!h) return;
                for (int ib = 1; ib <= h->GetNbinsX(); ib++)
                {
                    float v = h->GetBinContent(ib);
                    if (std::isfinite(v) && v > ymax_panel) ymax_panel = v;
                }
            };
            scan(r.inflow_rate);
            scan(r.outflow_rate);
            scan(r.inflow_rate_d);
            scan(r.outflow_rate_d);
            ymax_panel *= 1.4;

            TH1F *frame = new TH1F(Form("frame_bkg_flow_%d", ir),
                ";Reco #it{E}_{T} [GeV];Flow fraction",
                100, 8, 36);
            frame->GetYaxis()->SetRangeUser(0, ymax_panel);
            frame->Draw("axis");

            // Single interaction: filled markers
            r.inflow_rate->SetMarkerStyle(20);
            r.inflow_rate->SetMarkerSize(1.2);
            r.inflow_rate->SetMarkerColor(kBlack);
            r.inflow_rate->SetLineColor(kBlack);
            r.inflow_rate->Draw("E1 same");

            r.outflow_rate->SetMarkerStyle(21);
            r.outflow_rate->SetMarkerSize(1.2);
            r.outflow_rate->SetMarkerColor(kBlue + 1);
            r.outflow_rate->SetLineColor(kBlue + 1);
            r.outflow_rate->Draw("E1 same");

            // Double interaction: open markers, red
            if (r.inflow_rate_d)
            {
                r.inflow_rate_d->SetMarkerStyle(24);
                r.inflow_rate_d->SetMarkerSize(1.2);
                r.inflow_rate_d->SetMarkerColor(kRed + 1);
                r.inflow_rate_d->SetLineColor(kRed + 1);
                r.inflow_rate_d->Draw("E1 same");
            }
            if (r.outflow_rate_d)
            {
                r.outflow_rate_d->SetMarkerStyle(25);
                r.outflow_rate_d->SetMarkerSize(1.2);
                r.outflow_rate_d->SetMarkerColor(kOrange + 7);
                r.outflow_rate_d->SetLineColor(kOrange + 7);
                r.outflow_rate_d->Draw("E1 same");
            }

            myText(0.20, 0.90, 1, r.label, 0.045);
            myText(0.20, 0.83, 1, "Jet MC (bkg only, no-NPB)", 0.033);

            if (ir == 0)
            {
                TLegend *leg = new TLegend(0.42, 0.62, 0.93, 0.90);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.036);
                leg->SetHeader("Single (filled)  Double (open)");
                leg->AddEntry(r.inflow_rate,  "Inflow - single",  "lp");
                leg->AddEntry(r.outflow_rate, "Outflow - single", "lp");
                if (r.inflow_rate_d)  leg->AddEntry(r.inflow_rate_d,  "Inflow - double",  "lp");
                if (r.outflow_rate_d) leg->AddEntry(r.outflow_rate_d, "Outflow - double", "lp");
                leg->Draw();
            }
        }
        c11->SaveAs("figures/eta_migration_bkg_sideband_flow.pdf");
        std::cout << "[Plot 11] Saved figures/eta_migration_bkg_sideband_flow.pdf" << std::endl;
        delete c11;

        // (inflow - outflow) / native: net bias in the sideband count.
        TCanvas *c11b = new TCanvas("c_bkg_net_flow", "", 900, 650);
        TH1F *frame11b = new TH1F("frame_bkg_net_flow",
            ";Reco #it{E}_{T} [GeV];(Inflow #minus Outflow) / Native", 100, 8, 36);

        float ymax_net =  0.05;
        float ymin_net = -0.05;
        TH1D *net_hists_s[4] = {nullptr, nullptr, nullptr, nullptr};
        TH1D *net_hists_d[4] = {nullptr, nullptr, nullptr, nullptr};
        Color_t net_colors[4] = {kBlack, kBlue + 1, kGreen + 2, kOrange + 7};
        Style_t net_markers_s[4] = {20, 21, 22, 23};
        Style_t net_markers_d[4] = {24, 25, 26, 32};

        auto computeNet = [&](TH1D *rb_nat, TH1D *rb_in, TH1D *rb_out, const char *name) -> TH1D * {
            if (!rb_nat || !rb_in || !rb_out) return nullptr;
            TH1D *net = (TH1D *)rb_in->Clone(name);
            net->Add(rb_out, -1.0);
            net->Divide(rb_nat);
            return net;
        };

        for (int ir = 0; ir < 4; ir++)
        {
            RegionHists &r = regs[ir];
            net_hists_s[ir] = computeNet(r.native,   r.inflow,   r.outflow,   Form("net_s_%s", r.key));
            net_hists_d[ir] = computeNet(r.native_d, r.inflow_d, r.outflow_d, Form("net_d_%s", r.key));
            for (TH1D *h : {net_hists_s[ir], net_hists_d[ir]})
            {
                if (!h) continue;
                for (int ib = 1; ib <= h->GetNbinsX(); ib++)
                {
                    float v = h->GetBinContent(ib);
                    if (std::isfinite(v))
                    {
                        if (v > ymax_net) ymax_net = v;
                        if (v < ymin_net) ymin_net = v;
                    }
                }
            }
        }
        float span = std::max(std::abs(ymax_net), std::abs(ymin_net)) * 1.4;
        frame11b->GetYaxis()->SetRangeUser(-span, span);
        frame11b->Draw("axis");
        TLine *zero11 = new TLine(8, 0, 36, 0);
        zero11->SetLineStyle(7);
        zero11->SetLineColor(kGray + 1);
        zero11->Draw("same");

        TLegend *leg11b = new TLegend(0.18, 0.65, 0.55, 0.92);
        leg11b->SetBorderSize(0);
        leg11b->SetFillStyle(0);
        leg11b->SetTextSize(0.034);
        leg11b->SetHeader("Single (filled)  Double (open)");

        for (int ir = 0; ir < 4; ir++)
        {
            if (net_hists_s[ir])
            {
                net_hists_s[ir]->SetMarkerStyle(net_markers_s[ir]);
                net_hists_s[ir]->SetMarkerSize(1.2);
                net_hists_s[ir]->SetMarkerColor(net_colors[ir]);
                net_hists_s[ir]->SetLineColor(net_colors[ir]);
                net_hists_s[ir]->Draw("E1 same");
                leg11b->AddEntry(net_hists_s[ir], Form("%s - single", regs[ir].label), "lp");
            }
            if (net_hists_d[ir])
            {
                net_hists_d[ir]->SetMarkerStyle(net_markers_d[ir]);
                net_hists_d[ir]->SetMarkerSize(1.2);
                net_hists_d[ir]->SetMarkerColor(net_colors[ir]);
                net_hists_d[ir]->SetLineColor(net_colors[ir]);
                net_hists_d[ir]->Draw("E1 same");
                leg11b->AddEntry(net_hists_d[ir], Form("%s - double", regs[ir].label), "lp");
            }
        }
        leg11b->Draw();

        myText(0.60, 0.92, 1, strleg1.c_str(), 0.040);
        myText(0.60, 0.87, 1, strleg2.c_str(), 0.040);
        myText(0.60, 0.82, 1, "Jet MC (bkg only, no-NPB)", 0.034);

        c11b->SaveAs("figures/eta_migration_bkg_net_flow.pdf");
        std::cout << "[Plot 11b] Saved figures/eta_migration_bkg_net_flow.pdf" << std::endl;
        delete c11b;

        // R-factor bias estimate: R_meas / R_native, with R = B*C/D.
        TH1D *rR_s = nullptr, *rR_d = nullptr;
        auto makeRfactorRatio = [&](TFile *f, const char *tag) -> TH1D * {
            TH1D *hB_n = (TH1D *)f->Get("h_bkg_B_native_noNPB");
            TH1D *hC_n = (TH1D *)f->Get("h_bkg_C_native_noNPB");
            TH1D *hD_n = (TH1D *)f->Get("h_bkg_D_native_noNPB");
            TH1D *hB_i = (TH1D *)f->Get("h_bkg_B_inflow_noNPB");
            TH1D *hC_i = (TH1D *)f->Get("h_bkg_C_inflow_noNPB");
            TH1D *hD_i = (TH1D *)f->Get("h_bkg_D_inflow_noNPB");
            TH1D *hB_o = (TH1D *)f->Get("h_bkg_B_outflow_noNPB");
            TH1D *hC_o = (TH1D *)f->Get("h_bkg_C_outflow_noNPB");
            TH1D *hD_o = (TH1D *)f->Get("h_bkg_D_outflow_noNPB");
            if (!hB_n || !hC_n || !hD_n || !hB_i || !hC_i || !hD_i || !hB_o || !hC_o || !hD_o)
                return nullptr;

            TH1D *rB_n = rebinToPtBins(hB_n, Form("rB_n_Rf_%s", tag));
            TH1D *rC_n = rebinToPtBins(hC_n, Form("rC_n_Rf_%s", tag));
            TH1D *rD_n = rebinToPtBins(hD_n, Form("rD_n_Rf_%s", tag));
            TH1D *rB_i = rebinToPtBins(hB_i, Form("rB_i_Rf_%s", tag));
            TH1D *rC_i = rebinToPtBins(hC_i, Form("rC_i_Rf_%s", tag));
            TH1D *rD_i = rebinToPtBins(hD_i, Form("rD_i_Rf_%s", tag));
            TH1D *rB_o = rebinToPtBins(hB_o, Form("rB_o_Rf_%s", tag));
            TH1D *rC_o = rebinToPtBins(hC_o, Form("rC_o_Rf_%s", tag));
            TH1D *rD_o = rebinToPtBins(hD_o, Form("rD_o_Rf_%s", tag));

            // B_true = B_native, i.e. what the sideband SHOULD see without migration.
            // B_meas = B_native + B_inflow, i.e. what the analysis actually sees
            // (the sideband only counts reco-inside clusters). Outflow is the
            // counterfactual — clusters that would have contributed if reco were in-fid.
            TH1D *B_meas = (TH1D *)rB_n->Clone(Form("B_meas_%s", tag));
            B_meas->Add(rB_i);
            TH1D *C_meas = (TH1D *)rC_n->Clone(Form("C_meas_%s", tag));
            C_meas->Add(rC_i);
            TH1D *D_meas = (TH1D *)rD_n->Clone(Form("D_meas_%s", tag));
            D_meas->Add(rD_i);

            TH1D *R_native = (TH1D *)rB_n->Clone(Form("R_native_%s", tag));
            R_native->Multiply(rC_n);
            R_native->Divide(rD_n);

            TH1D *R_meas = (TH1D *)B_meas->Clone(Form("R_meas_%s", tag));
            R_meas->Multiply(C_meas);
            R_meas->Divide(D_meas);

            TH1D *rR = (TH1D *)R_meas->Clone(Form("R_ratio_%s", tag));
            rR->Divide(R_native);
            return rR;
        };

        rR_s = makeRfactorRatio(f_bkg, "s");
        if (have_bkg_double)
            rR_d = makeRfactorRatio(f_bkg_d, "d");

        if (rR_s)
        {
            TCanvas *c11c = new TCanvas("c_bkg_rfactor_bias", "", 900, 650);
            TH1F *frame11c = new TH1F("frame_bkg_rfactor",
                ";Reco #it{E}_{T} [GeV];R_{measured} / R_{native}", 100, 8, 36);
            float ymax_r = 1.05, ymin_r = 0.95;
            auto scanR = [&](TH1D *h) {
                if (!h) return;
                for (int ib = 1; ib <= h->GetNbinsX(); ib++)
                {
                    float v = h->GetBinContent(ib);
                    if (std::isfinite(v) && v > 0)
                    {
                        if (v > ymax_r) ymax_r = v;
                        if (v < ymin_r) ymin_r = v;
                    }
                }
            };
            scanR(rR_s);
            scanR(rR_d);
            frame11c->GetYaxis()->SetRangeUser(std::max(0.0f, ymin_r - 0.05f), ymax_r + 0.05f);
            frame11c->Draw("axis");

            TLine *unity11c = new TLine(8, 1.0, 36, 1.0);
            unity11c->SetLineStyle(7);
            unity11c->SetLineColor(kGray + 1);
            unity11c->Draw("same");

            rR_s->SetMarkerStyle(20);
            rR_s->SetMarkerSize(1.3);
            rR_s->SetMarkerColor(kBlack);
            rR_s->SetLineColor(kBlack);
            rR_s->Draw("E1 same");

            if (rR_d)
            {
                rR_d->SetMarkerStyle(24);
                rR_d->SetMarkerSize(1.3);
                rR_d->SetMarkerColor(kRed + 1);
                rR_d->SetLineColor(kRed + 1);
                rR_d->Draw("E1 same");
            }

            myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.80, 1, "Jet MC: R = B #times C / D", 0.036);

            TLegend *leg11c = new TLegend(0.55, 0.72, 0.93, 0.92);
            leg11c->SetBorderSize(0);
            leg11c->SetFillStyle(0);
            leg11c->SetTextSize(0.037);
            leg11c->AddEntry(rR_s, "Single interaction (jet20+jet30)", "lp");
            if (rR_d) leg11c->AddEntry(rR_d, "Double interaction (jet12_double)", "lp");
            leg11c->Draw();

            c11c->SaveAs("figures/eta_migration_bkg_rfactor_bias.pdf");
            std::cout << "[Plot 11c] Saved figures/eta_migration_bkg_rfactor_bias.pdf" << std::endl;
            delete c11c;
        }

        // --- 2D eta truth vs reco for background (single + double side by side) ---
        {
            auto drawMigrationMatrix = [&](TFile *f, const char *tag, const char *title_text,
                                           const char *outname)
            {
                TH2D *h2_bkg = (TH2D *)f->Get("h2_bkg_eta_truth_vs_reco");
                if (!h2_bkg) return;
                TCanvas *cc = new TCanvas(Form("c_bkg_mig_matrix_%s", tag), "", 900, 800);
                cc->SetRightMargin(0.15);
                cc->SetLogz();

                h2_bkg->SetStats(0);
                h2_bkg->SetTitle("");
                h2_bkg->GetXaxis()->SetTitle("Reco #eta");
                h2_bkg->GetYaxis()->SetTitle("Truth #eta (leading contributor)");
                h2_bkg->GetXaxis()->SetRangeUser(-1.2, 1.2);
                h2_bkg->GetYaxis()->SetRangeUser(-1.2, 1.2);
                h2_bkg->Draw("COLZ");
                drawEtaLines2D();

                myText(0.20, 0.92, 1, strleg1.c_str(), 0.04);
                myText(0.20, 0.87, 1, strleg2.c_str(), 0.04);
                myText(0.20, 0.82, 1, title_text, 0.035);

                cc->SaveAs(outname);
                std::cout << "[Plot 11d] Saved " << outname << std::endl;
                delete cc;
            };

            drawMigrationMatrix(f_bkg, "single", "Jet MC (bkg only, single int.)",
                                "figures/eta_migration_bkg_matrix.pdf");
            if (have_bkg_double)
                drawMigrationMatrix(f_bkg_d, "double", "Jet MC (bkg only, jet12_double)",
                                    "figures/eta_migration_bkg_matrix_double.pdf");
        }

        // Print numbers per pT bin for report
        std::cout << std::endl;
        std::cout << "============================================================" << std::endl;
        std::cout << "  Background sideband flow (weighted, per pT bin)" << std::endl;
        std::cout << "============================================================" << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        for (auto &r : regs)
        {
            if (!r.inflow_rate || !r.outflow_rate) continue;
            std::cout << "  Region " << r.key << " (" << r.label << ")" << std::endl;
            std::cout << "    pT bin         in_s    out_s";
            if (r.inflow_rate_d) std::cout << "    in_d     out_d";
            std::cout << std::endl;
            for (int ib = 1; ib <= r.inflow_rate->GetNbinsX(); ib++)
            {
                float lo = r.inflow_rate->GetBinLowEdge(ib);
                float hi = lo + r.inflow_rate->GetBinWidth(ib);
                std::cout << "    [" << std::setw(4) << lo << "," << std::setw(4) << hi << "]  "
                          << std::setw(9) << r.inflow_rate->GetBinContent(ib) << "  "
                          << std::setw(9) << r.outflow_rate->GetBinContent(ib);
                if (r.inflow_rate_d)
                    std::cout << "  " << std::setw(9) << r.inflow_rate_d->GetBinContent(ib)
                              << "  " << std::setw(9) << r.outflow_rate_d->GetBinContent(ib);
                std::cout << std::endl;
            }
        }
        if (rR_s)
        {
            std::cout << "  R-factor bias (R_meas/R_native) per pT bin:" << std::endl;
            for (int ib = 1; ib <= rR_s->GetNbinsX(); ib++)
            {
                float lo = rR_s->GetBinLowEdge(ib);
                float hi = lo + rR_s->GetBinWidth(ib);
                std::cout << "    [" << std::setw(4) << lo << "," << std::setw(4) << hi << "]  "
                          << "single: " << std::setw(7) << rR_s->GetBinContent(ib);
                if (rR_d)
                    std::cout << "  double: " << std::setw(7) << rR_d->GetBinContent(ib);
                std::cout << std::endl;
            }
        }
    }

    // =================================================================
    // Print summary of key numbers
    // =================================================================
    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "  Eta Migration Study -- Summary" << std::endl;
    std::cout << "============================================================" << std::endl;

    // Overall fake rate (single)
    TH1D *h_fake_num = (TH1D *)f_single->Get("h_reco_ET_truth_outside");
    TH1D *h_fake_den = (TH1D *)f_single->Get("h_reco_ET_matched_all");
    if (h_fake_num && h_fake_den)
    {
        double total_fakes = h_fake_num->Integral();
        double total_matched = h_fake_den->Integral();
        double overall_fake_rate = (total_matched > 0) ? total_fakes / total_matched : 0;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "  Overall inward migration fake rate (single): "
                  << overall_fake_rate << " (" << total_fakes << " / " << total_matched << ")" << std::endl;

        // Per pT bin
        std::cout << "  Fake rate by pT bin (single):" << std::endl;
        for (int i = 0; i < NptBins; i++)
        {
            int bin_lo = h_fake_num->FindBin(ptRanges[i] + 0.01);
            int bin_hi = h_fake_num->FindBin(ptRanges[i + 1] - 0.01);
            double num = h_fake_num->Integral(bin_lo, bin_hi);
            double den = h_fake_den->Integral(bin_lo, bin_hi);
            double rate = (den > 0) ? num / den : 0;
            std::cout << "    pT [" << std::setw(4) << ptRanges[i] << ", "
                      << std::setw(4) << ptRanges[i + 1] << "] GeV:  "
                      << rate << std::endl;
        }
    }

    // Overall fake rate (double)
    if (have_double)
    {
        TH1D *h_fake_num_d = (TH1D *)f_double->Get("h_reco_ET_truth_outside_double");
        TH1D *h_fake_den_d = (TH1D *)f_double->Get("h_reco_ET_matched_all_double");
        if (!h_fake_num_d) h_fake_num_d = (TH1D *)f_double->Get("h_reco_ET_truth_outside");
        if (!h_fake_den_d) h_fake_den_d = (TH1D *)f_double->Get("h_reco_ET_matched_all");
        if (h_fake_num_d && h_fake_den_d)
        {
            double total_fakes_d = h_fake_num_d->Integral();
            double total_matched_d = h_fake_den_d->Integral();
            double overall_fake_rate_d = (total_matched_d > 0) ? total_fakes_d / total_matched_d : 0;
            std::cout << "  Overall inward migration fake rate (double): "
                      << overall_fake_rate_d << " (" << total_fakes_d << " / " << total_matched_d << ")" << std::endl;

            // Single vs double ratio
            if (h_fake_num && h_fake_den)
            {
                double rate_s = h_fake_num->Integral() / h_fake_den->Integral();
                double ratio = (rate_s > 0) ? overall_fake_rate_d / rate_s : 0;
                std::cout << "  Double / Single fake rate ratio: " << ratio << std::endl;
            }
        }
    }

    // Response contamination fraction (single)
    TH2D *h2_contam = (TH2D *)f_single->Get("h2_response_contaminated");
    TH2D *h2_clean  = (TH2D *)f_single->Get("h2_response_clean");
    if (h2_contam && h2_clean)
    {
        double n_contam = h2_contam->Integral();
        double n_clean  = h2_clean->Integral();
        double frac = (n_contam + n_clean > 0) ? n_contam / (n_contam + n_clean) : 0;
        std::cout << "  Response contamination fraction (single): "
                  << frac << " (" << n_contam << " contam / " << n_contam + n_clean << " total)" << std::endl;
    }

    // Response contamination (double)
    if (have_double)
    {
        TH2D *h2_contam_d = (TH2D *)f_double->Get("h2_response_contaminated_double");
        TH2D *h2_clean_d  = (TH2D *)f_double->Get("h2_response_clean_double");
        if (!h2_contam_d) h2_contam_d = (TH2D *)f_double->Get("h2_response_contaminated");
        if (!h2_clean_d)  h2_clean_d  = (TH2D *)f_double->Get("h2_response_clean");
        if (h2_contam_d && h2_clean_d)
        {
            double n_contam_d = h2_contam_d->Integral();
            double n_clean_d  = h2_clean_d->Integral();
            double frac_d = (n_contam_d + n_clean_d > 0) ? n_contam_d / (n_contam_d + n_clean_d) : 0;
            std::cout << "  Response contamination fraction (double): " << frac_d << std::endl;
        }
    }

    // Loss rate (single)
    TH1D *h_loss_num = (TH1D *)f_single->Get("h_truth_pT_reco_outside");
    TH1D *h_loss_den = (TH1D *)f_single->Get("h_truth_pT_inside_all");
    if (h_loss_num && h_loss_den)
    {
        double total_loss = h_loss_num->Integral();
        double total_inside = h_loss_den->Integral();
        double overall_loss = (total_inside > 0) ? total_loss / total_inside : 0;
        std::cout << "  Overall outward migration loss rate (single): " << overall_loss << std::endl;
    }

    std::cout << "============================================================" << std::endl;
    std::cout << "  Plots saved to figures/eta_migration_*.pdf" << std::endl;
    std::cout << "============================================================" << std::endl;
}
