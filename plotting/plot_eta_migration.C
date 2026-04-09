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
    const std::string &double_file = "../efficiencytool/results/eta_migration_double_signal_nom.root")
{
    init_plot();
    gSystem->Exec("mkdir -p figures");

    TFile *f_single = TFile::Open(single_file.c_str());
    TFile *f_double = TFile::Open(double_file.c_str());

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
            frame->GetYaxis()->SetRangeUser(0, 0.15);
            frame->GetXaxis()->SetRangeUser(8, 36);
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
            frame->GetYaxis()->SetRangeUser(0, 0.15);
            frame->GetXaxis()->SetRangeUser(8, 36);
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
            frame->GetYaxis()->SetRangeUser(0, 0.20);
            frame->GetXaxis()->SetRangeUser(8, 36);
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
            frame->GetXaxis()->SetRangeUser(8, 36);
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
            frame->GetXaxis()->SetRangeUser(8, 36);
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
            h2_deta->GetXaxis()->SetRangeUser(8, 36);
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
            // Rebin to analysis pT bins for cleaner ratios
            const int nPt = NptBins;
            double edges[NptBins + 1];
            for (int i = 0; i <= NptBins; i++) edges[i] = ptRanges[i];

            auto rebinTo = [&](TH1D *h, const char *name) -> TH1D * {
                return (TH1D *)h->Rebin(nPt, name, edges);
            };

            // Denominator: fakes passing any ABCD region (A+B+C+D), not all fakes
            // This avoids the ~98% of double-interaction fakes that fail common cuts
            TH1D *rA_s = rebinTo(h_A_s, "rA_s");
            TH1D *rB_s = rebinTo(h_B_s, "rB_s");
            TH1D *rC_s = rebinTo(h_C_s, "rC_s");
            TH1D *rD_s = rebinTo(h_D_s, "rD_s");
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
                TH1D *rA_d = rebinTo(h_A_d, "rA_d");
                TH1D *rB_d = rebinTo(h_B_d, "rB_d");
                TH1D *rC_d = rebinTo(h_C_d, "rC_d");
                TH1D *rD_d = rebinTo(h_D_d, "rD_d");
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
