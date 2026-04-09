#include "plotcommon.h"

void plot_deltaR_single_study()
{
    init_plot();

    const std::string infile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/deltaR_single_study.root";
    const std::string savePath = "/sphenix/user/shuhangli/ppg12/plotting/figures/deltaR_study";
    gSystem->mkdir(savePath.c_str(), true);

    TFile *fin = TFile::Open(infile.c_str());
    if (!fin || fin->IsZombie())
    {
        std::cerr << "Cannot open " << infile << std::endl;
        return;
    }

    // Standard label positions (from plot_efficiency.C pattern)
    float xpos = 0.20, ypos = 0.885, dy = 0.054, fontsize = 0.046;

    // ===================================================================
    // Canvas 1: deltaR distribution
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_deltaR", "", 600, 600);

        TH1 *h_photon5  = (TH1 *)fin->Get("h_deltaR_photon5");
        TH1 *h_photon10 = (TH1 *)fin->Get("h_deltaR_photon10");
        TH1 *h_photon20 = (TH1 *)fin->Get("h_deltaR_photon20");
        TH1 *h_combined  = (TH1 *)fin->Get("h_deltaR_combined");

        if (!h_combined)
        {
            std::cerr << "Missing h_deltaR_combined" << std::endl;
            return;
        }

        // Normalize to unit area
        auto normHist = [](TH1 *h)
        {
            if (h && h->Integral() > 0)
                h->Scale(1.0 / h->Integral());
        };
        normHist(h_photon5);
        normHist(h_photon10);
        normHist(h_photon20);
        normHist(h_combined);

        TH1F *frame = new TH1F("frame_dR", ";#Delta#it{R}(cluster, truth);Normalized", 100, 0, 0.3);
        float ymax = h_combined->GetMaximum() * 5;
        frame->GetYaxis()->SetRangeUser(1e-4, ymax);
        frame->Draw("axis");

        h_combined->SetLineColor(kBlack);
        h_combined->SetLineWidth(3);
        h_combined->Draw("hist same");

        if (h_photon5)
        {
            h_photon5->SetLineColor(kBlue);
            h_photon5->SetLineWidth(2);
            h_photon5->SetLineStyle(1);
            h_photon5->Draw("hist same");
        }
        if (h_photon10)
        {
            h_photon10->SetLineColor(kRed);
            h_photon10->SetLineWidth(2);
            h_photon10->SetLineStyle(1);
            h_photon10->Draw("hist same");
        }
        if (h_photon20)
        {
            h_photon20->SetLineColor(kGreen + 2);
            h_photon20->SetLineWidth(2);
            h_photon20->SetLineStyle(1);
            h_photon20->Draw("hist same");
        }

        // Vertical dashed line at deltaR = 0.1
        TLine *lcut = new TLine(0.1, 1e-4, 0.1, ymax);
        lcut->SetLineStyle(7);
        lcut->SetLineColor(kGray + 2);
        lcut->SetLineWidth(2);
        lcut->Draw("same");

        TLegend *leg = new TLegend(0.50, 0.60, 0.88, 0.82);
        legStyle(leg, 0.20, 0.040);
        leg->AddEntry(h_combined, "Combined (weighted)", "l");
        if (h_photon5)  leg->AddEntry(h_photon5, "photon5 (0-14 GeV)", "l");
        if (h_photon10) leg->AddEntry(h_photon10, "photon10 (14-30 GeV)", "l");
        if (h_photon20) leg->AddEntry(h_photon20, "photon20 (30+ GeV)", "l");
        leg->Draw("same");

        myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
        myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
        myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 0);
        myText(xpos, ypos - 3 * dy, 1, "Single-interaction MC", fontsize, 0);

        gPad->SetLogy();
        c->SaveAs(Form("%s/deltaR_distribution.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Canvas 2: deltaR vs dz (2D + profile)
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_dR_dz", "", 1600, 700);
        c->Divide(2, 1);

        TH2 *h2 = (TH2 *)fin->Get("h_deltaR_vs_dz_combined");
        if (!h2)
        {
            std::cerr << "Missing h_deltaR_vs_dz_combined" << std::endl;
        }
        else
        {
            // Left panel: 2D colz
            c->cd(1);
            gPad->SetRightMargin(0.15);
            h2->SetTitle(";#it{dz} = vtx_{reco} #minus vtx_{truth} [cm];#Delta#it{R}(cluster, truth)");
            h2->Draw("colz");

            TLine *lh = new TLine(-30, 0.1, 30, 0.1);
            lh->SetLineStyle(7);
            lh->SetLineColor(kGray + 2);
            lh->SetLineWidth(2);
            lh->Draw("same");

            // Analytic curve: deltaR ~ |dz| / R_CEMC at eta=0
            // delta_eta = |dz| / (R_CEMC * cosh(0)) = |dz| / 93.5
            const double R_CEMC = 93.5;
            TGraph *g_analytic = new TGraph();
            for (int i = 0; i <= 60; ++i)
            {
                double dz = -30.0 + i;
                double dR_est = fabs(dz) / R_CEMC;
                g_analytic->SetPoint(i, dz, dR_est);
            }
            g_analytic->SetLineColor(kRed);
            g_analytic->SetLineWidth(2);
            g_analytic->SetLineStyle(2);
            g_analytic->Draw("L same");

            // Right panel: ProfileX
            c->cd(2);
            TProfile *prof = h2->ProfileX("prof_dR_dz");
            prof->SetTitle(";#it{dz} = vtx_{reco} #minus vtx_{truth} [cm];#LT#Delta#it{R}#GT");
            prof->SetMarkerStyle(20);
            prof->SetMarkerSize(0.8);
            prof->SetLineColor(kBlack);
            prof->Draw("e");

            TLine *lh2 = new TLine(-30, 0.1, 30, 0.1);
            lh2->SetLineStyle(7);
            lh2->SetLineColor(kGray + 2);
            lh2->SetLineWidth(2);
            lh2->Draw("same");

            g_analytic->Draw("L same");

            TLegend *leg = new TLegend(0.35, 0.75, 0.88, 0.88);
            legStyle(leg, 0.20, 0.040);
            leg->AddEntry(prof, "Profile (mean #Delta#it{R})", "pl");
            leg->AddEntry(g_analytic, "|#it{dz}| / #it{R}_{CEMC} (#eta = 0)", "l");
            leg->Draw("same");
        }

        c->cd(2);
        myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
        myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
        myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 0);

        c->SaveAs(Form("%s/deltaR_vs_dz.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Canvas 3: deltaR vs pT (2D + profile)
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_dR_pT", "", 1600, 700);
        c->Divide(2, 1);

        TH2 *h2 = (TH2 *)fin->Get("h_deltaR_vs_pT_combined");
        if (!h2)
        {
            std::cerr << "Missing h_deltaR_vs_pT_combined" << std::endl;
        }
        else
        {
            c->cd(1);
            gPad->SetRightMargin(0.15);
            h2->SetTitle(";#it{p}_{T}^{truth} [GeV];#Delta#it{R}(cluster, truth)");
            h2->Draw("colz");

            TLine *lh = new TLine(8, 0.1, 36, 0.1);
            lh->SetLineStyle(7);
            lh->SetLineColor(kGray + 2);
            lh->SetLineWidth(2);
            lh->Draw("same");

            c->cd(2);
            TProfile *prof = h2->ProfileX("prof_dR_pT");
            prof->SetTitle(";#it{p}_{T}^{truth} [GeV];#LT#Delta#it{R}#GT");
            prof->SetMarkerStyle(20);
            prof->SetMarkerSize(0.8);
            prof->SetLineColor(kBlack);
            prof->Draw("e");

            TLine *lh2 = new TLine(8, 0.1, 36, 0.1);
            lh2->SetLineStyle(7);
            lh2->SetLineColor(kGray + 2);
            lh2->SetLineWidth(2);
            lh2->Draw("same");
        }

        c->cd(2);
        myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
        myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
        myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 0);

        c->SaveAs(Form("%s/deltaR_vs_pT.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Canvas 4: deltaR vs eta (2D + profile)
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_dR_eta", "", 1600, 700);
        c->Divide(2, 1);

        TH2 *h2 = (TH2 *)fin->Get("h_deltaR_vs_eta_combined");
        if (!h2)
        {
            std::cerr << "Missing h_deltaR_vs_eta_combined" << std::endl;
        }
        else
        {
            c->cd(1);
            gPad->SetRightMargin(0.15);
            h2->SetTitle(";#it{#eta}^{truth};#Delta#it{R}(cluster, truth)");
            h2->Draw("colz");

            TLine *lh = new TLine(-0.7, 0.1, 0.7, 0.1);
            lh->SetLineStyle(7);
            lh->SetLineColor(kGray + 2);
            lh->SetLineWidth(2);
            lh->Draw("same");

            c->cd(2);
            TProfile *prof = h2->ProfileX("prof_dR_eta");
            prof->SetTitle(";#it{#eta}^{truth};#LT#Delta#it{R}#GT");
            prof->SetMarkerStyle(20);
            prof->SetMarkerSize(0.8);
            prof->SetLineColor(kBlack);
            prof->Draw("e");

            TLine *lh2 = new TLine(-0.7, 0.1, 0.7, 0.1);
            lh2->SetLineStyle(7);
            lh2->SetLineColor(kGray + 2);
            lh2->SetLineWidth(2);
            lh2->Draw("same");
        }

        c->cd(2);
        myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
        myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
        myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 0);

        c->SaveAs(Form("%s/deltaR_vs_eta.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Canvas 5: Efficiency threshold scan (KEY PLOT)
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_eff_scan", "", 600, 600);

        TH1F *frame = new TH1F("frame_eff_scan", ";#it{E}_{T}^{#gamma, truth} [GeV];Reconstruction Efficiency", 100, 8, 36);
        frame->GetYaxis()->SetRangeUser(0.80, 1.05);
        frame->Draw("axis");

        struct ThreshInfo
        {
            const char *name;
            const char *label;
            int color;
            int marker;
        };

        ThreshInfo thresholds[] = {
            {"eff_reco_dR0p05_combined", "#Delta#it{R} < 0.05", kRed + 1, 20},
            {"eff_reco_dR0p08_combined", "#Delta#it{R} < 0.08", kOrange + 1, 21},
            {"eff_reco_dR0p1_combined", "#Delta#it{R} < 0.1 (nominal)", kBlack, 33},
            {"eff_reco_dR0p15_combined", "#Delta#it{R} < 0.15", kGreen + 2, 34},
            {"eff_reco_dR0p2_combined", "#Delta#it{R} < 0.2", kBlue, 24},
            {"eff_reco_dR0p3_combined", "#Delta#it{R} < 0.3", kViolet + 1, 25},
            {"eff_reco_noDR_combined", "trkID only (no #Delta#it{R})", kGray + 2, 28}};
        int nThresh = sizeof(thresholds) / sizeof(thresholds[0]);

        TLegend *leg = new TLegend(0.42, 0.18, 0.88, 0.52);
        legStyle(leg, 0.20, 0.038);

        for (int i = 0; i < nThresh; ++i)
        {
            TEfficiency *eff = (TEfficiency *)fin->Get(thresholds[i].name);
            if (!eff)
            {
                std::cerr << "Missing " << thresholds[i].name << std::endl;
                continue;
            }
            eff->SetMarkerColor(thresholds[i].color);
            eff->SetLineColor(thresholds[i].color);
            eff->SetMarkerStyle(thresholds[i].marker);
            eff->SetMarkerSize(1.2);
            eff->Draw("same");
            leg->AddEntry(eff, thresholds[i].label, "pl");
        }
        leg->Draw("same");

        myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
        myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
        myText(xpos, ypos - 2 * dy, 1, "Single-interaction MC", fontsize, 0);
        myText(xpos, ypos - 3 * dy, 1, strleg3.c_str(), fontsize, 0);

        c->SaveAs(Form("%s/efficiency_threshold_scan.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Canvas 6: ET ratio for pass vs fail
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_ET_ratio", "", 600, 600);

        TH1 *h_pass = (TH1 *)fin->Get("h_ET_ratio_pass_combined");
        TH1 *h_fail = (TH1 *)fin->Get("h_ET_ratio_fail_combined");

        if (!h_pass || !h_fail)
        {
            std::cerr << "Missing h_ET_ratio_pass/fail_combined" << std::endl;
        }
        else
        {
            // Normalize to unit area
            if (h_pass->Integral() > 0)
                h_pass->Scale(1.0 / h_pass->Integral());
            if (h_fail->Integral() > 0)
                h_fail->Scale(1.0 / h_fail->Integral());

            TH1F *frame = new TH1F("frame_ET_ratio", ";#it{E}_{T}^{cluster} / #it{p}_{T}^{truth};Normalized", 100, 0, 2);
            float ymax = std::max(h_pass->GetMaximum(), h_fail->GetMaximum()) * 1.5;
            frame->GetYaxis()->SetRangeUser(0, ymax);
            frame->Draw("axis");

            h_pass->SetLineColor(kBlue);
            h_pass->SetLineWidth(2);
            h_pass->SetFillColorAlpha(kBlue, 0.15);
            h_pass->Draw("hist same");

            h_fail->SetLineColor(kRed);
            h_fail->SetLineWidth(2);
            h_fail->SetFillColorAlpha(kRed, 0.15);
            h_fail->Draw("hist same");

            TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.85);
            legStyle(leg, 0.20, 0.042);
            leg->AddEntry(h_pass, "#Delta#it{R} < 0.1", "lf");
            leg->AddEntry(h_fail, "#Delta#it{R} #geq 0.1", "lf");
            leg->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
            myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 0);
            myText(xpos, ypos - 3 * dy, 1, "Single-interaction MC", fontsize, 0);
        }

        c->SaveAs(Form("%s/ET_ratio_pass_fail.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Canvas 7: Fail fraction vs |dz|
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_fail_frac", "", 600, 600);

        TH1 *h_fail = (TH1 *)fin->Get("h_dz_abs_fail_combined");
        TH1 *h_all  = (TH1 *)fin->Get("h_dz_abs_all_combined");

        if (!h_fail || !h_all)
        {
            std::cerr << "Missing h_dz_abs_fail/all_combined" << std::endl;
        }
        else
        {
            TH1 *h_frac = (TH1 *)h_fail->Clone("h_fail_fraction");
            h_frac->Divide(h_fail, h_all, 1, 1, "B"); // binomial errors

            h_frac->SetTitle(";|#it{dz}| = |vtx_{reco} #minus vtx_{truth}| [cm];Fraction failing #Delta#it{R} > 0.1");
            h_frac->SetMarkerStyle(20);
            h_frac->SetMarkerSize(0.8);
            h_frac->SetLineColor(kBlack);
            h_frac->GetYaxis()->SetRangeUser(0, 1.0);
            h_frac->Draw("e");

            // Vertical line at critical threshold: dz = R_CEMC * deltaR_cut * cosh(0) = 93.5 * 0.1 = 9.35 cm
            TLine *lv = new TLine(9.35, 0, 9.35, 1.0);
            lv->SetLineStyle(7);
            lv->SetLineColor(kRed);
            lv->SetLineWidth(2);
            lv->Draw("same");

            TLegend *leg = new TLegend(0.45, 0.75, 0.88, 0.88);
            legStyle(leg, 0.20, 0.040);
            leg->AddEntry((TObject *)0, "Critical |#it{dz}| = 9.35 cm", "");
            leg->AddEntry((TObject *)0, "(#it{R}_{CEMC} #times #Delta#it{R}_{cut})", "");
            leg->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
            myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 0);
            myText(xpos, ypos - 3 * dy, 1, "Single-interaction MC", fontsize, 0);
        }

        c->SaveAs(Form("%s/fail_fraction_vs_dz.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Canvas 8: delta_eta vs dz (2D + analytic line)
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_deta_dz", "", 600, 600);
        gPad->SetRightMargin(0.15);

        TH2 *h2 = (TH2 *)fin->Get("h_deta_vs_dz_combined");
        if (!h2)
        {
            std::cerr << "Missing h_deta_vs_dz_combined" << std::endl;
        }
        else
        {
            h2->SetTitle(";#it{dz} = vtx_{reco} #minus vtx_{truth} [cm];#Delta#eta (cluster #minus truth)");
            h2->Draw("colz");

            // Analytic line: delta_eta = -dz / (R_CEMC * cosh(0)) = -dz / 93.5
            const double R_CEMC = 93.5;
            TGraph *g_analytic = new TGraph();
            for (int i = 0; i <= 60; ++i)
            {
                double dz = -30.0 + i;
                double deta = -dz / R_CEMC;
                g_analytic->SetPoint(i, dz, deta);
            }
            g_analytic->SetLineColor(kRed);
            g_analytic->SetLineWidth(2);
            g_analytic->SetLineStyle(2);
            g_analytic->Draw("L same");

            TLegend *leg = new TLegend(0.18, 0.75, 0.60, 0.88);
            legStyle(leg, 0.20, 0.040);
            leg->AddEntry(g_analytic, "#minus#it{dz} / #it{R}_{CEMC} (#eta = 0 approx.)", "l");
            leg->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
            myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 0);
        }

        c->SaveAs(Form("%s/deta_vs_dz.pdf", savePath.c_str()));
        delete c;
    }

    fin->Close();
    std::cout << "All plots saved to " << savePath << std::endl;
}
