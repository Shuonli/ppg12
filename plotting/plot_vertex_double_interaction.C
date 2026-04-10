#include "plotcommon.h"

void plot_vertex_double_interaction()
{
    init_plot();

    const std::string infile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/vertex_double_interaction_study.root";
    const std::string savePath = "/sphenix/user/shuhangli/ppg12/plotting/figures/vertex_study";
    const std::string savePath_dR = "/sphenix/user/shuhangli/ppg12/plotting/figures/deltaR_study";
    gSystem->mkdir(savePath.c_str(), true);
    gSystem->mkdir(savePath_dR.c_str(), true);

    TFile *fin = TFile::Open(infile.c_str());
    if (!fin || fin->IsZombie())
    {
        std::cerr << "Cannot open " << infile << std::endl;
        return;
    }

    float xpos = 0.20, ypos = 0.885, dy = 0.054, fontsize = 0.046;
    const double R_CEMC = 93.5;

    // ===================================================================
    // Figure 1: deta vs dz side-by-side (FOR DELTAR REPORT)
    //   Left: single (zoomed), Right: double (full range)
    //   Both with analytic line
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_deta_dz_comp", "", 1400, 600);
        c->Divide(2, 1);

        TH2 *h2_single = (TH2 *)fin->Get("h_deta_vs_dz_single");
        TH2 *h2_double = (TH2 *)fin->Get("h_deta_vs_dz_double");

        // Analytic: delta_eta = -dz / R_CEMC (at eta = 0)
        auto makeAnalytic = [&](double dz_lo, double dz_hi) -> TGraph * {
            TGraph *g = new TGraph();
            int np = 0;
            for (double dz = dz_lo; dz <= dz_hi; dz += 0.5)
            {
                g->SetPoint(np++, dz, -dz / R_CEMC);
            }
            g->SetLineColor(kRed);
            g->SetLineWidth(2);
            g->SetLineStyle(2);
            return g;
        };

        // Left: single interaction (zoomed)
        c->cd(1);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.14);
        if (h2_single)
        {
            h2_single->SetTitle(";#it{dz} = vtx_{reco} #minus vtx_{truth} [cm];#Delta#eta (cluster #minus truth)");
            h2_single->GetXaxis()->SetRangeUser(-30, 30);
            h2_single->GetYaxis()->SetRangeUser(-0.15, 0.15);
            h2_single->SetMinimum(1);
            h2_single->Draw("colz");
            gPad->SetLogz();

            TGraph *g1 = makeAnalytic(-30, 30);
            g1->Draw("L same");

            TLegend *leg1 = new TLegend(0.16, 0.72, 0.55, 0.82);
            legStyle(leg1, 0.20, 0.038);
            leg1->AddEntry(g1, "#minus#it{dz} / #it{R}_{CEMC}", "l");
            leg1->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "Single interaction", fontsize, 0);
            myText(xpos, ypos - 2 * dy, 1, "photon10 MC", fontsize, 0);
        }

        // Right: double interaction (full range)
        c->cd(2);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.14);
        if (h2_double)
        {
            h2_double->SetTitle(";#it{dz} = vtx_{reco} #minus vtx_{truth} [cm];#Delta#eta (cluster #minus truth)");
            h2_double->GetXaxis()->SetRangeUser(-120, 120);
            h2_double->GetYaxis()->SetRangeUser(-1.2, 1.2);
            h2_double->SetMinimum(1);
            h2_double->Draw("colz");
            gPad->SetLogz();

            TGraph *g2 = makeAnalytic(-120, 120);
            g2->Draw("L same");

            TLegend *leg2 = new TLegend(0.16, 0.72, 0.55, 0.82);
            legStyle(leg2, 0.20, 0.038);
            leg2->AddEntry(g2, "#minus#it{dz} / #it{R}_{CEMC}", "l");
            leg2->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "Double interaction", fontsize, 0);
            myText(xpos, ypos - 2 * dy, 1, "photon10_double MC", fontsize, 0);
        }

        c->SaveAs(Form("%s/deta_vs_dz_comparison.pdf", savePath_dR.c_str()));
        delete c;
    }

    // ===================================================================
    // Figure 2: Vertex distribution comparison (reco)
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_vtx_reco", "", 600, 600);

        TH1 *h_single = (TH1 *)fin->Get("h_vtxz_reco_single");
        TH1 *h_double = (TH1 *)fin->Get("h_vtxz_reco_double");

        if (h_single && h_double)
        {
            // Normalize to unit area
            if (h_single->Integral() > 0) h_single->Scale(1.0 / h_single->Integral());
            if (h_double->Integral() > 0) h_double->Scale(1.0 / h_double->Integral());

            TH1F *frame = new TH1F("frame_vtx", ";#it{z}_{vtx}^{reco} [cm];Normalized", 100, -120, 120);
            float ymax = std::max(h_single->GetMaximum(), h_double->GetMaximum()) * 1.5;
            frame->GetYaxis()->SetRangeUser(0, ymax);
            frame->Draw("axis");

            h_single->SetLineColor(kBlue);
            h_single->SetLineWidth(2);
            h_single->SetFillColorAlpha(kBlue, 0.10);
            h_single->Draw("hist same");

            h_double->SetLineColor(kRed);
            h_double->SetLineWidth(2);
            h_double->SetFillColorAlpha(kRed, 0.10);
            h_double->Draw("hist same");

            TLegend *leg = new TLegend(0.52, 0.60, 0.88, 0.73);
            legStyle(leg, 0.20, 0.042);
            leg->AddEntry(h_single, "Single interaction", "lf");
            leg->AddEntry(h_double, "Double interaction", "lf");
            leg->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "photon10 MC", fontsize, 0);
            myText(xpos, ypos - 2 * dy, 1, "Reconstructed vertex", fontsize, 0);

            float fs2 = 0.038;
            myText(0.52, 0.55, kBlue, Form("Single #sigma = %.1f cm", h_single->GetStdDev()), fs2, 0);
            myText(0.52, 0.50, kRed,  Form("Double #sigma = %.1f cm", h_double->GetStdDev()), fs2, 0);
        }

        c->SaveAs(Form("%s/vertex_reco_comparison.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Figure 3: Truth vertex distribution comparison
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_vtx_truth", "", 600, 600);

        TH1 *h_single = (TH1 *)fin->Get("h_vtxz_truth_single");
        TH1 *h_double = (TH1 *)fin->Get("h_vtxz_truth_double");

        if (h_single && h_double)
        {
            if (h_single->Integral() > 0) h_single->Scale(1.0 / h_single->Integral());
            if (h_double->Integral() > 0) h_double->Scale(1.0 / h_double->Integral());

            TH1F *frame = new TH1F("frame_vtx_truth", ";#it{z}_{vtx}^{truth} [cm];Normalized", 100, -120, 120);
            float ymax = std::max(h_single->GetMaximum(), h_double->GetMaximum()) * 1.5;
            frame->GetYaxis()->SetRangeUser(0, ymax);
            frame->Draw("axis");

            h_single->SetLineColor(kBlue);
            h_single->SetLineWidth(2);
            h_single->SetFillColorAlpha(kBlue, 0.10);
            h_single->Draw("hist same");

            h_double->SetLineColor(kRed);
            h_double->SetLineWidth(2);
            h_double->SetFillColorAlpha(kRed, 0.10);
            h_double->Draw("hist same");

            TLegend *leg = new TLegend(0.52, 0.60, 0.88, 0.73);
            legStyle(leg, 0.20, 0.042);
            leg->AddEntry(h_single, "Single interaction", "lf");
            leg->AddEntry(h_double, "Double interaction", "lf");
            leg->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "photon10 MC", fontsize, 0);
            myText(xpos, ypos - 2 * dy, 1, "Truth vertex (primary collision)", fontsize, 0);

            float fs2 = 0.038;
            myText(0.52, 0.55, kBlue, Form("Single #sigma = %.1f cm", h_single->GetStdDev()), fs2, 0);
            myText(0.52, 0.50, kRed,  Form("Double #sigma = %.1f cm", h_double->GetStdDev()), fs2, 0);
        }

        c->SaveAs(Form("%s/vertex_truth_comparison.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Figure 4: dz distribution comparison
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_dz", "", 600, 600);

        TH1 *h_single = (TH1 *)fin->Get("h_dz_event_single");
        TH1 *h_double = (TH1 *)fin->Get("h_dz_event_double");

        if (h_single && h_double)
        {
            if (h_single->Integral() > 0) h_single->Scale(1.0 / h_single->Integral());
            if (h_double->Integral() > 0) h_double->Scale(1.0 / h_double->Integral());

            TH1F *frame = new TH1F("frame_dz", ";#it{dz} = #it{z}_{vtx}^{reco} #minus #it{z}_{vtx}^{truth} [cm];Normalized", 100, -150, 150);
            float ymax = std::max(h_single->GetMaximum(), h_double->GetMaximum()) * 5;
            frame->GetYaxis()->SetRangeUser(1e-5, ymax);
            frame->Draw("axis");

            h_single->SetLineColor(kBlue);
            h_single->SetLineWidth(2);
            h_single->Draw("hist same");

            h_double->SetLineColor(kRed);
            h_double->SetLineWidth(2);
            h_double->Draw("hist same");

            // Critical threshold lines at +/- 9.35 cm
            TLine *l1 = new TLine(9.35, 1e-5, 9.35, ymax);
            l1->SetLineStyle(7); l1->SetLineColor(kGray + 2); l1->SetLineWidth(2);
            l1->Draw("same");
            TLine *l2 = new TLine(-9.35, 1e-5, -9.35, ymax);
            l2->SetLineStyle(7); l2->SetLineColor(kGray + 2); l2->SetLineWidth(2);
            l2->Draw("same");

            TLegend *leg = new TLegend(0.52, 0.68, 0.88, 0.85);
            legStyle(leg, 0.20, 0.040);
            leg->AddEntry(h_single, "Single interaction", "l");
            leg->AddEntry(h_double, "Double interaction", "l");
            leg->AddEntry(l1, "|#it{dz}| = 9.35 cm", "l");
            leg->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "photon10 MC", fontsize, 0);

            float fs2 = 0.038;
            myText(xpos, ypos - 3 * dy, kBlue, Form("Single #sigma_{dz} = %.1f cm", h_single->GetStdDev()), fs2, 0);
            myText(xpos, ypos - 4 * dy, kRed,  Form("Double #sigma_{dz} = %.1f cm", h_double->GetStdDev()), fs2, 0);

            gPad->SetLogy();
        }

        c->SaveAs(Form("%s/dz_distribution_comparison.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Figure 5: |dz| distribution with critical threshold
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_dz_abs", "", 600, 600);

        TH1 *h_single = (TH1 *)fin->Get("h_dz_abs_event_single");
        TH1 *h_double = (TH1 *)fin->Get("h_dz_abs_event_double");

        if (h_single && h_double)
        {
            if (h_single->Integral() > 0) h_single->Scale(1.0 / h_single->Integral());
            if (h_double->Integral() > 0) h_double->Scale(1.0 / h_double->Integral());

            TH1F *frame = new TH1F("frame_dz_abs", ";|#it{dz}| [cm];Normalized", 100, 0, 100);
            float ymax = std::max(h_single->GetMaximum(), h_double->GetMaximum()) * 5;
            frame->GetYaxis()->SetRangeUser(1e-5, ymax);
            frame->Draw("axis");

            h_single->SetLineColor(kBlue);
            h_single->SetLineWidth(2);
            h_single->Draw("hist same");

            h_double->SetLineColor(kRed);
            h_double->SetLineWidth(2);
            h_double->Draw("hist same");

            // Critical threshold
            TLine *lv = new TLine(9.35, 1e-5, 9.35, ymax);
            lv->SetLineStyle(7); lv->SetLineColor(kBlack); lv->SetLineWidth(2);
            lv->Draw("same");

            // Shade the region > 9.35 cm for both distributions
            TLegend *leg = new TLegend(0.45, 0.68, 0.88, 0.88);
            legStyle(leg, 0.20, 0.038);
            leg->AddEntry(h_single, "Single interaction", "l");
            leg->AddEntry(h_double, "Double interaction", "l");
            leg->AddEntry(lv, Form("|#it{dz}|_{crit} = 9.35 cm"), "l");
            leg->AddEntry((TObject *)0, Form("= #Delta#it{R}_{cut} #times #it{R}_{CEMC}"), "");
            leg->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "photon10 MC", fontsize, 0);

            float fs2 = 0.038;
            myText(xpos, ypos - 3 * dy, kBlue, Form("Single #sigma_{dz} = %.1f cm", h_single->GetStdDev()), fs2, 0);
            myText(xpos, ypos - 4 * dy, kRed,  Form("Double #sigma_{dz} = %.1f cm", h_double->GetStdDev()), fs2, 0);

            gPad->SetLogy();
        }

        c->SaveAs(Form("%s/dz_abs_comparison.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Figure 6: 2D reco vs truth vertex (side-by-side)
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_vtx_2d", "", 1400, 600);
        c->Divide(2, 1);

        TH2 *h2_single = (TH2 *)fin->Get("h_vtx_reco_vs_truth_single");
        TH2 *h2_double = (TH2 *)fin->Get("h_vtx_reco_vs_truth_double");

        TGraph *g_diag = new TGraph();
        g_diag->SetPoint(0, -120, -120);
        g_diag->SetPoint(1, 120, 120);
        g_diag->SetLineColor(kRed);
        g_diag->SetLineWidth(2);
        g_diag->SetLineStyle(2);

        c->cd(1);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.14);
        if (h2_single)
        {
            h2_single->SetTitle(";#it{z}_{vtx}^{truth} [cm];#it{z}_{vtx}^{reco} [cm]");
            h2_single->SetMinimum(1);
            h2_single->Draw("colz");
            gPad->SetLogz();
            g_diag->Draw("L same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "Single interaction", fontsize, 0);
        }

        c->cd(2);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.14);
        if (h2_double)
        {
            h2_double->SetTitle(";#it{z}_{vtx}^{truth} [cm];#it{z}_{vtx}^{reco} [cm]");
            h2_double->SetMinimum(1);
            h2_double->Draw("colz");
            gPad->SetLogz();
            g_diag->Draw("L same");

            TLegend *leg = new TLegend(0.16, 0.75, 0.50, 0.85);
            legStyle(leg, 0.20, 0.038);
            leg->AddEntry(g_diag, "diagonal", "l");
            leg->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "Double interaction", fontsize, 0);
        }

        c->SaveAs(Form("%s/vertex_reco_vs_truth.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Figure 7: deltaR distribution comparison
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_deltaR", "", 600, 600);

        TH1 *h_single = (TH1 *)fin->Get("h_deltaR_single");
        TH1 *h_double = (TH1 *)fin->Get("h_deltaR_double");

        if (h_single && h_double)
        {
            if (h_single->Integral() > 0) h_single->Scale(1.0 / h_single->Integral());
            if (h_double->Integral() > 0) h_double->Scale(1.0 / h_double->Integral());

            TH1F *frame = new TH1F("frame_deltaR", ";#Delta#it{R}(cluster, truth);Normalized", 100, 0, 1.5);
            float ymax = std::max(h_single->GetMaximum(), h_double->GetMaximum()) * 5;
            frame->GetYaxis()->SetRangeUser(1e-4, ymax);
            frame->Draw("axis");

            h_single->SetLineColor(kBlue);
            h_single->SetLineWidth(2);
            h_single->Draw("hist same");

            h_double->SetLineColor(kRed);
            h_double->SetLineWidth(2);
            h_double->Draw("hist same");

            TLine *lcut = new TLine(0.1, 1e-4, 0.1, ymax);
            lcut->SetLineStyle(7); lcut->SetLineColor(kGray + 2); lcut->SetLineWidth(2);
            lcut->Draw("same");

            TLegend *leg = new TLegend(0.45, 0.68, 0.88, 0.85);
            legStyle(leg, 0.20, 0.040);
            leg->AddEntry(h_single, "Single interaction", "l");
            leg->AddEntry(h_double, "Double interaction", "l");
            leg->AddEntry(lcut, "#Delta#it{R} = 0.1 cut", "l");
            leg->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "photon10 MC, trkID matched", fontsize, 0);

            gPad->SetLogy();
        }

        c->SaveAs(Form("%s/deltaR_comparison.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Figure 8: deta distribution comparison
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_deta", "", 600, 600);

        TH1 *h_single = (TH1 *)fin->Get("h_deta_single");
        TH1 *h_double = (TH1 *)fin->Get("h_deta_double");

        if (h_single && h_double)
        {
            if (h_single->Integral() > 0) h_single->Scale(1.0 / h_single->Integral());
            if (h_double->Integral() > 0) h_double->Scale(1.0 / h_double->Integral());

            TH1F *frame = new TH1F("frame_deta", ";#Delta#eta (cluster #minus truth);Normalized", 100, -1.0, 1.0);
            float ymax = std::max(h_single->GetMaximum(), h_double->GetMaximum()) * 5;
            frame->GetYaxis()->SetRangeUser(1e-4, ymax);
            frame->Draw("axis");

            h_single->SetLineColor(kBlue);
            h_single->SetLineWidth(2);
            h_single->Draw("hist same");

            h_double->SetLineColor(kRed);
            h_double->SetLineWidth(2);
            h_double->Draw("hist same");

            TLegend *leg = new TLegend(0.52, 0.72, 0.88, 0.85);
            legStyle(leg, 0.20, 0.042);
            leg->AddEntry(h_single, "Single interaction", "l");
            leg->AddEntry(h_double, "Double interaction", "l");
            leg->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "photon10 MC, trkID matched", fontsize, 0);

            gPad->SetLogy();
        }

        c->SaveAs(Form("%s/deta_comparison.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Figure 9: deta vs dz for double interaction ONLY (standalone view)
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_deta_dz_double", "", 700, 600);
        gPad->SetRightMargin(0.15);

        TH2 *h2 = (TH2 *)fin->Get("h_deta_vs_dz_double");
        if (h2)
        {
            h2->SetTitle(";#it{dz} = vtx_{reco} #minus vtx_{truth} [cm];#Delta#eta (cluster #minus truth)");
            h2->GetXaxis()->SetRangeUser(-120, 120);
            h2->GetYaxis()->SetRangeUser(-1.2, 1.2);
            h2->SetMinimum(1);
            h2->Draw("colz");
            gPad->SetLogz();

            TGraph *g = new TGraph();
            int np = 0;
            for (double dz = -120; dz <= 120; dz += 0.5)
                g->SetPoint(np++, dz, -dz / R_CEMC);
            g->SetLineColor(kRed);
            g->SetLineWidth(2);
            g->SetLineStyle(2);
            g->Draw("L same");

            // dR = 0.1 threshold lines (deta = +/- 0.1 at dphi ~ 0)
            TLine *lh1 = new TLine(-120, 0.1, 120, 0.1);
            lh1->SetLineStyle(3); lh1->SetLineColor(kGray + 1); lh1->SetLineWidth(1);
            lh1->Draw("same");
            TLine *lh2 = new TLine(-120, -0.1, 120, -0.1);
            lh2->SetLineStyle(3); lh2->SetLineColor(kGray + 1); lh2->SetLineWidth(1);
            lh2->Draw("same");

            TLegend *leg = new TLegend(0.16, 0.18, 0.55, 0.31);
            legStyle(leg, 0.20, 0.038);
            leg->AddEntry(g, "#minus#it{dz} / #it{R}_{CEMC}", "l");
            leg->AddEntry(lh1, "#Delta#eta = #pm 0.1", "l");
            leg->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "Double interaction MC", fontsize, 0);
            myText(xpos, ypos - 2 * dy, 1, "photon10_double, trkID matched", fontsize, 0);
        }

        c->SaveAs(Form("%s/deta_vs_dz_double.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Figure 10: deta vs eta side-by-side
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_deta_eta", "", 1400, 600);
        c->Divide(2, 1);

        TH2 *h2_single = (TH2 *)fin->Get("h_deta_vs_eta_single");
        TH2 *h2_double = (TH2 *)fin->Get("h_deta_vs_eta_double");

        c->cd(1);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.14);
        if (h2_single)
        {
            h2_single->SetTitle(";#eta_{truth};#Delta#eta (cluster #minus truth)");
            h2_single->GetYaxis()->SetRangeUser(-0.15, 0.15);
            h2_single->SetMinimum(1);
            h2_single->Draw("colz");
            gPad->SetLogz();
            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "Single interaction", fontsize, 0);
        }

        c->cd(2);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.14);
        if (h2_double)
        {
            h2_double->SetTitle(";#eta_{truth};#Delta#eta (cluster #minus truth)");
            h2_double->GetYaxis()->SetRangeUser(-1.0, 1.0);
            h2_double->SetMinimum(1);
            h2_double->Draw("colz");
            gPad->SetLogz();
            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "Double interaction", fontsize, 0);
        }

        c->SaveAs(Form("%s/deta_vs_eta_comparison.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Figure 11: ET ratio comparison
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_ET_ratio", "", 600, 600);

        TH1 *h_single = (TH1 *)fin->Get("h_ET_ratio_single");
        TH1 *h_double = (TH1 *)fin->Get("h_ET_ratio_double");

        if (h_single && h_double)
        {
            if (h_single->Integral() > 0) h_single->Scale(1.0 / h_single->Integral());
            if (h_double->Integral() > 0) h_double->Scale(1.0 / h_double->Integral());

            TH1F *frame = new TH1F("frame_ET", ";#it{E}_{T}^{cluster} / #it{p}_{T}^{truth};Normalized", 100, 0, 2);
            float ymax = std::max(h_single->GetMaximum(), h_double->GetMaximum()) * 1.5;
            frame->GetYaxis()->SetRangeUser(0, ymax);
            frame->Draw("axis");

            h_single->SetLineColor(kBlue);
            h_single->SetLineWidth(2);
            h_single->SetFillColorAlpha(kBlue, 0.10);
            h_single->Draw("hist same");

            h_double->SetLineColor(kRed);
            h_double->SetLineWidth(2);
            h_double->SetFillColorAlpha(kRed, 0.10);
            h_double->Draw("hist same");

            TLegend *leg = new TLegend(0.52, 0.72, 0.88, 0.85);
            legStyle(leg, 0.20, 0.042);
            leg->AddEntry(h_single, "Single interaction", "lf");
            leg->AddEntry(h_double, "Double interaction", "lf");
            leg->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "photon10 MC, trkID matched", fontsize, 0);
        }

        c->SaveAs(Form("%s/ET_ratio_comparison.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Figure 12: deltaR vs dz side-by-side
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_dR_dz", "", 1400, 600);
        c->Divide(2, 1);

        TH2 *h2_single = (TH2 *)fin->Get("h_deltaR_vs_dz_single");
        TH2 *h2_double = (TH2 *)fin->Get("h_deltaR_vs_dz_double");

        auto makeAnalyticDR = [&](double dz_lo, double dz_hi) -> TGraph * {
            TGraph *g = new TGraph();
            int np = 0;
            for (double dz = dz_lo; dz <= dz_hi; dz += 0.5)
                g->SetPoint(np++, dz, fabs(dz) / R_CEMC);
            g->SetLineColor(kRed);
            g->SetLineWidth(2);
            g->SetLineStyle(2);
            return g;
        };

        c->cd(1);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.14);
        if (h2_single)
        {
            h2_single->SetTitle(";#it{dz} [cm];#Delta#it{R}(cluster, truth)");
            h2_single->GetXaxis()->SetRangeUser(-30, 30);
            h2_single->GetYaxis()->SetRangeUser(0, 0.3);
            h2_single->SetMinimum(1);
            h2_single->Draw("colz");
            gPad->SetLogz();

            TGraph *g1 = makeAnalyticDR(-30, 30);
            g1->Draw("L same");

            TLine *lh = new TLine(-30, 0.1, 30, 0.1);
            lh->SetLineStyle(7); lh->SetLineColor(kGray + 2); lh->SetLineWidth(2);
            lh->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "Single interaction", fontsize, 0);
        }

        c->cd(2);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.14);
        if (h2_double)
        {
            h2_double->SetTitle(";#it{dz} [cm];#Delta#it{R}(cluster, truth)");
            h2_double->GetXaxis()->SetRangeUser(-120, 120);
            h2_double->GetYaxis()->SetRangeUser(0, 1.5);
            h2_double->SetMinimum(1);
            h2_double->Draw("colz");
            gPad->SetLogz();

            TGraph *g2 = makeAnalyticDR(-120, 120);
            g2->Draw("L same");

            TLine *lh = new TLine(-120, 0.1, 120, 0.1);
            lh->SetLineStyle(7); lh->SetLineColor(kGray + 2); lh->SetLineWidth(2);
            lh->Draw("same");

            TLegend *leg = new TLegend(0.16, 0.70, 0.60, 0.85);
            legStyle(leg, 0.20, 0.038);
            leg->AddEntry(g2, "|#it{dz}| / #it{R}_{CEMC}", "l");
            leg->AddEntry(lh, "#Delta#it{R} = 0.1", "l");
            leg->Draw("same");

            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, "Double interaction", fontsize, 0);
        }

        c->SaveAs(Form("%s/deltaR_vs_dz_comparison.pdf", savePath.c_str()));
        delete c;
    }

    fin->Close();
    std::cout << "\nAll plots saved to " << savePath << " and " << savePath_dR << std::endl;
}
