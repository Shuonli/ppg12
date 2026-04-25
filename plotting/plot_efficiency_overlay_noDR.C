#include "plotcommon.h"

// =====================================================================
// plot_efficiency_overlay_noDR.C
//
// Overlay single / double / mixed efficiency curves on the same canvas,
// using the noDR (no deltaR truth-matching cut) TEfficiency objects.
//
// Input:  efficiency_comparison_deltaR.root
// Output: 6 PDFs in figures/efficiency_overlay/
// =====================================================================

void plot_efficiency_overlay_noDR()
{
    init_plot();

    const std::string infile   = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/efficiency_comparison_deltaR.root";
    const std::string savePath = "/sphenix/user/shuhangli/ppg12/plotting/figures/efficiency_overlay";
    gSystem->mkdir(savePath.c_str(), true);

    TFile *fin = TFile::Open(infile.c_str());
    if (!fin || fin->IsZombie())
    {
        std::cerr << "Cannot open " << infile << std::endl;
        return;
    }

    // ---- Label positions ----
    float xpos = 0.20, ypos = 0.885, dy = 0.054, fontsize = 0.046;

    // ---- Stage definitions ----
    const char *stageNames[]  = {"reco", "iso", "id", "all"};
    const char *stageLabels[] = {"Reconstruction", "Isolation", "Photon ID (BDT)", "Overall"};
    const int nStages = 4;

    // ---- Sample set definitions ----
    const char *sampleSets[]   = {"single", "double", "mixed"};
    const char *sampleLabels[] = {"Single-interaction MC",
                                  "Double-interaction MC",
                                  "Mixed 0 mrad (77.6%+22.4%)"};
    const int nSamples = 3;

    // Colors and markers per sample set
    int sampleColors[]  = {kBlue, kRed + 1, kBlack};
    int sampleMarkers[] = {20, 21, 22};  // filled circle, square, triangle

    // Y-axis ranges per stage
    double ylo[] = {0.0, 0.45, 0.3, 0.0};
    double yhi[] = {1.2, 1.1, 1.15, 0.75};

    // ================================================================
    // Helper lambda: retrieve a noDR TEfficiency and style it
    // ================================================================
    auto getStyledEff = [&](const char *stage, int iSample) -> TEfficiency *
    {
        TEfficiency *eff = (TEfficiency *)fin->Get(
            Form("eff_%s_noDR_%s", stage, sampleSets[iSample]));
        if (!eff)
        {
            std::cerr << "Missing: eff_" << stage << "_noDR_"
                      << sampleSets[iSample] << std::endl;
            return nullptr;
        }
        eff->SetMarkerColor(sampleColors[iSample]);
        eff->SetLineColor(sampleColors[iSample]);
        eff->SetMarkerStyle(sampleMarkers[iSample]);
        eff->SetMarkerSize(1.2);
        eff->SetLineWidth(2);
        return eff;
    };

    // ================================================================
    // Helper lambda: draw one panel with 3 overlaid efficiency curves
    // Returns true if at least one curve was drawn
    // ================================================================
    auto drawPanel = [&](int iStage, bool drawLegend, bool drawWatermark,
                         float legX1, float legY1, float legX2, float legY2) -> bool
    {
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.16);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.06);

        TH1F *frame = new TH1F(
            Form("frame_overlay_%s_%d", stageNames[iStage], gPad->GetNumber()),
            Form(";#it{E}_{T}^{#gamma, truth} [GeV];%s Efficiency", stageLabels[iStage]),
            100, 8, 36);
        frame->GetYaxis()->SetRangeUser(ylo[iStage], yhi[iStage]);
        frame->GetXaxis()->SetTitleSize(0.055);
        frame->GetYaxis()->SetTitleSize(0.055);
        frame->GetXaxis()->SetLabelSize(0.045);
        frame->GetYaxis()->SetLabelSize(0.045);
        frame->GetYaxis()->SetTitleOffset(1.3);
        frame->Draw("axis");

        bool anyDrawn = false;
        TLegend *leg = nullptr;
        if (drawLegend)
        {
            leg = new TLegend(legX1, legY1, legX2, legY2);
            legStyle(leg, 0.20, 0.040);
        }

        for (int iS = 0; iS < nSamples; ++iS)
        {
            TEfficiency *eff = getStyledEff(stageNames[iStage], iS);
            if (!eff) continue;
            eff->Draw("same");
            anyDrawn = true;
            if (leg) leg->AddEntry(eff, sampleLabels[iS], "pl");
        }

        if (leg) leg->Draw("same");

        if (drawWatermark)
        {
            myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
            myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
            myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 0);
        }

        return anyDrawn;
    };

    // ================================================================
    // Canvas 1: 2x2 four-panel overlay (THE KEY PLOT)
    // ================================================================
    {
        TCanvas *c = new TCanvas("c_4panel", "", 1200, 1000);
        c->Divide(2, 2, 0.005, 0.005);

        for (int iStage = 0; iStage < nStages; ++iStage)
        {
            c->cd(iStage + 1);
            drawPanel(iStage,
                      /* drawLegend */ true,
                      /* drawWatermark */ true,
                      /* legend box */ 0.38, 0.18, 0.92, 0.38);
        }

        c->SaveAs(Form("%s/eff_overlay_noDR_4panel.pdf", savePath.c_str()));
        delete c;
    }

    // ================================================================
    // Canvases 2-5: single-panel per stage
    // ================================================================
    for (int iStage = 0; iStage < nStages; ++iStage)
    {
        TCanvas *c = new TCanvas(Form("c_%s", stageNames[iStage]), "", 600, 600);

        drawPanel(iStage,
                  /* drawLegend */ true,
                  /* drawWatermark */ true,
                  /* legend box */ 0.42, 0.20, 0.92, 0.40);

        c->SaveAs(Form("%s/eff_%s_overlay_noDR.pdf", savePath.c_str(), stageNames[iStage]));
        delete c;
    }

    // ================================================================
    // Canvas 6: BDT pass rate for "new clusters" (pass trkID, fail dR)
    // ================================================================
    {
        TCanvas *c = new TCanvas("c_newclusters", "", 600, 600);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.16);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.06);

        TH1F *frame = new TH1F("frame_newclusters",
            ";#it{p}_{T}^{truth} [GeV];BDT Efficiency (new clusters only)",
            100, 8, 36);
        frame->GetYaxis()->SetRangeUser(0.0, 1.05);
        frame->GetXaxis()->SetTitleSize(0.055);
        frame->GetYaxis()->SetTitleSize(0.050);
        frame->GetXaxis()->SetLabelSize(0.045);
        frame->GetYaxis()->SetLabelSize(0.045);
        frame->GetYaxis()->SetTitleOffset(1.3);
        frame->Draw("axis");

        TLegend *leg = new TLegend(0.42, 0.20, 0.92, 0.40);
        legStyle(leg, 0.20, 0.040);

        for (int iS = 0; iS < nSamples; ++iS)
        {
            TEfficiency *eff = (TEfficiency *)fin->Get(
                Form("eff_id_newclusters_%s", sampleSets[iS]));
            if (!eff)
            {
                std::cerr << "Missing: eff_id_newclusters_" << sampleSets[iS] << std::endl;
                continue;
            }
            eff->SetMarkerColor(sampleColors[iS]);
            eff->SetLineColor(sampleColors[iS]);
            eff->SetMarkerStyle(sampleMarkers[iS]);
            eff->SetMarkerSize(1.2);
            eff->SetLineWidth(2);
            eff->Draw("same");
            leg->AddEntry(eff, sampleLabels[iS], "pl");
        }
        leg->Draw("same");

        myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
        myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
        myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 0);

        c->SaveAs(Form("%s/eff_id_newclusters_overlay.pdf", savePath.c_str()));
        delete c;
    }

    fin->Close();

    std::cout << "\n=== All plots saved to " << savePath << " ===" << std::endl;
    std::cout << "  eff_overlay_noDR_4panel.pdf" << std::endl;
    std::cout << "  eff_reco_overlay_noDR.pdf" << std::endl;
    std::cout << "  eff_iso_overlay_noDR.pdf" << std::endl;
    std::cout << "  eff_id_overlay_noDR.pdf" << std::endl;
    std::cout << "  eff_all_overlay_noDR.pdf" << std::endl;
    std::cout << "  eff_id_newclusters_overlay.pdf" << std::endl;
}
