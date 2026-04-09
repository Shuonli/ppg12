#include "plotcommon.h"

// Helper: extract efficiency values from TEfficiency into arrays
// Returns number of bins filled
int extractEfficiency(TEfficiency *eff, double *xvals, double *yvals,
                      double *eylo, double *eyhi, int maxBins)
{
    if (!eff) return 0;
    int n = eff->GetTotalHistogram()->GetNbinsX();
    if (n > maxBins) n = maxBins;
    for (int i = 0; i < n; ++i)
    {
        xvals[i] = eff->GetTotalHistogram()->GetBinCenter(i + 1);
        yvals[i] = eff->GetEfficiency(i + 1);
        eylo[i]  = eff->GetEfficiencyErrorLow(i + 1);
        eyhi[i]  = eff->GetEfficiencyErrorUp(i + 1);
    }
    return n;
}

// Helper: build a ratio TGraphAsymmErrors from two TEfficiency objects (num / den)
TGraphAsymmErrors *makeRatioGraph(TEfficiency *effNum, TEfficiency *effDen,
                                  const char *name)
{
    const int kMax = 50;
    double xN[kMax], yN[kMax], eloN[kMax], ehiN[kMax];
    double xD[kMax], yD[kMax], eloD[kMax], ehiD[kMax];
    int nN = extractEfficiency(effNum, xN, yN, eloN, ehiN, kMax);
    int nD = extractEfficiency(effDen, xD, yD, eloD, ehiD, kMax);
    int n = std::min(nN, nD);

    TGraphAsymmErrors *gr = new TGraphAsymmErrors(n);
    gr->SetName(name);
    for (int i = 0; i < n; ++i)
    {
        double ratio = (yD[i] > 0) ? yN[i] / yD[i] : 0;
        // Propagate errors assuming uncorrelated
        double relEloN = (yN[i] > 0) ? eloN[i] / yN[i] : 0;
        double relEhiN = (yN[i] > 0) ? ehiN[i] / yN[i] : 0;
        double relEloD = (yD[i] > 0) ? eloD[i] / yD[i] : 0;
        double relEhiD = (yD[i] > 0) ? ehiD[i] / yD[i] : 0;
        double elo = ratio * std::sqrt(relEloN * relEloN + relEhiD * relEhiD);
        double ehi = ratio * std::sqrt(relEhiN * relEhiN + relEloD * relEloD);
        double ex = effNum->GetTotalHistogram()->GetBinWidth(i + 1) / 2.0;
        gr->SetPoint(i, xN[i], ratio);
        gr->SetPointError(i, ex, ex, elo, ehi);
    }
    return gr;
}

void plot_efficiency_comparison_deltaR()
{
    init_plot();

    const std::string infile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/efficiency_comparison_deltaR.root";
    const std::string savePath = "/sphenix/user/shuhangli/ppg12/plotting/figures/efficiency_comparison";
    gSystem->mkdir(savePath.c_str(), true);

    TFile *fin = TFile::Open(infile.c_str());
    if (!fin || fin->IsZombie())
    {
        std::cerr << "Cannot open " << infile << std::endl;
        return;
    }

    // Standard label positions
    float xpos = 0.20, ypos = 0.885, dy = 0.054, fontsize = 0.046;

    // Stage labels for axis titles and legends
    const char *stageNames[]  = {"reco", "iso", "id", "all"};
    const char *stageLabels[] = {"Reconstruction", "Isolation", "Photon ID (BDT)", "Overall"};
    const int nStages = 4;

    // Sample set info
    const char *sampleSets[]  = {"single", "double", "mixed"};
    const char *sampleLabels[] = {"Single-interaction MC", "Double-interaction MC",
                                   "Mixed 0 mrad (81.3% single + 18.7% double)"};
    const int nSamples = 3;

    // Y-axis ranges per stage
    double ylo_stage[] = {0.0, 0.45, 0.4, 0.0};
    double yhi_stage[] = {1.2, 1.1, 1.15, 0.75};

    // ===================================================================
    // Helper lambda: draw a 2x2 per-stage efficiency canvas for one sample set
    // ===================================================================
    auto drawPerStageCanvas = [&](int iSample, const char *pdfName)
    {
        TCanvas *c = new TCanvas(Form("c_perstage_%s", sampleSets[iSample]), "", 1200, 1000);
        c->Divide(2, 2, 0.005, 0.005);

        for (int iStage = 0; iStage < nStages; ++iStage)
        {
            c->cd(iStage + 1);
            gPad->SetLeftMargin(0.16);
            gPad->SetBottomMargin(0.14);
            gPad->SetRightMargin(0.04);
            gPad->SetTopMargin(0.06);

            TH1F *frame = new TH1F(Form("frame_%s_%s", stageNames[iStage], sampleSets[iSample]),
                                   Form(";#it{E}_{T}^{#gamma, truth} [GeV];%s Efficiency", stageLabels[iStage]),
                                   100, 8, 36);
            frame->GetYaxis()->SetRangeUser(ylo_stage[iStage], yhi_stage[iStage]);
            frame->GetXaxis()->SetTitleSize(0.055);
            frame->GetYaxis()->SetTitleSize(0.055);
            frame->GetXaxis()->SetLabelSize(0.045);
            frame->GetYaxis()->SetLabelSize(0.045);
            frame->GetYaxis()->SetTitleOffset(1.3);
            frame->Draw("axis");

            // withDR
            TEfficiency *effDR = (TEfficiency *)fin->Get(
                Form("eff_%s_withDR_%s", stageNames[iStage], sampleSets[iSample]));
            // noDR
            TEfficiency *effNoDR = (TEfficiency *)fin->Get(
                Form("eff_%s_noDR_%s", stageNames[iStage], sampleSets[iSample]));

            if (effDR)
            {
                effDR->SetMarkerColor(kRed + 1);
                effDR->SetLineColor(kRed + 1);
                effDR->SetMarkerStyle(24); // open circle
                effDR->SetMarkerSize(1.0);
                effDR->SetLineStyle(2);
                effDR->SetLineWidth(2);
                effDR->Draw("same");
            }
            else
            {
                std::cerr << "Missing eff_" << stageNames[iStage]
                          << "_withDR_" << sampleSets[iSample] << std::endl;
            }

            if (effNoDR)
            {
                effNoDR->SetMarkerColor(kBlue);
                effNoDR->SetLineColor(kBlue);
                effNoDR->SetMarkerStyle(20); // filled circle
                effNoDR->SetMarkerSize(1.0);
                effNoDR->SetLineStyle(1);
                effNoDR->SetLineWidth(2);
                effNoDR->Draw("same");
            }
            else
            {
                std::cerr << "Missing eff_" << stageNames[iStage]
                          << "_noDR_" << sampleSets[iSample] << std::endl;
            }

            // Legend
            TLegend *leg = new TLegend(0.45, 0.20, 0.92, 0.40);
            legStyle(leg, 0.20, 0.042);
            if (effNoDR) leg->AddEntry(effNoDR, "No #Delta#it{R} cut", "pl");
            if (effDR)   leg->AddEntry(effDR, "With #Delta#it{R} < 0.1", "pl");
            leg->Draw("same");

            // Labels only in top-left panel
            if (iStage == 0)
            {
                myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
                myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
                myText(xpos, ypos - 2 * dy, 1, sampleLabels[iSample], fontsize, 0);
                myText(xpos, ypos - 3 * dy, 1, strleg3.c_str(), fontsize, 0);
            }
        }

        c->SaveAs(Form("%s/%s", savePath.c_str(), pdfName));
        delete c;
    };

    // ===================================================================
    // Canvas 1-3: per-stage efficiency for single, double, mixed
    // ===================================================================
    drawPerStageCanvas(0, "efficiency_per_stage_single.pdf");
    drawPerStageCanvas(1, "efficiency_per_stage_double.pdf");
    drawPerStageCanvas(2, "efficiency_per_stage_mixed.pdf");

    // ===================================================================
    // Canvas 4: Ratio plot  eps(noDR) / eps(withDR)
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_ratio", "", 1200, 1000);
        c->Divide(2, 2, 0.005, 0.005);

        // Colors: single = blue, mixed = red
        int ratioColors[] = {kBlue, kRed + 1};
        int ratioMarkers[] = {20, 21};
        const char *ratioSamples[] = {"single", "mixed"};
        const char *ratioLabels[] = {"Single MC", "Mixed 0 mrad"};
        const int nRatio = 2;

        for (int iStage = 0; iStage < nStages; ++iStage)
        {
            c->cd(iStage + 1);
            gPad->SetLeftMargin(0.16);
            gPad->SetBottomMargin(0.14);
            gPad->SetRightMargin(0.04);
            gPad->SetTopMargin(0.06);

            TH1F *frame = new TH1F(
                Form("frame_ratio_%s", stageNames[iStage]),
                Form(";#it{E}_{T}^{#gamma, truth} [GeV];#varepsilon(no #Delta#it{R}) / #varepsilon(#Delta#it{R} < 0.1)  [%s]",
                     stageLabels[iStage]),
                100, 8, 36);
            frame->GetYaxis()->SetRangeUser(0.90, 1.25);
            frame->GetXaxis()->SetTitleSize(0.055);
            frame->GetYaxis()->SetTitleSize(0.045);
            frame->GetXaxis()->SetLabelSize(0.045);
            frame->GetYaxis()->SetLabelSize(0.045);
            frame->GetYaxis()->SetTitleOffset(1.4);
            frame->Draw("axis");

            // Horizontal line at 1.0
            TLine *lunit = new TLine(8, 1.0, 36, 1.0);
            lunit->SetLineStyle(7);
            lunit->SetLineColor(kGray + 2);
            lunit->SetLineWidth(2);
            lunit->Draw("same");

            TLegend *leg = new TLegend(0.45, 0.68, 0.92, 0.88);
            legStyle(leg, 0.20, 0.042);

            for (int iR = 0; iR < nRatio; ++iR)
            {
                TEfficiency *effNoDR = (TEfficiency *)fin->Get(
                    Form("eff_%s_noDR_%s", stageNames[iStage], ratioSamples[iR]));
                TEfficiency *effDR = (TEfficiency *)fin->Get(
                    Form("eff_%s_withDR_%s", stageNames[iStage], ratioSamples[iR]));
                if (!effNoDR || !effDR)
                {
                    std::cerr << "Missing TEfficiency for ratio: " << stageNames[iStage]
                              << " " << ratioSamples[iR] << std::endl;
                    continue;
                }

                TGraphAsymmErrors *gr = makeRatioGraph(effNoDR, effDR,
                    Form("gr_ratio_%s_%s", stageNames[iStage], ratioSamples[iR]));
                gr->SetMarkerColor(ratioColors[iR]);
                gr->SetLineColor(ratioColors[iR]);
                gr->SetMarkerStyle(ratioMarkers[iR]);
                gr->SetMarkerSize(1.0);
                gr->SetLineWidth(2);
                gr->Draw("P same");

                leg->AddEntry(gr, ratioLabels[iR], "pl");
            }
            leg->Draw("same");

            if (iStage == 0)
            {
                myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
                myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
                myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 0);
            }
        }

        c->SaveAs(Form("%s/efficiency_ratio.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Canvas 5: Photon ID efficiency focus (2 panels)
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_id_focus", "", 1400, 600);
        c->Divide(2, 1, 0.005, 0.005);

        // Colors per sample set: single = blue, double = red, mixed = black
        int sampleColors[] = {kBlue, kRed + 1, kBlack};
        // withDR = open markers, noDR = filled markers
        int openMarkers[]  = {24, 25, 26}; // circle, square, triangle
        int fillMarkers[]  = {20, 21, 22};

        // ---------------------------------------------------------------
        // Left panel: eff_id withDR vs noDR for all three sample sets
        // ---------------------------------------------------------------
        c->cd(1);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.14);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.06);

        TH1F *frameL = new TH1F("frame_id_left",
            ";#it{E}_{T}^{#gamma, truth} [GeV];Photon ID (BDT) Efficiency", 100, 8, 36);
        frameL->GetYaxis()->SetRangeUser(0.5, 1.05);
        frameL->GetXaxis()->SetTitleSize(0.055);
        frameL->GetYaxis()->SetTitleSize(0.055);
        frameL->GetXaxis()->SetLabelSize(0.045);
        frameL->GetYaxis()->SetLabelSize(0.045);
        frameL->GetYaxis()->SetTitleOffset(1.3);
        frameL->Draw("axis");

        TLegend *legL = new TLegend(0.35, 0.18, 0.92, 0.52);
        legStyle(legL, 0.22, 0.036);

        for (int iS = 0; iS < nSamples; ++iS)
        {
            // noDR (filled)
            TEfficiency *effNoDR = (TEfficiency *)fin->Get(
                Form("eff_id_noDR_%s", sampleSets[iS]));
            if (effNoDR)
            {
                effNoDR->SetMarkerColor(sampleColors[iS]);
                effNoDR->SetLineColor(sampleColors[iS]);
                effNoDR->SetMarkerStyle(fillMarkers[iS]);
                effNoDR->SetMarkerSize(1.0);
                effNoDR->SetLineWidth(2);
                effNoDR->SetLineStyle(1);
                effNoDR->Draw("same");
                legL->AddEntry(effNoDR,
                    Form("No #Delta#it{R}, %s", sampleSets[iS]), "pl");
            }

            // withDR (open)
            TEfficiency *effDR = (TEfficiency *)fin->Get(
                Form("eff_id_withDR_%s", sampleSets[iS]));
            if (effDR)
            {
                effDR->SetMarkerColor(sampleColors[iS]);
                effDR->SetLineColor(sampleColors[iS]);
                effDR->SetMarkerStyle(openMarkers[iS]);
                effDR->SetMarkerSize(1.0);
                effDR->SetLineWidth(2);
                effDR->SetLineStyle(2);
                effDR->Draw("same");
                legL->AddEntry(effDR,
                    Form("With #Delta#it{R}, %s", sampleSets[iS]), "pl");
            }
        }
        legL->Draw("same");

        myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
        myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
        myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 0);
        myText(xpos, ypos - 3 * dy, 1, "All sample sets", fontsize, 0);

        // ---------------------------------------------------------------
        // Right panel: BDT efficiency for "new clusters" (pass trkID, fail dR)
        // ---------------------------------------------------------------
        c->cd(2);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.14);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.06);

        TH1F *frameR = new TH1F("frame_id_right",
            ";#it{E}_{T}^{#gamma, truth} [GeV];BDT Efficiency (new clusters only)", 100, 8, 36);
        frameR->GetYaxis()->SetRangeUser(0.0, 1.05);
        frameR->GetXaxis()->SetTitleSize(0.055);
        frameR->GetYaxis()->SetTitleSize(0.050);
        frameR->GetXaxis()->SetLabelSize(0.045);
        frameR->GetYaxis()->SetLabelSize(0.045);
        frameR->GetYaxis()->SetTitleOffset(1.3);
        frameR->Draw("axis");

        TLegend *legR = new TLegend(0.42, 0.18, 0.92, 0.42);
        legStyle(legR, 0.22, 0.040);

        for (int iS = 0; iS < nSamples; ++iS)
        {
            TEfficiency *effNew = (TEfficiency *)fin->Get(
                Form("eff_id_newclusters_%s", sampleSets[iS]));
            if (effNew)
            {
                effNew->SetMarkerColor(sampleColors[iS]);
                effNew->SetLineColor(sampleColors[iS]);
                effNew->SetMarkerStyle(fillMarkers[iS]);
                effNew->SetMarkerSize(1.2);
                effNew->SetLineWidth(2);
                effNew->Draw("same");
                legR->AddEntry(effNew, sampleSets[iS], "pl");
            }
            else
            {
                std::cerr << "Missing eff_id_newclusters_" << sampleSets[iS] << std::endl;
            }
        }
        legR->Draw("same");

        myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
        myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
        myText(xpos, ypos - 2 * dy, 1, "Clusters: pass trkID, fail #Delta#it{R} < 0.1", 0.036, 0);

        c->SaveAs(Form("%s/photon_id_efficiency_focus.pdf", savePath.c_str()));
        delete c;
    }

    // ===================================================================
    // Canvas 6: Summary table (TLatex-based)
    // ===================================================================
    {
        TCanvas *c = new TCanvas("c_summary", "", 1000, 800);
        gPad->SetLeftMargin(0.02);
        gPad->SetRightMargin(0.02);
        gPad->SetTopMargin(0.06);
        gPad->SetBottomMargin(0.02);

        // Title
        TLatex title;
        title.SetNDC();
        title.SetTextSize(0.040);
        title.SetTextFont(62);
        title.DrawLatex(0.10, 0.94, "Efficiency Summary: with vs without #Delta#it{R} < 0.1 truth-matching cut");

        // Column headers
        TLatex tex;
        tex.SetNDC();
        tex.SetTextSize(0.030);
        tex.SetTextFont(42);

        // Table layout
        const float x0 = 0.04;  // row labels
        const float xCols[] = {0.22, 0.40, 0.58, 0.76}; // stage columns
        const float yHeader = 0.86;
        const float rowHeight = 0.042;

        // Header row: stage names
        tex.SetTextFont(62);
        tex.DrawLatex(x0, yHeader, "Sample / Matching");
        for (int iStage = 0; iStage < nStages; ++iStage)
            tex.DrawLatex(xCols[iStage], yHeader, stageLabels[iStage]);

        // Separator line
        TLine *sep = new TLine(0.03, yHeader - 0.015, 0.97, yHeader - 0.015);
        sep->SetNDC();
        sep->SetLineWidth(1);
        sep->Draw("same");

        tex.SetTextFont(42);
        int iRow = 0;

        for (int iS = 0; iS < nSamples; ++iS)
        {
            for (int iMatch = 0; iMatch < 2; ++iMatch)
            {
                const char *matching = (iMatch == 0) ? "noDR" : "withDR";
                const char *matchLabel = (iMatch == 0) ? "No #Delta#it{R}" : "#Delta#it{R} < 0.1";
                float yRow = yHeader - (iRow + 1) * rowHeight - 0.015;

                // Row label
                if (iMatch == 0)
                {
                    tex.SetTextFont(62);
                    tex.DrawLatex(x0, yRow, sampleSets[iS]);
                    tex.SetTextFont(42);
                }
                tex.DrawLatex(x0 + 0.07, yRow, matchLabel);

                // Efficiency values
                for (int iStage = 0; iStage < nStages; ++iStage)
                {
                    TEfficiency *eff = (TEfficiency *)fin->Get(
                        Form("eff_%s_%s_%s", stageNames[iStage], matching, sampleSets[iS]));
                    if (!eff) continue;

                    // Compute pT-averaged efficiency (all bins)
                    int nBins = eff->GetTotalHistogram()->GetNbinsX();
                    double sumPass = 0, sumTotal = 0;
                    for (int ib = 1; ib <= nBins; ++ib)
                    {
                        sumPass  += eff->GetPassedHistogram()->GetBinContent(ib);
                        sumTotal += eff->GetTotalHistogram()->GetBinContent(ib);
                    }
                    double avgEff = (sumTotal > 0) ? sumPass / sumTotal : 0;
                    double avgErr = (sumTotal > 0)
                        ? std::sqrt(avgEff * (1 - avgEff) / sumTotal) : 0;

                    tex.DrawLatex(xCols[iStage], yRow,
                        Form("%.3f #pm %.3f", avgEff, avgErr));
                }
                ++iRow;
            }

            // Ratio row
            {
                float yRow = yHeader - (iRow + 1) * rowHeight - 0.015;
                tex.SetTextColor(kBlue + 2);
                tex.DrawLatex(x0 + 0.07, yRow, "Ratio (no#Delta#it{R} / #Delta#it{R})");

                for (int iStage = 0; iStage < nStages; ++iStage)
                {
                    TEfficiency *effNoDR = (TEfficiency *)fin->Get(
                        Form("eff_%s_noDR_%s", stageNames[iStage], sampleSets[iS]));
                    TEfficiency *effDR = (TEfficiency *)fin->Get(
                        Form("eff_%s_withDR_%s", stageNames[iStage], sampleSets[iS]));
                    if (!effNoDR || !effDR) continue;

                    int nBins = effNoDR->GetTotalHistogram()->GetNbinsX();
                    double passN = 0, totN = 0, passD = 0, totD = 0;
                    for (int ib = 1; ib <= nBins; ++ib)
                    {
                        passN += effNoDR->GetPassedHistogram()->GetBinContent(ib);
                        totN  += effNoDR->GetTotalHistogram()->GetBinContent(ib);
                        passD += effDR->GetPassedHistogram()->GetBinContent(ib);
                        totD  += effDR->GetTotalHistogram()->GetBinContent(ib);
                    }
                    double eN = (totN > 0) ? passN / totN : 0;
                    double eD = (totD > 0) ? passD / totD : 0;
                    double ratio = (eD > 0) ? eN / eD : 0;

                    tex.DrawLatex(xCols[iStage], yRow, Form("%.4f", ratio));
                }
                tex.SetTextColor(kBlack);
                ++iRow;
            }

            // Blank row between sample sets
            ++iRow;
        }

        // Separator before new-clusters section
        float ySep2 = yHeader - iRow * rowHeight - 0.015;
        TLine *sep2 = new TLine(0.03, ySep2, 0.97, ySep2);
        sep2->SetNDC();
        sep2->SetLineWidth(1);
        sep2->Draw("same");

        // New clusters BDT efficiency
        tex.SetTextFont(62);
        float yNewHeader = ySep2 - rowHeight;
        tex.DrawLatex(x0, yNewHeader,
            "BDT pass rate for new clusters (pass trkID, fail #Delta#it{R})");
        tex.SetTextFont(42);

        for (int iS = 0; iS < nSamples; ++iS)
        {
            float yRow = yNewHeader - (iS + 1) * rowHeight;
            tex.DrawLatex(x0, yRow, sampleSets[iS]);

            TEfficiency *effNew = (TEfficiency *)fin->Get(
                Form("eff_id_newclusters_%s", sampleSets[iS]));
            if (!effNew) continue;

            int nBins = effNew->GetTotalHistogram()->GetNbinsX();
            double sumPass = 0, sumTotal = 0;
            for (int ib = 1; ib <= nBins; ++ib)
            {
                sumPass  += effNew->GetPassedHistogram()->GetBinContent(ib);
                sumTotal += effNew->GetTotalHistogram()->GetBinContent(ib);
            }
            double avgEff = (sumTotal > 0) ? sumPass / sumTotal : 0;
            double avgErr = (sumTotal > 0)
                ? std::sqrt(avgEff * (1 - avgEff) / sumTotal) : 0;

            tex.DrawLatex(xCols[0], yRow,
                Form("Integrated: %.3f #pm %.3f", avgEff, avgErr));
        }

        // sPHENIX watermark
        myText(0.70, 0.94, 1, strleg1.c_str(), 0.035, 0);

        c->SaveAs(Form("%s/efficiency_summary_table.pdf", savePath.c_str()));
        delete c;
    }

    fin->Close();
    std::cout << "All plots saved to " << savePath << std::endl;
}
