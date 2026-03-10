// FindBDTCut.C
// Reads 2D histograms from BDTScoreVsET output and computes the BDT score
// threshold corresponding to 60/70/80/90% signal efficiency as a function
// of cluster ET.  Results are fit with a linear function and plotted.
//
// Usage:
//   root -l -q 'FindBDTCut.C("results/BDTScoreVsET_bdt_nom.root","base_v3E",0)'

#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include "draw.C"

// ---------------------------------------------------------------------------
// findBDTCutoff
// Accumulates from the HIGH BDT end downward until the cumulative fraction
// reaches `eff`, then returns the lower edge of that BDT bin.
// ---------------------------------------------------------------------------
float findBDTCutoff(TH2D *h_signal, float eff, float pTlow, float pThigh)
{
    // X (ET) bins spanning [pTlow, pThigh]
    int xbinlow  = TMath::Max(1, h_signal->GetXaxis()->FindBin(pTlow));
    int xbinhigh = TMath::Min(h_signal->GetNbinsX(), h_signal->GetXaxis()->FindBin(pThigh));

    // Total entries in ET range (all BDT bins)
    float total = 0.0;
    for (int i = xbinlow; i <= xbinhigh; ++i)
        for (int j = 1; j <= h_signal->GetNbinsY(); ++j)
            total += h_signal->GetBinContent(i, j);

    if (total <= 0.0) return 0.0;

    float target     = total * eff;
    float currentsum = 0.0;

    // Scan from highest BDT bin downward
    for (int j = h_signal->GetNbinsY(); j >= 1; --j)
    {
        for (int i = xbinlow; i <= xbinhigh; ++i)
            currentsum += h_signal->GetBinContent(i, j);

        if (currentsum >= target)
            return h_signal->GetYaxis()->GetBinLowEdge(j);
    }

    // Fallback: lowest BDT bin
    return h_signal->GetYaxis()->GetBinLowEdge(1);
}

// ---------------------------------------------------------------------------
// FindBDTCut
// ---------------------------------------------------------------------------
void FindBDTCut(const std::string &inputfile  = "results/BDTScoreVsET_bdt_nom.root",
                const std::string &model_name = "base_v3E",
                int ieta = 0)
{
    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TFile *fin = TFile::Open(inputfile.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << inputfile << std::endl;
        return;
    }

    TH2D *h_signal = (TH2D*)fin->Get(
        Form("h_signal_bdt_%s_%d", model_name.c_str(), ieta));
    if (!h_signal)
    {
        std::cerr << "ERROR: h_signal_bdt_" << model_name << "_" << ieta
                  << " not found in " << inputfile << std::endl;
        fin->Close();
        return;
    }

    TH2D *h_incl = (TH2D*)fin->Get(
        Form("h_incl_bdt_%s_%d", model_name.c_str(), ieta));
    if (!h_incl)
    {
        std::cerr << "ERROR: h_incl_bdt_" << model_name << "_" << ieta
                  << " not found in " << inputfile << std::endl;
        fin->Close();
        return;
    }

    std::vector<int>         colors  = {kPink + 8, kSpring - 7, kAzure - 3, kViolet + 3};
    std::vector<int>         markers = {20, 21, 22, 23};
    std::vector<std::string> legends = {"90% Efficiency", "80% Efficiency",
                                        "70% Efficiency", "60% Efficiency"};

    // -----------------------------------------------------------------------
    // Helper lambda: fill bins, fit, draw, and print for one TH2D source
    // -----------------------------------------------------------------------
    auto runOnSource = [&](TH2D *h_src, const std::string &label, const std::string &outpdf)
    {
        // 13 equal-width bins from 10 to 36 GeV (same as FindETCut.C)
        const int nbins = 13;
        TH1F *h_eff90 = new TH1F(Form("h_eff90_%s", label.c_str()),
            "90% Eff BDT Cutoff; Cluster E_{T} [GeV]; BDT Score", nbins, 10, 36);
        TH1F *h_eff80 = new TH1F(Form("h_eff80_%s", label.c_str()),
            "80% Eff BDT Cutoff; Cluster E_{T} [GeV]; BDT Score", nbins, 10, 36);
        TH1F *h_eff70 = new TH1F(Form("h_eff70_%s", label.c_str()),
            "70% Eff BDT Cutoff; Cluster E_{T} [GeV]; BDT Score", nbins, 10, 36);
        TH1F *h_eff60 = new TH1F(Form("h_eff60_%s", label.c_str()),
            "60% Eff BDT Cutoff; Cluster E_{T} [GeV]; BDT Score", nbins, 10, 36);

        float ybinwidth = h_src->GetYaxis()->GetBinWidth(1);

        for (int i = 1; i <= nbins; ++i)
        {
            float pTlow  = h_eff90->GetXaxis()->GetBinLowEdge(i);
            float pThigh = h_eff90->GetXaxis()->GetBinUpEdge(i);

            h_eff90->SetBinContent(i, findBDTCutoff(h_src, 0.9, pTlow, pThigh));
            h_eff90->SetBinError  (i, ybinwidth / 2.0);

            h_eff80->SetBinContent(i, findBDTCutoff(h_src, 0.8, pTlow, pThigh));
            h_eff80->SetBinError  (i, ybinwidth / 2.0);

            h_eff70->SetBinContent(i, findBDTCutoff(h_src, 0.7, pTlow, pThigh));
            h_eff70->SetBinError  (i, ybinwidth / 2.0);

            h_eff60->SetBinContent(i, findBDTCutoff(h_src, 0.6, pTlow, pThigh));
            h_eff60->SetBinError  (i, ybinwidth / 2.0);
        }

        // Y-axis range from the threshold values (with 10% padding)
        float ymin = 1.0, ymax = 0.0;
        for (TH1F *h : {h_eff90, h_eff80, h_eff70, h_eff60})
        {
            for (int i = 1; i <= nbins; ++i)
            {
                float v = h->GetBinContent(i);
                if (v > 0 && v < ymin) ymin = v;
                if (v > ymax)          ymax = v;
            }
        }
        float ypad = (ymax - ymin) * 0.15;
        ymin = TMath::Max(0.0f, ymin - ypad);
        ymax = ymax + ypad;

        // Linear fits
        TF1 *fit90 = new TF1(Form("fit90_%s", label.c_str()), "pol1", 8, 40);
        fit90->SetLineColor(colors[0]);
        h_eff90->Fit(fit90, "RQ");

        TF1 *fit80 = new TF1(Form("fit80_%s", label.c_str()), "pol1", 8, 40);
        fit80->SetLineColor(colors[1]);
        h_eff80->Fit(fit80, "RQ");

        TF1 *fit70 = new TF1(Form("fit70_%s", label.c_str()), "pol1", 8, 40);
        fit70->SetLineColor(colors[2]);
        h_eff70->Fit(fit70, "RQ");

        TF1 *fit60 = new TF1(Form("fit60_%s", label.c_str()), "pol1", 8, 40);
        fit60->SetLineColor(colors[3]);
        h_eff60->Fit(fit60, "RQ");

        // Draw
        std::vector<std::string> fitTexts = {
            Form("Model: %s (%s)", model_name.c_str(), label.c_str()),
            Form("90%%: BDT = %.3f + %.4f #times E_{T}", fit90->GetParameter(0), fit90->GetParameter(1)),
            Form("80%%: BDT = %.3f + %.4f #times E_{T}", fit80->GetParameter(0), fit80->GetParameter(1)),
            Form("70%%: BDT = %.3f + %.4f #times E_{T}", fit70->GetParameter(0), fit70->GetParameter(1)),
            Form("60%%: BDT = %.3f + %.4f #times E_{T}", fit60->GetParameter(0), fit60->GetParameter(1)),
        };

        draw_1D_multiple_plot(
            {h_eff90, h_eff80, h_eff70, h_eff60}, colors, markers,
            false, 1, false,                                          // no rebin/normalise
            true,  10, 36, false,                                     // x range
            true,  ymin, ymax, false,                                 // y range
            true,  "Cluster E_{T} [GeV]", "BDT Score Cutoff",        // axis titles
            false, "",                                                 // not data
            true,  fitTexts, 0.15, 0.85, 0.030,                      // text annotations
            true,  legends,  0.55, 0.80, 0.030,                      // legend
            outpdf                                                     // output
        );

        // Print fit results
        std::cout << "\n=== BDT cut fit results: model=" << model_name
                  << ", ieta=" << ieta << ", source=" << label << " ===" << std::endl;
        std::cout << "90%: BDT_cut = " << fit90->GetParameter(0)
                  << " + " << fit90->GetParameter(1) << " * ET" << std::endl;
        std::cout << "80%: BDT_cut = " << fit80->GetParameter(0)
                  << " + " << fit80->GetParameter(1) << " * ET" << std::endl;
        std::cout << "70%: BDT_cut = " << fit70->GetParameter(0)
                  << " + " << fit70->GetParameter(1) << " * ET" << std::endl;
        std::cout << "60%: BDT_cut = " << fit60->GetParameter(0)
                  << " + " << fit60->GetParameter(1) << " * ET" << std::endl;
    };

    runOnSource(h_signal, "signal",
        Form("BDTCut_FitResults_%s_signal.pdf", model_name.c_str()));
    runOnSource(h_incl,   "incl",
        Form("BDTCut_FitResults_%s_incl.pdf",   model_name.c_str()));
}
