#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <string>

#include "plotcommon.h"

void plot_showershapes()
{
    init_plot();
    string savePath = "../PPG12-analysis-note/Figures/showershapes/";

    std::vector<double> eta_bins = {-0.7, 0.7};
    // std::vector<double> pT_bin_edges = {8, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50};
    std::vector<double> pT_bin_edges = {8, 10, 12, 15, 20, 25};
    // Number of bins is one less than the number of edges
    int nEtaBins = eta_bins.size() - 1;
    int nPtBins = pT_bin_edges.size() - 1;
    int nCuts = 2; // e.g. 0 or 1 (two cut options)

    //------------------------------------------------------------------------------
    // 2) Open your three ROOT files: data, signal, background
    //    (Update the filenames as needed.)
    //------------------------------------------------------------------------------
    TFile *f_data = TFile::Open("/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histoshower_shape_.root", "READ");
    TFile *f_sig  = TFile::Open("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_signal.root", "READ");
    TFile *f_bkg  = TFile::Open("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_jet.root", "READ");
    
    if (!f_data || f_data->IsZombie())
    {
        std::cerr << "Error: Could not open data file!" << std::endl;
        return;
    }
    if (!f_sig || f_sig->IsZombie())
    {
        std::cerr << "Error: Could not open signal file!" << std::endl;
        return;
    }
    if (!f_bkg || f_bkg->IsZombie())
    {
        std::cerr << "Error: Could not open background file!" << std::endl;
        return;
    }

    //------------------------------------------------------------------------------
    // 3) List of base histogram names (the keys in your map).
    //    They should match exactly what you used to create them (e.g. "h2d_prob").
    //------------------------------------------------------------------------------
    std::vector<std::string> histNames = {
        "h2d_prob",
        "h2d_CNN_prob",
        "h2d_e17_to_e77",
        "h2d_e37_to_e77",
        "h2d_e32_to_e35",
        "h2d_e33_to_e35",
        "h2d_e11_to_e33",
        "h2d_e11_to_E",
        "h2d_e33_to_E",
        "h2d_hcalet33_to_ettot",
        "h2d_ihcalet33_to_ettot",
        "h2d_ohcalet33_to_ettot",
        "h2d_hcalet22_to_ettot",
        "h2d_ihcalet22_to_ettot",
        "h2d_ohcalet22_to_ettot",
        //"h2d_detamax",
        //"h2d_dphimax",
        //"h2d_e1",
        //"h2d_e2",
        //"h2d_e3",
        //"h2d_e4",
        "h2d_et1",
        "h2d_et2",
        "h2d_et3",
        "h2d_et4",
        "h2d_weta",
        "h2d_wphi",
        "h2d_w32",
        "h2d_w52",
        "h2d_w72",
        "h2d_wr",
        "h2d_wrr",
        //"h2d_weta_cog",
        //"h2d_wphi_cog",
        "h2d_weta_cogx",
        "h2d_wphi_cogx"};

    //------------------------------------------------------------------------------
    // Helper lambda to scale a TH1D to unit area
    //------------------------------------------------------------------------------
    auto scaleToUnit = [](TH1D *h)
    {
        if (!h)
            return;
        double integral = h->Integral();
        if (integral > 1e-12)
            h->Scale(1.0 / integral);
    };

    //------------------------------------------------------------------------------
    // 4) Loop over each histogram base name, then over cuts, eta-bins, and pT-bins
    //------------------------------------------------------------------------------
    for (const auto &hbase : histNames)
    {
        //set

        for (int icut = 0; icut < nCuts; ++icut)
        {
            for (int ieta = 0; ieta < nEtaBins; ++ieta)
            {
                for (int ipt = 0; ipt < nPtBins; ++ipt)
                {
                    // x axis name is the base name after first 4 characters
                    TString xaxisname = hbase.substr(4, hbase.size() - 4);

                    float xaxismax = 1.0;
                    int nrebin = 4;

                    //if the first character of xaxisname is 'w', then the xaxismax is 2.0
                    if (xaxisname[0] == 'w')
                    {
                        xaxismax = 2.0;
                    }
                    //if first is h , then xaxismax is 0.4
                    if (xaxisname[0] == 'h')
                    {
                        xaxismax = 0.4;
                        nrebin = 1;
                    }



                    TString histNameFull = Form(
                        "%s_eta%d_pt%d_cut%d",
                        hbase.c_str(),
                        ieta,
                        ipt,
                        icut);

                    TString histNamesave = Form(
                        "%s_eta%d_pt%d_cut%d",
                        xaxisname.Data(),
                        ieta,
                        ipt,
                        icut);

                    // Retrieve TH2F from each file
                    TH2F *h2_data = dynamic_cast<TH2F *>(f_data->Get(histNameFull));
                    TH2F *h2_sig = dynamic_cast<TH2F *>(f_sig->Get(histNameFull));
                    TH2F *h2_bkg = dynamic_cast<TH2F *>(f_bkg->Get(histNameFull));
                    h2_bkg->RebinX(nrebin);
                    h2_sig->RebinX(nrebin);
                    h2_data->RebinX(nrebin);

                    h2_bkg->GetXaxis()->SetRangeUser(0, xaxismax);
                    h2_sig->GetXaxis()->SetRangeUser(0, xaxismax);
                    h2_data->GetXaxis()->SetRangeUser(0, xaxismax);

                    // Check if they exist (skip if not found)
                    if (!h2_data || !h2_sig || !h2_bkg)
                    {
                        // Possibly comment this out if you expect many missing
                        std::cerr << "Warning: Could not retrieve "
                                  << histNameFull << " from one or more files!\n";
                        continue;
                    }

                    TH1D *proj_data = h2_data->ProjectionX(
                        Form("%s_px_data", histNameFull.Data()));
                    TH1D *proj_sig = h2_sig->ProjectionX(
                        Form("%s_px_sig", histNameFull.Data()));
                    TH1D *proj_bkg = h2_bkg->ProjectionX(
                        Form("%s_px_bkg", histNameFull.Data()));

                    scaleToUnit(proj_data);
                    scaleToUnit(proj_sig);
                    scaleToUnit(proj_bkg);

                    float maxy = std::max({proj_data->GetMaximum(), proj_sig->GetMaximum(), proj_bkg->GetMaximum()});

                    TCanvas *c_proj = new TCanvas(
                        Form("c_proj_%s", histNameFull.Data()),
                        Form("ProjectionX - %s", histNameFull.Data()),
                        600, 600);
                    c_proj->cd();

                    // Set colors
                    proj_data->SetLineColor(kBlack);
                    proj_data->SetMarkerColor(kBlack);
                    proj_sig->SetLineColor(kRed);
                    proj_sig->SetMarkerColor(kRed);
                    proj_bkg->SetLineColor(kBlue);
                    proj_bkg->SetMarkerColor(kBlue);

                    // Draw them
                    proj_sig->SetYTitle("normalized counts");
                    proj_sig->SetXTitle(xaxisname.Data());
                    proj_sig->GetYaxis()->SetRangeUser(0, maxy * 1.2);
                    proj_sig->Draw("HIST");
                    proj_data->Draw("HIST SAME");

                    proj_bkg->Draw("HIST SAME");

                    myText(0.50, 0.90, 1, strleg1.c_str(), 0.04);
                    myText(0.50, 0.85, 1, strleg2.c_str(), 0.04);
                    myText(0.50, 0.80, 1, strleg3.c_str(), 0.04);
                    float pTlow = pT_bin_edges[ipt];
                    float pThigh = pT_bin_edges[ipt + 1];
                    float etalow = eta_bins[ieta];
                    float etahigh = eta_bins[ieta + 1];
                    std::string bgcut = (icut == 0) ? "w/o bg cut" : "w/ bg cut";
                    myText(0.35, 0.75, 1, Form("%.1f<#eta< %.1f, %.0f<p_{T}<%.0fGeV,%s", etalow, etahigh, pTlow, pThigh, bgcut.c_str()), 0.04);

                    myMarkerLineText(0.55, 0.70, 0, kBlack, 0, kBlack, 1,
                                     "Data", 0.05, true);
                    myMarkerLineText(0.55, 0.65, 0, kRed, 0, kRed, 1,
                                     "Signal", 0.05, true);
                    myMarkerLineText(0.55, 0.60, 0, kBlue, 0, kBlue, 1,
                                     "Background", 0.05, true);

                    c_proj->SaveAs(Form("%s/dis_%s.pdf", savePath.c_str(), histNamesave.Data()));

                    TProfile *pfx_data = h2_data->ProfileX(
                        Form("%s_pfx_data", histNameFull.Data()),
                        1,  // first y-bin
                        -1, // last y-bin (use all bins)
                        "s" // 's' option to store RMS in the bin error
                    );

                    TProfile *pfx_sig = h2_sig->ProfileX(
                        Form("%s_pfx_sig", histNameFull.Data()),
                        1, -1, "s");

                    TProfile *pfx_bkg = h2_bkg->ProfileX(
                        Form("%s_pfx_bkg", histNameFull.Data()),
                        1, -1, "");
                    TCanvas *c_prof = new TCanvas(
                        Form("c_prof_%s", histNameFull.Data()),
                        Form("ProfileX - %s", histNameFull.Data()),
                        600, 600);
                    c_prof->cd();

                    // Style
                    if (pfx_data)
                    {
                        pfx_data->SetLineColor(kBlack);
                        pfx_data->SetMarkerColor(kBlack);
                    }
                    if (pfx_sig)
                    {
                        pfx_sig->SetLineColor(kRed);
                        pfx_sig->SetMarkerColor(kRed);
                    }
                    if (pfx_bkg)
                    {
                        pfx_bkg->SetLineColor(kBlue);
                        pfx_bkg->SetMarkerColor(kBlue);
                    }

                    // Draw
                    if (pfx_data)
                        pfx_data->Draw("E");
                    if (pfx_sig)
                        pfx_sig->Draw("E SAME");
                    if (pfx_bkg){
                        pfx_bkg->Draw("hist SAME");
                        pfx_bkg->Draw("ex0 SAME");
                        pfx_bkg->SetMarkerSize(0);
                    }

                    myText(0.50, 0.90, 1, strleg1.c_str(), 0.04);
                    myText(0.50, 0.85, 1, strleg2.c_str(), 0.04);
                    myText(0.50, 0.80, 1, strleg3.c_str(), 0.04);

                    myMarkerLineText(0.55, 0.75, 0, kBlack, 0, kBlack, 1,
                                     "Data", 0.05, true);
                    myMarkerLineText(0.55, 0.70, 0, kRed, 0, kRed, 1,
                                     "Signal", 0.05, true);
                    myMarkerLineText(0.55, 0.65, 0, kBlue, 0, kBlue, 1,
                                     "Background", 0.05, true);

                    TCanvas *c_bkg = new TCanvas(
                        Form("c_bkg_%s", histNameFull.Data()),
                        Form("Bkg only - %s", histNameFull.Data()),
                        600, 600);
                    c_bkg->cd();

                    TH1D *proj_bkg_clone = (TH1D *)pfx_bkg->Clone(
                        Form("%s_px_bkgOnly", histNameFull.Data()));
                    proj_bkg_clone->SetLineColor(kBlue);

                    proj_bkg_clone->SetYTitle("<#it{E}_{T}^{iso} [GeV]>");
                    proj_bkg_clone->SetXTitle(xaxisname.Data());
                    float max = proj_bkg_clone->GetMaximum();
                    proj_bkg_clone->GetYaxis()->SetRangeUser(0,max*1.3);

                    proj_bkg_clone->Draw("HIST");
                    proj_bkg_clone->Draw("same ex0");
                    proj_bkg_clone->SetMarkerSize(0);

                    float corr = h2_bkg->GetCorrelationFactor();

                    myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                    myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                    myText(0.20, 0.80, 1, strleg3.c_str(), 0.04);
                    myText(0.55, 0.90, 1, Form("%.0f<p_{T}<%.0fGeV,%s", pTlow, pThigh, bgcut.c_str()), 0.04);
                    myText(0.55, 0.85, 1, "Background Only MC", 0.04);
                    myText(0.55, 0.80, 1, Form("Correlation: %.3f", corr), 0.04);

                    // Optionally save
                    c_bkg->SaveAs(Form("%s/pfx_%s.pdf", savePath.c_str(), histNamesave.Data()));

                } // end ipt
            } // end ieta
        } // end icut
    } // end loop over histNames

    //------------------------------------------------------------------------------
    // 5) Close files (if you want).
    //------------------------------------------------------------------------------
    f_data->Close();
    f_sig->Close();
    f_bkg->Close();
}
