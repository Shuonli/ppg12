#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TLine.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <yaml-cpp/yaml.h>
#include <TSystem.h>

#include "plotcommon.h"

void plot_showershapes_selections()
{
    init_plot();
    string savePath = "../PPG12-analysis-note/Figures/showershapes_selections/";

    // Create output directory if it doesn't exist
    gSystem->Exec(Form("mkdir -p %s", savePath.c_str()));

    // Load config for BDT bins
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node config = YAML::LoadFile("../efficiencytool/config_showershape.yaml");

    // Optional plotting toggles
    const bool makeDataMinusNPB = true;

    // BDT bins with default values if not specified in config
    std::vector<double> bdt_bins_default = {0.0, 0.3, 0.7, 1.0};
    std::vector<double> bdt_bins = config["analysis"]["bdt_bins"].as<std::vector<double>>(bdt_bins_default);
    std::vector<double> pT_bin_edges = config["analysis"]["pT_bins"].as<std::vector<double>>();
    std::vector<double> eta_bins = {-0.7, 0.7};

    int nBdtBins = bdt_bins.size() - 1;
    int nPtBins = pT_bin_edges.size() - 1;
    int nEtaBins = eta_bins.size() - 1;

    std::cout << "Using " << nBdtBins << " BDT bins from config file" << std::endl;
    std::cout << "Using " << nPtBins << " pT bins from config file" << std::endl;

    // NPB score ranges for Part 4
    std::vector<std::pair<double, double>> npb_score_ranges = {
        {0.0, 0.2},   // Very low NPB score
        {0.2, 0.4},   // Low NPB score
        {0.4, 0.6},   // Medium NPB score (spans cut threshold 0.5)
        {0.6, 0.8},   // High NPB score
        {0.8, 1.0}    // Very high NPB score
    };
    int nNpbBins = npb_score_ranges.size();

    // Color scheme for NPB score progression (low → high)
    std::vector<Color_t> npb_time_colors = {kBlue, kGreen+2, kOrange+7, kMagenta+1, kRed};

    // Open files (keep these paths in one place; allows swapping inclusive vs background-only inputs)
    const std::string dataFile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histoshower_shape_.root";
    const std::string sigFile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_signal.root";

    // Background MC: you may want BOTH
    // - inclusive: produced with doinclusive=true (default in ShowerShapeCheck)
    // - bkgOnly: produced with doinclusive=false (truth-matched background-only for jet samples)
    const std::string bkgInclusiveFile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_jet_inclusive.root";
    const std::string bkgOnlyFile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_jet.root"; // set this to a background-only file if you have it

    TFile *f_data = TFile::Open(dataFile.c_str(), "READ");
    TFile *f_sig = TFile::Open(sigFile.c_str(), "READ");
    TFile *f_bkg_inc = TFile::Open(bkgInclusiveFile.c_str(), "READ");

    bool hasBkgOnly = false;
    TFile *f_bkg_only = nullptr;
    if (!bkgOnlyFile.empty())
    {
        f_bkg_only = TFile::Open(bkgOnlyFile.c_str(), "READ");
        hasBkgOnly = (f_bkg_only && !f_bkg_only->IsZombie());
    }
    if (!hasBkgOnly)
    {
        f_bkg_only = f_bkg_inc; // fallback
    }

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
    if (!f_bkg_inc || f_bkg_inc->IsZombie())
    {
        std::cerr << "Error: Could not open inclusive background file!" << std::endl;
        return;
    }


    // Histogram base names for standard data/MC comparison (also used in Part 5)
    std::vector<std::string> histNames = {
        "h2d_weta_cogx", "h2d_wphi_cogx", "h2d_et1", "h2d_et2", "h2d_et3", "h2d_et4",
        "h2d_e11_to_e33", "h2d_e17_to_e77", "h2d_e32_to_e35", "h2d_bdt", "h2d_npb_score"
    };

    int nCuts = 4; // cut0, cut1, cut2 (tight), cut3 (non-tight)

    // Common color palette for pT bins (used in Parts 6, 7)
    std::vector<Color_t> pt_colors = {kBlue, kRed, kGreen+2, kMagenta+1, kOrange+7,
                                      kCyan+2, kViolet+1, kTeal+2, kPink+1, kAzure+2};

    // Common file/sample definitions (used in Parts 4, 6, 7)
    std::vector<std::tuple<std::string, TFile*, Color_t, std::string>> all_samples = {
        {"data", f_data, kBlack, "Data"},
        {"signal", f_sig, kRed, "Signal MC"},
        {"bkg_inclusive", f_bkg_inc, kBlue, "Inclusive MC"}
    };

    // Helper: get axis settings for a variable
    auto getAxisSettings = [](const TString& xaxisname, float& xmin, float& xmax, int& nrebin) {
        xmin = 0.0; xmax = 1.0; nrebin = 4;
        if (xaxisname[0] == 'w') { xmax = 2.0; }
        if (xaxisname[0] == 'h') { xmax = 0.3; nrebin = 1; }
        if (xaxisname == "et4") { xmax = 0.3; nrebin = 1; }
        if (xaxisname == "e32_to_e35") { xmin = 0.4; xmax = 1.0; nrebin = 1; }
        if (xaxisname == "et1") { xmin = 0.3; xmax = 1.0; nrebin = 1; }
        if (xaxisname == "bdt" || xaxisname == "npb_score") { xmin = 0.0; xmax = 1.0; nrebin = 2; }
    };

    // Helper lambda to scale to unit area
    auto scaleToUnit = [](TH1D *h)
    {
        if (!h)
            return;
        double integral = h->Integral();
        if (integral > 1e-12)
            h->Scale(1.0 / integral);
    };

    // Helper lambda to overlay vertical error bars for MC histograms:
    // draw a copy of the same histogram with error bars only (no markers/fill).
    auto overlayVerticalErrors = [](TH1 *h)
    {
        if (!h)
            return;
        TH1 *herr = dynamic_cast<TH1 *>(h->Clone(Form("%s_err", h->GetName())));
        if (!herr)
            return;
        herr->SetDirectory(nullptr);
        herr->SetFillStyle(0);
        herr->SetMarkerStyle(1);
        herr->SetMarkerSize(0);
        herr->SetLineColor(h->GetLineColor());
        herr->SetLineWidth(h->GetLineWidth());
        herr->SetMarkerColor(h->GetLineColor());
        // allow pad/canvas to clean this up
        herr->SetBit(kCanDelete);
        herr->Draw("E0 SAME");
    };

    //==============================================================================
    // Precompute NPB normalization factors per (eta, pT) from DATA:
    //
    // Goal (for cut0 overlays where we draw unit-area shapes):
    //   Draw NPB with the *fractional normalization* estimated from a tail region:
    //     weta_cogx > 1
    //
    // Method:
    //   Let D(x) be the cut0 data distribution (raw counts), N(x) be the cut4 NPB-tagged distribution (raw counts).
    //   Estimate NPB fraction in data as:
    //
    //     f_NPB = ( ∫_{weta>1} D / ∫ D ) / ( ∫_{weta>1} N / ∫ N )
    //
    //   Then in the cut0 plot (which normalizes shapes to unit area), we draw:
    //     f_NPB * N_norm(x)
    //
    // This makes the NPB curve match the data tail fraction assuming the tail is NPB-dominated.
    //==============================================================================
    std::vector<std::vector<double>> npbScale(nEtaBins, std::vector<double>(nPtBins, 0.0));
    for (int ieta = 0; ieta < nEtaBins; ++ieta)
    {
        for (int ipt = 0; ipt < nPtBins; ++ipt)
        {
            TString h_weta_all_name = Form("h2d_weta_cogx_eta%d_pt%d_cut0", ieta, ipt);
            TString h_weta_npb_name = Form("h2d_weta_cogx_eta%d_pt%d_cut4", ieta, ipt);
            TH2F *h2_weta_all = dynamic_cast<TH2F *>(f_data->Get(h_weta_all_name));
            TH2F *h2_weta_npb = dynamic_cast<TH2F *>(f_data->Get(h_weta_npb_name));
            if (!h2_weta_all || !h2_weta_npb)
            {
                npbScale[ieta][ipt] = 0.0;
                std::cerr << "[NPB scale] Missing hist(s): "
                          << h_weta_all_name << " or " << h_weta_npb_name << std::endl;
                continue;
            }

            TH1D *p_all = h2_weta_all->ProjectionX(Form("%s_px_all", h_weta_all_name.Data()));
            TH1D *p_npb = h2_weta_npb->ProjectionX(Form("%s_px_npb", h_weta_npb_name.Data()));
            if (!p_all || !p_npb)
            {
                npbScale[ieta][ipt] = 0.0;
                std::cerr << "[NPB scale] Failed ProjectionX for: "
                          << h_weta_all_name << " or " << h_weta_npb_name << std::endl;
                continue;
            }

            // Match Part 0 settings for weta-like vars
            p_all->Rebin(4);
            p_npb->Rebin(4);
            p_all->GetXaxis()->SetRangeUser(0.0, 2.0);
            p_npb->GetXaxis()->SetRangeUser(0.0, 2.0);

            const int binStart = std::max(1, p_all->FindBin(1.400001));
            const double all_tot = p_all->Integral(1, p_all->GetNbinsX());
            const double npb_tot = p_npb->Integral(1, p_npb->GetNbinsX());
            const double all_tail = p_all->Integral(binStart, p_all->GetNbinsX());
            const double npb_tail = p_npb->Integral(binStart, p_npb->GetNbinsX());

            if (all_tot <= 0 || npb_tot <= 0 || all_tail < 0 || npb_tail <= 0)
            {
                npbScale[ieta][ipt] = 0.0;
                std::cerr << "[NPB scale] Bad integrals for (eta,pt)=(" << ieta << "," << ipt << ") "
                          << " all_tot=" << all_tot << " all_tail=" << all_tail
                          << " npb_tot=" << npb_tot << " npb_tail=" << npb_tail << std::endl;
            }
            else
            {
                const double frac_all = all_tail / all_tot;
                const double frac_npb = npb_tail / npb_tot;
                const double f_npb = (frac_npb > 1e-12) ? (frac_all / frac_npb) : 0.0;

                // Guard against unphysical scaling (NPB fraction should be <= 1; allow small tolerance for stats)
                if (!(f_npb >= 0.0) || f_npb > 1.2)
                {
                    std::cerr << "[NPB scale] Unphysical f_NPB=" << f_npb
                              << " for (eta,pt)=(" << ieta << "," << ipt << ") "
                              << " frac_all=" << frac_all << " frac_npb=" << frac_npb
                              << " [weta>1]" << std::endl;
                    npbScale[ieta][ipt] = 0.0;
                }
                else
                {
                    npbScale[ieta][ipt] = f_npb;
                }
            }

            delete p_all;
            delete p_npb;
        }
    }

    //==============================================================================
    // Part 0: Standard Data/MC Comparison Plots (from original plot_showershapes.C)
    //==============================================================================
    std::cout << "Creating standard data/MC comparison plots..." << std::endl;

    for (const auto &hbase : histNames)
    {
        for (int icut = 0; icut < nCuts; ++icut)
        {
            for (int ieta = 0; ieta < nEtaBins; ++ieta)
            {
                for (int ipt = 0; ipt < nPtBins; ++ipt)
                {
                    // x axis name is the base name after first 4 characters
                    TString xaxisname = hbase.substr(4, hbase.size() - 4);
                    float xaxismin, xaxismax; int nrebin;
                    getAxisSettings(xaxisname, xaxismin, xaxismax, nrebin);

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
                    TH2F *h2_bkg = dynamic_cast<TH2F *>(f_bkg_inc->Get(histNameFull));

                    // Check if they exist
                    if (!h2_data || !h2_sig || !h2_bkg)
                    {
                        std::cerr << "Warning: Could not retrieve histograms for "
                                  << histNameFull << " from one or more files!\n";
                        continue;
                    }

                    h2_bkg->RebinX(nrebin);
                    h2_sig->RebinX(nrebin);
                    h2_data->RebinX(nrebin);

                    h2_bkg->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
                    h2_sig->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
                    h2_data->GetXaxis()->SetRangeUser(xaxismin, xaxismax);

                    // Create projection plots
                    TH1D *proj_data = h2_data->ProjectionX(
                        Form("%s_px_data", histNameFull.Data()));
                    TH1D *proj_sig = h2_sig->ProjectionX(
                        Form("%s_px_sig", histNameFull.Data()));
                    TH1D *proj_bkg = h2_bkg->ProjectionX(
                        Form("%s_px_bkg", histNameFull.Data()));

                    // Optional: DATA NPB overlay for cut0 (uses cut4 written by ShowerShapeCheck)
                    TH1D *proj_data_npb = nullptr;
                    if (icut == 0)
                    {
                        TString histNameNPBFull = Form("%s_eta%d_pt%d_cut4", hbase.c_str(), ieta, ipt);
                        TH2F *h2_data_npb = dynamic_cast<TH2F *>(f_data->Get(histNameNPBFull));
                        if (h2_data_npb)
                        {
                            // Match the same X rebinning and X range as the cut0 plots
                            h2_data_npb->RebinX(nrebin);
                            h2_data_npb->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
                            proj_data_npb = h2_data_npb->ProjectionX(Form("%s_px_data_npb", histNameNPBFull.Data()));
                        }
                        if(!proj_data_npb)
                        {
                            std::cerr << "Warning: Could not retrieve NPB histogram for " << histNameNPBFull << std::endl;
                        }
                    }



                    scaleToUnit(proj_data);
                    scaleToUnit(proj_sig);
                    scaleToUnit(proj_bkg);

                    // Set consistent styles early (used by both the cut0 overlay and the NPB-subtracted overlay)
                    proj_data->SetLineColor(kBlack);
                    proj_data->SetMarkerColor(kBlack);
                    proj_data->SetMarkerStyle(20);
                    proj_data->SetMarkerSize(1.0);
                    proj_data->SetStats(0);
                    proj_sig->SetLineColor(kRed);
                    proj_sig->SetMarkerColor(kRed);
                    proj_sig->SetLineWidth(2);
                    proj_sig->SetStats(0);
                    proj_bkg->SetLineColor(kBlue);
                    proj_bkg->SetMarkerColor(kBlue);
                    proj_bkg->SetLineWidth(2);
                    proj_bkg->SetStats(0);

                    // Prepare NPB overlay (normalize by weta>1 tail matching)
                    if (proj_data_npb)
                    {
                        proj_data_npb->Sumw2();
                        scaleToUnit(proj_data_npb);
                        const double s = npbScale[ieta][ipt];
                        if (s > 0.0)
                        {
                            proj_data_npb->Scale(s);
                        }
                        else
                        {
                            TString histNameNPBFull = Form("%s_eta%d_pt%d_cut4", hbase.c_str(), ieta, ipt);
                            std::cout << "Warning: No valid scale for NPB overlay for " << histNameNPBFull << std::endl;
                            std::cout << "s: " << s << std::endl;
                            // No valid scale; drop overlay to avoid confusion
                            delete proj_data_npb;
                            proj_data_npb = nullptr;
                        }
                    }

                    // Optional: plot Data - scaled NPB in the SAME normalization convention as the cut0 overlay.
                    // (i.e. both are unit-area shapes; NPB is then scaled by the tail-matched fraction f_NPB)
                    if (makeDataMinusNPB && icut == 0 && proj_data_npb)
                    {
                        TH1D *proj_sub = (TH1D *)proj_data->Clone(Form("%s_px_data_minus_npb", histNameFull.Data()));
                        proj_sub->SetDirectory(nullptr);
                        proj_sub->Sumw2();
                        proj_sub->Add(proj_data_npb, -1.0);

                        // Determine y-range: match cut0's yMax choice (max * 1.3), only extend below 0 if needed.
                        double subMin = 0.0;
                        double subMax = 0.0;
                        bool subInit = false;
                        for (int ib = 1; ib <= proj_sub->GetNbinsX(); ++ib)
                        {
                            const double v = proj_sub->GetBinContent(ib);
                            if (!subInit)
                            {
                                subMin = v;
                                subMax = v;
                                subInit = true;
                            }
                            else
                            {
                                if (v < subMin) subMin = v;
                                if (v > subMax) subMax = v;
                            }
                        }
                        if (!subInit)
                        {
                            subMin = -1e-6;
                            subMax =  1e-6;
                        }

                        const double maxy_cut0_like = std::max(
                            std::max((double)proj_data->GetMaximum(), std::max((double)proj_sig->GetMaximum(), (double)proj_bkg->GetMaximum())),
                            (proj_data_npb ? (double)proj_data_npb->GetMaximum() : 0.0));

                        double yMax = (maxy_cut0_like > 0.0) ? (1.3 * maxy_cut0_like) : 1e-6;
                        // If subtraction goes negative, give it room; otherwise keep same baseline (0) as cut0.
                        double yMin = 0.0;
                        if (subMin < 0.0)
                        {
                            // extend below zero a bit more than the min bin so error bars/markers aren't on the frame
                            yMin = 1.1 * subMin;
                        }
                        if (yMax <= yMin) yMax = yMin + 1e-6;

                        TCanvas *c_sub = new TCanvas(
                            Form("c_sub_%s", histNameFull.Data()),
                            Form("Data - scaled NPB - %s", histNameFull.Data()),
                            600, 600);
                        c_sub->cd();

                        // Draw like cut0: MC lines first, then (Data - scaled NPB) markers on top.
                        proj_sig->SetTitle("");
                        proj_sig->SetYTitle("normalized counts");
                        proj_sig->GetYaxis()->SetTitleOffset(1.5);
                        proj_sig->SetXTitle(xaxisname.Data());
                        proj_sig->GetXaxis()->SetNdivisions(505);
                        proj_sig->GetYaxis()->SetRangeUser(yMin, yMax);
                        proj_sig->Draw("HIST");
                        proj_bkg->Draw("HIST SAME");
                        overlayVerticalErrors(proj_sig);
                        overlayVerticalErrors(proj_bkg);

                        proj_sub->SetStats(0);
                        proj_sub->SetLineColor(kBlack);
                        proj_sub->SetMarkerColor(kBlack);
                        proj_sub->SetMarkerStyle(proj_data->GetMarkerStyle());
                        proj_sub->SetMarkerSize(proj_data->GetMarkerSize());
                        proj_sub->Draw("ex0 SAME");

                        TLine l0(xaxismin, 0.0, xaxismax, 0.0);
                        l0.SetLineStyle(2);
                        l0.SetLineColor(kGray + 2);
                        l0.Draw("SAME");

                        myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                        myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                        myText(0.20, 0.80, 1, strleg3.c_str(), 0.04);
                        float pTlow = pT_bin_edges[ipt];
                        float pThigh = pT_bin_edges[ipt + 1];
                        myText(0.20, 0.75, 1, Form("%.0f<p_{T}<%.0fGeV, cut0", pTlow, pThigh), 0.04);

                        myMarkerLineText(0.55, 0.90, 1.5, kBlack, 20, kBlack, 1,
                                         "Data - scaled NPB", 0.045, true);
                        myMarkerLineText(0.55, 0.85, 0, kRed, 0, kRed, 1,
                                         "Signal MC", 0.045, true);
                        myMarkerLineText(0.55, 0.80, 0, kBlue, 0, kBlue, 1,
                                         "Inclusive MC", 0.045, true);

                        c_sub->SaveAs(Form("%s/dis_subNPB_%s.pdf", savePath.c_str(), histNamesave.Data()));
                        delete c_sub;
                        delete proj_sub;
                    }

                    float maxy = std::max({proj_data->GetMaximum(), proj_sig->GetMaximum(), proj_bkg->GetMaximum()});
                    if (proj_data_npb)
                    {
                        maxy = std::max(maxy, (float)proj_data_npb->GetMaximum());
                    }

                    TCanvas *c_proj = new TCanvas(
                        Form("c_proj_%s", histNameFull.Data()),
                        Form("ProjectionX - %s", histNameFull.Data()),
                        600, 600);
                    c_proj->cd();

                    // Draw them
                    proj_sig->SetYTitle("normalized counts");
                    proj_sig->GetYaxis()->SetTitleOffset(1.5);
                    proj_sig->SetXTitle(xaxisname.Data());
                    proj_sig->GetYaxis()->SetRangeUser(0, maxy * 1.3);
                    proj_sig->GetXaxis()->SetNdivisions(505);
                    proj_sig->SetStats(0);
                    proj_sig->Draw("HIST");

                    proj_bkg->Draw("HIST SAME");
                    overlayVerticalErrors(proj_sig);
                    overlayVerticalErrors(proj_bkg);
                    proj_data->Draw("ex0 SAME");
                    if (proj_data_npb)
                    {
                        proj_data_npb->SetLineColor(kGreen + 2);
                        proj_data_npb->SetMarkerColor(kGreen + 2);
                        proj_data_npb->SetStats(0);
                        proj_data_npb->Draw("hist SAME");
                        overlayVerticalErrors(proj_data_npb);
                    }

                    myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                    myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                    myText(0.20, 0.80, 1, strleg3.c_str(), 0.04);
                    float pTlow = pT_bin_edges[ipt];
                    float pThigh = pT_bin_edges[ipt + 1];
                    std::string bgcut;
                    if (icut == 0)
                        bgcut = "w/o nbkg cut";
                    if (icut == 1)
                        bgcut = "w/  nbkg cut";
                    if (icut == 2)
                        bgcut = "w/ tight cut";
                    if (icut == 3)
                        bgcut = "w/ nontight cut";
                    myText(0.2, 0.75, 1, Form("%.0f<p_{T}<%.0fGeV,%s", pTlow, pThigh, bgcut.c_str()), 0.04);

                    myMarkerLineText(0.6, 0.90, 1.5, kBlack, 20, kBlack, 1,
                                     "Data", 0.05, true);
                    myMarkerLineText(0.6, 0.85, 0, kRed, 0, kRed, 1,
                                     "Signal MC", 0.05, true);
                    myMarkerLineText(0.6, 0.80, 0, kBlue, 0, kBlue, 1,
                                     "Inclusive MC", 0.05, true);
                    if (icut == 0)
                    {
                        myMarkerLineText(0.6, 0.75, 0, kGreen + 2, 0, kGreen + 2, 2,
                                         "Data NPB", 0.05, true);
                    }

                    c_proj->SaveAs(Form("%s/dis_%s.pdf", savePath.c_str(), histNamesave.Data()));

                    // Create profile plots for background only
                    TProfile *pfx_bkg = h2_bkg->ProfileX(
                        Form("%s_pfx_bkg", histNameFull.Data()),
                        1, -1, "");

                    TCanvas *c_bkg = new TCanvas(
                        Form("c_bkg_%s", histNameFull.Data()),
                        Form("Bkg only - %s", histNameFull.Data()),
                        600, 600);
                    c_bkg->cd();

                    TH1D *proj_bkg_clone = (TH1D *)pfx_bkg->Clone(
                        Form("%s_px_bkgOnly", histNameFull.Data()));
                    proj_bkg_clone->SetLineColor(kBlue);

                    proj_bkg_clone->SetYTitle("<#it{E}_{T}^{iso}> [GeV]");
                    proj_bkg_clone->SetXTitle(xaxisname.Data());
                    proj_bkg_clone->SetStats(0);
                    float max = proj_bkg_clone->GetMaximum();
                    proj_bkg_clone->GetYaxis()->SetRangeUser(0, max * 1.3);

                    proj_bkg_clone->Draw("HIST");
                    overlayVerticalErrors(proj_bkg_clone);

                    float corr = h2_bkg->GetCorrelationFactor();

                    myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                    myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                    myText(0.20, 0.80, 1, strleg3.c_str(), 0.04);
                    myText(0.55, 0.90, 1, Form("%.0f<p_{T}<%.0fGeV,%s", pTlow, pThigh, bgcut.c_str()), 0.04);
                    myText(0.55, 0.85, 1, "Background MC", 0.04);
                    myText(0.55, 0.80, 1, Form("Correlation: %.3f", corr), 0.04);

                    c_bkg->SaveAs(Form("%s/pfx_%s.pdf", savePath.c_str(), histNamesave.Data()));

                    delete c_proj;
                    delete c_bkg;
                    if (proj_data_npb) delete proj_data_npb;
                } // end ipt
            }     // end ieta
        }         // end icut
    }             // end loop over histNames

    //==============================================================================
    // Part 1: BDT Overlay Plots (using isolation ET histograms binned in BDT)
    //==============================================================================
    std::cout << "Creating BDT overlay plots..." << std::endl;

    // Check if BDT-binned histograms exist
    TH1D *test_hist = (TH1D *)f_sig->Get("h_reco_iso_eta0_pt0_bdt0");
    if (!test_hist)
    {
        std::cout << "Warning: BDT-binned histograms not found in signal file." << std::endl;
        std::cout << "Please run ShowerShapeCheck.C with updated code to generate these histograms." << std::endl;
        std::cout << "Skipping Part 1 (BDT overlay plots)..." << std::endl;
    }
    else
    {
        for (int ieta = 0; ieta < nEtaBins; ieta++)
        {
            for (int ipt = 0; ipt < nPtBins; ipt++)
            {
                TCanvas *c = new TCanvas(Form("c_bdt_iso_eta%d_pt%d", ieta, ipt),
                                         Form("BDT overlay - Iso ET"),
                                         600, 600);

                TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.88);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.04);

                // Overlay different BDT bins (robust for arbitrary nBdtBins)
                std::vector<Color_t> colors;
                colors.reserve(nBdtBins);
                Color_t baseColors[] = {kRed, kBlue, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet + 1, kTeal + 2};
                const int nBaseColors = sizeof(baseColors) / sizeof(baseColors[0]);
                for (int ibdt = 0; ibdt < nBdtBins; ++ibdt)
                {
                    colors.push_back(baseColors[ibdt % nBaseColors]);
                }
                std::vector<TH1D *> proj_bdt(nBdtBins, nullptr);

                float maxY = 0;
                int nFound = 0;

                for (int ibdt = 0; ibdt < nBdtBins; ibdt++)
                {
                    // Get 1D histogram for this BDT bin
                    TString histName = Form("h_reco_iso_eta%d_pt%d_bdt%d", ieta, ipt, ibdt);
                    proj_bdt[ibdt] = (TH1D *)f_sig->Get(histName);

                    if (!proj_bdt[ibdt])
                    {
                        continue;
                    }
                    nFound++;

                    // Clone and rebin
                    proj_bdt[ibdt] = (TH1D *)proj_bdt[ibdt]->Clone(Form("%s_clone", histName.Data()));
                    proj_bdt[ibdt]->Rebin(4);

                    // Normalize to unit integral
                    if (proj_bdt[ibdt]->Integral() > 0)
                    {
                        proj_bdt[ibdt]->Scale(1.0 / proj_bdt[ibdt]->Integral());
                        // Scale by bin width
                        proj_bdt[ibdt]->Scale(1.0 / proj_bdt[ibdt]->GetBinWidth(1));
                    }

                    proj_bdt[ibdt]->SetLineColor(colors[ibdt]);
                    proj_bdt[ibdt]->SetLineWidth(2);
                    proj_bdt[ibdt]->SetStats(0);

                    if (proj_bdt[ibdt]->GetMaximum() > maxY)
                        maxY = proj_bdt[ibdt]->GetMaximum();

                    if (ibdt == 0)
                    {
                        proj_bdt[ibdt]->SetTitle("");
                        proj_bdt[ibdt]->GetXaxis()->SetTitle("Iso E_{T} [GeV]");
                        proj_bdt[ibdt]->GetXaxis()->SetRangeUser(-3, 15);
                        proj_bdt[ibdt]->GetXaxis()->SetNdivisions(505);
                        proj_bdt[ibdt]->GetYaxis()->SetTitle("Normalized counts / GeV");
                        proj_bdt[ibdt]->GetYaxis()->SetTitleOffset(1.5);
                        proj_bdt[ibdt]->GetYaxis()->SetRangeUser(0, maxY * 1.3);
                        proj_bdt[ibdt]->Draw("HIST");
                        overlayVerticalErrors(proj_bdt[ibdt]);
                    }
                    else
                    {
                        proj_bdt[ibdt]->Draw("HIST SAME");
                        overlayVerticalErrors(proj_bdt[ibdt]);
                    }

                    leg->AddEntry(proj_bdt[ibdt],
                                  Form("%.2f < BDT < %.2f", bdt_bins[ibdt], bdt_bins[ibdt + 1]),
                                  "l");
                }

                if (nFound == 0)
                {
                    delete c;
                    continue;
                }

                // Update max Y
                if (proj_bdt[0])
                    proj_bdt[0]->GetYaxis()->SetRangeUser(0, maxY * 1.3);

                myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                myText(0.20, 0.80, 1, strleg3.c_str(), 0.04);
                myText(0.20, 0.75, 1, strSigMC.c_str(), 0.04);
                myText(0.20, 0.70, 1, "", 0.04);
                myText(0.20, 0.65, 1,
                       Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                       0.04);

                leg->Draw();
                c->SaveAs(Form("%s/isoET_reco_BDTbin_overlay_eta%d_pt%d.pdf", savePath.c_str(), ieta, ipt));
                delete c;
            }
        }
    }

    //==============================================================================
    // Part 2: Tight vs Non-tight (shower-shape sideband) Overlay Plots
    // 
    //==============================================================================
    std::cout << "Creating isolation region overlay plots..." << std::endl;

    std::vector<std::string> histNames2D = {
        //"h2d_weta_cogx",
        //"h2d_wphi_cogx",
        //"h2d_et1",
        //"h2d_prob",
        //"h2d_e32_to_e35",
        //"h2d_e11_to_e33",
        //"h2d_bdt",
        "h2d_npb_score"
    };

    for (const auto &hbase : histNames2D)
    {
        for (int ieta = 0; ieta < nEtaBins; ieta++)
        {
            for (int ipt = 0; ipt < nPtBins; ipt++)
            {
                TCanvas *c = new TCanvas(Form("c_iso_%s_eta%d_pt%d", hbase.c_str(), ieta, ipt),
                                         Form("Iso overlay - %s", hbase.c_str()),
                                         600, 600);

                // x axis name is the base name after first 4 characters
                TString xaxisname = hbase.substr(4, hbase.size() - 4);
                float xaxismin, xaxismax; int nrebin;
                getAxisSettings(xaxisname, xaxismin, xaxismax, nrebin);

                // Tight shower-shape selection - cut2
                TString histName_tight = Form("%s_eta%d_pt%d_cut2", hbase.c_str(), ieta, ipt);
                TH2F *h2_mc_tight = (TH2F *)f_bkg_inc->Get(histName_tight);
                TH2F *h2_data_tight = (TH2F *)f_data->Get(histName_tight);

                // Non-tight shower-shape sideband - cut3
                TString histName_nt = Form("%s_eta%d_pt%d_cut3", hbase.c_str(), ieta, ipt);
                TH2F *h2_mc_nt = (TH2F *)f_bkg_inc->Get(histName_nt);
                TH2F *h2_data_nt = (TH2F *)f_data->Get(histName_nt);

                if (!h2_mc_tight || !h2_mc_nt || !h2_data_tight || !h2_data_nt)
                {
                    std::cerr << "Warning: Could not retrieve histograms for " << hbase << std::endl;
                    delete c;
                    continue;
                }

                TH1D *proj_mc_tight = h2_mc_tight->ProjectionX(Form("%s_px_mc_tight", histName_tight.Data()));
                TH1D *proj_mc_nt = h2_mc_nt->ProjectionX(Form("%s_px_mc_nt", histName_nt.Data()));
                TH1D *proj_data_tight = h2_data_tight->ProjectionX(Form("%s_px_data_tight", histName_tight.Data()));
                TH1D *proj_data_nt = h2_data_nt->ProjectionX(Form("%s_px_data_nt", histName_nt.Data()));

                if (nrebin > 1)
                {
                    proj_mc_tight->Rebin(nrebin);
                    proj_mc_nt->Rebin(nrebin);
                    proj_data_tight->Rebin(nrebin);
                    proj_data_nt->Rebin(nrebin);
                }

                proj_mc_tight->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
                proj_mc_nt->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
                proj_data_tight->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
                proj_data_nt->GetXaxis()->SetRangeUser(xaxismin, xaxismax);

                // Normalize
                scaleToUnit(proj_mc_tight);
                scaleToUnit(proj_mc_nt);
                scaleToUnit(proj_data_tight);
                scaleToUnit(proj_data_nt);

                // Inclusive MC as lines
                proj_mc_tight->SetLineColor(kBlue);
                proj_mc_tight->SetLineWidth(2);
                proj_mc_tight->SetStats(0);

                proj_mc_nt->SetLineColor(kRed);
                proj_mc_nt->SetLineWidth(2);
                proj_mc_nt->SetStats(0);

                // Data as markers
                proj_data_tight->SetLineColor(kBlack);
                proj_data_tight->SetMarkerColor(kBlack);
                proj_data_tight->SetMarkerStyle(20);
                proj_data_tight->SetMarkerSize(0.9);
                proj_data_tight->SetStats(0);

                proj_data_nt->SetLineColor(kGray + 2);
                proj_data_nt->SetMarkerColor(kGray + 2);
                proj_data_nt->SetMarkerStyle(24);
                proj_data_nt->SetMarkerSize(0.9);
                proj_data_nt->SetStats(0);

                float maxY = std::max({proj_mc_tight->GetMaximum(),
                                       proj_mc_nt->GetMaximum(),
                                       proj_data_tight->GetMaximum(),
                                       proj_data_nt->GetMaximum()});

                proj_mc_tight->SetTitle("");
                proj_mc_tight->GetXaxis()->SetTitle(xaxisname.Data());
                proj_mc_tight->GetXaxis()->SetNdivisions(505);
                proj_mc_tight->GetYaxis()->SetTitle("Normalized counts");
                proj_mc_tight->GetYaxis()->SetTitleOffset(1.5);
                proj_mc_tight->GetYaxis()->SetRangeUser(0, maxY * 1.3);
                proj_mc_tight->Draw("HIST");
                proj_mc_nt->Draw("HIST SAME");
                overlayVerticalErrors(proj_mc_tight);
                overlayVerticalErrors(proj_mc_nt);
                proj_data_tight->Draw("E SAME");
                proj_data_nt->Draw("E SAME");

                myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                myText(0.20, 0.80, 1, "Data + Inclusive MC", 0.04);
                myText(0.20, 0.75, 1,
                       Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                       0.04);

                TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.88);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.04);
                leg->AddEntry(proj_mc_tight, "Inclusive MC (tight)", "l");
                leg->AddEntry(proj_mc_nt, "Inclusive MC (non-tight)", "l");
                leg->AddEntry(proj_data_tight, "Data (tight)", "lep");
                leg->AddEntry(proj_data_nt, "Data (non-tight)", "lep");
                leg->Draw();
                // Name explicitly what is being compared to avoid confusion (esp. when xaxisname == "bdt")
                c->SaveAs(Form("%s/isoRegionOverlay_%s_tight_vs_nontight_incMC_withData_eta%d_pt%d.pdf", savePath.c_str(), xaxisname.Data(), ieta, ipt));
                delete c;
            }
        }
    }

    //==============================================================================
    // Part 3: Data vs MC Isolation ET Comparison
    //==============================================================================
    std::cout << "Creating data vs MC isolation ET comparison plots..." << std::endl;

    for (int ieta = 0; ieta < nEtaBins; ieta++)
    {
        for (int ipt = 0; ipt < nPtBins; ipt++)
        {
            TCanvas *c = new TCanvas(Form("c_isoET_datamc_eta%d_pt%d", ieta, ipt),
                                     Form("Data vs MC Iso ET"),
                                     600, 600);

            // Data tight iso
            TH1D *h_data = (TH1D *)f_data->Get(Form("h_tight_reco_iso_eta%d_pt%d", ieta, ipt));
            // MC background tight iso
            TH1D *h_mc = (TH1D *)f_bkg_inc->Get(Form("h_tight_reco_iso_eta%d_pt%d", ieta, ipt));

            if (!h_data || !h_mc)
            {
                std::cerr << "Warning: Could not retrieve iso ET histograms for eta" << ieta << " pt" << ipt << std::endl;
                delete c;
                continue;
            }

            // Clone to avoid modifying originals
            TH1D *h_data_clone = (TH1D *)h_data->Clone(Form("h_data_clone_eta%d_pt%d", ieta, ipt));
            TH1D *h_mc_clone = (TH1D *)h_mc->Clone(Form("h_mc_clone_eta%d_pt%d", ieta, ipt));

            // Rebin for better statistics
            h_data_clone->Rebin(4);
            h_mc_clone->Rebin(4);

            // Normalize
            if (h_data_clone->Integral() > 0)
                h_data_clone->Scale(1.0 / h_data_clone->Integral());
            if (h_mc_clone->Integral() > 0)
                h_mc_clone->Scale(1.0 / h_mc_clone->Integral());

            // Scale by bin width
            h_data_clone->Scale(1.0 / h_data_clone->GetBinWidth(1));
            h_mc_clone->Scale(1.0 / h_mc_clone->GetBinWidth(1));

            h_data_clone->SetLineColor(kBlack);
            h_data_clone->SetMarkerColor(kBlack);
            h_data_clone->SetMarkerStyle(20);
            h_data_clone->SetMarkerSize(1.0);
            h_data_clone->SetStats(0);

            h_mc_clone->SetLineColor(kRed);
            h_mc_clone->SetLineWidth(2);
            h_mc_clone->SetMarkerColor(kRed);
            h_mc_clone->SetStats(0);

            float maxY = std::max(h_data_clone->GetMaximum(), h_mc_clone->GetMaximum());

            h_mc_clone->SetTitle("");
            h_mc_clone->GetXaxis()->SetTitle("Iso E_{T} [GeV]");
            h_mc_clone->GetXaxis()->SetNdivisions(505);
            h_mc_clone->GetYaxis()->SetTitle("Normalized counts / GeV");
            h_mc_clone->GetYaxis()->SetTitleOffset(1.5);
            h_mc_clone->GetYaxis()->SetRangeUser(0, maxY * 1.3);
            h_mc_clone->GetXaxis()->SetRangeUser(-3, 15);
            h_mc_clone->Draw("HIST");
            overlayVerticalErrors(h_mc_clone);
            h_data_clone->Draw("E SAME");

            myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.80, 1, strleg3.c_str(), 0.04);
            myText(0.20, 0.75, 1,
                   Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                   0.04);

            myMarkerLineText(0.60, 0.65, 1.5, kBlack, 20, kBlack, 1,
                           "Data", 0.04, true);
            myMarkerLineText(0.60, 0.60, 0, kRed, 0, kRed, 1,
                           "Inclusive MC", 0.04, true);

            c->SaveAs(Form("%s/isoET_datamc_eta%d_pt%d.pdf", savePath.c_str(), ieta, ipt));
            delete c;
        }
    }

    //==============================================================================
    // Part 4: Cluster-MBD Time Distribution in Different NPB Score Ranges
    //==============================================================================
    std::cout << "Creating NPB score vs time overlay plots..." << std::endl;

    std::vector<std::pair<std::string, std::string>> npb_time_hist_variants = {
        {"h_npb_score_vs_time", ""}, {"h_npb_score_vs_time_clean", "_clean"}, {"h_npb_score_vs_time_strict", "_strict"}
    };

    for (const auto& [hist_base, hist_suffix] : npb_time_hist_variants)
    {
        for (const auto& [file_tag, file_ptr, file_color, file_label] : all_samples)
        {
            if (!file_ptr || file_ptr->IsZombie()) continue;

            for (int ieta = 0; ieta < nEtaBins; ++ieta)
            {
                for (int ipt = 0; ipt < nPtBins; ++ipt)
                {
                    TString hist_name = Form("%s_eta%d_pt%d", hist_base.c_str(), ieta, ipt);
                    TH2D *h2_npb_time = dynamic_cast<TH2D*>(file_ptr->Get(hist_name));
                    if (!h2_npb_time) continue;

                    TCanvas *c = new TCanvas(Form("c_npb_time%s_%s_eta%d_pt%d", hist_suffix.c_str(), file_tag.c_str(), ieta, ipt),
                                             "", 600, 600);
                    TLegend *leg = new TLegend(0.65, 0.60, 0.88, 0.88);
                    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);

                    std::vector<TH1D*> time_projections;
                    float maxY = 0.0;

                    for (int inpb = 0; inpb < nNpbBins; ++inpb)
                    {
                        int ybin_min = h2_npb_time->GetYaxis()->FindBin(npb_score_ranges[inpb].first + 1e-6);
                        int ybin_max = h2_npb_time->GetYaxis()->FindBin(npb_score_ranges[inpb].second - 1e-6);
                        TH1D *proj = h2_npb_time->ProjectionX(Form("%s_time_npb%d", hist_name.Data(), inpb), ybin_min, ybin_max);
                        if (!proj || proj->Integral() < 1e-12) { if (proj) delete proj; continue; }

                        scaleToUnit(proj);
                        proj->SetLineColor(npb_time_colors[inpb]); proj->SetLineWidth(2); proj->SetStats(0);
                        time_projections.push_back(proj);
                        if (proj->GetMaximum() > maxY) maxY = proj->GetMaximum();
                        leg->AddEntry(proj, Form("%.1f < NPB < %.1f", npb_score_ranges[inpb].first, npb_score_ranges[inpb].second), "l");
                    }

                    if (time_projections.empty()) { delete c; delete leg; continue; }

                    time_projections[0]->SetTitle("");
                    time_projections[0]->GetXaxis()->SetTitle("Cluster-MBD Time [ns]");
                    time_projections[0]->GetXaxis()->SetRangeUser(hist_suffix == "_strict" ? -10 : -20, 20);
                    time_projections[0]->GetXaxis()->SetNdivisions(505);
                    time_projections[0]->GetYaxis()->SetTitle("Normalized counts");
                    time_projections[0]->GetYaxis()->SetTitleOffset(1.5);
                    time_projections[0]->GetYaxis()->SetRangeUser(0, maxY * 1.3);
                    time_projections[0]->Draw("HIST"); overlayVerticalErrors(time_projections[0]);
                    for (size_t i = 1; i < time_projections.size(); ++i) { time_projections[i]->Draw("HIST SAME"); overlayVerticalErrors(time_projections[i]); }

                    myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                    myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                    myText(0.20, 0.80, 1, strleg3.c_str(), 0.04);
                    myText(0.20, 0.75, 1, file_label.c_str(), 0.04);
                    myText(0.20, 0.70, 1, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), 0.04);
                    if (!hist_suffix.empty()) myText(0.20, 0.65, 1, Form("Selection: %s", hist_suffix.substr(1).c_str()), 0.04);
                    leg->Draw();

                    c->SaveAs(Form("%s/npb_time_overlay%s_%s_eta%d_pt%d.pdf", savePath.c_str(), hist_suffix.c_str(), file_tag.c_str(), ieta, ipt));
                    delete c; delete leg;
                }
            }
        }
    }

    // Part 4b: Alternative NPB score ranges (threshold-based)
    std::vector<std::pair<double, double>> npb_threshold_ranges = {
        {0.0, 0.2}, {0.2, 1.0}, {0.5, 1.0}, {0.8, 1.0}
    };
    std::vector<Color_t> npb_threshold_colors = {kRed, kBlue, kGreen+2, kMagenta+1};
    int nThresholdBins = npb_threshold_ranges.size();

    for (const auto& [hist_base, hist_suffix] : npb_time_hist_variants)
    {
        for (const auto& [file_tag, file_ptr, file_color, file_label] : all_samples)
        {
            if (!file_ptr || file_ptr->IsZombie()) continue;

            for (int ieta = 0; ieta < nEtaBins; ++ieta)
            {
                for (int ipt = 0; ipt < nPtBins; ++ipt)
                {
                    TString hist_name = Form("%s_eta%d_pt%d", hist_base.c_str(), ieta, ipt);
                    TH2D *h2_npb_time = dynamic_cast<TH2D*>(file_ptr->Get(hist_name));
                    if (!h2_npb_time) continue;

                    TCanvas *c = new TCanvas(Form("c_npb_time_thresh%s_%s_eta%d_pt%d", hist_suffix.c_str(), file_tag.c_str(), ieta, ipt),
                                             "", 600, 600);
                    TLegend *leg = new TLegend(0.65, 0.65, 0.88, 0.88);
                    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);

                    std::vector<TH1D*> time_projections;
                    float maxY = 0.0;

                    for (int inpb = 0; inpb < nThresholdBins; ++inpb)
                    {
                        int ybin_min = h2_npb_time->GetYaxis()->FindBin(npb_threshold_ranges[inpb].first + 1e-6);
                        int ybin_max = h2_npb_time->GetYaxis()->FindBin(npb_threshold_ranges[inpb].second - 1e-6);
                        TH1D *proj = h2_npb_time->ProjectionX(Form("%s_time_thresh%d", hist_name.Data(), inpb), ybin_min, ybin_max);
                        if (!proj || proj->Integral() < 1e-12) { if (proj) delete proj; continue; }

                        scaleToUnit(proj);
                        proj->SetLineColor(npb_threshold_colors[inpb]); proj->SetLineWidth(2); proj->SetStats(0);
                        time_projections.push_back(proj);
                        if (proj->GetMaximum() > maxY) maxY = proj->GetMaximum();
                        leg->AddEntry(proj, Form("%.1f < NPB < %.1f", npb_threshold_ranges[inpb].first, npb_threshold_ranges[inpb].second), "l");
                    }

                    if (time_projections.empty()) { delete c; delete leg; continue; }

                    time_projections[0]->SetTitle("");
                    time_projections[0]->GetXaxis()->SetTitle("Cluster-MBD Time [ns]");
                    time_projections[0]->GetXaxis()->SetRangeUser(hist_suffix == "_strict" ? -10 : -20, 20);
                    time_projections[0]->GetXaxis()->SetNdivisions(505);
                    time_projections[0]->GetYaxis()->SetTitle("Normalized counts");
                    time_projections[0]->GetYaxis()->SetTitleOffset(1.5);
                    time_projections[0]->GetYaxis()->SetRangeUser(0, maxY * 1.3);
                    time_projections[0]->Draw("HIST"); overlayVerticalErrors(time_projections[0]);
                    for (size_t i = 1; i < time_projections.size(); ++i) { time_projections[i]->Draw("HIST SAME"); overlayVerticalErrors(time_projections[i]); }

                    myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                    myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                    myText(0.20, 0.80, 1, strleg3.c_str(), 0.04);
                    myText(0.20, 0.75, 1, file_label.c_str(), 0.04);
                    myText(0.20, 0.70, 1, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), 0.04);
                    if (!hist_suffix.empty()) myText(0.20, 0.65, 1, Form("Selection: %s", hist_suffix.substr(1).c_str()), 0.04);
                    leg->Draw();

                    c->SaveAs(Form("%s/npb_time_threshold%s_%s_eta%d_pt%d.pdf", savePath.c_str(), hist_suffix.c_str(), file_tag.c_str(), ieta, ipt));
                    delete c; delete leg;
                }
            }
        }
    }

    std::cout << "Part 4 complete: NPB score vs time plots created" << std::endl;

    //==============================================================================
    // Part 5: Inclusive MC Shower Shape Distributions (cut0 vs cut1 comparison)
    //==============================================================================
    std::cout << "Creating cut0 vs cut1 comparison plots for inclusive MC..." << std::endl;

    for (const auto &hbase : histNames)  // Reuse histNames from Part 0
    {
        for (int ieta = 0; ieta < nEtaBins; ++ieta)
        {
            for (int ipt = 0; ipt < nPtBins; ++ipt)
            {
                TCanvas *c = new TCanvas(
                    Form("c_cutcomp_%s_eta%d_pt%d", hbase.c_str(), ieta, ipt),
                    Form("Cut comparison - %s", hbase.c_str()),
                    600, 600);

                TString xaxisname = hbase.substr(4, hbase.size() - 4);
                float xaxismin, xaxismax; int nrebin;
                getAxisSettings(xaxisname, xaxismin, xaxismax, nrebin);

                // Retrieve histograms for cut0 and cut1
                TString hist_cut0_name = Form("%s_eta%d_pt%d_cut0", hbase.c_str(), ieta, ipt);
                TString hist_cut1_name = Form("%s_eta%d_pt%d_cut1", hbase.c_str(), ieta, ipt);

                TH2F *h2_cut0 = dynamic_cast<TH2F*>(f_bkg_inc->Get(hist_cut0_name));
                TH2F *h2_cut1 = dynamic_cast<TH2F*>(f_bkg_inc->Get(hist_cut1_name));

                if (!h2_cut0 || !h2_cut1)
                {
                    std::cerr << "Warning: Could not retrieve histograms for "
                              << hbase << " eta" << ieta << " pt" << ipt << std::endl;
                    delete c;
                    continue;
                }

                // Rebin and set range
                h2_cut0->RebinX(nrebin);
                h2_cut1->RebinX(nrebin);
                h2_cut0->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
                h2_cut1->GetXaxis()->SetRangeUser(xaxismin, xaxismax);

                // Project to 1D
                TH1D *proj_cut0 = h2_cut0->ProjectionX(
                    Form("%s_px_cut0", hist_cut0_name.Data()));
                TH1D *proj_cut1 = h2_cut1->ProjectionX(
                    Form("%s_px_cut1", hist_cut1_name.Data()));

                // Normalize to unit area
                scaleToUnit(proj_cut0);
                scaleToUnit(proj_cut1);

                // Style cut0 (baseline - blue)
                proj_cut0->SetLineColor(kBlue);
                proj_cut0->SetLineWidth(2);
                proj_cut0->SetStats(0);

                // Style cut1 (common cuts - red)
                proj_cut1->SetLineColor(kRed);
                proj_cut1->SetLineWidth(2);
                proj_cut1->SetStats(0);

                // Determine Y-axis range
                float maxY = std::max(proj_cut0->GetMaximum(), proj_cut1->GetMaximum());

                // Draw
                proj_cut0->SetTitle("");
                proj_cut0->GetXaxis()->SetTitle(xaxisname.Data());
                proj_cut0->GetXaxis()->SetNdivisions(505);
                proj_cut0->GetYaxis()->SetTitle("Normalized counts");
                proj_cut0->GetYaxis()->SetTitleOffset(1.5);
                proj_cut0->GetYaxis()->SetRangeUser(0, maxY * 1.3);
                proj_cut0->Draw("HIST");
                overlayVerticalErrors(proj_cut0);

                proj_cut1->Draw("HIST SAME");
                overlayVerticalErrors(proj_cut1);

                // Add text labels
                myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                myText(0.20, 0.80, 1, "Inclusive MC", 0.04);
                myText(0.20, 0.75, 1,
                       Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]),
                       0.04);

                // Legend
                TLegend *leg = new TLegend(0.55, 0.75, 0.88, 0.88);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.04);
                leg->AddEntry(proj_cut0, "no cut", "l");
                leg->AddEntry(proj_cut1, "npb cut", "l");
                leg->Draw();

                c->SaveAs(Form("%s/cut_comparison_%s_eta%d_pt%d.pdf",
                              savePath.c_str(), xaxisname.Data(), ieta, ipt));

                delete c;
                delete leg;
            }
        }
    }

    std::cout << "Part 5 complete: Cut comparison plots created" << std::endl;

    //==============================================================================
    // Part 6: NPB Selection Efficiency vs Threshold (Signal MC & Inclusive MC)
    //==============================================================================
    std::cout << "Creating NPB selection efficiency vs threshold plots..." << std::endl;

    // MC samples only (exclude data)
    std::vector<std::tuple<std::string, TFile*, Color_t, std::string>> mc_samples = {
        {"inclusive", f_bkg_inc, kBlue, "Inclusive MC"}, {"signal", f_sig, kRed, "Signal MC"}
    };

    for (const auto& [mc_tag, mc_file, mc_color, mc_label] : mc_samples)
    {
        if (!mc_file || mc_file->IsZombie()) continue;

        // Individual plots per (eta, pt)
        for (int ieta = 0; ieta < nEtaBins; ++ieta)
        {
            for (int ipt = 0; ipt < nPtBins; ++ipt)
            {
                TString hist_npb_name = Form("h2d_npb_score_eta%d_pt%d_cut0", ieta, ipt);
                TH2F *h2_npb = dynamic_cast<TH2F*>(mc_file->Get(hist_npb_name));
                if (!h2_npb) continue;

                TH1D *h_npb = h2_npb->ProjectionX(Form("%s_px_%s", hist_npb_name.Data(), mc_tag.c_str()));
                if (!h_npb || h_npb->Integral() < 1e-12) { if (h_npb) delete h_npb; continue; }

                int nbins = h_npb->GetNbinsX();
                TH1D *h_eff = new TH1D(Form("h_npb_eff_%s_eta%d_pt%d", mc_tag.c_str(), ieta, ipt),
                                       "", nbins, h_npb->GetXaxis()->GetXmin(), h_npb->GetXaxis()->GetXmax());
                h_eff->SetDirectory(nullptr);

                double total = h_npb->Integral(1, nbins);
                for (int ibin = 1; ibin <= nbins; ++ibin)
                {
                    double eff = (total > 0) ? h_npb->Integral(ibin, nbins) / total : 0.0;
                    h_eff->SetBinContent(ibin, eff);
                    h_eff->SetBinError(ibin, sqrt(eff * (1 - eff) / total));
                }

                TCanvas *c = new TCanvas(Form("c_npb_eff_%s_eta%d_pt%d", mc_tag.c_str(), ieta, ipt),
                                         Form("NPB Selection Efficiency (%s)", mc_label.c_str()), 600, 600);

                h_eff->SetLineColor(mc_color); h_eff->SetMarkerColor(mc_color);
                h_eff->SetMarkerStyle(20); h_eff->SetMarkerSize(0.8);
                h_eff->SetLineWidth(1); h_eff->SetStats(0); h_eff->SetTitle("");
                h_eff->GetXaxis()->SetTitle("NPB Score Threshold"); h_eff->GetXaxis()->SetNdivisions(505);
                h_eff->GetYaxis()->SetTitle("Selection Efficiency"); h_eff->GetYaxis()->SetTitleOffset(1.5);
                h_eff->GetYaxis()->SetRangeUser(0.8, 1.1);
                h_eff->Draw("E P");

                myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                myText(0.20, 0.80, 1, mc_label.c_str(), 0.04);
                myText(0.20, 0.75, 1, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), 0.04);

                c->SaveAs(Form("%s/npb_efficiency_vs_threshold_%s_eta%d_pt%d.pdf", savePath.c_str(), mc_tag.c_str(), ieta, ipt));
                delete h_eff; delete h_npb; delete c;
            }
        }

        // Overlay plot showing all pT bins together
        TCanvas *c_overlay = new TCanvas(Form("c_npb_eff_overlay_%s", mc_tag.c_str()),
                                         Form("NPB Efficiency Overlay (%s)", mc_label.c_str()), 700, 600);
        TLegend *leg = new TLegend(0.55, 0.15, 0.88, 0.45);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);

        std::vector<TH1D*> eff_hists;
        bool first = true;
        int ieta = 0;

        for (int ipt = 0; ipt < nPtBins; ++ipt)
        {
            TString hist_npb_name = Form("h2d_npb_score_eta%d_pt%d_cut0", ieta, ipt);
            TH2F *h2_npb = dynamic_cast<TH2F*>(mc_file->Get(hist_npb_name));
            if (!h2_npb) continue;

            TH1D *h_npb = h2_npb->ProjectionX(Form("%s_px_%s_overlay", hist_npb_name.Data(), mc_tag.c_str()));
            if (!h_npb || h_npb->Integral() < 1e-12) { if (h_npb) delete h_npb; continue; }

            int nbins = h_npb->GetNbinsX();
            TH1D *h_eff = new TH1D(Form("h_npb_eff_%s_overlay_pt%d", mc_tag.c_str(), ipt),
                                   "", nbins, h_npb->GetXaxis()->GetXmin(), h_npb->GetXaxis()->GetXmax());
            h_eff->SetDirectory(nullptr);

            double total = h_npb->Integral(1, nbins);
            for (int ibin = 1; ibin <= nbins; ++ibin)
                h_eff->SetBinContent(ibin, (total > 0) ? h_npb->Integral(ibin, nbins) / total : 0.0);

            h_eff->SetLineColor(pt_colors[ipt % pt_colors.size()]);
            h_eff->SetMarkerColor(pt_colors[ipt % pt_colors.size()]);
            h_eff->SetMarkerStyle(20 + (ipt % 10)); h_eff->SetMarkerSize(0.8);
            h_eff->SetLineWidth(1); h_eff->SetStats(0);

            if (first)
            {
                h_eff->SetTitle(""); h_eff->GetXaxis()->SetTitle("NPB Score Threshold");
                h_eff->GetXaxis()->SetNdivisions(505); h_eff->GetYaxis()->SetTitle("Selection Efficiency");
                h_eff->GetYaxis()->SetTitleOffset(1.5); h_eff->GetYaxis()->SetRangeUser(0.5, 1.2);
                h_eff->Draw("E P"); first = false;
            }
            else h_eff->Draw("E P SAME");

            eff_hists.push_back(h_eff);
            leg->AddEntry(h_eff, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), "lep");
            delete h_npb;
        }

        if (!eff_hists.empty())
        {
            myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.80, 1, mc_label.c_str(), 0.04);
            leg->Draw();
            c_overlay->SaveAs(Form("%s/npb_efficiency_vs_threshold_%s_overlay.pdf", savePath.c_str(), mc_tag.c_str()));
        }

        for (auto h : eff_hists) delete h;
        delete leg; delete c_overlay;
    }

    std::cout << "Part 6 complete: NPB efficiency plots created" << std::endl;

    //==============================================================================
    // Part 7: NPB Score Ratio vs MBD Sigma (using h_mbd_avgsigma_vs_npb_score)
    //         Ratio = N(0.5 < NPB < 0.9) / N(0.9 < NPB < 1.0) as function of MBD sigma
    //==============================================================================
    std::cout << "Creating NPB score ratio vs MBD sigma plots..." << std::endl;

    for (const auto& [file_tag, file_ptr, file_color, file_label] : all_samples)
    {
        if (!file_ptr || file_ptr->IsZombie()) continue;

        for (int ieta = 0; ieta < nEtaBins; ++ieta)
        {
            for (int ipt = 0; ipt < nPtBins; ++ipt)
            {
                TString hist_name = Form("h_mbd_avgsigma_vs_npb_score_eta%d_pt%d", ieta, ipt);
                TH2D *h2 = dynamic_cast<TH2D*>(file_ptr->Get(hist_name));

                if (!h2)
                {
                    std::cerr << "Warning: Could not retrieve " << hist_name << " from " << file_tag << std::endl;
                    continue;
                }

                // Get bin ranges for NPB score (X-axis: 0 to 1, 100 bins)
                // [0.5, 0.9] and [0.9, 1.0]
                int xbin_05 = h2->GetXaxis()->FindBin(0.5 + 1e-6);
                int xbin_09 = h2->GetXaxis()->FindBin(0.9 - 1e-6);
                int xbin_09_start = h2->GetXaxis()->FindBin(0.9 + 1e-6);
                int xbin_10 = h2->GetXaxis()->FindBin(1.0 - 1e-6);

                // Create ratio histogram (same Y-axis binning as MBD sigma)
                int nbins_y = h2->GetNbinsY();
                double ymin = h2->GetYaxis()->GetXmin();
                double ymax = h2->GetYaxis()->GetXmax();

                TH1D *h_ratio = new TH1D(Form("h_npb_ratio_vs_mbd_sigma_%s_eta%d_pt%d", file_tag.c_str(), ieta, ipt),
                                         "", nbins_y, ymin, ymax);
                h_ratio->SetDirectory(nullptr);

                // For each MBD sigma bin (Y-axis), compute ratio
                for (int iy = 1; iy <= nbins_y; ++iy)
                {
                    double n_mid = 0.0;  // NPB in [0.5, 0.9]
                    double n_high = 0.0; // NPB in [0.9, 1.0]

                    for (int ix = xbin_05; ix <= xbin_09; ++ix)
                        n_mid += h2->GetBinContent(ix, iy);

                    for (int ix = xbin_09_start; ix <= xbin_10; ++ix)
                        n_high += h2->GetBinContent(ix, iy);

                    double ratio = 0.0;
                    double ratio_err = 0.0;
                    if (n_high > 0)
                    {
                        ratio = n_mid / n_high;
                        // Binomial-like error approximation
                        ratio_err = ratio * sqrt(1.0/n_mid + 1.0/n_high);
                        if (n_mid < 1) ratio_err = 0;
                    }

                    h_ratio->SetBinContent(iy, ratio);
                    h_ratio->SetBinError(iy, ratio_err);
                }

                // Create canvas
                TCanvas *c = new TCanvas(Form("c_npb_ratio_mbd_%s_eta%d_pt%d", file_tag.c_str(), ieta, ipt),
                                         "NPB Ratio vs MBD Sigma", 600, 600);

                h_ratio->SetLineColor(file_color);
                h_ratio->SetMarkerColor(file_color);
                h_ratio->SetMarkerStyle(20);
                h_ratio->SetMarkerSize(0.8);
                h_ratio->SetLineWidth(1);
                h_ratio->SetStats(0);
                h_ratio->SetTitle("");
                h_ratio->GetXaxis()->SetTitle("MBD Avg #sigma_{t} [ns]");
                h_ratio->GetXaxis()->SetNdivisions(505);
                h_ratio->GetYaxis()->SetTitle("N(0.5<NPB<0.9) / N(0.9<NPB<1.0)");
                h_ratio->GetYaxis()->SetTitleOffset(1.5);
                h_ratio->GetXaxis()->SetRangeUser(0, 5);
                // Auto-range Y based on content
                double max_ratio = h_ratio->GetMaximum();
                h_ratio->GetYaxis()->SetRangeUser(0, max_ratio > 0 ? max_ratio * 1.3 : 1.0);
                h_ratio->Draw("E P");

                myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                myText(0.20, 0.80, 1, file_label.c_str(), 0.04);
                myText(0.20, 0.75, 1, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), 0.04);

                c->SaveAs(Form("%s/npb_ratio_vs_mbd_sigma_%s_eta%d_pt%d.pdf", savePath.c_str(), file_tag.c_str(), ieta, ipt));

                delete h_ratio;
                delete c;
            }
        }

        // Overlay plot for all pT bins
        TCanvas *c_overlay = new TCanvas(Form("c_npb_ratio_mbd_overlay_%s", file_tag.c_str()),
                                         Form("NPB Ratio vs MBD Sigma Overlay (%s)", file_label.c_str()), 700, 600);
        TLegend *leg = new TLegend(0.55, 0.60, 0.88, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);

        std::vector<TH1D*> ratio_hists;
        bool first = true;
        double global_max = 0;
        int ieta = 0;

        for (int ipt = 0; ipt < nPtBins; ++ipt)
        {
            TString hist_name = Form("h_mbd_avgsigma_vs_npb_score_eta%d_pt%d", ieta, ipt);
            TH2D *h2 = dynamic_cast<TH2D*>(file_ptr->Get(hist_name));
            if (!h2) continue;

            int xbin_05 = h2->GetXaxis()->FindBin(0.5 + 1e-6);
            int xbin_09 = h2->GetXaxis()->FindBin(0.9 - 1e-6);
            int xbin_09_start = h2->GetXaxis()->FindBin(0.9 + 1e-6);
            int xbin_10 = h2->GetXaxis()->FindBin(1.0 - 1e-6);

            int nbins_y = h2->GetNbinsY();
            TH1D *h_ratio = new TH1D(Form("h_npb_ratio_overlay_%s_pt%d", file_tag.c_str(), ipt),
                                     "", nbins_y, h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax());
            h_ratio->SetDirectory(nullptr);

            for (int iy = 1; iy <= nbins_y; ++iy)
            {
                double n_mid = 0.0, n_high = 0.0;
                for (int ix = xbin_05; ix <= xbin_09; ++ix) n_mid += h2->GetBinContent(ix, iy);
                for (int ix = xbin_09_start; ix <= xbin_10; ++ix) n_high += h2->GetBinContent(ix, iy);

                double ratio = (n_high > 0) ? n_mid / n_high : 0.0;
                double ratio_err = (n_high > 0 && n_mid > 0) ? ratio * sqrt(1.0/n_mid + 1.0/n_high) : 0.0;
                h_ratio->SetBinContent(iy, ratio);
                h_ratio->SetBinError(iy, ratio_err);
            }

            h_ratio->SetLineColor(pt_colors[ipt % pt_colors.size()]);
            h_ratio->SetMarkerColor(pt_colors[ipt % pt_colors.size()]);
            h_ratio->SetMarkerStyle(20 + (ipt % 10));
            h_ratio->SetMarkerSize(0.8);
            h_ratio->SetLineWidth(1);
            h_ratio->SetStats(0);

            if (h_ratio->GetMaximum() > global_max) global_max = h_ratio->GetMaximum();

            if (first)
            {
                h_ratio->SetTitle("");
                h_ratio->GetXaxis()->SetTitle("MBD Avg #sigma_{t} [ns]");
                h_ratio->GetXaxis()->SetNdivisions(505);
                h_ratio->GetXaxis()->SetRangeUser(0, 5);
                h_ratio->GetYaxis()->SetTitle("N(0.5<NPB<0.9) / N(0.9<NPB<1.0)");
                h_ratio->GetYaxis()->SetTitleOffset(1.5);
                h_ratio->Draw("E P");
                first = false;
            }
            else h_ratio->Draw("E P SAME");

            ratio_hists.push_back(h_ratio);
            leg->AddEntry(h_ratio, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), "lep");
        }

        if (!ratio_hists.empty())
        {
            ratio_hists[0]->GetYaxis()->SetRangeUser(0, global_max * 1.3);
            myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
            myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
            myText(0.20, 0.80, 1, file_label.c_str(), 0.04);
            leg->Draw();
            c_overlay->SaveAs(Form("%s/npb_ratio_vs_mbd_sigma_%s_overlay.pdf", savePath.c_str(), file_tag.c_str()));
        }

        for (auto h : ratio_hists) delete h;
        delete leg;
        delete c_overlay;
    }

    std::cout << "Part 7 complete: NPB ratio vs MBD sigma plots created" << std::endl;

    //==============================================================================
    // Close files
    //==============================================================================
    f_data->Close();
    f_sig->Close();
    f_bkg_inc->Close();
    if (hasBkgOnly && f_bkg_only && f_bkg_only != f_bkg_inc)
    {
        f_bkg_only->Close();
    }

    std::cout << "All plots created successfully!" << std::endl;
    std::cout << "Output saved to: " << savePath << std::endl;
}
