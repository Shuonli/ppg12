// plot_showershapes_variations.C
// Loops over all config_showershape*.yaml files and produces per-variant plot subdirectories
// under plotting/figures/{suffix}/.
//
// Usage:
//   root -l -b -q 'plot_showershapes_variations.C'                       // all variants
//   root -l -b -q 'plot_showershapes_variations.C("showershape_nom")'    // one variant

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

// ---------------------------------------------------------------------------
// plotOneConfig: process a single config variant.
// configname is the bare filename (e.g. "config_showershape_nom.yaml").
// Returns false if any required input file is missing/zombie.
// ---------------------------------------------------------------------------
static bool plotOneConfig(const std::string &configname)
{
    // Derive suffix from configname: strip path, "config_" prefix, ".yaml" extension
    std::string config_suffix = configname;
    {
        size_t slash = config_suffix.rfind('/');
        if (slash != std::string::npos) config_suffix = config_suffix.substr(slash + 1);
        if (config_suffix.size() > 5 && config_suffix.substr(config_suffix.size() - 5) == ".yaml")
            config_suffix = config_suffix.substr(0, config_suffix.size() - 5);
        if (config_suffix.size() > 7 && config_suffix.substr(0, 7) == "config_")
            config_suffix = config_suffix.substr(7);
    }
    std::cout << "[plotOneConfig] Config suffix: " << config_suffix << std::endl;

    // Per-variant output directory
    string savePath = "figures/" + config_suffix;
    gSystem->Exec(Form("mkdir -p %s", savePath.c_str()));

    // Load config for pT/BDT bins
    YAML::Node config = YAML::LoadFile(("../efficiencytool/" + configname).c_str());

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

    // Open files
    const std::string dataFile        = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histoshower_shape_" + config_suffix + ".root";
    const std::string sigFile         = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_signal_" + config_suffix + ".root";
    const std::string bkgInclusiveFile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_jet_inclusive_" + config_suffix + ".root";
    const std::string bkgOnlyFile     = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_jet_" + config_suffix + ".root";

    TFile *f_data    = TFile::Open(dataFile.c_str(), "READ");
    TFile *f_sig     = TFile::Open(sigFile.c_str(), "READ");
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

    // Lambda to close and delete all open files before early return
    auto closeAll = [&]() {
        if (f_data)    { f_data->Close();    delete f_data;    }
        if (f_sig)     { f_sig->Close();     delete f_sig;     }
        if (f_bkg_inc) { f_bkg_inc->Close(); delete f_bkg_inc; }
        if (hasBkgOnly && f_bkg_only && f_bkg_only != f_bkg_inc)
            { f_bkg_only->Close(); delete f_bkg_only; }
    };

    if (!f_data || f_data->IsZombie())
    {
        std::cerr << "[SKIP] " << config_suffix << ": Could not open data file: " << dataFile << std::endl;
        closeAll();
        return false;
    }
    if (!f_sig || f_sig->IsZombie())
    {
        std::cerr << "[SKIP] " << config_suffix << ": Could not open signal file: " << sigFile << std::endl;
        closeAll();
        return false;
    }
    if (!f_bkg_inc || f_bkg_inc->IsZombie())
    {
        std::cerr << "[SKIP] " << config_suffix << ": Could not open inclusive background file: " << bkgInclusiveFile << std::endl;
        closeAll();
        return false;
    }

    // Histogram base names for standard data/MC comparison
    std::vector<std::string> histNames = {
        "h2d_weta_cogx", "h2d_wphi_cogx", "h2d_et1", "h2d_et2", "h2d_et3", "h2d_et4",
        "h2d_e11_to_e33", "h2d_e17_to_e77", "h2d_e32_to_e35", "h2d_bdt", "h2d_npb_score"
    };

    int nCuts = 4; // cut0, cut1, cut2 (tight), cut3 (non-tight)

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

    // Helper lambda to overlay vertical error bars for MC histograms
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
        herr->SetBit(kCanDelete);
        herr->Draw("E0 SAME");
    };

    //==============================================================================
    // Precompute NPB normalization factors per (eta, pT) from DATA
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

            p_all->Rebin(4);
            p_npb->Rebin(4);
            p_all->GetXaxis()->SetRangeUser(0.0, 2.0);
            p_npb->GetXaxis()->SetRangeUser(0.0, 2.0);

            const int binStart = std::max(1, p_all->FindBin(1.400001));
            const double all_tot  = p_all->Integral(1, p_all->GetNbinsX());
            const double npb_tot  = p_npb->Integral(1, p_npb->GetNbinsX());
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
    // Part 0: Standard Data/MC Comparison Plots
    //==============================================================================
    std::cout << "Creating standard data/MC comparison plots for " << config_suffix << " ..." << std::endl;

    for (const auto &hbase : histNames)
    {
        for (int icut = 0; icut < nCuts; ++icut)
        {
            for (int ieta = 0; ieta < nEtaBins; ++ieta)
            {
                for (int ipt = 0; ipt < nPtBins; ++ipt)
                {
                    TString xaxisname = hbase.substr(4, hbase.size() - 4);
                    float xaxismin, xaxismax; int nrebin;
                    getAxisSettings(xaxisname, xaxismin, xaxismax, nrebin);

                    TString histNameFull = Form(
                        "%s_eta%d_pt%d_cut%d",
                        hbase.c_str(), ieta, ipt, icut);

                    TString histNamesave = Form(
                        "%s_eta%d_pt%d_cut%d",
                        xaxisname.Data(), ieta, ipt, icut);

                    TH2F *h2_data = dynamic_cast<TH2F *>(f_data->Get(histNameFull));
                    TH2F *h2_sig  = dynamic_cast<TH2F *>(f_sig->Get(histNameFull));
                    TH2F *h2_bkg  = dynamic_cast<TH2F *>(f_bkg_inc->Get(histNameFull));

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

                    TH1D *proj_data = h2_data->ProjectionX(Form("%s_px_data", histNameFull.Data()));
                    TH1D *proj_sig  = h2_sig->ProjectionX(Form("%s_px_sig",  histNameFull.Data()));
                    TH1D *proj_bkg  = h2_bkg->ProjectionX(Form("%s_px_bkg",  histNameFull.Data()));

                    // Optional: DATA NPB overlay for cut0
                    TH1D *proj_data_npb = nullptr;
                    if (icut == 0)
                    {
                        TString histNameNPBFull = Form("%s_eta%d_pt%d_cut4", hbase.c_str(), ieta, ipt);
                        TH2F *h2_data_npb = dynamic_cast<TH2F *>(f_data->Get(histNameNPBFull));
                        if (h2_data_npb)
                        {
                            h2_data_npb->RebinX(nrebin);
                            h2_data_npb->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
                            proj_data_npb = h2_data_npb->ProjectionX(Form("%s_px_data_npb", histNameNPBFull.Data()));
                        }
                        if (!proj_data_npb)
                        {
                            std::cerr << "Warning: Could not retrieve NPB histogram for " << histNameNPBFull << std::endl;
                        }
                    }

                    scaleToUnit(proj_data);
                    scaleToUnit(proj_sig);
                    scaleToUnit(proj_bkg);

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

                    // Prepare NPB overlay
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
                            delete proj_data_npb;
                            proj_data_npb = nullptr;
                        }
                    }

                    // Optional: plot Data - scaled NPB
                    if (makeDataMinusNPB && icut == 0 && proj_data_npb)
                    {
                        TH1D *proj_sub = (TH1D *)proj_data->Clone(Form("%s_px_data_minus_npb", histNameFull.Data()));
                        proj_sub->SetDirectory(nullptr);
                        proj_sub->Sumw2();
                        proj_sub->Add(proj_data_npb, -1.0);

                        double subMin = 0.0, subMax = 0.0;
                        bool subInit = false;
                        for (int ib = 1; ib <= proj_sub->GetNbinsX(); ++ib)
                        {
                            const double v = proj_sub->GetBinContent(ib);
                            if (!subInit) { subMin = v; subMax = v; subInit = true; }
                            else { if (v < subMin) subMin = v; if (v > subMax) subMax = v; }
                        }
                        if (!subInit) { subMin = -1e-6; subMax = 1e-6; }

                        const double maxy_cut0_like = std::max(
                            std::max((double)proj_data->GetMaximum(), std::max((double)proj_sig->GetMaximum(), (double)proj_bkg->GetMaximum())),
                            (proj_data_npb ? (double)proj_data_npb->GetMaximum() : 0.0));

                        double yMax = (maxy_cut0_like > 0.0) ? (1.3 * maxy_cut0_like) : 1e-6;
                        double yMin = 0.0;
                        if (subMin < 0.0) yMin = 1.1 * subMin;
                        if (yMax <= yMin) yMax = yMin + 1e-6;

                        TCanvas *c_sub = new TCanvas(
                            Form("c_sub_%s", histNameFull.Data()),
                            Form("Data - scaled NPB - %s", histNameFull.Data()),
                            600, 600);
                        c_sub->cd();

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
                        maxy = std::max(maxy, (float)proj_data_npb->GetMaximum());

                    TCanvas *c_proj = new TCanvas(
                        Form("c_proj_%s", histNameFull.Data()),
                        Form("ProjectionX - %s", histNameFull.Data()),
                        600, 600);
                    c_proj->cd();

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
                    if (icut == 0) bgcut = "w/o nbkg cut";
                    if (icut == 1) bgcut = "w/  nbkg cut";
                    if (icut == 2) bgcut = "w/ tight cut";
                    if (icut == 3) bgcut = "w/ nontight cut";
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

                    // Profile plot for background only
                    TProfile *pfx_bkg = h2_bkg->ProfileX(
                        Form("%s_pfx_bkg", histNameFull.Data()), 1, -1, "");

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
                    float pmax = proj_bkg_clone->GetMaximum();
                    proj_bkg_clone->GetYaxis()->SetRangeUser(0, pmax * 1.3);
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

    // Close and clean up files
    f_data->Close();
    delete f_data;
    f_sig->Close();
    delete f_sig;
    if (hasBkgOnly && f_bkg_only && f_bkg_only != f_bkg_inc)
    {
        f_bkg_only->Close();
        delete f_bkg_only;
    }
    f_bkg_inc->Close();
    delete f_bkg_inc;

    std::cout << "[plotOneConfig] Done: " << config_suffix << " -> " << savePath << "/" << std::endl;
    return true;
}

// ---------------------------------------------------------------------------
// Main entry point: discover all config_showershape*.yaml and process each.
// Pass only_suffix (e.g. "showershape_nom") to process a single variant.
// ---------------------------------------------------------------------------
void plot_showershapes_variations(const std::string &only_suffix = "")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    init_plot();

    // Discover all showershape config files
    TString raw = gSystem->GetFromPipe(
        "ls ../efficiencytool/config_showershape*.yaml 2>/dev/null");
    if (raw.IsNull())
    {
        std::cerr << "No config_showershape*.yaml files found in ../efficiencytool/" << std::endl;
        return;
    }

    TObjArray *tokens = raw.Tokenize("\n");
    int nProcessed = 0, nSkipped = 0;

    for (int i = 0; i < tokens->GetEntries(); ++i)
    {
        TString configpath = ((TObjString *)tokens->At(i))->GetString().Strip(TString::kBoth);
        if (configpath.IsNull()) continue;

        // Extract basename (e.g. "config_showershape_nom.yaml")
        std::string configname = std::string(gSystem->BaseName(configpath.Data()));

        // Derive suffix for filtering (strip path, "config_", ".yaml")
        std::string suffix = configname;
        {
            size_t slash = suffix.rfind('/');
            if (slash != std::string::npos) suffix = suffix.substr(slash + 1);
            if (suffix.size() > 5 && suffix.substr(suffix.size() - 5) == ".yaml")
                suffix = suffix.substr(0, suffix.size() - 5);
            if (suffix.size() > 7 && suffix.substr(0, 7) == "config_")
                suffix = suffix.substr(7);
        }

        // Apply suffix filter if requested
        if (!only_suffix.empty() && suffix != only_suffix) continue;

        std::cout << "\n=== Processing: " << configname << " ===" << std::endl;
        bool ok = plotOneConfig(configname);
        if (ok) ++nProcessed;
        else    ++nSkipped;
    }

    delete tokens;
    std::cout << "\nDone. Processed " << nProcessed << " config(s), skipped " << nSkipped << "." << std::endl;
}
