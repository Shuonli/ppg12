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

void plot_showershapes_selections(const std::string &configname = "config_showershape_0rad.yaml",
                                  bool   combine_double     = true,
                                  double sig_double_weight  = 0.75,
                                  double bkg_double_weight  = 0.25)
{
    init_plot();
    string savePath = "../PPG12-analysis-note/Figures/showershapes_selections/";
    gSystem->Exec(Form("mkdir -p %s", savePath.c_str()));

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
    std::cout << "Config suffix: " << config_suffix << std::endl;

    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node config = YAML::LoadFile(("../efficiencytool/" + configname).c_str());

    const bool makeDataMinusNPB = true;

    std::vector<double> bdt_bins_default = {0.0, 0.3, 0.7, 1.0};
    std::vector<double> bdt_bins     = config["analysis"]["bdt_bins"].as<std::vector<double>>(bdt_bins_default);
    std::vector<double> pT_bin_edges = config["analysis"]["pT_bins"].as<std::vector<double>>();
    std::vector<double> eta_bins     = {-0.7, 0.7};

    int nBdtBins = bdt_bins.size() - 1;
    int nPtBins  = pT_bin_edges.size() - 1;
    int nEtaBins = eta_bins.size() - 1;

    std::cout << "Using " << nBdtBins << " BDT bins from config file" << std::endl;
    std::cout << "Using " << nPtBins  << " pT bins from config file"  << std::endl;

    // Open files
    const std::string resultsDir      = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/";
    const std::string dataFile        = resultsDir + "data_histoshower_shape_" + config_suffix + ".root";
    const std::string sigFile         = resultsDir + "MC_efficiencyshower_shape_signal_" + config_suffix + ".root";
    const std::string bkgInclusiveFile= resultsDir + "MC_efficiencyshower_shape_jet_inclusive_" + config_suffix + ".root";
    const std::string bkgOnlyFile     = resultsDir + "MC_efficiencyshower_shape_jet_" + config_suffix + ".root";
    const std::string sigDoubleFile   = resultsDir + "MC_efficiencyshower_shape_photon10_double_showershape.root";
    const std::string bkgIncDoubleFile= resultsDir + "MC_efficiencyshower_shape_jet12_double_inclusive_showershape.root";

    TFile *f_data    = TFile::Open(dataFile.c_str(),         "READ");
    TFile *f_sig     = TFile::Open(sigFile.c_str(),          "READ");
    TFile *f_bkg_inc = TFile::Open(bkgInclusiveFile.c_str(), "READ");

    TFile *f_bkg_only = TFile::Open(bkgOnlyFile.c_str(), "READ");
    bool hasBkgOnly   = (f_bkg_only && !f_bkg_only->IsZombie());
    if (!hasBkgOnly) f_bkg_only = f_bkg_inc;

    TFile *f_sig_double    = nullptr;
    TFile *f_bkginc_double = nullptr;
    if (combine_double)
    {
        f_sig_double    = TFile::Open(sigDoubleFile.c_str(),    "READ");
        f_bkginc_double = TFile::Open(bkgIncDoubleFile.c_str(), "READ");
        if (!f_sig_double    || f_sig_double->IsZombie())
        { std::cerr << "Warning: signal double file not found; combine_double disabled for signal.\n"; f_sig_double = nullptr; }
        if (!f_bkginc_double || f_bkginc_double->IsZombie())
        { std::cerr << "Warning: bkg double file not found; combine_double disabled for bkg.\n";    f_bkginc_double = nullptr; }
    }

    if (!f_data    || f_data->IsZombie())    { std::cerr << "Error: Could not open data file!\n";               return; }
    if (!f_sig     || f_sig->IsZombie())     { std::cerr << "Error: Could not open signal file!\n";             return; }
    if (!f_bkg_inc || f_bkg_inc->IsZombie()) { std::cerr << "Error: Could not open inclusive background file!\n"; return; }

    // Histogram base names for standard data/MC comparison (reused in Part 5)
    std::vector<std::string> histNames = {
        "h2d_weta_cogx", "h2d_wphi_cogx", "h2d_wr", "h2d_et1", "h2d_et2", "h2d_et3", "h2d_et4",
        "h2d_e11_to_e33", "h2d_e17_to_e77", "h2d_e32_to_e35", "h2d_bdt", "h2d_npb_score"
    };
    int nCuts = 4; // cut0, cut1, cut2 (tight), cut3 (non-tight)

    // Color palette for pT bins (Parts 6, 7)
    std::vector<Color_t> pt_colors = {kBlue, kRed, kGreen+2, kMagenta+1, kOrange+7,
                                      kCyan+2, kViolet+1, kTeal+2, kPink+1, kAzure+2};

    //--------------------------------------------------------------------------
    // Core helpers
    //--------------------------------------------------------------------------

    // getCombined: owned clone; (1-w)*nom + w*dbl when combine_double active.
    auto getCombined = [&](TFile* f_nominal, TFile* f_double, double w,
                           const char* histname) -> TH1* {
        TH1 *h_nom = dynamic_cast<TH1*>(f_nominal->Get(histname));
        if (!h_nom) return nullptr;
        TH1 *result = dynamic_cast<TH1*>(h_nom->Clone(Form("%s_comb", histname)));
        result->SetDirectory(nullptr);
        if (combine_double && f_double && w != 0.0) {
            TH1 *h_dbl = dynamic_cast<TH1*>(f_double->Get(histname));
            if (h_dbl) { result->Scale(1.0 - w); result->Add(h_dbl, w); }
        }
        return result;
    };

    // Fetchers: fetch_data → borrowed; fetch_sig / fetch_bkg_inc → owned clones.
    using Fetcher = std::function<TH1*(const char*)>;
    auto fetch_data    = [&](const char* n) -> TH1* { return dynamic_cast<TH1*>(f_data->Get(n)); };
    auto fetch_sig     = [&](const char* n) -> TH1* { return getCombined(f_sig,     f_sig_double,    sig_double_weight, n); };
    auto fetch_bkg_inc = [&](const char* n) -> TH1* { return getCombined(f_bkg_inc, f_bkginc_double, bkg_double_weight, n); };

    // all_samples: used in Parts 4, 7
    std::vector<std::tuple<std::string, Fetcher, Color_t, std::string>> all_samples = {
        {"data",          fetch_data,    kBlack, "Data"},
        {"signal",        fetch_sig,     kRed,   "Signal MC"},
        {"bkg_inclusive", fetch_bkg_inc, kBlue,  "Inclusive MC"}
    };

    // Axis range + rebin for each shower-shape variable
    auto getAxisSettings = [](const TString& xaxisname, float& xmin, float& xmax, int& nrebin) {
        xmin = 0.0; xmax = 1.0; nrebin = 4;
        if (xaxisname[0] == 'w') { xmax = 2.0; }
        if (xaxisname[0] == 'h') { xmax = 0.3; nrebin = 1; }
        if (xaxisname == "et4")          { xmax = 0.3; nrebin = 1; }
        if (xaxisname == "e32_to_e35")   { xmin = 0.4; xmax = 1.0; nrebin = 1; }
        if (xaxisname == "et1")          { xmin = 0.3; xmax = 1.0; nrebin = 1; }
        if (xaxisname == "bdt" || xaxisname == "npb_score") { xmin = 0.0; xmax = 1.0; nrebin = 2; }
    };

    // Scale histogram to unit area
    auto scaleToUnit = [](TH1D *h) {
        if (!h) return;
        double i = h->Integral();
        if (i > 1e-12) h->Scale(1.0 / i);
    };

    // Draw error bars on MC histogram (clone with kCanDelete so ROOT cleans it up)
    auto overlayVerticalErrors = [](TH1 *h) {
        if (!h) return;
        TH1 *herr = dynamic_cast<TH1 *>(h->Clone(Form("%s_err", h->GetName())));
        if (!herr) return;
        herr->SetDirectory(nullptr);
        herr->SetFillStyle(0);
        herr->SetMarkerStyle(1); herr->SetMarkerSize(0);
        herr->SetLineColor(h->GetLineColor()); herr->SetLineWidth(h->GetLineWidth());
        herr->SetMarkerColor(h->GetLineColor());
        herr->SetBit(kCanDelete);
        herr->Draw("E0 SAME");
    };

    //--------------------------------------------------------------------------
    // Convenience helpers (eliminate repetitive boilerplate)
    //--------------------------------------------------------------------------

    // Draw the standard 3-line sPHENIX / run header
    auto drawLabels = [&](float x, float y) {
        myText(x, y,       1, strleg1.c_str(), 0.04);
        myText(x, y-0.05f, 1, strleg2.c_str(), 0.04);
        myText(x, y-0.10f, 1, strleg3.c_str(), 0.04);
    };

    // Create a legend with standard no-border style
    auto makeLegend = [](float x1, float y1, float x2, float y2, float sz = 0.04) -> TLegend* {
        TLegend *l = new TLegend(x1, y1, x2, y2);
        l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(sz);
        return l;
    };

    // Style a histogram as an MC line
    auto styleMC = [](TH1 *h, Color_t col) {
        h->SetLineColor(col); h->SetMarkerColor(col); h->SetLineWidth(2); h->SetStats(0);
    };

    // Style a histogram as data markers
    auto styleData = [](TH1 *h, int mstyle = 20, float msize = 1.0f) {
        h->SetLineColor(kBlack); h->SetMarkerColor(kBlack);
        h->SetMarkerStyle(mstyle); h->SetMarkerSize(msize); h->SetStats(0);
    };

    // Compute efficiency-vs-threshold histogram from a 1D NPB score distribution.
    // eff(bin) = integral(bin → N) / total.  Error bars always set.
    auto computeEffHist = [](TH1D *h_npb, const char *name) -> TH1D* {
        int    nbins = h_npb->GetNbinsX();
        double total = h_npb->Integral(1, nbins);
        TH1D  *h_eff = new TH1D(name, "",
                                 nbins, h_npb->GetXaxis()->GetXmin(), h_npb->GetXaxis()->GetXmax());
        h_eff->SetDirectory(nullptr);
        for (int i = 1; i <= nbins; ++i) {
            double eff = (total > 0) ? h_npb->Integral(i, nbins) / total : 0.0;
            h_eff->SetBinContent(i, eff);
            h_eff->SetBinError(i, sqrt(eff * (1.0 - eff) / std::max(total, 1.0)));
        }
        return h_eff;
    };

    // Compute NPB ratio histogram N(0.5<NPB<0.9) / N(0.9<NPB<1.0) vs MBD sigma (Y-axis of h2).
    auto computeNPBRatioHist = [](TH2D *h2, const char *name) -> TH1D* {
        int xbin_lo   = h2->GetXaxis()->FindBin(0.5 + 1e-6);
        int xbin_mid  = h2->GetXaxis()->FindBin(0.9 - 1e-6);
        int xbin_hi   = h2->GetXaxis()->FindBin(0.9 + 1e-6);
        int xbin_max  = h2->GetXaxis()->FindBin(1.0 - 1e-6);
        int nbins_y   = h2->GetNbinsY();
        TH1D *h_ratio = new TH1D(name, "", nbins_y,
                                  h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax());
        h_ratio->SetDirectory(nullptr);
        for (int iy = 1; iy <= nbins_y; ++iy) {
            double n_mid = 0.0, n_high = 0.0;
            for (int ix = xbin_lo;  ix <= xbin_mid; ++ix) n_mid  += h2->GetBinContent(ix, iy);
            for (int ix = xbin_hi;  ix <= xbin_max; ++ix) n_high += h2->GetBinContent(ix, iy);
            double ratio = 0.0, ratio_err = 0.0;
            if (n_high > 0) {
                ratio = n_mid / n_high;
                ratio_err = (n_mid >= 1) ? ratio * sqrt(1.0/n_mid + 1.0/n_high) : 0.0;
            }
            h_ratio->SetBinContent(iy, ratio);
            h_ratio->SetBinError(iy, ratio_err);
        }
        return h_ratio;
    };

    //--------------------------------------------------------------------------
    // Precompute NPB normalization factors per (eta, pT) from DATA
    // f_NPB = (data tail fraction) / (NPB-tagged tail fraction), weta > 1
    //--------------------------------------------------------------------------
    std::vector<std::vector<double>> npbScale(nEtaBins, std::vector<double>(nPtBins, 0.0));
    for (int ieta = 0; ieta < nEtaBins; ++ieta)
    {
        for (int ipt = 0; ipt < nPtBins; ++ipt)
        {
            TString h_weta_all_name = Form("h2d_weta_cogx_eta%d_pt%d_cut0", ieta, ipt);
            TString h_weta_npb_name = Form("h2d_weta_cogx_eta%d_pt%d_cut4", ieta, ipt);
            TH2F *h2_weta_all = dynamic_cast<TH2F *>(f_data->Get(h_weta_all_name));
            TH2F *h2_weta_npb = dynamic_cast<TH2F *>(f_data->Get(h_weta_npb_name));
            if (!h2_weta_all || !h2_weta_npb) {
                std::cerr << "[NPB scale] Missing hist(s): "
                          << h_weta_all_name << " or " << h_weta_npb_name << std::endl;
                continue;
            }

            TH1D *p_all = h2_weta_all->ProjectionX(Form("%s_px_all", h_weta_all_name.Data()));
            TH1D *p_npb = h2_weta_npb->ProjectionX(Form("%s_px_npb", h_weta_npb_name.Data()));
            if (!p_all || !p_npb) {
                std::cerr << "[NPB scale] Failed ProjectionX for: "
                          << h_weta_all_name << " or " << h_weta_npb_name << std::endl;
                continue;
            }

            p_all->Rebin(4); p_npb->Rebin(4);
            p_all->GetXaxis()->SetRangeUser(0.0, 2.0);
            p_npb->GetXaxis()->SetRangeUser(0.0, 2.0);

            const int    binStart  = std::max(1, p_all->FindBin(1.400001));
            const double all_tot   = p_all->Integral(1, p_all->GetNbinsX());
            const double npb_tot   = p_npb->Integral(1, p_npb->GetNbinsX());
            const double all_tail  = p_all->Integral(binStart, p_all->GetNbinsX());
            const double npb_tail  = p_npb->Integral(binStart, p_npb->GetNbinsX());

            if (all_tot <= 0 || npb_tot <= 0 || all_tail < 0 || npb_tail <= 0) {
                std::cerr << "[NPB scale] Bad integrals for (eta,pt)=(" << ieta << "," << ipt << ")"
                          << " all_tot=" << all_tot << " all_tail=" << all_tail
                          << " npb_tot=" << npb_tot << " npb_tail=" << npb_tail << std::endl;
            } else {
                const double frac_all = all_tail / all_tot;
                const double frac_npb = npb_tail / npb_tot;
                const double f_npb    = (frac_npb > 1e-12) ? (frac_all / frac_npb) : 0.0;
                if (!(f_npb >= 0.0) || f_npb > 1.2)
                    std::cerr << "[NPB scale] Unphysical f_NPB=" << f_npb
                              << " for (eta,pt)=(" << ieta << "," << ipt << ")"
                              << " frac_all=" << frac_all << " frac_npb=" << frac_npb << std::endl;
                else
                    npbScale[ieta][ipt] = f_npb;
            }
            delete p_all; delete p_npb;
        }
    }

    // Cut labels used in Part 0
    static const std::vector<std::string> bgcut_labels = {
        "w/o nbkg cut", "w/  nbkg cut", "w/ tight cut", "w/ nontight cut"
    };

    //==========================================================================
    // Part 0: Standard Data/MC Comparison Plots
    //==========================================================================
    std::cout << "Creating standard data/MC comparison plots..." << std::endl;

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

                    TString histNameFull = Form("%s_eta%d_pt%d_cut%d", hbase.c_str(), ieta, ipt, icut);
                    TString histNamesave = Form("%s_eta%d_pt%d_cut%d", xaxisname.Data(), ieta, ipt, icut);

                    // h2_data is borrowed; h2_sig/h2_bkg are owned clones
                    TH2F *h2_data = dynamic_cast<TH2F *>(f_data->Get(histNameFull));
                    TH2F *h2_sig  = dynamic_cast<TH2F *>(getCombined(f_sig,     f_sig_double,    sig_double_weight, histNameFull));
                    TH2F *h2_bkg  = dynamic_cast<TH2F *>(getCombined(f_bkg_inc, f_bkginc_double, bkg_double_weight, histNameFull));

                    if (!h2_data || !h2_sig || !h2_bkg) {
                        std::cerr << "Warning: Could not retrieve histograms for "
                                  << histNameFull << " from one or more files!\n";
                        delete h2_sig; delete h2_bkg;
                        continue;
                    }

                    for (TH2F *h : {h2_bkg, h2_sig, h2_data}) {
                        h->RebinX(nrebin);
                        h->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
                    }

                    TH1D *proj_data = h2_data->ProjectionX(Form("%s_px_data", histNameFull.Data()));
                    TH1D *proj_sig  = h2_sig->ProjectionX( Form("%s_px_sig",  histNameFull.Data()));
                    TH1D *proj_bkg  = h2_bkg->ProjectionX( Form("%s_px_bkg",  histNameFull.Data()));

                    // Optional NPB overlay for cut0 (cut4 = NPB-tagged data)
                    TH1D *proj_data_npb = nullptr;
                    if (icut == 0) {
                        TString npbName = Form("%s_eta%d_pt%d_cut4", hbase.c_str(), ieta, ipt);
                        TH2F *h2_npb = dynamic_cast<TH2F *>(f_data->Get(npbName));
                        if (h2_npb) {
                            h2_npb->RebinX(nrebin);
                            h2_npb->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
                            proj_data_npb = h2_npb->ProjectionX(Form("%s_px_data_npb", npbName.Data()));
                        }
                        if (!proj_data_npb)
                            std::cerr << "Warning: Could not retrieve NPB histogram for " << npbName << std::endl;
                    }

                    scaleToUnit(proj_data); scaleToUnit(proj_sig); scaleToUnit(proj_bkg);
                    styleData(proj_data);
                    styleMC(proj_sig, kRed);
                    styleMC(proj_bkg, kBlue);

                    // Scale NPB overlay by tail-matched fraction
                    if (proj_data_npb) {
                        proj_data_npb->Sumw2();
                        scaleToUnit(proj_data_npb);
                        const double s = npbScale[ieta][ipt];
                        if (s > 0.0) {
                            proj_data_npb->Scale(s);
                        } else {
                            std::cout << "Warning: No valid scale for NPB overlay for "
                                      << Form("%s_eta%d_pt%d_cut4", hbase.c_str(), ieta, ipt)
                                      << " s=" << s << std::endl;
                            delete proj_data_npb; proj_data_npb = nullptr;
                        }
                    }

                    // Compute shared y-axis maximum (used by both c_sub and c_proj)
                    float maxy = std::max({proj_data->GetMaximum(), proj_sig->GetMaximum(), proj_bkg->GetMaximum()});
                    if (proj_data_npb) maxy = std::max(maxy, (float)proj_data_npb->GetMaximum());

                    float pTlow  = pT_bin_edges[ipt];
                    float pThigh = pT_bin_edges[ipt + 1];

                    // --- Data minus scaled-NPB canvas (cut0 only) ---
                    if (makeDataMinusNPB && icut == 0 && proj_data_npb)
                    {
                        TH1D *proj_sub = (TH1D *)proj_data->Clone(Form("%s_px_data_minus_npb", histNameFull.Data()));
                        proj_sub->SetDirectory(nullptr);
                        proj_sub->Sumw2();
                        proj_sub->Add(proj_data_npb, -1.0);

                        double subMin = proj_sub->GetMinimum();
                        double yMin   = (subMin < 0.0) ? 1.1 * subMin : 0.0;
                        double yMax   = std::max(1.3 * (double)maxy, yMin + 1e-6);

                        TCanvas *c_sub = new TCanvas(Form("c_sub_%s", histNameFull.Data()),
                                                     Form("Data - scaled NPB - %s", histNameFull.Data()), 600, 600);
                        c_sub->cd();
                        proj_sig->SetTitle(""); proj_sig->SetYTitle("normalized counts");
                        proj_sig->GetYaxis()->SetTitleOffset(1.5);
                        proj_sig->SetXTitle(xaxisname.Data());
                        proj_sig->GetXaxis()->SetNdivisions(505);
                        proj_sig->GetYaxis()->SetRangeUser(yMin, yMax);
                        proj_sig->Draw("HIST"); proj_bkg->Draw("HIST SAME");
                        overlayVerticalErrors(proj_sig); overlayVerticalErrors(proj_bkg);

                        proj_sub->SetStats(0);
                        proj_sub->SetLineColor(kBlack); proj_sub->SetMarkerColor(kBlack);
                        proj_sub->SetMarkerStyle(proj_data->GetMarkerStyle());
                        proj_sub->SetMarkerSize(proj_data->GetMarkerSize());
                        proj_sub->Draw("ex0 SAME");

                        TLine l0(xaxismin, 0.0, xaxismax, 0.0);
                        l0.SetLineStyle(2); l0.SetLineColor(kGray + 2); l0.Draw("SAME");

                        drawLabels(0.20, 0.90);
                        myText(0.20, 0.75, 1, Form("%.0f<p_{T}<%.0fGeV, cut0", pTlow, pThigh), 0.04);
                        myMarkerLineText(0.55, 0.90, 1.5, kBlack, 20, kBlack, 1, "Data - scaled NPB", 0.045, true);
                        myMarkerLineText(0.55, 0.85, 0,   kRed,    0, kRed,   1, "Signal MC",         0.045, true);
                        myMarkerLineText(0.55, 0.80, 0,   kBlue,   0, kBlue,  1, "Inclusive MC",      0.045, true);

                        c_sub->SaveAs(Form("%s/dis_subNPB_%s.pdf", savePath.c_str(), histNamesave.Data()));
                        delete c_sub; delete proj_sub;
                    }

                    // --- Main data/MC overlay canvas ---
                    TCanvas *c_proj = new TCanvas(Form("c_proj_%s", histNameFull.Data()),
                                                  Form("ProjectionX - %s", histNameFull.Data()), 600, 600);
                    c_proj->cd();
                    proj_sig->SetYTitle("normalized counts");
                    proj_sig->GetYaxis()->SetTitleOffset(1.5);
                    proj_sig->SetXTitle(xaxisname.Data());
                    proj_sig->GetYaxis()->SetRangeUser(0, maxy * 1.3);
                    proj_sig->GetXaxis()->SetNdivisions(505);
                    proj_sig->SetStats(0);
                    proj_sig->Draw("HIST"); proj_bkg->Draw("HIST SAME");
                    overlayVerticalErrors(proj_sig); overlayVerticalErrors(proj_bkg);
                    proj_data->Draw("ex0 SAME");
                    if (proj_data_npb) {
                        styleMC(proj_data_npb, kGreen + 2);
                        proj_data_npb->Draw("hist SAME");
                        overlayVerticalErrors(proj_data_npb);
                    }

                    drawLabels(0.20, 0.90);
                    myText(0.2, 0.75, 1, Form("%.0f<p_{T}<%.0fGeV,%s",
                                              pTlow, pThigh, bgcut_labels[icut].c_str()), 0.04);
                    myMarkerLineText(0.6, 0.90, 1.5, kBlack, 20, kBlack, 1, "Data",        0.05, true);
                    myMarkerLineText(0.6, 0.85, 0,   kRed,    0, kRed,   1, "Signal MC",   0.05, true);
                    myMarkerLineText(0.6, 0.80, 0,   kBlue,   0, kBlue,  1, "Inclusive MC",0.05, true);
                    if (icut == 0)
                        myMarkerLineText(0.6, 0.75, 0, kGreen+2, 0, kGreen+2, 2, "Data NPB", 0.05, true);

                    c_proj->SaveAs(Form("%s/dis_%s.pdf", savePath.c_str(), histNamesave.Data()));

                    // --- Background profile vs ET_iso ---
                    TProfile *pfx_bkg = h2_bkg->ProfileX(Form("%s_pfx_bkg", histNameFull.Data()), 1, -1, "");
                    TCanvas  *c_bkg   = new TCanvas(Form("c_bkg_%s", histNameFull.Data()),
                                                    Form("Bkg only - %s", histNameFull.Data()), 600, 600);
                    c_bkg->cd();
                    TH1D *h_bkg_profile = (TH1D *)pfx_bkg->Clone(Form("%s_px_bkgOnly", histNameFull.Data()));
                    h_bkg_profile->SetLineColor(kBlue);
                    h_bkg_profile->SetYTitle("<#it{E}_{T}^{iso}> [GeV]");
                    h_bkg_profile->SetXTitle(xaxisname.Data());
                    h_bkg_profile->SetStats(0);
                    h_bkg_profile->GetYaxis()->SetRangeUser(0, h_bkg_profile->GetMaximum() * 1.3);
                    h_bkg_profile->Draw("HIST");
                    overlayVerticalErrors(h_bkg_profile);

                    drawLabels(0.20, 0.90);
                    myText(0.55, 0.90, 1, Form("%.0f<p_{T}<%.0fGeV,%s", pTlow, pThigh, bgcut_labels[icut].c_str()), 0.04);
                    myText(0.55, 0.85, 1, "Background MC", 0.04);
                    myText(0.55, 0.80, 1, Form("Correlation: %.3f", h2_bkg->GetCorrelationFactor()), 0.04);

                    c_bkg->SaveAs(Form("%s/pfx_%s.pdf", savePath.c_str(), histNamesave.Data()));

                    delete h2_sig; delete h2_bkg;
                    delete pfx_bkg; delete h_bkg_profile;
                    delete c_proj; delete c_bkg;
                    if (proj_data_npb) delete proj_data_npb;
                } // ipt
            }     // ieta
        }         // icut
    }             // histNames

    //==========================================================================
    // Part 1: BDT Overlay Plots (isolation ET binned in BDT score)
    //==========================================================================
    std::cout << "Creating BDT overlay plots..." << std::endl;

    {
        TH1 *test = getCombined(f_sig, f_sig_double, sig_double_weight, "h_reco_iso_eta0_pt0_bdt0");
        bool hasBdtHists = (test != nullptr);
        delete test;
        if (!hasBdtHists) {
            std::cout << "Warning: BDT-binned histograms not found; skipping Part 1.\n";
        } else {
            Color_t baseColors[] = {kRed, kBlue, kGreen+2, kMagenta+1, kOrange+7, kCyan+2, kViolet+1, kTeal+2};
            const int nBaseColors = sizeof(baseColors) / sizeof(baseColors[0]);

            for (int ieta = 0; ieta < nEtaBins; ieta++)
            for (int ipt  = 0; ipt  < nPtBins;  ipt++)
            {
                std::vector<TH1D *> proj_bdt(nBdtBins, nullptr);
                float maxY = 0;
                int   nFound = 0;

                // Collect + style all BDT-bin histograms first, then draw
                for (int ibdt = 0; ibdt < nBdtBins; ibdt++) {
                    TString hname = Form("h_reco_iso_eta%d_pt%d_bdt%d", ieta, ipt, ibdt);
                    TH1 *tmp = getCombined(f_sig, f_sig_double, sig_double_weight, hname);
                    if (!tmp) continue;
                    proj_bdt[ibdt] = (TH1D *)tmp->Clone(Form("%s_clone", hname.Data()));
                    delete tmp;
                    proj_bdt[ibdt]->Rebin(4);
                    if (proj_bdt[ibdt]->Integral() > 0) {
                        proj_bdt[ibdt]->Scale(1.0 / proj_bdt[ibdt]->Integral());
                        proj_bdt[ibdt]->Scale(1.0 / proj_bdt[ibdt]->GetBinWidth(1));
                    }
                    styleMC(proj_bdt[ibdt], baseColors[ibdt % nBaseColors]);
                    maxY = std::max(maxY, (float)proj_bdt[ibdt]->GetMaximum());
                    ++nFound;
                }
                if (nFound == 0) continue;

                TCanvas *c = new TCanvas(Form("c_bdt_iso_eta%d_pt%d", ieta, ipt), "", 600, 600);
                TLegend *leg = makeLegend(0.55, 0.70, 0.88, 0.88);

                // Draw with correct max Y known up front
                bool first = true;
                for (int ibdt = 0; ibdt < nBdtBins; ibdt++) {
                    if (!proj_bdt[ibdt]) continue;
                    if (first) {
                        proj_bdt[ibdt]->SetTitle("");
                        proj_bdt[ibdt]->GetXaxis()->SetTitle("Iso E_{T} [GeV]");
                        proj_bdt[ibdt]->GetXaxis()->SetRangeUser(-3, 15);
                        proj_bdt[ibdt]->GetXaxis()->SetNdivisions(505);
                        proj_bdt[ibdt]->GetYaxis()->SetTitle("Normalized counts / GeV");
                        proj_bdt[ibdt]->GetYaxis()->SetTitleOffset(1.5);
                        proj_bdt[ibdt]->GetYaxis()->SetRangeUser(0, maxY * 1.3);
                        proj_bdt[ibdt]->Draw("HIST");
                        first = false;
                    } else {
                        proj_bdt[ibdt]->Draw("HIST SAME");
                    }
                    overlayVerticalErrors(proj_bdt[ibdt]);
                    leg->AddEntry(proj_bdt[ibdt], Form("%.2f < BDT < %.2f", bdt_bins[ibdt], bdt_bins[ibdt+1]), "l");
                }

                drawLabels(0.20, 0.90);
                myText(0.20, 0.75, 1, strSigMC.c_str(), 0.04);
                myText(0.20, 0.65, 1, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), 0.04);
                leg->Draw();

                c->SaveAs(Form("%s/isoET_reco_BDTbin_overlay_eta%d_pt%d.pdf", savePath.c_str(), ieta, ipt));
                for (auto h : proj_bdt) if (h) delete h;
                delete c; delete leg;
            }
        }
    }

    //==========================================================================
    // Part 2: Tight vs Non-tight (shower-shape sideband) Overlay Plots
    //==========================================================================
    std::cout << "Creating isolation region overlay plots..." << std::endl;

    std::vector<std::string> histNames2D = { "h2d_npb_score" };

    for (const auto &hbase : histNames2D)
    for (int ieta = 0; ieta < nEtaBins; ieta++)
    for (int ipt  = 0; ipt  < nPtBins;  ipt++)
    {
        TString xaxisname = hbase.substr(4, hbase.size() - 4);
        float xaxismin, xaxismax; int nrebin;
        getAxisSettings(xaxisname, xaxismin, xaxismax, nrebin);

        TString hname_tight = Form("%s_eta%d_pt%d_cut2", hbase.c_str(), ieta, ipt);
        TString hname_nt    = Form("%s_eta%d_pt%d_cut3", hbase.c_str(), ieta, ipt);

        TH2F *h2_mc_tight   = dynamic_cast<TH2F*>(getCombined(f_bkg_inc, f_bkginc_double, bkg_double_weight, hname_tight));
        TH2F *h2_mc_nt      = dynamic_cast<TH2F*>(getCombined(f_bkg_inc, f_bkginc_double, bkg_double_weight, hname_nt));
        TH2F *h2_data_tight = dynamic_cast<TH2F*>(f_data->Get(hname_tight));
        TH2F *h2_data_nt    = dynamic_cast<TH2F*>(f_data->Get(hname_nt));

        if (!h2_mc_tight || !h2_mc_nt || !h2_data_tight || !h2_data_nt) {
            std::cerr << "Warning: Could not retrieve histograms for " << hbase << std::endl;
            delete h2_mc_tight; delete h2_mc_nt;
            continue;
        }

        TH1D *proj_mc_tight   = h2_mc_tight->ProjectionX(Form("%s_px_mc_tight",   hname_tight.Data()));
        TH1D *proj_mc_nt      = h2_mc_nt->ProjectionX(   Form("%s_px_mc_nt",      hname_nt.Data()));
        TH1D *proj_data_tight = h2_data_tight->ProjectionX(Form("%s_px_data_tight",hname_tight.Data()));
        TH1D *proj_data_nt    = h2_data_nt->ProjectionX(  Form("%s_px_data_nt",    hname_nt.Data()));

        if (nrebin > 1) {
            proj_mc_tight->Rebin(nrebin); proj_mc_nt->Rebin(nrebin);
            proj_data_tight->Rebin(nrebin); proj_data_nt->Rebin(nrebin);
        }
        for (TH1D *h : {proj_mc_tight, proj_mc_nt, proj_data_tight, proj_data_nt})
            h->GetXaxis()->SetRangeUser(xaxismin, xaxismax);

        scaleToUnit(proj_mc_tight); scaleToUnit(proj_mc_nt);
        scaleToUnit(proj_data_tight); scaleToUnit(proj_data_nt);

        styleMC(proj_mc_tight, kBlue);
        styleMC(proj_mc_nt,    kRed);
        styleData(proj_data_tight, 20, 0.9f);
        proj_data_nt->SetLineColor(kGray+2); proj_data_nt->SetMarkerColor(kGray+2);
        proj_data_nt->SetMarkerStyle(24); proj_data_nt->SetMarkerSize(0.9); proj_data_nt->SetStats(0);

        float maxY = std::max({proj_mc_tight->GetMaximum(), proj_mc_nt->GetMaximum(),
                               proj_data_tight->GetMaximum(), proj_data_nt->GetMaximum()});

        TCanvas *c = new TCanvas(Form("c_iso_%s_eta%d_pt%d", hbase.c_str(), ieta, ipt), "", 600, 600);
        proj_mc_tight->SetTitle("");
        proj_mc_tight->GetXaxis()->SetTitle(xaxisname.Data());
        proj_mc_tight->GetXaxis()->SetNdivisions(505);
        proj_mc_tight->GetYaxis()->SetTitle("Normalized counts");
        proj_mc_tight->GetYaxis()->SetTitleOffset(1.5);
        proj_mc_tight->GetYaxis()->SetRangeUser(0, maxY * 1.3);
        proj_mc_tight->Draw("HIST"); proj_mc_nt->Draw("HIST SAME");
        overlayVerticalErrors(proj_mc_tight); overlayVerticalErrors(proj_mc_nt);
        proj_data_tight->Draw("E SAME"); proj_data_nt->Draw("E SAME");

        drawLabels(0.20, 0.90);
        myText(0.20, 0.80, 1, "Data + Inclusive MC", 0.04);
        myText(0.20, 0.75, 1, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), 0.04);

        TLegend *leg = makeLegend(0.55, 0.70, 0.88, 0.88);
        leg->AddEntry(proj_mc_tight,   "Inclusive MC (tight)",    "l");
        leg->AddEntry(proj_mc_nt,      "Inclusive MC (non-tight)","l");
        leg->AddEntry(proj_data_tight, "Data (tight)",            "lep");
        leg->AddEntry(proj_data_nt,    "Data (non-tight)",        "lep");
        leg->Draw();

        c->SaveAs(Form("%s/isoRegionOverlay_%s_tight_vs_nontight_incMC_withData_eta%d_pt%d.pdf",
                       savePath.c_str(), xaxisname.Data(), ieta, ipt));
        delete proj_mc_tight; delete proj_mc_nt; delete proj_data_tight; delete proj_data_nt;
        delete h2_mc_tight; delete h2_mc_nt; delete c; delete leg;
    }

    //==========================================================================
    // Part 3: Data vs MC Isolation ET Comparison
    //==========================================================================
    std::cout << "Creating data vs MC isolation ET comparison plots..." << std::endl;

    for (int ieta = 0; ieta < nEtaBins; ieta++)
    for (int ipt  = 0; ipt  < nPtBins;  ipt++)
    {
        TH1D *h_data = dynamic_cast<TH1D*>(f_data->Get(Form("h_tight_reco_iso_eta%d_pt%d", ieta, ipt)));
        TH1D *h_mc   = dynamic_cast<TH1D*>(getCombined(f_bkg_inc, f_bkginc_double, bkg_double_weight,
                                                        Form("h_tight_reco_iso_eta%d_pt%d", ieta, ipt)));
        if (!h_data || !h_mc) {
            std::cerr << "Warning: Could not retrieve iso ET histograms for eta" << ieta << " pt" << ipt << std::endl;
            delete h_mc; continue;
        }

        TH1D *h_data_c = (TH1D *)h_data->Clone(Form("h_data_clone_eta%d_pt%d", ieta, ipt));
        TH1D *h_mc_c   = (TH1D *)h_mc->Clone(  Form("h_mc_clone_eta%d_pt%d",   ieta, ipt));
        delete h_mc;

        h_data_c->Rebin(4); h_mc_c->Rebin(4);
        if (h_data_c->Integral() > 0) h_data_c->Scale(1.0 / h_data_c->Integral());
        if (h_mc_c->Integral()   > 0) h_mc_c->Scale(  1.0 / h_mc_c->Integral());
        h_data_c->Scale(1.0 / h_data_c->GetBinWidth(1));
        h_mc_c->Scale(  1.0 / h_mc_c->GetBinWidth(1));

        styleData(h_data_c);
        styleMC(h_mc_c, kRed);

        float maxY = std::max(h_data_c->GetMaximum(), h_mc_c->GetMaximum());

        TCanvas *c = new TCanvas(Form("c_isoET_datamc_eta%d_pt%d", ieta, ipt), "", 600, 600);
        h_mc_c->SetTitle("");
        h_mc_c->GetXaxis()->SetTitle("Iso E_{T} [GeV]");
        h_mc_c->GetXaxis()->SetNdivisions(505);
        h_mc_c->GetYaxis()->SetTitle("Normalized counts / GeV");
        h_mc_c->GetYaxis()->SetTitleOffset(1.5);
        h_mc_c->GetYaxis()->SetRangeUser(0, maxY * 1.3);
        h_mc_c->GetXaxis()->SetRangeUser(-3, 15);
        h_mc_c->Draw("HIST"); overlayVerticalErrors(h_mc_c);
        h_data_c->Draw("E SAME");

        drawLabels(0.20, 0.90);
        myText(0.20, 0.75, 1, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), 0.04);
        myMarkerLineText(0.60, 0.65, 1.5, kBlack, 20, kBlack, 1, "Data",        0.04, true);
        myMarkerLineText(0.60, 0.60, 0,   kRed,    0, kRed,   1, "Inclusive MC",0.04, true);

        c->SaveAs(Form("%s/isoET_datamc_eta%d_pt%d.pdf", savePath.c_str(), ieta, ipt));
        delete h_data_c; delete h_mc_c; delete c;
    }

    //==========================================================================
    // Part 4 & 4b: NPB Score vs Time Distribution in Sliced Score Ranges
    // (merged: two score-range definitions, same loop body)
    //==========================================================================
    std::cout << "Creating NPB score vs time overlay plots..." << std::endl;

    struct NPBTimeConfig {
        std::vector<std::pair<double,double>> ranges;
        std::vector<Color_t>                  colors;
        float                                 legend_y1;
        std::string                           output_label; // "overlay" or "threshold"
    };
    std::vector<NPBTimeConfig> npb_configs = {
        { {{0.0,0.2},{0.2,0.4},{0.4,0.6},{0.6,0.8},{0.8,1.0}},
          {kBlue, kGreen+2, kOrange+7, kMagenta+1, kRed},
          0.60f, "overlay" },
        { {{0.0,0.2},{0.2,1.0},{0.5,1.0},{0.8,1.0}},
          {kRed, kBlue, kGreen+2, kMagenta+1},
          0.65f, "threshold" }
    };

    std::vector<std::pair<std::string, std::string>> npb_time_hist_variants = {
        {"h_npb_score_vs_time", ""}, {"h_npb_score_vs_time_clean", "_clean"}, {"h_npb_score_vs_time_strict", "_strict"}
    };

    for (const auto &cfg : npb_configs)
    for (const auto &[hist_base, hist_suffix] : npb_time_hist_variants)
    for (const auto &[file_tag, file_fetcher, file_color, file_label] : all_samples)
    for (int ieta = 0; ieta < nEtaBins; ++ieta)
    for (int ipt  = 0; ipt  < nPtBins;  ++ipt)
    {
        TString hist_name = Form("%s_eta%d_pt%d", hist_base.c_str(), ieta, ipt);
        TH2D *h2 = dynamic_cast<TH2D*>(file_fetcher(hist_name));
        if (!h2) continue;

        std::vector<TH1D*> time_projs;
        float maxY = 0.0;

        TLegend *leg = makeLegend(0.65f, cfg.legend_y1, 0.88f, 0.88f, 0.032f);

        int nRanges = (int)cfg.ranges.size();
        for (int ir = 0; ir < nRanges; ++ir) {
            int ybin_min = h2->GetYaxis()->FindBin(cfg.ranges[ir].first  + 1e-6);
            int ybin_max = h2->GetYaxis()->FindBin(cfg.ranges[ir].second - 1e-6);
            TH1D *proj = h2->ProjectionX(Form("%s_%s_%d", hist_name.Data(), cfg.output_label.c_str(), ir),
                                         ybin_min, ybin_max);
            if (!proj || proj->Integral() < 1e-12) { if (proj) delete proj; continue; }
            scaleToUnit(proj);
            styleMC(proj, cfg.colors[ir]);
            time_projs.push_back(proj);
            maxY = std::max(maxY, (float)proj->GetMaximum());
            leg->AddEntry(proj, Form("%.1f < NPB < %.1f", cfg.ranges[ir].first, cfg.ranges[ir].second), "l");
        }
        if (file_tag != "data") delete h2;

        if (time_projs.empty()) { delete leg; continue; }

        TCanvas *c = new TCanvas(Form("c_npb_%s%s_%s_eta%d_pt%d",
                                      cfg.output_label.c_str(), hist_suffix.c_str(), file_tag.c_str(), ieta, ipt),
                                 "", 600, 600);
        float xlo = (hist_suffix == "_strict") ? -10.f : -20.f;
        time_projs[0]->SetTitle("");
        time_projs[0]->GetXaxis()->SetTitle("Cluster-MBD Time [ns]");
        time_projs[0]->GetXaxis()->SetRangeUser(xlo, 20);
        time_projs[0]->GetXaxis()->SetNdivisions(505);
        time_projs[0]->GetYaxis()->SetTitle("Normalized counts");
        time_projs[0]->GetYaxis()->SetTitleOffset(1.5);
        time_projs[0]->GetYaxis()->SetRangeUser(0, maxY * 1.3);
        time_projs[0]->Draw("HIST"); overlayVerticalErrors(time_projs[0]);
        for (size_t i = 1; i < time_projs.size(); ++i)
        { time_projs[i]->Draw("HIST SAME"); overlayVerticalErrors(time_projs[i]); }

        drawLabels(0.20, 0.90);
        myText(0.20, 0.75, 1, file_label.c_str(), 0.04);
        myText(0.20, 0.70, 1, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), 0.04);
        if (!hist_suffix.empty())
            myText(0.20, 0.65, 1, Form("Selection: %s", hist_suffix.substr(1).c_str()), 0.04);
        leg->Draw();

        c->SaveAs(Form("%s/npb_time_%s%s_%s_eta%d_pt%d.pdf",
                       savePath.c_str(), cfg.output_label.c_str(), hist_suffix.c_str(), file_tag.c_str(), ieta, ipt));
        delete c; delete leg;
    }
    std::cout << "Parts 4 & 4b complete: NPB score vs time plots created" << std::endl;

    //==========================================================================
    // Part 5: Inclusive MC cut0 vs cut1 (no-cut vs NPB-cut) Comparison
    //==========================================================================
    std::cout << "Creating cut0 vs cut1 comparison plots for inclusive MC..." << std::endl;

    for (const auto &hbase : histNames)
    for (int ieta = 0; ieta < nEtaBins; ++ieta)
    for (int ipt  = 0; ipt  < nPtBins;  ++ipt)
    {
        TString xaxisname = hbase.substr(4, hbase.size() - 4);
        float xaxismin, xaxismax; int nrebin;
        getAxisSettings(xaxisname, xaxismin, xaxismax, nrebin);

        TString hname_cut0 = Form("%s_eta%d_pt%d_cut0", hbase.c_str(), ieta, ipt);
        TString hname_cut1 = Form("%s_eta%d_pt%d_cut1", hbase.c_str(), ieta, ipt);

        TH2F *h2_cut0 = dynamic_cast<TH2F*>(getCombined(f_bkg_inc, f_bkginc_double, bkg_double_weight, hname_cut0));
        TH2F *h2_cut1 = dynamic_cast<TH2F*>(getCombined(f_bkg_inc, f_bkginc_double, bkg_double_weight, hname_cut1));

        if (!h2_cut0 || !h2_cut1) {
            std::cerr << "Warning: Could not retrieve histograms for "
                      << hbase << " eta" << ieta << " pt" << ipt << std::endl;
            delete h2_cut0; delete h2_cut1; continue;
        }

        h2_cut0->RebinX(nrebin); h2_cut1->RebinX(nrebin);
        h2_cut0->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
        h2_cut1->GetXaxis()->SetRangeUser(xaxismin, xaxismax);

        TH1D *proj_cut0 = h2_cut0->ProjectionX(Form("%s_px_cut0", hname_cut0.Data()));
        TH1D *proj_cut1 = h2_cut1->ProjectionX(Form("%s_px_cut1", hname_cut1.Data()));

        scaleToUnit(proj_cut0); scaleToUnit(proj_cut1);
        styleMC(proj_cut0, kBlue);
        styleMC(proj_cut1, kRed);

        float maxY = std::max(proj_cut0->GetMaximum(), proj_cut1->GetMaximum());

        TCanvas *c = new TCanvas(Form("c_cutcomp_%s_eta%d_pt%d", hbase.c_str(), ieta, ipt), "", 600, 600);
        proj_cut0->SetTitle("");
        proj_cut0->GetXaxis()->SetTitle(xaxisname.Data());
        proj_cut0->GetXaxis()->SetNdivisions(505);
        proj_cut0->GetYaxis()->SetTitle("Normalized counts");
        proj_cut0->GetYaxis()->SetTitleOffset(1.5);
        proj_cut0->GetYaxis()->SetRangeUser(0, maxY * 1.3);
        proj_cut0->Draw("HIST"); overlayVerticalErrors(proj_cut0);
        proj_cut1->Draw("HIST SAME"); overlayVerticalErrors(proj_cut1);

        drawLabels(0.20, 0.90);
        myText(0.20, 0.80, 1, "Inclusive MC", 0.04);
        myText(0.20, 0.75, 1, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), 0.04);

        TLegend *leg = makeLegend(0.55, 0.75, 0.88, 0.88);
        leg->AddEntry(proj_cut0, "no cut",  "l");
        leg->AddEntry(proj_cut1, "npb cut", "l");
        leg->Draw();

        c->SaveAs(Form("%s/cut_comparison_%s_eta%d_pt%d.pdf", savePath.c_str(), xaxisname.Data(), ieta, ipt));
        delete proj_cut0; delete proj_cut1;
        delete h2_cut0; delete h2_cut1; delete c; delete leg;
    }
    std::cout << "Part 5 complete: Cut comparison plots created" << std::endl;

    //==========================================================================
    // Part 6: NPB Selection Efficiency vs Threshold
    //==========================================================================
    std::cout << "Creating NPB selection efficiency vs threshold plots..." << std::endl;

    std::vector<std::tuple<std::string, Fetcher, Color_t, std::string>> mc_samples = {
        {"inclusive", fetch_bkg_inc, kBlue, "Inclusive MC"},
        {"signal",    fetch_sig,     kRed,  "Signal MC"}
    };

    for (const auto &[mc_tag, mc_fetcher, mc_color, mc_label] : mc_samples)
    {
        // Individual per-(eta, pt) efficiency plots
        for (int ieta = 0; ieta < nEtaBins; ++ieta)
        for (int ipt  = 0; ipt  < nPtBins;  ++ipt)
        {
            TString hname = Form("h2d_npb_score_eta%d_pt%d_cut0", ieta, ipt);
            TH2F  *h2_npb = dynamic_cast<TH2F*>(mc_fetcher(hname));
            if (!h2_npb) continue;

            TH1D *h_npb = h2_npb->ProjectionX(Form("%s_px_%s", hname.Data(), mc_tag.c_str()));
            if (mc_tag != "data") delete h2_npb;
            if (!h_npb || h_npb->Integral() < 1e-12) { if (h_npb) delete h_npb; continue; }

            TH1D *h_eff = computeEffHist(h_npb, Form("h_npb_eff_%s_eta%d_pt%d", mc_tag.c_str(), ieta, ipt));
            delete h_npb;

            TCanvas *c = new TCanvas(Form("c_npb_eff_%s_eta%d_pt%d", mc_tag.c_str(), ieta, ipt),
                                     Form("NPB Eff (%s)", mc_label.c_str()), 600, 600);
            h_eff->SetLineColor(mc_color); h_eff->SetMarkerColor(mc_color);
            h_eff->SetMarkerStyle(20); h_eff->SetMarkerSize(0.8);
            h_eff->SetLineWidth(1); h_eff->SetStats(0); h_eff->SetTitle("");
            h_eff->GetXaxis()->SetTitle("NPB Score Threshold"); h_eff->GetXaxis()->SetNdivisions(505);
            h_eff->GetYaxis()->SetTitle("Selection Efficiency"); h_eff->GetYaxis()->SetTitleOffset(1.5);
            h_eff->GetYaxis()->SetRangeUser(0.8, 1.1);
            h_eff->Draw("E P");

            drawLabels(0.20, 0.90);
            myText(0.20, 0.80, 1, mc_label.c_str(), 0.04);
            myText(0.20, 0.75, 1, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), 0.04);

            c->SaveAs(Form("%s/npb_efficiency_vs_threshold_%s_eta%d_pt%d.pdf", savePath.c_str(), mc_tag.c_str(), ieta, ipt));
            delete h_eff; delete c;
        }

        // Overlay: all pT bins, ieta=0
        {
            TCanvas *c_ov = new TCanvas(Form("c_npb_eff_overlay_%s", mc_tag.c_str()),
                                        Form("NPB Efficiency Overlay (%s)", mc_label.c_str()), 700, 600);
            TLegend *leg  = makeLegend(0.55, 0.15, 0.88, 0.45, 0.035f);
            std::vector<TH1D*> eff_hists;
            bool first = true;

            for (int ipt = 0; ipt < nPtBins; ++ipt) {
                TString hname = Form("h2d_npb_score_eta0_pt%d_cut0", ipt);
                TH2F  *h2_npb = dynamic_cast<TH2F*>(mc_fetcher(hname));
                if (!h2_npb) continue;

                TH1D *h_npb = h2_npb->ProjectionX(Form("%s_px_%s_ov", hname.Data(), mc_tag.c_str()));
                if (mc_tag != "data") delete h2_npb;
                if (!h_npb || h_npb->Integral() < 1e-12) { if (h_npb) delete h_npb; continue; }

                TH1D *h_eff = computeEffHist(h_npb, Form("h_npb_eff_%s_overlay_pt%d", mc_tag.c_str(), ipt));
                delete h_npb;

                h_eff->SetLineColor(pt_colors[ipt % pt_colors.size()]);
                h_eff->SetMarkerColor(pt_colors[ipt % pt_colors.size()]);
                h_eff->SetMarkerStyle(20 + (ipt % 10)); h_eff->SetMarkerSize(0.8);
                h_eff->SetLineWidth(1); h_eff->SetStats(0);

                if (first) {
                    h_eff->SetTitle(""); h_eff->GetXaxis()->SetTitle("NPB Score Threshold");
                    h_eff->GetXaxis()->SetNdivisions(505); h_eff->GetYaxis()->SetTitle("Selection Efficiency");
                    h_eff->GetYaxis()->SetTitleOffset(1.5); h_eff->GetYaxis()->SetRangeUser(0.5, 1.2);
                    h_eff->Draw("E P"); first = false;
                } else {
                    h_eff->Draw("E P SAME");
                }
                eff_hists.push_back(h_eff);
                leg->AddEntry(h_eff, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), "lep");
            }

            if (!eff_hists.empty()) {
                drawLabels(0.20, 0.90);
                myText(0.20, 0.80, 1, mc_label.c_str(), 0.04);
                leg->Draw();
                c_ov->SaveAs(Form("%s/npb_efficiency_vs_threshold_%s_overlay.pdf", savePath.c_str(), mc_tag.c_str()));
            }
            for (auto h : eff_hists) delete h;
            delete leg; delete c_ov;
        }
    }
    std::cout << "Part 6 complete: NPB efficiency plots created" << std::endl;

    //==========================================================================
    // Part 7: NPB Score Ratio N(0.5<NPB<0.9) / N(0.9<NPB<1.0) vs MBD Sigma
    //==========================================================================
    std::cout << "Creating NPB score ratio vs MBD sigma plots..." << std::endl;

    for (const auto &[file_tag, file_fetcher, file_color, file_label] : all_samples)
    {
        // Individual per-(eta, pt) ratio plots
        for (int ieta = 0; ieta < nEtaBins; ++ieta)
        for (int ipt  = 0; ipt  < nPtBins;  ++ipt)
        {
            TString hname = Form("h_mbd_avgsigma_vs_npb_score_eta%d_pt%d", ieta, ipt);
            TH2D *h2 = dynamic_cast<TH2D*>(file_fetcher(hname));
            if (!h2) { std::cerr << "Warning: Could not retrieve " << hname << " from " << file_tag << std::endl; continue; }

            TH1D *h_ratio = computeNPBRatioHist(h2, Form("h_npb_ratio_vs_mbd_%s_eta%d_pt%d", file_tag.c_str(), ieta, ipt));
            if (file_tag != "data") delete h2;

            TCanvas *c = new TCanvas(Form("c_npb_ratio_mbd_%s_eta%d_pt%d", file_tag.c_str(), ieta, ipt),
                                     "NPB Ratio vs MBD Sigma", 600, 600);
            h_ratio->SetLineColor(file_color); h_ratio->SetMarkerColor(file_color);
            h_ratio->SetMarkerStyle(20); h_ratio->SetMarkerSize(0.8);
            h_ratio->SetLineWidth(1); h_ratio->SetStats(0); h_ratio->SetTitle("");
            h_ratio->GetXaxis()->SetTitle("MBD Avg #sigma_{t} [ns]"); h_ratio->GetXaxis()->SetNdivisions(505);
            h_ratio->GetYaxis()->SetTitle("N(0.5<NPB<0.9) / N(0.9<NPB<1.0)"); h_ratio->GetYaxis()->SetTitleOffset(1.5);
            h_ratio->GetXaxis()->SetRangeUser(0, 5);
            double max_ratio = h_ratio->GetMaximum();
            h_ratio->GetYaxis()->SetRangeUser(0, max_ratio > 0 ? max_ratio * 1.3 : 1.0);
            h_ratio->Draw("E P");

            drawLabels(0.20, 0.90);
            myText(0.20, 0.80, 1, file_label.c_str(), 0.04);
            myText(0.20, 0.75, 1, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), 0.04);

            c->SaveAs(Form("%s/npb_ratio_vs_mbd_sigma_%s_eta%d_pt%d.pdf", savePath.c_str(), file_tag.c_str(), ieta, ipt));
            delete h_ratio; delete c;
        }

        // Overlay: all pT bins at ieta=0
        {
            TCanvas *c_ov = new TCanvas(Form("c_npb_ratio_mbd_overlay_%s", file_tag.c_str()),
                                        Form("NPB Ratio vs MBD Sigma Overlay (%s)", file_label.c_str()), 700, 600);
            TLegend *leg  = makeLegend(0.55, 0.60, 0.88, 0.88, 0.035f);
            std::vector<TH1D*> ratio_hists;
            double global_max = 0;
            bool first = true;

            for (int ipt = 0; ipt < nPtBins; ++ipt) {
                TString hname = Form("h_mbd_avgsigma_vs_npb_score_eta0_pt%d", ipt);
                TH2D *h2 = dynamic_cast<TH2D*>(file_fetcher(hname));
                if (!h2) continue;

                TH1D *h_ratio = computeNPBRatioHist(h2, Form("h_npb_ratio_overlay_%s_pt%d", file_tag.c_str(), ipt));
                if (file_tag != "data") delete h2;

                h_ratio->SetLineColor(pt_colors[ipt % pt_colors.size()]);
                h_ratio->SetMarkerColor(pt_colors[ipt % pt_colors.size()]);
                h_ratio->SetMarkerStyle(20 + (ipt % 10)); h_ratio->SetMarkerSize(0.8);
                h_ratio->SetLineWidth(1); h_ratio->SetStats(0);
                if (h_ratio->GetMaximum() > global_max) global_max = h_ratio->GetMaximum();

                if (first) {
                    h_ratio->SetTitle(""); h_ratio->GetXaxis()->SetTitle("MBD Avg #sigma_{t} [ns]");
                    h_ratio->GetXaxis()->SetNdivisions(505); h_ratio->GetXaxis()->SetRangeUser(0, 5);
                    h_ratio->GetYaxis()->SetTitle("N(0.5<NPB<0.9) / N(0.9<NPB<1.0)");
                    h_ratio->GetYaxis()->SetTitleOffset(1.5);
                    h_ratio->Draw("E P"); first = false;
                } else {
                    h_ratio->Draw("E P SAME");
                }
                ratio_hists.push_back(h_ratio);
                leg->AddEntry(h_ratio, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt+1]), "lep");
            }

            if (!ratio_hists.empty()) {
                // Fix BUG-1: update Y range and refresh canvas after all histograms are drawn
                ratio_hists[0]->GetYaxis()->SetRangeUser(0, global_max * 1.3);
                c_ov->Modified(); c_ov->Update();
                drawLabels(0.20, 0.90);
                myText(0.20, 0.80, 1, file_label.c_str(), 0.04);
                leg->Draw();
                c_ov->SaveAs(Form("%s/npb_ratio_vs_mbd_sigma_%s_overlay.pdf", savePath.c_str(), file_tag.c_str()));
            }
            for (auto h : ratio_hists) delete h;
            delete leg; delete c_ov;
        }
    }
    std::cout << "Part 7 complete: NPB ratio vs MBD sigma plots created" << std::endl;

    //==========================================================================
    // Close files
    //==========================================================================
    f_data->Close(); f_sig->Close(); f_bkg_inc->Close();
    if (hasBkgOnly && f_bkg_only && f_bkg_only != f_bkg_inc) f_bkg_only->Close();
    if (f_sig_double)    f_sig_double->Close();
    if (f_bkginc_double) f_bkginc_double->Close();

    std::cout << "All plots created successfully!" << std::endl;
    std::cout << "Output saved to: " << savePath << std::endl;
}
