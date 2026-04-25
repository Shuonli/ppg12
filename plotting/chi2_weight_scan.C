// chi2_weight_scan.C
// Scans the double-interaction weight w in combined = (1-w)*nominal + w*double
// and computes chi2 between data and combined MC for each shower-shape variable.
// Produces chi2 vs weight plots (total, per-ndf, per-shape) and a summary ROOT file.
//
// Usage:
//   root -l -b -q 'chi2_weight_scan.C'
//   root -l -b -q 'chi2_weight_scan.C("config_showershape_0rad.yaml", 21, 0.0, 1.0)'
//   root -l -b -q 'chi2_weight_scan.C("config_showershape_0rad.yaml", 21, 0.0, 1.0, true)'
//     (useBdtForBest=true → best weight from BDT score chi2 only)
//   root -l -b -q 'chi2_weight_scan.C("config_showershape_0rad.yaml", 21, 0.0, 1.0, true, 1, 3)'
//     (ptBinMin=1, ptBinMax=3 → chi2 uses only pT bins [1,3]; -1 means use all)

#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TLine.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <yaml-cpp/yaml.h>
#include <TSystem.h>

#include "plotcommon.h"

void chi2_weight_scan(
    const std::string &configname = "config_showershape_1p5mrad.yaml",
    int    nSteps = 201,
    double wMin   = 0.0,
    double wMax   = 1.0,
    bool   useBdtForBest = true,
    int    ptBinMin = 1,   // first pT bin index for chi2 (-1 = use all)
    int    ptBinMax = 3)   // last  pT bin index for chi2 (-1 = use all)
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    init_plot();
    gSystem->Exec("mkdir -p figures");

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

    // Load config for pT bins
    YAML::Node config = YAML::LoadFile(("../efficiencytool/" + configname).c_str());
    std::vector<double> pT_bin_edges = config["analysis"]["pT_bins"].as<std::vector<double>>();
    std::vector<double> eta_bins = {-0.7, 0.7};
    double vertex_cut = config["analysis"]["vertex_cut"].as<double>(60.0);

    int nPtBins  = pT_bin_edges.size() - 1;
    int nEtaBins = eta_bins.size() - 1;

    // Resolve pT bin range for chi2 (clamp to valid indices)
    int chi2PtMin = (ptBinMin < 0) ? 0 : std::max(0, ptBinMin);
    int chi2PtMax = (ptBinMax < 0) ? nPtBins - 1 : std::min(nPtBins - 1, ptBinMax);
    if (chi2PtMin > chi2PtMax) {
        std::cerr << "Error: ptBinMin (" << ptBinMin << ") > ptBinMax (" << ptBinMax
                  << ") after clamping. Using all pT bins." << std::endl;
        chi2PtMin = 0; chi2PtMax = nPtBins - 1;
    }
    std::cout << "Using " << nPtBins << " pT bins from config file" << std::endl;
    std::cout << "Chi2 computed over pT bins [" << chi2PtMin << ", " << chi2PtMax << "] = ["
              << pT_bin_edges[chi2PtMin] << ", " << pT_bin_edges[chi2PtMax + 1] << "] GeV" << std::endl;

    // Open files
    const std::string resultsDir       = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/";
    const std::string dataFile         = resultsDir + "data_histoshower_shape_" + config_suffix + ".root";
    const std::string bkgInclusiveFile = resultsDir + "MC_efficiencyshower_shape_jet12_nom_inclusive_" + config_suffix + ".root";
    const std::string bkgDoubleFile    = resultsDir + "MC_efficiencyshower_shape_jet12_double_inclusive_" + config_suffix + ".root";

    TFile *f_data    = TFile::Open(dataFile.c_str(), "READ");
    TFile *f_bkg_nom = TFile::Open(bkgInclusiveFile.c_str(), "READ");

    if (!f_data || f_data->IsZombie())
    {
        std::cerr << "Error: Could not open data file: " << dataFile << std::endl;
        return;
    }
    if (!f_bkg_nom || f_bkg_nom->IsZombie())
    {
        std::cerr << "Error: Could not open nominal bkg file: " << bkgInclusiveFile << std::endl;
        return;
    }

    TFile *f_bkg_dbl = TFile::Open(bkgDoubleFile.c_str(), "READ");
    bool hasDouble = (f_bkg_dbl && !f_bkg_dbl->IsZombie());
    if (!hasDouble)
    {
        std::cerr << "Warning: Could not open double-interaction bkg file: " << bkgDoubleFile << std::endl;
        std::cerr << "         Will scan with nominal only (all weights give same result)." << std::endl;
    }

    const std::string sigFile = resultsDir + "MC_efficiencyshower_shape_photon10_nom_" + config_suffix + ".root";
    TFile *f_sig = TFile::Open(sigFile.c_str(), "READ");
    bool hasSig = (f_sig && !f_sig->IsZombie());
    if (!hasSig)
        std::cerr << "Warning: Could not open signal file: " << sigFile << " — signal omitted from overlays." << std::endl;

    // Pre-scale double-interaction template to nominal using h_vertexz integral
    // in the configured vertex window, then apply scan weight on top.
    double doubleVertexNorm = 1.0;
    if (hasDouble)
    {
        TH1 *h_vtx_nom = dynamic_cast<TH1 *>(f_bkg_nom->Get("h_vertexz"));
        TH1 *h_vtx_dbl = dynamic_cast<TH1 *>(f_bkg_dbl->Get("h_vertexz"));
        if (!h_vtx_nom || !h_vtx_dbl)
        {
            std::cerr << "Warning: h_vertexz missing in nominal or double file. "
                      << "Skipping vertex integral pre-normalization (factor=1)." << std::endl;
        }
        else
        {
            const int binLoNom = h_vtx_nom->GetXaxis()->FindBin(-vertex_cut + 1e-6);
            const int binHiNom = h_vtx_nom->GetXaxis()->FindBin( vertex_cut - 1e-6);
            const int binLoDbl = h_vtx_dbl->GetXaxis()->FindBin(-vertex_cut + 1e-6);
            const int binHiDbl = h_vtx_dbl->GetXaxis()->FindBin( vertex_cut - 1e-6);

            const double intNom = h_vtx_nom->Integral(binLoNom, binHiNom);
            const double intDbl = h_vtx_dbl->Integral(binLoDbl, binHiDbl);
            if (intDbl > 1e-12)
            {
                doubleVertexNorm = intNom / intDbl;
                std::cout << "[VertexNorm] Using |z| < " << vertex_cut
                          << " cm: nominal integral = " << intNom
                          << ", double integral = " << intDbl
                          << ", double scale = " << doubleVertexNorm << std::endl;
            }
            else
            {
                std::cerr << "Warning: double h_vertexz integral in |z| < " << vertex_cut
                          << " is zero/non-positive. Using factor=1." << std::endl;
            }
        }
    }

    // Histogram base names
    std::vector<std::string> histNames = {
        "h2d_weta_cogx", "h2d_wphi_cogx", "h2d_et1", "h2d_et2", "h2d_et3", "h2d_et4",
        "h2d_e11_to_e33", "h2d_e17_to_e77", "h2d_e32_to_e35", "h2d_bdt", "h2d_npb_score"
    };
    int nShapes = histNames.size();

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
    auto scaleToUnit = [](TH1D *h) {
        if (!h) return;
        double integral = h->Integral();
        if (integral > 1e-12)
            h->Scale(1.0 / integral);
    };

    // =========================================================================
    // Pre-compute projections (outside weight loop)
    // =========================================================================
    // Indexed by [ishape][ieta][ipt]
    struct ProjSet {
        TH1D *data = nullptr;
        TH1D *nom  = nullptr;
        TH1D *dbl  = nullptr;
        TH1D *sig  = nullptr;
        float xmin = 0;
        float xmax = 1;
    };

    std::vector<std::vector<std::vector<ProjSet>>> projs(
        nShapes, std::vector<std::vector<ProjSet>>(
            nEtaBins, std::vector<ProjSet>(nPtBins)));

    for (int ishape = 0; ishape < nShapes; ++ishape)
    {
        const std::string &hbase = histNames[ishape];
        TString xaxisname = hbase.substr(4, hbase.size() - 4);
        float xmin, xmax;
        int nrebin;
        getAxisSettings(xaxisname, xmin, xmax, nrebin);

        for (int ieta = 0; ieta < nEtaBins; ++ieta)
        {
            for (int ipt = 0; ipt < nPtBins; ++ipt)
            {
                TString histNameFull = Form("%s_eta%d_pt%d_cut1",
                                            hbase.c_str(), ieta, ipt);

                // Data
                TH2F *h2_data = dynamic_cast<TH2F *>(f_data->Get(histNameFull));
                // Nominal bkg
                TH2F *h2_nom = dynamic_cast<TH2F *>(f_bkg_nom->Get(histNameFull));
                // Double bkg
                TH2F *h2_dbl = hasDouble ? dynamic_cast<TH2F *>(f_bkg_dbl->Get(histNameFull)) : nullptr;

                if (!h2_data || !h2_nom)
                {
                    std::cerr << "Warning: Missing histogram " << histNameFull
                              << " in data or nominal bkg file, skipping." << std::endl;
                    continue;
                }

                // Clone to avoid modifying originals, then rebin and set range
                TH2F *h2d_clone = (TH2F *)h2_data->Clone(Form("%s_clone_data", histNameFull.Data()));
                h2d_clone->SetDirectory(nullptr);
                h2d_clone->Sumw2();
                h2d_clone->RebinX(nrebin);
                h2d_clone->GetXaxis()->SetRangeUser(xmin, xmax);
                TH1D *proj_data = h2d_clone->ProjectionX(Form("%s_px_data_chi2", histNameFull.Data()));
                proj_data->SetDirectory(nullptr);
                proj_data->Sumw2();
                scaleToUnit(proj_data);
                delete h2d_clone;

                TH2F *h2n_clone = (TH2F *)h2_nom->Clone(Form("%s_clone_nom", histNameFull.Data()));
                h2n_clone->SetDirectory(nullptr);
                h2n_clone->Sumw2();
                h2n_clone->RebinX(nrebin);
                h2n_clone->GetXaxis()->SetRangeUser(xmin, xmax);
                TH1D *proj_nom = h2n_clone->ProjectionX(Form("%s_px_nom_chi2", histNameFull.Data()));
                proj_nom->SetDirectory(nullptr);
                proj_nom->Sumw2();
                delete h2n_clone;

                TH1D *proj_dbl = nullptr;
                if (h2_dbl)
                {
                    TH2F *h2db_clone = (TH2F *)h2_dbl->Clone(Form("%s_clone_dbl", histNameFull.Data()));
                    h2db_clone->SetDirectory(nullptr);
                    h2db_clone->Sumw2();
                    h2db_clone->RebinX(nrebin);
                    h2db_clone->GetXaxis()->SetRangeUser(xmin, xmax);
                    proj_dbl = h2db_clone->ProjectionX(Form("%s_px_dbl_chi2", histNameFull.Data()));
                    proj_dbl->SetDirectory(nullptr);
                    proj_dbl->Sumw2();
                    delete h2db_clone;
                }

                // Signal
                TH1D *proj_sig = nullptr;
                if (hasSig) {
                    TH2F *h2_sig = dynamic_cast<TH2F *>(f_sig->Get(histNameFull));
                    if (h2_sig) {
                        TH2F *h2s_clone = (TH2F *)h2_sig->Clone(Form("%s_clone_sig", histNameFull.Data()));
                        h2s_clone->SetDirectory(nullptr);
                        h2s_clone->Sumw2();
                        h2s_clone->RebinX(nrebin);
                        h2s_clone->GetXaxis()->SetRangeUser(xmin, xmax);
                        proj_sig = h2s_clone->ProjectionX(Form("%s_px_sig_chi2", histNameFull.Data()));
                        proj_sig->SetDirectory(nullptr);
                        proj_sig->Sumw2();
                        delete h2s_clone;
                    }
                }

                projs[ishape][ieta][ipt].data = proj_data;
                projs[ishape][ieta][ipt].nom  = proj_nom;
                projs[ishape][ieta][ipt].dbl  = proj_dbl;
                projs[ishape][ieta][ipt].sig  = proj_sig;
                projs[ishape][ieta][ipt].xmin = xmin;
                projs[ishape][ieta][ipt].xmax = xmax;
            }
        }
    }

    // =========================================================================
    // Weight scan
    // =========================================================================
    std::vector<double> weights(nSteps);
    for (int iw = 0; iw < nSteps; ++iw)
    {
        weights[iw] = (nSteps > 1) ? wMin + iw * (wMax - wMin) / (nSteps - 1) : wMin;
    }

    // chi2 accumulators: [iw] for total, [iw][ishape] for per-shape
    std::vector<double> chi2Total(nSteps, 0.0);
    std::vector<int>    ndofTotal(nSteps, 0);
    std::vector<std::vector<double>> chi2PerShape(nSteps, std::vector<double>(nShapes, 0.0));
    std::vector<std::vector<int>>    ndofPerShape(nSteps, std::vector<int>(nShapes, 0));

    for (int iw = 0; iw < nSteps; ++iw)
    {
        double w = weights[iw];

        for (int ishape = 0; ishape < nShapes; ++ishape)
        {
            for (int ieta = 0; ieta < nEtaBins; ++ieta)
            {
                for (int ipt = 0; ipt < nPtBins; ++ipt)
                {
                    if (ipt < chi2PtMin || ipt > chi2PtMax) continue;

                    const ProjSet &ps = projs[ishape][ieta][ipt];
                    if (!ps.data || !ps.nom) continue;

                    // Clone nominal, combine with double
                    TH1D *proj_comb = (TH1D *)ps.nom->Clone(
                        Form("comb_iw%d_s%d_e%d_p%d", iw, ishape, ieta, ipt));
                    proj_comb->SetDirectory(nullptr);

                    if (ps.dbl && w > 1e-9)
                    {
                        proj_comb->Scale(1.0 - w);
                        proj_comb->Add(ps.dbl, w * doubleVertexNorm);
                    }

                    scaleToUnit(proj_comb);

                    // Compute chi2 over bins in [xmin, xmax]
                    int nbins = ps.data->GetNbinsX();
                    for (int ib = 1; ib <= nbins; ++ib)
                    {
                        double bc = ps.data->GetBinCenter(ib);
                        if (bc < ps.xmin || bc > ps.xmax) continue;

                        double d  = ps.data->GetBinContent(ib);
                        double m  = proj_comb->GetBinContent(ib);
                        double ed = ps.data->GetBinError(ib);
                        double em = proj_comb->GetBinError(ib);
                        double denom = ed * ed + em * em;
                        if (denom < 1e-30) continue;

                        double chi2_bin = (d - m) * (d - m) / denom;
                        // Default: shower shapes only (exclude bdt and npb_score).
                        // useBdtForBest: use bdt score chi2 only to determine best weight.
                        bool inTotal = useBdtForBest
                            ? (histNames[ishape] == "h2d_bdt")
                            : (histNames[ishape] != "h2d_bdt" &&
                               histNames[ishape] != "h2d_npb_score");
                        if (inTotal) {
                            chi2Total[iw] += chi2_bin;
                            ndofTotal[iw] += 1;
                        }
                        chi2PerShape[iw][ishape] += chi2_bin;
                        ndofPerShape[iw][ishape] += 1;
                    }

                    delete proj_comb;
                }
            }
        }
    }

    // =========================================================================
    // Print results table
    // =========================================================================
    // Short names for table header
    std::vector<std::string> shortNames;
    for (const auto &hbase : histNames)
        shortNames.push_back(hbase.substr(4, hbase.size() - 4));

    std::cout << "\n============================================================" << std::endl;
    std::cout << "Chi2 Weight Scan Results (" << config_suffix << ")" << std::endl;
    std::cout << "============================================================" << std::endl;
    printf("%8s %12s %6s %10s", "weight", "chi2_tot", "ndof", "chi2/ndf");
    for (int is = 0; is < nShapes; ++is)
        printf(" %14s", shortNames[is].c_str());
    printf("\n");

    for (int iw = 0; iw < nSteps; ++iw)
    {
        double ndf = (ndofTotal[iw] > 0) ? (double)ndofTotal[iw] : 1.0;
        printf("%8.4f %12.2f %6d %10.4f", weights[iw], chi2Total[iw], ndofTotal[iw], chi2Total[iw] / ndf);
        for (int is = 0; is < nShapes; ++is)
            printf(" %14.2f", chi2PerShape[iw][is]);
        printf("\n");
    }

    // Find best weight (minimum chi2Total)
    int bestIdx = 0;
    for (int iw = 1; iw < nSteps; ++iw)
    {
        if (chi2Total[iw] < chi2Total[bestIdx])
            bestIdx = iw;
    }
    std::cout << "\nBest weight (" << (useBdtForBest ? "BDT score chi2" : "shower shapes chi2") << "): " << weights[bestIdx]
              << " (chi2 = " << chi2Total[bestIdx]
              << ", ndof = " << ndofTotal[bestIdx]
              << ", chi2/ndf = " << chi2Total[bestIdx] / std::max(1, ndofTotal[bestIdx])
              << ")" << std::endl;

    // Optimal weight per shower shape
    std::cout << "\n--- Optimal weight per shower shape ---" << std::endl;
    printf("%-22s  %8s  %12s  %6s  %10s\n", "shape", "best_w", "chi2_min", "ndof", "chi2/ndf");
    for (int is = 0; is < nShapes; ++is)
    {
        int bestShapeIdx = 0;
        for (int iw = 1; iw < nSteps; ++iw)
            if (chi2PerShape[iw][is] < chi2PerShape[bestShapeIdx][is])
                bestShapeIdx = iw;
        int ndf = std::max(1, ndofPerShape[bestShapeIdx][is]);
        printf("%-22s  %8.4f  %12.2f  %6d  %10.4f\n",
               shortNames[is].c_str(), weights[bestShapeIdx],
               chi2PerShape[bestShapeIdx][is], ndofPerShape[bestShapeIdx][is],
               (double)chi2PerShape[bestShapeIdx][is] / ndf);
    }
    std::cout << "---------------------------------------" << std::endl;

    // 1-sigma interval for one scan parameter from Delta(chi2)=1
    const double chi2Min = chi2Total[bestIdx];
    const double chi2At1Sigma = chi2Min + 1.0;
    bool hasLower1Sigma = false;
    bool hasUpper1Sigma = false;
    double w1SigmaLow = weights.front();
    double w1SigmaHigh = weights.back();

    // Left crossing: first crossing of chi2 = chi2Min+1 when moving left from best
    for (int i = bestIdx - 1; i >= 0; --i)
    {
        double x1 = weights[i];
        double x2 = weights[i + 1];
        double y1 = chi2Total[i];
        double y2 = chi2Total[i + 1];
        if ((y1 - chi2At1Sigma) * (y2 - chi2At1Sigma) <= 0.0)
        {
            if (std::abs(y2 - y1) < 1e-12) w1SigmaLow = 0.5 * (x1 + x2);
            else w1SigmaLow = x1 + (chi2At1Sigma - y1) * (x2 - x1) / (y2 - y1);
            hasLower1Sigma = true;
            break;
        }
    }

    // Right crossing: first crossing of chi2 = chi2Min+1 when moving right from best
    for (int i = bestIdx; i < nSteps - 1; ++i)
    {
        double x1 = weights[i];
        double x2 = weights[i + 1];
        double y1 = chi2Total[i];
        double y2 = chi2Total[i + 1];
        if ((y1 - chi2At1Sigma) * (y2 - chi2At1Sigma) <= 0.0)
        {
            if (std::abs(y2 - y1) < 1e-12) w1SigmaHigh = 0.5 * (x1 + x2);
            else w1SigmaHigh = x1 + (chi2At1Sigma - y1) * (x2 - x1) / (y2 - y1);
            hasUpper1Sigma = true;
            break;
        }
    }

    // =========================================================================
    // Output plots
    // =========================================================================
    // Build TGraphs
    TGraph *gr_chi2_total = new TGraph(nSteps);
    gr_chi2_total->SetName("gr_chi2_total");
    gr_chi2_total->SetTitle(";Double-interaction weight #it{w};#chi^{2}_{total}");
    gr_chi2_total->SetLineWidth(2);
    gr_chi2_total->SetMarkerStyle(20);
    gr_chi2_total->SetMarkerSize(0.8);

    TGraph *gr_chi2_ndf = new TGraph(nSteps);
    gr_chi2_ndf->SetName("gr_chi2_ndf");
    gr_chi2_ndf->SetTitle(";Double-interaction weight #it{w};#chi^{2}_{total}");
    gr_chi2_ndf->SetLineWidth(2);
    gr_chi2_ndf->SetMarkerStyle(20);
    gr_chi2_ndf->SetMarkerSize(0.8);

    for (int iw = 0; iw < nSteps; ++iw)
    {
        gr_chi2_total->SetPoint(iw, weights[iw], chi2Total[iw]);
        gr_chi2_ndf->SetPoint(iw, weights[iw], chi2Total[iw]);
    }

    std::vector<TGraph *> gr_per_shape(nShapes);
    std::vector<Color_t> colors = {kBlue, kRed, kGreen + 2, kMagenta + 1, kOrange + 7,
                                   kCyan + 2, kViolet + 1, kTeal + 2, kPink + 1,
                                   kAzure + 2, kYellow + 2};
    for (int is = 0; is < nShapes; ++is)
    {
        gr_per_shape[is] = new TGraph(nSteps);
        gr_per_shape[is]->SetName(Form("gr_chi2_%s", shortNames[is].c_str()));
        gr_per_shape[is]->SetTitle(shortNames[is].c_str());
        gr_per_shape[is]->SetLineWidth(2);
        gr_per_shape[is]->SetLineColor(colors[is % colors.size()]);
        gr_per_shape[is]->SetMarkerColor(colors[is % colors.size()]);
        gr_per_shape[is]->SetMarkerStyle(20 + is % 10);
        gr_per_shape[is]->SetMarkerSize(0.7);
        for (int iw = 0; iw < nSteps; ++iw)
            gr_per_shape[is]->SetPoint(iw, weights[iw], chi2PerShape[iw][is]);
    }

    // --- Plot 1: chi2_total vs weight ---
    {
        TCanvas *c = new TCanvas("c_chi2_total", "chi2 total vs weight", 700, 500);
        c->cd();
        gr_chi2_total->Draw("APL");

        double yMax_total = *std::max_element(chi2Total.begin(), chi2Total.end()) * 1.2;
        TLine *lbest = new TLine(weights[bestIdx], 0.0,
                                 weights[bestIdx], yMax_total);
        lbest->SetLineStyle(2);
        lbest->SetLineColor(kRed);
        lbest->SetLineWidth(2);
        lbest->Draw("SAME");

        myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
        myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
        myText(0.20, 0.80, 1, Form("Best w = %.3f", weights[bestIdx]), 0.035);

        c->SaveAs("figures/chi2_weight_scan_total.pdf");
        delete lbest;
        delete c;
    }

    // --- Plot 2: chi2 vs weight with 1-sigma interval ---
    {
        TCanvas *c = new TCanvas("c_chi2_ndf", "chi2 vs weight", 700, 500);
        c->cd();
        gr_chi2_ndf->Draw("APL");

        double yMax_ndf = *std::max_element(chi2Total.begin(), chi2Total.end()) * 1.2;
        TLine *lbest = new TLine(weights[bestIdx], 0.0,
                                 weights[bestIdx], yMax_ndf);
        lbest->SetLineStyle(2);
        lbest->SetLineColor(kRed);
        lbest->SetLineWidth(2);
        lbest->Draw("SAME");

        TLine *l1s = new TLine(weights.front(), chi2At1Sigma, weights.back(), chi2At1Sigma);
        l1s->SetLineStyle(2);
        l1s->SetLineColor(kGray + 2);
        l1s->SetLineWidth(2);
        l1s->Draw("SAME");

        TLine *llow = new TLine(w1SigmaLow, 0.0, w1SigmaLow, yMax_ndf);
        llow->SetLineStyle(3);
        llow->SetLineColor(kBlue + 1);
        llow->SetLineWidth(2);
        llow->Draw("SAME");

        TLine *lhigh = new TLine(w1SigmaHigh, 0.0, w1SigmaHigh, yMax_ndf);
        lhigh->SetLineStyle(3);
        lhigh->SetLineColor(kBlue + 1);
        lhigh->SetLineWidth(2);
        lhigh->Draw("SAME");

        myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
        myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
        myText(0.20, 0.80, 1, Form("Best w = %.3f", weights[bestIdx]), 0.035);
        myText(0.20, 0.75, 1, Form("1#sigma: [%.3f, %.3f] (#Delta#chi^{2}=1)", w1SigmaLow, w1SigmaHigh), 0.035);
        if (!hasLower1Sigma || !hasUpper1Sigma)
            myText(0.20, 0.70, 1, "Note: 1#sigma bound reaches scan edge", 0.03);

        c->SaveAs("figures/chi2_weight_scan_ndf.pdf");
        delete lhigh;
        delete llow;
        delete l1s;
        delete lbest;
        delete c;
    }

    // --- Plot 3: per-shape chi2 vs weight ---
    {
        TCanvas *c = new TCanvas("c_chi2_per_shape", "chi2 per shape vs weight", 800, 600);
        c->cd();

        TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle(";Double-interaction weight #it{w};#chi^{2} per shape");
        for (int is = 0; is < nShapes; ++is)
        {
            if (shortNames[is] == "npb_score") continue;
            mg->Add(gr_per_shape[is], "PL");
        }
        mg->Draw("A");

        TLegend *leg = new TLegend(0.55, 0.45, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.030);
        for (int is = 0; is < nShapes; ++is)
        {
            if (shortNames[is] == "npb_score") continue;
            leg->AddEntry(gr_per_shape[is], shortNames[is].c_str(), "lp");
        }
        leg->Draw();

        myText(0.20, 0.90, 1, strleg1.c_str(), 0.035);
        myText(0.20, 0.85, 1, strleg2.c_str(), 0.035);

        c->SaveAs("figures/chi2_weight_scan_per_shape.pdf");
        // delete c also deletes mg via pad primitive cleanup (fPrimitives->Delete())
        // Do NOT delete mg explicitly — it's a double-free
        delete c;
    }

    // --- Save all TGraphs to ROOT file ---
    {
        TFile *fout = TFile::Open("figures/chi2_weight_scan.root", "RECREATE");
        gr_chi2_total->Write();
        gr_chi2_ndf->Write();
        for (int is = 0; is < nShapes; ++is)
            gr_per_shape[is]->Write();
        fout->Close();
        delete fout;
    }

    // =========================================================================
    // Overlay plots: data vs combined_bkg_inc vs signal for selected weights
    // Histograms are summed over eta bins for each pT bin, normalized to unit area.
    // Selected weights: wMin, best, wMax (deduplicated).
    // =========================================================================
    {
        // Build selected weight indices (deduplicated, sorted)
        std::vector<int> selIdx;
        for (int iw : {0, bestIdx, nSteps-1}) {
            iw = std::max(0, std::min(nSteps-1, iw));
            if (std::find(selIdx.begin(), selIdx.end(), iw) == selIdx.end())
                selIdx.push_back(iw);
        }
        std::sort(selIdx.begin(), selIdx.end());

        std::vector<Color_t> wColors = {kBlue+1, kGreen+2, kRed+1, kOrange+7, kViolet+1};

        // Mirror reference overlayVerticalErrors: clone MC hist, draw "E0 SAME" for stat errors
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

        TCanvas *coverlay = new TCanvas("c_overlay", "overlay", 900, 700);

        for (int ipt = 0; ipt < nPtBins; ++ipt)
        {
            TString outPdf = Form("figures/chi2_weight_scan_overlays_cut1_%s_pt%d.pdf",
                                  config_suffix.c_str(), ipt);
            coverlay->Print(Form("%s[", outPdf.Data())); // open multi-page per pT bin

            for (int ishape = 0; ishape < nShapes; ++ishape)
            {
                const std::string &hbase = histNames[ishape];
                TString var = hbase.substr(4);

                // Sum data and signal projections over eta bins at fixed pT bin
                TH1D *h_data_sum = nullptr;
                TH1D *h_sig_sum  = nullptr;
                for (int ieta = 0; ieta < nEtaBins; ++ieta) {
                    const ProjSet &ps = projs[ishape][ieta][ipt];
                    if (ps.data) {
                        if (!h_data_sum) {
                            h_data_sum = (TH1D *)ps.data->Clone(Form("sum_data_%s_pt%d", var.Data(), ipt));
                            h_data_sum->SetDirectory(nullptr);
                        } else {
                            h_data_sum->Add(ps.data);
                        }
                    }
                    if (ps.sig) {
                        if (!h_sig_sum) {
                            h_sig_sum = (TH1D *)ps.sig->Clone(Form("sum_sig_%s_pt%d", var.Data(), ipt));
                            h_sig_sum->SetDirectory(nullptr);
                        } else {
                            h_sig_sum->Add(ps.sig);
                        }
                    }
                }
                if (!h_data_sum) continue;

                // Match reference style: normalize each summed distribution to unit area
                scaleToUnit(h_data_sum);
                if (h_sig_sum) scaleToUnit(h_sig_sum);

                // Style data
                h_data_sum->SetLineColor(kBlack); h_data_sum->SetMarkerColor(kBlack);
                h_data_sum->SetMarkerStyle(20);   h_data_sum->SetMarkerSize(0.8);
                h_data_sum->SetStats(0);
                h_data_sum->GetXaxis()->SetTitle(var.Data());
                h_data_sum->GetYaxis()->SetTitle("normalized counts");

                // Build combined bkg sums for each selected weight
                std::vector<TH1D *> h_bkg_sums(selIdx.size(), nullptr);
                for (int si = 0; si < (int)selIdx.size(); ++si) {
                    int iw = selIdx[si];
                    double w = weights[iw];
                    for (int ieta = 0; ieta < nEtaBins; ++ieta) {
                        const ProjSet &ps = projs[ishape][ieta][ipt];
                        if (!ps.nom) continue;
                        TH1D *pcomb = (TH1D *)ps.nom->Clone(
                            Form("ov_comb_s%d_w%d_e%d_p%d", ishape, iw, ieta, ipt));
                        pcomb->SetDirectory(nullptr);
                        if (ps.dbl && w > 1e-9) {
                            pcomb->Scale(1.0 - w);
                            pcomb->Add(ps.dbl, w * doubleVertexNorm);
                        }
                        if (!h_bkg_sums[si]) {
                            h_bkg_sums[si] = pcomb;
                        } else {
                            h_bkg_sums[si]->Add(pcomb);
                            delete pcomb;
                        }
                    }
                    if (h_bkg_sums[si]) {
                        scaleToUnit(h_bkg_sums[si]);
                        h_bkg_sums[si]->SetLineColor(wColors[si % wColors.size()]);
                        h_bkg_sums[si]->SetLineWidth(2);
                        h_bkg_sums[si]->SetMarkerSize(0);
                        h_bkg_sums[si]->SetStats(0);
                        h_bkg_sums[si]->GetXaxis()->SetTitle(var.Data());
                        h_bkg_sums[si]->GetYaxis()->SetTitle("normalized counts");
                    }
                }

                // Draw
                coverlay->cd();
                coverlay->Clear();

                // Set y-scale from Signal MC (fallback to data/MC max if missing)
                double ymax_base = 0.0;
                if (h_sig_sum) {
                    ymax_base = h_sig_sum->GetMaximum();
                } else {
                    ymax_base = h_data_sum->GetMaximum();
                    for (auto *h : h_bkg_sums)
                        if (h) ymax_base = std::max(ymax_base, (double)h->GetMaximum());
                }
                double ymax = std::max(1e-9, 1.4 * ymax_base);

                h_data_sum->GetYaxis()->SetRangeUser(0, ymax);
                if (h_sig_sum) h_sig_sum->GetYaxis()->SetRangeUser(0, ymax);
                for (auto *h : h_bkg_sums) if (h) h->GetYaxis()->SetRangeUser(0, ymax);
                // Draw MC first (HIST), then overlay their stat error bars, then data on top
                // Use first bkg as frame for axes
                bool firstDrawn = false;
                for (auto *h : h_bkg_sums) {
                    if (!h) continue;
                    if (!firstDrawn) { h->Draw("HIST"); firstDrawn = true; }
                    else              h->Draw("HIST SAME");
                    overlayVerticalErrors(h);
                }
                if (h_sig_sum) {
                    h_sig_sum->SetLineColor(kRed); h_sig_sum->SetLineWidth(2);
                    h_sig_sum->SetLineStyle(2);    h_sig_sum->SetMarkerSize(0);
                    h_sig_sum->GetXaxis()->SetTitle(var.Data());
                    h_sig_sum->GetYaxis()->SetTitle("normalized counts");
                    if (!firstDrawn) { h_sig_sum->Draw("HIST"); firstDrawn = true; }
                    else              h_sig_sum->Draw("HIST SAME");
                    overlayVerticalErrors(h_sig_sum);
                }
                if (!firstDrawn) h_data_sum->Draw("ex0");
                else             h_data_sum->Draw("ex0 SAME");

                // Legend
                TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.88);
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.04);
                leg->AddEntry(h_data_sum, "Data", "lep");
                for (int si = 0; si < (int)selIdx.size(); ++si)
                    if (h_bkg_sums[si])
                        leg->AddEntry(h_bkg_sums[si],
                            Form("Incl. MC w=%.2f%s", weights[selIdx[si]],
                                 selIdx[si]==bestIdx ? " (best)" : ""), "l");
                if (h_sig_sum) leg->AddEntry(h_sig_sum, "Signal MC", "l");
                leg->Draw();

                myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                myText(0.20, 0.80, 1, strleg3.c_str(), 0.04);
                myText(0.20, 0.75, 1, Form("%s, %.0f < p_{T}^{clus} < %.0f GeV, cut1",
                                            var.Data(), pT_bin_edges[ipt], pT_bin_edges[ipt + 1]), 0.04);
                myText(0.20, 0.70, 1, Form("Best scan weight: w = %.3f", weights[bestIdx]), 0.04);

                coverlay->Print(outPdf.Data());

                // Cleanup this shape's temporaries
                delete h_data_sum;
                if (h_sig_sum) delete h_sig_sum;
                for (auto *h : h_bkg_sums) if (h) delete h;
            }

            coverlay->Print(Form("%s]", outPdf.Data())); // close multi-page per pT bin
        }
        delete coverlay;
    }

    // =========================================================================
    // Clean up pre-computed projections
    // =========================================================================
    for (int ishape = 0; ishape < nShapes; ++ishape)
        for (int ieta = 0; ieta < nEtaBins; ++ieta)
            for (int ipt = 0; ipt < nPtBins; ++ipt)
            {
                if (projs[ishape][ieta][ipt].data) delete projs[ishape][ieta][ipt].data;
                if (projs[ishape][ieta][ipt].nom)  delete projs[ishape][ieta][ipt].nom;
                if (projs[ishape][ieta][ipt].dbl)  delete projs[ishape][ieta][ipt].dbl;
                if (projs[ishape][ieta][ipt].sig)  delete projs[ishape][ieta][ipt].sig;
            }

    // Clean up TGraphs
    delete gr_chi2_total;
    delete gr_chi2_ndf;
    for (int is = 0; is < nShapes; ++is)
        delete gr_per_shape[is];

    // Close files
    f_data->Close();
    delete f_data;
    f_bkg_nom->Close();
    delete f_bkg_nom;
    if (hasDouble)
    {
        f_bkg_dbl->Close();
        delete f_bkg_dbl;
    }
    if (hasSig)
    {
        f_sig->Close();
        delete f_sig;
    }

    std::cout << "\nDone. Output in figures/chi2_weight_scan_*.pdf and figures/chi2_weight_scan.root" << std::endl;
}
