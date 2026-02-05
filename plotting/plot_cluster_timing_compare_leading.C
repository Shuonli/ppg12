#include "plotcommon.h"
#include <algorithm>

// Compare cluster timing definitions in MC by overlaying cluster-MBD time distributions:
// - energy-weighted average tower time (default)
// - leading tower time
//
// Expected inputs are the MERGED MC outputs produced by efficiencytool/merge_cluster_time.C
// with use_leading_tower_time = false/true.
void plot_cluster_timing_compare_leading(
    const std::string &mc_inputfile_default = "../efficiencytool/results/cluster_time_analysis_photon.root",
    const std::string &mc_inputfile_leading = "../efficiencytool/results/cluster_time_analysis_photon_leadingTowerTime.root",
    const std::string &mc_inputfile_leading_corr = "../efficiencytool/results/cluster_time_analysis_photon_leadingTowerTimeCorr.root",
    const std::string &outputdir = "figures/")
{
    init_plot();

    std::cout << "==================================================" << std::endl;
    std::cout << "  Cluster Timing (MC) Compare: Default vs Leading" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "Default file: " << mc_inputfile_default << std::endl;
    std::cout << "Leading file: " << mc_inputfile_leading << std::endl;
    std::cout << "Leading+corr: " << mc_inputfile_leading_corr << std::endl;
    std::cout << "Output dir  : " << outputdir << std::endl;
    std::cout << std::endl;

    // Open files
    TFile *f_default = TFile::Open(mc_inputfile_default.c_str());
    if (!f_default || f_default->IsZombie())
    {
        std::cerr << "ERROR: Cannot open default input file: " << mc_inputfile_default << std::endl;
        return;
    }

    TFile *f_leading = TFile::Open(mc_inputfile_leading.c_str());
    if (!f_leading || f_leading->IsZombie())
    {
        std::cerr << "ERROR: Cannot open leading input file: " << mc_inputfile_leading << std::endl;
        f_default->Close();
        return;
    }

    TFile *f_leading_corr = TFile::Open(mc_inputfile_leading_corr.c_str());
    if (!f_leading_corr || f_leading_corr->IsZombie())
    {
        std::cerr << "ERROR: Cannot open leading+corr input file: " << mc_inputfile_leading_corr << std::endl;
        f_leading->Close();
        f_default->Close();
        return;
    }

    // Grab 2D histograms: eta vs (cluster - MBD time)
    TH2D *h_all_mbd_vs_eta_default = (TH2D *)f_default->Get("h_all_delta_t_mbd_vs_eta");
    TH2D *h_all_mbd_vs_eta_leading = (TH2D *)f_leading->Get("h_all_delta_t_mbd_vs_eta");
    TH2D *h_all_mbd_vs_eta_leading_corr = (TH2D *)f_leading_corr->Get("h_all_delta_t_mbd_vs_eta");

    if (!h_all_mbd_vs_eta_default || !h_all_mbd_vs_eta_leading || !h_all_mbd_vs_eta_leading_corr)
    {
        std::cerr << "ERROR: Missing h_all_delta_t_mbd_vs_eta in one or both files" << std::endl;
        if (!h_all_mbd_vs_eta_default)
            std::cerr << "  - Missing in default file" << std::endl;
        if (!h_all_mbd_vs_eta_leading)
            std::cerr << "  - Missing in leading file" << std::endl;
        if (!h_all_mbd_vs_eta_leading_corr)
            std::cerr << "  - Missing in leading+corr file" << std::endl;
        f_leading_corr->Close();
        f_leading->Close();
        f_default->Close();
        return;
    }

    // Project to 1D delta_t distributions (inclusive in eta)
    TH1D *h_dt_default = (TH1D *)h_all_mbd_vs_eta_default->ProjectionY("h_dt_mbd_default");
    TH1D *h_dt_leading = (TH1D *)h_all_mbd_vs_eta_leading->ProjectionY("h_dt_mbd_leading");
    TH1D *h_dt_leading_corr = (TH1D *)h_all_mbd_vs_eta_leading_corr->ProjectionY("h_dt_mbd_leading_corr");

    // Normalize
    h_dt_default->Rebin(1);
    h_dt_leading->Rebin(1);
    h_dt_leading_corr->Rebin(1);
    if (h_dt_default->Integral() > 0)
        h_dt_default->Scale(1.0 / h_dt_default->Integral());
    if (h_dt_leading->Integral() > 0)
        h_dt_leading->Scale(1.0 / h_dt_leading->Integral());
    if (h_dt_leading_corr->Integral() > 0)
        h_dt_leading_corr->Scale(1.0 / h_dt_leading_corr->Integral());

    const double xmax = 20.0;
    const double xmin = -20.0;
    h_dt_default->GetXaxis()->SetRangeUser(xmin, xmax);
    h_dt_leading->GetXaxis()->SetRangeUser(xmin, xmax);
    h_dt_leading_corr->GetXaxis()->SetRangeUser(xmin, xmax);

    const double ymax = std::max({h_dt_default->GetMaximum(), h_dt_leading->GetMaximum(), h_dt_leading_corr->GetMaximum()});

    TCanvas *c = new TCanvas("c_cluster_mbd_dt_compare", "", 650, 650);
    c->SetLogy();

    TH1D *h_frame = new TH1D("h_frame_cluster_mbd_dt_compare", "", 80, xmin, xmax);
    h_frame->SetXTitle("Cluster - MBD Time [ns]");
    h_frame->SetYTitle("Normalized Counts");
    h_frame->GetYaxis()->SetRangeUser(1e-5, 2.0 * std::max(1e-5, ymax));
    h_frame->Draw("axis");

    // Style + draw
    h_dt_default->SetLineColor(kRed);
    h_dt_default->SetMarkerColor(kRed);
    h_dt_default->SetLineWidth(2);
    h_dt_default->Draw("same hist");

    h_dt_leading->SetLineColor(kBlue);
    h_dt_leading->SetMarkerColor(kBlue);
    h_dt_leading->SetLineWidth(2);
    h_dt_leading->Draw("same hist");

    h_dt_leading_corr->SetLineColor(kGreen + 2);
    h_dt_leading_corr->SetMarkerColor(kGreen + 2);
    h_dt_leading_corr->SetLineWidth(2);
    h_dt_leading_corr->Draw("same hist");

    // Labels
    myText(0.20, 0.92, 1, strleg1.c_str(), 0.045);
    myText(0.20, 0.87, 1, strIncMC.c_str(), 0.045);
    myText(0.20, 0.82, 1, "All clusters", 0.045);

    myMarkerLineText(0.60, 0.88, 1, kRed, 20, kRed, 1, "E-weighted avg time", 0.04, true);
    myMarkerLineText(0.60, 0.83, 1, kBlue, 20, kBlue, 1, "Leading tower time", 0.04, true);
    myMarkerLineText(0.60, 0.78, 1, kGreen + 2, 20, kGreen + 2, 1, "Leading tower time (corr)", 0.04, true);
    myText(0.20, 0.78, 1, Form("RMS (E-weighted) = %.3f ns", h_dt_default->GetRMS()), 0.04);
    myText(0.20, 0.73, 1, Form("RMS (Leading)   = %.3f ns", h_dt_leading->GetRMS()), 0.04);
    myText(0.20, 0.68, 1, Form("RMS (Lead+corr) = %.3f ns", h_dt_leading_corr->GetRMS()), 0.04);

    c->SaveAs(Form("%s/cluster_mbd_delta_t_compare_leadingTowerTime.pdf", outputdir.c_str()));

    // Cleanup
    delete h_frame;
    delete c;
    delete h_dt_default;
    delete h_dt_leading;
    delete h_dt_leading_corr;

    f_leading_corr->Close();
    f_leading->Close();
    f_default->Close();
}


