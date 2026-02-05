#include "plotcommon.h"
#include <vector>
#include <algorithm>

const int col_sel[] = {kBlack, kRed, kBlue, kGreen + 2, kMagenta + 2};
const int col_pt[] = {kBlack, kRed, kBlue, kGreen + 2, kMagenta + 2, kOrange + 1};

void plot_cluster_timing(const std::string &inputfile = "../efficiencytool/results/cluster_time_analysis_data.root",
                          const std::string &outputdir = "figures/",
                          const std::string &mc_inputfile = "")
{
    init_plot();

    // Determine which file to use: MC if provided, otherwise data
    std::string file_to_use = inputfile;
    bool is_mc = false;

    if (!mc_inputfile.empty())
    {
        file_to_use = mc_inputfile;
        is_mc = true;
        std::cout << "Using inclusive MC file: " << file_to_use << std::endl;
    }
    else
    {
        std::cout << "Using data file: " << file_to_use << std::endl;
    }

    // Create output filename suffix
    std::string suffix = is_mc ? "_mc" : "";

    // Open the selected file
    TFile *fin = TFile::Open(file_to_use.c_str());
    if (!fin || fin->IsZombie())
    {
        std::cerr << "ERROR: Cannot open input file: " << file_to_use << std::endl;
        return;
    }

    std::cout << "Processing file: " << file_to_use << std::endl;

    // ========================================================================
    // Part 1: Cluster-Jet Delta T for each pT bin with different selections
    // ========================================================================

    // Get the 2D histograms for different selections
    TH2D *h_all_delta_t_jet = (TH2D *)fin->Get("h_all_delta_t_jet_all");
    TH2D *h_common_delta_t_jet = (TH2D *)fin->Get("h_common_delta_t_jet_all");
    TH2D *h_tight_iso_delta_t_jet = (TH2D *)fin->Get("h_tight_iso_delta_t_jet_all");
    TH2D *h_tight_noniso_delta_t_jet = (TH2D *)fin->Get("h_tight_noniso_delta_t_jet_all");
    TH2D *h_npb_delta_t_jet = (TH2D *)fin->Get("h_npb_delta_t_jet_all");

    if (!h_all_delta_t_jet || !h_common_delta_t_jet || !h_tight_iso_delta_t_jet ||
        !h_tight_noniso_delta_t_jet || !h_npb_delta_t_jet)
    {
        std::cerr << "ERROR: Cannot find required jet-time histograms" << std::endl;
        fin->Close();
        return;
    }

    // Define pT ranges for projection
    const int n_pt_ranges = 5;
    float pt_mins[n_pt_ranges] = {5, 10, 15, 20, 30};
    float pt_maxs[n_pt_ranges] = {10, 15, 20, 30, 50};
    std::string pt_labels[n_pt_ranges] = {
        "5<p_{T}^{#gamma}<10 GeV",
        "10<p_{T}^{#gamma}<15 GeV",
        "15<p_{T}^{#gamma}<20 GeV",
        "20<p_{T}^{#gamma}<30 GeV",
        "30<p_{T}^{#gamma}<50 GeV"
    };

    // Loop over pT ranges and create plots with different selections overlaid
    for (int irange = 0; irange < n_pt_ranges; irange++)
    {
        float pT_min = pt_mins[irange];
        float pT_max = pt_maxs[irange];

        // Set pT range for all histograms
        h_all_delta_t_jet->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_common_delta_t_jet->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_tight_iso_delta_t_jet->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_tight_noniso_delta_t_jet->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_npb_delta_t_jet->GetXaxis()->SetRangeUser(pT_min, pT_max);

        // Project and normalize
        auto h_proj_norm = [](TH2D *h2, const char *name) -> TH1D * {
            TH1D *h_proj = (TH1D *)h2->ProjectionY(name);
            h_proj->Rebin(2);
            double integral = h_proj->Integral();
            if (integral > 0)
                h_proj->Scale(1.0 / integral);
            return h_proj;
        };

        TH1D *h_all_proj = h_proj_norm(h_all_delta_t_jet, Form("h_all_proj_%d", irange));
        TH1D *h_common_proj = h_proj_norm(h_common_delta_t_jet, Form("h_common_proj_%d", irange));
        TH1D *h_tight_iso_proj = h_proj_norm(h_tight_iso_delta_t_jet, Form("h_tight_iso_proj_%d", irange));
        TH1D *h_tight_noniso_proj = h_proj_norm(h_tight_noniso_delta_t_jet, Form("h_tight_noniso_proj_%d", irange));
        TH1D *h_npb_proj = h_proj_norm(h_npb_delta_t_jet, Form("h_npb_proj_%d", irange));

        // Find y-range
        double ymax = std::max({
            h_all_proj->GetMaximum(),
            h_common_proj->GetMaximum(),
            h_tight_iso_proj->GetMaximum(),
            h_tight_noniso_proj->GetMaximum(),
            h_npb_proj->GetMaximum()
        });

        TCanvas *c1 = new TCanvas(Form("c_jet_time_%d", irange), "", 600, 600);
        c1->SetLogy();
        TH1D *h_frame = new TH1D(Form("h_frame_jet_%d", irange), "", 80, -30, 20);
        h_frame->GetYaxis()->SetRangeUser(1e-5, 5.0 * ymax);
        h_frame->SetXTitle("Cluster - Jet Time [ns]");
        h_frame->SetYTitle("Normalized Counts");
        h_frame->Draw("axis");

        // Draw histograms
        h_all_proj->SetLineColor(col_sel[0]);
        h_all_proj->SetMarkerColor(col_sel[0]);
        h_all_proj->SetLineWidth(2);
        h_all_proj->Draw("same hist");

        h_common_proj->SetLineColor(col_sel[1]);
        h_common_proj->SetMarkerColor(col_sel[1]);
        h_common_proj->SetLineWidth(2);
        h_common_proj->Draw("same hist");

        h_tight_iso_proj->SetLineColor(col_sel[2]);
        h_tight_iso_proj->SetMarkerColor(col_sel[2]);
        h_tight_iso_proj->SetLineWidth(2);
        h_tight_iso_proj->Draw("same hist");

        h_tight_noniso_proj->SetLineColor(col_sel[3]);
        h_tight_noniso_proj->SetMarkerColor(col_sel[3]);
        h_tight_noniso_proj->SetLineWidth(2);
        h_tight_noniso_proj->Draw("same hist");

        h_npb_proj->SetLineColor(col_sel[4]);
        h_npb_proj->SetMarkerColor(col_sel[4]);
        h_npb_proj->SetLineWidth(2);
        h_npb_proj->Draw("same hist");

        // Add labels
        myMarkerLineText(0.25, 0.88, 1, col_sel[0], 20, col_sel[0], 1, "All clusters", 0.045, true);
        myMarkerLineText(0.25, 0.83, 1, col_sel[1], 20, col_sel[1], 1, "Common selection", 0.045, true);
        myMarkerLineText(0.25, 0.78, 1, col_sel[2], 20, col_sel[2], 1, "Tight+Iso", 0.045, true);
        myMarkerLineText(0.25, 0.73, 1, col_sel[3], 20, col_sel[3], 1, "Tight+NonIso", 0.045, true);
        myMarkerLineText(0.25, 0.68, 1, col_sel[4], 20, col_sel[4], 1, "NPB", 0.045, true);
        myText(0.65, 0.88, 1, strleg1.c_str(), 0.045);
        if (is_mc)
        {
            myText(0.65, 0.83, 1, strIncMC.c_str(), 0.045);
            myText(0.65, 0.78, 1, "|#eta^{#gamma}|<0.7", 0.045);
            myText(0.65, 0.73, 1, pt_labels[irange].c_str(), 0.045);
        }
        else
        {
            myText(0.65, 0.78, 1, "|#eta^{#gamma}|<0.7", 0.045);
            myText(0.65, 0.73, 1, pt_labels[irange].c_str(), 0.045);
        }

        c1->SaveAs(Form("%s/cluster_jet_delta_t_pt_%.0f_%.0f%s.pdf", outputdir.c_str(), pT_min, pT_max, suffix.c_str()));

        delete h_frame;
        delete h_all_proj;
        delete h_common_proj;
        delete h_tight_iso_proj;
        delete h_tight_noniso_proj;
        delete h_npb_proj;
        delete c1;
    }

    std::cout << "Created " << n_pt_ranges << " cluster-jet delta T plots for different pT bins" << std::endl;

    // ========================================================================
    // Part 2: Cluster-Jet Delta T for different pT ranges (all clusters only)
    // ========================================================================

    TCanvas *c_pt_comparison = new TCanvas("c_pt_comparison", "", 600, 600);
    c_pt_comparison->SetLogy();

    std::vector<TH1D *> projections;
    double y_max_pt = 0;

    for (int irange = 0; irange < n_pt_ranges; irange++)
    {
        float pT_min = pt_mins[irange];
        float pT_max = pt_maxs[irange];

        h_all_delta_t_jet->GetXaxis()->SetRangeUser(pT_min, pT_max);
        TH1D *h_proj = (TH1D *)h_all_delta_t_jet->ProjectionY(Form("h_all_pt_proj_%d", irange));
        h_proj->Rebin(2);
        if (h_proj->Integral() > 0)
            h_proj->Scale(1.0 / h_proj->Integral());

        h_proj->SetLineColor(col_pt[irange]);
        h_proj->SetMarkerColor(col_pt[irange]);
        h_proj->SetLineWidth(2);

        if (h_proj->GetMaximum() > y_max_pt)
            y_max_pt = h_proj->GetMaximum();

        projections.push_back(h_proj);
    }

    TH1D *h_frame_pt = new TH1D("h_frame_pt", "", 80, -20, 20);
    h_frame_pt->GetYaxis()->SetRangeUser(1e-4, 1.2 * y_max_pt);
    h_frame_pt->SetXTitle("Cluster - Jet Time [ns]");
    h_frame_pt->SetYTitle("Normalized Counts");
    h_frame_pt->Draw("axis");

    for (size_t i = 0; i < projections.size(); i++)
    {
        projections[i]->Draw("same hist");
    }

    for (int i = 0; i < n_pt_ranges; i++)
    {
        myMarkerLineText(0.55, 0.88 - i * 0.05, 1, col_pt[i], 20, col_pt[i], 1,
                         pt_labels[i].c_str(), 0.04, true);
    }
    myText(0.25, 0.88, 1, strleg1.c_str(), 0.045);
    myText(0.25, 0.83, 1, "All clusters", 0.045);
    if (is_mc)
    {
        myText(0.25, 0.78, 1, strIncMC.c_str(), 0.045);
        myText(0.25, 0.73, 1, "|#eta^{#gamma}|<0.7", 0.045);
    }
    else
    {
        myText(0.25, 0.78, 1, "|#eta^{#gamma}|<0.7", 0.045);
    }

    c_pt_comparison->SaveAs(Form("%s/cluster_jet_delta_t_pt_comparison%s.pdf", outputdir.c_str(), suffix.c_str()));

    delete h_frame_pt;
    for (auto h : projections)
        delete h;
    delete c_pt_comparison;

    std::cout << "Created cluster-jet delta T pT comparison plot" << std::endl;

    // ========================================================================
    // Part 2b: Cluster-MBD Delta T vs. pT 2D plots (All, NPB) -> separate output files
    // ========================================================================

    TH2D *h_all_delta_t_mbd_all = (TH2D *)fin->Get("h_all_delta_t_mbd_all");
    TH2D *h_npb_delta_t_mbd_all = (TH2D *)fin->Get("h_npb_delta_t_mbd_all");

    if (!h_all_delta_t_mbd_all || !h_npb_delta_t_mbd_all)
    {
        std::cerr << "ERROR: Cannot find required MBD-time vs pT histograms (h_all_delta_t_mbd_all, h_npb_delta_t_mbd_all)" << std::endl;
        fin->Close();
        return;
    }

    {
        // --- All clusters ---
        TCanvas *c_mbd_pt_2d_all = new TCanvas("c_mbd_pt_2d_all", "", 700, 600);
        c_mbd_pt_2d_all->SetLeftMargin(0.12);
        c_mbd_pt_2d_all->SetRightMargin(0.16);
        c_mbd_pt_2d_all->SetBottomMargin(0.12);
        c_mbd_pt_2d_all->SetLogz();

        h_all_delta_t_mbd_all->SetTitle(";Cluster p_{T} [GeV];Cluster - MBD Time [ns]");
        h_all_delta_t_mbd_all->GetYaxis()->SetRangeUser(-20, 20);
        h_all_delta_t_mbd_all->GetXaxis()->SetTitleSize(0.045);
        h_all_delta_t_mbd_all->GetYaxis()->SetTitleSize(0.045);
        h_all_delta_t_mbd_all->GetYaxis()->SetTitleOffset(1.2);
        h_all_delta_t_mbd_all->Draw("colz");

        myText(0.16, 0.92, 1, strleg1.c_str(), 0.045);
        if (is_mc)
            myText(0.16, 0.86, 1, strIncMC.c_str(), 0.045);
        myText(0.16, 0.80, 1, "All clusters (p_{T}>10 GeV)", 0.045);
        myText(0.16, 0.74, 1, "|#eta^{#gamma}|<0.7", 0.045);

        c_mbd_pt_2d_all->SaveAs(Form("%s/cluster_mbd_time_vs_pt_2d_all%s.pdf", outputdir.c_str(), suffix.c_str()));
        c_mbd_pt_2d_all->SaveAs(Form("%s/cluster_mbd_time_vs_pt_2d_all%s.png", outputdir.c_str(), suffix.c_str()));
        delete c_mbd_pt_2d_all;

        // --- NPB selection ---
        TCanvas *c_mbd_pt_2d_npb = new TCanvas("c_mbd_pt_2d_npb", "", 700, 600);
        c_mbd_pt_2d_npb->SetLeftMargin(0.12);
        c_mbd_pt_2d_npb->SetRightMargin(0.16);
        c_mbd_pt_2d_npb->SetBottomMargin(0.12);
        c_mbd_pt_2d_npb->SetLogz();

        h_npb_delta_t_mbd_all->SetTitle(";Cluster p_{T} [GeV];Cluster - MBD Time [ns]");
        h_npb_delta_t_mbd_all->GetYaxis()->SetRangeUser(-20, 20);
        h_npb_delta_t_mbd_all->GetXaxis()->SetTitleSize(0.045);
        h_npb_delta_t_mbd_all->GetYaxis()->SetTitleSize(0.045);
        h_npb_delta_t_mbd_all->GetYaxis()->SetTitleOffset(1.2);
        h_npb_delta_t_mbd_all->Draw("colz");

        myText(0.16, 0.92, 1, strleg1.c_str(), 0.045);
        if (is_mc)
            myText(0.16, 0.86, 1, strIncMC.c_str(), 0.045);
        myText(0.16, 0.80, 1, "NPB selection (p_{T}>10 GeV)", 0.045);
        myText(0.16, 0.74, 1, "|#eta^{#gamma}|<0.7", 0.045);

        c_mbd_pt_2d_npb->SaveAs(Form("%s/cluster_mbd_time_vs_pt_2d_npb%s.pdf", outputdir.c_str(), suffix.c_str()));
        c_mbd_pt_2d_npb->SaveAs(Form("%s/cluster_mbd_time_vs_pt_2d_npb%s.png", outputdir.c_str(), suffix.c_str()));
        delete c_mbd_pt_2d_npb;

        std::cout << "Created cluster-MBD time vs pT 2D plots (separate files: all, npb)" << std::endl;
    }

    // ========================================================================
    // Part 3: Cluster-MBD Time vs. Eta 2D plot with profiles
    // ========================================================================

    TH2D *h_all_mbd_vs_eta = (TH2D *)fin->Get("h_all_delta_t_mbd_vs_eta");
    TH2D *h_npb_mbd_vs_eta = (TH2D *)fin->Get("h_npb_delta_t_mbd_vs_eta");

    if (!h_all_mbd_vs_eta || !h_npb_mbd_vs_eta)
    {
        std::cerr << "ERROR: Cannot find MBD-time vs eta histograms" << std::endl;
        fin->Close();
        return;
    }

    // ========================================================================
    // Part 3b: NPB selection - overlay delta_t(MBD) projections for eta slices
    // ========================================================================
    {
        // Define eta slices (adjust as needed)
        const int n_eta_slices = 4;
        const double eta_edges[n_eta_slices + 1] = {-1.0, -0.35, 0.35, 0.7, 1.0};
        const std::string eta_labels[n_eta_slices] = {
            "-1.0<#eta<-0.35",
            "-0.35<#eta<0.35",
            "0.35<#eta<0.7",
            "0.7<#eta<1.0"};

        TCanvas *c_eta_slices = new TCanvas("c_mbd_dt_eta_slices", "", 650, 650);
        //c_eta_slices->SetLogy();

        std::vector<TH1D *> h_slices;
        double ymax_eta = 0.0;

        for (int islice = 0; islice < n_eta_slices; ++islice)
        {
            const double eta_min = eta_edges[islice];
            const double eta_max = eta_edges[islice + 1];

            h_npb_mbd_vs_eta->GetXaxis()->SetRangeUser(eta_min, eta_max);
            TH1D *h_proj = (TH1D *)h_npb_mbd_vs_eta->ProjectionY(Form("h_npb_dt_mbd_eta_%d", islice));
            h_proj->Rebin(2);
            const double integral = h_proj->Integral();
            if (integral > 0)
            {
                h_proj->Scale(1.0 / integral);
            }

            const int col = col_pt[islice % (int)(sizeof(col_pt) / sizeof(col_pt[0]))];
            h_proj->SetLineColor(col);
            h_proj->SetMarkerColor(col);
            h_proj->SetLineWidth(2);

            ymax_eta = std::max(ymax_eta, h_proj->GetMaximum());
            h_slices.push_back(h_proj);
        }

        // Reset range so subsequent ProfileX uses full eta range
        h_npb_mbd_vs_eta->GetXaxis()->SetRange(0, 0);

        TH1D *h_frame_eta = new TH1D("h_frame_eta_slices", "", 80, -20, 20);
        h_frame_eta->SetXTitle("Cluster - MBD Time [ns]");
        h_frame_eta->SetYTitle("Normalized Counts");
        h_frame_eta->GetYaxis()->SetRangeUser(1e-5, 1.2 * std::max(1e-5, ymax_eta));
        h_frame_eta->Draw("axis");

        for (size_t i = 0; i < h_slices.size(); ++i)
        {
            h_slices[i]->Draw("same hist");
        }

        myText(0.20, 0.92, 1, strleg1.c_str(), 0.045);
        if (is_mc)
        {
            myText(0.20, 0.87, 1, strIncMC.c_str(), 0.045);
        }
        myText(0.20, 0.82, 1, "NPB selection", 0.045);

        for (int islice = 0; islice < n_eta_slices; ++islice)
        {
            const int col = col_pt[islice % (int)(sizeof(col_pt) / sizeof(col_pt[0]))];
            myMarkerLineText(0.65, 0.88 - islice * 0.05, 1, col, 20, col, 1, eta_labels[islice].c_str(), 0.04, true);
        }

        c_eta_slices->SaveAs(Form("%s/cluster_mbd_delta_t_eta_slices_npb%s.pdf", outputdir.c_str(), suffix.c_str()));

        delete h_frame_eta;
        for (auto h : h_slices)
            delete h;
        delete c_eta_slices;

        std::cout << "Created NPB cluster-MBD delta T overlay for eta slices" << std::endl;
    }

    TCanvas *c_mbd_2d = new TCanvas("c_mbd_2d", "", 800, 600);
    c_mbd_2d->SetLeftMargin(0.12);
    c_mbd_2d->SetRightMargin(0.15);
    c_mbd_2d->SetLogz();

    h_all_mbd_vs_eta->SetTitle(";Cluster #eta;Cluster - MBD Time [ns]");
    h_all_mbd_vs_eta->GetXaxis()->SetTitleSize(0.045);
    h_all_mbd_vs_eta->GetYaxis()->SetTitleSize(0.045);
    h_all_mbd_vs_eta->GetYaxis()->SetTitleOffset(1.2);
    h_all_mbd_vs_eta->GetYaxis()->SetRangeUser(-20, 20);
    h_all_mbd_vs_eta->GetZaxis()->SetRangeUser(100, 1.2 * h_all_mbd_vs_eta->GetMaximum());
    h_all_mbd_vs_eta->Draw("colz");

    // Create profiles
    TProfile *prof_all = h_all_mbd_vs_eta->ProfileX("prof_all_mbd", 1, -1, "");
    TProfile *prof_npb = h_npb_mbd_vs_eta->ProfileX("prof_npb_mbd", 1, -1, "");

    prof_all->SetLineColor(kRed);
    prof_all->SetLineWidth(3);
    prof_all->SetMarkerColor(kRed);
    prof_all->SetMarkerStyle(20);
    prof_all->SetMarkerSize(0.8);

    prof_npb->SetLineColor(kBlue);
    prof_npb->SetLineWidth(3);
    prof_npb->SetMarkerColor(kBlue);
    prof_npb->SetMarkerStyle(21);
    prof_npb->SetMarkerSize(0.8);

    prof_all->Draw("same");
    prof_npb->Draw("same");

    myMarkerLineText(0.18, 0.88, 1, kRed, 20, kRed, 1, "All Clusters (Mean)", 0.04, true);
    myMarkerLineText(0.18, 0.83, 1, kBlue, 21, kBlue, 1, "NPB Selection (Mean)", 0.04, true);
    myText(0.18, 0.78, 1, strleg1.c_str(), 0.04);
    if (is_mc)
    {
        myText(0.18, 0.73, 1, strIncMC.c_str(), 0.04);
        myText(0.18, 0.68, 1, "p_{T} > 10 GeV", 0.04);
    }
    else
    {
        myText(0.18, 0.73, 1, "p_{T} > 10 GeV", 0.04);
    }

    c_mbd_2d->SaveAs(Form("%s/cluster_mbd_time_vs_eta_2d%s.pdf", outputdir.c_str(), suffix.c_str()));
    c_mbd_2d->SaveAs(Form("%s/cluster_mbd_time_vs_eta_2d%s.png", outputdir.c_str(), suffix.c_str()));

    std::cout << "Created cluster-MBD time vs eta 2D plot" << std::endl;

    // ========================================================================
    // Part 4: Mean time profiles only
    // ========================================================================

    TCanvas *c_profiles = new TCanvas("c_profiles", "", 600, 600);

    TH1D *h_frame_prof = new TH1D("h_frame_prof", "", 50, -1.0, 1.0);
    h_frame_prof->SetTitle(";Cluster #eta;Mean Cluster - MBD Time [ns]");
    h_frame_prof->GetXaxis()->SetTitleSize(0.045);
    h_frame_prof->GetYaxis()->SetTitleSize(0.045);
    h_frame_prof->GetYaxis()->SetTitleOffset(1.4);
    h_frame_prof->SetMinimum(-10);
    h_frame_prof->SetMaximum(10);
    h_frame_prof->Draw("axis");

    prof_all->Draw("same");
    prof_npb->Draw("same");

    myMarkerLineText(0.25, 0.88, 1, kRed, 20, kRed, 1, "All Clusters (Mean)", 0.045, true);
    myMarkerLineText(0.25, 0.83, 1, kBlue, 21, kBlue, 1, "NPB Selection (Mean)", 0.045, true);
    myText(0.25, 0.78, 1, strleg1.c_str(), 0.045);
    if (is_mc)
    {
        myText(0.25, 0.73, 1, strIncMC.c_str(), 0.045);
        myText(0.25, 0.68, 1, "p_{T} > 10 GeV", 0.045);
    }
    else
    {
        myText(0.25, 0.73, 1, "p_{T} > 10 GeV", 0.045);
    }

    c_profiles->SaveAs(Form("%s/cluster_mbd_time_vs_eta_profiles%s.pdf", outputdir.c_str(), suffix.c_str()));
    c_profiles->SaveAs(Form("%s/cluster_mbd_time_vs_eta_profiles%s.png", outputdir.c_str(), suffix.c_str()));

    std::cout << "Created cluster-MBD time vs eta profiles plot" << std::endl;

    // ========================================================================
    // Summary
    // ========================================================================

    std::cout << "\n========================================" << std::endl;
    std::cout << "Summary:" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Created plots:" << std::endl;
    std::cout << "  - " << n_pt_ranges << " cluster-jet delta T plots (different selections per pT bin)" << std::endl;
    std::cout << "  - 1 cluster-jet delta T pT comparison plot (all clusters)" << std::endl;
    std::cout << "  - 1 cluster-MBD time vs eta 2D plot with profiles" << std::endl;
    std::cout << "  - 1 cluster-MBD time vs eta profiles only plot" << std::endl;
    std::cout << "========================================" << std::endl;

    delete h_frame_prof;
    delete c_profiles;
    delete c_mbd_2d;

    fin->Close();
}
