#include "plotcommon.h"
#include <vector>
#include <algorithm>

const int col_ratio[] = {kAzure + 2, kPink + 5, kSpring - 6, kRed - 4, kBlue - 4};

void plot_cluster_mbd_time()
{
    bool plot_sim = false;
    init_plot();

    TFile *f_in_time = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_showershape.root");
    if (plot_sim)
    {
        f_in_time = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_bdt_none.root");
    }
    // Histogram pointers
    TH2D* h_delta_t_tight_cluster = (TH2D*)f_in_time->Get("h_delta_t_tight_cluster_0");
    TH2D* h_delta_t_tight_iso_cluster = (TH2D*)f_in_time->Get("h_delta_t_tight_iso_cluster_0");
    TH2D* h_delta_t_tight_noniso_cluster = (TH2D*)f_in_time->Get("h_delta_t_tight_noniso_cluster_0");
    TH2D* h_delta_t_nontight_iso_cluster = (TH2D*)f_in_time->Get("h_delta_t_nontight_iso_cluster_0");
    TH2D* h_delta_t_nontight_noniso_cluster = (TH2D*)f_in_time->Get("h_delta_t_nontight_noniso_cluster_0");
    TH2D* h_delta_t_common_cluster = (TH2D*)f_in_time->Get("h_delta_t_common_cluster_0");
    TH2D* h_delta_t_all_cluster = (TH2D*)f_in_time->Get("h_delta_t_all_cluster_0");
    TH2D* h_delta_t_tight_cluster_b2bjet = (TH2D*)f_in_time->Get("h_delta_t_tight_cluster_b2bjet_0");
    TH2D* h_delta_t_nontight_iso_cluster_b2bjet = (TH2D*)f_in_time->Get("h_delta_t_nontight_iso_cluster_b2bjet_0");
    TH2D* h_delta_t_nontight_noniso_cluster_b2bjet = (TH2D*)f_in_time->Get("h_delta_t_nontight_noniso_cluster_b2bjet_0");
    TH2D* h_delta_t_npb_cluster = (TH2D*)f_in_time->Get("h_delta_t_npb_cluster_0");
    TH2D* h_mbd_t_vs_cluster_t = (TH2D*)f_in_time->Get("h_mbd_t_vs_cluster_t");

    // Define pT ranges
    const int n_ranges = 3;
    float pt_mins[n_ranges] = {10, 15, 20};
    float pt_maxs[n_ranges] = {15, 20, 30};
    std::string pt_labels[n_ranges] = {
        "10.0<p_{T}^{#gamma}<15.0GeV",
        "15.0<p_{T}^{#gamma}<20.0GeV",
        "20.0<p_{T}^{#gamma}<30.0GeV"
    };
    std::string pt_fnames[n_ranges] = {
        "cluster_mbd_time_10_15.pdf",
        "cluster_mbd_time_15_20.pdf",
        "cluster_mbd_time_20_30.pdf"
    };
    if (plot_sim)
    {
        pt_fnames[0] = "cluster_mbd_time_10_15_sim.pdf";
        pt_fnames[1] = "cluster_mbd_time_15_20_sim.pdf";
        pt_fnames[2] = "cluster_mbd_time_20_30_sim.pdf";
    }

    std::string pt_sideband_fnames[n_ranges] = {
        "cluster_mbd_time_10_15_sideband.pdf",
        "cluster_mbd_time_15_20_sideband.pdf",
        "cluster_mbd_time_20_30_sideband.pdf"
    };
    if (plot_sim)
    {
        pt_sideband_fnames[0] = "cluster_mbd_time_10_15_sideband_sim.pdf";
        pt_sideband_fnames[1] = "cluster_mbd_time_15_20_sideband_sim.pdf";
        pt_sideband_fnames[2] = "cluster_mbd_time_20_30_sideband_sim.pdf";
    }

    std::string pt_b2bjet_fnames[n_ranges] = {
        "cluster_mbd_time_10_15_b2bjet.pdf",
        "cluster_mbd_time_15_20_b2bjet.pdf",
        "cluster_mbd_time_20_30_b2bjet.pdf"
    };
    if (plot_sim)
    {
        pt_b2bjet_fnames[0] = "cluster_mbd_time_10_15_b2bjet_sim.pdf";
        pt_b2bjet_fnames[1] = "cluster_mbd_time_15_20_b2bjet_sim.pdf";
        pt_b2bjet_fnames[2] = "cluster_mbd_time_20_30_b2bjet_sim.pdf";
    }

    const std::string ana_label = plot_sim ? "inclusive MC" : "ana521";
    std::string mbd_time_fname = plot_sim ? "mbd_time_sim.pdf" : "mbd_time.pdf";

    for(int irange=0; irange<n_ranges; ++irange)
    {
        float pT_min = pt_mins[irange];
        float pT_max = pt_maxs[irange];

        // Set pT range for all used histograms
        h_delta_t_tight_cluster->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_delta_t_tight_iso_cluster->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_delta_t_tight_noniso_cluster->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_delta_t_nontight_iso_cluster->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_delta_t_nontight_noniso_cluster->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_delta_t_common_cluster->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_delta_t_all_cluster->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_delta_t_tight_cluster_b2bjet->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_delta_t_nontight_iso_cluster_b2bjet->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_delta_t_nontight_noniso_cluster_b2bjet->GetXaxis()->SetRangeUser(pT_min, pT_max);
        h_delta_t_npb_cluster->GetXaxis()->SetRangeUser(pT_min, pT_max);
        // Project and normalize
        auto h_proj_norm = [](TH2D* h2, const char* name) -> TH1D* {
            TH1D* h_proj = (TH1D*)h2->ProjectionY(name);
            h_proj->Rebin(2);
            double integral = h_proj->Integral();
            if(integral > 0) h_proj->Scale(1.0 / integral);
            return h_proj;
        };

        TH1D* h_delta_t_tight_cluster_proj = h_proj_norm(h_delta_t_tight_cluster, Form("h_delta_t_tight_cluster_proj_%d",irange));
        TH1D* h_delta_t_tight_iso_cluster_proj = h_proj_norm(h_delta_t_tight_iso_cluster, Form("h_delta_t_tight_iso_cluster_proj_%d",irange));
        TH1D* h_delta_t_tight_noniso_cluster_proj = h_proj_norm(h_delta_t_tight_noniso_cluster, Form("h_delta_t_tight_noniso_cluster_proj_%d",irange));
        TH1D* h_delta_t_nontight_iso_cluster_proj = h_proj_norm(h_delta_t_nontight_iso_cluster, Form("h_delta_t_nontight_iso_cluster_proj_%d",irange));
        TH1D* h_delta_t_nontight_noniso_cluster_proj = h_proj_norm(h_delta_t_nontight_noniso_cluster, Form("h_delta_t_nontight_noniso_cluster_proj_%d",irange));
        TH1D* h_delta_t_common_cluster_proj = h_proj_norm(h_delta_t_common_cluster, Form("h_delta_t_common_cluster_proj_%d",irange));
        TH1D* h_delta_t_all_cluster_proj = h_proj_norm(h_delta_t_all_cluster, Form("h_delta_t_all_cluster_proj_%d",irange));
        TH1D* h_delta_t_tight_cluster_b2bjet_proj = h_proj_norm(h_delta_t_tight_cluster_b2bjet, Form("h_delta_t_tight_cluster_b2bjet_proj_%d",irange));
        TH1D* h_delta_t_nontight_iso_cluster_b2bjet_proj = h_proj_norm(h_delta_t_nontight_iso_cluster_b2bjet, Form("h_delta_t_nontight_iso_cluster_b2bjet_proj_%d",irange));
        TH1D* h_delta_t_nontight_noniso_cluster_b2bjet_proj = h_proj_norm(h_delta_t_nontight_noniso_cluster_b2bjet, Form("h_delta_t_nontight_noniso_cluster_b2bjet_proj_%d",irange));
        TH1D* h_delta_t_npb_cluster_proj = h_proj_norm(h_delta_t_npb_cluster, Form("h_delta_t_npb_cluster_proj_%d",irange));
        // Prepare y-range
        double ymax = std::max({
            h_delta_t_tight_cluster_proj->GetMaximum(),
            h_delta_t_common_cluster_proj->GetMaximum(),
            h_delta_t_all_cluster_proj->GetMaximum(),
            h_delta_t_tight_cluster_b2bjet_proj->GetMaximum(),
            h_delta_t_npb_cluster_proj->GetMaximum()
        });

        TCanvas *c1 = new TCanvas(Form("c1_%d",irange), "", 600, 600);
        TH1D* h_frame = new TH1D(Form("h_frame_%d",irange), "h_frame", 160, -20, 10);
        h_frame->GetYaxis()->SetRangeUser(0, 1.2*ymax);
        h_frame->SetXTitle("Cluster Time - MBD Time [ns]");
        h_frame->Draw("axis");
        h_delta_t_all_cluster_proj->SetLineColor(kBlack);
        h_delta_t_all_cluster_proj->SetMarkerColor(kBlack);
        h_delta_t_all_cluster_proj->Draw("same");
        h_delta_t_common_cluster_proj->SetLineColor(kRed);
        h_delta_t_common_cluster_proj->SetMarkerColor(kRed);
        h_delta_t_common_cluster_proj->Draw("same");
        h_delta_t_tight_cluster_proj->SetLineColor(kBlue);
        h_delta_t_tight_cluster_proj->SetMarkerColor(kBlue);
        h_delta_t_tight_cluster_proj->Draw("same");
        h_delta_t_tight_cluster_b2bjet_proj->SetLineColor(kPink + 5);
        h_delta_t_tight_cluster_b2bjet_proj->SetMarkerColor(kPink + 5);
        h_delta_t_tight_cluster_b2bjet_proj->Draw("same");
        h_delta_t_npb_cluster_proj->SetLineColor(kSpring - 6);
        h_delta_t_npb_cluster_proj->SetMarkerColor(kSpring - 6);
        h_delta_t_npb_cluster_proj->Draw("same");

        myMarkerLineText(0.25, 0.9, 1, kBlack, 20, kBlack, 1, "All clusters", 0.05, true);
        myMarkerLineText(0.25, 0.85, 1, kRed, 20, kRed, 1, "pre-selected clusters", 0.05, true);
        myMarkerLineText(0.25, 0.8, 1, kBlue, 20, kBlue, 1, "Tight clusters", 0.05, true);
        myMarkerLineText(0.25, 0.75, 1, kPink + 5, 20, kPink + 5, 1, "Tight w/ B2BJet", 0.05, true);
        myMarkerLineText(0.25, 0.70, 1, kSpring - 6, 20, kSpring - 6, 1, "NPB clusters", 0.05, true);
        myText(0.25, 0.65, 1, "|#eta^{#gamma}|<0.7", 0.05);
        myText(0.25, 0.60, 1, "ana509", 0.05);
        myText(0.25, 0.55, 1, pt_labels[irange].c_str(), 0.05);

        c1->SaveAs(pt_fnames[irange].c_str());

        TCanvas *c2 = new TCanvas(Form("c2_%d",irange), "", 600, 600);
        TH1D* h_frame_sideband = new TH1D(Form("h_frame_sideband_%d",irange), "h_frame_sideband", 160, -20, 10);
        h_frame_sideband->GetYaxis()->SetRangeUser(0, 1.2*ymax);
        h_frame_sideband->SetXTitle("Cluster Time - MBD Time [ns]");
        h_frame_sideband->Draw("axis");
        h_delta_t_tight_iso_cluster_proj->SetLineColor(kBlack);
        h_delta_t_tight_iso_cluster_proj->SetMarkerColor(kBlack);
        h_delta_t_tight_iso_cluster_proj->Draw("same");
        h_delta_t_tight_noniso_cluster_proj->SetLineColor(kRed);
        h_delta_t_tight_noniso_cluster_proj->SetMarkerColor(kRed);
        h_delta_t_tight_noniso_cluster_proj->Draw("same");
        h_delta_t_nontight_iso_cluster_proj->SetLineColor(kBlue);
        h_delta_t_nontight_iso_cluster_proj->SetMarkerColor(kBlue);
        h_delta_t_nontight_iso_cluster_proj->Draw("same");
        h_delta_t_nontight_noniso_cluster_proj->SetLineColor(kPink + 5);
        h_delta_t_nontight_noniso_cluster_proj->SetMarkerColor(kPink + 5);
        h_delta_t_nontight_noniso_cluster_proj->Draw("same");

        myMarkerLineText(0.25, 0.9, 1, kBlack, 20, kBlack, 1, "Tight iso clusters", 0.05, true);
        myMarkerLineText(0.25, 0.85, 1, kRed, 20, kRed, 1, "Tight noniso clusters", 0.05, true);
        myMarkerLineText(0.25, 0.8, 1, kBlue, 20, kBlue, 1, "Nontight iso clusters", 0.05, true);
        myMarkerLineText(0.25, 0.75, 1, kPink + 5, 20, kPink + 5, 1, "Nontight noniso clusters", 0.05, true);
        myText(0.25, 0.70, 1, "|#eta^{#gamma}|<0.7", 0.05);
        myText(0.25, 0.65, 1, ana_label.c_str(), 0.05);
        myText(0.25, 0.60, 1, pt_labels[irange].c_str(), 0.05);

        c2->SaveAs(pt_sideband_fnames[irange].c_str());

        // Canvas 3: B2BJet comparison for nontight iso and noniso
        TCanvas *c3 = new TCanvas(Form("c3_%d",irange), "", 600, 600);
        TH1D* h_frame_b2bjet = new TH1D(Form("h_frame_b2bjet_%d",irange), "h_frame_b2bjet", 160, -20, 10);
        double ymax_b2bjet = std::max({
            h_delta_t_nontight_iso_cluster_proj->GetMaximum(),
            h_delta_t_nontight_iso_cluster_b2bjet_proj->GetMaximum(),
            h_delta_t_nontight_noniso_cluster_proj->GetMaximum(),
            h_delta_t_nontight_noniso_cluster_b2bjet_proj->GetMaximum()
        });
        h_frame_b2bjet->GetYaxis()->SetRangeUser(0, 1.2*ymax_b2bjet);
        h_frame_b2bjet->SetXTitle("Cluster Time - MBD Time [ns]");
        h_frame_b2bjet->Draw("axis");

        // Nontight iso without b2bjet
        h_delta_t_nontight_iso_cluster_proj->SetLineColor(col_ratio[4]);
        h_delta_t_nontight_iso_cluster_proj->SetMarkerColor(col_ratio[4]);
        h_delta_t_nontight_iso_cluster_proj->SetLineStyle(1);
        h_delta_t_nontight_iso_cluster_proj->Draw("same");

        // Nontight iso with b2bjet
        h_delta_t_nontight_iso_cluster_b2bjet_proj->SetLineColor(col_ratio[0]);
        h_delta_t_nontight_iso_cluster_b2bjet_proj->SetMarkerColor(col_ratio[0]);
        h_delta_t_nontight_iso_cluster_b2bjet_proj->SetLineStyle(2);
        h_delta_t_nontight_iso_cluster_b2bjet_proj->Draw("same");

        // Nontight noniso without b2bjet
        h_delta_t_nontight_noniso_cluster_proj->SetLineColor(col_ratio[3]);
        h_delta_t_nontight_noniso_cluster_proj->SetMarkerColor(col_ratio[3]);
        h_delta_t_nontight_noniso_cluster_proj->SetLineStyle(1);
        h_delta_t_nontight_noniso_cluster_proj->Draw("same");

        // Nontight noniso with b2bjet
        h_delta_t_nontight_noniso_cluster_b2bjet_proj->SetLineColor(col_ratio[1]);
        h_delta_t_nontight_noniso_cluster_b2bjet_proj->SetMarkerColor(col_ratio[1]);
        h_delta_t_nontight_noniso_cluster_b2bjet_proj->SetLineStyle(2);
        h_delta_t_nontight_noniso_cluster_b2bjet_proj->Draw("same");

        myMarkerLineText(0.25, 0.9, 1, col_ratio[4], 20, col_ratio[4], 1, "Nontight iso", 0.05, true);
        myMarkerLineText(0.25, 0.85, 1, col_ratio[0], 20, col_ratio[0], 2, "Nontight iso w/ B2BJet", 0.05, true);
        myMarkerLineText(0.25, 0.80, 1, col_ratio[3], 20, col_ratio[3], 1, "Nontight noniso", 0.05, true);
        myMarkerLineText(0.25, 0.75, 1, col_ratio[1], 20, col_ratio[1], 2, "Nontight noniso w/ B2BJet", 0.05, true);
        myText(0.25, 0.70, 1, "|#eta^{#gamma}|<0.7", 0.05);
        myText(0.25, 0.65, 1, ana_label.c_str(), 0.05);
        myText(0.25, 0.60, 1, pt_labels[irange].c_str(), 0.05);

        c3->SaveAs(pt_b2bjet_fnames[irange].c_str());

        // Optionally, clean up histograms to avoid leaks in long runs:
        delete h_frame;
        delete h_frame_sideband;
        delete h_frame_b2bjet;
        delete h_delta_t_tight_cluster_proj;
        delete h_delta_t_tight_iso_cluster_proj;
        delete h_delta_t_tight_noniso_cluster_proj;
        delete h_delta_t_nontight_iso_cluster_proj;
        delete h_delta_t_nontight_noniso_cluster_proj;
        delete h_delta_t_common_cluster_proj;
        delete h_delta_t_all_cluster_proj;
        delete h_delta_t_tight_cluster_b2bjet_proj;
        delete h_delta_t_nontight_iso_cluster_b2bjet_proj;
        delete h_delta_t_nontight_noniso_cluster_b2bjet_proj;
        delete h_delta_t_npb_cluster_proj;
        delete c1;
        delete c2;
        delete c3;
    }

    if (h_mbd_t_vs_cluster_t)
    {
        TH1D* h_mbd_time_proj = (TH1D*)h_mbd_t_vs_cluster_t->ProjectionX("h_mbd_time_proj");
        h_mbd_time_proj->Rebin(2);
        double integral_mbd = h_mbd_time_proj->Integral();
        if (integral_mbd > 0)
        {
            h_mbd_time_proj->Scale(1.0 / integral_mbd);
        }
        double ymax_mbd = 1.2 * h_mbd_time_proj->GetMaximum();
        if (ymax_mbd <= 0)
        {
            ymax_mbd = 1.0;
        }

        TCanvas *c_mbd = new TCanvas("c_mbd", "", 600, 600);
        h_mbd_time_proj->GetXaxis()->SetRangeUser(-10, 10);
        h_mbd_time_proj->GetYaxis()->SetRangeUser(0, ymax_mbd);
        h_mbd_time_proj->SetLineColor(kBlack);
        h_mbd_time_proj->SetMarkerColor(kBlack);
        h_mbd_time_proj->SetLineWidth(2);
        h_mbd_time_proj->GetXaxis()->SetTitle("MBD Time [ns]");
        h_mbd_time_proj->GetYaxis()->SetTitle("Probability");
        h_mbd_time_proj->Draw("hist");

        myText(0.25, 0.85, 1, ana_label.c_str(), 0.05);

        c_mbd->SaveAs(mbd_time_fname.c_str());

        delete h_mbd_time_proj;
        delete c_mbd;
    }
}