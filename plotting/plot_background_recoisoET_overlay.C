#include "plotcommon.h"
#include <yaml-cpp/yaml.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TSystem.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

void plot_background_recoisoET_overlay(const std::string &configname = "config_bdt_nom.yaml",
                                       const std::string &infile = "",
                                       bool normalize = true,
                                       int rebin_iso = 40)
{
    init_plot();

    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node config = YAML::LoadFile(("../efficiencytool/" + configname).c_str());

    const std::vector<double> pT_bin_edges = config["analysis"]["pT_bins"].as<std::vector<double>>();
    const int n_pT_bins = static_cast<int>(pT_bin_edges.size()) - 1;

    const std::vector<double> eta_bins = config["analysis"]["eta_bins"].as<std::vector<double>>();
    const int n_eta_bins = static_cast<int>(eta_bins.size()) - 1;

    const std::string eff_outfile = config["output"]["eff_outfile"].as<std::string>();
    const std::string var_type = config["output"]["var_type"].as<std::string>();
    const std::string input_file = infile.empty() ? (eff_outfile + "_jet_" + var_type + ".root") : infile;

    TFile *fin = TFile::Open(input_file.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "Error: cannot open input file: " << input_file << std::endl;
        return;
    }

    const std::string savePath = "figures/background_recoisoET_overlay";
    gSystem->Exec(Form("mkdir -p %s", savePath.c_str()));

    for (int ieta = 0; ieta < n_eta_bins; ++ieta)
    {
        TH2D *h2_tight = dynamic_cast<TH2D *>(fin->Get(Form("h_tight_recoisoET_background_%d", ieta)));
        TH2D *h2_nontight = dynamic_cast<TH2D *>(fin->Get(Form("h_nontight_recoisoET_background_%d", ieta)));

        if (!h2_tight || !h2_nontight)
        {
            std::cerr << "Warning: missing background recoisoET TH2Ds for eta bin " << ieta << std::endl;
            continue;
        }

        for (int ipt = 0; ipt < n_pT_bins; ++ipt)
        {
            const double pt_lo = pT_bin_edges[ipt];
            const double pt_hi = pT_bin_edges[ipt + 1];

            const int xbin_lo = h2_tight->GetXaxis()->FindBin(pt_lo + 1e-6);
            const int xbin_hi = h2_tight->GetXaxis()->FindBin(pt_hi - 1e-6);

            TH1D *h_tight = h2_tight->ProjectionY(
                Form("h_tight_bkg_recoiso_eta%d_pt%d", ieta, ipt), xbin_lo, xbin_hi);
            TH1D *h_nontight = h2_nontight->ProjectionY(
                Form("h_nontight_bkg_recoiso_eta%d_pt%d", ieta, ipt), xbin_lo, xbin_hi);
            h_tight->SetDirectory(nullptr);
            h_nontight->SetDirectory(nullptr);

            if (rebin_iso > 1)
            {
                h_tight->Rebin(rebin_iso);
                h_nontight->Rebin(rebin_iso);
            }

            const double int_tight = h_tight->Integral("width");
            const double int_nontight = h_nontight->Integral("width");
            if (int_tight <= 0 && int_nontight <= 0)
            {
                delete h_tight;
                delete h_nontight;
                continue;
            }

            if (normalize)
            {
                if (int_tight > 0) h_tight->Scale(1.0 / int_tight);
                if (int_nontight > 0) h_nontight->Scale(1.0 / int_nontight);
            }

            h_tight->SetLineColor(kBlue + 1);
            h_tight->SetLineWidth(2);
            h_tight->SetMarkerColor(kBlue + 1);
            h_tight->SetMarkerStyle(20);
            h_tight->SetMarkerSize(0.8);

            h_nontight->SetLineColor(kRed + 1);
            h_nontight->SetLineWidth(2);
            h_nontight->SetMarkerColor(kRed + 1);
            h_nontight->SetMarkerStyle(24);
            h_nontight->SetMarkerSize(0.8);

            const double ymax = 1.35 * std::max(h_tight->GetMaximum(), h_nontight->GetMaximum());

            TCanvas *c = new TCanvas(Form("c_bkg_recoiso_eta%d_pt%d", ieta, ipt), "", 650, 600);
            c->SetLeftMargin(0.15);
            c->SetBottomMargin(0.13);
            c->SetTopMargin(0.07);
            c->SetRightMargin(0.05);

            TH1F *fr = dynamic_cast<TH1F *>(frame_isoET->Clone(Form("fr_bkg_recoiso_eta%d_pt%d", ieta, ipt)));
            fr->GetXaxis()->SetRangeUser(-5, 30);
            fr->GetYaxis()->SetRangeUser(0, std::max(1e-8, ymax));
            fr->SetYTitle(normalize ? "Probability density / GeV" : "Counts");
            fr->Draw("axis");

            h_tight->Draw("same hist");
            h_nontight->Draw("same hist");
            h_tight->Draw("same e0");
            h_nontight->Draw("same e0");

            myText(0.18, 0.90, 1, strleg1.c_str(), 0.042);
            myText(0.18, 0.85, 1, strleg2.c_str(), 0.042);
            myText(0.18, 0.80, 1, strIncMC.c_str(), 0.040);
            myText(0.18, 0.75, 1, Form("%.1f < #it{#eta} < %.1f", eta_bins[ieta], eta_bins[ieta + 1]), 0.040);
            myText(0.18, 0.70, 1, Form("%.0f < #it{E}_{T}^{#gamma,rec} < %.0f GeV", pt_lo, pt_hi), 0.040);

            TLegend *leg = new TLegend(0.55, 0.75, 0.90, 0.90);
            legStyle(leg, 0.17, 0.038);
            leg->AddEntry(h_tight, "Background tight", "lep");
            leg->AddEntry(h_nontight, "Background non-tight", "lep");
            leg->Draw();

            c->SaveAs(Form("%s/bkg_recoiso_overlay_eta%d_pt%d.pdf", savePath.c_str(), ieta, ipt));

            delete fr;
            delete leg;
            delete c;
            delete h_tight;
            delete h_nontight;
        }
    }

    fin->Close();
    std::cout << "Saved overlays to " << savePath << std::endl;
}
