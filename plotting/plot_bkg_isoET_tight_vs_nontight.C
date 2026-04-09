// Background-only isoET: tight vs non-tight BDT (truth-tagged, no signal)
// Uses h_tight_recoisoET_background_0 / h_nontight_recoisoET_background_0 (TH2D)
// projected per pT bin. Pure background — signal photons removed by truth matching.
//
// Usage: root -l -b -q 'plot_bkg_isoET_tight_vs_nontight.C("config_bdt_nom.yaml")'

#include "plotcommon.h"
#include <yaml-cpp/yaml.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TSystem.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

void plot_bkg_isoET_tight_vs_nontight(const std::string &configname = "config_bdt_nom.yaml",
                                       int rebin_iso = 40)
{
    init_plot();

    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node config = YAML::LoadFile(("../efficiencytool/" + configname).c_str());

    const std::vector<double> pT_bin_edges = config["analysis"]["pT_bins"].as<std::vector<double>>();
    const int n_pT_bins = static_cast<int>(pT_bin_edges.size()) - 1;

    const std::string eff_outfile = config["output"]["eff_outfile"].as<std::string>();
    const std::string var_type = config["output"]["var_type"].as<std::string>();
    const std::string jet_file = eff_outfile + "_jet_" + var_type + ".root";

    TFile *fin = TFile::Open(jet_file.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "Error: cannot open " << jet_file << std::endl;
        return;
    }

    gSystem->Exec("mkdir -p figures");

    TH2D *h2_tight = dynamic_cast<TH2D *>(fin->Get("h_tight_recoisoET_background_0"));
    TH2D *h2_nontight = dynamic_cast<TH2D *>(fin->Get("h_nontight_recoisoET_background_0"));

    if (!h2_tight || !h2_nontight)
    {
        std::cerr << "Error: missing background recoisoET TH2Ds" << std::endl;
        return;
    }

    // ── Multi-panel 4x3 ──
    const int ncols = 3, nrows = 4;
    TCanvas *c1 = new TCanvas("c_bkg_iso", "", 550 * ncols, 450 * nrows);
    c1->Divide(ncols, nrows, 0.001, 0.001);

    for (int ipt = 0; ipt < n_pT_bins; ++ipt)
    {
        const double pt_lo = pT_bin_edges[ipt];
        const double pt_hi = pT_bin_edges[ipt + 1];

        const int xbin_lo = h2_tight->GetXaxis()->FindBin(pt_lo + 1e-6);
        const int xbin_hi = h2_tight->GetXaxis()->FindBin(pt_hi - 1e-6);

        TH1D *ht = h2_tight->ProjectionY(Form("ht_bkg_%d", ipt), xbin_lo, xbin_hi);
        TH1D *hnt = h2_nontight->ProjectionY(Form("hnt_bkg_%d", ipt), xbin_lo, xbin_hi);
        ht->SetDirectory(nullptr);
        hnt->SetDirectory(nullptr);

        if (rebin_iso > 1)
        {
            ht->Rebin(rebin_iso);
            hnt->Rebin(rebin_iso);
        }

        const double int_t = ht->Integral("width");
        const double int_nt = hnt->Integral("width");
        if (int_t > 0) ht->Scale(1.0 / int_t);
        if (int_nt > 0) hnt->Scale(1.0 / int_nt);

        ht->SetLineColor(kBlue + 1);
        ht->SetLineWidth(2);
        ht->SetMarkerColor(kBlue + 1);
        ht->SetMarkerStyle(20);
        ht->SetMarkerSize(0.7);

        hnt->SetLineColor(kRed + 1);
        hnt->SetLineWidth(2);
        hnt->SetMarkerColor(kRed + 1);
        hnt->SetMarkerStyle(24);
        hnt->SetMarkerSize(0.7);

        TH1D *hratio = dynamic_cast<TH1D *>(ht->Clone(Form("hr_bkg_%d", ipt)));
        hratio->Divide(hnt);
        hratio->SetLineColor(kBlack);
        hratio->SetMarkerColor(kBlack);
        hratio->SetMarkerStyle(20);
        hratio->SetMarkerSize(0.5);

        c1->cd(ipt + 1);
        TPad *pm = new TPad(Form("pm_%d", ipt), "", 0, 0.3, 1, 1);
        pm->SetLeftMargin(0.15);
        pm->SetBottomMargin(0.005);
        pm->SetTopMargin(0.08);
        pm->SetRightMargin(0.05);
        pm->Draw();

        TPad *pr = new TPad(Form("pr_%d", ipt), "", 0, 0, 1, 0.3);
        pr->SetLeftMargin(0.15);
        pr->SetBottomMargin(0.30);
        pr->SetTopMargin(0.02);
        pr->SetRightMargin(0.05);
        pr->Draw();

        pm->cd();
        const double ymax = 1.5 * std::max(ht->GetMaximum(), hnt->GetMaximum());

        TH1F *fr = dynamic_cast<TH1F *>(frame_isoET->Clone(Form("fr_bkg_%d", ipt)));
        fr->GetYaxis()->SetRangeUser(0, std::max(1e-8, ymax));
        fr->SetYTitle("Prob. density / GeV");
        fr->GetYaxis()->SetTitleSize(0.055);
        fr->GetYaxis()->SetLabelSize(0.050);
        fr->GetYaxis()->SetTitleOffset(1.3);
        fr->GetXaxis()->SetRangeUser(-3, 15);
        fr->Draw("axis");

        ht->Draw("same hist");
        hnt->Draw("same hist");
        ht->Draw("same e0");
        hnt->Draw("same e0");

        myText(0.19, 0.88, 1, Form("%.0f < #it{p}_{T} < %.0f GeV", pt_lo, pt_hi), 0.055);
        myText(0.19, 0.78, 1, strleg1.c_str(), 0.050);
        myText(0.19, 0.70, 1, "PYTHIA8 Jet MC (bkg only)", 0.042);

        TLegend *leg = new TLegend(0.48, 0.65, 0.92, 0.88);
        legStyle(leg, 0.17, 0.045);
        leg->AddEntry(ht, "Background tight", "lep");
        leg->AddEntry(hnt, "Background non-tight", "lep");
        leg->Draw();

        pr->cd();
        TH1F *fr_r = new TH1F(Form("frr_bkg_%d", ipt), "", 100, -3, 15);
        fr_r->SetXTitle("#it{E}_{T}^{iso} [GeV]");
        fr_r->SetYTitle("T / NT");
        fr_r->GetXaxis()->SetTitleSize(0.11);
        fr_r->GetXaxis()->SetLabelSize(0.10);
        fr_r->GetXaxis()->SetTitleOffset(1.0);
        fr_r->GetYaxis()->SetTitleSize(0.10);
        fr_r->GetYaxis()->SetLabelSize(0.09);
        fr_r->GetYaxis()->SetTitleOffset(0.55);
        fr_r->GetYaxis()->SetNdivisions(505);
        fr_r->GetYaxis()->SetRangeUser(0, 3.5);
        fr_r->GetXaxis()->SetRangeUser(-3, 15);
        fr_r->Draw("axis");

        TLine *l1 = new TLine(-3, 1, 15, 1);
        l1->SetLineColor(kGray + 1);
        l1->SetLineStyle(7);
        l1->Draw();

        hratio->Draw("same e0");
    }

    c1->SaveAs(Form("figures/bkg_isoET_tight_vs_nontight_jet_%s.pdf", var_type.c_str()));

    // ── Individual per-pT canvases ──
    for (int ipt = 0; ipt < n_pT_bins; ++ipt)
    {
        const double pt_lo = pT_bin_edges[ipt];
        const double pt_hi = pT_bin_edges[ipt + 1];

        const int xbin_lo = h2_tight->GetXaxis()->FindBin(pt_lo + 1e-6);
        const int xbin_hi = h2_tight->GetXaxis()->FindBin(pt_hi - 1e-6);

        TH1D *ht = h2_tight->ProjectionY(Form("ht2_bkg_%d", ipt), xbin_lo, xbin_hi);
        TH1D *hnt = h2_nontight->ProjectionY(Form("hnt2_bkg_%d", ipt), xbin_lo, xbin_hi);
        ht->SetDirectory(nullptr);
        hnt->SetDirectory(nullptr);

        if (rebin_iso > 1) { ht->Rebin(rebin_iso); hnt->Rebin(rebin_iso); }

        const double int_t = ht->Integral("width");
        const double int_nt = hnt->Integral("width");
        if (int_t <= 0 && int_nt <= 0) continue;
        if (int_t > 0) ht->Scale(1.0 / int_t);
        if (int_nt > 0) hnt->Scale(1.0 / int_nt);

        ht->SetLineColor(kBlue + 1); ht->SetLineWidth(2);
        ht->SetMarkerColor(kBlue + 1); ht->SetMarkerStyle(20); ht->SetMarkerSize(0.8);
        hnt->SetLineColor(kRed + 1); hnt->SetLineWidth(2);
        hnt->SetMarkerColor(kRed + 1); hnt->SetMarkerStyle(24); hnt->SetMarkerSize(0.8);

        TH1D *hratio = dynamic_cast<TH1D *>(ht->Clone(Form("hr2_bkg_%d", ipt)));
        hratio->Divide(hnt);
        hratio->SetLineColor(kBlack); hratio->SetMarkerColor(kBlack);
        hratio->SetMarkerStyle(20); hratio->SetMarkerSize(0.6);

        TCanvas *c2 = new TCanvas(Form("c_bkg_%d", ipt), "", 650, 750);
        TPad *p1 = new TPad("p1", "", 0, 0.35, 1, 1);
        p1->SetLeftMargin(0.15); p1->SetBottomMargin(0.005);
        p1->SetTopMargin(0.07); p1->SetRightMargin(0.05);
        p1->Draw();
        TPad *p2 = new TPad("p2", "", 0, 0, 1, 0.35);
        p2->SetLeftMargin(0.15); p2->SetBottomMargin(0.25);
        p2->SetTopMargin(0.02); p2->SetRightMargin(0.05);
        p2->Draw();

        p1->cd();
        const double ymax = 1.45 * std::max(ht->GetMaximum(), hnt->GetMaximum());
        TH1F *fr = dynamic_cast<TH1F *>(frame_isoET->Clone(Form("fr2b_%d", ipt)));
        fr->GetXaxis()->SetRangeUser(-3, 15);
        fr->GetYaxis()->SetRangeUser(0, std::max(1e-8, ymax));
        fr->SetYTitle("Probability density / GeV");
        fr->Draw("axis");

        ht->Draw("same hist"); hnt->Draw("same hist");
        ht->Draw("same e0"); hnt->Draw("same e0");

        myText(0.18, 0.90, 1, strleg1.c_str(), 0.045);
        myText(0.18, 0.84, 1, strleg2.c_str(), 0.042);
        myText(0.18, 0.78, 1, "PYTHIA8 Jet MC (bkg only)", 0.038);
        myText(0.18, 0.72, 1, strleg3.c_str(), 0.040);
        myText(0.18, 0.66, 1, Form("%.0f < #it{E}_{T}^{#gamma,rec} < %.0f GeV", pt_lo, pt_hi), 0.040);

        TLegend *leg = new TLegend(0.50, 0.72, 0.92, 0.90);
        legStyle(leg, 0.17, 0.040);
        leg->AddEntry(ht, "Background tight", "lep");
        leg->AddEntry(hnt, "Background non-tight", "lep");
        leg->Draw();

        p2->cd();
        TH1F *fr_r = new TH1F(Form("frr2b_%d", ipt), "", 100, -3, 15);
        fr_r->SetXTitle("#it{E}_{T}^{iso} [GeV]");
        fr_r->SetYTitle("Tight / Non-tight");
        fr_r->GetXaxis()->SetTitleSize(0.09); fr_r->GetXaxis()->SetLabelSize(0.08);
        fr_r->GetXaxis()->SetTitleOffset(1.0);
        fr_r->GetYaxis()->SetTitleSize(0.08); fr_r->GetYaxis()->SetLabelSize(0.07);
        fr_r->GetYaxis()->SetTitleOffset(0.65); fr_r->GetYaxis()->SetNdivisions(505);
        fr_r->GetYaxis()->SetRangeUser(0, 3.5);
        fr_r->GetXaxis()->SetRangeUser(-3, 15);
        fr_r->Draw("axis");

        TLine *l1 = new TLine(-3, 1, 15, 1);
        l1->SetLineColor(kGray + 1); l1->SetLineStyle(7); l1->Draw();
        hratio->Draw("same e0");

        c2->SaveAs(Form("figures/bkg_isoET_tight_vs_nontight_jet_%s_pt%d.pdf",
                         var_type.c_str(), ipt));
        delete c2;
    }

    fin->Close();
    std::cout << "Saved to figures/bkg_isoET_tight_vs_nontight_jet_" << var_type << "*.pdf" << std::endl;
}
