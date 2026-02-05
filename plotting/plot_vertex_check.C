#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGraph.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <string>

#include "plotcommon.h"

// Compute efficiency vs threshold: fraction of events with score > threshold
TGraph* computeEfficiency(TH1D* h, const std::string& name) {
    if (!h || h->Integral() == 0) return nullptr;
    int nbins = h->GetNbinsX();
    double total = h->Integral(1, nbins);
    TGraph* g = new TGraph();
    g->SetName(name.c_str());
    for (int i = 1; i <= nbins; i++) {
        g->SetPoint(g->GetN(), h->GetBinLowEdge(i), h->Integral(i, nbins) / total);
    }
    g->SetPoint(g->GetN(), h->GetBinLowEdge(nbins + 1), 0.0);
    return g;
}

// Style helper for original/shifted graphs
void styleGraphPair(TGraph* g_orig, TGraph* g_shift) {
    g_orig->SetLineColor(kBlue);   g_orig->SetLineWidth(2);
    g_shift->SetLineColor(kRed);   g_shift->SetLineWidth(2); g_shift->SetLineStyle(2);
}

// Draw legend for original/shifted
TLegend* drawLegend(TGraph* g_orig, TGraph* g_shift, double x1=0.55, double y1=0.70) {
    TLegend* leg = new TLegend(x1, y1, x1+0.33, y1+0.18);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(g_orig, "Original Vertex", "l");
    leg->AddEntry(g_shift, "Shifted Vertex", "l");
    leg->Draw();
    return leg;
}

void plot_vertex_check(const std::string &inputfile_signal = "../efficiencytool/results/MC_vertex_check_photon.root",
                       const std::string &inputfile_bkg = "../efficiencytool/results/MC_vertex_check_jet.root",
                       const std::string &savePath = "../PPG12-analysis-note/Figures/")
{
    init_plot();
    gSystem->Exec(Form("mkdir -p %s", savePath.c_str()));

    TFile *fin_sig = TFile::Open(inputfile_signal.c_str(), "READ");
    TFile *fin_bkg = TFile::Open(inputfile_bkg.c_str(), "READ");
    if (!fin_sig || fin_sig->IsZombie()) { std::cerr << "Cannot open signal file\n"; return; }
    bool hasBkg = fin_bkg && !fin_bkg->IsZombie();

    std::vector<double> pT_bins = {10, 14, 18, 22, 26, 35, 50};
    int nPtBins = pT_bins.size() - 1;

    auto saveCanvas = [&](TCanvas* c, const std::string& name) {
        c->SaveAs((savePath + name + ".pdf").c_str());
    };

    // Config: NPB uses inclusive MC, BDT uses signal MC
    struct PlotConfig { const char* type; TFile* file; const char* mcLabel; };
    std::vector<PlotConfig> configs;
    if (hasBkg) configs.push_back({"npb", fin_bkg, "Inclusive MC"});
    configs.push_back({"bdt", fin_sig, "Signal MC"});

    // =========================================================================
    // 2D Correlation plots - pT inclusive
    // =========================================================================
    for (const auto& cfg : configs) {
        TH2D* h2d = (TH2D*)cfg.file->Get(Form("h_%s_2d_eta0_inclusive", cfg.type));
        if (!h2d) continue;

        TCanvas* c = new TCanvas(Form("c_%s_2d", cfg.type), "", 800, 700);
        c->SetRightMargin(0.15);
        h2d->SetTitle("");
        h2d->GetXaxis()->SetTitle(Form("Original %s Score", cfg.type == std::string("npb") ? "NPB" : "BDT"));
        h2d->GetYaxis()->SetTitle(Form("Shifted %s Score", cfg.type == std::string("npb") ? "NPB" : "BDT"));
        h2d->Draw("COLZ");

        TLine* diag = new TLine(0, 0, 1, 1);
        diag->SetLineColor(kBlack); diag->SetLineStyle(2); diag->SetLineWidth(2);
        diag->Draw("same");

        TLatex latex; latex.SetNDC(); latex.SetTextSize(0.04);
        latex.DrawLatex(0.18, 0.92, strleg1.c_str());
        latex.DrawLatex(0.18, 0.87, (std::string(cfg.mcLabel) + " (PYTHIA8)").c_str());
        latex.DrawLatex(0.18, 0.82, strleg3.c_str());

        saveCanvas(c, std::string(cfg.type) + "_score_2d_correlation");
        delete c;
    }

    // =========================================================================
    // 1D Overlay plots - pT inclusive
    // =========================================================================
    for (const auto& cfg : configs) {
        TH2D* h2d = (TH2D*)cfg.file->Get(Form("h_%s_2d_eta0_inclusive", cfg.type));
        if (!h2d) continue;

        TCanvas* c = new TCanvas(Form("c_%s_1d", cfg.type), "", 800, 600);
        c->SetLogy();

        TH1D* h_orig = h2d->ProjectionX(Form("h_%s_orig", cfg.type));
        TH1D* h_shift = h2d->ProjectionY(Form("h_%s_shift", cfg.type));
        if (h_orig->Integral() > 0) h_orig->Scale(1.0 / h_orig->Integral());
        if (h_shift->Integral() > 0) h_shift->Scale(1.0 / h_shift->Integral());

        h_orig->SetTitle("");
        h_orig->GetXaxis()->SetTitle(Form("%s Score", cfg.type == std::string("npb") ? "NPB" : "BDT"));
        h_orig->GetYaxis()->SetTitle("Normalized Counts");
        h_orig->SetLineColor(kBlue); h_orig->SetLineWidth(2);
        h_orig->SetMaximum(h_orig->GetMaximum() * 5); h_orig->SetMinimum(1e-4);
        h_shift->SetLineColor(kRed); h_shift->SetLineWidth(2); h_shift->SetLineStyle(2);

        h_orig->Draw("HIST"); h_shift->Draw("HIST same");

        TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->AddEntry(h_orig, "Original Vertex", "l");
        leg->AddEntry(h_shift, "Shifted Vertex", "l");
        leg->Draw();

        TLatex latex; latex.SetNDC(); latex.SetTextSize(0.04);
        latex.DrawLatex(0.18, 0.92, strleg1.c_str());
        latex.DrawLatex(0.18, 0.87, (std::string(cfg.mcLabel) + " (PYTHIA8)").c_str());
        latex.DrawLatex(0.18, 0.82, strleg3.c_str());

        saveCanvas(c, std::string(cfg.type) + "_score_1d_overlay");
        delete c;
    }

    // =========================================================================
    // 2D per pT bin - individual plots
    // =========================================================================
    for (const auto& cfg : configs) {
        for (int ipt = 0; ipt < nPtBins; ipt++) {
            TH2D* h = (TH2D*)cfg.file->Get(Form("h_%s_2d_eta0_pt%d", cfg.type, ipt));
            if (!h) continue;

            TCanvas* c = new TCanvas(Form("c_%s_2d_pt%d", cfg.type, ipt), "", 800, 700);
            c->SetRightMargin(0.15);
            h->SetTitle("");
            h->GetXaxis()->SetTitle(Form("Original %s Score", cfg.type == std::string("npb") ? "NPB" : "BDT"));
            h->GetYaxis()->SetTitle(Form("Shifted %s Score", cfg.type == std::string("npb") ? "NPB" : "BDT"));
            h->Draw("COLZ");

            TLine* diag = new TLine(0, 0, 1, 1);
            diag->SetLineColor(kBlack); diag->SetLineStyle(2); diag->SetLineWidth(2);
            diag->Draw("same");

            TLatex latex; latex.SetNDC(); latex.SetTextSize(0.04);
            latex.DrawLatex(0.18, 0.92, strleg1.c_str());
            latex.DrawLatex(0.18, 0.87, (std::string(cfg.mcLabel) + " (PYTHIA8)").c_str());
            latex.DrawLatex(0.18, 0.82, Form("%.0f < p_{T} < %.0f GeV", pT_bins[ipt], pT_bins[ipt + 1]));

            saveCanvas(c, Form("%s_score_2d_pt%d", cfg.type, ipt));
            delete c;
        }
    }

    // =========================================================================
    // 1D per pT bin - individual plots
    // =========================================================================
    for (const auto& cfg : configs) {
        for (int ipt = 0; ipt < nPtBins; ipt++) {
            TH2D* h2d = (TH2D*)cfg.file->Get(Form("h_%s_2d_eta0_pt%d", cfg.type, ipt));
            if (!h2d) continue;

            TCanvas* c = new TCanvas(Form("c_%s_1d_pt%d", cfg.type, ipt), "", 800, 600);
            c->SetLogy();

            TH1D* h_orig = h2d->ProjectionX(Form("h_%s_orig_pt%d", cfg.type, ipt));
            TH1D* h_shift = h2d->ProjectionY(Form("h_%s_shift_pt%d", cfg.type, ipt));
            if (h_orig->Integral() > 0) h_orig->Scale(1.0 / h_orig->Integral());
            if (h_shift->Integral() > 0) h_shift->Scale(1.0 / h_shift->Integral());

            h_orig->SetTitle("");
            h_orig->GetXaxis()->SetTitle(Form("%s Score", cfg.type == std::string("npb") ? "NPB" : "BDT"));
            h_orig->GetYaxis()->SetTitle("Normalized Counts");
            h_orig->SetLineColor(kBlue); h_orig->SetLineWidth(2);
            h_orig->SetMaximum(h_orig->GetMaximum() * 5); h_orig->SetMinimum(1e-4);
            h_shift->SetLineColor(kRed); h_shift->SetLineWidth(2); h_shift->SetLineStyle(2);

            h_orig->Draw("HIST"); h_shift->Draw("HIST same");

            TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
            leg->SetBorderSize(0); leg->SetFillStyle(0);
            leg->AddEntry(h_orig, "Original Vertex", "l");
            leg->AddEntry(h_shift, "Shifted Vertex", "l");
            leg->Draw();

            TLatex latex; latex.SetNDC(); latex.SetTextSize(0.04);
            latex.DrawLatex(0.18, 0.92, strleg1.c_str());
            latex.DrawLatex(0.18, 0.87, (std::string(cfg.mcLabel) + " (PYTHIA8)").c_str());
            latex.DrawLatex(0.18, 0.82, Form("%.0f < p_{T} < %.0f GeV", pT_bins[ipt], pT_bins[ipt + 1]));

            saveCanvas(c, Form("%s_score_1d_pt%d", cfg.type, ipt));
            delete c;
        }
    }

    // =========================================================================
    // Efficiency vs Threshold - pT inclusive
    // =========================================================================
    for (const auto& cfg : configs) {
        TH2D* h2d = (TH2D*)cfg.file->Get(Form("h_%s_2d_eta0_inclusive", cfg.type));
        if (!h2d) continue;

        TCanvas* c = new TCanvas(Form("c_%s_eff", cfg.type), "", 800, 600);
        TH1D* h_orig = h2d->ProjectionX(Form("h_%s_eff_orig", cfg.type));
        TH1D* h_shift = h2d->ProjectionY(Form("h_%s_eff_shift", cfg.type));

        TGraph* g_orig = computeEfficiency(h_orig, Form("g_%s_eff_orig", cfg.type));
        TGraph* g_shift = computeEfficiency(h_shift, Form("g_%s_eff_shift", cfg.type));
        if (!g_orig || !g_shift) { delete c; continue; }

        styleGraphPair(g_orig, g_shift);

        TH1F* frame = new TH1F(Form("frame_%s_eff", cfg.type), "", 100, 0, 1);
        frame->SetTitle("");
        frame->GetXaxis()->SetTitle(Form("%s Score Threshold", cfg.type == std::string("npb") ? "NPB" : "BDT"));
        frame->GetYaxis()->SetTitle("Selection Efficiency (score > threshold)");
        frame->GetYaxis()->SetRangeUser(0, 1.05);
        frame->Draw();

        g_orig->Draw("L same"); g_shift->Draw("L same");
        drawLegend(g_orig, g_shift);

        TLatex latex; latex.SetNDC(); latex.SetTextSize(0.04);
        latex.DrawLatex(0.18, 0.92, strleg1.c_str());
        latex.DrawLatex(0.18, 0.87, (std::string(cfg.mcLabel) + " (PYTHIA8)").c_str());
        latex.DrawLatex(0.18, 0.82, strleg3.c_str());

        saveCanvas(c, std::string(cfg.type) + "_efficiency_vs_threshold");
        delete c;
    }

    // =========================================================================
    // Efficiency per pT bin - individual plots
    // =========================================================================
    for (const auto& cfg : configs) {
        for (int ipt = 0; ipt < nPtBins; ipt++) {
            TH2D* h2d = (TH2D*)cfg.file->Get(Form("h_%s_2d_eta0_pt%d", cfg.type, ipt));
            if (!h2d) continue;

            TCanvas* c = new TCanvas(Form("c_%s_eff_pt%d", cfg.type, ipt), "", 800, 600);

            TH1D* h_orig = h2d->ProjectionX(Form("h_%s_eff_orig_pt%d", cfg.type, ipt));
            TH1D* h_shift = h2d->ProjectionY(Form("h_%s_eff_shift_pt%d", cfg.type, ipt));
            TGraph* g_orig = computeEfficiency(h_orig, Form("g_%s_eff_orig_pt%d", cfg.type, ipt));
            TGraph* g_shift = computeEfficiency(h_shift, Form("g_%s_eff_shift_pt%d", cfg.type, ipt));
            if (!g_orig || !g_shift) { delete c; continue; }

            styleGraphPair(g_orig, g_shift);

            TH1F* frame = new TH1F(Form("frame_%s_eff_pt%d", cfg.type, ipt), "", 100, 0, 1);
            frame->GetXaxis()->SetTitle(Form("%s Score Threshold", cfg.type == std::string("npb") ? "NPB" : "BDT"));
            frame->GetYaxis()->SetTitle("Selection Efficiency (score > threshold)");
            frame->GetYaxis()->SetRangeUser(0, 1.05);
            frame->Draw();

            g_orig->Draw("L same"); g_shift->Draw("L same");
            drawLegend(g_orig, g_shift);

            TLatex latex; latex.SetNDC(); latex.SetTextSize(0.04);
            latex.DrawLatex(0.18, 0.92, strleg1.c_str());
            latex.DrawLatex(0.18, 0.87, (std::string(cfg.mcLabel) + " (PYTHIA8)").c_str());
            latex.DrawLatex(0.18, 0.82, Form("%.0f < p_{T} < %.0f GeV", pT_bins[ipt], pT_bins[ipt + 1]));

            saveCanvas(c, Form("%s_efficiency_pt%d", cfg.type, ipt));
            delete c;
        }
    }

    fin_sig->Close();
    if (hasBkg) fin_bkg->Close();
    std::cout << "Done. Plots saved to " << savePath << std::endl;
}
