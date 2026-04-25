// TruthVertexRecoCheckDraw.C
// Re-draw the canvases from results/truth_vertex_reco_check.root without re-processing.
//
// Uses the official sPHENIX style loaded from cvmfs.

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include <TLine.h>
#include <iostream>
#include <vector>
#include <string>

namespace {
struct Sample { std::string key, label; int cR, cN; };
const std::vector<Sample> kSamples = {
    {"photon10",        "PYTHIA #gamma+jet (single)",     kBlue+1, kRed+1},
    {"photon10_double", "PYTHIA #gamma+jet (double)",     kBlue+1, kRed+1},
    {"jet12",           "PYTHIA inclusive jet (single)",  kBlue+1, kRed+1},
    {"jet12_double",    "PYTHIA inclusive jet (double)",  kBlue+1, kRed+1},
};

// Load the official sPHENIX style from cvmfs and layer on the multi-panel
// tweaks we need (wider L/R margins for 2x2 layouts, A4 landscape canvas).
void applyStyle() {
    const char* style_path =
        "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/"
        "new.20/rootmacros/sPhenixStyle.C";
    gROOT->LoadMacro(style_path);
    gROOT->ProcessLine("SetsPhenixStyle();");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPaperSize(29.7, 21.0);
    // Tighten margins slightly for dense multi-panel grids.
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetPadBottomMargin(0.14);
}

// sPHENIX legend text, kept in sync with plotting/plotcommon.h
const char* kLegSP  = "#bf{#it{sPHENIX}} Internal";
const char* kLegPP  = "#it{p}+#it{p} #kern[-0.1]{#sqrt{#it{s}} = 200 GeV}";

// Draw a compact header block in the top-left of the current pad.
// Text is top-anchored (SetTextAlign 13) so `y` is the TOP of the first
// line -- keeps the tallest ascenders below the frame top for any
// topMargin >= 0.08.
void header(double x, double y, const std::string& lab, const std::string& sub="") {
    TLatex t; t.SetNDC();
    t.SetTextAlign(13); // left, top
    t.SetTextFont(62); t.SetTextSize(0.058);
    t.DrawLatex(x, y, kLegSP);
    t.SetTextFont(42); t.SetTextSize(0.044);
    t.DrawLatex(x, y-0.062, kLegPP);
    t.DrawLatex(x, y-0.114, lab.c_str());
    if (!sub.empty()) t.DrawLatex(x, y-0.166, sub.c_str());
}

void set_per_pad_margins() {
    gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.09);
    gPad->SetBottomMargin(0.15);
}
}

void TruthVertexRecoCheckDraw() {
    applyStyle();
    TFile* fin = TFile::Open("results/truth_vertex_reco_check.root", "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "cannot open input\n"; return; }

    // ---- Canvas 1: 2x2 shape overlay (aspect matches A4 landscape 29.7/21.0)
    TCanvas* c1 = new TCanvas("c_shape_draw", "shape", 1584, 1120);
    c1->Divide(2, 2, 0.003, 0.003);
    for (size_t i = 0; i < kSamples.size(); ++i) {
        c1->cd(i+1);
        set_per_pad_margins();

        TH1D* hR = (TH1D*)fin->Get(("h_vtxtruth_reco_"   + kSamples[i].key).c_str());
        TH1D* hN = (TH1D*)fin->Get(("h_vtxtruth_noreco_" + kSamples[i].key).c_str());
        if (!hR || !hN) { std::cerr << "missing histos for " << kSamples[i].key << "\n"; continue; }

        long long nR  = (long long)hR->Integral(0, hR->GetNbinsX()+1);
        long long nN  = (long long)hN->Integral(0, hN->GetNbinsX()+1);
        long long nTot= nR + nN;

        TH1D* hRn = (TH1D*)hR->Clone(); hRn->SetDirectory(nullptr);
        TH1D* hNn = (TH1D*)hN->Clone(); hNn->SetDirectory(nullptr);
        if (hRn->Integral() > 0) hRn->Scale(1.0 / hRn->Integral());
        if (hNn->Integral() > 0) hNn->Scale(1.0 / hNn->Integral());

        hRn->SetLineColor(kSamples[i].cR); hRn->SetLineWidth(3);
        hNn->SetLineColor(kSamples[i].cN); hNn->SetLineWidth(3); hNn->SetLineStyle(2);

        hRn->GetXaxis()->SetTitle("truth vertex #it{z} [cm]");
        hRn->GetYaxis()->SetTitle("normalized counts");
        hRn->GetYaxis()->SetTitleOffset(1.55);
        hRn->GetXaxis()->SetTitleOffset(1.10);
        double ymax = std::max(hRn->GetMaximum(), hNn->GetMaximum());
        hRn->SetMaximum(1.75 * ymax);
        hRn->SetMinimum(0.0);
        hRn->Draw("HIST");
        hNn->Draw("HIST SAME");

        TLine* l1 = new TLine(-30, 0, -30, 1.75*ymax);
        TLine* l2 = new TLine( 30, 0,  30, 1.75*ymax);
        l1->SetLineStyle(3); l1->SetLineColor(kGray+2); l1->Draw();
        l2->SetLineStyle(3); l2->SetLineColor(kGray+2); l2->Draw();

        double fR = nTot>0 ? 100.0*nR/(double)nTot : 0.0;
        double fN = nTot>0 ? 100.0*nN/(double)nTot : 0.0;
        TLegend* leg = new TLegend(0.50, 0.66, 0.96, 0.82);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->SetTextFont(42); leg->SetTextSize(0.040);
        leg->AddEntry(hRn, Form("reco vertex (%.1f%%)",    fR), "l");
        leg->AddEntry(hNn, Form("no reco vertex (%.1f%%)", fN), "l");
        leg->Draw();

        header(0.20, 0.90, kSamples[i].label);
    }
    c1->SaveAs("reports/figures/truth_vertex_reco_check.pdf");

    // ---- Canvas 2: 2x2 ratio  (noReco-shape / reco-shape)
    TCanvas* c2 = new TCanvas("c_ratio_draw", "ratio", 1584, 1120);
    c2->Divide(2, 2, 0.003, 0.003);
    for (size_t i = 0; i < kSamples.size(); ++i) {
        c2->cd(i+1);
        set_per_pad_margins();

        TH1D* hR = (TH1D*)fin->Get(("h_vtxtruth_reco_"   + kSamples[i].key).c_str());
        TH1D* hN = (TH1D*)fin->Get(("h_vtxtruth_noreco_" + kSamples[i].key).c_str());
        if (!hR || !hN) continue;
        TH1D* hRn = (TH1D*)hR->Clone(); hRn->SetDirectory(nullptr);
        TH1D* hNn = (TH1D*)hN->Clone(); hNn->SetDirectory(nullptr);
        if (hRn->Integral() > 0) hRn->Scale(1.0 / hRn->Integral());
        if (hNn->Integral() > 0) hNn->Scale(1.0 / hNn->Integral());
        TH1D* hRatio = (TH1D*)hNn->Clone(("h_ratio_" + kSamples[i].key).c_str());
        hRatio->SetDirectory(nullptr);
        hRatio->Divide(hRn);
        hRatio->SetMarkerStyle(20);
        hRatio->SetMarkerColor(kBlack);
        hRatio->SetLineColor(kBlack);
        hRatio->GetXaxis()->SetTitle("truth vertex #it{z} [cm]");
        hRatio->GetYaxis()->SetTitle("(no reco) / (reco)   shape");
        hRatio->GetYaxis()->SetTitleOffset(1.55);
        hRatio->GetXaxis()->SetTitleOffset(1.10);
        hRatio->SetMinimum(0.75);
        hRatio->SetMaximum(1.25);
        hRatio->GetXaxis()->SetRangeUser(-120, 120);
        hRatio->Draw("E1");
        TLine* l0 = new TLine(-120, 1.0, 120, 1.0);
        l0->SetLineStyle(2); l0->SetLineColor(kGray+2); l0->Draw();
        TLine* la = new TLine(-30, 0.75, -30, 1.25);
        TLine* lb = new TLine( 30, 0.75,  30, 1.25);
        la->SetLineStyle(3); la->SetLineColor(kGray+2); la->Draw();
        lb->SetLineStyle(3); lb->SetLineColor(kGray+2); lb->Draw();
        header(0.20, 0.90, kSamples[i].label);
    }
    c2->SaveAs("reports/figures/truth_vertex_reco_check_ratio.pdf");

    fin->Close();
    std::cout << "Drew reports/figures/truth_vertex_reco_check{,_ratio}.pdf\n";
}
