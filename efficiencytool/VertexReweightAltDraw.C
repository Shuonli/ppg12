// VertexReweightAltDraw.C
// Draw plots for the alternative truth-based vertex reweighting.
// Reads:  results/vertex_reweight_alt_jet.root
// Writes: reports/figures/vertex_reweight_alt_shapes.pdf
//         reports/figures/vertex_reweight_alt_weights.pdf
//         reports/figures/vertex_reweight_alt_closure.pdf
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
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <string>

namespace {

struct Variant {
    std::string key;          // single_0mrad, mixed_1p5mrad, ...
    std::string label;        // display label
    std::string crossing;     // "0 mrad" / "1.5 mrad"
    std::string mctype;       // "single" / "double" / "mixed"
};

const std::vector<Variant> kVariants = {
    {"single_0mrad",  "jet12 (single)",      "0 mrad",  "single"},
    {"double_0mrad",  "jet12 (double)",      "0 mrad",  "double"},
    {"mixed_0mrad",   "jet12 (mixed 22.4%)", "0 mrad",  "mixed"},
    {"single_1p5mrad","jet12 (single)",      "1.5 mrad","single"},
    {"double_1p5mrad","jet12 (double)",      "1.5 mrad","double"},
    {"mixed_1p5mrad", "jet12 (mixed 7.9%)",  "1.5 mrad","mixed"},
};

// Load the official sPHENIX style from cvmfs and layer on the multi-panel
// tweaks we need (wider L/R margins for 2x3 layouts, A4 landscape canvas).
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

TH1D* getHist(TFile* f, const std::string& name) {
    TH1D* h = (TH1D*)f->Get(name.c_str());
    if (!h) {
        std::cerr << "[draw] missing " << name << "\n";
        return nullptr;
    }
    h = (TH1D*)h->Clone((name+"_c").c_str());
    h->SetDirectory(nullptr);
    return h;
}

TH1D* normalize(TH1D* h) {
    if (!h) return nullptr;
    double I = h->Integral();
    if (I > 0) h->Scale(1.0 / I);
    return h;
}

void set_per_pad_margins() {
    // Wide left margin so the rotated y-axis title on the leftmost column
    // of a TCanvas::Divide grid has clearance from the canvas edge.
    gPad->SetLeftMargin(0.22);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.09);
    gPad->SetBottomMargin(0.17);
}

// TH1D axis style for subpads in a 3x2 (or 2x2) grid. The sPHENIX default
// title/label sizes are in canvas-relative NDC, so they render too large
// in divided subpads. Override with explicit sizes, and keep the y-axis
// title close to the frame so it doesn't fall outside the canvas for the
// leftmost column when TCanvas::Divide is used.
void set_subpad_axes(TH1D* h, double tsize=0.055, double lsize=0.050,
                     double title_off_y=1.30, double title_off_x=1.10) {
    h->GetXaxis()->SetTitleSize(tsize);
    h->GetXaxis()->SetLabelSize(lsize);
    h->GetXaxis()->SetTitleOffset(title_off_x);
    h->GetYaxis()->SetTitleSize(tsize);
    h->GetYaxis()->SetLabelSize(lsize);
    h->GetYaxis()->SetTitleOffset(title_off_y);
}

} // namespace

void VertexReweightAltDraw() {
    applyStyle();

    TFile* fin = TFile::Open("results/vertex_reweight_alt_jet.root", "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "[draw] cannot open results/vertex_reweight_alt_jet.root\n";
        return;
    }

    gSystem->mkdir("reports/figures", true);

    // ============== Figure 1: 2x3 shape overlay D/M_reco/T_truth ===========
    // 1500x1050 matches A4 landscape aspect (29.7/21.0 = 1.414) so nothing
    // clips when TCanvas::Print("*.pdf") maps to the paper.
    TCanvas* c1 = new TCanvas("c_shapes", "shapes", 1500, 1050);
    c1->Divide(3, 2, 0.003, 0.003);
    for (size_t i = 0; i < kVariants.size(); ++i) {
        c1->cd(i+1);
        set_per_pad_margins();

        TH1D* hD = normalize(getHist(fin, "hD_"+kVariants[i].key));
        TH1D* hM = normalize(getHist(fin, "hM_"+kVariants[i].key));
        TH1D* hT = normalize(getHist(fin, "hT_"+kVariants[i].key));
        if (!hD || !hM || !hT) continue;

        hD->SetLineColor(kBlack);      hD->SetLineWidth(3);
        hD->SetMarkerColor(kBlack);    hD->SetMarkerStyle(20); hD->SetMarkerSize(0.8);
        hM->SetLineColor(kBlue+1);     hM->SetLineWidth(3); hM->SetLineStyle(1);
        hT->SetLineColor(kRed+1);      hT->SetLineWidth(3); hT->SetLineStyle(2);

        hD->GetXaxis()->SetTitle("vertex #it{z} [cm]");
        hD->GetYaxis()->SetTitle("normalized counts");
        set_subpad_axes(hD);
        hD->GetXaxis()->SetRangeUser(-120, 120);

        double ymax = 0;
        for (auto* h : {hD, hM, hT}) if (h->GetMaximum() > ymax) ymax = h->GetMaximum();
        hD->SetMaximum(1.75 * ymax);
        hD->SetMinimum(0.0);

        hD->Draw("E1");
        hM->Draw("HIST SAME");
        hT->Draw("HIST SAME");
        hD->Draw("E1 SAME"); // bring markers to top

        // Legend upper-right, fully clear of the peak at z=0.
        TLegend* leg = new TLegend(0.54, 0.56, 0.95, 0.78);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->SetTextFont(42); leg->SetTextSize(0.044);
        leg->AddEntry(hD, "Data (#it{z}_{reco})",  "lep");
        leg->AddEntry(hM, "MC (#it{z}_{reco})",    "l");
        leg->AddEntry(hT, "MC (#it{z}_{truth})",   "l");
        leg->Draw();

        header(0.28, 0.90, kVariants[i].label, kVariants[i].crossing);
    }
    c1->SaveAs("reports/figures/vertex_reweight_alt_shapes.pdf");

    // ============== Figure 2: 2x3 w1 vs w2 vs w_iter overlay ==============
    TCanvas* c2 = new TCanvas("c_weights", "weights", 1500, 1050);
    c2->Divide(3, 2, 0.003, 0.003);
    for (size_t i = 0; i < kVariants.size(); ++i) {
        c2->cd(i+1);
        set_per_pad_margins();

        TH1D* hw1 = getHist(fin, "h_w1_"+kVariants[i].key);
        TH1D* hw2 = getHist(fin, "h_w2_"+kVariants[i].key);
        TH1D* hwi = getHist(fin, "h_w_iter_"+kVariants[i].key);
        if (!hw1 || !hw2) continue;

        hw1->SetLineColor(kBlue+1);   hw1->SetLineWidth(3); hw1->SetLineStyle(1);
        hw1->SetMarkerColor(kBlue+1); hw1->SetMarkerStyle(20); hw1->SetMarkerSize(0.8);
        hw2->SetLineColor(kRed+1);    hw2->SetLineWidth(3); hw2->SetLineStyle(2);
        hw2->SetMarkerColor(kRed+1);  hw2->SetMarkerStyle(24); hw2->SetMarkerSize(0.8);
        if (hwi) {
            hwi->SetLineColor(kGreen+2);   hwi->SetLineWidth(3); hwi->SetLineStyle(9);
            hwi->SetMarkerColor(kGreen+2); hwi->SetMarkerStyle(22); hwi->SetMarkerSize(0.8);
        }

        hw1->GetXaxis()->SetTitle("truth vertex #it{z} [cm]");
        hw1->GetYaxis()->SetTitle("weight  Data/MC");
        set_subpad_axes(hw1);
        hw1->GetXaxis()->SetRangeUser(-120, 120);

        // Per-row y-axis: 0 mrad peaks are ~1.3 (double edges up to 1.6),
        // 1.5 mrad peaks are ~3.8. Use generous headroom so the legend at
        // y-NDC 0.68-0.90 stays above the highest curve.
        const bool is_1p5 = (kVariants[i].crossing == "1.5 mrad");
        const double ymax = is_1p5 ? 6.2 : 2.5;
        hw1->SetMinimum(0.0);
        hw1->SetMaximum(ymax);
        hw1->Draw("E1");
        hw2->Draw("E1 SAME");
        if (hwi) hwi->Draw("E1 SAME");

        TLine* l1 = new TLine(-120, 1.0, 120, 1.0);
        l1->SetLineStyle(3); l1->SetLineColor(kGray+2); l1->Draw();
        TLine* la = new TLine(-60, 0.0, -60, ymax);
        TLine* lb = new TLine( 60, 0.0,  60, ymax);
        la->SetLineStyle(3); la->SetLineColor(kGray+2); la->Draw();
        lb->SetLineStyle(3); lb->SetLineColor(kGray+2); lb->Draw();

        // Legend in the upper-right strip, above the highest peak.
        // With ymax=6.0 the 1.5 mrad peak at 3.7 is at y-NDC=0.62, so the
        // legend sitting at y=[0.70, 0.92] clears it.
        TLegend* leg = new TLegend(0.52, 0.68, 0.96, 0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->SetTextFont(42); leg->SetTextSize(0.038);
        leg->AddEntry(hw1, "#it{w}_{1}: #it{D}/#it{T} direct",       "l");
        leg->AddEntry(hw2, "#it{w}_{2}: two-step via #it{f}(#it{z}_{r})", "l");
        if (hwi) leg->AddEntry(hwi, "#it{w}_{iter}: D'Agostini unfold", "l");
        leg->Draw();

        header(0.28, 0.90, kVariants[i].label, kVariants[i].crossing);
    }
    c2->SaveAs("reports/figures/vertex_reweight_alt_weights.pdf");

    // ============== Figure 3: 2x2 closure (mixed + double) ==============
    TCanvas* c3 = new TCanvas("c_closure", "closure", 1500, 1050);
    c3->Divide(2, 2, 0.003, 0.003);
    std::vector<std::string> closureKeys = {
        "mixed_0mrad", "mixed_1p5mrad",
        "double_0mrad", "double_1p5mrad"
    };
    std::vector<std::string> closureLabels = {
        "jet12 (mixed 22.4%)", "jet12 (mixed 7.9%)",
        "jet12 (double only)", "jet12 (double only)"
    };
    std::vector<std::string> closureCross = {
        "0 mrad", "1.5 mrad",
        "0 mrad", "1.5 mrad"
    };
    for (size_t i = 0; i < closureKeys.size(); ++i) {
        c3->cd(i+1);
        set_per_pad_margins();

        TH1D* hD = normalize(getHist(fin, "hD_"+closureKeys[i]));
        TH1D* hM = normalize(getHist(fin, "hM_"+closureKeys[i]));
        TH1D* hCw1 = normalize(getHist(fin, "hM_closure_w1_"+closureKeys[i]));
        TH1D* hCw2 = normalize(getHist(fin, "hM_closure_w2_"+closureKeys[i]));
        TH1D* hCit = normalize(getHist(fin, "hM_closure_iter_"+closureKeys[i]));
        if (!hD || !hM || !hCw1 || !hCw2) continue;

        hD->SetLineColor(kBlack);     hD->SetLineWidth(3);
        hD->SetMarkerColor(kBlack);   hD->SetMarkerStyle(20); hD->SetMarkerSize(0.9);
        hM->SetLineColor(kGray+2);    hM->SetLineWidth(2);
        hCw1->SetLineColor(kBlue+1);  hCw1->SetLineWidth(3); hCw1->SetLineStyle(1);
        hCw2->SetLineColor(kRed+1);   hCw2->SetLineWidth(3); hCw2->SetLineStyle(2);
        if (hCit) { hCit->SetLineColor(kGreen+2); hCit->SetLineWidth(3); hCit->SetLineStyle(9); }

        hD->GetXaxis()->SetTitle("reco vertex #it{z} [cm]");
        hD->GetYaxis()->SetTitle("normalized counts");
        set_subpad_axes(hD);
        hD->GetXaxis()->SetRangeUser(-120, 120);

        double ymax = 0;
        for (auto* h : {hD, hM, hCw1, hCw2}) if (h->GetMaximum() > ymax) ymax = h->GetMaximum();
        if (hCit && hCit->GetMaximum() > ymax) ymax = hCit->GetMaximum();
        hD->SetMaximum(1.80 * ymax);
        hD->SetMinimum(0.0);

        hD->Draw("E1");
        hM->Draw("HIST SAME");
        hCw1->Draw("HIST SAME");
        hCw2->Draw("HIST SAME");
        if (hCit) hCit->Draw("HIST SAME");
        hD->Draw("E1 SAME");

        // Legend upper-right, away from the sharp 1.5 mrad peak at z=0.
        TLegend* leg = new TLegend(0.55, 0.50, 0.96, 0.80);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->SetTextFont(42); leg->SetTextSize(0.040);
        leg->AddEntry(hD,   "Data",                                 "lep");
        leg->AddEntry(hM,   "MC raw",                               "l");
        leg->AddEntry(hCw1, "MC #times #it{w}_{1}(#it{z}_{truth})", "l");
        leg->AddEntry(hCw2, "MC #times #it{w}_{2}(#it{z}_{truth})", "l");
        if (hCit) leg->AddEntry(hCit, "MC #times #it{w}_{iter}(#it{z}_{truth})", "l");
        leg->Draw();

        header(0.28, 0.90, closureLabels[i], closureCross[i]);
    }
    c3->SaveAs("reports/figures/vertex_reweight_alt_closure.pdf");

    fin->Close();
    std::cout << "[draw] wrote:\n"
              << "  reports/figures/vertex_reweight_alt_shapes.pdf\n"
              << "  reports/figures/vertex_reweight_alt_weights.pdf\n"
              << "  reports/figures/vertex_reweight_alt_closure.pdf\n";
}
