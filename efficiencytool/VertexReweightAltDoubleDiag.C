// VertexReweightAltDoubleDiag.C
// Diagnostic plot for the double_1p5mrad closure failure.
// Shows the conditional reco kernel P(z_r | z_t=0) from the double MC
// response matrix, overlaid with the 1.5 mrad data reco shape. The
// kernel width (~34 cm) is wider than the data width (~22 cm), which
// geometrically bounds the achievable closure.
//
// Uses the official sPHENIX style loaded from cvmfs.

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include <TLine.h>
#include <TSystem.h>
#include <algorithm>
#include <iostream>
#include <string>

namespace {

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
    gPad->SetLeftMargin(0.20);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.20);
}

// Extract column j (fixed z_t bin) of TH2D as a TH1D along z_r
TH1D* slice_column(TH2D* h2, int iy, const std::string& name) {
    const int nx = h2->GetNbinsX();
    TH1D* h = new TH1D(name.c_str(), "", nx,
                       h2->GetXaxis()->GetXmin(),
                       h2->GetXaxis()->GetXmax());
    h->SetDirectory(nullptr);
    double s = 0.0;
    for (int i = 1; i <= nx; ++i) {
        const double v = h2->GetBinContent(i, iy);
        h->SetBinContent(i, v);
        s += v;
    }
    if (s > 0) h->Scale(1.0 / s);
    return h;
}
} // namespace

void VertexReweightAltDoubleDiag() {
    applyStyle();
    gSystem->mkdir("reports/figures", true);

    TFile* f = TFile::Open("results/vertex_reweight_alt_jet.root", "READ");
    if (!f || f->IsZombie()) { std::cerr << "cannot open input\n"; return; }

    // Get response matrix for double_1p5mrad and extract the z_t=0 column
    TH2D* hR = (TH2D*)f->Get("hR_double_1p5mrad");
    if (!hR) { std::cerr << "no hR_double_1p5mrad\n"; return; }
    const int iy_zero = hR->GetYaxis()->FindBin(0.0);
    TH1D* hKernel = slice_column(hR, iy_zero, "hKernel_zt0_double_1p5mrad");

    // Extract kernel at z_t = -50 and +50 too for visual comparison
    TH1D* hKernel_m50 = slice_column(hR, hR->GetYaxis()->FindBin(-50.0), "hKernel_ztm50");
    TH1D* hKernel_p50 = slice_column(hR, hR->GetYaxis()->FindBin(+50.0), "hKernel_ztp50");

    // 1.5 mrad data reco shape (normalized)
    TH1D* hD = (TH1D*)f->Get("hD_mixed_1p5mrad");
    TH1D* hD_n = (TH1D*)hD->Clone("hD_n");
    hD_n->SetDirectory(nullptr);
    if (hD_n->Integral() > 0) hD_n->Scale(1.0/hD_n->Integral());

    // Single MC reco (for comparison of kernel)
    TH1D* hM_sig = (TH1D*)f->Get("hM_single_1p5mrad");
    TH1D* hM_sig_n = (TH1D*)hM_sig->Clone("hM_sig_n");
    hM_sig_n->SetDirectory(nullptr);
    if (hM_sig_n->Integral() > 0) hM_sig_n->Scale(1.0/hM_sig_n->Integral());

    // Closure iter for double_1p5mrad
    TH1D* hC = (TH1D*)f->Get("hM_closure_iter_double_1p5mrad");
    TH1D* hC_n = (TH1D*)hC->Clone("hC_n");
    hC_n->SetDirectory(nullptr);
    if (hC_n->Integral() > 0) hC_n->Scale(1.0/hC_n->Integral());

    TCanvas* c = new TCanvas("c_diag", "diag", 1400, 1100);
    c->Divide(2, 2, 0.002, 0.002);

    // Panel 1: conditional kernel at z_t = 0 vs 1.5 mrad data
    c->cd(1);
    set_per_pad_margins();
    hKernel->SetLineColor(kBlue+1); hKernel->SetLineWidth(3); hKernel->SetLineStyle(1);
    hKernel->SetMarkerColor(kBlue+1); hKernel->SetMarkerStyle(20); hKernel->SetMarkerSize(0.8);
    hD_n->SetLineColor(kBlack);     hD_n->SetLineWidth(3);
    hD_n->SetMarkerColor(kBlack);   hD_n->SetMarkerStyle(20); hD_n->SetMarkerSize(0.9);

    hKernel->GetXaxis()->SetTitle("reco vertex z [cm]");
    hKernel->GetYaxis()->SetTitle("normalized counts");
    hKernel->GetXaxis()->SetTitleSize(0.055); hKernel->GetXaxis()->SetLabelSize(0.050);
    hKernel->GetYaxis()->SetTitleSize(0.055); hKernel->GetYaxis()->SetLabelSize(0.050);
    hKernel->GetYaxis()->SetTitleOffset(1.30);
    hKernel->GetXaxis()->SetTitleOffset(1.10);
    hKernel->GetXaxis()->SetRangeUser(-120, 120);
    hKernel->SetMinimum(0.0);
    double ymax1 = std::max(hKernel->GetMaximum(), hD_n->GetMaximum());
    hKernel->SetMaximum(1.6 * ymax1);
    hKernel->Draw("HIST");
    hD_n->Draw("E1 SAME");

    TLegend* leg1 = new TLegend(0.55, 0.55, 0.96, 0.80);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0);
    leg1->SetTextFont(42); leg1->SetTextSize(0.040);
    leg1->AddEntry(hKernel, "double MC: #it{P}(#it{z}_{r} | #it{z}_{t}=0)",  "l");
    leg1->AddEntry(hD_n,   "data mixed 1.5 mrad (#it{z}_{r})", "lep");
    leg1->Draw();
    header(0.24, 0.90, "Conditional kernel vs data", "double MC 1.5 mrad");

    // Panel 2: compare kernels at z_t = -50, 0, +50
    c->cd(2);
    set_per_pad_margins();
    hKernel->Draw("HIST");
    hKernel_m50->SetLineColor(kRed+1); hKernel_m50->SetLineWidth(3); hKernel_m50->SetLineStyle(2);
    hKernel_p50->SetLineColor(kGreen+2); hKernel_p50->SetLineWidth(3); hKernel_p50->SetLineStyle(9);
    double ymax2 = std::max({hKernel->GetMaximum(),
                             hKernel_m50->GetMaximum(),
                             hKernel_p50->GetMaximum()});
    hKernel->SetMaximum(1.6 * ymax2);
    hKernel_m50->Draw("HIST SAME");
    hKernel_p50->Draw("HIST SAME");
    TLegend* leg2 = new TLegend(0.55, 0.55, 0.96, 0.80);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0);
    leg2->SetTextFont(42); leg2->SetTextSize(0.040);
    leg2->AddEntry(hKernel,     "#it{z}_{t} = 0",   "l");
    leg2->AddEntry(hKernel_m50, "#it{z}_{t} = #minus50 cm", "l");
    leg2->AddEntry(hKernel_p50, "#it{z}_{t} = +50 cm",     "l");
    leg2->Draw();
    header(0.24, 0.90, "Kernel vs truth position", "double MC 1.5 mrad");

    // Panel 3: closure_iter vs data shows that the minimum reco width
    // after unfolding is close to the kernel width at z_t=0
    c->cd(3);
    set_per_pad_margins();
    TH1D* hDD = (TH1D*)f->Get("hD_double_1p5mrad");
    TH1D* hDD_n = (TH1D*)hDD->Clone("hDD_n");
    hDD_n->SetDirectory(nullptr);
    if (hDD_n->Integral() > 0) hDD_n->Scale(1.0/hDD_n->Integral());
    hDD_n->SetLineColor(kBlack);     hDD_n->SetLineWidth(3);
    hDD_n->SetMarkerColor(kBlack);   hDD_n->SetMarkerStyle(20); hDD_n->SetMarkerSize(0.9);
    hC_n->SetLineColor(kGreen+2);    hC_n->SetLineWidth(3); hC_n->SetLineStyle(9);
    hKernel->SetLineColor(kBlue+1);  hKernel->SetLineWidth(2); hKernel->SetLineStyle(2);

    hDD_n->GetXaxis()->SetTitle("reco vertex z [cm]");
    hDD_n->GetYaxis()->SetTitle("normalized counts");
    hDD_n->GetXaxis()->SetTitleSize(0.055); hDD_n->GetXaxis()->SetLabelSize(0.050);
    hDD_n->GetYaxis()->SetTitleSize(0.055); hDD_n->GetYaxis()->SetLabelSize(0.050);
    hDD_n->GetYaxis()->SetTitleOffset(1.30);
    hDD_n->GetXaxis()->SetTitleOffset(1.10);
    hDD_n->GetXaxis()->SetRangeUser(-120, 120);
    hDD_n->SetMinimum(0.0);
    double ymax3 = std::max(hDD_n->GetMaximum(), std::max(hC_n->GetMaximum(), hKernel->GetMaximum()));
    hDD_n->SetMaximum(1.8 * ymax3);

    hDD_n->Draw("E1");
    hC_n->Draw("HIST SAME");
    hKernel->Draw("HIST SAME");
    hDD_n->Draw("E1 SAME");

    TLegend* leg3 = new TLegend(0.55, 0.55, 0.96, 0.80);
    leg3->SetBorderSize(0); leg3->SetFillStyle(0);
    leg3->SetTextFont(42); leg3->SetTextSize(0.040);
    leg3->AddEntry(hDD_n,   "data (#it{z}_{r})",                 "lep");
    leg3->AddEntry(hC_n,   "MC #times w_{iter} (best unfold)", "l");
    leg3->AddEntry(hKernel, "kernel floor (truth #rightarrow #delta)", "l");
    leg3->Draw();
    header(0.24, 0.90, "Double MC 1.5 mrad", "unfolding bounded by kernel");

    // Panel 4: narrative -- widths summary as a TLatex table.
    // NOTE: no header() call here -- this is a text-only pad. Drawing the
    // sPHENIX header on top of a dense text stack makes the body illegible.
    // Instead we prepend two bold branding lines at the top of the stack.
    c->cd(4);
    gPad->SetLeftMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.05);
    TLatex tx;
    tx.SetNDC();

    // Branding lines (replace the missing header() call).
    tx.SetTextFont(62); tx.SetTextSize(0.06);
    tx.DrawLatex(0.05, 0.92, "#bf{#it{sPHENIX}} Internal");
    tx.SetTextFont(42); tx.SetTextSize(0.05);
    tx.DrawLatex(0.05, 0.85, "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, double MC 1.5 mrad");

    // Width-budget bullets, starting at y=0.76 with dy=0.06.
    tx.SetTextFont(42); tx.SetTextSize(0.048);
    double y = 0.76;
    const double dy = 0.06;
    tx.DrawLatex(0.05, y, "#bullet  1.5 mrad data reco width:");             tx.DrawLatex(0.78, y, "22.3");
    y -= dy;
    tx.DrawLatex(0.05, y, "#bullet  single MC truth width:");                tx.DrawLatex(0.78, y, "57.6");
    y -= dy;
    tx.DrawLatex(0.05, y, "#bullet  double MC truth width:");                tx.DrawLatex(0.78, y, "62.7");
    y -= dy;
    tx.DrawLatex(0.05, y, "#bullet  double MC reco width:");                 tx.DrawLatex(0.78, y, "48.6");
    y -= dy;
    tx.DrawLatex(0.05, y, "#bullet  conditional kernel at #it{z}_{t}=0:");   tx.DrawLatex(0.78, y, "34.1");
    y -= dy;
    tx.DrawLatex(0.05, y, "#color[4]{Gap between kernel floor and data width:}");
    tx.DrawLatex(0.78, y, "#color[4]{+11.8}");
    y -= dy;
    tx.DrawLatex(0.05, y, "#color[2]{No truth reweight can narrow double MC reco below 34 cm.}");
    y -= dy;
    tx.DrawLatex(0.05, y, "Expected double-MC reco if generation were correct:");
    y -= dy;
    tx.DrawLatex(0.05, y, "   #sigma_{data} / #sqrt{2} = 22.3/#sqrt{2} = 15.75 cm");
    y -= dy;
    tx.DrawLatex(0.05, y, "   Actual MC kernel is 3.1#times too wide.");

    c->SaveAs("reports/figures/vertex_reweight_alt_double_diag.pdf");
    std::cout << "Wrote reports/figures/vertex_reweight_alt_double_diag.pdf\n";
}
