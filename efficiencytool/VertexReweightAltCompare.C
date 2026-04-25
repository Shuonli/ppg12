// VertexReweightAltCompare.C
// Overlay jet12 and photon10 alternative weights + closure.
// Inputs:
//   results/vertex_reweight_alt_jet.root
//   results/vertex_reweight_alt_photon.root
// Outputs:
//   reports/figures/vertex_reweight_alt_jet_vs_photon_weights.pdf
//   reports/figures/vertex_reweight_alt_jet_vs_photon_closure.pdf

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

struct Panel {
    std::string crossing;   // "0 mrad" / "1.5 mrad"
    std::string variant;    // "mixed_0mrad" / "mixed_1p5mrad"
};
const std::vector<Panel> kPanels = {
    {"0 mrad",   "mixed_0mrad"},
    {"1.5 mrad", "mixed_1p5mrad"},
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
    gPad->SetTopMargin(0.09);
    gPad->SetBottomMargin(0.18);
}

// Axis style helper for subpads in a 2x2 Divide grid.
void set_subpad_axes(TH1D* h, double tsize=0.055, double lsize=0.050,
                     double title_off_y=1.30, double title_off_x=1.10) {
    h->GetXaxis()->SetTitleSize(tsize);
    h->GetXaxis()->SetLabelSize(lsize);
    h->GetXaxis()->SetTitleOffset(title_off_x);
    h->GetYaxis()->SetTitleSize(tsize);
    h->GetYaxis()->SetLabelSize(lsize);
    h->GetYaxis()->SetTitleOffset(title_off_y);
}
TH1D* getH(TFile* f, const std::string& name) {
    TH1D* h = (TH1D*)f->Get(name.c_str());
    if (!h) { std::cerr << "missing " << name << "\n"; return nullptr; }
    h = (TH1D*)h->Clone((name+"_c").c_str());
    h->SetDirectory(nullptr);
    return h;
}
TH1D* normed(TH1D* h) {
    if (!h) return nullptr;
    double I = h->Integral();
    if (I > 0) h->Scale(1.0 / I);
    return h;
}
} // namespace

void VertexReweightAltCompare() {
    applyStyle();
    gSystem->mkdir("reports/figures", true);

    TFile* fJ = TFile::Open("results/vertex_reweight_alt_jet.root", "READ");
    TFile* fP = TFile::Open("results/vertex_reweight_alt_photon.root", "READ");
    if (!fJ || fJ->IsZombie() || !fP || fP->IsZombie()) {
        std::cerr << "[compare] missing input file(s)\n"; return;
    }

    // ---------- Figure A: weights overlay (2 columns = 0,1.5 mrad; 2 rows = w1, w2)
    TCanvas* c1 = new TCanvas("c_cmp_w", "cmp_w", 1400, 1100);
    c1->Divide(2, 2, 0.002, 0.002);
    for (size_t col = 0; col < kPanels.size(); ++col) {
        const auto& p = kPanels[col];
        for (int row = 0; row < 2; ++row) {
            // row 0 = w1; row 1 = w2
            int pad = row*2 + col + 1;
            c1->cd(pad);
            set_per_pad_margins();

            const char* wname = (row == 0) ? "h_w1_" : "h_w2_";
            TH1D* hJ = getH(fJ, std::string(wname) + p.variant);
            TH1D* hP = getH(fP, std::string(wname) + p.variant);
            if (!hJ || !hP) continue;

            hJ->SetLineColor(kBlue+1);   hJ->SetLineWidth(3);
            hJ->SetMarkerColor(kBlue+1); hJ->SetMarkerStyle(20); hJ->SetMarkerSize(0.7);
            hP->SetLineColor(kRed+1);    hP->SetLineWidth(3); hP->SetLineStyle(2);
            hP->SetMarkerColor(kRed+1);  hP->SetMarkerStyle(24); hP->SetMarkerSize(0.7);

            hJ->GetXaxis()->SetTitle("truth vertex #it{z} [cm]");
            hJ->GetYaxis()->SetTitle((std::string("weight ") + (row == 0 ? "#it{w}_{1}" : "#it{w}_{2}")).c_str());
            set_subpad_axes(hJ);
            hJ->GetXaxis()->SetRangeUser(-120, 120);
            double ymax = (p.variant == "mixed_1p5mrad") ? 4.5 : 2.0;
            hJ->SetMinimum(0.0);
            hJ->SetMaximum(ymax);
            hJ->Draw("E1");
            hP->Draw("E1 SAME");

            TLine* l0 = new TLine(-120, 1.0, 120, 1.0);
            l0->SetLineStyle(3); l0->SetLineColor(kGray+2); l0->Draw();

            TLegend* leg = new TLegend(0.58, 0.54, 0.96, 0.70);
            leg->SetBorderSize(0); leg->SetFillStyle(0);
            leg->SetTextFont(42); leg->SetTextSize(0.044);
            leg->AddEntry(hJ, "jet12",    "l");
            leg->AddEntry(hP, "photon10", "l");
            leg->Draw();
            std::string sub = (row == 0 ? "mixed, #it{w}_{1}" : "mixed, #it{w}_{2}");
            header(0.26, 0.90, p.crossing, sub);
        }
    }
    c1->SaveAs("reports/figures/vertex_reweight_alt_jet_vs_photon_weights.pdf");

    // ---------- Figure B: closure overlay per sample (2x2 = sample x crossing)
    TCanvas* c2 = new TCanvas("c_cmp_cl", "cmp_cl", 1400, 1100);
    c2->Divide(2, 2, 0.002, 0.002);
    struct CPanel { const char* label; const char* sub; TFile* f; const char* key; };
    CPanel cps[4] = {
        {"jet12 mixed",    "0 mrad",   fJ, "mixed_0mrad"},
        {"jet12 mixed",    "1.5 mrad", fJ, "mixed_1p5mrad"},
        {"photon10 mixed", "0 mrad",   fP, "mixed_0mrad"},
        {"photon10 mixed", "1.5 mrad", fP, "mixed_1p5mrad"},
    };
    for (int i = 0; i < 4; ++i) {
        c2->cd(i+1);
        set_per_pad_margins();

        TH1D* hD  = normed(getH(cps[i].f, std::string("hD_")  + cps[i].key));
        TH1D* hM  = normed(getH(cps[i].f, std::string("hM_")  + cps[i].key));
        TH1D* hC1 = normed(getH(cps[i].f, std::string("hM_closure_w1_") + cps[i].key));
        TH1D* hC2 = normed(getH(cps[i].f, std::string("hM_closure_w2_") + cps[i].key));
        if (!hD || !hM || !hC1 || !hC2) continue;

        hD->SetLineColor(kBlack);   hD->SetLineWidth(3);
        hD->SetMarkerColor(kBlack); hD->SetMarkerStyle(20); hD->SetMarkerSize(0.8);
        hM->SetLineColor(kGray+2);  hM->SetLineWidth(2);
        hC1->SetLineColor(kBlue+1); hC1->SetLineWidth(3);
        hC2->SetLineColor(kRed+1);  hC2->SetLineWidth(3); hC2->SetLineStyle(2);

        hD->GetXaxis()->SetTitle("reco vertex #it{z} [cm]");
        hD->GetYaxis()->SetTitle("normalized counts");
        set_subpad_axes(hD);
        hD->GetXaxis()->SetRangeUser(-120, 120);
        double ymax = 0;
        for (auto* h : {hD, hM, hC1, hC2}) if (h->GetMaximum() > ymax) ymax = h->GetMaximum();
        hD->SetMaximum(1.80 * ymax);
        hD->SetMinimum(0.0);

        hD->Draw("E1");
        hM->Draw("HIST SAME");
        hC1->Draw("HIST SAME");
        hC2->Draw("HIST SAME");
        hD->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.56, 0.42, 0.96, 0.70);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->SetTextFont(42); leg->SetTextSize(0.040);
        leg->AddEntry(hD,   "Data",                                    "lep");
        leg->AddEntry(hM,   "MC raw",                                  "l");
        leg->AddEntry(hC1, "MC #times #it{w}_{1}(#it{z}_{truth})",     "l");
        leg->AddEntry(hC2, "MC #times #it{w}_{2}(#it{z}_{truth})",     "l");
        leg->Draw();

        header(0.26, 0.90, cps[i].label, cps[i].sub);
    }
    c2->SaveAs("reports/figures/vertex_reweight_alt_jet_vs_photon_closure.pdf");

    fJ->Close(); fP->Close();
    std::cout << "[compare] wrote:\n"
              << "  reports/figures/vertex_reweight_alt_jet_vs_photon_weights.pdf\n"
              << "  reports/figures/vertex_reweight_alt_jet_vs_photon_closure.pdf\n";
}
