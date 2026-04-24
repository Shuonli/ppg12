#include "plotcommon.h"

// Cross-section ratio plot: b2bjet variants / nominal
// Variants: pt5, pt7, pt9, pt5_npb0, pt7_npb0. deta cancels in ratio.

void plot_b2bjet_over_nom()
{
    init_plot();

    const char *base = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";

    TFile *fnom   = TFile::Open(Form("%s/Photon_final_bdt_nom.root",             base));
    TFile *f5     = TFile::Open(Form("%s/Photon_final_bdt_b2bjet_pt5.root",      base));
    TFile *f7     = TFile::Open(Form("%s/Photon_final_bdt_b2bjet_pt7.root",      base));
    TFile *f9     = TFile::Open(Form("%s/Photon_final_bdt_b2bjet_pt9.root",      base));
    TFile *f5np   = TFile::Open(Form("%s/Photon_final_bdt_b2bjet_pt5_npb0.root", base));
    TFile *f7np   = TFile::Open(Form("%s/Photon_final_bdt_b2bjet_pt7_npb0.root", base));

    TH1F *hnom = (TH1F*)fnom->Get("h_unfold_sub_result");
    TH1F *h5   = (TH1F*)f5->Get("h_unfold_sub_result");
    TH1F *h7   = (TH1F*)f7->Get("h_unfold_sub_result");
    TH1F *h9   = (TH1F*)f9->Get("h_unfold_sub_result");
    TH1F *h5np = (TH1F*)f5np->Get("h_unfold_sub_result");
    TH1F *h7np = (TH1F*)f7np->Get("h_unfold_sub_result");

    TH1F *r5   = (TH1F*)h5 ->Clone("r5");   r5 ->Divide(hnom);
    TH1F *r7   = (TH1F*)h7 ->Clone("r7");   r7 ->Divide(hnom);
    TH1F *r9   = (TH1F*)h9 ->Clone("r9");   r9 ->Divide(hnom);
    TH1F *r5np = (TH1F*)h5np->Clone("r5np"); r5np->Divide(hnom);
    TH1F *r7np = (TH1F*)h7np->Clone("r7np"); r7np->Divide(hnom);

    TCanvas *c = new TCanvas("c_b2b_over_nom", "", 800, 640);
    c->cd();
    c->SetTickx();
    c->SetTicky();

    TH2F *frame = new TH2F("frame", ";#it{E}_{T}^{#gamma} [GeV];d#sigma(b2b requirement) / d#sigma(nominal)",
                           100, 10, 32, 100, 0, 1.4);
    frame->GetYaxis()->SetTitleOffset(1.2);
    frame->Draw();

    TLine *one = new TLine(10, 1.0, 32, 1.0);
    one->SetLineStyle(2);
    one->SetLineColor(kGray + 1);
    one->Draw("same");

    auto style = [](TH1F *h, int col, int mk) {
        h->SetLineColor(col); h->SetMarkerColor(col);
        h->SetMarkerStyle(mk); h->SetMarkerSize(1.2); h->SetLineWidth(2);
    };
    style(r5,   kAzure + 2, 20);
    style(r7,   kSpring - 5, 21);
    style(r9,   kRed - 4, 22);
    style(r5np, kAzure + 2, 24);   // open circle
    style(r7np, kSpring - 5, 25);  // open square

    r5  ->Draw("HIST P same");
    r7  ->Draw("HIST P same");
    r9  ->Draw("HIST P same");
    r5np->Draw("HIST P same");
    r7np->Draw("HIST P same");

    TLegend *leg = new TLegend(0.20, 0.55, 0.55, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->SetTextSize(0.032);
    leg->AddEntry(r5,   "b2bjet #it{p}_{T}#geq5 GeV", "pl");
    leg->AddEntry(r5np, "b2bjet #it{p}_{T}#geq5 GeV, no NPB", "pl");
    leg->AddEntry(r7,   "b2bjet #it{p}_{T}#geq7 GeV", "pl");
    leg->AddEntry(r7np, "b2bjet #it{p}_{T}#geq7 GeV, no NPB", "pl");
    leg->AddEntry(r9,   "b2bjet #it{p}_{T}#geq9 GeV", "pl");
    leg->Draw();

    myText(0.60, 0.30, 1, strleg1.c_str(), 0.045);
    myText(0.60, 0.25, 1, strleg2_1.c_str(), 0.035);
    myText(0.60, 0.21, 1, strleg3.c_str(), 0.035);

    c->SaveAs("/sphenix/user/shuhangli/ppg12/plotting/figures/b2bjet_over_nom.pdf");
    c->SaveAs("/sphenix/user/shuhangli/ppg12/plotting/figures/b2bjet_over_nom.png");
}
