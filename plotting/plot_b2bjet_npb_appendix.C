#include "plotcommon.h"

// Two-panel cross-check plot for the NPB-BDT b2b-jet validation.
// Top panel: d^2 sigma / d eta d E_T^gamma vs. E_T^gamma (log y) for nominal
//            and four b2bjet variants (pt5, pt5_npb0, pt7, pt7_npb0).
// Bottom panel: variant / nominal ratio.
// Saves to PPG12-analysis-note/Figures/analysis/b2bjet_npb_xcheck.pdf

void plot_b2bjet_npb_appendix()
{
    init_plot();

    const char *base    = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    const char *outpath = "/sphenix/user/shuhangli/ppg12/PPG12-analysis-note/Figures/analysis";
    const float deta    = 1.4;

    TFile *fnom   = TFile::Open(Form("%s/Photon_final_bdt_nom.root",             base));
    TFile *f5     = TFile::Open(Form("%s/Photon_final_bdt_b2bjet_pt5.root",      base));
    TFile *f5np   = TFile::Open(Form("%s/Photon_final_bdt_b2bjet_pt5_npb0.root", base));
    TFile *f7     = TFile::Open(Form("%s/Photon_final_bdt_b2bjet_pt7.root",      base));
    TFile *f7np   = TFile::Open(Form("%s/Photon_final_bdt_b2bjet_pt7_npb0.root", base));

    auto getXsec = [&](TFile *f, const char *cname) -> TH1F * {
        TH1F *h = (TH1F *)f->Get("h_unfold_sub_result");
        TH1F *c = (TH1F *)h->Clone(cname);
        c->Scale(1.0 / deta);
        c->SetDirectory(0);
        return c;
    };

    TH1F *hnom = getXsec(fnom, "h_xsec_nom");
    TH1F *h5   = getXsec(f5,   "h_xsec_b2b_pt5");
    TH1F *h5np = getXsec(f5np, "h_xsec_b2b_pt5_npb0");
    TH1F *h7   = getXsec(f7,   "h_xsec_b2b_pt7");
    TH1F *h7np = getXsec(f7np, "h_xsec_b2b_pt7_npb0");

    TH1F *r5   = (TH1F *)h5  ->Clone("r_b2b_pt5");        r5  ->Divide(hnom);
    TH1F *r5np = (TH1F *)h5np->Clone("r_b2b_pt5_npb0");   r5np->Divide(hnom);
    TH1F *r7   = (TH1F *)h7  ->Clone("r_b2b_pt7");        r7  ->Divide(hnom);
    TH1F *r7np = (TH1F *)h7np->Clone("r_b2b_pt7_npb0");   r7np->Divide(hnom);

    auto styleAbs = [](TH1F *h, int col, int mk) {
        h->SetLineColor(col);   h->SetMarkerColor(col);
        h->SetMarkerStyle(mk);  h->SetMarkerSize(1.3);  h->SetLineWidth(2);
    };
    const int colNom = kBlack;
    const int col5   = kAzure + 2;
    const int col7   = kRed - 4;
    styleAbs(hnom, colNom, 20);
    styleAbs(h5,   col5,   20);
    styleAbs(h5np, col5,   24);   // open marker = NPB off
    styleAbs(h7,   col7,   21);
    styleAbs(h7np, col7,   25);
    styleAbs(r5,   col5,   20);
    styleAbs(r5np, col5,   24);
    styleAbs(r7,   col7,   21);
    styleAbs(r7np, col7,   25);

    TCanvas *c = new TCanvas("c_b2b_npb_xcheck", "", 720, 720);
    c->cd();
    TPad *pTop = new TPad("pTop", "", 0.0, 0.36, 1.0, 1.0);
    TPad *pBot = new TPad("pBot", "", 0.0, 0.00, 1.0, 0.36);
    pTop->SetBottomMargin(0.02);
    pTop->SetLeftMargin(0.15);
    pTop->SetRightMargin(0.04);
    pTop->SetTopMargin(0.06);
    pBot->SetTopMargin(0.02);
    pBot->SetBottomMargin(0.30);
    pBot->SetLeftMargin(0.15);
    pBot->SetRightMargin(0.04);
    pTop->SetTickx(); pTop->SetTicky();
    pBot->SetTickx(); pBot->SetTicky();
    pTop->Draw();
    pBot->Draw();

    // Top panel
    pTop->cd();
    pTop->SetLogy();
    TH2F *frTop = new TH2F("frTop", ";;d^{2}#it{#sigma}/d#it{#eta}d#it{E}_{T}^{#gamma} [pb/GeV]",
                            100, 12, 32, 100, 1.0e-1, 5.0e3);
    frTop->GetYaxis()->SetTitleOffset(1.2);
    frTop->GetYaxis()->SetTitleSize(0.06);
    frTop->GetYaxis()->SetLabelSize(0.05);
    frTop->GetXaxis()->SetLabelSize(0.0);
    frTop->Draw();

    hnom->Draw("E1 same");
    h5  ->Draw("E1 same");
    h5np->Draw("E1 same");
    h7  ->Draw("E1 same");
    h7np->Draw("E1 same");

    TLegend *leg = new TLegend(0.50, 0.50, 0.94, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.040);
    leg->AddEntry(hnom, "Nominal selection",                          "pl");
    leg->AddEntry(h5,   "b2b jet #it{p}_{T}^{jet}#geq5 GeV",            "pl");
    leg->AddEntry(h5np, "b2b jet #it{p}_{T}^{jet}#geq5 GeV, no NPB cut", "pl");
    leg->AddEntry(h7,   "b2b jet #it{p}_{T}^{jet}#geq7 GeV",            "pl");
    leg->AddEntry(h7np, "b2b jet #it{p}_{T}^{jet}#geq7 GeV, no NPB cut", "pl");
    leg->Draw();

    myText(0.20, 0.88, 1, strleg1.c_str(),   0.052);
    myText(0.20, 0.82, 1, strleg2_1.c_str(), 0.044);
    myText(0.20, 0.77, 1, strleg3.c_str(),   0.044);

    // Bottom panel — ratio
    pBot->cd();
    TH2F *frBot = new TH2F("frBot", ";#it{E}_{T}^{#gamma} [GeV];Variant / Nominal",
                            100, 12, 32, 100, 0.4, 1.6);
    frBot->GetYaxis()->SetTitleOffset(0.65);
    frBot->GetYaxis()->SetTitleSize(0.10);
    frBot->GetXaxis()->SetTitleSize(0.10);
    frBot->GetYaxis()->SetLabelSize(0.085);
    frBot->GetXaxis()->SetLabelSize(0.085);
    frBot->GetYaxis()->SetNdivisions(505);
    frBot->Draw();

    TLine *one = new TLine(12, 1.0, 32, 1.0);
    one->SetLineStyle(2);
    one->SetLineColor(kGray + 1);
    one->Draw("same");

    r5  ->Draw("E1 same");
    r5np->Draw("E1 same");
    r7  ->Draw("E1 same");
    r7np->Draw("E1 same");

    c->cd();
    c->SaveAs(Form("%s/b2bjet_npb_xcheck.pdf", outpath));

    fnom->Close(); f5->Close(); f5np->Close(); f7->Close(); f7np->Close();
}
