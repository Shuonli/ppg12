// plot_xsec_period_ratio.C
//
// Cross-check: overlay the cross-sections for the 60cm fiducial per-period
// pipelines (bdt_0rad, bdt_nom = 1.5 mrad, bdt_all = 0mrad + 1.5mrad merged)
// and show the ratios σ_period/σ_all vs ET.
//
// Physics expectation: if the lumi-weighted merge is self-consistent, the
// all-range cross-section should equal the lumi-weighted average of the two
// per-period results. At the per-pT-bin level σ_0mrad and σ_1.5mrad can
// differ from σ_all by statistical fluctuations (each covers a fraction of
// the total lumi), but the two ratios should bracket 1 with per-period error
// bars consistent with statistics.
//
// Usage: root -l -b -q 'plot_xsec_period_ratio.C'
// Output: plotting/figures/xsec_period_ratio.pdf

#include "plotcommon.h"

namespace {

TH1F *load_xsec(const char *tune, const char *alias)
{
    TString path = Form("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_%s.root", tune);
    TFile *f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "[ERROR] cannot open " << path << std::endl;
        return nullptr;
    }
    TH1F *h = (TH1F *)f->Get("h_unfold_sub_result");
    if (!h) {
        std::cerr << "[ERROR] h_unfold_sub_result missing in " << path << std::endl;
        return nullptr;
    }
    TH1F *hc = (TH1F *)h->Clone(Form("h_%s", alias));
    hc->SetDirectory(nullptr);
    hc->Scale(1.0 / 1.4);  // divide by delta eta = 1.4 to match plot_final_selection.C
    return hc;
}

} // namespace

void plot_xsec_period_ratio()
{
    init_plot();

    TH1F *h_0rad = load_xsec("bdt_0rad", "0rad");
    TH1F *h_nom  = load_xsec("bdt_nom",  "nom");
    TH1F *h_all  = load_xsec("bdt_all",  "all");
    if (!h_0rad || !h_nom || !h_all) return;

    TH1F *h_r_0rad = (TH1F *)h_0rad->Clone("h_r_0rad");
    h_r_0rad->Divide(h_all);
    TH1F *h_r_nom = (TH1F *)h_nom->Clone("h_r_nom");
    h_r_nom->Divide(h_all);
    // Direct period-to-period ratio σ_0mrad / σ_1.5mrad (doesn't involve σ_all)
    TH1F *h_r_0rad_over_nom = (TH1F *)h_0rad->Clone("h_r_0rad_over_nom");
    h_r_0rad_over_nom->Divide(h_nom);

    // Colors: 0mrad = red-ish, 1.5mrad = blue-ish, all = black, 0/1.5 = green
    const int col_0rad     = kOrange + 7;
    const int col_nom      = kAzure + 2;
    const int col_all      = kBlack;
    const int col_0_over_1 = kGreen + 3;

    auto style = [](TH1 *h, int c, int m) {
        h->SetLineColor(c);
        h->SetMarkerColor(c);
        h->SetMarkerStyle(m);
        h->SetMarkerSize(1.2);
        h->SetLineWidth(2);
    };
    style(h_0rad, col_0rad, 21);
    style(h_nom,  col_nom,  20);
    style(h_all,  col_all,  33);
    style(h_r_0rad, col_0rad, 21);
    style(h_r_nom,  col_nom,  20);
    style(h_r_0rad_over_nom, col_0_over_1, 22);

    TCanvas *c1 = new TCanvas("c1", "xsec period ratio", 600, 700);
    c1->cd();

    TPad *pad_up = new TPad("pad_up", "", 0.0, 0.36, 1.0, 1.0);
    pad_up->SetBottomMargin(0.02);
    pad_up->SetLeftMargin(0.16);
    pad_up->SetRightMargin(0.04);
    pad_up->SetTopMargin(0.08);
    pad_up->SetLogy();
    pad_up->Draw();

    TPad *pad_dn = new TPad("pad_dn", "", 0.0, 0.0, 1.0, 0.36);
    pad_dn->SetTopMargin(0.02);
    pad_dn->SetBottomMargin(0.34);
    pad_dn->SetLeftMargin(0.16);
    pad_dn->SetRightMargin(0.04);
    pad_dn->Draw();

    // ---- Top: overlay cross-sections ----
    pad_up->cd();
    TH1F *fr = (TH1F *)frame_et_rec->Clone("fr_top");
    fr->GetYaxis()->SetTitle("d#sigma/d#it{p}_{T} d#eta [pb/GeV]");
    fr->GetYaxis()->SetRangeUser(0.2, 5e4);
    fr->GetXaxis()->SetLabelSize(0);
    fr->GetXaxis()->SetTitleSize(0);
    fr->GetYaxis()->SetLabelSize(0.050);
    fr->GetYaxis()->SetTitleSize(0.060);
    fr->GetYaxis()->SetTitleOffset(1.15);
    fr->Draw();

    h_all->Draw("E1 X0 same");
    h_0rad->Draw("E1 X0 same");
    h_nom->Draw("E1 X0 same");

    TLegend *leg = new TLegend(0.50, 0.64, 0.92, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.042);
    leg->AddEntry((TObject *)nullptr, strleg1.c_str(), "");
    leg->AddEntry((TObject *)nullptr, "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV", "");
    leg->AddEntry(h_all,  "All runs (48.9 pb^{-1})",   "pl");
    leg->AddEntry(h_0rad, "0 mrad (32.7 pb^{-1})",     "pl");
    leg->AddEntry(h_nom,  "1.5 mrad (16.3 pb^{-1})",   "pl");
    leg->Draw();

    // ---- Bottom: ratio σ_period / σ_all ----
    pad_dn->cd();
    TH1F *fr_r = (TH1F *)frame_et_rec->Clone("fr_bot");
    fr_r->GetYaxis()->SetTitle("ratio");
    fr_r->GetYaxis()->SetRangeUser(0.3, 2.4);
    fr_r->GetYaxis()->SetNdivisions(505);
    fr_r->GetXaxis()->SetTitle("#it{p}_{T}^{#gamma} [GeV]");
    fr_r->GetXaxis()->SetLabelSize(0.089);
    fr_r->GetXaxis()->SetTitleSize(0.106);
    fr_r->GetXaxis()->SetTitleOffset(1.35);
    fr_r->GetYaxis()->SetLabelSize(0.089);
    fr_r->GetYaxis()->SetTitleSize(0.106);
    fr_r->GetYaxis()->SetTitleOffset(0.68);
    fr_r->Draw();

    TLine *l1 = new TLine(8, 1.0, 36, 1.0);
    l1->SetLineColor(kGray + 1);
    l1->SetLineStyle(2);
    l1->Draw();

    h_r_0rad->Draw("E1 X0 same");
    h_r_nom->Draw("E1 X0 same");
    h_r_0rad_over_nom->Draw("E1 X0 same");

    // 3-entry legend along the top of the ratio pad; y-range above (2.0-2.4)
    // keeps data points (max ~1.7 with error bars to ~2.0) clear of the legend.
    TLegend *leg_r = new TLegend(0.18, 0.86, 0.98, 0.98);
    leg_r->SetBorderSize(0);
    leg_r->SetFillStyle(0);
    leg_r->SetTextSize(0.065);
    leg_r->SetNColumns(3);
    leg_r->AddEntry(h_r_0rad,          "0mrad/all",      "pl");
    leg_r->AddEntry(h_r_nom,           "1.5mrad/all",    "pl");
    leg_r->AddEntry(h_r_0rad_over_nom, "0mrad/1.5mrad",  "pl");
    leg_r->Draw();

    c1->SaveAs("/sphenix/user/shuhangli/ppg12/plotting/figures/xsec_period_ratio.pdf");

    // Print numerics for the per-bin ratios
    std::cout << "\n=== per-bin ratios ===\n";
    std::cout << "pT bin   pT_ctr    σ_0rad/σ_all       σ_1.5mrad/σ_all     σ_0rad/σ_1.5mrad\n";
    for (int i = 1; i <= h_r_0rad->GetNbinsX(); ++i) {
        double pt = h_r_0rad->GetBinCenter(i);
        double r0 = h_r_0rad->GetBinContent(i);
        double e0 = h_r_0rad->GetBinError(i);
        double r1 = h_r_nom->GetBinContent(i);
        double e1 = h_r_nom->GetBinError(i);
        double r01 = h_r_0rad_over_nom->GetBinContent(i);
        double e01 = h_r_0rad_over_nom->GetBinError(i);
        std::cout << Form(" %2d    %6.2f   %.3f ± %.3f     %.3f ± %.3f     %.3f ± %.3f\n",
                          i, pt, r0, e0, r1, e1, r01, e01);
    }
}
