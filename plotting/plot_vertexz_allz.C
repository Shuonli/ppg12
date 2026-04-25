// Quick diagnostic: reco vertex-z distribution, data vs MC, for the allz
// cross-check. Overlays area-normalized shapes so the reweight quality +
// fiducial vs beam-delivered effects are visible at a glance.

#include "plotcommon.h"

void plot_vertexz_allz()
{
    init_plot();

    auto load = [](const char* fpath, const char* label) {
        TFile* f = TFile::Open(fpath, "READ");
        if (!f || f->IsZombie()) { std::cerr << "[ERR] cannot open " << fpath << "\n"; return (TH1*)nullptr; }
        TH1* h = (TH1*)f->Get("h_vertexz");
        if (!h) { std::cerr << "[ERR] no h_vertexz in " << fpath << "\n"; return (TH1*)nullptr; }
        TH1* hc = (TH1*)h->Clone(Form("h_vtx_%s", label));
        hc->SetDirectory(nullptr);
        f->Close();
        return hc;
    };

    TH1* h_mc_allz = load("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_allz.root", "mc_allz");
    TH1* h_mc_nom  = load("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root",  "mc_nom");
    TH1* h_dt_allz = load("/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_allz.root",    "dt_allz");
    TH1* h_dt_nom  = load("/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom.root",     "dt_nom");
    if (!h_mc_allz || !h_mc_nom || !h_dt_allz || !h_dt_nom) return;

    // Area-normalize to 1 for shape comparison
    for (TH1* h : {h_mc_allz, h_mc_nom, h_dt_allz, h_dt_nom}) {
        if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
    }

    TCanvas* c = new TCanvas("c_vtx_allz", "", 800, 600);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.08);
    c->SetBottomMargin(0.12);

    double ymax = 0;
    for (TH1* h : {h_mc_allz, h_mc_nom, h_dt_allz, h_dt_nom}) ymax = std::max(ymax, h->GetMaximum());

    TH1F* fr = new TH1F("fr_vtx", "", 100, -100, 100);
    fr->GetXaxis()->SetTitle("reco vertex z [cm]");
    fr->GetYaxis()->SetTitle("fraction / 1 cm");
    fr->GetXaxis()->SetTitleSize(0.05);
    fr->GetYaxis()->SetTitleSize(0.05);
    fr->GetXaxis()->SetLabelSize(0.045);
    fr->GetYaxis()->SetLabelSize(0.045);
    fr->GetYaxis()->SetTitleOffset(1.30);
    fr->SetMinimum(0);
    fr->SetMaximum(ymax * 1.25);
    fr->Draw("axis");

    auto style = [](TH1* h, int c, int m, int ls) {
        h->SetLineColor(c);
        h->SetMarkerColor(c);
        h->SetMarkerStyle(m);
        h->SetMarkerSize(0.9);
        h->SetLineWidth(2);
        h->SetLineStyle(ls);
    };
    style(h_dt_nom,  kBlack,      20, 1);
    style(h_dt_allz, kGray + 2,   24, 1);
    style(h_mc_nom,  kAzure + 2,  21, 1);
    style(h_mc_allz, kOrange + 7, 25, 2);

    h_dt_nom->Draw("E1 X0 same");
    h_dt_allz->Draw("E1 X0 same");
    h_mc_nom->Draw("HIST same");
    h_mc_allz->Draw("HIST same");

    // Fiducial |z|<60 reference lines
    TLine* lp = new TLine(60, 0, 60, ymax * 1.05);   lp->SetLineColor(kGray + 1); lp->SetLineStyle(2); lp->Draw();
    TLine* lm = new TLine(-60, 0, -60, ymax * 1.05); lm->SetLineColor(kGray + 1); lm->SetLineStyle(2); lm->Draw();

    TLegend* leg = new TLegend(0.54, 0.58, 0.96, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.033);
    leg->AddEntry((TObject*)nullptr, strleg1.c_str(), "");
    leg->AddEntry((TObject*)nullptr, "#it{p}+#it{p} #sqrt{#it{s}}=200 GeV, run 47289-54000", "");
    leg->AddEntry(h_dt_nom,  "data (nom, |z_{reco}|<60)",   "pl");
    leg->AddEntry(h_dt_allz, "data (allz, |z_{reco}|<60)",  "pl");
    leg->AddEntry(h_mc_nom,  "MC nom: w_{0}/w_{1}=2.01 (60cm)", "l");
    leg->AddEntry(h_mc_allz, "MC allz: w_{0}/w_{1}=2.75 (allz)", "l");
    leg->Draw();

    c->SaveAs("/sphenix/user/shuhangli/ppg12/plotting/figures/vertexz_allz_vs_data.pdf");

    // Also print a summary
    std::cout << "\n=== normalized h_vertexz integrals in |z|<60 ===\n";
    int b_m60 = h_mc_nom->FindBin(-60.0);
    int b_p60 = h_mc_nom->FindBin( 60.0) - 1;
    for (auto [h, lab] : std::vector<std::pair<TH1*, const char*>>{
         {h_dt_nom, "data nom"}, {h_dt_allz, "data allz"},
         {h_mc_nom, "MC nom"}, {h_mc_allz, "MC allz"} }) {
        std::cout << "  " << lab << ": in-60 fraction = " << h->Integral(b_m60, b_p60) << "\n";
    }
}
