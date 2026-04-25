// Per-period reco vertex-z data vs MC, fiducial |z|<60.
// If MC matches data per period, the truth-vertex reweight is fit correctly
// and any merged-plot mismatch (plot_vertexz_allz.C) is purely from merge-
// weight construction, not reweight-fit quality.

#include "plotcommon.h"

namespace {
TH1* load(const char* fpath, const char* alias) {
    TFile* f = TFile::Open(fpath, "READ");
    if (!f || f->IsZombie()) { std::cerr << "[ERR] " << fpath << "\n"; return nullptr; }
    TH1* h = (TH1*)f->Get("h_vertexz");
    if (!h) { std::cerr << "[ERR] no h_vertexz in " << fpath << "\n"; return nullptr; }
    TH1* hc = (TH1*)h->Clone(alias);
    hc->SetDirectory(nullptr);
    f->Close();
    return hc;
}
}

void plot_vertexz_per_period()
{
    init_plot();

    struct Per { const char* period; double lumi; TH1* data; TH1* mc; };
    std::vector<Per> entries = {
        { "0rad",    32.6574, nullptr, nullptr },
        { "1p5mrad", 16.2735, nullptr, nullptr },
    };

    for (auto& e : entries) {
        TString dpath = Form("/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom_%s.root", e.period);
        TString mpath = Form("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom_%s.root", e.period);
        e.data = load(dpath, Form("h_dt_%s", e.period));
        e.mc   = load(mpath, Form("h_mc_%s", e.period));
        if (!e.data || !e.mc) return;
        if (e.data->Integral() > 0) e.data->Scale(1.0 / e.data->Integral());
        if (e.mc->Integral() > 0)   e.mc->Scale(1.0 / e.mc->Integral());
    }

    TCanvas* c = new TCanvas("c_vtx_per_period", "", 1200, 500);
    c->Divide(2, 1);

    for (size_t ip = 0; ip < entries.size(); ++ip) {
        auto& e = entries[ip];
        c->cd(ip + 1);
        gPad->SetLeftMargin(0.13);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.08);
        gPad->SetBottomMargin(0.13);

        double ymax = std::max(e.data->GetMaximum(), e.mc->GetMaximum());

        TH1F* fr = new TH1F(Form("fr_%s", e.period), "", 100, -100, 100);
        fr->GetXaxis()->SetTitle("reco vertex z [cm]");
        fr->GetYaxis()->SetTitle("fraction / 1 cm");
        fr->GetXaxis()->SetTitleSize(0.050);
        fr->GetYaxis()->SetTitleSize(0.050);
        fr->GetXaxis()->SetLabelSize(0.045);
        fr->GetYaxis()->SetLabelSize(0.045);
        fr->GetYaxis()->SetTitleOffset(1.25);
        fr->SetMinimum(0);
        fr->SetMaximum(ymax * 1.35);
        fr->Draw("axis");

        // Fiducial boundary
        TLine* lp = new TLine(60, 0, 60, ymax * 1.1);   lp->SetLineColor(kGray + 1); lp->SetLineStyle(2); lp->Draw();
        TLine* lm = new TLine(-60, 0, -60, ymax * 1.1); lm->SetLineColor(kGray + 1); lm->SetLineStyle(2); lm->Draw();

        e.data->SetLineColor(kBlack);
        e.data->SetMarkerColor(kBlack);
        e.data->SetMarkerStyle(20);
        e.data->SetMarkerSize(0.9);
        e.data->SetLineWidth(1);
        e.mc->SetLineColor(kAzure + 2);
        e.mc->SetLineWidth(2);

        e.mc->Draw("HIST same");
        e.data->Draw("E1 X0 same");

        TLegend* leg = new TLegend(0.16, 0.72, 0.56, 0.89);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.045);
        leg->AddEntry((TObject*)nullptr, strleg1.c_str(), "");
        leg->AddEntry((TObject*)nullptr, Form("%s, %.2f pb^{-1}",
                      (strcmp(e.period, "0rad") == 0) ? "0 mrad" : "1.5 mrad", e.lumi), "");
        leg->AddEntry(e.data, "data", "pl");
        leg->AddEntry(e.mc,   "MC (reweighted)", "l");
        leg->Draw();

        // Print per-bin residual
        std::cout << "\n=== " << e.period << " (lumi " << e.lumi << " pb^{-1}) ===\n";
        std::cout << "z-bin   data     MC      data/MC\n";
        int n = e.data->GetNbinsX();
        for (int i = 1; i <= n; i += 10) {
            double z = e.data->GetBinCenter(i);
            double d = e.data->GetBinContent(i);
            double m = e.mc->GetBinContent(i);
            std::cout << Form("%+5.0f   %.4f   %.4f   %.3f\n", z, d, m, (m>0 ? d/m : 0));
        }
    }

    c->SaveAs("/sphenix/user/shuhangli/ppg12/plotting/figures/vertexz_per_period.pdf");
}
