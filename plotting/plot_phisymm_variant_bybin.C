#include "plotcommon.h"

// Per-ET-bin σ_mask/σ_nom for the 3 phi-symm mask variants (one per
// selection level). All three overlaid on a single panel.

void plot_phisymm_variant_bybin()
{
    init_plot();
    const char *RES = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    auto read = [&](const char *fn) -> TH1F * {
        TFile *f = TFile::Open(Form("%s/%s", RES, fn), "READ");
        TH1F *h = (TH1F *) f->Get("h_unfold_sub_result");
        h->SetDirectory(nullptr);
        return h;
    };
    TH1F *h_nom = read("Photon_final_bdt_nom.root");
    struct Var { const char *label; const char *file; int color; int style; };
    std::vector<Var> V = {
        {"preselect", "Photon_final_bdt_mask_phisymm_preselect.root", kAzure - 3, 20},
        {"common",    "Photon_final_bdt_mask_phisymm_common.root",    kGreen + 3, 21},
        {"tight",     "Photon_final_bdt_mask_phisymm_tight.root",     kOrange + 7, 22},
        {"OR (common #cup tight)", "Photon_final_bdt_mask_phisymm_or.root", kViolet + 1, 23},
    };
    std::vector<TH1F *> rats;
    for (const auto &v : V) {
        TH1F *h = read(v.file);
        TH1F *r = (TH1F *) h->Clone(Form("r_%s", v.label));
        r->Divide(h_nom);
        r->SetMarkerStyle(v.style);
        r->SetMarkerColor(v.color);
        r->SetLineColor(v.color);
        r->SetMarkerSize(1.5);
        r->SetLineWidth(2);
        rats.push_back(r);
    }

    TCanvas *c = new TCanvas("c_bybin_phi", "", 900, 700);
    c->SetLeftMargin(0.14); c->SetRightMargin(0.05);
    c->SetTopMargin(0.10);  c->SetBottomMargin(0.14);

    TH1F *frame = new TH1F("frame_phi", "", NptBins, ptRanges);
    frame->GetYaxis()->SetRangeUser(0.82, 1.10);
    frame->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");
    frame->SetYTitle("#sigma_{masked} / #sigma_{nominal}");
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetRangeUser(10, 36);
    frame->Draw("AXIS");

    TLine *ln = new TLine(10, 1.0, 36, 1.0);
    ln->SetLineStyle(2); ln->SetLineColor(kGray + 1); ln->SetLineWidth(2);
    ln->Draw();

    for (auto *r : rats) r->Draw("P same E1");

    TLatex lx; lx.SetNDC();
    lx.SetTextSize(0.04);
    lx.DrawLatex(0.17, 0.935, strleg1.c_str());
    lx.SetTextSize(0.032);
    lx.DrawLatex(0.17, 0.895, "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, 64.4 pb^{-1} (all-z)");
    lx.DrawLatex(0.17, 0.86,  "phi-symmetry data-driven mask, z<-2");

    TLegend *leg = new TLegend(0.60, 0.18, 0.93, 0.38);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.035);
    for (size_t i = 0; i < V.size(); ++i)
        leg->AddEntry(rats[i], Form("%s mask", V[i].label), "LP");
    leg->Draw();

    c->SaveAs("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures/phisymm_variant_bybin.pdf");
    c->SaveAs("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures/phisymm_variant_bybin.png");
}
