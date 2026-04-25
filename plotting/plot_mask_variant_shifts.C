#include "plotcommon.h"

// σ(masked)/σ(nominal) for the 9 tower-mask variants vs mask severity, in
// three families: common-level, tight-level, OR(common, tight). σ integrated
// over pT 10-36 GeV from h_unfold_sub_result.
// Nominal is the new allz-based nom (lumi=64.3718, vertex_cut_truth=9999).

void plot_mask_variant_shifts()
{
    init_plot();
    const char *RES = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    const double pt_lo = 10.0, pt_hi = 36.0;

    auto integrate = [&](const char *fn) -> double {
        TFile *f = TFile::Open(Form("%s/%s", RES, fn), "READ");
        if (!f || f->IsZombie()) return -1;
        TH1F *h = (TH1F *) f->Get("h_unfold_sub_result");
        int b_lo = h->FindBin(pt_lo + 1e-6), b_hi = h->FindBin(pt_hi - 1e-6);
        double s = 0;
        for (int i = b_lo; i <= b_hi; ++i) s += h->GetBinContent(i) * h->GetBinWidth(i);
        f->Close();
        return s;
    };

    double sig_nom = integrate("Photon_final_bdt_nom.root");
    std::cout << "σ_nom (allz) = " << sig_nom << std::endl;

    struct Family {
        const char *label;
        int color;
        int marker;
        std::array<const char *, 3> files;
    };
    std::vector<Family> F = {
        {"common-level mask",     kAzure - 3, 20,
         {"Photon_final_bdt_mask_common_harddead.root",
          "Photon_final_bdt_mask_common_harddead_zlt5.root",
          "Photon_final_bdt_mask_common_harddead_zlt2.root"}},
        {"tight-level mask",      kOrange + 7, 21,
         {"Photon_final_bdt_mask_tight_harddead.root",
          "Photon_final_bdt_mask_tight_harddead_zlt5.root",
          "Photon_final_bdt_mask_tight_harddead_zlt2.root"}},
        {"OR(common, tight) mask", kGreen + 3, 22,
         {"Photon_final_bdt_mask_or_harddead.root",
          "Photon_final_bdt_mask_or_harddead_zlt5.root",
          "Photon_final_bdt_mask_or_harddead_zlt2.root"}},
    };

    std::vector<TGraphErrors *> gs;
    for (const auto &fam : F) {
        TGraphErrors *g = new TGraphErrors();
        for (int i = 0; i < 3; ++i) {
            double r = integrate(fam.files[i]) / sig_nom;
            g->SetPoint(i, i + 1, r);
            g->SetPointError(i, 0, 0);
            std::cout << "  " << fam.label << " sev " << (i + 1) << ": " << r << std::endl;
        }
        g->SetMarkerStyle(fam.marker);
        g->SetMarkerColor(fam.color);
        g->SetLineColor(fam.color);
        g->SetMarkerSize(1.6);
        g->SetLineWidth(2);
        gs.push_back(g);
    }

    TGraph *g_nom_ref = new TGraph();
    g_nom_ref->SetPoint(0, 0.5, 1.0);
    g_nom_ref->SetPoint(1, 3.5, 1.0);
    g_nom_ref->SetLineStyle(2); g_nom_ref->SetLineColor(kGray + 1); g_nom_ref->SetLineWidth(2);

    TCanvas *c = new TCanvas("c_mask_shift", "", 900, 650);
    c->SetLeftMargin(0.14); c->SetRightMargin(0.05);
    c->SetTopMargin(0.12);  c->SetBottomMargin(0.16);

    TH1F *frame = new TH1F("frame_mask", "", 3, 0.5, 3.5);
    frame->GetYaxis()->SetRangeUser(0.99, 1.09);
    frame->SetXTitle("mask severity (1 = hard-dead,  2 = + z<-5,  3 = + z<-2)");
    frame->SetYTitle("#sigma_{masked} / #sigma_{nominal (allz)}");
    frame->GetXaxis()->SetNdivisions(3, kFALSE);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->Draw("AXIS");

    g_nom_ref->Draw("L");
    for (auto *g : gs) g->Draw("LP same");

    TLatex lx; lx.SetNDC();
    lx.SetTextSize(0.035);
    lx.DrawLatex(0.16, 0.955, strleg1.c_str());
    lx.SetTextSize(0.028);
    lx.DrawLatex(0.16, 0.915, "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, 64.4 pb^{-1} (all-z)");

    TLegend *leg = new TLegend(0.18, 0.58, 0.55, 0.82);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.032);
    for (size_t i = 0; i < gs.size(); ++i) leg->AddEntry(gs[i], F[i].label, "LP");
    leg->AddEntry(g_nom_ref, "nominal reference", "L");
    leg->Draw();

    const char *outdir = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures";
    c->SaveAs(Form("%s/mask_variant_shifts.pdf", outdir));
    c->SaveAs(Form("%s/mask_variant_shifts.png", outdir));
}
