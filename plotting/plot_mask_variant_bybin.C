#include "plotcommon.h"

// Per-pT-bin σ_mask/σ_nom for the 9 tower-mask variants. Three panels,
// one per mask family (common / tight / OR), each showing the 3 severities
// vs unfolded ET.

void plot_mask_variant_bybin()
{
    init_plot();
    const char *RES = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results";

    auto read_spec = [&](const char *fn) -> TH1F * {
        TFile *f = TFile::Open(Form("%s/%s", RES, fn), "READ");
        if (!f || f->IsZombie()) return nullptr;
        TH1F *h = (TH1F *) f->Get("h_unfold_sub_result");
        h->SetDirectory(nullptr);
        return h;
    };

    TH1F *h_nom = read_spec("Photon_final_bdt_nom.root");

    struct Var { const char *label; const char *file; int color; int style; };
    struct Family { const char *name; std::vector<Var> vars; };

    std::vector<Family> F = {
        {"common-level mask", {
            {"hard-dead only", "Photon_final_bdt_mask_common_harddead.root",       kAzure - 3, 20},
            {"+ z<-5",         "Photon_final_bdt_mask_common_harddead_zlt5.root",  kTeal + 2,  21},
            {"+ z<-2",         "Photon_final_bdt_mask_common_harddead_zlt2.root",  kBlue + 2,  22},
        }},
        {"tight-level mask", {
            {"hard-dead only", "Photon_final_bdt_mask_tight_harddead.root",        kOrange + 7, 20},
            {"+ z<-5",         "Photon_final_bdt_mask_tight_harddead_zlt5.root",   kRed - 3,    21},
            {"+ z<-2",         "Photon_final_bdt_mask_tight_harddead_zlt2.root",   kRed + 2,    22},
        }},
        {"OR(common, tight) mask", {
            {"hard-dead only", "Photon_final_bdt_mask_or_harddead.root",           kGreen + 1,  20},
            {"+ z<-5",         "Photon_final_bdt_mask_or_harddead_zlt5.root",      kGreen + 3,  21},
            {"+ z<-2",         "Photon_final_bdt_mask_or_harddead_zlt2.root",      kSpring - 7, 22},
        }},
    };

    auto make_ratio = [&](const Var &v) -> TH1F * {
        TH1F *h = read_spec(v.file);
        TH1F *r = (TH1F *) h->Clone(Form("ratio_%s", v.file));
        r->Divide(h_nom);
        r->SetMarkerStyle(v.style);
        r->SetMarkerColor(v.color);
        r->SetLineColor(v.color);
        r->SetMarkerSize(1.3);
        r->SetLineWidth(2);
        return r;
    };

    TCanvas *c = new TCanvas("c_bybin", "", 900, 1300);
    c->Divide(1, 3, 0, 0.01);

    auto draw_panel = [&](TVirtualPad *p, const Family &fam,
                          int ipanel, int npanel)
    {
        bool is_top    = (ipanel == 1);
        bool is_bottom = (ipanel == npanel);
        p->SetLeftMargin(0.14);
        p->SetRightMargin(0.06);
        p->SetTopMargin(is_top ? 0.12 : 0.02);
        p->SetBottomMargin(is_bottom ? 0.18 : 0.02);
        p->cd();

        TH1F *frame = new TH1F(Form("frame_p%d", ipanel), "", NptBins, ptRanges);
        frame->GetYaxis()->SetRangeUser(0.70, 1.22);
        frame->SetYTitle("#sigma_{masked} / #sigma_{nominal}");
        frame->GetYaxis()->SetTitleSize(0.058);
        frame->GetYaxis()->SetTitleOffset(1.02);
        frame->GetYaxis()->SetLabelSize(0.055);
        frame->GetYaxis()->SetNdivisions(505);
        if (!is_bottom) {
            frame->GetXaxis()->SetLabelSize(0);
            frame->GetXaxis()->SetTitleSize(0);
        } else {
            frame->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");
            frame->GetXaxis()->SetTitleSize(0.058);
            frame->GetXaxis()->SetLabelSize(0.055);
        }
        frame->GetXaxis()->SetRangeUser(10, 36);
        frame->Draw("AXIS");

        TLine *ln = new TLine(10, 1.0, 36, 1.0);
        ln->SetLineStyle(2); ln->SetLineColor(kGray + 1); ln->SetLineWidth(2);
        ln->Draw();

        TLegend *leg = new TLegend(0.50, 0.54, 0.94, 0.85);
        leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.048);
        leg->SetHeader(fam.name, "L");
        for (const auto &v : fam.vars) {
            TH1F *r = make_ratio(v);
            r->Draw("P same E1");
            leg->AddEntry(r, v.label, "LP");
        }
        leg->Draw();

        if (is_top) {
            TLatex lx; lx.SetNDC();
            lx.SetTextSize(0.065);
            lx.DrawLatex(0.16, 0.93, strleg1.c_str());
            lx.SetTextSize(0.045);
            lx.DrawLatex(0.16, 0.88, "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, 64.4 pb^{-1} (all-z)");
        }
    };

    draw_panel(c->cd(1), F[0], 1, 3);
    draw_panel(c->cd(2), F[1], 2, 3);
    draw_panel(c->cd(3), F[2], 3, 3);

    const char *outdir = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures";
    c->SaveAs(Form("%s/mask_variant_bybin.pdf", outdir));
    c->SaveAs(Form("%s/mask_variant_bybin.png", outdir));
}
