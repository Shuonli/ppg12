#include "plotcommon.h"

// =====================================================================
// plot_energy_response.C
//
// Plot photon energy response (ET_reco / pT_truth) across interaction
// categories and vertex-reweighting variants from energy_response_comparison.root.
//
// Output (sPHENIX style, PDF):
//   plotting/figures/energy_response/
//     resp_mean_reco.pdf          — mean response vs pT, level=reco
//     resp_sigma_reco.pdf         — sigma response vs pT, level=reco
//     resp_mean_matched.pdf       — mean response vs pT, level=matched
//     resp_sigma_matched.pdf      — sigma response vs pT, level=matched
//     resp_summary_4panel.pdf     — 2x2 panel of mean/sigma × matched/reco
//     resp_projection_pt_<low/mid/high>.pdf
//                                 — 1D projections at 3 representative pT bins
//     resp_vtxrw_diagnostics.pdf  — h_vtxz per sample + reweighting curves
//     resp_pileup_only.pdf        — single/double/0mrad/1p5mrad no-reweight
// =====================================================================

namespace {

struct SetStyle {
    const char *name;
    const char *label;
    int         color;
    int         marker;
    int         linestyle;
};

static const SetStyle kSets[] = {
    {"single",               "Single-interaction MC",             kBlue,       20, 1},
    {"double",               "Double-interaction MC",             kRed+1,      21, 1},
    {"mixed_0mrad",          "Mixed 0 mrad (22.4% DI)",           kBlack,      22, 1},
    {"mixed_0mrad_vtxrw",    "Mixed 0 mrad + vtx reweight",       kGreen+2,    23, 2},
    {"mixed_1p5mrad",        "Mixed 1.5 mrad (7.9% DI)",          kOrange+1,   33, 1},
    {"mixed_1p5mrad_vtxrw",  "Mixed 1.5 mrad + vtx reweight",     kMagenta+1,  34, 2},
};
static const int kNsets = sizeof(kSets) / sizeof(kSets[0]);

// Subset useful for the pileup-only overlay (no vtx rw)
static const int kPileupIdx[] = {0, 1, 2, 4};
static const int kNpileupIdx = 4;

// Subset useful for the vertex-reweight study (paired no-rw vs rw, 0 and 1.5 mrad)
static const int kVtxRwIdx[] = {2, 3, 4, 5};
static const int kNvtxRwIdx  = 4;

static void style_h(TH1 *h, const SetStyle &s, double mrk_size = 1.1)
{
    if (!h) return;
    h->SetMarkerColor(s.color);
    h->SetLineColor(s.color);
    h->SetMarkerStyle(s.marker);
    h->SetMarkerSize(mrk_size);
    h->SetLineWidth(2);
    h->SetLineStyle(s.linestyle);
}

// Draw sPHENIX labels in the upper left
static void draw_labels(float x = 0.22, float y0 = 0.88, float dy = 0.055, float sz = 0.045,
                        const char *extra = nullptr)
{
    TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(sz);
    l.DrawLatex(x, y0,            strleg1.c_str());
    l.DrawLatex(x, y0 -   dy,     strleg2.c_str());
    l.DrawLatex(x, y0 - 2*dy,     strleg3.c_str());
    if (extra) l.DrawLatex(x, y0 - 3*dy, extra);
}

// Single-panel overlay of a set of TH1 summary histograms vs truth pT
//   leg_loc: "lower-left" (default), "upper-right", "upper-left", "lower-right"
static void overlay_summary(TFile *fin,
                             const char *hname_fmt,
                             const char *level,
                             const int *idx, int nidx,
                             double ylo, double yhi,
                             const char *ytitle,
                             const char *extra_label,
                             const char *outpdf,
                             const char *leg_loc = "lower-left")
{
    TCanvas c("c_summary", "", 800, 650);
    gPad->SetLeftMargin(0.17);
    gPad->SetBottomMargin(0.17);

    TH1F *frame = new TH1F("fr", "", 1, 7, 40);
    frame->SetXTitle("#it{p}_{T}^{#gamma, truth} [GeV]");
    frame->SetYTitle(ytitle);
    frame->GetYaxis()->SetRangeUser(ylo, yhi);
    frame->GetXaxis()->SetRangeUser(7, 40);
    frame->Draw();

    // Pick legend coordinates based on requested corner (NDC)
    double lx1 = 0.22, ly1 = 0.18, lx2 = 0.62, ly2 = 0.18 + 0.05 * nidx;
    float label_x = 0.22, label_y0 = 0.88;
    if (std::string(leg_loc) == "upper-right") {
        // narrow column on the right; sPHENIX labels stay upper-left
        lx1 = 0.48; lx2 = 0.93;
        ly2 = 0.93; ly1 = ly2 - 0.045 * nidx;
    } else if (std::string(leg_loc) == "upper-left") {
        lx1 = 0.22; lx2 = 0.62;
        ly2 = 0.88; ly1 = ly2 - 0.05 * nidx;
        // push sPHENIX labels down so they don't collide
        label_y0 = 0.55;
    } else if (std::string(leg_loc) == "lower-right") {
        lx1 = 0.52; lx2 = 0.93;
        ly1 = 0.18; ly2 = 0.18 + 0.05 * nidx;
    }

    TLegend leg(lx1, ly1, lx2, ly2);
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextFont(42); leg.SetTextSize(0.038);

    std::vector<TH1*> keep;
    for (int ii = 0; ii < nidx; ii++) {
        int i = idx[ii];
        TH1 *h = (TH1*)fin->Get(Form(hname_fmt, kSets[i].name, level));
        if (!h) { std::cerr << "MISSING: " << Form(hname_fmt, kSets[i].name, level) << std::endl; continue; }
        style_h(h, kSets[i]);
        h->Draw("E1 SAME");
        leg.AddEntry(h, kSets[i].label, "lp");
        keep.push_back(h);
    }
    leg.Draw();
    draw_labels(label_x, label_y0, 0.055, 0.042, extra_label);
    c.SaveAs(outpdf);
}

// 1D projection of the 2D response at a given pT bin
static void overlay_projection(TFile *fin,
                                const char *obs_key,        // respET or respE
                                const char *level,
                                const int *idx, int nidx,
                                int ipt,                     // 1-indexed pT bin
                                const char *pt_label,
                                const char *obs_ytitle,
                                const char *outpdf)
{
    TCanvas c("c_proj", "", 800, 650);
    gPad->SetLeftMargin(0.17);
    gPad->SetBottomMargin(0.17);
    gPad->SetLogy(0);

    TH1F *frame = new TH1F("frpr", "", 250, 0, 2.5);
    frame->SetXTitle(obs_ytitle);
    frame->SetYTitle("Normalized entries");
    frame->GetYaxis()->SetRangeUser(0, 0.06);
    frame->GetXaxis()->SetRangeUser(0.2, 1.6);
    frame->Draw();

    TLegend leg(0.58, 0.58, 0.92, 0.58 + 0.05 * nidx);
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextFont(42); leg.SetTextSize(0.035);

    double ymax = 0;
    std::vector<TH1*> keep;
    for (int ii = 0; ii < nidx; ii++) {
        int i = idx[ii];
        TH2F *h2 = (TH2F*)fin->Get(Form("h_%s_%s_%s", obs_key, kSets[i].name, level));
        if (!h2) { std::cerr << "MISSING 2D: h_" << obs_key << "_" << kSets[i].name << "_" << level << std::endl; continue; }
        TH1D *hp = h2->ProjectionY(Form("_p_%s_%s_%s_%d", obs_key, kSets[i].name, level, ipt), ipt, ipt);
        if (hp->Integral() <= 0) { delete hp; continue; }
        hp->Scale(1.0 / hp->Integral("width"));
        hp->Scale(hp->GetBinWidth(1));  // back to normalized entries per bin
        style_h(hp, kSets[i]);
        hp->Draw("HIST SAME");
        leg.AddEntry(hp, kSets[i].label, "l");
        ymax = std::max(ymax, hp->GetMaximum());
        keep.push_back(hp);
    }
    frame->GetYaxis()->SetRangeUser(0, std::max(0.04, 1.3 * ymax));
    leg.Draw();

    draw_labels(0.22, 0.88, 0.055, 0.042, pt_label);

    TLine line; line.SetLineStyle(7); line.SetLineColor(kGray+2);
    line.DrawLine(1.0, 0, 1.0, frame->GetMaximum());

    c.SaveAs(outpdf);
}

}  // anonymous namespace

void plot_energy_response()
{
    init_plot();

    const std::string infile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/energy_response_comparison.root";
    const std::string outDir = "/sphenix/user/shuhangli/ppg12/plotting/figures/energy_response";
    gSystem->mkdir(outDir.c_str(), true);

    TFile *fin = TFile::Open(infile.c_str());
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open " << infile << std::endl;
        return;
    }

    // -----------------------------------------------------------------
    // All-6 overlay — mean and sigma summaries at "reco" level
    // -----------------------------------------------------------------
    int kAllIdx[kNsets]; for (int i = 0; i < kNsets; i++) kAllIdx[i] = i;

    overlay_summary(fin, "h_respET_mean_%s_%s",  "reco",
        kAllIdx, kNsets, 0.90, 1.05,
        "#LTE_{T}^{reco} / #it{p}_{T}^{truth}#GT",
        "post common cuts (reco)",
        (outDir + "/resp_mean_reco.pdf").c_str());
    overlay_summary(fin, "h_respET_sigma_%s_%s", "reco",
        kAllIdx, kNsets, 0.02, 0.12,
        "#sigma(E_{T}^{reco} / #it{p}_{T}^{truth})",
        "post common cuts (reco)",
        (outDir + "/resp_sigma_reco.pdf").c_str());

    overlay_summary(fin, "h_respET_mean_%s_%s",  "matched",
        kAllIdx, kNsets, 0.82, 1.05,
        "#LTE_{T}^{reco} / #it{p}_{T}^{truth}#GT",
        "trkID-matched only (matched)",
        (outDir + "/resp_mean_matched.pdf").c_str());
    overlay_summary(fin, "h_respET_sigma_%s_%s", "matched",
        kAllIdx, kNsets, 0.02, 0.20,
        "#sigma(E_{T}^{reco} / #it{p}_{T}^{truth})",
        "trkID-matched only (matched)",
        (outDir + "/resp_sigma_matched.pdf").c_str());

    // -----------------------------------------------------------------
    // Pileup-only overlay (no vtx reweighting) at reco level
    // -----------------------------------------------------------------
    overlay_summary(fin, "h_respET_mean_%s_%s",  "reco",
        kPileupIdx, kNpileupIdx, 0.92, 1.03,
        "#LTE_{T}^{reco} / #it{p}_{T}^{truth}#GT",
        "post common cuts (reco)",
        (outDir + "/resp_mean_reco_pileup.pdf").c_str());
    overlay_summary(fin, "h_respET_sigma_%s_%s", "reco",
        kPileupIdx, kNpileupIdx, 0.02, 0.10,
        "#sigma(E_{T}^{reco} / #it{p}_{T}^{truth})",
        "post common cuts (reco)",
        (outDir + "/resp_sigma_reco_pileup.pdf").c_str());

    // -----------------------------------------------------------------
    // Vertex-reweighting overlay (0 mrad / 1.5 mrad, w/wo rw) at reco level
    // -----------------------------------------------------------------
    // Data slopes down from upper-left to lower-right; move legend to upper-right
    // where the panel is empty (previously lower-left overlapped the data curves).
    overlay_summary(fin, "h_respET_mean_%s_%s",  "reco",
        kVtxRwIdx, kNvtxRwIdx, 0.94, 1.02,
        "#LTE_{T}^{reco} / #it{p}_{T}^{truth}#GT",
        "effect of z_{vtx} reweighting",
        (outDir + "/resp_mean_reco_vtxrw.pdf").c_str(),
        "upper-right");
    // Data values for mixed_0mrad variants extend up to ~0.093 at low pT;
    // widen yhi so they are not clipped (previously clipped at 0.09, making
    // the 0 mrad curves disappear for pT < 20 GeV).
    overlay_summary(fin, "h_respET_sigma_%s_%s", "reco",
        kVtxRwIdx, kNvtxRwIdx, 0.03, 0.11,
        "#sigma(E_{T}^{reco} / #it{p}_{T}^{truth})",
        "effect of z_{vtx} reweighting",
        (outDir + "/resp_sigma_reco_vtxrw.pdf").c_str(),
        "lower-right");

    // -----------------------------------------------------------------
    // Response projections at 3 representative pT bins (matched level)
    //   pT bins (pT_bins_truth from config): 7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45
    //   low = bin 3 (10-12 GeV), mid = bin 9 (22-24), high = bin 12 (28-32)
    // -----------------------------------------------------------------
    overlay_projection(fin, "respET", "matched", kPileupIdx, kNpileupIdx, 3,
        "10 < #it{p}_{T}^{truth} < 12 GeV",
        "E_{T}^{reco} / #it{p}_{T}^{truth}",
        (outDir + "/resp_projection_matched_pt10_12.pdf").c_str());
    overlay_projection(fin, "respET", "matched", kPileupIdx, kNpileupIdx, 9,
        "22 < #it{p}_{T}^{truth} < 24 GeV",
        "E_{T}^{reco} / #it{p}_{T}^{truth}",
        (outDir + "/resp_projection_matched_pt22_24.pdf").c_str());
    overlay_projection(fin, "respET", "matched", kPileupIdx, kNpileupIdx, 12,
        "28 < #it{p}_{T}^{truth} < 32 GeV",
        "E_{T}^{reco} / #it{p}_{T}^{truth}",
        (outDir + "/resp_projection_matched_pt28_32.pdf").c_str());
    overlay_projection(fin, "respET", "reco",    kPileupIdx, kNpileupIdx, 9,
        "22 < #it{p}_{T}^{truth} < 24 GeV (post common)",
        "E_{T}^{reco} / #it{p}_{T}^{truth}",
        (outDir + "/resp_projection_reco_pt22_24.pdf").c_str());

    // -----------------------------------------------------------------
    // 4-panel summary figure: mean/sigma × matched/reco (all 6 sets)
    // -----------------------------------------------------------------
    {
        TCanvas c("c4", "", 1400, 1100);
        c.Divide(2, 2, 0.002, 0.002);

        auto drawPanel = [&](int ipad, const char *hfmt, const char *lvl,
                             double ylo, double yhi, const char *yt,
                             const char *lvl_label, const char *leg_loc)
        {
            c.cd(ipad);
            gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.17);
            TH1F *fr = new TH1F(Form("fr_%d", ipad), "", 1, 7, 40);
            fr->SetXTitle("#it{p}_{T}^{#gamma, truth} [GeV]");
            fr->SetYTitle(yt);
            fr->GetYaxis()->SetRangeUser(ylo, yhi);
            fr->GetXaxis()->SetRangeUser(7, 40);
            fr->Draw();
            for (int i = 0; i < kNsets; i++) {
                TH1 *h = (TH1*)fin->Get(Form(hfmt, kSets[i].name, lvl));
                if (!h) continue;
                style_h(h, kSets[i]);
                h->Draw("E1 SAME");
            }
            TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.040);
            l.DrawLatex(0.22, 0.90, strleg1.c_str());
            l.DrawLatex(0.22, 0.85, strleg2.c_str());
            l.DrawLatex(0.22, 0.80, strleg3.c_str());
            l.DrawLatex(0.22, 0.75, lvl_label);

            // Per-panel series legend. Abbreviated 6-entry labels in a 2-column
            // layout fit into the upper-right strip that is the only region
            // consistently empty across all 4 panels (sPHENIX header at
            // upper-left; data on mean panels ends at NDC y=0.69, on sigma
            // panels the double curve ends at NDC y=0.80).
            static const char *kShortLabels[6] = {
                "Single",           "Double",
                "Mix 0 mrad",       "Mix 0 mrad +vtxrw",
                "Mix 1.5 mrad",     "Mix 1.5 mrad +vtxrw",
            };
            TLegend *leg = new TLegend(0.44, 0.74, 0.93, 0.93);
            leg->SetFillStyle(0); leg->SetBorderSize(0);
            leg->SetTextFont(42); leg->SetTextSize(0.026);
            leg->SetNColumns(2);
            for (int i = 0; i < kNsets; i++) {
                TH1 *h = (TH1*)fin->Get(Form(hfmt, kSets[i].name, lvl));
                if (!h) continue;
                leg->AddEntry(h, kShortLabels[i], "lp");
            }
            leg->Draw();
        };

        // All 4 panels carry the same 6-entry series legend (upper-right,
        // 2-column) so the figure is readable in isolation without needing
        // a separate legend key.
        drawPanel(1, "h_respET_mean_%s_%s",  "matched", 0.86, 1.02,
                  "#LTE_{T}^{reco} / #it{p}_{T}^{truth}#GT", "matched", "upper-right");
        drawPanel(2, "h_respET_mean_%s_%s",  "reco",    0.88, 1.05,
                  "#LTE_{T}^{reco} / #it{p}_{T}^{truth}#GT", "reco", "upper-right");
        drawPanel(3, "h_respET_sigma_%s_%s", "matched", 0.02, 0.22,
                  "#sigma(E_{T}^{reco} / #it{p}_{T}^{truth})", "matched", "upper-right");
        drawPanel(4, "h_respET_sigma_%s_%s", "reco",    0.02, 0.16,
                  "#sigma(E_{T}^{reco} / #it{p}_{T}^{truth})", "reco", "upper-right");

        c.SaveAs((outDir + "/resp_summary_4panel.pdf").c_str());
    }

    // -----------------------------------------------------------------
    // Diagnostics: vertex distributions and reweighting curves
    // -----------------------------------------------------------------
    {
        TCanvas c("cvtx", "", 1200, 550);
        c.Divide(2, 1, 0.002, 0.002);

        // TLatex reused for both panels: sPHENIX header + panel title
        TLatex ltx; ltx.SetNDC(); ltx.SetTextFont(42);

        // Heap-allocate legends so they outlive the if-block scope (otherwise
        // the stack TLegend is destroyed before c.SaveAs runs, leaving the pad
        // with no legend rendered in the final PDF).
        TLegend *leg1 = nullptr;
        TLegend *leg2 = nullptr;

        c.cd(1);
        gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.06);
        TH1 *h_single = (TH1*)fin->Get("h_vtxz_single");
        TH1 *h_double = (TH1*)fin->Get("h_vtxz_double");
        if (h_single && h_double) {
            h_single->Scale(1.0 / std::max(1.0, h_single->Integral()));
            h_double->Scale(1.0 / std::max(1.0, h_double->Integral()));
            h_single->SetLineColor(kBlue);    h_single->SetLineWidth(2);
            h_double->SetLineColor(kRed+1);   h_double->SetLineWidth(2);
            h_single->SetTitle(";z_{vtx} [cm];Normalized entries");
            h_single->GetYaxis()->SetRangeUser(0, 1.35 * std::max(h_single->GetMaximum(),
                                                                   h_double->GetMaximum()));
            h_single->Draw("HIST");
            h_double->Draw("HIST SAME");
            // sPHENIX labels + panel title (upper-left)
            ltx.SetTextSize(0.045);
            ltx.DrawLatex(0.20, 0.90, "MC Vertex Distributions");
            ltx.SetTextSize(0.040);
            ltx.DrawLatex(0.20, 0.85, strleg1.c_str());
            ltx.DrawLatex(0.20, 0.80, strleg2.c_str());
            // Legend: middle-right (the peaks sit near z=0 at the top-center,
            // and the tails drop off outside |z|>60 cm, leaving room on the right)
            leg1 = new TLegend(0.52, 0.55, 0.93, 0.70);
            leg1->SetFillStyle(0); leg1->SetBorderSize(0);
            leg1->SetTextFont(42); leg1->SetTextSize(0.040);
            leg1->AddEntry(h_single, "photon10 (single) MC", "l");
            leg1->AddEntry(h_double, "photon10_double MC",   "l");
            leg1->Draw();
        }

        c.cd(2);
        gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.06);
        TH1 *w0  = (TH1*)fin->Get("h_vtx_weight_mixed_0mrad_vtxrw");
        TH1 *w15 = (TH1*)fin->Get("h_vtx_weight_mixed_1p5mrad_vtxrw");
        if (w0 && w15) {
            w0 ->SetLineColor(kGreen+2);  w0 ->SetLineWidth(2);
            w15->SetLineColor(kMagenta+1); w15->SetLineWidth(2);
            w0->SetTitle(";z_{vtx} [cm];data / sim ratio");
            double m = std::max(w0->GetMaximum(), w15->GetMaximum());
            w0->GetYaxis()->SetRangeUser(0, 1.45 * m);
            w0->GetXaxis()->SetRangeUser(-60, 60);
            w0 ->Draw("HIST");
            w15->Draw("HIST SAME");
            // sPHENIX labels + panel title (upper-left)
            ltx.SetTextSize(0.045);
            ltx.DrawLatex(0.20, 0.90, "Data / MC Vertex Reweighting");
            ltx.SetTextSize(0.040);
            ltx.DrawLatex(0.20, 0.85, strleg1.c_str());
            ltx.DrawLatex(0.20, 0.80, strleg2.c_str());
            // Legend: middle-right (1.5 mrad curve peaks at z~0; 0 mrad curve
            // is flat-ish so left edges are busy — right side is cleaner)
            leg2 = new TLegend(0.52, 0.60, 0.93, 0.74);
            leg2->SetFillStyle(0); leg2->SetBorderSize(0);
            leg2->SetTextFont(42); leg2->SetTextSize(0.038);
            leg2->AddEntry(w0,  "0 mrad (data / sim)",   "l");
            leg2->AddEntry(w15, "1.5 mrad (data / sim)", "l");
            leg2->Draw();
        }

        c.SaveAs((outDir + "/resp_vtxrw_diagnostics.pdf").c_str());
        delete leg1;
        delete leg2;
    }

    std::cout << "[plot_energy_response] done: figures in " << outDir << std::endl;
    fin->Close();
    delete fin;
}
