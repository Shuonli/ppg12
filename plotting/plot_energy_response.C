#include "plotcommon.h"

// Double-sided Crystal-Ball PDF (same definition as the producer).
// par[0]=N par[1]=mu par[2]=sigma par[3]=alpha_L par[4]=n_L par[5]=alpha_H par[6]=n_H
static Double_t dscb_fn(Double_t *x, Double_t *par)
{
    Double_t N  = par[0];
    Double_t mu = par[1];
    Double_t sg = par[2];
    Double_t aL = par[3], nL = par[4];
    Double_t aH = par[5], nH = par[6];
    if (sg <= 0) return 0;
    Double_t t = (x[0] - mu) / sg;
    if (t > -aL && t < aH)
        return N * TMath::Exp(-0.5 * t * t);
    if (t <= -aL) {
        Double_t A = TMath::Power(nL / TMath::Abs(aL), nL) * TMath::Exp(-0.5 * aL * aL);
        Double_t B = nL / TMath::Abs(aL) - TMath::Abs(aL);
        return N * A * TMath::Power(B - t, -nL);
    }
    Double_t A = TMath::Power(nH / TMath::Abs(aH), nH) * TMath::Exp(-0.5 * aH * aH);
    Double_t B = nH / TMath::Abs(aH) - TMath::Abs(aH);
    return N * A * TMath::Power(B + t, -nH);
}

// =====================================================================
// plot_energy_response.C
//
// Plot photon energy response (ET_reco / pT_truth) across interaction
// categories and vertex-reweighting variants from energy_response_comparison.root.
//
// Level naming convention (producer: compare_energy_response.C):
//   reco       : trkID match only, no extra cuts
//   common_cut : reco + common cuts (cluster probability, e11/e33, omega_r,
//                omega_eta, NPB if enabled); working point that enters ABCD
//   tight_iso  : common_cut + isolation + ET-dependent tight BDT (signal
//                region A of ABCD)
//
// Estimators stored by the producer (both available for every set/level):
//   (Gauss)  h_respET_{mean,sigma}_{set}_{level}      — truncated Gaussian on [0.3, 1.7]
//   (CB)     h_respET_cb_{mean,sigma}_{set}_{level}   — Crystal-Ball on [0.3, 2.0]
//
// Output (sPHENIX style, PDF):
//   plotting/figures/energy_response/
//     resp_mean_reco.pdf                  — Gauss mean vs pT, reco
//     resp_sigma_reco.pdf                 — Gauss sigma vs pT, reco
//     resp_mean_common_cut.pdf            — Gauss mean vs pT, common_cut
//     resp_sigma_common_cut.pdf           — Gauss sigma vs pT, common_cut
//     resp_mean_cb_reco.pdf               — CB peak vs pT, reco
//     resp_sigma_cb_reco.pdf              — CB sigma vs pT, reco
//     resp_mean_cb_common_cut.pdf         — CB peak vs pT, common_cut
//     resp_sigma_cb_common_cut.pdf        — CB sigma vs pT, common_cut
//     resp_summary_4panel.pdf             — 2x2 panel of Gauss mean/sigma × reco/common_cut
//     resp_summary_4panel_cb.pdf          — 2x2 panel of CB peak/sigma × reco/common_cut
//     resp_projection_*_pt_*.pdf          — 1D projections at representative pT bins
//     resp_projection_fits_{single,double}_{reco,common_cut}_pt22_24.pdf
//                                         — 1D projection with Gauss + CB fit curves overlaid
//     resp_fit_qc_grid_{set}_{reco,common_cut}.pdf
//                                         — 4x3 QC grid: projection + Gauss/CB fits per pT bin
//     resp_fit_qc_grid_dscb_{single,double,mixed_1p5mrad}_tight_iso.pdf
//                                         — 4x3 QC grid at tight_iso with Gauss + CB + DSCB
//     resp_vtxrw_diagnostics.pdf          — h_vtxz per sample + reweighting curves
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
    frame->GetXaxis()->SetRangeUser(10, 40);
    frame->Draw();

    // Pick legend coordinates based on requested corner (NDC)
    double lx1 = 0.22, ly1 = 0.18, lx2 = 0.62, ly2 = 0.18 + 0.05 * nidx;
    float label_x = 0.22, label_y0 = 0.88;
    if (std::string(leg_loc) == "upper-right") {
        // narrow column on the right; sPHENIX labels stay upper-left. lx1
        // bumped from 0.48 -> 0.54 so the legend does not overlap the
        // right edge of the sPHENIX/lumi/|eta| header block (which ends
        // near NDC x ~ 0.50 at label text size 0.042).
        lx1 = 0.54; lx2 = 0.93;
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

// 1D projection of ONE set at one pT bin, overlaid with both Gaussian
// [0.3, 1.7] and Crystal-Ball [0.3, 2.0] fits. Diagnostic: shows how much
// of the low-side tail the Gauss misses and the CB absorbs.
static void overlay_projection_with_fits(TFile *fin,
                                          const char *obs_key,
                                          int set_idx,
                                          const char *level,
                                          int ipt,
                                          const char *pt_label,
                                          const char *obs_ytitle,
                                          const char *outpdf)
{
    TCanvas c("c_projfit", "", 900, 700);
    gPad->SetLeftMargin(0.17);
    gPad->SetBottomMargin(0.17);

    TH2F *h2 = (TH2F*)fin->Get(Form("h_%s_%s_%s", obs_key, kSets[set_idx].name, level));
    if (!h2) { std::cerr << "MISSING 2D: h_" << obs_key << "_"
                          << kSets[set_idx].name << "_" << level << std::endl; return; }

    TH1D *hp = h2->ProjectionY(Form("_pfit_%s_%s_%s_%d", obs_key, kSets[set_idx].name, level, ipt),
                                ipt, ipt);
    if (hp->Integral() <= 0) { delete hp; return; }

    // Keep the histogram as raw counts so the fit norm constant reads naturally.
    hp->SetMarkerColor(kBlack); hp->SetLineColor(kBlack);
    hp->SetMarkerStyle(20); hp->SetMarkerSize(0.9);
    hp->SetLineWidth(1);
    hp->GetXaxis()->SetTitle(obs_ytitle);
    hp->GetYaxis()->SetTitle("entries");
    hp->GetXaxis()->SetRangeUser(0.2, 1.7);
    double ymax = hp->GetMaximum() * 1.30;
    hp->GetYaxis()->SetRangeUser(0, ymax);
    hp->Draw("E1");

    // Gaussian fit
    TF1 *fG = new TF1(Form("fG_%s_%s_%s_%d", obs_key, kSets[set_idx].name, level, ipt),
                      "gaus", 0.3, 1.7);
    fG->SetParameter(0, hp->GetMaximum());
    fG->SetParameter(1, hp->GetMean());
    fG->SetParameter(2, std::max(0.05, (double)hp->GetStdDev()));
    fG->SetLineColor(kBlue+1); fG->SetLineStyle(2); fG->SetLineWidth(2);
    hp->Fit(fG, "RQN");
    fG->SetRange(0.2, 1.7);
    fG->Draw("SAME");

    // Crystal-Ball fit
    TF1 *fCB = new TF1(Form("fCB_%s_%s_%s_%d", obs_key, kSets[set_idx].name, level, ipt),
                       "crystalball", 0.3, 2.0);
    fCB->SetParameter(0, fG->GetParameter(0));
    fCB->SetParameter(1, fG->GetParameter(1));
    fCB->SetParameter(2, std::max(0.03, (double)fG->GetParameter(2)));
    fCB->SetParameter(3, 1.0);
    fCB->SetParameter(4, 5.0);
    fCB->SetParLimits(2, 0.005, 0.4);
    fCB->SetParLimits(3, 0.1, 10.0);
    fCB->SetParLimits(4, 1.05, 100.0);
    fCB->SetLineColor(kRed+1); fCB->SetLineStyle(1); fCB->SetLineWidth(2);
    hp->Fit(fCB, "RQN");
    fCB->SetRange(0.2, 2.0);
    fCB->Draw("SAME");

    draw_labels(0.22, 0.88, 0.055, 0.040, pt_label);

    TLatex lbl; lbl.SetNDC(); lbl.SetTextFont(42); lbl.SetTextSize(0.038);
    lbl.DrawLatex(0.22, 0.64, kSets[set_idx].label);
    lbl.DrawLatex(0.22, 0.59, Form("level: %s", level));

    // Wider legend + smaller text so the full #mu=X.XXX, #sigma=X.XXX
    // fits inside the canvas right edge.
    TLegend *leg = new TLegend(0.43, 0.72, 0.95, 0.90);
    leg->SetFillStyle(0); leg->SetBorderSize(0);
    leg->SetTextFont(42); leg->SetTextSize(0.030);
    leg->AddEntry(hp,  "MC projection", "lp");
    leg->AddEntry(fG,  Form("Gauss:  #mu=%.3f, #sigma=%.3f",
                            fG->GetParameter(1), fG->GetParameter(2)), "l");
    leg->AddEntry(fCB, Form("CB:  #mu=%.3f, #sigma=%.3f, #alpha=%.2f",
                             fCB->GetParameter(1), fCB->GetParameter(2),
                             fCB->GetParameter(3)), "l");
    leg->Draw();

    // Reference line at r=1
    TLine line; line.SetLineStyle(7); line.SetLineColor(kGray+2);
    line.DrawLine(1.0, 0, 1.0, ymax);

    c.SaveAs(outpdf);

    delete fG; delete fCB; delete leg; delete hp;
}

// 12-panel QC grid (4 cols × 3 rows): 1D R_ET projection with both Gauss
// and CB fits overlaid at each pT bin in the analysis range (8-36 GeV).
// Used for the report appendix to document fit convergence quality.
//
// pT bin mapping (pT_bins_truth = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45}):
//   ipt=2..12  →  (10-12), (12-14), ..., (32-36) GeV  (11 panels)
static void draw_fit_qc_grid(TFile *fin,
                              const char *obs_key,
                              int set_idx,
                              const char *level,
                              const char *outpdf)
{
    static const int ipt_range[11] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    static const double pT_edges[14] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45};

    // Extra top margin on the canvas so the figure-level title above the
    // pad grid does not overlap the top-row y-axis tick labels.
    TCanvas c("c_qcgrid", "", 1600, 1180);
    TPad *pad_title = new TPad("pad_title", "pad_title", 0.0, 0.955, 1.0, 1.0);
    pad_title->SetFillStyle(0); pad_title->SetBorderMode(0);
    pad_title->Draw();
    TPad *pad_grid  = new TPad("pad_grid",  "pad_grid",  0.0, 0.0,   1.0, 0.955);
    pad_grid->SetFillStyle(0); pad_grid->SetBorderMode(0);
    pad_grid->Draw();
    pad_grid->cd();
    pad_grid->Divide(4, 3, 0.003, 0.003);

    TH2F *h2 = (TH2F*)fin->Get(Form("h_%s_%s_%s", obs_key, kSets[set_idx].name, level));
    if (!h2) {
        std::cerr << "[qc_grid] MISSING 2D: h_" << obs_key << "_"
                  << kSets[set_idx].name << "_" << level << std::endl;
        return;
    }

    std::vector<TF1*> keepFits;
    std::vector<TH1D*> keepProj;

    for (int k = 0; k < 11; k++) {
        int ipt = ipt_range[k];
        pad_grid->cd(k + 1);
        gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.17); gPad->SetTopMargin(0.04); gPad->SetRightMargin(0.03);

        TH1D *hp = h2->ProjectionY(Form("_qc_%s_%s_%s_%d", obs_key, kSets[set_idx].name, level, ipt),
                                    ipt, ipt);
        if (hp->Integral() <= 0) { delete hp; continue; }

        hp->SetMarkerColor(kBlack); hp->SetLineColor(kBlack);
        hp->SetMarkerStyle(20); hp->SetMarkerSize(0.7);
        hp->GetXaxis()->SetTitle("E_{T}^{reco} / p_{T}^{truth}");
        hp->GetYaxis()->SetTitle("entries");
        hp->GetXaxis()->SetRangeUser(0.2, 1.7);
        double ymax = hp->GetMaximum() * 1.55;
        hp->GetYaxis()->SetRangeUser(0, ymax);
        hp->GetXaxis()->SetTitleSize(0.06);
        hp->GetYaxis()->SetTitleSize(0.06);
        hp->GetXaxis()->SetLabelSize(0.055);
        hp->GetYaxis()->SetLabelSize(0.055);
        hp->GetYaxis()->SetTitleOffset(1.1);
        hp->Draw("E1");

        // Thinner fit lines (width 1) so they don't overlay the data markers.
        TF1 *fG = new TF1(Form("fG_qc_%s_%s_%s_%d", obs_key, kSets[set_idx].name, level, ipt),
                          "gaus", 0.3, 1.7);
        fG->SetParameter(0, hp->GetMaximum());
        fG->SetParameter(1, hp->GetMean());
        fG->SetParameter(2, std::max(0.05, (double)hp->GetStdDev()));
        fG->SetLineColor(kBlue+1); fG->SetLineStyle(2); fG->SetLineWidth(1);
        hp->Fit(fG, "RQN");
        fG->Draw("SAME");

        TF1 *fCB = new TF1(Form("fCB_qc_%s_%s_%s_%d", obs_key, kSets[set_idx].name, level, ipt),
                           "crystalball", 0.3, 2.0);
        fCB->SetParameter(0, fG->GetParameter(0));
        fCB->SetParameter(1, fG->GetParameter(1));
        fCB->SetParameter(2, std::max(0.03, (double)fG->GetParameter(2)));
        fCB->SetParameter(3, 1.0);
        fCB->SetParameter(4, 5.0);
        fCB->SetParLimits(2, 0.005, 0.4);
        fCB->SetParLimits(3, 0.1, 10.0);
        fCB->SetParLimits(4, 1.05, 100.0);
        fCB->SetLineColor(kRed+1); fCB->SetLineStyle(1); fCB->SetLineWidth(1);
        hp->Fit(fCB, "RQN");
        fCB->Draw("SAME");

        TLatex lbl; lbl.SetNDC(); lbl.SetTextFont(42);
        lbl.SetTextSize(0.065);
        lbl.DrawLatex(0.20, 0.90, Form("%g < p_{T}^{truth} < %g GeV",
                                        pT_edges[ipt-1], pT_edges[ipt]));
        lbl.SetTextSize(0.055);
        lbl.SetTextColor(kBlue+1);
        lbl.DrawLatex(0.20, 0.82, Form("G: #mu=%.3f, #sigma=%.3f",
                                        fG->GetParameter(1), fG->GetParameter(2)));
        lbl.SetTextColor(kRed+1);
        lbl.DrawLatex(0.20, 0.75, Form("CB: #mu=%.3f, #sigma=%.3f, #alpha=%.2f",
                                        fCB->GetParameter(1), fCB->GetParameter(2),
                                        fCB->GetParameter(3)));
        lbl.SetTextColor(kBlack);

        // Reference line at R_ET = 1
        TLine *line = new TLine(1.0, 0, 1.0, ymax);
        line->SetLineStyle(7); line->SetLineColor(kGray+2);
        line->Draw();

        keepFits.push_back(fG);
        keepFits.push_back(fCB);
        keepProj.push_back(hp);
    }

    // Figure-wide title drawn inside the dedicated top pad (ASCII-safe: the
    // earlier em-dash rendered as mojibake on some backends).
    pad_title->cd();
    TLatex top; top.SetNDC(); top.SetTextFont(42); top.SetTextSize(0.55);
    top.SetTextAlign(22);  // centered
    top.DrawLatex(0.5, 0.5, Form("#bf{#it{sPHENIX}} Internal -- %s (%s) -- Gauss (blue dashed) vs Crystal-Ball (red) fits",
                                   kSets[set_idx].label, level));

    c.SaveAs(outpdf);

    for (auto *f : keepFits) delete f;
    for (auto *h : keepProj) delete h;
}

// Like draw_fit_qc_grid but also overlays a double-sided CB fit (green
// solid). Used for tight_iso level where the one-sided CB misses the
// residual high-side tail from isolation leakage.
static void draw_fit_qc_grid_dscb(TFile *fin,
                                   const char *obs_key,
                                   int set_idx,
                                   const char *level,
                                   const char *outpdf)
{
    static const int ipt_range[11] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    static const double pT_edges[14] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45};

    TCanvas c("c_qcgrid_dscb", "", 1600, 1180);
    TPad *pad_title = new TPad("pad_title_d", "", 0.0, 0.955, 1.0, 1.0);
    pad_title->SetFillStyle(0); pad_title->SetBorderMode(0);
    pad_title->Draw();
    TPad *pad_grid  = new TPad("pad_grid_d",  "", 0.0, 0.0, 1.0, 0.955);
    pad_grid->SetFillStyle(0); pad_grid->SetBorderMode(0);
    pad_grid->Draw();
    pad_grid->cd();
    pad_grid->Divide(4, 3, 0.003, 0.003);

    TH2F *h2 = (TH2F*)fin->Get(Form("h_%s_%s_%s", obs_key, kSets[set_idx].name, level));
    if (!h2) {
        std::cerr << "[qc_grid_dscb] MISSING 2D: h_" << obs_key << "_"
                  << kSets[set_idx].name << "_" << level << std::endl;
        return;
    }

    std::vector<TF1*>  keepFits;
    std::vector<TH1D*> keepProj;

    for (int k = 0; k < 11; k++) {
        int ipt = ipt_range[k];
        pad_grid->cd(k + 1);
        gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.17); gPad->SetTopMargin(0.04); gPad->SetRightMargin(0.03);

        TH1D *hp = h2->ProjectionY(Form("_qcD_%s_%s_%s_%d", obs_key, kSets[set_idx].name, level, ipt),
                                    ipt, ipt);
        if (hp->Integral() <= 0) { delete hp; continue; }

        hp->SetMarkerColor(kBlack); hp->SetLineColor(kBlack);
        hp->SetMarkerStyle(20); hp->SetMarkerSize(0.7);
        hp->GetXaxis()->SetTitle("E_{T}^{reco} / p_{T}^{truth}");
        hp->GetYaxis()->SetTitle("entries");
        hp->GetXaxis()->SetRangeUser(0.2, 1.7);
        double ymax = hp->GetMaximum() * 1.75;
        hp->GetYaxis()->SetRangeUser(0, ymax);
        hp->GetXaxis()->SetTitleSize(0.06);
        hp->GetYaxis()->SetTitleSize(0.06);
        hp->GetXaxis()->SetLabelSize(0.055);
        hp->GetYaxis()->SetLabelSize(0.055);
        hp->GetYaxis()->SetTitleOffset(1.1);
        hp->Draw("E1");

        // Gauss fit
        TF1 *fG = new TF1(Form("fG_qcD_%s_%s_%s_%d", obs_key, kSets[set_idx].name, level, ipt),
                          "gaus", 0.3, 1.7);
        fG->SetParameter(0, hp->GetMaximum());
        fG->SetParameter(1, hp->GetMean());
        fG->SetParameter(2, std::max(0.05, (double)hp->GetStdDev()));
        fG->SetLineColor(kBlue+1); fG->SetLineStyle(2); fG->SetLineWidth(1);
        hp->Fit(fG, "RQN");
        fG->Draw("SAME");

        // One-sided CB fit
        TF1 *fCB = new TF1(Form("fCB_qcD_%s_%s_%s_%d", obs_key, kSets[set_idx].name, level, ipt),
                           "crystalball", 0.3, 2.0);
        fCB->SetParameter(0, fG->GetParameter(0));
        fCB->SetParameter(1, fG->GetParameter(1));
        fCB->SetParameter(2, std::max(0.03, (double)fG->GetParameter(2)));
        fCB->SetParameter(3, 1.0);
        fCB->SetParameter(4, 5.0);
        fCB->SetParLimits(2, 0.005, 0.4);
        fCB->SetParLimits(3, 0.1, 10.0);
        fCB->SetParLimits(4, 1.05, 100.0);
        fCB->SetLineColor(kRed+1); fCB->SetLineStyle(3); fCB->SetLineWidth(1);
        hp->Fit(fCB, "RQN");
        fCB->Draw("SAME");

        // Double-sided CB fit
        TF1 *fD = new TF1(Form("fD_qcD_%s_%s_%s_%d", obs_key, kSets[set_idx].name, level, ipt),
                          dscb_fn, 0.3, 2.0, 7);
        fD->SetParNames("N", "mu", "sigma", "alphaL", "nL", "alphaH", "nH");
        fD->SetParameter(0, fCB->GetParameter(0));
        fD->SetParameter(1, fCB->GetParameter(1));
        fD->SetParameter(2, std::max(0.03, (double)fCB->GetParameter(2)));
        fD->SetParameter(3, 1.0);
        fD->SetParameter(4, 5.0);
        fD->SetParameter(5, 1.5);
        fD->SetParameter(6, 5.0);
        fD->SetParLimits(2, 0.005, 0.4);
        fD->SetParLimits(3, 0.1, 10.0);
        fD->SetParLimits(4, 1.05, 100.0);
        fD->SetParLimits(5, 0.1, 10.0);
        fD->SetParLimits(6, 1.05, 100.0);
        fD->SetLineColor(kGreen+2); fD->SetLineStyle(1); fD->SetLineWidth(1);
        hp->Fit(fD, "RQN");
        fD->Draw("SAME");

        TLatex lbl; lbl.SetNDC(); lbl.SetTextFont(42);
        lbl.SetTextSize(0.060);
        lbl.DrawLatex(0.20, 0.92, Form("%g < p_{T}^{truth} < %g GeV",
                                        pT_edges[ipt-1], pT_edges[ipt]));
        lbl.SetTextSize(0.048);
        lbl.SetTextColor(kBlue+1);
        lbl.DrawLatex(0.20, 0.85, Form("G: #mu=%.3f, #sigma=%.3f",
                                        fG->GetParameter(1), fG->GetParameter(2)));
        lbl.SetTextColor(kRed+1);
        lbl.DrawLatex(0.20, 0.79, Form("CB: #mu=%.3f, #sigma=%.3f, #alpha=%.2f",
                                        fCB->GetParameter(1), fCB->GetParameter(2),
                                        fCB->GetParameter(3)));
        lbl.SetTextColor(kGreen+2);
        lbl.DrawLatex(0.20, 0.73, Form("DSCB: #mu=%.3f, #sigma=%.3f",
                                        fD->GetParameter(1), fD->GetParameter(2)));
        lbl.DrawLatex(0.20, 0.67, Form("    #alpha_{L}=%.2f, #alpha_{H}=%.2f",
                                        fD->GetParameter(3), fD->GetParameter(5)));
        lbl.SetTextColor(kBlack);

        TLine *line = new TLine(1.0, 0, 1.0, ymax);
        line->SetLineStyle(7); line->SetLineColor(kGray+2);
        line->Draw();

        keepFits.push_back(fG);
        keepFits.push_back(fCB);
        keepFits.push_back(fD);
        keepProj.push_back(hp);
    }

    pad_title->cd();
    TLatex top; top.SetNDC(); top.SetTextFont(42); top.SetTextSize(0.55);
    top.SetTextAlign(22);
    top.DrawLatex(0.5, 0.5, Form("#bf{#it{sPHENIX}} Internal -- %s (%s) -- Gauss (blue) / CB (red) / DSCB (green)",
                                  kSets[set_idx].label, level));

    c.SaveAs(outpdf);

    for (auto *f : keepFits) delete f;
    for (auto *h : keepProj) delete h;
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
    // All-6 overlay — mean and sigma summaries at "common_cut" level
    // -----------------------------------------------------------------
    int kAllIdx[kNsets]; for (int i = 0; i < kNsets; i++) kAllIdx[i] = i;

    overlay_summary(fin, "h_respET_mean_%s_%s",  "common_cut",
        kAllIdx, kNsets, 0.90, 1.05,
        "#LTE_{T}^{reco} / #it{p}_{T}^{truth}#GT",
        "common cut",
        (outDir + "/resp_mean_common_cut.pdf").c_str());
    overlay_summary(fin, "h_respET_sigma_%s_%s", "common_cut",
        kAllIdx, kNsets, 0.02, 0.12,
        "#sigma(E_{T}^{reco} / #it{p}_{T}^{truth})",
        "common cut",
        (outDir + "/resp_sigma_common_cut.pdf").c_str());

    overlay_summary(fin, "h_respET_mean_%s_%s",  "reco",
        kAllIdx, kNsets, 0.82, 1.05,
        "#LTE_{T}^{reco} / #it{p}_{T}^{truth}#GT",
        "reco (trkID match)",
        (outDir + "/resp_mean_reco.pdf").c_str());
    overlay_summary(fin, "h_respET_sigma_%s_%s", "reco",
        kAllIdx, kNsets, 0.02, 0.20,
        "#sigma(E_{T}^{reco} / #it{p}_{T}^{truth})",
        "reco (trkID match)",
        (outDir + "/resp_sigma_reco.pdf").c_str());

    // -----------------------------------------------------------------
    // All-6 overlay — Crystal-Ball summaries (peak mu / sigma) at reco
    // and common_cut. CB absorbs the low-side brem/leakage tail and
    // gives a less-biased peak estimator, especially for double-interaction.
    // -----------------------------------------------------------------
    // CB curves cluster near y=1 at mean panels and descend near data markers
    // at sigma panels — lower-left legend collides with data. Use lower-right
    // for mean (upper-right is occupied by Single/Mixed at 1.00), and upper-
    // right for sigma (data descends down from upper-left there).
    overlay_summary(fin, "h_respET_cb_mean_%s_%s",  "reco",
        kAllIdx, kNsets, 0.86, 1.04,
        "CB peak  #mu(E_{T}^{reco} / #it{p}_{T}^{truth})",
        "reco (CB fit)",
        (outDir + "/resp_mean_cb_reco.pdf").c_str(),
        "lower-right");
    overlay_summary(fin, "h_respET_cb_sigma_%s_%s", "reco",
        kAllIdx, kNsets, 0.02, 0.18,
        "CB  #sigma(E_{T}^{reco} / #it{p}_{T}^{truth})",
        "reco (CB fit)",
        (outDir + "/resp_sigma_cb_reco.pdf").c_str(),
        "upper-right");
    overlay_summary(fin, "h_respET_cb_mean_%s_%s",  "common_cut",
        kAllIdx, kNsets, 0.92, 1.04,
        "CB peak  #mu(E_{T}^{reco} / #it{p}_{T}^{truth})",
        "common cut (CB fit)",
        (outDir + "/resp_mean_cb_common_cut.pdf").c_str(),
        "lower-right");
    overlay_summary(fin, "h_respET_cb_sigma_%s_%s", "common_cut",
        kAllIdx, kNsets, 0.02, 0.16,
        "CB  #sigma(E_{T}^{reco} / #it{p}_{T}^{truth})",
        "common cut (CB fit)",
        (outDir + "/resp_sigma_cb_common_cut.pdf").c_str(),
        "upper-right");

    // -----------------------------------------------------------------
    // Pileup-only overlay (no vtx reweighting) at reco level
    // -----------------------------------------------------------------
    overlay_summary(fin, "h_respET_mean_%s_%s",  "common_cut",
        kPileupIdx, kNpileupIdx, 0.92, 1.03,
        "#LTE_{T}^{reco} / #it{p}_{T}^{truth}#GT",
        "common cut",
        (outDir + "/resp_mean_common_cut_pileup.pdf").c_str());
    overlay_summary(fin, "h_respET_sigma_%s_%s", "common_cut",
        kPileupIdx, kNpileupIdx, 0.02, 0.10,
        "#sigma(E_{T}^{reco} / #it{p}_{T}^{truth})",
        "common cut",
        (outDir + "/resp_sigma_common_cut_pileup.pdf").c_str());

    // -----------------------------------------------------------------
    // Vertex-reweighting overlay (0 mrad / 1.5 mrad, w/wo rw) at reco level
    // -----------------------------------------------------------------
    // Data slopes down from upper-left to lower-right; move legend to upper-right
    // where the panel is empty (previously lower-left overlapped the data curves).
    overlay_summary(fin, "h_respET_mean_%s_%s",  "common_cut",
        kVtxRwIdx, kNvtxRwIdx, 0.94, 1.02,
        "#LTE_{T}^{reco} / #it{p}_{T}^{truth}#GT",
        "effect of z_{vtx} reweighting",
        (outDir + "/resp_mean_common_cut_vtxrw.pdf").c_str(),
        "upper-right");
    // Data values for mixed_0mrad variants extend up to ~0.093 at low pT;
    // widen yhi so they are not clipped (previously clipped at 0.09, making
    // the 0 mrad curves disappear for pT < 20 GeV).
    overlay_summary(fin, "h_respET_sigma_%s_%s", "common_cut",
        kVtxRwIdx, kNvtxRwIdx, 0.03, 0.11,
        "#sigma(E_{T}^{reco} / #it{p}_{T}^{truth})",
        "effect of z_{vtx} reweighting",
        (outDir + "/resp_sigma_common_cut_vtxrw.pdf").c_str(),
        "lower-right");

    // -----------------------------------------------------------------
    // Response projections at 3 representative pT bins (matched level)
    //   pT bins (pT_bins_truth from config): 7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45
    //   low = bin 3 (10-12 GeV), mid = bin 9 (22-24), high = bin 12 (28-32)
    // -----------------------------------------------------------------
    overlay_projection(fin, "respET", "reco", kPileupIdx, kNpileupIdx, 3,
        "10 < #it{p}_{T}^{truth} < 12 GeV",
        "E_{T}^{reco} / #it{p}_{T}^{truth}",
        (outDir + "/resp_projection_reco_pt10_12.pdf").c_str());
    overlay_projection(fin, "respET", "reco", kPileupIdx, kNpileupIdx, 9,
        "22 < #it{p}_{T}^{truth} < 24 GeV",
        "E_{T}^{reco} / #it{p}_{T}^{truth}",
        (outDir + "/resp_projection_reco_pt22_24.pdf").c_str());
    overlay_projection(fin, "respET", "reco", kPileupIdx, kNpileupIdx, 12,
        "28 < #it{p}_{T}^{truth} < 32 GeV",
        "E_{T}^{reco} / #it{p}_{T}^{truth}",
        (outDir + "/resp_projection_reco_pt28_32.pdf").c_str());
    overlay_projection(fin, "respET", "common_cut",    kPileupIdx, kNpileupIdx, 9,
        "22 < #it{p}_{T}^{truth} < 24 GeV (common cut)",
        "E_{T}^{reco} / #it{p}_{T}^{truth}",
        (outDir + "/resp_projection_common_cut_pt22_24.pdf").c_str());

    // -----------------------------------------------------------------
    // Fit overlays: single-set projections with Gauss + CB fit curves
    // at 22-24 GeV (index 9). Shows how much of the low-side tail the
    // Gaussian misses and the CB absorbs — especially for double MC.
    // -----------------------------------------------------------------
    overlay_projection_with_fits(fin, "respET", 0 /*single*/, "reco", 9,
        "22 < #it{p}_{T}^{truth} < 24 GeV",
        "E_{T}^{reco} / #it{p}_{T}^{truth}",
        (outDir + "/resp_projection_fits_single_reco_pt22_24.pdf").c_str());
    overlay_projection_with_fits(fin, "respET", 1 /*double*/, "reco", 9,
        "22 < #it{p}_{T}^{truth} < 24 GeV",
        "E_{T}^{reco} / #it{p}_{T}^{truth}",
        (outDir + "/resp_projection_fits_double_reco_pt22_24.pdf").c_str());
    overlay_projection_with_fits(fin, "respET", 0 /*single*/, "common_cut", 9,
        "22 < #it{p}_{T}^{truth} < 24 GeV",
        "E_{T}^{reco} / #it{p}_{T}^{truth}",
        (outDir + "/resp_projection_fits_single_common_cut_pt22_24.pdf").c_str());
    overlay_projection_with_fits(fin, "respET", 1 /*double*/, "common_cut", 9,
        "22 < #it{p}_{T}^{truth} < 24 GeV",
        "E_{T}^{reco} / #it{p}_{T}^{truth}",
        (outDir + "/resp_projection_fits_double_common_cut_pt22_24.pdf").c_str());

    // -----------------------------------------------------------------
    // Fit QC grids for the report appendix: 12-panel projection + fits
    // across analysis pT range for single, double, and mixed_1p5mrad MC.
    // -----------------------------------------------------------------
    draw_fit_qc_grid(fin, "respET", 0 /*single*/,         "reco",
                     (outDir + "/resp_fit_qc_grid_single_reco.pdf").c_str());
    draw_fit_qc_grid(fin, "respET", 1 /*double*/,         "reco",
                     (outDir + "/resp_fit_qc_grid_double_reco.pdf").c_str());
    draw_fit_qc_grid(fin, "respET", 4 /*mixed_1p5mrad*/,  "reco",
                     (outDir + "/resp_fit_qc_grid_mixed_1p5mrad_reco.pdf").c_str());
    draw_fit_qc_grid(fin, "respET", 2 /*mixed_0mrad*/,    "reco",
                     (outDir + "/resp_fit_qc_grid_mixed_0mrad_reco.pdf").c_str());
    draw_fit_qc_grid(fin, "respET", 0 /*single*/,         "common_cut",
                     (outDir + "/resp_fit_qc_grid_single_common_cut.pdf").c_str());
    draw_fit_qc_grid(fin, "respET", 1 /*double*/,         "common_cut",
                     (outDir + "/resp_fit_qc_grid_double_common_cut.pdf").c_str());

    // -----------------------------------------------------------------
    // DSCB-overlay QC grids for tight_iso level: the residual R_ET
    // distribution after tight+iso is narrow and roughly symmetric, so a
    // one-sided CB misses the high-side tail from isolation leakage.
    // Green (DSCB) tracks both sides.
    // -----------------------------------------------------------------
    draw_fit_qc_grid_dscb(fin, "respET", 0 /*single*/,        "tight_iso",
                          (outDir + "/resp_fit_qc_grid_dscb_single_tight_iso.pdf").c_str());
    draw_fit_qc_grid_dscb(fin, "respET", 1 /*double*/,        "tight_iso",
                          (outDir + "/resp_fit_qc_grid_dscb_double_tight_iso.pdf").c_str());
    draw_fit_qc_grid_dscb(fin, "respET", 4 /*mixed_1p5mrad*/, "tight_iso",
                          (outDir + "/resp_fit_qc_grid_dscb_mixed_1p5mrad_tight_iso.pdf").c_str());
    draw_fit_qc_grid_dscb(fin, "respET", 2 /*mixed_0mrad*/,   "tight_iso",
                          (outDir + "/resp_fit_qc_grid_dscb_mixed_0mrad_tight_iso.pdf").c_str());

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
            fr->GetXaxis()->SetRangeUser(10, 40);
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
            // Append level tag onto strleg3 line position but on a new line
            // below the sPHENIX header block, in a smaller font and in the
            // bottom-left of the pad where there is no data or legend.
            l.SetTextSize(0.032);
            l.DrawLatex(0.22, 0.23, Form("level: %s", lvl_label));

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
            // Honour leg_loc so mean panels (data near top) can send the
            // legend to lower-right; sigma panels (data near bottom) keep
            // it upper-right. For upper-right, bump lx1 from 0.44 -> 0.52
            // so the legend's "Mix 1.5 mrad" row does not clip the
            // "reco (CB)" / "common cut (CB)" level label at x ~ 0.22-0.42.
            double lx1 = 0.52, ly1 = 0.74, lx2 = 0.93, ly2 = 0.93;
            if (std::string(leg_loc) == "lower-right") {
                lx1 = 0.44; ly1 = 0.18; lx2 = 0.93; ly2 = 0.37;
            }
            TLegend *leg = new TLegend(lx1, ly1, lx2, ly2);
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
        // All 4 panels use lower-right for the legend: mean panels have
        // data near top (NDC > 0.5); sigma panels have data that descends
        // with pT — at the high-pT end values are ≲ 0.11 (NDC < 0.45 at
        // range 0.02-0.22), so the 0.18-0.37 NDC band is below the data.
        drawPanel(1, "h_respET_mean_%s_%s",  "reco",       0.86, 1.02,
                  "#LTE_{T}^{reco} / #it{p}_{T}^{truth}#GT", "reco",       "lower-right");
        drawPanel(2, "h_respET_mean_%s_%s",  "common_cut", 0.88, 1.05,
                  "#LTE_{T}^{reco} / #it{p}_{T}^{truth}#GT", "common cut", "lower-right");
        drawPanel(3, "h_respET_sigma_%s_%s", "reco",       0.02, 0.22,
                  "#sigma(E_{T}^{reco} / #it{p}_{T}^{truth})", "reco",       "lower-right");
        drawPanel(4, "h_respET_sigma_%s_%s", "common_cut", 0.02, 0.16,
                  "#sigma(E_{T}^{reco} / #it{p}_{T}^{truth})", "common cut", "lower-right");

        c.SaveAs((outDir + "/resp_summary_4panel.pdf").c_str());
    }

    // -----------------------------------------------------------------
    // 4-panel CB summary: CB peak/sigma × reco/common_cut (all 6 sets)
    // CB absorbs the low-side tail; peaks cluster at unity more tightly
    // and sigmas drop to the core resolution.
    // -----------------------------------------------------------------
    {
        TCanvas c("c4cb", "", 1400, 1100);
        c.Divide(2, 2, 0.002, 0.002);

        auto drawPanel = [&](int ipad, const char *hfmt, const char *lvl,
                             double ylo, double yhi, const char *yt,
                             const char *lvl_label, const char *leg_loc)
        {
            c.cd(ipad);
            gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.17);
            TH1F *fr = new TH1F(Form("frcb_%d", ipad), "", 1, 7, 40);
            fr->SetXTitle("#it{p}_{T}^{#gamma, truth} [GeV]");
            fr->SetYTitle(yt);
            fr->GetYaxis()->SetRangeUser(ylo, yhi);
            fr->GetXaxis()->SetRangeUser(10, 40);
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
            // Append level tag onto strleg3 line position but on a new line
            // below the sPHENIX header block, in a smaller font and in the
            // bottom-left of the pad where there is no data or legend.
            l.SetTextSize(0.032);
            l.DrawLatex(0.22, 0.23, Form("level: %s", lvl_label));

            static const char *kShortLabels[6] = {
                "Single",           "Double",
                "Mix 0 mrad",       "Mix 0 mrad +vtxrw",
                "Mix 1.5 mrad",     "Mix 1.5 mrad +vtxrw",
            };
            // Honour leg_loc so mean panels (data near top) can send the
            // legend to lower-right; sigma panels (data near bottom) keep
            // it upper-right. For upper-right, bump lx1 from 0.44 -> 0.52
            // so the legend's "Mix 1.5 mrad" row does not clip the
            // "reco (CB)" / "common cut (CB)" level label at x ~ 0.22-0.42.
            double lx1 = 0.52, ly1 = 0.74, lx2 = 0.93, ly2 = 0.93;
            if (std::string(leg_loc) == "lower-right") {
                lx1 = 0.44; ly1 = 0.18; lx2 = 0.93; ly2 = 0.37;
            }
            TLegend *leg = new TLegend(lx1, ly1, lx2, ly2);
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

        drawPanel(1, "h_respET_cb_mean_%s_%s",  "reco",       0.90, 1.03,
                  "CB peak  #mu(E_{T}^{reco} / #it{p}_{T}^{truth})", "reco (CB)",       "lower-right");
        drawPanel(2, "h_respET_cb_mean_%s_%s",  "common_cut", 0.94, 1.03,
                  "CB peak  #mu(E_{T}^{reco} / #it{p}_{T}^{truth})", "common cut (CB)", "lower-right");
        drawPanel(3, "h_respET_cb_sigma_%s_%s", "reco",       0.02, 0.16,
                  "CB  #sigma(E_{T}^{reco} / #it{p}_{T}^{truth})", "reco (CB)",       "lower-right");
        drawPanel(4, "h_respET_cb_sigma_%s_%s", "common_cut", 0.02, 0.14,
                  "CB  #sigma(E_{T}^{reco} / #it{p}_{T}^{truth})", "common cut (CB)", "lower-right");

        c.SaveAs((outDir + "/resp_summary_4panel_cb.pdf").c_str());
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
            // Legend: lower-right (curve tails drop to ~0 outside |z|>60 cm,
            // leaving this corner empty; also clear of the panel title block
            // at upper-left).
            leg1 = new TLegend(0.60, 0.22, 0.88, 0.36);
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
            // Legend: lower-right (the x-axis x=0 peak is above, and the
            // upper portion holds the panel-title/sPHENIX block).  Narrow
            // and compact so the 0-mrad data/sim ratio (which rises into
            // the right half of the frame) does not clip the entries.
            leg2 = new TLegend(0.56, 0.62, 0.88, 0.76);
            leg2->SetFillStyle(0); leg2->SetBorderSize(0);
            leg2->SetTextFont(42); leg2->SetTextSize(0.034);
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
