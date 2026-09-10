// plot_eff_compare_herwig.C
//
// Compare iso, id, and MBD-vertex efficiencies for the PYTHIA-vs-HERWIG
// 1.5 mrad cross-check. Each output is a 2-pad plot (top: efficiency curves,
// bottom: HERWIG / PYTHIA ratio) in the main-analysis plotting style.
//
// Style follows plot_efficiency.C + plot_unfolded_woeff_compare.C:
//   - frame_et_truth from plotcommon.h
//   - PYTHIA = kAzure+2 / marker 20, HERWIG = kSpring-7 / marker 22
//   - Header strings strleg1..4 + strMC
//
// Inputs: MC_efficiency_*.root (eff_iso_eta_0, eff_id_eta_0)
//         Photon_final_*.root  (g_mbd_eff)
//
// Outputs:
//   figures/eff_iso_compare_1p5mrad.pdf
//   figures/eff_id_compare_1p5mrad.pdf
//   figures/eff_mbd_compare_1p5mrad.pdf

#include "plotcommon.h"

namespace
{

// Convert a TEfficiency into a TGraphAsymmErrors (central value + asymmetric
// Clopper-Pearson errors). Skips bins outside [pT_lo, pT_hi].
TGraphAsymmErrors *toGraph(TEfficiency *eff, const char *name,
                           double pT_lo = 8., double pT_hi = 36.)
{
    const TH1 *htot = eff->GetTotalHistogram();
    int n = htot->GetNbinsX();
    TGraphAsymmErrors *g = new TGraphAsymmErrors();
    g->SetName(name);
    int k = 0;
    for (int i = 1; i <= n; ++i) {
        double x = htot->GetBinCenter(i);
        if (x < pT_lo || x > pT_hi) continue;
        double y = eff->GetEfficiency(i);
        double el = eff->GetEfficiencyErrorLow(i);
        double eh = eff->GetEfficiencyErrorUp(i);
        g->SetPoint(k, x, y);
        g->SetPointError(k, htot->GetBinWidth(i)/2., htot->GetBinWidth(i)/2., el, eh);
        ++k;
    }
    return g;
}

// Bin-by-bin ratio gH / gP. Assumes the two graphs share the same x-binning;
// matches points by x within tol. Errors combine in quadrature (independent
// MC samples).
TGraphAsymmErrors *makeRatio(TGraphAsymmErrors *gH, TGraphAsymmErrors *gP,
                             const char *name, double tol = 1e-3)
{
    TGraphAsymmErrors *r = new TGraphAsymmErrors();
    r->SetName(name);
    int kr = 0;
    for (int i = 0; i < gH->GetN(); ++i) {
        double xh = gH->GetX()[i];
        double yh = gH->GetY()[i];
        double ehl = gH->GetErrorYlow(i);
        double ehu = gH->GetErrorYhigh(i);
        // find matching point in gP
        int jmatch = -1;
        for (int j = 0; j < gP->GetN(); ++j) {
            if (std::fabs(gP->GetX()[j] - xh) < tol) { jmatch = j; break; }
        }
        if (jmatch < 0) continue;
        double yp = gP->GetY()[jmatch];
        double epl = gP->GetErrorYlow(jmatch);
        double epu = gP->GetErrorYhigh(jmatch);
        if (yp <= 0 || yh <= 0) continue;
        double ratio = yh / yp;
        // average asymmetric errors -> symmetrize for ratio (small errors)
        double rel_h = 0.5 * (ehl + ehu) / yh;
        double rel_p = 0.5 * (epl + epu) / yp;
        double rel_r = std::sqrt(rel_h*rel_h + rel_p*rel_p);
        r->SetPoint(kr, xh, ratio);
        double exl = gH->GetErrorXlow(i);
        double exh = gH->GetErrorXhigh(i);
        r->SetPointError(kr, exl, exh, rel_r * ratio, rel_r * ratio);
        ++kr;
    }
    return r;
}

// Two-pad efficiency vs ratio plot.
void drawTwoPad(TGraphAsymmErrors *gP, TGraphAsymmErrors *gH,
                const std::string &ytitle_top, double ymin_top, double ymax_top,
                const std::string &out_pdf,
                int col_p, int col_h, int mk_p, int mk_h,
                const std::string &leg_label = "")
{
    TGraphAsymmErrors *r = makeRatio(gH, gP, "r_HoverP");

    TCanvas *c = new TCanvas(("can_" + out_pdf).c_str(), "", 800, 900);
    c->Divide(1, 2);

    // -- top pad
    TPad *pad_1 = (TPad *)c->cd(1);
    pad_1->SetPad(0, 0.36, 1, 1);
    pad_1->SetTopMargin(0.06);
    pad_1->SetLeftMargin(0.13);
    pad_1->SetBottomMargin(0.002);
    pad_1->SetRightMargin(0.05);

    frame_et_truth->SetYTitle(ytitle_top.c_str());
    frame_et_truth->SetXTitle("");
    frame_et_truth->GetXaxis()->SetRangeUser(10, 36);
    frame_et_truth->GetYaxis()->SetRangeUser(ymin_top, ymax_top);
    frame_et_truth->GetXaxis()->SetTitleOffset(1.05);
    frame_et_truth->GetYaxis()->SetTitleOffset(1.10);
    frame_et_truth->GetYaxis()->SetTitleSize(0.053);
    frame_et_truth->GetXaxis()->SetLabelSize(0.050);
    frame_et_truth->GetYaxis()->SetLabelSize(0.050);
    frame_et_truth->GetXaxis()->SetLabelOffset(2);
    frame_et_truth->GetXaxis()->SetNdivisions(505);
    frame_et_truth->Draw("axis");

    gP->SetMarkerStyle(mk_p);
    gP->SetMarkerSize(1.4);
    gP->SetMarkerColor(col_p);
    gP->SetLineColor(col_p);
    gP->SetLineWidth(2);
    gP->Draw("p same");

    gH->SetMarkerStyle(mk_h);
    gH->SetMarkerSize(1.4);
    gH->SetMarkerColor(col_h);
    gH->SetLineColor(col_h);
    gH->SetLineWidth(2);
    gH->Draw("p same");

    // header text
    myText(0.50, 0.88, 1, strleg1.c_str(), 0.05);
    myText(0.50, 0.82, 1, "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, 17.16 pb^{-1}  (1.5 mrad)", 0.04);
    myText(0.50, 0.77, 1, strleg3.c_str(), 0.04);
    myText(0.50, 0.72, 1, strleg4.c_str(), 0.04);
    if (!leg_label.empty()) myText(0.50, 0.65, 1, leg_label.c_str(), 0.04);

    myMarkerLineText(0.18, 0.30, 1, col_p, mk_p, col_p, 1, "PYTHIA", 0.040, true);
    myMarkerLineText(0.18, 0.24, 1, col_h, mk_h, col_h, 1, "HERWIG", 0.040, true);

    // -- bottom pad
    TPad *pad_2 = (TPad *)c->cd(2);
    pad_2->SetPad(0, 0, 1, 0.36);
    pad_2->SetTopMargin(0.02);
    pad_2->SetLeftMargin(0.13);
    pad_2->SetBottomMargin(0.30);
    pad_2->SetRightMargin(0.05);

    // need a separate frame for bottom pad — clone
    TH1F *frame_bot = (TH1F *)frame_et_truth->Clone("frame_bot");
    frame_bot->SetYTitle("HERWIG / PYTHIA");
    frame_bot->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
    frame_bot->GetYaxis()->SetNdivisions(506);
    frame_bot->GetYaxis()->SetRangeUser(0.85, 1.15);
    frame_bot->GetXaxis()->SetRangeUser(10, 36);
    frame_bot->GetXaxis()->SetTitleOffset(frame_et_truth->GetXaxis()->GetTitleOffset() * 4 / 6.* 1.4);
    frame_bot->GetYaxis()->SetTitleOffset(frame_et_truth->GetYaxis()->GetTitleOffset() * 4 / 6.);
    frame_bot->GetXaxis()->SetLabelSize(frame_et_truth->GetXaxis()->GetLabelSize() * 6 / 4.);
    frame_bot->GetYaxis()->SetLabelSize(frame_et_truth->GetYaxis()->GetLabelSize() * 6 / 4.);
    frame_bot->GetXaxis()->SetTitleSize(frame_et_truth->GetXaxis()->GetTitleSize() * 6 / 4. * 1.2);
    frame_bot->GetYaxis()->SetTitleSize(frame_et_truth->GetYaxis()->GetTitleSize() * 6 / 4.);
    frame_bot->GetXaxis()->SetLabelOffset(0.005);
    frame_bot->GetXaxis()->SetNdivisions(505);
    frame_bot->Draw("axis");

    r->SetMarkerStyle(mk_h);
    r->SetMarkerSize(1.4);
    r->SetMarkerColor(col_h);
    r->SetLineColor(col_h);
    r->SetLineWidth(2);
    r->Draw("p same");

    lineone->Draw("L");

    c->SaveAs(("figures/" + out_pdf).c_str());
    std::cout << "Saved figures/" << out_pdf << std::endl;
}

}  // namespace

void plot_eff_compare_herwig()
{
    init_plot();

    const std::string base = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/";

    // MC efficiency files for iso & id
    TFile *fmc_p = TFile::Open((base + "MC_efficiency_bdt_nom_1p5mrad_newrewt.root").c_str(), "READ");
    TFile *fmc_h = TFile::Open((base + "MC_efficiency_bdt_herwig_1p5mrad_common.root").c_str(), "READ");
    // Final yield files for MBD-vertex efficiency
    TFile *fdat_p = TFile::Open((base + "Photon_final_bdt_nom_1p5mrad_newrewt.root").c_str(), "READ");
    TFile *fdat_h = TFile::Open((base + "Photon_final_bdt_herwig_1p5mrad_common.root").c_str(), "READ");

    if (!fmc_p || fmc_p->IsZombie() || !fmc_h || fmc_h->IsZombie() ||
        !fdat_p || fdat_p->IsZombie() || !fdat_h || fdat_h->IsZombie()) {
        std::cerr << "ERROR: cannot open inputs" << std::endl;
        return;
    }

    const int col_p = kAzure + 2;
    const int col_h = kSpring - 7;
    const int mk_p  = 20;
    const int mk_h  = 22;

    // -- iso eff
    TEfficiency *eff_iso_p = (TEfficiency *)fmc_p->Get("eff_iso_eta_0");
    TEfficiency *eff_iso_h = (TEfficiency *)fmc_h->Get("eff_iso_eta_0");
    TGraphAsymmErrors *g_iso_p = toGraph(eff_iso_p, "g_iso_p");
    TGraphAsymmErrors *g_iso_h = toGraph(eff_iso_h, "g_iso_h");
    drawTwoPad(g_iso_p, g_iso_h,
               "Isolation efficiency #varepsilon_{iso}",
               0.85, 1.05,
               "eff_iso_compare_1p5mrad.pdf",
               col_p, col_h, mk_p, mk_h);

    // -- id eff
    TEfficiency *eff_id_p = (TEfficiency *)fmc_p->Get("eff_id_eta_0");
    TEfficiency *eff_id_h = (TEfficiency *)fmc_h->Get("eff_id_eta_0");
    TGraphAsymmErrors *g_id_p = toGraph(eff_id_p, "g_id_p");
    TGraphAsymmErrors *g_id_h = toGraph(eff_id_h, "g_id_h");
    drawTwoPad(g_id_p, g_id_h,
               "Identification efficiency #varepsilon_{ID}",
               0.45, 1.05,
               "eff_id_compare_1p5mrad.pdf",
               col_p, col_h, mk_p, mk_h);

    // -- MBD eff
    TGraphAsymmErrors *g_mbd_p = (TGraphAsymmErrors *)fdat_p->Get("g_mbd_eff");
    TGraphAsymmErrors *g_mbd_h = (TGraphAsymmErrors *)fdat_h->Get("g_mbd_eff");
    g_mbd_p = (TGraphAsymmErrors *)g_mbd_p->Clone("g_mbd_p");
    g_mbd_h = (TGraphAsymmErrors *)g_mbd_h->Clone("g_mbd_h");
    drawTwoPad(g_mbd_p, g_mbd_h,
               "MBD-vertex efficiency #varepsilon_{MBD}",
               0.45, 1.05,
               "eff_mbd_compare_1p5mrad.pdf",
               col_p, col_h, mk_p, mk_h);
}
