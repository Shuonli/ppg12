// Closure check: weighted MC reco-vertex distribution (cross-section x
// lumi_weight x truth-vertex-reweight applied per event, hadd'd across the
// two crossing-angle periods) vs data reco-vertex distribution, both after
// the nominal cross-section pipeline.
//
// Inputs: efficiencytool/results/{data_histo,MC_efficiency}_bdt_nom*.root
// (h_vertexz; 200 bins, [-100,100] cm).
//
// Output: plotting/figures/truth_vertex_closure_nominal.pdf

#include "plotcommon.h"

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TLine.h"
#include "TMath.h"

#include <iostream>

namespace {
constexpr const char *kData0rad   = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom_0rad.root";
constexpr const char *kData1p5rad = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom_1p5mrad.root";
constexpr const char *kDataAll    = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom.root";
constexpr const char *kMC0rad     = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom_0rad.root";
constexpr const char *kMC1p5rad   = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom_1p5mrad.root";
constexpr const char *kMCAll      = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root";

TH1F *grab(const char *fn, const char *hn, const char *newname) {
  TFile *f = TFile::Open(fn, "READ");
  if (!f || f->IsZombie()) { std::cerr << "cannot open " << fn << std::endl; return nullptr; }
  TH1F *h = dynamic_cast<TH1F *>(f->Get(hn));
  if (!h) { std::cerr << hn << " not found in " << fn << std::endl; f->Close(); return nullptr; }
  TH1F *out = (TH1F *)h->Clone(newname);
  out->SetDirectory(0);
  f->Close();
  return out;
}

void stylePanel(TH1F *h, int colour, int marker) {
  h->SetLineColor(colour);
  h->SetMarkerColor(colour);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(0.6);
  h->SetLineWidth(2);
}

// chi2 / ndf and KS on the SHAPE (both unit-integral normalised)
// Limit to the fiducial |z| < 60 cm window (analysis vertex cut).
void computeStats(TH1F *hData, TH1F *hMC, double zlo, double zhi,
                  double &chi2ndf, double &kolS, int &ndf) {
  int bLo = hData->FindBin(zlo + 1e-6);
  int bHi = hData->FindBin(zhi - 1e-6);
  double chi2 = 0;
  ndf = 0;
  for (int b = bLo; b <= bHi; ++b) {
    double d = hData->GetBinContent(b);
    double m = hMC->GetBinContent(b);
    double ed = hData->GetBinError(b);
    double em = hMC->GetBinError(b);
    double var = ed * ed + em * em;
    if (var <= 0) continue;
    chi2 += (d - m) * (d - m) / var;
    ++ndf;
  }
  chi2ndf = (ndf > 0) ? chi2 / ndf : -1;
  kolS = hData->KolmogorovTest(hMC, "MN");
}

void drawPanel(TPad *padTop, TPad *padBot,
               TH1F *hData, TH1F *hMC, const char *periodLabel,
               double chi2ndf, double kolS, int ndf) {
  // normalise to unit integral in the |z| < 60 cm window
  const double zlo = -60, zhi = 60;
  auto integrate = [&](TH1F *h) {
    int bLo = h->FindBin(zlo + 1e-6);
    int bHi = h->FindBin(zhi - 1e-6);
    return h->Integral(bLo, bHi);
  };
  double idata = integrate(hData);
  double imc   = integrate(hMC);
  if (idata > 0) hData->Scale(1.0 / idata);
  if (imc   > 0) hMC->Scale(1.0 / imc);

  stylePanel(hData, kBlack,       20);
  stylePanel(hMC,   kAzure + 2,    0);  // line only
  hMC->SetMarkerSize(0);
  hMC->SetLineStyle(1);

  // Top pad: overlay
  padTop->cd();
  padTop->SetTopMargin(0.08);
  padTop->SetBottomMargin(0.02);
  padTop->SetLeftMargin(0.16);
  padTop->SetRightMargin(0.04);

  double yMax = std::max(hData->GetMaximum(), hMC->GetMaximum()) * 1.5;
  hData->GetXaxis()->SetRangeUser(-80, 80);
  hData->GetYaxis()->SetRangeUser(0, yMax);
  hData->GetYaxis()->SetTitle("normalised / 1 cm");
  hData->GetYaxis()->SetTitleSize(0.06);
  hData->GetYaxis()->SetTitleOffset(1.2);
  hData->GetYaxis()->SetLabelSize(0.05);
  hData->GetXaxis()->SetLabelSize(0);
  hData->GetXaxis()->SetTitleSize(0);
  hData->SetTitle("");
  hData->Draw("E");
  hMC->Draw("HIST SAME");

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.055);
  lat.DrawLatex(0.20, 0.86, strleg1.c_str());
  lat.SetTextSize(0.045);
  lat.DrawLatex(0.20, 0.80, strleg2_1.c_str());
  lat.DrawLatex(0.20, 0.74, periodLabel);
  lat.DrawLatex(0.20, 0.68, Form("#chi^{2}/ndf = %.2f  (ndf=%d)", chi2ndf, ndf));
  lat.DrawLatex(0.20, 0.62, Form("KS (shape) = %.3f", kolS));

  TLegend *leg = new TLegend(0.62, 0.72, 0.92, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->AddEntry(hData, "Data", "pe");
  leg->AddEntry(hMC,   "MC (weighted)", "l");
  leg->Draw();

  // Bottom pad: ratio
  padBot->cd();
  padBot->SetTopMargin(0.02);
  padBot->SetBottomMargin(0.30);
  padBot->SetLeftMargin(0.16);
  padBot->SetRightMargin(0.04);

  TH1F *hRatio = (TH1F *)hData->Clone(Form("ratio_%s", periodLabel));
  hRatio->Divide(hMC);
  hRatio->SetTitle("");
  hRatio->GetYaxis()->SetTitle("Data / MC");
  hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetTitleSize(0.12);
  hRatio->GetYaxis()->SetTitleOffset(0.55);
  hRatio->GetYaxis()->SetLabelSize(0.10);
  hRatio->GetXaxis()->SetTitle("reco vertex z [cm]");
  hRatio->GetXaxis()->SetTitleSize(0.13);
  hRatio->GetXaxis()->SetTitleOffset(1.0);
  hRatio->GetXaxis()->SetLabelSize(0.10);
  hRatio->GetXaxis()->SetRangeUser(-80, 80);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(0.6);
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerColor(kBlack);
  hRatio->Draw("E");

  TLine *l1 = new TLine(-80, 1, 80, 1);
  l1->SetLineStyle(2);
  l1->SetLineColor(kGray + 2);
  l1->Draw();
  TLine *lzLo = new TLine(-60, 0.5, -60, 1.5);
  TLine *lzHi = new TLine( 60, 0.5,  60, 1.5);
  for (TLine *L : {lzLo, lzHi}) {
    L->SetLineStyle(3);
    L->SetLineColor(kGray + 2);
    L->Draw();
  }
}
} // anonymous namespace

void plot_truth_vertex_closure()
{
  init_plot();

  // --- grab the six histograms -----------------------------------------
  TH1F *hData_0   = grab(kData0rad,   "h_vertexz", "hData_0rad");
  TH1F *hData_1p5 = grab(kData1p5rad, "h_vertexz", "hData_1p5mrad");
  TH1F *hData_all = grab(kDataAll,    "h_vertexz", "hData_all");
  TH1F *hMC_0     = grab(kMC0rad,     "h_vertexz", "hMC_0rad");
  TH1F *hMC_1p5   = grab(kMC1p5rad,   "h_vertexz", "hMC_1p5mrad");
  TH1F *hMC_all   = grab(kMCAll,      "h_vertexz", "hMC_all");

  if (!hData_0 || !hData_1p5 || !hData_all ||
      !hMC_0   || !hMC_1p5   || !hMC_all) {
    std::cerr << "Missing input histograms; aborting." << std::endl;
    return;
  }

  // Data-vs-MC statistics on SHAPE (unit-normalised in |z| < 60 cm).
  // Do stats first on *copies*, because drawPanel rescales the inputs.
  auto copy = [](TH1F *h, const char *n) {
    TH1F *c = (TH1F *)h->Clone(n);
    c->SetDirectory(0);
    return c;
  };
  TH1F *cD0   = copy(hData_0,   "cD0");
  TH1F *cM0   = copy(hMC_0,     "cM0");
  TH1F *cD1p5 = copy(hData_1p5, "cD1p5");
  TH1F *cM1p5 = copy(hMC_1p5,   "cM1p5");
  TH1F *cDA   = copy(hData_all, "cDA");
  TH1F *cMA   = copy(hMC_all,   "cMA");

  // Normalise the copies too so the chi2 is on shape.
  const double zlo = -60, zhi = 60;
  auto normShape = [&](TH1F *h) {
    int bLo = h->FindBin(zlo + 1e-6);
    int bHi = h->FindBin(zhi - 1e-6);
    double I = h->Integral(bLo, bHi);
    if (I > 0) h->Scale(1.0 / I);
  };
  for (TH1F *h : {cD0, cM0, cD1p5, cM1p5, cDA, cMA}) normShape(h);

  double c2_0, ks_0, c2_1p5, ks_1p5, c2_A, ks_A;
  int    n_0, n_1p5, n_A;
  computeStats(cD0,   cM0,   zlo, zhi, c2_0,   ks_0,   n_0);
  computeStats(cD1p5, cM1p5, zlo, zhi, c2_1p5, ks_1p5, n_1p5);
  computeStats(cDA,   cMA,   zlo, zhi, c2_A,   ks_A,   n_A);

  std::cout << "=== truth-vertex reweight closure (shape, |z|<60 cm) ===" << std::endl;
  std::cout << "0 mrad:   chi2/ndf = " << c2_0   << "  (ndf=" << n_0   << "),  KS = " << ks_0   << std::endl;
  std::cout << "1.5 mrad: chi2/ndf = " << c2_1p5 << "  (ndf=" << n_1p5 << "),  KS = " << ks_1p5 << std::endl;
  std::cout << "all:      chi2/ndf = " << c2_A   << "  (ndf=" << n_A   << "),  KS = " << ks_A   << std::endl;

  // --- 3-panel canvas ---------------------------------------------------
  TCanvas *c = new TCanvas("cVtxClosure", "", 1800, 700);
  c->Divide(3, 1, 0.002, 0.002);

  auto buildSubpanel = [&](int idx, TH1F *hD, TH1F *hM, const char *lbl,
                           double c2, double ks, int ndf) {
    c->cd(idx);
    TPad *pParent = (TPad *)gPad;
    pParent->SetLeftMargin(0);
    pParent->SetRightMargin(0);
    pParent->SetTopMargin(0);
    pParent->SetBottomMargin(0);

    TPad *pTop = new TPad(Form("pTop_%d", idx), "", 0, 0.32, 1, 1);
    TPad *pBot = new TPad(Form("pBot_%d", idx), "", 0, 0.00, 1, 0.32);
    pTop->Draw();
    pBot->Draw();
    drawPanel(pTop, pBot, hD, hM, lbl, c2, ks, ndf);
  };

  buildSubpanel(1, hData_0,   hMC_0,   "0 mrad  (47289-51274)",   c2_0,   ks_0,   n_0);
  buildSubpanel(2, hData_1p5, hMC_1p5, "1.5 mrad  (51274-54000)", c2_1p5, ks_1p5, n_1p5);
  buildSubpanel(3, hData_all, hMC_all, "all runs (hadd lumi-weighted)", c2_A, ks_A, n_A);

  c->SaveAs("/sphenix/user/shuhangli/ppg12/plotting/figures/truth_vertex_closure_nominal.pdf");

  std::cout << "Saved: /sphenix/user/shuhangli/ppg12/plotting/figures/truth_vertex_closure_nominal.pdf" << std::endl;
}
