// VertexReweighting.C
// Goal:
//  - Read vertex z distributions (histogram name: h_vertexz) from:
//      data: results/data_histo_bdt_none.root
//      MC  : results/MC_efficiency_bdt_none.root
//      MCj : results/MC_efficiency_jet_bdt_none.root
//  - Normalize to unit area (shape-only)
//  - Produce ratio histograms Data/MC for future vertex reweighting
//  - Save outputs to a ROOT file (and optional PDF)
//
// Usage (in an environment with ROOT available):
//   root -l -b -q 'VertexReweighting.C("results/data_histo_bdt_none.root", \
//                                     "results/MC_efficiency_bdt_none.root", \
//                                     "results/MC_efficiency_jet_bdt_none.root", \
//                                     "results/vertex_reweight_bdt_none.root", \
//                                     1, 0, 0)'
//
// Args:
//   dataFile, mcFile, mcJetFile: input ROOT files
//   outFile: output ROOT file to write ratio hists into
//   normalize: 1 -> normalize each input hist to unit integral (recommended)
//   rebin: rebin factor (0 or 1 means no rebin)
//   makePdf: 1 -> also write results/vertex_reweighting.pdf
//

#include "TFile.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TClass.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include <iostream>
#include <string>

static TH1* FindHistRecursive(TDirectory* dir, const std::string& histName)
{
  if (!dir) return nullptr;

  TIter nextkey(dir->GetListOfKeys());
  while (TKey* key = (TKey*)nextkey()) {
    TObject* obj = key->ReadObj();
    if (!obj) continue;

    // Direct match
    if (obj->InheritsFrom(TH1::Class())) {
      if (histName == obj->GetName()) {
        return (TH1*)obj; // owned by file/dir; caller should Clone if needed
      }
    }

    // Recurse into subdirectories
    if (obj->InheritsFrom(TDirectory::Class())) {
      TH1* found = FindHistRecursive((TDirectory*)obj, histName);
      if (found) return found;
    }
  }

  return nullptr;
}

static TH1D* LoadVertexHist(const std::string& fileName,
                            const std::string& histName = "h_vertexz",
                            int rebin = 0,
                            int normalize = 1)
{
  TFile* f = TFile::Open(fileName.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[VertexReweighting] ERROR: cannot open file: " << fileName << std::endl;
    return nullptr;
  }

  TH1* h = nullptr;
  // First try top-level lookup
  h = dynamic_cast<TH1*>(f->Get(histName.c_str()));
  // Fallback: recursive search
  if (!h) h = FindHistRecursive(f, histName);

  if (!h) {
    std::cerr << "[VertexReweighting] ERROR: could not find histogram '" << histName
              << "' in file " << fileName << std::endl;
    f->Close();
    delete f;
    return nullptr;
  }

  TH1D* hc = dynamic_cast<TH1D*>(h->Clone((std::string(h->GetName()) + "_clone").c_str()));
  if (!hc) {
    // If original is TH1F/TH1 etc, clone into TH1D explicitly
    hc = new TH1D((std::string(h->GetName()) + "_clone").c_str(),
                  h->GetTitle(),
                  h->GetNbinsX(),
                  h->GetXaxis()->GetXmin(),
                  h->GetXaxis()->GetXmax());
    hc->Sumw2();
    for (int i = 0; i <= h->GetNbinsX() + 1; ++i) {
      hc->SetBinContent(i, h->GetBinContent(i));
      hc->SetBinError(i, h->GetBinError(i));
    }
  }
  hc->SetDirectory(nullptr);
  hc->Sumw2();

  if (rebin > 1) hc->Rebin(rebin);

  if (normalize) {
    const double integral = hc->Integral(0, hc->GetNbinsX() + 1);
    if (integral > 0) hc->Scale(1.0 / integral);
    else std::cerr << "[VertexReweighting] WARNING: zero integral for " << fileName << std::endl;
  }

  f->Close();
  delete f;
  return hc;
}

static TH1D* MakeRatioSafe(const TH1D* hData, const TH1D* hMc, const std::string& outName)
{
  if (!hData || !hMc) return nullptr;
  if (hData->GetNbinsX() != hMc->GetNbinsX()) {
    std::cerr << "[VertexReweighting] ERROR: bin mismatch for ratio " << outName
              << " (data bins " << hData->GetNbinsX() << ", mc bins " << hMc->GetNbinsX() << ")"
              << std::endl;
    return nullptr;
  }

  TH1D* hRatio = (TH1D*)hData->Clone(outName.c_str());
  hRatio->SetDirectory(nullptr);
  hRatio->Sumw2();

  // Use ROOT's error propagation for division, but guard MC=0 bins.
  for (int i = 0; i <= hRatio->GetNbinsX() + 1; ++i) {
    const double mc = hMc->GetBinContent(i);
    if (mc <= 0) {
      hRatio->SetBinContent(i, 0.0);
      hRatio->SetBinError(i, 0.0);
    } else {
      const double d = hData->GetBinContent(i);
      const double ed = hData->GetBinError(i);
      const double em = hMc->GetBinError(i);
      const double r = d / mc;
      // standard propagation: (d/m) * sqrt( (ed/d)^2 + (em/m)^2 ), protect d=0
      double er = 0.0;
      if (d > 0) er = r * std::sqrt((ed / d) * (ed / d) + (em / mc) * (em / mc));
      else er = r * (em / mc);
      hRatio->SetBinContent(i, r);
      hRatio->SetBinError(i, er);
    }
  }

  hRatio->SetTitle("Vertex reweighting: Data/MC;vtx z (cm);Data/MC");
  return hRatio;
}

void VertexReweighting(const char* dataFile = "results/data_histo_bdt_none.root",
                       const char* mcFile = "results/MC_efficiency_bdt_none.root",
                       const char* mcJetFile = "results/MC_efficiency_jet_bdt_none.root",
                       const char* outFile = "results/vertex_reweight_bdt_none.root",
                       int normalize = 1,
                       int rebin = 0,
                       int makePdf = 1)
{
  gStyle->SetOptStat(0);

  // Load and normalize vertex-z hists
  TH1D* hData = LoadVertexHist(dataFile, "h_vertexz", rebin, normalize);
  TH1D* hMc = LoadVertexHist(mcFile, "h_vertexz", rebin, normalize);
  TH1D* hMcJet = LoadVertexHist(mcJetFile, "h_vertexz", rebin, normalize);
  if (!hData || !hMc || !hMcJet) {
    std::cerr << "[VertexReweighting] ERROR: failed to load one or more input hists" << std::endl;
    return;
  }

  // Build a combined MC shape (simple sum of normalized shapes then renormalize)
  TH1D* hMcCombined = (TH1D*)hMc->Clone("h_vertexz_mc_combined");
  hMcCombined->SetDirectory(nullptr);
  hMcCombined->Add(hMcJet);
  const double combInt = hMcCombined->Integral(0, hMcCombined->GetNbinsX() + 1);
  if (combInt > 0) hMcCombined->Scale(1.0 / combInt);

  // Ratios
  TH1D* hRatioMc = MakeRatioSafe(hData, hMc, "h_vertexz_ratio_data_over_mc");
  TH1D* hRatioMcJet = MakeRatioSafe(hData, hMcJet, "h_vertexz_ratio_data_over_mcjet");
  TH1D* hRatioCombined = MakeRatioSafe(hData, hMcCombined, "h_vertexz_ratio_data_over_mccombined");
  if (!hRatioMc || !hRatioMcJet || !hRatioCombined) {
    std::cerr << "[VertexReweighting] ERROR: failed to compute ratio hist(s)" << std::endl;
    return;
  }

  // Save ROOT outputs
  TFile* fout = TFile::Open(outFile, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "[VertexReweighting] ERROR: cannot open output file: " << outFile << std::endl;
    return;
  }
  fout->cd();
  hData->SetName("h_vertexz_data");
  hMc->SetName("h_vertexz_mc");
  hMcJet->SetName("h_vertexz_mcjet");
  hData->Write();
  hMc->Write();
  hMcJet->Write();
  hMcCombined->Write();
  hRatioMc->Write();
  hRatioMcJet->Write();
  hRatioCombined->Write();
  fout->Close();
  delete fout;

  // Optional PDF quick-look
  if (makePdf) {
    TCanvas* c = new TCanvas("c_vertexz", "vertexz", 900, 900);
    c->Divide(1, 2);

    // Top: shapes
    c->cd(1);
    gPad->SetBottomMargin(0.12);
    hData->SetLineColor(kBlack);
    hData->SetMarkerColor(kBlack);
    hData->SetMarkerStyle(20);
    hData->SetTitle("Vertex z (normalized);vtx z (cm);Arb. norm.");

    hMc->SetLineColor(kRed + 1);
    hMcJet->SetLineColor(kBlue + 1);
    hMcCombined->SetLineColor(kGreen + 2);
    hMc->SetLineWidth(2);
    hMcJet->SetLineWidth(2);
    hMcCombined->SetLineWidth(2);

    double maxy = hData->GetMaximum();
    if (hMc->GetMaximum() > maxy) maxy = hMc->GetMaximum();
    if (hMcJet->GetMaximum() > maxy) maxy = hMcJet->GetMaximum();
    if (hMcCombined->GetMaximum() > maxy) maxy = hMcCombined->GetMaximum();
    hData->SetMaximum(1.25 * maxy);

    hData->Draw("E1");
    hMc->Draw("HIST SAME");
    hMcJet->Draw("HIST SAME");
    hMcCombined->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.60, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hData, "Data", "lep");
    leg->AddEntry(hMc, "MC (photon)", "l");
    leg->AddEntry(hMcJet, "MC (jet)", "l");
    leg->AddEntry(hMcCombined, "MC combined", "l");
    leg->Draw();

    // Bottom: ratio (combined)
    c->cd(2);
    gPad->SetTopMargin(0.06);
    gPad->SetBottomMargin(0.14);
    hRatioCombined->SetLineColor(kGreen + 2);
    hRatioCombined->SetMarkerColor(kGreen + 2);
    hRatioCombined->SetMarkerStyle(20);
    hRatioCombined->GetYaxis()->SetRangeUser(0.0, 2.0);
    hRatioCombined->Draw("E1");
    TLine* line = new TLine(hRatioCombined->GetXaxis()->GetXmin(), 1.0,
                            hRatioCombined->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->Draw("SAME");

    c->SaveAs("results/vertex_reweighting.pdf");
    delete c;
  }

  std::cout << "[VertexReweighting] Wrote: " << outFile << std::endl;
  if (makePdf) std::cout << "[VertexReweighting] Wrote: results/vertex_reweighting.pdf" << std::endl;
}











