// plot_paper_iso_template.C -- paper Fig.~\ref{fig:isoET}
//
// Lifts the h1D_iso_<tune>_<lo>_<hi>.pdf canvas (L1-164) from
// plotting/CONF_plots.C and pins it to the (lowbin=3, highbin=5) cell shown
// in the journal draft, which corresponds to the 14 < ETg < 20 GeV slice
// (ptRanges = [10,12,14,16,18,20,22,...]).
//
// The fit panel from CONF_plots.C (L165-228, h1D_iso_fit_*.pdf) is not used
// in the paper and is intentionally NOT carried over.
//
// Inputs:
//   data_histo_<tune>.root      -- per-pT TH1Ds h_tight_isoET_0_<i>, h_nontight_isoET_0_<i>
//   MC_efficiency_<tune>.root   -- per-pT TH1Ds h_tight_isoET_0_<i> (signal MC)
//
// Output: PPG12-Paper/figures/h1D_iso_nom_3_5.pdf  (filename retained for
// stability with main.tex; the underlying data is the bdt_nom analysis.)

#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <string>

#include "paper_style.h"

namespace
{
constexpr int kLowBin  = 3;   // ptRanges[3] = 16
constexpr int kHighBin = 5;   // ptRanges[6] = 22  -> slice 16-22 GeV
// NOTE: the paper figure file is named "..._3_5.pdf" but the inclusive sum
// runs from bin 3 through bin 5 inclusive => upper edge ptRanges[5+1].
}  // namespace

void plot_paper_iso_template(const std::string &tune = "bdt_nom")
{
    paper_init();

    const std::string resdir = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    TFile *fdata = TFile::Open(Form("%s/data_histo_%s.root",   resdir.c_str(), tune.c_str()));
    TFile *fmc   = TFile::Open(Form("%s/MC_efficiency_%s.root", resdir.c_str(), tune.c_str()));
    if (!fdata || fdata->IsZombie() || !fmc || fmc->IsZombie())
    {
        std::cerr << "[plot_paper_iso_template] missing data_histo/MC_efficiency for tune=" << tune << std::endl;
        return;
    }

    // Per-pT input histograms.
    TH1D *h_tight_isoET_pt[NptBins];
    TH1D *h_tight_isoET_mcSig_pt[NptBins];
    TH1D *h_nontight_isoET_pt[NptBins];
    for (int ipt = 0; ipt < NptBins; ++ipt)
    {
        h_tight_isoET_pt[ipt]       = (TH1D *)fdata->Get(Form("h_tight_isoET_0_%d",    ipt));
        h_tight_isoET_mcSig_pt[ipt] = (TH1D *)fmc  ->Get(Form("h_tight_isoET_0_%d",    ipt));
        h_nontight_isoET_pt[ipt]    = (TH1D *)fdata->Get(Form("h_nontight_isoET_0_%d", ipt));
    }

    // Sum the requested pT slice (lo .. hi inclusive) into a single distribution.
    TH1D *h_tight_isoET       = (TH1D *)h_tight_isoET_pt      [kLowBin]->Clone("h_tight_isoET");
    TH1D *h_tight_isoET_mcSig = (TH1D *)h_tight_isoET_mcSig_pt[kLowBin]->Clone("h_tight_isoET_mcSig");
    TH1D *h_nontight_isoET    = (TH1D *)h_nontight_isoET_pt   [kLowBin]->Clone("h_nontight_isoET_pt");
    for (int ib = kLowBin + 1; ib <= kHighBin; ++ib)
    {
        h_tight_isoET      ->Add(h_tight_isoET_pt[ib]);
        h_tight_isoET_mcSig->Add(h_tight_isoET_mcSig_pt[ib]);
        h_nontight_isoET   ->Add(h_nontight_isoET_pt[ib]);
    }

    // Variable rebinning: 1-GeV bins below 2.5 GeV, 5-GeV bins above.
    int nbinsOrig = h_tight_isoET->GetNbinsX();
    std::vector<double> newEdges;
    newEdges.push_back(h_tight_isoET->GetBinLowEdge(1));
    int i = 1;
    while (i <= nbinsOrig)
    {
        int groupSize = (h_tight_isoET->GetBinLowEdge(i) < 2.5) ? 1 : 5;
        int last = i + groupSize - 1;
        if (last > nbinsOrig) last = nbinsOrig;
        double upperEdge = h_tight_isoET->GetBinLowEdge(last) + h_tight_isoET->GetBinWidth(last);
        newEdges.push_back(upperEdge);
        i += groupSize;
    }
    int Nnew = newEdges.size() - 1;
    h_tight_isoET       = (TH1D *)h_tight_isoET      ->Rebin(Nnew, "h_tight_isoET_rebin",       newEdges.data());
    h_tight_isoET_mcSig = (TH1D *)h_tight_isoET_mcSig->Rebin(Nnew, "h_tight_isoET_mcSig_rebin", newEdges.data());
    h_nontight_isoET    = (TH1D *)h_nontight_isoET   ->Rebin(Nnew, "h_nontight_isoET_rebin",    newEdges.data());
    h_tight_isoET      ->Scale(1., "width");
    h_tight_isoET_mcSig->Scale(1., "width");
    h_nontight_isoET   ->Scale(1., "width");

    // Normalise non-isolated background and signal MC into the data shape.
    const float ptcutnorm = 6;
    int   ptB         = h_tight_isoET->FindBin(ptcutnorm);
    int   ptBs        = h_tight_isoET->GetNbinsX();
    float normTight   = h_tight_isoET   ->Integral(ptB, ptBs);
    float normNontight= h_nontight_isoET->Integral(ptB, ptBs);
    h_nontight_isoET->Scale(normTight / normNontight);
    float dataSig = h_tight_isoET->Integral() - h_nontight_isoET->Integral();
    h_tight_isoET_mcSig->Scale(dataSig / h_tight_isoET_mcSig->Integral());

    TCanvas *c1 = new TCanvas("c1_paper_isoET", "", 600, 560);
    h_tight_isoET->Draw("ex0");
    h_tight_isoET->GetXaxis()->SetRangeUser(-1, 15);

    h_tight_isoET_mcSig->SetFillColorAlpha(kBlue, 0.3);
    h_tight_isoET_mcSig->SetFillStyle(1001);
    h_tight_isoET_mcSig->SetLineColor(kBlue - 1);
    h_tight_isoET_mcSig->SetLineWidth(2);

    h_nontight_isoET->SetFillColorAlpha(kRed, 0.3);
    h_nontight_isoET->SetFillStyle(1001);
    h_nontight_isoET->SetLineColor(kRed - 1);
    h_nontight_isoET->SetLineWidth(2);

    THStack *hs = new THStack("hs_paper_isoET", "");
    hs->Add(h_nontight_isoET);
    hs->Add(h_tight_isoET_mcSig);
    hs->Draw("hist same");
    h_tight_isoET->Draw("same ex0");
    h_tight_isoET->Draw("same axis");

    h_tight_isoET->SetYTitle("Counts / Bin Width");
    h_tight_isoET->SetXTitle("#it{E}_{T}^{iso} [GeV]");
    h_tight_isoET->GetXaxis()->SetTitleOffset(1.2);

    std::string st_etbin = Form("%0.0f < #kern[-0.25]{#it{E}_{T}^{#gamma}} < %0.0f GeV",
                                ptRanges[kLowBin], ptRanges[kHighBin + 1]);
    float xpos(0.915), ypos(0.875), dy(0.054), dy1(0.06), fontsize(0.046), fontsize1(0.048);
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),     fontsize1, 1);
    myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(),     fontsize,  1);
    myText(xpos, ypos - 2 * dy, 1, st_etbin.c_str(),    fontsize,  1);
    myText(xpos, ypos - 3 * dy, 1, strleg3.c_str(),     fontsize,  1);

    int nLegend(3), legOffset(3);
    TLegend *l1 = new TLegend(0.51, ypos - (nLegend + legOffset) * dy1, xpos, ypos - legOffset * dy + 0.03);
    legStyle(l1, 0.20, fontsize);
    l1->AddEntry(h_tight_isoET,       "Data (Signal)",     "pe");
    l1->AddEntry(h_nontight_isoET,    "Data (Background)", "f");
    l1->AddEntry(h_tight_isoET_mcSig, "Signal MC",         "f");
    l1->Draw("same");

    // Filename retained as h1D_iso_nom_<lo>_<hi>.pdf for main.tex stability;
    // see header comment.
    c1->SaveAs(Form("%s/h1D_iso_nom_%d_%d.pdf",
                    paper_savepath().c_str(), kLowBin, kHighBin));
}
