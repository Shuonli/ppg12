// plot_paper_purity.C -- paper Fig.~\ref{fig:photon_purity}
//
// Lifts the data-only purity panel from plotting/plot_purity_selection.C
// (the plotMC_truth = false branch). Drops the inclusive-MC truth overlay
// and the leak-vs-no-leak split so the paper figure shows only the leak-
// corrected purity with its 68% CI band and parametric fit. The "no-leak"
// markers are kept so the reader can see the size of the signal-leakage
// correction.
//
// Inputs:
//   Photon_final_<tune>.root      -- gpurity, gpurity_leak, grFineConf_leak,
//                                    f_purity_leak_fit
//
// Output: PPG12-Paper/figures/purity_<tune>.pdf  (canonical name; main.tex
// includes purity_nom.pdf -- a symlink/copy is made by the driver script.)

#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <string>

#include "paper_style.h"

void plot_paper_purity(const std::string &tune = "bdt_nom")
{
    paper_init();

    const std::string resdir = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    TFile *fdata = TFile::Open(Form("%s/Photon_final_%s.root", resdir.c_str(), tune.c_str()));
    if (!fdata || fdata->IsZombie())
    {
        std::cerr << "[plot_paper_purity] could not open Photon_final_" << tune << ".root" << std::endl;
        return;
    }

    TGraphErrors *gpurity         = (TGraphErrors *)fdata->Get("gpurity");
    TGraphErrors *gpurity_leak    = (TGraphErrors *)fdata->Get("gpurity_leak");
    TGraphErrors *grFineConf_leak = (TGraphErrors *)fdata->Get("grFineConf_leak");
    TF1          *f_purity_fit    = (TF1 *)         fdata->Get("f_purity_leak_fit");

    if (!gpurity || !gpurity_leak || !grFineConf_leak || !f_purity_fit)
    {
        std::cerr << "[plot_paper_purity] one or more graphs missing in Photon_final_" << tune << ".root" << std::endl;
        return;
    }

    TCanvas *c1 = new TCanvas("c1_paper_purity", "", 600, 600);
    frame_et_rec->SetTitle(";#it{E}_{T}^{#gamma} [GeV];Purity");
    frame_et_rec->GetYaxis()->SetRangeUser(0.0, 1.1);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 35);
    frame_et_rec->Draw("axis");

    gpurity->SetMarkerColor(kBlack);
    gpurity->SetMarkerStyle(20);
    gpurity->SetMarkerSize(1.5);
    gpurity->SetLineColor(kBlack);
    gpurity->Draw("P same");

    gpurity_leak->SetMarkerColor(kBlue);
    gpurity_leak->SetMarkerStyle(24);
    gpurity_leak->SetMarkerSize(1.5);
    gpurity_leak->SetLineColor(kBlue);
    gpurity_leak->Draw("P same");

    grFineConf_leak->SetMarkerColor(kAzure + 2);
    grFineConf_leak->SetLineColor(kAzure + 2);
    grFineConf_leak->SetFillColorAlpha(kAzure + 2, 0.2);
    grFineConf_leak->Draw("e3 same");

    f_purity_fit->SetLineColor(kAzure + 2);
    f_purity_fit->Draw("same");

    float xpos(0.2), xpos2(0.915), ypos(0.88), ypos2(0.19), dy(0.056), dy1(0.075), fontsize(0.048);
    myText(xpos,  ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
    myText(xpos,  ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);
    myText(xpos2, ypos - 0 * dy, 1, strleg3.c_str(), fontsize, 1);

    int nEntry = 3;
    TLegend *l1 = new TLegend(xpos, ypos2, xpos2, ypos2 + nEntry * dy1);
    legStyle(l1, 0.14, fontsize);
    l1->AddEntry(gpurity_leak,    "w/ sig. leak. corr.",  "pl");
    l1->AddEntry(gpurity,         "w/o sig. leak. corr.", "pl");
    l1->AddEntry(grFineConf_leak, "fit w/ 68% C.L.",      "fl");
    l1->Draw("same");

    // Canonical paper name (no tune suffix) for main.tex stability.
    c1->SaveAs(Form("%s/purity_nom.pdf", paper_savepath().c_str()));
}
