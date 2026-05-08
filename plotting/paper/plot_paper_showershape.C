// plot_paper_showershape.C -- paper Fig.~\ref{fig:dis_showershape}
//
// Single-pad data/MC overlay for the journal draft. Cosmetics differ from
// the analysis-note version (plot_showershapes_selections.C):
//   * single-pad layout (no data-incMC diff bottom pad)
//   * no chi^2/ndf or p-value annotation
//   * eta label is |eta|<0.7 (no gamma superscript)
//   * bin label uses ET instead of pT, and "w/ preselection" not "w/ nbkg cut"
//
// Pinned to:
//   variables : weta_cogx, e32_to_e35
//   eta bin   : 0  (|eta| < 0.7, the analysis fiducial)
//   pT bin    : 1  (14 < ETg < 18 GeV; pT_bins = [10,14,18,22,28,30])
//   cut       : 1  (preselection applied)
//
// Outputs:
//   PPG12-Paper/figures/dis_weta_cogx_eta0_pt1_cut1.pdf
//   PPG12-Paper/figures/dis_e32_to_e35_eta0_pt1_cut1.pdf
//
// Inputs: DI-blended showershape ROOT files from the single-pass
// truth-vertex pipeline (signal_combined_<suffix>,
// jet_inclusive_combined_<suffix>). Default suffix = "showershape".

#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "paper_style.h"

namespace
{
const std::vector<std::string> kPaperVars = {"h2d_weta_cogx", "h2d_e32_to_e35"};
constexpr int kEtaBin = 0;
constexpr int kPtBin  = 1;     // pT_bins = [10,14,18,22,28,30]; pt1 = 14-18 GeV
constexpr int kCutBin = 1;     // 1 = preselection applied
constexpr float kPtLow  = 14.f;
constexpr float kPtHigh = 18.f;
const std::string kCutLabel = "w/ preselection";

// Paper-specific eta label: drop the gamma superscript present in
// plotcommon.h's strleg3 ("|#it{#eta^{#gamma}}| < 0.7").
const std::string kPaperEtaLabel = "|#it{#eta}| < 0.7";
}  // namespace

void plot_paper_showershape(const std::string &configsuffix = "showershape")
{
    paper_init();

    // DI-blended outputs (single-pass truth-vertex pipeline). When
    // combine_double=false the parent macro just clones the nominal SI
    // file, so we mirror that behaviour by reading the combined files
    // directly.
    const std::string resdir = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    TFile *f_data = TFile::Open(Form("%s/data_histoshower_shape_%s.root",                          resdir.c_str(), configsuffix.c_str()), "READ");
    TFile *f_sig  = TFile::Open(Form("%s/MC_efficiencyshower_shape_signal_combined_%s.root",       resdir.c_str(), configsuffix.c_str()), "READ");
    TFile *f_bkg  = TFile::Open(Form("%s/MC_efficiencyshower_shape_jet_inclusive_combined_%s.root", resdir.c_str(), configsuffix.c_str()), "READ");

    if (!f_data || f_data->IsZombie() ||
        !f_sig  || f_sig->IsZombie()  ||
        !f_bkg  || f_bkg->IsZombie())
    {
        std::cerr << "[plot_paper_showershape] missing one of the input files in " << resdir << std::endl;
        return;
    }

    auto scaleToUnit = [](TH1D *h)
    {
        if (!h) return;
        double integral = h->Integral();
        if (integral > 1e-12) h->Scale(1.0 / integral);
    };
    auto styleMC = [](TH1 *h, Color_t col)
    {
        h->SetLineColor(col); h->SetMarkerColor(col);
        h->SetLineWidth(2);   h->SetStats(0);
    };
    auto styleData = [](TH1 *h)
    {
        h->SetLineColor(kBlack); h->SetMarkerColor(kBlack);
        h->SetMarkerStyle(20);   h->SetMarkerSize(1.0f);
        h->SetStats(0);
    };
    auto overlayVerticalErrors = [](TH1 *h)
    {
        if (!h) return;
        TH1 *herr = dynamic_cast<TH1 *>(h->Clone(Form("%s_err", h->GetName())));
        if (!herr) return;
        herr->SetDirectory(nullptr);
        herr->SetFillStyle(0);
        herr->SetMarkerStyle(1); herr->SetMarkerSize(0);
        herr->SetLineColor(h->GetLineColor());
        herr->SetLineWidth(h->GetLineWidth());
        herr->SetMarkerColor(h->GetLineColor());
        herr->SetBit(kCanDelete);
        herr->Draw("E0 SAME");
    };
    auto drawLabels = [](float x, float y)
    {
        myText(x, y,        1, strleg1.c_str(),         0.04);
        myText(x, y - 0.05, 1, strleg2.c_str(),         0.04);
        myText(x, y - 0.10, 1, kPaperEtaLabel.c_str(),  0.04);
    };

    for (const auto &hbase : kPaperVars)
    {
        TString xaxisname = hbase.substr(4, hbase.size() - 4);

        float xaxismin = 0.f, xaxismax = 1.f;
        int   nrebin   = 4;
        if (xaxisname[0] == 'w')                       xaxismax = 2.f;
        if (xaxisname.CompareTo("e32_to_e35") == 0)  { xaxismin = 0.4f; xaxismax = 1.f; nrebin = 1; }

        TString histNameFull = Form("%s_eta%d_pt%d_cut%d",
                                    hbase.c_str(), kEtaBin, kPtBin, kCutBin);
        TString histNamesave = Form("%s_eta%d_pt%d_cut%d",
                                    xaxisname.Data(), kEtaBin, kPtBin, kCutBin);

        TH2F *h2_data = dynamic_cast<TH2F *>(f_data->Get(histNameFull));
        TH2F *h2_sig  = dynamic_cast<TH2F *>(f_sig ->Get(histNameFull));
        TH2F *h2_bkg  = dynamic_cast<TH2F *>(f_bkg ->Get(histNameFull));
        if (!h2_data || !h2_sig || !h2_bkg)
        {
            std::cerr << "[plot_paper_showershape] missing " << histNameFull << " in one or more files; skipping." << std::endl;
            continue;
        }
        for (TH2F *h : {h2_bkg, h2_sig, h2_data})
        {
            h->RebinX(nrebin);
            h->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
        }

        TH1D *proj_data = h2_data->ProjectionX(Form("%s_px_data", histNameFull.Data()));
        TH1D *proj_sig  = h2_sig ->ProjectionX(Form("%s_px_sig",  histNameFull.Data()));
        TH1D *proj_bkg  = h2_bkg ->ProjectionX(Form("%s_px_bkg",  histNameFull.Data()));

        scaleToUnit(proj_data);
        scaleToUnit(proj_sig);
        scaleToUnit(proj_bkg);
        styleData(proj_data);
        styleMC(proj_sig, kRed);
        styleMC(proj_bkg, kBlue);

        float maxy = std::max({proj_data->GetMaximum(),
                               proj_sig ->GetMaximum(),
                               proj_bkg ->GetMaximum()});

        // Single-pad canvas (paper draft -- no diff pad).
        TCanvas *c_proj = new TCanvas(Form("c_paper_%s", histNameFull.Data()),
                                      Form("ProjectionX - %s", histNameFull.Data()),
                                      600, 600);
        c_proj->cd();

        proj_sig->SetYTitle("normalized counts");
        proj_sig->GetYaxis()->SetTitleOffset(1.5);
        proj_sig->SetXTitle(xaxisname.Data());
        proj_sig->GetYaxis()->SetRangeUser(0, maxy * 1.3);
        proj_sig->GetXaxis()->SetNdivisions(505);
        proj_sig->SetStats(0);
        proj_sig->Draw("HIST");
        proj_bkg->Draw("HIST SAME");
        overlayVerticalErrors(proj_sig);
        overlayVerticalErrors(proj_bkg);
        proj_data->Draw("ex0 SAME");

        drawLabels(0.20, 0.90);
        // Two-line bin label so it stays clear of the legend on the right.
        myText(0.20, 0.75, 1,
               Form("%.0f < #it{E}_{T} < %.0f GeV", kPtLow, kPtHigh),
               0.04);
        myText(0.20, 0.70, 1, kCutLabel.c_str(), 0.04);
        {
            const double y_top = 0.90;
            const double y_bot = 0.72;
            TLegend *leg = new TLegend(0.62, y_bot, 0.92, y_top);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.045);
            leg->AddEntry(proj_data, "Data",         "lep");
            leg->AddEntry(proj_sig,  "Signal MC",    "l");
            leg->AddEntry(proj_bkg,  "Inclusive MC", "l");
            leg->Draw();
        }

        c_proj->SaveAs(Form("%s/dis_%s.pdf", paper_savepath().c_str(), histNamesave.Data()));
    }

    f_data->Close();
    f_sig ->Close();
    f_bkg ->Close();
}
