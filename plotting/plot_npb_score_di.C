#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>
#include <vector>
#include <string>

#include "plotcommon.h"

// Sum h_npb_score_eta0_pt{ipt}_bdt{ibdt} across all (ipt, ibdt) bins from
// the per-sample showershape files of one mix flavour ("nom" = SI,
// "double" = DI), normalised at the end. Each per-sample file already
// carries cross-section x mix-weight in the bin contents, so a plain Add
// across samples gives the inclusive-MC shape.
TH1D *buildInclusiveNPBScore(const std::vector<std::string> &samples,
                             const std::string &mix,
                             const std::string &period,
                             const char *outName,
                             int nPtBins,
                             int nBdtBins)
{
    const std::string base =
        "/sphenix/user/shuhangli/ppg12/efficiencytool/results/";

    TH1D *h_sum = nullptr;
    for (const auto &s : samples) {
        const std::string fn = base + Form(
            "MC_efficiencyshower_shape_%s_%s_inclusive_showershape_%s.root",
            s.c_str(), mix.c_str(), period.c_str());
        TFile *f = TFile::Open(fn.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "[NPB-DI] Skipping missing file: " << fn << std::endl;
            if (f) delete f;
            continue;
        }
        for (int ipt = 0; ipt < nPtBins; ++ipt) {
            for (int ibdt = 0; ibdt < nBdtBins; ++ibdt) {
                TString hname = Form("h_npb_score_eta0_pt%d_bdt%d", ipt, ibdt);
                TH1D *h = dynamic_cast<TH1D *>(f->Get(hname));
                if (!h) continue;
                if (!h_sum) {
                    h_sum = (TH1D *)h->Clone(outName);
                    h_sum->SetDirectory(nullptr);
                } else {
                    h_sum->Add(h);
                }
            }
        }
        f->Close();
        delete f;
    }
    if (!h_sum) {
        std::cerr << "[NPB-DI] Empty sum for mix=" << mix << std::endl;
        return nullptr;
    }
    if (h_sum->Integral() > 0) h_sum->Scale(1.0 / h_sum->Integral());
    return h_sum;
}

void plot_npb_score_di(const std::string &period = "0rad",
                       const std::string &savePath =
                           "../PPG12-analysis-note/Figures/")
{
    init_plot();
    gSystem->Exec(Form("mkdir -p %s", savePath.c_str()));

    // Background MC channels with full GEANT DI partners (matches CLAUDE.md
    // pipeline doc: jet8 / jet12 / jet20 / jet30 / jet40 each have a `_double`
    // companion). jet5 has no DI partner so it is excluded here.
    std::vector<std::string> samples = {"jet8", "jet12", "jet20", "jet30", "jet40"};
    const int nPtBins = 5;   // matches config_showershape pT_bins [10,14,18,22,28,30]
    const int nBdtBins = 3;  // matches config bdt_bins [0, 0.3, 0.7, 1.0]

    TH1D *h_si = buildInclusiveNPBScore(samples, "nom",    period, "h_npb_si", nPtBins, nBdtBins);
    TH1D *h_di = buildInclusiveNPBScore(samples, "double", period, "h_npb_di", nPtBins, nBdtBins);
    if (!h_si || !h_di) {
        std::cerr << "[NPB-DI] Cannot build SI or DI shape — abort." << std::endl;
        return;
    }

    h_si->SetLineColor(kBlue);
    h_si->SetLineWidth(2);
    h_si->SetLineStyle(1);
    h_di->SetLineColor(kRed);
    h_di->SetLineWidth(2);
    h_di->SetLineStyle(2);

    h_si->SetTitle("");
    h_si->GetXaxis()->SetTitle("NPB Score");
    h_si->GetYaxis()->SetTitle("Normalized Counts");
    h_si->GetYaxis()->SetTitleOffset(1.4);
    const double ymax = std::max(h_si->GetMaximum(), h_di->GetMaximum()) * 5.0;
    h_si->SetMaximum(ymax);
    h_si->SetMinimum(1e-4);

    TCanvas *c = new TCanvas("c_npb_si_vs_di", "", 800, 600);
    c->SetLogy();
    h_si->Draw("HIST");
    h_di->Draw("HIST SAME");

    // Integrated fraction below NPB < 0.5 for SI and DI. The histograms are
    // already unit-area normalised, so this is a per-mille-level migration
    // metric for the analysis NPB > 0.5 cut.
    const int    bin_lo = 1;
    const int    bin_hi = h_si->FindBin(0.5) - 1;  // strict NPB < 0.5
    const double frac_si_lo = h_si->Integral(bin_lo, bin_hi);
    const double frac_di_lo = h_di->Integral(bin_lo, bin_hi);
    std::cout << "[NPB-DI] integrated fraction with NPB<0.5: SI=" << frac_si_lo
              << "  DI=" << frac_di_lo
              << "  delta=" << (frac_di_lo - frac_si_lo) << std::endl;

    TLegend *leg = new TLegend(0.55, 0.66, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->AddEntry(h_si, Form("Single Interaction (NPB<0.5: %.2g)", frac_si_lo), "l");
    leg->AddEntry(h_di, Form("Double Interaction (NPB<0.5: %.2g)", frac_di_lo), "l");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.18, 0.92, strleg1.c_str());
    latex.DrawLatex(0.18, 0.87, "Inclusive MC (PYTHIA8)");
    latex.DrawLatex(0.18, 0.82, strleg3.c_str());
    const char *xangle = (period == "0rad") ? "0 mrad crossing"
                                            : "1.5 mrad crossing";
    latex.DrawLatex(0.18, 0.77, xangle);
    latex.DrawLatex(0.18, 0.72,
        Form("#Delta(DI#minusSI) below NPB<0.5: %.1g", frac_di_lo - frac_si_lo));

    c->SaveAs((savePath + "npb_score_1d_overlay_di.pdf").c_str());
    delete c;
    delete leg;
    delete h_si;
    delete h_di;

    std::cout << "Done. Saved npb_score_1d_overlay_di.pdf to " << savePath
              << std::endl;
}
