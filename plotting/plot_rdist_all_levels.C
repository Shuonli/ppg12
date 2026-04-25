#include "plotcommon.h"

// 4-level panel of R = p_data / p_MC distributions with core-Gaussian
// fit overlay for each. Panels: preselect, common, tight, tight_iso.

static std::pair<double, double> fit_gauss_3sigma(const std::vector<double> &X,
                                                  double *out_median = nullptr)
{
    // 3-sigma iterative Gaussian fit. Works well when input X is nearly
    // Gaussian (we feed it log(R), not R, so this is satisfied).
    if (X.empty()) return {0.0, 1.0};
    std::vector<double> sorted = X;
    std::sort(sorted.begin(), sorted.end());
    if (out_median) *out_median = sorted[sorted.size()/2];
    double mean = 0, sigma = 0;
    for (double x : X) mean += x; mean /= X.size();
    for (double x : X) sigma += (x - mean) * (x - mean);
    sigma = std::sqrt(sigma / X.size());
    for (int it = 0; it < 3; ++it) {
        double m2 = 0, s2 = 0; int nk = 0;
        for (double x : X) if (std::fabs(x - mean) < 3.0 * sigma) { m2 += x; nk++; }
        if (nk == 0) break;
        mean = m2 / nk;
        for (double x : X) if (std::fabs(x - mean) < 3.0 * sigma) s2 += (x - mean) * (x - mean);
        sigma = std::sqrt(s2 / nk);
    }
    return {mean, sigma};
}

void plot_rdist_all_levels()
{
    init_plot();
    const char *infile = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root";
    TFile *f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) { std::cerr << "cannot open " << infile << std::endl; return; }

    std::vector<std::string> levels = {"preselect", "common", "tight", "tight_iso"};

    TCanvas *c = new TCanvas("c_rdist", "", 1400, 1000);
    c->Divide(2, 2, 0.02, 0.02);

    for (size_t i = 0; i < levels.size(); ++i) {
        const std::string &lvl = levels[i];
        TH2F *h_mc = (TH2F *) f->Get(Form("h_etaphi_tower_%s_mc_inclusive", lvl.c_str()));
        TH2F *h_da = (TH2F *) f->Get(Form("h_etaphi_tower_%s_data", lvl.c_str()));
        if (!h_mc || !h_da) continue;

        double Nm = h_mc->Integral(), Nd = h_da->Integral();
        std::vector<double> logR;
        for (int ix = 1; ix <= h_mc->GetNbinsX(); ++ix)
          for (int iy = 1; iy <= h_mc->GetNbinsY(); ++iy) {
            double nm = h_mc->GetBinContent(ix,iy), nd = h_da->GetBinContent(ix,iy);
            if (nm <= 0 || nd <= 0) continue;
            double r = (nd/Nd) / (nm/Nm);
            if (r <= 0) continue;
            logR.push_back(std::log(r));
          }
        double median;
        auto stats = fit_gauss_3sigma(logR, &median);
        double mean = stats.first, sigma = stats.second;

        TH1F *h = new TH1F(Form("h_logr_%s", lvl.c_str()),
                           Form(";log(R) = log(p_{data}/p_{MC});bins / 0.05"), 160, -4, 4);
        for (double x : logR) h->Fill(x);

        c->cd(i + 1);
        gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.11); gPad->SetBottomMargin(0.16);
        gPad->SetLogy();
        h->SetLineColor(kBlack); h->SetLineWidth(2);
        h->SetXTitle("log(R) = log(p_{data}/p_{MC})");
        h->SetYTitle("bins / 0.05");
        h->GetXaxis()->SetTitleSize(0.045);
        h->GetYaxis()->SetTitleSize(0.045);
        h->GetXaxis()->SetLabelSize(0.04);
        h->GetYaxis()->SetLabelSize(0.04);
        h->GetYaxis()->SetTitleOffset(1.3);
        h->GetXaxis()->SetTitleOffset(1.2);
        h->GetXaxis()->SetRangeUser(-4, 4);
        h->Draw("HIST");

        TF1 *g = new TF1(Form("g_%s", lvl.c_str()),
                         "[0] * TMath::Gaus(x, [1], [2])", -4, 4);
        double norm = h->Integral() * 0.05 / (sigma * std::sqrt(2 * M_PI));
        g->SetParameters(norm, mean, sigma);
        g->SetLineColor(kRed); g->SetLineStyle(7); g->SetLineWidth(2);
        g->Draw("same");

        TLatex lx; lx.SetNDC();
        lx.SetTextSize(0.04);
        lx.DrawLatex(0.17, 0.935,
                     Form("#bf{#it{sPHENIX}} Internal -- %s", lvl.c_str()));
        lx.SetTextSize(0.035);
        lx.DrawLatex(0.55, 0.84, Form("#mu_{log R} = %.3f", mean));
        lx.DrawLatex(0.55, 0.78, Form("#sigma_{log R} = %.3f", sigma));
        lx.DrawLatex(0.55, 0.72, Form("#LT R #GT = e^{#mu} = %.3f", std::exp(mean)));
        lx.DrawLatex(0.55, 0.66, Form("N_{bins} = %zu", logR.size()));

        TLine *ln = new TLine(mean, h->GetMinimum() + 1, mean, h->GetMaximum());
        ln->SetLineStyle(2); ln->SetLineColor(kRed); ln->SetLineWidth(1);
        ln->Draw();
    }

    const char *outdir = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures";
    c->SaveAs(Form("%s/rdist_all_levels.pdf", outdir));
    c->SaveAs(Form("%s/rdist_all_levels.png", outdir));
    f->Close();
}
