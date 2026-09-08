// plot_phisymm_band_example.C
//
// Example per-phi cluster-count distribution for a single 2-eta-row band
// at the `common` selection level, with the iterative-3-sigma Gaussian
// fit, the +/- sigma band, and the -2 sigma flagging threshold drawn
// on top. Illustrates the phi-symmetry dead-tower flagging procedure
// described in Section 12 of the PPG12 analysis note (reviewer comment
// 35).
//
// Output: figures/phisymm_band_example.pdf

#include "plotcommon.h"

namespace {

// 3-sigma iterative clipped Gaussian: same algorithm as
// make_tower_masks.C:fit_gauss_3sigma. Returned as (mu, sigma).
std::pair<double, double> fit_gauss_3sigma(const std::vector<double> &x)
{
    if (x.empty()) return {0.0, 1.0};
    double mean = 0, sigma = 0;
    for (double v : x) mean += v;
    mean /= x.size();
    for (double v : x) sigma += (v - mean) * (v - mean);
    sigma = std::sqrt(sigma / x.size());
    for (int it = 0; it < 3; ++it) {
        double m2 = 0, s2 = 0; int nk = 0;
        for (double v : x) if (std::fabs(v - mean) < 3.0 * sigma) { m2 += v; nk++; }
        if (nk == 0) break;
        mean = m2 / nk;
        for (double v : x) if (std::fabs(v - mean) < 3.0 * sigma) s2 += (v - mean) * (v - mean);
        sigma = std::sqrt(s2 / nk);
    }
    return {mean, sigma};
}

}  // namespace

// ieta_lo/ieta_hi are 0-indexed inclusive tower indices (e.g. 64-67 for
// bands 32 and 33 combined). Sums all rows in [ieta_lo, ieta_hi] into a
// single per-phi distribution and fits a 3-sigma-clipped Gaussian.
void plot_phisymm_band_example(int ieta_lo = 64, int ieta_hi = 67)
{
    init_plot();

    const std::string inFile  =
        "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/"
        "results/Photon_final_bdt_nom.root";
    const std::string outFile = "figures/phisymm_band_example.pdf";

    TFile *f = TFile::Open(inFile.c_str(), "READ");
    if (!f || f->IsZombie()) { printf("cannot open %s\n", inFile.c_str()); return; }

    TH2F *h2 = dynamic_cast<TH2F*>(f->Get("h_etaphi_tower_common_data"));
    if (!h2) { printf("missing h_etaphi_tower_common_data\n"); return; }

    const int nx = h2->GetNbinsX();   // 96 eta towers
    const int ny = h2->GetNbinsY();   // 256 phi towers
    if (ieta_lo < 0 || ieta_hi >= nx || ieta_lo > ieta_hi) {
        printf("ieta range [%d,%d] out of bounds [0,%d)\n",
               ieta_lo, ieta_hi, nx);
        return;
    }
    const int n_rows = ieta_hi - ieta_lo + 1;

    // Sum the requested eta rows into a per-phi count vector
    std::vector<double> band(ny, 0.0);
    for (int de = 0; de < n_rows; ++de) {
        int ix = ieta_lo + de + 1;     // ROOT 1-indexed
        for (int iy = 1; iy <= ny; ++iy)
            band[iy - 1] += h2->GetBinContent(ix, iy);
    }

    auto stats = fit_gauss_3sigma(band);
    const double mu  = stats.first;
    const double sig = stats.second;
    const double thr = mu - 2.0 * sig;

    // Build a TH1F of the per-phi counts for plotting
    TH1F *h_band = new TH1F(Form("h_phi_band_%d_%d", ieta_lo, ieta_hi),
                            "", ny, 0, ny);
    h_band->SetDirectory(nullptr);
    int n_flagged = 0;
    for (int iy = 0; iy < ny; ++iy) {
        h_band->SetBinContent(iy + 1, band[iy]);
        h_band->SetBinError(iy + 1, std::sqrt(std::max(band[iy], 0.0)));
        const double z = (sig > 0) ? (band[iy] - mu) / sig : 0.0;
        if (z < -2.0) ++n_flagged;
    }
    h_band->SetMarkerStyle(20);
    h_band->SetMarkerColor(kBlack);
    h_band->SetLineColor(kBlack);
    h_band->SetMarkerSize(0.6);

    TCanvas *c = new TCanvas("c_phisymm_band", "", 800, 600);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.13);

    h_band->SetStats(0);
    h_band->GetXaxis()->SetTitle("cluster i#it{#phi} (tower index)");
    h_band->GetYaxis()->SetTitle("cluster count in band");
    h_band->GetYaxis()->SetTitleOffset(1.30);
    h_band->GetXaxis()->SetTitleOffset(1.10);
    const double ymax_data = h_band->GetMaximum();
    h_band->GetYaxis()->SetRangeUser(0, 1.4 * ymax_data);
    h_band->Draw("p e");

    // Mu line + +/- sigma band drawn as horizontal lines using TLine
    auto drawHLine = [&](double y, int color, int style, int width) {
        TLine *l = new TLine(0, y, ny, y);
        l->SetLineColor(color);
        l->SetLineStyle(style);
        l->SetLineWidth(width);
        l->Draw("same");
        return l;
    };
    drawHLine(mu,             kRed + 1,    1, 2);
    drawHLine(mu + sig,        kAzure + 2,  2, 1);
    drawHLine(mu - sig,        kAzure + 2,  2, 1);
    drawHLine(thr,            kOrange + 7, 1, 2);

    // Annotation
    const double eta_lo_val = -1.1 + (ieta_lo) * (2.2 / 96.0);
    const double eta_hi_val = -1.1 + (ieta_hi + 1) * (2.2 / 96.0);

    myText(0.18, 0.88, 1, strleg1.c_str(),   0.040, 0);
    myText(0.18, 0.83, 1, strleg2_1.c_str(), 0.040, 0);
    myText(0.18, 0.78, 1, "Data, common pre-selection", 0.040, 0);
    myText(0.18, 0.73, 1,
           Form("i#it{#eta} %d#minus%d (#it{#eta}: %.2f to %.2f)",
                ieta_lo, ieta_hi, eta_lo_val, eta_hi_val), 0.040, 0);
    myText(0.18, 0.68, 1, Form("flagged towers: %d / %d", n_flagged, ny), 0.040, 0);

    TLegend *leg = new TLegend(0.62, 0.66, 0.94, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.038);
    leg->AddEntry(h_band, "Per-#it{#phi} count", "lep");
    {
        TLine *l_mu = new TLine();
        l_mu->SetLineColor(kRed + 1);     l_mu->SetLineStyle(1); l_mu->SetLineWidth(2);
        leg->AddEntry(l_mu, Form("Gauss fit #it{#mu} = %.0f", mu), "l");
        TLine *l_sig = new TLine();
        l_sig->SetLineColor(kAzure + 2);  l_sig->SetLineStyle(2); l_sig->SetLineWidth(1);
        leg->AddEntry(l_sig, Form("#it{#mu} #pm #it{#sigma} (#sigma = %.0f)", sig), "l");
        TLine *l_th = new TLine();
        l_th->SetLineColor(kOrange + 7); l_th->SetLineStyle(1); l_th->SetLineWidth(2);
        leg->AddEntry(l_th, "#it{#mu} #minus 2#it{#sigma} flag threshold", "l");
    }
    leg->Draw();

    c->SaveAs(outFile.c_str());
    printf("ieta=%d-%d  mu=%.2f sigma=%.2f  threshold=mu-2sigma=%.2f  flagged=%d\n",
           ieta_lo, ieta_hi, mu, sig, thr, n_flagged);
    printf("saved %s\n", outFile.c_str());
}
