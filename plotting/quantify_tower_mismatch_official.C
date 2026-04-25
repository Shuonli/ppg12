#include "plotcommon.h"

// Quantify the data vs inclusive-MC tower acceptance mismatch at each of the
// 4 selection levels using the R-distribution z-score (empirical mean/sigma
// from the core of R = p_data/p_MC, iterative core-Gauss fit). Per-bin
// Poisson is NOT used because the lumi-weighted combined MC violates it.

static std::pair<double, double> fit_logR(TH2F *h_mc, TH2F *h_da, const char *lbl)
{
    // log(R) z-score. Per-bin count ratios are log-normal, not normal.
    const double Nm = h_mc->Integral(), Nd = h_da->Integral();
    std::vector<double> logR;
    for (int ix = 1; ix <= h_mc->GetNbinsX(); ++ix)
      for (int iy = 1; iy <= h_mc->GetNbinsY(); ++iy) {
        double nm = h_mc->GetBinContent(ix,iy), nd = h_da->GetBinContent(ix,iy);
        if (nm <= 0 || nd <= 0) continue;
        double r = (nd/Nd) / (nm/Nm);
        if (r <= 0) continue;
        logR.push_back(std::log(r));
      }
    if (logR.empty()) return {0.0, 1.0};
    double mean = 0, sigma = 0;
    for (double x : logR) mean += x; mean /= logR.size();
    for (double x : logR) sigma += (x - mean) * (x - mean);
    sigma = std::sqrt(sigma / logR.size());
    for (int it = 0; it < 3; ++it) {
        double m2 = 0, s2 = 0; int nk = 0;
        for (double x : logR) if (std::fabs(x - mean) < 3.0 * sigma) { m2 += x; nk++; }
        if (nk == 0) break;
        mean = m2 / nk;
        for (double x : logR) if (std::fabs(x - mean) < 3.0 * sigma) s2 += (x - mean) * (x - mean);
        sigma = std::sqrt(s2 / nk);
    }
    std::cout << "[" << lbl << "] log(R) fit: mean=" << mean << " sigma=" << sigma
              << " (<R>=" << std::exp(mean) << ")" << std::endl;
    return {mean, sigma};
}

static void analyse_level(TH2F *h_mc, TH2F *h_da, const std::string &lvl,
                          const char *outdir)
{
    if (!h_mc || !h_da) { return; }
    auto stats = fit_logR(h_mc, h_da, lvl.c_str());
    double mean_lR = stats.first, sigma_lR = stats.second;

    const int nx = h_mc->GetNbinsX();
    const int ny = h_mc->GetNbinsY();
    double N_mc = h_mc->Integral(), N_da = h_da->Integral();
    if (N_mc <= 0 || N_da <= 0) return;

    TH2F *h_ratio = (TH2F *) h_mc->Clone(Form("h_ratio_%s", lvl.c_str()));
    TH2F *h_z     = (TH2F *) h_mc->Clone(Form("h_z_%s",     lvl.c_str()));
    h_ratio->Reset(); h_z->Reset();
    h_ratio->SetZTitle("p_{data}/p_{MC}");
    h_z->SetZTitle("z-score (R-stat)");
    h_ratio->SetTitle(Form("Data/inclusive-MC ratio (%s)", lvl.c_str()));
    h_z->SetTitle(Form("R-stat z-score (%s)", lvl.c_str()));

    int n_fid = 0, n_zlo = 0, n_zhi = 0, n_zvlo = 0, n_zvhi = 0;
    int n_empty_mc = 0, n_empty_da = 0;
    double def_zlt_m2 = 0, def_zlt_m5 = 0;

    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            double nm = h_mc->GetBinContent(ix, iy);
            double nd = h_da->GetBinContent(ix, iy);
            n_fid++;
            if (nm <= 0) { n_empty_mc++; continue; }
            if (nd <= 0)   n_empty_da++;
            double p_m = nm / N_mc;
            double p_d = nd / N_da;
            double R   = p_d / p_m;
            double z   = (R > 0) ? (std::log(R) - mean_lR) / sigma_lR : -999;
            h_ratio->SetBinContent(ix, iy, R);
            h_z->SetBinContent(ix, iy, z);
            if (z < -2) {
                n_zlo++;
                double mu = nm * (N_da / N_mc);
                def_zlt_m2 += std::max(0.0, mu - nd);
            }
            if (z >  2) n_zhi++;
            if (z < -5) {
                n_zvlo++;
                double mu = nm * (N_da / N_mc);
                def_zlt_m5 += std::max(0.0, mu - nd);
            }
            if (z >  5) n_zvhi++;
        }
    }

    double f_anom = (double)(n_zlo + n_zhi) / n_fid;
    std::cout << "---------------------------------------------------------\n"
              << "Level: " << lvl << " (log(R)-stat z)\n"
              << "  log(R) fit: mean = " << mean_lR << "  sigma = " << sigma_lR << "\n"
              << "  N_mc_incl = " << N_mc << "  N_data = " << N_da << "\n"
              << "  |z|>2 total         : " << (n_zlo + n_zhi)
              << " (" << 100.0 * f_anom << " %)\n"
              << "     z<-2 (data<MC)   : " << n_zlo
              << " (" << 100.0 * n_zlo / (double)n_fid << " %)\n"
              << "     z>+2 (MC<data)   : " << n_zhi
              << " (" << 100.0 * n_zhi / (double)n_fid << " %)\n"
              << "  |z|>5 total         : " << (n_zvlo + n_zvhi) << "\n"
              << "Data deficit z<-2 : " << def_zlt_m2
              << " clusters (" << 100.0 * def_zlt_m2 / N_da << " % of data)\n"
              << "Data deficit z<-5 : " << def_zlt_m5
              << " clusters (" << 100.0 * def_zlt_m5 / N_da << " % of data)\n"
              << "---------------------------------------------------------" << std::endl;

    auto draw_one = [&](TH2F *h, const char *title, const char *subtitle,
                        const char *savename, bool logz,
                        double xmin, double xmax, double zmin, double zmax)
    {
        TCanvas *cc = new TCanvas(Form("c_%s", savename), "", 1000, 750);
        cc->SetRightMargin(0.17); cc->SetLeftMargin(0.12); cc->SetTopMargin(0.11);
        if (logz) cc->SetLogz();
        h->GetXaxis()->SetRangeUser(xmin, xmax);
        if (zmin != zmax) h->GetZaxis()->SetRangeUser(zmin, zmax);
        h->GetXaxis()->SetTitle("cluster i#it{#eta} (tower)");
        h->GetYaxis()->SetTitle("cluster i#it{#phi} (tower)");
        h->Draw("COLZ");
        TLatex lx; lx.SetNDC();
        lx.SetTextSize(0.034);
        lx.DrawLatex(0.13, 0.955, title);
        lx.SetTextSize(0.025);
        lx.DrawLatex(0.13, 0.915, subtitle);
        cc->SaveAs(Form("%s/official_tower_%s.pdf", outdir, savename));
        cc->SaveAs(Form("%s/official_tower_%s.png", outdir, savename));
    };

    draw_one(h_mc,
             Form("#bf{#it{sPHENIX}} Internal -- inclusive MC (%s)", lvl.c_str()),
             "16.6 pb^{-1} (1.5 mrad), all-range merge",
             (lvl + "_mc").c_str(), true, 0, 96, 0, 0);
    draw_one(h_da,
             Form("#bf{#it{sPHENIX}} Internal -- Data (%s)", lvl.c_str()),
             "Photon-4-GeV trigger, all-range",
             (lvl + "_data").c_str(), true, 0, 96, 0, 0);
    draw_one(h_ratio,
             Form("Data/MC ratio (%s)", lvl.c_str()),
             Form("log(R): #mu=%.3f #sigma=%.3f; z<-2: %.2f%%, z>+2: %.2f%%",
                  mean_lR, sigma_lR, 100.0 * n_zlo / (double)n_fid,
                  100.0 * n_zhi / (double)n_fid),
             (lvl + "_ratio").c_str(), false, 0, 96, 0.0, 2.0);
    draw_one(h_z,
             Form("log(R) z-score (%s)", lvl.c_str()),
             Form("log(R) #mu=%.3f #sigma=%.3f; def(z<-2)=%.2f%% of data",
                  mean_lR, sigma_lR, 100.0 * def_zlt_m2 / N_da),
             (lvl + "_zscore").c_str(), false, 0, 96, -6.0, 6.0);
}

void quantify_tower_mismatch_official()
{
    init_plot();
    const char *infile = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root";
    TFile *f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) { std::cerr << "cannot open " << infile << std::endl; return; }

    const std::vector<std::string> levels = {"preselect", "common", "tight", "tight_iso"};
    const char *outdir = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures";

    for (const auto &lvl : levels)
    {
        TH2F *h_mc = (TH2F *) f->Get(Form("h_etaphi_tower_%s_mc_inclusive", lvl.c_str()));
        TH2F *h_da = (TH2F *) f->Get(Form("h_etaphi_tower_%s_data", lvl.c_str()));
        analyse_level(h_mc, h_da, lvl, outdir);
    }
}
