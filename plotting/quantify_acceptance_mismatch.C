#include "plotcommon.h"

// Quantify the data/MC acceptance mismatch in the fiducial |eta|<0.7 region.
// Strategy:
//   1. Load the two eta-phi TH2 histograms already produced by
//      plot_photon10_etaphi_acceptance.C and plot_data_etaphi_acceptance.C.
//   2. Restrict to |eta| < 0.7 (the analysis fiducial region).
//   3. Normalize each to unit integral over the fiducial region so that
//      physics shape differences are removed to first order.
//   4. For each bin compute R = p_data / p_MC and
//      sigma(R)/R = sqrt(1/n_data + 1/n_MC) (Poisson in both).
//   5. z[i,j] = (R - 1) / (sigma(R)/R * R) = (R - 1) / sigma(R)
//   6. Count bins with z < -2 (dead in data relative to MC) and z > +2
//      (dead in MC relative to data), report fractional sky coverage.
//   7. Save a 2x2 summary panel PDF.

void quantify_acceptance_mismatch()
{
    init_plot();

    const char *fmc   = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures/photon10_etaphi_acceptance.root";
    const char *fdata = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures/data_etaphi_acceptance.root";

    TFile *f_mc = TFile::Open(fmc,   "READ");
    TFile *f_da = TFile::Open(fdata, "READ");
    TH2F *h_mc = (TH2F *) f_mc->Get("h_etaphi");
    TH2F *h_da = (TH2F *) f_da->Get("h_etaphi");
    if (!h_mc || !h_da) { std::cerr << "missing histos\n"; return; }

    const int nx = h_mc->GetNbinsX();
    const int ny = h_mc->GetNbinsY();

    // Integrate over fiducial |eta|<0.7 to normalize
    // eta axis: 96 bins over [-1.1, 1.1] -> each bin is 2.2/96 ~ 0.02292 wide
    // |eta|<0.7 => bins with centre in [-0.7, 0.7]
    double N_mc_fid = 0, N_da_fid = 0;
    for (int ix = 1; ix <= nx; ++ix) {
        double eta = h_mc->GetXaxis()->GetBinCenter(ix);
        if (std::fabs(eta) > 0.7) continue;
        for (int iy = 1; iy <= ny; ++iy) {
            N_mc_fid += h_mc->GetBinContent(ix, iy);
            N_da_fid += h_da->GetBinContent(ix, iy);
        }
    }
    std::cout << "N_mc_fid  = " << N_mc_fid  << "\n"
              << "N_da_fid  = " << N_da_fid  << std::endl;

    // Build ratio and z-score maps
    TH2F *h_ratio = (TH2F *) h_mc->Clone("h_ratio");
    TH2F *h_z     = (TH2F *) h_mc->Clone("h_z");
    h_ratio->Reset(); h_z->Reset();
    h_ratio->SetTitle("");
    h_ratio->SetXTitle("#it{#eta^{cluster}}");
    h_ratio->SetYTitle("#it{#phi^{cluster}} [rad]");
    h_ratio->SetZTitle("p_{data}/p_{MC}");
    h_z->SetTitle("");
    h_z->SetXTitle("#it{#eta^{cluster}}");
    h_z->SetYTitle("#it{#phi^{cluster}} [rad]");
    h_z->SetZTitle("z-score");

    TH1F *h_zdist = new TH1F("h_zdist", ";z-score;bins", 120, -12, 12);

    int n_bins_fid      = 0;
    int n_bins_z_low    = 0;  // z < -2: dead in data vs MC
    int n_bins_z_high   = 0;  // z > +2: dead in MC vs data
    int n_bins_z_vlow   = 0;  // z < -5
    int n_bins_z_vhigh  = 0;  // z > +5
    int n_bins_empty_da = 0;
    int n_bins_empty_mc = 0;
    // yield-deficit accounting:
    // expected data count in bin i, given MC shape and data total: mu_i = n_mc * (N_da/N_mc)
    // deficit_i = max(0, mu_i - n_da_i)  for z<-2 bins only.
    // total missing data clusters from data-only dead regions divided by N_da_fid
    double deficit_z_lt_m2 = 0;
    double deficit_z_lt_m5 = 0;

    for (int ix = 1; ix <= nx; ++ix) {
        double eta = h_mc->GetXaxis()->GetBinCenter(ix);
        if (std::fabs(eta) > 0.7) continue;
        for (int iy = 1; iy <= ny; ++iy) {
            double n_mc = h_mc->GetBinContent(ix, iy);
            double n_da = h_da->GetBinContent(ix, iy);
            n_bins_fid++;

            if (n_mc <= 0) { n_bins_empty_mc++; continue; }
            if (n_da <= 0) { n_bins_empty_da++; }

            double p_mc = n_mc / N_mc_fid;
            double p_da = n_da / N_da_fid;
            double R    = p_da / p_mc;
            // relative error: sigma_R/R = sqrt(1/n_d + 1/n_m)
            double relerr = std::sqrt(1.0 / std::max(n_da, 1.0) + 1.0 / n_mc);
            double sig_R  = R * relerr;
            double z      = (R - 1.0) / sig_R;

            h_ratio->SetBinContent(ix, iy, R);
            h_z    ->SetBinContent(ix, iy, z);
            h_zdist->Fill(z);

            if (z < -2) {
                n_bins_z_low++;
                double mu = n_mc * (N_da_fid / N_mc_fid);
                deficit_z_lt_m2 += std::max(0.0, mu - n_da);
            }
            if (z >  2) n_bins_z_high++;
            if (z < -5) {
                n_bins_z_vlow++;
                double mu = n_mc * (N_da_fid / N_mc_fid);
                deficit_z_lt_m5 += std::max(0.0, mu - n_da);
            }
            if (z >  5) n_bins_z_vhigh++;
        }
    }

    int n_bins_anomalous = n_bins_z_low + n_bins_z_high;
    double frac_anomalous  = (double) n_bins_anomalous / n_bins_fid;
    double frac_dead_data  = (double) n_bins_z_low     / n_bins_fid;
    double frac_dead_mc    = (double) n_bins_z_high    / n_bins_fid;

    std::cout << "-------------------------------------------------\n"
              << "Fiducial |eta|<0.7 summary\n"
              << "  total fiducial bins : " << n_bins_fid        << "\n"
              << "  empty MC bins       : " << n_bins_empty_mc   << " (skipped)\n"
              << "  empty data bins     : " << n_bins_empty_da   << "\n"
              << "  |z|>2 total         : " << n_bins_anomalous  << " (" << 100.0 * frac_anomalous << " %)\n"
              << "     z<-2 (data<MC)   : " << n_bins_z_low      << " (" << 100.0 * frac_dead_data << " %)\n"
              << "     z>+2 (MC<data)   : " << n_bins_z_high     << " (" << 100.0 * frac_dead_mc   << " %)\n"
              << "  |z|>5 total         : " << (n_bins_z_vlow + n_bins_z_vhigh)
              << " ("<< 100.0 * (n_bins_z_vlow + n_bins_z_vhigh) / n_bins_fid <<" %)\n"
              << "     z<-5             : " << n_bins_z_vlow     << "\n"
              << "     z>+5             : " << n_bins_z_vhigh    << "\n"
              << "-------------------------------------------------\n";

    // Expected |z|>2 fraction under pure Gaussian stats: 4.55%. Anything beyond
    // this is acceptance-mismatch / physics-shape excess.
    double expected_z2_frac = 0.0455;
    std::cout << "Expected |z|>2 under Poisson (Gaussian approx): "
              << 100.0*expected_z2_frac << " %\n"
              << "Excess over Gaussian: " << 100.0 * (frac_anomalous - expected_z2_frac) << " %\n";
    // Data yield deficit from data-only dead regions (z<-2):
    // The numerator is "extra clusters that MC predicts in bins where data is
    // anomalously low"; divided by data-fiducial total gives the "rough lower
    // bound on fractional acceptance correction missing from the cross-section".
    // Deficit for z<-5: only the clearly-dead bins; safer underestimate.
    std::cout << "Data yield deficit in z<-2 dead regions: "
              << deficit_z_lt_m2 << " clusters, "
              << 100.0 * deficit_z_lt_m2 / N_da_fid << " % of fiducial data count\n"
              << "Data yield deficit in z<-5 dead regions: "
              << deficit_z_lt_m5 << " clusters, "
              << 100.0 * deficit_z_lt_m5 / N_da_fid << " % of fiducial data count\n";

    // -----------------------------------------------------------
    // Individual panels (full-size) for each derived map
    // -----------------------------------------------------------
    const char *outdir = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures";

    auto draw_single = [&](TH2F *h, const char *title, const char *subtitle,
                           const char *savename, bool logz, double zmin, double zmax)
    {
        TCanvas *cc = new TCanvas(Form("c_%s", savename), "", 1000, 750);
        cc->SetRightMargin(0.17);
        cc->SetLeftMargin(0.12);
        cc->SetTopMargin(0.09);
        if (logz) cc->SetLogz();
        h->GetXaxis()->SetRangeUser(-0.7, 0.7);
        if (zmin != zmax) h->GetZaxis()->SetRangeUser(zmin, zmax);
        h->Draw("COLZ");
        TLatex lx; lx.SetNDC();
        lx.SetTextSize(0.04);
        lx.DrawLatex(0.13, 0.935, title);
        lx.SetTextSize(0.03);
        lx.DrawLatex(0.13, 0.895, subtitle);
        cc->SaveAs(Form("%s/acc_mismatch_%s.pdf", outdir, savename));
        cc->SaveAs(Form("%s/acc_mismatch_%s.png", outdir, savename));
    };

    draw_single(h_mc, "#bf{#it{sPHENIX}} Internal -- MC photon10",
                "|z_{vtx}|<10 cm, #it{E}_{T}^{cluster}>10 GeV, |#it{#eta}|<0.7",
                "mc", true, 0, 0);
    draw_single(h_da, "#bf{#it{sPHENIX}} Internal -- Data (Photon 4 GeV trig)",
                "|z_{vtx}|<10 cm, #it{E}_{T}^{cluster}>10 GeV, |#it{#eta}|<0.7",
                "data", true, 0, 0);
    draw_single(h_ratio, "Normalized ratio  (data / MC)",
                Form("z<-2: %.2f%%  z>+2: %.2f%%  total |z|>2: %.2f%%",
                     100.0*frac_dead_data, 100.0*frac_dead_mc, 100.0*frac_anomalous),
                "ratio", false, 0.0, 2.0);
    draw_single(h_z, "z-score (data - MC, Poisson)",
                Form("deficit(z<-2) = %.1f%% of data;  deficit(z<-5) = %.1f%% of data",
                     100.0*deficit_z_lt_m2/N_da_fid, 100.0*deficit_z_lt_m5/N_da_fid),
                "zscore", false, -6.0, 6.0);

    // -----------------------------------------------------------
    // z-score 1D distribution (vs unit gaussian reference)
    // -----------------------------------------------------------
    TCanvas *c2 = new TCanvas("c2", "", 900, 700);
    c2->SetLogy();
    h_zdist->SetXTitle("z-score (data vs MC, Poisson)");
    h_zdist->SetYTitle("bins / 0.2");
    h_zdist->SetLineColor(kBlack);
    h_zdist->SetLineWidth(2);
    h_zdist->Draw("HIST");
    // overlay expected gaussian
    double norm = h_zdist->Integral() * 0.2;
    TF1 *gaus = new TF1("gaus_ref", "[0]*TMath::Gaus(x,0,1,true)", -12, 12);
    gaus->SetParameter(0, norm);
    gaus->SetLineColor(kRed);
    gaus->SetLineWidth(2);
    gaus->SetLineStyle(7);
    gaus->Draw("same");
    TLatex lx; lx.SetNDC(); lx.SetTextSize(0.035);
    lx.DrawLatex(0.15, 0.93, strleg1.c_str());
    lx.DrawLatex(0.15, 0.89, "Red dashed: unit Gaussian expectation");
    c2->SaveAs(Form("%s/acceptance_mismatch_zdist.pdf", outdir));
    c2->SaveAs(Form("%s/acceptance_mismatch_zdist.png", outdir));

    // also save the derived TH2s for later use
    TFile *fo = TFile::Open(Form("%s/acceptance_mismatch.root", outdir), "RECREATE");
    h_ratio->Write("h_ratio");
    h_z    ->Write("h_z");
    h_zdist->Write("h_zdist");
    fo->Close();
}
