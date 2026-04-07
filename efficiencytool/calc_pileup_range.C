#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

// Compute aggregate pile-up statistics over a run range.
//
// Reads per-run ROOT files produced by analyze.C:
//   output/added_509/triggeroutput_nblair_{rn}_clt.root  (TTree "outt")
//
// Aggregate pile-up is computed as a time-weighted mean:
//   mu_corr = Σ(nmbdc) / (rB * Σ(ttseg))
// NOT the arithmetic mean of per-run mu, because runs have different durations.
//
// rB = 111 * 78200 Hz = 8.682 MHz  [Run-3 standard; revisit for other run periods]
//
// Usage:
//   root -b -q -l 'calc_pileup_range.C(47289, 53862)'
//   root -b -q -l 'calc_pileup_range.C(47289, 51274, "goodrunlist.list", "output/added_509", false)'

int calc_pileup_range(
    int rn_min,
    int rn_max,
    const char* goodrunlist = "/sphenix/user/jocl/projects/analysis/LuminosityCounterGoodRuns/run/fullgoodrunlist.list",
    const char* rootdir     = "/sphenix/user/jocl/projects/analysis/LuminosityCounterGoodRuns/run/output/added_509",
    bool make_plots         = true)
{
    // --- Physics constants (Run-3 Au+Au) ---
    const double nBunch    = 111.0;         // filled bunch crossings per revolution
    const double f_rev     = 78200.0;       // Hz, RHIC revolution frequency per bunch
    const double rB        = nBunch * f_rev; // = 8.682e6 Hz, beam crossing rate
    const double sigma_MB  = 25.2e9;        // MBD cross-section in units giving pb^-1 (= 25.2 pb * 1e9)

    // --- Validate range ---
    if (rn_min > rn_max)
    {
        std::cerr << "ERROR: rn_min (" << rn_min << ") > rn_max (" << rn_max << ")" << std::endl;
        return 1;
    }

    // --- Load good run list, filter to [rn_min, rn_max] ---
    std::vector<int> runs;
    {
        std::ifstream rnfile(goodrunlist);
        if (!rnfile.is_open())
        {
            std::cerr << "ERROR: cannot open good run list: " << goodrunlist << std::endl;
            return 1;
        }
        std::string line;
        while (std::getline(rnfile, line))
        {
            if (line.empty()) continue;
            if (!std::isdigit((unsigned char)line[0])) continue;
            int rn = std::stoi(line);
            if (rn >= rn_min && rn <= rn_max)
                runs.push_back(rn);
        }
    }
    if (runs.empty())
    {
        std::cerr << "WARNING: no runs found in [" << rn_min << ", " << rn_max << "] in " << goodrunlist << std::endl;
        return 2;
    }

    // --- Accumulators ---
    double total_nmbdc   = 0.0;
    double total_nocorN  = 0.0;
    double total_rawmbd  = 0.0;  // Σ sumgoodraw[10]
    double total_ttseg   = 0.0;  // seconds
    ULong64_t total_nevt = 0;
    int run_count        = 0;
    int skipped_missing  = 0;
    int skipped_quality  = 0;

    // Luminosity-weighted accumulators for nonlinear pile-up quantities
    // (Jensen's inequality: <f(mu)> != f(<mu>), so must compute per-run then average)
    double wsum_p_pileup = 0.0;
    double wsum_f_single = 0.0;
    double wsum_f_double = 0.0;
    double wsum_f_ge3    = 0.0;

    // --- TGraph arrays for optional plots ---
    int nruns = (int)runs.size();
    std::vector<double> g_rn, g_mu_raw, g_mu_eff, g_mu_corr, g_p_pileup, g_corr_fac;
    if (make_plots)
    {
        g_rn.reserve(nruns);
        g_mu_raw.reserve(nruns);
        g_mu_eff.reserve(nruns);
        g_mu_corr.reserve(nruns);
        g_p_pileup.reserve(nruns);
        g_corr_fac.reserve(nruns);
    }

    // --- Per-run table header ---
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    printf("\n");
    printf("%-9s  %10s  %9s  %9s  %9s  %9s  %9s  %12s\n",
           "RN", "ttseg[s]", "mu_raw", "mu_eff", "mu_corr", "corr_fac", "p_pileup", "L_run[pb-1]");
    std::string sep(100, '-');
    std::cout << sep << std::endl;

    // --- Branch variables (types must match analyze.C TTree declarations) ---
    float     nmbdc_v, nocorN_v, ttseg_v;
    ULong64_t sumgoodraw[64];
    ULong64_t nevt_v;

    // --- Main loop ---
    for (int rn : runs)
    {
        TString filepath = TString::Format("%s/triggeroutput_nblair_%d_clt.root", rootdir, rn);
        TFile* file = TFile::Open(filepath, "READ");
        if (!file || file->IsZombie())
        {
            std::cerr << "WARN: missing file for run " << rn << ", skipping." << std::endl;
            ++skipped_missing;
            if (file) file->Close();
            continue;
        }

        TTree* tree = (TTree*)file->Get("outt");
        if (!tree)
        {
            std::cerr << "WARN: no TTree 'outt' in file for run " << rn << ", skipping." << std::endl;
            file->Close();
            ++skipped_missing;
            continue;
        }

        tree->SetBranchAddress("nmbdc",       &nmbdc_v);
        tree->SetBranchAddress("nocorN",      &nocorN_v);
        tree->SetBranchAddress("ttseg",       &ttseg_v);
        tree->SetBranchAddress("sumgoodraw",   sumgoodraw);
        tree->SetBranchAddress("nevt",        &nevt_v);
        tree->GetEntry(0);
        file->Close();  // safe: branch data already copied to local vars

        // Quality cuts (matching get_luminosity_182630.C:98)
        if (ttseg_v <= 0.0f)
        {
            std::cerr << "WARN: run " << rn << " ttseg=0, skipping." << std::endl;
            ++skipped_quality;
            continue;
        }
        if (nocorN_v <= 0.0f)
        {
            std::cerr << "WARN: run " << rn << " nocorN<=0, skipping." << std::endl;
            ++skipped_quality;
            continue;
        }
        if ((double)nmbdc_v / (double)nocorN_v < 1.0)
        {
            std::cerr << "WARN: run " << rn << " nmbdc/nocorN=" << (double)nmbdc_v/(double)nocorN_v
                      << " < 1 (unphysical), skipping." << std::endl;
            ++skipped_quality;
            continue;
        }

        // Per-run quantities
        double ttseg_d   = (double)ttseg_v;
        double nmbdc_d   = (double)nmbdc_v;
        double nocorN_d  = (double)nocorN_v;
        double rawmbd_d  = (double)sumgoodraw[10];

        double mu_raw    = rawmbd_d  / (rB * ttseg_d);
        double mu_eff    = nocorN_d  / (rB * ttseg_d);  // after MBD efficiency map
        double mu_corr   = nmbdc_d   / (rB * ttseg_d);  // full Poisson pile-up correction
        double corr_fac  = nmbdc_d   / nocorN_d;
        double L_run     = nmbdc_d   / sigma_MB;

        // Pile-up event fraction: P(>=2 collisions | >=1 collision)
        double p_pileup  = (1.0 - std::exp(-mu_corr) * (1.0 + mu_corr)) / (1.0 - std::exp(-mu_corr));

        // Per-run interaction-mix fractions, conditioned on >=1 collision (triggered events)
        // P(n|n>=1) = P(n) / P(n>=1),  where P(n) = Poisson(n; mu)
        double p_ge1_run    = 1.0 - std::exp(-mu_corr);                                // P(n>=1)
        double f_single_run = mu_corr * std::exp(-mu_corr) / p_ge1_run;                // P(n=1|n>=1)
        double f_double_run = 0.5 * mu_corr * mu_corr * std::exp(-mu_corr) / p_ge1_run; // P(n=2|n>=1)
        double f_ge3_run    = 1.0 - f_single_run - f_double_run;                        // P(n>=3|n>=1)

        printf("%-9d  %10.1f  %9.5f  %9.5f  %9.5f  %9.5f  %9.5f  %12.4e\n",
               rn, ttseg_d, mu_raw, mu_eff, mu_corr, corr_fac, p_pileup, L_run);

        // Accumulate
        total_nmbdc  += nmbdc_d;
        total_nocorN += nocorN_d;
        total_rawmbd += rawmbd_d;
        total_ttseg  += ttseg_d;
        total_nevt   += nevt_v;
        ++run_count;

        // Luminosity-weight: L_i ∝ nmbdc_i (sigma_MB cancels in ratio)
        wsum_p_pileup += nmbdc_d * p_pileup;
        wsum_f_single += nmbdc_d * f_single_run;
        wsum_f_double += nmbdc_d * f_double_run;
        wsum_f_ge3    += nmbdc_d * f_ge3_run;

        if (make_plots)
        {
            g_rn.push_back((double)rn);
            g_mu_raw.push_back(mu_raw);
            g_mu_eff.push_back(mu_eff);
            g_mu_corr.push_back(mu_corr);
            g_corr_fac.push_back(corr_fac);
            g_p_pileup.push_back(p_pileup);
        }
    }

    // --- Guard against empty result ---
    if (run_count == 0 || total_ttseg <= 0.0)
    {
        std::cerr << "ERROR: no valid runs found after quality cuts." << std::endl;
        return 3;
    }

    // --- Aggregate quantities ---
    double total_crossings = rB * total_ttseg;
    double mu_raw_agg   = total_rawmbd  / total_crossings;
    double mu_eff_agg   = total_nocorN  / total_crossings;
    double mu_corr_agg  = total_nmbdc   / total_crossings;
    double corr_fac_agg = total_nmbdc   / total_nocorN;
    double L_int        = total_nmbdc   / sigma_MB;

    // --- Luminosity-weighted interaction-mix fractions ---
    // Nonlinear functions of mu (pile-up fraction, f_n) must be computed per-run
    // and then luminosity-weighted, because <f(mu)> != f(<mu>) (Jensen's inequality).
    // Weight = L_i ∝ nmbdc_i; sigma_MB cancels in the ratio.
    double p_pileup_agg = wsum_p_pileup / total_nmbdc;
    double f_single     = wsum_f_single / total_nmbdc;
    double f_double     = wsum_f_double / total_nmbdc;
    double f_ge3        = wsum_f_ge3    / total_nmbdc;

    // --- Summary ---
    std::cout << sep << std::endl;
    printf("\n=== SUMMARY: Run range [%d, %d] ===\n", rn_min, rn_max);
    printf("Runs processed : %d  (skipped: %d missing, %d quality-cut)\n",
           run_count, skipped_missing, skipped_quality);
    printf("Total time     : %.1f s  (%.2f hours)\n", total_ttseg, total_ttseg / 3600.0);
    printf("rB (assumed)   : %.1f Hz  (%.0f bunches x %.0f Hz) [Run-3 standard]\n", rB, nBunch, f_rev);
    printf("\nAggregate pile-up:\n");
    printf("  mu_raw   = sum(raw MBD counts) / (rB * sum_ttseg) = %.5f\n", mu_raw_agg);
    printf("  mu_eff   = sum(nocorN)         / (rB * sum_ttseg) = %.5f  [after MBD eff map]\n", mu_eff_agg);
    printf("  mu_corr  = sum(nmbdc)          / (rB * sum_ttseg) = %.5f  [full Poisson correction]\n", mu_corr_agg);
    printf("\nOverall correction factor  nmbdc/nocorN  = %.5f\n", corr_fac_agg);
    printf("Pile-up event fraction P(>=2 | >=1)     = %.5f  [lumi-weighted]\n", p_pileup_agg);
    printf("\nTriggered interaction-mix fractions  P(n|n>=1)  (lumi-weighted avg of per-run):\n");
    printf("  f_single  [n=1]  = <mu*e^-mu / (1-e^-mu)>_L      = %.6f\n", f_single);
    printf("  f_double  [n=2]  = <mu^2/2*e^-mu / (1-e^-mu)>_L  = %.6f\n", f_double);
    printf("  f_triple+ [n>=3] = <1 - f1 - f2>_L                = %.6f\n", f_ge3);
    printf("  (sum check: f_single + f_double + f_ge3 = %.6f)\n", f_single + f_double + f_ge3);
    double p_ge1_agg = 1.0 - std::exp(-mu_corr_agg);
    printf("  (cf. from aggregate mu: f1(mu_agg) = %.6f, f2(mu_agg) = %.6f)\n",
           mu_corr_agg * std::exp(-mu_corr_agg) / p_ge1_agg,
           0.5 * mu_corr_agg * mu_corr_agg * std::exp(-mu_corr_agg) / p_ge1_agg);

    // Cluster-weighted mixing fractions: a double-interaction event contributes
    // ~2x as many clusters, triple ~3x, etc.  w(n) = n*f(n) / sum(k*f(k)).
    // Triple+ is folded into double (only double-interaction MC is available).
    double denom_cw = 1.0 * f_single + 2.0 * f_double + 3.0 * f_ge3;
    double cw_single = (1.0 * f_single) / denom_cw;
    double cw_double = (2.0 * f_double + 3.0 * f_ge3) / denom_cw;
    printf("\nCluster-weighted mixing fractions  w(n) = n*f(n)/sum(k*f(k))  (triple+ folded into double):\n");
    printf("  cw_single = %.6f\n", cw_single);
    printf("  cw_double = %.6f\n", cw_double);

    printf("\nIntegrated luminosity (MBD corrected)   = %.4f pb^-1\n", L_int);
    printf("Total events                            = %llu\n\n", (unsigned long long)total_nevt);

    // --- Optional plots ---
    if (make_plots && run_count > 0)
    {
        int n = (int)g_rn.size();

        // Canvas 1: mu vs run number
        TCanvas* c1 = new TCanvas("c_mu","mu vs run",1400,600);
        c1->SetLeftMargin(0.08);
        c1->SetRightMargin(0.04);

        TGraph* gr_raw  = new TGraph(n, g_rn.data(), g_mu_raw.data());
        TGraph* gr_eff  = new TGraph(n, g_rn.data(), g_mu_eff.data());
        TGraph* gr_corr = new TGraph(n, g_rn.data(), g_mu_corr.data());

        gr_raw ->SetMarkerStyle(24); gr_raw ->SetMarkerColor(kBlue);    gr_raw ->SetMarkerSize(0.7);
        gr_eff ->SetMarkerStyle(25); gr_eff ->SetMarkerColor(kGreen+2); gr_eff ->SetMarkerSize(0.7);
        gr_corr->SetMarkerStyle(20); gr_corr->SetMarkerColor(kRed);     gr_corr->SetMarkerSize(0.7);
        gr_raw ->SetLineColor(kBlue);    gr_raw ->SetLineWidth(0);
        gr_eff ->SetLineColor(kGreen+2); gr_eff ->SetLineWidth(0);
        gr_corr->SetLineColor(kRed);     gr_corr->SetLineWidth(0);

        // Draw using gr_corr as frame (largest values)
        gr_corr->SetTitle(TString::Format(";Run Number;#mu (mean collisions/crossing)  [%d-%d]", rn_min, rn_max));
        gr_corr->Draw("AP");
        gr_eff ->Draw("P SAME");
        gr_raw ->Draw("P SAME");

        TLegend* leg1 = new TLegend(0.65, 0.72, 0.94, 0.90);
        leg1->SetBorderSize(0); leg1->SetFillStyle(0);
        leg1->AddEntry(gr_corr, "#mu_{corr} (Poisson corrected)", "p");
        leg1->AddEntry(gr_eff,  "#mu_{eff}  (MBD eff. map)",      "p");
        leg1->AddEntry(gr_raw,  "#mu_{raw}  (MBD scaler only)",   "p");
        leg1->Draw();

        TString outname1 = TString::Format("pileup_mu_vs_rn_%d_%d.png", rn_min, rn_max);
        c1->SaveAs(outname1);
        std::cout << "Saved: " << outname1 << std::endl;

        // Canvas 2: pile-up event fraction vs run number
        TCanvas* c2 = new TCanvas("c_pp","pileup fraction vs run",1400,600);
        c2->SetLeftMargin(0.08);
        c2->SetRightMargin(0.04);

        TGraph* gr_pp = new TGraph(n, g_rn.data(), g_p_pileup.data());
        gr_pp->SetMarkerStyle(20); gr_pp->SetMarkerColor(kRed); gr_pp->SetMarkerSize(0.7);
        gr_pp->SetTitle(TString::Format(";Run Number;P(#geq2 | #geq1 collision)  [%d-%d]", rn_min, rn_max));
        gr_pp->Draw("AP");

        // Reference line at aggregate value
        TLine* hline = new TLine(g_rn.front(), p_pileup_agg, g_rn.back(), p_pileup_agg);
        hline->SetLineColor(kBlue); hline->SetLineWidth(2); hline->SetLineStyle(2);
        hline->Draw();

        TLegend* leg2 = new TLegend(0.65, 0.78, 0.94, 0.88);
        leg2->SetBorderSize(0); leg2->SetFillStyle(0);
        leg2->AddEntry(gr_pp, "Per-run P(#geq2|#geq1)", "p");
        leg2->AddEntry(hline, TString::Format("Aggregate = %.4f", p_pileup_agg), "l");
        leg2->Draw();

        TString outname2 = TString::Format("pileup_fraction_vs_rn_%d_%d.png", rn_min, rn_max);
        c2->SaveAs(outname2);
        std::cout << "Saved: " << outname2 << std::endl;
    }

    return 0;
}
