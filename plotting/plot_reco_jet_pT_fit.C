#include "plotcommon.h"

void plot_reco_jet_pT_fit(const std::string &infile =
    "../efficiencytool/results/MC_efficiency_bdt_nom.root")
{
    init_plot();

    TFile *fin = TFile::Open(infile.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "Cannot open " << infile << std::endl;
        return;
    }

    TH1F *h = dynamic_cast<TH1F *>(fin->Get("h_max_reco_jet_pT"));
    if (!h)
    {
        std::cerr << "Missing h_max_reco_jet_pT" << std::endl;
        return;
    }
    h = (TH1F *) h->Clone("h_max_reco_jet_pT_clone");

    // Rebin 0.1 GeV -> 1 GeV for fitting & display.
    h->Rebin(10);
    h->Scale(1.0, "width");
    h->SetLineColor(kBlack);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.6);
    h->SetMarkerColor(kBlack);

    // ------------------------------------------------------------------
    // Fit: "Modified Hagedorn" (same form as plotting/plot_combine.C:68
    // used for the combined photon5+10+20 sum). Power law with log- and
    // linear-in-pT corrections to the exponent:
    //   f(pT) = A * (p0/pT)^( n + c1*ln(pT/p0) + c2*pT )
    // Fit 10-50 GeV. Double-fit for convergence (matches plot_combine.C).
    // ------------------------------------------------------------------
    const double fit_lo = 10.0;
    const double fit_hi = 50.0;
    TF1 *f_pow = new TF1("f_pow",
        "[0]*TMath::Power([1]/x, [2] + [3]*TMath::Log(x/[1]) + [4]*x)",
        5, 80);
    f_pow->SetParameters(1e16, 1.0, 1.0, 2.0, 0.01);
    f_pow->SetParNames("A", "p0", "n", "c1", "c2");
    f_pow->SetLineColor(kRed + 1);
    f_pow->SetLineWidth(3);
    f_pow->SetLineStyle(2);

    h->Fit(f_pow, "REMQN", "", fit_lo, fit_hi);
    TFitResultPtr r = h->Fit(f_pow, "RSQN0", "", fit_lo, fit_hi);
    double A  = f_pow->GetParameter(0);
    double p0 = f_pow->GetParameter(1);
    double n  = f_pow->GetParameter(2);
    double c1 = f_pow->GetParameter(3);
    double c2 = f_pow->GetParameter(4);
    double chi2 = f_pow->GetChisquare();
    int    ndf  = f_pow->GetNDF();
    // Local power-law index n_eff(x) = - d ln f / d ln x. For
    //   ln f = ln A - (n + c1 L + c2 x) * L,  L = ln(x/p0),
    //   n_eff = n + 2 c1 L + c2 x (1 + L).
    auto nEff = [&](double x) {
        double L = std::log(x/p0);
        return n + 2*c1*L + c2*x*(1.0 + L);
    };
    double n20 = nEff(20.0);
    double n30 = nEff(30.0);
    std::cout << "=== Modified Hagedorn  A*(p0/x)^(n + c1*ln(x/p0) + c2*x)"
              << ", range " << fit_lo << "-" << fit_hi << " GeV ===" << std::endl;
    std::cout << "  A  = " << A  << std::endl;
    std::cout << "  p0 = " << p0 << std::endl;
    std::cout << "  n  = " << n  << std::endl;
    std::cout << "  c1 = " << c1 << std::endl;
    std::cout << "  c2 = " << c2 << std::endl;
    std::cout << "  Local power-law index at pT=20 GeV: n_eff = " << n20 << std::endl;
    std::cout << "  Local power-law index at pT=30 GeV: n_eff = " << n30 << std::endl;
    std::cout << "  chi2/ndf = " << chi2 << " / " << ndf
              << " = " << (ndf>0 ? chi2/ndf : 0) << std::endl;

    // ------------------------------------------------------------------
    // Build data/fit ratio over the full spectrum; scan for the widest
    // contiguous range where |ratio - 1| <= eff_tol AND rel stat err
    // <= stat_tol. That's the "MC efficient" range.
    // ------------------------------------------------------------------
    const double eff_tol  = 0.20;   // allow 20% deviation from power law
    const double stat_tol = 0.30;   // require bin rel_err <= 30%

    TH1F *h_ratio = (TH1F *) h->Clone("h_ratio");
    h_ratio->Reset();
    int nb = h->GetNbinsX();
    std::vector<bool> good(nb + 2, false);
    for (int i = 1; i <= nb; ++i)
    {
        double x  = h->GetBinCenter(i);
        double v  = h->GetBinContent(i);
        double e  = h->GetBinError(i);
        double ff = f_pow->Eval(x);
        if (v <= 0 || ff <= 0) continue;
        double r   = v / ff;
        double er  = e / ff;
        h_ratio->SetBinContent(i, r);
        h_ratio->SetBinError(i, er);
        bool ok_eff  = std::fabs(r - 1.0) <= eff_tol;
        bool ok_stat = (e / v) <= stat_tol;
        good[i] = ok_eff && ok_stat;
    }
    // Find the widest contiguous good band.
    int best_lo = -1, best_hi = -1, best_w = 0;
    int cur_lo = -1;
    for (int i = 1; i <= nb; ++i)
    {
        if (good[i])
        {
            if (cur_lo < 0) cur_lo = i;
            int w = i - cur_lo + 1;
            if (w > best_w) { best_w = w; best_lo = cur_lo; best_hi = i; }
        }
        else
        {
            cur_lo = -1;
        }
    }
    double eff_lo_pT = (best_lo > 0) ? h->GetBinLowEdge(best_lo) : 0;
    double eff_hi_pT = (best_hi > 0) ? h->GetBinLowEdge(best_hi) + h->GetBinWidth(best_hi) : 0;
    std::cout << "\n=== MC-efficient range (|ratio-1| <= " << eff_tol
              << " AND rel_err <= " << stat_tol << ") ===" << std::endl;
    std::cout << "  pT = " << eff_lo_pT << " - " << eff_hi_pT
              << " GeV  (width = " << best_w << " bins)" << std::endl;

    // ------------------------------------------------------------------
    // Layout: upper pad = spectrum + fit; lower pad = data/fit ratio.
    // ------------------------------------------------------------------
    TCanvas *c = new TCanvas("c_reco_jet_fit", "reco jet fit", 800, 800);
    TPad *pU = new TPad("pU", "", 0, 0.32, 1, 1);
    TPad *pL = new TPad("pL", "", 0, 0.00, 1, 0.33);
    pU->SetBottomMargin(0.02);
    pU->SetTickx(); pU->SetTicky();
    pL->SetTopMargin(0.02);
    pL->SetBottomMargin(0.30);
    pL->SetTickx(); pL->SetTicky();
    pU->Draw(); pL->Draw();

    pU->cd();
    pU->SetLogy();
    TH1F *frameU = new TH1F("frameU_reco_jet_pT", "", 100, 0, 60);
    frameU->SetXTitle("");
    frameU->SetYTitle("Weighted 1/#it{p}_{T} d#it{N}/d#it{p}_{T}");
    frameU->GetXaxis()->SetLabelSize(0);
    frameU->GetXaxis()->SetTitleSize(0);
    frameU->GetYaxis()->SetTitleSize(0.050);
    frameU->GetYaxis()->SetLabelSize(0.040);
    frameU->GetYaxis()->SetTitleOffset(1.05);
    frameU->GetXaxis()->SetRangeUser(0, 60);
    double ymax = h->GetMaximum() * 50;
    double ymin = 1e-1;
    frameU->GetYaxis()->SetRangeUser(ymin, ymax);
    frameU->Draw("AXIS");
    h    ->Draw("E SAME");
    f_pow->Draw("SAME");

    // Fit range indicators
    TLine *lnLo = new TLine(fit_lo, ymin, fit_lo, ymax);
    TLine *lnHi = new TLine(fit_hi, ymin, fit_hi, ymax);
    lnLo->SetLineColor(kGray + 2); lnLo->SetLineStyle(3);
    lnHi->SetLineColor(kGray + 2); lnHi->SetLineStyle(3);
    lnLo->Draw(); lnHi->Draw();

    // Efficient-range band
    TBox *effBox = new TBox(eff_lo_pT, ymin, eff_hi_pT, ymax);
    effBox->SetFillColorAlpha(kGreen - 9, 0.25);
    effBox->SetLineColor(0);
    effBox->Draw("SAME");
    h    ->Draw("E SAME");
    f_pow->Draw("SAME");
    frameU->Draw("AXIS SAME");

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.045);
    tex.DrawLatex(0.55, 0.87, strleg1.c_str());
    tex.SetTextSize(0.038);
    tex.DrawLatex(0.55, 0.82,
                  "#it{p}+#it{p} #kern[-0.05]{#sqrt{#it{s}} = 200 GeV, 48.9 pb^{-1}}");
    tex.DrawLatex(0.55, 0.77,
                  (strSigMC + ", leading jet/evt").c_str());
    tex.DrawLatex(0.55, 0.72, "anti-#it{k}_{T} #it{R}=0.4, unsubtracted");

    TLegend *leg = new TLegend(0.53, 0.46, 0.94, 0.68);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
    leg->AddEntry(h,     "h_max_reco_jet_pT (1 GeV)", "lep");
    leg->AddEntry(f_pow, Form("Modified Hagedorn, #chi^{2}/ndf=%.1f", ndf>0?chi2/ndf:0), "l");
    leg->AddEntry(effBox, Form("Efficient: %.0f #minus %.0f GeV",
                               eff_lo_pT, eff_hi_pT), "f");
    leg->Draw();

    tex.SetTextSize(0.034);
    tex.DrawLatex(0.20, 0.30,
                  Form("Fit range: %.0f #minus %.0f GeV", fit_lo, fit_hi));
    tex.DrawLatex(0.20, 0.25,
                  Form("Local n(20) = %.1f, n(30) = %.1f", n20, n30));
    tex.DrawLatex(0.20, 0.20, Form("#chi^{2}/ndf = %.1f", ndf>0?chi2/ndf:0));

    pL->cd();
    TH1F *frameL = new TH1F("frameL_reco_jet_pT", "", 100, 0, 60);
    frameL->SetXTitle("#it{p}_{T}^{jet,rec} [GeV]");
    frameL->SetYTitle("data / fit");
    frameL->GetXaxis()->SetRangeUser(0, 60);
    frameL->GetYaxis()->SetRangeUser(0.45, 1.45);
    frameL->GetYaxis()->SetNdivisions(505);
    frameL->GetXaxis()->SetLabelSize(0.09);
    frameL->GetXaxis()->SetTitleSize(0.10);
    frameL->GetXaxis()->SetTitleOffset(1.15);
    frameL->GetYaxis()->SetLabelSize(0.085);
    frameL->GetYaxis()->SetTitleSize(0.095);
    frameL->GetYaxis()->SetTitleOffset(0.55);
    frameL->Draw("AXIS");

    TBox *effBoxL = new TBox(eff_lo_pT, 0.45, eff_hi_pT, 1.45);
    effBoxL->SetFillColorAlpha(kGreen - 9, 0.25);
    effBoxL->SetLineColor(0);
    effBoxL->Draw("SAME");

    TLine *l1 = new TLine(0, 1.0, 60, 1.0);
    l1->SetLineColor(kRed + 1); l1->SetLineStyle(2); l1->SetLineWidth(2);
    l1->Draw();
    TLine *lHi = new TLine(0, 1 + eff_tol, 60, 1 + eff_tol);
    TLine *lLo = new TLine(0, 1 - eff_tol, 60, 1 - eff_tol);
    lHi->SetLineColor(kGray + 1); lHi->SetLineStyle(3); lHi->SetLineWidth(2); lHi->Draw();
    lLo->SetLineColor(kGray + 1); lLo->SetLineStyle(3); lLo->SetLineWidth(2); lLo->Draw();

    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerSize(0.6);
    h_ratio->SetMarkerColor(kBlack);
    h_ratio->SetLineColor(kBlack);
    h_ratio->Draw("E SAME");
    frameL->Draw("AXIS SAME");

    gSystem->Exec("mkdir -p figures");
    c->SaveAs("figures/reco_jet_pT_signal_MC_fit.pdf");
    std::cout << "\nSaved figures/reco_jet_pT_signal_MC_fit.pdf" << std::endl;
}
