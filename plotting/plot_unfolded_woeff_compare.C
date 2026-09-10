// plot_unfolded_woeff_compare.C
//
// Compare the unfolded BG-subtracted yield BEFORE efficiency division
// (h_unfold_sub_result_woeff) for three pipelines in the 1.5 mrad analysis:
//   PYTHIA  with the proper 1.5-mrad-derived reweight (apples-to-apples)
//   HERWIG  with its own iter=1-derived reweight ("proper")
//   HERWIG  reweighted to PYTHIA's iter=1 ("common target") — this isolates
//           the truly irreducible response-matrix difference (UE/FSR).
//
// Style follows the main-analysis plotting convention via plotcommon.h.
//
// Output: plotting/figures/unfolded_woeff_compare_1p5mrad.pdf

#include "plotcommon.h"

void plot_unfolded_woeff_compare()
{
    init_plot();

    const std::string base = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/";

    TFile *f_pythia = new TFile((base + "Photon_final_bdt_nom_1p5mrad_newrewt.root").c_str(), "READ");
    TFile *f_h_common = new TFile((base + "Photon_final_bdt_herwig_1p5mrad_common.root").c_str(), "READ");

    if (!f_pythia || f_pythia->IsZombie() || !f_h_common || f_h_common->IsZombie()) {
        std::cerr << "ERROR: cannot open inputs" << std::endl;
        return;
    }

    TH1F *h_pyt = (TH1F *)f_pythia->Get("h_unfold_sub_result_woeff");
    TH1F *h_hc  = (TH1F *)f_h_common->Get("h_unfold_sub_result_woeff");

    h_pyt = (TH1F *)h_pyt->Clone("h_pyt");
    h_hc  = (TH1F *)h_hc ->Clone("h_hc");

    // 1/Δη normalization (same as main analysis)
    const float deta = 1.4;
    h_pyt->Scale(1.0 / deta);
    h_hc ->Scale(1.0 / deta);

    // Ratio HERWIG/PYTHIA (woeff). The DATA input to both pipelines is bit-
    // identical (h_data_sub matches across runs), so the data-statistical
    // fluctuation cancels in the ratio. The naive Divide() error treats the
    // two as independent and inflates the ratio uncertainty by ~10x. We
    // override the ratio errors with the proper Gaussian MC-stat from the
    // response-matrix truth-marginal Sumw2: for each truth bin t the MC fills
    // are weighted (cross_section / N_gen) and combine across multiple samples
    // (photon10 + photon20 + jets), so per-bin variance = sum_i w_i^2 (not
    // Poisson). We pull rel_err = bin_err / bin_content from the truth-marginal
    // projection of h_response_full_0 in each pipeline and combine in
    // quadrature for the ratio.
    TH1F *r_hc = (TH1F *)h_hc->Clone("r_hc");
    r_hc->Divide(h_pyt);

    TH2F *Rp = (TH2F *)f_pythia->Get("h_response_full_0");
    TH2F *Rh = (TH2F *)f_h_common->Get("h_response_full_0");
    TH1D *mp = Rp->ProjectionY("mp_truth_marg");
    TH1D *mh = Rh->ProjectionY("mh_truth_marg");

    for (int i = 1; i <= r_hc->GetNbinsX(); ++i) {
        double v  = r_hc->GetBinContent(i);
        if (v <= 0) continue;
        double pt_c = r_hc->GetBinCenter(i);
        // Look up the truth-marginal bin in the response matrix
        int jp = mp->FindBin(pt_c);
        int jh = mh->FindBin(pt_c);
        double Np = mp->GetBinContent(jp), Ep = mp->GetBinError(jp);
        double Nh = mh->GetBinContent(jh), Eh = mh->GetBinError(jh);
        if (Np <= 0 || Nh <= 0) continue;
        double rel_p = Ep / Np;
        double rel_h = Eh / Nh;
        double rel_ratio = std::sqrt(rel_p * rel_p + rel_h * rel_h);
        r_hc->SetBinError(i, rel_ratio * v);
    }

    // Color palette (same convention as plot_final_selection.C)
    const int col_p   = kAzure + 2;       // PYTHIA blue
    const int col_hc  = kSpring - 7;      // HERWIG common-target — green

    // ============================================================
    // Canvas with 2 pads (top: spectra log-y, bottom: ratio H/P)
    // ============================================================
    TCanvas *c1 = new TCanvas("can_woeff_cmp", "", 800, 900);
    c1->Divide(1, 2);

    // -- top pad: woeff spectra
    TPad *pad_1 = (TPad *)c1->cd(1);
    pad_1->SetPad(0, 0.36, 1, 1);
    pad_1->SetTopMargin(0.06);
    pad_1->SetLeftMargin(0.13);
    pad_1->SetBottomMargin(0.002);
    pad_1->SetRightMargin(0.05);
    pad_1->SetLogy();

    frame_et_rec->SetYTitle("d^{2}#sigma/dp_{T}d#eta  (woeff) [pb / GeV]");
    frame_et_rec->GetXaxis()->SetRangeUser(10, 36);
    frame_et_rec->GetYaxis()->SetRangeUser(5e-2, 1e3);
    frame_et_rec->GetXaxis()->SetTitleOffset(1.05);
    frame_et_rec->GetYaxis()->SetTitleOffset(1.10);
    frame_et_rec->GetYaxis()->SetTitleSize(0.053);
    frame_et_rec->GetXaxis()->SetLabelSize(0.050);
    frame_et_rec->GetYaxis()->SetLabelSize(0.050);
    frame_et_rec->GetXaxis()->SetLabelOffset(2);
    frame_et_rec->GetXaxis()->SetNdivisions(505);
    frame_et_rec->Draw("axis");

    // PYTHIA
    h_pyt->SetMarkerStyle(20);
    h_pyt->SetMarkerSize(1.4);
    h_pyt->SetMarkerColor(col_p);
    h_pyt->SetLineColor(col_p);
    h_pyt->SetLineWidth(2);
    h_pyt->Draw("PE same");

    // HERWIG common-target
    h_hc->SetMarkerStyle(22);
    h_hc->SetMarkerSize(1.4);
    h_hc->SetMarkerColor(col_hc);
    h_hc->SetLineColor(col_hc);
    h_hc->SetLineWidth(2);
    h_hc->Draw("PE same");

    // header text
    myText(0.50, 0.88, 1, strleg1.c_str(), 0.05);
    myText(0.50, 0.82, 1, "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, 17.16 pb^{-1}  (1.5 mrad)", 0.04);
    myText(0.50, 0.77, 1, strleg3.c_str(), 0.04);
    myText(0.50, 0.72, 1, strleg4.c_str(), 0.04);
    myText(0.50, 0.65, 1, "Unfolded yield, eff-uncorrected", 0.04);

    myMarkerLineText(0.18, 0.26, 1, col_p,  20, col_p,  1, "PYTHIA",  0.040, true);
    myMarkerLineText(0.18, 0.20, 1, col_hc, 22, col_hc, 1, "HERWIG  (common-target rewt)",  0.040, true);

    // -- bottom pad: ratio
    TPad *pad_2 = (TPad *)c1->cd(2);
    pad_2->SetPad(0, 0, 1, 0.36);
    pad_2->SetTopMargin(0.02);
    pad_2->SetLeftMargin(0.13);
    pad_2->SetBottomMargin(0.30);
    pad_2->SetRightMargin(0.05);

    frame_et_truth->SetYTitle("HERWIG / PYTHIA");
    frame_et_truth->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
    frame_et_truth->GetYaxis()->SetNdivisions(506);
    frame_et_truth->GetYaxis()->SetRangeUser(0.7, 1.3);
    frame_et_truth->GetXaxis()->SetRangeUser(10, 36);
    frame_et_truth->GetXaxis()->SetTitleOffset(frame_et_rec->GetXaxis()->GetTitleOffset() * 4 / 6.* 1.4);
    frame_et_truth->GetYaxis()->SetTitleOffset(frame_et_rec->GetYaxis()->GetTitleOffset() * 4 / 6.);
    frame_et_truth->GetXaxis()->SetLabelSize(frame_et_rec->GetXaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetYaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetTitleSize(frame_et_rec->GetXaxis()->GetTitleSize() * 6 / 4. * 1.2);
    frame_et_truth->GetYaxis()->SetTitleSize(frame_et_rec->GetYaxis()->GetTitleSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetNdivisions(505);
    frame_et_truth->Draw("axis");

    r_hc->SetMarkerStyle(22);
    r_hc->SetMarkerSize(1.4);
    r_hc->SetMarkerColor(col_hc);
    r_hc->SetLineColor(col_hc);
    r_hc->SetLineWidth(2);
    r_hc->Draw("PE same");

    lineone->Draw("L");

    // Note about correlated data error
    TLatex lat_note;
    lat_note.SetNDC();
    lat_note.SetTextSize(0.045);
    lat_note.SetTextColor(kGray + 2);
    lat_note.DrawLatex(0.16, 0.88, "Errors: Gaussian MC-stat (Sumw2) only \x96 data corr cancels in ratio");

    c1->SaveAs("figures/unfolded_woeff_compare_1p5mrad.pdf");
    std::cout << "Saved figures/unfolded_woeff_compare_1p5mrad.pdf" << std::endl;
}
