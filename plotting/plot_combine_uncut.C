#include "plotcommon.h"

// Regenerates the photon-sample turn-on Fig 3 (`10_over_05.pdf`,
// `20_over_10.pdf`) and stitched-spectrum Fig 4 (`combine.pdf`) using
// UNRESTRICTED leading-truth-photon-pT distributions (i.e., without
// the per-sample event cut applied in the analysis pipeline). This
// recovers the natural turn-on shape per sample in their overlap
// regions and yields a clean stitched spectrum for the Hagedorn fit.
//
// Inputs:
//   plotting/photon_max_pT_uncut.root  (built by
//   plotting/build_h_max_photon_pT_uncut.py)
//
// Outputs:
//   plotting/figures/10_over_05.pdf
//   plotting/figures/20_over_10.pdf
//   plotting/figures/combine.pdf

void plot_combine_uncut()
{
    init_plot();

    TFile *fin = TFile::Open(
        "/sphenix/user/shuhangli/ppg12/plotting/photon_max_pT_uncut.root",
        "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open photon_max_pT_uncut.root --"
                  << " run plotting/build_h_max_photon_pT_uncut.py first."
                  << std::endl;
        return;
    }

    TH1D *h5  = (TH1D *)fin->Get("h_max_photon_pT_photon5");
    TH1D *h10 = (TH1D *)fin->Get("h_max_photon_pT_photon10");
    TH1D *h20 = (TH1D *)fin->Get("h_max_photon_pT_photon20");
    if (!h5 || !h10 || !h20) {
        std::cerr << "ERROR: input histograms missing.\n";
        return;
    }
    h5->Sumw2();  h10->Sumw2();  h20->Sumw2();

    // --------- ratios (Fig 3) ---------
    TH1D *h10_over_05 = (TH1D *)h10->Clone("h10_over_05");
    h10_over_05->Divide(h5);
    TH1D *h20_over_10 = (TH1D *)h20->Clone("h20_over_10");
    h20_over_10->Divide(h10);

    {
        TCanvas *c = new TCanvas("c_10_over_05", "", 600, 600);
        h10_over_05->GetXaxis()->SetRangeUser(10, 32);
        h10_over_05->GetYaxis()->SetRangeUser(0.6, 1.4);
        h10_over_05->SetYTitle("photon10/photon5");
        h10_over_05->SetXTitle("Leading #it{E}_{T}^{#gamma} [GeV]");
        h10_over_05->SetMarkerStyle(20);
        h10_over_05->SetMarkerColor(kAzure + 7);
        h10_over_05->SetLineColor(kAzure + 7);
        h10_over_05->Draw("PE");
        lineone->Draw("L");
        myText(0.20, 0.88, 1, strleg1.c_str(), 0.040, 0);
        myText(0.20, 0.83, 1, strleg2.c_str(), 0.040, 0);
        myText(0.88, 0.88, 1, strSigMC.c_str(), 0.040, 1);
        c->SaveAs("figures/10_over_05.pdf");
    }
    {
        TCanvas *c = new TCanvas("c_20_over_10", "", 600, 600);
        h20_over_10->GetXaxis()->SetRangeUser(20, 40);
        h20_over_10->GetYaxis()->SetRangeUser(0.6, 1.4);
        h20_over_10->SetYTitle("photon20/photon10");
        h20_over_10->SetXTitle("Leading #it{E}_{T}^{#gamma} [GeV]");
        h20_over_10->SetMarkerStyle(20);
        h20_over_10->SetMarkerColor(kPink + 8);
        h20_over_10->SetLineColor(kPink + 8);
        h20_over_10->Draw("PE");
        lineone->Draw("L");
        myText(0.20, 0.88, 1, strleg1.c_str(), 0.040, 0);
        myText(0.20, 0.83, 1, strleg2.c_str(), 0.040, 0);
        myText(0.88, 0.88, 1, strSigMC.c_str(), 0.040, 1);
        c->SaveAs("figures/20_over_10.pdf");
    }

    // --------- stitched spectrum + Hagedorn fit (Fig 4) ---------
    // Use each sample only in its CANONICAL pT window so the stitched
    // sum reproduces the physical spectrum without double counting:
    //   photon5  : pT in [0, 14]
    //   photon10 : pT in [14, 22]
    //   photon20 : pT in [22, infty)
    TH1D *h5_w  = (TH1D *)h5->Clone("h5_w");
    TH1D *h10_w = (TH1D *)h10->Clone("h10_w");
    TH1D *h20_w = (TH1D *)h20->Clone("h20_w");
    for (int i = 1; i <= h5_w->GetNbinsX(); ++i) {
        double c = h5_w->GetBinCenter(i);
        if (c >= 14.0) { h5_w->SetBinContent(i, 0); h5_w->SetBinError(i, 0); }
    }
    for (int i = 1; i <= h10_w->GetNbinsX(); ++i) {
        double c = h10_w->GetBinCenter(i);
        if (c < 14.0 || c >= 22.0) { h10_w->SetBinContent(i, 0); h10_w->SetBinError(i, 0); }
    }
    for (int i = 1; i <= h20_w->GetNbinsX(); ++i) {
        double c = h20_w->GetBinCenter(i);
        if (c < 22.0) { h20_w->SetBinContent(i, 0); h20_w->SetBinError(i, 0); }
    }
    TH1D *h_sum = (TH1D *)h5_w->Clone("h_sum_stitched");
    h_sum->Add(h10_w);
    h_sum->Add(h20_w);

    // Hagedorn-style modified power-law fit
    TF1 *f1 = new TF1("f1", "[0]*pow([1]/x,[2]+[3]*log(x/[1])+ [4]*x)", 10, 50);
    f1->SetParameters(2.09375e+9, 1.0, 1.0, 2.0, 0.01);
    h_sum->Fit(f1, "REMNQ", "", 10, 36);
    h_sum->Fit(f1, "REMNQ", "", 10, 36);

    TH1D *h_ratio = (TH1D *)h_sum->Clone("h_ratio");
    h_ratio->Divide(f1);

    TCanvas *c4 = new TCanvas("can_combine", "", 800, 889);
    c4->Divide(1, 2);

    TPad *p1 = (TPad *)c4->cd(1);
    p1->SetPad(0, 0.4, 1, 1);
    p1->SetTopMargin(0.05);
    p1->SetLeftMargin(0.13);
    p1->SetBottomMargin(0.002);
    p1->SetRightMargin(0.08);
    p1->SetLogy();

    // Auto-range Y based on data: max from h5 low-pT tail, min from h20 high-pT tail.
    double y_max = h5_w->GetMaximum();
    double y_min = h20_w->GetBinContent(h20_w->FindBin(35.0));
    if (y_min <= 0) y_min = 1e-4;
    frame_et_rec->SetYTitle("d#sigma / d#it{E}_{T}^{#gamma} [pb / GeV]");
    frame_et_rec->GetYaxis()->SetRangeUser(y_min * 0.3, y_max * 5);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 40);
    frame_et_rec->GetXaxis()->SetTitleOffset(1.05);
    frame_et_rec->GetYaxis()->SetTitleOffset(1.05);
    frame_et_rec->GetYaxis()->SetTitleSize(0.053);
    frame_et_rec->GetXaxis()->SetLabelSize(0.050);
    frame_et_rec->GetYaxis()->SetLabelSize(0.050);
    frame_et_rec->GetXaxis()->SetLabelOffset(2);
    frame_et_rec->GetXaxis()->SetNdivisions(505);
    frame_et_rec->Draw("axis");

    h5_w->SetMarkerStyle(20);
    h5_w->SetMarkerColor(kPink + 5);
    h5_w->SetLineColor(kPink + 5);
    h5_w->Draw("P same");

    h10_w->SetMarkerStyle(20);
    h10_w->SetMarkerColor(kGreen - 2);
    h10_w->SetLineColor(kGreen - 2);
    h10_w->Draw("P same");

    h20_w->SetMarkerStyle(20);
    h20_w->SetMarkerColor(kAzure + 7);
    h20_w->SetLineColor(kAzure + 7);
    h20_w->Draw("P same");

    f1->SetLineColor(kRed - 4);
    f1->SetLineWidth(2);
    f1->Draw("same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, strMC.c_str(), 0.05);

    myMarkerLineText(0.6, 0.75, 1, kPink + 5, 20, kPink + 5, 1, "photon5",  0.05, true);
    myMarkerLineText(0.6, 0.70, 1, kGreen - 2, 20, kGreen - 2, 1, "photon10", 0.05, true);
    myMarkerLineText(0.6, 0.65, 1, kAzure + 7, 20, kAzure + 7, 1, "photon20", 0.05, true);

    TPad *p2 = (TPad *)c4->cd(2);
    p2->SetPad(0, 0, 1, 0.4);
    p2->SetTopMargin(0.023);
    p2->SetLeftMargin(0.13);
    p2->SetBottomMargin(0.25);
    p2->SetRightMargin(0.08);

    frame_et_truth->SetYTitle("Data / Fit");
    frame_et_truth->SetXTitle("Leading #it{E}_{T}^{#gamma} [GeV]");
    frame_et_truth->GetYaxis()->SetNdivisions(506);
    frame_et_truth->GetYaxis()->SetRangeUser(0.85, 1.15);
    frame_et_truth->GetXaxis()->SetRangeUser(10, 40);
    frame_et_truth->GetXaxis()->SetTitleOffset(frame_et_rec->GetXaxis()->GetTitleOffset() * 4 / 6.* 1.4);
    frame_et_truth->GetYaxis()->SetTitleOffset(frame_et_rec->GetYaxis()->GetTitleOffset() * 4 / 6.);
    frame_et_truth->GetYaxis()->SetLabelOffset(frame_et_rec->GetYaxis()->GetLabelOffset() * 4 / 6.);
    frame_et_truth->GetXaxis()->SetLabelSize(frame_et_rec->GetXaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetYaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetTitleSize(frame_et_rec->GetXaxis()->GetTitleSize() * 6 / 4. * 1.2);
    frame_et_truth->GetYaxis()->SetTitleSize(frame_et_rec->GetYaxis()->GetTitleSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetNdivisions(505);
    frame_et_truth->Draw("axis");

    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerColor(kBlack);
    h_ratio->SetLineColor(kBlack);
    h_ratio->Draw("PE same");
    lineone->Draw("L");

    c4->SaveAs("figures/combine.pdf");
    fin->Close();
}
