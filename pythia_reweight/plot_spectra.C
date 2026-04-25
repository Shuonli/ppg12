// plot_spectra.C
// ROOT macro: read output_reweighted.root and output_bins.root,
// compute dσ/dpT cross sections, and produce a 2×2 canvas saved as spectra.pdf
//
// Usage: root -l 'plot_spectra.C'

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

void plot_spectra()
{
    // ---------------------------------------------------------------
    // pT binning for histograms
    // ---------------------------------------------------------------
    const int    nPtBins  = 20;
    const double ptBinLo  = 5.;
    const double ptBinHi  = 55.;

    // ---------------------------------------------------------------
    // Helper: normalise histogram in-place to dσ/dpT [mb/GeV]
    //   scale_factor = sigma_gen / weight_sum
    //   each bin: content *= scale_factor / bin_width
    // ---------------------------------------------------------------
    auto normalise = [](TH1D* h, double sigma, double wsum) {
        if (wsum <= 0.) return;
        double sf = sigma / wsum;
        for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
            double bw  = h->GetBinWidth(ib);
            double val = h->GetBinContent(ib);
            double err = h->GetBinError(ib);
            h->SetBinContent(ib, val * sf / bw);
            h->SetBinError  (ib, err * sf / bw);
        }
    };

    // ---------------------------------------------------------------
    // Pass 1: output_reweighted.root
    // ---------------------------------------------------------------
    TFile* f1 = TFile::Open("output_reweighted.root", "READ");
    if (!f1 || f1->IsZombie()) {
        std::cerr << "ERROR: cannot open output_reweighted.root" << std::endl;
        return;
    }

    // Read meta: sigma_gen and weight_sum summed over both sub-runs
    // (weight for histograms is per-event; normalisation uses total sigma/wsum per sub-run)
    TTree* meta1 = (TTree*)f1->Get("meta");
    if (!meta1) { std::cerr << "ERROR: meta tree not found in pass1 file" << std::endl; return; }

    Int_t    m1_bin_id;
    Double_t m1_sigma, m1_wsum;
    meta1->SetBranchAddress("bin_id",     &m1_bin_id);
    meta1->SetBranchAddress("sigma_gen",  &m1_sigma);
    meta1->SetBranchAddress("weight_sum", &m1_wsum);

    // Store per-sub-run normalisation factors
    std::vector<double> sf1(2, 0.);   // sigma/wsum for bin 0 and bin 1
    for (Long64_t ie = 0; ie < meta1->GetEntries(); ++ie) {
        meta1->GetEntry(ie);
        if (m1_wsum > 0.) sf1[m1_bin_id] = m1_sigma / m1_wsum;
        std::cout << "Pass1 sub-run " << m1_bin_id
                  << ": sigma=" << m1_sigma << " mb, wsum=" << m1_wsum
                  << ", sf=" << sf1[m1_bin_id] << std::endl;
    }

    // Weighted fill: accumulate (event_weight * sf_bin) per pT bin,
    // then divide by bin width in normalise step is already handled above.
    // But since each sub-run has its own sf, we fill two intermediate hists
    // and sum them.
    TH1D* h_jet_p1_s0 = new TH1D("h_jet_p1_s0","",nPtBins,ptBinLo,ptBinHi);
    TH1D* h_jet_p1_s1 = new TH1D("h_jet_p1_s1","",nPtBins,ptBinLo,ptBinHi);
    TH1D* h_pho_p1_s0 = new TH1D("h_pho_p1_s0","",nPtBins,ptBinLo,ptBinHi);
    TH1D* h_pho_p1_s1 = new TH1D("h_pho_p1_s1","",nPtBins,ptBinLo,ptBinHi);
    h_jet_p1_s0->Sumw2(); h_jet_p1_s1->Sumw2();
    h_pho_p1_s0->Sumw2(); h_pho_p1_s1->Sumw2();

    TTree* ev1 = (TTree*)f1->Get("events");
    if (!ev1) { std::cerr << "ERROR: events tree not found" << std::endl; return; }

    Int_t    ev1_bin_id;
    Double_t ev1_wt;
    std::vector<float>* ev1_jet_pt  = nullptr;
    std::vector<float>* ev1_pho_pt  = nullptr;

    ev1->SetBranchAddress("bin_id",       &ev1_bin_id);
    ev1->SetBranchAddress("event_weight", &ev1_wt);
    ev1->SetBranchAddress("jet_pt",       &ev1_jet_pt);
    ev1->SetBranchAddress("photon_pt",    &ev1_pho_pt);

    for (Long64_t ie = 0; ie < ev1->GetEntries(); ++ie) {
        ev1->GetEntry(ie);
        TH1D* hj = (ev1_bin_id == 0) ? h_jet_p1_s0 : h_jet_p1_s1;
        TH1D* hp = (ev1_bin_id == 0) ? h_pho_p1_s0 : h_pho_p1_s1;
        for (float pt : *ev1_jet_pt) hj->Fill(pt, ev1_wt);
        for (float pt : *ev1_pho_pt) hp->Fill(pt, ev1_wt);
    }

    // Normalise each sub-run histogram then sum
    for (int ib = 0; ib < 2; ++ib) {
        TH1D* hj = (ib == 0) ? h_jet_p1_s0 : h_jet_p1_s1;
        TH1D* hp = (ib == 0) ? h_pho_p1_s0 : h_pho_p1_s1;
        // Scale by sf (sigma/wsum) and divide by bin width
        for (int k = 1; k <= hj->GetNbinsX(); ++k) {
            double bw = hj->GetBinWidth(k);
            hj->SetBinContent(k, hj->GetBinContent(k) * sf1[ib] / bw);
            hj->SetBinError  (k, hj->GetBinError(k)   * sf1[ib] / bw);
            hp->SetBinContent(k, hp->GetBinContent(k) * sf1[ib] / bw);
            hp->SetBinError  (k, hp->GetBinError(k)   * sf1[ib] / bw);
        }
    }

    TH1D* h_jet_pass1 = (TH1D*)h_jet_p1_s0->Clone("h_jet_xsec_pass1");
    h_jet_pass1->SetTitle("Jet p_{T} cross section (Pass 1);p_{T} [GeV];d#sigma/dp_{T} [mb/GeV]");
    h_jet_pass1->Add(h_jet_p1_s1);

    TH1D* h_pho_pass1 = (TH1D*)h_pho_p1_s0->Clone("h_photon_xsec_pass1");
    h_pho_pass1->SetTitle("Photon p_{T} cross section (Pass 1);p_{T} [GeV];d#sigma/dp_{T} [mb/GeV]");
    h_pho_pass1->Add(h_pho_p1_s1);

    f1->Close();

    // ---------------------------------------------------------------
    // Pass 2: output_bins.root
    // ---------------------------------------------------------------
    TFile* f2 = TFile::Open("output_bins.root", "READ");
    if (!f2 || f2->IsZombie()) {
        std::cerr << "ERROR: cannot open output_bins.root" << std::endl;
        return;
    }

    TTree* meta2 = (TTree*)f2->Get("meta");
    if (!meta2) { std::cerr << "ERROR: meta tree not found in pass2 file" << std::endl; return; }

    Double_t m2_binlo, m2_binhi, m2_sigma, m2_wsum;
    meta2->SetBranchAddress("bin_low",    &m2_binlo);
    meta2->SetBranchAddress("bin_high",   &m2_binhi);
    meta2->SetBranchAddress("sigma_gen",  &m2_sigma);
    meta2->SetBranchAddress("weight_sum", &m2_wsum);

    int nBins2 = (int)meta2->GetEntries();

    TH1D* h_jet_pass2 = new TH1D("h_jet_xsec_pass2",
        "Jet p_{T} cross section (Pass 2);p_{T} [GeV];d#sigma/dp_{T} [mb/GeV]",
        nPtBins, ptBinLo, ptBinHi);
    TH1D* h_pho_pass2 = new TH1D("h_photon_xsec_pass2",
        "Photon p_{T} cross section (Pass 2);p_{T} [GeV];d#sigma/dp_{T} [mb/GeV]",
        nPtBins, ptBinLo, ptBinHi);
    h_jet_pass2->Sumw2();
    h_pho_pass2->Sumw2();

    for (int iBin = 0; iBin < nBins2; ++iBin) {
        meta2->GetEntry(iBin);
        double sf2 = (m2_wsum > 0.) ? m2_sigma / m2_wsum : 0.;
        std::cout << "Pass2 bin " << iBin
                  << " [" << m2_binlo << "-" << m2_binhi << "] GeV"
                  << ": sigma=" << m2_sigma << " mb, wsum=" << m2_wsum
                  << ", sf=" << sf2 << std::endl;

        std::string tname = "events_bin" + std::to_string(iBin);
        TTree* ev2 = (TTree*)f2->Get(tname.c_str());
        if (!ev2) {
            std::cerr << "WARNING: tree " << tname << " not found, skipping" << std::endl;
            continue;
        }

        std::vector<float>* ev2_jet_pt = nullptr;
        std::vector<float>* ev2_pho_pt = nullptr;
        ev2->SetBranchAddress("jet_pt",    &ev2_jet_pt);
        ev2->SetBranchAddress("photon_pt", &ev2_pho_pt);

        // Temporary histograms for this bin
        TH1D* htj = new TH1D(Form("htj%d",iBin),"",nPtBins,ptBinLo,ptBinHi);
        TH1D* htp = new TH1D(Form("htp%d",iBin),"",nPtBins,ptBinLo,ptBinHi);
        htj->Sumw2(); htp->Sumw2();

        for (Long64_t ie = 0; ie < ev2->GetEntries(); ++ie) {
            ev2->GetEntry(ie);
            for (float pt : *ev2_jet_pt) htj->Fill(pt);
            for (float pt : *ev2_pho_pt) htp->Fill(pt);
        }

        // Normalise by sf2 / bin_width and add to total
        for (int k = 1; k <= htj->GetNbinsX(); ++k) {
            double bw = htj->GetBinWidth(k);
            h_jet_pass2->SetBinContent(k, h_jet_pass2->GetBinContent(k) + htj->GetBinContent(k)*sf2/bw);
            h_jet_pass2->SetBinError  (k, std::hypot(h_jet_pass2->GetBinError(k), htj->GetBinError(k)*sf2/bw));
            h_pho_pass2->SetBinContent(k, h_pho_pass2->GetBinContent(k) + htp->GetBinContent(k)*sf2/bw);
            h_pho_pass2->SetBinError  (k, std::hypot(h_pho_pass2->GetBinError(k), htp->GetBinError(k)*sf2/bw));
        }

        delete htj;
        delete htp;
    }

    f2->Close();

    // ---------------------------------------------------------------
    // Style
    // ---------------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);

    auto styleHist = [](TH1D* h, Color_t col, Style_t mst) {
        h->SetLineColor(col);
        h->SetMarkerColor(col);
        h->SetMarkerStyle(mst);
        h->SetMarkerSize(0.8);
        h->GetXaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetTitleSize(0.05);
        h->GetXaxis()->SetLabelSize(0.04);
        h->GetYaxis()->SetLabelSize(0.04);
    };

    styleHist(h_jet_pass1,  kBlue+1,   20);
    styleHist(h_jet_pass2,  kRed+1,    21);
    styleHist(h_pho_pass1,  kGreen+2,  22);
    styleHist(h_pho_pass2,  kMagenta+1,23);

    // ---------------------------------------------------------------
    // Draw 2×2 canvas
    // ---------------------------------------------------------------
    TCanvas* cv = new TCanvas("cv", "Cross Section Spectra", 1200, 1000);
    cv->Divide(2, 2, 0.01, 0.01);

    // --- pad 1: jet pass1 ---
    cv->cd(1);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    if (h_jet_pass1->GetMaximum() > 0) {
        h_jet_pass1->SetMinimum(1e-12);
        h_jet_pass1->Draw("E");
    }

    // --- pad 2: jet pass2 ---
    cv->cd(2);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    if (h_jet_pass2->GetMaximum() > 0) {
        h_jet_pass2->SetMinimum(1e-12);
        h_jet_pass2->Draw("E");
    }

    // --- pad 3: photon pass1 ---
    cv->cd(3);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    if (h_pho_pass1->GetMaximum() > 0) {
        h_pho_pass1->SetMinimum(1e-12);
        h_pho_pass1->Draw("E");
    }

    // --- pad 4: photon pass2 ---
    cv->cd(4);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    if (h_pho_pass2->GetMaximum() > 0) {
        h_pho_pass2->SetMinimum(1e-12);
        h_pho_pass2->Draw("E");
    }

    cv->SaveAs("spectra.pdf");
    std::cout << "Saved spectra.pdf" << std::endl;

    // ---------------------------------------------------------------
    // Also draw pass1 vs pass2 overlaid for comparison
    // ---------------------------------------------------------------
    TCanvas* cv2 = new TCanvas("cv2","Comparison Pass1 vs Pass2", 1200, 500);
    cv2->Divide(2, 1, 0.01, 0.01);

    cv2->cd(1);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    double jetMax = std::max(h_jet_pass1->GetMaximum(), h_jet_pass2->GetMaximum());
    double jetMin = 1e-12;
    if (h_jet_pass1->GetMaximum() > 0) {
        h_jet_pass1->SetMinimum(jetMin);
        h_jet_pass1->SetMaximum(jetMax * 10.);
        h_jet_pass1->Draw("E");
    }
    if (h_jet_pass2->GetMaximum() > 0) {
        h_jet_pass2->SetMinimum(jetMin);
        h_jet_pass2->Draw("E SAME");
    }
    TLegend* leg1 = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg1->AddEntry(h_jet_pass1, "Pass 1 (reweighted)", "lpe");
    leg1->AddEntry(h_jet_pass2, "Pass 2 (bins)",       "lpe");
    leg1->SetBorderSize(0);
    leg1->Draw();

    cv2->cd(2);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    double phoMax = std::max(h_pho_pass1->GetMaximum(), h_pho_pass2->GetMaximum());
    double phoMin = 1e-12;
    if (h_pho_pass1->GetMaximum() > 0) {
        h_pho_pass1->SetMinimum(phoMin);
        h_pho_pass1->SetMaximum(phoMax * 10.);
        h_pho_pass1->Draw("E");
    }
    if (h_pho_pass2->GetMaximum() > 0) {
        h_pho_pass2->SetMinimum(phoMin);
        h_pho_pass2->Draw("E SAME");
    }
    TLegend* leg2 = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg2->AddEntry(h_pho_pass1, "Pass 1 (reweighted)", "lpe");
    leg2->AddEntry(h_pho_pass2, "Pass 2 (bins)",       "lpe");
    leg2->SetBorderSize(0);
    leg2->Draw();

    cv2->SaveAs("spectra_comparison.pdf");
    std::cout << "Saved spectra_comparison.pdf" << std::endl;
}
