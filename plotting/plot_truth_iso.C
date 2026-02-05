#include "plotcommon.h"

void plot_truth_iso()
{
    init_plot();

    string savePath = "figures/";

    // Open input files
    TFile* f_sig = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_signal.root", "READ");
    TFile* f_bkg = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_jet.root", "READ");

    if (!f_sig || f_sig->IsZombie()) {
        std::cout << "Error: Could not open signal file!" << std::endl;
        return;
    }
    if (!f_bkg || f_bkg->IsZombie()) {
        std::cout << "Error: Could not open background file!" << std::endl;
        return;
    }

    // Arrays to store histograms for each pT bin
    TH1D* h_sig_truth_pt[NptBins];
    TH1D* h_bkg_truth_pt[NptBins];
    TH1D* h_sig_reco_pt[NptBins];
    TH1D* h_bkg_reco_pt[NptBins];

    // Load and process histograms for all pT bins
    for(int ipt=0; ipt<NptBins; ++ipt) {
        // Load histograms
        h_sig_truth_pt[ipt] = (TH1D*) f_sig->Get(Form("h_signal_truth_iso_eta0_pt%d", ipt));
        h_bkg_truth_pt[ipt] = (TH1D*) f_bkg->Get(Form("h_background_truth_iso_eta0_pt%d", ipt));
        h_sig_reco_pt[ipt] = (TH1D*) f_sig->Get(Form("h_tight_reco_iso_eta0_pt%d", ipt));
        h_bkg_reco_pt[ipt] = (TH1D*) f_bkg->Get(Form("h_tight_reco_iso_eta0_pt%d", ipt));

        // Check if histograms exist
        if (!h_sig_truth_pt[ipt]) {
            std::cout << "Warning: Could not load h_signal_truth_iso_eta0_pt" << ipt << std::endl;
            continue;
        }
        if (!h_bkg_truth_pt[ipt]) {
            std::cout << "Warning: Could not load h_background_truth_iso_eta0_pt" << ipt << std::endl;
            continue;
        }
        if (!h_sig_reco_pt[ipt]) {
            std::cout << "Warning: Could not load h_tight_reco_iso_eta0_pt from signal" << ipt << std::endl;
            continue;
        }
        if (!h_bkg_reco_pt[ipt]) {
            std::cout << "Warning: Could not load h_tight_reco_iso_eta0_pt from background" << ipt << std::endl;
            continue;
        }

        // Rebin to improve statistics
        h_sig_truth_pt[ipt]->Rebin(4);
        h_bkg_truth_pt[ipt]->Rebin(4);
        h_sig_reco_pt[ipt]->Rebin(4);
        h_bkg_reco_pt[ipt]->Rebin(4);

        // Normalize by integral (area = 1)
        h_sig_truth_pt[ipt]->Scale(1.0 / h_sig_truth_pt[ipt]->Integral());
        h_bkg_truth_pt[ipt]->Scale(1.0 / h_bkg_truth_pt[ipt]->Integral());
        h_sig_reco_pt[ipt]->Scale(1.0 / h_sig_reco_pt[ipt]->Integral());
        h_bkg_reco_pt[ipt]->Scale(1.0 / h_bkg_reco_pt[ipt]->Integral());

        // Scale by bin width (convert to probability density)
        h_sig_truth_pt[ipt]->Scale(1.0 / h_sig_truth_pt[ipt]->GetBinWidth(1));
        h_bkg_truth_pt[ipt]->Scale(1.0 / h_bkg_truth_pt[ipt]->GetBinWidth(1));
        h_sig_reco_pt[ipt]->Scale(1.0 / h_sig_reco_pt[ipt]->GetBinWidth(1));
        h_bkg_reco_pt[ipt]->Scale(1.0 / h_bkg_reco_pt[ipt]->GetBinWidth(1));
    }

    // Create plots for each pT bin
    for(int ipt=0; ipt<NptBins; ++ipt) {
        if (!h_sig_truth_pt[ipt] || !h_bkg_truth_pt[ipt] || !h_sig_reco_pt[ipt] || !h_bkg_reco_pt[ipt]) continue;

        // ============================================================
        // Plot A: Signal truth iso vs Background truth iso
        // ============================================================
        TCanvas* c1 = new TCanvas(Form("c1_pt%d", ipt), Form("c1_pt%d", ipt), 600, 600);
        frame_isoET->Draw("axis");

        float max1 = std::max(h_sig_truth_pt[ipt]->GetMaximum(),
                              h_bkg_truth_pt[ipt]->GetMaximum());
	frame_isoET->SetXTitle("#it{E}_{T}^{iso, truth} [GeV]");
        frame_isoET->GetYaxis()->SetRangeUser(0, max1 * 1.3);
        frame_isoET->SetYTitle("Arb. Unit");

        // Draw histograms
        h_bkg_truth_pt[ipt]->SetLineColor(kRed);
        h_bkg_truth_pt[ipt]->SetMarkerColor(kRed);
        h_bkg_truth_pt[ipt]->SetLineWidth(2);
        h_bkg_truth_pt[ipt]->Draw("same hist");
        h_bkg_truth_pt[ipt]->Draw("same ex0");
        h_bkg_truth_pt[ipt]->SetMarkerSize(0);

        h_sig_truth_pt[ipt]->SetLineColor(kBlue);
        h_sig_truth_pt[ipt]->SetMarkerColor(kBlue);
        h_sig_truth_pt[ipt]->SetLineWidth(2);
        h_sig_truth_pt[ipt]->Draw("same hist");

        // Add labels
        myText(0.2, 0.90, 1, Form("%0.0f < #it{E}_{T}^{#gamma,rec} < %0.0f [GeV]",
               ptRanges[ipt], ptRanges[ipt+1]), 0.04);

        myText(0.50, 0.90-0.05, 1, strleg1.c_str(), 0.04);
        myText(0.50, 0.85-0.05, 1, strleg2.c_str(), 0.04);
        myText(0.50, 0.80-0.05, 1, strleg3.c_str(), 0.04);

        myMarkerLineText(0.55, 0.75-0.05, 0, kBlue, 0, kBlue, 1,
                        "Tight Signal", 0.05, true);
        myMarkerLineText(0.55, 0.70-0.05, 0, kRed, 0, kRed, 1,
                        "Tight Background", 0.05, true);

        c1->SaveAs(Form("%s/truth_iso_sig_vs_bkg_pt%d.pdf", savePath.c_str(), ipt));
        delete c1;

        // ============================================================
        // Plot B: Signal truth iso vs Signal reco iso
        // ============================================================
        TCanvas* c2 = new TCanvas(Form("c2_pt%d", ipt), Form("c2_pt%d", ipt), 600, 600);
        frame_isoET->Draw("axis");

        float max2 = std::max(h_sig_truth_pt[ipt]->GetMaximum(),
                              h_sig_reco_pt[ipt]->GetMaximum());
        frame_isoET->GetYaxis()->SetRangeUser(0, max2 * 1.3);
        frame_isoET->SetYTitle("Probability / GeV");

        // Draw histograms
        h_sig_reco_pt[ipt]->SetLineColor(kRed);
        h_sig_reco_pt[ipt]->SetMarkerColor(kRed);
        h_sig_reco_pt[ipt]->SetLineWidth(2);
        h_sig_reco_pt[ipt]->Draw("same hist");
        h_sig_reco_pt[ipt]->Draw("same ex0");
        h_sig_reco_pt[ipt]->SetMarkerSize(0);

        h_sig_truth_pt[ipt]->SetLineColor(kBlue);
        h_sig_truth_pt[ipt]->SetMarkerColor(kBlue);
        h_sig_truth_pt[ipt]->SetLineWidth(2);
        h_sig_truth_pt[ipt]->Draw("same hist");

        // Add labels
        myText(0.2, 0.90, 1, Form("%0.0f < #it{E}_{T}^{#gamma,rec} < %0.0f [GeV]",
               ptRanges[ipt], ptRanges[ipt+1]), 0.04);

        myText(0.50, 0.90-0.05, 1, strleg1.c_str(), 0.04);
        myText(0.50, 0.85-0.05, 1, strleg2.c_str(), 0.04);
        myText(0.50, 0.80-0.05, 1, strleg3.c_str(), 0.04);
        myText(0.50, 0.75-0.05, 1, strSigMC.c_str(), 0.04);

        myMarkerLineText(0.55, 0.70-0.05, 0, kBlue, 0, kBlue, 1,
                        "Truth isolation", 0.05, true);
        myMarkerLineText(0.55, 0.65-0.05, 0, kRed, 0, kRed, 1,
                        "Reco isolation", 0.05, true);

        c2->SaveAs(Form("%s/truth_iso_sig_truth_vs_reco_pt%d.pdf", savePath.c_str(), ipt));
        delete c2;

        // ============================================================
        // Plot C: Signal reco iso vs Background reco iso
        // ============================================================
        TCanvas* c3 = new TCanvas(Form("c3_pt%d", ipt), Form("c3_pt%d", ipt), 600, 600);
        frame_isoET->Draw("axis");

        float max3 = std::max(h_sig_reco_pt[ipt]->GetMaximum(),
                              h_bkg_reco_pt[ipt]->GetMaximum());
        frame_isoET->SetXTitle("#it{E}_{T}^{iso, reco} [GeV]");
        frame_isoET->GetYaxis()->SetRangeUser(0, max3 * 1.3);
        frame_isoET->SetYTitle("Arb. Unit");

        // Draw histograms
        h_bkg_reco_pt[ipt]->SetLineColor(kRed);
        h_bkg_reco_pt[ipt]->SetMarkerColor(kRed);
        h_bkg_reco_pt[ipt]->SetLineWidth(2);
        h_bkg_reco_pt[ipt]->Draw("same hist");
        h_bkg_reco_pt[ipt]->Draw("same ex0");
        h_bkg_reco_pt[ipt]->SetMarkerSize(0);

        h_sig_reco_pt[ipt]->SetLineColor(kBlue);
        h_sig_reco_pt[ipt]->SetMarkerColor(kBlue);
        h_sig_reco_pt[ipt]->SetLineWidth(2);
        h_sig_reco_pt[ipt]->Draw("same hist");

        // Add labels
        myText(0.2, 0.90, 1, Form("%0.0f < #it{E}_{T}^{#gamma,rec} < %0.0f [GeV]",
               ptRanges[ipt], ptRanges[ipt+1]), 0.04);

        myText(0.50, 0.90-0.05, 1, strleg1.c_str(), 0.04);
        myText(0.50, 0.85-0.05, 1, strleg2.c_str(), 0.04);
        myText(0.50, 0.80-0.05, 1, strleg3.c_str(), 0.04);

        myMarkerLineText(0.55, 0.75-0.05, 0, kBlue, 0, kBlue, 1,
                        "Tight Signal", 0.05, true);
        myMarkerLineText(0.55, 0.70-0.05, 0, kRed, 0, kRed, 1,
                        "Tight Background", 0.05, true);

        c3->SaveAs(Form("%s/reco_iso_sig_vs_bkg_pt%d.pdf", savePath.c_str(), ipt));
        delete c3;
    }

    f_sig->Close();
    f_bkg->Close();

    std::cout << "Plots saved to " << savePath << std::endl;
}
