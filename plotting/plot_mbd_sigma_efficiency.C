#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <iostream>

#include "plotcommon.h"

void plot_mbd_sigma_efficiency()
{
    init_plot();

    // Open data file
    TFile *f_data = TFile::Open("/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histoshower_shape_.root", "READ");
    if (!f_data || f_data->IsZombie())
    {
        std::cerr << "Error: Could not open data file!" << std::endl;
        return;
    }

    // pT bin edges (combine all bins)
    std::vector<double> pT_bin_edges = {10, 15, 20, 30};
    int nPtBins = pT_bin_edges.size() - 1;
    int ieta = 0;

    // Combine MBD sigma distributions from all pT bins
    TH1D *h_mbd_combined = nullptr;

    for (int ipt = 0; ipt < nPtBins; ++ipt)
    {
        TString hist_name = Form("h_mbd_avgsigma_vs_npb_score_eta%d_pt%d", ieta, ipt);
        TH2D *h2 = dynamic_cast<TH2D*>(f_data->Get(hist_name));

        if (!h2)
        {
            std::cerr << "Warning: Could not retrieve " << hist_name << std::endl;
            continue;
        }

        // Project onto Y-axis (MBD sigma)
        TH1D *h_mbd = h2->ProjectionY(Form("h_mbd_sigma_pt%d", ipt));
        if (!h_mbd || h_mbd->Integral() < 1)
        {
            std::cerr << "Warning: Empty histogram for pt bin " << ipt << std::endl;
            continue;
        }

        if (!h_mbd_combined)
        {
            h_mbd_combined = (TH1D*)h_mbd->Clone("h_mbd_combined");
            h_mbd_combined->SetDirectory(nullptr);
        }
        else
        {
            h_mbd_combined->Add(h_mbd);
        }
    }

    if (!h_mbd_combined || h_mbd_combined->Integral() < 1)
    {
        std::cerr << "Error: No data to plot!" << std::endl;
        f_data->Close();
        return;
    }

    // Create efficiency histogram
    int nbins = h_mbd_combined->GetNbinsX();
    double xmin = h_mbd_combined->GetXaxis()->GetXmin();
    double xmax = h_mbd_combined->GetXaxis()->GetXmax();

    TH1D *h_eff = new TH1D("h_eff_mbd_sigma", "", nbins, xmin, xmax);
    h_eff->SetDirectory(nullptr);

    // Calculate efficiency: fraction passing cut (MBD sigma < threshold)
    double total = h_mbd_combined->Integral(1, nbins);
    for (int ibin = 1; ibin <= nbins; ++ibin)
    {
        double passing = h_mbd_combined->Integral(1, ibin);
        double eff = (total > 0) ? passing / total : 0.0;
        double eff_err = (total > 0) ? sqrt(eff * (1 - eff) / total) : 0.0;
        h_eff->SetBinContent(ibin, eff);
        h_eff->SetBinError(ibin, eff_err);
    }

    // Create canvas
    TCanvas *c = new TCanvas("c_mbd_sigma_eff", "MBD Sigma Selection Efficiency", 600, 600);

    // Style
    h_eff->SetLineColor(kBlack);
    h_eff->SetMarkerColor(kBlack);
    h_eff->SetMarkerStyle(20);
    h_eff->SetMarkerSize(1.0);
    h_eff->SetLineWidth(2);
    h_eff->SetStats(0);
    h_eff->SetTitle("");

    h_eff->GetXaxis()->SetTitle("MBD Avg #sigma_{t} Cut [ns]");
    h_eff->GetXaxis()->SetNdivisions(505);
    h_eff->GetXaxis()->SetRangeUser(0, 5);
    h_eff->GetYaxis()->SetTitle("Selection Efficiency");
    h_eff->GetYaxis()->SetTitleOffset(1.5);
    h_eff->GetYaxis()->SetRangeUser(0, 1.1);

    h_eff->Draw("L P");

    // Draw reference lines at common cut values
    TLine *line1 = new TLine(0.5, 0, 0.5, 1.0);
    line1->SetLineStyle(2);
    line1->SetLineColor(kGray+2);
    line1->Draw("SAME");

    TLine *line2 = new TLine(2.0, 0, 2.0, 1.0);
    line2->SetLineStyle(2);
    line2->SetLineColor(kGray+2);
    line2->Draw("SAME");

    // Labels
    myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
    myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.20, 0.80, 1, strleg3.c_str(), 0.04);
    myText(0.20, 0.75, 1, "Data", 0.04);
    myText(0.20, 0.70, 1, "10 < #it{E}_{T} < 30 GeV", 0.04);

    myMarkerLineText(0.55, 0.30, 1.5, kBlack, 20, kBlack, 2,
                     "Eff = N(#sigma_{t} < cut) / N_{total}", 0.035, true);

    // Save
    c->SaveAs("mbd_sigma_selection_efficiency.pdf");

    std::cout << "Plot saved to mbd_sigma_selection_efficiency.pdf" << std::endl;

    // Cleanup
    delete h_eff;
    delete h_mbd_combined;
    f_data->Close();
}
