#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <iostream>

#include "plotcommon.h"

// Default suffix = "showershape" (all-range bare config). Pass e.g. "showershape_0rad"
// to restrict to the 0 mrad period.
void plot_mbd_sigma_efficiency(const std::string &config_suffix = "showershape")
{
    init_plot();

    const std::string base = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/";
    const std::string dataFile = base + "data_histoshower_shape_" + config_suffix + ".root";
    const std::string mcFile   = base + "MC_efficiencyshower_shape_jet12_inclusive_" + config_suffix + ".root";

    // Open data file
    TFile *f_data = TFile::Open(dataFile.c_str(), "READ");
    if (!f_data || f_data->IsZombie())
    {
        std::cerr << "Error: Could not open data file: " << dataFile << std::endl;
        return;
    }

    // Open inclusive MC file (jet12 sample as single-sample proxy for MC sigma_t shape)
    TFile *f_mc = TFile::Open(mcFile.c_str(), "READ");
    if (!f_mc || f_mc->IsZombie())
    {
        std::cerr << "Error: Could not open MC file: " << mcFile << std::endl;
        f_data->Close();
        return;
    }

    // pT bin edges (combine all bins)
    std::vector<double> pT_bin_edges = {10, 15, 20, 30};
    int nPtBins = pT_bin_edges.size() - 1;
    int ieta = 0;

    // Helper lambda: combine MBD sigma projections across pT bins from a given file
    auto buildCombined = [&](TFile *f, const char *tag) -> TH1D* {
        TH1D *h_combined = nullptr;
        for (int ipt = 0; ipt < nPtBins; ++ipt)
        {
            TString hist_name = Form("h_mbd_avgsigma_vs_npb_score_eta%d_pt%d", ieta, ipt);
            TH2D *h2 = dynamic_cast<TH2D*>(f->Get(hist_name));
            if (!h2)
            {
                std::cerr << "Warning: Could not retrieve " << hist_name << " from " << tag << std::endl;
                continue;
            }
            TH1D *h_mbd = h2->ProjectionY(Form("h_mbd_sigma_%s_pt%d", tag, ipt));
            if (!h_mbd || h_mbd->Integral() < 1) continue;
            if (!h_combined)
            {
                h_combined = (TH1D*)h_mbd->Clone(Form("h_mbd_combined_%s", tag));
                h_combined->SetDirectory(nullptr);
            }
            else
            {
                h_combined->Add(h_mbd);
            }
        }
        return h_combined;
    };

    TH1D *h_mbd_data = buildCombined(f_data, "data");
    TH1D *h_mbd_mc   = buildCombined(f_mc,   "mc");

    if (!h_mbd_data || h_mbd_data->Integral() < 1)
    {
        std::cerr << "Error: No data to plot!" << std::endl;
        f_data->Close(); f_mc->Close();
        return;
    }
    if (!h_mbd_mc || h_mbd_mc->Integral() < 1)
    {
        std::cerr << "Error: No MC to plot!" << std::endl;
        f_data->Close(); f_mc->Close();
        return;
    }

    // Build cumulative fraction histogram from a raw distribution
    auto buildFraction = [](TH1D *h_raw, const char *name) -> TH1D* {
        int nbins = h_raw->GetNbinsX();
        double xmin = h_raw->GetXaxis()->GetXmin();
        double xmax = h_raw->GetXaxis()->GetXmax();
        TH1D *h_frac = new TH1D(name, "", nbins, xmin, xmax);
        h_frac->SetDirectory(nullptr);
        double total = h_raw->Integral(1, nbins);
        for (int ibin = 1; ibin <= nbins; ++ibin)
        {
            double passing = h_raw->Integral(1, ibin);
            double frac = (total > 0) ? passing / total : 0.0;
            double err  = (total > 0) ? sqrt(frac * (1 - frac) / total) : 0.0;
            h_frac->SetBinContent(ibin, frac);
            h_frac->SetBinError(ibin, err);
        }
        return h_frac;
    };

    TH1D *h_frac_data = buildFraction(h_mbd_data, "h_frac_data");
    TH1D *h_frac_mc   = buildFraction(h_mbd_mc,   "h_frac_mc");

    // Create canvas
    TCanvas *c = new TCanvas("c_mbd_sigma_eff", "MBD Sigma Fraction of Events", 600, 600);

    // Style — data
    h_frac_data->SetLineColor(kBlack);
    h_frac_data->SetMarkerColor(kBlack);
    h_frac_data->SetMarkerStyle(20);
    h_frac_data->SetMarkerSize(1.0);
    h_frac_data->SetLineWidth(2);
    h_frac_data->SetStats(0);
    h_frac_data->SetTitle("");

    h_frac_data->GetXaxis()->SetTitle("MBD Avg #sigma_{t} Cut [ns]");
    h_frac_data->GetXaxis()->SetNdivisions(505);
    h_frac_data->GetXaxis()->SetRangeUser(0, 5);
    h_frac_data->GetYaxis()->SetTitle("Fraction of Events");
    h_frac_data->GetYaxis()->SetTitleOffset(1.5);
    h_frac_data->GetYaxis()->SetRangeUser(0.5, 1.3);

    // Style — MC
    h_frac_mc->SetLineColor(kRed+1);
    h_frac_mc->SetMarkerColor(kRed+1);
    h_frac_mc->SetMarkerStyle(24);
    h_frac_mc->SetMarkerSize(1.0);
    h_frac_mc->SetLineWidth(2);
    h_frac_mc->SetStats(0);

    h_frac_data->Draw("L P");
    h_frac_mc->Draw("L P SAME");

    // Draw reference lines at common cut values
    TLine *line1 = new TLine(0.5, 0.5, 0.5, 1.0);
    line1->SetLineStyle(2);
    line1->SetLineColor(kGray+2);
    line1->Draw("SAME");

    TLine *line2 = new TLine(2.0, 0.5, 2.0, 1.0);
    line2->SetLineStyle(2);
    line2->SetLineColor(kGray+2);
    line2->Draw("SAME");

    // Labels — top left
    myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
    myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.20, 0.80, 1, strleg3.c_str(), 0.04);
    myText(0.20, 0.75, 1, "10 < #it{E}_{T} < 30 GeV", 0.04);

    // Legend — top right
    myMarkerLineText(0.60, 0.90, 1.5, kBlack, 20, kBlack, 2, "Data", 0.035, true);
    myMarkerLineText(0.60, 0.83, 1.5, kRed+1, 24, kRed+1, 2, "Inclusive MC", 0.035, true);

    // Save
    const std::string outPdf = "mbd_sigma_selection_efficiency_" + config_suffix + ".pdf";
    c->SaveAs(outPdf.c_str());

    std::cout << "Plot saved to " << outPdf << std::endl;

    // Cleanup
    delete h_frac_data;
    delete h_frac_mc;
    delete h_mbd_data;
    delete h_mbd_mc;
    f_data->Close();
    f_mc->Close();
}
