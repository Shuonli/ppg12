#include "plotcommon.h"
#include <vector>
#include <algorithm>

const int col_ratio[] = {kAzure + 2, kPink + 5, kSpring - 6, kRed - 4, kBlue - 4};

void plot_time_energy_corr(const char* infile = "/sphenix/user/shuhangli/ppg12/efficiencytool/time_energy_corr.root")
{
    init_plot();

    std::string base_dir = "/sphenix/user/shuhangli/ppg12/efficiencytool/";

    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie())
    {
        std::cerr << "Failed to open file: " << infile << std::endl;
        return;
    }

    // Prefer TH3D h_time_ratio_truth: axes (time, ratio, truth pT)
    TH3D *h_time_ratio_truth = dynamic_cast<TH3D*>(f->Get("h_time_ratio_truth"));
    TH3D *h_time_reco_truth = h_time_ratio_truth ? nullptr : dynamic_cast<TH3D*>(f->Get("h_time_reco_truth"));
    if (!h_time_ratio_truth && !h_time_reco_truth)
    {
        std::cerr << "Neither h_time_ratio_truth (TH3) nor h_time_reco_truth (TH3) found in file: " << infile << std::endl;
        return;
    }

    // Axes pointers are filled depending on which histogram we use
    TAxis *xTime = nullptr;
    TAxis *yReco = nullptr;   // only for fallback TH3 path
    TAxis *zTruth = nullptr;  // for both TH3 paths
    TAxis *yRatio = nullptr;  // for TH3 ratio path
    if (h_time_ratio_truth)
    {
        xTime = h_time_ratio_truth->GetXaxis();
        yRatio = h_time_ratio_truth->GetYaxis();
        zTruth = h_time_ratio_truth->GetZaxis();
    }


    const int nTimeBinsToPlot = 5;
    const double tmin_all = xTime->GetXmin();
    const double tmax_all = xTime->GetXmax();
    double time_edges[nTimeBinsToPlot + 1] = {0, 0.5, 1, 1.5, 2, 2.5};

    // Build validated time edges vector (monotonic). If invalid, fallback to equal-width bins
    std::vector<double> edges(nTimeBinsToPlot + 1);
    for (int i = 0; i <= nTimeBinsToPlot; ++i) edges[i] = time_edges[i];


    std::vector<TH1D*> h_ratio; h_ratio.reserve(nTimeBinsToPlot);
    std::vector<TString> time_bin_labels; time_bin_labels.reserve(nTimeBinsToPlot);
    double xlow = 0.6, xhigh = 1.2;
    // Create ratio histograms for each time slice
    for (int i = 0; i < nTimeBinsToPlot; ++i)
    {
        const double tmin = edges[i];
        const double tmax = edges[i + 1];
        // map edges to x-axis bins
        const int binStart = xTime->FindBin(tmin + 1e-6);
        const int binEnd = xTime->FindBin(tmax - 1e-6);
        time_bin_labels.push_back(Form("%.1f < t < %.1f", tmin, tmax));

        TH1D *h = nullptr;
        if (h_time_ratio_truth)
        {
            // Project ratio distribution Y over full truth pT (Z) in this time range
            h = h_time_ratio_truth->ProjectionY(Form("h_ratio_timebin_%d", i), binStart, binEnd, 1, zTruth->GetNbins());
        }
        
        h->SetDirectory(nullptr);
        h->SetLineColor(col_ratio[i]);
        h->SetMarkerColor(col_ratio[i]);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(1.2);

        // Normalize to unit area if non-empty
        const double integral = h->Integral(h->GetXaxis()->FindBin(xlow), h->GetXaxis()->FindBin(xhigh));
        if (integral > 0) h->Scale(1.0 / integral);

        h_ratio.push_back(h);
    }

    // Determine y-axis max
    double ymax = 0.0;
    for (auto *h : h_ratio)
        ymax = std::max(ymax, h->GetMaximum());

    // Draw
    TCanvas *c = new TCanvas("c_time_ratio", "Reco/Truth ET ratio by time bins", 700, 700);


    TH1F *frame = new TH1F("frame_time_ratio", "", 20, xlow, xhigh);
    frame->GetXaxis()->SetTitle("Reco E_{T} / Truth p_{T}");
    frame->GetYaxis()->SetTitle("Normalized counts");
    frame->GetYaxis()->SetRangeUser(0.0, ymax * 1.7);
    frame->Draw("axis");

    for (size_t i = 0; i < h_ratio.size(); ++i)
    {
        h_ratio[i]->Draw(i == 0 ? "E0 X0 same" : "E0 X0 same");
    }

    float xpos(0.2), xpos2(0.915), ypos(0.885), ypos2(0.65), dy(0.054), dy1(0.055), fontsize(0.046), fontsize1(0.048);
    myText(xpos2, ypos - 0 * dy, 1, strMC.c_str(), fontsize, 1);

    TLegend *leg = new TLegend(0.45, ypos2, xpos2, ypos2 + nTimeBinsToPlot * dy1);
    legStyle(leg, 0.20, 0.06);
    for (int i = 0; i < nTimeBinsToPlot; ++i)
    {
        leg->AddEntry(h_ratio[i], time_bin_labels[i], "l");
    }
    leg->Draw("same");

    c->Update();

    c->SaveAs("figures/time_energy_corr.pdf");

    TCanvas *c2 = new TCanvas("c_time", "1D time distribution", 700, 600);
    TH1D *h_time = h_time_ratio_truth->ProjectionX("h_time", 1, xTime->GetNbins());
    h_time->Draw();
    h_time->GetXaxis()->SetTitle("Time[nsamples]");
    h_time->GetXaxis()->SetRangeUser(-0.5, 3.5);
    h_time->GetYaxis()->SetRangeUser(0, h_time->GetMaximum() * 1.25);
    myText(0.8, 0.8, 1, strMC.c_str(), fontsize, 1);
    c2->Update();
    c2->SaveAs("figures/time_distribution.pdf");
}
