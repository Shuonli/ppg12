#include <iostream>
#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH1D.h>

void CheckPhotonPtRange(const std::string &inputfile = "results/truth_photon_tower_analysis.root",
                        const std::string &outputfile = "results/projection_check.root",
                        double eta_lo = -0.35,
                        double eta_hi = 0.35)
{
    std::cout << "==================================================" << std::endl;
    std::cout << "  Check Histogram Ranges" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "Input file: " << inputfile << std::endl;
    std::cout << "Output file: " << outputfile << std::endl;
    std::cout << "Eta range for 3D projection: [" << eta_lo << ", " << eta_hi << "]" << std::endl;
    std::cout << std::endl;

    // Create output file for saving projections
    TFile *fout = new TFile(outputfile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie())
    {
        std::cerr << "ERROR: Cannot create output file: " << outputfile << std::endl;
        return;
    }

    // Open input file
    TFile *fin = TFile::Open(inputfile.c_str());
    if (!fin || fin->IsZombie())
    {
        std::cerr << "ERROR: Cannot open input file: " << inputfile << std::endl;
        return;
    }

    // Get 2D histograms
    TH2D *h_tower_e_vs_time = (TH2D *)fin->Get("h_tower_e_vs_time");
    TH2D *h_leading_tower_e_vs_time = (TH2D *)fin->Get("h_leading_tower_e_vs_time");
    TH2D *h_nonleading_tower_e_vs_time = (TH2D *)fin->Get("h_nonleading_tower_e_vs_time");

    // Get 3D histograms
    TH3D *h3_tower_e_vs_time_vs_eta = (TH3D *)fin->Get("h3_tower_e_vs_time_vs_eta");
    TH3D *h3_leading_tower_e_vs_time_vs_eta = (TH3D *)fin->Get("h3_leading_tower_e_vs_time_vs_eta");
    TH3D *h3_leading_tower_adc_vs_time_vs_eta = (TH3D *)fin->Get("h3_leading_tower_adc_vs_time_vs_eta");

    if (!h_tower_e_vs_time)
    {
        std::cerr << "ERROR: Cannot find h_tower_e_vs_time" << std::endl;
        fin->Close();
        return;
    }

    std::cout << "==================================================" << std::endl;
    std::cout << "  All Towers (h_tower_e_vs_time)" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "Total entries: " << h_tower_e_vs_time->GetEntries() << std::endl;
    std::cout << "Histogram range (Y-axis, Energy): [" << h_tower_e_vs_time->GetYaxis()->GetXmin()
              << ", " << h_tower_e_vs_time->GetYaxis()->GetXmax() << "] GeV" << std::endl;

    // Find actual filled range
    TH1D *proj_e_all = (TH1D *)h_tower_e_vs_time->ProjectionY("proj_e_all");
    if (proj_e_all->GetEntries() > 0)
    {
        int first_bin = proj_e_all->FindFirstBinAbove(0);
        int last_bin = proj_e_all->FindLastBinAbove(0);
        double min_e = proj_e_all->GetBinLowEdge(first_bin);
        double max_e = proj_e_all->GetBinLowEdge(last_bin) + proj_e_all->GetBinWidth(last_bin);
        std::cout << "Actual filled range: [" << min_e << ", " << max_e << "] GeV" << std::endl;
        std::cout << "Mean tower energy: " << proj_e_all->GetMean() << " GeV" << std::endl;
        std::cout << "RMS: " << proj_e_all->GetRMS() << " GeV" << std::endl;

        // Check entries above certain thresholds
        std::cout << "\nEntries above energy thresholds:" << std::endl;
        for (double threshold : {5.0, 10.0, 12.0, 15.0, 18.0, 20.0})
        {
            int bin = proj_e_all->FindBin(threshold);
            double entries = proj_e_all->Integral(bin, proj_e_all->GetNbinsX());
            std::cout << "  E > " << threshold << " GeV: " << entries << " entries ("
                      << (100.0 * entries / proj_e_all->GetEntries()) << "%)" << std::endl;
        }
    }
    fout->cd();
    proj_e_all->Write();
    delete proj_e_all;

    if (h_leading_tower_e_vs_time)
    {
        std::cout << "\n==================================================" << std::endl;
        std::cout << "  Leading Towers (h_leading_tower_e_vs_time)" << std::endl;
        std::cout << "==================================================" << std::endl;
        std::cout << "Total entries: " << h_leading_tower_e_vs_time->GetEntries() << std::endl;
        std::cout << "Histogram range (Y-axis, Energy): [" << h_leading_tower_e_vs_time->GetYaxis()->GetXmin()
                  << ", " << h_leading_tower_e_vs_time->GetYaxis()->GetXmax() << "] GeV" << std::endl;

        TH1D *proj_e_lead = (TH1D *)h_leading_tower_e_vs_time->ProjectionY("proj_e_lead");
        if (proj_e_lead->GetEntries() > 0)
        {
            int first_bin = proj_e_lead->FindFirstBinAbove(0);
            int last_bin = proj_e_lead->FindLastBinAbove(0);
            double min_e = proj_e_lead->GetBinLowEdge(first_bin);
            double max_e = proj_e_lead->GetBinLowEdge(last_bin) + proj_e_lead->GetBinWidth(last_bin);
            std::cout << "Actual filled range: [" << min_e << ", " << max_e << "] GeV" << std::endl;
            std::cout << "Mean leading tower energy: " << proj_e_lead->GetMean() << " GeV" << std::endl;
            std::cout << "RMS: " << proj_e_lead->GetRMS() << " GeV" << std::endl;

            std::cout << "\nEntries above energy thresholds:" << std::endl;
            for (double threshold : {5.0, 10.0, 12.0, 15.0, 18.0, 20.0})
            {
                int bin = proj_e_lead->FindBin(threshold);
                double entries = proj_e_lead->Integral(bin, proj_e_lead->GetNbinsX());
                std::cout << "  E > " << threshold << " GeV: " << entries << " entries ("
                          << (100.0 * entries / proj_e_lead->GetEntries()) << "%)" << std::endl;
            }
        }
        fout->cd();
        proj_e_lead->Write();
        delete proj_e_lead;
    }

    if (h_nonleading_tower_e_vs_time)
    {
        std::cout << "\n==================================================" << std::endl;
        std::cout << "  Non-Leading Towers (h_nonleading_tower_e_vs_time)" << std::endl;
        std::cout << "==================================================" << std::endl;
        std::cout << "Total entries: " << h_nonleading_tower_e_vs_time->GetEntries() << std::endl;

        TH1D *proj_e_nonlead = (TH1D *)h_nonleading_tower_e_vs_time->ProjectionY("proj_e_nonlead");
        if (proj_e_nonlead->GetEntries() > 0)
        {
            int first_bin = proj_e_nonlead->FindFirstBinAbove(0);
            int last_bin = proj_e_nonlead->FindLastBinAbove(0);
            double min_e = proj_e_nonlead->GetBinLowEdge(first_bin);
            double max_e = proj_e_nonlead->GetBinLowEdge(last_bin) + proj_e_nonlead->GetBinWidth(last_bin);
            std::cout << "Actual filled range: [" << min_e << ", " << max_e << "] GeV" << std::endl;
            std::cout << "Mean non-leading tower energy: " << proj_e_nonlead->GetMean() << " GeV" << std::endl;
            std::cout << "RMS: " << proj_e_nonlead->GetRMS() << " GeV" << std::endl;
        }
        fout->cd();
        proj_e_nonlead->Write();
        delete proj_e_nonlead;
    }

    // Check 3D histograms
    if (h3_leading_tower_adc_vs_time_vs_eta)
    {
        std::cout << "\n==================================================" << std::endl;
        std::cout << "  3D Leading Tower ADC (h3_leading_tower_adc_vs_time_vs_eta)" << std::endl;
        std::cout << "==================================================" << std::endl;
        std::cout << "Total entries: " << h3_leading_tower_adc_vs_time_vs_eta->GetEntries() << std::endl;
        std::cout << "Axis ranges:" << std::endl;
        std::cout << "  X (Time): [" << h3_leading_tower_adc_vs_time_vs_eta->GetXaxis()->GetXmin()
                  << ", " << h3_leading_tower_adc_vs_time_vs_eta->GetXaxis()->GetXmax() << "] ns" << std::endl;
        std::cout << "  Y (ADC): [" << h3_leading_tower_adc_vs_time_vs_eta->GetYaxis()->GetXmin()
                  << ", " << h3_leading_tower_adc_vs_time_vs_eta->GetYaxis()->GetXmax() << "]" << std::endl;
        std::cout << "  Z (Eta): [" << h3_leading_tower_adc_vs_time_vs_eta->GetZaxis()->GetXmin()
                  << ", " << h3_leading_tower_adc_vs_time_vs_eta->GetZaxis()->GetXmax() << "]" << std::endl;
        std::cout << "  Z axis bins: " << h3_leading_tower_adc_vs_time_vs_eta->GetZaxis()->GetNbins() << std::endl;

        // Manual projection in eta range
        std::cout << "\nProjecting to XY in eta range [" << eta_lo << ", " << eta_hi << "]..." << std::endl;
        TH3D *h3_clone = (TH3D *)h3_leading_tower_adc_vs_time_vs_eta->Clone("h3_clone");
        h3_clone->SetDirectory(nullptr);

        const double eps = 1e-6;
        int zlo = h3_clone->GetZaxis()->FindBin(eta_lo + eps);
        int zhi = h3_clone->GetZaxis()->FindBin(eta_hi - eps);
        std::cout << "  Z-axis bin range: [" << zlo << ", " << zhi << "]" << std::endl;
        std::cout << "  Corresponding eta values: [" << h3_clone->GetZaxis()->GetBinLowEdge(zlo)
                  << ", " << h3_clone->GetZaxis()->GetBinUpEdge(zhi) << "]" << std::endl;

        h3_clone->GetZaxis()->SetRange(zlo, zhi);
        TH2D *h2_proj = (TH2D *)h3_clone->Project3D("xy");  // Fixed: use "xy" not "yx"
        h2_proj->SetName("h2_leading_adc_projected");
        h2_proj->SetDirectory(nullptr);

        std::cout << "\nProjected 2D histogram (Time vs ADC in eta range):" << std::endl;
        std::cout << "  Entries: " << h2_proj->GetEntries() << std::endl;

        TH1D *proj_adc = (TH1D *)h2_proj->ProjectionY("proj_adc_from3d");  // Now projects ADC (Y-axis)
        if (proj_adc->GetEntries() > 0)
        {
            int first_bin = proj_adc->FindFirstBinAbove(0);
            int last_bin = proj_adc->FindLastBinAbove(0);
            double min_adc = proj_adc->GetBinLowEdge(first_bin);
            double max_adc = proj_adc->GetBinLowEdge(last_bin) + proj_adc->GetBinWidth(last_bin);
            std::cout << "  ADC filled range: [" << min_adc << ", " << max_adc << "]" << std::endl;
            std::cout << "  Mean ADC: " << proj_adc->GetMean() << std::endl;
        }

        fout->cd();
        h2_proj->Write();
        proj_adc->Write();
        delete proj_adc;
        delete h2_proj;
        delete h3_clone;
    }

    if (h3_leading_tower_e_vs_time_vs_eta)
    {
        std::cout << "\n==================================================" << std::endl;
        std::cout << "  3D Leading Tower Energy (h3_leading_tower_e_vs_time_vs_eta)" << std::endl;
        std::cout << "==================================================" << std::endl;
        std::cout << "Total entries: " << h3_leading_tower_e_vs_time_vs_eta->GetEntries() << std::endl;
        std::cout << "Axis ranges:" << std::endl;
        std::cout << "  X (Time): [" << h3_leading_tower_e_vs_time_vs_eta->GetXaxis()->GetXmin()
                  << ", " << h3_leading_tower_e_vs_time_vs_eta->GetXaxis()->GetXmax() << "] ns" << std::endl;
        std::cout << "  Y (Energy): [" << h3_leading_tower_e_vs_time_vs_eta->GetYaxis()->GetXmin()
                  << ", " << h3_leading_tower_e_vs_time_vs_eta->GetYaxis()->GetXmax() << "] GeV" << std::endl;
        std::cout << "  Z (Eta): [" << h3_leading_tower_e_vs_time_vs_eta->GetZaxis()->GetXmin()
                  << ", " << h3_leading_tower_e_vs_time_vs_eta->GetZaxis()->GetXmax() << "]" << std::endl;
        std::cout << "  Z axis bins: " << h3_leading_tower_e_vs_time_vs_eta->GetZaxis()->GetNbins() << std::endl;

        // Manual projection in eta range
        std::cout << "\nProjecting to XY in eta range [" << eta_lo << ", " << eta_hi << "]..." << std::endl;
        TH3D *h3_clone = (TH3D *)h3_leading_tower_e_vs_time_vs_eta->Clone("h3_clone_e");
        h3_clone->SetDirectory(nullptr);

        const double eps = 1e-6;
        int zlo = h3_clone->GetZaxis()->FindBin(eta_lo + eps);
        int zhi = h3_clone->GetZaxis()->FindBin(eta_hi - eps);
        std::cout << "  Z-axis bin range: [" << zlo << ", " << zhi << "]" << std::endl;
        std::cout << "  Corresponding eta values: [" << h3_clone->GetZaxis()->GetBinLowEdge(zlo)
                  << ", " << h3_clone->GetZaxis()->GetBinUpEdge(zhi) << "]" << std::endl;

        h3_clone->GetZaxis()->SetRange(zlo, zhi);
        TH2D *h2_proj = (TH2D *)h3_clone->Project3D("xy");  // Fixed: use "xy" not "yx"
        h2_proj->SetName("h2_leading_e_projected");
        h2_proj->SetDirectory(nullptr);

        std::cout << "\nProjected 2D histogram (Time vs Energy in eta range):" << std::endl;
        std::cout << "  Entries: " << h2_proj->GetEntries() << std::endl;

        TH1D *proj_e = (TH1D *)h2_proj->ProjectionY("proj_e_lead_from3d");  // Now projects Energy (Y-axis)
        if (proj_e->GetEntries() > 0)
        {
            int first_bin = proj_e->FindFirstBinAbove(0);
            int last_bin = proj_e->FindLastBinAbove(0);
            double min_e = proj_e->GetBinLowEdge(first_bin);
            double max_e = proj_e->GetBinLowEdge(last_bin) + proj_e->GetBinWidth(last_bin);
            std::cout << "  Energy filled range: [" << min_e << ", " << max_e << "] GeV" << std::endl;
            std::cout << "  Mean energy: " << proj_e->GetMean() << " GeV" << std::endl;
            std::cout << "  RMS: " << proj_e->GetRMS() << " GeV" << std::endl;

            std::cout << "\nEntries above energy thresholds:" << std::endl;
            for (double threshold : {5.0, 10.0, 12.0, 15.0, 18.0, 20.0})
            {
                int bin = proj_e->FindBin(threshold);
                double entries = proj_e->Integral(bin, proj_e->GetNbinsX());
                std::cout << "  E > " << threshold << " GeV: " << entries << " entries ("
                          << (100.0 * entries / proj_e->GetEntries()) << "%)" << std::endl;
            }
        }

        fout->cd();
        h2_proj->Write();
        proj_e->Write();
        delete proj_e;
        delete h2_proj;
        delete h3_clone;
    }

    std::cout << "\n==================================================" << std::endl;
    std::cout << "Done!" << std::endl;
    std::cout << "==================================================" << std::endl;

    // Close and save output file
    fout->cd();
    fout->Close();
    std::cout << "\nProjection histograms saved to: " << outputfile << std::endl;
    std::cout << "Histograms in output file:" << std::endl;
    std::cout << "  - proj_e_all (from 2D h_tower_e_vs_time)" << std::endl;
    std::cout << "  - proj_e_lead (from 2D h_leading_tower_e_vs_time)" << std::endl;
    std::cout << "  - proj_e_nonlead (from 2D h_nonleading_tower_e_vs_time)" << std::endl;
    std::cout << "  - h2_leading_adc_projected (2D projection from 3D)" << std::endl;
    std::cout << "  - proj_adc_from3d (1D projection from 3D)" << std::endl;
    std::cout << "  - h2_leading_e_projected (2D projection from 3D)" << std::endl;
    std::cout << "  - proj_e_lead_from3d (1D projection from 3D)" << std::endl;
    std::cout << "==================================================" << std::endl;

    fin->Close();
    delete fin;
    delete fout;
}
