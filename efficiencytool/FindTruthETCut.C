#include <iostream>
#include <TFile.h>
#include <TH2.h>
#include <TMath.h>
#include "draw.C"

float findEffCutoff(TH2D *h_input, float eff, float pTlow, float pThigh, float isolow)
{
    // Determine x bins for pT range
    int xbinlow = h_input->GetXaxis()->FindBin(pTlow);
    int xbinhigh = h_input->GetXaxis()->FindBin(pThigh);

    // Clamp x bins to valid range
    xbinlow = TMath::Max(1, xbinlow);
    xbinhigh = TMath::Min(h_input->GetNbinsX(), xbinhigh);

    // Calculate total sum in the pT range and all isoET bins
    float sum = 0.0;
    for (int i = xbinlow; i <= xbinhigh; ++i)
    {
        for (int j = 1; j <= h_input->GetNbinsY(); ++j)
        {
            sum += h_input->GetBinContent(i, j);
        }
    }

    if (sum <= 0.0)
    {
        return 0.0; // Handle cases with no entries
    }

    float targetsum = sum * eff;

    // Determine starting y bin
    int ybinlow = h_input->GetYaxis()->FindBin(isolow);
    // Clamp ybinlow to valid range
    ybinlow = TMath::Max(1, ybinlow);
    ybinlow = TMath::Min(h_input->GetNbinsY(), ybinlow);

    int ybinhigh = ybinlow;
    float currentsum = 0.0;
    bool found = false;

    // Iterate over y bins from low to high isoET
    for (int j = ybinlow; j <= h_input->GetNbinsY(); ++j)
    {
        // Sum over all pT bins in the current y bin
        for (int i = xbinlow; i <= xbinhigh; ++i)
        {
            currentsum += h_input->GetBinContent(i, j);
        }
        if (currentsum >= targetsum)
        {
            ybinhigh = j;
            found = true;
            break;
        }
    }
    

    if (!found)
    {
        ybinhigh = h_input->GetNbinsY();
    }

    // Return the upper edge of the y bin where the target was met
    float isohigh_eff = h_input->GetYaxis()->GetBinUpEdge(ybinhigh);
    return isohigh_eff;
}

float findEfficiency(TH2D* h_input, float pTlow, float pThigh, float isolow, float isohigh)
{
    // Determine x bins for pT range
    int xbinlow  = h_input->GetXaxis()->FindBin(pTlow);
    int xbinhigh = h_input->GetXaxis()->FindBin(pThigh);

    // Clamp pT bins to valid range
    xbinlow  = TMath::Max(1, xbinlow);
    xbinhigh = TMath::Min(h_input->GetNbinsX(), xbinhigh);

    // Calculate total sum in the pT range (all iso bins)
    float totalSum = 0.0;
    for (int i = xbinlow; i <= xbinhigh; ++i)
    {
        for (int j = 1; j <= h_input->GetNbinsY(); ++j)
        {
            totalSum += h_input->GetBinContent(i, j);
        }
    }

    // If totalSum is zero, avoid division by zero
    if (totalSum <= 0.0)
    {
        return 0.0;
    }

    // Determine y bins for isolation range
    int ybinlow  = h_input->GetYaxis()->FindBin(isolow);
    int ybinhigh = h_input->GetYaxis()->FindBin(isohigh);

    // Clamp iso bins to valid range
    ybinlow  = TMath::Max(1, ybinlow);
    ybinhigh = TMath::Min(h_input->GetNbinsY(), ybinhigh);

    // Calculate sum in the specified iso range and pT range
    float sumInRange = 0.0;
    for (int i = xbinlow; i <= xbinhigh; ++i)
    {
        for (int j = ybinlow; j <= ybinhigh; ++j)
        {
            sumInRange += h_input->GetBinContent(i, j);
        }
    }

    // Calculate and return the efficiency
    float efficiency = sumInRange / totalSum;
    return efficiency;
}

/*
void FindTruthETCut()
{

    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TFile *fin = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_merge.root", "READ");

    bool x_is_truth = false;

    TH2D *h_eff = (TH2D *)fin->Get("h_singal_truth_isoET_0");

    if(!x_is_truth){
        h_eff = (TH2D *)fin->Get("h_singal_reco_isoET_0");
    }

    
    const int nbins = 13; // 8-40 GeV in 1 GeV bins
    TH1F *h_eff90 = new TH1F("h_eff90", "90% Efficiency Cutoff; pT [GeV]; isoET [GeV]", nbins, 10, 36);
    TH1F *h_eff80 = new TH1F("h_eff80", "80% Efficiency Cutoff; pT [GeV]; isoET [GeV]", nbins, 10, 36);
    TH1F *h_eff70 = new TH1F("h_eff70", "70% Efficiency Cutoff; pT [GeV]; isoET [GeV]", nbins, 10, 36);
    


    // Populate histograms
    float ybinwidth = h_eff->GetYaxis()->GetBinWidth(1);
    for (int i = 1; i <= nbins; ++i)
    {
        float pTlow = h_eff90->GetXaxis()->GetBinLowEdge(i);
        float pThigh = h_eff90->GetXaxis()->GetBinUpEdge(i);


        
        h_eff90->SetBinError(i, ybinwidth / 2.0);
        h_eff90->SetBinContent(i, findEffCutoff(h_eff, 0.9, pTlow, pThigh, -1.0));
        h_eff80->SetBinError(i, ybinwidth / 2.0);
        h_eff80->SetBinContent(i, findEffCutoff(h_eff, 0.8, pTlow, pThigh, -1.0));
        h_eff70->SetBinError(i, ybinwidth / 2.0);
        h_eff70->SetBinContent(i, findEffCutoff(h_eff, 0.7, pTlow, pThigh, -1.0));
        
    }

    std::vector<int> colors = {kPink + 8, kSpring - 7, kAzure - 3, kViolet + 3, kOrange + 10};
    // Perform linear fits
    TF1 *fit90 = new TF1("fit90", "pol1", 8, 40);
    fit90->SetLineColor(colors[0]);
    h_eff90->Fit(fit90, "RQ");
    TF1 *fit80 = new TF1("fit80", "pol1", 8, 40);
    fit80->SetLineColor(colors[1]);
    h_eff80->Fit(fit80, "RQ");
    TF1 *fit70 = new TF1("fit70", "pol1", 8, 40);
    fit70->SetLineColor(colors[2]);
    h_eff70->Fit(fit70, "RQ");
    

    // Prepare fit result text
    std::vector<std::string> fitTexts = {
        "Pythia, #sqrt{s} = 200 GeV",
        Form("90%%: E^{iso}_{T} = %.3f  + %.3fp_{T}",
             fit90->GetParameter(0),
             fit90->GetParameter(1)),
        Form("80%%: E^{iso}_{T} = %.3f  + %.3fp_{T}",
             fit80->GetParameter(0),
             fit80->GetParameter(1)),
        Form("70%%: E^{iso}_{T} = %.3f  + %.3fp_{T}",
             fit70->GetParameter(0),
             fit70->GetParameter(1)),
        "vtx |z| < 30 cm, |#eta^{#gamma}|<0.7"};

    // Plotting parameters

    std::vector<int> markers = {20, 21, 22};
    std::vector<std::string> legends = {"90% Efficiency", "80% Efficiency", "70% Efficiency"};

    std::string xtitle = "Cluster p_{T} [GeV]";

    if (x_is_truth)
    {
        xtitle = "Truth #gamma p_{T} [GeV]";
    }

    // Draw plots with fit results
    draw_1D_multiple_plot(
        {h_eff90, h_eff80, h_eff70}, colors, markers,
        false, 1, false,                                   // Rebin/Normalize
        true, 10, 3, false,                                // X-axis
        true, 0.5, 3.5, false,                                 // Y-axis
        true, xtitle, "E^{iso}_{T} Cutoff [GeV]", // Titles
        false, "",                                         // Data/Run status
        true, fitTexts, 0.15, 0.85, 0.035,                 // Text annotations
        true, legends, 0.6, 0.80, 0.035,                   // Legend
        "plots/ETCut_FitResults.pdf"                             // Output
    );

    //print the fit result
    std::cout << "90% Efficiency Fit: EisoT = " << fit90->GetParameter(0) << " + " << fit90->GetParameter(1) << "pT" << std::endl;
    std::cout << "80% Efficiency Fit: EisoT = " << fit80->GetParameter(0) << " + " << fit80->GetParameter(1) << "pT" << std::endl;
    std::cout << "70% Efficiency Fit: EisoT = " << fit70->GetParameter(0) << " + " << fit70->GetParameter(1) << "pT" << std::endl;

    // Cleanup
    // fin->Close();
    // delete fin;
}
*/
void FindTruthETCut()
{
    // Style options.
    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // Open the input file.
    TFile *fin = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_noiso.root", "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }

    // Retrieve the two 2D histograms.
    // For this analysis:
    // - h_direct is used for direct photons.
    // - h_frag is used for fragmentation photons.
    TH2D *h_direct = (TH2D*)fin->Get("h_direct_pT_truth_isoET_0");
    TH2D *h_frag   = (TH2D*)fin->Get("h_frag_pT_truth_isoET_0");
    if (!h_direct || !h_frag) {
        std::cerr << "Error: one or both histograms not found!" << std::endl;
        return;
    }

    // Define three pₜ ranges (in GeV) for the efficiency study.
    std::vector< std::pair<float, float> > pT_ranges;
    pT_ranges.push_back(std::make_pair(10.0, 15.0));
    pT_ranges.push_back(std::make_pair(15.0, 20.0));
    pT_ranges.push_back(std::make_pair(25.0, 30.0));

    //--- Efficiency versus isoET cutoff ---
    // We scan the isoET cutoff from 1 to 20 GeV in bins of width 2 GeV.
    // (A histogram defined from 1 to 21 GeV yields 10 bins.)
    int nbinsCut = 20;
    float cutoff_low = 0.0;
    float cutoff_high = 20.0;

    // Create efficiency histograms (one per pₜ range) for direct photons.
    std::vector<TH1F*> h_eff_direct;
    std::vector<std::string> legend_direct;
    for (size_t i = 0; i < pT_ranges.size(); ++i) {
         TString name = Form("h_eff_direct_%zu", i);
         TString title = Form("Direct Efficiency (p_{T}: %.0f-%.0f GeV); IsoET Cutoff [GeV]; Efficiency",
                              pT_ranges[i].first, pT_ranges[i].second);
         TH1F* hist = new TH1F(name, title, nbinsCut, cutoff_low, cutoff_high);
         h_eff_direct.push_back(hist);
         legend_direct.push_back(Form("Direct #gamma p_{T} [%.0f, %.0f]", pT_ranges[i].first, pT_ranges[i].second));
    }
    
    // Create efficiency histograms for fragmentation photons.
    std::vector<TH1F*> h_eff_frag;
    std::vector<std::string> legend_frag;
    for (size_t i = 0; i < pT_ranges.size(); ++i) {
         TString name = Form("h_eff_frag_%zu", i);
         TString title = Form("Frag Efficiency (p_{T}: %.0f-%.0f GeV); IsoET Cutoff [GeV]; Efficiency",
                              pT_ranges[i].first, pT_ranges[i].second);
         TH1F* hist = new TH1F(name, title, nbinsCut, cutoff_low, cutoff_high);
         h_eff_frag.push_back(hist);
         legend_frag.push_back(Form("Fragmentation #gamma p_{T} [%.0f, %.0f]", pT_ranges[i].first, pT_ranges[i].second));
    }

    // Loop over the cutoff bins.
    for (int ibin = 1; ibin <= nbinsCut; ++ibin) {
         float cutoff_val = h_eff_direct[0]->GetXaxis()->GetBinUpEdge(ibin);
         // For each pₜ range, calculate the efficiency using findEfficiency.
         // We use isolow = -1.0 to include all isoET values below the cutoff.
         for (size_t i = 0; i < pT_ranges.size(); ++i) {
              float pTlow  = pT_ranges[i].first;
              float pThigh = pT_ranges[i].second;
              float eff_direct = findEfficiency(h_direct, pTlow, pThigh, -1.0, cutoff_val);
              h_eff_direct[i]->SetBinContent(ibin, eff_direct);
              h_eff_direct[i]->SetBinError(ibin, 0.0);
              float eff_frag = findEfficiency(h_frag, pTlow, pThigh, -1.0, cutoff_val);
              h_eff_frag[i]->SetBinContent(ibin, eff_frag);
              h_eff_frag[i]->SetBinError(ibin, 0.0);
         }
    }

    //--- Define color arrays using the provided values.
    // Provided colors:
    // { kPink+8, kSpring-7, kAzure-3, kViolet+3, kOrange+10 }
    // For the efficiency curves we will use the first three colors for direct photons and
    // a different set for fragmentation photons.
    std::vector<int> colors = {kPink + 8, kSpring - 7, kAzure - 3, kViolet + 3, kOrange + 10};
    std::vector<int> direct_colors = { kPink+8, kSpring-7, kAzure-3 };
    std::vector<int> frag_colors   = { kAzure-3, kViolet+3, kOrange+10 };
    std::vector<int> markers = {20, 21, 22};

    // Some text annotations.
    std::vector<std::string> texts;
    texts.push_back("Pythia, #sqrt{s}=200 GeV");
    texts.push_back("vtx |z| < 30 cm, |#eta^{#gamma}| < 0.7");

    //--- Draw Direct Photon Efficiency vs. IsoET Cutoff ---
    draw_1D_multiple_plot(h_eff_direct, colors, markers,
                          false, 1, false,        // no rebin, no normalization
                          false, 0, 0, false,       // no forced x-range
                          true, 0.9, 1.05, false,       // no forced y-range
                          true, "E_{T}^{iso} Cutoff [GeV]", "Fraction of Events",
                          false, "Photon+Jet Samples", // run_status text
                          true, texts, 0.15, 0.85, 0.035,
                          true, legend_direct, 0.6, 0.80, 0.035,
                          "plots/Direct_Efficiency_vs_IsoET_Cutoff.pdf");

    //--- Draw Fragmentation Photon Efficiency vs. IsoET Cutoff ---
    draw_1D_multiple_plot(h_eff_frag, colors, markers,
                          false, 1, false,
                          false, 0, 0, false,
                          true, 0.8, 1.1, false,
                          true, "E_{T}^{iso} Cutoff [GeV]", "Fraction of Events",
                          false, "Photon+Jet Samples",
                          true, texts, 0.15, 0.85, 0.035,
                          true, legend_frag, 0.6, 0.80, 0.035,
                          "plots/Frag_Efficiency_vs_IsoET_Cutoff.pdf");

    //--- pₜ Spectrum with a 4 GeV isoET Cut ---
    // Project a pₜ spectrum from each 2D histogram by summing over isoET bins with upper edge ≤ 4 GeV.
    int nptbins = h_direct->GetNbinsX();
    float pt_min = h_direct->GetXaxis()->GetBinLowEdge(1);
    float pt_max = h_direct->GetXaxis()->GetBinUpEdge(nptbins);
    
    TH1F *h_spec_direct = new TH1F("h_spec_direct", "Direct p_{T} Spectrum (IsoET < 4 GeV); p_{T} [GeV]; Counts", nptbins, pt_min, pt_max);
    TH1F *h_spec_frag   = new TH1F("h_spec_frag",   "Frag p_{T} Spectrum (IsoET < 4 GeV); p_{T} [GeV]; Counts", nptbins, pt_min, pt_max);
    
    for (int ix = 1; ix <= nptbins; ++ix) {
         float sum_direct = 0.0;
         float sum_frag   = 0.0;
         for (int iy = 1; iy <= h_direct->GetNbinsY(); ++iy) {
              float isoET_up = h_direct->GetYaxis()->GetBinUpEdge(iy);
              if (isoET_up <= 4.0) {
                  sum_direct += h_direct->GetBinContent(ix, iy);
                  sum_frag   += h_frag->GetBinContent(ix, iy);
              }
         }
         h_spec_direct->SetBinContent(ix, sum_direct);
         h_spec_frag->SetBinContent(ix, sum_frag);
    }
    h_spec_direct->Rebin(10);
    h_spec_frag->Rebin(10);

    TH1F* h_spec_total = (TH1F*)h_spec_direct->Clone("h_spec_total");
    h_spec_total->Add(h_spec_frag);
    
    
    // For the pₜ spectrum, we compare direct and fragmentation photons.
    std::vector<TH1F*> h_spec;
    h_spec.push_back(h_spec_total);
    h_spec.push_back(h_spec_direct);
    h_spec.push_back(h_spec_frag);
    // For the pₜ spectrum, we use two colors from the provided array.
    std::vector<int> spec_markers = {20, 21, 22};
    std::vector<std::string> spec_legend;
    spec_legend.push_back("Total");
    spec_legend.push_back("Direct");
    spec_legend.push_back("Fragmentation");

    texts.clear();
    texts.push_back("Pythia, #sqrt{s}=200 GeV");
    texts.push_back("|#eta^{#gamma}| < 0.7");
    //texts.push_back("R = 0.3, E_{T}^{iso} < 4 GeV");
    texts.push_back("no E_{T}^{iso} requirement");
    
    draw_1D_multiple_plot_ratio(h_spec, colors, spec_markers,
                          false, 1, false,
                          true, 10, 35, false,
                          false, 0, 0, true,
                          true, 0.5, 1,
                          true, "p_{T} [GeV]", "Counts", "Direct/Total",
                          false, "Photon Jet Samples",
                          true, texts, 0.25, 0.85, 0.035,
                          true, spec_legend, 0.6, 0.80, 0.035,
                          "plots/pT_Spectrum_IsoET_200GeV_Cut.pdf");

    // Optionally, close the file.
    // fin->Close();
}