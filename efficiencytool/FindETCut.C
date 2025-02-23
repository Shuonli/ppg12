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

void FindETCut()
{

    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TFile *fin = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_nom.root", "READ");

    bool x_is_truth = false;

    TH2D *h_eff = (TH2D *)fin->Get("h_singal_truth_isoET_0");

    if(!x_is_truth){
        h_eff = (TH2D *)fin->Get("h_singal_reco_isoET_0");
    }

    
    const int nbins = 13; // 8-40 GeV in 1 GeV bins
    TH1F *h_eff90 = new TH1F("h_eff90", "90% Efficiency Cutoff; pT [GeV]; isoET [GeV]", nbins, 10, 36);
    TH1F *h_eff80 = new TH1F("h_eff80", "80% Efficiency Cutoff; pT [GeV]; isoET [GeV]", nbins, 10, 36);
    TH1F *h_eff70 = new TH1F("h_eff70", "70% Efficiency Cutoff; pT [GeV]; isoET [GeV]", nbins, 10, 36);
    

   /*
    const int nbins = 5; // 8-40 GeV in 1 GeV bins
    TH1F *h_eff90 = new TH1F("h_eff90", "70% Efficiency Cutoff; pT [GeV]; isoET [GeV]", nbins, 10, 15);
    TH1F *h_eff80 = new TH1F("h_eff80", "50% Efficiency Cutoff; pT [GeV]; isoET [GeV]", nbins, 10, 15);
    TH1F *h_eff70 = new TH1F("h_eff70", "30% Efficiency Cutoff; pT [GeV]; isoET [GeV]", nbins, 10, 15);
    */

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
        
        /*
        h_eff90->SetBinError(i, ybinwidth / 2.0);
        h_eff90->SetBinContent(i, findEffCutoff(h_eff, 0.6, pTlow, pThigh, -1.0));
        h_eff80->SetBinError(i, ybinwidth / 2.0);
        h_eff80->SetBinContent(i, findEffCutoff(h_eff, 0.5, pTlow, pThigh, -1.0));
        h_eff70->SetBinError(i, ybinwidth / 2.0);
        h_eff70->SetBinContent(i, findEffCutoff(h_eff, 0.4, pTlow, pThigh, -1.0));
        */
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
        "ETCut_FitResults.pdf"                             // Output
    );

    //print the fit result
    std::cout << "90% Efficiency Fit: EisoT = " << fit90->GetParameter(0) << " + " << fit90->GetParameter(1) << "pT" << std::endl;
    std::cout << "80% Efficiency Fit: EisoT = " << fit80->GetParameter(0) << " + " << fit80->GetParameter(1) << "pT" << std::endl;
    std::cout << "70% Efficiency Fit: EisoT = " << fit70->GetParameter(0) << " + " << fit70->GetParameter(1) << "pT" << std::endl;

    // Cleanup
    // fin->Close();
    // delete fin;
}
