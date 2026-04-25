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

void FindETCut(const char* var_type, double target_eff = 0.9)
{

    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    std::string infile = std::string("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_") + var_type + ".root";
    TFile *fin = new TFile(infile.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << infile << std::endl;
        return;
    }

    TH2D *h_eff_input = (TH2D *)fin->Get("h_singal_reco_isoET_0");
    if (!h_eff_input)
    {
        std::cerr << "ERROR: histogram h_singal_reco_isoET_0 missing in " << infile << std::endl;
        return;
    }

    const int nbins = 13; // 8-40 GeV in 1 GeV bins
    std::string htitle = Form("%d%% Efficiency Cutoff; pT [GeV]; isoET [GeV]", (int)(100 * target_eff));
    TH1F *h_eff = new TH1F("h_eff", htitle.c_str(), nbins, 10, 36);

    // Populate histogram
    float ybinwidth = h_eff_input->GetYaxis()->GetBinWidth(1);
    for (int i = 1; i <= nbins; ++i)
    {
        float pTlow = h_eff->GetXaxis()->GetBinLowEdge(i);
        float pThigh = h_eff->GetXaxis()->GetBinUpEdge(i);

        h_eff->SetBinError(i, ybinwidth / 2.0);
        h_eff->SetBinContent(i, findEffCutoff(h_eff_input, target_eff, pTlow, pThigh, -1.0));
    }

    std::vector<int> colors = {kPink + 8, kPink + 8};
    // Perform linear fit over the matched pT range [10, 36]
    TF1 *fit = new TF1("fit", "pol1", 10, 36);
    fit->SetLineColor(colors[0]);
    h_eff->Fit(fit, "RQ");

    double intercept = fit->GetParameter(0);
    double slope = fit->GetParameter(1);

    // Prepare fit result text
    std::vector<std::string> fitTexts = {
        "Pythia, #sqrt{s} = 200 GeV",
        Form("%d%%: E^{iso}_{T} = %.3f  + %.3fp_{T}",
             (int)(100 * target_eff), intercept, slope),
        "vtx |z| < 30 cm, |#eta^{#gamma}|<0.7"};

    std::vector<int> markers = {20, 20};
    std::vector<std::string> legends = {Form("%d%% Efficiency", (int)(100 * target_eff)), ""};

    std::string xtitle = "Cluster p_{T} [GeV]";

    std::string outpdf = std::string("ETCut_FitResults_") + var_type + "_eff" + std::to_string((int)(100 * target_eff)) + ".pdf";

    // draw_1D_multiple_plot requires >=2 histograms; duplicate for visual (fit is already shown)
    TH1F *h_eff_dup = (TH1F*)h_eff->Clone("h_eff_dup");

    // Draw plots with fit results
    draw_1D_multiple_plot(
        {h_eff, h_eff_dup}, colors, markers,
        false, 1, false,                                   // Rebin/Normalize
        true, 10, 3, false,                                // X-axis
        true, -0.5, 5.0, false,                                 // Y-axis
        true, xtitle, "E^{iso}_{T} Cutoff [GeV]", // Titles
        false, "",                                         // Data/Run status
        true, fitTexts, 0.15, 0.85, 0.035,                 // Text annotations
        true, legends, 0.6, 0.80, 0.035,                   // Legend
        outpdf.c_str()                                      // Output
    );

    // Machine-parseable line for driver scripts
    std::cout.precision(17);
    std::cout << "DERIVED_ISO_CUT var_type=" << var_type
              << " target_eff=" << target_eff
              << " intercept=" << intercept
              << " slope=" << slope << std::endl;
}

void FindETCut()
{
    // Preserves original no-arg behaviour which opened MC_efficiency_bdt_nom.root
    FindETCut("nom", 0.9);
}
