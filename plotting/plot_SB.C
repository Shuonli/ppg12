#include "plotcommon.h"

void plot_SB()
{
    init_plot();

    // Signal and background both come from the inclusive jet MC. It already
    // contains prompt-photon production at the PYTHIA rate, so adding the
    // photon MC on top would count the prompt photons twice.
    TFile *fin = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_jet_inclusive_combined_showershape.root", "READ");

    // h_ET_isoET: all truth-matched clusters. h_ET_isoET_signal: the subset
    // matched to a truth signal photon (direct or fragmentation, truth iso
    // below the fiducial cut). Background is everything else.
    TH2D* h_all = (TH2D *)fin->Get("h_ET_isoET_eta0");
    TH2D* h_sig = (TH2D *)fin->Get("h_ET_isoET_signal_eta0");
    if (!h_all || !h_sig)
    {
        std::cerr << "plot_SB: h_ET_isoET_eta0 or h_ET_isoET_signal_eta0 missing in " << fin->GetName() << std::endl;
        return;
    }

    int rebinx = 16;
    h_all->RebinX(rebinx);
    h_sig->RebinX(rebinx);

    // No isolation or identification selection: project over the full iso-ET axis.
    TH1D *h_all_proj = h_all->ProjectionX("h_all_proj");
    TH1D *h_sig_proj = h_sig->ProjectionX("h_sig_proj");

    // background = all - signal. The signal clusters are a subset of all
    // with the same weights, so the background sum of squared weights is
    // the difference as well.
    TH1D *h_bg_proj = (TH1D *)h_all_proj->Clone("h_bg_proj");
    for (int ix = 1; ix <= h_bg_proj->GetNbinsX(); ++ix) {
        double err2 = std::pow(h_all_proj->GetBinError(ix), 2) - std::pow(h_sig_proj->GetBinError(ix), 2);
        h_bg_proj->SetBinContent(ix, h_all_proj->GetBinContent(ix) - h_sig_proj->GetBinContent(ix));
        h_bg_proj->SetBinError(ix, std::sqrt(std::max(err2, 0.0)));
    }

    //calculate s/b
    TH1D* h_sb = (TH1D *)h_sig_proj->Clone("h_sb");

    h_sb->Divide(h_bg_proj);

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    frame_et_rec->SetYTitle("S/B");
    double y_max = 1.2 * h_sb->GetMaximum();
    frame_et_rec->GetYaxis()->SetRangeUser(0.0, y_max < 1.0 ? 1.0 : y_max);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 32);
    frame_et_rec->Draw("axis");

    h_sb->SetMarkerColor(kBlack);
    h_sb->SetMarkerStyle(20);
    h_sb->SetMarkerSize(1.5);
    h_sb->SetLineColor(kBlack);
    h_sb->Draw("P same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5, 0.80, 1, strMC.c_str(), 0.04);

    for (int ix = h_sb->FindBin(10.01); ix <= h_sb->FindBin(31.99); ++ix)
        std::cout << "S/B " << h_sb->GetXaxis()->GetBinLowEdge(ix) << "-" << h_sb->GetXaxis()->GetBinUpEdge(ix)
                  << " GeV: " << h_sb->GetBinContent(ix) << " +- " << h_sb->GetBinError(ix) << std::endl;

    c1->SaveAs("figures/SB.pdf");
}
