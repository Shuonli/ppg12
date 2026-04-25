#include "plotcommon.h"

void plot_reco_jet_pT(const std::string &infile =
    "../efficiencytool/results/MC_efficiency_bdt_nom.root")
{
    init_plot();

    TFile *fin = TFile::Open(infile.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "Cannot open " << infile << std::endl;
        return;
    }

    TH1F *h_reco_jet_pT     = dynamic_cast<TH1F *>(fin->Get("h_reco_jet_pT"));
    TH1F *h_max_reco_jet_pT = dynamic_cast<TH1F *>(fin->Get("h_max_reco_jet_pT"));
    if (!h_reco_jet_pT || !h_max_reco_jet_pT)
    {
        std::cerr << "Missing reco-jet-pT histograms in " << infile << std::endl;
        return;
    }

    // Rebin from 1000×0.1 GeV bins to 100×1 GeV for display.
    const int rebin = 10;
    h_reco_jet_pT    ->Rebin(rebin);
    h_max_reco_jet_pT->Rebin(rebin);

    h_reco_jet_pT    ->SetLineColor(kAzure + 2);
    h_reco_jet_pT    ->SetLineWidth(3);
    h_reco_jet_pT    ->SetMarkerStyle(20);
    h_reco_jet_pT    ->SetMarkerColor(kAzure + 2);
    h_reco_jet_pT    ->SetMarkerSize(0.7);

    h_max_reco_jet_pT->SetLineColor(kRed + 1);
    h_max_reco_jet_pT->SetLineWidth(3);
    h_max_reco_jet_pT->SetLineStyle(2);
    h_max_reco_jet_pT->SetMarkerStyle(24);
    h_max_reco_jet_pT->SetMarkerColor(kRed + 1);
    h_max_reco_jet_pT->SetMarkerSize(0.7);

    TH1F *frame = new TH1F("frame_reco_jet_pT", "", 100, 0, 60);
    frame->SetXTitle("#it{p}_{T}^{jet,rec} [GeV]");
    frame->SetYTitle("Weighted counts / 1 GeV");
    frame->GetXaxis()->SetRangeUser(0, 60);
    double ymin = 1e2;
    double ymax = std::max(h_reco_jet_pT->GetMaximum(),
                           h_max_reco_jet_pT->GetMaximum()) * 50;
    frame->GetYaxis()->SetRangeUser(ymin, ymax);

    TCanvas *c = new TCanvas("c_reco_jet_pT", "reco jet pT", 800, 700);
    c->SetLogy();
    c->SetTickx();
    c->SetTicky();
    frame->Draw("AXIS");
    h_reco_jet_pT    ->Draw("HIST E SAME");
    h_max_reco_jet_pT->Draw("HIST E SAME");
    frame->Draw("AXIS SAME");

    TLegend *leg = new TLegend(0.58, 0.20, 0.92, 0.35);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(h_reco_jet_pT,     "All reco jets",            "l");
    leg->AddEntry(h_max_reco_jet_pT, "Leading reco jet / event", "l");
    leg->Draw();

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.040);
    tex.DrawLatex(0.20, 0.86, strleg1.c_str());
    tex.SetTextSize(0.034);
    tex.DrawLatex(0.20, 0.81,
                  "#it{p}+#it{p} #kern[-0.05]{#sqrt{#it{s}} = 200 GeV, 48.9 pb^{-1}}");
    tex.DrawLatex(0.20, 0.76, (strSigMC + " (photon5+10+20, all runs)").c_str());
    tex.DrawLatex(0.20, 0.71, "anti-#it{k}_{T} #it{R}=0.4, unsubtracted");

    gSystem->Exec("mkdir -p figures");
    c->SaveAs("figures/reco_jet_pT_signal_MC.pdf");
    std::cout << "Saved figures/reco_jet_pT_signal_MC.pdf" << std::endl;
}
