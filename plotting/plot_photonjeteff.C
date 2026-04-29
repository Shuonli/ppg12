#include "plotcommon.h"

const int col[] = {kAzure + 2, kPink + 5, kAzure+3, kTeal -6,         kSpring - 6, kRed - 4, kBlack, kBlue - 3, kPink - 5, kGreen + 3, kBlue - 3};
const int mkcol[] = {kAzure+2, kPink + 5, kAzure+3, kTeal-6,           kSpring - 6, kRed - 4, kBlack, kBlue - 3, kPink - 5, kGreen + 3, kBlue - 3};
const int linecol[] = {kAzure+2, kPink + 5, kAzure+3, kTeal-6,           kSpring - 6, kRed - 4, kBlack, kBlue - 3, kPink - 5, kGreen + 3, kBlue - 3};


void plot_photonjeteff()
{
    init_plot();

    TFile *fmc = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root");

    TH2D *h_all_cluster_Et_max_b2bjet = (TH2D *)fmc->Get("h_all_cluster_Et_max_b2bjet_0");
    TH1D *h_all_cluster_Et = (TH1D *)fmc->Get("h_all_cluster_signal_0");

    assert(h_all_cluster_Et);
    assert(h_all_cluster_Et_max_b2bjet);


    std::vector<float> jet_pT_bins = {5, 7, 9};
    std::vector<TH1D*> h_eff = {};
    float max_jet_pT = 100;
    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    frame_et_rec->GetYaxis()->SetRangeUser(0, 1);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 36);
    frame_et_rec->SetYTitle("Probability of recoil jet");
    frame_et_rec->Draw("axis");
    for (int i = 0; i < (int)jet_pT_bins.size(); i++)
    {
        h_all_cluster_Et_max_b2bjet->GetYaxis()->SetRangeUser(jet_pT_bins[i], max_jet_pT);
        TH1D *h_jet_eff = h_all_cluster_Et_max_b2bjet->ProjectionX(Form("h_jet_eff_%d", i));
        h_jet_eff->Divide(h_all_cluster_Et);
        h_jet_eff->SetLineColor(col[i]);
        h_jet_eff->SetMarkerColor(col[i]);
        h_jet_eff->SetMarkerStyle(20);
        h_jet_eff->SetMarkerSize(1.5);
        h_jet_eff->Draw("same");
        h_eff.push_back(h_jet_eff);
    }

    float xpos(0.2), xpos2(0.915), ypos(0.885), ypos2(0.22), dy(0.054), dy1(0.07), fontsize(0.046);
    myText(xpos2,ypos-0*dy,1,strMC.c_str(),fontsize,1);
    myText(xpos2,ypos-1*dy,1,strleg3.c_str(),fontsize,1);
    myText(xpos2,ypos-2*dy,1,"|#it{#eta}^{jet}| < 0.6",fontsize,1);
    myText(xpos2,ypos-3*dy,1,strdijet.c_str(),fontsize,1);
    myText(xpos,ypos-0*dy,1,strleg1.c_str(),fontsize,0);
    myText(xpos,ypos-1*dy,1,strleg2.c_str(),fontsize,0);

    TLegend *leg = new TLegend(0.50, ypos2, xpos2, ypos2+(int)jet_pT_bins.size()*dy1);
    legStyle(leg, 0.20, 0.05);
    for (int i = 0; i < (int)jet_pT_bins.size(); i++)
    {
        leg->AddEntry(h_eff[i], Form("#it{p}_{T}^{jet} > %d GeV", (int)jet_pT_bins[i]), "lp");
    }
    leg->Draw("same");

    c1->SaveAs("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/PPG12-analysis-note/Figures/analysis/dijet_partner_probability.pdf");
}