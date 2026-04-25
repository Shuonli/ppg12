#include "/sphenix/u/shuhang98/AtlasUtils.C"
#include "/sphenix/u/shuhang98/AtlasStyle.C"

void PlotSubstruction(){
    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TFile* fin = new TFile("/sphenix/user/shuhangli/pi0pythiasim/macro/demo_photon_yield_CLUSTERINFO_CEMC_CNN_single_2_5_8_10.root");

    TGraphErrors* gyield = (TGraphErrors*)fin->Get("photonyield");
    TGraphErrors* grawyield = (TGraphErrors*)fin->Get("rawphotonyield");

    TCanvas* c1 = new TCanvas("c1", "Photon Yield Overlay", 800, 800);

    gyield->GetXaxis()->SetTitle("Raw E_{T} [GeV]");
    gyield->GetYaxis()->SetTitle("Raw Photon Counts");
    gyield->GetYaxis()->SetRangeUser(15, 3900);
    gyield->GetYaxis()->SetTitleOffset(1.7);
    
    gyield->SetLineColor(kRed);
    gyield->SetMarkerStyle(20);
    gyield->SetMarkerColor(kRed);

    
    gyield->Draw("AP");

    grawyield->SetLineColor(kBlue);
    grawyield->SetMarkerStyle(21);
    grawyield->SetMarkerColor(kBlue);
    grawyield->Draw("P SAME");

    myMarkerLineText(0.5, 0.85, 1, kBlue, 21, kBlue, 1, "w/o. sideband substruction");
    myMarkerLineText(0.5, 0.80, 1, kRed, 20, kRed, 1, "w/ sideband substruction");



}