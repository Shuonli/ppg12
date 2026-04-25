#include <TFile.h>
#include <TH2.h>

#include "/sphenix/u/bseidlitz/plotstyle/AtlasUtils.C"
#include "/sphenix/u/shuhang98/AtlasStyle.C"


void PlotIsoCorrelation(){
    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    TFile *fin = new TFile("MC_showershape_treeana.root", "READ");
    std::vector<TH2F*> hrecoiso_corrs;
    std::vector<std::string> hrecoiso_names = {"rhad22", "rhad33", "reta77", "reta", "rphi77", "rphi", "reta55", "rphi55", "weta", "wphi", "et1", "et2", "et3", "et4", "prob", "w72"};
    std::vector<TCanvas*> canvases;
    for (auto name : hrecoiso_names){
        TH2F *hrecoiso_corr = (TH2F*)fin->Get(Form("hrecoiso_%s", name.c_str()));
        hrecoiso_corrs.push_back(hrecoiso_corr);
        TCanvas *c = new TCanvas(Form("c_%s", name.c_str()), Form("c_%s", name.c_str()), 800, 600);
        canvases.push_back(c);
    }

    float iso_min = -1;
    float iso_max = 20;

    for (int i = 0; i < hrecoiso_corrs.size(); i++){
        canvases[i]->cd();
        canvases[i]->SetTopMargin(0.1);
        canvases[i]->SetRightMargin(0.15);
        //logz
        canvases[i]->SetLogz();
        hrecoiso_corrs[i]->GetXaxis()->SetRangeUser(iso_min, iso_max);
        hrecoiso_corrs[i]->SetXTitle("Cluster Iso_{ET} [GeV]");
        hrecoiso_corrs[i]->SetYTitle(hrecoiso_names[i].c_str());
        TProfile *prof = hrecoiso_corrs[i]->ProfileX();
        hrecoiso_corrs[i]->Draw("colz");
        prof->Draw("same");

        myText(0.1, 0.95, kBlack, "#bf{#it{sPHENIX}} Internal");
        float correlction = hrecoiso_corrs[i]->GetCorrelationFactor();
        myText(0.4, 0.95, kBlack, Form("Correlation: %.2f", correlction));

        

    }


}