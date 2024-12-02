// this macro simulate pi0 decay to 2 photons
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMath.h>
#include <TRandom3.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "/sphenix/u/bseidlitz/plotstyle/AtlasUtils.C"
#include "/sphenix/u/shuhang98/AtlasStyle.C"
// #include<cmath>


R__LOAD_LIBRARY(/sphenix/user/egm2153/calib_study/JetValidation/analysis/roounfold/libRooUnfold.so)
void resolutionsim()
{
    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    
    // output TFile
    TFile *f = new TFile("resosimout.root", "RECREATE");
    float xlow = 0;
    float xup = 50;
    int nbins = 10;
    float tau = 3.5;// Example decay constant, adjust as needed
    std::vector<std::pair<float, float>> reso = {{0.1, 0.1}
    , {0.1, 0.11}, {0.1, 0.09}, {0.05, 0.06}, {0.05, 0.04}
    };

    const int npairs = (int)reso.size();
    TH1D *hBias[npairs];
    TCanvas *cBias[npairs];
    int iplot = 0;
    for(auto mc_data : reso)
    {
        float reso_mc = mc_data.first;
        float reso_data = mc_data.second;
        //response matrix
        RooUnfoldResponse *response = new RooUnfoldResponse(nbins, xlow, xup);

    
    // for loop 
    TRandom3 *r = new TRandom3();
    for (int i = 0; i < 10000000; i++)
    {
        //trainning
        double photonpT = r->Exp(tau);
        double rerocorrection = r->Gaus(1, reso_mc);
        double recophotonpT = photonpT * rerocorrection;
        if(recophotonpT > xup){
            response->Miss(photonpT);
        }
        response->Fill(recophotonpT, photonpT);

    

 
    }
    //data pass
    TH1D *hMeasPT = new TH1D("hMeasPT", "", nbins, xlow, xup);
    TH1D *hTruthPT = new TH1D("hTruthPT", "", nbins, xlow, xup);
    
    for (int i = 0; i < 10000000; i++)
    {
        double photonpT = r->Exp(tau);
        double rerocorrection = r->Gaus(1, reso_data);
        double recophotonpT = photonpT * rerocorrection;
        hMeasPT->Fill(recophotonpT);
        hTruthPT->Fill(photonpT);
    }
    RooUnfoldBayes * unfold = new RooUnfoldBayes(response, hMeasPT, 2);
    TH1D *hRecoPT = (TH1D *)unfold->Hreco();
    hRecoPT->SetName("hRecoPT");

    //bias
    std::string biasname = "bias" + std::to_string(iplot);
    hBias[iplot] = (TH1D*)hRecoPT->Clone("hBias");
    hBias[iplot]->Add(hTruthPT, -1);
    hBias[iplot]->Divide(hTruthPT);
    hBias[iplot]->SetXTitle("p_{T}^{truth}");
    hBias[iplot]->SetYTitle("Bias");

    cBias[iplot] = new TCanvas(biasname.c_str(), biasname.c_str(), 800, 600);
    hBias[iplot]->GetYaxis()->SetRangeUser(-0.1, 0.1);
    hBias[iplot]->Draw();
/*
void myText(Double_t x,Double_t y,Color_t color, const char *text, Double_t tsize) {

  TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}

*/


    myText(0.2, 0.85, kBlack, Form("MC resolution: %.2f", reso_mc), 0.04);
    myText(0.2, 0.80, kBlack, Form("Data resolution: %.2f", reso_data), 0.04);
    
    


    iplot++;

    }
    f->cd();
    for(int i = 0; i < npairs; i++)
    {
        //hBias[i]->Write();
        cBias[i]->Write();
    }
    f->Write();
    f->Close();
}
