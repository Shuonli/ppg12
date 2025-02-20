// author: Yeonju Go 
// date: 2024 Dec. 9th 
// 
// v4: streak cluster rejection
// v5: skim tree for data w/ jets 
// v6: various streak event rejection purpose shower shapes 
#include "../include/yjUtility.h"
// #include "setBranch_slimtree_v3.h"
// #include "setBranch_skim_v1.h"
#include "setBranch_skim_v2.h"

void makeHist_showerShapes_v6(
    // string infile = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run22/jet10/condorout_waveform/caloana.root" 
    // , string cap = "jet10_v6"
    // string infile = "/direct/sphenix+tg+tg01/jets/yeonjugo/ppg12/skimtree/data_photon_skim_eventSelection.root"
    string infile = "/direct/sphenix+tg+tg01/jets/yeonjugo/ppg12/skimtree/data462_2024p010.root"
    , string cap = "data_v7"
    , int IDcut = 6 // 0: no ID cut, 1: Prelim. ID, 2: CNN prob, 3: CNN prob + Prelim. ID, 4: streak event removal, 5: only streak events 
    , bool fillPart = 0
    , bool doJetCut = 0
    , bool doVZreweight = 1
    , bool doMBDreweight = 0
    , bool doIso = 0
    ){
  cout << " :::::: makeHist_showerShapes_v6.C :::::: " << endl;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  SetyjPadStyle();

  bool isBkg = false;
  if(infile.find("jet")!=string::npos) isBkg = true; 

  bool isData = false;
  if(cap.find("data")!=string::npos) isData = true; 

  if(isData && (IDcut > 0)) cap += Form("_IDcut%d", IDcut);

  if(fillPart) cap += "_partial";

  if(isData){
    // if(!doJetCut) cap += "_requireSameSideJet";
    if(doJetCut) cap += "_requireOppositeSideJet";
  }
  if(!isData && doVZreweight) cap += "_vzreweighted";
  if(!isData && doMBDreweight) cap += "_mbdreweighted";
  if(!isData && doIso) cap += "_iso";

  TFile* f1 = new TFile(infile.c_str());
  TTree* t1 = (TTree*) f1->Get("slimtree");
  // TChain* t1 = new TChain("slimtree");
  // t1->Add("/sphenix/tg/tg01/commissioning/CaloCalibWG/sli/ppg12/ana450/ana450/*.root");
  InitBranches(t1, isData);

  ///////////////////////////////////////////////////

  for(int ivz = 0; ivz< nVZ; ivz++){
    vzStr[ivz] = Form("Vz%dto%d", (int)vzBins[ivz], (int)vzBins[ivz+1]);
  }
  vzStr[nVZ] = Form("Vz%dto%d", (int)vzBins[0], (int)vzBins[nVZ]);

  for(int ipt = 0; ipt< nPT; ipt++){
    ptStr[ipt] = Form("Pt%dto%d", (int)ptBins[ipt], (int)ptBins[ipt+1]);
  }
  ptStr[nPT] = Form("Pt%dto%d", (int)ptBins[0], (int)ptBins[nPT]);

  for(int ieta = 0; ieta< nETA; ieta++){
    etaStr[ieta] = Form("Eta%dto%d", (int)(etaBins[ieta]*10), (int)(etaBins[ieta+1]*10));
  }
  etaStr[nETA] = Form("Eta%dto%d", (int)(etaBins[0]*10), (int)(etaBins[nETA]*10));

  if(isBkg) 
    classStr[nCLASS] = "NonpromptPho"; 
  else
    classStr[nCLASS] = "PromptPho"; 

  if(isData) 
    classStr[nCLASS] = "AllPho"; 

  TH1D* h1D_vz = new TH1D("h1D_vz", "", 100, -40, 40);
  TH2D* h2D_mbdnorthhit_vs_mbdsouthhit = new TH2D("h2D_mbdnorthhit_vs_mbdsouthhit", "", 20, 0, 20, 20, 0, 20);

  TH1D* h1D_pt[nVZ+1][nETA+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_eta[nVZ+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_phi[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];

  TH1D* h1D_clus_E[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_Et[nVZ+1][nETA+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_Eta[nVZ+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_Phi[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];

  TH1D* h1D_clus_prob[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_CNN_prob[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_e17_to_e77[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_e37_to_e77[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_e33_to_e37[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_e32_to_e35[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_e33_to_e35[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_e11_to_e33[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_e11_to_E[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_e33_to_E[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_hcalet33_to_ettot[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_ihcalet33_to_ettot[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_ohcalet33_to_ettot[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_hcalet22_to_ettot[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_ihcalet22_to_ettot[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_ohcalet22_to_ettot[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_detamax[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_dphimax[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_e1[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_e2[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_e3[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_e4[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_et1[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_et2[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_et3[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_et4[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_weta[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wphi[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_w32[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_w52[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_w72[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH2D* h2D_clus_w72_vs_wphi[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_iso2[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_iso3[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_iso4[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_iso4_emcal[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_iso4_hcalin[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_iso4_hcalout[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH2D* h2D_clus_Eta_vs_Phi[nVZ+1][nPT+1][nCLASS+1][nCONVERSION+1];

  ///////////////////////////////////////////
  // test variables
  TH1D* h1D_clus_wphi72[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wphi73[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_weta72[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_weta73[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wr72n[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wr73n[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wr72[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wr73[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wr[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wrr72n[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wrr73n[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wrr72[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wrr73[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wrr[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];

  TH1D* h1D_clus_wetan[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wphin[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_weta_cog[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wphi_cog[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_weta_cogx[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wphi_cogx[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wrn[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wr_cog[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH1D* h1D_clus_wr_cogx[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  TH2D* h2D_clus_weta_cogx_vs_wr_cogx[nVZ+1][nETA+1][nPT+1][nCLASS+1][nCONVERSION+1];
  

  const int nbins_ratio = 100;
  const float min_ratio = -0.1;
  const float max_ratio = 1.05;
  const float max_ratio2 = 0.2;
  for(int ivz = 0; ivz< nVZ+1; ivz++){
    for(int ieta = 0; ieta< nETA+1; ieta++){
      for(int ipt = 0; ipt< nPT+1; ipt++){
        for(int icl = 0; icl< nCLASS+1; icl++){
          for(int iconv = 0; iconv< nCONVERSION+1; iconv++){
            string tempSt = Form("%s_%s_%s_%s_%s", vzStr[ivz].data(), etaStr[ieta].data(), ptStr[ipt].data(), classStr[icl].data(), convStr[iconv].data());
            h1D_clus_E[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_E_%s", tempSt.data()), ";Cluster E [GeV];", 80, 0, 80);  
            if(ipt==0){
              string tempSt2 = Form("%s_%s_%s_%s", vzStr[ivz].data(), etaStr[ieta].data(), classStr[icl].data(), convStr[iconv].data());
              h1D_clus_Et[ivz][ieta][icl][iconv] = new TH1D(Form("h1D_clus_Et_%s", tempSt2.data()), ";Cluster E_{T} [GeV];", 80, 0, 80);  
            }
            if(ieta==0){
              string tempSt2 = Form("%s_%s_%s_%s", vzStr[ivz].data(), ptStr[ipt].data(), classStr[icl].data(), convStr[iconv].data());
              h1D_clus_Eta[ivz][ipt][icl][iconv] = new TH1D(Form("h1D_clus_Eta_%s", tempSt2.data()), ";Cluster #eta;", 80, -1, 1);  
              h2D_clus_Eta_vs_Phi[ivz][ipt][icl][iconv] = new TH2D(Form("h2D_clus_Eta_vs_Phi_%s", tempSt2.data()), ";Cluster #eta;Cluster #phi", 80, -1, 1, 100, -1*TMath::Pi(), TMath::Pi());  
            }
            h1D_clus_Phi[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_Phi_%s", tempSt.data()), ";Cluster #phi;", 80, -1*TMath::Pi(), TMath::Pi());  

            /////////////////////////////////////
            h1D_clus_wphi72[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wphi72_%s", tempSt.data()), ";wphi72;",  120,0,2);  
            h1D_clus_wphi73[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wphi73_%s", tempSt.data()), ";wphi73;",  120,0,2);  
            h1D_clus_weta72[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_weta72_%s", tempSt.data()), ";weta72;",  120,0,2.5);  
            h1D_clus_weta73[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_weta73_%s", tempSt.data()), ";weta73;",  120,0,2.5);  
            h1D_clus_wr72n[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wr72n_%s", tempSt.data()), ";wr72n;",  120,0,5);  
            h1D_clus_wr73n[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wr73n_%s", tempSt.data()), ";wr73n;",  120,0,5);  
            h1D_clus_wr72[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wr72_%s", tempSt.data()), ";wr72;",  120,0,5);  
            h1D_clus_wr73[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wr73_%s", tempSt.data()), ";wr73;",  120,0,5);  
            h1D_clus_wr[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wr_%s", tempSt.data()), ";wr;",  200,0,10);  
            h1D_clus_wrr72n[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wrr72n_%s", tempSt.data()), ";wrr72n;",  120,0,10);  
            h1D_clus_wrr73n[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wrr73n_%s", tempSt.data()), ";wrr73n;",  120,0,10);  
            h1D_clus_wrr72[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wrr72_%s", tempSt.data()), ";wrr72;",  120,0,10);  
            h1D_clus_wrr73[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wrr73_%s", tempSt.data()), ";wrr73;",  120,0,10);  
            h1D_clus_wrr[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wrr_%s", tempSt.data()), ";wrr;",  120,0,10);  

            h1D_clus_wetan[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wetan_%s", tempSt.data()), ";wetan;",  100,0,2.5);  
            h1D_clus_wphin[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wphin_%s", tempSt.data()), ";wphin;",  100,0,2);  
            h1D_clus_weta_cog[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_weta_cog_%s", tempSt.data()), ";weta_cog;",  100,0,2.5);  
            h1D_clus_wphi_cog[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wphi_cog_%s", tempSt.data()), ";wphi_cog;",  100,0,2);  
            h1D_clus_weta_cogx[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_weta_cogx_%s", tempSt.data()), ";weta_cogx;",  100,0,2.5);  
            h1D_clus_wphi_cogx[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wphi_cogx_%s", tempSt.data()), ";wphi_cogx;",  100,0,2);  
            h1D_clus_wrn[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wrn_%s", tempSt.data()), ";wrn;",  100,0,5);  
            h1D_clus_wr_cog[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wr_cog_%s", tempSt.data()), ";wr_cog;",  100,0,5);  
            h1D_clus_wr_cogx[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wr_cogx_%s", tempSt.data()), ";wr_cogx;",  100,0,5);  
            h2D_clus_weta_cogx_vs_wr_cogx[ivz][ieta][ipt][icl][iconv] = new TH2D(Form("h2D_clus_weta_cogx_vs_wr_cogx_%s", tempSt.data()), ";weta_cogx;wr_cogx",  100,0,3, 100, 0, 5);  

            /////////////////////////////////////
            h1D_clus_prob[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_prob_%s", tempSt.data()), ";Cluster Probability;", 100, 0, 1);  
            h1D_clus_CNN_prob[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_CNN_prob_%s", tempSt.data()), ";Cluster CNN Probability;", 100, 0, 1);  

            h1D_clus_e17_to_e77[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_e17_to_e77_%s", tempSt.data()), ";E1x7/E7x7;", nbins_ratio, 0, max_ratio);  
            h1D_clus_e37_to_e77[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_e37_to_e77_%s", tempSt.data()), ";E3x7/E7x7;", nbins_ratio, 0.5, max_ratio);  
            h1D_clus_e33_to_e37[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_e33_to_e37_%s", tempSt.data()), ";E3x3/E3x7;", nbins_ratio, 0.6, max_ratio);  
            h1D_clus_e32_to_e35[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_e32_to_e35_%s", tempSt.data()), ";E3x2/E3x5;", nbins_ratio, 0.6, max_ratio);  
            h1D_clus_e33_to_e35[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_e33_to_e35_%s", tempSt.data()), ";E3x3/E3x5;", nbins_ratio, 0.7, max_ratio);  
            h1D_clus_e11_to_e33[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_e11_to_e33_%s", tempSt.data()), ";E1x1/E3x3;", nbins_ratio, 0,   max_ratio);  
            h1D_clus_e11_to_E[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_e11_to_E_%s", tempSt.data()), ";E1x1/E;", nbins_ratio, 0, max_ratio);  
            h1D_clus_e33_to_E[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_e33_to_E_%s", tempSt.data()), ";E3x3/E;", nbins_ratio, 0.5, max_ratio);  
            h1D_clus_hcalet33_to_ettot[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_hcalet33_to_ettot_%s", tempSt.data()), ";E_{T}^{HCal,3x3}/E_{T};", nbins_ratio, min_ratio, 0.1);  
            h1D_clus_hcalet22_to_ettot[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_hcalet22_to_ettot_%s", tempSt.data()), ";E_{T}^{HCal,2x2}/E_{T};", nbins_ratio, min_ratio, 0.1);  
            h1D_clus_ihcalet33_to_ettot[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_ihcalet33_to_ettot_%s", tempSt.data()), ";E_{T}^{iHCal,3x3}/E_{T};", nbins_ratio, min_ratio, 0.1);  
            h1D_clus_ihcalet22_to_ettot[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_ihcalet22_to_ettot_%s", tempSt.data()), ";E_{T}^{iHCal,2x2}/E_{T};", nbins_ratio, min_ratio, 0.1);  
            h1D_clus_ohcalet33_to_ettot[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_ohcalet33_to_ettot_%s", tempSt.data()), ";E_{T}^{oHCal,3x3}/E_{T};", nbins_ratio, min_ratio, 0.1);  
            h1D_clus_ohcalet22_to_ettot[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_ohcalet22_to_ettot_%s", tempSt.data()), ";E_{T}^{oHCal,2x2}/E_{T};", nbins_ratio, min_ratio, 0.1);  
            h1D_clus_detamax[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_detamax_%s", tempSt.data()), ";#Deltai#eta^{max};", 20,0,10);  
            h1D_clus_dphimax[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_dphimax_%s", tempSt.data()), ";#Deltai#phi^{max};", 20,0,10);  
            h1D_clus_e1[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_e1_%s", tempSt.data()), ";e1;", 80,0,20);  
            h1D_clus_e2[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_e2_%s", tempSt.data()), ";e2;", 80,0,20);  
            h1D_clus_e3[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_e3_%s", tempSt.data()), ";e3;", 80,0,20);  
            h1D_clus_e4[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_e4_%s", tempSt.data()), ";e4;", 80,0,20);  
            h1D_clus_et1[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_et1_%s", tempSt.data()), ";et1;", 80,0,1.1);  
            h1D_clus_et2[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_et2_%s", tempSt.data()), ";et2;", 80,0,1.1);  
            h1D_clus_et3[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_et3_%s", tempSt.data()), ";et3",  80,0,1.1);  
            h1D_clus_et4[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_et4_%s", tempSt.data()), ";et4",  80,0,0.6);  
            h1D_clus_weta[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_weta_%s", tempSt.data()), ";weta",  80,0,2);  
            h1D_clus_wphi[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_wphi_%s", tempSt.data()), ";wphi",  80,0,1.2);  
            h1D_clus_w32[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_w32_%s", tempSt.data()), ";w32",  80,0,1);  
            h1D_clus_w52[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_w52_%s", tempSt.data()), ";w52",  80,0,1.5);  
            h1D_clus_w72[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_w72_%s", tempSt.data()), ";w72",  80,0,2);  
            h2D_clus_w72_vs_wphi[ivz][ieta][ipt][icl][iconv] = new TH2D(Form("h2D_clus_w72_vs_wphi_%s", tempSt.data()), ";w72;wphi",  80,0,2, 80, 0, 1.2);  
            h1D_clus_iso2[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_iso2_%s", tempSt.data()), ";E_{T}^{iso, #DeltaR<0.2} [GeV]",  400,-15,15);  
            h1D_clus_iso3[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_iso3_%s", tempSt.data()), ";E_{T}^{iso, #DeltaR<0.3} [GeV]",  400,-15,15);  
            h1D_clus_iso4[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_iso4_%s", tempSt.data()), ";E_{T}^{iso, #DeltaR<0.4} [GeV]",  400,-15,15);  
            h1D_clus_iso4_emcal[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_iso4_emcal_%s", tempSt.data()), ";EMCal E_{T}^{iso, #DeltaR<0.4} [GeV]",  400,-15,15);  
            h1D_clus_iso4_hcalin[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_iso4_hcalin_%s", tempSt.data()), ";iHCal E_{T}^{iso, #DeltaR<0.4} [GeV]",  400,-15,15);  
            h1D_clus_iso4_hcalout[ivz][ieta][ipt][icl][iconv] = new TH1D(Form("h1D_clus_iso4_hcalout_%s", tempSt.data()), ";oHCal E_{T}^{iso, #DeltaR<0.4} [GeV]",  400,-15,15);  
          }
        }
      }
    }
  }

  double weightFactor = 1.;

  TF1* fvz = new TF1("fvz", "pol4", -100, 100); // vertex data/mc fit 
  fvz->SetParameters(0.733211, -0.00185783, 0.000765095, 1.25162e-06, 1.45691e-06);

  TF1* fmbd = new TF1("fmbd", "pol4", 0, 100); // mbd vertex reweight
  fmbd->SetParameters(0.623345, -0.0164323, 0.05306, -0.00665635, 0.000429544);//south, jet10
  ///////////////////////////////////////////////////
  // event loop
  int nTotEntries = t1->GetEntries();
  // if(!isData) nTotEntries = nTotEntries/2;
  // else nTotEntries = nTotEntries/2;
  int nTriggeredEvents = 0;
  // for(int ievent=0; ievent<10000; ievent++){
  int nNonIso = 0;
  for(int ievent=0; ievent<nTotEntries; ievent++){
    if(ievent%100000==0) cout << " events " << ievent << "/" << nTotEntries << " = " << (float)ievent/nTotEntries*100. << "% done..." << endl;
    // if(isData && runnumber==51576 && eventnumber==2261368) continue;
    // if(isData && runnumber==51909 && eventnumber==11450906) continue;

    t1->GetEntry(ievent);
    if(ncluster_CLUSTERINFO_CEMC_NO_SPLIT>5) continue;


    ///////////////////////////////////
    // trigger 
    if(isData){
      if( !(scaledtrigger[24] || scaledtrigger[25] ||scaledtrigger[26] || scaledtrigger[27]) ) continue;
    }

    ///////////////////////////////////
    // event selection: vertex cut
    if(abs(vertexz)>vzBins[nVZ]) continue;

    weightFactor = 1.;
    if(!isData && doVZreweight){ 
      weightFactor = weightFactor*fvz->Eval(vertexz); 
    }

    h1D_vz->Fill(vertexz, weightFactor);

    ///////////////////////////////////
    // event selection: MBD cut
    // int nmbdhitnorth = 0;
    // int nmbdhitsouth = 0;
    // for(int imbd=0;imbd<64;imbd++){
    //   if(mbdnorthq[imbd]>mbdqcut_onehit) nmbdhitnorth++;
    //   if(mbdsouthq[imbd]>mbdqcut_onehit) nmbdhitsouth++;
    // }
    // if(!(nmbdhitnorth >=1 && nmbdhitsouth>=1)) continue;
    if(!(mbdnorthhit >=1 && mbdsouthhit>=1)) continue;

    if(!isData && doMBDreweight){ 
      weightFactor = weightFactor*fmbd->Eval(mbdsouthhit); 
    }

    h2D_mbdnorthhit_vs_mbdsouthhit->Fill(mbdsouthhit, mbdnorthhit, weightFactor);

    nTriggeredEvents++;


    ////////////////////////////////////////
    // in MC, leading truth photon selection
    int phoIndex = -1; 
    int phoPt = -99; 
    if(!isData){
      for(int ip=0; ip<nparticles; ip++){
        if(particle_E[ip]<10) continue;
        if(isBkg && (particle_photonclass[ip]==1 || particle_photonclass[ip]==2)) continue;
        if(!isBkg && particle_pid[ip]!=22) continue;
        if(!isBkg && !(particle_photonclass[ip]==1 || particle_photonclass[ip]==2)) continue;

        if(particle_Pt[ip]>phoPt){
          phoIndex = ip;
          phoPt = particle_Pt[ip];
        }
      }
      if(phoIndex==-1) continue;
    }

    ////////////////////////////
    // cluster loop only for MC to find the matching 
    int leadingClusIndex = -1;
    float leadingClusEt = -99;
    if(!isData){
      for(int clusIndex=0; clusIndex<ncluster_CLUSTERINFO_CEMC_NO_SPLIT; clusIndex++){
        // MC particle - cluster matching
        // if(cluster_truthtrkID_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]==particle_trkid[phoIndex]){}
        if(getDR(particle_Eta[phoIndex], particle_Phi[phoIndex], cluster_Eta_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], cluster_Phi_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]) < truthRecoMatchingR) {
          if(leadingClusEt < cluster_Eta_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex])
          {
            leadingClusEt = cluster_Eta_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex];
            leadingClusIndex = clusIndex;
            // break;
          }
        } 
      }
    }
    if(!isData && leadingClusIndex==-1) continue;


    //////////////////////////////////////
    // cluster loop for both data and MC
    for(int clusIndex=0; clusIndex<ncluster_CLUSTERINFO_CEMC_NO_SPLIT; clusIndex++){
      if(!isData && clusIndex!=leadingClusIndex) continue;
      if(cluster_Et_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<minClusEtCut) continue;

      //////////////////////////////////////
      // jet cut 
      if(isData){
        bool passJetCut = false;
        for(int ij = 0; ij < njet; ij++){
          if(jet_Et[ij]<5) continue;
          if(getDPHI(cluster_Phi_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], jet_Phi[ij]) > 3*TMath::Pi()/4.){
            passJetCut = true; 
            break;
          } 
        }
        if(doJetCut){
          if(!passJetCut) continue; 
        }
      }
      //////////////////////////////////////
      // cout <<"335" << endl;

      int vzPos = findBinPos(fabs(vertexz), vzBins, nVZ);
      if(vzPos == -1) continue;

      int etaPos = findBinPos(fabs(cluster_Eta_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]), etaBins, nETA);
      // int etaPos = isData ? (findBinPos(fabs(cluster_Eta_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]), etaBins, nETA)) : (findBinPos(fabs(particle_Eta[phoIndex]), etaBins, nETA));
      if(etaPos == -1) continue;

      int ptPos = findBinPos(cluster_Et_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], ptBins, nPT);
      // int ptPos = isData ? (findBinPos(cluster_Et_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], ptBins, nPT)) : (findBinPos(particle_Pt[phoIndex], ptBins, nPT));
      if(ptPos == -1) continue;

      int clPos = isData ? 0 : particle_photonclass[phoIndex];
      int convPos = isData ? 0 : particle_converted[phoIndex];

      //////////////////////////////////////
      ///////////////////////////////////////
      // shower shape calculations
      TowerSum towerSum = calculateTowerInfo(cluster_e_array_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], cluster_e_array_idx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], cluster_status_array_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], cluster_ietacent_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], cluster_iphicent_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], 0.1);
      int nti_phi = towerSum.nearByClusIndex_phi;
      int nti_eta = towerSum.nearByClusIndex_eta;

      float weta = towerSum.weta_clus;
      float wphi = towerSum.wphi_clus;

      //w72 numerators 
      float wphi72n = (towerSum.esum_tw_phi7_clusc[1] + towerSum.esum_tw_phi7_clusc[nti_eta]);
      float wphi73n = (towerSum.esum_tw_phi7_clusc[1] + towerSum.esum_tw_phi7_clusc[0] + towerSum.esum_tw_phi7_clusc[2]);
      float weta72n = (towerSum.esum_tw_eta7_clusc[1] + towerSum.esum_tw_eta7_clusc[nti_phi]);
      float weta73n = (towerSum.esum_tw_eta7_clusc[1] + towerSum.esum_tw_eta7_clusc[0] + towerSum.esum_tw_eta7_clusc[2]);

      float esum_phi72 = (towerSum.esum_phi7_clusc[1] + towerSum.esum_phi7_clusc[nti_eta]);
      float esum_phi73 = (towerSum.esum_phi7_clusc[1] + towerSum.esum_phi7_clusc[0] + towerSum.esum_phi7_clusc[2]);
      float esum_eta72 = (towerSum.esum_eta7_clusc[1] + towerSum.esum_eta7_clusc[nti_phi]);
      float esum_eta73 = (towerSum.esum_eta7_clusc[1] + towerSum.esum_eta7_clusc[0] + towerSum.esum_eta7_clusc[2]);

      float wphi72 = wphi72n/esum_phi72;
      float wphi73 = wphi73n/esum_phi73;
      float weta72 = weta72n/esum_eta72;
      float weta73 = weta73n/esum_eta73;

      float wr72n = wphi72n/weta72n; 
      float wr73n = wphi73n/weta73n; 
      float wr72 = wphi72/weta72; 
      float wr73 = wphi73/weta73; 
      float wr = wphi/weta; 

      float wrr72n = weta72n/wphi72n; 
      float wrr73n = weta73n/wphi73n; 
      float wrr72 = weta72/wphi72; 
      float wrr73 = weta73/wphi73; 
      float wrr = weta/wphi; 
      float wr_cogx = cluster_wphi_cogx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_weta_cogx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex];

      float clus_et = cluster_E_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/TMath::CosH(cluster_Eta_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]); 
      float hcalet33 = cluster_ihcal_et33_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex] + cluster_ohcal_et33_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex];
      float hcalet22 = cluster_ihcal_et22_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex] + cluster_ohcal_et22_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex];

      // if(!fillPart){
        if(isData && IDcut>0){
          if(IDcut==1){
            if(cluster_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.02) continue;
            if(wr_cogx<0.4 && cluster_weta_cogx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>1) continue; // streak event removal
            if(cluster_e37_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_e77_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.9) continue;
            if(hcalet33/(hcalet33+clus_et)>0.1) continue;
            if(cluster_w72_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>0.75) continue;
            if(cluster_e11_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_E_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>0.98) continue;
          } else if(IDcut==2){
            if(cluster_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.02) continue;
            // if(cluster_CNN_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.8) continue;
          } else if(IDcut==3){
            if(cluster_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.02) continue;
            if(wr_cogx<0.4 && cluster_weta_cogx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>1) continue; // streak event removal
            if(cluster_e37_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_e77_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.9) continue;
            if(hcalet33/(hcalet33+clus_et)>0.1) continue;
            if(cluster_w72_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>0.75) continue;
            if(cluster_iso_03_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>1) continue;
            if(cluster_e11_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_E_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>0.98) continue;
          } else if(IDcut==4){
            if(wr_cogx<0.4 && cluster_weta_cogx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>1 && cluster_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.02) continue; // streak event removal
          } else if(IDcut==5){
            if(!(wr_cogx<0.4 && cluster_weta_cogx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>1 && cluster_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.02)) continue; // only streak events..
          } else if(IDcut==6){
            if(cluster_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.02) continue;
            if(wr_cogx<0.4 && cluster_weta_cogx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>1) continue; // streak event removal
            if(cluster_e37_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_e77_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.9) continue;
            if(hcalet33/(hcalet33+clus_et)>0.1) continue;
            if(cluster_w72_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>0.75) continue;
            if(cluster_iso_03_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>1) continue;
            if(cluster_e11_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_E_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>0.98) continue;
            if(cluster_CNN_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.8) continue;
          } else if(IDcut==7){
            if(cluster_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.02) continue;
            if(wr_cogx<0.4 && cluster_weta_cogx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>1) continue; // streak event removal
            if(cluster_e37_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_e77_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.9) continue;
            if(hcalet33/(hcalet33+clus_et)>0.1) continue;
            if(cluster_w72_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>0.75) continue;
            if(cluster_iso_03_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>1 || cluster_iso_03_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<-2) continue;
            if(cluster_e11_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_E_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>0.98) continue;
            if(cluster_CNN_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.8) continue;
          } else if(IDcut==8){
            if(cluster_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.4) continue;
            if(wr_cogx<0.4 && cluster_weta_cogx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>1) continue; // streak event removal
            if(cluster_e37_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_e77_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.9) continue;
            if(hcalet33/(hcalet33+clus_et)>0.1) continue;
            if(cluster_w72_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>0.75) continue;
            if(cluster_e11_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_E_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>0.98) continue;
            if(cluster_CNN_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.8) continue;

            if(cluster_iso_03_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>1 || cluster_iso_03_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<-2) continue;

          } else if(IDcut==9){
            if(cluster_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.4) continue;
            if(wr_cogx<0.4 && cluster_weta_cogx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>1) continue; // streak event removal
            if(cluster_e37_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_e77_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.9) continue;
            if(hcalet33/(hcalet33+clus_et)>0.1) continue;
            if(cluster_w72_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>0.75) continue;
            if(cluster_e11_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_E_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>0.98) continue;
            if(cluster_CNN_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<0.8) continue;

            if(cluster_iso_03_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<3) continue;
          }
        }

        if(!isData && doIso){ 
          if(cluster_iso_03_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]>1 || cluster_iso_03_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]<-2){
            nNonIso++;
            if(nNonIso<50)
            cout << "nonIso MC :: " << ievent << endl; 
            continue;
          }
        }
      // }

      h1D_clus_wphi72[vzPos][etaPos][ptPos][clPos][convPos]->Fill(wphi72, weightFactor);
      h1D_clus_wphi73[vzPos][etaPos][ptPos][clPos][convPos]->Fill(wphi73, weightFactor);
      h1D_clus_weta72[vzPos][etaPos][ptPos][clPos][convPos]->Fill(weta72, weightFactor);
      h1D_clus_weta73[vzPos][etaPos][ptPos][clPos][convPos]->Fill(weta73, weightFactor);
      h1D_clus_wr72n[vzPos][etaPos][ptPos][clPos][convPos]->Fill(wr72n, weightFactor);
      h1D_clus_wr73n[vzPos][etaPos][ptPos][clPos][convPos]->Fill(wr73n, weightFactor);
      h1D_clus_wr72[vzPos][etaPos][ptPos][clPos][convPos]->Fill(wr72, weightFactor);
      h1D_clus_wr73[vzPos][etaPos][ptPos][clPos][convPos]->Fill(wr73, weightFactor);
      h1D_clus_wr[vzPos][etaPos][ptPos][clPos][convPos]->Fill(wr, weightFactor);
      h1D_clus_wrr72n[vzPos][etaPos][ptPos][clPos][convPos]->Fill(wrr72n, weightFactor);
      h1D_clus_wrr73n[vzPos][etaPos][ptPos][clPos][convPos]->Fill(wrr73n, weightFactor);
      h1D_clus_wrr72[vzPos][etaPos][ptPos][clPos][convPos]->Fill(wrr72, weightFactor);
      h1D_clus_wrr73[vzPos][etaPos][ptPos][clPos][convPos]->Fill(wrr73, weightFactor);
      h1D_clus_wrr[vzPos][etaPos][ptPos][clPos][convPos]->Fill(wrr, weightFactor);

      h1D_clus_wetan[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_weta_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
      h1D_clus_wphin[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_wphi_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
      h1D_clus_weta_cog[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_weta_cog_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
      h1D_clus_wphi_cog[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_wphi_cog_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
      h1D_clus_weta_cogx[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_weta_cogx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
      h1D_clus_wphi_cogx[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_wphi_cogx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
      h1D_clus_wrn[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_wphi_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_weta_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
      h1D_clus_wr_cog[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_wphi_cog_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_weta_cog_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
      h1D_clus_wr_cogx[vzPos][etaPos][ptPos][clPos][convPos]->Fill(wr_cogx, weightFactor);
      h2D_clus_weta_cogx_vs_wr_cogx[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_weta_cogx_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], wr_cogx, weightFactor);


      if(!fillPart){
        // cout << "Pos : "<< vzPos << ", " << etaPos << ", " << ptPos << ", " << clPos << ", " << convPos << std::endl;

        h1D_clus_E[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_E_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_Et[vzPos][etaPos][clPos][convPos]->Fill(cluster_Et_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_Eta[vzPos][ptPos][clPos][convPos]->Fill(cluster_Eta_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h2D_clus_Eta_vs_Phi[vzPos][ptPos][clPos][convPos]->Fill(cluster_Eta_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], cluster_Phi_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_Phi[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_Phi_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_prob[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_CNN_prob[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_CNN_prob_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_e17_to_e77[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_e17_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_e77_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_e37_to_e77[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_e37_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_e77_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_e33_to_e37[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_e33_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_e37_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_e32_to_e35[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_e32_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_e35_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_e33_to_e35[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_e33_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_e35_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_e11_to_e33[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_e11_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_e33_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_e11_to_E[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_e11_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_E_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_e33_to_E[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_e33_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/cluster_E_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);

        h1D_clus_hcalet33_to_ettot[vzPos][etaPos][ptPos][clPos][convPos]->Fill(hcalet33/(hcalet33+clus_et), weightFactor);
        h1D_clus_hcalet22_to_ettot[vzPos][etaPos][ptPos][clPos][convPos]->Fill(hcalet22/(hcalet33+clus_et), weightFactor);
        h1D_clus_ihcalet33_to_ettot[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_ihcal_et33_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/(cluster_ihcal_et33_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]+clus_et), weightFactor);
        h1D_clus_ihcalet22_to_ettot[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_ihcal_et22_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/(cluster_ihcal_et22_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]+clus_et), weightFactor);
        h1D_clus_ohcalet33_to_ettot[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_ohcal_et33_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/(cluster_ohcal_et33_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]+clus_et), weightFactor);
        h1D_clus_ohcalet22_to_ettot[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_ohcal_et22_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]/(cluster_ohcal_et22_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex]+clus_et), weightFactor);

        h1D_clus_detamax[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_detamax_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_dphimax[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_dphimax_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);

        h1D_clus_e1[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_e1_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_e2[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_e2_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_e3[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_e3_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_e4[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_e4_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_et1[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_et1_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_et2[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_et2_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_et3[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_et3_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_et4[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_et4_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_weta[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_weta_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_wphi[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_wphi_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_w32[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_w32_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_w52[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_w52_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_w72[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_w72_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_iso2[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_iso_02_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_iso3[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_iso_03_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_iso4[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_iso_04_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_iso4_emcal[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_iso_04_emcal_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_iso4_hcalin[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_iso_04_hcalin_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);
        h1D_clus_iso4_hcalout[vzPos][etaPos][ptPos][clPos][convPos]->Fill(cluster_iso_04_hcalout_CLUSTERINFO_CEMC_NO_SPLIT[clusIndex], weightFactor);

      }
      h2D_clus_w72_vs_wphi[vzPos][etaPos][ptPos][clPos][convPos]->Fill(weta72, wphi, weightFactor);
    }//cluster for loop
  }

  cout << " >>>>> Triggered events " << nTriggeredEvents << " out of " << nTotEntries << endl; 

  /////////////////////////////////////////////////////////
  // merge bins
  cout << "Merge in CONVERSION " << endl;
  for(int ivz = 0; ivz< nVZ; ivz++){
    for(int ieta = 0; ieta< nETA; ieta++){
      for(int ipt = 0; ipt< nPT; ipt++){
        for(int icl = 0; icl< nCLASS; icl++){
          for(int iconv = 0; iconv< nCONVERSION; iconv++){
            h1D_clus_E[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_E[ivz][ieta][ipt][icl][iconv]); 
            if(ipt==0)
              h1D_clus_Et[ivz][ieta][icl][nCONVERSION] ->Add( h1D_clus_Et[ivz][ieta][icl][iconv]); 
            if(ieta==0){
              h1D_clus_Eta[ivz][ipt][icl][nCONVERSION] ->Add( h1D_clus_Eta[ivz][ipt][icl][iconv]); 
              h2D_clus_Eta_vs_Phi[ivz][ipt][icl][nCONVERSION] ->Add( h2D_clus_Eta_vs_Phi[ivz][ipt][icl][iconv]); 
            }
            /////////////////////
            h1D_clus_wphi72[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wphi72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi73[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wphi73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta72[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_weta72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta73[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_weta73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr72n[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wr72n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr73n[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wr73n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr72[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wr72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr73[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wr73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wr[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr72n[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wrr72n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr73n[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wrr73n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr72[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wrr72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr73[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wrr73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wrr[ivz][ieta][ipt][icl][iconv]); 

            h1D_clus_wetan[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wetan[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphin[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wphin[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta_cog[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_weta_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi_cog[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wphi_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta_cogx[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_weta_cogx[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi_cogx[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wphi_cogx[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrn[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wrn[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr_cog[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wr_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr_cogx[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_wr_cogx[ivz][ieta][ipt][icl][iconv]); 
            h2D_clus_weta_cogx_vs_wr_cogx[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h2D_clus_weta_cogx_vs_wr_cogx[ivz][ieta][ipt][icl][iconv]); 
            /////////////////////

            h1D_clus_Phi[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_Phi[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_prob[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_prob[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_CNN_prob[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_CNN_prob[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_e17_to_e77[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_e17_to_e77[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e37_to_e77[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_e37_to_e77[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e33_to_e37[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_e33_to_e37[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e32_to_e35[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_e32_to_e35[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e33_to_e35[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_e33_to_e35[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e11_to_e33[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_e11_to_e33[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e11_to_E[ivz][ieta][ipt][icl][nCONVERSION] ->Add(   h1D_clus_e11_to_E[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_e33_to_E[ivz][ieta][ipt][icl][nCONVERSION] ->Add(   h1D_clus_e33_to_E[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_hcalet33_to_ettot[ivz][ieta][ipt][icl][nCONVERSION] ->Add(  h1D_clus_hcalet33_to_ettot[ivz][ieta][ipt][icl][iconv] );
            h1D_clus_hcalet22_to_ettot[ivz][ieta][ipt][icl][nCONVERSION] ->Add(  h1D_clus_hcalet22_to_ettot[ivz][ieta][ipt][icl][iconv] );
            h1D_clus_ihcalet33_to_ettot[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_ihcalet33_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ihcalet22_to_ettot[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_ihcalet22_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ohcalet33_to_ettot[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_ohcalet33_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ohcalet22_to_ettot[ivz][ieta][ipt][icl][nCONVERSION] ->Add( h1D_clus_ohcalet22_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_detamax[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_detamax[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_dphimax[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_dphimax[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e1[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_e1[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e2[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_e2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e3[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_e3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e4[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_e4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et1[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_et1[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et2[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_et2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et3[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_et3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et4[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_et4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_weta[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_wphi[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w32[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_w32[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w52[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_w52[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w72[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_w72[ivz][ieta][ipt][icl][iconv]); 
            h2D_clus_w72_vs_wphi[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h2D_clus_w72_vs_wphi[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso2[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_iso2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso3[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_iso3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_iso4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_emcal[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_iso4_emcal[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_hcalin[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_iso4_hcalin[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_hcalout[ivz][ieta][ipt][icl][nCONVERSION] ->Add(            h1D_clus_iso4_hcalout[ivz][ieta][ipt][icl][iconv]); 
          }
        }
      }
    }
  }

  cout << "Merge in CLASS " << endl;
  for(int ivz = 0; ivz< nVZ; ivz++){
    for(int ieta = 0; ieta< nETA; ieta++){
      for(int ipt = 0; ipt< nPT; ipt++){
        for(int iconv = 0; iconv< nCONVERSION+1; iconv++){
          for(int icl = 0; icl< nCLASS; icl++){
            h1D_clus_E[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_E[ivz][ieta][ipt][icl][iconv]); 
            if(ipt==0)
              h1D_clus_Et[ivz][ieta][nCLASS][iconv] ->Add( h1D_clus_Et[ivz][ieta][icl][iconv]); 
            if(ieta==0){
              h1D_clus_Eta[ivz][ipt][nCLASS][iconv] ->Add( h1D_clus_Eta[ivz][ipt][icl][iconv]); 
              h2D_clus_Eta_vs_Phi[ivz][ipt][nCLASS][iconv] ->Add( h2D_clus_Eta_vs_Phi[ivz][ipt][icl][iconv]); 
            }
            h1D_clus_Phi[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_Phi[ivz][ieta][ipt][icl][iconv]); 

            /////////////////////
            h1D_clus_wphi72[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wphi72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi73[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wphi73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta72[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_weta72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta73[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_weta73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr72n[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wr72n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr73n[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wr73n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr72[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wr72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr73[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wr73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wr[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr72n[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wrr72n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr73n[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wrr73n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr72[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wrr72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr73[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wrr73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wrr[ivz][ieta][ipt][icl][iconv]); 

            h1D_clus_wetan[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wetan[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphin[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wphin[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta_cog[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_weta_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi_cog[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wphi_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta_cogx[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_weta_cogx[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi_cogx[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wphi_cogx[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrn[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wrn[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr_cog[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wr_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr_cogx[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_wr_cogx[ivz][ieta][ipt][icl][iconv]); 
            h2D_clus_weta_cogx_vs_wr_cogx[ivz][ieta][ipt][nCLASS][iconv] ->Add( h2D_clus_weta_cogx_vs_wr_cogx[ivz][ieta][ipt][icl][iconv]); 
            /////////////////////

            h1D_clus_prob[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_prob[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_CNN_prob[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_CNN_prob[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_e17_to_e77[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_e17_to_e77[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e37_to_e77[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_e37_to_e77[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e33_to_e37[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_e33_to_e37[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e32_to_e35[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_e32_to_e35[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e33_to_e35[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_e33_to_e35[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e11_to_e33[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_e11_to_e33[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e11_to_E[ivz][ieta][ipt][nCLASS][iconv] ->Add(   h1D_clus_e11_to_E[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_e33_to_E[ivz][ieta][ipt][nCLASS][iconv] ->Add(   h1D_clus_e33_to_E[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_hcalet33_to_ettot[ivz][ieta][ipt][nCLASS][iconv] ->Add(  h1D_clus_hcalet33_to_ettot[ivz][ieta][ipt][icl][iconv] );
            h1D_clus_hcalet22_to_ettot[ivz][ieta][ipt][nCLASS][iconv] ->Add(  h1D_clus_hcalet22_to_ettot[ivz][ieta][ipt][icl][iconv] );
            h1D_clus_ihcalet33_to_ettot[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_ihcalet33_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ihcalet22_to_ettot[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_ihcalet22_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ohcalet33_to_ettot[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_ohcalet33_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ohcalet22_to_ettot[ivz][ieta][ipt][nCLASS][iconv] ->Add( h1D_clus_ohcalet22_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_detamax[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_detamax[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_dphimax[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_dphimax[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e1[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_e1[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e2[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_e2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e3[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_e3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e4[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_e4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et1[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_et1[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et2[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_et2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et3[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_et3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et4[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_et4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_weta[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_wphi[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w32[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_w32[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w52[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_w52[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w72[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_w72[ivz][ieta][ipt][icl][iconv]); 
            h2D_clus_w72_vs_wphi[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h2D_clus_w72_vs_wphi[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso2[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_iso2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso3[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_iso3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_iso4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_emcal[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_iso4_emcal[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_hcalin[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_iso4_hcalin[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_hcalout[ivz][ieta][ipt][nCLASS][iconv] ->Add(            h1D_clus_iso4_hcalout[ivz][ieta][ipt][icl][iconv]); 
          }
        }
      }
    }
  }

  cout << "Merge in PT " << endl;
  for(int ivz = 0; ivz< nVZ; ivz++){
    for(int ieta = 0; ieta< nETA; ieta++){
      for(int icl = 0; icl< nCLASS+1; icl++){
        for(int iconv = 0; iconv< nCONVERSION+1; iconv++){
          for(int ipt = 0; ipt< nPT; ipt++){
            h1D_clus_E[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_E[ivz][ieta][ipt][icl][iconv]); 
            if(ieta==0){
              h1D_clus_Eta[ivz][nPT][icl][iconv] ->Add( h1D_clus_Eta[ivz][ipt][icl][iconv]); 
              h2D_clus_Eta_vs_Phi[ivz][nPT][icl][iconv] ->Add( h2D_clus_Eta_vs_Phi[ivz][ipt][icl][iconv]); 
            }
            h1D_clus_Phi[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_Phi[ivz][ieta][ipt][icl][iconv]); 

            /////////////////////
            h1D_clus_wphi72[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wphi72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi73[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wphi73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta72[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_weta72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta73[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_weta73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr72n[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wr72n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr73n[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wr73n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr72[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wr72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr73[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wr73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wr[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr72n[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wrr72n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr73n[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wrr73n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr72[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wrr72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr73[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wrr73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wrr[ivz][ieta][ipt][icl][iconv]); 

            h1D_clus_wetan[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wetan[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphin[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wphin[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta_cog[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_weta_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi_cog[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wphi_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta_cogx[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_weta_cogx[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi_cogx[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wphi_cogx[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrn[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wrn[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr_cog[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wr_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr_cogx[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_wr_cogx[ivz][ieta][ipt][icl][iconv]); 
            h2D_clus_weta_cogx_vs_wr_cogx[ivz][ieta][nPT][icl][iconv] ->Add( h2D_clus_weta_cogx_vs_wr_cogx[ivz][ieta][ipt][icl][iconv]); 
            /////////////////////

            h1D_clus_prob[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_prob[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_CNN_prob[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_CNN_prob[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_e17_to_e77[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_e17_to_e77[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e37_to_e77[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_e37_to_e77[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e33_to_e37[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_e33_to_e37[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e32_to_e35[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_e32_to_e35[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e33_to_e35[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_e33_to_e35[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e11_to_e33[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_e11_to_e33[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e11_to_E[ivz][ieta][nPT][icl][iconv] ->Add(   h1D_clus_e11_to_E[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_e33_to_E[ivz][ieta][nPT][icl][iconv] ->Add(   h1D_clus_e33_to_E[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_hcalet33_to_ettot[ivz][ieta][nPT][icl][iconv] ->Add(  h1D_clus_hcalet33_to_ettot[ivz][ieta][ipt][icl][iconv] );
            h1D_clus_hcalet22_to_ettot[ivz][ieta][nPT][icl][iconv] ->Add(  h1D_clus_hcalet22_to_ettot[ivz][ieta][ipt][icl][iconv] );
            h1D_clus_ihcalet33_to_ettot[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_ihcalet33_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ihcalet22_to_ettot[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_ihcalet22_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ohcalet33_to_ettot[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_ohcalet33_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ohcalet22_to_ettot[ivz][ieta][nPT][icl][iconv] ->Add( h1D_clus_ohcalet22_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_detamax[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_detamax[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_dphimax[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_dphimax[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e1[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_e1[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e2[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_e2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e3[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_e3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e4[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_e4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et1[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_et1[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et2[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_et2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et3[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_et3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et4[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_et4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_weta[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_wphi[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w32[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_w32[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w52[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_w52[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w72[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_w72[ivz][ieta][ipt][icl][iconv]); 
            h2D_clus_w72_vs_wphi[ivz][ieta][nPT][icl][iconv] ->Add(            h2D_clus_w72_vs_wphi[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso2[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_iso2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso3[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_iso3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_iso4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_emcal[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_iso4_emcal[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_hcalin[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_iso4_hcalin[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_hcalout[ivz][ieta][nPT][icl][iconv] ->Add(            h1D_clus_iso4_hcalout[ivz][ieta][ipt][icl][iconv]); 
          }
        }
      }
    }
  }

  cout << "Merge in ETA " << endl;
  for(int ivz = 0; ivz< nVZ; ivz++){
    for(int ipt = 0; ipt< nPT+1; ipt++){
      for(int icl = 0; icl< nCLASS+1; icl++){
        for(int iconv = 0; iconv< nCONVERSION+1; iconv++){
          for(int ieta = 0; ieta< nETA; ieta++){
            h1D_clus_E[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_E[ivz][ieta][ipt][icl][iconv]); 
            if(ipt==0)
              h1D_clus_Et[ivz][nETA][icl][iconv] ->Add( h1D_clus_Et[ivz][ieta][icl][iconv]); 
            h1D_clus_Phi[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_Phi[ivz][ieta][ipt][icl][iconv]); 

            /////////////////////
            h1D_clus_wphi72[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wphi72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi73[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wphi73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta72[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_weta72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta73[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_weta73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr72n[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wr72n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr73n[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wr73n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr72[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wr72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr73[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wr73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wr[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr72n[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wrr72n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr73n[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wrr73n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr72[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wrr72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr73[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wrr73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wrr[ivz][ieta][ipt][icl][iconv]); 

            h1D_clus_wetan[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wetan[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphin[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wphin[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta_cog[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_weta_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi_cog[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wphi_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta_cogx[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_weta_cogx[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi_cogx[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wphi_cogx[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrn[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wrn[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr_cog[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wr_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr_cogx[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_wr_cogx[ivz][ieta][ipt][icl][iconv]); 
            h2D_clus_weta_cogx_vs_wr_cogx[ivz][nETA][ipt][icl][iconv] ->Add( h2D_clus_weta_cogx_vs_wr_cogx[ivz][ieta][ipt][icl][iconv]); 
            /////////////////////

            h1D_clus_prob[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_prob[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_CNN_prob[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_CNN_prob[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_e17_to_e77[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_e17_to_e77[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e37_to_e77[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_e37_to_e77[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e33_to_e37[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_e33_to_e37[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e32_to_e35[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_e32_to_e35[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e33_to_e35[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_e33_to_e35[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e11_to_e33[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_e11_to_e33[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e11_to_E[ivz][nETA][ipt][icl][iconv] ->Add(   h1D_clus_e11_to_E[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_e33_to_E[ivz][nETA][ipt][icl][iconv] ->Add(   h1D_clus_e33_to_E[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_hcalet33_to_ettot[ivz][nETA][ipt][icl][iconv] ->Add(  h1D_clus_hcalet33_to_ettot[ivz][ieta][ipt][icl][iconv] );
            h1D_clus_hcalet22_to_ettot[ivz][nETA][ipt][icl][iconv] ->Add(  h1D_clus_hcalet22_to_ettot[ivz][ieta][ipt][icl][iconv] );
            h1D_clus_ihcalet33_to_ettot[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_ihcalet33_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ihcalet22_to_ettot[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_ihcalet22_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ohcalet33_to_ettot[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_ohcalet33_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ohcalet22_to_ettot[ivz][nETA][ipt][icl][iconv] ->Add( h1D_clus_ohcalet22_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_detamax[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_detamax[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_dphimax[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_dphimax[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e1[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_e1[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e2[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_e2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e3[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_e3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e4[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_e4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et1[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_et1[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et2[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_et2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et3[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_et3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et4[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_et4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_weta[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_wphi[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w32[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_w32[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w52[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_w52[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w72[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_w72[ivz][ieta][ipt][icl][iconv]); 
            h2D_clus_w72_vs_wphi[ivz][nETA][ipt][icl][iconv] ->Add(            h2D_clus_w72_vs_wphi[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso2[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_iso2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso3[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_iso3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_iso4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_emcal[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_iso4_emcal[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_hcalin[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_iso4_hcalin[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_hcalout[ivz][nETA][ipt][icl][iconv] ->Add(            h1D_clus_iso4_hcalout[ivz][ieta][ipt][icl][iconv]); 
          }
        }
      }
    }
  }

  cout << "Merge in VZ " << endl;
  for(int ieta = 0; ieta< nETA+1; ieta++){
    for(int ipt = 0; ipt< nPT+1; ipt++){
      for(int icl = 0; icl< nCLASS+1; icl++){
        for(int iconv = 0; iconv< nCONVERSION+1; iconv++){
          for(int ivz = 0; ivz< nVZ; ivz++){
            h1D_clus_E[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_E[ivz][ieta][ipt][icl][iconv]); 
            if(ipt==0)
              h1D_clus_Et[nVZ][ieta][icl][iconv] ->Add( h1D_clus_Et[ivz][ieta][icl][iconv]); 
            if(ieta==0){
              h1D_clus_Eta[nVZ][ipt][icl][iconv] ->Add( h1D_clus_Eta[ivz][ipt][icl][iconv]); 
              h2D_clus_Eta_vs_Phi[nVZ][ipt][icl][iconv] ->Add( h2D_clus_Eta_vs_Phi[ivz][ipt][icl][iconv]); 
            }
            h1D_clus_Phi[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_Phi[ivz][ieta][ipt][icl][iconv]); 

            /////////////////////
            h1D_clus_wphi72[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wphi72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi73[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wphi73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta72[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_weta72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta73[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_weta73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr72n[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wr72n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr73n[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wr73n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr72[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wr72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr73[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wr73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wr[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr72n[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wrr72n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr73n[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wrr73n[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr72[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wrr72[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr73[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wrr73[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrr[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wrr[ivz][ieta][ipt][icl][iconv]); 

            h1D_clus_wetan[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wetan[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphin[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wphin[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta_cog[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_weta_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi_cog[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wphi_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta_cogx[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_weta_cogx[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi_cogx[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wphi_cogx[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wrn[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wrn[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr_cog[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wr_cog[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wr_cogx[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_wr_cogx[ivz][ieta][ipt][icl][iconv]); 
            h2D_clus_weta_cogx_vs_wr_cogx[nVZ][ieta][ipt][icl][iconv] ->Add( h2D_clus_weta_cogx_vs_wr_cogx[ivz][ieta][ipt][icl][iconv]); 
            /////////////////////

            h1D_clus_prob[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_prob[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_CNN_prob[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_CNN_prob[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_e17_to_e77[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_e17_to_e77[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e37_to_e77[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_e37_to_e77[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e33_to_e37[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_e33_to_e37[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e32_to_e35[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_e32_to_e35[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e33_to_e35[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_e33_to_e35[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e11_to_e33[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_e11_to_e33[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e11_to_E[nVZ][ieta][ipt][icl][iconv] ->Add(   h1D_clus_e11_to_E[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_e33_to_E[nVZ][ieta][ipt][icl][iconv] ->Add(   h1D_clus_e33_to_E[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_hcalet33_to_ettot[nVZ][ieta][ipt][icl][iconv] ->Add(  h1D_clus_hcalet33_to_ettot[ivz][ieta][ipt][icl][iconv] );
            h1D_clus_hcalet22_to_ettot[nVZ][ieta][ipt][icl][iconv] ->Add(  h1D_clus_hcalet22_to_ettot[ivz][ieta][ipt][icl][iconv] );
            h1D_clus_ihcalet33_to_ettot[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_ihcalet33_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ihcalet22_to_ettot[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_ihcalet22_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ohcalet33_to_ettot[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_ohcalet33_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_ohcalet22_to_ettot[nVZ][ieta][ipt][icl][iconv] ->Add( h1D_clus_ohcalet22_to_ettot[ivz][ieta][ipt][icl][iconv]);
            h1D_clus_detamax[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_detamax[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_dphimax[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_dphimax[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e1[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_e1[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e2[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_e2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e3[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_e3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_e4[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_e4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et1[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_et1[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et2[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_et2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et3[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_et3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_et4[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_et4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_weta[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_weta[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_wphi[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_wphi[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w32[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_w32[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w52[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_w52[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_w72[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_w72[ivz][ieta][ipt][icl][iconv]); 
            h2D_clus_w72_vs_wphi[nVZ][ieta][ipt][icl][iconv] ->Add(            h2D_clus_w72_vs_wphi[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso2[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_iso2[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso3[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_iso3[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_iso4[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_emcal[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_iso4_emcal[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_hcalin[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_iso4_hcalin[ivz][ieta][ipt][icl][iconv]); 
            h1D_clus_iso4_hcalout[nVZ][ieta][ipt][icl][iconv] ->Add(            h1D_clus_iso4_hcalout[ivz][ieta][ipt][icl][iconv]); 
          }
        }
      }
    }
  }

  cout << "Save files " << endl;

  //////////////////////////////////////
  // save
  string fout_name = Form("/sphenix/u/yeonjugo/SE/sPHENIX_PPG12/showerShape/hist/showershape_hist_%s.root",cap.data());
  TFile* fout = new TFile(fout_name.c_str(),"recreate");
  h1D_vz->Write();
  h2D_mbdnorthhit_vs_mbdsouthhit->Write();
  for(int ivz = 0; ivz< nVZ+1; ivz++){
    for(int ieta = 0; ieta< nETA+1; ieta++){
      for(int ipt = 0; ipt< nPT+1; ipt++){
        for(int icl = 0; icl< nCLASS+1; icl++){
          for(int iconv = 0; iconv< nCONVERSION+1; iconv++){
            h1D_clus_wphi72[ivz][ieta][ipt][icl][iconv]->Write();
            h1D_clus_wphi73[ivz][ieta][ipt][icl][iconv]->Write();
            h1D_clus_weta72[ivz][ieta][ipt][icl][iconv]->Write();
            h1D_clus_weta73[ivz][ieta][ipt][icl][iconv]->Write();
            h1D_clus_wr72n[ivz][ieta][ipt][icl][iconv]->Write();
            h1D_clus_wr73n[ivz][ieta][ipt][icl][iconv]->Write();
            h1D_clus_wr72[ivz][ieta][ipt][icl][iconv]->Write();
            h1D_clus_wr73[ivz][ieta][ipt][icl][iconv]->Write();
            h1D_clus_wr[ivz][ieta][ipt][icl][iconv]->Write();
            h1D_clus_wrr72n[ivz][ieta][ipt][icl][iconv]->Write();
            h1D_clus_wrr73n[ivz][ieta][ipt][icl][iconv]->Write();
            h1D_clus_wrr72[ivz][ieta][ipt][icl][iconv]->Write();
            h1D_clus_wrr73[ivz][ieta][ipt][icl][iconv]->Write();
            h1D_clus_wrr[ivz][ieta][ipt][icl][iconv]->Write();

            h1D_clus_wetan[ivz][ieta][ipt][icl][iconv] ->Write();
            h1D_clus_wphin[ivz][ieta][ipt][icl][iconv] ->Write();
            h1D_clus_weta_cog[ivz][ieta][ipt][icl][iconv] ->Write();
            h1D_clus_wphi_cog[ivz][ieta][ipt][icl][iconv] ->Write();
            h1D_clus_weta_cogx[ivz][ieta][ipt][icl][iconv] ->Write();
            h1D_clus_wphi_cogx[ivz][ieta][ipt][icl][iconv] ->Write();
            h1D_clus_wrn[ivz][ieta][ipt][icl][iconv] ->Write();
            h1D_clus_wr_cog[ivz][ieta][ipt][icl][iconv] ->Write();
            h1D_clus_wr_cogx[ivz][ieta][ipt][icl][iconv] ->Write();
            h2D_clus_weta_cogx_vs_wr_cogx[ivz][ieta][ipt][icl][iconv] ->Write();

            if(!fillPart){
              h1D_clus_E[ivz][ieta][ipt][icl][iconv]->Write();
              if(ipt==0)
                h1D_clus_Et[ivz][ieta][icl][iconv]->Write();
              if(ieta==0){
                h1D_clus_Eta[ivz][ipt][icl][iconv]->Write();
                h2D_clus_Eta_vs_Phi[ivz][ipt][icl][iconv]->Write();
              }
              h1D_clus_Phi[ivz][ieta][ipt][icl][iconv]->Write();
              h1D_clus_prob[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_CNN_prob[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_e17_to_e77[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_e37_to_e77[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_e33_to_e37[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_e32_to_e35[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_e33_to_e35[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_e11_to_e33[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_e11_to_E[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_e33_to_E[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_hcalet33_to_ettot[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_hcalet22_to_ettot[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_ihcalet33_to_ettot[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_ihcalet22_to_ettot[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_ohcalet33_to_ettot[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_ohcalet22_to_ettot[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_detamax[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_dphimax[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_e1[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_e2[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_e3[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_e4[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_et1[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_et2[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_et3[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_et4[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_weta[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_wphi[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_w32[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_w52[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_w72[ivz][ieta][ipt][icl][iconv] ->Write();
              h2D_clus_w72_vs_wphi[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_iso2[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_iso3[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_iso4[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_iso4_emcal[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_iso4_hcalin[ivz][ieta][ipt][icl][iconv] ->Write();
              h1D_clus_iso4_hcalout[ivz][ieta][ipt][icl][iconv] ->Write();
            }
          }
        }
      }
    }
  }

  cout << " :::::: DONE:::: makeHist_showerShapes_v6.C for " << cap << " :::::: " << endl;
}
