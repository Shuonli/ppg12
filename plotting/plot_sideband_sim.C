#include"plotcommon.h"


void plot_sideband_sim(){

  init_plot();

  string savePath="/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/PPG12-analysis-note/Figures/";


  TFile* fdata = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_nom.root");
  //TFile* fdata = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_nom.root");
  TFile* fmc   = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_nom.root");

  TH1D* h_tight_iso_cluster       = (TH1D*) fdata->Get("h_tight_iso_cluster_0");
  TH1D* h_tight_noniso_cluster    = (TH1D*) fdata->Get("h_tight_noniso_cluster_0");
  TH1D* h_nontight_iso_cluster    = (TH1D*) fdata->Get("h_nontight_iso_cluster_0");
  TH1D* h_nontight_noniso_cluster = (TH1D*) fdata->Get("h_nontight_noniso_cluster_0");

  TH1D* h_tight_iso_cluster_mcSig       = (TH1D*) fmc->Get("h_tight_iso_cluster_0");
  TH1D* h_tight_noniso_cluster_mcSig    = (TH1D*) fmc->Get("h_tight_noniso_cluster_0");
  TH1D* h_nontight_iso_cluster_mcSig    = (TH1D*) fmc->Get("h_nontight_iso_cluster_0");
  TH1D* h_nontight_noniso_cluster_mcSig = (TH1D*) fmc->Get("h_nontight_noniso_cluster_0");

  TH1D* h_BoverA = (TH1D*) h_tight_noniso_cluster->Clone("h_BoverA");
  h_BoverA->Divide(h_tight_iso_cluster);
  TH1D* h_CoverA = (TH1D*) h_nontight_iso_cluster->Clone("h_CoverA");
  h_CoverA->Divide(h_tight_iso_cluster);
  TH1D* h_DoverA = (TH1D*) h_nontight_noniso_cluster->Clone("h_DoverA");
  h_DoverA->Divide(h_tight_iso_cluster);

  TH1D* h_BoverA_mcSig = (TH1D*) h_tight_noniso_cluster_mcSig->Clone("h_BoverA_mcSig");
  h_BoverA_mcSig->Divide(h_tight_iso_cluster_mcSig);
  TH1D* h_CoverA_mcSig = (TH1D*) h_nontight_iso_cluster_mcSig->Clone("h_CoverA_mcSig");
  h_CoverA_mcSig->Divide(h_tight_iso_cluster_mcSig);
  TH1D* h_DoverA_mcSig = (TH1D*) h_nontight_noniso_cluster_mcSig->Clone("h_DoverA_mcSig");
  h_DoverA_mcSig->Divide(h_tight_iso_cluster_mcSig);

  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  frame_et_rec->Draw("axis");
  frame_et_rec->GetXaxis()->SetRangeUser(10,30);
  frame_et_rec->GetYaxis()->SetRangeUser(500, 5e6);
  h_tight_iso_cluster->Draw("same hist");

  h_tight_noniso_cluster->Draw("same hist");
  h_tight_noniso_cluster->SetLineColor(kRed);

  h_nontight_iso_cluster->Draw("same hist");
  h_nontight_iso_cluster->SetLineColor(kBlue);

  h_nontight_noniso_cluster->Draw("same hist");
  h_nontight_noniso_cluster->SetLineColor(kMagenta);

  myText(0.5,0.9 ,1,strleg1.c_str(),0.04);
  myText(0.5,0.85,1,strleg2.c_str(),0.04);
  myText(0.5,0.80,1,strleg3.c_str(),0.04);
  myMarkerLineText(0.55,0.75, 0, kBlack, 0, kBlack, 1,"A: tight iso", 0.05, true);
  myMarkerLineText(0.55,0.70, 0, kRed, 0, kRed, 1,"B: tight noniso", 0.05, true);
  myMarkerLineText(0.55,0.65, 0, kBlue, 0, kBlue, 1,"C: nontight iso", 0.05, true);
  myMarkerLineText(0.55,0.60, 0, kMagenta, 0, kMagenta, 1,"D: nontight noniso", 0.05, true);

  gPad->SetLogy();
  c1->SaveAs(Form("%s/et_sbs_sim.pdf",savePath.c_str()));

    //unset logy
  gPad->SetLogy(0);
  TCanvas* c3 = new TCanvas("c3","c3",600,600);
  frame_et_rec->Draw("axis");
  frame_et_rec->GetYaxis()->SetRangeUser(7e-3,2);
  frame_et_rec->GetXaxis()->SetRangeUser(10,30);
  frame_et_rec->GetYaxis()->SetTitle("Ratio");
  h_BoverA->Draw("same hist");

  h_CoverA->Draw("same hist");
  h_CoverA->SetLineColor(kRed);

  h_DoverA->Draw("same hist");
  h_DoverA->SetLineColor(kBlue);

  myText(0.18,0.9 ,1,strleg1.c_str(),0.04);
  myText(0.18,0.85,1,strleg2.c_str(),0.04);
  myText(0.18,0.80,1,strleg3.c_str(),0.04);
  myMarkerLineText(0.55,0.75+0.15, 0, kBlack, 0, kBlack,1,"B/A: tight noniso"   , 0.05, true);
  myMarkerLineText(0.55,0.70+0.15, 0, kRed  , 0, kRed  ,1,"C/A: nontight iso", 0.05, true);
  myMarkerLineText(0.55,0.65+0.15, 0, kBlue , 0, kBlue ,1,"D/A: nontight noniso", 0.05, true);

  c3->SaveAs(Form("%s/et_sbs_ratio_sim.pdf",savePath.c_str()));


  //////////////////////////////////////////
  // Leakage fraction

  TCanvas* c4 = new TCanvas("c4","c4",600,600);
  frame_et_rec->Draw("axis");
  frame_et_rec->GetYaxis()->SetRangeUser(1e-3,1.3);
  frame_et_rec->GetYaxis()->SetTitle("Signal leakage");

  h_BoverA_mcSig->Draw("same hist");
  h_CoverA_mcSig->Draw("same hist");
  h_CoverA_mcSig->SetLineColor(kRed);
  h_DoverA_mcSig->Draw("same hist");
  h_DoverA_mcSig->SetLineColor(kBlue);

  myText(0.18,0.9 ,1,strleg1.c_str(),0.04);
  myText(0.18,0.85,1,strleg2.c_str(),0.04);
  myText(0.18,0.80,1,strleg3.c_str(),0.04);
  myText(0.18,0.75,1,strSigMC.c_str(),0.04);
  myMarkerLineText(0.53,0.75+0.15, 0, kBlack, 0, kBlack,1,"#it{N}^{sig}_{B}/#it{N}^{sig}_{A} tight noniso"   , 0.05, true);
  myMarkerLineText(0.53,0.68+0.15, 0, kRed  , 0, kRed  ,1,"#it{N}^{sig}_{C}/#it{N}^{sig}_{A} nontight iso", 0.05, true);
  myMarkerLineText(0.53,0.61+0.15, 0, kBlue , 0, kBlue ,1,"#it{N}^{sig}_{D}/#it{N}^{sig}_{A} nontight noniso", 0.05, true);

  //gPad->SetLogy();
  c4->SaveAs(Form("%s/leakage_fraction_et_sim.pdf",savePath.c_str()));



  ////////////////////////////////////////////////
  // Isolation ET

  TH1D* h_tight_isoET_pt[NptBins];
  TH1D* h_tight_isoET_mcSig_pt[NptBins];
  TH1D* h_nontight_isoET_pt[NptBins];
  for(int ipt=0; ipt<NptBins; ++ipt){
    h_tight_isoET_pt[ipt] = (TH1D*) fdata->Get(Form("h_tight_isoET_0_%d",ipt));
    h_tight_isoET_mcSig_pt[ipt] = (TH1D*) fmc->Get(Form("h_tight_isoET_0_%d",ipt));
    h_nontight_isoET_pt[ipt] = (TH1D*) fdata->Get(Form("h_nontight_isoET_0_%d",ipt));

    h_tight_isoET_pt[ipt]      ->Scale(1./h_tight_isoET_pt[ipt]->Integral());
    h_tight_isoET_mcSig_pt[ipt]->Scale(1./h_tight_isoET_mcSig_pt[ipt]->Integral());
    h_nontight_isoET_pt[ipt]   ->Scale(1./h_nontight_isoET_pt[ipt]->Integral());

    h_tight_isoET_pt[ipt]      ->Rebin(8);
    h_tight_isoET_mcSig_pt[ipt]->Rebin(4);
    h_nontight_isoET_pt[ipt]   ->Rebin(8);

    //scale by bin width
    h_tight_isoET_pt[ipt]->Scale(1./h_tight_isoET_pt[ipt]->GetBinWidth(1));
    h_tight_isoET_mcSig_pt[ipt]->Scale(1./h_tight_isoET_mcSig_pt[ipt]->GetBinWidth(1));
    h_nontight_isoET_pt[ipt]->Scale(1./h_nontight_isoET_pt[ipt]->GetBinWidth(1));

    //[ipt]->Scale(1./[ipt]->Integral());
  }

  for(int ipt=0; ipt<NptBins; ++ipt){

    TCanvas* c2 = new TCanvas("c2","c2",600,600);
    frame_isoET->Draw("axis");
    float max = h_tight_isoET_mcSig_pt[ipt]->GetMaximum();
    frame_isoET->GetYaxis()->SetRangeUser(0,max*1.3);

    h_tight_isoET_pt[ipt]   ->Draw("same hist");

    h_nontight_isoET_pt[ipt]->Draw("same hist");
    h_nontight_isoET_pt[ipt]->SetLineColor(kRed);

    h_tight_isoET_mcSig_pt[ipt]->Draw("same hist");
    h_tight_isoET_mcSig_pt[ipt]->SetLineColor(kBlue);

    myText(0.2,0.9,1,Form("%0.0f < #it{E}_{T}^{#gamma,rec} < %0.0f [GeV]",ptRanges[ipt],ptRanges[ipt+1]),0.04);

    myText          (0.50,0.9 -0.1,1,strleg1.c_str(),0.04);
    myText          (0.50,0.85-0.1,1,strleg2.c_str(),0.04);
    myText          (0.50,0.80-0.1,1,strleg3.c_str(),0.04);
    myMarkerLineText(0.55,0.75-0.1, 0, kBlack, 0, kBlack, 1,"tight data", 0.05, true);
    myMarkerLineText(0.55,0.70-0.1, 0, kRed, 0, kRed, 1,"nontight data", 0.05, true);
    myMarkerLineText(0.55,0.65-0.1, 0, kBlue, 0, kBlue, 1,"tight signal MC", 0.05, true);
    //myMarkerLineText(0.55,0.60-0.1, 0, kMagenta, 0, kMagenta, 1,"", 0.05, true);

    c2->SaveAs(Form("%s/iso_ET_sim_pt%d.pdf",savePath.c_str(),ipt));

  }

}
