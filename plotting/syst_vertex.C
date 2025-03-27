#include "plotcommon.h"

void syst_vertex()
{
  init_plot();
  string varStr = "vertex";
  string savePath = "figures";

  TFile *f1 = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nom.root");
  TFile *f2 = new TFile(Form("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_vertex.root"));

  string legf1 = "nominal";
  string legf2 = "vertex var.";

  /////////////////////////////
  // plotting
  TH1F* h_res_1 = (TH1F *)f1->Get("h_unfold_sub_result");
  TH1F* h_res_2 = (TH1F *)f2->Get("h_unfold_sub_result");

  std::pair<TH1F *, TH1F *> h_dev_pair = calcDelta(h_res_2, h_res_1, "h_dev");

  TH1F *h_dev = h_dev_pair.first;
  TH1F *h_dev_rel = h_dev_pair.second;

  TCanvas *c1 = new TCanvas("can", "", 800, 889);
  c1->Divide(1, 2);
  gStyle->SetPalette(57);
  TPad *pad_1 = (TPad *)c1->cd(1);
  pad_1->SetPad(0, 0.4, 1, 1);
  pad_1->SetTopMargin(0.05);
  pad_1->SetLeftMargin(0.13);
  pad_1->SetBottomMargin(0.03);
  pad_1->SetRightMargin(0.08);
  pad_1->SetLogy();

  frame_et_rec->SetYTitle("d#sigma/dE_{T} [pb/GeV]");
  frame_et_rec->GetYaxis()->SetRangeUser(0.01, 1e3);
  frame_et_rec->GetXaxis()->SetRangeUser(10, 30);

  frame_et_rec->GetXaxis()->SetTitleOffset(0.98);
  frame_et_rec->GetYaxis()->SetTitleOffset(1.15);
  frame_et_rec->GetXaxis()->SetLabelSize(0.045);
  frame_et_rec->GetYaxis()->SetLabelSize(0.045);
  frame_et_rec->GetXaxis()->SetLabelOffset(2);
  frame_et_rec->GetXaxis()->CenterTitle();
  frame_et_rec->GetYaxis()->CenterTitle();
  frame_et_rec->GetXaxis()->SetNdivisions(505);

  frame_et_rec->Draw("axis");

  h_res_1->SetMarkerStyle(20);
  h_res_1->SetMarkerColor(kBlack);
  h_res_1->SetLineColor(kBlack);
  h_res_1->Draw("same");

  h_res_2->SetMarkerStyle(20);
  h_res_2->SetMarkerColor(kBlue);
  h_res_2->SetLineColor(kBlue);
  h_res_2->Draw("same");

  myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
  myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);

  myMarkerLineText(0.25, 0.25, 1, kBlack, 20, kBlack, 1, legf1.c_str(), 0.05, true);
  myMarkerLineText(0.25, 0.20, 1, kBlue, 20, kBlue, 1, legf2.c_str(), 0.05, true);

  TPad *pad_2 = (TPad *)c1->cd(2);
  pad_2->SetPad(0, 0, 1, 0.4);
  pad_2->SetTopMargin(0.02);
  pad_2->SetLeftMargin(0.13);
  pad_2->SetBottomMargin(0.25);
  pad_2->SetRightMargin(0.08);

  frame_et_truth->SetYTitle("Relative difference");
  frame_et_truth->GetYaxis()->SetNdivisions(506);
  frame_et_truth->GetYaxis()->SetRangeUser(-0.5, 0.5);
  frame_et_truth->GetXaxis()->SetRangeUser(10, 30);
  frame_et_truth->GetYaxis()->SetTitleOffset(frame_et_rec->GetYaxis()->GetTitleOffset()*4/6.);
  frame_et_truth->GetYaxis()->SetLabelOffset(frame_et_rec->GetYaxis()->GetLabelOffset()*4/6.);
  frame_et_truth->GetXaxis()->SetLabelSize(frame_et_rec->GetXaxis()->GetLabelSize()*6/4.);
  frame_et_truth->GetYaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize()*6/4.);
  frame_et_truth->GetXaxis()->SetTitleSize(frame_et_rec->GetXaxis()->GetTitleSize()*6/4.);
  frame_et_truth->GetYaxis()->SetTitleSize(frame_et_rec->GetYaxis()->GetTitleSize()*6/4.);
  frame_et_truth->GetXaxis()->SetNdivisions(505);


  frame_et_truth->Draw("axis");

  h_dev_rel->SetMarkerStyle(20);
  h_dev_rel->SetMarkerColor(kBlue);
  h_dev_rel->SetLineColor(kBlue);

  h_dev_rel->Draw("same hist p");

  linezero->Draw("L");

  c1->SaveAs(Form("%s/syst_spec_%s.pdf", savePath.c_str(), varStr.c_str()));


  TCanvas *c2 = new TCanvas("can", "", 900, 600);
  init_plot();
  frame_et_truth->SetYTitle("Relative difference");
  frame_et_truth->GetYaxis()->SetRangeUser(-0.1, 0.15);
  frame_et_truth->GetXaxis()->SetRangeUser(pTmin, pTmax);
  frame_et_truth->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");

  frame_et_truth->Draw("axis");
  linezero->Draw("L");

  h_dev_rel->SetMarkerStyle(20);
  h_dev_rel->SetMarkerColor(kBlue);
  h_dev_rel->SetLineColor(kBlue);

  h_dev_rel->Draw("same p hist");

  myText(0.6, 0.9, 1, strleg1.c_str(), 0.05);
  myText(0.6, 0.85, 1, strleg2.c_str(), 0.05);
  myText(0.6, 0.80, 1, legf2.c_str(), 0.05);

  c2->SaveAs(Form("%s/syst_rel_%s.pdf", savePath.c_str(), varStr.c_str()));

  TH1F* h_dev_low = (TH1F *)h_dev->Clone("h_dev_low");
  TH1F* h_dev_high = (TH1F *)h_dev->Clone("h_dev_high");

  TH1F* h_dev_rel_low = (TH1F *)h_dev_rel->Clone("h_dev_rel_low");
  TH1F* h_dev_rel_high = (TH1F *)h_dev_rel->Clone("h_dev_rel_high");

  for (int i = 1; i <= h_dev->GetNbinsX(); i++) {
    h_dev_low->SetBinContent(i, abs(h_dev->GetBinContent(i)));
    h_dev_high->SetBinContent(i, abs(h_dev->GetBinContent(i)));

    h_dev_rel_low->SetBinContent(i, abs(h_dev_rel->GetBinContent(i)));
    h_dev_rel_high->SetBinContent(i, abs(h_dev_rel->GetBinContent(i))); 
  }

  TFile *systOut = new TFile(Form("rootFiles/syst_%s.root", varStr.c_str()), "RECREATE");
  h_dev->Write();
  h_dev_rel->Write();
  h_dev_low->Write();
  h_dev_high->Write();
  h_dev_rel_low->Write();
  h_dev_rel_high->Write();
  systOut->Close();
}
