#include "plotcommon.h"

void CONF_plots(string tune = "nom", int lowbin = 6, int highbin = 8)
{



  init_plot();
  string savePath = "./figures";
  TFile *fdata = new TFile(
    Form("/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_%s.root", tune.data()));
  TFile *fmc = new TFile(
    Form("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_%s.root", tune.data()));

  /////////////
  // get plots

  TH1D *h_tight_isoET_pt[NptBins];
  TH1D *h_tight_isoET_mcSig_pt[NptBins];
  TH1D *h_nontight_isoET_pt[NptBins];
  for (int ipt = 0; ipt < NptBins; ++ipt)
  {
    h_tight_isoET_pt[ipt] = (TH1D *)fdata->Get(Form("h_tight_isoET_0_%d", ipt));
    h_tight_isoET_mcSig_pt[ipt] = (TH1D *)fmc->Get(Form("h_tight_isoET_0_%d", ipt));
    h_nontight_isoET_pt[ipt] = (TH1D *)fdata->Get(Form("h_nontight_isoET_0_%d", ipt));
  }

  ///////
  // rebin in pt

  int sumBlow = lowbin;
  int sumBhigh = highbin;
  TH1D *h_tight_isoET = (TH1D *)h_tight_isoET_pt[sumBlow]->Clone("h_tight_isoET");
  TH1D *h_tight_isoET_mcSig = (TH1D *)h_tight_isoET_mcSig_pt[sumBlow]->Clone("h_tight_isoET_mcSig");
  TH1D *h_nontight_isoET = (TH1D *)h_nontight_isoET_pt[sumBlow]->Clone("h_nontight_isoET_pt");

  for (int ib = sumBlow + 1; ib <= sumBhigh; ib++)
  {
    h_tight_isoET->Add(h_tight_isoET_pt[ib]);
    h_tight_isoET_mcSig->Add(h_tight_isoET_mcSig_pt[ib]);
    h_nontight_isoET->Add(h_nontight_isoET_pt[ib]);
  }

  int nbinsOrig = h_tight_isoET->GetNbinsX();
  std::vector<double> newEdges;
  newEdges.push_back(h_tight_isoET->GetBinLowEdge(1)); // start edge

  int i = 1;
  while (i <= nbinsOrig)
  {
    // Decide the grouping size based on the lower edge of the current bin.
    int groupSize = (h_tight_isoET->GetBinLowEdge(i) < 2.5) ? 2 : 5;
    int last = i + groupSize - 1;
    if (last > nbinsOrig)
      last = nbinsOrig; // Do not exceed available bins

    // Get the upper edge of the last bin in the group.
    double upperEdge = h_tight_isoET->GetBinLowEdge(last) + h_tight_isoET->GetBinWidth(last);
    newEdges.push_back(upperEdge);
    i += groupSize;
  }

  int Nnew = newEdges.size() - 1; // new number of bins = number of edges - 1

  // Convert newEdges vector to an array.
  std::vector<double> edges = newEdges; // if needed, you can use &edges[0] directly
  TH1D *h_tight_isoET_rebinned = (TH1D *)h_tight_isoET->Rebin(Nnew, "h_tight_isoET_rebinned", &edges[0]);
  TH1D *h_tight_isoET_mcSig_rebinned = (TH1D *)h_tight_isoET_mcSig->Rebin(Nnew, "h_tight_isoET_mcSig_rebinned", &edges[0]);
  TH1D *h_nontight_isoET_rebinned = (TH1D *)h_nontight_isoET->Rebin(Nnew, "h_tight_isoET_rebinned", &edges[0]);

  h_tight_isoET = h_tight_isoET_rebinned;
  h_tight_isoET_mcSig = h_tight_isoET_mcSig_rebinned;
  h_nontight_isoET = h_nontight_isoET_rebinned;

  h_tight_isoET->Scale(1., "width");
  h_tight_isoET_mcSig->Scale(1., "width");
  h_nontight_isoET->Scale(1., "width");

  int rbf = 2;

  ////////////
  // normalize
  float ptcutnorm = 4;
  int ptB = h_tight_isoET->FindBin(ptcutnorm);
  int ptBs = h_tight_isoET->GetNbinsX();
  float normTight = h_tight_isoET->Integral(ptB, ptBs);
  float normNontight = h_nontight_isoET->Integral(ptB, ptBs);
  h_nontight_isoET->Scale(normTight / normNontight);
  float dataSig = h_tight_isoET->Integral() - h_nontight_isoET->Integral();
  h_tight_isoET_mcSig->Scale(dataSig / h_tight_isoET_mcSig->Integral());

  TH1D *h_tight_isoET_sub = (TH1D *)h_tight_isoET->Clone("h_tight_isoET_sub");
  h_tight_isoET_sub->Add(h_nontight_isoET, -1);

  //////////
  // plot

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 560);
  h_tight_isoET->Draw("ex0");
  h_tight_isoET->GetXaxis()->SetRangeUser(-1, 15);
  /*
    h_tight_isoET_mcSig->SetFillColorAlpha(kBlue, 0.8);
    h_tight_isoET_mcSig->SetFillStyle(1001);
    h_tight_isoET_mcSig->SetLineColor(kBlue-1);
    h_tight_isoET_mcSig->SetLineWidth(2);
    h_tight_isoET_mcSig->Draw("same hist F");

    h_nontight_isoET->SetFillColorAlpha(kRed, 0.8);
    h_nontight_isoET->SetFillStyle(1001);
    h_nontight_isoET->SetLineColor(kRed-1);
    h_nontight_isoET->SetLineWidth(2);
    h_nontight_isoET->Draw("same hist F");
  */

  h_tight_isoET_mcSig->SetFillColorAlpha(kBlue, 0.3);
  h_tight_isoET_mcSig->SetFillStyle(1001);
  h_tight_isoET_mcSig->SetLineColor(kBlue - 1);
  h_tight_isoET_mcSig->SetLineWidth(2);

  h_nontight_isoET->SetFillColorAlpha(kRed, 0.3);
  h_nontight_isoET->SetFillStyle(1001);
  h_nontight_isoET->SetLineColor(kRed - 1);
  h_nontight_isoET->SetLineWidth(2);

  THStack *hs = new THStack("hs", "");
  hs->Add(h_nontight_isoET);
  hs->Add(h_tight_isoET_mcSig);
  hs->Draw("hist same");

  h_tight_isoET->Draw("same ex0");

  h_tight_isoET->Draw("same axis");
  float binWidth = h_tight_isoET->GetBinWidth(1);
  // h_tight_isoET->SetYTitle(Form("Counts / %.1f", binWidth));
  h_tight_isoET->SetYTitle("Counts / Bin Width");
  h_tight_isoET->SetXTitle("#it{E}_{T}^{iso} [GeV]");
  h_tight_isoET->GetXaxis()->SetTitleOffset(1.2);

  string st_etbin = Form("%0.0f < #kern[-0.25]{#it{E}_{T}^{#gamma}} < %0.0f GeV", ptRanges[sumBlow], ptRanges[sumBhigh + 1]);
  float shift(0.05), xpos(0.915), ypos(0.875), dy(0.054), dy1(0.06), fontsize(0.046), fontsize1(0.048);
  myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize1, 1);
  myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 1);
  myText(xpos, ypos - 2 * dy, 1, Form("%s", st_etbin.data()), fontsize, 1);
  myText(xpos, ypos - 3 * dy, 1, Form("%s", strleg3.c_str()), fontsize, 1);
  int nLegend(3), legOffset(3);
  TLegend *l1 = new TLegend(0.51, ypos - (nLegend + legOffset) * dy1, xpos, ypos - legOffset * dy + 0.03);
  legStyle(l1, 0.20, fontsize);
  l1->AddEntry(h_tight_isoET, "Data (Signal)", "pe");
  l1->AddEntry(h_nontight_isoET, "Data (Background)", "f");
  l1->AddEntry(h_tight_isoET_mcSig, "Signal MC", "f");
  l1->Draw("same");

  /*
 float shift(0.05), xpos(0.88), ypos(0.9), dy(0.05), fontsize(0.04);
  myText          (0.35+shift,ypos-0*dy,1,strleg1.c_str(),fontsize,0);
  myText          (0.35+shift,ypos-1*dy,1,strleg2.c_str(),fontsize,0);
  myText          (0.35+shift,ypos-2*dy,1,strleg3.c_str(),fontsize,0);
  myText          (0.50+shift,ypos-3*dy,1,Form("%0.0f < #it{E}_{T}^{#gamma} < %0.0f GeV",ptRanges[sumBlow],ptRanges[sumBhigh+1]),0.04);
  myMarkerLineText(0.38+shift,ypos-4*dy, 1.5, kBlack, 20, kBlack, 1,"Data Signal Sel.", 0.05, true);
  myOnlyBoxText   (0.38+shift,ypos-5*dy, 1, kRed, kRed-1, 0,"Data Background Sel.", 0.05, 1001);
  myOnlyBoxText   (0.38+shift,ypos-6*dy, 1, kBlue, kBlue-1, 1,"Signal MC Signal Sel.", 0.05, 1001);
  */

  c1->SaveAs(Form("%s/h1D_iso_%s_%d_%d.pdf", savePath.c_str(), tune.c_str(), sumBlow, sumBhigh));

  TCanvas *c2 = new TCanvas("c2", "c2", 600, 560);
  // fit h_tight_isoET_sub and h_tight_isoET_mcSig
  /*
  TF1 *f1 = new TF1("f1", "gaus", 0.2, 1.5);

  */
  TF1 *f1 = new TF1("f1", "crystalball", 0.2, 4);

  f1->SetParameter(0, 10000); // Constant
  f1->SetParameter(1, 1);
  f1->SetParLimits(1, 0, 1.3); // Mean (mu) limits
  f1->SetParLimits(2, 0.01, 1.0); // Standard deviation (sigma)
  f1->SetParLimits(3, -3, 0);      // Alpha
  f1->SetParameter(4, 1.0);       // n (exponential tail power)

  h_tight_isoET_sub->Fit(f1, "REM");
  h_tight_isoET_sub->Fit(f1, "REM");

  TF1 *f2 = new TF1("f2", "crystalball", 0.2, 4);

  f2->SetParameter(0, 10000); // Constant
  f2->SetParameter(1, 1);
  f2->SetParLimits(1, 0, 1.3); // Mean (mu) limits
  f2->SetParLimits(2, 0.01, 1.0); // Standard deviation (sigma)
  f2->SetParLimits(3, -3, 0);      // Alpha
  f2->SetParameter(4, 1.0);       // n (exponential tail power)

  h_tight_isoET_mcSig->Fit(f2, "REM");
  h_tight_isoET_mcSig->Fit(f2, "REM");

  h_tight_isoET_mcSig->SetLineColor(kBlue);
  h_tight_isoET_mcSig->SetMarkerColor(kBlue);

  f2->SetLineColor(kBlue);

  h_tight_isoET_mcSig->GetXaxis()->SetRangeUser(-0.5, 7);
  h_tight_isoET_mcSig->SetXTitle("#it{E}_{T}^{iso} [GeV]");
  h_tight_isoET_mcSig->Draw("ex0");
  f2->Draw("same");

  f1->SetLineColor(kBlack);
  h_tight_isoET_sub->Draw("same ex0");
  f1->Draw("same");

  myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize1, 1);
  myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 1);
  myText(xpos, ypos - 2 * dy, 1, Form("%s", st_etbin.data()), fontsize, 1);

  float fit1mean = f1->GetParameter(1);
  float fit1sigma = f1->GetParameter(2);

  float fit2mean = f2->GetParameter(1);
  float fit2sigma = f2->GetParameter(2);

  nLegend = 2;
  l1 = new TLegend(0.36, ypos - (nLegend + legOffset) * dy1, xpos - 0.05, ypos - legOffset * dy + 0.01);
  legStyle(l1, 0.20, fontsize); 
  l1->AddEntry(h_tight_isoET_sub, Form("Data Signal #mu= %.2f, #sigma= %.2f", fit1mean, fit1sigma), "pe");
  l1->AddEntry(h_tight_isoET_mcSig, Form("Signal MC #mu= %.2f, #sigma = %.2f", fit2mean, fit2sigma), "pe");
  l1->Draw("same");

  c2->SaveAs(Form("%s/h1D_iso_fit_%s_%d_%d.pdf", savePath.c_str(), tune.c_str(), sumBlow, sumBhigh));
}