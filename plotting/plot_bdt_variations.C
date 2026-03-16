#include "plotcommon.h"

void plot_bdt_variations(bool useMC = true)
{
  init_plot();

  string basepath = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/";
  string savePath  = "figures";
  string mc_suffix = useMC ? "_mc" : "";

  // -----------------------------------------------------------------------
  // Variation table: {suffix, display label, color, marker style}
  // -----------------------------------------------------------------------
  struct VarEntry
  {
    string suffix;
    string label;
    int    color;
    int    marker;
  };

  vector<VarEntry> variations = {
    // No BDT
    //{"none",          "no BDT",             kBlack,     20},
    // BDT model variant
    {"modelbase",     "model base",         kViolet,    34},
    {"modelbase_v0",  "model base v0",      kViolet-1,  20},
    {"modelbase_v1",  "model base v1",      kViolet-2,  21},
    {"modelbase_v2",  "model base v2",      kViolet-3,  22},
    {"modelbase_E",   "model base E",       kMagenta,   23},
    {"modelbase_v0E", "model base v0E",     kMagenta+1, 29},
    {"modelbase_v1E", "model base v1E",     kMagenta+2, 30},
    {"modelbase_v2E", "model base v2E",     kMagenta+3, 33},
    {"modelbase_v3E", "model base v3E",     kViolet+2,  24},
  };

  // -----------------------------------------------------------------------
  // Load nominal
  // -----------------------------------------------------------------------
  TFile *fnom = new TFile((basepath + "Photon_final_bdt_nom" + mc_suffix + ".root").c_str());
  if (!fnom || fnom->IsZombie())
  {
    cout << "Error: cannot open nominal file" << endl;
    return;
  }
  TH1F *h_nom = (TH1F *)fnom->Get("h_unfold_sub_result");
  if (!h_nom)
  {
    cout << "Error: h_unfold_sub_result not found in nominal file" << endl;
    return;
  }

  // -----------------------------------------------------------------------
  // Loop over variations: open file, store spectrum and compute ratio var/nom
  // -----------------------------------------------------------------------
  vector<TH1F *> h_var_vec;
  vector<TH1F *> h_ratio_vec;

  for (auto &v : variations)
  {
    string fname = basepath + "Photon_final_bdt_" + v.suffix + mc_suffix + ".root";
    TFile *fvar  = new TFile(fname.c_str());
    if (!fvar || fvar->IsZombie())
    {
      cout << "Warning: cannot open " << fname << ", skipping" << endl;
      h_var_vec.push_back(nullptr);
      h_ratio_vec.push_back(nullptr);
      continue;
    }
    TH1F *h_var = (TH1F *)fvar->Get("h_unfold_sub_result");
    if (!h_var)
    {
      cout << "Warning: h_unfold_sub_result not found in " << fname << ", skipping" << endl;
      h_var_vec.push_back(nullptr);
      h_ratio_vec.push_back(nullptr);
      continue;
    }
    h_var_vec.push_back((TH1F *)h_var->Clone(("h_var_" + v.suffix).c_str()));
    TH1F *h_ratio = (TH1F *)h_var->Clone(("h_ratio_" + v.suffix).c_str());
    h_ratio->Divide(h_nom);
    h_ratio_vec.push_back(h_ratio);
  }

  // -----------------------------------------------------------------------
  // Single-pad canvas
  // -----------------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1", "", 1000, 700);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.07);
  c1->SetBottomMargin(0.15);

  // Configure frame axes
  frame_et_truth->SetYTitle("var / nom");
  frame_et_truth->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");
  frame_et_truth->GetYaxis()->SetRangeUser(0.8, 1.5);
  frame_et_truth->GetXaxis()->SetRangeUser(pTmin, pTmax);
  frame_et_truth->GetYaxis()->SetNdivisions(506);
  frame_et_truth->GetXaxis()->SetNdivisions(505);
  frame_et_truth->GetXaxis()->CenterTitle();
  frame_et_truth->GetYaxis()->CenterTitle();

  frame_et_truth->Draw("axis");
  TGraph *lineone = new TGraph();
  lineone->SetPoint(0, -1e6, 1);
  lineone->SetPoint(1, 1e6, 1);
  lineone->SetLineColor(kBlack);
  lineone->SetLineStyle(7);
  lineone->Draw("L");

  // -----------------------------------------------------------------------
  // Draw variation histograms
  // -----------------------------------------------------------------------
  for (int i = 0; i < (int)variations.size(); i++)
  {
    if (!h_ratio_vec[i]) continue;
    h_ratio_vec[i]->SetMarkerStyle(variations[i].marker);
    h_ratio_vec[i]->SetMarkerColor(variations[i].color);
    h_ratio_vec[i]->SetLineColor(variations[i].color);
    h_ratio_vec[i]->SetMarkerSize(1.2);
    h_ratio_vec[i]->Draw("same hist p");
  }

  // -----------------------------------------------------------------------
  // Legend (3 columns)
  // -----------------------------------------------------------------------
  TLegend *leg = new TLegend(0.14, 0.64, 0.95, 0.92);
  leg->SetNColumns(3);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.033);

  for (int i = 0; i < (int)variations.size(); i++)
  {
    if (!h_ratio_vec[i]) continue;
    leg->AddEntry(h_ratio_vec[i], variations[i].label.c_str(), "lp");
  }
  leg->Draw();

  // -----------------------------------------------------------------------
  // sPHENIX labels
  // -----------------------------------------------------------------------
  myText(0.25, 0.25, 1, strleg1.c_str(), 0.040);
  myText(0.25, 0.20, 1, strleg2.c_str(), 0.040);
  myText(0.25, 0.15, 1, useMC ? "inclusive MC" : "data", 0.040);

  c1->SaveAs(Form("%s/bdt_variations_rel%s.pdf", savePath.c_str(), mc_suffix.c_str()));

  // -----------------------------------------------------------------------
  // Spectrum canvas: nom + all variations
  // -----------------------------------------------------------------------
  TCanvas *c2 = new TCanvas("c2", "", 1000, 700);
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.05);
  c2->SetTopMargin(0.07);
  c2->SetBottomMargin(0.15);
  c2->SetLogy();

  frame_et_truth->SetYTitle("d#sigma/d#it{E}_{T}^{#gamma} [pb/GeV]");
  frame_et_truth->GetYaxis()->SetRangeUser(1e0, 1e4);
  frame_et_truth->Draw("axis");

  h_nom->SetMarkerStyle(20);
  h_nom->SetMarkerColor(kBlack);
  h_nom->SetLineColor(kBlack);
  h_nom->SetMarkerSize(1.2);
  h_nom->Draw("same hist p");

  for (int i = 0; i < (int)variations.size(); i++)
  {
    if (!h_var_vec[i]) continue;
    h_var_vec[i]->SetMarkerStyle(variations[i].marker);
    h_var_vec[i]->SetMarkerColor(variations[i].color);
    h_var_vec[i]->SetLineColor(variations[i].color);
    h_var_vec[i]->SetMarkerSize(1.2);
    h_var_vec[i]->Draw("same hist p");
  }

  TLegend *leg2 = new TLegend(0.14, 0.64, 0.95, 0.92);
  leg2->SetNColumns(3);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.033);
  leg2->AddEntry(h_nom, "nominal", "lp");
  for (int i = 0; i < (int)variations.size(); i++)
  {
    if (!h_var_vec[i]) continue;
    leg2->AddEntry(h_var_vec[i], variations[i].label.c_str(), "lp");
  }
  leg2->Draw();

  myText(0.25, 0.25, 1, strleg1.c_str(), 0.040);
  myText(0.25, 0.20, 1, strleg2.c_str(), 0.040);
  myText(0.25, 0.15, 1, useMC ? "inclusive MC" : "data", 0.040);

  c2->SaveAs(Form("%s/bdt_variations_spectra%s.pdf", savePath.c_str(), mc_suffix.c_str()));
}
