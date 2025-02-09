#include "/sphenix/user/hanpuj/plotstyle/AtlasStyle.C"
#include "/sphenix/user/hanpuj/plotstyle/AtlasUtils.C"

void draw_1D_single_plot(TH1F *h_input,
                         bool do_rebin, int rebin_factor, bool do_normalize,
                         bool set_xrange, float xlow, float xhigh, bool set_logx,
                         bool set_yrange, float ylow, float yhigh, bool set_logy,
                         bool set_title, std::string xtitle, std::string ytitle,
                         bool do_data, std::string run_status,
                         bool set_text, std::vector<std::string> text, float xstart, float ystart, float size,
                         std::string output_name) {
  TH1F *h = (TH1F*)h_input->Clone("h");
  TCanvas *can = new TCanvas("can", "", 800, 600);
  gStyle->SetPalette(57);
  can->SetTopMargin(0.12);
  can->SetLeftMargin(0.13);
  can->SetBottomMargin(0.12);
  can->SetRightMargin(0.08);
  if (do_rebin) h->Rebin(rebin_factor);
  if (do_normalize) h->Scale(1.0/(float)h->GetBinWidth(1));
  if (set_xrange) h->GetXaxis()->SetRangeUser(xlow, xhigh);
  if (set_logx) can->SetLogx();
  if (set_yrange) h->GetYaxis()->SetRangeUser(ylow, yhigh);
  if (set_logy) can->SetLogy();
  if (set_title) {
    h->GetXaxis()->SetTitle(xtitle.c_str());
    h->GetYaxis()->SetTitle(ytitle.c_str());
  }
  h->GetXaxis()->SetTitleOffset(0.98);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->SetMarkerStyle(20);
  h->SetMarkerColor(kBlack);
  h->SetLineColor(kBlack);
  h->Draw();
  if (do_data) myText(0.065, 0.97, 1, "#bf{#it{sPHENIX}} Data", 0.04);
  else myText(0.065, 0.97, 1, "#bf{#it{sPHENIX}} Simulation", 0.04);
  myText(0.065, 0.92, 1, run_status.c_str(), 0.04);
  if (set_text) {
    for (int i = 0; i < text.size(); i++) {
      myText(xstart, ystart-i*0.05, 1, text[i].c_str(), size);
    }
  }
  can->SaveAs(output_name.c_str());
  delete can;
}

void draw_1D_multiple_plot_ratio(std::vector<TH1F*> h_input, std::vector<int> color, std::vector<int> markerstyle,
                                 bool do_rebin, int rebin_factor, bool do_normalize,
                                 bool set_xrange, float xlow, float xhigh, bool set_logx,
                                 bool set_yrange, float ylow, float yhigh, bool set_logy,
                                 bool set_ratioyrange, float ratiolow, float ratiohigh,
                                 bool set_title, std::string xtitle, std::string ytitle, std::string ratio_title,
                                 bool do_data, std::string run_status,
                                 bool set_text, std::vector<std::string> text, float xstart, float ystart, float size,
                                 bool set_legend, std::vector<std::string> legend, float xstart_legend, float ystart_legend, float size_legend,
                                 std::string output_name) {
  if (h_input.size() <= 1) {
    std::cout << "Error: no enough input histograms for draw_1D_multiple_plot_ratio" << std::endl;
    return;
  }
  std::vector<TH1F*> h;
  for (int i = 0; i < h_input.size(); ++i) {
    h.push_back((TH1F*)h_input[i]->Clone(Form("h_%d", i)));
  }
  TCanvas *can = new TCanvas("can", "", 800, 889);
  can->Divide(1, 2);
  gStyle->SetPalette(57);
  TPad *pad_1 = (TPad*)can->cd(1);
  pad_1->SetPad(0, 0.4, 1, 1);
  pad_1->SetTopMargin(0.12);
  pad_1->SetLeftMargin(0.13);
  pad_1->SetBottomMargin(0.035);
  pad_1->SetRightMargin(0.08);
  for (int i = 0; i < h.size(); ++i) {
    h.at(i)->SetMarkerStyle(markerstyle.at(i));
    h.at(i)->SetMarkerColor(color.at(i));
    h.at(i)->SetLineColor(color.at(i));
    if (do_rebin) h.at(i)->Rebin(rebin_factor);
    if (do_normalize) h.at(i)->Scale(1.0/(float)h.at(i)->GetBinWidth(1));
  }
  if (set_xrange) h.at(0)->GetXaxis()->SetRangeUser(xlow, xhigh);
  if (set_logx) pad_1->SetLogx();
  if (set_yrange) h.at(0)->GetYaxis()->SetRangeUser(ylow, yhigh);
  if (set_logy) pad_1->SetLogy();
  if (set_title) {
    h.at(0)->GetXaxis()->SetTitle(xtitle.c_str());
    h.at(0)->GetYaxis()->SetTitle(ytitle.c_str());
  }
  h.at(0)->GetXaxis()->SetTitleOffset(0.98);
  h.at(0)->GetYaxis()->SetTitleOffset(1.15);
  h.at(0)->GetXaxis()->SetLabelSize(0.045);
  h.at(0)->GetYaxis()->SetLabelSize(0.045);
  h.at(0)->GetXaxis()->SetLabelOffset(2);
  h.at(0)->GetXaxis()->CenterTitle();
  h.at(0)->GetYaxis()->CenterTitle();

  h.at(0)->Draw("p");
  for (int i = 0; i < h.size(); ++i) {
    if (i == 0) continue;
    h.at(i)->Draw("p same");
  }

  if (do_data) myText(0.065, 0.97, 1, "#bf{#it{sPHENIX}} Data", 0.04);
  else myText(0.065, 0.97, 1, "#bf{#it{sPHENIX}} Simulation", 0.04);
  myText(0.065, 0.92, 1, run_status.c_str(), 0.04);
  if (set_text) {
    for (int i = 0; i < text.size(); i++) {
      myText(xstart, ystart-i*0.05, 1, text[i].c_str(), size);
    }
  }
  if (set_legend) {
    for (int i = 0; i < legend.size(); i++) {
      myMarkerLineText(xstart_legend, ystart_legend-i*0.05, 1, color.at(i), markerstyle.at(i), color.at(i), 1, legend.at(i).c_str(), size_legend, true);
    }
  }

  TPad *pad_2 = (TPad*)can->cd(2);
  pad_2->SetPad(0, 0, 1, 0.4);
  pad_2->SetTopMargin(0.02);
  pad_2->SetLeftMargin(0.13);
  pad_2->SetBottomMargin(0.25);
  pad_2->SetRightMargin(0.08);
  std::vector<TH1F*> h_ratio;
  for (int i = 0; i < h.size()-1; ++i) {
    TH1F *h_temp = (TH1F*)h.at(i+1)->Clone(Form("h_ratio_%d", i));
    h_temp->Divide(h.at(0));
    h_ratio.push_back(h_temp);
    h_ratio.at(i)->SetMarkerStyle(markerstyle.at(i+1));
    h_ratio.at(i)->SetMarkerColor(color.at(i+1));
    h_ratio.at(i)->SetLineColor(color.at(i+1));
  }
  h.at(0)->Draw();
  if (set_title) {
    h_ratio.at(0)->GetXaxis()->SetTitle(xtitle.c_str());
    h_ratio.at(0)->GetYaxis()->SetTitle(ratio_title.c_str());
  } else {
    h_ratio.at(0)->GetXaxis()->SetTitle(h.at(0)->GetXaxis()->GetTitle());
    h_ratio.at(0)->GetYaxis()->SetTitle("Ratio");
  }
  if (set_xrange) h_ratio.at(0)->GetXaxis()->SetRangeUser(xlow, xhigh);
  if (set_ratioyrange) h_ratio.at(0)->GetYaxis()->SetRangeUser(ratiolow, ratiohigh);
  h_ratio.at(0)->GetXaxis()->CenterTitle();
  h_ratio.at(0)->GetYaxis()->CenterTitle();
  h_ratio.at(0)->GetYaxis()->SetTitleOffset(h.at(0)->GetYaxis()->GetTitleOffset()*4/6.);
  h_ratio.at(0)->GetYaxis()->SetLabelOffset(h.at(0)->GetYaxis()->GetLabelOffset()*4/6.);
  h_ratio.at(0)->GetXaxis()->SetLabelSize(h.at(0)->GetXaxis()->GetLabelSize()*6/4.);
  h_ratio.at(0)->GetYaxis()->SetLabelSize(h.at(0)->GetYaxis()->GetLabelSize()*6/4.);
  h_ratio.at(0)->GetXaxis()->SetTitleSize(h.at(0)->GetXaxis()->GetTitleSize()*6/4.);
  h_ratio.at(0)->GetYaxis()->SetTitleSize(h.at(0)->GetYaxis()->GetTitleSize()*6/4.);
  h_ratio.at(0)->Draw();
  for (int i = 0; i < h_ratio.size(); ++i) {
    if (i == 0) continue;
    h_ratio.at(i)->Draw("EX0 same");
  }

  TLine *line;
  if (set_xrange) line = new TLine(xlow, 1, xhigh, 1);
  else line = new TLine(h.at(0)->GetXaxis()->GetBinLowEdge(1), 1, h.at(0)->GetXaxis()->GetBinUpEdge(h.at(0)->GetNbinsX()), 1);
  line->SetLineColor(kBlack);
  line->SetLineStyle(3);
  line->Draw("same");

  can->SaveAs(output_name.c_str());
  delete can;
  delete line;
}

void draw_1D_multiple_plot(std::vector<TH1F*> h_input, std::vector<int> color, std::vector<int> markerstyle,
                           bool do_rebin, int rebin_factor, bool do_normalize,
                           bool set_xrange, float xlow, float xhigh, bool set_logx,
                           bool set_yrange, float ylow, float yhigh, bool set_logy,
                           bool set_title, std::string xtitle, std::string ytitle,
                           bool do_data, std::string run_status,
                           bool set_text, std::vector<std::string> text, float xstart, float ystart, float size,
                           bool set_legend, std::vector<std::string> legend, float xstart_legend, float ystart_legend, float size_legend,
                           std::string output_name) {
  if (h_input.size() <= 1) {
    std::cout << "Error: no enough input histograms for draw_1D_multiple_plot" << std::endl;
    return;
  }
  std::vector<TH1F*> h;
  for (int i = 0; i < h_input.size(); ++i) {
    h.push_back((TH1F*)h_input[i]->Clone(Form("h_%d", i)));
  }
  TCanvas *can = new TCanvas("can", "", 800, 600);
  gStyle->SetPalette(57);
  can->SetTopMargin(0.12);
  can->SetLeftMargin(0.13);
  can->SetBottomMargin(0.12);
  can->SetRightMargin(0.08);
  for (int i = 0; i < h.size(); ++i) {
    h.at(i)->SetMarkerStyle(markerstyle.at(i));
    h.at(i)->SetMarkerColor(color.at(i));
    h.at(i)->SetLineColor(color.at(i));
    if (do_rebin) h.at(i)->Rebin(rebin_factor);
    if (do_normalize) h.at(i)->Scale(1.0/(float)h.at(i)->GetBinWidth(1));
  }
  if (set_xrange) h.at(0)->GetXaxis()->SetRangeUser(xlow, xhigh);
  if (set_logx) can->SetLogx();
  if (set_yrange) h.at(0)->GetYaxis()->SetRangeUser(ylow, yhigh);
  if (set_logy) can->SetLogy();
  if (set_title) {
    h.at(0)->GetXaxis()->SetTitle(xtitle.c_str());
    h.at(0)->GetYaxis()->SetTitle(ytitle.c_str());
  }
  h.at(0)->GetXaxis()->SetTitleOffset(0.98);
  h.at(0)->GetYaxis()->SetTitleOffset(1.15);
  h.at(0)->GetXaxis()->SetLabelSize(0.045);
  h.at(0)->GetYaxis()->SetLabelSize(0.045);
  h.at(0)->GetXaxis()->CenterTitle();
  h.at(0)->GetYaxis()->CenterTitle();

  h.at(0)->Draw("p");
  for (int i = 0; i < h.size(); ++i) {
    if (i == 0) continue;
    h.at(i)->Draw("p same");
  }

  if (do_data) myText(0.065, 0.97, 1, "#bf{#it{sPHENIX}} Data", 0.04);
  else myText(0.065, 0.97, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
  myText(0.065, 0.92, 1, run_status.c_str(), 0.04);
  if (set_text) {
    for (int i = 0; i < text.size(); i++) {
      myText(xstart, ystart-i*0.05, 1, text[i].c_str(), size);
    }
  }
  if (set_legend) {
    for (int i = 0; i < legend.size(); i++) {
      myMarkerLineText(xstart_legend, ystart_legend-i*0.05, 1, color.at(i), markerstyle.at(i), color.at(i), 1, legend.at(i).c_str(), size_legend, true);
    }
  }

  can->SaveAs(output_name.c_str());
  delete can;
}

void draw() {
  SetAtlasStyle();
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TFile *fin_Jet10GeV = new TFile("output_Jet10GeV_backup3.root", "READ");
  //TFile *fin_Jet10GeV = new TFile("output_unfolded.root", "READ");
  if (!fin_Jet10GeV) {
    std::cout << "Error: cannot open input file" << std::endl;
    return;
  }

  TH1F *h_truth_test = (TH1F*)fin_Jet10GeV->Get("h_truth_test");

  TH1F *h_measure_reco_test = (TH1F*)fin_Jet10GeV->Get("h_measure_reco_test");
  TH1F *h_unfold_reco_half_2 = (TH1F*)fin_Jet10GeV->Get("h_unfold_reco_half_2");

  TH1F *h_measure_calib_test = (TH1F*)fin_Jet10GeV->Get("h_measure_calib_test");
  TH1F *h_unfold_calib_half_2 = (TH1F*)fin_Jet10GeV->Get("h_unfold_calib_half_2");

  TH2F *h_respmatrix_reco_all = (TH2F*)fin_Jet10GeV->Get("h_respmatrix_reco_all");
  TH2F *h_respmatrix_calib_all = (TH2F*)fin_Jet10GeV->Get("h_respmatrix_calib_all");

  std::vector<TH1F*> h_input;
  std::vector<int> color;
  std::vector<int> markerstyle;
  std::vector<std::string> text;
  std::vector<std::string> legend;

  h_input.push_back((TH1F*)h_truth_test);
  h_input.push_back((TH1F*)h_measure_reco_test);
  h_input.push_back((TH1F*)h_unfold_reco_half_2);
  color.push_back(kRed);
  color.push_back(kBlack);
  color.push_back(kGreen);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("|z_{vertex}| < 30 cm");
  text.push_back("#eta^{jet} within detector range");
  text.push_back("E^{jet} > 0 GeV");
  legend.push_back("Truth jet spectrum");
  legend.push_back("Reco jet spectrum");
  legend.push_back("Unfolded jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              false, 0, 40, false,
                              false, 0, 0.5, true,
                              true, 0, 2,
                              true, "p_{T}^{jet} [GeV]", "dN_{jet}/dp_{T}^{jet}", "Ratio",
                              false, "Run 21 Jet10GeV simulation data set",
                              true, text, 0.55, 0.81, 0.04,
                              true, legend, 0.58, 0.66, 0.04,
                              "figure/h_unfold_reco_half_2.png");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_truth_test);
  h_input.push_back((TH1F*)h_measure_calib_test);
  h_input.push_back((TH1F*)h_unfold_calib_half_2);
  color.push_back(kRed);
  color.push_back(kBlack);
  color.push_back(kGreen);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("|z_{vertex}| < 30 cm");
  text.push_back("#eta^{jet} within detector range");
  text.push_back("E^{jet} > 0 GeV");
  legend.push_back("Truth jet spectrum");
  legend.push_back("Calib jet spectrum");
  legend.push_back("Unfolded jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              false, 0, 40, false,
                              false, 0, 0.5, true,
                              true, 0, 2,
                              true, "p_{T}^{jet} [GeV]", "dN_{jet}/dp_{T}^{jet}", "Ratio",
                              false, "Run 21 Jet10GeV simulation data set",
                              true, text, 0.55, 0.81, 0.04,
                              true, legend, 0.58, 0.66, 0.04,
                              "figure/h_unfold_calib_half_2.png");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();
}