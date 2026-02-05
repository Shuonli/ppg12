#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "plotcommon.h"

// Helper functions
void styleDataHist(TH1D *h) {
    h->SetLineColor(kBlack);
    h->SetMarkerColor(kBlack);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.0);
    h->SetStats(0);
    h->SetTitle("");
}

void styleBkgHist(TH1D *h, Color_t color = kRed, bool fill = false) {
    h->SetLineColor(color);
    h->SetLineWidth(2);
    h->SetStats(0);
    if (fill) {
        h->SetFillColor(color);
        h->SetFillStyle(3354);
    }
}

void setupTimeAxis(TH1 *h, double ymax, double ymin = 0.0, const char *ytitle = "Counts") {
    h->GetXaxis()->SetTitle("Cluster-MBD Time [ns]");
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetTitle(ytitle);
    h->GetYaxis()->SetTitleOffset(1.5);
    h->GetYaxis()->SetRangeUser(ymin, ymax);
}

void drawLabels(double x, double y, int ipt, const std::vector<double> &pT_bin_edges,
                double size = 0.04, double spacing = 0.05) {
    myText(x, y, 1, strleg1.c_str(), size);
    myText(x, y - spacing, 1, strleg2.c_str(), size);
    myText(x, y - 2*spacing, 1, strleg3.c_str(), size);
    myText(x, y - 3*spacing, 1, Form("%.0f < E_{T} < %.0f GeV", pT_bin_edges[ipt], pT_bin_edges[ipt + 1]), size);
}

TLegend* makeLegend(double x1 = 0.55, double y1 = 0.75, double x2 = 0.88, double y2 = 0.88) {
    TLegend *leg = new TLegend(x1, y1, x2, y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    return leg;
}

void drawTailCutLine(double xcut, double ymax) {
    TLine *lcut = new TLine(xcut, 0.0, xcut, ymax);
    lcut->SetLineStyle(2);
    lcut->SetLineColor(kGray + 2);
    lcut->Draw("SAME");
}

void plot_npb_time_bkgsub()
{
    init_plot();

    const std::string savePath = "../PPG12-analysis-note/Figures/npb_time_bkgsub/";
    gSystem->Exec(Form("mkdir -p %s", savePath.c_str()));

    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node config = YAML::LoadFile("../efficiencytool/config_showershape.yaml");

    std::vector<double> pT_bin_edges = config["analysis"]["pT_bins"].as<std::vector<double>>();
    const int nPtBins = pT_bin_edges.size() - 1;
    const int nEtaBins = 1;

    const std::string dataFile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histoshower_shape_.root";
    TFile *f_data = TFile::Open(dataFile.c_str(), "READ");
    if (!f_data || f_data->IsZombie()) {
        std::cerr << "Error: Could not open data file!" << std::endl;
        return;
    }

    const double tail_time_cut = -7.0;
    const int rebin_time = 2;
    const double npb_bkg_max = 0.2;
    std::vector<double> npb_thresholds = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    for (int ieta = 0; ieta < nEtaBins; ++ieta)
    {
        for (int ipt = 0; ipt < nPtBins; ++ipt)
        {
            // ========== Original plots using h_npb_score_vs_time ==========
            TString hist_name = Form("h_npb_score_vs_time_eta%d_pt%d", ieta, ipt);
            TH2D *h2_time = dynamic_cast<TH2D *>(f_data->Get(hist_name));
            if (!h2_time) { std::cerr << "Warning: Missing " << hist_name << std::endl; continue; }

            auto getYBin = [&](double val) { return h2_time->GetYaxis()->FindBin(val); };
            
            TH1D *h_time_all = h2_time->ProjectionX(Form("%s_all", hist_name.Data()), 1, -1);
            TH1D *h_time_bkg = h2_time->ProjectionX(Form("%s_bkg", hist_name.Data()), 1, getYBin(0.2 - 1e-6));
            TH1D *h_time_bkg_high = h2_time->ProjectionX(Form("%s_high", hist_name.Data()), getYBin(0.5 + 1e-6), getYBin(1.0 - 1e-6));
            
            if (!h_time_all || !h_time_bkg || !h_time_bkg_high) continue;
            h_time_all->Rebin(rebin_time);
            h_time_bkg->Rebin(rebin_time);
            h_time_bkg_high->Rebin(rebin_time);

            // 2D correlation plot
            {
                TH2D *h2_norm = (TH2D *)h2_time->Clone(Form("%s_norm", hist_name.Data()));
                h2_norm->SetDirectory(nullptr);
                h2_norm->RebinY(4);
                for (int iy = 1; iy <= h2_norm->GetNbinsY(); ++iy) {
                    double row_sum = 0.0;
                    for (int ix = 1; ix <= h2_norm->GetNbinsX(); ++ix) row_sum += h2_norm->GetBinContent(ix, iy);
                    if (row_sum > 0.0)
                        for (int ix = 1; ix <= h2_norm->GetNbinsX(); ++ix)
                            h2_norm->SetBinContent(ix, iy, h2_norm->GetBinContent(ix, iy) / row_sum);
                }
                TCanvas *c = new TCanvas(Form("c_2d_eta%d_pt%d", ieta, ipt), "", 700, 600);
                c->SetRightMargin(0.15);
                h2_norm->GetXaxis()->SetRangeUser(-15, 10);
                h2_norm->SetStats(0); h2_norm->SetTitle("");
                h2_norm->GetXaxis()->SetTitle("Cluster-MBD Time [ns]"); h2_norm->GetXaxis()->SetNdivisions(505);
                h2_norm->GetYaxis()->SetTitle("NPB Score"); h2_norm->GetYaxis()->SetTitleOffset(1.2);
                h2_norm->GetZaxis()->SetTitleOffset(1.2); h2_norm->GetZaxis()->SetRangeUser(0.02, 0.2);
                h2_norm->Draw("COLZ");
                drawLabels(0.18, 0.90, ipt, pT_bin_edges, 0.035);
                myText(0.18, 0.70, 1, "Normalized per NPB bin", 0.035);
                c->SaveAs(Form("%s/npb_time_2d_eta%d_pt%d.pdf", savePath.c_str(), ieta, ipt));
                delete h2_norm; delete c;
            }

            int bin_tail_max = std::min(h_time_all->GetNbinsX(), h_time_all->GetXaxis()->FindBin(tail_time_cut - 1e-6));
            double data_tail = h_time_all->Integral(1, bin_tail_max);
            double bkg_tail = h_time_bkg->Integral(1, bin_tail_max);
            double bkg_high_tail = h_time_bkg_high->Integral(1, bin_tail_max);

            if (bkg_tail <= 0.0) { std::cerr << "Warning: zero bkg tail for " << hist_name << std::endl; continue; }

            TH1D *h_bkg_scaled = (TH1D *)h_time_bkg->Clone(Form("%s_scaled", hist_name.Data()));
            h_bkg_scaled->Scale(data_tail / bkg_tail);
            TH1D *h_bkg_low_scaled = (TH1D *)h_time_bkg->Clone(Form("%s_low_scaled", hist_name.Data()));
            h_bkg_low_scaled->Scale(bkg_high_tail / bkg_tail);

            // Overlay plots
            auto makeOverlay = [&](TH1D *hData, TH1D *hBkg, const char *suffix, const char *legData, const char *legBkg, const char *extra = "") {
                styleDataHist(hData); styleBkgHist(hBkg);
                double ymax = std::max(hData->GetMaximum(), hBkg->GetMaximum()) * 1.3;
                TCanvas *c = new TCanvas(Form("c_%s_eta%d_pt%d", suffix, ieta, ipt), "", 600, 600);
                setupTimeAxis(hData, ymax);
                hData->Draw("E"); hBkg->Draw("HIST SAME");
                drawTailCutLine(tail_time_cut, ymax);
                drawLabels(0.20, 0.90, ipt, pT_bin_edges);
                if (strlen(extra) > 0) myText(0.20, 0.70, 1, extra, 0.04);
                TLegend *leg = makeLegend();
                leg->AddEntry(hData, legData, "lep"); leg->AddEntry(hBkg, legBkg, "l");
                leg->Draw();
                c->SaveAs(Form("%s/%s_eta%d_pt%d.pdf", savePath.c_str(), suffix, ieta, ipt));
                delete leg; delete c;
            };

            makeOverlay(h_time_all, h_bkg_scaled, "npb_time_bkg_overlay", "All data", "NPB 0-0.2",
                        Form("NPB 0-0.2 scaled to t < %.0f ns", tail_time_cut));
            makeOverlay(h_time_bkg_high, h_bkg_low_scaled, "npb_time_low_vs_high", "NPB 0.5-1.0", "NPB 0-0.2 (scaled)",
                        Form("Scaled with t < %.0f ns", tail_time_cut));

            // Subtraction plots
            auto makeSubPlot = [&](TH1D *hData, TH1D *hBkg, const char *suffix, const char *ytitle, const char *extra) {
                TH1D *h_sub = (TH1D *)hData->Clone(Form("%s_%s", hist_name.Data(), suffix));
                h_sub->SetDirectory(nullptr);
                h_sub->Add(hBkg, -1.0);
                double ymin = (h_sub->GetMinimum() < 0) ? 1.1 * h_sub->GetMinimum() : 0.0;
                double ymax = (h_sub->GetMaximum() > 0) ? 1.3 * h_sub->GetMaximum() : 1.0;
                TCanvas *c = new TCanvas(Form("c_%s_eta%d_pt%d", suffix, ieta, ipt), "", 600, 600);
                styleDataHist(h_sub);
                setupTimeAxis(h_sub, ymax, ymin, ytitle);
                h_sub->Draw("E");
                TLine l0(h_sub->GetXaxis()->GetXmin(), 0.0, h_sub->GetXaxis()->GetXmax(), 0.0);
                l0.SetLineStyle(2); l0.SetLineColor(kGray + 2); l0.Draw("SAME");
                drawLabels(0.20, 0.90, ipt, pT_bin_edges);
                myText(0.20, 0.70, 1, extra, 0.04);
                c->SaveAs(Form("%s/%s_eta%d_pt%d.pdf", savePath.c_str(), suffix, ieta, ipt));
                delete h_sub; delete c;
            };

            makeSubPlot(h_time_all, h_bkg_scaled, "npb_time_bkg_sub", "Counts (Data - scaled npb)", "Data - scaled NPB 0-0.2");
            makeSubPlot(h_time_bkg_high, h_bkg_low_scaled, "npb_time_high_minus_low", "Counts (NPB 0.5-1.0 - scaled 0-0.2)", "NPB 0.5-1.0 minus scaled 0-0.2");

            delete h_bkg_scaled; delete h_bkg_low_scaled;
            delete h_time_all; delete h_time_bkg; delete h_time_bkg_high;

            // ========== Purity analysis using h_npb_score_vs_time_clean ==========
            TString hist_name_clean = Form("h_npb_score_vs_time_clean_eta%d_pt%d", ieta, ipt);
            TH2D *h2_clean = dynamic_cast<TH2D *>(f_data->Get(hist_name_clean));
            if (!h2_clean) { std::cerr << "Warning: Missing " << hist_name_clean << std::endl; continue; }

            TH1D *h_bkg_template = h2_clean->ProjectionX(Form("%s_bkg", hist_name_clean.Data()), 1, h2_clean->GetYaxis()->FindBin(npb_bkg_max - 1e-6));
            h_bkg_template->Rebin(rebin_time);
            
            int bin_tail_max_clean = std::min(h_bkg_template->GetNbinsX(), h_bkg_template->GetXaxis()->FindBin(tail_time_cut - 1e-6));
            double bkg_template_tail = h_bkg_template->Integral(1, bin_tail_max_clean);
            if (bkg_template_tail <= 0.0) { delete h_bkg_template; continue; }

            std::vector<double> thresh_vals, purity_vals, purity_errs, total_vals, bkg_vals, sig_vals;

            for (double npb_thresh : npb_thresholds) {
                TH1D *h_sel = h2_clean->ProjectionX(Form("%s_sel%.0f", hist_name_clean.Data(), npb_thresh*10),
                    h2_clean->GetYaxis()->FindBin(npb_thresh + 1e-6), h2_clean->GetYaxis()->FindBin(1.0 - 1e-6));
                h_sel->Rebin(rebin_time);

                double sel_tail = h_sel->Integral(1, bin_tail_max_clean);
                double scale = sel_tail / bkg_template_tail;
                double total = h_sel->Integral();
                double bkg = h_bkg_template->Integral() * scale;
                double sig = std::max(0.0, total - bkg);
                double purity = (total > 0) ? sig / total : 0.0;
                double purity_err = (total > 0) ? std::sqrt(purity * (1 - purity) / total) : 0.0;

                thresh_vals.push_back(npb_thresh);
                purity_vals.push_back(purity);
                purity_errs.push_back(purity_err);
                total_vals.push_back(total);
                bkg_vals.push_back(bkg);
                sig_vals.push_back(sig);

                std::cout << Form("eta%d pt%d NPB>%.1f: total=%.0f, bkg=%.0f, sig=%.0f, purity=%.3f", ieta, ipt, npb_thresh, total, bkg, sig, purity) << std::endl;

                // Overlay plot for this threshold
                TH1D *h_bkg_sc = (TH1D *)h_bkg_template->Clone(Form("%s_bkg_sc%.0f", hist_name_clean.Data(), npb_thresh*10));
                h_bkg_sc->Scale(scale);
                styleDataHist(h_sel); styleBkgHist(h_bkg_sc, kRed, true);
                double ymax = std::max(h_sel->GetMaximum(), h_bkg_sc->GetMaximum()) * 1.4;
                TCanvas *c = new TCanvas(Form("c_pur_ov_eta%d_pt%d_npb%.0f", ieta, ipt, npb_thresh*10), "", 600, 600);
                setupTimeAxis(h_sel, ymax);
                h_sel->Draw("E"); h_bkg_sc->Draw("HIST SAME");
                drawTailCutLine(tail_time_cut, ymax);
                drawLabels(0.20, 0.90, ipt, pT_bin_edges, 0.035, 0.04);
                myText(0.20, 0.74, 1, Form("NPB > %.1f (clean)", npb_thresh), 0.035);
                myText(0.20, 0.70, 1, Form("Purity = %.1f%%", purity * 100.0), 0.035);
                TLegend *leg = makeLegend();
                leg->AddEntry(h_sel, Form("NPB %.1f-1.0", npb_thresh), "lep");
                leg->AddEntry(h_bkg_sc, "Bkg (NPB<0.2)", "lf");
                leg->Draw();
                c->SaveAs(Form("%s/purity_overlay_eta%d_pt%d_npb%d.pdf", savePath.c_str(), ieta, ipt, (int)(npb_thresh*10)));
                delete leg; delete h_bkg_sc; delete h_sel; delete c;
            }

            // Purity vs threshold plot
            {
                TCanvas *c = new TCanvas(Form("c_purity_eta%d_pt%d", ieta, ipt), "", 600, 600);
                TGraphErrors *g = new TGraphErrors(thresh_vals.size());
                for (size_t i = 0; i < thresh_vals.size(); ++i) { g->SetPoint(i, thresh_vals[i], purity_vals[i]); g->SetPointError(i, 0, purity_errs[i]); }
                g->SetMarkerStyle(20); g->SetMarkerSize(1.2); g->SetMarkerColor(kBlue+1); g->SetLineColor(kBlue+1); g->SetLineWidth(2);
                TH1F *fr = new TH1F(Form("fr_pur_%d_%d", ieta, ipt), "", 10, 0.05, 0.95);
                fr->SetStats(0); fr->GetXaxis()->SetTitle("NPB Score Threshold"); fr->GetYaxis()->SetTitle("Purity");
                fr->GetYaxis()->SetTitleOffset(1.4); fr->GetYaxis()->SetRangeUser(0.8, 1.1); fr->Draw();
                g->Draw("PL SAME");
                drawLabels(0.20, 0.90, ipt, pT_bin_edges);
                myText(0.20, 0.70, 1, "Bkg template: NPB < 0.2", 0.04);
                myText(0.20, 0.65, 1, Form("Tail fit: t < %.0f ns", tail_time_cut), 0.04);
                c->SaveAs(Form("%s/purity_vs_npb_eta%d_pt%d.pdf", savePath.c_str(), ieta, ipt));
                delete g; delete fr; delete c;
            }

            // Yield vs threshold plot
            {
                TCanvas *c = new TCanvas(Form("c_yield_eta%d_pt%d", ieta, ipt), "", 600, 600);
                TGraph *gT = new TGraph(thresh_vals.size()), *gB = new TGraph(thresh_vals.size()), *gS = new TGraph(thresh_vals.size());
                for (size_t i = 0; i < thresh_vals.size(); ++i) { gT->SetPoint(i, thresh_vals[i], total_vals[i]); gB->SetPoint(i, thresh_vals[i], bkg_vals[i]); gS->SetPoint(i, thresh_vals[i], sig_vals[i]); }
                auto styleG = [](TGraph *g, int m, Color_t col) { g->SetMarkerStyle(m); g->SetMarkerSize(1.2); g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(2); };
                styleG(gT, 20, kBlack); styleG(gB, 21, kRed); styleG(gS, 22, kBlue+1);
                double maxY = *std::max_element(total_vals.begin(), total_vals.end()) * 1.3;
                TH1F *fr = new TH1F(Form("fr_yld_%d_%d", ieta, ipt), "", 10, 0.15, 0.95);
                fr->SetStats(0); fr->GetXaxis()->SetTitle("NPB Score Threshold"); fr->GetYaxis()->SetTitle("Yield");
                fr->GetYaxis()->SetTitleOffset(1.5); fr->GetYaxis()->SetRangeUser(0, maxY); fr->Draw();
                gT->Draw("PL SAME"); gB->Draw("PL SAME"); gS->Draw("PL SAME");
                drawLabels(0.20, 0.90, ipt, pT_bin_edges, 0.035);
                TLegend *leg = makeLegend(0.55, 0.72, 0.88, 0.88);
                leg->AddEntry(gT, "Total", "pl"); leg->AddEntry(gB, "Background", "pl"); leg->AddEntry(gS, "Signal", "pl");
                leg->Draw();
                c->SaveAs(Form("%s/yield_vs_npb_eta%d_pt%d.pdf", savePath.c_str(), ieta, ipt));
                delete gT; delete gB; delete gS; delete fr; delete leg; delete c;
            }

            delete h_bkg_template;
        }
    }

    f_data->Close();
    std::cout << "Done. Output saved to: " << savePath << std::endl;
}
