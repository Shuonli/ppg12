// PlotROC.C
// Reads merged roc_isoET ROOT file and produces ROC curves per cluster ET bin
// comparing 7 isolation energy definitions (iso70, iso60, iso120, isotopo03, isotopo04, isotoposoft03, isotoposoft04).
// Also produces signal-comparison plots (direct vs fragmentation) per ET bin.
//
// Usage: root -l -b -q 'PlotROC.C("results/roc_isoET_merged_isoroc.root","roc_plots")'

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TLine.h>

#include "../plotting/plotcommon.h"

// Trapezoidal AUC for a TGraph (x = bg eff, y = sig eff).
double ComputeAUC(TGraph *g)
{
    if (!g || g->GetN() < 2) return 0.;
    double *x = g->GetX();
    double *y = g->GetY();
    double auc = 0.;
    for (int i = 0; i < g->GetN() - 1; i++)
        auc += 0.5 * (y[i] + y[i + 1]) * std::abs(x[i + 1] - x[i]);
    return auc;
}

// Build ROC TGraph from signal and background TH1D.
// Scans the isoET cut from low to high (inclusive cut: isoET < cut).
// Returns TGraph of (bg_efficiency, sig_efficiency).
TGraph *BuildROC(TH1D *h_sig, TH1D *h_bg)
{
    if (!h_sig || !h_bg) return nullptr;
    double total_sig = h_sig->Integral(0, h_sig->GetNbinsX() + 1);
    double total_bg  = h_bg ->Integral(0, h_bg ->GetNbinsX() + 1);
    if (total_sig <= 0 || total_bg <= 0) return nullptr;

    int nbins = h_sig->GetNbinsX();
    std::vector<double> bg_eff, sig_eff;

    for (int j = 0; j <= nbins + 1; j++)
    {
        double cs = h_sig->Integral(0, j);
        double cb = h_bg ->Integral(0, j);
        sig_eff.push_back(cs / total_sig);
        bg_eff.push_back(cb / total_bg);
    }

    TGraph *g = new TGraph((int)sig_eff.size(), bg_eff.data(), sig_eff.data());
    return g;
}

void PlotROC(const std::string &infile  = "results/roc_isoET_merged_isoroc.root",
             const std::string &outdir  = "roc_plots")
{
    gSystem->mkdir(outdir.c_str(), kTRUE);

    TFile *fin = TFile::Open(infile.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << infile << std::endl;
        return;
    }

    // Read pT bin info from stored histogram
    TH1D *h_ptbins = dynamic_cast<TH1D*>(fin->Get("h_pT_bins"));
    if (!h_ptbins)
    {
        std::cerr << "ERROR: h_pT_bins not found in " << infile << std::endl;
        return;
    }
    int n_pT_bins = h_ptbins->GetNbinsX();
    std::vector<float> pT_bins(n_pT_bins + 1);
    for (int i = 0; i <= n_pT_bins; i++) pT_bins[i] = h_ptbins->GetXaxis()->GetBinLowEdge(i + 1);

    const char *iso_labels[7]  = {"iso70",  "iso60",  "iso120",  "isotopo03",    "isotopo04",    "isotoposoft03",    "isotoposoft04"};
    const char *iso_titles[7]  = {"70 MeV", "60 MeV", "120 MeV", "topo R=0.3",  "topo R=0.4",  "topo soft R=0.3", "topo soft R=0.4"};
    const int   iso_colors[7]  = {kBlue, kRed, kGreen + 2, kOrange + 1, kViolet + 1, kOrange + 1, kViolet + 1};
    const int   iso_styles[7]  = {1,     2,    3,          1,           1,            2,            2};

    SetsPhenixStyle();

    const char *strInternal = "#bf{#it{sPHENIX}} Internal";
    const char *strPP       = "#it{p}+#it{p} #kern[-0.1]{#sqrt{#it{s}} = 200 GeV}";
    const char *strEta      = "|#it{#eta}^{cluster}| < 0.7";
    const float tsize       = 0.036;
    const float tsize_leg   = 0.032;

    // -----------------------------------------------------------------------
    // 1) Per-ET-bin ROC: compare 5 iso definitions for combined signal + AUC
    // -----------------------------------------------------------------------
    for (int ipt = 0; ipt < n_pT_bins; ipt++)
    {
        float pt_lo = pT_bins[ipt];
        float pt_hi = pT_bins[ipt + 1];

        TCanvas *c = new TCanvas(Form("c_roc_pt%d", ipt),
                                 Form("ROC %.0f-%.0f GeV", pt_lo, pt_hi), 600, 600);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.05);
        c->SetTopMargin(0.07);
        c->SetBottomMargin(0.14);
        c->SetGrid();

        // Diagonal reference line
        TGraph *gdiag = new TGraph(2);
        gdiag->SetPoint(0, 0, 0);
        gdiag->SetPoint(1, 1, 1);
        gdiag->SetLineStyle(kDashed);
        gdiag->SetLineColor(kGray + 1);
        gdiag->SetLineWidth(1);

        // Legend: lower-right
        TLegend *leg = new TLegend(0.38, 0.18, 0.93, 0.58);
        legStyle(leg, 0.12, tsize_leg, 42,
                 Form("%.0f < #it{E}_{T}^{cluster} < %.0f GeV", pt_lo, pt_hi));

        bool first = true;
        for (int iiso = 0; iiso < 7; iiso++)
        {
            TH1D *h_sig = dynamic_cast<TH1D*>(fin->Get(Form("h_signal_%s_pt%d", iso_labels[iiso], ipt)));
            TH1D *h_bg  = dynamic_cast<TH1D*>(fin->Get(Form("h_bg_%s_pt%d",     iso_labels[iiso], ipt)));
            TGraph *g = BuildROC(h_sig, h_bg);
            if (!g) continue;

            double auc = ComputeAUC(g);

            g->SetLineColor(iso_colors[iiso]);
            g->SetLineStyle(iso_styles[iiso]);
            g->SetLineWidth(2);
            g->SetMarkerStyle(0);

            if (first)
            {
                g->GetXaxis()->SetTitle("Background efficiency");
                g->GetYaxis()->SetTitle("Signal efficiency");
                g->GetXaxis()->SetTitleSize(tsize);
                g->GetYaxis()->SetTitleSize(tsize);
                g->GetXaxis()->SetLabelSize(tsize);
                g->GetYaxis()->SetLabelSize(tsize);
                g->GetXaxis()->SetTitleOffset(1.1);
                g->GetYaxis()->SetTitleOffset(1.3);
                g->GetXaxis()->SetLimits(0, 1);
                g->GetHistogram()->SetMinimum(0);
                g->GetHistogram()->SetMaximum(1.05);
                g->Draw("AL");
                gdiag->Draw("L same");
                first = false;
            }
            else
            {
                g->Draw("L same");
            }
            leg->AddEntry(g, Form("%s  (AUC = %.3f)", iso_titles[iiso], auc), "L");
        }
        leg->Draw();

        myText(0.93, 0.93, 1, strInternal,                    tsize, true);
        myText(0.93, 0.88, 1, strPP,                          tsize, true);
        myText(0.93, 0.83, 1, strEta,                         tsize, true);
        myText(0.93, 0.78, 1, "Direct + frag #gamma signal",  tsize, true);

        c->SaveAs(Form("%s/roc_combined_pt%d_%.0f_%.0f.pdf",  outdir.c_str(), ipt, pt_lo, pt_hi));
        c->SaveAs(Form("%s/roc_combined_pt%d_%.0f_%.0f.root", outdir.c_str(), ipt, pt_lo, pt_hi));
        delete c;
        delete gdiag;
        delete leg;
    }

    // -----------------------------------------------------------------------
    // 2) Per-ET-bin: direct vs frag signal comparison (one canvas per iso def)
    // -----------------------------------------------------------------------
    const char *sig_labels[2]  = {"direct", "frag"};
    const char *sig_titles[2]  = {"Direct #gamma signal", "Frag #gamma signal"};
    const int   sig_colors[2]  = {kBlue, kRed};
    const int   sig_styles[2]  = {1, 2};

    for (int iiso = 0; iiso < 7; iiso++)
    {
        for (int ipt = 0; ipt < n_pT_bins; ipt++)
        {
            float pt_lo = pT_bins[ipt];
            float pt_hi = pT_bins[ipt + 1];

            TH1D *h_bg = dynamic_cast<TH1D*>(fin->Get(Form("h_bg_%s_pt%d", iso_labels[iiso], ipt)));
            if (!h_bg) continue;

            TCanvas *c = new TCanvas(Form("c_roc_%s_sigcomp_pt%d", iso_labels[iiso], ipt),
                                     Form("ROC %s %.0f-%.0f GeV", iso_titles[iiso], pt_lo, pt_hi), 600, 600);
            c->SetLeftMargin(0.15);
            c->SetRightMargin(0.05);
            c->SetTopMargin(0.07);
            c->SetBottomMargin(0.14);
            c->SetGrid();

            TLegend *leg = new TLegend(0.38, 0.10, 0.93, 0.42);
            legStyle(leg, 0.12, tsize_leg, 42,
                     Form("%s, %.0f < #it{E}_{T} < %.0f GeV", iso_titles[iiso], pt_lo, pt_hi));

            TGraph *gdiag = new TGraph(2);
            gdiag->SetPoint(0, 0, 0);
            gdiag->SetPoint(1, 1, 1);
            gdiag->SetLineStyle(kDashed);
            gdiag->SetLineColor(kGray + 1);
            gdiag->SetLineWidth(1);

            bool first = true;
            for (int isig = 0; isig < 2; isig++)
            {
                TH1D *h_sig = dynamic_cast<TH1D*>(fin->Get(
                    Form("h_%s_%s_pt%d", sig_labels[isig], iso_labels[iiso], ipt)));
                TGraph *g = BuildROC(h_sig, h_bg);
                if (!g) continue;

                double auc = ComputeAUC(g);
                g->SetLineColor(sig_colors[isig]);
                g->SetLineStyle(sig_styles[isig]);
                g->SetLineWidth(2);
                g->SetMarkerStyle(0);

                if (first)
                {
                    g->GetXaxis()->SetTitle("Background efficiency");
                    g->GetYaxis()->SetTitle("Signal efficiency");
                    g->GetXaxis()->SetTitleSize(tsize);
                    g->GetYaxis()->SetTitleSize(tsize);
                    g->GetXaxis()->SetLabelSize(tsize);
                    g->GetYaxis()->SetLabelSize(tsize);
                    g->GetXaxis()->SetTitleOffset(1.1);
                    g->GetYaxis()->SetTitleOffset(1.3);
                    g->GetXaxis()->SetLimits(0, 1);
                    g->GetHistogram()->SetMinimum(0);
                    g->GetHistogram()->SetMaximum(1.05);
                    g->Draw("AL");
                    gdiag->Draw("L same");
                    first = false;
                }
                else
                {
                    g->Draw("L same");
                }
                leg->AddEntry(g, Form("%s  (AUC = %.3f)", sig_titles[isig], auc), "L");
            }
            leg->Draw();

            myText(0.93, 0.93, 1, strInternal, tsize, true);
            myText(0.93, 0.88, 1, strPP,       tsize, true);
            myText(0.93, 0.83, 1, strEta,      tsize, true);

            c->SaveAs(Form("%s/roc_%s_sigcomp_pt%d_%.0f_%.0f.pdf",
                           outdir.c_str(), iso_labels[iiso], ipt, pt_lo, pt_hi));
            delete c;
            delete gdiag;
            delete leg;
        }
    }

    // -----------------------------------------------------------------------
    // 3) Summary: overlay ET bins for each iso definition (combined signal)
    // -----------------------------------------------------------------------
    const int pt_colors[] = {kBlack, kBlue, kRed, kGreen+2, kMagenta, kOrange+1,
                              kCyan+1, kViolet, kTeal, kPink+1};
    const int pt_ncolors = sizeof(pt_colors) / sizeof(pt_colors[0]);

    for (int iiso = 0; iiso < 7; iiso++)
    {
        TCanvas *c = new TCanvas(Form("c_roc_summary_%s", iso_labels[iiso]),
                                 Form("ROC summary %s", iso_titles[iiso]), 600, 600);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.05);
        c->SetTopMargin(0.07);
        c->SetBottomMargin(0.14);
        c->SetGrid();

        TLegend *leg = new TLegend(0.38, 0.10, 0.93, 0.50);
        legStyle(leg, 0.12, tsize_leg, 42,
                 Form("Isolation: %s", iso_titles[iiso]));

        TGraph *gdiag = new TGraph(2);
        gdiag->SetPoint(0, 0, 0);
        gdiag->SetPoint(1, 1, 1);
        gdiag->SetLineStyle(kDashed);
        gdiag->SetLineColor(kGray + 1);
        gdiag->SetLineWidth(1);

        bool first = true;
        for (int ipt = 0; ipt < n_pT_bins; ipt++)
        {
            TH1D *h_sig = dynamic_cast<TH1D*>(fin->Get(Form("h_signal_%s_pt%d", iso_labels[iiso], ipt)));
            TH1D *h_bg  = dynamic_cast<TH1D*>(fin->Get(Form("h_bg_%s_pt%d",     iso_labels[iiso], ipt)));
            TGraph *g = BuildROC(h_sig, h_bg);
            if (!g) continue;

            int col = pt_colors[ipt % pt_ncolors];
            g->SetLineColor(col);
            g->SetLineWidth(2);
            g->SetMarkerStyle(0);

            float pt_lo = pT_bins[ipt];
            float pt_hi = pT_bins[ipt + 1];

            if (first)
            {
                g->GetXaxis()->SetTitle("Background efficiency");
                g->GetYaxis()->SetTitle("Signal efficiency");
                g->GetXaxis()->SetTitleSize(tsize);
                g->GetYaxis()->SetTitleSize(tsize);
                g->GetXaxis()->SetLabelSize(tsize);
                g->GetYaxis()->SetLabelSize(tsize);
                g->GetXaxis()->SetTitleOffset(1.1);
                g->GetYaxis()->SetTitleOffset(1.3);
                g->GetXaxis()->SetLimits(0, 1);
                g->GetHistogram()->SetMinimum(0);
                g->GetHistogram()->SetMaximum(1.05);
                g->Draw("AL");
                gdiag->Draw("L same");
                first = false;
            }
            else
            {
                g->Draw("L same");
            }
            leg->AddEntry(g, Form("%.0f#font[122]{-}%.0f GeV", pt_lo, pt_hi), "L");
        }
        leg->Draw();

        myText(0.93, 0.93, 1, strInternal, tsize, true);
        myText(0.93, 0.88, 1, strPP,       tsize, true);
        myText(0.93, 0.83, 1, strEta,      tsize, true);

        c->SaveAs(Form("%s/roc_summary_%s.pdf",  outdir.c_str(), iso_labels[iiso]));
        c->SaveAs(Form("%s/roc_summary_%s.root", outdir.c_str(), iso_labels[iiso]));
        delete c;
        delete gdiag;
        delete leg;
    }

    // -----------------------------------------------------------------------
    // 4) Summary: overlay 5 iso definitions for each ET bin (combined signal)
    //    (same as plot 1 but all ET bins on separate pages of a single PDF)
    // -----------------------------------------------------------------------
    {
        TCanvas *c = new TCanvas("c_roc_all", "ROC all ET bins", 600, 600);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.05);
        c->SetTopMargin(0.07);
        c->SetBottomMargin(0.14);
        c->Print(Form("%s/roc_all_etbins.pdf[", outdir.c_str())); // open PDF

        for (int ipt = 0; ipt < n_pT_bins; ipt++)
        {
            float pt_lo = pT_bins[ipt];
            float pt_hi = pT_bins[ipt + 1];

            c->cd();
            c->Clear();
            c->SetGrid();

            TLegend *leg = new TLegend(0.38, 0.18, 0.93, 0.58);
            legStyle(leg, 0.12, tsize_leg, 42,
                     Form("%.0f < #it{E}_{T}^{cluster} < %.0f GeV", pt_lo, pt_hi));

            TGraph *gdiag = new TGraph(2);
            gdiag->SetPoint(0, 0, 0);
            gdiag->SetPoint(1, 1, 1);
            gdiag->SetLineStyle(kDashed);
            gdiag->SetLineColor(kGray + 1);
            gdiag->SetLineWidth(1);

            bool first = true;
            for (int iiso = 0; iiso < 7; iiso++)
            {
                TH1D *h_sig = dynamic_cast<TH1D*>(fin->Get(Form("h_signal_%s_pt%d", iso_labels[iiso], ipt)));
                TH1D *h_bg  = dynamic_cast<TH1D*>(fin->Get(Form("h_bg_%s_pt%d",     iso_labels[iiso], ipt)));
                TGraph *g = BuildROC(h_sig, h_bg);
                if (!g) continue;

                double auc = ComputeAUC(g);
                g->SetLineColor(iso_colors[iiso]);
                g->SetLineStyle(iso_styles[iiso]);
                g->SetLineWidth(2);
                g->SetMarkerStyle(0);

                if (first)
                {
                    g->GetXaxis()->SetTitle("Background efficiency");
                    g->GetYaxis()->SetTitle("Signal efficiency");
                    g->GetXaxis()->SetTitleSize(tsize);
                    g->GetYaxis()->SetTitleSize(tsize);
                    g->GetXaxis()->SetLabelSize(tsize);
                    g->GetYaxis()->SetLabelSize(tsize);
                    g->GetXaxis()->SetTitleOffset(1.1);
                    g->GetYaxis()->SetTitleOffset(1.3);
                    g->GetXaxis()->SetLimits(0, 1);
                    g->GetHistogram()->SetMinimum(0);
                    g->GetHistogram()->SetMaximum(1.05);
                    g->Draw("AL");
                    gdiag->Draw("L same");
                    first = false;
                }
                else
                {
                    g->Draw("L same");
                }
                leg->AddEntry(g, Form("%s  (AUC = %.3f)", iso_titles[iiso], auc), "L");
            }
            leg->Draw();

            myText(0.93, 0.93, 1, strInternal,                   tsize, true);
            myText(0.93, 0.88, 1, strPP,                         tsize, true);
            myText(0.93, 0.83, 1, strEta,                        tsize, true);
            myText(0.93, 0.78, 1, "Direct + frag #gamma signal", tsize, true);

            c->Print(Form("%s/roc_all_etbins.pdf", outdir.c_str()));

            delete gdiag;
            delete leg;
        }
        c->Print(Form("%s/roc_all_etbins.pdf]", outdir.c_str())); // close PDF
        delete c;
    }

    fin->Close();
    std::cout << "ROC plots saved to: " << outdir << "/" << std::endl;
}
