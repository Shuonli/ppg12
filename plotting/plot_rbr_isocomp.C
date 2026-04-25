// plot_rbr_isocomp.C
// Run-by-run comparison plots for three isolation ET definitions
// (60 MeV tower, 120 MeV tower, topo R=0.4) in sPHENIX publication style.
//
// Usage: root -l -b -q 'plot_rbr_isocomp.C'

#include "plotcommon.h"
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TLine.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

struct IsoDef
{
    std::string file;
    std::string label;
    int         color;
    int         marker;
};

// Compute y-range from multiple TGraphErrors (ignoring zeros)
void getYRange(const std::vector<TGraphErrors *> &graphs,
               double &ymin, double &ymax)
{
    ymin = 1e30;
    ymax = -1e30;
    for (auto *g : graphs)
    {
        if (!g) continue;
        for (int i = 0; i < g->GetN(); i++)
        {
            double x, y;
            g->GetPoint(i, x, y);
            if (y <= 0) continue;
            if (y < ymin) ymin = y;
            if (y > ymax) ymax = y;
        }
    }
    if (ymin > ymax) { ymin = 0; ymax = 1; }
}

// Build ratio TGraphErrors from two TGraphErrors (C/A)
TGraphErrors *buildRatio(TGraphErrors *grNum, TGraphErrors *grDen)
{
    if (!grNum || !grDen) return nullptr;
    std::vector<double> x_r, y_r, ex_r, ey_r;
    for (int i = 0; i < std::min(grNum->GetN(), grDen->GetN()); i++)
    {
        double xN, yN, xD, yD;
        grNum->GetPoint(i, xN, yN);
        grDen->GetPoint(i, xD, yD);
        if (std::abs(xN - xD) > 1) continue;
        if (yD <= 0 || yN <= 0) continue;
        double eN = grNum->GetErrorY(i);
        double eD = grDen->GetErrorY(i);
        double ratio = yN / yD;
        double err = ratio * std::sqrt((eN / yN) * (eN / yN) + (eD / yD) * (eD / yD));
        x_r.push_back(xN);
        y_r.push_back(ratio);
        ex_r.push_back(0);
        ey_r.push_back(err);
    }
    if (x_r.empty()) return nullptr;
    return new TGraphErrors((int)x_r.size(), x_r.data(), y_r.data(),
                            ex_r.data(), ey_r.data());
}

void plot_rbr_isocomp()
{
    init_plot();
    gSystem->Exec("mkdir -p figures/rbr_isocomp");

    const std::string savePath = "figures/rbr_isocomp";

    // Crossing angle boundary
    const double xing_run = 51274;

    // Isolation definitions
    std::vector<IsoDef> defs = {
        {"../efficiencytool/rbrQA_60MeV.root",  "60 MeV (R=0.3)",  kRed,       20},
        {"../efficiencytool/rbrQA_120MeV.root", "120 MeV (R=0.3)", kGreen + 2, 21},
        {"../efficiencytool/rbrQA_topo.root",   "topo R=0.4",      kBlue,      22},
    };

    // Open files
    std::vector<TFile *> files;
    for (auto &d : defs)
    {
        TFile *f = TFile::Open(d.file.c_str(), "READ");
        if (!f || f->IsZombie())
        {
            std::cerr << "Cannot open " << d.file << std::endl;
            files.push_back(nullptr);
        }
        else
        {
            files.push_back(f);
        }
    }

    // Helper: get TGraphErrors from file
    auto getGr = [&](int idef, const char *name) -> TGraphErrors *
    {
        if (!files[idef]) return nullptr;
        return dynamic_cast<TGraphErrors *>(files[idef]->Get(name));
    };

    // Helper: draw multi-definition overlay
    auto drawOverlay = [&](const char *grName,
                           const char *ytitle,
                           const char *canvasName,
                           const char *extra_label = nullptr,
                           double forced_ymin = -1,
                           double forced_ymax = -1)
    {
        std::vector<TGraphErrors *> graphs;
        for (int i = 0; i < (int)defs.size(); i++)
            graphs.push_back(getGr(i, grName));

        double ymin, ymax;
        getYRange(graphs, ymin, ymax);
        double ypad = 0.15 * (ymax - ymin);
        if (forced_ymin >= 0) ymin = forced_ymin;
        else ymin = std::max(0.0, ymin - ypad);
        if (forced_ymax >= 0) ymax = forced_ymax;
        else ymax = ymax + ypad;

        // Get x-range
        double xmin = 1e30, xmax = -1e30;
        for (auto *g : graphs)
        {
            if (!g) continue;
            for (int i = 0; i < g->GetN(); i++)
            {
                double x, y;
                g->GetPoint(i, x, y);
                if (x < xmin) xmin = x;
                if (x > xmax) xmax = x;
            }
        }
        double xpad = 0.01 * (xmax - xmin);

        TCanvas *c = new TCanvas(canvasName, "", 1200, 500);
        c->SetLeftMargin(0.10);
        c->SetRightMargin(0.03);
        c->SetTopMargin(0.06);
        c->SetBottomMargin(0.14);

        TH1F *frame = new TH1F(Form("frame_%s", canvasName), "",
                                100, xmin - xpad, xmax + xpad);
        frame->GetYaxis()->SetRangeUser(ymin, ymax);
        frame->GetXaxis()->SetTitle("Run Number");
        frame->GetYaxis()->SetTitle(ytitle);
        frame->GetXaxis()->SetTitleSize(0.05);
        frame->GetYaxis()->SetTitleSize(0.05);
        frame->GetXaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetTitleOffset(0.8);
        frame->SetStats(0);
        frame->Draw("axis");

        // Crossing angle line
        TLine *xline = new TLine(xing_run, ymin, xing_run, ymax);
        xline->SetLineStyle(2);
        xline->SetLineColor(kGray + 1);
        xline->SetLineWidth(1);
        xline->Draw("same");

        float ly = 0.88;
        for (int i = 0; i < (int)defs.size(); i++)
        {
            if (!graphs[i]) continue;
            graphs[i]->SetMarkerStyle(defs[i].marker);
            graphs[i]->SetMarkerSize(0.4);
            graphs[i]->SetMarkerColor(defs[i].color);
            graphs[i]->SetLineColor(defs[i].color);
            graphs[i]->SetLineWidth(1);
            graphs[i]->Draw("p same");
            myMarkerLineText(0.62, ly, 0.4, defs[i].color, defs[i].marker,
                             defs[i].color, 1,
                             defs[i].label.c_str(), 0.035, true);
            ly -= 0.055;
        }

        myText(0.12, 0.92, 1, strleg1.c_str(), 0.04);
        myText(0.12, 0.87, 1, strleg2_1.c_str(), 0.035);
        myText(0.12, 0.82, 1, strleg3.c_str(), 0.035);
        if (extra_label)
            myText(0.12, 0.77, 1, extra_label, 0.035);

        // 0/1.5 mrad label
        myText(0.42, 0.92, kGray + 1, "0 mrad", 0.03);
        myText(0.56, 0.92, kGray + 1, "1.5 mrad", 0.03);

        c->SaveAs(Form("%s/%s.pdf", savePath.c_str(), canvasName));
        delete c;
    };

    // ---------------------------------------------------------------
    // Plot 1: Mean isoET vs run
    // ---------------------------------------------------------------
    drawOverlay("gr_avg_isoET",
                "#LT#it{E}_{T}^{iso}#GT [GeV]",
                "rbr_mean_isoET",
                "Pre-selected clusters");

    // ---------------------------------------------------------------
    // Plot 2: Region A (tight+iso) yield vs run
    // ---------------------------------------------------------------
    drawOverlay("gr_tight_iso",
                "Tight+Iso Clusters / Lumi [mb^{-1}]",
                "rbr_yield_A",
                "Region A (tight + iso)");

    // ---------------------------------------------------------------
    // Plot 3: Region C (nontight+iso) yield vs run
    // ---------------------------------------------------------------
    drawOverlay("gr_nontight_iso",
                "Nontight+Iso Clusters / Lumi [mb^{-1}]",
                "rbr_yield_C",
                "Region C (nontight + iso)");

    // ---------------------------------------------------------------
    // Plot 4: Region B (tight+noniso) yield vs run
    // ---------------------------------------------------------------
    drawOverlay("gr_tight_noniso",
                "Tight+Noniso Clusters / Lumi [mb^{-1}]",
                "rbr_yield_B",
                "Region B (tight + noniso)");

    // ---------------------------------------------------------------
    // Plot 5: Region D (nontight+noniso) yield vs run
    // ---------------------------------------------------------------
    drawOverlay("gr_nontight_noniso",
                "Nontight+Noniso Clusters / Lumi [mb^{-1}]",
                "rbr_yield_D",
                "Region D (nontight + noniso)");

    // ---------------------------------------------------------------
    // Plot 6: C/A ratio vs run
    // ---------------------------------------------------------------
    {
        TCanvas *c = new TCanvas("c_ratio_CA", "", 1200, 500);
        c->SetLeftMargin(0.10);
        c->SetRightMargin(0.03);
        c->SetTopMargin(0.06);
        c->SetBottomMargin(0.14);

        TH1F *frame = new TH1F("frame_ratio_CA", "",
                                100, 46500, 54500);
        frame->GetYaxis()->SetRangeUser(0, 4.5);
        frame->GetXaxis()->SetTitle("Run Number");
        frame->GetYaxis()->SetTitle("C / A  (nontight iso / tight iso)");
        frame->GetXaxis()->SetTitleSize(0.05);
        frame->GetYaxis()->SetTitleSize(0.05);
        frame->GetXaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetTitleOffset(0.8);
        frame->SetStats(0);
        frame->Draw("axis");

        TLine *xline = new TLine(xing_run, 0, xing_run, 4.5);
        xline->SetLineStyle(2);
        xline->SetLineColor(kGray + 1);
        xline->Draw("same");

        float ly = 0.88;
        for (int i = 0; i < (int)defs.size(); i++)
        {
            TGraphErrors *grA = getGr(i, "gr_tight_iso");
            TGraphErrors *grC = getGr(i, "gr_nontight_iso");
            TGraphErrors *grRatio = buildRatio(grC, grA);
            if (!grRatio) continue;

            grRatio->SetMarkerStyle(defs[i].marker);
            grRatio->SetMarkerSize(0.4);
            grRatio->SetMarkerColor(defs[i].color);
            grRatio->SetLineColor(defs[i].color);
            grRatio->SetLineWidth(1);
            grRatio->Draw("p same");

            myMarkerLineText(0.62, ly, 0.4, defs[i].color, defs[i].marker,
                             defs[i].color, 1,
                             defs[i].label.c_str(), 0.035, true);
            ly -= 0.055;
        }

        myText(0.12, 0.92, 1, strleg1.c_str(), 0.04);
        myText(0.12, 0.87, 1, strleg2_1.c_str(), 0.035);
        myText(0.12, 0.82, 1, strleg3.c_str(), 0.035);
        myText(0.42, 0.92, kGray + 1, "0 mrad", 0.03);
        myText(0.56, 0.92, kGray + 1, "1.5 mrad", 0.03);

        c->SaveAs(Form("%s/rbr_ratio_CA.pdf", savePath.c_str()));
        delete c;
    }

    // ---------------------------------------------------------------
    // Plot 7: D/B ratio vs run
    // ---------------------------------------------------------------
    {
        TCanvas *c = new TCanvas("c_ratio_DB", "", 1200, 500);
        c->SetLeftMargin(0.10);
        c->SetRightMargin(0.03);
        c->SetTopMargin(0.06);
        c->SetBottomMargin(0.14);

        TH1F *frame = new TH1F("frame_ratio_DB", "",
                                100, 46500, 54500);
        frame->GetYaxis()->SetRangeUser(0, 4.5);
        frame->GetXaxis()->SetTitle("Run Number");
        frame->GetYaxis()->SetTitle("D / B  (nontight noniso / tight noniso)");
        frame->GetXaxis()->SetTitleSize(0.05);
        frame->GetYaxis()->SetTitleSize(0.05);
        frame->GetXaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetTitleOffset(0.8);
        frame->SetStats(0);
        frame->Draw("axis");

        TLine *xline = new TLine(xing_run, 0, xing_run, 4.5);
        xline->SetLineStyle(2);
        xline->SetLineColor(kGray + 1);
        xline->Draw("same");

        float ly = 0.88;
        for (int i = 0; i < (int)defs.size(); i++)
        {
            TGraphErrors *grB = getGr(i, "gr_tight_noniso");
            TGraphErrors *grD = getGr(i, "gr_nontight_noniso");
            TGraphErrors *grRatio = buildRatio(grD, grB);
            if (!grRatio) continue;

            grRatio->SetMarkerStyle(defs[i].marker);
            grRatio->SetMarkerSize(0.4);
            grRatio->SetMarkerColor(defs[i].color);
            grRatio->SetLineColor(defs[i].color);
            grRatio->SetLineWidth(1);
            grRatio->Draw("p same");

            myMarkerLineText(0.62, ly, 0.4, defs[i].color, defs[i].marker,
                             defs[i].color, 1,
                             defs[i].label.c_str(), 0.035, true);
            ly -= 0.055;
        }

        myText(0.12, 0.92, 1, strleg1.c_str(), 0.04);
        myText(0.12, 0.87, 1, strleg2_1.c_str(), 0.035);
        myText(0.12, 0.82, 1, strleg3.c_str(), 0.035);
        myText(0.42, 0.92, kGray + 1, "0 mrad", 0.03);
        myText(0.56, 0.92, kGray + 1, "1.5 mrad", 0.03);

        c->SaveAs(Form("%s/rbr_ratio_DB.pdf", savePath.c_str()));
        delete c;
    }

    // Close files
    for (auto *f : files)
        if (f) f->Close();

    std::cout << "Plots saved to " << savePath << "/" << std::endl;
}
