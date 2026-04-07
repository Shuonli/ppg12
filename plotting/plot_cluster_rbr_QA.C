#include "plotcommon.h"
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TSystem.h>
#include <vector>
#include <string>
#include <algorithm>

// -----------------------------------------------------------------------
// plot_cluster_rbr_QA.C
// Plots run-by-run QA quantities from the output of Cluster_rbr.C
// (rbrQA.root).  Produces PDF files in the current directory.
//
// Usage (ROOT prompt):
//   .x plot_cluster_rbr_QA.C("path/to/rbrQA.root")
// -----------------------------------------------------------------------

// Helper: draw a TGraphErrors vs run number on a fresh canvas with a
// floating y-axis.  Returns the canvas (caller owns it).
TCanvas *drawRunGraph(TGraphErrors *gr,
                      const char *cname,
                      const char *ytitle,
                      int color = kBlack,
                      int marker = 20)
{
    TCanvas *c = new TCanvas(cname, "", 800, 500);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.14);
    c->SetTopMargin(0.07);
    c->SetRightMargin(0.05);

    gr->SetMarkerStyle(marker);
    gr->SetMarkerSize(0.5);
    gr->SetMarkerColor(color);
    gr->SetLineColor(color);
    gr->GetXaxis()->SetTitle("Run Number");
    gr->GetYaxis()->SetTitle(ytitle);
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleSize(0.05);
    gr->GetXaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetTitleOffset(1.1);
    gr->Draw("ap");
    return c;
}

// Helper: overlay multiple TGraphErrors on a shared canvas / frame.
// grlist = {graph, legend_label, color, marker}
struct GrEntry
{
    TGraphErrors *gr;
    std::string   label;
    int           color;
    int           marker;
};

TCanvas *drawMultiRunGraph(const std::vector<GrEntry> &entries,
                           const char *cname,
                           const char *ytitle,
                           float leg_x = 0.55,
                           float leg_y_top = 0.88)
{
    TCanvas *c = new TCanvas(cname, "", 800, 500);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.14);
    c->SetTopMargin(0.07);
    c->SetRightMargin(0.05);

    // Determine global y-range
    double ymin =  1e30;
    double ymax = -1e30;
    double xmin =  1e30;
    double xmax = -1e30;
    for (const auto &e : entries)
    {
        if (!e.gr || e.gr->GetN() == 0) continue;
        e.gr->ComputeRange(xmin, ymin, xmax, ymax);
    }
    if (ymin > ymax) { ymin = 0; ymax = 1; }
    double yrange = ymax - ymin;
    if (yrange <= 0) yrange = 1;
    double ypad = 0.12 * yrange;

    // Draw a blank frame to set axes
    TH1F *frame = new TH1F(Form("frame_%s", cname), "",
                           100,
                           xmin - 0.01 * (xmax - xmin),
                           xmax + 0.01 * (xmax - xmin));
    frame->GetYaxis()->SetRangeUser(std::max(0.0, ymin - ypad), ymax + ypad);
    frame->GetXaxis()->SetTitle("Run Number");
    frame->GetYaxis()->SetTitle(ytitle);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetTitleOffset(1.1);
    frame->SetStats(0);
    frame->Draw("axis");

    float leg_dy = 0.06f;
    float ly = leg_y_top;
    for (const auto &e : entries)
    {
        if (!e.gr || e.gr->GetN() == 0) continue;
        e.gr->SetMarkerStyle(e.marker);
        e.gr->SetMarkerSize(0.5);
        e.gr->SetMarkerColor(e.color);
        e.gr->SetLineColor(e.color);
        e.gr->Draw("p same");
        myMarkerLineText(leg_x, ly, 0.5, e.color, e.marker, e.color, 1,
                         e.label.c_str(), 0.040, true);
        ly -= leg_dy;
    }
    return c;
}

void plot_cluster_rbr_QA(const std::string &infile = "../efficiencytool/rbrQA_topo.root")
{
    init_plot();

    TFile *fin = TFile::Open(infile.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "Cannot open " << infile << std::endl;
        return;
    }

    auto getGrE = [&](const char *name) -> TGraphErrors *
    {
        TGraphErrors *g = dynamic_cast<TGraphErrors *>(fin->Get(name));
        if (!g) std::cerr << "WARNING: " << name << " not found in file\n";
        return g;
    };
    auto getGr = [&](const char *name) -> TGraph *
    {
        TGraph *g = dynamic_cast<TGraph *>(fin->Get(name));
        if (!g) std::cerr << "WARNING: " << name << " not found in file\n";
        return g;
    };

    TGraphErrors *gr_event          = getGrE("gr_event");
    TGraphErrors *gr_common         = getGrE("gr_common");
    TGraphErrors *gr_tight          = getGrE("gr_tight");
    TGraphErrors *gr_nontight       = getGrE("gr_nontight");
    TGraphErrors *gr_tight_iso      = getGrE("gr_tight_iso");
    TGraphErrors *gr_tight_noniso   = getGrE("gr_tight_noniso");
    TGraphErrors *gr_nontight_iso   = getGrE("gr_nontight_iso");
    TGraphErrors *gr_nontight_noniso= getGrE("gr_nontight_noniso");
    TGraphErrors *gr_corr           = getGrE("gr_corr");
    TGraphErrors *gr_avg_isoET      = getGrE("gr_avg_isoET");
    TGraph       *gr_min_scaler     = getGr ("gr_min_scaler");

    gSystem->Exec("mkdir -p figures/rbr_QA");

    const std::string savePath = "figures/rbr_QA";

    // Standard labels
    const char *sph_label = "#bf{#it{sPHENIX}} Internal";
    const char *eta_label = "|#it{#eta}^{#gamma}| < 0.7";

    // ---------------------------------------------------------------
    // Plot 1: Event rate vs run number
    // ---------------------------------------------------------------
    if (gr_event)
    {
        TCanvas *c = drawRunGraph(gr_event, "c_event",
                                  "Events / Lumi [mb]",
                                  kBlack, 20);
        myText(0.15, 0.88, 1, sph_label, 0.04);
        myText(0.15, 0.83, 1, eta_label, 0.04);
        c->SaveAs(Form("%s/rbr_event_rate.pdf", savePath.c_str()));
        delete c;
    }

    // ---------------------------------------------------------------
    // Plot 2: Min scaler vs run number
    // ---------------------------------------------------------------
    if (gr_min_scaler)
    {
        TCanvas *c = new TCanvas("c_scaler", "", 800, 500);
        c->SetLeftMargin(0.12);
        c->SetBottomMargin(0.14);
        c->SetTopMargin(0.07);
        c->SetRightMargin(0.05);
        gr_min_scaler->SetMarkerStyle(20);
        gr_min_scaler->SetMarkerSize(0.5);
        gr_min_scaler->SetMarkerColor(kBlack);
        gr_min_scaler->SetLineColor(kBlack);
        gr_min_scaler->GetXaxis()->SetTitle("Run Number");
        gr_min_scaler->GetYaxis()->SetTitle("Min Scaler (scaled trigger)");
        gr_min_scaler->GetXaxis()->SetTitleSize(0.05);
        gr_min_scaler->GetYaxis()->SetTitleSize(0.05);
        gr_min_scaler->GetXaxis()->SetLabelSize(0.04);
        gr_min_scaler->GetYaxis()->SetLabelSize(0.04);
        gr_min_scaler->GetYaxis()->SetTitleOffset(1.1);
        gr_min_scaler->Draw("ap");
        myText(0.15, 0.88, 1, sph_label, 0.04);
        c->SaveAs(Form("%s/rbr_min_scaler.pdf", savePath.c_str()));
        delete c;
    }

    // ---------------------------------------------------------------
    // Plot 3: Cluster rates — common, tight, nontight
    // ---------------------------------------------------------------
    if (gr_common && gr_tight && gr_nontight)
    {
        std::vector<GrEntry> entries = {
            {gr_common,   "Pre-selected", kBlack,    20},
            {gr_tight,    "Tight",        kBlue,     21},
            {gr_nontight, "Non-tight",    kRed,      22},
        };
        TCanvas *c = drawMultiRunGraph(entries, "c_cluster_rates",
                                       "Clusters / Lumi [mb^{-1}]",
                                       0.55, 0.88);
        myText(0.15, 0.88, 1, sph_label, 0.04);
        myText(0.15, 0.83, 1, eta_label, 0.04);
        c->SaveAs(Form("%s/rbr_cluster_rates.pdf", savePath.c_str()));
        delete c;
    }

    // ---------------------------------------------------------------
    // Plot 4: ABCD region cluster rates — one canvas per region
    // ---------------------------------------------------------------
    struct ABCDEntry
    {
        TGraphErrors *gr;
        const char   *tag;
        const char   *label;
        int           color;
    };
    std::vector<ABCDEntry> abcd = {
        {gr_tight_iso,       "A_tight_iso",       "A: tight iso",       kBlack  },
        {gr_tight_noniso,    "B_tight_noniso",    "B: tight noniso",    kRed    },
        {gr_nontight_iso,    "C_nontight_iso",    "C: nontight iso",    kBlue   },
        {gr_nontight_noniso, "D_nontight_noniso", "D: nontight noniso", kMagenta},
    };
    for (const auto &e : abcd)
    {
        if (!e.gr) continue;
        TCanvas *c = drawRunGraph(e.gr,
                                  Form("c_%s", e.tag),
                                  "Clusters / Lumi [mb^{-1}]",
                                  e.color, 20);
        myText(0.15, 0.88, 1, sph_label, 0.04);
        myText(0.15, 0.83, 1, eta_label, 0.04);
        myText(0.15, 0.78, 1, e.label,   0.04);
        c->SaveAs(Form("%s/rbr_rate_%s.pdf", savePath.c_str(), e.tag));
        delete c;
    }

    // ---------------------------------------------------------------
    // Plot 5: Average isolation ET vs run number
    // ---------------------------------------------------------------
    if (gr_avg_isoET)
    {
        TCanvas *c = drawRunGraph(gr_avg_isoET, "c_avg_isoET",
                                  "#LT#it{E}_{T}^{iso}#GT [GeV]",
                                  kBlack, 20);
        myText(0.15, 0.88, 1, sph_label, 0.04);
        myText(0.15, 0.83, 1, eta_label, 0.04);
        myText(0.15, 0.78, 1, "Pre-selected clusters", 0.04);
        c->SaveAs(Form("%s/rbr_avg_isoET.pdf", savePath.c_str()));
        delete c;
    }

    // ---------------------------------------------------------------
    // Plot 6: Lumi correction factor vs run number
    // ---------------------------------------------------------------
    if (gr_corr)
    {
        TCanvas *c = drawRunGraph(gr_corr, "c_corr",
                                  "Lumi correction (corrected / nominal)",
                                  kBlack, 20);
        myText(0.15, 0.88, 1, sph_label, 0.04);
        c->SaveAs(Form("%s/rbr_lumi_correction.pdf", savePath.c_str()));
        delete c;
    }

    // ---------------------------------------------------------------
    // Plot 7: Ratio C/A (nontight_iso / tight_iso) vs run number
    //         Key stability check for the sideband method
    // ---------------------------------------------------------------
    if (gr_tight_iso && gr_nontight_iso &&
        gr_tight_iso->GetN() > 0 && gr_nontight_iso->GetN() > 0)
    {
        // Build ratio graph point-by-point (graphs share the same run ordering
        // since both come from the same map iteration in Cluster_rbr.C)
        int nA = gr_tight_iso->GetN();
        int nC = gr_nontight_iso->GetN();
        int npts = std::min(nA, nC);

        std::vector<double> x_r, y_r, ex_r, ey_r;
        for (int i = 0; i < npts; i++)
        {
            double xA, yA, xC, yC;
            gr_tight_iso->GetPoint(i, xA, yA);
            gr_nontight_iso->GetPoint(i, xC, yC);
            if (xA != xC) continue;   // mismatched run numbers — skip
            if (yA <= 0)  continue;
            double eA = gr_tight_iso->GetErrorY(i);
            double eC = gr_nontight_iso->GetErrorY(i);
            double ratio = yC / yA;
            double err   = ratio * std::sqrt((eA/yA)*(eA/yA) + (eC/yC > 0 ? (eC/yC)*(eC/yC) : 0));
            x_r.push_back(xA);
            y_r.push_back(ratio);
            ex_r.push_back(0);
            ey_r.push_back(err);
        }

        if (!x_r.empty())
        {
            TGraphErrors *gr_ratio_CA = new TGraphErrors(
                (int)x_r.size(), x_r.data(), y_r.data(), ex_r.data(), ey_r.data());
            gr_ratio_CA->SetName("gr_ratio_CA");

            TCanvas *c = drawRunGraph(gr_ratio_CA, "c_ratio_CA",
                                      "C / A  (nontight iso / tight iso)",
                                      kBlack, 20);
            myText(0.15, 0.88, 1, sph_label, 0.04);
            myText(0.15, 0.83, 1, eta_label, 0.04);
            c->SaveAs(Form("%s/rbr_ratio_CA.pdf", savePath.c_str()));
            delete c;
            delete gr_ratio_CA;
        }
    }

    // ---------------------------------------------------------------
    // Plot 8: Ratio D/B (nontight_noniso / tight_noniso) vs run number
    // ---------------------------------------------------------------
    if (gr_tight_noniso && gr_nontight_noniso &&
        gr_tight_noniso->GetN() > 0 && gr_nontight_noniso->GetN() > 0)
    {
        int nB = gr_tight_noniso->GetN();
        int nD = gr_nontight_noniso->GetN();
        int npts = std::min(nB, nD);

        std::vector<double> x_r, y_r, ex_r, ey_r;
        for (int i = 0; i < npts; i++)
        {
            double xB, yB, xD, yD;
            gr_tight_noniso->GetPoint(i, xB, yB);
            gr_nontight_noniso->GetPoint(i, xD, yD);
            if (xB != xD) continue;
            if (yB <= 0)  continue;
            double eB = gr_tight_noniso->GetErrorY(i);
            double eD = gr_nontight_noniso->GetErrorY(i);
            double ratio = yD / yB;
            double err   = ratio * std::sqrt((eB/yB)*(eB/yB) + (eD/yD > 0 ? (eD/yD)*(eD/yD) : 0));
            x_r.push_back(xB);
            y_r.push_back(ratio);
            ex_r.push_back(0);
            ey_r.push_back(err);
        }

        if (!x_r.empty())
        {
            TGraphErrors *gr_ratio_DB = new TGraphErrors(
                (int)x_r.size(), x_r.data(), y_r.data(), ex_r.data(), ey_r.data());
            gr_ratio_DB->SetName("gr_ratio_DB");

            TCanvas *c = drawRunGraph(gr_ratio_DB, "c_ratio_DB",
                                      "D / B  (nontight noniso / tight noniso)",
                                      kBlack, 20);
            myText(0.15, 0.88, 1, sph_label, 0.04);
            myText(0.15, 0.83, 1, eta_label, 0.04);
            c->SaveAs(Form("%s/rbr_ratio_DB.pdf", savePath.c_str()));
            delete c;
            delete gr_ratio_DB;
        }
    }

    // ---------------------------------------------------------------
    // Plot 9: Summary — tight_iso and nontight_iso on the same canvas
    //         (signal region and sideband iso side by side)
    // ---------------------------------------------------------------
    if (gr_tight_iso && gr_nontight_iso)
    {
        std::vector<GrEntry> entries = {
            {gr_tight_iso,    "A: tight iso",    kBlue,  20},
            {gr_nontight_iso, "C: nontight iso", kRed,   21},
        };
        TCanvas *c = drawMultiRunGraph(entries, "c_iso_comparison",
                                       "Clusters / Lumi [mb^{-1}]",
                                       0.55, 0.88);
        myText(0.15, 0.88, 1, sph_label, 0.04);
        myText(0.15, 0.83, 1, eta_label, 0.04);
        c->SaveAs(Form("%s/rbr_iso_comparison.pdf", savePath.c_str()));
        delete c;
    }

    fin->Close();
    std::cout << "Done. Plots saved to " << savePath << "/" << std::endl;
}
