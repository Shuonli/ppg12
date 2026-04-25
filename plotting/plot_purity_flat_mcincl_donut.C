#include "plotcommon.h"

// =====================================================================
// plot_purity_flat_mcincl_donut.C
//
// Donut-style single-panel overlays for the inclusive-MC (photon+jet
// cocktail) closure test. Reads g_purity_truth and g_mc_purity_fit_ratio
// from Photon_final_bdt_*_mcincl.root.
//
// Produces two PDFs:
//   purity_sim_flat_mcincl_scan.pdf     — truth-matched signal purity
//   purity_closure_flat_mcincl_scan.pdf — ABCD / truth ratio
//
// ntbdt variants are absent because ntbdtmin05/10 crash the Brent fit
// in the inclusive-MC run (see purity_nonclosure_ntbdt.tex §5 caveats).
// =====================================================================

namespace
{
void drawPuritySim(const std::vector<std::string> &files,
                   const std::vector<std::string> &labels,
                   const std::vector<int> &colors,
                   const std::vector<int> &markers,
                   const std::string &resultsDir,
                   const std::string &saveName,
                   const std::string &canvasName,
                   const std::string &headerExtra,
                   double ymin, double ymax,
                   const std::string &ytitle)
{
    TCanvas *c = new TCanvas(canvasName.c_str(), "", 700, 600);
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.06);

    frame_et_rec->GetYaxis()->SetRangeUser(ymin, ymax);
    frame_et_rec->SetYTitle(ytitle.c_str());
    frame_et_rec->Draw("axis");

    if (ytitle.find("ratio") != std::string::npos ||
        ytitle.find("/") != std::string::npos)
    {
        TLine *unity = new TLine(8, 1.0, 40, 1.0);
        unity->SetLineColor(kGray + 2);
        unity->SetLineStyle(2);
        unity->SetLineWidth(1);
        unity->Draw("same");
    }

    TLegend *leg = new TLegend(0.55, 0.22, 0.92, 0.42);
    legStyle(leg, 0.18, 0.036);

    std::vector<TGraphAsymmErrors *> keep_asym;
    std::vector<TGraphErrors *>      keep_sym;
    int nDrawn = 0;

    const bool is_ratio = ytitle.find("ratio") != std::string::npos ||
                          ytitle.find("/") != std::string::npos;

    for (size_t i = 0; i < files.size(); ++i)
    {
        const std::string path = resultsDir + "/" + files[i];
        TFile *fin = TFile::Open(path.c_str(), "READ");
        if (!fin || fin->IsZombie())
        {
            std::cerr << "WARN: cannot open " << path << std::endl;
            if (fin) { fin->Close(); delete fin; }
            continue;
        }

        if (is_ratio)
        {
            // ratio uses g_mc_purity_fit_ratio (TGraphErrors)
            TGraphErrors *g = (TGraphErrors *)fin->Get("g_mc_purity_fit_ratio");
            if (!g)
            {
                std::cerr << "WARN: no g_mc_purity_fit_ratio in " << path << std::endl;
                fin->Close(); delete fin; continue;
            }
            TGraphErrors *gc = (TGraphErrors *)g->Clone(
                Form("gratio_%zu_%s", i, canvasName.c_str()));
            gc->SetMarkerColor(colors[i]); gc->SetLineColor(colors[i]);
            gc->SetMarkerStyle(markers[i]); gc->SetMarkerSize(1.2);
            gc->SetLineWidth(2);
            gc->Draw("PZ same");
            leg->AddEntry(gc, labels[i].c_str(), "pl");
            keep_sym.push_back(gc);
        }
        else
        {
            // sim uses g_purity_truth (TGraphAsymmErrors)
            TGraphAsymmErrors *g = (TGraphAsymmErrors *)fin->Get("g_purity_truth");
            if (!g)
            {
                std::cerr << "WARN: no g_purity_truth in " << path << std::endl;
                fin->Close(); delete fin; continue;
            }
            TGraphAsymmErrors *gc = (TGraphAsymmErrors *)g->Clone(
                Form("gpsim_%zu_%s", i, canvasName.c_str()));
            gc->SetMarkerColor(colors[i]); gc->SetLineColor(colors[i]);
            gc->SetMarkerStyle(markers[i]); gc->SetMarkerSize(1.2);
            gc->SetLineWidth(2);
            gc->Draw("PZ same");
            leg->AddEntry(gc, labels[i].c_str(), "pl");
            keep_asym.push_back(gc);
        }
        ++nDrawn;
        fin->Close(); delete fin;
    }
    leg->Draw("same");

    const float xpos = 0.20, ypos = 0.885, dy = 0.054, fs = 0.040;
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),    fs, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(),  fs, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),    fs, 0);
    myText(xpos, ypos - 3 * dy, 1, headerExtra.c_str(), fs, 0);

    c->SaveAs(saveName.c_str());
    std::cout << saveName << "  (" << nDrawn << "/" << files.size()
              << " curves drawn)" << std::endl;
    delete c;
}
} // anonymous namespace

void plot_purity_flat_mcincl_donut()
{
    init_plot();
    const std::string resultsDir = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    const std::string figDir = "/sphenix/user/shuhangli/ppg12/plotting/figures";
    gSystem->mkdir(figDir.c_str(), true);

    const std::vector<std::string> files = {
        "Photon_final_bdt_nom_mcincl.root",
        "Photon_final_bdt_flat_t85_nt50_mcincl.root",
        "Photon_final_bdt_flat_t90_nt50_mcincl.root",
        "Photon_final_bdt_flat_t90_nt70_mcincl.root",
        "Photon_final_bdt_flat_t95_nt50_mcincl.root",
    };
    const std::vector<std::string> labels = {
        "parametric nom (ref)",
        "flat t=[0.85,1], nt=[0.50,0.85]",
        "flat t=[0.90,1], nt=[0.50,0.90]",
        "flat t=[0.90,1], nt=[0.70,0.90]",
        "flat t=[0.95,1], nt=[0.50,0.95]",
    };
    const std::vector<int> colors  = {kGray + 2, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 2};
    const std::vector<int> markers = {24, 20, 21, 22, 23};

    drawPuritySim(files, labels, colors, markers, resultsDir,
                  figDir + "/purity_sim_flat_mcincl_scan.pdf",
                  "c_psim_flat_mcincl",
                  "Flat BDT scan  (inclusive MC: #gamma+jet)",
                  0.0, 1.1, "Truth-matched signal purity");

    drawPuritySim(files, labels, colors, markers, resultsDir,
                  figDir + "/purity_closure_flat_mcincl_scan.pdf",
                  "c_close_flat_mcincl",
                  "Flat BDT scan  (inclusive MC: #gamma+jet)",
                  -0.5, 2.0, "ABCD purity / truth purity");
}
