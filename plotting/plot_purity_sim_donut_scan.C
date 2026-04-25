#include "plotcommon.h"

// =====================================================================
// plot_purity_sim_donut_scan.C
//
// Overlay truth-matched signal purity (g_purity_truth) vs cluster pT
// from donut-iso variant Photon_final_*.root files on two canvases:
//   1) donutFull: R_inner = 0.05 / 0.075 / 0.20
//   2) donutExcl: R_inner = 0.05 / 0.075 / 0.10 / 0.20
// Missing files are skipped with a warning — partial PDFs are still
// produced so the macro can be run while Pass 2 is still in flight.
// =====================================================================

namespace
{
void drawDonutScan(const std::vector<std::string> &files,
                   const std::vector<std::string> &labels,
                   const std::vector<int> &colors,
                   const std::vector<int> &markers,
                   const std::string &resultsDir,
                   const std::string &saveName,
                   const std::string &canvasName)
{
    TCanvas *c = new TCanvas(canvasName.c_str(), "", 700, 600);
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.06);

    frame_et_rec->GetYaxis()->SetRangeUser(0.0, 1.1);
    frame_et_rec->SetYTitle("Truth-matched signal purity");
    frame_et_rec->Draw("axis");

    TLegend *leg = new TLegend(0.55, 0.22, 0.92, 0.42);
    legStyle(leg, 0.18, 0.038);

    std::vector<TGraphAsymmErrors *> keep; // keep clones alive
    int nDrawn = 0;

    for (size_t i = 0; i < files.size(); ++i)
    {
        const std::string path = resultsDir + "/" + files[i];
        TFile *fin = TFile::Open(path.c_str(), "READ");
        if (!fin || fin->IsZombie())
        {
            std::cerr << "[plot_purity_sim_donut_scan] WARNING: cannot open "
                      << path << " — skipping" << std::endl;
            if (fin) { fin->Close(); delete fin; }
            continue;
        }
        TGraphAsymmErrors *g = (TGraphAsymmErrors *)fin->Get("g_purity_truth");
        if (!g)
        {
            std::cerr << "[plot_purity_sim_donut_scan] WARNING: g_purity_truth "
                      << "not found in " << path << " — skipping" << std::endl;
            fin->Close();
            delete fin;
            continue;
        }
        TGraphAsymmErrors *gc = (TGraphAsymmErrors *)g->Clone(
            Form("g_purity_truth_clone_%zu_%s", i, canvasName.c_str()));
        gc->SetMarkerColor(colors[i]);
        gc->SetLineColor(colors[i]);
        gc->SetMarkerStyle(markers[i]);
        gc->SetMarkerSize(1.2);
        gc->SetLineWidth(2);
        gc->Draw("PZ same");
        leg->AddEntry(gc, labels[i].c_str(), "pl");
        keep.push_back(gc);
        ++nDrawn;

        fin->Close();
        delete fin;
    }

    leg->Draw("same");

    // Watermark: 4 lines in header area (top-left)
    const float xpos = 0.20, ypos = 0.885, dy = 0.054, fs = 0.040;
    const std::string strHeader =
        "Donut iso (all-calo, 120 MeV): #it{R}_{outer} = 0.4";
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),   fs, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(), fs, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),   fs, 0);
    myText(xpos, ypos - 3 * dy, 1, strHeader.c_str(), fs, 0);

    c->SaveAs(saveName.c_str());

    std::cout << "[plot_purity_sim_donut_scan] " << saveName
              << "  (" << nDrawn << "/" << files.size() << " curves drawn)"
              << std::endl;

    delete c;
}
} // anonymous namespace

void plot_purity_sim_donut_scan()
{
    init_plot();

    const std::string resultsDir =
        "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    const std::string figDir =
        "/sphenix/user/shuhangli/ppg12/plotting/figures";
    gSystem->mkdir(figDir.c_str(), true);

    // -------- donutFull: R = 0.05, 0.075, 0.20 (3 curves) --------
    {
        std::vector<std::string> files = {
            "Photon_final_bdt_donutFull_005_mc.root",
            "Photon_final_bdt_donutFull_0075_mc.root",
            "Photon_final_bdt_donutFull_02_mc.root",
        };
        std::vector<std::string> labels = {
            "#it{R}_{inner} = 0.05",
            "#it{R}_{inner} = 0.075",
            "#it{R}_{inner} = 0.20",
        };
        std::vector<int> colors  = {kBlue + 1, kRed + 1, kGreen + 2};
        std::vector<int> markers = {20, 21, 22};

        drawDonutScan(files, labels, colors, markers, resultsDir,
                      figDir + "/purity_sim_donutFull_scan.pdf",
                      "c_purity_donutFull");
    }

    // -------- donutExcl: R = 0.05, 0.075, 0.10, 0.20 (4 curves) --------
    {
        std::vector<std::string> files = {
            "Photon_final_bdt_donutExcl_005_mc.root",
            "Photon_final_bdt_donutExcl_0075_mc.root",
            "Photon_final_bdt_donutExcl_01_mc.root",
            "Photon_final_bdt_donutExcl_02_mc.root",
        };
        std::vector<std::string> labels = {
            "#it{R}_{inner} = 0.05",
            "#it{R}_{inner} = 0.075",
            "#it{R}_{inner} = 0.10",
            "#it{R}_{inner} = 0.20",
        };
        std::vector<int> colors  = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 2};
        std::vector<int> markers = {20, 21, 22, 23};

        drawDonutScan(files, labels, colors, markers, resultsDir,
                      figDir + "/purity_sim_donutExcl_scan.pdf",
                      "c_purity_donutExcl");
    }
}
