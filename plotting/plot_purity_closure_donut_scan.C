#include "plotcommon.h"

// plot_purity_closure_donut_scan.C
// Overlay ABCD-on-MC / truth purity ratio (g_mc_purity_fit_ratio) vs cluster pT
// for donut variants. Closure = 1.0 is perfect.

namespace
{
void drawClosureScan(const std::vector<std::string> &files,
                     const std::vector<std::string> &labels,
                     const std::vector<int> &colors,
                     const std::vector<int> &markers,
                     const std::string &resultsDir,
                     const std::string &saveName,
                     const std::string &canvasName,
                     const std::string &headerExtra)
{
    TCanvas *c = new TCanvas(canvasName.c_str(), "", 700, 600);
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.06);

    frame_et_rec->GetYaxis()->SetRangeUser(0.3, 2.0);
    frame_et_rec->SetYTitle("ABCD purity / truth purity");
    frame_et_rec->Draw("axis");

    TLine *unity = new TLine(8, 1.0, 40, 1.0);
    unity->SetLineColor(kGray + 2);
    unity->SetLineStyle(2);
    unity->SetLineWidth(1);
    unity->Draw("same");

    TLegend *leg = new TLegend(0.55, 0.68, 0.92, 0.88);
    legStyle(leg, 0.18, 0.038);

    std::vector<TGraphErrors *> keep;
    int nDrawn = 0;

    for (size_t i = 0; i < files.size(); ++i)
    {
        const std::string path = resultsDir + "/" + files[i];
        TFile *fin = TFile::Open(path.c_str(), "READ");
        if (!fin || fin->IsZombie())
        {
            std::cerr << "[closure] WARNING: cannot open " << path << std::endl;
            if (fin) { fin->Close(); delete fin; }
            continue;
        }
        TGraphErrors *g = (TGraphErrors *)fin->Get("g_mc_purity_fit_ratio");
        if (!g)
        {
            std::cerr << "[closure] WARNING: no g_mc_purity_fit_ratio in " << path << std::endl;
            fin->Close();
            delete fin;
            continue;
        }
        TGraphErrors *gc = (TGraphErrors *)g->Clone(
            Form("g_closure_clone_%zu_%s", i, canvasName.c_str()));
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

    const float xpos = 0.20, ypos = 0.885, dy = 0.054, fs = 0.040;
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),   fs, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(), fs, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),   fs, 0);
    myText(xpos, ypos - 3 * dy, 1, headerExtra.c_str(), fs, 0);

    c->SaveAs(saveName.c_str());
    std::cout << "[closure] " << saveName << "  (" << nDrawn << "/" << files.size()
              << " curves drawn)" << std::endl;
}
}

void plot_purity_closure_donut_scan()
{
    init_plot();

    const std::string resultsDir = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    const std::string figDir = "/sphenix/user/shuhangli/ppg12/plotting/figures";
    gSystem->mkdir(figDir.c_str(), true);

    const std::string hdrFull = "Donut iso (all-calo, 120 MeV), #it{R}_{outer} = 0.4";
    const std::string hdrExcl = "Donut iso excl. EMCal-owned towers, #it{R}_{outer} = 0.4";

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

        drawClosureScan(files, labels, colors, markers, resultsDir,
                        figDir + "/purity_closure_donutFull_scan.pdf",
                        "c_closure_donutFull", hdrFull);
    }

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

        drawClosureScan(files, labels, colors, markers, resultsDir,
                        figDir + "/purity_closure_donutExcl_scan.pdf",
                        "c_closure_donutExcl", hdrExcl);
    }
}
