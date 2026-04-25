#include "plotcommon.h"

// =====================================================================
// plot_purity_closure_ntbdt_flat_scan.C
//
// Donut-style single-panel overlays of ABCD/truth purity ratio
// (g_mc_purity_fit_ratio) applied to the BDT-scan variants in
// purity_nonclosure_ntbdt:
//   - ntbdt scan: 3 variants, truth and region A identical, so the ratio
//     curves differ only in the ABCD numerator
//   - flat BDT partitions: 5 variants (nom ref + 4 flat), each with its
//     own truth denominator
//
// Matches plot_purity_closure_donut_scan.C aesthetic.
//
// Outputs:
//   purity_closure_ntbdt_scan.pdf
//   purity_closure_flat_scan.pdf
// =====================================================================

namespace
{
void drawClosure(const std::vector<std::string> &files,
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

    TLegend *leg = new TLegend(0.55, 0.66, 0.92, 0.90);
    legStyle(leg, 0.18, 0.036);

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
            std::cerr << "[closure] WARNING: no g_mc_purity_fit_ratio in "
                      << path << std::endl;
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
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),    fs, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(),  fs, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),    fs, 0);
    myText(xpos, ypos - 3 * dy, 1, headerExtra.c_str(), fs, 0);

    c->SaveAs(saveName.c_str());
    std::cout << "[closure] " << saveName << "  (" << nDrawn << "/" << files.size()
              << " curves drawn)" << std::endl;

    delete c;
}
} // anonymous namespace

void plot_purity_closure_ntbdt_flat_scan()
{
    init_plot();

    const std::string resultsDir =
        "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    const std::string figDir =
        "/sphenix/user/shuhangli/ppg12/plotting/figures";
    gSystem->mkdir(figDir.c_str(), true);

    // -------- ntbdt scan: 3 variants (nt_bdt_min in {0.02, 0.05, 0.10}) --------
    {
        std::vector<std::string> files = {
            "Photon_final_bdt_nom_mc.root",
            "Photon_final_bdt_ntbdtmin05_mc.root",
            "Photon_final_bdt_ntbdtmin10_mc.root",
        };
        std::vector<std::string> labels = {
            "nt_bdt_min = 0.02 (nom)",
            "nt_bdt_min = 0.05",
            "nt_bdt_min = 0.10",
        };
        std::vector<int> colors  = {kBlue + 1, kRed + 1, kGreen + 2};
        std::vector<int> markers = {20, 21, 22};

        drawClosure(files, labels, colors, markers, resultsDir,
                    figDir + "/purity_closure_ntbdt_scan.pdf",
                    "c_close_ntbdt",
                    "Parametric scan: nt_bdt_min lower bound");
    }

    // -------- flat BDT: 5 variants --------
    {
        std::vector<std::string> files = {
            "Photon_final_bdt_nom_mc.root",
            "Photon_final_bdt_flat_t85_nt50_mc.root",
            "Photon_final_bdt_flat_t90_nt50_mc.root",
            "Photon_final_bdt_flat_t90_nt70_mc.root",
            "Photon_final_bdt_flat_t95_nt50_mc.root",
        };
        std::vector<std::string> labels = {
            "parametric nom (ref)",
            "flat t=[0.85,1], nt=[0.50,0.85]",
            "flat t=[0.90,1], nt=[0.50,0.90]",
            "flat t=[0.90,1], nt=[0.70,0.90]",
            "flat t=[0.95,1], nt=[0.50,0.95]",
        };
        std::vector<int> colors  = {kGray + 2, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 2};
        std::vector<int> markers = {24, 20, 21, 22, 23};

        drawClosure(files, labels, colors, markers, resultsDir,
                    figDir + "/purity_closure_flat_scan.pdf",
                    "c_close_flat",
                    "Flat BDT partition scan");
    }
}
