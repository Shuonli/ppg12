#include "plotcommon.h"

// Overlay ratio of donut cross-section to nominal for both families.
// Cross section proxy: h_unfold_sub_result (TH1F).

namespace
{
void drawRatio(const std::vector<std::string> &files,
               const std::vector<std::string> &labels,
               const std::vector<int> &colors,
               const std::vector<int> &markers,
               const std::string &resultsDir,
               const std::string &nomPath,
               const std::string &saveName,
               const std::string &canvasName,
               const std::string &headerExtra)
{
    TCanvas *c = new TCanvas(canvasName.c_str(), "", 700, 600);
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.06);

    frame_et_rec->GetYaxis()->SetRangeUser(0.6, 1.4);
    frame_et_rec->SetYTitle("Variant / nominal cross section");
    frame_et_rec->Draw("axis");

    TLine *unity = new TLine(8, 1.0, 40, 1.0);
    unity->SetLineColor(kGray + 2);
    unity->SetLineStyle(2);
    unity->SetLineWidth(1);
    unity->Draw("same");

    TFile *fnom = TFile::Open(nomPath.c_str(), "READ");
    if (!fnom || fnom->IsZombie())
    {
        std::cerr << "[ratio] cannot open nominal " << nomPath << std::endl;
        return;
    }
    TH1F *hnom = (TH1F*)fnom->Get("h_unfold_sub_result");
    if (!hnom)
    {
        std::cerr << "[ratio] no h_unfold_sub_result in nominal" << std::endl;
        fnom->Close(); delete fnom; return;
    }

    TLegend *leg = new TLegend(0.55, 0.68, 0.92, 0.88);
    legStyle(leg, 0.18, 0.038);

    std::vector<TH1F *> keep;
    for (size_t i = 0; i < files.size(); ++i)
    {
        const std::string path = resultsDir + "/" + files[i];
        TFile *fin = TFile::Open(path.c_str(), "READ");
        if (!fin || fin->IsZombie()) { if (fin) delete fin; continue; }
        TH1F *hv = (TH1F*)fin->Get("h_unfold_sub_result");
        if (!hv) { fin->Close(); delete fin; continue; }
        TH1F *r = (TH1F*)hv->Clone(Form("ratio_%zu_%s", i, canvasName.c_str()));
        r->SetDirectory(0);
        r->Divide(hnom);
        r->SetMarkerColor(colors[i]);
        r->SetLineColor(colors[i]);
        r->SetMarkerStyle(markers[i]);
        r->SetMarkerSize(1.2);
        r->SetLineWidth(2);
        r->Draw("E1 same");
        leg->AddEntry(r, labels[i].c_str(), "pl");
        keep.push_back(r);
        fin->Close(); delete fin;
    }
    leg->Draw("same");

    const float xpos = 0.20, ypos = 0.885, dy = 0.054, fs = 0.040;
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),   fs, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(), fs, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),   fs, 0);
    myText(xpos, ypos - 3 * dy, 1, headerExtra.c_str(), fs, 0);

    c->SaveAs(saveName.c_str());
    std::cout << "[ratio] " << saveName << std::endl;
    fnom->Close(); delete fnom;
}
}

void plot_xsec_ratio_donut_scan()
{
    init_plot();
    const std::string resultsDir = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    const std::string figDir = "/sphenix/user/shuhangli/ppg12/plotting/figures";
    const std::string nomPath = resultsDir + "/Photon_final_bdt_nom.root";

    {
        std::vector<std::string> files = {
            "Photon_final_bdt_donutFull_005.root",
            "Photon_final_bdt_donutFull_0075.root",
            "Photon_final_bdt_donutFull_02.root",
        };
        std::vector<std::string> labels = {"#it{R}_{inner}=0.05","#it{R}_{inner}=0.075","#it{R}_{inner}=0.20"};
        std::vector<int> colors = {kBlue+1, kRed+1, kGreen+2};
        std::vector<int> markers = {20,21,22};
        drawRatio(files, labels, colors, markers, resultsDir, nomPath,
                  figDir + "/xsec_ratio_donutFull_scan.pdf", "c_xsec_donutFull",
                  "donutFull / nom cross section");
    }
    {
        std::vector<std::string> files = {
            "Photon_final_bdt_donutExcl_005.root",
            "Photon_final_bdt_donutExcl_0075.root",
            "Photon_final_bdt_donutExcl_01.root",
            "Photon_final_bdt_donutExcl_02.root",
        };
        std::vector<std::string> labels = {"#it{R}_{inner}=0.05","#it{R}_{inner}=0.075","#it{R}_{inner}=0.10","#it{R}_{inner}=0.20"};
        std::vector<int> colors = {kBlue+1, kRed+1, kGreen+2, kMagenta+2};
        std::vector<int> markers = {20,21,22,23};
        drawRatio(files, labels, colors, markers, resultsDir, nomPath,
                  figDir + "/xsec_ratio_donutExcl_scan.pdf", "c_xsec_donutExcl",
                  "donutExcl / nom cross section");
    }
}
