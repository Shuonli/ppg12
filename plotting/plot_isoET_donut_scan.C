#include "plotcommon.h"

// plot_isoET_donut_scan.C
// Overlay isoET distributions (h_tight_isoET_0_{ipt}) across donut variants for
// two representative pT bins. Normalized to unit area.

namespace {
void drawIsoETpanel(const std::vector<std::string> &files,
                    const std::vector<std::string> &labels,
                    const std::vector<int> &colors,
                    const std::vector<int> &markers,
                    const std::string &resdir,
                    int ipt,
                    const std::string &ptLabel,
                    const std::string &saveName,
                    const std::string &canvasName,
                    const std::string &headerExtra)
{
    TCanvas *c = new TCanvas(canvasName.c_str(), "", 700, 600);
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.06);

    gPad->SetLogy();
    TH1F *frame = new TH1F(Form("frame_%s", canvasName.c_str()), "", 1, -2, 10);
    frame->SetMinimum(5e-4);
    frame->SetMaximum(3.0);
    frame->SetXTitle("E_{T}^{iso} [GeV]");
    frame->SetYTitle("Normalized entries / bin");
    frame->Draw("axis");

    TLegend *leg = new TLegend(0.60, 0.35, 0.93, 0.58);
    legStyle(leg, 0.18, 0.034);

    std::vector<TH1D *> keep;
    for (size_t i = 0; i < files.size(); ++i) {
        const std::string path = resdir + "/data_histo_" + files[i] + ".root";
        TFile *fin = TFile::Open(path.c_str(), "READ");
        if (!fin || fin->IsZombie()) {
            std::cerr << "[isoET] cannot open " << path << std::endl;
            if (fin) delete fin; continue;
        }
        TH1D *h = (TH1D*)fin->Get(Form("h_tight_isoET_0_%d", ipt));
        if (!h) { fin->Close(); delete fin; continue; }
        TH1D *hc = (TH1D*)h->Clone(Form("isoet_%zu_%s", i, canvasName.c_str()));
        hc->SetDirectory(0);
        hc->Rebin(2);
        if (hc->Integral() > 0) hc->Scale(1.0/hc->Integral());
        hc->SetMarkerColor(colors[i]);
        hc->SetLineColor(colors[i]);
        hc->SetMarkerStyle(markers[i]);
        hc->SetMarkerSize(0.9);
        hc->SetLineWidth(2);
        hc->Draw("E1 same");
        leg->AddEntry(hc, labels[i].c_str(), "pl");
        keep.push_back(hc);
        fin->Close(); delete fin;
    }
    leg->Draw("same");

    const float xpos = 0.20, ypos = 0.885, dy = 0.050, fs = 0.036;
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),   fs, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(), fs, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),   fs, 0);
    myText(xpos, ypos - 3 * dy, 1, headerExtra.c_str(), fs, 0);
    myText(xpos, ypos - 4 * dy, 1, ("Tight data, " + ptLabel).c_str(), fs, 0);

    c->SaveAs(saveName.c_str());
    std::cout << "[isoET] " << saveName << std::endl;
}
}

void plot_isoET_donut_scan()
{
    init_plot();
    const std::string resdir = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    const std::string figdir = "/sphenix/user/shuhangli/ppg12/plotting/figures";

    std::vector<std::string> fullFiles  = {"bdt_donutFull_005","bdt_donutFull_0075","bdt_donutFull_02"};
    std::vector<std::string> fullLabels = {"#it{R}_{inner}=0.05","#it{R}_{inner}=0.075","#it{R}_{inner}=0.20"};
    std::vector<int>         fullColors = {kBlue+1, kRed+1, kGreen+2};
    std::vector<int>         fullMarks  = {20, 21, 22};

    std::vector<std::string> exclFiles  = {"bdt_donutExcl_005","bdt_donutExcl_0075","bdt_donutExcl_01","bdt_donutExcl_02"};
    std::vector<std::string> exclLabels = {"#it{R}_{inner}=0.05","#it{R}_{inner}=0.075","#it{R}_{inner}=0.10","#it{R}_{inner}=0.20"};
    std::vector<int>         exclColors = {kBlue+1, kRed+1, kGreen+2, kMagenta+2};
    std::vector<int>         exclMarks  = {20, 21, 22, 23};

    const std::string hdrFull = "donutFull, R_{outer}=0.4";
    const std::string hdrExcl = "donutExcl, R_{outer}=0.4";

    // Low pT: ipt=1 (10-12 GeV). High pT: ipt=7 (22-24 GeV).
    drawIsoETpanel(fullFiles, fullLabels, fullColors, fullMarks, resdir,
                   1, "10 < p_{T} < 12 GeV",
                   figdir + "/isoET_donutFull_ptlow.pdf",  "c_iso_full_lo", hdrFull);
    drawIsoETpanel(fullFiles, fullLabels, fullColors, fullMarks, resdir,
                   7, "22 < p_{T} < 24 GeV",
                   figdir + "/isoET_donutFull_pthigh.pdf", "c_iso_full_hi", hdrFull);
    drawIsoETpanel(exclFiles, exclLabels, exclColors, exclMarks, resdir,
                   1, "10 < p_{T} < 12 GeV",
                   figdir + "/isoET_donutExcl_ptlow.pdf",  "c_iso_excl_lo", hdrExcl);
    drawIsoETpanel(exclFiles, exclLabels, exclColors, exclMarks, resdir,
                   7, "22 < p_{T} < 24 GeV",
                   figdir + "/isoET_donutExcl_pthigh.pdf", "c_iso_excl_hi", hdrExcl);
}
