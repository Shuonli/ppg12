#include "plotcommon.h"

// plot_isoET_donut_mc_scan.C
// Same layout as plot_isoET_donut_scan.C (data) but for:
//   sigMC  = MC_efficiency_bdt_<var>.root  (merged photon: h_tight_isoET_0_{ipt})
//   incMC  = sigMC + MC_efficiency_jet_bdt_<var>.root (merged jet added)
// Both normalized to unit area. Per donut family (Full, Excl), at low/high pT bins.

namespace {
TH1D* loadIso(const std::string &path, int ipt)
{
    TFile *f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) { if (f) delete f; return nullptr; }
    TH1D *h = (TH1D*)f->Get(Form("h_tight_isoET_0_%d", ipt));
    if (!h) { f->Close(); delete f; return nullptr; }
    TH1D *hc = (TH1D*)h->Clone(Form("iso_loader_%p", (void*)f));
    hc->SetDirectory(0);
    f->Close(); delete f;
    return hc;
}

void drawIsoETmc(const std::vector<std::string> &files,
                 const std::vector<std::string> &labels,
                 const std::vector<int> &colors,
                 const std::vector<int> &markers,
                 const std::string &resdir,
                 int ipt,
                 const std::string &ptLabel,
                 const std::string &mcMode,          // "sig" or "inc"
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
        std::string path = (mcMode == "inc")
            ? resdir + "/MC_efficiency_jet_" + files[i] + ".root"
            : resdir + "/MC_efficiency_"     + files[i] + ".root";
        TH1D *h = loadIso(path, ipt);
        if (!h) { std::cerr << "[isoET-mc] miss: " << path << std::endl; continue; }
        if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
        h->SetMarkerColor(colors[i]);
        h->SetLineColor(colors[i]);
        h->SetMarkerStyle(markers[i]);
        h->SetMarkerSize(0.9);
        h->SetLineWidth(2);
        h->Draw("E1 same");
        leg->AddEntry(h, labels[i].c_str(), "pl");
        keep.push_back(h);
    }
    leg->Draw("same");

    const std::string mcLabel = (mcMode == "sig") ? "Truth-matched MC signal" : "Inclusive MC (jet)";
    const float xpos = 0.20, ypos = 0.885, dy = 0.050, fs = 0.036;
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),   fs, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(), fs, 0);
    myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),   fs, 0);
    myText(xpos, ypos - 3 * dy, 1, headerExtra.c_str(), fs, 0);
    myText(xpos, ypos - 4 * dy, 1, (mcLabel + ", tight, " + ptLabel).c_str(), fs, 0);

    c->SaveAs(saveName.c_str());
    std::cout << "[isoET-mc] " << saveName << std::endl;
}
}

void plot_isoET_donut_mc_scan()
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

    // sigMC (signal-only photon MC)
    drawIsoETmc(fullFiles, fullLabels, fullColors, fullMarks, resdir, 1, "10 < p_{T} < 12 GeV", "sig",
                figdir + "/isoET_donutFull_ptlow_sigMC.pdf",  "c_iso_full_lo_sig", hdrFull);
    drawIsoETmc(fullFiles, fullLabels, fullColors, fullMarks, resdir, 7, "22 < p_{T} < 24 GeV", "sig",
                figdir + "/isoET_donutFull_pthigh_sigMC.pdf", "c_iso_full_hi_sig", hdrFull);
    drawIsoETmc(exclFiles, exclLabels, exclColors, exclMarks, resdir, 1, "10 < p_{T} < 12 GeV", "sig",
                figdir + "/isoET_donutExcl_ptlow_sigMC.pdf",  "c_iso_excl_lo_sig", hdrExcl);
    drawIsoETmc(exclFiles, exclLabels, exclColors, exclMarks, resdir, 7, "22 < p_{T} < 24 GeV", "sig",
                figdir + "/isoET_donutExcl_pthigh_sigMC.pdf", "c_iso_excl_hi_sig", hdrExcl);

    // incMC (photon + jet)
    drawIsoETmc(fullFiles, fullLabels, fullColors, fullMarks, resdir, 1, "10 < p_{T} < 12 GeV", "inc",
                figdir + "/isoET_donutFull_ptlow_incMC.pdf",  "c_iso_full_lo_inc", hdrFull);
    drawIsoETmc(fullFiles, fullLabels, fullColors, fullMarks, resdir, 7, "22 < p_{T} < 24 GeV", "inc",
                figdir + "/isoET_donutFull_pthigh_incMC.pdf", "c_iso_full_hi_inc", hdrFull);
    drawIsoETmc(exclFiles, exclLabels, exclColors, exclMarks, resdir, 1, "10 < p_{T} < 12 GeV", "inc",
                figdir + "/isoET_donutExcl_ptlow_incMC.pdf",  "c_iso_excl_lo_inc", hdrExcl);
    drawIsoETmc(exclFiles, exclLabels, exclColors, exclMarks, resdir, 7, "22 < p_{T} < 24 GeV", "inc",
                figdir + "/isoET_donutExcl_pthigh_incMC.pdf", "c_iso_excl_hi_inc", hdrExcl);
}
