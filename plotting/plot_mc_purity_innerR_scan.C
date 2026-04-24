#include "plotcommon.h"

// Overlay MC purity closure ratio (g_mc_purity_fit_ratio = f_fit / f_truth)
// for 4 inner-R values x 2 reweight states (nominal + no-reweight).
// Produces a 2-panel PDF: left = reweight on (nominal), right = reweight off.
//
// Inputs: /sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_innerR_{suffix}_mc.root
// Output: plotting/figures/mc_purity_innerR_scan.pdf

namespace {
struct Variant
{
    const char *suffix;
    const char *label;
    int color;
    int style;
};

// Keep TFile + TGraph alive across panels (avoids closed-file data deallocation)
std::vector<TFile *> g_files;
std::vector<TGraphErrors *> g_graphs;

void draw_panel(const Variant *vars, int n,
                const char *stateLabel,
                const char *resultsDir)
{
    frame_et_rec->SetTitle(";#it{E}_{T}^{#gamma} [GeV];f_{purity}^{fit} / f_{purity}^{truth}");
    frame_et_rec->GetXaxis()->SetRangeUser(10, 35);
    frame_et_rec->GetYaxis()->SetRangeUser(0.6, 1.4);
    frame_et_rec->Draw("axis");

    lineone->SetLineColor(kGray + 2);
    lineone->SetLineStyle(2);
    lineone->Draw("L same");

    TLegend *leg = new TLegend(0.54, 0.72, 0.94, 0.92);
    legStyle(leg, 0.12, 0.038);

    for (int i = 0; i < n; ++i)
    {
        const Variant &v = vars[i];
        TString path = Form("%s/Photon_final_bdt_innerR_%s_mc.root", resultsDir, v.suffix);
        TFile *f = TFile::Open(path, "READ");
        if (!f || f->IsZombie())
        {
            std::cerr << "WARN: cannot open " << path << std::endl;
            continue;
        }
        TGraphErrors *g = (TGraphErrors *)f->Get("g_mc_purity_fit_ratio");
        if (!g)
        {
            std::cerr << "WARN: g_mc_purity_fit_ratio missing in " << path << std::endl;
            f->Close();
            continue;
        }
        g->SetMarkerStyle(v.style);
        g->SetMarkerColor(v.color);
        g->SetLineColor(v.color);
        g->SetLineWidth(2);
        g->SetMarkerSize(1.3);
        g->Draw("P same");
        leg->AddEntry(g, v.label, "pl");
        // Keep file + graph alive for SaveAs
        g_files.push_back(f);
        g_graphs.push_back(g);
    }
    leg->Draw("same");

    myText(0.20, 0.90, 1, strleg1.c_str(), 0.045, 0);
    myText(0.20, 0.85, 1, strleg3.c_str(), 0.042, 0);
    myText(0.20, 0.30, 1, stateLabel, 0.040, 0);
}
} // namespace

void plot_mc_purity_innerR_scan()
{
    init_plot();

    const char *resultsDir = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";

    Variant reweight_on[4] = {
        {"005",  "inner R = 0.05",  kBlack,   20},
        {"0075", "inner R = 0.075", kBlue + 1, 21},
        {"01",   "inner R = 0.10",  kGreen + 2, 22},
        {"02",   "inner R = 0.20",  kRed + 1, 33},
    };
    Variant reweight_off[4] = {
        {"005_nore",  "inner R = 0.05",  kBlack,   20},
        {"0075_nore", "inner R = 0.075", kBlue + 1, 21},
        {"01_nore",   "inner R = 0.10",  kGreen + 2, 22},
        {"02_nore",   "inner R = 0.20",  kRed + 1, 33},
    };

    TCanvas *c = new TCanvas("c_innerR_scan", "", 1400, 600);
    c->Divide(2, 1, 0.001, 0.001);

    c->cd(1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.14);
    draw_panel(reweight_on, 4, "Reweight: ON (nominal)", resultsDir);

    c->cd(2);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.14);
    draw_panel(reweight_off, 4, "Reweight: OFF", resultsDir);

    c->SaveAs("figures/mc_purity_innerR_scan.pdf");

    std::cout << "\nWrote figures/mc_purity_innerR_scan.pdf with "
              << g_graphs.size() << " graphs total\n";

    // Only NOW close the files
    for (auto *f : g_files)
        if (f)
            f->Close();
}
