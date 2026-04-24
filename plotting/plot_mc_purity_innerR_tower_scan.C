#include "plotcommon.h"

// Overlay MC purity closure ratio (g_mc_purity_fit_ratio = f_fit / f_truth)
// for the TOWER-based inner-R scan: 70 MeV outer cone R=0.3 (emcal+hcalin+hcalout)
// minus EMCal-only inner ring at 70 MeV. 4 inner-R values x 2 reweight states.
// Includes the nominal topo R=0.4 closure (no inner ring) as a reference.
//
// Inputs: efficiencytool/results/Photon_final_bdt_{nom, innerRtower_*{,_nore}}_mc.root
// Output: plotting/figures/mc_purity_innerR_tower_scan.pdf

namespace {
struct Variant
{
    const char *suffix;
    const char *label;
    int color;
    int style;
};

std::vector<TFile *> g_files;
std::vector<TGraphErrors *> g_graphs;

void draw_panel(const Variant *vars, int n,
                const char *stateLabel,
                const char *resultsDir)
{
    frame_et_rec->SetTitle(";#it{E}_{T}^{#gamma} [GeV];f_{purity}^{fit} / f_{purity}^{truth}");
    frame_et_rec->GetXaxis()->SetRangeUser(10, 35);
    frame_et_rec->GetYaxis()->SetRangeUser(0.6, 1.45);
    frame_et_rec->Draw("axis");

    lineone->SetLineColor(kGray + 2);
    lineone->SetLineStyle(2);
    lineone->Draw("L same");

    TLegend *leg = new TLegend(0.56, 0.66, 0.94, 0.92);
    legStyle(leg, 0.12, 0.034);

    for (int i = 0; i < n; ++i)
    {
        const Variant &v = vars[i];
        TString path = Form("%s/Photon_final_bdt_%s_mc.root", resultsDir, v.suffix);
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
        g_files.push_back(f);
        g_graphs.push_back(g);
    }
    leg->Draw("same");

    myText(0.18, 0.90, 1, strleg1.c_str(), 0.045, 0);
    myText(0.18, 0.85, 1, strleg3.c_str(), 0.042, 0);
    myText(0.18, 0.30, 1, stateLabel, 0.040, 0);
    myText(0.18, 0.25, 1, "Tower 70 MeV donut, R_{out}=0.3", 0.036, 0);
}
} // namespace

void plot_mc_purity_innerR_tower_scan()
{
    init_plot();

    const char *resultsDir = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";

    Variant reweight_on[5] = {
        {"nom",              "nominal topo R = 0.4",  kGray + 2,   24},
        {"innerRtower_005",  "tower inner R = 0.05",  kBlack,      20},
        {"innerRtower_0075", "tower inner R = 0.075", kBlue + 1,   21},
        {"innerRtower_01",   "tower inner R = 0.10",  kGreen + 2,  22},
        {"innerRtower_02",   "tower inner R = 0.20",  kRed + 1,    33},
    };
    Variant reweight_off[5] = {
        {"no_unfolding_reweighting", "nominal topo R = 0.4", kGray + 2,  24},
        {"innerRtower_005_nore",  "tower inner R = 0.05",  kBlack,      20},
        {"innerRtower_0075_nore", "tower inner R = 0.075", kBlue + 1,   21},
        {"innerRtower_01_nore",   "tower inner R = 0.10",  kGreen + 2,  22},
        {"innerRtower_02_nore",   "tower inner R = 0.20",  kRed + 1,    33},
    };

    TCanvas *c = new TCanvas("c_innerRtower_scan", "", 1400, 600);
    c->Divide(2, 1, 0.001, 0.001);

    c->cd(1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.14);
    draw_panel(reweight_on, 5, "Reweight: ON (nominal)", resultsDir);

    c->cd(2);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.14);
    draw_panel(reweight_off, 5, "Reweight: OFF", resultsDir);

    c->SaveAs("figures/mc_purity_innerR_tower_scan.pdf");

    std::cout << "\nWrote figures/mc_purity_innerR_tower_scan.pdf with "
              << g_graphs.size() << " graphs total\n";

    for (auto *f : g_files)
        if (f)
            f->Close();
}
