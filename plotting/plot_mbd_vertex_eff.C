// plot_mbd_vertex_eff.C
//
// Plotting macro for StudyMBDVertexEff.C output.
// Reads efficiencytool/results/mbd_vertex_eff_study.root and produces:
//   plotting/figures/mbd_vertex_eff_study.pdf      (2x2 panels)
//   plotting/figures/mbd_vertex_eff_summary.pdf    (loose, photon vs jet, no reweight)
//
// Run:
//   cd plotting
//   root -l -b -q plot_mbd_vertex_eff.C+

#include "plotcommon.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>

#include <string>
#include <vector>


namespace {

const char* kInPath =
    "/sphenix/user/shuhangli/ppg12/efficiencytool/results/mbd_vertex_eff_study.root";
const char* kFigDir = "/sphenix/user/shuhangli/ppg12/plotting/figures";


TGraphAsymmErrors* GetGraph(TFile* f, const std::string& name)
{
    auto* g = dynamic_cast<TGraphAsymmErrors*>(f->Get(name.c_str()));
    if (!g) Printf("  ERROR: %s missing", name.c_str());
    return g;
}

void StyleGraph(TGraphAsymmErrors* g, int color, int marker)
{
    if (!g) return;
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.1);
    g->SetLineWidth(2);
}

}  // namespace


void plot_mbd_vertex_eff()
{
    init_plot();
    gStyle->SetOptStat(0);

    TFile* fin = TFile::Open(kInPath, "READ");
    if (!fin || fin->IsZombie()) {
        Printf("ERROR: cannot open %s", kInPath);
        return;
    }

    auto* g_p_l_nr = GetGraph(fin, "g_eff_photon_loose_noreweight");
    auto* g_p_l_0  = GetGraph(fin, "g_eff_photon_loose_0mrad");
    auto* g_p_l_15 = GetGraph(fin, "g_eff_photon_loose_1p5mrad");
    auto* g_p_b_nr = GetGraph(fin, "g_eff_photon_bundled_noreweight");
    auto* g_p_b_0  = GetGraph(fin, "g_eff_photon_bundled_0mrad");
    auto* g_p_b_15 = GetGraph(fin, "g_eff_photon_bundled_1p5mrad");

    auto* g_j_l_nr = GetGraph(fin, "g_eff_jet_loose_noreweight");
    auto* g_j_l_0  = GetGraph(fin, "g_eff_jet_loose_0mrad");
    auto* g_j_l_15 = GetGraph(fin, "g_eff_jet_loose_1p5mrad");
    auto* g_j_b_nr = GetGraph(fin, "g_eff_jet_bundled_noreweight");
    auto* g_j_b_0  = GetGraph(fin, "g_eff_jet_bundled_0mrad");
    auto* g_j_b_15 = GetGraph(fin, "g_eff_jet_bundled_1p5mrad");

    StyleGraph(g_p_l_nr, kBlack,     20);
    StyleGraph(g_p_l_0,  kRed + 1,   21);
    StyleGraph(g_p_l_15, kAzure + 2, 22);
    StyleGraph(g_p_b_nr, kBlack,     20);
    StyleGraph(g_p_b_0,  kRed + 1,   21);
    StyleGraph(g_p_b_15, kAzure + 2, 22);
    StyleGraph(g_j_l_nr, kBlack,     20);
    StyleGraph(g_j_l_0,  kRed + 1,   21);
    StyleGraph(g_j_l_15, kAzure + 2, 22);
    StyleGraph(g_j_b_nr, kBlack,     20);
    StyleGraph(g_j_b_0,  kRed + 1,   21);
    StyleGraph(g_j_b_15, kAzure + 2, 22);

    // ---- Figure 1: 2x2 (channel x definition) ----
    TCanvas* c1 = new TCanvas("c_mbd_vertex_eff", "MBD vertex efficiency",
                              1400, 1200);
    c1->Divide(2, 2, 0.001, 0.001);

    auto draw_panel = [&](int pad, const char* xtitle, const char* hdr,
                          const char* def_line,
                          const std::vector<TGraphAsymmErrors*>& gs,
                          const std::vector<std::string>& labels) {
        c1->cd(pad);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.07);

        auto* frame = new TH1F(Form("frame_p%d", pad), "", 100, 3.0, 36.0);
        frame->SetMinimum(0.0);
        frame->SetMaximum(1.15);
        frame->GetXaxis()->SetTitle(xtitle);
        frame->GetYaxis()->SetTitle("MBD vertex efficiency");
        frame->Draw();

        for (auto* g : gs) {
            if (g) g->Draw("PZ same");
        }

        TLatex tl;
        tl.SetNDC();
        tl.SetTextFont(42);
        tl.SetTextSize(0.040);
        tl.DrawLatex(0.20, 0.88, strleg1.c_str());
        tl.SetTextSize(0.034);
        tl.DrawLatex(0.20, 0.83, strMC.c_str());
        tl.DrawLatex(0.20, 0.78, hdr);
        tl.DrawLatex(0.20, 0.73, def_line);

        auto* leg = new TLegend(0.55, 0.18, 0.93, 0.42);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.032);
        leg->SetHeader("truth-z reweight:");
        for (size_t k = 0; k < gs.size(); ++k) {
            if (gs[k]) leg->AddEntry(gs[k], labels[k].c_str(), "lp");
        }
        leg->Draw();
    };

    draw_panel(1,
               "leading truth #it{p}_{T}^{#gamma} [GeV]",
               "Photon channel (photon5/10/20)",
               "loose: vz_{reco} #neq -9999",
               {g_p_l_nr, g_p_l_0, g_p_l_15},
               {"none (intrinsic)", "0 mrad", "1.5 mrad"});

    draw_panel(2,
               "leading truth #it{p}_{T}^{jet} [GeV]",
               "Jet channel (jet5/8/12/20/30/40)",
               "loose: vz_{reco} #neq -9999",
               {g_j_l_nr, g_j_l_0, g_j_l_15},
               {"none (intrinsic)", "0 mrad", "1.5 mrad"});

    draw_panel(3,
               "leading truth #it{p}_{T}^{#gamma} [GeV]",
               "Photon channel (photon5/10/20)",
               "bundled: |vz_{reco}|<60, MBD N&S #geq 1",
               {g_p_b_nr, g_p_b_0, g_p_b_15},
               {"none (intrinsic)", "0 mrad", "1.5 mrad"});

    draw_panel(4,
               "leading truth #it{p}_{T}^{jet} [GeV]",
               "Jet channel (jet5/8/12/20/30/40)",
               "bundled: |vz_{reco}|<60, MBD N&S #geq 1",
               {g_j_b_nr, g_j_b_0, g_j_b_15},
               {"none (intrinsic)", "0 mrad", "1.5 mrad"});

    c1->SaveAs(Form("%s/mbd_vertex_eff_study.pdf", kFigDir));

    // ---- Figure 2: summary overlay (loose, no reweight, photon vs jet) ----
    TCanvas* c2 = new TCanvas("c_mbd_vertex_eff_summary",
                              "MBD vertex efficiency (loose, intrinsic)",
                              900, 700);
    c2->cd();
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.06);

    auto* frame2 = new TH1F("frame_summary", "", 100, 3.0, 36.0);
    frame2->SetMinimum(0.0);
    frame2->SetMaximum(1.15);
    frame2->GetXaxis()->SetTitle("leading truth #it{p}_{T} [GeV]");
    frame2->GetYaxis()->SetTitle("MBD vertex efficiency");
    frame2->Draw();

    StyleGraph(g_p_l_nr, kBlack,   20);
    StyleGraph(g_j_l_nr, kRed + 1, 21);
    if (g_p_l_nr) g_p_l_nr->Draw("PZ same");
    if (g_j_l_nr) g_j_l_nr->Draw("PZ same");

    TLatex tl2;
    tl2.SetNDC();
    tl2.SetTextFont(42);
    tl2.SetTextSize(0.045);
    tl2.DrawLatex(0.20, 0.88, strleg1.c_str());
    tl2.SetTextSize(0.038);
    tl2.DrawLatex(0.20, 0.83, strMC.c_str());
    tl2.DrawLatex(0.20, 0.78, "loose: vz_{reco} #neq -9999");
    tl2.DrawLatex(0.20, 0.73, "no truth-z reweight (intrinsic)");

    auto* leg2 = new TLegend(0.55, 0.20, 0.93, 0.36);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.038);
    if (g_p_l_nr)
        leg2->AddEntry(g_p_l_nr,
                       "Photon (leading truth #gamma)", "lp");
    if (g_j_l_nr)
        leg2->AddEntry(g_j_l_nr,
                       "Jet (leading R=0.4 truth jet)", "lp");
    leg2->Draw();

    c2->SaveAs(Form("%s/mbd_vertex_eff_summary.pdf", kFigDir));

    fin->Close();
    Printf("Wrote:");
    Printf("  %s/mbd_vertex_eff_study.pdf",   kFigDir);
    Printf("  %s/mbd_vertex_eff_summary.pdf", kFigDir);
}
