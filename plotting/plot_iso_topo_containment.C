#include "plotcommon.h"

// Read /sphenix/user/shuhangli/ppg12/efficiencytool/results/iso_topo_containment.root
// Produce three PDFs under plotting/figures/:
//   iso_topo_containment_vs_R.pdf      — f_contain vs R, inclusive + 4 pT bins, pre/post common
//   iso_topo_dist_data_vs_mc.pdf       — 6-panel (1 per R), data (marker) vs inclusive MC (line), post-common
//   iso_topo_dist_overlay.pdf          — 3-panel (truth-matched / inclusive MC / data), 6 R overlaid, post-common

namespace {
const char *INPUT = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/iso_topo_containment.root";
const char *FIGDIR = "figures";

// R set and palette
struct Rentry { const char *rs; double r; int color; int style; };
const std::vector<Rentry> RADII = {
    {"005",  0.05,  kBlack,     20},
    {"0075", 0.075, kBlue + 1,  21},
    {"01",   0.10,  kGreen + 2, 22},
    {"02",   0.20,  kMagenta + 1, 33},
    {"03",   0.30,  kOrange + 7, 34},
    {"04",   0.40,  kRed + 1,   29},
};

// -------------------------------------------------------------------------
// Figure 1: photon containment <f_contain(R)> vs R, for inclusive + 4 pT bins
void plot_containment_vs_R(TFile *fin)
{
    TTree *t = (TTree *)fin->Get("fcont_summary");
    if (!t) { std::cerr << "fcont_summary tree missing\n"; return; }
    // Branches: cut (string), pt_lo (float), pt_hi (float), r (float), mean (float), rms (float)
    char cut[16];
    Long64_t pt_lo, pt_hi;
    Double_t rv;
    Float_t mean, rms;
    t->SetBranchAddress("cut",    cut);
    t->SetBranchAddress("pt_lo",  &pt_lo);
    t->SetBranchAddress("pt_hi",  &pt_hi);
    t->SetBranchAddress("r",      &rv);
    t->SetBranchAddress("mean",   &mean);
    t->SetBranchAddress("rms",    &rms);

    struct Bin { std::string label; float lo, hi; int color; int style; };
    std::vector<Bin> bins = {
        {"inclusive",  8, 40, kBlack,     20},
        {"p_{T}#in[8,14] GeV",  8, 14, kBlue + 1, 21},
        {"p_{T}#in[14,20] GeV", 14, 20, kGreen + 2, 22},
        {"p_{T}#in[20,28] GeV", 20, 28, kOrange + 7, 33},
        {"p_{T}#in[28,40] GeV", 28, 40, kRed + 1, 34},
    };

    // collect (cut, bin) -> vector of (R, mean, rms)
    std::map<std::pair<std::string, std::pair<int,int>>, std::vector<std::tuple<double,double,double>>> data;
    int N = t->GetEntries();
    for (int i = 0; i < N; ++i) {
        t->GetEntry(i);
        auto key = std::make_pair(std::string(cut),
                                  std::make_pair((int)pt_lo, (int)pt_hi));
        data[key].push_back(std::make_tuple(rv, mean, rms));
    }

    TCanvas *c = new TCanvas("c_contain", "", 1400, 600);
    c->Divide(2, 1, 0.001, 0.001);

    auto draw_side = [&](const char *cutname, const char *title){
        TH2F *frame = new TH2F(Form("frame_%s", cutname), "", 10, 0, 0.45, 10, 0.3, 1.3);
        frame->SetTitle(";Inner cone radius R;#LT(iso_{topo}^{R}+E_{T}^{cluster})/p_{T}^{#gamma, truth}#GT");
        frame->Draw("axis");

        TLine *lu = new TLine(0, 1.0, 0.45, 1.0);
        lu->SetLineStyle(2); lu->SetLineColor(kGray + 2); lu->Draw("L same");

        TLegend *leg = new TLegend(0.18, 0.18, 0.55, 0.42);
        legStyle(leg, 0.12, 0.032);

        for (auto &b : bins) {
            auto key = std::make_pair(std::string(cutname),
                                      std::make_pair((int)b.lo, (int)b.hi));
            auto it = data.find(key);
            if (it == data.end()) continue;
            auto &pts = it->second;
            int n = pts.size();
            TGraphErrors *g = new TGraphErrors(n);
            for (int i = 0; i < n; ++i) {
                g->SetPoint(i, std::get<0>(pts[i]), std::get<1>(pts[i]));
                g->SetPointError(i, 0, 0);  // no RMS band, keep clean
            }
            g->SetMarkerStyle(b.style);
            g->SetMarkerColor(b.color);
            g->SetLineColor(b.color);
            g->SetMarkerSize(1.3);
            g->SetLineWidth(2);
            g->Draw("PL same");
            leg->AddEntry(g, b.label.c_str(), "pl");
        }
        leg->Draw("same");

        myText(0.22, 0.90, 1, strleg1.c_str(), 0.045, 0);
        myText(0.22, 0.85, 1, strleg3.c_str(), 0.042, 0);
        myText(0.22, 0.80, 1, title, 0.040, 0);
    };

    c->cd(1);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    draw_side("pre",  "Pre-common (truth-matched direct photons)");
    c->cd(2);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    draw_side("post", "Post-common (npb_score > 0.5)");

    TString out = Form("%s/iso_topo_containment_vs_R.pdf", FIGDIR);
    c->SaveAs(out);
    std::cout << "wrote " << out << "\n";
}

// -------------------------------------------------------------------------
// Figure 2: 6-panel iso_topo distribution per R, data vs inclusive MC (post-common)
void plot_iso_vs_R_data_vs_mc(TFile *fin)
{
    TCanvas *c = new TCanvas("c_iso6", "", 1500, 900);
    c->Divide(3, 2, 0.005, 0.005);
    for (size_t i = 0; i < RADII.size(); ++i) {
        c->cd(i + 1);
        gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.17);
        gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.05);
        gPad->SetLogy();
        const auto &re = RADII[i];
        TH1F *hd  = (TH1F *)fin->Get(Form("h_iso_d_%s_post", re.rs));
        TH1F *hm  = (TH1F *)fin->Get(Form("h_iso_mc_%s_post", re.rs));
        if (!hd || !hm) continue;
        // MC line
        hm->SetLineColor(kBlue + 1); hm->SetLineWidth(2);
        hm->SetFillStyle(0);
        hm->SetTitle(";iso_{topo}^{R} [GeV];Normalised");
        hm->GetXaxis()->SetTitleSize(0.06);
        hm->GetXaxis()->SetLabelSize(0.05);
        hm->GetYaxis()->SetTitleSize(0.06);
        hm->GetYaxis()->SetLabelSize(0.05);
        hm->GetYaxis()->SetRangeUser(1e-4, 1);
        hm->Draw("HIST");
        // Data markers — explicit, no fill
        hd->SetMarkerStyle(20); hd->SetMarkerSize(0.8);
        hd->SetMarkerColor(kBlack); hd->SetLineColor(kBlack);
        hd->SetFillStyle(0); hd->SetLineWidth(1);
        hd->Draw("PE SAME");

        // Per-panel R label (top-right)
        myText(0.70, 0.87, 1, Form("R = %g", re.r), 0.07, 0);

        if (i == 0) {
            TLegend *leg = new TLegend(0.55, 0.62, 0.94, 0.80);
            legStyle(leg, 0.15, 0.05);
            leg->AddEntry(hd, "Data", "lep");
            leg->AddEntry(hm, "Inclusive MC", "l");
            leg->Draw("same");
            myText(0.21, 0.90, 1, strleg1.c_str(), 0.06, 0);
            myText(0.21, 0.83, 1, "post-common", 0.05, 0);
        }
    }
    TString out = Form("%s/iso_topo_dist_data_vs_mc.pdf", FIGDIR);
    c->SaveAs(out);
    std::cout << "wrote " << out << "\n";
}

// -------------------------------------------------------------------------
// Figure 3: 3-panel overlay (truth-MC / inclusive-MC / data), 6 R overlaid
void plot_overlay_per_sample(TFile *fin)
{
    TCanvas *c = new TCanvas("c_overlay", "", 1800, 600);
    c->Divide(3, 1, 0.003, 0.005);
    const char *samples[3]   = {"tr", "mc", "d"};
    const char *titles[3]    = {"MC truth-matched direct #gamma", "Inclusive MC", "Data"};
    for (int si = 0; si < 3; ++si) {
        c->cd(si + 1);
        gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.17);
        gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.05);
        gPad->SetLogy();
        // Legend top-right, 2 columns; opaque white to cover curves passing through
        TLegend *leg = new TLegend(0.60, 0.55, 0.94, 0.92);
        legStyle(leg, 0.12, 0.038);
        leg->SetNColumns(2);
        leg->SetFillStyle(1001); leg->SetFillColor(kWhite);
        bool first = true;
        for (size_t i = 0; i < RADII.size(); ++i) {
            const auto &re = RADII[i];
            TH1F *h = (TH1F *)fin->Get(Form("h_iso_%s_%s_post", samples[si], re.rs));
            if (!h) continue;
            h->SetLineColor(re.color); h->SetLineWidth(2);
            h->SetFillStyle(0);
            if (first) {
                h->SetTitle(";iso_{topo}^{R} [GeV];Normalised");
                h->GetYaxis()->SetRangeUser(1e-4, 1);
                h->GetXaxis()->SetTitleSize(0.06);
                h->GetXaxis()->SetLabelSize(0.05);
                h->GetYaxis()->SetTitleSize(0.06);
                h->GetYaxis()->SetLabelSize(0.05);
                h->Draw("HIST");
                first = false;
            } else {
                h->Draw("HIST SAME");
            }
            leg->AddEntry(h, Form("R = %g", re.r), "l");
        }
        leg->Draw("same");
        // sPHENIX top-left line 1, sample title line 2, post-common line 3
        if (si == 0) myText(0.22, 0.90, 1, strleg1.c_str(), 0.06, 0);
        myText(0.22, (si == 0 ? 0.83 : 0.90), 1, titles[si], 0.048, 0);
        myText(0.22, (si == 0 ? 0.78 : 0.85), 1, "post-common", 0.040, 0);
    }
    TString out = Form("%s/iso_topo_dist_overlay.pdf", FIGDIR);
    c->SaveAs(out);
    std::cout << "wrote " << out << "\n";
}
} // namespace

void plot_iso_topo_containment()
{
    init_plot();

    TFile *fin = TFile::Open(INPUT, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "cannot open " << INPUT << std::endl;
        return;
    }

    gSystem->mkdir(FIGDIR, true);

    plot_containment_vs_R(fin);
    plot_iso_vs_R_data_vs_mc(fin);
    plot_overlay_per_sample(fin);

    fin->Close();
}
