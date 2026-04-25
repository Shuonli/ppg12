#include "plotcommon.h"

// Tower-based iso at R=0.3, threshold scan (EMCal+HCalIn+HCalOut).
// Compares no-threshold, 60 MeV, 70 MeV, and 120 MeV per-tower thresholds.
// 4 pages (one per cluster-E_T bin), 3 panels (truth-matched / inclusive MC / data).
//
// Output: figures/iso_tower_thresh_r03_per_clEt.pdf

namespace {
const char *INPUT  = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/iso_topo_containment.root";
const char *FIGDIR = "figures";

struct Ventry { const char *vkey; const char *label; int color; };
const std::vector<Ventry> VARS = {
    {"thr_nothr_03", "no threshold",  kBlack},
    {"thr_60_03",    "60 MeV",        kBlue + 1},
    {"thr_70_03",    "70 MeV",        kOrange + 7},
    {"thr_120_03",   "120 MeV",       kRed + 1},
};

struct ClEtBin { int lo, hi; };
const std::vector<ClEtBin> CLET_BINS = {{8, 14}, {14, 20}, {20, 28}, {28, 40}};
} // namespace

void plot_iso_tower_thresh_r03_per_clusterEt()
{
    init_plot();

    TFile *fin = TFile::Open(INPUT, "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "cannot open " << INPUT << "\n"; return; }
    gSystem->mkdir(FIGDIR, true);

    TString outpdf = Form("%s/iso_tower_thresh_r03_per_clEt.pdf", FIGDIR);
    TCanvas *c = new TCanvas("c_iso_tw_thr_r03", "", 1800, 600);
    c->Print(outpdf + "[");

    const char *samples[3] = {"tr", "mc", "d"};
    const char *titles[3]  = {"Truth-matched #gamma MC", "Inclusive MC", "Data"};

    for (const auto &eb : CLET_BINS) {
        c->Clear();
        c->Divide(3, 1, 0.003, 0.005);
        for (int si = 0; si < 3; ++si) {
            c->cd(si + 1);
            gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.17);
            gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.05);

            TLegend *leg = new TLegend(0.62, 0.62, 0.94, 0.92);
            legStyle(leg, 0.12, 0.035);
            leg->SetFillStyle(1001); leg->SetFillColor(kWhite);
            leg->SetBorderSize(0);

            bool first = true; double ymax = 0;
            std::vector<TH1F *> hists;
            for (const auto &v : VARS) {
                TString key = Form("h_iso_tw_%s_%s_post_clEt_%d_%d",
                                   v.vkey, samples[si], eb.lo, eb.hi);
                TH1F *h = (TH1F *)fin->Get(key);
                if (!h) { std::cerr << "missing " << key << "\n"; hists.push_back(nullptr); continue; }
                h->SetLineColor(v.color); h->SetLineWidth(2); h->SetFillStyle(0);
                if (h->GetMaximum() > ymax) ymax = h->GetMaximum();
                hists.push_back(h);
            }
            for (size_t i = 0; i < hists.size(); ++i) {
                if (!hists[i]) continue;
                TH1F *h = hists[i];
                if (first) {
                    h->SetTitle(";iso^{R=0.3}_{tower} (EMCAL+HCALIN+HCALOUT) [GeV];Normalised");
                    h->GetYaxis()->SetRangeUser(0, ymax * 1.35);
                    h->GetXaxis()->SetTitleSize(0.055);
                    h->GetXaxis()->SetLabelSize(0.05);
                    h->GetYaxis()->SetTitleSize(0.06);
                    h->GetYaxis()->SetLabelSize(0.05);
                    h->Draw("HIST");
                    first = false;
                } else {
                    h->Draw("HIST SAME");
                }
                leg->AddEntry(h, VARS[i].label, "l");
            }
            leg->Draw("same");

            if (si == 0) myText(0.22, 0.90, 1, strleg1.c_str(), 0.055, 0);
            myText(0.22, (si == 0 ? 0.83 : 0.90), 1, titles[si], 0.047, 0);
            myText(0.22, (si == 0 ? 0.77 : 0.84), 1,
                   Form("E_{T} #in [%d, %d] GeV, R=0.3", eb.lo, eb.hi), 0.038, 0);
        }
        c->Print(outpdf);
    }
    c->Print(outpdf + "]");
    std::cout << "wrote " << outpdf << " (" << CLET_BINS.size() << " pages)\n";
    fin->Close();
}
