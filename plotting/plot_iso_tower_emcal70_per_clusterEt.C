#include "plotcommon.h"

// Tower-based iso, EMCal-only with 70 MeV per-tower threshold.
// Scans R in {0.05, 0.075, 0.10, 0.20, 0.30}.
// 4 pages (one per cluster-E_T bin), 3 panels (truth-matched / inclusive MC / data).
//
// Output: figures/iso_tower_emcal70_per_clEt.pdf

namespace {
const char *INPUT  = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/iso_topo_containment.root";
const char *FIGDIR = "figures";

struct Rentry { const char *rs; double r; int color; };
const std::vector<Rentry> RADII_EMCAL70 = {
    {"005",  0.05,  kBlack},
    {"0075", 0.075, kBlue + 1},
    {"01",   0.10,  kGreen + 2},
    {"02",   0.20,  kMagenta + 1},
    {"03",   0.30,  kOrange + 7},
};

struct ClEtBin { int lo, hi; };
const std::vector<ClEtBin> CLET_BINS = {{8, 14}, {14, 20}, {20, 28}, {28, 40}};
} // namespace

void plot_iso_tower_emcal70_per_clusterEt()
{
    init_plot();

    TFile *fin = TFile::Open(INPUT, "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "cannot open " << INPUT << "\n"; return; }
    gSystem->mkdir(FIGDIR, true);

    TString outpdf = Form("%s/iso_tower_emcal70_per_clEt.pdf", FIGDIR);
    TCanvas *c = new TCanvas("c_iso_tw_emcal70", "", 1800, 600);
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

            TLegend *leg = new TLegend(0.67, 0.58, 0.94, 0.92);
            legStyle(leg, 0.12, 0.035);
            leg->SetNColumns(2);
            leg->SetFillStyle(1001); leg->SetFillColor(kWhite);
            leg->SetBorderSize(0);

            bool first = true; double ymax = 0;
            std::vector<TH1F *> hists;
            for (const auto &re : RADII_EMCAL70) {
                TString key = Form("h_iso_tw_emcal70_%s_%s_post_clEt_%d_%d",
                                   re.rs, samples[si], eb.lo, eb.hi);
                TH1F *h = (TH1F *)fin->Get(key);
                if (!h) { std::cerr << "missing " << key << "\n"; hists.push_back(nullptr); continue; }
                h->SetLineColor(re.color); h->SetLineWidth(2); h->SetFillStyle(0);
                if (h->GetMaximum() > ymax) ymax = h->GetMaximum();
                hists.push_back(h);
            }
            for (size_t i = 0; i < hists.size(); ++i) {
                if (!hists[i]) continue;
                TH1F *h = hists[i];
                if (first) {
                    h->SetTitle(";iso^{R}_{tower} (EMCAL only, 70 MeV thr.) [GeV];Normalised");
                    h->GetYaxis()->SetRangeUser(0, ymax * 1.25);
                    h->GetXaxis()->SetTitleSize(0.055);
                    h->GetXaxis()->SetLabelSize(0.05);
                    h->GetYaxis()->SetTitleSize(0.06);
                    h->GetYaxis()->SetLabelSize(0.05);
                    h->Draw("HIST");
                    first = false;
                } else {
                    h->Draw("HIST SAME");
                }
                leg->AddEntry(h, Form("R = %g", RADII_EMCAL70[i].r), "l");
            }
            leg->Draw("same");

            if (si == 0) myText(0.22, 0.90, 1, strleg1.c_str(), 0.055, 0);
            myText(0.22, (si == 0 ? 0.83 : 0.90), 1, titles[si], 0.047, 0);
            myText(0.22, (si == 0 ? 0.77 : 0.84), 1,
                   Form("E_{T} #in [%d, %d] GeV, EMCal 70 MeV", eb.lo, eb.hi), 0.038, 0);
        }
        c->Print(outpdf);
    }
    c->Print(outpdf + "]");
    std::cout << "wrote " << outpdf << " (" << CLET_BINS.size() << " pages)\n";
    fin->Close();
}
