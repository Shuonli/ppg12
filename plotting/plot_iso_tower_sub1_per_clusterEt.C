#include "plotcommon.h"

// SUB1 (background-subtracted) tower iso vs plain tower iso at R=0.3 and R=0.4.
// Both variants sum EMCal + HCalIn + HCalOut.
// 4 pages (one per cluster-E_T bin), 3 panels (truth-matched / inclusive MC / data).
//
// Output: figures/iso_tower_sub1_per_clEt.pdf

namespace {
const char *INPUT  = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/iso_topo_containment.root";
const char *FIGDIR = "figures";

struct Ventry { const char *vkey; bool is_plain; const char *label; int color; int style; };
// plain_{03,04} use existing h_iso_tw_{sample}_{R}_{cut}_{bin}  (sample before R)
// sub1_{03,04} use new-variant h_iso_tw_sub1_{R}_{sample}_{cut}_{bin}  (variant before sample)
const std::vector<Ventry> VARS = {
    {"03",      true,  "R=0.3  plain",        kBlack,    1},
    {"sub1_03", false, "R=0.3  SUB1-subtr.",  kBlue + 1, 2},
    {"04",      true,  "R=0.4  plain",        kRed + 1,  1},
    {"sub1_04", false, "R=0.4  SUB1-subtr.",  kOrange + 7, 2},
};

struct ClEtBin { int lo, hi; };
const std::vector<ClEtBin> CLET_BINS = {{8, 14}, {14, 20}, {20, 28}, {28, 40}};

TH1F *get_h(TFile *fin, const Ventry &v, const char *sample, const ClEtBin &eb)
{
    TString key;
    if (v.is_plain) {
        // existing naming: sample before R
        key = Form("h_iso_tw_%s_%s_post_clEt_%d_%d", sample, v.vkey, eb.lo, eb.hi);
    } else {
        // new-variant naming: variant before sample
        key = Form("h_iso_tw_%s_%s_post_clEt_%d_%d", v.vkey, sample, eb.lo, eb.hi);
    }
    return (TH1F *)fin->Get(key);
}
} // namespace

void plot_iso_tower_sub1_per_clusterEt()
{
    init_plot();

    TFile *fin = TFile::Open(INPUT, "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "cannot open " << INPUT << "\n"; return; }
    gSystem->mkdir(FIGDIR, true);

    TString outpdf = Form("%s/iso_tower_sub1_per_clEt.pdf", FIGDIR);
    TCanvas *c = new TCanvas("c_iso_tw_sub1", "", 1800, 600);
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

            TLegend *leg = new TLegend(0.55, 0.62, 0.94, 0.92);
            legStyle(leg, 0.12, 0.034);
            leg->SetFillStyle(1001); leg->SetFillColor(kWhite);
            leg->SetBorderSize(0);

            bool first = true; double ymax = 0;
            std::vector<TH1F *> hists;
            for (const auto &v : VARS) {
                TH1F *h = get_h(fin, v, samples[si], eb);
                if (!h) { std::cerr << "missing " << v.vkey << " " << samples[si] << "\n";
                         hists.push_back(nullptr); continue; }
                h->SetLineColor(v.color); h->SetLineStyle(v.style); h->SetLineWidth(2);
                h->SetFillStyle(0);
                if (h->GetMaximum() > ymax) ymax = h->GetMaximum();
                hists.push_back(h);
            }
            for (size_t i = 0; i < hists.size(); ++i) {
                if (!hists[i]) continue;
                TH1F *h = hists[i];
                if (first) {
                    h->SetTitle(";iso_{tower} (EMCAL+HCALIN+HCALOUT) [GeV];Normalised");
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
                   Form("E_{T} #in [%d, %d] GeV", eb.lo, eb.hi), 0.038, 0);
        }
        c->Print(outpdf);
    }
    c->Print(outpdf + "]");
    std::cout << "wrote " << outpdf << " (" << CLET_BINS.size() << " pages)\n";
    fin->Close();
}
