#include "plotcommon.h"

// =====================================================================
// plot_xsec_flat_scan.C
//
// Donut-style single-panel overlays of the extracted isolated-photon
// cross-section ("h_unfold_sub_result") for the flat-BDT-partition scan
// in the purity_nonclosure_ntbdt study.
//
// Input source: Photon_final_bdt_{nom,flat_t85_nt50,flat_t90_nt50,
// flat_t90_nt70,flat_t95_nt50}_mc.root. The "_mc" suffix means the
// data side of CalculatePhotonYield is the inclusive PYTHIA MC (the
// "jet" merge is already inclusive — HardQCD:all + PromptPhoton:all),
// so this is the MC-extraction-equivalent cross-section. The data-side
// re-run for the flat variants is pending (data_histo_bdt_flat*.root
// was filled before the njet-branch fix, all zero); this plot is the
// closest-to-nominal pipeline output available today.
//
// Outputs:
//   xsec_flat_scan.pdf        absolute d^2N_iso/dE_T dyc (log-y)
//   xsec_ratio_flat_scan.pdf  variant / nominal ratio (linear-y)
// =====================================================================

namespace
{
struct Variant
{
    std::string tag;
    std::string label;
    int color;
    int marker;
};

const std::vector<Variant> SCAN = {
    {"nom",           "parametric nom (ref)",               kGray + 2,    24},
    {"flat_t85_nt50", "flat t=[0.85,1], nt=[0.50,0.85]",    kRed + 1,     20},
    {"flat_t90_nt50", "flat t=[0.90,1], nt=[0.50,0.90]",    kBlue + 1,    21},
    {"flat_t90_nt70", "flat t=[0.90,1], nt=[0.70,0.90]",    kGreen + 2,   22},
    {"flat_t95_nt50", "flat t=[0.95,1], nt=[0.50,0.95]",    kMagenta + 2, 23},
};

TH1F *loadXsec(const std::string &tag, const std::string &resultsDir,
               const std::string &cloneName)
{
    const std::string path = resultsDir + "/Photon_final_bdt_" + tag + "_mc.root";
    TFile *fin = TFile::Open(path.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "WARN cannot open " << path << std::endl;
        if (fin) { fin->Close(); delete fin; }
        return nullptr;
    }
    TH1F *h = (TH1F *)fin->Get("h_unfold_sub_result");
    if (!h)
    {
        std::cerr << "WARN no h_unfold_sub_result in " << path << std::endl;
        fin->Close();
        delete fin;
        return nullptr;
    }
    TH1F *clone = (TH1F *)h->Clone(cloneName.c_str());
    clone->SetDirectory(0);
    fin->Close();
    delete fin;
    return clone;
}
} // anonymous namespace

void plot_xsec_flat_scan()
{
    init_plot();

    const std::string resultsDir =
        "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    const std::string figDir =
        "/sphenix/user/shuhangli/ppg12/plotting/figures";
    gSystem->mkdir(figDir.c_str(), true);

    // ---------- Absolute cross-section (log y) ----------
    {
        TCanvas *c = new TCanvas("c_xsec_flat", "", 700, 600);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.14);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.06);
        gPad->SetLogy();

        frame_et_rec->GetYaxis()->SetRangeUser(1.0, 5e5);
        frame_et_rec->SetYTitle("Extracted isolated-photon yield (a.u.)");
        frame_et_rec->Draw("axis");

        TLegend *leg = new TLegend(0.55, 0.60, 0.94, 0.90);
        legStyle(leg, 0.18, 0.034);

        std::vector<TH1F *> keep;
        for (size_t i = 0; i < SCAN.size(); ++i)
        {
            TH1F *h = loadXsec(SCAN[i].tag, resultsDir,
                               Form("hxsec_%zu", i));
            if (!h) continue;
            h->SetMarkerColor(SCAN[i].color);
            h->SetLineColor(SCAN[i].color);
            h->SetMarkerStyle(SCAN[i].marker);
            h->SetMarkerSize(1.2);
            h->SetLineWidth(2);
            h->Draw("E1 same");
            leg->AddEntry(h, SCAN[i].label.c_str(), "pl");
            keep.push_back(h);
        }
        leg->Draw("same");

        const float xpos = 0.20, ypos = 0.885, dy = 0.054, fs = 0.040;
        myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),    fs, 0);
        myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(),  fs, 0);
        myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),    fs, 0);
        myText(xpos, ypos - 3 * dy, 1,
               "Flat BDT scan (inclusive PYTHIA MC extraction)", fs, 0);

        c->SaveAs((figDir + "/xsec_flat_scan.pdf").c_str());
        std::cout << "[xsec] wrote " << figDir << "/xsec_flat_scan.pdf ("
                  << keep.size() << "/" << SCAN.size() << ")" << std::endl;
        delete c;
    }

    // ---------- Ratio to nominal (linear y) ----------
    {
        TH1F *hnom = loadXsec("nom", resultsDir, "hxsec_nom_for_ratio");
        if (!hnom)
        {
            std::cerr << "cannot build ratio without nominal — skipping\n";
            return;
        }

        TCanvas *c = new TCanvas("c_xsec_ratio_flat", "", 700, 600);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.14);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.06);

        frame_et_rec->GetYaxis()->SetRangeUser(0.5, 3.0);
        frame_et_rec->SetYTitle("Variant / nominal extracted cross-section");
        frame_et_rec->Draw("axis");

        TLine *unity = new TLine(8, 1.0, 40, 1.0);
        unity->SetLineColor(kGray + 2);
        unity->SetLineStyle(2);
        unity->SetLineWidth(1);
        unity->Draw("same");

        TLegend *leg = new TLegend(0.40, 0.22, 0.92, 0.44);
        legStyle(leg, 0.18, 0.034);

        std::vector<TH1F *> keep;
        for (size_t i = 0; i < SCAN.size(); ++i)
        {
            if (SCAN[i].tag == "nom") continue;   // ratio-to-itself
            TH1F *hv = loadXsec(SCAN[i].tag, resultsDir,
                                Form("hxsec_r_%zu", i));
            if (!hv) continue;
            hv->Divide(hnom);
            hv->SetMarkerColor(SCAN[i].color);
            hv->SetLineColor(SCAN[i].color);
            hv->SetMarkerStyle(SCAN[i].marker);
            hv->SetMarkerSize(1.2);
            hv->SetLineWidth(2);
            hv->Draw("E1 same");
            leg->AddEntry(hv, SCAN[i].label.c_str(), "pl");
            keep.push_back(hv);
        }
        leg->Draw("same");

        const float xpos = 0.20, ypos = 0.885, dy = 0.054, fs = 0.040;
        myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(),    fs, 0);
        myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(),  fs, 0);
        myText(xpos, ypos - 2 * dy, 1, strleg3.c_str(),    fs, 0);
        myText(xpos, ypos - 3 * dy, 1,
               "Flat BDT scan / nominal  (inclusive MC extraction)",
               fs, 0);

        c->SaveAs((figDir + "/xsec_ratio_flat_scan.pdf").c_str());
        std::cout << "[xsec ratio] wrote " << figDir
                  << "/xsec_ratio_flat_scan.pdf" << std::endl;
        delete c;
    }
}
