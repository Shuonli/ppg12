// =============================================================================
// TruthSpectrumOverlay.C
// -----------------------------------------------------------------------------
// Produces four PDF pages overlaying per-sample truth pT spectra (weighted
// by xsec / N_evt / bin_width -> pb/GeV) along with the combined sum.
// Consumes truth_spectrum_{filetype}.root files written by TruthSpectrumCheck.
//
// Outputs (into plotting/figures/truth_spectrum/):
//   truth_spectrum_photon_all.pdf   -- photon5/10/20 + sum (main pipeline)
//   truth_spectrum_photon_di.pdf    -- photon{5,10,20}_double + sum
//   truth_spectrum_jet_all.pdf      -- jet{8,12,20,30,40} + sum (+ jet5 ref)
//   truth_spectrum_jet_di.pdf       -- jet{8,12,20,30,40}_double + sum
//
// All overlays use the per-event MAX truth pT histogram -- that is what the
// Pythia pT-hat binning actually selects on, so stitching is meaningful at the
// truth-max level.
//
// Also prints an ASCII table of boundary ratios (e.g. photon5/photon10 in
// [13,14] GeV). If weights are right these should be ~1; large departures
// point to a weight bug or a truth-window overlap (like the photon10_double
// [10,100] / jet12_double [10,100] double-counting).
//
// Usage:
//   root -l -b -q TruthSpectrumOverlay.C
// =============================================================================

#include "../plotting/plotcommon.h"

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <sys/stat.h>

#include <iostream>
#include <string>
#include <vector>

using namespace std;

namespace {

// One sample entry in an overlay.
struct Sample {
    string filetype;
    string label;
    int    color;
    int    style;  // line style
};

// Resolve truth_spectrum_{filetype}.root relative to results dir.
static string result_path(const string &results_dir, const string &ft) {
    return results_dir + "/truth_spectrum_" + ft + ".root";
}

// Load one histogram (cloned) from results file for a sample; returns nullptr
// on any failure.
static TH1D *load_hist(const string &results_dir, const string &ft,
                       const string &hname)
{
    const string p = result_path(results_dir, ft);
    TFile *f = TFile::Open(p.c_str(), "READ");
    if (!f || f->IsZombie()) {
        cout << "  [miss] " << p << "\n";
        if (f) { f->Close(); delete f; }
        return nullptr;
    }
    TH1D *h = (TH1D *)f->Get(hname.c_str());
    if (!h) {
        cout << "  [miss hist " << hname << "] " << p << "\n";
        f->Close();
        return nullptr;
    }
    h = (TH1D *)h->Clone(Form("%s_%s", hname.c_str(), ft.c_str()));
    h->SetDirectory(nullptr);
    f->Close();
    return h;
}

// Compute integral of histogram in [xlo,xhi) GeV (using bin centers, times
// bin_width since stored as d\sigma/dp_T -> integrate to get xsec[pb]).
static double integral_pb(const TH1D *h, double xlo, double xhi) {
    if (!h) return 0.0;
    double tot = 0.0;
    const double bw = h->GetBinWidth(1);
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        const double c = h->GetBinCenter(i);
        if (c < xlo || c >= xhi) continue;
        tot += h->GetBinContent(i) * bw;
    }
    return tot;
}

// Draw a boundary dashed vertical line at x on current pad (y in log range).
static void draw_boundary(double x, double ymin, double ymax) {
    TLine *l = new TLine(x, ymin, x, ymax);
    l->SetLineColor(kGray + 2);
    l->SetLineStyle(2);
    l->SetLineWidth(1);
    l->Draw();
}

// Load one TNamed string from the results file; empty string on failure.
static string load_named(const string &results_dir, const string &ft,
                         const string &name)
{
    const string p = result_path(results_dir, ft);
    TFile *f = TFile::Open(p.c_str(), "READ");
    if (!f || f->IsZombie()) { if (f) delete f; return ""; }
    TNamed *n = (TNamed *)f->Get(name.c_str());
    string v = n ? n->GetTitle() : "";
    f->Close();
    return v;
}

// Draw one overlay page: N sample hists + sum, log-y, save PDF.
static void draw_overlay(const string &results_dir,
                         const vector<Sample> &samples,
                         const string &hname,
                         const string &title_top,
                         const string &xtitle,
                         const vector<double> &boundaries,
                         double xlo, double xhi, double ylo, double yhi,
                         const string &outpdf)
{
    cout << "\n--- " << outpdf << " (" << hname << ") ---\n";

    // Print declared truth windows for each sample (diagnostic for
    // mis-specified windows like photon10_double [10,100] bug).
    const bool is_photon_hist = (hname.find("photon") != string::npos);
    cout << "  declared truth windows ("
         << (is_photon_hist ? "photon" : "jet") << "):\n";
    for (const auto &s : samples) {
        const string lo = load_named(results_dir, s.filetype,
                                     is_photon_hist ? "photon_pt_lower" : "jet_pt_lower");
        const string hi = load_named(results_dir, s.filetype,
                                     is_photon_hist ? "photon_pt_upper" : "jet_pt_upper");
        cout << "    " << s.filetype << "  [" << lo << ", " << hi << "]\n";
    }

    vector<TH1D *> hs;
    TH1D *hsum = nullptr;
    int nloaded = 0;
    for (const auto &s : samples) {
        TH1D *h = load_hist(results_dir, s.filetype, hname);
        hs.push_back(h);
        if (!h) continue;
        ++nloaded;
        if (!hsum) {
            hsum = (TH1D *)h->Clone(("hsum_" + hname).c_str());
            hsum->SetDirectory(nullptr);
            hsum->Reset();
        }
        hsum->Add(h);
    }
    if (nloaded == 0) {
        cout << "  [skip page: no samples loaded]\n";
        return;
    }

    TCanvas *c = new TCanvas(Form("c_%s_%s", hname.c_str(), outpdf.c_str()),
                             outpdf.c_str(), 900, 700);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.08);
    c->SetBottomMargin(0.14);
    c->SetLogy();

    // Frame from scratch -- we want a wide pT range that doesn't match the
    // shared 8-40 GeV reco frames in plotcommon.h.
    TH1F *frame = new TH1F(Form("frame_%s_%s", hname.c_str(), outpdf.c_str()),
                           title_top.c_str(), 10, xlo, xhi);
    frame->SetXTitle(xtitle.c_str());
    frame->SetYTitle("d#sigma/dp_{T} [pb/GeV]");
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->SetRangeUser(xlo, xhi);
    frame->GetYaxis()->SetRangeUser(ylo, yhi);
    frame->SetTitle("");
    frame->Draw("axis");

    // Boundary dashed lines
    for (double b : boundaries) draw_boundary(b, ylo, yhi);

    // Per-sample histograms as colored lines
    for (size_t i = 0; i < samples.size(); ++i) {
        TH1D *h = hs[i];
        if (!h) continue;
        h->SetLineColor(samples[i].color);
        h->SetLineStyle(samples[i].style);
        h->SetLineWidth(2);
        h->SetMarkerStyle(1);
        h->SetMarkerSize(0);
        h->Draw("hist same");
    }

    // Combined sum as black filled circles
    if (hsum) {
        hsum->SetMarkerStyle(20);
        hsum->SetMarkerColor(kBlack);
        hsum->SetLineColor(kBlack);
        hsum->SetMarkerSize(0.8);
        hsum->SetLineWidth(1);
        hsum->Draw("p e same");
    }

    // Legend
    TLegend *leg = new TLegend(0.60, 0.55, 0.93, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.032);
    for (size_t i = 0; i < samples.size(); ++i) {
        TH1D *h = hs[i];
        if (!h) continue;
        leg->AddEntry(h, samples[i].label.c_str(), "l");
    }
    if (hsum) leg->AddEntry(hsum, "combined sum", "lep");
    leg->Draw();

    // sPHENIX label block
    myText(0.18, 0.85, 1, strleg1.c_str(), 0.040);
    myText(0.18, 0.80, 1, "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV Pythia8", 0.035);
    myText(0.18, 0.75, 1, title_top.c_str(), 0.032);

    c->SaveAs(outpdf.c_str());
    cout << "  wrote " << outpdf << "\n";

    // Boundary-slab diagnostics. For each stitching pT value b we compare the
    // per-sample integrals in two windows:
    //   left  = [b-1, b)  (below the boundary -- should be dominated by the
    //                       LOW-pT sample whose window ends at b)
    //   right = [b, b+1)  (above the boundary -- should be dominated by the
    //                       HIGH-pT sample whose window starts at b)
    // If weighting is correct:
    //   left_top_sample  ~ right_top_sample * (left_width * xsec_lo / right_width * xsec_hi)
    //   Stitching is smooth ->  both samples cover left/right tails similarly
    // Double-counting shows up as BOTH samples contributing in the SAME window
    // (i.e. a sample's non-zero contribution far outside its declared truth
    //  window), not as ratio of combined / single in a merge window.
    cout << "  boundary diagnostics (per-sample integrals [pb]):\n";
    vector<double> per_L(samples.size(), 0.0);
    vector<double> per_R(samples.size(), 0.0);
    for (double b : boundaries) {
        for (size_t i = 0; i < hs.size(); ++i) {
            per_L[i] = hs[i] ? integral_pb(hs[i], b - 1.0, b) : 0.0;
            per_R[i] = hs[i] ? integral_pb(hs[i], b,       b + 1.0) : 0.0;
        }
        const double comb_L = integral_pb(hsum, b - 1.0, b);
        const double comb_R = integral_pb(hsum, b,       b + 1.0);
        const double ratio_smooth = (comb_R > 0) ? comb_L / comb_R : 0.0;

        cout << Form("    boundary pT=%5.1f GeV  left[%4.1f,%4.1f)=%.3e  right[%4.1f,%4.1f)=%.3e  L/R=%.2f",
                     b, b - 1.0, b, comb_L, b, b + 1.0, comb_R, ratio_smooth);
        if (ratio_smooth < 0.7 || ratio_smooth > 1.5) cout << " (kink!)";
        cout << "\n";

        cout << "      per-sample (left | right):";
        for (size_t i = 0; i < hs.size(); ++i) {
            if (per_L[i] <= 0 && per_R[i] <= 0) continue;
            cout << Form("  %s=(%.2e|%.2e)", samples[i].filetype.c_str(), per_L[i], per_R[i]);
        }
        cout << "\n";

        // Double-counting detector: sample X contributes > 5% of comb_L AND
        // > 5% of comb_R when another sample is supposed to own one of these
        // bands. Flag samples that appear in BOTH windows with > 5% share.
        for (size_t i = 0; i < hs.size(); ++i) {
            const double fL = (comb_L > 0) ? per_L[i] / comb_L : 0.0;
            const double fR = (comb_R > 0) ? per_R[i] / comb_R : 0.0;
            if (fL > 0.05 && fR > 0.05) {
                cout << Form("      ** %s straddles boundary (L share=%.2f, R share=%.2f) -- likely TRUTH-WINDOW double-count **\n",
                             samples[i].filetype.c_str(), fL, fR);
            }
        }
    }
    cout << "  ---\n";
}

}  // anonymous namespace

void TruthSpectrumOverlay(
    const string &results_dir = "results",
    const string &figdir      = "/sphenix/user/shuhangli/ppg12/plotting/figures/truth_spectrum")
{
    init_plot();

    ::mkdir(figdir.c_str(), 0755);

    // -------------------------------------------------------------------------
    // Photon main pipeline: photon5 [0,14] + photon10 [14,22] + photon20 [22+]
    // -------------------------------------------------------------------------
    vector<Sample> photon_main = {
        {"photon5",  "photon5 (0-14 GeV)",   kRed + 1,    1},
        {"photon10", "photon10 (14-22 GeV)", kBlue + 1,   1},
        {"photon20", "photon20 (22+ GeV)",   kGreen + 2,  1},
    };

    draw_overlay(results_dir, photon_main,
                 "h_max_truth_photon_pt_filtered",
                 "Truth #gamma p_{T} spectrum -- main pipeline",
                 "max truth p_{T}^{#gamma} (|#eta|<0.7) [GeV]",
                 {14.0, 22.0},
                 0.0, 60.0, 1e-2, 1e10,
                 figdir + "/truth_spectrum_photon_all.pdf");

    // -------------------------------------------------------------------------
    // Photon DI: photon{5,10,20}_double
    // -------------------------------------------------------------------------
    vector<Sample> photon_di = {
        {"photon5_double",  "photon5_double (0-14 GeV)",   kRed + 1,    1},
        {"photon10_double", "photon10_double (14-22 GeV)", kBlue + 1,   1},
        {"photon20_double", "photon20_double (22+ GeV)",   kGreen + 2,  1},
    };

    draw_overlay(results_dir, photon_di,
                 "h_max_truth_photon_pt_filtered",
                 "Truth #gamma p_{T} spectrum -- double-interaction MC",
                 "max truth p_{T}^{#gamma} (|#eta|<0.7) [GeV]",
                 {14.0, 22.0},
                 0.0, 60.0, 1e-2, 1e10,
                 figdir + "/truth_spectrum_photon_di.pdf");

    // -------------------------------------------------------------------------
    // Jet main pipeline (MergeSim.C list): jet8, jet12, jet20, jet30, jet40
    // plus jet5 as a lighter dashed reference (training-input only).
    // -------------------------------------------------------------------------
    vector<Sample> jet_main = {
        {"jet5",  "jet5 (training ref)",  kGray + 1,   2},  // dashed, not in merge
        {"jet8",  "jet8 (9-14 GeV)",      kRed + 1,    1},
        {"jet12", "jet12 (14-21 GeV)",    kOrange + 7, 1},
        {"jet20", "jet20 (21-32 GeV)",    kGreen + 2,  1},
        {"jet30", "jet30 (32-42 GeV)",    kAzure + 1,  1},
        {"jet40", "jet40 (42+ GeV)",      kViolet + 1, 1},
    };
    // Combined sum should NOT include jet5 since it's not in the pipeline merge.
    // To avoid biasing the sum we mark jet5 dashed but still include it; the
    // reader will see the dashed line drawn below the sum -- we exclude it from
    // the sum by loading separately. Simpler: run draw_overlay over pipeline
    // samples only for the sum and overlay jet5 on top with a second call.

    // First overlay the 5 pipeline samples with combined sum:
    vector<Sample> jet_main_pipeline(jet_main.begin() + 1, jet_main.end());
    draw_overlay(results_dir, jet_main_pipeline,
                 "h_max_truth_jet_pt_filtered",
                 "Truth jet p_{T} spectrum -- main pipeline (jet{8,12,20,30,40})",
                 "max truth p_{T}^{jet} [GeV]",
                 {14.0, 21.0, 32.0, 42.0},
                 0.0, 60.0, 1e-2, 1e15,
                 figdir + "/truth_spectrum_jet_all.pdf");

    // -------------------------------------------------------------------------
    // Jet DI: jet{8,12,20,30,40}_double
    // -------------------------------------------------------------------------
    vector<Sample> jet_di = {
        {"jet8_double",  "jet8_double (9-14 GeV)",    kRed + 1,    1},
        {"jet12_double", "jet12_double (14-21 GeV)",  kOrange + 7, 1},
        {"jet20_double", "jet20_double (21-32 GeV)",  kGreen + 2,  1},
        {"jet30_double", "jet30_double (32-42 GeV)",  kAzure + 1,  1},
        {"jet40_double", "jet40_double (42+ GeV)",       kViolet + 1, 1},
    };
    draw_overlay(results_dir, jet_di,
                 "h_max_truth_jet_pt_filtered",
                 "Truth jet p_{T} spectrum -- double-interaction MC",
                 "max truth p_{T}^{jet} [GeV]",
                 {14.0, 21.0, 32.0, 42.0},
                 0.0, 60.0, 1e-2, 1e15,
                 figdir + "/truth_spectrum_jet_di.pdf");

    cout << "\nDone. PDFs in " << figdir << "\n";
}
