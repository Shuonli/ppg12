// plot_paper_final.C -- paper Fig.~\ref{fig:final}, Fig.~\ref{fig:final_phenix},
// Fig.~\ref{fig:final_sphenix}.
//
// Lifted verbatim from plotting/plot_final_selection.C with three changes:
//   1. Function renamed plot_final_selection -> plot_paper_final, and the
//      output paths redirected to PPG12-Paper/figures/<canonical name>.pdf.
//   2. The intermediate diagnostic canvases (final_common_cluster,
//      final_tight_iso_cluster, final_all, final_phenix_fit) are kept in
//      this file because their data products (h_*_data, h_unfold_ratio)
//      feed the cross-section panel; only their SaveAs lines are redirected
//      to a PPG12-Paper/figures/_unused/ scratch directory so the paper-
//      figure folder stays uncluttered.
//   3. A new "sPHENIX-only" canvas is appended at the bottom (lifted and
//      modernized from plot_final_backup250327_v1.C, L877-933) producing
//      final_sphenix.pdf -- our data + total systematic band + lumi
//      legend, no PHENIX overlay.
//
// Cosmetic iteration on the paper figures lives entirely in this file.
// Edit colours, fonts, legends, ranges, axis labels here without touching
// the parent plotting/plot_final_selection.C (which feeds the analysis
// note + comparison reports).

#include "paper_style.h"
#include <yaml-cpp/yaml.h>
#include <TSystem.h>

namespace {
// 0: our data, 1: PHENIX, 2: JETPHOX, 3: PYTHIA, 4: Werner
const int col[] = {kAzure + 2, kPink + 5, kPink + 5, kOrange + 7, kYellow + 2, kSpring - 7,kCyan + 4, kRed - 4, kBlack, kBlue - 3, kPink - 5, kGreen + 3, kBlue - 3};
const int mkcol[] = {kAzure + 2, kPink + 5, kPink + 5, kOrange + 7, kTeal+4, kSpring +4, kSpring - 7,kCyan + 4, kRed - 4, kBlack, kBlue - 3, kPink - 5, kGreen + 3, kBlue - 3};
const float trans[] = {0.35, 0.5, 0.40, 1.0, 0.0};
const float boxlinewidth[] = {0., 0., 0., 0.0, 1.0};
const int mkStyle[] = {20, 21, 25, 47, 27, 33, 25, 27, 28, 24, 29, 28, 22};
const float mkSize[] = {1.4, 1.4, 1.2, 1.5, 2.0, 1, 1, 1, 1, 1, 1, 1, 1};
const int lineWidth[] = {2, 2, 2, 2, 2};
const int fillStyle[] = {1, 1, 1001, 1, 0};
// const int fillStyle[] = {1, 1, 1001, 1, 1001};
// 0: CT18, 1: NNPDF, 2: CTEQ, 3: MSHT, 4: CT14
const float mkSizepdf[] = {1.4, 1.4, 1.7, 1.5, 1.0, 1.0, 1.0};
const float mkcolpdf[] = {kPink + 5, kRed - 4, kSpring -7, kBlue+1, kOrange + 7, kAzure - 4};
}  // namespace

void plot_paper_final_yj(string tune = "bdt_nom")
{
    paper_init();

    // Output paths. Stored as locals (not as repeated paper_savepath().c_str()
    // calls inline) because Cling occasionally hands back a dangling/stale
    // c_str() when the function-local static is read mid-function -- both
    // paths and the canvas titles share the same overflow string buffer.
    const std::string paperdir = paper_savepath();
    const std::string scratch  = paperdir + "/_unused";
    gSystem->mkdir(scratch.c_str(), kTRUE);

    // Read luminosity from the matching analysis config
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    std::string configpath = Form("/sphenix/user/shuhangli/ppg12/efficiencytool/config_%s.yaml", tune.data());
    float datalumi = 49.562; // fallback
    try {
        YAML::Node cfg = YAML::LoadFile(configpath);
        datalumi = cfg["analysis"]["lumi"].as<float>();
        std::cout << "[plot] Loaded lumi = " << datalumi << " pb^-1 from " << configpath << std::endl;
    } catch (...) {
        std::cerr << "[plot] WARNING: could not read lumi from " << configpath
                  << ", using fallback " << datalumi << " pb^-1" << std::endl;
    }

    // Build luminosity legend strings from the loaded value (#1: fixed broken
    // 2-line lumi label — leading-space hack collapsed and made "= 64.4 pb^{-1}"
    // look like a syntactically broken assignment).
    std::string strleg_lumi = Form("#it{p}+#it{p} #kern[-0.05]{#sqrt{#it{s}} = 200 GeV, %.1f pb^{-1}}", datalumi);
    std::string strleg_lumi_line2 = Form("#it{L} = %.1f pb^{-1}", datalumi);
    // std::string strleg_lumi_line2 = Form("#int #it{L} d#it{t} = %.1f pb^{-1}", datalumi);

    std::string MCstring = "JETPHOX";

    std::string reweightedstring = "";

    std::string datastring = "data";
    std::string bg_MCstring = "Inclusive Sim";

    float deta = 1.4;

    float lowery = 0.2;
    float lowerx = 12;
    float upperx = 32;

    TFile *fin_data = new TFile(Form("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_%s.root", tune.data()));

    TFile *fin_syst = new TFile("/sphenix/user/shuhangli/ppg12/plotting/rootFiles/syst_sum.root");
    TFile *fin_NLO = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_ct18_10_chunked.root");
    TFile *fin_NLO_up = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_ct18_05_chunked.root");
    TFile *fin_NLO_down = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_ct18_20_chunked.root");
    // Alternative-PDF central-scale JETPHOX productions for the PDF
    // comparison panel (Lower-most pad of the final figure).
    TFile *fin_NLO_ct14  = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_nlo_10_chunked.root");
    TFile *fin_NLO_nnpdf = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_nnpdf4_10_chunked.root");
    TFile *fin_NLO_cteq  = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_cteq_10_chunked.root");
    TFile *fin_NLO_msht  = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_msht_10_chunked.root");
    // Scale variations for the additional PDFs (used by the bottom-pad
    // per-PDF systematic bands).
    TFile *fin_NLO_nnpdf_up = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_nnpdf_05_chunked.root");
    TFile *fin_NLO_nnpdf_dn = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_nnpdf_20_chunked.root");
    TFile *fin_NLO_cteq_up  = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_cteq_05_chunked.root");
    TFile *fin_NLO_cteq_dn  = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_cteq_20_chunked.root");
    TFile *fin_NLO_msht_up  = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_msht_05_chunked.root");
    TFile *fin_NLO_msht_dn  = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_msht_20_chunked.root");

    // Per-PDF Hessian uncertainty bands (used in the bottom pad).
    TFile *fin_pdfunc = new TFile("/sphenix/user/shuhangli/ppg12/NLO/rootFiles/jetPHOX_pdfunc.root");
    TGraphAsymmErrors *g_pdf_band_ct18_abs   = (TGraphAsymmErrors*)fin_pdfunc->Get("g_pdf_band_ct18");
    TGraphAsymmErrors *g_pdf_band_nnpdf4_abs = (TGraphAsymmErrors*)fin_pdfunc->Get("g_pdf_band_nnpdf4");
    TGraphAsymmErrors *g_pdf_band_cteq_abs   = (TGraphAsymmErrors*)fin_pdfunc->Get("g_pdf_band_cteq");
    TGraphAsymmErrors *g_pdf_band_msht_abs   = (TGraphAsymmErrors*)fin_pdfunc->Get("g_pdf_band_msht");

    TFile *fin_mc = new TFile(Form("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_%s_mc.root", tune.data()));

    TH1F *h_data = (TH1F *)fin_data->Get("h_unfold_sub_result");
    h_data->Scale(1.0 / deta);
    TH1F *h_data_cp = (TH1F *)h_data->Clone("h_data_cp");
    TH1F *h_pythia = (TH1F *)fin_data->Get("h_truth_pT_novtx_0");
    h_pythia->Scale(1.0 / deta);
    TH1F *h_common_cluster_data = (TH1F *)fin_data->Get("h_common_cluster_0");
    h_common_cluster_data->Scale(1.0 / deta);
    TH1F *h_common_cluster_mc = (TH1F *)fin_mc->Get("h_common_cluster_0");
    h_common_cluster_mc->Scale(1.0 / deta);
    TH1F *h_tight_iso_cluster_data = (TH1F *)fin_data->Get("h_tight_iso_cluster_0");
    h_tight_iso_cluster_data->Scale(1.0 / deta);
    TH1F *h_sub_data = (TH1F *)fin_data->Get("h_data_sub_leak");
    h_sub_data->Scale(1.0 / deta);
    TH1F *h_sub_data_unfold = (TH1F *)fin_data->Get("h_unfold_sub_result_woeff");
    h_sub_data_unfold->Scale(1.0 / deta);
    TH1F *h_tight_iso_cluster_mc = (TH1F *)fin_mc->Get("h_tight_iso_cluster_0");
    h_tight_iso_cluster_mc->Scale(1.0 / deta);

    TH1F *h_NLO = (TH1F *)fin_NLO->Get("h_truth_pT");
    h_NLO->Scale(1.0 / deta);
    TH1F *h_NLO_up = (TH1F *)fin_NLO_up->Get("h_truth_pT");
    h_NLO_up->Scale(1.0 / deta);
    TH1F *h_NLO_down = (TH1F *)fin_NLO_down->Get("h_truth_pT");
    h_NLO_down->Scale(1.0 / deta);

    TH1F *h_data_NLO = (TH1F *)h_data->Clone("h_data_NLO");
    h_data_NLO->Divide(h_NLO);

    TH1F *h_NLO_data = (TH1F *)h_NLO->Clone("h_NLO_data");
    h_NLO_data->Divide(h_data);
    TH1F *h_NLO_data_up = (TH1F *)h_NLO_up->Clone("h_NLO_data_up");
    h_NLO_data_up->Divide(h_data);
    TH1F *h_NLO_data_down = (TH1F *)h_NLO_down->Clone("h_NLO_data_down");
    h_NLO_data_down->Divide(h_data);

    // Alternative-PDF JETPHOX/data ratios for the PDF-comparison panel.
    TH1F *h_NLO_ct14 = (TH1F *)fin_NLO_ct14->Get("h_truth_pT");
    h_NLO_ct14->Scale(1.0 / deta);
    TH1F *h_NLO_nnpdf = (TH1F *)fin_NLO_nnpdf->Get("h_truth_pT");
    h_NLO_nnpdf->Scale(1.0 / deta);
    TH1F *h_NLO_data_ct14  = (TH1F *)h_NLO_ct14 ->Clone("h_NLO_data_ct14");
    h_NLO_data_ct14 ->Divide(h_data);
    TH1F *h_NLO_data_nnpdf = (TH1F *)h_NLO_nnpdf->Clone("h_NLO_data_nnpdf");
    h_NLO_data_nnpdf->Divide(h_data);
    TH1F *h_NLO_cteq = (TH1F *)fin_NLO_cteq->Get("h_truth_pT");
    h_NLO_cteq->Scale(1.0 / deta);
    TH1F *h_NLO_msht = (TH1F *)fin_NLO_msht->Get("h_truth_pT");
    h_NLO_msht->Scale(1.0 / deta);
    TH1F *h_NLO_data_cteq = (TH1F *)h_NLO_cteq->Clone("h_NLO_data_cteq");
    h_NLO_data_cteq->Divide(h_data);
    TH1F *h_NLO_data_msht = (TH1F *)h_NLO_msht->Clone("h_NLO_data_msht");
    h_NLO_data_msht->Divide(h_data);

    // Build per-PDF Hessian uncertainty bands in JETPHOX/Data form. The
    // bands in fin_pdfunc are cross-section absolute values; the relative
    // uncertainty (eyhigh/y, eylow/y) is applied to each PDF's nominal
    // ratio (h_NLO_data_<pdf>) so the band sits around each central line
    // in the bottom-pad ratio plot.
    // Horizontal errors set to 0 so the filled band edges interpolate
    // between bin centers, matching the central-line draws (which also
    // pass through bin centers). Restricted to the reported analysis
    // range 12 < ETg < 32 GeV so the band starts/ends at the same bin
    // centers as the central lines.
    auto build_pdf_band_ratio = [&](TGraphAsymmErrors *g_band_abs, TH1F *h_ratio) {
        TGraphAsymmErrors *g_out = new TGraphAsymmErrors();
        int n = g_band_abs->GetN();
        for (int i = 0; i < n; i++) {
            double x = g_band_abs->GetX()[i];
            if (x < 12.0 || x > 32.0) continue;
            double y_abs = g_band_abs->GetY()[i];
            if (y_abs <= 0) continue;
            double eyhi_rel  = g_band_abs->GetErrorYhigh(i) / y_abs;
            double eylow_rel = g_band_abs->GetErrorYlow(i)  / y_abs;
            int bin = h_ratio->FindBin(x);
            double ratio = h_ratio->GetBinContent(bin);
            if (ratio <= 0) continue;
            int p = g_out->GetN();
            g_out->SetPoint(p, x, ratio);
            g_out->SetPointError(p, 0.0, 0.0, ratio * eylow_rel, ratio * eyhi_rel);
        }
        return g_out;
    };

    TGraphAsymmErrors *g_pdfband_ct18  = build_pdf_band_ratio(g_pdf_band_ct18_abs,   h_NLO_data);
    TGraphAsymmErrors *g_pdfband_nnpdf = build_pdf_band_ratio(g_pdf_band_nnpdf4_abs, h_NLO_data_nnpdf);
    TGraphAsymmErrors *g_pdfband_cteq  = build_pdf_band_ratio(g_pdf_band_cteq_abs,   h_NLO_data_cteq);
    TGraphAsymmErrors *g_pdfband_msht  = build_pdf_band_ratio(g_pdf_band_msht_abs,   h_NLO_data_msht);

    // Up/down scale variations -> theory/data ratios for NNPDF/CTEQ/MSHT.
    auto load_NLO = [&](TFile *f) {
        TH1F *h = (TH1F *)f->Get("h_truth_pT");
        h->Scale(1.0 / deta);
        return h;
    };
    TH1F *h_NLO_nnpdf_up = load_NLO(fin_NLO_nnpdf_up);
    TH1F *h_NLO_nnpdf_dn = load_NLO(fin_NLO_nnpdf_dn);
    TH1F *h_NLO_cteq_up  = load_NLO(fin_NLO_cteq_up);
    TH1F *h_NLO_cteq_dn  = load_NLO(fin_NLO_cteq_dn);
    TH1F *h_NLO_msht_up  = load_NLO(fin_NLO_msht_up);
    TH1F *h_NLO_msht_dn  = load_NLO(fin_NLO_msht_dn);

    TH1F *h_NLO_data_nnpdf_up = (TH1F *)h_NLO_nnpdf_up->Clone("h_NLO_data_nnpdf_up");
    h_NLO_data_nnpdf_up->Divide(h_data);
    TH1F *h_NLO_data_nnpdf_dn = (TH1F *)h_NLO_nnpdf_dn->Clone("h_NLO_data_nnpdf_dn");
    h_NLO_data_nnpdf_dn->Divide(h_data);
    TH1F *h_NLO_data_cteq_up = (TH1F *)h_NLO_cteq_up->Clone("h_NLO_data_cteq_up");
    h_NLO_data_cteq_up->Divide(h_data);
    TH1F *h_NLO_data_cteq_dn = (TH1F *)h_NLO_cteq_dn->Clone("h_NLO_data_cteq_dn");
    h_NLO_data_cteq_dn->Divide(h_data);
    TH1F *h_NLO_data_msht_up = (TH1F *)h_NLO_msht_up->Clone("h_NLO_data_msht_up");
    h_NLO_data_msht_up->Divide(h_data);
    TH1F *h_NLO_data_msht_dn = (TH1F *)h_NLO_msht_dn->Clone("h_NLO_data_msht_dn");
    h_NLO_data_msht_dn->Divide(h_data);

    TGraphAsymmErrors *g_syst_rel_nnpdf = new TGraphAsymmErrors(h_NLO_data_nnpdf);
    TGraphAsymmErrors *g_syst_rel_cteq  = new TGraphAsymmErrors(h_NLO_data_cteq);
    TGraphAsymmErrors *g_syst_rel_msht  = new TGraphAsymmErrors(h_NLO_data_msht);
    for (int i = 0; i < h_NLO_nnpdf->GetNbinsX(); i++)
    {
        float xlow = g_syst_rel_nnpdf->GetErrorXlow(i);
        float xup  = g_syst_rel_nnpdf->GetErrorXhigh(i);
        g_syst_rel_nnpdf->SetPointError(i, xlow, xup,
            h_NLO_data_nnpdf->GetBinContent(i + 1) - h_NLO_data_nnpdf_dn->GetBinContent(i + 1),
            h_NLO_data_nnpdf_up->GetBinContent(i + 1) - h_NLO_data_nnpdf->GetBinContent(i + 1));
        g_syst_rel_cteq->SetPointError(i, xlow, xup,
            h_NLO_data_cteq->GetBinContent(i + 1) - h_NLO_data_cteq_dn->GetBinContent(i + 1),
            h_NLO_data_cteq_up->GetBinContent(i + 1) - h_NLO_data_cteq->GetBinContent(i + 1));
        g_syst_rel_msht->SetPointError(i, xlow, xup,
            h_NLO_data_msht->GetBinContent(i + 1) - h_NLO_data_msht_dn->GetBinContent(i + 1),
            h_NLO_data_msht_up->GetBinContent(i + 1) - h_NLO_data_msht->GetBinContent(i + 1));
    }

    TH1F *h_pythia_data = (TH1F *)h_pythia->Clone("h_pythia_data");
    h_pythia_data->Divide(h_data);

    TH1F *h_syst_low = (TH1F *)fin_syst->Get("h_sum_low");
    h_syst_low->Scale(1.0 / deta);
    TH1F *h_syst_high = (TH1F *)fin_syst->Get("h_sum_high");
    h_syst_high->Scale(1.0 / deta);

    TH1F *h_syst_rel_low = (TH1F *)fin_syst->Get("h_sum_rel_low");
    TH1F *h_syst_rel_high = (TH1F *)fin_syst->Get("h_sum_rel_high");

    TGraphAsymmErrors *g_syst = new TGraphAsymmErrors(h_data);
    TGraphAsymmErrors *g_syst_rel_pythia = new TGraphAsymmErrors(h_pythia_data);
    TGraphAsymmErrors *g_syst_rel = new TGraphAsymmErrors(h_data);

    for (int i = 0; i < h_data->GetNbinsX(); i++)
    {
        float xlowerror = g_syst->GetErrorXlow(i);
        float xuperror = g_syst->GetErrorXhigh(i);

        g_syst->SetPointError(i, xlowerror, xuperror, h_syst_low->GetBinContent(i + 1), h_syst_high->GetBinContent(i + 1));

        // float y = g_syst_rel_NLO->GetY()[i];
        // std::cout << "y: " << y << std::endl;
        // why this is here?

        // g_syst_rel_NLO->SetPointError(i, xlowerror, xuperror, y * h_syst_rel_high->GetBinContent(i + 1), y * h_syst_rel_low->GetBinContent(i + 1));
        // std::cout<<h_data->GetBinCenter(i+1)<<std::endl;
        // std::cout << "h_syst_rel_low->GetBinContent(i + 1): " << h_syst_rel_low->GetBinContent(i + 1) << std::endl;
        // std::cout << "h_syst_rel_high->GetBinContent(i + 1): " << h_syst_rel_high->GetBinContent(i + 1) << std::endl;
        // float y_pythia = g_syst_rel_pythia->GetY()[i];

        std::cout << "xlowerror" << xlowerror << std::endl;
        // g_syst_rel_pythia->SetPointError(i, xlowerror, xuperror, y_pythia * h_syst_rel_high->GetBinContent(i + 1), y_pythia * h_syst_rel_low->GetBinContent(i + 1));

        // float x_data = h_data->GetBinContent(i + 1);
        g_syst_rel->SetPoint(i, g_syst_rel->GetX()[i], 1);
        g_syst_rel->SetPointError(i, xlowerror, xuperror, abs(h_syst_rel_low->GetBinContent(i + 1)), h_syst_rel_high->GetBinContent(i + 1));
    }

    TGraphAsymmErrors *g_syst_NLO = new TGraphAsymmErrors(h_NLO);
    TGraphAsymmErrors *g_syst_rel_NLO = new TGraphAsymmErrors(h_NLO_data);

    for (int i = 0; i < h_NLO->GetNbinsX(); i++)
    {
        float xlowerror = g_syst_NLO->GetErrorXlow(i);
        float xuperror = g_syst_NLO->GetErrorXhigh(i);

        g_syst_NLO->SetPointError(i, xlowerror, xuperror, h_NLO->GetBinContent(i + 1) - h_NLO_down->GetBinContent(i + 1), h_NLO_up->GetBinContent(i + 1) - h_NLO->GetBinContent(i + 1));

        g_syst_rel_NLO->SetPointError(i, xlowerror, xuperror, h_NLO_data->GetBinContent(i + 1) - h_NLO_data_down->GetBinContent(i + 1), h_NLO_data_up->GetBinContent(i + 1) - h_NLO_data->GetBinContent(i + 1));
    }

    TCanvas *c1 = new TCanvas("can", "", 800, 1110);
    c1->Divide(1, 3);

    TPad *pad_1 = (TPad *)c1->cd(1);
    pad_1->SetPad(0, 0.5, 1, 1);
    pad_1->SetTopMargin(0.05);
    pad_1->SetLeftMargin(0.13);
    pad_1->SetBottomMargin(0.002);
    pad_1->SetRightMargin(0.08);
    pad_1->SetLogy();

    frame_et_rec->SetYTitle("d^{2}#it{#sigma}/(d#it{#eta} d#it{E}_{T}^{#gamma}) [pb/GeV]");
    // frame_et_rec->SetYTitle("d#sigma/d#eta/dE_{T} [pb/GeV]");
    frame_et_rec->GetYaxis()->SetRangeUser(lowery, 500);
    // frame_et_rec->GetYaxis()->SetRangeUser(lowery, 1500);
    // frame_et_rec->GetYaxis()->SetRangeUser(0.2, 4e3);
    frame_et_rec->GetXaxis()->SetRangeUser(lowerx, upperx);

    frame_et_rec->GetXaxis()->SetTitleOffset(1.05);
    frame_et_rec->GetYaxis()->SetTitleOffset(1.05);
    frame_et_rec->GetYaxis()->SetTitleSize(0.053);
    frame_et_rec->GetXaxis()->SetLabelSize(0.05);  // top-panel x-labels hidden (ratio panel below has them) — #6
    // frame_et_rec->GetXaxis()->SetLabelSize(0);  // top-panel x-labels hidden (ratio panel below has them) — #6
    frame_et_rec->GetYaxis()->SetLabelSize(0.050);
    // frame_et_rec->GetXaxis()->CenterTitle();
    // frame_et_rec->GetYaxis()->CenterTitle();
    frame_et_rec->GetXaxis()->SetNdivisions(505);

    frame_et_rec->Draw("axis");

    // PHENIX data
    //  Number of data points

    const Int_t n = 18;

    // Central pT
    Double_t x[n] = {
        5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75,
        11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0};

    // pT low edges
    Double_t xLow[n] = {
        5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
        10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0};

    // pT high edges
    Double_t xHigh[n] = {
        5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
        12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0};

    // Cross section values, E d^3sigma/dp^3 [pb * GeV^-2 * c^3]
    Double_t y[n] = {
        1140, 613, 348, 231, 136, 93, 67, 48.3, 32.1, 20.4,
        9.8, 2.97, 1.06, 0.34, 0.173, 0.088, 0.042, 0.029};

    // Stat errors (asymmetric: + and -)
    Double_t statUp[n] = {
        30, 19, 13, 9, 6, 4, 3.2, 2.5, 1.9, 1.5,
        0.4, 0.19, 0.1, 0.06, 0.034, 0.02, 0.014, 0.009};

    Double_t statDown[n] = {
        30, 19, 13, 9, 6, 4, 3.2, 2.5, 1.9, 1.5,
        0.4, 0.19, 0.1, 0.06, 0.034, 0.021, 0.015, 0.014};

    // Sys errors (asymmetric: + and -)
    Double_t sysUp[n] = {
        478, 221, 101, 62, 31, 20, 13.4, 9.2, 6.1, 3.7,
        1.7, 0.48, 0.17, 0.05, 0.028, 0.015, 0.007, 0.004};

    Double_t sysDown[n] = {
        478, 221, 101, 62, 31, 20, 13.4, 9.2, 6.1, 3.7,
        1.7, 0.48, 0.17, 0.05, 0.028, 0.015, 0.007, 0.004};

    // Create the TGraphAsymmErrors for STAT
    TGraphAsymmErrors *gStat_PHENIX = new TGraphAsymmErrors(n);
    gStat_PHENIX->SetName("gStat");
    gStat_PHENIX->SetTitle("Cross Section (scaled) with Stat. Errors");

    // Create the TGraphAsymmErrors for SYS
    TGraphAsymmErrors *gSys_PHENIX = new TGraphAsymmErrors(n);
    gSys_PHENIX->SetName("gSys");
    gSys_PHENIX->SetTitle("Cross Section (scaled) with Syst. Errors");

    // Common factor that does NOT include pT, since pT changes per bin
    // Overall factor: 2*pi * 2 * deta, with deta=0.7
    // Double_t factorCommon = 2.0 * TMath::Pi() * 2.0 * 0.7; // ~ 8.7964
    Double_t factorCommon = 2.0 * TMath::Pi();

    // Cosmetic-only x half-width for the published PHENIX sys box; small
    // enough not to suggest a bin-width meaning (PHENIX has bin-shifted
    // these to the bin centre themselves), large enough that the box
    // renders. Smaller than every visible PHENIX bin's half-width.
    const Double_t kPhenixSysHalfX = 0.5;
    for (Int_t i = 0; i < n; i++)
    {
        Double_t scaleFactor = x[i] * factorCommon;
        Double_t yScaled = y[i] * scaleFactor;
        Double_t statUpScaled = statUp[i] * scaleFactor;
        Double_t statDownScaled = statDown[i] * scaleFactor;
        Double_t sysUpScaled = sysUp[i] * scaleFactor;
        Double_t sysDownScaled = sysDown[i] * scaleFactor;

        gStat_PHENIX->SetPoint(i, x[i], yScaled);
        gStat_PHENIX->SetPointError(i, 0.0, 0.0, statDownScaled, statUpScaled);

        gSys_PHENIX->SetPoint(i, x[i], yScaled);
        gSys_PHENIX->SetPointError(i, kPhenixSysHalfX, kPhenixSysHalfX,
                                   sysDownScaled, sysUpScaled);
    }

    // ---------------------------------------------------------------
    // PHENIX corrected to |eta|<0.7: undo bin-width via modified power
    // law fit, then rescale by 1/R where R = (dN/deta)_{|eta|<0.25} /
    // (dN/deta)_{|eta|<0.7} from truth_eta_ratio_inclusive.root
    // (Pythia truth, NO isolation — matches PHENIX inclusive fiducial).
    // Fit form A*(1+pT^2/b)^c (PHENIX 1405.3940) lifted to differential
    // d^2sigma/(deta dpT) by the 2*pi*pT factor. Fit window 10-26 GeV
    // to focus on the pT region overlapping PPG12.
    // ---------------------------------------------------------------
    TGraphErrors *gFit_PHENIX = new TGraphErrors(n);
    for (Int_t i = 0; i < n; ++i)
    {
        double pT = x[i];
        double sf = pT * factorCommon;
        double y_d  = y[i] * sf;
        double s_d  = 0.5 * (statUp[i] + statDown[i]) * sf;
        gFit_PHENIX->SetPoint(i, pT, y_d);
        gFit_PHENIX->SetPointError(i, 0.0, s_d);
    }
    TF1 *f_phenix_mpl = new TF1("f_phenix_mpl",
        "2.0*TMath::Pi()*x*[0]*pow(1.0 + x*x/[1], [2])",
        xLow[0], xHigh[n - 1]);
    f_phenix_mpl->SetParameters(8.3e-3, 2.26, -3.45);
    f_phenix_mpl->SetParLimits(1, 0.05, 1e3);
    gFit_PHENIX->Fit(f_phenix_mpl, "QRN");

    TFile *f_eta_ratio = TFile::Open(
        "/sphenix/user/shuhangli/ppg12/efficiencytool/truth_eta_ratio_inclusive.root");
    TH1D *h_eta_ratio = (TH1D *) f_eta_ratio->Get("h_ratio_central_over_full");

    TGraphAsymmErrors *gStat_PHENIX_corr = new TGraphAsymmErrors(n);
    gStat_PHENIX_corr->SetName("gStat_PHENIX_corr");
    TGraphAsymmErrors *gSys_PHENIX_corr = new TGraphAsymmErrors(n);
    gSys_PHENIX_corr->SetName("gSys_PHENIX_corr");

    for (Int_t i = 0; i < n; ++i)
    {
        double pT = x[i];
        double dpT = xHigh[i] - xLow[i];
        double y_binavg = f_phenix_mpl->Integral(xLow[i], xHigh[i]) / dpT;
        double f_center = f_phenix_mpl->Eval(pT);
        double bin_shape = (f_center > 0.0) ? y_binavg / f_center : 1.0;

        int rb = h_eta_ratio->FindBin(pT);
        if (rb < 1) rb = 1;
        if (rb > h_eta_ratio->GetNbinsX()) rb = h_eta_ratio->GetNbinsX();
        double R = h_eta_ratio->GetBinContent(rb);
        if (R <= 0.0) R = 1.0;

        double sf = pT * factorCommon;
        double y_pub = y[i] * sf;
        double total_scale = bin_shape / R;
        double y_corr = y_pub * total_scale;

        double exl = pT - xLow[i];
        double exh = xHigh[i] - pT;
        gStat_PHENIX_corr->SetPoint(i, pT, y_corr);
        gStat_PHENIX_corr->SetPointError(i, exl, exh,
            statDown[i] * sf * total_scale, statUp[i] * sf * total_scale);
        gSys_PHENIX_corr->SetPoint(i, pT, y_corr);
        gSys_PHENIX_corr->SetPointError(i, exl, exh,
            sysDown[i] * sf * total_scale, sysUp[i] * sf * total_scale);
    }
    f_eta_ratio->Close();

    // getting a different NLO
    ifstream myfile;
    TGraph *tphoton = new TGraph();
    TGraph *tphoton05 = new TGraph();
    TGraph *tphoton02 = new TGraph();
    // Werner Vogelsang NLO inputs: refreshed sPhenix-sc{05,1,2}.dat tables
    // (8-30 GeV, 1 GeV steps, 23 rows). Same 4-column format as the legacy
    // photons_newphenix_sc*.dat: pT, direct, fragmentation, total (last col
    // ignored — direct + frag is recomputed).
    const int n_werner = 23;
    myfile.open("/sphenix/user/shuhangli/ppg12/sPhenix-sc1.dat");
    if (!myfile)
    {
        cout << "No file found!" << endl;
        return;
    }

    for (int index = 0; index < n_werner; index++)
    {
        double pt;
        double yield;
        double fragyield;
        double foo;
        myfile >> pt >> yield >> fragyield >> foo;
        // rescale points
        yield = (yield + fragyield) * factorCommon * pt;
        tphoton->SetPoint(index, pt, yield);
    }
    myfile.close();
    TF1 *fphoton = new TF1("f", [&](double *x, double *)
                           { return tphoton->Eval(x[0]); }, lowerx, upperx, 0);

    myfile.open("/sphenix/user/shuhangli/ppg12/sPhenix-sc05.dat");
    if (!myfile)
    {
        cout << "No file found!" << endl;
        return;
    }

    for (int index = 0; index < n_werner; index++)
    {
        double pt;
        double yield;
        double fragyield;
        double foo;
        myfile >> pt >> yield >> fragyield >> foo;
        // rescale points
        yield = (yield + fragyield) * factorCommon * pt;
        tphoton05->SetPoint(index, pt, yield);
    }
    myfile.close();
    TF1 *fphoton05 = new TF1("f05", [&](double *x, double *)
                             { return tphoton05->Eval(x[0]); }, lowerx, upperx, 0);

    myfile.open("/sphenix/user/shuhangli/ppg12/sPhenix-sc2.dat");
    if (!myfile)
    {
        cout << "No file found!" << endl;
        return;
    }

    for (int index = 0; index < n_werner; index++)
    {
        double pt;
        double yield;
        double fragyield;
        double foo;
        myfile >> pt >> yield >> fragyield >> foo;
        // rescale points
        yield = (yield + fragyield) * factorCommon * pt;
        tphoton02->SetPoint(index, pt, yield);
    }
    myfile.close();
    TF1 *fphoton02 = new TF1("f02", [&](double *x, double *)
                             { return tphoton02->Eval(x[0]); }, lowerx, upperx, 0);

    TGraphAsymmErrors *g_syst_NLO_werner = (TGraphAsymmErrors *)g_syst_NLO->Clone("g_syst_NLO_werner");
    TGraphAsymmErrors *g_NLO_werner = (TGraphAsymmErrors *)g_syst_NLO->Clone("g_NLO_werner");
    TGraphAsymmErrors *g_syst_rel_NLO_werner = (TGraphAsymmErrors *)g_syst_rel_NLO->Clone("g_syst_rel_NLO_werner");
    TGraphAsymmErrors *g_rel_NLO_werner = (TGraphAsymmErrors *)g_syst_rel_NLO->Clone("g_rel_NLO_werner");

    // loop over bins
    for (int i = 0; i < g_rel_NLO_werner->GetN(); i++)
    {
        // get the bin center
        double x = g_rel_NLO_werner->GetX()[i];

        double xlow = g_rel_NLO_werner->GetErrorXlow(i);
        double xup = g_rel_NLO_werner->GetErrorXhigh(i);
        double xmin = x - xlow;
        double xmax = x + xup;

        float y_yield = fphoton->Integral(xmin, xmax) / (xmax - xmin);
        float y_yield05 = fphoton05->Integral(xmin, xmax) / (xmax - xmin);
        float y_yield02 = fphoton02->Integral(xmin, xmax) / (xmax - xmin);

        int binxdata = h_data->FindBin(x);
        float y_data = h_data->GetBinContent(binxdata);
        float y_data_err = h_data->GetBinError(binxdata);
        float data_rel_err = y_data_err / y_data;

        float ratio = y_yield / y_data;
        float ratio_err = ratio * data_rel_err;
        float ratio05 = y_yield05 / y_data;
        float ratio20 = y_yield02 / y_data;

        g_NLO_werner->SetPoint(i, x, y_yield);
        g_NLO_werner->SetPointError(i, xlow, xup, 0, 0);

        g_syst_NLO_werner->SetPoint(i, x, y_yield);
        g_syst_NLO_werner->SetPointError(i, xlow, xup, y_yield - y_yield02, y_yield05 - y_yield);

        g_rel_NLO_werner->SetPoint(i, x, ratio);
        g_rel_NLO_werner->SetPointError(i, xlow, xup, ratio_err, ratio_err);

        g_syst_rel_NLO_werner->SetPoint(i, x, ratio);
        g_syst_rel_NLO_werner->SetPointError(i, xlow, xup, ratio - ratio20, ratio05 - ratio);
    }

    tphoton->SetMarkerColor(mkcol[4]);
    tphoton->SetLineColor(mkcol[4]);
    // tphoton->Draw(" l,same");

    tphoton05->SetMarkerColor(mkcol[4]);
    tphoton05->SetLineColor(mkcol[4]);
    tphoton05->SetLineStyle(7);
    // tphoton05->Draw("l,same");

    tphoton02->SetMarkerColor(mkcol[4]);
    tphoton02->SetLineColor(mkcol[4]);
    tphoton02->SetLineStyle(7);
    // tphoton02->Draw("l,same");

    // gStat_PHENIX->Draw("P");
    // gSys_PHENIX->Draw("2 same");

    // Werner
    g_syst_NLO_werner->SetMarkerStyle(mkStyle[4]);
    g_syst_NLO_werner->SetMarkerSize(mkSize[4]);
    g_syst_NLO_werner->SetMarkerColor(mkcol[4]);
    g_syst_NLO_werner->SetLineColor(mkcol[4]);
    g_syst_NLO_werner->SetFillColorAlpha(col[4], trans[4]);
    g_syst_NLO_werner->Draw("5 same");

    // Data syst
    g_syst->SetMarkerStyle(mkStyle[0]);
    g_syst->SetMarkerColor(col[0]);
    g_syst->SetLineColor(col[0]);
    g_syst->SetFillColorAlpha(col[0], trans[0]);
    g_syst->Draw("2 same");

    // JETPHOX syst
    g_syst_NLO->SetMarkerStyle(mkStyle[2]);
    g_syst_NLO->SetMarkerColor(col[2]);
    g_syst_NLO->SetLineColor(col[2]);
    g_syst_NLO->SetFillColorAlpha(col[2], trans[2]);
    g_syst_NLO->Draw("2 same");

    // Werner points
    g_NLO_werner->SetMarkerStyle(mkStyle[4]);
    g_NLO_werner->SetMarkerSize(mkSize[4]);
    g_NLO_werner->SetMarkerColor(mkcol[4]);
    g_NLO_werner->SetLineColor(mkcol[4]);
    g_NLO_werner->SetLineWidth(lineWidth[4]);
    g_NLO_werner->Draw("p same");

    // JETPHOX
    h_NLO->SetMarkerStyle(mkStyle[2]);
    h_NLO->SetMarkerSize(mkSize[2]);
    h_NLO->SetMarkerColor(mkcol[2]);
    h_NLO->SetLineColor(mkcol[2]);
    h_NLO->SetFillColorAlpha(col[2], trans[2]);
    h_NLO->SetLineWidth(lineWidth[2]);
    h_NLO->Draw("same");

    h_pythia->SetMarkerStyle(mkStyle[3]);
    h_pythia->SetMarkerSize(mkSize[3]);
    h_pythia->SetMarkerColor(col[3]);
    h_pythia->SetLineColor(col[3]);
    h_pythia->SetLineWidth(lineWidth[3]);
    h_pythia->Draw("same");

    h_data->SetMarkerStyle(mkStyle[0]);
    h_data->SetMarkerSize(mkSize[0]);
    h_data->SetMarkerColor(col[0]);
    h_data->SetLineColor(col[0]);
    h_data->SetLineWidth(lineWidth[0]);
    h_data->Draw("same");

    TH1F *htemp_data = (TH1F *)h_data->Clone("htemp");
    htemp_data->SetFillColorAlpha(col[0], trans[0]);
    TH1F *htemp_NLO = (TH1F *)h_NLO->Clone("htemp_NLO");
    htemp_NLO->SetFillColorAlpha(col[2], trans[2]);

    //--------------------------------------------------lower panel

    float xpos(0.15), xpos2(0.875), ypos(0.87), ypos2(0.1), dy(0.065), dy1(0.078), fontsize(0.052), fontsize1(0.055);
    myText(xpos2, ypos - 0 * dy, 1, strleg1.c_str(), fontsize1, 1);
    myText(xpos2, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 1);
    myText(xpos2, ypos - 2 * dy, 1, strleg_lumi_line2.c_str(), fontsize, 1);
    myText(xpos2, ypos - 3 * dy, 1, strleg3.c_str(), fontsize, 1);
    myText(xpos2, ypos - 4 * dy, 1, strleg4.c_str(), fontsize, 1);
    // myText(xpos2,ypos-1*dy,1,strleg2_1.c_str(),fontsize,1);
    // myText(xpos2,ypos-2*dy,1,strleg3.c_str(),fontsize,1);
    // myText(xpos2,ypos-3*dy,1,strleg4.c_str(),fontsize,1);

    int nEntry = 6;
    TLegend *l1 = new TLegend(xpos, ypos2, 0.5, ypos2 + nEntry * dy1);
    legStyle(l1, 0.20, fontsize);
    l1->AddEntry(htemp_data, "Data", "fpl");
    l1->AddEntry(h_pythia, "PYTHIA8", "pl");
    l1->AddEntry(htemp_NLO, "NLO pQCD JETPHOX", "fpl");

    // #2: shared scale-choice annotation for ALL theory predictions (was
    // previously duplicated in two places — once after JETPHOX line, once
    // after Vogelsang). Keep PDF/FF info on the JETPHOX line; combine the
    // common scale choice into a single sub-caption below the legend.
    string st_thScale = "#kern[-0.55]{#it{#mu}_{f}} = #kern[-0.55]{#it{#mu}_{F}} = #kern[-0.55]{#it{#mu}_{R}} = #kern[-0.55]{#it{E}_{T}^{#gamma}}";
    l1->AddEntry((TObject *)0, "#scale[0.93]{BFG II FF}", "");
    // l1->AddEntry((TObject *)0, "#scale[0.93]{CT18NLO PDF / BFG II FF}", "");
    l1->AddEntry(g_syst_NLO_werner, "NLO pQCD by W. Vogelsang", "fp");
    // l1->AddEntry(g_syst_NLO_werner, "NLO pQCD by W. Vogelsang", "fpl");
    l1->AddEntry((TObject *)0, "#scale[0.93]{GRV FF}", "");
    l1->Draw("same");
    // Single shared scale-choice line below the legend (replaces the two
    // floating myText calls that previously sat between the legend rows).
    myText(xpos, ypos2 - 0.06, 1, Form("#scale[0.93]{All NLO pQCD: CT18NLO PDF / %s}", st_thScale.data()), fontsize, 0);

    TPad *pad_2 = (TPad *)c1->cd(2);
    pad_2->SetPad(0, 0.25, 1, 0.5);
    pad_2->SetTopMargin(0.023);
    pad_2->SetLeftMargin(0.13);
    pad_2->SetBottomMargin(0.002);
    pad_2->SetRightMargin(0.08);

    frame_et_truth->SetYTitle("Theory / Data");
    // x-axis labels are hidden here (the third pad below carries them).
    frame_et_truth->SetXTitle("");
    frame_et_truth->GetXaxis()->SetLabelSize(0);
    frame_et_truth->GetYaxis()->SetNdivisions(509);
    frame_et_truth->GetYaxis()->SetRangeUser(0.61, 1.79);  // #8: widen so band doesn't clip; 0.5 lower edge avoids "0.4" label clipping at the pad seam
    // frame_et_truth->GetYaxis()->SetRangeUser(0.5, 2.2);  // #8: widen so band doesn't clip; 0.5 lower edge avoids "0.4" label clipping at the pad seam
    frame_et_truth->GetXaxis()->SetRangeUser(lowerx, upperx);
    frame_et_truth->GetXaxis()->SetTitleOffset(frame_et_rec->GetXaxis()->GetTitleOffset() * 4 / 6. * 1.4);
    frame_et_truth->GetYaxis()->SetTitleOffset(frame_et_rec->GetYaxis()->GetTitleOffset() * 4 / 6. * 0.7);
    frame_et_truth->GetYaxis()->SetLabelOffset(frame_et_rec->GetYaxis()->GetLabelOffset() * 4 / 6.);
    // Derive the ratio-panel x-label size from the *y*-axis label of the top
    // panel (which is preserved at 0.050). The top-panel x-label is hidden
    // (size=0) per the two-pad layout, so reading from it would yield 0.
    frame_et_truth->GetXaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetYaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize() * 6 / 4. * 1.1);
    frame_et_truth->GetXaxis()->SetTitleSize(frame_et_rec->GetXaxis()->GetTitleSize() * 6 / 4. * 1.25);
    frame_et_truth->GetYaxis()->SetTitleSize(frame_et_rec->GetYaxis()->GetTitleSize() * 6 / 4. * 1.25);
    frame_et_truth->GetXaxis()->SetNdivisions(505);  // 5 primary divisions, optimized — matches top panel; gives 15/20/25/30 labels
    frame_et_truth->Draw("axis");

    g_syst_rel->SetMarkerStyle(mkStyle[0]);
    g_syst_rel->SetMarkerColor(col[0]);
    g_syst_rel->SetLineColor(col[0]);
    g_syst_rel->SetFillColorAlpha(col[0], trans[0]);

    g_syst_rel->Draw("2 same");
    // #9: unity reference drawn AFTER the band so it is not overpainted; also
    // bumped to width=2 dashed black so it stays the dominant eyeline.
    lineone->SetLineColor(kBlack);
    lineone->SetLineStyle(2);
    lineone->SetLineWidth(2);
    lineone->Draw("L same");

    // Werner syst
    g_syst_rel_NLO_werner->SetMarkerStyle(mkStyle[4]);
    g_syst_rel_NLO_werner->SetMarkerColor(mkcol[4]);
    g_syst_rel_NLO_werner->SetLineColor(mkcol[4]);
    g_syst_rel_NLO_werner->SetLineWidth(lineWidth[4]);
    g_syst_rel_NLO_werner->SetFillColorAlpha(col[4], trans[4]);
    g_syst_rel_NLO_werner->Draw("5 same");

    // JETPHOX syst
    g_syst_rel_NLO->SetMarkerStyle(mkStyle[2]);
    g_syst_rel_NLO->SetMarkerColor(col[2]);
    g_syst_rel_NLO->SetLineColor(col[2]);
    g_syst_rel_NLO->SetLineWidth(lineWidth[2]);
    g_syst_rel_NLO->SetFillColorAlpha(col[2], trans[2]);
    g_syst_rel_NLO->Draw("2 same");


    // Werner/Data
    g_rel_NLO_werner->SetMarkerStyle(mkStyle[4]);
    g_rel_NLO_werner->SetMarkerSize(mkSize[4]);
    g_rel_NLO_werner->SetMarkerColor(mkcol[4]);
    g_rel_NLO_werner->SetLineColor(mkcol[4]);
    g_rel_NLO_werner->SetLineWidth(lineWidth[4]);
    g_rel_NLO_werner->Draw("p same");

    // Pythia/Data
    h_pythia_data->SetMarkerStyle(mkStyle[3]);
    h_pythia_data->SetMarkerColor(col[3]);
    h_pythia_data->SetLineColor(col[3]);
    h_pythia_data->SetMarkerSize(mkSize[3]);
    h_pythia_data->SetLineWidth(lineWidth[3]);
    h_pythia_data->Draw("same");
    // JETPHOX/Data
    h_NLO_data->SetMarkerStyle(mkStyle[2]);
    h_NLO_data->SetMarkerSize(mkSize[2]);
    h_NLO_data->SetMarkerColor(mkcol[2]);
    h_NLO_data->SetLineColor(mkcol[2]);
    h_NLO_data->SetLineWidth(lineWidth[2]);
    // h_NLO_data->SetFillColorAlpha(col[2], trans[2]);
    h_NLO_data->Draw("same");

    // ----------------------------------------------------------------
    // Third pad: JETPHOX/Data ratio for three NLO PDF sets
    // (CT18NLO, CT14NLO, NNPDF3.1) at the central scale mu=pT.
    // ----------------------------------------------------------------
    TPad *pad_3 = (TPad *)c1->cd(3);
    pad_3->SetPad(0, 0, 1, 0.25);
    pad_3->SetTopMargin(0.023);
    pad_3->SetLeftMargin(0.13);
    pad_3->SetBottomMargin(0.30);
    pad_3->SetRightMargin(0.08);

    TH1F *frame_pdf = (TH1F *)frame_et_truth->Clone("frame_pdf");
    frame_pdf->SetYTitle("JETPHOX / Data");
    frame_pdf->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");
    frame_pdf->GetYaxis()->SetRangeUser(0.71, 1.6);
    // frame_pdf->GetYaxis()->SetRangeUser(0.4, 2.0);
    frame_pdf->GetYaxis()->SetNdivisions(505);
    frame_pdf->GetYaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize() * 6.0 / 4.0 * 1.1);
    frame_pdf->GetXaxis()->SetLabelSize(frame_et_rec->GetXaxis()->GetLabelSize() * 6.0 / 4.0 * 1.1);
    frame_pdf->GetXaxis()->SetTitleSize(frame_et_rec->GetXaxis()->GetTitleSize() * 6.0 / 4.0 * 1.4);
    frame_pdf->GetXaxis()->SetTitleOffset(frame_et_rec->GetXaxis()->GetTitleOffset() * 4.0 / 6.0 * 1.2);
    frame_pdf->Draw("axis");

    // const Color_t kCTEQcol  = kOrange + 7;
    // const Color_t kMSHTcol  = kAzure - 4;
    // const Color_t kCT14col  = kSpring - 7;
    // const Color_t kCT18col  = mkcol[2];
    // const Color_t kNNPDFcol = kRed - 4;

    // 0: CT18, 1: NNPDF, 2: CTEQ, 3: MSHT, 4: CT14
    h_NLO_data_ct14->SetMarkerStyle(24);
    h_NLO_data_ct14->SetMarkerSize(mkSizepdf[4]);
    h_NLO_data_ct14->SetMarkerColor(mkcolpdf[4]);
    h_NLO_data_ct14->SetLineColor(mkcolpdf[4]);
    h_NLO_data_ct14->SetLineWidth(lineWidth[2]);

    h_NLO_data->SetMarkerStyle(mkStyle[2]);
    h_NLO_data->SetMarkerSize(mkSizepdf[0]);
    h_NLO_data->SetMarkerColor(mkcolpdf[0]);
    h_NLO_data->SetLineColor(mkcolpdf[0]);
    h_NLO_data->SetLineWidth(lineWidth[2]);

    h_NLO_data_nnpdf->SetMarkerStyle(24);
    h_NLO_data_nnpdf->SetMarkerSize(mkSizepdf[1]);
    h_NLO_data_nnpdf->SetMarkerColor(mkcolpdf[1]);
    h_NLO_data_nnpdf->SetLineColor(mkcolpdf[1]);
    h_NLO_data_nnpdf->SetLineWidth(lineWidth[2]);

    h_NLO_data_cteq->SetMarkerStyle(27);
    h_NLO_data_cteq->SetMarkerSize(mkSizepdf[2]);
    h_NLO_data_cteq->SetMarkerColor(mkcolpdf[2]);
    h_NLO_data_cteq->SetLineColor(mkcolpdf[2]);
    h_NLO_data_cteq->SetLineWidth(lineWidth[2]);

    h_NLO_data_msht->SetMarkerStyle(28);
    h_NLO_data_msht->SetMarkerSize(mkSizepdf[3]);
    h_NLO_data_msht->SetMarkerColor(mkcolpdf[3]);
    h_NLO_data_msht->SetLineColor(mkcolpdf[3]);
    h_NLO_data_msht->SetLineWidth(lineWidth[2]);

    // Bands: data systematic centred on y=1 (azure), plus per-PDF JETPHOX
    // scale-uncertainty envelopes (pink/red/orange/blue) tracking each
    // PDF's nominal ratio. Drawn first so the PDF markers render on top.
    const float kPdfBandAlpha = 0.25;
    g_syst_rel       ->Draw("2 same");
    g_syst_rel_NLO   ->SetFillColorAlpha(col[2],   kPdfBandAlpha);
    // g_syst_rel_NLO   ->Draw("2 same");
    g_syst_rel_nnpdf ->SetFillColorAlpha(mkcolpdf[1], kPdfBandAlpha);
    g_syst_rel_nnpdf ->SetLineColor(mkcolpdf[1]);
    // g_syst_rel_nnpdf ->Draw("2 same");
    g_syst_rel_cteq  ->SetFillColorAlpha(mkcolpdf[2], kPdfBandAlpha);
    g_syst_rel_cteq  ->SetLineColor(mkcolpdf[2]);
    // g_syst_rel_cteq  ->Draw("2 same");
    g_syst_rel_msht  ->SetFillColorAlpha(mkcolpdf[3], kPdfBandAlpha);
    g_syst_rel_msht  ->SetLineColor(mkcolpdf[3]);
    // g_syst_rel_msht  ->Draw("2 same");

    // Per-PDF Hessian uncertainty bands, drawn first so the central lines
    // render on top. Filled box per ET bin in each PDF's color, low alpha
    // so the four bands remain readable when they overlap.
    const float kPdfHessAlpha = 0.30;
    g_pdfband_ct18 ->SetFillColorAlpha(mkcolpdf[0], kPdfHessAlpha);
    g_pdfband_ct18 ->SetLineColor(mkcolpdf[0]);
    g_pdfband_ct18 ->SetLineWidth(2);
    g_pdfband_ct18 ->Draw("3 same");
    g_pdfband_nnpdf->SetFillColorAlpha(mkcolpdf[1], kPdfHessAlpha);
    g_pdfband_nnpdf->SetLineColor(mkcolpdf[1]);
    g_pdfband_nnpdf->SetLineWidth(2);
    g_pdfband_nnpdf->Draw("3 same");
    g_pdfband_cteq ->SetFillColorAlpha(mkcolpdf[2], kPdfHessAlpha);
    g_pdfband_cteq ->SetLineColor(mkcolpdf[2]);
    g_pdfband_cteq ->SetLineWidth(2);
    g_pdfband_cteq ->Draw("3 same");
    g_pdfband_msht ->SetFillColorAlpha(mkcolpdf[3], kPdfHessAlpha);
    g_pdfband_msht ->SetLineColor(mkcolpdf[3]);
    g_pdfband_msht ->SetLineWidth(2);
    g_pdfband_msht ->Draw("3 same");

    lineone->SetLineColor(kBlack);
    lineone->SetLineStyle(2);
    lineone->SetLineWidth(2);
    lineone->Draw("L same");

    h_NLO_data      ->Draw("C hist same");
    h_NLO_data_nnpdf->Draw("C hist same");
    h_NLO_data_cteq ->Draw("C hist same");
    h_NLO_data_msht ->Draw("C hist same");

    TLegend *l_pdf = new TLegend(0.14, 0.85, 0.92, 0.96);
    l_pdf->SetBorderSize(0);
    l_pdf->SetFillStyle(0);
    l_pdf->SetTextFont(42);
    l_pdf->SetTextSize(0.085);
    l_pdf->SetNColumns(4);
    // Use the band+line composite (filled rectangle around the central
    // line) so the legend marker conveys the Hessian uncertainty band too.
    l_pdf->AddEntry(g_pdfband_ct18,  "CT18NLO",  "lf");
    l_pdf->AddEntry(g_pdfband_nnpdf, "NNPDF4.0", "lf");
    l_pdf->AddEntry(g_pdfband_cteq,  "CTEQ6.6",  "lf");
    l_pdf->AddEntry(g_pdfband_msht,  "MSHT20NLO","lf");
    l_pdf->Draw("same");

    // Paper Fig. final.pdf -- 3-pad cross-section + theory ratios.
    std::string outputname = Form("%s/final.pdf", paperdir.c_str());
    c1->SaveAs(outputname.c_str());

    pad_1->cd();
    frame_et_rec->GetYaxis()->SetRangeUser(1, 2e4);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 30);
    frame_et_rec->Draw("axis");

    h_common_cluster_data->SetMarkerStyle(20);
    h_common_cluster_data->SetMarkerColor(kBlack);
    h_common_cluster_data->SetLineColor(kBlack);

    h_common_cluster_data->Draw("same");

    h_common_cluster_mc->SetMarkerStyle(24);
    h_common_cluster_mc->SetMarkerColor(kPink + 8);
    h_common_cluster_mc->SetLineColor(kPink + 8);

    h_common_cluster_mc->Draw("same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, "Clusters pass preliminary cuts", 0.05);

    myMarkerLineText(0.25, 0.25, 1, kBlack, 20, kBlack, 1, Form("%s", datastring.c_str()), 0.05, true);
    myMarkerLineText(0.25, 0.20, 1, kPink + 8, 24, kPink + 8, 1, bg_MCstring.c_str(), 0.05, true);

    pad_2->cd();
    frame_et_truth->GetYaxis()->SetRangeUser(0.5, 3.0);
    frame_et_truth->GetXaxis()->SetRangeUser(10, 30);
    frame_et_truth->SetYTitle("MC/data");
    frame_et_truth->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
    frame_et_truth->Draw("axis");

    lineone->Draw("L");

    TH1F *h_common_cluster_data_ratio = (TH1F *)h_common_cluster_mc->Clone("h_common_cluster_data_ratio");
    h_common_cluster_data_ratio->Divide(h_common_cluster_data);

    h_common_cluster_data_ratio->SetMarkerStyle(20);
    h_common_cluster_data_ratio->SetMarkerColor(kBlack);
    h_common_cluster_data_ratio->SetLineColor(kBlack);

    h_common_cluster_data_ratio->Draw("same");

    c1->SaveAs(Form("%s/final_common_cluster_%s.pdf", scratch.c_str(), tune.data()));

    pad_1->cd();
    frame_et_rec->GetYaxis()->SetRangeUser(0.1, 1e3);
    frame_et_rec->GetXaxis()->SetRangeUser(10, 30);
    frame_et_rec->Draw("axis");

    h_tight_iso_cluster_data->SetMarkerStyle(20);
    h_tight_iso_cluster_data->SetMarkerColor(kBlack);
    h_tight_iso_cluster_data->SetLineColor(kBlack);

    h_tight_iso_cluster_data->Draw("same");

    h_tight_iso_cluster_mc->SetMarkerStyle(24);
    h_tight_iso_cluster_mc->SetMarkerColor(kPink + 8);
    h_tight_iso_cluster_mc->SetLineColor(kPink + 8);

    h_tight_iso_cluster_mc->Draw("same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, "Tight iso clusters", 0.05);

    myMarkerLineText(0.25, 0.25, 1, kBlack, 20, kBlack, 1, Form("%s", datastring.c_str()), 0.05, true);
    myMarkerLineText(0.25, 0.20, 1, kPink + 8, 24, kPink + 8, 1, bg_MCstring.c_str(), 0.05, true);

    pad_2->cd();
    frame_et_truth->GetYaxis()->SetRangeUser(0.2, 3.0);
    frame_et_truth->GetXaxis()->SetRangeUser(10, 30);
    frame_et_truth->SetYTitle("MC/data");
    frame_et_truth->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
    frame_et_truth->Draw("axis");

    lineone->Draw("L");

    TH1F *h_tight_iso_cluster_data_ratio = (TH1F *)h_tight_iso_cluster_mc->Clone("h_tight_iso_cluster_data_ratio");
    h_tight_iso_cluster_data_ratio->Divide(h_tight_iso_cluster_data);

    h_tight_iso_cluster_data_ratio->SetMarkerStyle(20);
    h_tight_iso_cluster_data_ratio->SetMarkerColor(kBlack);
    h_tight_iso_cluster_data_ratio->SetLineColor(kBlack);

    h_tight_iso_cluster_data_ratio->Draw("same");

    c1->SaveAs(Form("%s/final_tight_iso_cluster_%s.pdf", scratch.c_str(), tune.data()));

    /*
    TCanvas *c2 = new TCanvas("can", "", 800, 900);
    //log y
    c2->SetLogy();
    TH2F *frame = new TH2F("frame", "", 1, 10, 30, 1, 0.1, 1e5);
    frame->SetXTitle("E_{T}^{#gamma} [GeV]");
    frame->SetYTitle("d#sigma/dE_{T} [pb/GeV]");
    frame->GetYaxis()->SetRangeUser(0.5, 5e4);
    frame->GetXaxis()->SetRangeUser(10, 30);

    frame->Draw("axis");
    */
    pad_1->cd();
    // common, tight iso, sub, sub with unfold and final result
    frame_et_rec->GetYaxis()->SetRangeUser(1, 1e5);
    frame_et_rec->GetXaxis()->SetRangeUser(12, 32);
    frame_et_rec->SetYTitle("dN/dE_{T} [GeV^{-1}]");
    frame_et_rec->Draw("axis");

    std::vector<int> marker_styles = {20, 24, 25, 27, 28};
    std::vector<float> marker_sizes = {1, 1, 1.0, 1.5, 1.5};
    std::vector<int> colors = {kPink + 8, kSpring - 7, kAzure - 3, kOrange + 10, kViolet + 3, kYellow + 4};

    // scale things up by lumi
    h_common_cluster_data->Scale(datalumi);
    h_tight_iso_cluster_data->Scale(datalumi);
    h_sub_data->Scale(datalumi);
    h_sub_data_unfold->Scale(datalumi);
    h_data->Scale(datalumi);

    h_common_cluster_data->SetMarkerStyle(marker_styles[0]);
    h_common_cluster_data->SetMarkerColor(colors[0]);
    h_common_cluster_data->SetLineColor(colors[0]);
    h_common_cluster_data->SetMarkerSize(marker_sizes[0]);

    h_common_cluster_data->Draw("same");

    h_tight_iso_cluster_data->SetMarkerStyle(marker_styles[1]);
    h_tight_iso_cluster_data->SetMarkerColor(colors[1]);
    h_tight_iso_cluster_data->SetLineColor(colors[1]);
    h_tight_iso_cluster_data->SetMarkerSize(marker_sizes[1]);
    h_tight_iso_cluster_data->Draw("same");

    h_sub_data->SetMarkerStyle(marker_styles[2]);
    h_sub_data->SetMarkerColor(colors[2]);
    h_sub_data->SetLineColor(colors[2]);
    h_sub_data->SetMarkerSize(marker_sizes[2]);
    h_sub_data->Draw("same");

    h_sub_data_unfold->SetMarkerStyle(marker_styles[3]);
    h_sub_data_unfold->SetMarkerColor(colors[3]);
    h_sub_data_unfold->SetLineColor(colors[3]);
    h_sub_data_unfold->SetMarkerSize(marker_sizes[3]);
    h_sub_data_unfold->Draw("same");

    h_data->SetMarkerStyle(marker_styles[4]);
    h_data->SetMarkerColor(colors[4]);
    h_data->SetLineColor(colors[4]);
    h_data->SetMarkerSize(marker_sizes[4]);
    h_data->Draw("same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, "|#eta^{#gamma}|<0.7", 0.05);

    myMarkerLineText(0.6, 0.25 + 0.5, 1, colors[0], marker_styles[0], colors[0], 1, "w/ ps cut", 0.05, true);
    myMarkerLineText(0.6, 0.20 + 0.5, 1, colors[1], marker_styles[1], colors[1], 1, "tight iso", 0.05, true);
    myMarkerLineText(0.6, 0.15 + 0.5, 1, colors[2], marker_styles[2], colors[2], 1, "purity corrected", 0.05, true);
    myMarkerLineText(0.6, 0.10 + 0.5, 1.5, colors[3], marker_styles[3], colors[3], 1, "unfolded", 0.05, true);
    myMarkerLineText(0.6, 0.05 + 0.5, 1.5, colors[4], marker_styles[4], colors[4], 1, "efficiency corrected", 0.05, true);

    pad_2->cd();
    frame_et_truth->GetYaxis()->SetRangeUser(0.8, 1.05);
    frame_et_truth->GetXaxis()->SetRangeUser(12, 32);
    frame_et_truth->SetYTitle("after/before unfolding");
    frame_et_truth->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");
    frame_et_truth->Draw("axis");

    lineone->Draw("L");

    TH1F *h_unfold_ratio = (TH1F *)h_sub_data_unfold->Clone("h_unfold_ratio");
    for (int i = 2; i <= h_unfold_ratio->GetNbinsX(); i++)
    {
        float centerX = h_unfold_ratio->GetBinCenter(i);
        int bin = h_sub_data->FindBin(centerX);
        float value = h_sub_data->GetBinContent(bin);
        float error = h_sub_data->GetBinError(bin);

        if (value != 0)
        {
            h_unfold_ratio->SetBinContent(i, h_sub_data_unfold->GetBinContent(i) / value);
            std::cout << "bincenter: " << centerX << " " << h_sub_data_unfold->GetBinContent(i) << " " << value << std::endl;
            h_unfold_ratio->SetBinError(i, 0);
        }
    }

    h_unfold_ratio->SetMarkerStyle(marker_styles[3]);
    h_unfold_ratio->SetMarkerColor(colors[3]);
    h_unfold_ratio->SetLineColor(colors[3]);
    h_unfold_ratio->SetMarkerSize(marker_sizes[3]);
    h_unfold_ratio->Draw("same");

    c1->SaveAs(Form("%s/final_all_%s.pdf", scratch.c_str(), tune.data()));

    //-----------------------------------------------------------------
    // only compare to PHENIX
    init_plot();
    h_data = h_data_cp;
    // Two-pad layout mirrors c1 (60/40 split) so font scaling is reusable.
    TCanvas *c2 = new TCanvas("can2", "", 800, 889);

    TPad *pad2_top = new TPad("pad2_top", "", 0, 0.4, 1, 1);
    pad2_top->SetTopMargin(0.05);
    pad2_top->SetBottomMargin(0.002);
    pad2_top->SetLeftMargin(0.13);
    pad2_top->SetRightMargin(0.08);
    pad2_top->SetLogy();
    pad2_top->Draw();

    TPad *pad2_bot = new TPad("pad2_bot", "", 0, 0, 1, 0.4);
    pad2_bot->SetTopMargin(0.023);
    pad2_bot->SetBottomMargin(0.25);
    pad2_bot->SetLeftMargin(0.13);
    pad2_bot->SetRightMargin(0.08);
    pad2_bot->Draw();

    pad2_top->cd();
    frame_et_rec->SetTitle(";;d^{2}#sigma/d#it{#eta}d#it{E}_{T}^{#gamma} [pb/GeV]");
    frame_et_rec->GetYaxis()->SetRangeUser(lowery, 500);
    frame_et_rec->GetXaxis()->SetRangeUser(lowerx, upperx);
    frame_et_rec->GetXaxis()->SetLabelSize(0);
    frame_et_rec->GetYaxis()->SetTitleSize(0.069);
    frame_et_rec->GetYaxis()->SetLabelSize(0.055);
    frame_et_rec->GetYaxis()->SetTitleOffset(0.85);

    frame_et_rec->Draw("axis");

    gStat_PHENIX->SetMarkerStyle(mkStyle[1]);
    gStat_PHENIX->SetMarkerSize(mkSize[1]);
    gStat_PHENIX->SetMarkerColor(col[1]);
    gStat_PHENIX->SetLineColor(col[1]);
    gStat_PHENIX->Draw("P");

    gSys_PHENIX->SetMarkerStyle(mkStyle[1]);
    gSys_PHENIX->SetMarkerSize(mkSize[1]);
    gSys_PHENIX->SetMarkerColor(col[1]);
    gSys_PHENIX->SetLineColor(col[1]);
    gSys_PHENIX->SetFillColorAlpha(col[1], 0.35);
    gSys_PHENIX->Draw("2 same");

    const Color_t kCorrColor = kViolet + 1;
    int PHENIXcorrMarker = 54; 
    gSys_PHENIX_corr->SetMarkerStyle(PHENIXcorrMarker);
    gSys_PHENIX_corr->SetMarkerSize(mkSize[1]);
    gSys_PHENIX_corr->SetMarkerColor(kCorrColor);
    gSys_PHENIX_corr->SetLineColor(kCorrColor);
    gSys_PHENIX_corr->SetFillColorAlpha(kCorrColor, 0.30);
    gSys_PHENIX_corr->Draw("2 same");

    gStat_PHENIX_corr->SetMarkerStyle(25);
    gStat_PHENIX_corr->SetMarkerSize(mkSize[1]);
    gStat_PHENIX_corr->SetMarkerColor(kCorrColor);
    gStat_PHENIX_corr->SetLineColor(kCorrColor);
    gStat_PHENIX_corr->Draw("P same");

    g_syst->SetMarkerStyle(mkStyle[0]);
    g_syst->SetMarkerColor(col[0]);
    g_syst->SetLineColor(col[0]);
    g_syst->SetFillColorAlpha(col[0], 0.35);

    g_syst->Draw("2 same");

    h_data->SetMarkerStyle(mkStyle[0]);
    h_data->SetMarkerSize(mkSize[0]);
    h_data->SetMarkerColor(col[0]);
    h_data->SetLineColor(col[0]);

    h_data->Draw("same");

    TH1F *htemp_PHENIX = (TH1F *)h_data->Clone("htemp_PHENIX");
    htemp_PHENIX->SetMarkerStyle(mkStyle[1]);
    htemp_PHENIX->SetMarkerSize(mkSize[1]);
    htemp_PHENIX->SetMarkerColor(col[1]);
    htemp_PHENIX->SetLineColor(col[1]);
    htemp_PHENIX->SetLineWidth(2);
    htemp_PHENIX->SetFillColorAlpha(col[1], 0.25);

    TH1F *htemp_PHENIX_corr = (TH1F *)h_data->Clone("htemp_PHENIX_corr");
    htemp_PHENIX_corr->SetMarkerStyle(PHENIXcorrMarker);
    htemp_PHENIX_corr->SetMarkerSize(mkSize[1]);
    htemp_PHENIX_corr->SetMarkerColor(kViolet + 1);
    htemp_PHENIX_corr->SetLineColor(kViolet + 1);
    htemp_PHENIX_corr->SetLineWidth(2);
    htemp_PHENIX_corr->SetFillColorAlpha(kViolet + 1, 0.30);

    // float xpos(0.15), xpos2(0.875), ypos(0.87), ypos2(0.1), dy(0.065), dy1(0.078), fontsize(0.052), fontsize1(0.055);
    xpos2 = 0.87;
    fontsize = 0.053;
    fontsize1 = 0.047;
    dy = 0.065;
    xpos = 0.16;
    ypos2 = 0.09;
    myText(xpos2, ypos - 0 * dy, 1, strleg1.c_str(), fontsize1, 1);
    myText(xpos2, ypos - 1 * dy, 1, strleg_lumi.c_str(), fontsize, 1);
    myText(xpos2, ypos - 2 * dy, 1, strleg3.c_str(), fontsize, 1);
    myText(xpos2, ypos - 3 * dy, 1, strleg4.c_str(), fontsize, 1);
    // myText(xpos2,ypos-1*dy,1,strleg2.c_str(),fontsize,1);
    // myText(xpos2,ypos-2*dy,1,strleg5.c_str(),fontsize,1);
    // myText(xpos2,ypos-3*dy,1,strleg3.c_str(),fontsize,1);
    // myText(xpos2,ypos-4*dy,1,strleg4.c_str(),fontsize,1);

    dy1 = 0.083;
    nEntry = 4;
    TLegend *l2 = new TLegend(xpos, ypos2, 0.6, ypos2 + nEntry * dy1);
    legStyle(l2, 0.21, fontsize);
    l2->AddEntry(htemp_data, "Data", "fpl");
    l2->AddEntry(htemp_PHENIX, "#scale[0.93]{PHENIX |#eta^{#gamma}|<0.25}", "fpl");
    // l2->AddEntry((TObject*)0, "#scale[0.93]{#it{PRD 86 072008}}", "");
    l2->AddEntry(htemp_PHENIX_corr, "#scale[0.93]{PHENIX}", "fpl");
    // l2->AddEntry(htemp_PHENIX_corr, "#scale[0.93]{PHENIX, (corrected for bin-avg, |#eta^{#gamma}|<0.7)}", "fpl");
    l2->AddEntry((TObject*)0, "#lower[-0.35]{#scale[0.93]{corrected for bin-avg and |#eta^{#gamma}|<0.7}}", "");
    // l2->AddEntry((TObject*)0, "#scale[0.93]{(no #kern[-0.2]{#it{E}_{T}^{iso}} requirement)}", "");
    l2->Draw("same");
    myText(xpos+0.02, ypos2 - 0.04, 1, "#scale[0.93]{(PHENIX #kern[-0.05]{#it{PRD 86 072008}}: no #kern[-0.2]{#it{E}_{T}^{iso}} requirement)}", fontsize, 0);
    // myText(xpos + 0.08, ypos2 - 0.03, 1, "#scale[0.93]{(PHENIX: no #kern[-0.2]{#it{E}_{T}^{iso}} requirement)}", fontsize, 0);

    // -----------------------------------------------------------------
    // Bottom pad: PHENIX_{|eta|<0.7-corrected} / data ratio.
    // Mirrors the c1 (theory/data) convention: per-point markers with
    // combined stat as vertical bars and bin-width horizontals, the
    // numerator (PHENIX-corr) systematic as filled boxes around the
    // markers, and the denominator (data) systematic as a band at y=1.
    // -----------------------------------------------------------------
    pad2_bot->cd();

    TH1F *frame_ratio_phenix = new TH1F("frame_ratio_phenix", "", 1, lowerx, upperx);
    frame_ratio_phenix->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");
    frame_ratio_phenix->SetYTitle("PHENIX / Data");
    frame_ratio_phenix->GetYaxis()->SetRangeUser(0.6, 2.3);
    frame_ratio_phenix->GetYaxis()->SetNdivisions(505);
    // Match the c1 lower-panel font scaling (6/4 = 1.5) for visual parity.
    const double k_p2bot_scale = 6.0 / 4.0;
    frame_ratio_phenix->GetXaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize() * k_p2bot_scale * 1.0);
    frame_ratio_phenix->GetYaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize() * k_p2bot_scale * 1.0);
    frame_ratio_phenix->GetXaxis()->SetTitleSize(frame_et_rec->GetYaxis()->GetTitleSize() * k_p2bot_scale * 1.0);
    frame_ratio_phenix->GetYaxis()->SetTitleSize(frame_et_rec->GetYaxis()->GetTitleSize() * k_p2bot_scale);
    frame_ratio_phenix->GetXaxis()->SetTitleOffset(frame_et_rec->GetXaxis()->GetTitleOffset() * 4.0 / 6.0 * 1.0);
    frame_ratio_phenix->GetYaxis()->SetTitleOffset(frame_et_rec->GetYaxis()->GetTitleOffset() * 4.0 / 6.0);
    frame_ratio_phenix->GetXaxis()->SetNdivisions(505);
    frame_ratio_phenix->Draw("axis");

    TGraphAsymmErrors *g_ratio_phenix       = new TGraphAsymmErrors();
    TGraphAsymmErrors *g_ratio_phenix_pstat = new TGraphAsymmErrors();
    TGraphAsymmErrors *g_ratio_phenix_dstat = new TGraphAsymmErrors();
    TGraphAsymmErrors *g_ratio_phenix_psys  = new TGraphAsymmErrors();
    TGraphAsymmErrors *g_data_sys_band      = new TGraphAsymmErrors();

    const double kStatOffset = 0.25;  // small horizontal offset to separate the two stat bars

    int rp_idx = 0;
    for (Int_t i = 1; i <= h_data->GetNbinsX(); ++i)
    {
        double pT_c = h_data->GetBinCenter(i);
        if (pT_c < lowerx || pT_c > 26.0) continue;
        double y_d = h_data->GetBinContent(i);
        double e_d = h_data->GetBinError(i);
        if (y_d <= 0.0) continue;

        // Match PHENIX-corrected bin
        double x_p = 0.0, y_p = 0.0;
        int p_idx = -1;
        for (Int_t j = 0; j < gStat_PHENIX_corr->GetN(); ++j)
        {
            gStat_PHENIX_corr->GetPoint(j, x_p, y_p);
            if (std::abs(x_p - pT_c) < 0.1) { p_idx = j; break; }
        }
        if (p_idx < 0 || y_p <= 0.0) continue;

        // Data sys from g_syst (same binning as h_data)
        double sys_d_lo = 0.0, sys_d_hi = 0.0;
        for (Int_t j = 0; j < g_syst->GetN(); ++j)
        {
            double x_s = 0.0, y_s = 0.0;
            g_syst->GetPoint(j, x_s, y_s);
            if (std::abs(x_s - pT_c) < 0.1) {
                sys_d_lo = g_syst->GetErrorYlow(j);
                sys_d_hi = g_syst->GetErrorYhigh(j);
                break;
            }
        }
        // PHENIX-corrected sys (same index as the matched stat point)
        double sys_p_lo = gSys_PHENIX_corr->GetErrorYlow(p_idx);
        double sys_p_hi = gSys_PHENIX_corr->GetErrorYhigh(p_idx);
        // PHENIX-corrected stat (asymmetric on gStat_PHENIX_corr)
        double stat_p_lo = gStat_PHENIX_corr->GetErrorYlow(p_idx);
        double stat_p_hi = gStat_PHENIX_corr->GetErrorYhigh(p_idx);

        double r          = y_p / y_d;
        // Per-source stat propagated to the ratio
        double r_pstat_lo = r * (stat_p_lo / y_p);
        double r_pstat_hi = r * (stat_p_hi / y_p);
        double r_dstat    = r * (e_d / y_d);
        double r_psys_lo  = sys_p_lo / y_d;
        double r_psys_hi  = sys_p_hi / y_d;
        double r_dsys_lo  = sys_d_lo / y_d;
        double r_dsys_hi  = sys_d_hi / y_d;
        double bin_w      = h_data->GetBinWidth(i);

        // Marker at the bin centre with bin-width horizontal bars only
        g_ratio_phenix->SetPoint(rp_idx, pT_c, r);
        g_ratio_phenix->SetPointError(rp_idx, bin_w / 2.0, bin_w / 2.0, 0.0, 0.0);

        // PHENIX stat: vertical bar at slight left offset
        g_ratio_phenix_pstat->SetPoint(rp_idx, pT_c - kStatOffset, r);
        g_ratio_phenix_pstat->SetPointError(rp_idx, 0.0, 0.0, r_pstat_lo, r_pstat_hi);

        // sPHENIX stat: vertical bar at slight right offset
        g_ratio_phenix_dstat->SetPoint(rp_idx, pT_c + kStatOffset, r);
        g_ratio_phenix_dstat->SetPointError(rp_idx, 0.0, 0.0, r_dstat, r_dstat);

        g_ratio_phenix_psys->SetPoint(rp_idx, pT_c, r);
        g_ratio_phenix_psys->SetPointError(rp_idx, bin_w / 2.0, bin_w / 2.0, r_psys_lo, r_psys_hi);

        g_data_sys_band->SetPoint(rp_idx, pT_c, 1.0);
        g_data_sys_band->SetPointError(rp_idx, bin_w / 2.0, bin_w / 2.0, r_dsys_lo, r_dsys_hi);

        ++rp_idx;
    }

    // Layer 1 (back): data sys band centred on y=1
    g_data_sys_band->SetFillColorAlpha(col[0], trans[0]);
    g_data_sys_band->SetLineColor(col[0]);
    g_data_sys_band->Draw("2 same");

    // Unity reference line
    lineone->SetLineColor(kBlack);
    lineone->SetLineStyle(2);
    lineone->SetLineWidth(2);
    lineone->Draw("L same");

    // Layer 2: PHENIX-corrected sys boxes around the ratio markers
    g_ratio_phenix_psys->SetMarkerStyle(PHENIXcorrMarker);
    g_ratio_phenix_psys->SetMarkerColor(kViolet + 1);
    g_ratio_phenix_psys->SetLineColor(kViolet + 1);
    g_ratio_phenix_psys->SetFillColorAlpha(kViolet + 1, 0.30);
    g_ratio_phenix_psys->Draw("2 same");

    // Layer 3a: PHENIX stat bar (violet, left of marker)
    g_ratio_phenix_pstat->SetMarkerStyle(0);
    g_ratio_phenix_pstat->SetLineColor(kViolet + 1);
    g_ratio_phenix_pstat->SetLineWidth(2);
    g_ratio_phenix_pstat->Draw("Z same");

    // Layer 3b: sPHENIX stat bar (azure, right of marker)
    g_ratio_phenix_dstat->SetMarkerStyle(0);
    g_ratio_phenix_dstat->SetLineColor(col[0]);
    g_ratio_phenix_dstat->SetLineWidth(2);
    g_ratio_phenix_dstat->Draw("Z same");

    // Layer 4 (front): markers with bin-width horizontal bars only
    g_ratio_phenix->SetMarkerStyle(PHENIXcorrMarker);
    g_ratio_phenix->SetMarkerSize(mkSize[1]);
    g_ratio_phenix->SetMarkerColor(kViolet + 1);
    g_ratio_phenix->SetLineColor(kViolet + 1);
    g_ratio_phenix->Draw("P same");

    // Paper Fig. final_phenix.pdf -- 2-pad PHENIX overlay + ratio.
    c2->SaveAs(Form("%s/final_phenix.pdf", paperdir.c_str()));

    //-----------------------------------------------------------------
    // PHENIX fit + pull diagnostic
    //-----------------------------------------------------------------
    TCanvas *c3 = new TCanvas("c3_phenix_fit", "", 800, 800);

    TPad *p3_top = new TPad("p3_top", "", 0, 0.32, 1, 1);
    p3_top->SetTopMargin(0.06);
    p3_top->SetBottomMargin(0.02);
    p3_top->SetLeftMargin(0.13);
    p3_top->SetRightMargin(0.05);
    p3_top->SetLogy();
    p3_top->Draw();

    TPad *p3_bot = new TPad("p3_bot", "", 0, 0, 1, 0.32);
    p3_bot->SetTopMargin(0.02);
    p3_bot->SetBottomMargin(0.30);
    p3_bot->SetLeftMargin(0.13);
    p3_bot->SetRightMargin(0.05);
    p3_bot->Draw();

    p3_top->cd();
    TH1F *frame_fit_top = new TH1F("frame_fit_top", "", 1, xLow[0], xHigh[n - 1]);
    frame_fit_top->SetYTitle("d^{2}#sigma/d#it{#eta}d#it{E}_{T}^{#gamma} [pb/GeV]");
    frame_fit_top->GetYaxis()->SetRangeUser(0.05, 1e5);
    frame_fit_top->GetYaxis()->SetTitleSize(0.055);
    frame_fit_top->GetYaxis()->SetTitleOffset(1.05);
    frame_fit_top->GetYaxis()->SetLabelSize(0.050);
    frame_fit_top->GetXaxis()->SetLabelSize(0);
    frame_fit_top->GetXaxis()->SetTickLength(0.03);
    frame_fit_top->Draw("axis");

    gFit_PHENIX->SetMarkerStyle(mkStyle[1]);
    gFit_PHENIX->SetMarkerSize(mkSize[1]);
    gFit_PHENIX->SetMarkerColor(col[1]);
    gFit_PHENIX->SetLineColor(col[1]);
    gFit_PHENIX->Draw("P same");

    f_phenix_mpl->SetLineColor(kBlack);
    f_phenix_mpl->SetLineWidth(2);
    f_phenix_mpl->Draw("L same");

    myText(0.55, 0.87, 1, strleg1.c_str(), 0.046, 0);
    myText(0.55, 0.81, 1, "PHENIX PRD 86 072008, |#eta^{#gamma}| < 0.25", 0.038, 0);
    myText(0.55, 0.76, 1,
        Form("Fit: A(1 + #it{p}_{T}^{2}/b)^{c}, range %.0f#minus%.0f GeV",
             xLow[0], xHigh[n - 1]),
        0.034, 0);
    myText(0.55, 0.71, 1,
        Form("A = %.3g, b = %.2f, c = %.3f",
             f_phenix_mpl->GetParameter(0),
             f_phenix_mpl->GetParameter(1),
             f_phenix_mpl->GetParameter(2)),
        0.034, 0);
    myText(0.55, 0.66, 1,
        Form("#chi^{2}/ndf = %.2f / %d",
             f_phenix_mpl->GetChisquare(),
             f_phenix_mpl->GetNDF()),
        0.034, 0);

    p3_bot->cd();
    TH1F *frame_fit_bot = new TH1F("frame_fit_bot", "", 1, xLow[0], xHigh[n - 1]);
    frame_fit_bot->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");
    frame_fit_bot->SetYTitle("(data #minus fit) / #sigma_{stat}");
    frame_fit_bot->GetYaxis()->SetRangeUser(-3.5, 3.5);
    frame_fit_bot->GetYaxis()->SetNdivisions(505);
    frame_fit_bot->GetYaxis()->SetTitleSize(0.085);
    frame_fit_bot->GetYaxis()->SetTitleOffset(0.65);
    frame_fit_bot->GetYaxis()->SetLabelSize(0.080);
    frame_fit_bot->GetXaxis()->SetTitleSize(0.095);
    frame_fit_bot->GetXaxis()->SetTitleOffset(1.20);
    frame_fit_bot->GetXaxis()->SetLabelSize(0.080);
    frame_fit_bot->Draw("axis");

    TGraph *g_pull = new TGraph(n);
    for (Int_t i = 0; i < n; ++i)
    {
        double pT = x[i];
        double sf = pT * factorCommon;
        double y_d = y[i] * sf;
        double s_d = 0.5 * (statUp[i] + statDown[i]) * sf;
        double y_binavg_p = f_phenix_mpl->Integral(xLow[i], xHigh[i]) / (xHigh[i] - xLow[i]);
        double pull = (s_d > 0.0) ? (y_d - y_binavg_p) / s_d : 0.0;
        g_pull->SetPoint(i, pT, pull);
    }
    g_pull->SetMarkerStyle(mkStyle[1]);
    g_pull->SetMarkerSize(mkSize[1]);
    g_pull->SetMarkerColor(col[1]);
    g_pull->SetLineColor(col[1]);
    g_pull->Draw("P same");

    TLine *l_pull0 = new TLine(xLow[0], 0.0, xHigh[n - 1], 0.0);
    l_pull0->SetLineStyle(2);
    l_pull0->SetLineColor(kGray + 2);
    l_pull0->Draw();

    c3->SaveAs(Form("%s/final_phenix_fit_%s.pdf", scratch.c_str(), tune.data()));

    //-----------------------------------------------------------------
    // Paper Fig. final_sphenix.pdf -- single-pad sPHENIX-only cross-section.
    //
    // Modernized from plotting/plot_final_backup250327_v1.C L877-933:
    //   * uses the current 64.37 pb^-1 lumi (from strleg_lumi)
    //   * drops the PHENIX overlay (data-only)
    //   * keeps the data systematic band + markers
    //
    // The PHENIX overlay is intentionally removed -- final_phenix.pdf
    // already serves the cross-experiment comparison.
    //-----------------------------------------------------------------
    init_plot();
    h_data = h_data_cp;

    TCanvas *c4 = new TCanvas("c4_paper_sphenix", "", 800, 700);
    c4->SetLogy();

    frame_et_rec->SetTitle(";#it{E}_{T}^{#gamma} [GeV];d^{2}#it{#sigma}/(d#it{#eta} d#it{E}_{T}^{#gamma}) [pb/GeV]");
    frame_et_rec->GetYaxis()->SetRangeUser(lowery, 500);
    frame_et_rec->GetXaxis()->SetRangeUser(lowerx, upperx);
    frame_et_rec->GetXaxis()->SetTitleOffset(1.0);
    frame_et_rec->GetYaxis()->SetTitleOffset(1.0);
    frame_et_rec->GetYaxis()->SetTitleSize(0.05);
    frame_et_rec->GetYaxis()->SetLabelSize(0.045);
    frame_et_rec->GetXaxis()->SetLabelSize(0.045);
    frame_et_rec->GetXaxis()->SetTitleSize(0.05);
    frame_et_rec->Draw("axis");

    g_syst->SetMarkerStyle(mkStyle[0]);
    g_syst->SetMarkerColor(col[0]);
    g_syst->SetLineColor(col[0]);
    g_syst->SetFillColorAlpha(col[0], trans[0]);
    g_syst->Draw("2 same");

    h_data->SetMarkerStyle(mkStyle[0]);
    h_data->SetMarkerSize(mkSize[0]);
    h_data->SetMarkerColor(col[0]);
    h_data->SetLineColor(col[0]);
    h_data->Draw("same");

    {
        const float xpos_s = 0.20, xpos2_s = 0.91, ypos_s = 0.87,
                    dy_s = 0.055, fs1 = 0.047, fs = 0.043;
        myText(xpos2_s, ypos_s - 0 * dy_s, 1, strleg1.c_str(),     fs1, 1);
        myText(xpos2_s, ypos_s - 1 * dy_s, 1, strleg_lumi.c_str(), fs,  1);
        myText(xpos2_s, ypos_s - 2 * dy_s, 1, strleg3.c_str(),     fs,  1);
        myText(xpos2_s, ypos_s - 3 * dy_s, 1, strleg4.c_str(),     fs,  1);

        TH1F *htemp_data = (TH1F *)h_data->Clone("htemp_data_sphenix");
        htemp_data->SetMarkerStyle(mkStyle[0]);
        htemp_data->SetMarkerSize(mkSize[0]);
        htemp_data->SetMarkerColor(col[0]);
        htemp_data->SetLineColor(col[0]);
        htemp_data->SetFillColorAlpha(col[0], trans[0]);

        TLegend *l_s = new TLegend(xpos_s, 0.25, 0.55, 0.34);
        legStyle(l_s, 0.21, fs);
        l_s->AddEntry(htemp_data, "Data", "fpl");
        l_s->Draw("same");
    }

    c4->SaveAs(Form("%s/final_sphenix.pdf", paperdir.c_str()));
}
