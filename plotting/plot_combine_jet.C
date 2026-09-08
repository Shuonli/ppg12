// plot_combine_jet.C
//
// Companion to plot_combine.C — same combining strategy applied to jet
// MC samples. Reads h_max_truth_jet_pT from the 5 jet samples used in
// the nominal SI/DI pipeline (jet8/12/20/30/40), each pre-weighted by
// effective cross section at fill time, and overlays them with their
// stitched sum and a smooth power-law fit.
//
// Output: figures/combine_jet.pdf

#include "plotcommon.h"

void plot_combine_jet()
{
    init_plot();

    const std::string base =
        "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/";

    // Per-sample truth-leading-jet pT window (matches CrossSectionWeights.h).
    // jet8 and jet30 have a sub-permille float-precision leak one fine bin
    // above their upper edge that we trim before plotting; otherwise the
    // leaked bin shows up as a stray low-statistics point on the log canvas
    // 5 orders of magnitude below the curve. The cross-section weights are
    // applied at fill time in RecoEffCalculator_TTreeReader.C; this plotter
    // just stacks the per-sample histograms unchanged.
    // June 2026 production: per-sample MC exists as the four merge-feeder
    // outputs (SI "_nom" and DI "_double", 0 mrad and 1.5 mrad). Each is
    // pre-scaled by lumi_weight and mix_weight at fill time, so the plain sum
    // of the four is the all-range SI/DI-blended MC used in the analysis.
    struct Sample {
        std::string name;
        int color;
        float pt_lo;
        float pt_hi;       // events with ptj >= pt_hi rejected
    };
    std::vector<Sample> samples = {
        {"jet8",  kPink + 5,     9.0f, 14.0f},
        {"jet12", kGreen - 2,   14.0f, 21.0f},
        {"jet20", kAzure + 7,   21.0f, 32.0f},
        {"jet30", kOrange + 7,  32.0f, 42.0f},
        {"jet40", kMagenta + 1, 42.0f, 1e6f},
    };
    const std::vector<std::string> parts = {"_nom_bdt_nom_0rad", "_nom_bdt_nom_1p5mrad",
                                           "_double_bdt_nom_0rad", "_double_bdt_nom_1p5mrad"};

    int rebinx = 10;

    std::vector<TH1F*> h_per_sample;
    for (auto &s : samples) {
        TH1F *h = nullptr;
        for (const auto &part : parts) {
            const std::string file = base + "MC_efficiency_" + s.name + part + ".root";
            TFile *f = TFile::Open(file.c_str(), "READ");
            if (!f || f->IsZombie()) { printf("missing %s\n", file.c_str()); return; }
            TH1F *hp = dynamic_cast<TH1F*>(f->Get("h_max_truth_jet_pT"));
            if (!hp) { printf("missing h_max_truth_jet_pT in %s\n", file.c_str()); return; }
            if (!h) { h = (TH1F*) hp->Clone(("h_max_truth_jet_pT_" + s.name).c_str()); h->SetDirectory(nullptr); }
            else h->Add(hp);
        }

        // Zero out bins whose center falls outside the declared truth-pT window
        // (handles float-precision leak at the upper boundary).
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            const double xc = h->GetXaxis()->GetBinCenter(i);
            if (xc < s.pt_lo || xc >= s.pt_hi) {
                h->SetBinContent(i, 0.0);
                h->SetBinError(i, 0.0);
            }
        }

        h->Rebin(rebinx);
        h->SetDirectory(nullptr);
        h_per_sample.push_back(h);
    }

    TH1F *h_sum = (TH1F*) h_per_sample[0]->Clone("h_max_truth_jet_pT_sum");
    for (size_t i = 1; i < h_per_sample.size(); ++i) h_sum->Add(h_per_sample[i]);

    // Power-law fit on the sum, range matched to the analysis truth window
    float xlower = 9.0f;
    float xupper = 60.0f;
    TF1 *f_fit = new TF1("f_fit_jet",
                         "[0]*pow([1]/x,[2]+[3]*log(x/[1])+ [4]*x)",
                         xlower, xupper);
    f_fit->SetParameters(1e9, 1.0, 1.0, 2.0, 0.01);
    h_sum->Fit(f_fit, "REMN", "", xlower, xupper);
    h_sum->Fit(f_fit, "REMN", "", xlower, xupper);

    TH1F *h_ratio = (TH1F*) h_sum->Clone("h_max_truth_jet_pT_sum_overfit");
    h_ratio->Divide(f_fit);

    // 2-pad canvas (top: stack + fit, bottom: sim/fit ratio)
    TCanvas *c = new TCanvas("c_combine_jet", "", 800, 889);
    c->Divide(1, 2);

    TPad *pad_1 = (TPad *) c->cd(1);
    pad_1->SetPad(0, 0.4, 1, 1);
    pad_1->SetTopMargin(0.05);
    pad_1->SetLeftMargin(0.13);
    pad_1->SetBottomMargin(0.002);
    pad_1->SetRightMargin(0.08);
    pad_1->SetLogy();

    frame_et_rec->SetYTitle("counts");
    frame_et_rec->GetYaxis()->SetRangeUser(5e3, 1e12);
    frame_et_rec->GetXaxis()->SetRangeUser(8, 60);
    frame_et_rec->GetXaxis()->SetTitleOffset(1.05);
    frame_et_rec->GetYaxis()->SetTitleOffset(1.05);
    frame_et_rec->GetYaxis()->SetTitleSize(0.053);
    frame_et_rec->GetXaxis()->SetLabelSize(0.050);
    frame_et_rec->GetYaxis()->SetLabelSize(0.050);
    frame_et_rec->GetXaxis()->SetLabelOffset(2);
    frame_et_rec->GetXaxis()->SetNdivisions(505);
    frame_et_rec->Draw("axis");

    for (size_t i = 0; i < samples.size(); ++i) {
        h_per_sample[i]->SetMarkerStyle(20);
        h_per_sample[i]->SetMarkerColor(samples[i].color);
        h_per_sample[i]->SetLineColor(samples[i].color);
        h_per_sample[i]->Draw("same p");
    }

    f_fit->SetLineColor(kRed);
    f_fit->Draw("same");

    myText(0.5, 0.90, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, strMC.c_str(),   0.05);

    {
        TLegend *leg = new TLegend(0.55, 0.45, 0.92, 0.78);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.045);
        for (size_t i = 0; i < samples.size(); ++i)
            leg->AddEntry(h_per_sample[i], samples[i].name.c_str(), "lep");
        leg->Draw();
    }

    TPad *pad_2 = (TPad *) c->cd(2);
    pad_2->SetPad(0, 0, 1, 0.4);
    pad_2->SetTopMargin(0.023);
    pad_2->SetLeftMargin(0.13);
    pad_2->SetBottomMargin(0.25);
    pad_2->SetRightMargin(0.08);

    frame_et_truth->SetYTitle("MC / Fit");
    frame_et_truth->SetXTitle("Leading #it{p}_{T}^{jet} [GeV]");
    frame_et_truth->GetYaxis()->SetNdivisions(506);
    frame_et_truth->GetYaxis()->SetRangeUser(0.85, 1.15);
    frame_et_truth->GetXaxis()->SetRangeUser(8, 60);
    frame_et_truth->GetXaxis()->SetTitleOffset(frame_et_rec->GetXaxis()->GetTitleOffset() * 4 / 6.* 1.4);
    frame_et_truth->GetYaxis()->SetTitleOffset(frame_et_rec->GetYaxis()->GetTitleOffset() * 4 / 6.);
    frame_et_truth->GetYaxis()->SetLabelOffset(frame_et_rec->GetYaxis()->GetLabelOffset() * 4 / 6.);
    frame_et_truth->GetXaxis()->SetLabelSize(frame_et_rec->GetXaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetYaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetTitleSize(frame_et_rec->GetXaxis()->GetTitleSize() * 6 / 4. * 1.2);
    frame_et_truth->GetYaxis()->SetTitleSize(frame_et_rec->GetYaxis()->GetTitleSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetNdivisions(505);
    frame_et_truth->Draw("axis");

    h_ratio->Draw("same");
    lineone->Draw("L");

    c->SaveAs("figures/combine_jet.pdf");
    printf("saved figures/combine_jet.pdf\n");
}
