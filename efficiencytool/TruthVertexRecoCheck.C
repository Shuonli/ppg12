// TruthVertexRecoCheck.C
//
// Truth vertex-z distribution split by reco-vertex status
// (vertexz == -9999 sentinel => no MBD vertex reconstructed).
//
// Samples: photon10, photon10_double, jet12, jet12_double
// Output : results/truth_vertex_reco_check.root
//          reports/figures/truth_vertex_reco_check.pdf
//
// Usage: root -l -b -q efficiencytool/TruthVertexRecoCheck.C+

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPad.h>
#include <TLine.h>
#include <TPaveText.h>
#include <iostream>
#include <vector>
#include <string>

namespace {

struct Sample {
    std::string key;      // photon10, photon10_double, ...
    std::string label;    // display label
    std::string file;     // combined.root path
    int colorReco;
    int colorNoReco;
};

const std::vector<Sample> kSamples = {
    {"photon10",        "PYTHIA #gamma+jet (single)",
     "anatreemaker/macro_maketree/sim/run28/photon10/condorout/combined.root",
     kBlue+1,  kRed+1},
    {"photon10_double", "PYTHIA #gamma+jet (double)",
     "anatreemaker/macro_maketree/sim/run28/photon10_double/condorout/combined.root",
     kBlue+1,  kRed+1},
    {"jet12",           "PYTHIA inclusive jet (single)",
     "anatreemaker/macro_maketree/sim/run28/jet12/condorout/combined.root",
     kBlue+1,  kRed+1},
    {"jet12_double",    "PYTHIA inclusive jet (double)",
     "anatreemaker/macro_maketree/sim/run28/jet12_double/condorout/combined.root",
     kBlue+1,  kRed+1},
};

// Minimal sPHENIX-like style (since sPhenixStyle.C is not available locally).
void applyStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetLegendFont(42);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetEndErrorSize(0);
}

void drawHeader(double x, double y, const std::string& sampleLabel) {
    TLatex t;
    t.SetNDC();
    t.SetTextFont(62);
    t.SetTextSize(0.055);
    t.DrawLatex(x, y,        "#bf{#it{sPHENIX}} Internal");
    t.SetTextFont(42);
    t.SetTextSize(0.045);
    t.DrawLatex(x, y-0.065,  "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV, PYTHIA8");
    t.DrawLatex(x, y-0.12,   sampleLabel.c_str());
}

} // namespace

void TruthVertexRecoCheck() {
    applyStyle();

    const int nbins = 80;
    const double xmin = -200.0;
    const double xmax =  200.0;

    TFile* fout = TFile::Open("results/truth_vertex_reco_check.root", "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: cannot open results/truth_vertex_reco_check.root\n";
        return;
    }

    struct HistPair {
        TH1D* hReco;     // truth vz for events WITH reco vertex
        TH1D* hNoReco;   // truth vz for events WITHOUT reco vertex
        long long nReco;
        long long nNoReco;
    };
    std::vector<HistPair> hpairs(kSamples.size());

    for (size_t i = 0; i < kSamples.size(); ++i) {
        const auto& s = kSamples[i];
        std::cout << "[TruthVertexRecoCheck] processing " << s.key << "\n";

        TChain chain("slimtree");
        int added = chain.Add(s.file.c_str());
        if (added <= 0) {
            std::cerr << "  WARN: could not add " << s.file << " (added=" << added << ")\n";
            continue;
        }

        std::string hnRec = "h_vtxtruth_reco_"   + s.key;
        std::string hnNoR = "h_vtxtruth_noreco_" + s.key;
        hpairs[i].hReco   = new TH1D(hnRec.c_str(), "", nbins, xmin, xmax);
        hpairs[i].hNoReco = new TH1D(hnNoR.c_str(), "", nbins, xmin, xmax);
        hpairs[i].hReco->Sumw2();
        hpairs[i].hNoReco->Sumw2();

        // Fast branch-enable
        chain.SetBranchStatus("*", 0);
        chain.SetBranchStatus("vertexz", 1);
        chain.SetBranchStatus("vertexz_truth", 1);

        float vz = 0.f, vzt = 0.f;
        chain.SetBranchAddress("vertexz",       &vz);
        chain.SetBranchAddress("vertexz_truth", &vzt);

        const Long64_t N = chain.GetEntries();
        std::cout << "    entries = " << N << "\n";

        long long nR = 0, nNR = 0;
        const Long64_t report_every = std::max<Long64_t>(1, N/10);
        for (Long64_t ie = 0; ie < N; ++ie) {
            chain.GetEntry(ie);
            if (ie % report_every == 0) {
                std::cout << "    " << ie << " / " << N << "\n";
            }
            // "no reco vertex" is defined as the sentinel -9999 written in CaloAna24.cc
            // when no MBD-type GlobalVertex is found.
            const bool has_reco = (vz > -999.f);
            if (has_reco) {
                hpairs[i].hReco->Fill(vzt);
                ++nR;
            } else {
                hpairs[i].hNoReco->Fill(vzt);
                ++nNR;
            }
        }
        hpairs[i].nReco   = nR;
        hpairs[i].nNoReco = nNR;
        std::cout << "    with reco vtx: " << nR  << "  (" << (100.0*nR /(double)(nR+nNR)) << "%)\n";
        std::cout << "    no   reco vtx: " << nNR << "  (" << (100.0*nNR/(double)(nR+nNR)) << "%)\n";

        fout->cd();
        hpairs[i].hReco->Write();
        hpairs[i].hNoReco->Write();
    }

    // ---- Canvas 1: 2x2 shape comparison (area-normalized truth vz)
    TCanvas* c1 = new TCanvas("c_shape", "shape", 1400, 1100);
    c1->Divide(2, 2, 0.001, 0.001);

    for (size_t i = 0; i < kSamples.size(); ++i) {
        c1->cd(i+1);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.07);
        gPad->SetBottomMargin(0.13);

        auto* hR = hpairs[i].hReco;
        auto* hN = hpairs[i].hNoReco;
        if (!hR || !hN) continue;

        // Clone and area-normalize for shape comparison
        TH1D* hRn = (TH1D*)hR->Clone((std::string(hR->GetName())+"_norm").c_str());
        TH1D* hNn = (TH1D*)hN->Clone((std::string(hN->GetName())+"_norm").c_str());
        hRn->SetDirectory(nullptr);
        hNn->SetDirectory(nullptr);
        if (hRn->Integral() > 0) hRn->Scale(1.0 / hRn->Integral());
        if (hNn->Integral() > 0) hNn->Scale(1.0 / hNn->Integral());

        hRn->SetLineColor(kSamples[i].colorReco);
        hRn->SetLineWidth(3);
        hRn->SetMarkerStyle(20);
        hRn->SetMarkerColor(kSamples[i].colorReco);
        hRn->SetMarkerSize(0.8);

        hNn->SetLineColor(kSamples[i].colorNoReco);
        hNn->SetLineWidth(3);
        hNn->SetLineStyle(2);
        hNn->SetMarkerStyle(24);
        hNn->SetMarkerColor(kSamples[i].colorNoReco);
        hNn->SetMarkerSize(0.8);

        hRn->GetXaxis()->SetTitle("truth vertex z [cm]");
        hRn->GetYaxis()->SetTitle("normalized counts");
        hRn->GetYaxis()->SetTitleOffset(1.35);
        double ymax = std::max(hRn->GetMaximum(), hNn->GetMaximum());
        hRn->SetMaximum(1.45 * ymax);
        hRn->SetMinimum(0.0);

        hRn->Draw("HIST");
        hNn->Draw("HIST SAME");

        // 30 cm analysis cut reference lines
        TLine* l1 = new TLine(-30, 0, -30, 1.45*ymax);
        TLine* l2 = new TLine( 30, 0,  30, 1.45*ymax);
        l1->SetLineStyle(3); l1->SetLineColor(kGray+2); l1->Draw();
        l2->SetLineStyle(3); l2->SetLineColor(kGray+2); l2->Draw();

        // Legend with fractions
        long long ntot = hpairs[i].nReco + hpairs[i].nNoReco;
        double fR = ntot > 0 ? 100.0*hpairs[i].nReco/(double)ntot : 0.0;
        double fN = ntot > 0 ? 100.0*hpairs[i].nNoReco/(double)ntot : 0.0;

        TLegend* leg = new TLegend(0.50, 0.72, 0.94, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.040);
        leg->AddEntry(hRn, Form("reco vertex (%.1f%%)",    fR), "l");
        leg->AddEntry(hNn, Form("no reco vertex (%.1f%%)", fN), "l");
        leg->Draw();

        drawHeader(0.16, 0.88, kSamples[i].label);
    }

    c1->SaveAs("reports/figures/truth_vertex_reco_check.pdf");
    fout->cd();
    c1->Write("c_shape");

    // ---- Canvas 2: 2x2 RATIO (no-reco / reco, both area-normalized)
    TCanvas* c2 = new TCanvas("c_ratio", "ratio", 1400, 1100);
    c2->Divide(2, 2, 0.001, 0.001);
    for (size_t i = 0; i < kSamples.size(); ++i) {
        c2->cd(i+1);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.07);
        gPad->SetBottomMargin(0.13);

        auto* hR = hpairs[i].hReco;
        auto* hN = hpairs[i].hNoReco;
        if (!hR || !hN) continue;

        TH1D* hRn = (TH1D*)hR->Clone((std::string(hR->GetName())+"_rn").c_str());
        TH1D* hNn = (TH1D*)hN->Clone((std::string(hN->GetName())+"_rn").c_str());
        hRn->SetDirectory(nullptr);
        hNn->SetDirectory(nullptr);
        if (hRn->Integral() > 0) hRn->Scale(1.0 / hRn->Integral());
        if (hNn->Integral() > 0) hNn->Scale(1.0 / hNn->Integral());

        TH1D* hRatio = (TH1D*)hNn->Clone(("h_ratio_" + kSamples[i].key).c_str());
        hRatio->SetDirectory(nullptr);
        hRatio->Divide(hRn);
        hRatio->SetLineColor(kBlack);
        hRatio->SetMarkerColor(kBlack);
        hRatio->SetMarkerStyle(20);
        hRatio->SetMarkerSize(0.9);
        hRatio->GetXaxis()->SetTitle("truth vertex z [cm]");
        hRatio->GetYaxis()->SetTitle("(no reco) / (reco)   shape");
        hRatio->GetYaxis()->SetTitleOffset(1.35);
        hRatio->SetMinimum(0.0);
        hRatio->SetMaximum(3.5);
        hRatio->Draw("E1");

        TLine* l1 = new TLine(xmin, 1.0, xmax, 1.0);
        l1->SetLineStyle(2); l1->SetLineColor(kGray+2); l1->Draw();

        drawHeader(0.16, 0.88, kSamples[i].label);
    }
    c2->SaveAs("reports/figures/truth_vertex_reco_check_ratio.pdf");
    fout->cd();
    c2->Write("c_ratio");

    // ---- Summary text file: per-sample yields and fractions
    std::cout << "\n====== summary ======\n";
    std::cout << "sample                N_total       N_reco       N_noReco    frac_noReco\n";
    for (size_t i = 0; i < kSamples.size(); ++i) {
        long long ntot = hpairs[i].nReco + hpairs[i].nNoReco;
        double f = ntot > 0 ? 100.0*hpairs[i].nNoReco/(double)ntot : 0.0;
        printf("%-22s %10lld  %10lld  %10lld    %6.2f%%\n",
               kSamples[i].key.c_str(), ntot, hpairs[i].nReco, hpairs[i].nNoReco, f);
    }

    fout->Close();
    delete fout;
    std::cout << "[TruthVertexRecoCheck] wrote results/truth_vertex_reco_check.root\n";
    std::cout << "[TruthVertexRecoCheck] wrote reports/figures/truth_vertex_reco_check.pdf\n";
    std::cout << "[TruthVertexRecoCheck] wrote reports/figures/truth_vertex_reco_check_ratio.pdf\n";
}
