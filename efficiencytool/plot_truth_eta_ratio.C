// plot_truth_eta_ratio.C
//
// Pythia truth-level eta-density ratio for prompt photons:
//
//   R(pT) = ( N(|eta|<0.25) / 0.5 ) / ( N(|eta|<0.7) / 1.4 )
//
// Fills numerator/denominator from the same events on photon5/10/20, so the
// xs weight cancels per pT bin and the ratio is purely the Pythia eta shape.
// Numerator window |eta|<0.25 is contained in the denominator window |eta|<0.7,
// so bin-by-bin Sumw2 errors are correlated; treat the error bars as upper bounds.
//
// Usage (from efficiencytool/):
//   root -l -b -q 'plot_truth_eta_ratio.C("config_bdt_bdtmodel_v0_0rad.yaml")'
//   root -l -b -q 'plot_truth_eta_ratio.C("config_bdt_bdtmodel_v0_0rad.yaml", 4.0, 3, "truth_eta_ratio_iso4")'
//   root -l -b -q 'plot_truth_eta_ratio.C("config_bdt_bdtmodel_v0_0rad.yaml", -1.0, 3, "truth_eta_ratio_inclusive")'
//
// Args:
//   configname       : YAML config (used only for input paths + tree name)
//   truth_iso_max    : truth iso (R=0.3) ET cut in GeV. Set < 0 to disable.
//                      Default 4.0 matches the analysis fiducial cross section.
//   prompt_class_max : keep particles with photonclass < this value.
//                      3 = direct + fragmentation (default); 2 = direct only.
//   outname          : prefix for .pdf/.png/.root output

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <yaml-cpp/yaml.h>

#include "CrossSectionWeights.h"

void plot_truth_eta_ratio(const std::string &configname = "config_bdt_bdtmodel_v0_0rad.yaml",
                          float truth_iso_max  = 4.0f,
                          int   prompt_class_max = 3,
                          const std::string &outname = "truth_eta_ratio")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    const std::string root_dir   = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
    const std::string branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
    const std::string treename   = configYaml["input"]["tree"].as<std::string>();

    const float dEta_central = 2.0f * 0.25f;  // 0.5
    const float dEta_full    = 2.0f * 0.7f;   // 1.4

    // Truth pT binning. Matches analysis pT_bins_truth from config.
    std::vector<double> pt_edges = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45};
    const int nbins = (int)pt_edges.size() - 1;

    auto make_hist = [&](const char *name, const char *title) {
        TH1D *h = new TH1D(name, title, nbins, pt_edges.data());
        h->Sumw2();
        return h;
    };
    TH1D *h_central = make_hist("h_pt_eta025",
        ";truth photon p_{T} [GeV];weighted N (|#eta|<0.25)");
    TH1D *h_full    = make_hist("h_pt_eta07",
        ";truth photon p_{T} [GeV];weighted N (|#eta|<0.7)");

    // Diagnostic per-sample copies (for verification)
    std::vector<std::pair<std::string, std::pair<float,float>>> samples = {
        {"photon5",  {0.0f,  14.0f}},
        {"photon10", {14.0f, 22.0f}},
        {"photon20", {22.0f, 1e6f}},
    };

    for (const auto &sp : samples) {
        const std::string &sample = sp.first;
        const float pt_lo = sp.second.first;
        const float pt_hi = sp.second.second;

        const std::string fname = root_dir + sample + branch_dir;
        std::cout << "[truth_eta_ratio] Loading " << fname << std::endl;

        TChain chain(treename.c_str());
        if (chain.Add(fname.c_str()) == 0) {
            std::cerr << "  no entries; skipping " << sample << std::endl;
            continue;
        }

        const PPG12::SampleConfig sc = PPG12::GetSampleConfig(sample);
        const float xs_weight = sc.weight;

        TTreeReader reader(&chain);
        TTreeReaderValue<int>     nparticles(reader, "nparticles");
        TTreeReaderArray<float>   particle_Pt(reader, "particle_Pt");
        TTreeReaderArray<float>   particle_Eta(reader, "particle_Eta");
        TTreeReaderArray<int>     particle_pid(reader, "particle_pid");
        TTreeReaderArray<int>     particle_photonclass(reader, "particle_photonclass");
        TTreeReaderArray<float>   particle_truth_iso_03(reader, "particle_truth_iso_03");

        const Long64_t nentries = chain.GetEntries();
        std::cout << "  entries: " << nentries
                  << "  xs_weight: " << xs_weight
                  << "  truth pT window: [" << pt_lo << ", " << pt_hi << ")" << std::endl;

        Long64_t ient = 0;
        Long64_t accepted_events = 0;
        while (reader.Next()) {
            if ((++ient) % 200000 == 0) {
                std::cout << "    " << sample << "  " << ient << " / " << nentries << std::endl;
            }

            // Gate event into this sample's truth pT window using leading photon pT
            // (matches RecoEffCalculator_TTreeReader.C convention; prevents double-
            // counting between photon5/10/20 in their canonical windows).
            float max_phpT = 0.0f;
            for (int ip = 0; ip < *nparticles; ++ip) {
                if (particle_pid[ip] == 22 && particle_Pt[ip] > max_phpT) {
                    max_phpT = particle_Pt[ip];
                }
            }
            if (max_phpT < pt_lo || max_phpT >= pt_hi) continue;
            ++accepted_events;

            for (int ip = 0; ip < *nparticles; ++ip) {
                if (particle_pid[ip] != 22) continue;
                if (particle_photonclass[ip] >= prompt_class_max) continue;
                if (truth_iso_max >= 0.0f && particle_truth_iso_03[ip] >= truth_iso_max) continue;

                const float pT   = particle_Pt[ip];
                const float aeta = std::fabs(particle_Eta[ip]);

                if (aeta < 0.7f)  h_full   ->Fill(pT, xs_weight);
                if (aeta < 0.25f) h_central->Fill(pT, xs_weight);
            }
        }
        std::cout << "  accepted events in pT window: " << accepted_events << std::endl;
    }

    // Density per dEta and ratio
    TH1D *h_central_density = (TH1D*)h_central->Clone("h_pt_eta025_density");
    TH1D *h_full_density    = (TH1D*)h_full   ->Clone("h_pt_eta07_density");
    h_central_density->Scale(1.0 / dEta_central);
    h_full_density   ->Scale(1.0 / dEta_full);

    TH1D *h_ratio = (TH1D*)h_central_density->Clone("h_ratio_central_over_full");
    h_ratio->Divide(h_full_density);
    h_ratio->SetTitle(";truth photon p_{T} [GeV];"
        "(dN/d#eta)_{|#eta|<0.25} / (dN/d#eta)_{|#eta|<0.7}");

    // Plot
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c_truth_eta_ratio", "", 800, 700);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.13);
    c->SetTopMargin(0.06);
    c->SetRightMargin(0.05);

    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerSize(1.0);
    h_ratio->SetMarkerColor(kBlack);
    h_ratio->SetLineColor(kBlack);
    h_ratio->GetYaxis()->SetTitleOffset(1.5);
    h_ratio->GetXaxis()->SetTitleOffset(1.1);
    h_ratio->GetYaxis()->SetRangeUser(0.7, 1.3);
    h_ratio->GetXaxis()->SetRangeUser(8, 45);
    h_ratio->Draw("E");

    TLine *l1 = new TLine(8, 1.0, 45, 1.0);
    l1->SetLineColor(kGray + 2);
    l1->SetLineStyle(2);
    l1->Draw();

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.034);
    tex.DrawLatex(0.18, 0.88, "#bf{#it{sPHENIX}} Internal");
    tex.DrawLatex(0.18, 0.84, "Pythia 8 truth, prompt #gamma");
    if (prompt_class_max == 2) {
        tex.DrawLatex(0.18, 0.80, "(direct only)");
    } else {
        tex.DrawLatex(0.18, 0.80, "(direct + fragmentation)");
    }
    if (truth_iso_max >= 0.0f) {
        tex.DrawLatex(0.18, 0.76, Form("Truth iso (R=0.3) < %.1f GeV", truth_iso_max));
    } else {
        tex.DrawLatex(0.18, 0.76, "No isolation requirement");
    }

    const std::string pdf  = outname + ".pdf";
    const std::string png  = outname + ".png";
    const std::string root = outname + ".root";
    c->SaveAs(pdf.c_str());
    c->SaveAs(png.c_str());

    TFile fout(root.c_str(), "RECREATE");
    h_central->Write();
    h_full->Write();
    h_central_density->Write();
    h_full_density->Write();
    h_ratio->Write();
    fout.Close();

    std::cout << "[truth_eta_ratio] saved: " << pdf << ", " << png << ", " << root << std::endl;

    // Print bin-by-bin numbers for sanity
    std::cout << "\npT [GeV)        N(|eta|<0.25)    N(|eta|<0.7)     ratio" << std::endl;
    for (int ib = 1; ib <= h_ratio->GetNbinsX(); ++ib) {
        const double xlo = h_ratio->GetXaxis()->GetBinLowEdge(ib);
        const double xhi = h_ratio->GetXaxis()->GetBinUpEdge(ib);
        const double nC = h_central->GetBinContent(ib);
        const double nF = h_full   ->GetBinContent(ib);
        const double R  = h_ratio  ->GetBinContent(ib);
        const double Re = h_ratio  ->GetBinError(ib);
        printf("[%5.1f, %5.1f)  %14.6g  %14.6g  %6.3f +/- %5.3f\n",
               xlo, xhi, nC, nF, R, Re);
    }
}
