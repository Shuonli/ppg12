#include <iostream>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TSystem.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <yaml-cpp/yaml.h>
// unfolding
#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>

#include "draw.C"

void Closure(const std::string &configname = "config.yaml")
{
    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);
    
    // std::string infilename = configYaml["input"]["photon_jet_file"].as<std::string>();

    std::string infilename = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_response_photon10.root";

    std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();

    TFile *fresin = new TFile(infilename.c_str(), "READ");

    // unfold response matrix
    std::vector<RooUnfoldResponse *> responses_full;
    std::vector<RooUnfoldResponse *> responses_half;
    // id histogram for unfolding
    std::vector<TH1D *> h_pT_truth_response;
    std::vector<TH1D *> h_pT_reco_response;

    std::vector<TH1D *> h_pT_truth_half_response;
    std::vector<TH1D *> h_pT_reco_half_response;

    std::vector<TH1D *> h_pT_truth_secondhalf_response;
    std::vector<TH1D *> h_pT_reco_secondhalf_response;

    std::vector<int> colors = {kPink + 8, kSpring - 7, kAzure - 3, kViolet + 3, kOrange + 10};

    std::vector<int> markerstyle = {20, 20, 20, 20, 20};

    for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
    {
        responses_full.push_back((RooUnfoldResponse *)fresin->Get(Form("response_matrix_full_%d", ieta)));
        responses_half.push_back((RooUnfoldResponse *)fresin->Get(Form("response_matrix_half_%d", ieta)));

        h_pT_truth_response.push_back((TH1D *)fresin->Get(Form("h_pT_truth_response_%d", ieta)));
        h_pT_reco_response.push_back((TH1D *)fresin->Get(Form("h_pT_reco_response_%d", ieta)));

        h_pT_truth_half_response.push_back((TH1D *)fresin->Get(Form("h_pT_truth_half_response_%d", ieta)));
        h_pT_reco_half_response.push_back((TH1D *)fresin->Get(Form("h_pT_reco_half_response_%d", ieta)));

        h_pT_truth_secondhalf_response.push_back((TH1D *)fresin->Get(Form("h_pT_truth_secondhalf_response_%d", ieta)));
        h_pT_reco_secondhalf_response.push_back((TH1D *)fresin->Get(Form("h_pT_reco_secondhalf_response_%d", ieta)));

        // full closure test
        int niterations = 1;
        for (int iit = 0; iit < niterations; iit++)
        {
            RooUnfoldBayes *full_unfold = new RooUnfoldBayes(responses_full[ieta], h_pT_reco_response[ieta], iit + 1);
            TH1D *hRecoPT = (TH1D *)full_unfold->Hunfold();

            hRecoPT->SetName("hRecoPT");

            std::cout << "full closure test iteration: " << iit << std::endl;

            std::cout << "hRecoPT: " << hRecoPT->GetBinContent(2) << std::endl;

            std::vector<TH1F *> h_input;

            h_input.push_back((TH1F *)h_pT_truth_response[ieta]);
            h_input.push_back((TH1F *)h_pT_reco_response[ieta]);
            h_input.push_back((TH1F *)hRecoPT);

            std::vector<std::string> text;
            std::vector<std::string> legend;

            text.push_back("Full closure test");
            text.push_back(Form("Iteration: %d", iit + 1));
            text.push_back(Form("%.1f < #eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]));
            text.push_back("|vtxz|<30cm, tight iso truth #gamma");

            legend.push_back("Truth photon spectrum");
            legend.push_back("Reco photon spectrum");
            legend.push_back("Unfolded photon spectrum");

            draw_1D_multiple_plot_ratio(h_input, colors, markerstyle,
                                        false, 10, true,
                                        false, 0, 40, false,
                                        false, 0, 0.5, true,
                                        true, 0, 2,
                                        true, "p_{T}^{jet} [GeV]", "dN_{#gamma}/dp_{T}^{#gamma}", "Ratio",
                                        false, "Run 22 Photon10GeV sim",
                                        true, text, 0.25, 0.41, 0.04,
                                        true, legend, 0.58, 0.76, 0.04,
                                        "figure/h_unfold_reco_full_1.png");
        }

        //half closure test
    }
}