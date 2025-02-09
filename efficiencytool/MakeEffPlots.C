#include "/sphenix/u/bseidlitz/plotstyle/AtlasUtils.C"
#include "/sphenix/u/shuhang98/AtlasStyle.C"
#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TH2.h>
#include <yaml-cpp/yaml.h>

void MakeEffPlots(const std::string &configname = "config.yaml")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    std::string infilename = configYaml["output"]["eff_outfile"].as<std::string>();

    float truthisocut = configYaml["analysis"]["truth_iso_max"].as<float>();

    float recoiso_min = configYaml["analysis"]["reco_iso_min"].as<float>();
    float recoiso_max = configYaml["analysis"]["reco_iso_max"].as<float>();

    float vertexcut = configYaml["analysis"]["vertex_cut"].as<float>();
    std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();

    int conesize = configYaml["analysis"]["cone_size"].as<int>();

    float reco_min_ET = configYaml["analysis"]["reco_min_ET"].as<float>();

    float eff_dR = configYaml["analysis"]["eff_dR"].as<float>();

    // getting cuts from the config file
    float tight_reta77_min = configYaml["analysis"]["tight"]["reta77_min"].as<float>();
    float tight_reta77_max = configYaml["analysis"]["tight"]["reta77_max"].as<float>();

    float tight_rhad33_max = configYaml["analysis"]["tight"]["rhad33_max"].as<float>();
    float tight_rhad33_min = configYaml["analysis"]["tight"]["rhad33_min"].as<float>();

    float tight_w72_max = configYaml["analysis"]["tight"]["w72_max"].as<float>();
    float tight_w72_min = configYaml["analysis"]["tight"]["w72_min"].as<float>();

    float tight_re11_E_max = configYaml["analysis"]["tight"]["re11_E_max"].as<float>();
    float tight_re11_E_min = configYaml["analysis"]["tight"]["re11_E_min"].as<float>();

    TFile *fin = new TFile(infilename.c_str(), "READ");

    std::vector<TEfficiency *> eff_reco_eta;
    std::vector<TEfficiency *> eff_id_eta;
    std::vector<TEfficiency *> eff_converts_eta;
    std::vector<TEfficiency *> eff_iso_eta;


    for (int ieta = 0; ieta < (int)eta_bins.size() - 4; ieta++)
    {
        //get file
        eff_reco_eta.push_back((TEfficiency *)fin->Get(Form("eff_reco_eta_%d", ieta)));
        eff_id_eta.push_back((TEfficiency *)fin->Get(Form("eff_id_eta_%d", ieta)));
        eff_converts_eta.push_back((TEfficiency *)fin->Get(Form("eff_converts_eta_%d", ieta)));
        eff_iso_eta.push_back((TEfficiency *)fin->Get(Form("eff_iso_eta_%d", ieta)));

    }

    std::vector<int> colors = {kPink+8, kSpring-7, kAzure-3, kViolet+3, kOrange +10, kYellow + 4};

    TCanvas *ceff_reco = new TCanvas("ceff_reco", "", 800, 750);
    ceff_reco->SetTopMargin(0.1);
    TH2D* hrecoframe = new TH2D("hrecoframe", "", 10, 10, 36, 1, 0.8, 1.02);
    hrecoframe->GetXaxis()->SetTitle("Truth #gamma p_{T} [GeV]");
    hrecoframe->GetYaxis()->SetTitle("Reco Efficiency");
    hrecoframe->Draw();
    for(int ieta = 0; ieta < (int)eta_bins.size() - 4; ieta++){
        eff_reco_eta[ieta]->SetMarkerStyle(20);
        eff_reco_eta[ieta]->SetMarkerSize(1);
        eff_reco_eta[ieta]->SetMarkerColor(colors[ieta]);
        eff_reco_eta[ieta]->SetLineColor(colors[ieta]);
        eff_reco_eta[ieta]->Draw("same");
        myMarkerText(0.24, 0.5 - 0.05 * ieta, colors[ieta], 20, Form("%.1f < #eta_{#gamma} < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),1,0.04);
    }
    myText(0.2, 0.95, kBlack, "#bf{#it{sPHENIX}} Simulation", 0.04);
    myText(0.6, 0.95, kBlack, "Pythia 200 GeV", 0.04);
    myText(0.2, 0.65, kBlack, Form("|vtx_z| < %.1f cm", vertexcut), 0.04);
    
    //myText(0.2, 0.60, kBlack, "Iso dR < 0.4", 0.04);
    //myText(0.2, 0.55, kBlack, Form("%.1f GeV < Iso E_{T}^{reco} < %.1f GeV", recoiso_min, recoiso_max), 0.04);
    
    TCanvas *ceff_id = new TCanvas("ceff_id", "", 800, 750);
    ceff_id->SetTopMargin(0.1);
    TH2D* hidframe = new TH2D("hidframe", "", 10, 10, 36, 1, 0.6, 1.02);
    hidframe->GetXaxis()->SetTitle("Truth #gamma p_{T} [GeV]");
    hidframe->GetYaxis()->SetTitle("ID Efficiency");
    hidframe->Draw();
    for(int ieta = 0; ieta < (int)eta_bins.size() - 4; ieta++){
        eff_id_eta[ieta]->SetMarkerStyle(20);
        eff_id_eta[ieta]->SetMarkerSize(1);
        eff_id_eta[ieta]->SetMarkerColor(colors[ieta]);
        eff_id_eta[ieta]->SetLineColor(colors[ieta]);
        eff_id_eta[ieta]->Draw("same");
        myMarkerText(0.24, 0.5 - 0.05 * ieta, colors[ieta], 20, Form("%.1f < #eta_{#gamma} < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),1,0.04);
    }
    myText(0.2, 0.95, kBlack, "#bf{#it{sPHENIX}} Simulation", 0.04);
    myText(0.6, 0.95, kBlack, "Pythia 200 GeV", 0.04);
    myText(0.2, 0.55, kBlack, Form("|vtx_z| < %.1f cm", vertexcut), 0.04);

    TCanvas *ceffiso = new TCanvas("ceffiso", "", 800, 750);
    ceffiso->SetTopMargin(0.1);
    TH2D* hisoframe = new TH2D("hisoframe", "", 10, 10, 36, 1, 0.5, 0.95);
    hisoframe->GetXaxis()->SetTitle("Truth #gamma p_{T} [GeV]");
    hisoframe->GetYaxis()->SetTitle("Iso Efficiency");
    hisoframe->Draw();
    for(int ieta = 0; ieta < (int)eta_bins.size() - 4; ieta++){
        eff_iso_eta[ieta]->SetMarkerStyle(20);
        eff_iso_eta[ieta]->SetMarkerSize(1);
        eff_iso_eta[ieta]->SetMarkerColor(colors[ieta]);
        eff_iso_eta[ieta]->SetLineColor(colors[ieta]);
        eff_iso_eta[ieta]->Draw("same");
        myMarkerText(0.24, 0.5 - 0.05 * ieta, colors[ieta], 20, Form("%.1f < #eta_{#gamma} < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),1,0.04);
    }
    myText(0.2, 0.65, kBlack, Form("|vtx_z| < %.1f cm", vertexcut), 0.04);
    myText(0.2, 0.55, kBlack, Form("%.1f GeV < Iso E_{T}^{reco} < %.1f GeV", recoiso_min, recoiso_max), 0.04);
    myText(0.2, 0.60, kBlack, Form("Iso dR < 0.%.1d", conesize), 0.04);
}
