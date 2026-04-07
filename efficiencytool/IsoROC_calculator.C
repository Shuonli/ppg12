// IsoROC_calculator.C
// Computes isoET distributions for truth-classified signal (direct+frag photons)
// and background clusters in ET bins for 4 isolation definitions:
//   iso70:   cluster_iso_03_70_emcal  + cluster_iso_03_70_hcalin  + cluster_iso_03_70_hcalout
//   iso60:   cluster_iso_03_60_emcal  + cluster_iso_03_60_hcalin  + cluster_iso_03_60_hcalout
//   iso120:  cluster_iso_03_120_emcal + cluster_iso_03_120_hcalin + cluster_iso_03_120_hcalout
//   isotopo: cluster_iso_topo_03_<clusternodename>
// Run on photon5/10/20 and jet5/12/20/30/40 samples, then hadd and plot ROC curves.
// No photon ID (shower-shape) cuts applied. Raw isoET (no MC scale/shift).

#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <vector>
#include <sstream>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TSystem.h>
#include <TObjString.h>
#include <yaml-cpp/yaml.h>
#include "CrossSectionWeights.h"
using namespace PPG12;

void IsoROC_calculator(const std::string &configname = "config_bdt_isoroc.yaml",
                       const std::string filetype = "jet5")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    bool issim = true;
    if (filetype == "data") issim = false;

    static const bool isbackground = filetype.find("jet") != std::string::npos;

    float max_photon_lower = 0;
    float max_photon_upper = 200;
    float max_jet_lower    = 0;
    float max_jet_upper    = 100;
    float cluster_ET_upper = 100;
    float cross_weight     = 1.0;

    if      (filetype == "photon5")  { max_photon_lower = 0;  max_photon_upper = 14;  cross_weight = photon5cross  / photon20cross; }
    else if (filetype == "photon10") { max_photon_lower = 14; max_photon_upper = 30;  cross_weight = photon10cross / photon20cross; }
    else if (filetype == "photon20") { max_photon_lower = 30; max_photon_upper = 200; cross_weight = 1.0; }
    else if (filetype == "jet5")     { max_jet_lower = 7;  max_jet_upper = 14;  cluster_ET_upper = 14;  cross_weight = jet5cross  / jet50cross; }
    else if (filetype == "jet12")    { max_jet_lower = 14; max_jet_upper = 21;  cluster_ET_upper = 23;  cross_weight = jet12cross / jet50cross; }
    else if (filetype == "jet20")    { max_jet_lower = 21; max_jet_upper = 32;  cluster_ET_upper = 34;  cross_weight = jet20cross / jet50cross; }
    else if (filetype == "jet30")    { max_jet_lower = 32; max_jet_upper = 42;  cluster_ET_upper = 44;  cross_weight = jet30cross / jet50cross; }
    else if (filetype == "jet40")    { max_jet_lower = 42; max_jet_upper = 100; cluster_ET_upper = 100; cross_weight = jet40cross / jet50cross; }
    else if (filetype == "jet50")    { max_jet_lower = 52; max_jet_upper = 100; cluster_ET_upper = 100; cross_weight = 1.0; }

    // Vertex reweighting
    TH1 *h_vertex_reweight = nullptr;
    int vertex_reweight_on = 1;
    if (issim)
    {
        vertex_reweight_on = configYaml["analysis"]["vertex_reweight_on"].as<int>(1);
        std::string vertex_reweight_file =
            configYaml["analysis"]["vertex_reweight_file"].as<std::string>("results/vertex_reweight.root");
        if (vertex_reweight_on)
        {
            TFile *fvtx = TFile::Open(vertex_reweight_file.c_str(), "READ");
            if (!fvtx || fvtx->IsZombie())
            {
                std::cerr << "[VertexReweight] ERROR: cannot open " << vertex_reweight_file << std::endl;
                return;
            }
            std::string vtx_histname = "h_vertexz_ratio_data_over_mccombined";
            h_vertex_reweight = dynamic_cast<TH1 *>(fvtx->Get(vtx_histname.c_str()));
            if (!h_vertex_reweight)
            {
                std::cerr << "[VertexReweight] ERROR: histogram not found" << std::endl;
                return;
            }
            h_vertex_reweight->SetDirectory(nullptr);
            fvtx->Close();
        }
    }

    // Config parameters
    std::string infilename_root_dir  = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
    std::string infilename           = infilename_root_dir + filetype + infilename_branch_dir;
    std::string treename             = configYaml["input"]["tree"].as<std::string>();
    std::string clusternodename      = configYaml["input"]["cluster_node_name"].as<std::string>();
    std::string bdt_model_name       = configYaml["input"]["bdt_model_name"].as<std::string>("base");
    std::string var_type             = configYaml["output"]["var_type"].as<std::string>();
    std::string outdir               = configYaml["output"]["eff_outfile"].as<std::string>();
    // Strip trailing filename from eff_outfile to get directory prefix
    // Output: results/roc_isoET_<filetype>_<var_type>.root
    {
        size_t slash = outdir.rfind('/');
        if (slash != std::string::npos) outdir = outdir.substr(0, slash);
    }
    std::string outfilename = outdir + "/roc_isoET_" + filetype + "_" + var_type + ".root";

    float vertexcut    = configYaml["analysis"]["vertex_cut"].as<float>();
    float truthisocut  = configYaml["analysis"]["truth_iso_max"].as<float>();
    float reco_min_ET  = configYaml["analysis"]["reco_min_ET"].as<float>();
    float eff_dR       = configYaml["analysis"]["eff_dR"].as<float>();

    std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();
    std::vector<float> pT_bins  = configYaml["analysis"]["pT_bins"].as<std::vector<float>>();
    int n_pT_bins = (int)pT_bins.size() - 1;
    double pT_bin_edges[pT_bins.size()];
    for (int i = 0; i < (int)pT_bins.size(); i++) pT_bin_edges[i] = pT_bins[i];

    std::cout << "infilename:  " << infilename  << std::endl;
    std::cout << "outfilename: " << outfilename << std::endl;

    // Build chain
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());

    // TTreeReader
    TTreeReader reader(&chain);

    // Event-level branches
    TTreeReaderValue<int>   nparticles(reader, "nparticles");
    TTreeReaderValue<int>   ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
    TTreeReaderValue<float> vertexz(reader, "vertexz");
    TTreeReaderValue<float> vertexz_truth(reader, "vertexz_truth");

    // Truth jet (needed for max_jet cut)
    TTreeReaderValue<int>         njet_truth(reader, "njet_truth");
    TTreeReaderArray<float>       jet_truth_Pt(reader, "jet_truth_Pt");

    // Particle arrays
    TTreeReaderArray<float> particle_Pt(reader, "particle_Pt");
    TTreeReaderArray<float> particle_Eta(reader, "particle_Eta");
    TTreeReaderArray<float> particle_Phi(reader, "particle_Phi");
    TTreeReaderArray<float> particle_truth_iso_03(reader, "particle_truth_iso_03");
    TTreeReaderArray<int>   particle_pid(reader, "particle_pid");
    TTreeReaderArray<int>   particle_trkid(reader, "particle_trkid");
    TTreeReaderArray<int>   particle_photonclass(reader, "particle_photonclass");

    // Cluster kinematics
    TTreeReaderArray<float> cluster_Et(reader,  Form("cluster_Et_%s",  clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusternodename.c_str()));
    TTreeReaderArray<int>   cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusternodename.c_str()));

    // Isolation branches (4 definitions)
    TTreeReaderArray<float> cluster_iso_03_70_emcal(reader,  Form("cluster_iso_03_70_emcal_%s",  clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalin(reader, Form("cluster_iso_03_70_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalout(reader,Form("cluster_iso_03_70_hcalout_%s",clusternodename.c_str()));

    TTreeReaderArray<float> cluster_iso_03_60_emcal(reader,  Form("cluster_iso_03_60_emcal_%s",  clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_hcalin(reader, Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_hcalout(reader,Form("cluster_iso_03_60_hcalout_%s",clusternodename.c_str()));

    TTreeReaderArray<float> cluster_iso_03_120_emcal(reader,  Form("cluster_iso_03_120_emcal_%s",  clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_120_hcalin(reader, Form("cluster_iso_03_120_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_120_hcalout(reader,Form("cluster_iso_03_120_hcalout_%s",clusternodename.c_str()));

    TTreeReaderArray<float> cluster_iso_topo_03(reader, Form("cluster_iso_topo_03_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_topo_04(reader, Form("cluster_iso_topo_04_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_topo_soft_03(reader, Form("cluster_iso_topo_soft_03_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_topo_soft_04(reader, Form("cluster_iso_topo_soft_04_%s", clusternodename.c_str()));

    // Output file and histograms
    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    TH1::SetDefaultSumw2(kTRUE);

    // iso definition labels
    const char *iso_labels[7]   = {"iso70", "iso60", "iso120", "isotopo03", "isotopo04", "isotoposoft03", "isotoposoft04"};
    const char *iso_titles[7]   = {"70 MeV threshold", "60 MeV threshold", "120 MeV threshold", "topo R=0.3", "topo R=0.4", "topo soft R=0.3", "topo soft R=0.4"};

    // isoET histogram range: -5 to 20 GeV, 500 bins (0.05 GeV resolution)
    const int iso_nbins = 500;
    const float iso_lo  = -5.0;
    const float iso_hi  = 20.0;

    // Histograms per [iso_def][pT_bin]
    // Categories: direct signal, frag signal, signal (direct+frag combined), background
    std::vector<std::vector<TH1D*>> h_direct(7);
    std::vector<std::vector<TH1D*>> h_frag(7);
    std::vector<std::vector<TH1D*>> h_signal(7);
    std::vector<std::vector<TH1D*>> h_bg(7);

    for (int iiso = 0; iiso < 7; iiso++)
    {
        for (int ipt = 0; ipt < n_pT_bins; ipt++)
        {
            float pt_lo = pT_bins[ipt];
            float pt_hi = pT_bins[ipt + 1];

            h_direct[iiso].push_back(new TH1D(
                Form("h_direct_%s_pt%d", iso_labels[iiso], ipt),
                Form("Direct photon isoET (%s) %.0f<ET<%.0f GeV; isoET [GeV]; Weighted counts",
                     iso_titles[iiso], pt_lo, pt_hi),
                iso_nbins, iso_lo, iso_hi));

            h_frag[iiso].push_back(new TH1D(
                Form("h_frag_%s_pt%d", iso_labels[iiso], ipt),
                Form("Frag photon isoET (%s) %.0f<ET<%.0f GeV; isoET [GeV]; Weighted counts",
                     iso_titles[iiso], pt_lo, pt_hi),
                iso_nbins, iso_lo, iso_hi));

            h_signal[iiso].push_back(new TH1D(
                Form("h_signal_%s_pt%d", iso_labels[iiso], ipt),
                Form("Signal isoET (%s) %.0f<ET<%.0f GeV; isoET [GeV]; Weighted counts",
                     iso_titles[iiso], pt_lo, pt_hi),
                iso_nbins, iso_lo, iso_hi));

            h_bg[iiso].push_back(new TH1D(
                Form("h_bg_%s_pt%d", iso_labels[iiso], ipt),
                Form("Background isoET (%s) %.0f<ET<%.0f GeV; isoET [GeV]; Weighted counts",
                     iso_titles[iiso], pt_lo, pt_hi),
                iso_nbins, iso_lo, iso_hi));
        }
    }

    // Also store pT bin edges in output for reference
    TH1D *h_pT_bins = new TH1D("h_pT_bins", "pT bin edges", n_pT_bins, pT_bin_edges);

    int nentries = chain.GetEntries();
    int ientry   = 0;

    while (reader.Next())
    {
        if (ientry % 50000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;

        // Per-event weight: start from cross-section weight
        float weight       = cross_weight;
        float vertex_weight = 1.0;

        // Vertex reweighting
        if (issim && vertex_reweight_on && h_vertex_reweight)
        {
            int bin = h_vertex_reweight->FindBin(*vertexz);
            if (bin < 1) bin = 1;
            if (bin > h_vertex_reweight->GetNbinsX()) bin = h_vertex_reweight->GetNbinsX();
            vertex_weight = h_vertex_reweight->GetBinContent(bin);
            if (!std::isfinite(vertex_weight) || vertex_weight <= 0.0) vertex_weight = 1.0;
            weight *= vertex_weight;
        }

        // --- Particle loop: build truth classification maps ---
        // photon_map: iparticle → photonclass (1=direct, 2=frag) for signal photons
        std::map<int, int> particle_trkidmap;       // trkid → iparticle index
        std::map<int, int> photon_map;              // iparticle → photonclass (1 or 2)

        if (issim)
        {
            float maxphotonpT = 0;
            for (int ipart = 0; ipart < *nparticles; ipart++)
            {
                particle_trkidmap[particle_trkid[ipart]] = ipart;

                if (particle_pid[ipart] == 22)
                {
                    if (particle_Pt[ipart] > maxphotonpT) maxphotonpT = particle_Pt[ipart];

                    if (particle_photonclass[ipart] == 1 || particle_photonclass[ipart] == 2)
                    {
                        float truthiso = particle_truth_iso_03[ipart];
                        if (truthiso < truthisocut)
                        {
                            photon_map[ipart] = particle_photonclass[ipart];
                        }
                    }
                }
            }

            // Apply photon sample max pT filter
            if (!isbackground)
            {
                if (maxphotonpT > max_photon_upper || maxphotonpT < max_photon_lower)
                {
                    ientry++;
                    continue;
                }
            }
        }

        // Apply jet sample max jet pT filter
        if (isbackground && issim)
        {
            float maxjetpT = 0;
            for (int ijet = 0; ijet < *njet_truth; ijet++)
            {
                if (jet_truth_Pt[ijet] > maxjetpT) maxjetpT = jet_truth_Pt[ijet];
            }
            if (maxjetpT == 0 || maxjetpT > max_jet_upper || maxjetpT < max_jet_lower)
            {
                ientry++;
                continue;
            }
        }

        // Vertex cut
        if (std::abs(*vertexz) > vertexcut)
        {
            ientry++;
            continue;
        }

        // --- Cluster loop ---
        for (int icluster = 0; icluster < *ncluster; icluster++)
        {
            float clusterET  = cluster_Et[icluster];
            float clusterEta = cluster_Eta[icluster];

            // Basic ET cut
            if (clusterET < reco_min_ET) continue;

            // ET upper cut for jet (bg) samples
            if (isbackground && clusterET > cluster_ET_upper) continue;

            // Eta acceptance
            int etabin = -1;
            for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
            {
                if (clusterEta > eta_bins[ieta] && clusterEta < eta_bins[ieta + 1])
                {
                    etabin = ieta;
                    break;
                }
            }
            if (etabin == -1) continue;

            // ET bin
            int pTbin = -1;
            for (int ipt = 0; ipt < n_pT_bins; ipt++)
            {
                if (clusterET > pT_bins[ipt] && clusterET < pT_bins[ipt + 1])
                {
                    pTbin = ipt;
                    break;
                }
            }
            if (pTbin == -1) continue;

            // Compute the 4 raw isoET values (no MC scale/shift applied)
            float iso70      = cluster_iso_03_70_emcal[icluster]  + cluster_iso_03_70_hcalin[icluster]  + cluster_iso_03_70_hcalout[icluster];
            float iso60      = cluster_iso_03_60_emcal[icluster]  + cluster_iso_03_60_hcalin[icluster]  + cluster_iso_03_60_hcalout[icluster];
            float iso120     = cluster_iso_03_120_emcal[icluster] + cluster_iso_03_120_hcalin[icluster] + cluster_iso_03_120_hcalout[icluster];
            float isotopo03      = cluster_iso_topo_03[icluster];
            float isotopo04      = cluster_iso_topo_04[icluster];
            float isotoposoft03  = cluster_iso_topo_soft_03[icluster];
            float isotoposoft04  = cluster_iso_topo_soft_04[icluster];
            float iso_vals[7] = {iso70, iso60, iso120, isotopo03, isotopo04, isotoposoft03, isotoposoft04};

            // Truth classification (simulation only)
            int cluster_class = 0; // 0=background, 1=direct signal, 2=frag signal
            if (issim)
            {
                auto it = particle_trkidmap.find(cluster_truthtrkID[icluster]);
                if (it != particle_trkidmap.end())
                {
                    int ipart = it->second;

                    // Check dR between cluster and matched truth particle
                    float deta = cluster_Eta[icluster] - particle_Eta[ipart];
                    float dphi = cluster_Phi[icluster] - particle_Phi[ipart];
                    while (dphi >  M_PI) dphi -= 2 * M_PI;
                    while (dphi < -M_PI) dphi += 2 * M_PI;
                    float dR = std::sqrt(deta * deta + dphi * dphi);

                    if (dR < eff_dR)
                    {
                        auto jt = photon_map.find(ipart);
                        if (jt != photon_map.end())
                        {
                            cluster_class = jt->second; // 1 or 2
                        }
                    }
                }
            }

            // Fill histograms
            for (int iiso = 0; iiso < 7; iiso++)
            {
                float iso_val = iso_vals[iiso];
                if (cluster_class == 1)
                {
                    h_direct[iiso][pTbin]->Fill(iso_val, weight);
                    h_signal[iiso][pTbin]->Fill(iso_val, weight);
                }
                else if (cluster_class == 2)
                {
                    h_frag[iiso][pTbin]->Fill(iso_val, weight);
                    h_signal[iiso][pTbin]->Fill(iso_val, weight);
                }
                else
                {
                    h_bg[iiso][pTbin]->Fill(iso_val, weight);
                }
            }
        } // end cluster loop

        ientry++;
    } // end event loop

    fout->cd();
    for (int iiso = 0; iiso < 7; iiso++)
    {
        for (int ipt = 0; ipt < n_pT_bins; ipt++)
        {
            h_direct[iiso][ipt]->Write();
            h_frag[iiso][ipt]->Write();
            h_signal[iiso][ipt]->Write();
            h_bg[iiso][ipt]->Write();
        }
    }
    h_pT_bins->Write();
    fout->Close();

    std::cout << "Output written to: " << outfilename << std::endl;
}
