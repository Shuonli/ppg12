#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <yaml-cpp/yaml.h>

// Keep the same time sample as in RecoEffCalculator_TTreeReader.C
const float TIME_SAMPLE_NS = 16.67;

static inline float wrapDeltaPhi(float dphi)
{
    while (dphi > M_PI) dphi -= 2 * M_PI;
    while (dphi < -M_PI) dphi += 2 * M_PI;
    return dphi;
}

void time_energy_corr(const std::string &configname = "config_bdt_none.yaml",
                      const std::string &filetype = "photon20")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    bool issim = true;
    if (filetype == "data")
    {
        issim = false;
    }

    // Build input path
    std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
    std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;
    //overwrite the input file name here:
    infilename = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/jet50/condor_0/timetest.root";
    if (!issim)
    {
        infilename = configYaml["input"]["data_file"].as<std::string>();
    }

    std::string treename = configYaml["input"]["tree"].as<std::string>();
    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();
    float eff_dR = configYaml["analysis"]["eff_dR"].as<float>();
    float reco_min_ET = configYaml["analysis"]["reco_min_ET"].as<float>();

    // Set up chain and reader
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());
    TTreeReader reader(&chain);

    // Event-level
    TTreeReaderValue<int> nparticles(reader, "nparticles");

    // Truth particle arrays
    TTreeReaderArray<int> particle_trkid(reader, "particle_trkid");
    TTreeReaderArray<int> particle_pid(reader, "particle_pid");
    TTreeReaderArray<float> particle_Pt(reader, "particle_Pt");
    TTreeReaderArray<float> particle_Eta(reader, "particle_Eta");
    TTreeReaderArray<float> particle_Phi(reader, "particle_Phi");

    // Cluster arrays (dynamic names)
    TTreeReaderValue<int> ncluster(reader, ("ncluster_" + clusternodename).c_str());
    TTreeReaderArray<float> cluster_Et(reader, ("cluster_Et_" + clusternodename).c_str());
    TTreeReaderArray<float> cluster_Eta(reader, ("cluster_Eta_" + clusternodename).c_str());
    TTreeReaderArray<float> cluster_Phi(reader, ("cluster_Phi_" + clusternodename).c_str());
    TTreeReaderArray<int> cluster_truthtrkID(reader, ("cluster_truthtrkID_" + clusternodename).c_str());

    // 2D arrays flattened to 1D (49 cells per cluster)
    TTreeReaderArray<float> cluster_e_array(reader, ("cluster_e_array_" + clusternodename).c_str());
    TTreeReaderArray<float> cluster_time_array(reader, ("cluster_time_array_" + clusternodename).c_str());
    TTreeReaderArray<int> cluster_ownership_array(reader, ("cluster_ownership_array_" + clusternodename).c_str());


    TFile *fout = new TFile("time_energy_corr_0.root", "RECREATE");
    TH3D *h_time_reco_truth = new TH3D("h_time_reco_truth", "h_time_reco_truth", 100, -3, 7, 100, 0, 100, 100, 0, 100);
    h_time_reco_truth->GetXaxis()->SetTitle("Time (sample)");
    h_time_reco_truth->GetYaxis()->SetTitle("Reco Cluster Et(GeV)");
    h_time_reco_truth->GetZaxis()->SetTitle("Truth Particle Pt (GeV)");
    TH3D *h_time_ratio_truth = new TH3D("h_time_ratio_truth", "h_time_ratio_truth", 100, -3, 7, 1000, 0, 2, 100, 0, 100);
    h_time_ratio_truth->GetXaxis()->SetTitle("Time (sample)");
    h_time_ratio_truth->GetYaxis()->SetTitle("Reco/Truth Et ratio");
    h_time_ratio_truth->GetZaxis()->SetTitle("Truth Particle Pt (GeV)");
    // Event loop
    int event_count = 0;
    while (reader.Next())
    {
        if(event_count%10000==0)
        {
            std::cout << "Processing event " << event_count << std::endl;
        }
        event_count++;
        // Build track id -> particle index map per event (for truth matching)
        std::map<int, int> particle_trkidmap;
        if (issim)
        {
            for (int ip = 0; ip < *nparticles; ++ip)
            {
                particle_trkidmap[particle_trkid[ip]] = ip;
            }
        }

        // Loop clusters
        for (int icluster = 0; icluster < *ncluster; ++icluster)
        {
            if (cluster_Et[icluster] < 10)
            {
                continue;
            }
            if (cluster_Eta[icluster] > 0.7 || cluster_Eta[icluster] < -0.7)
            {
                continue;
            }

            // Compute energy-weighted average time in ns
            float sum_e_t = 0;
            float sum_e = 0;
            for (int i = 0; i < 49; ++i)
            {
                if (cluster_ownership_array[icluster * 49 + i] == 1)
                {
                    const float e = cluster_e_array[icluster * 49 + i];
                    sum_e += e;
                    sum_e_t += e * cluster_time_array[icluster * 49 + i];
                }
            }
            float avg_time_samples = (sum_e > 0) ? (sum_e_t / sum_e) : 0.0f;
            float avg_time_ns = avg_time_samples * TIME_SAMPLE_NS;



            // Truth matching (if simulation)
            if (issim)
            {
                auto it = particle_trkidmap.find(cluster_truthtrkID[icluster]);
                if (it != particle_trkidmap.end())
                {
                    int ip = it->second;
                    float deta = cluster_Eta[icluster] - particle_Eta[ip];
                    float dphi = cluster_Phi[icluster] - particle_Phi[ip];
                    
                    if (particle_Pt[ip] < 10)
                    {
                        continue;
                    }
                    //require photon type 1 and 2
                    if (particle_pid[ip]!=22 && particle_pid[ip]!=111)
                    {
                        continue;
                    }
                    if (dphi > M_PI)
                    {
                        dphi = 2 * M_PI - dphi;
                    }
                    float dR = std::sqrt(deta * deta + dphi * dphi);
                    if (dR < eff_dR)
                    {
                        h_time_reco_truth->Fill(avg_time_samples, cluster_Et[icluster], particle_Pt[ip]);
                        h_time_ratio_truth->Fill(avg_time_samples, cluster_Et[icluster] / particle_Pt[ip], particle_Pt[ip]);
                    }


                }
            }

        }
    }

    fout->Write();
    fout->Close();
}


