#include <iostream>
#include <string>
#include <fstream>
#include <set>
#include <sstream>
#include <map>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TSystem.h>
#include <yaml-cpp/yaml.h>
#include <cmath>

const float TIME_SAMPLE_NS = 17.6;

void plot_cluster_jet_time(const std::string &configname = "config_bdt_none.yaml", const std::string filetype = "data")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    bool issim = (filetype != "data");

    std::string infilename;
    if (!issim)
    {
        infilename = configYaml["input"]["data_file"].as<std::string>();
    }
    else
    {
        std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
        std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
        infilename = infilename_root_dir + filetype + infilename_branch_dir;
    }

    std::cout << "Input file: " << infilename << std::endl;

    // Build input chain
    std::string treename = configYaml["input"]["tree"].as<std::string>();
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());

    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();

    // Read analysis parameters from config
    float vertexcut = configYaml["analysis"]["vertex_cut"].as<float>();
    float reco_min_ET = configYaml["analysis"]["reco_min_ET"].as<float>();

    // Get pT bins
    std::vector<float> pT_bins = configYaml["analysis"]["pT_bins"].as<std::vector<float>>();
    int n_pT_bins = pT_bins.size() - 1;
    double pT_bin_edges[n_pT_bins + 1];
    std::copy(pT_bins.begin(), pT_bins.end(), pT_bin_edges);

    // Common selection parameters
    float common_prob_max = configYaml["analysis"]["common"]["prob_max"].as<float>();
    float common_prob_min = configYaml["analysis"]["common"]["prob_min"].as<float>();
    float common_e11_over_e33_max = configYaml["analysis"]["common"]["e11_over_e33_max"].as<float>();
    float common_e11_over_e33_min = configYaml["analysis"]["common"]["e11_over_e33_min"].as<float>();
    float common_cluster_weta_cogx_bound = configYaml["analysis"]["common"]["cluster_weta_cogx_bound"].as<float>();
    int common_b2bjet_cut = configYaml["analysis"]["common_b2bjet_cut"].as<int>(0);
    float common_b2bjet_pt_min = configYaml["analysis"]["common_b2bjet_pt_min"].as<float>(7.0);

    // Tight selection parameters
    float tight_weta_cogx_min = configYaml["analysis"]["tight"]["weta_cogx_min"].as<float>();
    float tight_weta_cogx_max_b = configYaml["analysis"]["tight"]["weta_cogx_max_b"].as<float>();
    float tight_weta_cogx_max_s = configYaml["analysis"]["tight"]["weta_cogx_max_s"].as<float>();
    float tight_wphi_cogx_min = configYaml["analysis"]["tight"]["wphi_cogx_min"].as<float>();
    float tight_wphi_cogx_max_b = configYaml["analysis"]["tight"]["wphi_cogx_max_b"].as<float>();
    float tight_wphi_cogx_max_s = configYaml["analysis"]["tight"]["wphi_cogx_max_s"].as<float>();
    float tight_et1_min_b = configYaml["analysis"]["tight"]["et1_min_b"].as<float>();
    float tight_et1_min_s = configYaml["analysis"]["tight"]["et1_min_s"].as<float>();
    float tight_et1_max = configYaml["analysis"]["tight"]["et1_max"].as<float>();
    float tight_et2_min = configYaml["analysis"]["tight"]["et2_min"].as<float>(0.0);
    float tight_et2_max = configYaml["analysis"]["tight"]["et2_max"].as<float>(1.0);
    float tight_et3_min = configYaml["analysis"]["tight"]["et3_min"].as<float>(0.0);
    float tight_et3_max = configYaml["analysis"]["tight"]["et3_max"].as<float>(1.0);
    float tight_et4_min = configYaml["analysis"]["tight"]["et4_min"].as<float>();
    float tight_et4_max = configYaml["analysis"]["tight"]["et4_max"].as<float>();
    float tight_e11_over_e33_min = configYaml["analysis"]["tight"]["e11_over_e33_min"].as<float>();
    float tight_e11_over_e33_max = configYaml["analysis"]["tight"]["e11_over_e33_max"].as<float>();
    float tight_e32_over_e35_min = configYaml["analysis"]["tight"]["e32_over_e35_min"].as<float>();
    float tight_e32_over_e35_max = configYaml["analysis"]["tight"]["e32_over_e35_max"].as<float>();
    float tight_prob_min = configYaml["analysis"]["tight"]["prob_min"].as<float>();
    float tight_prob_max = configYaml["analysis"]["tight"]["prob_max"].as<float>();
    float tight_bdt_min = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0);
    float tight_bdt_max = configYaml["analysis"]["tight"]["bdt_max"].as<float>(1);

    // Isolation parameters
    float recoiso_min = configYaml["analysis"]["reco_iso_min"].as<float>();
    float recoiso_max_b = configYaml["analysis"]["reco_iso_max_b"].as<float>();
    float recoiso_max_s = configYaml["analysis"]["reco_iso_max_s"].as<float>();
    float recononiso_min_shift = configYaml["analysis"]["reco_noniso_min_shift"].as<float>();
    float recononiso_max = configYaml["analysis"]["reco_noniso_max"].as<float>();

    int conesize = configYaml["analysis"]["cone_size"].as<int>();
    std::string bdt_model_name = configYaml["input"]["bdt_model_name"].as<std::string>("base");

    // TTreeReader setup
    TTreeReader reader(&chain);
    TTreeReaderValue<int> mbdnorthhit(reader, "mbdnorthhit");
    TTreeReaderValue<int> mbdsouthhit(reader, "mbdsouthhit");
    TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
    TTreeReaderValue<int> runnumber(reader, "runnumber");
    TTreeReaderValue<float> vertexz(reader, "vertexz");

    // Jet arrays including jet_time
    TTreeReaderValue<int> njet(reader, "njet");
    TTreeReaderArray<float> jet_Pt(reader, "jet_Pt");
    TTreeReaderArray<float> jet_Eta(reader, "jet_Eta");
    TTreeReaderArray<float> jet_Phi(reader, "jet_Phi");
    TTreeReaderArray<float> jet_time(reader, "jet_time");  // NEW BRANCH

    // Cluster arrays
    TTreeReaderArray<float> cluster_E(reader, Form("cluster_E_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_prob(reader, Form("cluster_prob_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e11(reader, Form("cluster_e11_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e33(reader, Form("cluster_e33_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e35(reader, Form("cluster_e35_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e32(reader, Form("cluster_e32_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_weta_cogx(reader, Form("cluster_weta_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_wphi_cogx(reader, Form("cluster_wphi_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et1(reader, Form("cluster_et1_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et2(reader, Form("cluster_et2_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et3(reader, Form("cluster_et3_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et4(reader, Form("cluster_et4_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_nsaturated(reader, Form("cluster_nsaturated_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03(reader, Form("cluster_iso_03_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_bdt(reader, Form("cluster_bdt_%s_%s", clusternodename.c_str(), bdt_model_name.c_str()));

    // Cluster time arrays
    TTreeReaderArray<float> cluster_e_array(reader, Form("cluster_e_array_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_time_array(reader, Form("cluster_time_array_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_status_array(reader, Form("cluster_status_array_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ownership_array(reader, Form("cluster_ownership_array_%s", clusternodename.c_str()));

    // Create output file
    std::string outfilename = "cluster_jet_time_" + filetype + ".root";
    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");

    // Create histograms for different selections and pT bins
    // For each selection: common, tight_iso, tight_noniso, NPB
    std::vector<TH2D*> h_common_delta_t_jet;
    std::vector<TH2D*> h_tight_iso_delta_t_jet;
    std::vector<TH2D*> h_tight_noniso_delta_t_jet;
    std::vector<TH2D*> h_npb_delta_t_jet;

    for (int ipt = 0; ipt < n_pT_bins; ipt++)
    {
        h_common_delta_t_jet.push_back(new TH2D(
            Form("h_common_delta_t_jet_pt%d", ipt),
            Form("Common Selection: Cluster-Jet Time vs Cluster pT (%.1f < pT < %.1f GeV);Cluster pT (GeV);Cluster-Jet Time (ns)",
                 pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
            100, pT_bin_edges[ipt], pT_bin_edges[ipt + 1], 160, -40, 40));

        h_tight_iso_delta_t_jet.push_back(new TH2D(
            Form("h_tight_iso_delta_t_jet_pt%d", ipt),
            Form("Tight+Iso Selection: Cluster-Jet Time vs Cluster pT (%.1f < pT < %.1f GeV);Cluster pT (GeV);Cluster-Jet Time (ns)",
                 pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
            100, pT_bin_edges[ipt], pT_bin_edges[ipt + 1], 160, -40, 40));

        h_tight_noniso_delta_t_jet.push_back(new TH2D(
            Form("h_tight_noniso_delta_t_jet_pt%d", ipt),
            Form("Tight+NonIso Selection: Cluster-Jet Time vs Cluster pT (%.1f < pT < %.1f GeV);Cluster pT (GeV);Cluster-Jet Time (ns)",
                 pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
            100, pT_bin_edges[ipt], pT_bin_edges[ipt + 1], 160, -40, 40));

        h_npb_delta_t_jet.push_back(new TH2D(
            Form("h_npb_delta_t_jet_pt%d", ipt),
            Form("NPB Selection: Cluster-Jet Time vs Cluster pT (%.1f < pT < %.1f GeV);Cluster pT (GeV);Cluster-Jet Time (ns)",
                 pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
            100, pT_bin_edges[ipt], pT_bin_edges[ipt + 1], 160, -40, 40));
    }

    // Also create inclusive histograms (all pT bins combined)
    TH2D *h_common_delta_t_jet_all = new TH2D("h_common_delta_t_jet_all",
                                               "Common Selection: Cluster-Jet Time vs Cluster pT;Cluster pT (GeV);Cluster-Jet Time (ns)",
                                               n_pT_bins, pT_bin_edges, 160, -40, 40);
    TH2D *h_tight_iso_delta_t_jet_all = new TH2D("h_tight_iso_delta_t_jet_all",
                                                  "Tight+Iso Selection: Cluster-Jet Time vs Cluster pT;Cluster pT (GeV);Cluster-Jet Time (ns)",
                                                  n_pT_bins, pT_bin_edges, 160, -40, 40);
    TH2D *h_tight_noniso_delta_t_jet_all = new TH2D("h_tight_noniso_delta_t_jet_all",
                                                     "Tight+NonIso Selection: Cluster-Jet Time vs Cluster pT;Cluster pT (GeV);Cluster-Jet Time (ns)",
                                                     n_pT_bins, pT_bin_edges, 160, -40, 40);
    TH2D *h_npb_delta_t_jet_all = new TH2D("h_npb_delta_t_jet_all",
                                            "NPB Selection: Cluster-Jet Time vs Cluster pT;Cluster pT (GeV);Cluster-Jet Time (ns)",
                                            n_pT_bins, pT_bin_edges, 160, -40, 40);

    // Skip bad runs
    std::set<int> skiprunnumbers = {47698, 51489, 51721, 51725, 53284};

    int nentries = chain.GetEntries();
    int ientry = 0;

    std::cout << "Processing " << nentries << " entries..." << std::endl;

    while (reader.Next())
    {
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;

        // Skip bad runs for data
        if (!issim && skiprunnumbers.find(*runnumber) != skiprunnumbers.end())
        {
            ientry++;
            continue;
        }

        // Vertex cut
        if (std::abs(*vertexz) > vertexcut)
        {
            ientry++;
            continue;
        }

        // MBD hit requirement
        if (!(*mbdnorthhit >= 1 && *mbdsouthhit >= 1))
        {
            ientry++;
            continue;
        }

        // Find the closest jet (or leading jet) for timing reference
        // For now, let's use the leading jet
        float leading_jet_time = 0;
        float leading_jet_pt = -1;
        for (int ijet = 0; ijet < *njet; ijet++)
        {
            if (jet_Pt[ijet] > leading_jet_pt)
            {
                leading_jet_pt = jet_Pt[ijet];
                leading_jet_time = jet_time[ijet];
            }
        }

        // If no jets, skip this event
        if (leading_jet_pt < 0)
        {
            ientry++;
            continue;
        }

        // Loop over clusters
        for (int icluster = 0; icluster < *ncluster; icluster++)
        {
            // Minimum ET cut
            if (cluster_Et[icluster] < reco_min_ET)
                continue;

            // Calculate cluster time
            float clusteravgtime = 0;
            float cluster_total_e = 0;
            for (int i = 0; i < 49; i++)
            {
                if (cluster_ownership_array[icluster * 49 + i] == 1)
                {
                    int status = cluster_status_array[icluster * 49 + i];
                    if (status & (1 << 5))
                    {
                        continue;
                    }
                    clusteravgtime += cluster_time_array[icluster * 49 + i] * cluster_e_array[icluster * 49 + i];
                    cluster_total_e += cluster_e_array[icluster * 49 + i];
                }
            }
            clusteravgtime = cluster_total_e > 0 ? clusteravgtime / cluster_total_e : 0;
            clusteravgtime = clusteravgtime * TIME_SAMPLE_NS;

            // Calculate delta_t between cluster and leading jet
            float delta_t_jet = clusteravgtime - leading_jet_time;

            float cluster_eta = cluster_Eta[icluster];
            float clusterET = cluster_Et[icluster];

            // Calculate e11_over_e33 and e32_over_e35
            float e11_over_e33 = cluster_e11[icluster] / cluster_e33[icluster];
            float e32_over_e35 = cluster_e32[icluster] / cluster_e35[icluster];

            // Find pT bin
            int pTbin = -1;
            for (int ipt = 0; ipt < n_pT_bins; ipt++)
            {
                if (clusterET > pT_bins[ipt] && clusterET < pT_bins[ipt + 1])
                {
                    pTbin = ipt;
                    break;
                }
            }

            // Calculate isolation
            float recoisoET = cluster_iso_03[icluster];

            // Check for b2b jet and otherside_jet
            bool otherside_jet = false;
            float max_b2bjet_pT = -1;
            float b2bjet_dphi = 3 * M_PI / 4;
            float jet_eta_cut = 0.6;

            for (int ijet = 0; ijet < *njet; ijet++)
            {
                float dphi = cluster_Phi[icluster] - jet_Phi[ijet];
                while (dphi > M_PI)
                    dphi = dphi - 2 * M_PI;
                while (dphi < -M_PI)
                    dphi = dphi + 2 * M_PI;

                if (std::abs(dphi) > (M_PI / 2))
                {
                    otherside_jet = true;
                }

                if (std::abs(jet_Eta[ijet]) < jet_eta_cut)
                {
                    if (std::abs(dphi) > b2bjet_dphi)
                    {
                        if (jet_Pt[ijet] > max_b2bjet_pT && jet_Pt[ijet] > 5)
                        {
                            max_b2bjet_pT = jet_Pt[ijet];
                        }
                    }
                }
            }

            bool passes_common_b2bjet = true;
            if (common_b2bjet_cut)
            {
                passes_common_b2bjet = (max_b2bjet_pT >= common_b2bjet_pt_min);
            }

            // Common selection
            bool passes_common_shape =
                cluster_prob[icluster] > common_prob_min &&
                cluster_prob[icluster] < common_prob_max &&
                e11_over_e33 > common_e11_over_e33_min &&
                e11_over_e33 < common_e11_over_e33_max &&
                (cluster_weta_cogx[icluster] < common_cluster_weta_cogx_bound);

            bool common_pass = passes_common_shape && passes_common_b2bjet;

            // Tight selection
            float tight_weta_cogx_max = tight_weta_cogx_max_b + tight_weta_cogx_max_s * clusterET;
            float tight_wphi_cogx_max = tight_wphi_cogx_max_b + tight_wphi_cogx_max_s * clusterET;
            float tight_et1_min = tight_et1_min_b + tight_et1_min_s * clusterET;

            bool is_cluster_weta_cogx_tight =
                (cluster_weta_cogx[icluster] > tight_weta_cogx_min) &&
                (cluster_weta_cogx[icluster] < tight_weta_cogx_max);

            bool is_cluster_wphi_cogx_tight =
                (cluster_wphi_cogx[icluster] > tight_wphi_cogx_min) &&
                (cluster_wphi_cogx[icluster] < tight_wphi_cogx_max);

            bool is_cluster_et1_tight =
                (cluster_et1[icluster] > tight_et1_min) &&
                (cluster_et1[icluster] < tight_et1_max);

            bool is_cluster_et2_tight =
                (cluster_et2[icluster] > tight_et2_min) &&
                (cluster_et2[icluster] < tight_et2_max);

            bool is_cluster_et3_tight =
                (cluster_et3[icluster] > tight_et3_min) &&
                (cluster_et3[icluster] < tight_et3_max);

            bool is_e11_over_e33_tight =
                (e11_over_e33 > tight_e11_over_e33_min) &&
                (e11_over_e33 < tight_e11_over_e33_max);

            bool is_e32_over_e35_tight =
                (e32_over_e35 > tight_e32_over_e35_min) &&
                (e32_over_e35 < tight_e32_over_e35_max);

            bool is_cluster_et4_tight =
                (cluster_et4[icluster] > tight_et4_min) &&
                (cluster_et4[icluster] < tight_et4_max);

            bool is_cluster_prob_tight =
                (cluster_prob[icluster] > tight_prob_min) &&
                (cluster_prob[icluster] < tight_prob_max);

            bool is_bdt_tight =
                (cluster_bdt[icluster] > tight_bdt_min) &&
                (cluster_bdt[icluster] < tight_bdt_max);

            bool tight = false;
            if (common_pass &&
                is_cluster_weta_cogx_tight &&
                is_cluster_wphi_cogx_tight &&
                is_cluster_et1_tight &&
                is_cluster_et2_tight &&
                is_cluster_et3_tight &&
                is_e11_over_e33_tight &&
                is_e32_over_e35_tight &&
                is_cluster_et4_tight &&
                is_cluster_prob_tight &&
                is_bdt_tight)
            {
                tight = true;
            }

            // Isolation cuts
            float recoiso_max = recoiso_max_b + recoiso_max_s * clusterET;
            float recononiso_min = recoiso_max + recononiso_min_shift;

            bool iso = false;
            bool noniso = false;

            if (recoisoET > recoiso_min && recoisoET < recoiso_max)
            {
                iso = true;
            }
            if (recoisoET > recononiso_min && recoisoET < recononiso_max)
            {
                noniso = true;
            }

            // NPB selection
            float bg_timing_cut = -0.5;
            float npb_weta_min = 0.4;
            bool badtime = clusteravgtime < bg_timing_cut;
            bool isnpb = badtime && !otherside_jet && (cluster_weta_cogx[icluster] >= npb_weta_min);

            // Fill histograms
            if (common_pass)
            {
                h_common_delta_t_jet_all->Fill(clusterET, delta_t_jet);
                if (pTbin >= 0)
                {
                    h_common_delta_t_jet[pTbin]->Fill(clusterET, delta_t_jet);
                }
            }

            if (tight && iso)
            {
                h_tight_iso_delta_t_jet_all->Fill(clusterET, delta_t_jet);
                if (pTbin >= 0)
                {
                    h_tight_iso_delta_t_jet[pTbin]->Fill(clusterET, delta_t_jet);
                }
            }

            if (tight && noniso)
            {
                h_tight_noniso_delta_t_jet_all->Fill(clusterET, delta_t_jet);
                if (pTbin >= 0)
                {
                    h_tight_noniso_delta_t_jet[pTbin]->Fill(clusterET, delta_t_jet);
                }
            }

            if (isnpb)
            {
                h_npb_delta_t_jet_all->Fill(clusterET, delta_t_jet);
                if (pTbin >= 0)
                {
                    h_npb_delta_t_jet[pTbin]->Fill(clusterET, delta_t_jet);
                }
            }
        }

        ientry++;
    }

    // Write histograms
    fout->cd();
    h_common_delta_t_jet_all->Write();
    h_tight_iso_delta_t_jet_all->Write();
    h_tight_noniso_delta_t_jet_all->Write();
    h_npb_delta_t_jet_all->Write();

    for (int ipt = 0; ipt < n_pT_bins; ipt++)
    {
        h_common_delta_t_jet[ipt]->Write();
        h_tight_iso_delta_t_jet[ipt]->Write();
        h_tight_noniso_delta_t_jet[ipt]->Write();
        h_npb_delta_t_jet[ipt]->Write();
    }

    fout->Close();
    delete fout;

    std::cout << "Output written to: " << outfilename << std::endl;
    std::cout << "Created histograms for " << n_pT_bins << " pT bins" << std::endl;
    std::cout << "Selections: common, tight+iso, tight+noniso, NPB" << std::endl;
}
