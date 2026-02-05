#include <iostream>
#include <string>
#include <fstream>
#include <set>
#include <sstream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TF1.h>
#include <TSystem.h>
#include <yaml-cpp/yaml.h>
#include <cmath>

const float TIME_SAMPLE_NS = 17.6;

void plot_cluster_mbd_time_eta(const std::string &configname = "config_bdt_none.yaml", const std::string filetype = "data")
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

    // Common selection parameters
    float common_prob_max = configYaml["analysis"]["common"]["prob_max"].as<float>();
    float common_prob_min = configYaml["analysis"]["common"]["prob_min"].as<float>();
    float common_e11_over_e33_max = configYaml["analysis"]["common"]["e11_over_e33_max"].as<float>();
    float common_e11_over_e33_min = configYaml["analysis"]["common"]["e11_over_e33_min"].as<float>();
    float common_cluster_weta_cogx_bound = configYaml["analysis"]["common"]["cluster_weta_cogx_bound"].as<float>();
    int common_b2bjet_cut = configYaml["analysis"]["common_b2bjet_cut"].as<int>(0);
    float common_b2bjet_pt_min = configYaml["analysis"]["common_b2bjet_pt_min"].as<float>(7.0);

    // TTreeReader setup
    TTreeReader reader(&chain);
    TTreeReaderValue<int> mbdnorthhit(reader, "mbdnorthhit");
    TTreeReaderValue<int> mbdsouthhit(reader, "mbdsouthhit");
    TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
    TTreeReaderValue<int> runnumber(reader, "runnumber");
    TTreeReaderValue<float> vertexz(reader, "vertexz");
    TTreeReaderValue<float> mbd_time(reader, "mbd_time");

    // Cluster arrays
    TTreeReaderArray<float> cluster_E(reader, Form("cluster_E_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_prob(reader, Form("cluster_prob_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e11(reader, Form("cluster_e11_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e33(reader, Form("cluster_e33_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_weta_cogx(reader, Form("cluster_weta_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_nsaturated(reader, Form("cluster_nsaturated_%s", clusternodename.c_str()));

    // Cluster time arrays
    TTreeReaderArray<float> cluster_e_array(reader, Form("cluster_e_array_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_time_array(reader, Form("cluster_time_array_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_status_array(reader, Form("cluster_status_array_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ownership_array(reader, Form("cluster_ownership_array_%s", clusternodename.c_str()));

    // Jet arrays for b2b cut
    TTreeReaderValue<int> njet(reader, "njet");
    TTreeReaderArray<float> jet_Pt(reader, "jet_Pt");
    TTreeReaderArray<float> jet_Eta(reader, "jet_Eta");
    TTreeReaderArray<float> jet_Phi(reader, "jet_Phi");

    // Load MBD t0 correction
    std::cout << "Loading MBD t0 correction" << std::endl;
    std::map<int, float> mbd_t0_correction;
    std::ifstream file("/sphenix/user/shuhangli/ppg12/efficiencytool/MbdOut.corr");
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        int runnumber_corr;
        float t0;
        iss >> runnumber_corr >> t0;
        mbd_t0_correction[runnumber_corr] = t0;
    }
    std::cout << "Loaded MBD t0 correction" << std::endl;

    // Create output file and histograms
    std::string outfilename = "cluster_mbd_time_eta_" + filetype + ".root";
    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");

    // Histogram 1: (cluster - MBD time) vs. cluster eta for clusters with common selection and pT>10
    TH2D *h_common_delta_t_vs_eta = new TH2D("h_common_delta_t_vs_eta",
                                             "Cluster-MBD Time vs Eta (Common Selection, pT>10 GeV);Cluster #eta;Cluster-MBD Time (ns)",
                                             50, -1.0, 1.0, 160, -40, 40);

    // Histogram 2: (cluster - MBD time) vs. cluster eta for clusters with weta_cogx>1.0 and pT>10
    TH2D *h_weta1p0_delta_t_vs_eta = new TH2D("h_weta1p0_delta_t_vs_eta",
                                              "Cluster-MBD Time vs Eta (weta_cogx>1.0, pT>10 GeV);Cluster #eta;Cluster-MBD Time (ns)",
                                              50, -1.0, 1.0, 160, -40, 40);

    // Histogram 3: (cluster - MBD time) vs. cluster eta for clusters with NPB selection and pT>10
    TH2D *h_npb_delta_t_vs_eta = new TH2D("h_npb_delta_t_vs_eta",
                                          "Cluster-MBD Time vs Eta (NPB Selection, pT>10 GeV);Cluster #eta;Cluster-MBD Time (ns)",
                                          50, -1.0, 1.0, 160, -40, 40);

    // Histogram 4: cluster time vs. cluster eta for clusters with common selection and pT>10
    TH2D *h_common_cluster_t_vs_eta = new TH2D("h_common_cluster_t_vs_eta",
                                               "Cluster Time vs Eta (Common Selection, pT>10 GeV);Cluster #eta;Cluster Time (ns)",
                                               50, -1.0, 1.0, 160, -40, 40);

    // Histogram 5: cluster time vs. cluster eta for clusters with weta_cogx>1.0 and pT>10
    TH2D *h_weta1p0_cluster_t_vs_eta = new TH2D("h_weta1p0_cluster_t_vs_eta",
                                                "Cluster Time vs Eta (weta_cogx>1.0, pT>10 GeV);Cluster #eta;Cluster Time (ns)",
                                                50, -1.0, 1.0, 160, -40, 40);

    // Histogram 6: cluster time vs. cluster eta for clusters with NPB selection and pT>10
    TH2D *h_npb_cluster_t_vs_eta = new TH2D("h_npb_cluster_t_vs_eta",
                                            "Cluster Time vs Eta (NPB Selection, pT>10 GeV);Cluster #eta;Cluster Time (ns)",
                                            50, -1.0, 1.0, 160, -40, 40);

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

        // Calculate MBD mean time with correction
        float mbd_mean_time = *mbd_time;
        float mbdoffset = 0;
        if (mbd_t0_correction.find(*runnumber) != mbd_t0_correction.end())
        {
            mbdoffset = mbd_t0_correction[*runnumber];
        }
        if (issim)
        {
            mbdoffset = 0.0;
        }
        mbd_mean_time = mbd_mean_time - mbdoffset;

        // Loop over clusters
        for (int icluster = 0; icluster < *ncluster; icluster++)
        {
            // pT > 10 GeV requirement
            if (cluster_Et[icluster] < 10.0)
                continue;

            // Calculate cluster time
            float clusteravgtime = 0;
            float cluster_total_e = 0;
            for (int i = 0; i < 49; i++)
            {
                if (cluster_ownership_array[icluster * 49 + i] == 1)
                {
                    int status = cluster_status_array[icluster * 49 + i];
                    // Check bit 5 is set then it is ZS tower
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

            float delta_t = clusteravgtime - mbd_mean_time;
            float cluster_eta = cluster_Eta[icluster];

            // Calculate e11_over_e33
            float e11_over_e33 = cluster_e11[icluster] / cluster_e33[icluster];

            // Check for b2b jet if required, and also check for otherside_jet for NPB selection
            bool passes_common_b2bjet = true;
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

                // Check for opposite side jet (for NPB selection)
                if (std::abs(dphi) > (M_PI / 2))
                {
                    otherside_jet = true;
                }

                // Check for b2b jet (for common selection)
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

            // Fill histograms 1 & 4: common selection with pT>10
            if (common_pass)
            {
                h_common_delta_t_vs_eta->Fill(cluster_eta, delta_t);
                h_common_cluster_t_vs_eta->Fill(cluster_eta, clusteravgtime);
            }

            // Fill histograms 2 & 5: weta_cogx>1.0 with pT>10
            if (cluster_weta_cogx[icluster] > 1.0)
            {
                h_weta1p0_delta_t_vs_eta->Fill(cluster_eta, delta_t);
                h_weta1p0_cluster_t_vs_eta->Fill(cluster_eta, clusteravgtime);
            }

            // NPB selection: no opposite jet + weta_cogx >= 0.6
            float npb_weta_min = 0.6;
            //bool badtime = True;
            bool isnpb = !otherside_jet && (cluster_weta_cogx[icluster] >= npb_weta_min);

            // Fill histograms 3 & 6: NPB selection with pT>10
            if (isnpb)
            {
                h_npb_delta_t_vs_eta->Fill(cluster_eta, delta_t);
                h_npb_cluster_t_vs_eta->Fill(cluster_eta, clusteravgtime);
            }
        }

        ientry++;
    }

    // Get entries before closing file (histograms are deleted after file close)
    Long64_t entries_common = h_common_delta_t_vs_eta->GetEntries();
    Long64_t entries_weta = h_weta1p0_delta_t_vs_eta->GetEntries();
    Long64_t entries_npb = h_npb_delta_t_vs_eta->GetEntries();
    Long64_t entries_common_cluster_t = h_common_cluster_t_vs_eta->GetEntries();
    Long64_t entries_weta_cluster_t = h_weta1p0_cluster_t_vs_eta->GetEntries();
    Long64_t entries_npb_cluster_t = h_npb_cluster_t_vs_eta->GetEntries();

    // Write and close
    fout->cd();
    h_common_delta_t_vs_eta->Write();
    h_weta1p0_delta_t_vs_eta->Write();
    h_npb_delta_t_vs_eta->Write();
    h_common_cluster_t_vs_eta->Write();
    h_weta1p0_cluster_t_vs_eta->Write();
    h_npb_cluster_t_vs_eta->Write();
    fout->Close();
    delete fout;

    std::cout << "Output written to: " << outfilename << std::endl;
    std::cout << "Histogram 1: Common selection (Cluster-MBD time, pT>10) - " << entries_common << " entries" << std::endl;
    std::cout << "Histogram 2: weta_cogx>1.0 (Cluster-MBD time, pT>10) - " << entries_weta << " entries" << std::endl;
    std::cout << "Histogram 3: NPB selection (Cluster-MBD time, pT>10) - " << entries_npb << " entries" << std::endl;
    std::cout << "Histogram 4: Common selection (Cluster time, pT>10) - " << entries_common_cluster_t << " entries" << std::endl;
    std::cout << "Histogram 5: weta_cogx>1.0 (Cluster time, pT>10) - " << entries_weta_cluster_t << " entries" << std::endl;
    std::cout << "Histogram 6: NPB selection (Cluster time, pT>10) - " << entries_npb_cluster_t << " entries" << std::endl;
}
