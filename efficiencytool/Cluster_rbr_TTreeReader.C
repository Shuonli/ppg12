#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <yaml-cpp/yaml.h>

//void Cluster_rbr_TTreeReader(const std::string &configname = "config_bdt_none_split.yaml", const std::string filetype = "data")
void Cluster_rbr_TTreeReader(const std::string &configname = "config_nom468.yaml", const std::string filetype = "data")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    bool issim = true;
    if (filetype == "data")
    {
        issim = false;
    }

    std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
    std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;
    if (!issim)
    {
        infilename = configYaml["input"]["data_file"].as<std::string>();
    }

    // read in lumi file
    std::string lumifilename = "/sphenix/user/jocl/projects/LumiList/list_468_RN_18_26_30_2.list";

    std::map<int, float> lumivals;
    std::map<int, float> lumivals_nc;
    std:map<int, float> event_counts;
    std::map<int, float> all_clusters;
    std::map<int, float> common_clusters;
    std::map<int, float> tight_iso_clusters;
    std::map<int, float> tight_noniso_clusters;
    std::map<int, float> nontight_iso_clusters;
    std::map<int, float> nontight_noniso_clusters;

    std::ifstream lumifile(lumifilename);
    if (lumifile.is_open())
    {
        std::string line;
        while (std::getline(lumifile, line))
        {
            std::istringstream iss(line);
            int runnumber;
            float lumi;
            float lumi_nc;
            float bla;
            if (!(iss >> runnumber >> bla >> bla >> lumi >> lumi_nc >> bla >> bla))
            {
                continue;
            }
            lumivals[runnumber] = lumi;
            lumivals_nc[runnumber] = lumi_nc;
            std::cout << "runnumber: " << runnumber << " lumi: " << lumi << std::endl;
        }
    }

    // Build input chain and connect reader
    std::string treename = configYaml["input"]["tree"].as<std::string>();
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());

    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();

    int iso_threshold = configYaml["analysis"]["iso_threshold"].as<int>(0);
    int n_nt_fail = configYaml["analysis"]["n_nt_fail"].as<int>(1);
    int weta_fail = configYaml["analysis"]["weta_fail"].as<int>(0);
    int wphi_fail = configYaml["analysis"]["wphi_fail"].as<int>(0);
    int e11_to_e33_fail = configYaml["analysis"]["e11_to_e33_fail"].as<int>(0);
    int e32_to_e35_fail = configYaml["analysis"]["e32_to_e35_fail"].as<int>(0);
    int et1_fail = configYaml["analysis"]["et1_fail"].as<int>(0);

    int weta_on = configYaml["analysis"]["weta_on"].as<int>(1);
    int wphi_on = configYaml["analysis"]["wphi_on"].as<int>(1);
    int e11_to_e33_on = configYaml["analysis"]["e11_to_e33_on"].as<int>(1);
    int e32_to_e35_on = configYaml["analysis"]["e32_to_e35_on"].as<int>(1);
    int et1_on = configYaml["analysis"]["et1_on"].as<int>(1);
    int et2_on = configYaml["analysis"]["et2_on"].as<int>(1);
    int et3_on = configYaml["analysis"]["et3_on"].as<int>(1);
    int et4_on = configYaml["analysis"]["et3_on"].as<int>(1);
    int nosat = configYaml["analysis"]["nosat"].as<int>(0);

    float truthisocut = configYaml["analysis"]["truth_iso_max"].as<float>();
    float recoiso_min = configYaml["analysis"]["reco_iso_min"].as<float>();
    float recoiso_max_b = configYaml["analysis"]["reco_iso_max_b"].as<float>();
    float recoiso_max_s = configYaml["analysis"]["reco_iso_max_s"].as<float>();
    float recononiso_min_shift = configYaml["analysis"]["reco_noniso_min_shift"].as<float>();
    float recononiso_max = configYaml["analysis"]["reco_noniso_max"].as<float>();
    float vertexcut = configYaml["analysis"]["vertex_cut"].as<float>();
    std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();
    std::vector<float> pT_bins = configYaml["analysis"]["pT_bins"].as<std::vector<float>>();
    int n_pT_bins = pT_bins.size() - 1;
    double pT_bin_edges[n_pT_bins + 1];
    std::copy(pT_bins.begin(), pT_bins.end(), pT_bin_edges);
    std::vector<float> pT_bins_truth = configYaml["analysis"]["pT_bins_truth"].as<std::vector<float>>();
    int n_pT_bins_truth = pT_bins_truth.size() - 1;
    double pT_bin_edges_truth[n_pT_bins_truth + 1];
    double pTmin_truth = pT_bins_truth[0];
    double pTmax_truth = pT_bins_truth[n_pT_bins_truth];
    std::cout << "n_pT_bins_truth: " << n_pT_bins_truth << std::endl;
    for (int i = 0; i < n_pT_bins_truth + 1; i++)
    {
        pT_bin_edges_truth[i] = pT_bins_truth[i];
        std::cout << "pT_bin_edges_truth: " << pT_bin_edges_truth[i] << std::endl;
    }

    int conesize = configYaml["analysis"]["cone_size"].as<int>();
    float reco_min_ET = configYaml["analysis"]["reco_min_ET"].as<float>();
    float eff_dR = configYaml["analysis"]["eff_dR"].as<float>();
    int trigger_used = configYaml["analysis"]["trigger_used"].as<int>();

    // tight cuts
    float tight_weta_cogx_max = configYaml["analysis"]["tight"]["weta_cogx_max"].as<float>();
    float tight_weta_cogx_min = configYaml["analysis"]["tight"]["weta_cogx_min"].as<float>();
    float tight_weta_cogx_max_b = configYaml["analysis"]["tight"]["weta_cogx_max_b"].as<float>();
    float tight_weta_cogx_max_s = configYaml["analysis"]["tight"]["weta_cogx_max_s"].as<float>();
    float tight_wphi_cogx_max = configYaml["analysis"]["tight"]["wphi_cogx_max"].as<float>();
    float tight_wphi_cogx_min = configYaml["analysis"]["tight"]["wphi_cogx_min"].as<float>();
    float tight_wphi_cogx_max_b = configYaml["analysis"]["tight"]["wphi_cogx_max_b"].as<float>();
    float tight_wphi_cogx_max_s = configYaml["analysis"]["tight"]["wphi_cogx_max_s"].as<float>();
    float tight_e11_over_e33_max = configYaml["analysis"]["tight"]["e11_over_e33_max"].as<float>();
    float tight_e11_over_e33_min = configYaml["analysis"]["tight"]["e11_over_e33_min"].as<float>();
    float tight_et1_max = configYaml["analysis"]["tight"]["et1_max"].as<float>();
    float tight_et1_min = configYaml["analysis"]["tight"]["et1_min"].as<float>();
    float tight_et1_min_b = configYaml["analysis"]["tight"]["et1_min_b"].as<float>();
    float tight_et1_min_s = configYaml["analysis"]["tight"]["et1_min_s"].as<float>();
    float tight_et2_max = configYaml["analysis"]["tight"]["et2_max"].as<float>(1.0);
    float tight_et2_min = configYaml["analysis"]["tight"]["et2_min"].as<float>(0.0);
    float tight_et3_max = configYaml["analysis"]["tight"]["et3_max"].as<float>(1.0);
    float tight_et3_min = configYaml["analysis"]["tight"]["et3_min"].as<float>(0.0);
    float tight_e32_over_e35_max = configYaml["analysis"]["tight"]["e32_over_e35_max"].as<float>();
    float tight_e32_over_e35_min = configYaml["analysis"]["tight"]["e32_over_e35_min"].as<float>();
    float tight_prob_max = configYaml["analysis"]["tight"]["prob_max"].as<float>();
    float tight_prob_min = configYaml["analysis"]["tight"]["prob_min"].as<float>();
    float tight_et4_max = configYaml["analysis"]["tight"]["et4_max"].as<float>();
    float tight_et4_min = configYaml["analysis"]["tight"]["et4_min"].as<float>();

    // non tight cuts
    float non_tight_weta_cogx_max = configYaml["analysis"]["non_tight"]["weta_cogx_max"].as<float>();
    float non_tight_weta_cogx_min = configYaml["analysis"]["non_tight"]["weta_cogx_min"].as<float>();
    float non_tight_wphi_cogx_max = configYaml["analysis"]["non_tight"]["wphi_cogx_max"].as<float>();
    float non_tight_wphi_cogx_min = configYaml["analysis"]["non_tight"]["wphi_cogx_min"].as<float>();
    float non_tight_prob_max = configYaml["analysis"]["non_tight"]["prob_max"].as<float>();
    float non_tight_prob_min = configYaml["analysis"]["non_tight"]["prob_min"].as<float>();
    float non_tight_et1_max = configYaml["analysis"]["non_tight"]["et1_max"].as<float>();
    float non_tight_et1_min = configYaml["analysis"]["non_tight"]["et1_min"].as<float>();
    float non_tight_e11_over_e33_max = configYaml["analysis"]["non_tight"]["e11_over_e33_max"].as<float>();
    float non_tight_e11_over_e33_min = configYaml["analysis"]["non_tight"]["e11_over_e33_min"].as<float>();
    float non_tight_e32_over_e35_max = configYaml["analysis"]["non_tight"]["e32_over_e35_max"].as<float>();
    float non_tight_e32_over_e35_min = configYaml["analysis"]["non_tight"]["e32_over_e35_min"].as<float>();
    float non_tight_et4_max = configYaml["analysis"]["non_tight"]["et4_max"].as<float>();
    float non_tight_et4_min = configYaml["analysis"]["non_tight"]["et4_min"].as<float>();

    // common cuts
    float common_prob_max = configYaml["analysis"]["common"]["prob_max"].as<float>();
    float common_prob_min = configYaml["analysis"]["common"]["prob_min"].as<float>();
    float common_e11_over_e33_max = configYaml["analysis"]["common"]["e11_over_e33_max"].as<float>();
    float common_e11_over_e33_min = configYaml["analysis"]["common"]["e11_over_e33_min"].as<float>();
    float common_wr_cogx_bound = configYaml["analysis"]["common"]["wr_cogx_bound"].as<float>();
    float common_cluster_weta_cogx_bound = configYaml["analysis"]["common"]["cluster_weta_cogx_bound"].as<float>();

    // Reader setup
    TTreeReader reader(&chain);
    TTreeReaderValue<int> mbdnorthhit(reader, "mbdnorthhit");
    TTreeReaderValue<int> mbdsouthhit(reader, "mbdsouthhit");
    TTreeReaderValue<float> vertexz(reader, "vertexz");
    TTreeReaderArray<Bool_t> scaledtrigger(reader, "scaledtrigger");
    TTreeReaderValue<int> runnumber(reader, "runnumber");
    TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));

    TTreeReaderArray<float> cluster_E(reader, Form("cluster_E_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_prob(reader, Form("cluster_prob_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_02(reader, Form("cluster_iso_02_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03(reader, Form("cluster_iso_03_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_04(reader, Form("cluster_iso_04_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et1(reader, Form("cluster_et1_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et2(reader, Form("cluster_et2_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et3(reader, Form("cluster_et3_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et4(reader, Form("cluster_et4_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_weta_cogx(reader, Form("cluster_weta_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_wphi_cogx(reader, Form("cluster_wphi_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e11(reader, Form("cluster_e11_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e33(reader, Form("cluster_e33_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e35(reader, Form("cluster_e35_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e32(reader, Form("cluster_e32_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_emcal(reader, Form("cluster_iso_03_60_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_hcalin(reader, Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_hcalout(reader, Form("cluster_iso_03_60_hcalout_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_nsaturated(reader, Form("cluster_nsaturated_%s", clusternodename.c_str()));

    std::set<int> skiprunnumbers = {47698, 51489, 51721, 51725, 53284};

    TH1D *h_allclusterET = new TH1D("h_allclusterET", "h_allclusterET", 100, 0, 100);
    TH1D *h_commonclusterET = new TH1D("h_commonclusterET", "h_commonclusterET", 100, 0, 100);
    TH1D *h_vertexz = new TH1D("h_vertexz", "h_vertexz", 100, -100, 100);
    TH1D *h_n_cluster = new TH1D("h_n_cluster", "h_n_cluster", 100, 0, 100);

    int nentries = chain.GetEntries();
    int ientry = 0;
    while (reader.Next())
    {
        if (ientry % 100000 == 0)
        {
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        }

        if (!issim)
        {
            if (skiprunnumbers.find(*runnumber) != skiprunnumbers.end())
            {
                ientry++;
                continue;
            }
            if (scaledtrigger.GetSize() <= trigger_used || scaledtrigger[trigger_used] == 0)
            {
                ientry++;
                continue;
            }
        }

        if (std::abs(*vertexz) > vertexcut)
        {
            ientry++;
            continue;
        }
        if (!(*mbdnorthhit >= 1 && *mbdsouthhit >= 1))
        {
            ientry++;
            continue;
        }
        // if n cluster is 0, skip
        if (*ncluster == 0)
        {
            ientry++;
            continue;
        }
        event_counts[*runnumber] += 1;

        h_vertexz->Fill(*vertexz);
        h_n_cluster->Fill(*ncluster);
        for (int icluster = 0; icluster < *ncluster; icluster++)
        {
            if (cluster_Et[icluster] < 10)
                continue;
            if (cluster_Et[icluster] > 26)
                continue;
            if (std::abs(cluster_Eta[icluster]) > 0.7)
                continue;

            h_allclusterET->Fill(cluster_Et[icluster]);

            float e11_over_e33 = cluster_e11[icluster] / cluster_e33[icluster];
            float e32_over_e35 = cluster_e32[icluster] / cluster_e35[icluster];

            float recoisoET = -999;
            if (conesize == 4)
            {
                recoisoET = cluster_iso_04[icluster];
            }
            else if (conesize == 3)
            {
                recoisoET = cluster_iso_03[icluster];
            }
            else if (conesize == 2)
            {
                recoisoET = cluster_iso_02[icluster];
            }
            else
            {
                std::cout << "Error: conesize not supported" << std::endl;
                continue;
            }

            if (iso_threshold)
            {
                recoisoET = cluster_iso_03_60_emcal[icluster] + cluster_iso_03_60_hcalin[icluster] + cluster_iso_03_60_hcalout[icluster];
            }

            bool common_pass = false;
            bool tight = false;
            bool nontight = false;
            bool iso = false;
            bool noniso = false;

            float clusterET = cluster_Et[icluster];
            float recoiso_max = recoiso_max_b + recoiso_max_s * clusterET;
            float recononiso_min = recoiso_max + recononiso_min_shift;
            if (recoisoET > recoiso_min && recoisoET < recoiso_max)
            {
                iso = true;
            }
            if (recoisoET > recononiso_min && recoisoET < recononiso_max)
            {
                noniso = true;
            }

            if (cluster_prob[icluster] > common_prob_min &&
                cluster_prob[icluster] < common_prob_max &&
                e11_over_e33 > common_e11_over_e33_min &&
                e11_over_e33 < common_e11_over_e33_max &&
                (cluster_weta_cogx[icluster] < common_cluster_weta_cogx_bound))
            {
                common_pass = true;
            }

            if (common_pass)
            {
                tight_weta_cogx_max = tight_weta_cogx_max_b + tight_weta_cogx_max_s * clusterET;
                tight_wphi_cogx_max = tight_wphi_cogx_max_b + tight_wphi_cogx_max_s * clusterET;
                tight_et1_min = tight_et1_min_b + tight_et1_min_s * clusterET;

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

                if (is_cluster_weta_cogx_tight &&
                    is_cluster_wphi_cogx_tight &&
                    is_cluster_et1_tight &&
                    is_cluster_et2_tight &&
                    is_cluster_et3_tight &&
                    is_e11_over_e33_tight &&
                    is_e32_over_e35_tight &&
                    is_cluster_et4_tight &&
                    is_cluster_prob_tight)
                {
                    tight = true;
                }

                if (
                    cluster_weta_cogx[icluster] > non_tight_weta_cogx_min &&
                    cluster_weta_cogx[icluster] < non_tight_weta_cogx_max &&
                    cluster_wphi_cogx[icluster] > non_tight_wphi_cogx_min &&
                    cluster_wphi_cogx[icluster] < non_tight_wphi_cogx_max &&
                    cluster_prob[icluster] > non_tight_prob_min &&
                    cluster_prob[icluster] < non_tight_prob_max &&
                    e11_over_e33 > non_tight_e11_over_e33_min &&
                    e11_over_e33 < non_tight_e11_over_e33_max &&
                    e32_over_e35 > non_tight_e32_over_e35_min &&
                    e32_over_e35 < non_tight_e32_over_e35_max &&
                    cluster_et1[icluster] > non_tight_et1_min &&
                    cluster_et1[icluster] < non_tight_et1_max &&
                    cluster_et4[icluster] > non_tight_et4_min &&
                    cluster_et4[icluster] < non_tight_et4_max)
                {
                    int nfail = 0;
                    if (!is_cluster_weta_cogx_tight) nfail += weta_on;
                    if (!is_cluster_wphi_cogx_tight) nfail += wphi_on;
                    if (!is_cluster_et1_tight) nfail += et1_on;
                    if (!is_cluster_et2_tight) nfail += et2_on;
                    if (!is_cluster_et3_tight) nfail += et3_on;
                    if (!is_e11_over_e33_tight) nfail += e11_to_e33_on;
                    if (!is_e32_over_e35_tight) nfail += e32_to_e35_on;
                    if (!is_cluster_et4_tight) nfail += et4_on;
                    if (!is_cluster_prob_tight) nfail += 1;

                    if ((nfail > n_nt_fail))
                    {
                        bool all_flags_fail = true;
                        if (weta_fail && is_cluster_weta_cogx_tight) all_flags_fail = false;
                        if (wphi_fail && is_cluster_wphi_cogx_tight) all_flags_fail = false;
                        if (et1_fail && is_cluster_et1_tight) all_flags_fail = false;
                        if (e11_to_e33_fail && is_e11_over_e33_tight) all_flags_fail = false;
                        if (e32_to_e35_fail && is_e32_over_e35_tight) all_flags_fail = false;

                        if (all_flags_fail)
                        {
                            nontight = true;
                        }
                    }
                }

                common_clusters[*runnumber] += 1;
                h_commonclusterET->Fill(clusterET);
                if (tight && iso)
                {
                    tight_iso_clusters[*runnumber] += 1;
                }
                if (tight && noniso)
                {
                    tight_noniso_clusters[*runnumber] += 1;
                }
                if (nontight && iso)
                {
                    nontight_iso_clusters[*runnumber] += 1;
                }
                if (nontight && noniso)
                {
                    nontight_noniso_clusters[*runnumber] += 1;
                }
            }
        }

        ientry++;
    }

    // print out the results and write graphs
    std::vector<double> x_common, x_common_error, y_common, y_common_error;
    std::vector<double> x_tight_iso, x_tight_iso_error, y_tight_iso, y_tight_iso_error;
    std::vector<double> x_tight_noniso, x_tight_noniso_error, y_tight_noniso, y_tight_noniso_error;
    std::vector<double> x_nontight_iso, x_nontight_iso_error, y_nontight_iso, y_nontight_iso_error;
    std::vector<double> x_nontight_noniso, x_nontight_noniso_error, y_nontight_noniso, y_nontight_noniso_error;
    std::vector<double> x_corr, y_corr, x_corr_error, y_corr_error;
    std::vector<double> x_common_cluster_event, y_common_cluster_event, x_common_cluster_event_error, y_common_cluster_event_error;

    TH1D* h_n_common_cluster_event = new TH1D("h_n_common_cluster_event", "Number of common clusters per event", 100, 0, 0.5);

    float total_lumi = 0;
    for (auto const &entry : common_clusters)
    {
        float lumi = lumivals[entry.first];
        std::cout << entry.first << " " << lumi << std::endl;

        if ((entry.second < 10) || lumi == 0)
        {
            if ((entry.second > 0) && lumi == 0)
            {
                std::cout << entry.first << std::endl;
            }
            common_clusters[entry.first] = 0;
            tight_iso_clusters[entry.first] = 0;
            tight_noniso_clusters[entry.first] = 0;
            nontight_iso_clusters[entry.first] = 0;
            nontight_noniso_clusters[entry.first] = 0;
        }
        else
        {
            h_n_common_cluster_event->Fill(common_clusters[entry.first]/event_counts[entry.first]);
            x_common_cluster_event.push_back(entry.first);
            x_common_cluster_event_error.push_back(0);
            y_common_cluster_event.push_back(common_clusters[entry.first]/event_counts[entry.first]);
            y_common_cluster_event_error.push_back(0);
            total_lumi += lumi;

            float common_clusters_error = std::sqrt(common_clusters[entry.first]);
            common_clusters[entry.first] = common_clusters[entry.first] / lumi;
            common_clusters_error = common_clusters_error / lumi;
            x_common.push_back(entry.first);
            x_common_error.push_back(0);
            y_common.push_back(common_clusters[entry.first]);
            y_common_error.push_back(common_clusters_error);

            float tight_iso_clusters_error = std::sqrt(tight_iso_clusters[entry.first]);
            tight_iso_clusters[entry.first] = tight_iso_clusters[entry.first] / lumi;
            tight_iso_clusters_error = tight_iso_clusters_error / lumi;
            x_tight_iso.push_back(entry.first);
            x_tight_iso_error.push_back(0);
            y_tight_iso.push_back(tight_iso_clusters[entry.first]);
            y_tight_iso_error.push_back(tight_iso_clusters_error);

            float tight_noniso_clusters_error = std::sqrt(tight_noniso_clusters[entry.first]);
            tight_noniso_clusters[entry.first] = tight_noniso_clusters[entry.first] / lumi;
            tight_noniso_clusters_error = tight_noniso_clusters_error / lumi;
            x_tight_noniso.push_back(entry.first);
            x_tight_noniso_error.push_back(0);
            y_tight_noniso.push_back(tight_noniso_clusters[entry.first]);
            y_tight_noniso_error.push_back(tight_noniso_clusters_error);

            float nontight_iso_clusters_error = std::sqrt(nontight_iso_clusters[entry.first]);
            nontight_iso_clusters[entry.first] = nontight_iso_clusters[entry.first] / lumi;
            nontight_iso_clusters_error = nontight_iso_clusters_error / lumi;
            x_nontight_iso.push_back(entry.first);
            x_nontight_iso_error.push_back(0);
            y_nontight_iso.push_back(nontight_iso_clusters[entry.first]);
            y_nontight_iso_error.push_back(nontight_iso_clusters_error);

            float nontight_noniso_clusters_error = std::sqrt(nontight_noniso_clusters[entry.first]);
            nontight_noniso_clusters[entry.first] = nontight_noniso_clusters[entry.first] / lumi;
            nontight_noniso_clusters_error = nontight_noniso_clusters_error / lumi;
            x_nontight_noniso.push_back(entry.first);
            x_nontight_noniso_error.push_back(0);
            y_nontight_noniso.push_back(nontight_noniso_clusters[entry.first]);
            y_nontight_noniso_error.push_back(nontight_noniso_clusters_error);

            float correction = lumivals[entry.first] / lumivals_nc[entry.first];
            x_corr.push_back(entry.first);
            x_corr_error.push_back(0);
            y_corr.push_back(correction);
            y_corr_error.push_back(0);
        }
    }

    std::cout << "Total lumi: " << total_lumi << std::endl;

    TFile *fout = new TFile("rbrQA.root", "RECREATE");

    TGraphErrors *gr_common = new TGraphErrors(x_common.size(), &x_common[0], &y_common[0], &x_common_error[0], &y_common_error[0]);
    gr_common->SetName("gr_common");

    TGraphErrors *gr_tight_iso = new TGraphErrors(x_tight_iso.size(), &x_tight_iso[0], &y_tight_iso[0], &x_tight_iso_error[0], &y_tight_iso_error[0]);
    gr_tight_iso->SetName("gr_tight_iso");

    TGraphErrors *gr_tight_noniso = new TGraphErrors(x_tight_noniso.size(), &x_tight_noniso[0], &y_tight_noniso[0], &x_tight_noniso_error[0], &y_tight_noniso_error[0]);
    gr_tight_noniso->SetName("gr_tight_noniso");

    TGraphErrors *gr_nontight_iso = new TGraphErrors(x_nontight_iso.size(), &x_nontight_iso[0], &y_nontight_iso[0], &x_nontight_iso_error[0], &y_nontight_iso_error[0]);
    gr_nontight_iso->SetName("gr_nontight_iso");

    TGraphErrors *gr_nontight_noniso = new TGraphErrors(x_nontight_noniso.size(), &x_nontight_noniso[0], &y_nontight_noniso[0], &x_nontight_noniso_error[0], &y_nontight_noniso_error[0]);
    gr_nontight_noniso->SetName("gr_nontight_noniso");

    TGraphErrors *gr_corr = new TGraphErrors(x_corr.size(), &x_corr[0], &y_corr[0], &x_corr_error[0], &y_corr_error[0]);
    gr_corr->SetName("gr_corr");

    TGraphErrors *gr_common_cluster_event = new TGraphErrors(x_common_cluster_event.size(), &x_common_cluster_event[0], &y_common_cluster_event[0], &x_common_cluster_event_error[0], &y_common_cluster_event_error[0]);
    gr_common_cluster_event->SetName("gr_common_cluster_event");

    gr_common->Draw("AP");
    //gr_common_cluster_event->Draw("P same");

    fout->cd();
    gr_common->Write();
    gr_tight_iso->Write();
    gr_tight_noniso->Write();
    gr_nontight_iso->Write();
    gr_nontight_noniso->Write();
    gr_corr->Write();
    h_allclusterET->Write();
    h_commonclusterET->Write();
    h_vertexz->Write();
    h_n_common_cluster_event->Write();
    gr_common_cluster_event->Write();
    h_n_cluster->Write();
    fout->Write();
    fout->Close();
}


