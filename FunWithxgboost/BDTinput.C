#include <yaml-cpp/yaml.h>
#include "../efficiencytool/CrossSectionWeights.h"
using namespace PPG12;

// STL
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// ROOT
#include <TChain.h>
#include <TString.h> // Form
#include <TSystem.h> // gSystem
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

// math
#include <cmath>

// EMCAL time sample period used for NPB timing calculation
const float TIME_SAMPLE_NS = 17.6;

void BDTinput(const std::string &configname = "config_nom.yaml", const std::string filetype = "data")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    bool issim = true;

    if (filetype == "data")
    {
        issim = false;
    }
    bool isbackground = false;
    float max_photon_lower = 0;
    float max_photon_upper = 100;
    float max_jet_lower = 0;
    float max_jet_upper = 100;

    float energy_scale_lower = 0;
    float energy_scale_upper = 100;

    float cluster_ET_upper = 100;

    float weight = 1.0;
    float vertex_weight = 1.0;
    float cross_weight = 1.0;

    if (filetype == "photon5")
    {
        max_photon_lower = 0;
        // max_photon_upper = 14;
        //  max_photon_upper = 200;
        weight = photon5cross / photon20cross;
    }
    else if (filetype == "photon10")
    {
        max_photon_lower = 14;
        // max_photon_upper = 30;

        // max_photon_lower = 0;
        // max_photon_upper = 200;
        weight = photon10cross / photon20cross;
    }
    else if (filetype == "photon20")
    {
        max_photon_lower = 30;
        // max_photon_lower = 0;
        max_photon_upper = 200;
        weight = 1.0;
    }
    else if (filetype == "jet10")
    {
        max_jet_lower = 10;
        // max_jet_upper = 19;
        energy_scale_lower = 10;
        energy_scale_upper = 16;
        cluster_ET_upper = 25;
        weight = jet10cross / jet30cross;
        isbackground = true;
    }
    else if (filetype == "jet15")
    {
        max_jet_lower = 19;
        // max_jet_upper = 23;
        energy_scale_lower = 16;
        energy_scale_upper = 20;
        cluster_ET_upper = 25;
        weight = jet15cross / jet30cross;
        isbackground = true;
    }
    else if (filetype == "jet20")
    {
        max_jet_lower = 23;
        // max_jet_upper = 30;
        energy_scale_lower = 20;
        energy_scale_upper = 30;
        weight = jet20cross / jet30cross;
        isbackground = true;
    }
    else if (filetype == "jet30")
    {
        max_jet_lower = 30;
        max_jet_upper = 100;
        energy_scale_lower = 30;
        energy_scale_upper = 100;
        weight = 1.0;
        isbackground = true;
    }
    else if (filetype == "jet50")
    {
        max_jet_lower = 50;
        max_jet_upper = 70;
        weight = jet50cross / jet30cross;
        isbackground = true;
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
    // Build input chain (same pattern as RecoEffCalculator_TTreeReader.C)
    std::string treename = configYaml["input"]["tree"].as<std::string>();
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());
    TTreeReader reader(&chain);

    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();

    int iso_threshold = configYaml["analysis"]["iso_threshold"].as<int>(0);
    // If enabled, write MC clusters to the txt output even when they fail `common_pass`.
    // This is useful for NPB-score / preselection studies without photon-ID bias.
    int write_outside_common_pass_mc = configYaml["analysis"]["write_outside_common_pass_mc"].as<int>(0);

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
    int et4_on = configYaml["analysis"]["et4_on"].as<int>(1);
    int nosat = configYaml["analysis"]["nosat"].as<int>(0);

    // NPB selection parameters (data only)
    int npb_cut_on = 1;
    float npb_weta_min = 0.4;
    float npb_delta_t_cut = -5.0;
    int mbd_t0_correction_on = 0;
    std::string mbd_t0_correction_file = "";
    std::map<int, float> mbd_t0_correction;

    if (!issim) {
        npb_cut_on = configYaml["analysis"]["npb_cut_on"].as<int>(0);
        npb_weta_min = configYaml["analysis"]["npb_weta_min"].as<float>(0.4);
        npb_delta_t_cut = configYaml["analysis"]["npb_delta_t_cut"].as<float>(-5.0);
        mbd_t0_correction_on = configYaml["analysis"]["mbd_t0_correction_on"].as<int>(0);
        mbd_t0_correction_file = configYaml["analysis"]["mbd_t0_correction_file"]
            .as<std::string>("/sphenix/user/shuhangli/ppg12/efficiencytool/MbdOut.corr");
    }

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
    double pTmin = pT_bins[0];
    double pTmax = pT_bins[n_pT_bins];

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

    std::copy(pT_bins.begin(), pT_bins.end(), pT_bin_edges);

    int conesize = configYaml["analysis"]["cone_size"].as<int>();

    float reco_min_ET = configYaml["analysis"]["reco_min_ET"].as<float>();

    float eff_dR = configYaml["analysis"]["eff_dR"].as<float>();

    int trigger_used = configYaml["analysis"]["trigger_used"].as<int>();

    // getting cuts from the config file
    std::cout << "tight cuts" << std::endl;
    float tight_reta77_min = configYaml["analysis"]["tight"]["reta77_min"].as<float>();
    float tight_reta77_max = configYaml["analysis"]["tight"]["reta77_max"].as<float>();

    float tight_rhad33_max = configYaml["analysis"]["tight"]["rhad33_max"].as<float>();
    float tight_rhad33_min = configYaml["analysis"]["tight"]["rhad33_min"].as<float>();

    float tight_w72_max = configYaml["analysis"]["tight"]["w72_max"].as<float>();
    float tight_w72_min = configYaml["analysis"]["tight"]["w72_min"].as<float>();

    float tight_re11_E_max = configYaml["analysis"]["tight"]["re11_E_max"].as<float>();
    float tight_re11_E_min = configYaml["analysis"]["tight"]["re11_E_min"].as<float>();

    float tight_CNN_min = configYaml["analysis"]["tight"]["CNN_min"].as<float>();
    float tight_CNN_max = configYaml["analysis"]["tight"]["CNN_max"].as<float>();

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

    float tight_w32_max = configYaml["analysis"]["tight"]["w32_max"].as<float>();
    float tight_w32_min = configYaml["analysis"]["tight"]["w32_min"].as<float>();

    // non tight cuts
    std::cout << "non tight cuts" << std::endl;
    float non_tight_reta77_min = configYaml["analysis"]["non_tight"]["reta77_min"].as<float>();
    float non_tight_reta77_max = configYaml["analysis"]["non_tight"]["reta77_max"].as<float>();

    float non_tight_rhad33_max = configYaml["analysis"]["non_tight"]["rhad33_max"].as<float>();
    float non_tight_rhad33_min = configYaml["analysis"]["non_tight"]["rhad33_min"].as<float>();

    float non_tight_w72_max = configYaml["analysis"]["non_tight"]["w72_max"].as<float>();
    float non_tight_w72_min = configYaml["analysis"]["non_tight"]["w72_min"].as<float>();

    float non_tight_re11_E_max = configYaml["analysis"]["non_tight"]["re11_E_max"].as<float>();
    float non_tight_re11_E_min = configYaml["analysis"]["non_tight"]["re11_E_min"].as<float>();

    float non_tight_CNN_min = configYaml["analysis"]["non_tight"]["CNN_min"].as<float>();
    float non_tight_CNN_max = configYaml["analysis"]["non_tight"]["CNN_max"].as<float>();

    float non_tight_weta_cogx_max = configYaml["analysis"]["non_tight"]["weta_cogx_max"].as<float>();
    float non_tight_weta_cogx_min = configYaml["analysis"]["non_tight"]["weta_cogx_min"].as<float>();
    float non_tight_weta_cogx_max_b = configYaml["analysis"]["non_tight"]["weta_cogx_max_b"].as<float>();
    float non_tight_weta_cogx_max_s = configYaml["analysis"]["non_tight"]["weta_cogx_max_s"].as<float>();

    float non_tight_wphi_cogx_max = configYaml["analysis"]["non_tight"]["wphi_cogx_max"].as<float>();
    float non_tight_wphi_cogx_min = configYaml["analysis"]["non_tight"]["wphi_cogx_min"].as<float>();
    float non_tight_wphi_cogx_max_b = configYaml["analysis"]["non_tight"]["wphi_cogx_max_b"].as<float>();
    float non_tight_wphi_cogx_max_s = configYaml["analysis"]["non_tight"]["wphi_cogx_max_s"].as<float>();

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

    float non_tight_w32_max = configYaml["analysis"]["non_tight"]["w32_max"].as<float>();
    float non_tight_w32_min = configYaml["analysis"]["non_tight"]["w32_min"].as<float>();

    // common cuts for both tight and non tight

    float common_prob_max = configYaml["analysis"]["common"]["prob_max"].as<float>();
    float common_prob_min = configYaml["analysis"]["common"]["prob_min"].as<float>();

    float common_e11_over_e33_max = configYaml["analysis"]["common"]["e11_over_e33_max"].as<float>();
    float common_e11_over_e33_min = configYaml["analysis"]["common"]["e11_over_e33_min"].as<float>();

    float common_wr_cogx_bound = configYaml["analysis"]["common"]["wr_cogx_bound"].as<float>();
    float common_cluster_weta_cogx_bound = configYaml["analysis"]["common"]["cluster_weta_cogx_bound"].as<float>();

    int reweight = configYaml["analysis"]["unfold"]["reweight"].as<int>(); // 0 for no reweighting, 1 for reweighting

    float clusterescale = configYaml["analysis"]["cluster_escale"].as<float>(1.0);
    float clustereres = configYaml["analysis"]["cluster_eres"].as<float>(0.0);

    // Load run-by-run MBD t0 corrections for data NPB selection
    if (!issim && npb_cut_on && mbd_t0_correction_on) {
        std::ifstream file(mbd_t0_correction_file);
        if (!file.is_open()) {
            std::cerr << "[MBD t0] WARNING: cannot open " << mbd_t0_correction_file << std::endl;
        } else {
            std::string line;
            while (std::getline(file, line)) {
                if (line.empty()) continue;
                std::istringstream iss(line);
                int runnumber_temp = 0;
                float t0 = 0.0;
                if (!(iss >> runnumber_temp >> t0)) continue;
                mbd_t0_correction[runnumber_temp] = t0;
            }
            std::cout << "[MBD t0] Loaded " << mbd_t0_correction.size()
                      << " corrections from " << mbd_t0_correction_file << std::endl;
        }
    }

    // ------------------------------------------------------------------
    // TTreeReader bindings (like RecoEffCalculator_TTreeReader.C)
    // ------------------------------------------------------------------
    TTreeReaderValue<int> mbdnorthhit(reader, "mbdnorthhit");
    TTreeReaderValue<int> mbdsouthhit(reader, "mbdsouthhit");
    TTreeReaderValue<float> vertexz(reader, "vertexz");
    TTreeReaderValue<float> energy_scale(reader, "energy_scale");
    TTreeReaderValue<int> runnumber(reader, "runnumber");
    TTreeReaderArray<Bool_t> scaledtrigger(reader, "scaledtrigger");
    TTreeReaderArray<Bool_t> livetrigger(reader, "livetrigger");

    // MC truth particles (used only for MC pid labeling)
    TTreeReaderValue<int> nparticles(reader, "nparticles");
    TTreeReaderArray<int> particle_trkid(reader, "particle_trkid");
    TTreeReaderArray<int> particle_photonclass(reader, "particle_photonclass");

    // Clusters
    TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_prob(reader, Form("cluster_prob_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_02(reader, Form("cluster_iso_02_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03(reader, Form("cluster_iso_03_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_04(reader, Form("cluster_iso_04_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et1(reader, Form("cluster_et1_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et2(reader, Form("cluster_et2_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et3(reader, Form("cluster_et3_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et4(reader, Form("cluster_et4_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_weta_cogx(reader, Form("cluster_weta_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_wphi_cogx(reader, Form("cluster_wphi_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_nsaturated(reader, Form("cluster_nsaturated_%s", clusternodename.c_str()));

    // Shower-shape ingredients
    TTreeReaderArray<float> cluster_e11(reader, Form("cluster_e11_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e22(reader, Form("cluster_e22_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e13(reader, Form("cluster_e13_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e15(reader, Form("cluster_e15_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e17(reader, Form("cluster_e17_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e31(reader, Form("cluster_e31_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e51(reader, Form("cluster_e51_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e71(reader, Form("cluster_e71_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e33(reader, Form("cluster_e33_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e35(reader, Form("cluster_e35_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e37(reader, Form("cluster_e37_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e53(reader, Form("cluster_e53_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e73(reader, Form("cluster_e73_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e77(reader, Form("cluster_e77_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w32(reader, Form("cluster_w32_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w52(reader, Form("cluster_w52_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w72(reader, Form("cluster_w72_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e32(reader, Form("cluster_e32_%s", clusternodename.c_str()));

    // Isolation (used when iso_threshold is enabled)
    TTreeReaderArray<float> cluster_iso_03_60_emcal(reader, Form("cluster_iso_03_60_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_hcalin(reader, Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_hcalout(reader, Form("cluster_iso_03_60_hcalout_%s", clusternodename.c_str()));

    // Truth jets (MC)
    TTreeReaderValue<int> njet_truth(reader, "njet_truth");
    TTreeReaderArray<float> jet_truth_Pt(reader, "jet_truth_Pt");

    // Optional NPB inputs (data only when enabled)
    std::unique_ptr<TTreeReaderValue<float>> mbd_time;
    std::unique_ptr<TTreeReaderValue<int>> njet;
    std::unique_ptr<TTreeReaderArray<float>> jet_Pt;
    std::unique_ptr<TTreeReaderArray<float>> jet_Phi;
    std::unique_ptr<TTreeReaderArray<float>> cluster_time_array;
    std::unique_ptr<TTreeReaderArray<float>> cluster_e_array;
    std::unique_ptr<TTreeReaderArray<int>> cluster_ownership_array;

    if (!issim && npb_cut_on)
    {
        mbd_time = std::make_unique<TTreeReaderValue<float>>(reader, "mbd_time");
        njet = std::make_unique<TTreeReaderValue<int>>(reader, "njet");
        jet_Pt = std::make_unique<TTreeReaderArray<float>>(reader, "jet_Pt");
        jet_Phi = std::make_unique<TTreeReaderArray<float>>(reader, "jet_Phi");
        cluster_time_array = std::make_unique<TTreeReaderArray<float>>(reader, Form("cluster_time_array_%s", clusternodename.c_str()));
        cluster_e_array = std::make_unique<TTreeReaderArray<float>>(reader, Form("cluster_e_array_%s", clusternodename.c_str()));
        cluster_ownership_array = std::make_unique<TTreeReaderArray<int>>(reader, Form("cluster_ownership_array_%s", clusternodename.c_str()));
    }

    std::set<int> skiprunnumbers = {47698, 51489, 51721, 51725, 53284};
    // output text file for the showershapes
    std::ofstream outputfile;
    std::string outputname = "shapes_split_" + filetype + ".txt";

    // Override output name for data NPB
    if (!issim && npb_cut_on) {
        outputname = "shapes_split_data_npb.txt";
    }

    // Open output file for MC or data NPB
    if (issim || (!issim && npb_cut_on))
    {
        outputfile.open(outputname);

        // Write header info
        outputfile << "cluster_Et cluster_Eta cluster_Phi vertexz "
                   << "e11_over_e33 e32_over_e35 e11_over_e22 e11_over_e13 "
                   << "e11_over_e15 e11_over_e17 e11_over_e31 "
                   << "e11_over_e51 e11_over_e71 e22_over_e33 "
                   << "e22_over_e35 e22_over_e37 e22_over_e53 "
                   << "cluster_prob cluster_weta_cogx cluster_wphi_cogx "
                   << "cluster_et1 cluster_et2 cluster_et3 cluster_et4 "
                   << "cluster_w32 cluster_w52 cluster_w72 "
                   << "recoisoET "
                   << "is_tight pid " << std::endl;
    }

    const Long64_t nentries = chain.GetEntries();
    Long64_t ientry = 0;
    while (reader.Next())
    {

        if (ientry == 16152886)
        {
            ientry++;
            continue;
        }
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        std::map<int, int> particle_trkidmap;
        if (!issim)
        {
            if (skiprunnumbers.find(*runnumber) != skiprunnumbers.end())
            {
                ientry++;
                continue;
            }

            const auto ntrig = scaledtrigger.GetSize();
            if (trigger_used < 0 || (unsigned int)trigger_used >= ntrig)
            {
                ientry++;
                continue;
            }
            if (scaledtrigger[trigger_used] == 0)
            {
                ientry++;
                continue;
            }
        }
        else
        {
            for (int iparticle = 0; iparticle < *nparticles; iparticle++)
            {
                particle_trkidmap[particle_trkid[iparticle]] = iparticle;
            }
            if (isbackground)
            {
                float maxjetpT = 0;
                for (int ijet = 0; ijet < *njet_truth; ijet++)
                {
                    if (jet_truth_Pt[ijet] > maxjetpT)
                    {
                        maxjetpT = jet_truth_Pt[ijet];
                    }
                }
                /*
                if ((maxjetpT > max_jet_upper) || (maxjetpT < max_jet_lower))
                {
                    continue;
                }
                */

                // energy scale cut for now
                /*
                if ((energy_scale > energy_scale_upper) || (energy_scale < energy_scale_lower))
                {
                    continue;
                }
                */
            }
        }

        // check vertex cut
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

        for (int icluster = 0; icluster < *ncluster; icluster++)
        {
            // hard code the cut here for now :(
            if (cluster_Et[icluster] < 6)
                continue;
            if (cluster_Et[icluster] > 40)
                continue;
            if (std::abs(cluster_Eta[icluster]) > 0.7)
                continue;

            // NPB selection for data
            bool isnpb = false;
            int pid = 0;  // Default for MC

            if (!issim && npb_cut_on) {
                if (!mbd_time || !njet || !jet_Pt || !jet_Phi || !cluster_time_array || !cluster_e_array || !cluster_ownership_array)
                {
                    // missing optional branches
                    continue;
                }
                // Calculate energy-weighted cluster average time
                float clusteravgtime = 0;
                float cluster_total_e = 0;
                for (int i = 0; i < 49; i++) {
                    const int idx = icluster * 49 + i;
                    if ((*cluster_ownership_array)[idx] == 1) {
                        clusteravgtime += (*cluster_time_array)[idx] * (*cluster_e_array)[idx];
                        cluster_total_e += (*cluster_e_array)[idx];
                    }
                }
                clusteravgtime = cluster_total_e > 0 ? clusteravgtime / cluster_total_e : 0;
                clusteravgtime *= TIME_SAMPLE_NS;  // Convert to ns

                // Apply run-by-run MBD t0 correction
                float mbd_mean_time = **mbd_time;
                float mbdoffset = 0.0;
                if (mbd_t0_correction_on) {
                    auto it = mbd_t0_correction.find(*runnumber);
                    if (it != mbd_t0_correction.end()) {
                        mbdoffset = it->second;
                    }
                }
                mbd_mean_time = mbd_mean_time - mbdoffset;

                // Check "bad time" criterion
                const float delta_t_cluster_mbd = clusteravgtime - mbd_mean_time;
                bool badtime = (delta_t_cluster_mbd < npb_delta_t_cut);

                // Check for back-to-back jets (veto if present)
                bool otherside_jet = false;
                for (int ijet = 0; ijet < **njet; ijet++) {
                    if ((*jet_Pt)[ijet] < 5) continue;  // Require jet pT > 5 GeV

                    float dphi = cluster_Phi[icluster] - (*jet_Phi)[ijet];
                    while (dphi > M_PI) dphi -= 2 * M_PI;
                    while (dphi < -M_PI) dphi += 2 * M_PI;

                    if (std::abs(dphi) > (M_PI / 2)) {  // Back-to-back: Δφ > 90°
                        otherside_jet = true;
                        break;
                    }
                }

                // NPB tag
                isnpb = badtime && !otherside_jet && (cluster_weta_cogx[icluster] >= npb_weta_min);

                // Only output NPB-tagged clusters for data
                if (!isnpb) continue;

                pid = -1;  // NPB label
            }

            // Protect all ratio features against divide-by-zero to avoid inf/NaN in the txt output.
            // This is important for XGBoost training (it errors out on inf inputs).
            auto safe_div = [](float num, float den) -> float
            {
                return (std::isfinite(num) && std::isfinite(den) && den != 0.0f) ? (num / den) : 0.0f;
            };

            float e11_over_e33 = safe_div(cluster_e11[icluster], cluster_e33[icluster]);
            float e32_over_e35 = safe_div(cluster_e32[icluster], cluster_e35[icluster]);
            float e11_over_e22 = safe_div(cluster_e11[icluster], cluster_e22[icluster]);
            float e11_over_e13 = safe_div(cluster_e11[icluster], cluster_e13[icluster]);
            float e11_over_e15 = safe_div(cluster_e11[icluster], cluster_e15[icluster]);
            float e11_over_e17 = safe_div(cluster_e11[icluster], cluster_e17[icluster]);
            float e11_over_e31 = safe_div(cluster_e11[icluster], cluster_e31[icluster]);
            float e11_over_e51 = safe_div(cluster_e11[icluster], cluster_e51[icluster]);
            float e11_over_e71 = safe_div(cluster_e11[icluster], cluster_e71[icluster]);
            float e22_over_e33 = safe_div(cluster_e22[icluster], cluster_e33[icluster]);
            float e22_over_e35 = safe_div(cluster_e22[icluster], cluster_e35[icluster]);
            float e22_over_e37 = safe_div(cluster_e22[icluster], cluster_e37[icluster]);
            float e22_over_e53 = safe_div(cluster_e22[icluster], cluster_e53[icluster]);

            // Assign pid for MC as early as possible so it is available even if `common_pass` fails.
            // For data NPB, pid is already set to -1 earlier.
            if (issim)
            {
                auto it = particle_trkidmap.find(cluster_truthtrkID[icluster]);
                if (it != particle_trkidmap.end())
                {
                    pid = particle_photonclass[it->second];
                }
            }

            // reco cut
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
            // std::cout<<recoiso_max<<std::endl;
            float recononiso_min = recoiso_max + recononiso_min_shift;
            if (recoisoET > recoiso_min && recoisoET < recoiso_max)
            {
                iso = true;
            }
            if (recoisoET > recononiso_min && recoisoET < recononiso_max)
            {
                noniso = true;
            }

            // common cuts
            if (cluster_prob[icluster] > common_prob_min &&
                cluster_prob[icluster] < common_prob_max &&
                e11_over_e33 > common_e11_over_e33_min &&
                e11_over_e33 < common_e11_over_e33_max &&
                //(!(wr_cogx < common_wr_cogx_bound && cluster_weta_cogx[icluster] > common_cluster_weta_cogx_bound))
                (cluster_weta_cogx[icluster] < common_cluster_weta_cogx_bound))
            {
                common_pass = true;
            }

            if (common_pass)
            {
                {
                    tight_weta_cogx_max = tight_weta_cogx_max_b + tight_weta_cogx_max_s * clusterET;
                    tight_wphi_cogx_max = tight_wphi_cogx_max_b + tight_wphi_cogx_max_s * clusterET;
                    tight_et1_min = tight_et1_min_b + tight_et1_min_s * clusterET;
                }
                // need to update to a function to find tight and non tight
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

                // Combined condition
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
                    // fail at least one of the tight cuts with small correlation
                    int nfail = 0;
                    if (!is_cluster_weta_cogx_tight)
                    {
                        nfail += weta_on;
                    }
                    if (!is_cluster_wphi_cogx_tight)
                    {
                        nfail += wphi_on;
                    }
                    if (!is_cluster_et1_tight)
                    {
                        nfail += et1_on;
                    }
                    if (!is_cluster_et2_tight)
                    {
                        nfail += et2_on;
                    }
                    if (!is_cluster_et3_tight)
                    {
                        nfail += et3_on;
                    }
                    if (!is_e11_over_e33_tight)
                    {
                        nfail += e11_to_e33_on;
                    }
                    if (!is_e32_over_e35_tight)
                    {
                        nfail += e32_to_e35_on;
                    }
                    if (!is_cluster_et4_tight)
                    {
                        nfail += et4_on;
                    }
                    if (!is_cluster_prob_tight)
                    {
                        nfail++;
                    }

                    if ((nfail > n_nt_fail))
                    {
                        bool all_flags_fail = true;
                        if (weta_fail)
                        {
                            if (is_cluster_weta_cogx_tight)
                                all_flags_fail = false;
                        }
                        if (wphi_fail)
                        {
                            if (is_cluster_wphi_cogx_tight)
                                all_flags_fail = false;
                        }
                        if (et1_fail)
                        {
                            if (is_cluster_et1_tight)
                                all_flags_fail = false;
                        }
                        if (e11_to_e33_fail)
                        {
                            if (is_e11_over_e33_tight)
                                all_flags_fail = false;
                        }
                        if (e32_to_e35_fail)
                        {
                            if (is_e32_over_e35_tight)
                                all_flags_fail = false;
                        }

                        if (all_flags_fail)
                        {
                            nontight = true;
                        }
                    }
                    /*
                    if (
                        //!(e11_over_e33 > tight_e11_over_e33_min && e11_over_e33 < tight_e11_over_e33_max) ||
                        //!(cluster_et4[icluster] > tight_et4_min && cluster_et4[icluster] < tight_et4_max) ||
                        //!(cluster_w32[icluster] > tight_w32_min && cluster_w32[icluster] < tight_w32_max)
                        !(cluster_wphi_cogx[icluster] > tight_wphi_cogx_min && cluster_wphi_cogx[icluster] < tight_wphi_cogx_max) || !(cluster_weta_cogx[icluster] > tight_weta_cogx_min && cluster_weta_cogx[icluster] < tight_weta_cogx_max) || !(cluster_et1[icluster] > tight_et1_min && cluster_et1[icluster] < tight_et1_max))
                    {

                        nontight = true;
                    }
                    */
                }
                common_clusters[*runnumber] += 1;
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

                // Assign pid for MC (for data with NPB, pid is already set to -1)
                // (pid assignment moved earlier so it also exists when common_pass fails)

                /*
                    outputfile << "cluster_Et cluster_Eta cluster_Phi vertexz "
               << "e11_over_e33 e32_over_e35 e11_over_e22 e11_over_e13 "
               << "e11_over_e15 e11_over_e17 e11_over_e31 "
               << "e11_over_e51 e11_over_e71 e22_over_e33 "
               << "e22_over_e35 e22_over_e37 e22_over_e53 "
               << "cluster_prob cluster_weta_cogx cluster_wphi_cogx "
               << "cluster_et1 cluster_et2 cluster_et3 cluster_et4 "
               << "cluster_w32 cluster_w52 cluster_w72 "
               << "recoisoET "
               << "is_tight pid " << std::endl;

                */

                // if (tight || nontight)
                {
                    // save to output file
                    outputfile << cluster_Et[icluster] << " "
                               << cluster_Eta[icluster] << " "
                               << cluster_Phi[icluster] << " "
                               << *vertexz << " "
                               << e11_over_e33 << " "
                               << e32_over_e35 << " "
                               << e11_over_e22 << " "
                               << e11_over_e13 << " "
                               << e11_over_e15 << " "
                               << e11_over_e17 << " "
                               << e11_over_e31 << " "
                               << e11_over_e51 << " "
                               << e11_over_e71 << " "
                               << e22_over_e33 << " "
                               << e22_over_e35 << " "
                               << e22_over_e37 << " "
                               << e22_over_e53 << " "
                               << cluster_prob[icluster] << " "
                               << cluster_weta_cogx[icluster] << " "
                               << cluster_wphi_cogx[icluster] << " "
                               << cluster_et1[icluster] << " "
                               << cluster_et2[icluster] << " "
                               << cluster_et3[icluster] << " "
                               << cluster_et4[icluster] << " "
                               << cluster_w32[icluster] << " "
                               << cluster_w52[icluster] << " "
                               << cluster_w72[icluster] << " "
                               << recoisoET << " "
                               << (tight ? 1 : 0) << " "
                               << pid << std::endl;
                }
            }
            else if (!issim && npb_cut_on)
            {
                // DATA NPB MODE:
                // Write out NPB-tagged clusters even if they fail `common_pass`.
                // `tight` is only evaluated inside the `common_pass` block, so here we label is_tight=0.
                outputfile << cluster_Et[icluster] << " "
                           << cluster_Eta[icluster] << " "
                           << cluster_Phi[icluster] << " "
                           << *vertexz << " "
                           << e11_over_e33 << " "
                           << e32_over_e35 << " "
                           << e11_over_e22 << " "
                           << e11_over_e13 << " "
                           << e11_over_e15 << " "
                           << e11_over_e17 << " "
                           << e11_over_e31 << " "
                           << e11_over_e51 << " "
                           << e11_over_e71 << " "
                           << e22_over_e33 << " "
                           << e22_over_e35 << " "
                           << e22_over_e37 << " "
                           << e22_over_e53 << " "
                           << cluster_prob[icluster] << " "
                           << cluster_weta_cogx[icluster] << " "
                           << cluster_wphi_cogx[icluster] << " "
                           << cluster_et1[icluster] << " "
                           << cluster_et2[icluster] << " "
                           << cluster_et3[icluster] << " "
                           << cluster_et4[icluster] << " "
                           << cluster_w32[icluster] << " "
                           << cluster_w52[icluster] << " "
                           << cluster_w72[icluster] << " "
                           << recoisoET << " "
                           << 0 << " "
                           << pid << std::endl;
            }
            else if (issim && write_outside_common_pass_mc)
            {
                // MC mode option:
                // Write out clusters even if they fail `common_pass` (preselection-free output).
                // `tight` is only evaluated inside the `common_pass` block, so here we label is_tight=0.
                outputfile << cluster_Et[icluster] << " "
                           << cluster_Eta[icluster] << " "
                           << cluster_Phi[icluster] << " "
                           << *vertexz << " "
                           << e11_over_e33 << " "
                           << e32_over_e35 << " "
                           << e11_over_e22 << " "
                           << e11_over_e13 << " "
                           << e11_over_e15 << " "
                           << e11_over_e17 << " "
                           << e11_over_e31 << " "
                           << e11_over_e51 << " "
                           << e11_over_e71 << " "
                           << e22_over_e33 << " "
                           << e22_over_e35 << " "
                           << e22_over_e37 << " "
                           << e22_over_e53 << " "
                           << cluster_prob[icluster] << " "
                           << cluster_weta_cogx[icluster] << " "
                           << cluster_wphi_cogx[icluster] << " "
                           << cluster_et1[icluster] << " "
                           << cluster_et2[icluster] << " "
                           << cluster_et3[icluster] << " "
                           << cluster_et4[icluster] << " "
                           << cluster_w32[icluster] << " "
                           << cluster_w52[icluster] << " "
                           << cluster_w72[icluster] << " "
                           << recoisoET << " "
                           << 0 << " "
                           << pid << std::endl;
            }
        }
        ientry++;
    }
}
