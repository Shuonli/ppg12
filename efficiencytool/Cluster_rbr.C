#include <yaml-cpp/yaml.h>

void Cluster_rbr(const std::string &configname = "config_nom468.yaml", const std::string filetype = "data")
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
    std::string lumifilename = "/sphenix/user/jocl/projects/LumiList/list_RN_18_26_30.list";

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
            if (!(iss >> runnumber >>bla>>bla>> lumi >> lumi_nc>>bla>>bla))
            {
                continue;
            }
            lumivals[runnumber] = lumi;
            lumivals_nc[runnumber] = lumi_nc;
            std::cout << "runnumber: " << runnumber << " lumi: " << lumi << std::endl;
        }
    }
    TFile *ftreein = new TFile(infilename.c_str(), "READ");
    TTree *slimtree = (TTree *)ftreein->Get(configYaml["input"]["tree"].as<std::string>().c_str());

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

    int mbdnorthhit, mbdsouthhit;
    int pythiaid, nparticles;
    int ncluster;
    int runnumber;
    Bool_t scaledtrigger[32] = {0};
    Bool_t livetrigger[32] = {0};

    float energy_scale;

    float vertexz;
    float vertexz_truth;

    float trigger_prescale[32] = {0};

    static const int nparticlesmax = 100;
    static const int nclustercontainermx = 50;

    float particle_E[nparticlesmax], particle_Pt[nparticlesmax], particle_Eta[nparticlesmax], particle_Phi[nparticlesmax], particle_truth_iso_02[nparticlesmax], particle_truth_iso_03[nparticlesmax], particle_truth_iso_04[nparticlesmax];
    int particle_pid[nparticlesmax], particle_trkid[nparticlesmax], particle_photonclass[nparticlesmax], particle_converted[nparticlesmax];
    float cluster_E[nclustercontainermx], cluster_Et[nclustercontainermx], cluster_Eta[nclustercontainermx], cluster_Phi[nclustercontainermx], cluster_prob[nclustercontainermx], cluster_iso_02[nclustercontainermx], cluster_iso_03[nclustercontainermx], cluster_iso_04[nclustercontainermx], cluster_e1[nclustercontainermx], cluster_e2[nclustercontainermx], cluster_e3[nclustercontainermx], cluster_e4[nclustercontainermx], cluster_et1[nclustercontainermx], cluster_et2[nclustercontainermx], cluster_et3[nclustercontainermx], cluster_et4[nclustercontainermx], cluster_weta[nclustercontainermx], cluster_wphi[nclustercontainermx], cluster_ietacent[nclustercontainermx], cluster_iphicent[nclustercontainermx], cluster_e11[nclustercontainermx], cluster_e22[nclustercontainermx], cluster_e13[nclustercontainermx], cluster_e15[nclustercontainermx], cluster_e17[nclustercontainermx], cluster_e31[nclustercontainermx], cluster_e51[nclustercontainermx], cluster_e71[nclustercontainermx], cluster_e33[nclustercontainermx], cluster_e35[nclustercontainermx], cluster_e37[nclustercontainermx], cluster_e53[nclustercontainermx], cluster_e73[nclustercontainermx], cluster_e55[nclustercontainermx], cluster_e57[nclustercontainermx], cluster_e75[nclustercontainermx], cluster_e77[nclustercontainermx], cluster_w32[nclustercontainermx], cluster_e32[nclustercontainermx], cluster_w72[nclustercontainermx], cluster_e72[nclustercontainermx], cluster_ihcal_et[nclustercontainermx], cluster_ohcal_et[nclustercontainermx], cluster_ihcal_et22[nclustercontainermx], cluster_ohcal_et22[nclustercontainermx], cluster_ihcal_et33[nclustercontainermx], cluster_ohcal_et33[nclustercontainermx];
    float cluster_w52[nclustercontainermx];
    int cluster_truthtrkID[nclustercontainermx], cluster_pid[nclustercontainermx];

    int cluster_detamax[nclustercontainermx], cluster_dphimax[nclustercontainermx], cluster_ihcal_ieta[nclustercontainermx], cluster_ihcal_iphi[nclustercontainermx], cluster_ohcal_ieta[nclustercontainermx], cluster_ohcal_iphi[nclustercontainermx];
    float cluster_CNN_prob[nclustercontainermx];

    float cluster_weta_cogx[nclustercontainermx], cluster_wphi_cogx[nclustercontainermx];

    static const int arraysize = 49;

    int cluster_ownership_array[nclustercontainermx][arraysize] = {0};

    float cluster_time_array[nclustercontainermx][arraysize] = {0};

    float cluster_e_array[nclustercontainermx][arraysize] = {0};

    float cluster_adc_array[nclustercontainermx][arraysize] = {0};

    int cluster_e_array_idx[nclustercontainermx][arraysize] = {0};

    int cluster_status_array[nclustercontainermx][arraysize] = {0};

    float cluster_iso_03_emcal[nclustercontainermx], cluster_iso_03_hcalin[nclustercontainermx], cluster_iso_03_hcalout[nclustercontainermx];
    float cluster_iso_03_60_emcal[nclustercontainermx], cluster_iso_03_60_hcalin[nclustercontainermx], cluster_iso_03_60_hcalout[nclustercontainermx];
    float cluster_iso_03_120_emcal[nclustercontainermx], cluster_iso_03_120_hcalin[nclustercontainermx], cluster_iso_03_120_hcalout[nclustercontainermx];
    int cluster_nsaturated[nclustercontainermx];

    static const int njettruthmax = 100;
    int njet_truth;
    float jet_truth_E[njettruthmax], jet_truth_Pt[njettruthmax], jet_truth_Eta[njettruthmax], jet_truth_Phi[njettruthmax];

    int njet;

    static const int njetmax = 100;

    float jet_E[njetmax], jet_Pt[njetmax], jet_Eta[njetmax], jet_Phi[njetmax];

    slimtree->SetBranchStatus("*", 0);

    slimtree->SetBranchStatus("mbdnorthhit", 1);
    slimtree->SetBranchAddress("mbdnorthhit", &mbdnorthhit);
    slimtree->SetBranchStatus("mbdsouthhit", 1);
    slimtree->SetBranchAddress("mbdsouthhit", &mbdsouthhit);
    slimtree->SetBranchStatus("vertexz", 1);
    slimtree->SetBranchAddress("vertexz", &vertexz);
    slimtree->SetBranchStatus("energy_scale", 1);
    slimtree->SetBranchAddress("energy_scale", &energy_scale);
    slimtree->SetBranchStatus("scaledtrigger", 1);
    slimtree->SetBranchAddress("scaledtrigger", scaledtrigger);
    slimtree->SetBranchStatus("livetrigger", 1);
    slimtree->SetBranchAddress("livetrigger", livetrigger);
    slimtree->SetBranchStatus("runnumber", 1);
    slimtree->SetBranchAddress("runnumber", &runnumber);
    slimtree->SetBranchStatus(Form("ncluster_%s", clusternodename.c_str()));
    slimtree->SetBranchAddress(Form("ncluster_%s", clusternodename.c_str()), &ncluster);

    slimtree->SetBranchStatus(Form("cluster_E_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_E_%s", clusternodename.c_str()), &cluster_E);
    slimtree->SetBranchStatus(Form("cluster_Et_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_Et_%s", clusternodename.c_str()), &cluster_Et);
    slimtree->SetBranchStatus(Form("cluster_Eta_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_Eta_%s", clusternodename.c_str()), &cluster_Eta);
    slimtree->SetBranchStatus(Form("cluster_Phi_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_Phi_%s", clusternodename.c_str()), &cluster_Phi);
    slimtree->SetBranchStatus(Form("cluster_prob_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_prob_%s", clusternodename.c_str()), &cluster_prob);
    slimtree->SetBranchStatus(Form("cluster_CNN_prob_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_CNN_prob_%s", clusternodename.c_str()), &cluster_CNN_prob);
    slimtree->SetBranchStatus(Form("cluster_pid_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_pid_%s", clusternodename.c_str()), &cluster_pid);
    slimtree->SetBranchStatus(Form("cluster_iso_02_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_02_%s", clusternodename.c_str()), &cluster_iso_02);
    slimtree->SetBranchStatus(Form("cluster_iso_03_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_%s", clusternodename.c_str()), &cluster_iso_03);
    slimtree->SetBranchStatus(Form("cluster_iso_04_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_04_%s", clusternodename.c_str()), &cluster_iso_04);
    slimtree->SetBranchStatus(Form("cluster_et1_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_et1_%s", clusternodename.c_str()), &cluster_et1);
    slimtree->SetBranchStatus(Form("cluster_et2_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_et2_%s", clusternodename.c_str()), &cluster_et2);
    slimtree->SetBranchStatus(Form("cluster_et3_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_et3_%s", clusternodename.c_str()), &cluster_et3);
    slimtree->SetBranchStatus(Form("cluster_et4_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_et4_%s", clusternodename.c_str()), &cluster_et4);
    slimtree->SetBranchStatus(Form("cluster_weta_cogx_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_weta_cogx_%s", clusternodename.c_str()), &cluster_weta_cogx);
    slimtree->SetBranchStatus(Form("cluster_wphi_cogx_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_wphi_cogx_%s", clusternodename.c_str()), &cluster_wphi_cogx);
    slimtree->SetBranchStatus(Form("cluster_nsaturated_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_nsaturated_%s", clusternodename.c_str()), &cluster_nsaturated);
    slimtree->SetBranchStatus(Form("cluster_e11_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e11_%s", clusternodename.c_str()), &cluster_e11);
    slimtree->SetBranchStatus(Form("cluster_e22_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e22_%s", clusternodename.c_str()), &cluster_e22);
    slimtree->SetBranchStatus(Form("cluster_e13_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e13_%s", clusternodename.c_str()), &cluster_e13);
    slimtree->SetBranchStatus(Form("cluster_e15_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e15_%s", clusternodename.c_str()), &cluster_e15);
    slimtree->SetBranchStatus(Form("cluster_e17_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e17_%s", clusternodename.c_str()), &cluster_e17);
    slimtree->SetBranchStatus(Form("cluster_e31_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e31_%s", clusternodename.c_str()), &cluster_e31);
    slimtree->SetBranchStatus(Form("cluster_e51_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e51_%s", clusternodename.c_str()), &cluster_e51);
    slimtree->SetBranchStatus(Form("cluster_e71_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e71_%s", clusternodename.c_str()), &cluster_e71);
    slimtree->SetBranchStatus(Form("cluster_e33_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e33_%s", clusternodename.c_str()), &cluster_e33);
    slimtree->SetBranchStatus(Form("cluster_e35_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e35_%s", clusternodename.c_str()), &cluster_e35);
    slimtree->SetBranchStatus(Form("cluster_e37_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e37_%s", clusternodename.c_str()), &cluster_e37);
    slimtree->SetBranchStatus(Form("cluster_e53_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e53_%s", clusternodename.c_str()), &cluster_e53);
    slimtree->SetBranchStatus(Form("cluster_e73_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e73_%s", clusternodename.c_str()), &cluster_e73);
    slimtree->SetBranchStatus(Form("cluster_e55_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e55_%s", clusternodename.c_str()), &cluster_e55);
    slimtree->SetBranchStatus(Form("cluster_e57_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e57_%s", clusternodename.c_str()), &cluster_e57);
    slimtree->SetBranchStatus(Form("cluster_e75_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e75_%s", clusternodename.c_str()), &cluster_e75);
    slimtree->SetBranchStatus(Form("cluster_e77_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e77_%s", clusternodename.c_str()), &cluster_e77);
    slimtree->SetBranchStatus(Form("cluster_w32_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_w32_%s", clusternodename.c_str()), &cluster_w32);
    slimtree->SetBranchStatus(Form("cluster_e32_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e32_%s", clusternodename.c_str()), &cluster_e32);
    slimtree->SetBranchStatus(Form("cluster_w72_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_w72_%s", clusternodename.c_str()), &cluster_w72);
    slimtree->SetBranchStatus(Form("cluster_e72_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e72_%s", clusternodename.c_str()), &cluster_e72);
    slimtree->SetBranchStatus(Form("cluster_w52_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_w52_%s", clusternodename.c_str()), &cluster_w52);
    slimtree->SetBranchStatus(Form("cluster_iso_03_emcal_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_emcal_%s", clusternodename.c_str()), &cluster_iso_03_emcal);
    slimtree->SetBranchStatus(Form("cluster_iso_03_hcalin_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_hcalin_%s", clusternodename.c_str()), &cluster_iso_03_hcalin);
    slimtree->SetBranchStatus(Form("cluster_iso_03_hcalout_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_hcalout_%s", clusternodename.c_str()), &cluster_iso_03_hcalout);
    slimtree->SetBranchStatus(Form("cluster_iso_03_60_emcal_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_60_emcal_%s", clusternodename.c_str()), &cluster_iso_03_60_emcal);
    slimtree->SetBranchStatus(Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()), &cluster_iso_03_60_hcalin);
    slimtree->SetBranchStatus(Form("cluster_iso_03_60_hcalout_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_60_hcalout_%s", clusternodename.c_str()), &cluster_iso_03_60_hcalout);
    slimtree->SetBranchStatus(Form("cluster_iso_03_120_emcal_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_120_emcal_%s", clusternodename.c_str()), &cluster_iso_03_120_emcal);
    slimtree->SetBranchStatus(Form("cluster_iso_03_120_hcalin_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_120_hcalin_%s", clusternodename.c_str()), &cluster_iso_03_120_hcalin);
    slimtree->SetBranchStatus(Form("cluster_iso_03_120_hcalout_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_120_hcalout_%s", clusternodename.c_str()), &cluster_iso_03_120_hcalout);

    std::set<int> skiprunnumbers = {47698, 51489, 51721, 51725, 53284};
    int nentries = slimtree->GetEntries();
    for (int ientry = 0; ientry < nentries; ientry++)
    {

	 
        if (ientry == 16152886)
            continue;
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        slimtree->GetEntry(ientry);
        if (!issim)
        {
            if (skiprunnumbers.find(runnumber) != skiprunnumbers.end())
            {
                continue;
            }

            if (scaledtrigger[trigger_used] == 0)
                continue;
        }

        // check vertex cut
        if (abs(vertexz) > vertexcut)
            continue;

        if (!(mbdnorthhit >= 1 && mbdsouthhit >= 1))
            continue;

        for (int icluster = 0; icluster < ncluster; icluster++)
        {
            // hard code the cut here for now :(
            if (cluster_Et[icluster] < 7)
                continue;
            if (cluster_Et[icluster] > 26)
                continue;
            if (abs(cluster_Eta[icluster]) > 0.7)
                continue;

            float e11_over_e33 = cluster_e11[icluster] / cluster_e33[icluster];

            float e32_over_e35 = cluster_e32[icluster] / cluster_e35[icluster];

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
                common_clusters[runnumber] += 1;
                if (tight && iso)
                {
                    tight_iso_clusters[runnumber] += 1;
                }
                if (tight && noniso)
                {
                    tight_noniso_clusters[runnumber] += 1;
                }
                if (nontight && iso)
                {
                    nontight_iso_clusters[runnumber] += 1;
                }
                if (nontight && noniso)
                {
                    nontight_noniso_clusters[runnumber] += 1;
                }
            }
        }
    }
    // print out the results
    std::vector<double> x_common, x_common_error, y_common, y_common_error;
    std::vector<double> x_tight_iso, x_tight_iso_error, y_tight_iso, y_tight_iso_error;
    std::vector<double> x_tight_noniso, x_tight_noniso_error, y_tight_noniso, y_tight_noniso_error;
    std::vector<double> x_nontight_iso, x_nontight_iso_error, y_nontight_iso, y_nontight_iso_error;
    std::vector<double> x_nontight_noniso, x_nontight_noniso_error, y_nontight_noniso, y_nontight_noniso_error;
    std::vector<double> x_corr, y_corr, x_corr_error, y_corr_error;
    float total_lumi = 0;
    for (auto const &runnumber : common_clusters)
    {
        float lumi = lumivals[runnumber.first];
        std::cout << runnumber.first << " " << lumi << std::endl;


        if ((runnumber.second < 10) || lumi == 0)
        {
            // set to 0
            common_clusters[runnumber.first] = 0;
            tight_iso_clusters[runnumber.first] = 0;
            tight_noniso_clusters[runnumber.first] = 0;
            nontight_iso_clusters[runnumber.first] = 0;
            nontight_noniso_clusters[runnumber.first] = 0;
        }
        else
        {
            total_lumi += lumi;
            // scale every thing by luminosity
            float common_clusters_error = sqrt(common_clusters[runnumber.first]);
            common_clusters[runnumber.first] = common_clusters[runnumber.first] / lumi;
            common_clusters_error = common_clusters_error / lumi;
            x_common.push_back(runnumber.first);
            x_common_error.push_back(0);
            y_common.push_back(common_clusters[runnumber.first]);
            y_common_error.push_back(common_clusters_error);

            float tight_iso_clusters_error = sqrt(tight_iso_clusters[runnumber.first]);
            tight_iso_clusters[runnumber.first] = tight_iso_clusters[runnumber.first] / lumi;
            tight_iso_clusters_error = tight_iso_clusters_error / lumi;
            x_tight_iso.push_back(runnumber.first);
            x_tight_iso_error.push_back(0);
            y_tight_iso.push_back(tight_iso_clusters[runnumber.first]);
            y_tight_iso_error.push_back(tight_iso_clusters_error);

            float tight_noniso_clusters_error = sqrt(tight_noniso_clusters[runnumber.first]);
            tight_noniso_clusters[runnumber.first] = tight_noniso_clusters[runnumber.first] / lumi;
            tight_noniso_clusters_error = tight_noniso_clusters_error / lumi;
            x_tight_noniso.push_back(runnumber.first);
            x_tight_noniso_error.push_back(0);
            y_tight_noniso.push_back(tight_noniso_clusters[runnumber.first]);
            y_tight_noniso_error.push_back(tight_noniso_clusters_error);

            float nontight_iso_clusters_error = sqrt(nontight_iso_clusters[runnumber.first]);
            nontight_iso_clusters[runnumber.first] = nontight_iso_clusters[runnumber.first] / lumi;
            nontight_iso_clusters_error = nontight_iso_clusters_error / lumi;
            x_nontight_iso.push_back(runnumber.first);
            x_nontight_iso_error.push_back(0);
            y_nontight_iso.push_back(nontight_iso_clusters[runnumber.first]);
            y_nontight_iso_error.push_back(nontight_iso_clusters_error);

            float nontight_noniso_clusters_error = sqrt(nontight_noniso_clusters[runnumber.first]);
            nontight_noniso_clusters[runnumber.first] = nontight_noniso_clusters[runnumber.first] / lumi;
            nontight_noniso_clusters_error = nontight_noniso_clusters_error / lumi;
            x_nontight_noniso.push_back(runnumber.first);
            x_nontight_noniso_error.push_back(0);
            y_nontight_noniso.push_back(nontight_noniso_clusters[runnumber.first]);
            y_nontight_noniso_error.push_back(nontight_noniso_clusters_error);

            float correction = lumivals[runnumber.first] / lumivals_nc[runnumber.first];
            x_corr.push_back(runnumber.first);
            x_corr_error.push_back(0);
            y_corr.push_back(correction);
            y_corr_error.push_back(0);
        }
    }

    std::cout << "Total lumi: " << total_lumi << std::endl;

    TFile* fout = new TFile("rbrQA.root", "RECREATE");

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

    gr_common->Draw("AP");

    fout->cd();
    gr_common->Write();
    gr_tight_iso->Write();
    gr_tight_noniso->Write();
    gr_nontight_iso->Write();
    gr_nontight_noniso->Write();
    gr_corr->Write();

    fout->Write();
    fout->Close();

}
