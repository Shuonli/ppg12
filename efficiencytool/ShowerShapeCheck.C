#include <yaml-cpp/yaml.h>

void ShowerShapeCheck(const std::string &configname = "config.yaml", const std::string filetype = "data", bool doinclusive = false)
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    bool issim = true;

    bool isbackground = false;


    if (filetype == "data")
    {
        issim = false;
    }

    std::string incusive_str = doinclusive ? "_inclusive" : "";

    std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();

    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();

    std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;

    if (!issim)
    {
        infilename = configYaml["input"]["data_file"].as<std::string>();
    }

    std::cout << "infilename: " << infilename << std::endl;

    float bg_timing_cut = -0.5;
    float bg_timing_cut_more = -0.4;
    float eta_bg_lower = 0.0;

    float max_photon_lower = 0;
    float max_photon_upper = 100;

    // unit in pb
    const float photon5cross = 2.017e+08 * 0.000442571;
    const float photon10cross = 3.690e+07 * 0.000181474;
    const float photon20cross = 1.571e+05 * 0.000673448;

    // Hanpu uses unit in b
    const float jet10cross = 3.646e-6;
    const float jet15cross = 36864930.0 * 0.011059973 * 1e-12;
    const float jet20cross = 1392140.9 * 0.042 * 1e-12;
    const float jet30cross = 2.505e-9;

    float max_jet_lower = 0;
    float max_jet_upper = 100;
    float weight = 1.0;

    float energy_scale_lower = 0;
    float energy_scale_upper = 100;

    float cluster_ET_upper = 100;

    if (filetype == "photon5")
    {
        max_photon_lower = 5;
        max_photon_upper = 14;
        weight = photon5cross / photon20cross;
    }
    else if (filetype == "photon10")
    {
        max_photon_lower = 14;
        max_photon_upper = 30;
        weight = photon10cross / photon20cross;
    }
    else if (filetype == "photon20")
    {
        max_photon_lower = 30;
        max_photon_upper = 100;
        weight = 1.0;
    }
    else if (filetype == "jet10")
    {
        max_jet_lower = 10;
        max_jet_upper = 17;
        energy_scale_lower = 10;
        energy_scale_upper = 16;
        cluster_ET_upper = 17;
        weight = jet10cross / jet30cross;
        isbackground = true;
    }
    else if (filetype == "jet15")
    {
        max_jet_lower = 17;
        max_jet_upper = 23;
        energy_scale_lower = 16;
        energy_scale_upper = 20;
        cluster_ET_upper = 23;
        weight = jet15cross / jet30cross;
        isbackground = true;
    }
    else if (filetype == "jet20")
    {
        max_jet_lower = 23;
        max_jet_upper = 30;
        energy_scale_lower = 20;
        energy_scale_upper = 30;
        cluster_ET_upper = 30;
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
    TFile *ftreein = new TFile(infilename.c_str(), "READ");
    std::string var_type = configYaml["output"]["var_type"].as<std::string>();

    std::string outfilename = configYaml["output"]["eff_outfile"].as<std::string>() + "shower_shape" + "_" + filetype + incusive_str + ".root";

    std::string treestring = configYaml["output"]["eff_outfile"].as<std::string>() + "_tree" + "_" + filetype + incusive_str;

    std::string responsefilename = "bla";

    if (!issim)
    {
        outfilename = configYaml["output"]["data_outfile"].as<std::string>() + "shower_shape" + "_" + ".root";
        // unfolding is only for sim
        responsefilename = "bla";
    }

    std::cout << "outfilename: " << outfilename << std::endl;
    std::cout << "responsefilename: " << responsefilename << std::endl;

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

    float common_e32_over_e35_max = configYaml["analysis"]["common"]["e32_over_e35_max"].as<float>(1.0);
    float common_e32_over_e35_min = configYaml["analysis"]["common"]["e32_over_e35_min"].as<float>(0.8);

    float common_et1_min = configYaml["analysis"]["common"]["et1_min"].as<float>(0.6);
    float common_et1_max = configYaml["analysis"]["common"]["et1_max"].as<float>(1.0);

    float common_wr_cogx_bound = configYaml["analysis"]["common"]["wr_cogx_bound"].as<float>();
    float common_cluster_weta_cogx_bound = configYaml["analysis"]["common"]["cluster_weta_cogx_bound"].as<float>();

    int reweight = configYaml["analysis"]["unfold"]["reweight"].as<int>(); // 0 for no reweighting, 1 for reweighting

    float clusterescale = configYaml["analysis"]["cluster_escale"].as<float>(1.0);

    int mbdnorthhit, mbdsouthhit;
    int pythiaid, nparticles;
    int ncluster;
    Bool_t scaledtrigger[32] = {0};
    Bool_t livetrigger[32] = {0};

    float energy_scale;

    float vertexz;
    float vertexz_truth;

    static const int nparticlesmax = 1000;
    static const int nclustercontainermx = 100;

    static const int arraysize = 49;

    float particle_E[nparticlesmax], particle_Pt[nparticlesmax], particle_Eta[nparticlesmax], particle_Phi[nparticlesmax], particle_truth_iso_02[nparticlesmax], particle_truth_iso_03[nparticlesmax], particle_truth_iso_04[nparticlesmax];
    int particle_pid[nparticlesmax], particle_trkid[nparticlesmax], particle_photonclass[nparticlesmax], particle_converted[nparticlesmax];
    float cluster_E[nclustercontainermx], cluster_Et[nclustercontainermx], cluster_Eta[nclustercontainermx], cluster_Phi[nclustercontainermx], cluster_prob[nclustercontainermx], cluster_iso_02[nclustercontainermx], cluster_iso_03[nclustercontainermx], cluster_iso_04[nclustercontainermx], cluster_e1[nclustercontainermx], cluster_e2[nclustercontainermx], cluster_e3[nclustercontainermx], cluster_e4[nclustercontainermx], cluster_et1[nclustercontainermx], cluster_et2[nclustercontainermx], cluster_et3[nclustercontainermx], cluster_et4[nclustercontainermx], cluster_weta[nclustercontainermx], cluster_wphi[nclustercontainermx], cluster_ietacent[nclustercontainermx], cluster_iphicent[nclustercontainermx], cluster_e11[nclustercontainermx], cluster_e22[nclustercontainermx], cluster_e13[nclustercontainermx], cluster_e15[nclustercontainermx], cluster_e17[nclustercontainermx], cluster_e31[nclustercontainermx], cluster_e51[nclustercontainermx], cluster_e71[nclustercontainermx], cluster_e33[nclustercontainermx], cluster_e35[nclustercontainermx], cluster_e37[nclustercontainermx], cluster_e53[nclustercontainermx], cluster_e73[nclustercontainermx], cluster_e55[nclustercontainermx], cluster_e57[nclustercontainermx], cluster_e75[nclustercontainermx], cluster_e77[nclustercontainermx], cluster_w32[nclustercontainermx], cluster_e32[nclustercontainermx], cluster_w72[nclustercontainermx], cluster_e72[nclustercontainermx], cluster_ihcal_et[nclustercontainermx], cluster_ohcal_et[nclustercontainermx], cluster_ihcal_et22[nclustercontainermx], cluster_ohcal_et22[nclustercontainermx], cluster_ihcal_et33[nclustercontainermx], cluster_ohcal_et33[nclustercontainermx];
    float cluster_w52[nclustercontainermx];
    int cluster_truthtrkID[nclustercontainermx], cluster_pid[nclustercontainermx];

    int cluster_detamax[nclustercontainermx], cluster_dphimax[nclustercontainermx], cluster_ihcal_ieta[nclustercontainermx], cluster_ihcal_iphi[nclustercontainermx], cluster_ohcal_ieta[nclustercontainermx], cluster_ohcal_iphi[nclustercontainermx];
    float cluster_CNN_prob[nclustercontainermx];

    float cluster_weta_cogx[nclustercontainermx], cluster_wphi_cogx[nclustercontainermx], cluster_weta_cog[nclustercontainermx], cluster_wphi_cog[nclustercontainermx];

    int cluster_ownership_array[nclustercontainermx][arraysize] = {0};

    float cluster_time_array[nclustercontainermx][arraysize] = {0};

    float cluster_e_array[nclustercontainermx][arraysize] = {0};

    float cluster_adc_array[nclustercontainermx][arraysize] = {0};

    int cluster_e_array_idx[nclustercontainermx][arraysize] = {0};

    int cluster_status_array[nclustercontainermx][arraysize] = {0};

    float cluster_iso_03_emcal[nclustercontainermx], cluster_iso_03_hcalin[nclustercontainermx], cluster_iso_03_hcalout[nclustercontainermx];
    float cluster_iso_03_60_emcal[nclustercontainermx], cluster_iso_03_60_hcalin[nclustercontainermx], cluster_iso_03_60_hcalout[nclustercontainermx];
    float cluster_iso_03_120_emcal[nclustercontainermx], cluster_iso_03_120_hcalin[nclustercontainermx], cluster_iso_03_120_hcalout[nclustercontainermx];

    static const int njettruthmax = 100;
    int njet_truth;
    float jet_truth_E[njettruthmax], jet_truth_Pt[njettruthmax], jet_truth_Eta[njettruthmax], jet_truth_Phi[njettruthmax];

    int njet;

    static const int njetmax = 100;

    float jet_E[njetmax], jet_Pt[njetmax], jet_Eta[njetmax], jet_Phi[njetmax];

    slimtree->SetBranchAddress("mbdnorthhit", &mbdnorthhit);
    slimtree->SetBranchAddress("mbdsouthhit", &mbdsouthhit);
    slimtree->SetBranchAddress("pythiaid", &pythiaid);
    slimtree->SetBranchAddress("nparticles", &nparticles);
    slimtree->SetBranchAddress(Form("ncluster_%s", clusternodename.c_str()), &ncluster);
    slimtree->SetBranchAddress("vertexz", &vertexz);
    slimtree->SetBranchAddress("vertexz_truth", &vertexz_truth);
    slimtree->SetBranchAddress("energy_scale", &energy_scale);
    slimtree->SetBranchAddress("scaledtrigger", scaledtrigger);
    slimtree->SetBranchAddress("livetrigger", livetrigger);

    slimtree->SetBranchAddress("particle_E", &particle_E);
    slimtree->SetBranchAddress("particle_Pt", &particle_Pt);
    slimtree->SetBranchAddress("particle_Eta", &particle_Eta);
    slimtree->SetBranchAddress("particle_Phi", &particle_Phi);
    slimtree->SetBranchAddress("particle_truth_iso_02", &particle_truth_iso_02);
    slimtree->SetBranchAddress("particle_truth_iso_03", &particle_truth_iso_03);
    slimtree->SetBranchAddress("particle_truth_iso_04", &particle_truth_iso_04);
    slimtree->SetBranchAddress("particle_pid", &particle_pid);
    slimtree->SetBranchAddress("particle_trkid", &particle_trkid);
    slimtree->SetBranchAddress("particle_photonclass", &particle_photonclass);
    slimtree->SetBranchAddress("particle_converted", &particle_converted);

    slimtree->SetBranchAddress(Form("cluster_E_%s", clusternodename.c_str()), &cluster_E);
    slimtree->SetBranchAddress(Form("cluster_Et_%s", clusternodename.c_str()), &cluster_Et);
    slimtree->SetBranchAddress(Form("cluster_Eta_%s", clusternodename.c_str()), &cluster_Eta);
    slimtree->SetBranchAddress(Form("cluster_Phi_%s", clusternodename.c_str()), &cluster_Phi);
    slimtree->SetBranchAddress(Form("cluster_prob_%s", clusternodename.c_str()), &cluster_prob);
    slimtree->SetBranchAddress(Form("cluster_CNN_prob_%s", clusternodename.c_str()), &cluster_CNN_prob);
    slimtree->SetBranchAddress(Form("cluster_truthtrkID_%s", clusternodename.c_str()), &cluster_truthtrkID);
    slimtree->SetBranchAddress(Form("cluster_pid_%s", clusternodename.c_str()), &cluster_pid);
    slimtree->SetBranchAddress(Form("cluster_iso_02_%s", clusternodename.c_str()), &cluster_iso_02);
    slimtree->SetBranchAddress(Form("cluster_iso_03_%s", clusternodename.c_str()), &cluster_iso_03);
    slimtree->SetBranchAddress(Form("cluster_iso_04_%s", clusternodename.c_str()), &cluster_iso_04);
    slimtree->SetBranchAddress(Form("cluster_e1_%s", clusternodename.c_str()), &cluster_e1);
    slimtree->SetBranchAddress(Form("cluster_e2_%s", clusternodename.c_str()), &cluster_e2);
    slimtree->SetBranchAddress(Form("cluster_e3_%s", clusternodename.c_str()), &cluster_e3);
    slimtree->SetBranchAddress(Form("cluster_e4_%s", clusternodename.c_str()), &cluster_e4);
    slimtree->SetBranchAddress(Form("cluster_et1_%s", clusternodename.c_str()), &cluster_et1);
    slimtree->SetBranchAddress(Form("cluster_et2_%s", clusternodename.c_str()), &cluster_et2);
    slimtree->SetBranchAddress(Form("cluster_et3_%s", clusternodename.c_str()), &cluster_et3);
    slimtree->SetBranchAddress(Form("cluster_et4_%s", clusternodename.c_str()), &cluster_et4);
    slimtree->SetBranchAddress(Form("cluster_weta_%s", clusternodename.c_str()), &cluster_weta);
    slimtree->SetBranchAddress(Form("cluster_wphi_%s", clusternodename.c_str()), &cluster_wphi);
    slimtree->SetBranchAddress(Form("cluster_ietacent_%s", clusternodename.c_str()), &cluster_ietacent);
    slimtree->SetBranchAddress(Form("cluster_iphicent_%s", clusternodename.c_str()), &cluster_iphicent);
    slimtree->SetBranchAddress(Form("cluster_detamax_%s", clusternodename.c_str()), &cluster_detamax);
    slimtree->SetBranchAddress(Form("cluster_dphimax_%s", clusternodename.c_str()), &cluster_dphimax);
    slimtree->SetBranchAddress(Form("cluster_weta_cogx_%s", clusternodename.c_str()), &cluster_weta_cogx);
    slimtree->SetBranchAddress(Form("cluster_wphi_cogx_%s", clusternodename.c_str()), &cluster_wphi_cogx);
    slimtree->SetBranchAddress(Form("cluster_weta_cog_%s", clusternodename.c_str()), &cluster_weta_cog);
    slimtree->SetBranchAddress(Form("cluster_wphi_cog_%s", clusternodename.c_str()), &cluster_wphi_cog);

    slimtree->SetBranchAddress(Form("cluster_e11_%s", clusternodename.c_str()), &cluster_e11);
    slimtree->SetBranchAddress(Form("cluster_e22_%s", clusternodename.c_str()), &cluster_e22);
    slimtree->SetBranchAddress(Form("cluster_e13_%s", clusternodename.c_str()), &cluster_e13);
    slimtree->SetBranchAddress(Form("cluster_e15_%s", clusternodename.c_str()), &cluster_e15);
    slimtree->SetBranchAddress(Form("cluster_e17_%s", clusternodename.c_str()), &cluster_e17);
    slimtree->SetBranchAddress(Form("cluster_e31_%s", clusternodename.c_str()), &cluster_e31);
    slimtree->SetBranchAddress(Form("cluster_e51_%s", clusternodename.c_str()), &cluster_e51);
    slimtree->SetBranchAddress(Form("cluster_e71_%s", clusternodename.c_str()), &cluster_e71);
    slimtree->SetBranchAddress(Form("cluster_e33_%s", clusternodename.c_str()), &cluster_e33);
    slimtree->SetBranchAddress(Form("cluster_e35_%s", clusternodename.c_str()), &cluster_e35);
    slimtree->SetBranchAddress(Form("cluster_e37_%s", clusternodename.c_str()), &cluster_e37);
    slimtree->SetBranchAddress(Form("cluster_e53_%s", clusternodename.c_str()), &cluster_e53);
    slimtree->SetBranchAddress(Form("cluster_e73_%s", clusternodename.c_str()), &cluster_e73);
    slimtree->SetBranchAddress(Form("cluster_e55_%s", clusternodename.c_str()), &cluster_e55);
    slimtree->SetBranchAddress(Form("cluster_e57_%s", clusternodename.c_str()), &cluster_e57);
    slimtree->SetBranchAddress(Form("cluster_e75_%s", clusternodename.c_str()), &cluster_e75);
    slimtree->SetBranchAddress(Form("cluster_e77_%s", clusternodename.c_str()), &cluster_e77);
    slimtree->SetBranchAddress(Form("cluster_w32_%s", clusternodename.c_str()), &cluster_w32);
    slimtree->SetBranchAddress(Form("cluster_e32_%s", clusternodename.c_str()), &cluster_e32);
    slimtree->SetBranchAddress(Form("cluster_w72_%s", clusternodename.c_str()), &cluster_w72);
    slimtree->SetBranchAddress(Form("cluster_e72_%s", clusternodename.c_str()), &cluster_e72);
    slimtree->SetBranchAddress(Form("cluster_w52_%s", clusternodename.c_str()), &cluster_w52);

    slimtree->SetBranchAddress(Form("cluster_e_array_%s", clusternodename.c_str()), cluster_e_array);
    // slimtree->SetBranchAddress(Form("cluster_adc_array_%s", clusternodename.c_str()), &cluster_adc_array);
    // slimtree->SetBranchAddress(Form("cluster_e_array_idx_%s", clusternodename.c_str()), &cluster_e_array_idx);
    // slimtree->SetBranchAddress(Form("cluster_status_array_%s", clusternodename.c_str()), &cluster_status_array);
    slimtree->SetBranchAddress(Form("cluster_time_array_%s", clusternodename.c_str()), &cluster_time_array);
    slimtree->SetBranchAddress(Form("cluster_ownership_array_%s", clusternodename.c_str()), &cluster_ownership_array);

    slimtree->SetBranchAddress(Form("cluster_iso_03_emcal_%s", clusternodename.c_str()), &cluster_iso_03_emcal);
    slimtree->SetBranchAddress(Form("cluster_iso_03_hcalin_%s", clusternodename.c_str()), &cluster_iso_03_hcalin);
    slimtree->SetBranchAddress(Form("cluster_iso_03_hcalout_%s", clusternodename.c_str()), &cluster_iso_03_hcalout);
    slimtree->SetBranchAddress(Form("cluster_iso_03_60_emcal_%s", clusternodename.c_str()), &cluster_iso_03_60_emcal);
    slimtree->SetBranchAddress(Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()), &cluster_iso_03_60_hcalin);
    slimtree->SetBranchAddress(Form("cluster_iso_03_60_hcalout_%s", clusternodename.c_str()), &cluster_iso_03_60_hcalout);
    slimtree->SetBranchAddress(Form("cluster_iso_03_120_emcal_%s", clusternodename.c_str()), &cluster_iso_03_120_emcal);
    slimtree->SetBranchAddress(Form("cluster_iso_03_120_hcalin_%s", clusternodename.c_str()), &cluster_iso_03_120_hcalin);
    slimtree->SetBranchAddress(Form("cluster_iso_03_120_hcalout_%s", clusternodename.c_str()), &cluster_iso_03_120_hcalout);

    slimtree->SetBranchAddress(Form("cluster_ihcal_et_%s", clusternodename.c_str()), &cluster_ihcal_et);
    slimtree->SetBranchAddress(Form("cluster_ohcal_et_%s", clusternodename.c_str()), &cluster_ohcal_et);
    slimtree->SetBranchAddress(Form("cluster_ihcal_et22_%s", clusternodename.c_str()), &cluster_ihcal_et22);
    slimtree->SetBranchAddress(Form("cluster_ohcal_et22_%s", clusternodename.c_str()), &cluster_ohcal_et22);
    slimtree->SetBranchAddress(Form("cluster_ihcal_et33_%s", clusternodename.c_str()), &cluster_ihcal_et33);
    slimtree->SetBranchAddress(Form("cluster_ohcal_et33_%s", clusternodename.c_str()), &cluster_ohcal_et33);
    slimtree->SetBranchAddress(Form("cluster_ihcal_ieta_%s", clusternodename.c_str()), &cluster_ihcal_ieta);
    slimtree->SetBranchAddress(Form("cluster_ihcal_iphi_%s", clusternodename.c_str()), &cluster_ihcal_iphi);
    slimtree->SetBranchAddress(Form("cluster_ohcal_ieta_%s", clusternodename.c_str()), &cluster_ohcal_ieta);
    slimtree->SetBranchAddress(Form("cluster_ohcal_iphi_%s", clusternodename.c_str()), &cluster_ohcal_iphi);

    slimtree->SetBranchAddress("njet_truth", &njet_truth);
    slimtree->SetBranchAddress("jet_truth_E", &jet_truth_E);
    slimtree->SetBranchAddress("jet_truth_Pt", &jet_truth_Pt);
    slimtree->SetBranchAddress("jet_truth_Eta", &jet_truth_Eta);
    slimtree->SetBranchAddress("jet_truth_Phi", &jet_truth_Phi);

    slimtree->SetBranchAddress("njet", &njet);
    slimtree->SetBranchAddress("jet_E", &jet_E);
    slimtree->SetBranchAddress("jet_Pt", &jet_Pt);
    slimtree->SetBranchAddress("jet_Eta", &jet_Eta);
    slimtree->SetBranchAddress("jet_Phi", &jet_Phi);

    // ttree for TMVA training
    const static int nETbins = 1;
    double ETbins[nETbins + 1] = {5, 100};

    TFile *foutTMVA[nETbins];
    TTree *treeTMVA[nETbins];

    float prob_TMVA;
    float CNN_prob_TMVA;
    float e17_to_e77_TMVA;
    float e37_to_e77_TMVA;
    float e32_to_e35_TMVA;
    float e33_to_e35_TMVA;
    float e11_to_e33_TMVA;
    float e11_to_E_TMVA;
    float e33_to_E_TMVA;
    float hcalet33_to_ettot_TMVA;
    float ihcalet33_to_ettot_TMVA;
    float ohcalet33_to_ettot_TMVA;
    float hcalet22_to_ettot_TMVA;
    float ihcalet22_to_ettot_TMVA;
    float ohcalet22_to_ettot_TMVA;
    float detamax_TMVA;
    float dphimax_TMVA;
    float e1_TMVA;
    float e2_TMVA;
    float e3_TMVA;
    float e4_TMVA;
    float et1_TMVA;
    float et2_TMVA;
    float et3_TMVA;
    float et4_TMVA;
    float weta_TMVA;
    float wphi_TMVA;
    float w32_TMVA;
    float w52_TMVA;
    float w72_TMVA;
    float wr_TMVA;
    float wrr_TMVA;
    float weta_cog_TMVA;
    float wphi_cog_TMVA;
    float weta_cogx_TMVA;
    float wphi_cogx_TMVA;
    float vertex_z_TMVA;
    float cluster_eta_TMVA;
    float cluster_ET_TMVA;
    float cluster_iso_02_TMVA;
    float cluster_iso_03_TMVA;
    float cluster_iso_04_TMVA;
    float cluster_iso_03_emcal_TMVA;
    float cluster_iso_03_hcalin_TMVA;
    float cluster_iso_03_hcalout_TMVA;
    float cluster_iso_03_60_emcal_TMVA;
    float cluster_iso_03_60_hcalin_TMVA;
    float cluster_iso_03_60_hcalout_TMVA;
    float cluster_iso_03_120_emcal_TMVA;
    float cluster_iso_03_120_hcalin_TMVA;
    float cluster_iso_03_120_hcalout_TMVA;
    int label_TMVA;

    for (int iET = 0; iET < nETbins; iET++)
    {

        foutTMVA[iET] = new TFile(Form("%s_%d.root", treestring.c_str(), iET), "RECREATE");
        treeTMVA[iET] = new TTree("tree", "tree");

        treeTMVA[iET]->Branch("prob", &prob_TMVA, "prob/F");
        treeTMVA[iET]->Branch("CNN_prob", &CNN_prob_TMVA, "CNN_prob/F");
        treeTMVA[iET]->Branch("e17_to_e77", &e17_to_e77_TMVA, "e17_to_e77/F");
        treeTMVA[iET]->Branch("e37_to_e77", &e37_to_e77_TMVA, "e37_to_e77/F");
        treeTMVA[iET]->Branch("e32_to_e35", &e32_to_e35_TMVA, "e32_to_e35/F");
        treeTMVA[iET]->Branch("e33_to_e35", &e33_to_e35_TMVA, "e33_to_e35/F");
        treeTMVA[iET]->Branch("e11_to_e33", &e11_to_e33_TMVA, "e11_to_e33/F");
        treeTMVA[iET]->Branch("e11_to_E", &e11_to_E_TMVA, "e11_to_E/F");
        treeTMVA[iET]->Branch("e33_to_E", &e33_to_E_TMVA, "e33_to_E/F");
        treeTMVA[iET]->Branch("hcalet33_to_ettot", &hcalet33_to_ettot_TMVA, "hcalet33_to_ettot/F");
        treeTMVA[iET]->Branch("ihcalet33_to_ettot", &ihcalet33_to_ettot_TMVA, "ihcalet33_to_ettot/F");
        treeTMVA[iET]->Branch("ohcalet33_to_ettot", &ohcalet33_to_ettot_TMVA, "ohcalet33_to_ettot/F");
        treeTMVA[iET]->Branch("hcalet22_to_ettot", &hcalet22_to_ettot_TMVA, "hcalet22_to_ettot/F");
        treeTMVA[iET]->Branch("ihcalet22_to_ettot", &ihcalet22_to_ettot_TMVA, "ihcalet22_to_ettot/F");
        treeTMVA[iET]->Branch("ohcalet22_to_ettot", &ohcalet22_to_ettot_TMVA, "ohcalet22_to_ettot/F");
        treeTMVA[iET]->Branch("detamax", &detamax_TMVA, "detamax/F");
        treeTMVA[iET]->Branch("dphimax", &dphimax_TMVA, "dphimax/F");
        treeTMVA[iET]->Branch("e1", &e1_TMVA, "e1/F");
        treeTMVA[iET]->Branch("e2", &e2_TMVA, "e2/F");
        treeTMVA[iET]->Branch("e3", &e3_TMVA, "e3/F");
        treeTMVA[iET]->Branch("e4", &e4_TMVA, "e4/F");
        treeTMVA[iET]->Branch("et1", &et1_TMVA, "et1/F");
        treeTMVA[iET]->Branch("et2", &et2_TMVA, "et2/F");
        treeTMVA[iET]->Branch("et3", &et3_TMVA, "et3/F");
        treeTMVA[iET]->Branch("et4", &et4_TMVA, "et4/F");
        treeTMVA[iET]->Branch("weta", &weta_TMVA, "weta/F");
        treeTMVA[iET]->Branch("wphi", &wphi_TMVA, "wphi/F");
        treeTMVA[iET]->Branch("w32", &w32_TMVA, "w32/F");
        treeTMVA[iET]->Branch("w52", &w52_TMVA, "w52/F");
        treeTMVA[iET]->Branch("w72", &w72_TMVA, "w72/F");
        treeTMVA[iET]->Branch("wr", &wr_TMVA, "wr/F");
        treeTMVA[iET]->Branch("wrr", &wrr_TMVA, "wrr/F");
        treeTMVA[iET]->Branch("weta_cog", &weta_cog_TMVA, "weta_cog/F");
        treeTMVA[iET]->Branch("wphi_cog", &wphi_cog_TMVA, "wphi_cog/F");
        treeTMVA[iET]->Branch("weta_cogx", &weta_cogx_TMVA, "weta_cogx/F");
        treeTMVA[iET]->Branch("wphi_cogx", &wphi_cogx_TMVA, "wphi_cogx/F");
        treeTMVA[iET]->Branch("vertex_z", &vertex_z_TMVA, "vertex_z/F");
        treeTMVA[iET]->Branch("cluster_eta", &cluster_eta_TMVA, "cluster_eta/F");
        treeTMVA[iET]->Branch("cluster_Et", &cluster_ET_TMVA, "cluster_ET/F");
        treeTMVA[iET]->Branch("cluster_iso_02", &cluster_iso_02_TMVA, "cluster_iso_02/F");
        treeTMVA[iET]->Branch("cluster_iso_03", &cluster_iso_03_TMVA, "cluster_iso_03/F");
        treeTMVA[iET]->Branch("cluster_iso_04", &cluster_iso_04_TMVA, "cluster_iso_04/F");
        treeTMVA[iET]->Branch("cluster_iso_03_emcal", &cluster_iso_03_emcal_TMVA, "cluster_iso_03_emcal/F");
        treeTMVA[iET]->Branch("cluster_iso_03_hcalin", &cluster_iso_03_hcalin_TMVA, "cluster_iso_03_hcalin/F");
        treeTMVA[iET]->Branch("cluster_iso_03_hcalout", &cluster_iso_03_hcalout_TMVA, "cluster_iso_03_hcalout/F");
        treeTMVA[iET]->Branch("cluster_iso_03_60_emcal", &cluster_iso_03_60_emcal_TMVA, "cluster_iso_03_60_emcal/F");
        treeTMVA[iET]->Branch("cluster_iso_03_60_hcalin", &cluster_iso_03_60_hcalin_TMVA, "cluster_iso_03_60_hcalin/F");
        treeTMVA[iET]->Branch("cluster_iso_03_60_hcalout", &cluster_iso_03_60_hcalout_TMVA, "cluster_iso_03_60_hcalout/F");
        treeTMVA[iET]->Branch("cluster_iso_03_120_emcal", &cluster_iso_03_120_emcal_TMVA, "cluster_iso_03_120_emcal/F");
        treeTMVA[iET]->Branch("cluster_iso_03_120_hcalin", &cluster_iso_03_120_hcalin_TMVA, "cluster_iso_03_120_hcalin/F");
        treeTMVA[iET]->Branch("cluster_iso_03_120_hcalout", &cluster_iso_03_120_hcalout_TMVA, "cluster_iso_03_120_hcalout/F");
        treeTMVA[iET]->Branch("label", &label_TMVA, "label/I");
    }

    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    // TH2D *
    std::vector<std::vector<TH1D *>> h_tight_reco_iso;
    h_tight_reco_iso.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH1D *>> h_nt_reco_iso;
    h_nt_reco_iso.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH1D *>> h_signal_truth_iso;
    h_signal_truth_iso.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH1D *>> h_background_truth_iso;
    h_background_truth_iso.resize(eta_bins.size() - 1);
    std::vector<TH2D *> h_ET_isoET;

    for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
    {
        h_ET_isoET.push_back(new TH2D(Form("h_ET_isoET_eta%d", ieta),
                                      Form("ET vs isoET %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                      400, 0, 50, 4400, -5, 50));

        for (int ipt = 0; ipt < n_pT_bins; ipt++)
        {
            h_signal_truth_iso[ieta].push_back(new TH1D(Form("h_signal_truth_iso_eta%d_pt%d", ieta, ipt),
                                                        Form("Signal truth iso %.1f < eta < %.1f, %.1f < pT < %.1f", eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                        400, -10, 30));

            h_background_truth_iso[ieta].push_back(new TH1D(Form("h_background_truth_iso_eta%d_pt%d", ieta, ipt),
                                                            Form("Background truth iso %.1f < eta < %.1f, %.1f < pT < %.1f", eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                            400, -10, 30));

            h_tight_reco_iso[ieta].push_back(new TH1D(Form("h_tight_reco_iso_eta%d_pt%d", ieta, ipt),
                                                      Form("Tight Iso ET %.1f < eta < %.1f, %.1f < pT < %.1f", eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                      400, -10, 30));

            h_nt_reco_iso[ieta].push_back(new TH1D(Form("h_nt_reco_iso_eta%d_pt%d", ieta, ipt),
                                                   Form("Non Tight Iso ET %.1f < eta < %.1f, %.1f < pT < %.1f", eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                   400, -10, 30));
        }
    }

    // shower shape correlation with isoET before prelim cut
    // a master vector of 3D histograms
    std::map<std::string, std::vector<std::vector<std::vector<TH2F *>>>> h2d_all;

    // 2) Define your "histogram names" (the old map keys become string keys).
    //    For instance, you can directly populate them:

    static const int ncuts = 7;

    h2d_all["h2d_prob"];
    h2d_all["h2d_CNN_prob"];
    h2d_all["h2d_e17_to_e77"];
    h2d_all["h2d_e37_to_e77"];
    h2d_all["h2d_e32_to_e35"];
    h2d_all["h2d_e33_to_e35"];
    h2d_all["h2d_e11_to_e33"];
    h2d_all["h2d_e11_to_E"];
    h2d_all["h2d_e33_to_E"];
    h2d_all["h2d_hcalet33_to_ettot"];
    h2d_all["h2d_ihcalet33_to_ettot"];
    h2d_all["h2d_ohcalet33_to_ettot"];
    h2d_all["h2d_hcalet22_to_ettot"];
    h2d_all["h2d_ihcalet22_to_ettot"];
    h2d_all["h2d_ohcalet22_to_ettot"];
    h2d_all["h2d_detamax"];
    h2d_all["h2d_dphimax"];
    h2d_all["h2d_e1"];
    h2d_all["h2d_e2"];
    h2d_all["h2d_e3"];
    h2d_all["h2d_e4"];
    h2d_all["h2d_et1"];
    h2d_all["h2d_et2"];
    h2d_all["h2d_et3"];
    h2d_all["h2d_et4"];
    h2d_all["h2d_weta"];
    h2d_all["h2d_wphi"];
    h2d_all["h2d_w32"];
    h2d_all["h2d_w52"];
    h2d_all["h2d_w72"];
    h2d_all["h2d_wr"];
    h2d_all["h2d_wrr"];
    h2d_all["h2d_weta_cog"];
    h2d_all["h2d_wphi_cog"];
    h2d_all["h2d_weta_cogx"];
    h2d_all["h2d_wphi_cogx"];

    std::vector<std::vector<TH2D *>> h_ET_eta;
    h_ET_eta.resize(ncuts);

    for (int icut = 0; icut < ncuts; icut++)
    {
        for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
        {
            h_ET_eta[icut].push_back(new TH2D(Form("h_ET_eta_cut%d_eta%d", icut, ieta),
                                              Form("ET vs eta cut %d %.1f < eta < %.1f", icut, eta_bins[ieta], eta_bins[ieta + 1]),
                                              400, 0, 50, 100, -1, 1));
        }
    }

    for (auto &kv : h2d_all)
    {
        // kv is a reference to std::pair<const std::string, std::vector<std::vector<std::vector<TH2F*>>>>
        const std::string &basename = kv.first; // the key (string)
        auto &histVector3D = kv.second;         // the 3D vector of TH2F*, which we can now modify

        for (int icut = 0; icut < ncuts; icut++)
        {
            std::vector<std::vector<TH2F *>> h2d_eta;
            for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
            {
                //
                std::vector<TH2F *> h2d_pt;
                for (int ipt = 0; ipt < n_pT_bins; ipt++)
                {
                    TH2F *h2D = new TH2F(
                        Form("%s_eta%d_pt%d_cut%d", basename.c_str(), ieta, ipt, icut),
                        Form("Tight Iso ET %.1f < eta < %.1f, %.1f < pT < %.1f",
                             eta_bins[ieta],
                             eta_bins[ieta + 1],
                             pT_bin_edges[ipt],
                             pT_bin_edges[ipt + 1]),
                        300, -1, 2,  // X bins
                        200, -10, 30 // Y bins
                    );
                    h2D->GetXaxis()->SetTitle(basename.c_str());
                    h2D->GetYaxis()->SetTitle("Iso ET [GeV]");
                    h2d_pt.push_back(h2D);
                }
                h2d_eta.push_back(h2d_pt);
            }

            // Now we can push_back into histVector3D (the map value).
            histVector3D.push_back(h2d_eta);
        }
    }

    int nentries = slimtree->GetEntries();
    for (int ientry = 0; ientry < nentries; ientry++)
    {
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        slimtree->GetEntry(ientry);
        if (!issim)
        {

            if (scaledtrigger[trigger_used] == 0)
                continue;
        }

        std::set<int> signal_set;
        std::set<int> background_set;
        std::map<int, int> particle_trkidmap;
        if (issim)
        {
            float maxphotonpT = 0;
            int maxphotonclass = 0;

            for (int iparticle = 0; iparticle < nparticles; iparticle++)
            {
                int bg_particle = false;
                particle_trkidmap[particle_trkid[iparticle]] = iparticle;
                float truthisoET = 0;
                if (conesize == 4)
                {
                    truthisoET = particle_truth_iso_04[iparticle];
                }
                else if (conesize == 3)
                {
                    truthisoET = particle_truth_iso_03[iparticle];
                }
                else if (conesize == 2)
                {
                    truthisoET = particle_truth_iso_02[iparticle];
                }
                else
                {
                    std::cout << "Error: conesize not supported" << std::endl;
                    continue;
                }

                // find eta and pT bins
                float particle_eta = particle_Eta[iparticle];
                int etabin = -1;
                for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
                {
                    if (particle_eta > eta_bins[ieta] && particle_eta < eta_bins[ieta + 1])
                    {
                        etabin = ieta;
                        break;
                    }
                }

                float pTbin = -1;
                float particlePT = particle_Pt[iparticle];
                for (int ipt = 0; ipt < n_pT_bins; ipt++)
                {
                    if (particlePT > pT_bins[ipt] && particlePT < pT_bins[ipt + 1])
                    {
                        pTbin = ipt;
                        break;
                    }
                }

                if (particle_pid[iparticle] == 22)
                {
                    if (particle_Pt[iparticle] > maxphotonpT)
                    {
                        maxphotonpT = particle_Pt[iparticle];
                        maxphotonclass = particle_photonclass[iparticle];
                    }

                    if (particle_photonclass[iparticle] < 3) // direct or fragmentation
                    {
                        if (truthisoET < truthisocut)
                        {
                            signal_set.insert(iparticle);
                            if (etabin != -1 && pTbin != -1)
                            {
                                h_signal_truth_iso[etabin][pTbin]->Fill(truthisoET, weight);
                            }
                        }
                    }
                    else
                    {
                        background_set.insert(iparticle);
                        bg_particle = true;
                    }
                    if (bg_particle)
                    {
                        if (etabin != -1 && pTbin != -1)
                        {
                            h_background_truth_iso[etabin][pTbin]->Fill(truthisoET, weight);
                        }
                    }
                }
                else
                {
                    background_set.insert(iparticle);
                    bg_particle = true;
                }
            }
            if (!isbackground)
            {

                if ((maxphotonpT > max_photon_upper) || (maxphotonpT < max_photon_lower))
                {
                    continue;
                }
            }
            // loop over truth_jets for the background sample
            if (isbackground)
            {
                float maxjetpT = 0;
                for (int ijet = 0; ijet < njet_truth; ijet++)
                {
                    if (jet_truth_Pt[ijet] > maxjetpT)
                    {
                        maxjetpT = jet_truth_Pt[ijet];
                    }
                }

                if ((maxjetpT > max_jet_upper) || (maxjetpT < max_jet_lower))
                {
                    continue;
                }

                /*
                 if ((energy_scale > energy_scale_upper) || (energy_scale < energy_scale_lower))
                 {
                     continue;
                 }
                 */
            }
        }

        std::vector<float> jetphi;

        for (int ijet = 0; ijet < njet; ijet++)
        {
            if (jet_Pt[ijet] < 10)
                continue;
            jetphi.push_back(jet_Phi[ijet]);
        }
        // loop over clusters
        for (int icluster = 0; icluster < ncluster; icluster++)
        {
            // need ET > 10 GeV
            if (cluster_Et[icluster] < reco_min_ET)
                continue;

            // for jet 10 event we want to remove some high ET clusters to reduce the fluctuation
            if (isbackground)
            {
                if (cluster_Et[icluster] > cluster_ET_upper)
                {
                    continue;
                }
            }

            float rhad22 = (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]) / (cluster_Et[icluster] + (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]));
            float rhad33 = (cluster_ihcal_et33[icluster] + cluster_ohcal_et33[icluster]) / (cluster_Et[icluster] + (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]));

            float reta77 = cluster_e37[icluster] / cluster_e77[icluster];
            float rphi77 = cluster_e73[icluster] / cluster_e77[icluster];

            float reta55 = cluster_e35[icluster] / cluster_e55[icluster];
            float rphi55 = cluster_e53[icluster] / cluster_e55[icluster];

            float reta = cluster_e33[icluster] / cluster_e73[icluster];

            float rphi = cluster_e33[icluster] / cluster_e37[icluster];

            float re11_E = cluster_e11[icluster] / cluster_E[icluster];

            float wr_cogx = cluster_wphi_cogx[icluster] / cluster_weta_cogx[icluster];

            float hcalet33 = cluster_ihcal_et33[icluster] + cluster_ohcal_et33[icluster];
            float hcalet22 = cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster];

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

            //
            if (issim)
            {
                if (particle_trkidmap.find(cluster_truthtrkID[icluster]) == particle_trkidmap.end())
                {
                    // std::cout<<"trackid: "<<cluster_truthtrkID[icluster]<<std::endl;
                    // std::cout << "Error: cluster_truthtrkID not found in particle_trkidmap" << std::endl;
                    continue;
                }
                int iparticle = particle_trkidmap[cluster_truthtrkID[icluster]];
                if (!isbackground)
                {
                    if (signal_set.find(iparticle) == signal_set.end())
                    {
                        continue;
                    }
                }
                else
                {
                    if (!doinclusive)
                    {
                        if (background_set.find(iparticle) == background_set.end())
                        {
                            continue;
                        }
                    }
                }
            }

            prob_TMVA = cluster_prob[icluster];
            CNN_prob_TMVA = cluster_CNN_prob[icluster];
            e17_to_e77_TMVA = cluster_e17[icluster] / cluster_e77[icluster];
            e37_to_e77_TMVA = cluster_e37[icluster] / cluster_e77[icluster];
            e32_to_e35_TMVA = cluster_e32[icluster] / cluster_e35[icluster];
            e33_to_e35_TMVA = cluster_e33[icluster] / cluster_e35[icluster];
            e11_to_e33_TMVA = cluster_e11[icluster] / cluster_e33[icluster];
            e11_to_E_TMVA = cluster_e11[icluster] / cluster_E[icluster];
            e33_to_E_TMVA = cluster_e33[icluster] / cluster_E[icluster];
            hcalet33_to_ettot_TMVA = hcalet33 / (cluster_Et[icluster] + hcalet33);
            ihcalet33_to_ettot_TMVA = cluster_ihcal_et33[icluster] / (cluster_Et[icluster] + hcalet33);
            ohcalet33_to_ettot_TMVA = cluster_ohcal_et33[icluster] / (cluster_Et[icluster] + hcalet33);
            hcalet22_to_ettot_TMVA = hcalet22 / (cluster_Et[icluster] + hcalet22);
            ihcalet22_to_ettot_TMVA = cluster_ihcal_et22[icluster] / (cluster_Et[icluster] + hcalet22);
            ohcalet22_to_ettot_TMVA = cluster_ohcal_et22[icluster] / (cluster_Et[icluster] + hcalet22);
            detamax_TMVA = cluster_detamax[icluster] / 10.0;
            dphimax_TMVA = cluster_dphimax[icluster] / 10.0;
            e1_TMVA = cluster_e1[icluster];
            e2_TMVA = cluster_e2[icluster];
            e3_TMVA = cluster_e3[icluster];
            e4_TMVA = cluster_e4[icluster];
            et1_TMVA = cluster_et1[icluster];
            et2_TMVA = cluster_et2[icluster];
            et3_TMVA = cluster_et3[icluster];
            et4_TMVA = cluster_et4[icluster];
            weta_TMVA = cluster_weta[icluster];
            wphi_TMVA = cluster_wphi[icluster];
            w32_TMVA = cluster_w32[icluster];
            w52_TMVA = cluster_w52[icluster];
            w72_TMVA = cluster_w72[icluster];
            wr_TMVA = cluster_wphi[icluster] / cluster_weta[icluster];
            wrr_TMVA = cluster_weta[icluster] / cluster_wphi[icluster];
            weta_cog_TMVA = cluster_weta_cog[icluster];
            wphi_cog_TMVA = cluster_wphi_cog[icluster];
            weta_cogx_TMVA = cluster_weta_cogx[icluster];
            wphi_cogx_TMVA = cluster_wphi_cogx[icluster];
            vertex_z_TMVA = vertexz;
            cluster_eta_TMVA = cluster_Eta[icluster];
            cluster_ET_TMVA = cluster_Et[icluster];
            cluster_iso_02_TMVA = cluster_iso_02[icluster];
            cluster_iso_03_TMVA = cluster_iso_03[icluster];
            cluster_iso_04_TMVA = cluster_iso_04[icluster];
            cluster_iso_03_emcal_TMVA = cluster_iso_03_emcal[icluster];
            cluster_iso_03_hcalin_TMVA = cluster_iso_03_hcalin[icluster];
            cluster_iso_03_hcalout_TMVA = cluster_iso_03_hcalout[icluster];
            cluster_iso_03_60_emcal_TMVA = cluster_iso_03_60_emcal[icluster];
            cluster_iso_03_60_hcalin_TMVA = cluster_iso_03_60_hcalin[icluster];
            cluster_iso_03_60_hcalout_TMVA = cluster_iso_03_60_hcalout[icluster];
            cluster_iso_03_120_emcal_TMVA = cluster_iso_03_120_emcal[icluster];
            cluster_iso_03_120_hcalin_TMVA = cluster_iso_03_120_hcalin[icluster];
            cluster_iso_03_120_hcalout_TMVA = cluster_iso_03_120_hcalout[icluster];
            label_TMVA = isbackground ? 0 : 1;

            int ETbin = -1;
            for (int iET = 0; iET < nETbins; iET++)
            {
                if (cluster_Et[icluster] > ETbins[iET] && cluster_Et[icluster] < ETbins[iET + 1])
                {
                    ETbin = iET;
                    break;
                }
            }
            if (ETbin != -1 && issim)
            {
                treeTMVA[ETbin]->Fill();
            }

            if (abs(vertexz) > vertexcut)
                continue;

            if (!(mbdnorthhit >= 1 && mbdsouthhit >= 1))
                continue;

            float cluster_eta = cluster_Eta[icluster];
            int etabin = -1;
            for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
            {
                if (cluster_eta > eta_bins[ieta] && cluster_eta < eta_bins[ieta + 1])
                {
                    etabin = ieta;
                    break;
                }
            }
            if (etabin == -1)
            {
                continue;
            }

            float pTbin = -1;
            float clusterET = cluster_Et[icluster];
            for (int ipt = 0; ipt < n_pT_bins; ipt++)
            {
                if (clusterET > pT_bins[ipt] && clusterET < pT_bins[ipt + 1])
                {
                    pTbin = ipt;
                    break;
                }
            }
            if (pTbin == -1)
            {
                continue;
            }

            // fill iso ET and et histogram
            h_ET_isoET[etabin]->Fill(cluster_Et[icluster], recoisoET, weight);

            // now fill the histograms
            auto fillAllHists = [&](int idx, size_t icl)
            {
                h2d_all["h2d_prob"][idx][etabin][pTbin]->Fill(cluster_prob[icl], recoisoET, weight);
                h2d_all["h2d_CNN_prob"][idx][etabin][pTbin]->Fill(cluster_CNN_prob[icl], recoisoET, weight);
                h2d_all["h2d_e17_to_e77"][idx][etabin][pTbin]->Fill(cluster_e17[icl] / cluster_e77[icl], recoisoET, weight);
                h2d_all["h2d_e37_to_e77"][idx][etabin][pTbin]->Fill(cluster_e37[icl] / cluster_e77[icl], recoisoET, weight);
                h2d_all["h2d_e32_to_e35"][idx][etabin][pTbin]->Fill(cluster_e32[icl] / cluster_e35[icl], recoisoET, weight);
                h2d_all["h2d_e33_to_e35"][idx][etabin][pTbin]->Fill(cluster_e33[icl] / cluster_e35[icl], recoisoET, weight);
                h2d_all["h2d_e11_to_e33"][idx][etabin][pTbin]->Fill(cluster_e11[icl] / cluster_e33[icl], recoisoET, weight);
                h2d_all["h2d_e11_to_E"][idx][etabin][pTbin]->Fill(cluster_e11[icl] / cluster_E[icl], recoisoET, weight);
                h2d_all["h2d_e33_to_E"][idx][etabin][pTbin]->Fill(cluster_e33[icl] / cluster_E[icl], recoisoET, weight);
                h2d_all["h2d_hcalet33_to_ettot"][idx][etabin][pTbin]->Fill(hcalet33 / (cluster_Et[icl] + hcalet33), recoisoET, weight);
                h2d_all["h2d_ihcalet33_to_ettot"][idx][etabin][pTbin]->Fill(cluster_ihcal_et33[icl] / (cluster_Et[icl] + hcalet33), recoisoET, weight);
                h2d_all["h2d_ohcalet33_to_ettot"][idx][etabin][pTbin]->Fill(cluster_ohcal_et33[icl] / (cluster_Et[icl] + hcalet33), recoisoET, weight);
                h2d_all["h2d_hcalet22_to_ettot"][idx][etabin][pTbin]->Fill(hcalet22 / (cluster_Et[icl] + hcalet22), recoisoET, weight);
                h2d_all["h2d_ihcalet22_to_ettot"][idx][etabin][pTbin]->Fill(cluster_ihcal_et22[icl] / (cluster_Et[icl] + hcalet22), recoisoET, weight);
                h2d_all["h2d_ohcalet22_to_ettot"][idx][etabin][pTbin]->Fill(cluster_ohcal_et22[icl] / (cluster_Et[icl] + hcalet22), recoisoET, weight);
                h2d_all["h2d_detamax"][idx][etabin][pTbin]->Fill(cluster_detamax[icl] / 10.0, recoisoET, weight);
                h2d_all["h2d_dphimax"][idx][etabin][pTbin]->Fill(cluster_dphimax[icl] / 10.0, recoisoET, weight);
                h2d_all["h2d_e1"][idx][etabin][pTbin]->Fill(cluster_e1[icl], recoisoET, weight);
                h2d_all["h2d_e2"][idx][etabin][pTbin]->Fill(cluster_e2[icl], recoisoET, weight);
                h2d_all["h2d_e3"][idx][etabin][pTbin]->Fill(cluster_e3[icl], recoisoET, weight);
                h2d_all["h2d_e4"][idx][etabin][pTbin]->Fill(cluster_e4[icl], recoisoET, weight);
                h2d_all["h2d_et1"][idx][etabin][pTbin]->Fill(cluster_et1[icl], recoisoET, weight);
                h2d_all["h2d_et2"][idx][etabin][pTbin]->Fill(cluster_et2[icl], recoisoET, weight);
                h2d_all["h2d_et3"][idx][etabin][pTbin]->Fill(cluster_et3[icl], recoisoET, weight);
                h2d_all["h2d_et4"][idx][etabin][pTbin]->Fill(cluster_et4[icl], recoisoET, weight);
                h2d_all["h2d_weta"][idx][etabin][pTbin]->Fill(cluster_weta[icl], recoisoET, weight);
                h2d_all["h2d_wphi"][idx][etabin][pTbin]->Fill(cluster_wphi[icl], recoisoET, weight);
                h2d_all["h2d_w32"][idx][etabin][pTbin]->Fill(cluster_w32[icl], recoisoET, weight);
                h2d_all["h2d_w52"][idx][etabin][pTbin]->Fill(cluster_w52[icl], recoisoET, weight);
                h2d_all["h2d_w72"][idx][etabin][pTbin]->Fill(cluster_w72[icl], recoisoET, weight);
                h2d_all["h2d_wr"][idx][etabin][pTbin]->Fill(cluster_wphi[icl] / cluster_weta[icl], recoisoET, weight);
                h2d_all["h2d_wrr"][idx][etabin][pTbin]->Fill(cluster_weta[icl] / cluster_wphi[icl], recoisoET, weight);
                h2d_all["h2d_weta_cog"][idx][etabin][pTbin]->Fill(cluster_weta_cog[icl], recoisoET, weight);
                h2d_all["h2d_wphi_cog"][idx][etabin][pTbin]->Fill(cluster_wphi_cog[icl], recoisoET, weight);
                h2d_all["h2d_weta_cogx"][idx][etabin][pTbin]->Fill(cluster_weta_cogx[icl], recoisoET, weight);
                h2d_all["h2d_wphi_cogx"][idx][etabin][pTbin]->Fill(cluster_wphi_cogx[icl], recoisoET, weight);

                h_ET_eta[idx][etabin]->Fill(cluster_Et[icl], cluster_Eta[icl], weight);
            };

            fillAllHists(0, icluster);

            // calculate cluster time
            float clusteravgtime = 0;
            float cluster_total_e = 0;
            for (int i = 0; i < 49; i++)
            {

                if (cluster_ownership_array[icluster][i] == 1)
                {
                    clusteravgtime += cluster_time_array[icluster][i] * cluster_e_array[icluster][i];
                    // std::cout<<"cluster_time_array[icluster][i]: "<<cluster_time_array[icluster][i]<<std::endl;
                    cluster_total_e += cluster_e_array[icluster][i];
                }
            }
            clusteravgtime = cluster_total_e > 0 ? clusteravgtime / cluster_total_e : 0;

            bool badtime = clusteravgtime < bg_timing_cut;

            bool morebadtime = clusteravgtime < bg_timing_cut_more;

            bool otherside_jet = false;

            for (int ijet = 0; ijet < (int)jetphi.size(); ijet++)
            {
                float dphi = cluster_Phi[icluster] - jetphi[ijet];

                while (dphi > M_PI)
                    dphi = dphi - 2 * M_PI;
                while (dphi < -M_PI)
                    dphi = dphi + 2 * M_PI;

                if (abs(dphi) > (M_PI / 2))
                {
                    otherside_jet = true;
                    break;
                }
            }

            if (badtime && !otherside_jet)
            {
                fillAllHists(4, icluster);
            }

            if (morebadtime && !otherside_jet)
            {
                fillAllHists(5, icluster);
            }

            // streak event removal cut
            // if (wr_cogx < 0.4 && cluster_weta_cogx[icluster] > 1)
            //    continue;

            // common cuts
            bool common_pass = false;
            bool tight = false;
            bool nontight = false;
            bool pscut = false;
            if (cluster_prob[icluster] > common_prob_min &&
                cluster_prob[icluster] < common_prob_max &&
                e11_over_e33 > common_e11_over_e33_min &&
                e11_over_e33 < common_e11_over_e33_max &&
                //(!(wr_cogx < common_wr_cogx_bound && cluster_weta_cogx[icluster] > common_cluster_weta_cogx_bound))
                (cluster_weta_cogx[icluster] < 0.8))
            {
                common_pass = true;
                if(e32_over_e35 > common_e32_over_e35_min &&
                e32_over_e35 < common_e32_over_e35_max &&
                cluster_et1[icluster] > common_et1_min &&
                cluster_et1[icluster] < common_et1_max )
                {
                    pscut = true;
                }
            }
            if (!common_pass)
                continue;
            // fill the samething again with 1

            fillAllHists(1, icluster);

            if(pscut)
            {
                fillAllHists(6, icluster);
            }

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

            if (tight)
            {
                fillAllHists(2, icluster);
                h_tight_reco_iso[etabin][pTbin]->Fill(recoisoET, weight);
            }
            if (nontight)
            {
                fillAllHists(3, icluster);
                h_nt_reco_iso[etabin][pTbin]->Fill(recoisoET, weight);
            }
        }
    }
    fout->cd();
    fout->Write();
    fout->Close();

    for (int iET = 0; iET < nETbins; iET++)
    {
        // save the tree
        foutTMVA[iET]->cd();
        treeTMVA[iET]->Write();
        foutTMVA[iET]->Close();
    }
}
