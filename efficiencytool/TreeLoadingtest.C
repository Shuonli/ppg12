#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TF1.h>
#include <TSystem.h>
#include <TEfficiency.h>
#include <TObjString.h>
#include <TRandom3.h>
#include <yaml-cpp/yaml.h>
// unfolding
#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>

// R__LOAD_LIBRARY(/sphenix/user/egm2153/calib_study/JetValidation/analysis/roounfold/libRooUnfold.so)

void SaveYamlToRoot(TFile *f, const char *yaml_filename)
{
    // Read YAML file into a string
    std::ifstream yaml_file(yaml_filename);
    std::string yaml_content((std::istreambuf_iterator<char>(yaml_file)),
                             std::istreambuf_iterator<char>());

    // Create a ROOT file and save the YAML string
    TObjString yaml_obj(yaml_content.c_str());
    f->cd();
    yaml_obj.Write("config");
}

void TreeLoadingtest(const std::string &configname = "config_bdt_test.yaml", const std::string filetype = "jet10")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    bool issim = true;

    if (filetype == "data")
    {
        issim = false;
    }

    // bool isbackground = (bool)configYaml["input"]["isbackground"].as<int>();
    static const bool isbackground = filetype.find("jet")!=string::npos ? true : false;
    // hardcode here for now
    float bg_timing_cut = -0.5;
    float npb_weta_min = 0.4;

    std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();

    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();

    std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;

    if (!issim)
    {
        infilename = configYaml["input"]["data_file"].as<std::string>();
    }

    TFile *corrFile = new TFile("/sphenix/user/hanpuj/JES_MC_Calibration/offline/output_forshuhang.root", "READ");
    TF1 *f_corr = (TF1 *)corrFile->Get("f_corr_run21_r04_z0_eta0123");


    /*
 1  Constant     1.07368e+00   1.10770e-03  -6.32308e-06  -1.15696e-02
2  Mean        -1.73712e+00   8.31735e-02  -1.32021e-04  -1.25750e-03
3  Sigma        4.47289e+01   2.37499e-01   2.37499e-01  -5.45865e-02
*/
    // gaussian vertex reweight with a gaussian function
    TF1 *f_vertex_reweight = new TF1("f_vertex_reweight", "gaus", -50, 50);
    f_vertex_reweight->SetParameters(1.07368e+00, -1.73712e+00, 4.47289e+01);

    // std::string infilename = "/sphenix/tg/tg01/commissioning/CaloCalibWG/sli/ppg12/ana450/condorout/combine.root";
    TFile *ftreein = new TFile(infilename.c_str(), "READ");
    std::string var_type = configYaml["output"]["var_type"].as<std::string>();

    std::string outfilename = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + filetype + "_" + var_type + ".root";

    std::string responsefilename = configYaml["output"]["response_outfile"].as<std::string>() + "_" + filetype + "_" + var_type + ".root";

    if (!issim)
    {
        outfilename = configYaml["output"]["data_outfile"].as<std::string>() + "_" + var_type + ".root";
        // unfolding is only for sim
        responsefilename = "bla.root";
    }

    std::cout << "outfilename: " << outfilename << std::endl;
    std::cout << "responsefilename: " << responsefilename << std::endl;

    TTree *slimtree = (TTree *)ftreein->Get(configYaml["input"]["tree"].as<std::string>().c_str());

    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();

    std::string bdt_model_name = configYaml["input"]["bdt_model_name"].as<std::string>("base");



    // polynomial 3 for the reweighting
    // TF1 *f_reweight = new TF1("f_reweight", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0, 100);
    TF1 *f_reweight = new TF1("f_reweight", "([0] + [1]*x + [3]*x*x) / (1 + [2]*x + [4]*x*x)", 0, 100);
    // need to make this into the config file in the future!!!
    // f_reweight->SetParameters(1.04713, 0.00623875, -0.00106856, 2.64199e-06);
    // f_reweight->SetParameters(1.04713, 0.00623875, -0.00106856, 2.64199e-06);
    f_reweight->SetParameters(0.714962, -0.0856443, -0.125383, 0.00345831, 0.00462972);

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
    float cluster_bdt[nclustercontainermx];

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
    slimtree->SetBranchStatus("pythiaid", 1);
    slimtree->SetBranchAddress("pythiaid", &pythiaid);
    slimtree->SetBranchStatus("nparticles", 1);
    slimtree->SetBranchAddress("nparticles", &nparticles);
    slimtree->SetBranchStatus("vertexz", 1);
    slimtree->SetBranchAddress("vertexz", &vertexz);
    slimtree->SetBranchStatus("vertexz_truth", 1);
    slimtree->SetBranchAddress("vertexz_truth", &vertexz_truth);
    slimtree->SetBranchStatus("energy_scale", 1);
    slimtree->SetBranchAddress("energy_scale", &energy_scale);
    slimtree->SetBranchStatus("scaledtrigger", 1);
    slimtree->SetBranchAddress("scaledtrigger", scaledtrigger);
    slimtree->SetBranchStatus("livetrigger", 1);
    slimtree->SetBranchAddress("livetrigger", livetrigger);
    slimtree->SetBranchStatus("trigger_prescale", 1);
    slimtree->SetBranchAddress("trigger_prescale", trigger_prescale);
    slimtree->SetBranchStatus("runnumber", 1);
    slimtree->SetBranchAddress("runnumber", &runnumber);
    slimtree->SetBranchStatus(Form("ncluster_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("ncluster_%s", clusternodename.c_str()), &ncluster);

    slimtree->SetBranchStatus("particle_E", 1);
    slimtree->SetBranchAddress("particle_E", particle_E);
    slimtree->SetBranchStatus("particle_Pt", 1);
    slimtree->SetBranchAddress("particle_Pt", particle_Pt);
    slimtree->SetBranchStatus("particle_Eta", 1);
    slimtree->SetBranchAddress("particle_Eta", particle_Eta);
    slimtree->SetBranchStatus("particle_Phi", 1);
    slimtree->SetBranchAddress("particle_Phi", particle_Phi);
    slimtree->SetBranchStatus("particle_truth_iso_02", 1);
    slimtree->SetBranchAddress("particle_truth_iso_02", particle_truth_iso_02);
    slimtree->SetBranchStatus("particle_truth_iso_03", 1);
    slimtree->SetBranchAddress("particle_truth_iso_03", particle_truth_iso_03);
    slimtree->SetBranchStatus("particle_truth_iso_04", 1);
    slimtree->SetBranchAddress("particle_truth_iso_04", particle_truth_iso_04);
    slimtree->SetBranchStatus("particle_pid", 1);
    slimtree->SetBranchAddress("particle_pid", particle_pid);
    slimtree->SetBranchStatus("particle_trkid", 1);
    slimtree->SetBranchAddress("particle_trkid", particle_trkid);
    slimtree->SetBranchStatus("particle_photonclass", 1);
    slimtree->SetBranchAddress("particle_photonclass", particle_photonclass);
    slimtree->SetBranchStatus("particle_converted", 1);
    slimtree->SetBranchAddress("particle_converted", particle_converted);

    slimtree->SetBranchStatus(Form("cluster_E_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_E_%s", clusternodename.c_str()), cluster_E);
    slimtree->SetBranchStatus(Form("cluster_Et_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_Et_%s", clusternodename.c_str()), cluster_Et);
    slimtree->SetBranchStatus(Form("cluster_Eta_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_Eta_%s", clusternodename.c_str()), cluster_Eta);
    slimtree->SetBranchStatus(Form("cluster_Phi_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_Phi_%s", clusternodename.c_str()), cluster_Phi);
    slimtree->SetBranchStatus(Form("cluster_prob_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_prob_%s", clusternodename.c_str()), cluster_prob);
    slimtree->SetBranchStatus(Form("cluster_CNN_prob_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_CNN_prob_%s", clusternodename.c_str()), cluster_CNN_prob);
    slimtree->SetBranchStatus(Form("cluster_truthtrkID_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_truthtrkID_%s", clusternodename.c_str()), cluster_truthtrkID);
    slimtree->SetBranchStatus(Form("cluster_pid_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_pid_%s", clusternodename.c_str()), cluster_pid);
    slimtree->SetBranchStatus(Form("cluster_iso_02_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_02_%s", clusternodename.c_str()), cluster_iso_02);
    slimtree->SetBranchStatus(Form("cluster_iso_03_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_%s", clusternodename.c_str()), cluster_iso_03);
    slimtree->SetBranchStatus(Form("cluster_iso_04_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_04_%s", clusternodename.c_str()), cluster_iso_04);
    slimtree->SetBranchStatus(Form("cluster_e1_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e1_%s", clusternodename.c_str()), cluster_e1);
    slimtree->SetBranchStatus(Form("cluster_e2_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e2_%s", clusternodename.c_str()), cluster_e2);
    slimtree->SetBranchStatus(Form("cluster_e3_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e3_%s", clusternodename.c_str()), cluster_e3);
    slimtree->SetBranchStatus(Form("cluster_e4_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e4_%s", clusternodename.c_str()), cluster_e4);
    slimtree->SetBranchStatus(Form("cluster_et1_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_et1_%s", clusternodename.c_str()), cluster_et1);
    slimtree->SetBranchStatus(Form("cluster_et2_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_et2_%s", clusternodename.c_str()), cluster_et2);
    slimtree->SetBranchStatus(Form("cluster_et3_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_et3_%s", clusternodename.c_str()), cluster_et3);
    slimtree->SetBranchStatus(Form("cluster_et4_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_et4_%s", clusternodename.c_str()), cluster_et4);
    slimtree->SetBranchStatus(Form("cluster_weta_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_weta_%s", clusternodename.c_str()), cluster_weta);
    slimtree->SetBranchStatus(Form("cluster_wphi_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_wphi_%s", clusternodename.c_str()), cluster_wphi);
    slimtree->SetBranchStatus(Form("cluster_ietacent_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_ietacent_%s", clusternodename.c_str()), cluster_ietacent);
    slimtree->SetBranchStatus(Form("cluster_iphicent_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iphicent_%s", clusternodename.c_str()), cluster_iphicent);
    slimtree->SetBranchStatus(Form("cluster_detamax_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_detamax_%s", clusternodename.c_str()), cluster_detamax);
    slimtree->SetBranchStatus(Form("cluster_dphimax_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_dphimax_%s", clusternodename.c_str()), cluster_dphimax);
    slimtree->SetBranchStatus(Form("cluster_weta_cogx_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_weta_cogx_%s", clusternodename.c_str()), cluster_weta_cogx);
    slimtree->SetBranchStatus(Form("cluster_wphi_cogx_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_wphi_cogx_%s", clusternodename.c_str()), cluster_wphi_cogx);
    slimtree->SetBranchStatus(Form("cluster_nsaturated_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_nsaturated_%s", clusternodename.c_str()), cluster_nsaturated);

    slimtree->SetBranchStatus(Form("cluster_e11_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e11_%s", clusternodename.c_str()), cluster_e11);
    slimtree->SetBranchStatus(Form("cluster_e22_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e22_%s", clusternodename.c_str()), cluster_e22);
    slimtree->SetBranchStatus(Form("cluster_e13_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e13_%s", clusternodename.c_str()), cluster_e13);
    slimtree->SetBranchStatus(Form("cluster_e15_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e15_%s", clusternodename.c_str()), cluster_e15);
    slimtree->SetBranchStatus(Form("cluster_e17_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e17_%s", clusternodename.c_str()), cluster_e17);
    slimtree->SetBranchStatus(Form("cluster_e31_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e31_%s", clusternodename.c_str()), cluster_e31);
    slimtree->SetBranchStatus(Form("cluster_e51_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e51_%s", clusternodename.c_str()), cluster_e51);
    slimtree->SetBranchStatus(Form("cluster_e71_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e71_%s", clusternodename.c_str()), cluster_e71);
    slimtree->SetBranchStatus(Form("cluster_e33_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e33_%s", clusternodename.c_str()), cluster_e33);
    slimtree->SetBranchStatus(Form("cluster_e35_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e35_%s", clusternodename.c_str()), cluster_e35);
    slimtree->SetBranchStatus(Form("cluster_e37_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e37_%s", clusternodename.c_str()), cluster_e37);
    slimtree->SetBranchStatus(Form("cluster_e53_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e53_%s", clusternodename.c_str()), cluster_e53);
    slimtree->SetBranchStatus(Form("cluster_e73_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e73_%s", clusternodename.c_str()), cluster_e73);
    slimtree->SetBranchStatus(Form("cluster_e55_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e55_%s", clusternodename.c_str()), cluster_e55);
    slimtree->SetBranchStatus(Form("cluster_e57_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e57_%s", clusternodename.c_str()), cluster_e57);
    slimtree->SetBranchStatus(Form("cluster_e75_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e75_%s", clusternodename.c_str()), cluster_e75);
    slimtree->SetBranchStatus(Form("cluster_e77_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e77_%s", clusternodename.c_str()), cluster_e77);
    slimtree->SetBranchStatus(Form("cluster_w32_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_w32_%s", clusternodename.c_str()), cluster_w32);
    slimtree->SetBranchStatus(Form("cluster_e32_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e32_%s", clusternodename.c_str()), cluster_e32);
    slimtree->SetBranchStatus(Form("cluster_w72_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_w72_%s", clusternodename.c_str()), cluster_w72);
    slimtree->SetBranchStatus(Form("cluster_e72_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e72_%s", clusternodename.c_str()), cluster_e72);
    slimtree->SetBranchStatus(Form("cluster_w52_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_w52_%s", clusternodename.c_str()), cluster_w52);

    slimtree->SetBranchStatus(Form("cluster_e_array_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e_array_%s", clusternodename.c_str()), cluster_e_array);
    // slimtree->SetBranchAddress(Form("cluster_adc_array_%s", clusternodename.c_str()), &cluster_adc_array);
    // slimtree->SetBranchAddress(Form("cluster_e_array_idx_%s", clusternodename.c_str()), &cluster_e_array_idx);
    // slimtree->SetBranchAddress(Form("cluster_status_array_%s", clusternodename.c_str()), &cluster_status_array);
    slimtree->SetBranchStatus(Form("cluster_time_array_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_time_array_%s", clusternodename.c_str()), cluster_time_array);
    slimtree->SetBranchStatus(Form("cluster_ownership_array_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_ownership_array_%s", clusternodename.c_str()), cluster_ownership_array);

    slimtree->SetBranchStatus(Form("cluster_iso_03_emcal_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_emcal_%s", clusternodename.c_str()), cluster_iso_03_emcal);
    slimtree->SetBranchStatus(Form("cluster_iso_03_hcalin_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_hcalin_%s", clusternodename.c_str()), cluster_iso_03_hcalin);
    slimtree->SetBranchStatus(Form("cluster_iso_03_hcalout_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_hcalout_%s", clusternodename.c_str()), cluster_iso_03_hcalout);
    slimtree->SetBranchStatus(Form("cluster_iso_03_60_emcal_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_60_emcal_%s", clusternodename.c_str()), cluster_iso_03_60_emcal);
    slimtree->SetBranchStatus(Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()), cluster_iso_03_60_hcalin);
    slimtree->SetBranchStatus(Form("cluster_iso_03_60_hcalout_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_60_hcalout_%s", clusternodename.c_str()), cluster_iso_03_60_hcalout);
    slimtree->SetBranchStatus(Form("cluster_iso_03_120_emcal_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_120_emcal_%s", clusternodename.c_str()), cluster_iso_03_120_emcal);
    slimtree->SetBranchStatus(Form("cluster_iso_03_120_hcalin_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_120_hcalin_%s", clusternodename.c_str()), cluster_iso_03_120_hcalin);
    slimtree->SetBranchStatus(Form("cluster_iso_03_120_hcalout_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_120_hcalout_%s", clusternodename.c_str()), cluster_iso_03_120_hcalout);

    slimtree->SetBranchStatus(Form("cluster_ihcal_et_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_ihcal_et_%s", clusternodename.c_str()), cluster_ihcal_et);
    slimtree->SetBranchStatus(Form("cluster_ohcal_et_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_ohcal_et_%s", clusternodename.c_str()), cluster_ohcal_et);
    slimtree->SetBranchStatus(Form("cluster_ihcal_et22_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_ihcal_et22_%s", clusternodename.c_str()), cluster_ihcal_et22);
    slimtree->SetBranchStatus(Form("cluster_ohcal_et22_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_ohcal_et22_%s", clusternodename.c_str()), cluster_ohcal_et22);
    slimtree->SetBranchStatus(Form("cluster_ihcal_et33_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_ihcal_et33_%s", clusternodename.c_str()), cluster_ihcal_et33);
    slimtree->SetBranchStatus(Form("cluster_ohcal_et33_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_ohcal_et33_%s", clusternodename.c_str()), cluster_ohcal_et33);
    slimtree->SetBranchStatus(Form("cluster_ihcal_ieta_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_ihcal_ieta_%s", clusternodename.c_str()), cluster_ihcal_ieta);
    slimtree->SetBranchStatus(Form("cluster_ihcal_iphi_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_ihcal_iphi_%s", clusternodename.c_str()), cluster_ihcal_iphi);
    slimtree->SetBranchStatus(Form("cluster_ohcal_ieta_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_ohcal_ieta_%s", clusternodename.c_str()), cluster_ohcal_ieta);
    slimtree->SetBranchStatus(Form("cluster_ohcal_iphi_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_ohcal_iphi_%s", clusternodename.c_str()), cluster_ohcal_iphi);
    //bdt score
    slimtree->SetBranchStatus(Form("cluster_bdt_%s_%s", clusternodename.c_str(), bdt_model_name.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_bdt_%s_%s", clusternodename.c_str(), bdt_model_name.c_str()), cluster_bdt);

    slimtree->SetBranchStatus("njet_truth", 1);
    slimtree->SetBranchAddress("njet_truth", &njet_truth);
    slimtree->SetBranchStatus("jet_truth_E", 1);
    slimtree->SetBranchAddress("jet_truth_E", jet_truth_E);
    slimtree->SetBranchStatus("jet_truth_Pt", 1);
    slimtree->SetBranchAddress("jet_truth_Pt", jet_truth_Pt);
    slimtree->SetBranchStatus("jet_truth_Eta", 1);
    slimtree->SetBranchAddress("jet_truth_Eta", jet_truth_Eta);
    slimtree->SetBranchStatus("jet_truth_Phi", 1);
    slimtree->SetBranchAddress("jet_truth_Phi", jet_truth_Phi);

    slimtree->SetBranchStatus("njet", 1);
    slimtree->SetBranchAddress("njet", &njet);
    slimtree->SetBranchStatus("jet_E", 1);
    slimtree->SetBranchAddress("jet_E", jet_E);
    slimtree->SetBranchStatus("jet_Pt", 1);
    slimtree->SetBranchAddress("jet_Pt", jet_Pt);
    slimtree->SetBranchStatus("jet_Eta", 1);
    slimtree->SetBranchAddress("jet_Eta", jet_Eta);
    slimtree->SetBranchStatus("jet_Phi", 1);
    slimtree->SetBranchAddress("jet_Phi", jet_Phi);

    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");



    TRandom3 *rand = new TRandom3(0);
    std::set<int> skiprunnumbers = {47698, 51489, 51721, 51725, 53284};
    int nentries = slimtree->GetEntries();
    for (int ientry = 0; ientry < nentries; ientry++)
    {

        if (ientry == 16152886)
            continue;
        // new calib entry
        // if (ientry == 16095884)
        //     continue;

        // if (ientry == 22626876)
        //     continue;

        // std::cout<<ientry<<std::endl;
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        Long64_t bytesRead = slimtree->GetEntry(ientry);
        if (bytesRead <= 0) {
            std::cout << "Error reading entry " << ientry << std::endl;
            continue;
        }
        std::cout<<"vertexz: "<<vertexz<<std::endl;
        
    }


    if (!issim)
    {
        SaveYamlToRoot(fout, configname.c_str());
    }

}
