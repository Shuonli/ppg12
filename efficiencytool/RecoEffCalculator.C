#include <iostream>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TSystem.h>
#include <TEfficiency.h>
#include <yaml-cpp/yaml.h>
// unfolding
#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>

// R__LOAD_LIBRARY(/sphenix/user/egm2153/calib_study/JetValidation/analysis/roounfold/libRooUnfold.so)

void RecoEffCalculator(const std::string &configname = "config.yaml", const std::string filetype = "data")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    bool issim = true;

    if(filetype == "data")
    {
        issim = false;
    }

    // bool isbackground = (bool)configYaml["input"]["isbackground"].as<int>();
    bool isbackground = false;

    std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();

    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();

    std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;

    if(!issim)
    {
        infilename = configYaml["input"]["data_file"].as<std::string>();
    }

    std::cout<<"infilename: "<<infilename<<std::endl;


    float max_photon_lower = 0;
    float max_photon_upper = 100;

    const float photon5cross = 2.017e+08 * 0.000442571;
    const float photon10cross = 3.690e+07 * 0.000181474;
    const float photon20cross = 1.571e+05 * 0.000673448;

    float weight = 1.0;

    if (filetype == "photon5")
    {
        max_photon_lower = 5;
        max_photon_upper = 12;
        weight = photon5cross / photon20cross;
    }
    else if (filetype == "photon10")
    {
        max_photon_lower = 12;
        max_photon_upper = 25;
        weight = photon10cross / photon20cross;
    }
    else if (filetype == "photon20")
    {
        max_photon_lower = 25;
        max_photon_upper = 100;
        weight = 1.0;
    }
    

    // std::string infilename = "/sphenix/tg/tg01/commissioning/CaloCalibWG/sli/ppg12/ana450/condorout/combine.root";
    TFile *ftreein = new TFile(infilename.c_str(), "READ");
    std::string var_type = configYaml["output"]["var_type"].as<std::string>();

    std::string outfilename = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + filetype + "_" + var_type + ".root";

    std::string responsefilename = configYaml["output"]["response_outfile"].as<std::string>() + "_" + filetype + "_" + var_type + ".root";

    if(!issim)
    {
        outfilename = configYaml["output"]["data_outfile"].as<std::string>() + "_" + var_type + ".root";
        //unfolding is only for sim
        responsefilename = "bla.root";
    }

    std::cout<<"outfilename: "<<outfilename<<std::endl;
    std::cout<<"responsefilename: "<<responsefilename<<std::endl;

    TTree *slimtree = (TTree *)ftreein->Get(configYaml["input"]["tree"].as<std::string>().c_str());

    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();

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

    std::vector<float> pT_bins_truth = configYaml["analysis"]["pT_bins_truth"].as<std::vector<float>>();
    int n_pT_bins_truth = pT_bins_truth.size() - 1;
    double pT_bin_edges_truth[n_pT_bins_truth + 1];

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

    // getting cuts from the config file
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

    // non tight cuts
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

    int mbdnorthhit, mbdsouthhit;
    int pythiaid, nparticles;
    int ncluster;

    float energy_scale;

    float vertexz;
    float vertexz_truth;

    static const int nparticlesmax = 10000;
    static const int nclustercontainermx = 10000;

    float particle_E[nparticlesmax], particle_Pt[nparticlesmax], particle_Eta[nparticlesmax], particle_Phi[nparticlesmax], particle_truth_iso_02[nparticlesmax], particle_truth_iso_03[nparticlesmax], particle_truth_iso_04[nparticlesmax];
    int particle_pid[nparticlesmax], particle_trkid[nparticlesmax], particle_photonclass[nparticlesmax], particle_converted[nparticlesmax];
    float cluster_E[nclustercontainermx], cluster_Et[nclustercontainermx], cluster_Eta[nclustercontainermx], cluster_Phi[nclustercontainermx], cluster_prob[nclustercontainermx], cluster_iso_02[nclustercontainermx], cluster_iso_03[nclustercontainermx], cluster_iso_04[nclustercontainermx], cluster_e1[nclustercontainermx], cluster_e2[nclustercontainermx], cluster_e3[nclustercontainermx], cluster_e4[nclustercontainermx], cluster_et1[nclustercontainermx], cluster_et2[nclustercontainermx], cluster_et3[nclustercontainermx], cluster_et4[nclustercontainermx], cluster_weta[nclustercontainermx], cluster_wphi[nclustercontainermx], cluster_ietacent[nclustercontainermx], cluster_iphicent[nclustercontainermx], cluster_e11[nclustercontainermx], cluster_e22[nclustercontainermx], cluster_e13[nclustercontainermx], cluster_e15[nclustercontainermx], cluster_e17[nclustercontainermx], cluster_e31[nclustercontainermx], cluster_e51[nclustercontainermx], cluster_e71[nclustercontainermx], cluster_e33[nclustercontainermx], cluster_e35[nclustercontainermx], cluster_e37[nclustercontainermx], cluster_e53[nclustercontainermx], cluster_e73[nclustercontainermx], cluster_e55[nclustercontainermx], cluster_e57[nclustercontainermx], cluster_e75[nclustercontainermx], cluster_e77[nclustercontainermx], cluster_w32[nclustercontainermx], cluster_e32[nclustercontainermx], cluster_w72[nclustercontainermx], cluster_e72[nclustercontainermx], cluster_ihcal_et[nclustercontainermx], cluster_ohcal_et[nclustercontainermx], cluster_ihcal_et22[nclustercontainermx], cluster_ohcal_et22[nclustercontainermx], cluster_ihcal_et33[nclustercontainermx], cluster_ohcal_et33[nclustercontainermx];
    float cluster_w52[nclustercontainermx];
    int cluster_truthtrkID[nclustercontainermx], cluster_pid[nclustercontainermx];

    int cluster_detamax[nclustercontainermx], cluster_dphimax[nclustercontainermx], cluster_ihcal_ieta[nclustercontainermx], cluster_ihcal_iphi[nclustercontainermx], cluster_ohcal_ieta[nclustercontainermx], cluster_ohcal_iphi[nclustercontainermx];
    float cluster_CNN_prob[nclustercontainermx];

    slimtree->SetBranchAddress("mbdnorthhit", &mbdnorthhit);
    slimtree->SetBranchAddress("mbdsouthhit", &mbdsouthhit);
    slimtree->SetBranchAddress("pythiaid", &pythiaid);
    slimtree->SetBranchAddress("nparticles", &nparticles);
    slimtree->SetBranchAddress(Form("ncluster_%s", clusternodename.c_str()), &ncluster);
    slimtree->SetBranchAddress("vertexz", &vertexz);
    slimtree->SetBranchAddress("vertexz_truth", &vertexz_truth);
    slimtree->SetBranchAddress("energy_scale", &energy_scale);

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

    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    TH1F *h_max_photon_pT = new TH1F("h_max_photon_pT", "Max Photon pT", 1000, 0, 100);
    TH1F *h_photon_pT = new TH1F("h_photon_pT", "Photon pT", 1000, 0, 100);
    TH1F *h_max_direct_pT = new TH1F("h_max_direct_pT", "Max Direct Photon pT", 1000, 0, 100);
    TH1F *h_direct_pT = new TH1F("h_direct_pT", "Direct Photon pT", 1000, 0, 100);
    TH1F *h_max_frag_pT = new TH1F("h_max_frag_pT", "Max Fragmentation Photon pT", 1000, 0, 100);
    TH1F *h_frag_pT = new TH1F("h_frag_pT", "Fragmentation Photon pT", 1000, 0, 100);
    TH1F *h_max_decay_pT = new TH1F("h_max_decay_pT", "Max Decay Photon pT", 1000, 0, 100);
    TH1F *h_decay_photon_pT = new TH1F("h_decay_photon_pT", "Decay Photon pT", 1000, 0, 100);

    // TEfficiency for conversion and reco
    TEfficiency::EStatOption effopt = TEfficiency::kBUniform;
    TEfficiency *eff_reco = new TEfficiency("eff_reco", "Reco Efficiency", 40, 10, 50, 50, -1, 1);
    eff_reco->SetStatisticOption(effopt);
    TEfficiency *eff_id = new TEfficiency("eff_id", "ID Efficiency", 40, 10, 50, 50, -1, 1);
    eff_id->SetStatisticOption(effopt);
    TEfficiency *eff_converts = new TEfficiency("eff_converts", "Conversion Prob", 40, 10, 50, 50, -1, 1);
    eff_converts->SetStatisticOption(effopt);

    std::vector<TEfficiency *> eff_reco_eta;
    std::vector<TEfficiency *> eff_iso_eta;
    std::vector<TEfficiency *> eff_id_eta;
    std::vector<TEfficiency *> eff_converts_eta;

    std::vector<TH1D *> h_tight_iso_cluster_signal;
    std::vector<TH1D *> h_tight_noniso_cluster_signal;
    std::vector<TH1D *> h_nontight_iso_cluster_signal;
    std::vector<TH1D *> h_nontight_noniso_cluster_signal;

    std::vector<TH2D *> h_singal_reco_isoET;
    std::vector<TH2D *> h_singal_truth_isoET;

    // here are for the plots we gonna make for both data and simulation
    std::vector<TH1D *> h_tight_iso_cluster;
    std::vector<TH1D *> h_tight_noniso_cluster;
    std::vector<TH1D *> h_nontight_iso_cluster;
    std::vector<TH1D *> h_nontight_noniso_cluster;

    // unfold response matrix
    std::vector<RooUnfoldResponse *> responses_full;
    std::vector<RooUnfoldResponse *> responses_half;
    // id histogram for unfolding
    std::vector<TH1D *> h_pT_truth_response;
    std::vector<TH1D *> h_pT_reco_response;

    std::vector<TH1D *> h_pT_reco_fake;

    std::vector<TH1D *> h_pT_truth_half_response;
    std::vector<TH1D *> h_pT_reco_half_response;

    std::vector<TH1D *> h_pT_truth_secondhalf_response;
    std::vector<TH1D *> h_pT_reco_secondhalf_response;

    // direct and fragmentation photon pT vs truth iso ET
    std::vector<TH2D *> h_direct_pT_truth_isoET;
    std::vector<TH2D *> h_frag_pT_truth_isoET;

    // n cluster per photon
    std::vector<TH2D *> h_ncluster_truth;

    // energy resolution
    std::vector<TH2D *> h_pT_truth_reco;

    // isolation profile for debugging reasons
    std::vector<std::vector<TH1D *>> h_tight_cluster_pT;
    h_tight_cluster_pT.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH1D *>> h_nontight_cluster_pT;
    h_nontight_cluster_pT.resize(eta_bins.size() - 1);
    // truth iso vs reco iso for different pT bins and eta bins
    std::vector<std::vector<TH2D *>> h_iso_truth_reco;
    h_iso_truth_reco.resize(eta_bins.size() - 1);
    // response vs. isoET
    std::vector<std::vector<TH2D *>> h_response_isoET;
    h_response_isoET.resize(eta_bins.size() - 1);

    for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
    {
        eff_reco_eta.push_back(new TEfficiency(Form("eff_reco_eta_%d", ieta),
                                               Form("Reco Efficiency %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                               n_pT_bins, pT_bin_edges));
        eff_reco_eta[ieta]->SetStatisticOption(effopt);
        eff_reco_eta[ieta]->SetWeight(weight);

        eff_id_eta.push_back(new TEfficiency(Form("eff_id_eta_%d", ieta),
                                             Form("ID Efficiency %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                             n_pT_bins, pT_bin_edges));
        eff_id_eta[ieta]->SetStatisticOption(effopt);
        eff_id_eta[ieta]->SetWeight(weight);

        eff_converts_eta.push_back(new TEfficiency(Form("eff_converts_eta_%d", ieta),
                                                   Form("Conversion Prob %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                   n_pT_bins, pT_bin_edges));
        eff_converts_eta[ieta]->SetStatisticOption(effopt);
        eff_converts_eta[ieta]->SetWeight(weight);

        eff_iso_eta.push_back(new TEfficiency(Form("eff_iso_eta_%d", ieta),
                                              Form("Iso Efficiency %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                              n_pT_bins, pT_bin_edges));
        eff_iso_eta[ieta]->SetStatisticOption(effopt);
        eff_iso_eta[ieta]->SetWeight(weight);

        h_tight_iso_cluster_signal.push_back(new TH1D(Form("h_tight_iso_cluster_signal_%d", ieta),
                                                      Form("Tight Iso Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                      n_pT_bins, pT_bin_edges));
        h_tight_noniso_cluster_signal.push_back(new TH1D(Form("h_tight_noniso_cluster_signal_%d", ieta),
                                                         Form("Tight Non-Iso Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                         n_pT_bins, pT_bin_edges));
        h_nontight_iso_cluster_signal.push_back(new TH1D(Form("h_nontight_iso_cluster_signal_%d", ieta),
                                                         Form("Non-Tight Iso Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                         n_pT_bins, pT_bin_edges));
        h_nontight_noniso_cluster_signal.push_back(new TH1D(Form("h_nontight_noniso_cluster_signal_%d", ieta),
                                                            Form("Non-Tight Non-Iso Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                            n_pT_bins, pT_bin_edges));

        h_singal_reco_isoET.push_back(new TH2D(Form("h_singal_reco_isoET_%d", ieta),
                                               Form("Signal Reco Iso ET %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                               400, 0, 50, 4400, -5, 50));

        h_singal_truth_isoET.push_back(new TH2D(Form("h_singal_truth_isoET_%d", ieta),
                                                Form("Signal Truth Iso ET %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                400, 0, 50, 4400, -5, 50));

        // plot for both data and simulation
        h_tight_iso_cluster.push_back(new TH1D(Form("h_tight_iso_cluster_%d", ieta),
                                               Form("Tight Iso Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                               n_pT_bins, pT_bin_edges));
        h_tight_noniso_cluster.push_back(new TH1D(Form("h_tight_noniso_cluster_%d", ieta),
                                                  Form("Tight Non-Iso Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                  n_pT_bins, pT_bin_edges));
        h_nontight_iso_cluster.push_back(new TH1D(Form("h_nontight_iso_cluster_%d", ieta),
                                                  Form("Non-Tight Iso Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                  n_pT_bins, pT_bin_edges));
        h_nontight_noniso_cluster.push_back(new TH1D(Form("h_nontight_noniso_cluster_%d", ieta),
                                                     Form("Non-Tight Non-Iso Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                     n_pT_bins, pT_bin_edges));

        // unfold histograms
        h_pT_truth_response.push_back(new TH1D(Form("h_pT_truth_response_%d", ieta),
                                               Form("Truth pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                               n_pT_bins_truth, pT_bin_edges_truth));

        h_pT_reco_response.push_back(new TH1D(Form("h_pT_reco_response_%d", ieta),
                                              Form("Reco pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                              n_pT_bins, pT_bin_edges));

        h_pT_reco_fake.push_back(new TH1D(Form("h_pT_reco_fake_%d", ieta),
                                          Form("Reco Fake pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                          n_pT_bins, pT_bin_edges));

        std::cout << h_pT_reco_response[ieta]->GetNbinsX() << std::endl;

        TH2D *h_response_full = new TH2D(Form("h_response_full_%d", ieta),
                                         Form("Response Matrix %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                         n_pT_bins, pT_bin_edges, n_pT_bins_truth, pT_bin_edges_truth);

        responses_full.push_back(new RooUnfoldResponse((const TH1 *)h_pT_reco_response[ieta], (const TH1 *)h_pT_truth_response[ieta], h_response_full, Form("response_matrix_full_%d", ieta), "", false));

        h_pT_truth_half_response.push_back(new TH1D(Form("h_pT_truth_half_response_%d", ieta),
                                                    Form("Truth pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                    n_pT_bins_truth, pT_bin_edges_truth));

        h_pT_reco_half_response.push_back(new TH1D(Form("h_pT_reco_half_response_%d", ieta),
                                                   Form("Reco pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                   n_pT_bins, pT_bin_edges));

        h_pT_reco_secondhalf_response.push_back(new TH1D(Form("h_pT_reco_secondhalf_response_%d", ieta),
                                                         Form("Reco pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                         n_pT_bins, pT_bin_edges));

        h_pT_truth_secondhalf_response.push_back(new TH1D(Form("h_pT_truth_secondhalf_response_%d", ieta),
                                                          Form("Truth pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                          n_pT_bins_truth, pT_bin_edges_truth));

        TH2D *h_response_half = new TH2D(Form("h_response_half_%d", ieta),
                                         Form("Response Matrix %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                         n_pT_bins, pT_bin_edges, n_pT_bins_truth, pT_bin_edges_truth);

        responses_half.push_back(new RooUnfoldResponse(h_pT_reco_half_response[ieta], h_pT_truth_half_response[ieta], h_response_half, Form("response_matrix_half_%d", ieta), ""));

        h_ncluster_truth.push_back(new TH2D(Form("h_ncluster_truth_%d", ieta),
                                            Form("N Cluster From Truth %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                            400, 0, 40, 10, 0, 10));

        h_pT_truth_reco.push_back(new TH2D(Form("h_pT_truth_reco_%d", ieta),
                                           Form("Truth Reco %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                           400, 5, 40, 150, 0, 1.5));

        h_direct_pT_truth_isoET.push_back(new TH2D(Form("h_direct_pT_truth_isoET_%d", ieta),
                                                   Form("Direct Photon pT vs Truth Iso ET %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                   400, 0, 40, 400, 0, 40));

        h_frag_pT_truth_isoET.push_back(new TH2D(Form("h_frag_pT_truth_isoET_%d", ieta),
                                                 Form("Fragmentation Photon pT vs Truth Iso ET %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                 400, 0, 40, 400, 0, 40));

        for (int ipt = 0; ipt < n_pT_bins; ipt++)
        {
            h_tight_cluster_pT[ieta].push_back(new TH1D(Form("h_tight_isoET_%d_%d", ieta, ipt),
                                                        Form("Tight Iso ET %.1f < eta < %.1f %.1f < pT < %.1f", eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                        400, -10, 30));
            h_nontight_cluster_pT[ieta].push_back(new TH1D(Form("h_nontight_isoET_%d_%d", ieta, ipt),
                                                           Form("Non-Tight Iso ET %.1f < eta < %.1f %.1f < pT < %.1f", eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                           400, -10, 30));
            h_iso_truth_reco[ieta].push_back(new TH2D(Form("h_iso_truth_reco_%d_%d", ieta, ipt),
                                                      Form("Iso Truth Reco %.1f < eta < %.1f %.1f < pT < %.1f", eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                      300, 0, 30, 400, -10, 30));

            h_response_isoET[ieta].push_back(new TH2D(Form("h_response_isoET_%d_%d", ieta, ipt),
                                                      Form("Response Iso ET %.1f < eta < %.1f %.1f < pT < %.1f", eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                      150, 0, 1.5, 400, -10, 30));
        }
    }

    int nentries = slimtree->GetEntries();
    for (int ientry = 0; ientry < nentries; ientry++)
    {
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;

        slimtree->GetEntry(ientry);
        // check vertex cut
        if (abs(vertexz) > vertexcut)
            continue;
        std::map<int, int> particle_trkidmap;
        // map for photon reco eff
        std::map<int, bool> photon_converts;
        // map for photon reco eff
        std::map<int, bool> photon_reco;
        // map for iso eff
        std::map<int, bool> photon_iso;
        std::map<int, float> photon_iso_ET;
        // map for id eff
        std::map<int, bool> photon_id;
        // map for n cluster per truth photon
        std::map<int, int> photon_ncluster;
        if (issim)
        {
            float maxphotonpT = 0;
            int maxphotonclass = 0;
            for (int iparticle = 0; iparticle < nparticles; iparticle++)
            {
                particle_trkidmap[particle_trkid[iparticle]] = iparticle;
                if (!isbackground)
                {
                    if (particle_pid[iparticle] == 22)
                    {
                        if (particle_Pt[iparticle] > maxphotonpT)
                        {
                            maxphotonpT = particle_Pt[iparticle];
                            maxphotonclass = particle_photonclass[iparticle];
                        }

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
                        if (particle_photonclass[iparticle] < 3) // direct or fragmentation
                        {
                            if (truthisoET < truthisocut)
                            {
                                // we check for conversion here
                                photon_converts[iparticle] = (particle_converted[iparticle] > 0);
                                // initialize reco to false
                                photon_reco[iparticle] = false;
                                // initialize id to false
                                photon_id[iparticle] = false;
                                // initialize iso to false
                                photon_iso[iparticle] = false;

                                photon_ncluster[iparticle] = 0;
                            }
                            photon_iso_ET[iparticle] = truthisoET;
                        }
                    }
                }
                // if we are calculating things for the background sample
                else
                {
                    // place holder to make the code work
                    if ((particle_pid[iparticle] != 22) || (particle_photonclass[iparticle] == 3))
                    {
                        photon_converts[iparticle] = (particle_converted[iparticle] > 0);
                        // initialize reco to false
                        photon_reco[iparticle] = false;
                        // initialize id to false
                        photon_id[iparticle] = false;
                        // initialize iso to false
                        photon_iso[iparticle] = false;

                        photon_ncluster[iparticle] = 0;
                    }
                }
            }
            if ((maxphotonpT > max_photon_upper) || (maxphotonpT < max_photon_lower))
            {
                continue;
            }
            // if (abs(particle_Eta[iparticle]) < 0.7)
            {
                h_max_photon_pT->Fill(maxphotonpT, weight);
                if (maxphotonclass == 1)
                {
                    h_max_direct_pT->Fill(maxphotonpT, weight);
                }
                else if (maxphotonclass == 2)
                {
                    h_max_frag_pT->Fill(maxphotonpT, weight);
                }
                else if (maxphotonclass == 3)
                {
                    h_max_decay_pT->Fill(maxphotonpT, weight);
                }
            }

            for (int iparticle = 0; iparticle < nparticles; iparticle++)
            {
                particle_trkidmap[particle_trkid[iparticle]] = iparticle;
                if (!isbackground)
                {
                    if (particle_pid[iparticle] == 22)
                    {
                        if (particle_Pt[iparticle] > maxphotonpT)
                        {
                            maxphotonpT = particle_Pt[iparticle];
                            maxphotonclass = particle_photonclass[iparticle];
                        }

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

                        if (abs(particle_Eta[iparticle]) < 0.7)
                        {
                            if (particle_photonclass[iparticle] == 1)
                            {
                                h_direct_pT_truth_isoET[0]->Fill(particle_Pt[iparticle], truthisoET, weight);
                                h_direct_pT->Fill(particle_Pt[iparticle], weight);
                            }
                            else if (particle_photonclass[iparticle] == 2)
                            {
                                h_frag_pT_truth_isoET[0]->Fill(particle_Pt[iparticle], truthisoET, weight);
                                h_frag_pT->Fill(particle_Pt[iparticle], weight);
                            }
                            else if (particle_photonclass[iparticle] == 3)
                            {
                                h_decay_photon_pT->Fill(particle_Pt[iparticle], weight);
                            }

                            h_photon_pT->Fill(particle_Pt[iparticle], weight);
                        }
                    }
                }
            }
        }

        for (int icluster = 0; icluster < ncluster; icluster++)
        {
            // need ET > 10 GeV
            if (cluster_Et[icluster] < reco_min_ET)
                continue;

            float rhad22 = (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]) / (cluster_Et[icluster] + (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]));
            float rhad33 = (cluster_ihcal_et33[icluster] + cluster_ohcal_et33[icluster]) / (cluster_Et[icluster] + (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]));

            float reta77 = cluster_e37[icluster] / cluster_e77[icluster];
            float rphi77 = cluster_e73[icluster] / cluster_e77[icluster];

            float reta55 = cluster_e35[icluster] / cluster_e55[icluster];
            float rphi55 = cluster_e53[icluster] / cluster_e55[icluster];

            float reta = cluster_e33[icluster] / cluster_e73[icluster];

            float rphi = cluster_e33[icluster] / cluster_e37[icluster];

            float re11_E = cluster_e11[icluster] / cluster_E[icluster];

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

            bool tight = false;
            bool nontight = false;
            bool iso = false;
            bool noniso = false;
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

            // need to update to a function to find tight and non tight
            if (reta77 > tight_reta77_min &&
                reta77 < tight_reta77_max &&
                rhad33 > tight_rhad33_min &&
                rhad33 < tight_rhad33_max &&
                cluster_w72[icluster] > tight_w72_min &&
                cluster_w72[icluster] < tight_w72_max &&
                re11_E > tight_re11_E_min &&
                re11_E < tight_re11_E_max &&
                cluster_CNN_prob[icluster] > tight_CNN_min &&
                cluster_CNN_prob[icluster] < tight_CNN_max)
            {

                tight = true;
            }
            if (reta77 > non_tight_reta77_min &&
                     reta77 < non_tight_reta77_max &&
                     rhad33 > non_tight_rhad33_min &&
                     rhad33 < non_tight_rhad33_max &&
                     cluster_w72[icluster] > non_tight_w72_min &&
                     cluster_w72[icluster] < non_tight_w72_max &&
                     re11_E > non_tight_re11_E_min &&
                     re11_E < non_tight_re11_E_max &&
                     cluster_CNN_prob[icluster] > non_tight_CNN_min &&
                     cluster_CNN_prob[icluster] < non_tight_CNN_max)
            {
                nontight = true;
            }

            if (issim)
            {
                if (particle_trkidmap.find(cluster_truthtrkID[icluster]) == particle_trkidmap.end())
                {
                    // std::cout<<"trackid: "<<cluster_truthtrkID[icluster]<<std::endl;
                    // std::cout << "Error: cluster_truthtrkID not found in particle_trkidmap" << std::endl;
                    continue;
                }
                int iparticle = particle_trkidmap[cluster_truthtrkID[icluster]];

                // delta R cut
                float deta = cluster_Eta[icluster] - particle_Eta[iparticle];
                float dphi = cluster_Phi[icluster] - particle_Phi[iparticle];
                if (dphi > M_PI)
                {
                    dphi = 2 * M_PI - dphi;
                }
                float dR = sqrt(deta * deta + dphi * dphi);

                // truth iso vs. reco iso

                if (photon_iso_ET.find(iparticle) != photon_iso_ET.end())
                {
                    h_iso_truth_reco[etabin][pTbin]->Fill(photon_iso_ET[iparticle], recoisoET, weight);

                    if (photon_reco.find(iparticle) == photon_reco.end())
                    {
                        // then it is non truth signal, if it pass the reco, iso, and tight cuts, then it is a fake
                        if (iso && tight && (dR < eff_dR))
                        {
                            h_pT_reco_fake[etabin]->Fill(cluster_Et[icluster], weight);
                        }
                    }
                }

                // iparticle has to be in the photon map
                if (photon_reco.find(iparticle) == photon_reco.end())
                {
                    continue;
                }

                photon_ncluster[iparticle]++;

                if (dR < eff_dR)
                // if ( (dR < eff_dR) && ( (cluster_Et[icluster] / particle_Pt[iparticle]) > 0.8) )
                {
                    
                    // if(photon_converts[iparticle]) continue;

                    photon_reco[iparticle] = true;

                    h_pT_truth_reco[etabin]->Fill(particle_Pt[iparticle], cluster_Et[icluster] / particle_Pt[iparticle], weight);

                    h_response_isoET[etabin][pTbin]->Fill(cluster_Et[icluster] / particle_Pt[iparticle], recoisoET, weight);

                    if (iso)
                    {
                        photon_iso[iparticle] = true;
                    }

                    if (tight)
                    {
                        photon_id[iparticle] = true;
                    }

                    if (tight && iso)
                    {
                        h_tight_iso_cluster_signal[etabin]->Fill(cluster_Et[icluster], weight);
                        // fill the response matrix
                        h_pT_truth_response[etabin]->Fill(particle_Pt[iparticle], weight);
                        h_pT_reco_response[etabin]->Fill(cluster_Et[icluster], weight);
                        responses_full[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight);
                        if (ientry < (nentries / 2))
                        {
                            h_pT_truth_half_response[etabin]->Fill(particle_Pt[iparticle], weight);
                            h_pT_reco_half_response[etabin]->Fill(cluster_Et[icluster], weight);
                            responses_half[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight);
                        }
                        else
                        {
                            h_pT_truth_secondhalf_response[etabin]->Fill(particle_Pt[iparticle], weight);
                            h_pT_reco_secondhalf_response[etabin]->Fill(cluster_Et[icluster], weight);
                        }
                    }
                    if (tight && noniso)
                    {
                        h_tight_noniso_cluster_signal[etabin]->Fill(cluster_Et[icluster], weight);
                    }
                    if (nontight && iso)
                    {
                        h_nontight_iso_cluster_signal[etabin]->Fill(cluster_Et[icluster], weight);
                    }
                    if (nontight && noniso)
                    {
                        h_nontight_noniso_cluster_signal[etabin]->Fill(cluster_Et[icluster], weight);
                    }

                    h_singal_reco_isoET[etabin]->Fill(cluster_Et[icluster], recoisoET, weight);
                    h_singal_truth_isoET[etabin]->Fill(particle_Pt[iparticle], recoisoET, weight);
                }
            }
            if (tight && iso)
            {
                h_tight_iso_cluster[etabin]->Fill(cluster_Et[icluster], weight);
                // h_pT_reco_response[etabin]->Fill(cluster_Et[icluster]);
            }
            if (tight && noniso)
            {
                h_tight_noniso_cluster[etabin]->Fill(cluster_Et[icluster], weight);
            }
            if (nontight && iso)
            {
                h_nontight_iso_cluster[etabin]->Fill(cluster_Et[icluster], weight);
            }
            if (nontight && noniso)
            {
                h_nontight_noniso_cluster[etabin]->Fill(cluster_Et[icluster], weight);
            }

            if (tight)
            {
                h_tight_cluster_pT[etabin][pTbin]->Fill(recoisoET, weight);
            }
            if (nontight)
            {
                h_nontight_cluster_pT[etabin][pTbin]->Fill(recoisoET, weight);
            }
        }

        // go over the map and fill the TEfficiency
        for (auto it = photon_reco.begin(); it != photon_reco.end(); ++it)
        {
            float photon_pT = particle_Pt[it->first];
            float photon_eta = particle_Eta[it->first];
            int etabin = -1;
            for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
            {
                if (photon_eta > eta_bins[ieta] && photon_eta < eta_bins[ieta + 1])
                {
                    etabin = ieta;
                    break;
                }
            }
            if (etabin == -1)
            {
                continue;
            }
            eff_converts_eta[etabin]->Fill(photon_converts[it->first], photon_pT);

            h_ncluster_truth[etabin]->Fill(photon_pT, photon_ncluster[it->first], weight);

            // std::cout<<"photon_pT: "<<photon_pT<<" photon_eta: "<<photon_eta<<" photon_converts: "<<photon_converts[it->first]<<std::endl;
            // std::cout<<"effieicncy: "<<eff_converts->GetEfficiency(0)<<std::endl;
            eff_reco_eta[etabin]->Fill(photon_reco[it->first], photon_pT);

            if (photon_reco[it->first])
            {
                // if(!photon_converts[it->first])
                {
                    eff_iso_eta[etabin]->Fill(photon_iso[it->first], photon_pT);
                }

                if (photon_iso[it->first])
                {
                    eff_id_eta[etabin]->Fill(photon_id[it->first], photon_pT);
                }
            }
        }
    }
    TFile *fresponse = new TFile(responsefilename.c_str(), "RECREATE");
    for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
    {

        responses_full[ieta]->Write();

        responses_half[ieta]->Write();

        h_pT_truth_response[ieta]->Write();
        h_pT_reco_response[ieta]->Write();

        h_pT_truth_half_response[ieta]->Write();
        h_pT_reco_half_response[ieta]->Write();

        h_pT_truth_secondhalf_response[ieta]->Write();
        h_pT_reco_secondhalf_response[ieta]->Write();
    }

    fout->Write();
    fout->Close();

    fresponse->Write();
    fresponse->Close();
}
