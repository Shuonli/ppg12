#include <yaml-cpp/yaml.h>
#include <TMVA/RBDT.hxx>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TChain.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <cmath>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>

#include "MbdPileupHelper.h"

const float TIME_SAMPLE_NS = 17.6;

// TODO: fill in actual run ranges for crossing angle classification
bool is_0mrad(int run) { return (run >= 47000 && run <= 48000); }
bool is_1p5mrad(int run) { return (run >= 49000 && run <= 53500); }

void DoubleInteractionCheck(
    const std::string &configname = "config_showershape.yaml",
    const std::string filetype = "data",
    bool doinclusive = true)
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    bool issim = true;
    bool isbackground = false;

    if (filetype == "data")
    {
        issim = false;
    }

    std::string inclusive_str = doinclusive ? "_inclusive" : "";

    std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
    std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;

    if (!issim)
    {
        infilename = configYaml["input"]["data_file"].as<std::string>();
    }

    std::cout << "infilename: " << infilename << std::endl;

    // ---------------------------------------------------------------
    // NPB score configuration
    // ---------------------------------------------------------------
    float npb_score_cut = configYaml["analysis"]["npb_score_cut"].as<float>(0.5);

    int mbd_t0_correction_on = configYaml["analysis"]["mbd_t0_correction_on"].as<int>(1);
    std::string mbd_t0_correction_file =
        configYaml["analysis"]["mbd_t0_correction_file"].as<std::string>("/sphenix/user/shuhangli/ppg12/efficiencytool/MbdOut.corr");
    std::map<int, float> mbd_t0_correction;
    if (mbd_t0_correction_on && !issim)
    {
        std::ifstream file(mbd_t0_correction_file);
        if (!file.is_open())
        {
            std::cerr << "[MBD t0] WARNING: cannot open correction file: " << mbd_t0_correction_file << std::endl;
        }
        else
        {
            std::string line;
            while (std::getline(file, line))
            {
                if (line.empty()) continue;
                std::istringstream iss(line);
                int runnumber = 0;
                float t0 = 0.0;
                if (!(iss >> runnumber >> t0)) continue;
                mbd_t0_correction[runnumber] = t0;
            }
            std::cout << "[MBD t0] Loaded " << mbd_t0_correction.size() << " corrections" << std::endl;
        }
    }

    // ---------------------------------------------------------------
    // Cross-section weights (same as ShowerShapeCheck.C)
    // ---------------------------------------------------------------
    float max_photon_lower = 0;
    float max_photon_upper = 100;

    const float photon5cross = 146359.3;
    const float photon10cross = 6944.675;
    const float photon20cross = 130.4461;

    const float jet10cross = 3.997e+06;
    const float jet15cross = 4.073e+05;
    const float jet20cross = 6.218e+04;
    const float jet30cross = 2.502e+03;
    const float jet50cross = 7.2695;

    float max_jet_lower = 0;
    float max_jet_upper = 100;
    float weight = 1.0;
    float vertex_weight = 1.0;
    float cross_weight = 1.0;
    float cluster_ET_upper = 100;

    if (filetype == "photon5")
    {
        max_photon_lower = 0;
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
        max_photon_upper = 200;
        weight = 1.0;
    }
    else if (filetype == "jet10")
    {
        max_jet_lower = 10;
        max_jet_upper = 15;
        cluster_ET_upper = 18;
        weight = jet10cross / jet50cross;
        isbackground = true;
    }
    else if (filetype == "jet15")
    {
        max_jet_lower = 15;
        max_jet_upper = 20;
        cluster_ET_upper = 23;
        weight = jet15cross / jet50cross;
        isbackground = true;
    }
    else if (filetype == "jet20")
    {
        max_jet_lower = 20;
        max_jet_upper = 30;
        cluster_ET_upper = 33;
        weight = jet20cross / jet50cross;
        isbackground = true;
    }
    else if (filetype == "jet30")
    {
        max_jet_lower = 30;
        max_jet_upper = 50;
        cluster_ET_upper = 45;
        weight = jet30cross / jet50cross;
        isbackground = true;
    }
    else if (filetype == "jet50")
    {
        max_jet_lower = 50;
        max_jet_upper = 100;
        weight = jet50cross / jet50cross;
        isbackground = true;
    }
    cross_weight = weight;

    // ---------------------------------------------------------------
    // Vertex reweighting (sim)
    // ---------------------------------------------------------------
    TH1 *h_vertex_reweight = nullptr;
    int vertex_reweight_on = 1;
    std::string vertex_reweight_file = "results/vertex_reweight.root";

    TFile *fvtx = nullptr;
    TH1D *h_vertexz_data = nullptr;

    if (issim)
    {
        vertex_reweight_on = configYaml["analysis"]["vertex_reweight_on"].as<int>(1);
        vertex_reweight_file =
            configYaml["analysis"]["vertex_reweight_file"].as<std::string>("results/vertex_reweight.root");

        fvtx = TFile::Open(vertex_reweight_file.c_str(), "READ");
        if (!fvtx || fvtx->IsZombie())
        {
            std::cerr << "[VertexReweight] ERROR: cannot open " << vertex_reweight_file << std::endl;
            return;
        }

        if (vertex_reweight_on)
        {
            TH1 *htmp = dynamic_cast<TH1 *>(fvtx->Get("h_vertexz_ratio_data_over_mccombined"));
            if (!htmp)
            {
                std::cerr << "[VertexReweight] ERROR: histogram not found" << std::endl;
                fvtx->Close();
                return;
            }
            h_vertex_reweight = dynamic_cast<TH1 *>(htmp->Clone("h_vertexz_ratio_clone"));
            h_vertex_reweight->SetDirectory(nullptr);
            h_vertex_reweight->Sumw2();
            std::cout << "[VertexReweight] Loaded from " << vertex_reweight_file << std::endl;
        }

        // Data vertex distribution for toy double interaction
        TH1D *htmp_data = dynamic_cast<TH1D *>(fvtx->Get("h_vertexz_data"));
        if (!htmp_data)
        {
            std::cerr << "[VertexData] ERROR: h_vertexz_data not found in " << vertex_reweight_file << std::endl;
            fvtx->Close();
            return;
        }
        h_vertexz_data = dynamic_cast<TH1D *>(htmp_data->Clone("h_vertexz_data_clone"));
        h_vertexz_data->SetDirectory(nullptr);
        std::cout << "[VertexData] Loaded data vertex distribution (" << h_vertexz_data->GetEntries() << " entries)" << std::endl;
    }

    TRandom3 rng(42);

    // ---------------------------------------------------------------
    // NPB TMVA model (sim only, for rescoring)
    // ---------------------------------------------------------------
    TMVA::Experimental::RBDT *npb_bdt_ptr = nullptr;
    TMVA::Experimental::RBDT *ss_bdt_ptr = nullptr;
    std::string ss_bdt_model_name = "base_v3E";

    if (issim)
    {
        std::string npb_tmva_file = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/npb_models/npb_score_tmva.root";
        if (configYaml["analysis"] && configYaml["analysis"]["npb_tmva_file"])
        {
            npb_tmva_file = configYaml["analysis"]["npb_tmva_file"].as<std::string>();
        }
        if (gSystem->AccessPathName(npb_tmva_file.c_str()))
        {
            std::cerr << "ERROR: NPB TMVA model not found at '" << npb_tmva_file << "'" << std::endl;
            return;
        }
        std::cout << "[NPB] Loading TMVA model from " << npb_tmva_file << std::endl;
        npb_bdt_ptr = new TMVA::Experimental::RBDT("myBDT", npb_tmva_file);

        // Showershape BDT TMVA model
        if (configYaml["input"] && configYaml["input"]["bdt_model_name"])
        {
            ss_bdt_model_name = configYaml["input"]["bdt_model_name"].as<std::string>();
        }
        std::string ss_bdt_file = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_"
                                  + ss_bdt_model_name + "_single_tmva.root";
        if (gSystem->AccessPathName(ss_bdt_file.c_str()))
        {
            std::cerr << "ERROR: SS BDT model not found at '" << ss_bdt_file << "'" << std::endl;
            return;
        }
        std::cout << "[SS-BDT] Loading model (" << ss_bdt_model_name << ") from " << ss_bdt_file << std::endl;
        ss_bdt_ptr = new TMVA::Experimental::RBDT("myBDT", ss_bdt_file);
    }

    // ---------------------------------------------------------------
    // Config parsing
    // ---------------------------------------------------------------
    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();
    std::string bdt_model_name = configYaml["input"]["bdt_model_name"].as<std::string>("base");

    int iso_threshold = configYaml["analysis"]["iso_threshold"].as<int>(0);
    int iso_hcalonly = configYaml["analysis"]["iso_hcalonly"].as<int>(0);
    float iso_emcalinnerr = configYaml["analysis"]["iso_emcalinnerr"].as<float>(0.0);

    int n_nt_fail = configYaml["analysis"]["n_nt_fail"].as<int>(1);
    int weta_fail = configYaml["analysis"]["weta_fail"].as<int>(0);
    int wphi_fail = configYaml["analysis"]["wphi_fail"].as<int>(0);
    int e11_to_e33_fail = configYaml["analysis"]["e11_to_e33_fail"].as<int>(0);
    int e32_to_e35_fail = configYaml["analysis"]["e32_to_e35_fail"].as<int>(0);
    int et1_fail = configYaml["analysis"]["et1_fail"].as<int>(0);
    int bdt_fail = configYaml["analysis"]["bdt_fail"].as<int>(0);

    int weta_on = configYaml["analysis"]["weta_on"].as<int>(1);
    int wphi_on = configYaml["analysis"]["wphi_on"].as<int>(1);
    int e11_to_e33_on = configYaml["analysis"]["e11_to_e33_on"].as<int>(1);
    int e32_to_e35_on = configYaml["analysis"]["e32_to_e35_on"].as<int>(1);
    int et1_on = configYaml["analysis"]["et1_on"].as<int>(1);
    int et2_on = configYaml["analysis"]["et2_on"].as<int>(1);
    int et3_on = configYaml["analysis"]["et3_on"].as<int>(1);
    int et4_on = configYaml["analysis"]["et4_on"].as<int>(1);
    int bdt_on = configYaml["analysis"]["bdt_on"].as<int>(1);

    float truthisocut = configYaml["analysis"]["truth_iso_max"].as<float>();
    float recoiso_min = configYaml["analysis"]["reco_iso_min"].as<float>();
    float recoiso_max_b = configYaml["analysis"]["reco_iso_max_b"].as<float>();
    float recoiso_max_s = configYaml["analysis"]["reco_iso_max_s"].as<float>();
    float recononiso_min_shift = configYaml["analysis"]["reco_noniso_min_shift"].as<float>();
    float recononiso_max = configYaml["analysis"]["reco_noniso_max"].as<float>();

    float vertexcut = configYaml["analysis"]["vertex_cut"].as<float>();
    std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();
    int n_eta_bins = (int)eta_bins.size() - 1;

    std::vector<float> pT_bins = configYaml["analysis"]["pT_bins"].as<std::vector<float>>();
    int n_pT_bins = pT_bins.size() - 1;
    double pT_bin_edges[n_pT_bins + 1];
    std::copy(pT_bins.begin(), pT_bins.end(), pT_bin_edges);

    int conesize = configYaml["analysis"]["cone_size"].as<int>();
    float reco_min_ET = configYaml["analysis"]["reco_min_ET"].as<float>();

    // Trigger
    std::vector<int> trigger_used;
    {
        YAML::Node trigNode = configYaml["analysis"]["trigger_used"];
        if (trigNode && trigNode.IsSequence())
        {
            trigger_used = trigNode.as<std::vector<int>>();
        }
        else
        {
            trigger_used.push_back(configYaml["analysis"]["trigger_used"].as<int>());
        }
    }

    float mc_iso_shift = configYaml["analysis"]["mc_iso_shift"].as<float>(0.0);
    float mc_iso_scale = configYaml["analysis"]["mc_iso_scale"].as<float>(1.2);

    // ---------------------------------------------------------------
    // Tight cuts
    // ---------------------------------------------------------------
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

    float tight_bdt_max = configYaml["analysis"]["tight"]["bdt_max"].as<float>(1.0);
    float tight_bdt_min = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0.0);

    // ---------------------------------------------------------------
    // Non-tight cuts
    // ---------------------------------------------------------------
    float non_tight_weta_cogx_max = configYaml["analysis"]["non_tight"]["weta_cogx_max"].as<float>();
    float non_tight_weta_cogx_min = configYaml["analysis"]["non_tight"]["weta_cogx_min"].as<float>();
    float non_tight_wphi_cogx_max = configYaml["analysis"]["non_tight"]["wphi_cogx_max"].as<float>();
    float non_tight_wphi_cogx_min = configYaml["analysis"]["non_tight"]["wphi_cogx_min"].as<float>();
    float non_tight_prob_max = configYaml["analysis"]["non_tight"]["prob_max"].as<float>();
    float non_tight_prob_min = configYaml["analysis"]["non_tight"]["prob_min"].as<float>();
    float non_tight_e11_over_e33_max = configYaml["analysis"]["non_tight"]["e11_over_e33_max"].as<float>();
    float non_tight_e11_over_e33_min = configYaml["analysis"]["non_tight"]["e11_over_e33_min"].as<float>();
    float non_tight_e32_over_e35_max = configYaml["analysis"]["non_tight"]["e32_over_e35_max"].as<float>();
    float non_tight_e32_over_e35_min = configYaml["analysis"]["non_tight"]["e32_over_e35_min"].as<float>();
    float non_tight_et1_max = configYaml["analysis"]["non_tight"]["et1_max"].as<float>();
    float non_tight_et1_min = configYaml["analysis"]["non_tight"]["et1_min"].as<float>();
    float non_tight_et4_max = configYaml["analysis"]["non_tight"]["et4_max"].as<float>();
    float non_tight_et4_min = configYaml["analysis"]["non_tight"]["et4_min"].as<float>();
    float non_tight_bdt_max = configYaml["analysis"]["non_tight"]["bdt_max"].as<float>(1.0);
    float non_tight_bdt_min = configYaml["analysis"]["non_tight"]["bdt_min"].as<float>(0.0);

    // ---------------------------------------------------------------
    // Common cuts
    // ---------------------------------------------------------------
    float common_prob_max = configYaml["analysis"]["common"]["prob_max"].as<float>();
    float common_prob_min = configYaml["analysis"]["common"]["prob_min"].as<float>();
    float common_e11_over_e33_max = configYaml["analysis"]["common"]["e11_over_e33_max"].as<float>();
    float common_e11_over_e33_min = configYaml["analysis"]["common"]["e11_over_e33_min"].as<float>();
    float common_cluster_weta_cogx_bound = configYaml["analysis"]["common"]["cluster_weta_cogx_bound"].as<float>(0.8);

    float clusterescale = configYaml["analysis"]["cluster_escale"].as<float>(1.0);

    // ---------------------------------------------------------------
    // TChain + TTreeReader
    // ---------------------------------------------------------------
    std::string treename = configYaml["input"]["tree"].as<std::string>();
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());

    // NPB score branch detection
    std::string npb_score_branch = Form("cluster_npb_score_%s", clusternodename.c_str());
    if (!chain.GetBranch(npb_score_branch.c_str()))
    {
        if (chain.GetBranch("cluster_npb_score"))
        {
            npb_score_branch = "cluster_npb_score";
        }
        else
        {
            npb_score_branch.clear();
            std::cerr << "[NPBScore] WARNING: cannot find NPB score branch" << std::endl;
            return;
        }
    }

    TTreeReader reader(&chain);

    // Event-level
    TTreeReaderValue<int> mbdnorthhit(reader, "mbdnorthhit");
    TTreeReaderValue<int> mbdsouthhit(reader, "mbdsouthhit");
    TTreeReaderValue<int> runnumber(reader, "runnumber");
    TTreeReaderValue<float> mbd_time(reader, "mbd_time");
    TTreeReaderArray<float> mbd_north_time(reader, "mbdnortht");
    TTreeReaderArray<float> mbd_south_time(reader, "mbdsoutht");
    TTreeReaderArray<float> mbd_north_charge(reader, "mbdnorthq");
    TTreeReaderArray<float> mbd_south_charge(reader, "mbdsouthq");
    TTreeReaderValue<int> nparticles(reader, "nparticles");
    TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
    TTreeReaderValue<float> vertexz(reader, "vertexz");
    TTreeReaderArray<Bool_t> scaledtrigger(reader, "scaledtrigger");

    // Particle arrays
    TTreeReaderArray<float> particle_Pt(reader, "particle_Pt");
    TTreeReaderArray<float> particle_Eta(reader, "particle_Eta");
    TTreeReaderArray<float> particle_truth_iso_02(reader, "particle_truth_iso_02");
    TTreeReaderArray<float> particle_truth_iso_03(reader, "particle_truth_iso_03");
    TTreeReaderArray<float> particle_truth_iso_04(reader, "particle_truth_iso_04");
    TTreeReaderArray<int> particle_pid(reader, "particle_pid");
    TTreeReaderArray<int> particle_trkid(reader, "particle_trkid");
    TTreeReaderArray<int> particle_photonclass(reader, "particle_photonclass");

    // Cluster arrays
    TTreeReaderArray<float> cluster_E(reader, Form("cluster_E_%s", clusternodename.c_str()));
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

    // Cluster energy ratio arrays
    TTreeReaderArray<float> cluster_e11(reader, Form("cluster_e11_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e22(reader, Form("cluster_e22_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e13(reader, Form("cluster_e13_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e15(reader, Form("cluster_e15_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e17(reader, Form("cluster_e17_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e31(reader, Form("cluster_e31_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e51(reader, Form("cluster_e51_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e71(reader, Form("cluster_e71_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e32(reader, Form("cluster_e32_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e33(reader, Form("cluster_e33_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e35(reader, Form("cluster_e35_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e37(reader, Form("cluster_e37_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e53(reader, Form("cluster_e53_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w32(reader, Form("cluster_w32_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w52(reader, Form("cluster_w52_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w72(reader, Form("cluster_w72_%s", clusternodename.c_str()));

    // Isolation with threshold
    TTreeReaderArray<float> cluster_iso_03_70_emcal(reader, Form("cluster_iso_03_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalin(reader, Form("cluster_iso_03_70_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalout(reader, Form("cluster_iso_03_70_hcalout_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_005_70_emcal(reader, Form("cluster_iso_005_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_01_70_emcal(reader, Form("cluster_iso_01_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_02_70_emcal(reader, Form("cluster_iso_02_70_emcal_%s", clusternodename.c_str()));

    // BDT score
    TTreeReaderArray<float> cluster_bdt(reader, Form("cluster_bdt_%s_%s", clusternodename.c_str(), bdt_model_name.c_str()));

    // NPB score
    std::unique_ptr<TTreeReaderArray<float>> cluster_npb_score;
    if (!npb_score_branch.empty())
    {
        cluster_npb_score = std::make_unique<TTreeReaderArray<float>>(reader, npb_score_branch.c_str());
    }

    // Truth jets
    TTreeReaderValue<int> njet_truth(reader, "njet_truth");
    TTreeReaderArray<float> jet_truth_Pt(reader, "jet_truth_Pt");

    // ---------------------------------------------------------------
    // Output file
    // ---------------------------------------------------------------
    std::string outfilename = configYaml["output"]["eff_outfile"].as<std::string>()
                              + "_double_interaction_check_" + filetype + inclusive_str + ".root";
    if (!issim)
    {
        outfilename = configYaml["output"]["data_outfile"].as<std::string>()
                      + "_double_interaction_check.root";
    }
    std::cout << "outfilename: " << outfilename << std::endl;

    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    fout->cd();

    // ---------------------------------------------------------------
    // Histogram booking
    // ---------------------------------------------------------------

    // Task 1: ABCD yields as 2D (ET vs NPB score)
    std::vector<TH2D *> h_tight_iso_npb, h_tight_noniso_npb, h_nontight_iso_npb, h_nontight_noniso_npb;

    // Task 2: MBD avg sigma vs run number (data only)
    TProfile *h_mbd_avgsigma_vs_run_all = nullptr;
    TProfile *h_mbd_avgsigma_vs_run_0mrad = nullptr;
    TProfile *h_mbd_avgsigma_vs_run_1p5mrad = nullptr;
    TH2D *h2_mbd_avgsigma_vs_run = nullptr;

    if (!issim)
    {
        h_mbd_avgsigma_vs_run_all = new TProfile("h_mbd_avgsigma_vs_run_all",
            "MBD Avg #sigma_{t} vs Run Number;Run Number;MBD Avg #sigma_{t} [ns]",
            7000, 46000, 53000);
        h_mbd_avgsigma_vs_run_0mrad = new TProfile("h_mbd_avgsigma_vs_run_0mrad",
            "MBD Avg #sigma_{t} vs Run (0 mrad);Run Number;MBD Avg #sigma_{t} [ns]",
            7000, 46000, 53000);
        h_mbd_avgsigma_vs_run_1p5mrad = new TProfile("h_mbd_avgsigma_vs_run_1p5mrad",
            "MBD Avg #sigma_{t} vs Run (1.5 mrad);Run Number;MBD Avg #sigma_{t} [ns]",
            7000, 46000, 53000);
        h2_mbd_avgsigma_vs_run = new TH2D("h2_mbd_avgsigma_vs_run",
            "MBD Avg #sigma_{t} vs Run Number;Run Number;MBD Avg #sigma_{t} [ns]",
            7000, 46000, 53000, 100, 0, 5.0);
    }

    // Task 3: ABCD yields as 2D (ET vs MBD avg sigma)
    std::vector<TH2D *> h_tight_iso_mbdsigma, h_tight_noniso_mbdsigma, h_nontight_iso_mbdsigma, h_nontight_noniso_mbdsigma;

    // Task 4: Single vs double interaction ABCD (sim only)
    std::vector<TH1D *> h_tight_iso_single, h_tight_noniso_single, h_nontight_iso_single, h_nontight_noniso_single;
    std::vector<TH1D *> h_tight_iso_single_signal, h_tight_noniso_single_signal, h_nontight_iso_single_signal, h_nontight_noniso_single_signal;
    std::vector<TH1D *> h_tight_iso_double, h_tight_noniso_double, h_nontight_iso_double, h_nontight_noniso_double;
    std::vector<TH1D *> h_tight_iso_double_signal, h_tight_noniso_double_signal, h_nontight_iso_double_signal, h_nontight_noniso_double_signal;

    // Inclusive signal histograms: signal tagged for all MC (including jet)
    std::vector<TH1D *> h_tight_iso_single_incsig, h_tight_noniso_single_incsig, h_nontight_iso_single_incsig, h_nontight_noniso_single_incsig;
    std::vector<TH1D *> h_tight_iso_double_incsig, h_tight_noniso_double_incsig, h_nontight_iso_double_incsig, h_nontight_noniso_double_incsig;

    for (int ieta = 0; ieta < n_eta_bins; ieta++)
    {
        // Task 1: NPB score vs ET (2D)
        h_tight_iso_npb.push_back(new TH2D(Form("h_tight_iso_cluster_npb_%d", ieta),
            Form("Tight Iso NPB %.1f<#eta<%.1f;Cluster E_{T} [GeV];NPB Score", eta_bins[ieta], eta_bins[ieta + 1]),
            n_pT_bins, pT_bin_edges, 500, 0, 1));
        h_tight_noniso_npb.push_back(new TH2D(Form("h_tight_noniso_cluster_npb_%d", ieta),
            Form("Tight NonIso NPB %.1f<#eta<%.1f;Cluster E_{T} [GeV];NPB Score", eta_bins[ieta], eta_bins[ieta + 1]),
            n_pT_bins, pT_bin_edges, 500, 0, 1));
        h_nontight_iso_npb.push_back(new TH2D(Form("h_nontight_iso_cluster_npb_%d", ieta),
            Form("NonTight Iso NPB %.1f<#eta<%.1f;Cluster E_{T} [GeV];NPB Score", eta_bins[ieta], eta_bins[ieta + 1]),
            n_pT_bins, pT_bin_edges, 500, 0, 1));
        h_nontight_noniso_npb.push_back(new TH2D(Form("h_nontight_noniso_cluster_npb_%d", ieta),
            Form("NonTight NonIso NPB %.1f<#eta<%.1f;Cluster E_{T} [GeV];NPB Score", eta_bins[ieta], eta_bins[ieta + 1]),
            n_pT_bins, pT_bin_edges, 500, 0, 1));

        // Task 3: MBD avg sigma vs ET (2D, data only)
        if (!issim)
        {
            h_tight_iso_mbdsigma.push_back(new TH2D(Form("h_tight_iso_cluster_mbdsigma_%d", ieta),
                Form("Tight Iso MBD#sigma %.1f<#eta<%.1f;Cluster E_{T} [GeV];MBD Avg #sigma_{t} [ns]", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges, 50, 0, 5.0));
            h_tight_noniso_mbdsigma.push_back(new TH2D(Form("h_tight_noniso_cluster_mbdsigma_%d", ieta),
                Form("Tight NonIso MBD#sigma %.1f<#eta<%.1f;Cluster E_{T} [GeV];MBD Avg #sigma_{t} [ns]", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges, 50, 0, 5.0));
            h_nontight_iso_mbdsigma.push_back(new TH2D(Form("h_nontight_iso_cluster_mbdsigma_%d", ieta),
                Form("NonTight Iso MBD#sigma %.1f<#eta<%.1f;Cluster E_{T} [GeV];MBD Avg #sigma_{t} [ns]", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges, 50, 0, 5.0));
            h_nontight_noniso_mbdsigma.push_back(new TH2D(Form("h_nontight_noniso_cluster_mbdsigma_%d", ieta),
                Form("NonTight NonIso MBD#sigma %.1f<#eta<%.1f;Cluster E_{T} [GeV];MBD Avg #sigma_{t} [ns]", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges, 50, 0, 5.0));
        }

        // Task 4: Single/double interaction (sim)
        if (issim)
        {
            h_tight_iso_single.push_back(new TH1D(Form("h_tight_iso_cluster_single_%d", ieta),
                Form("Tight Iso Single %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_tight_noniso_single.push_back(new TH1D(Form("h_tight_noniso_cluster_single_%d", ieta),
                Form("Tight NonIso Single %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_nontight_iso_single.push_back(new TH1D(Form("h_nontight_iso_cluster_single_%d", ieta),
                Form("NonTight Iso Single %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_nontight_noniso_single.push_back(new TH1D(Form("h_nontight_noniso_cluster_single_%d", ieta),
                Form("NonTight NonIso Single %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));

            h_tight_iso_single_signal.push_back(new TH1D(Form("h_tight_iso_cluster_single_signal_%d", ieta),
                Form("Tight Iso Single Signal %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_tight_noniso_single_signal.push_back(new TH1D(Form("h_tight_noniso_cluster_single_signal_%d", ieta),
                Form("Tight NonIso Single Signal %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_nontight_iso_single_signal.push_back(new TH1D(Form("h_nontight_iso_cluster_single_signal_%d", ieta),
                Form("NonTight Iso Single Signal %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_nontight_noniso_single_signal.push_back(new TH1D(Form("h_nontight_noniso_cluster_single_signal_%d", ieta),
                Form("NonTight NonIso Single Signal %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));

            h_tight_iso_double.push_back(new TH1D(Form("h_tight_iso_cluster_double_%d", ieta),
                Form("Tight Iso Double %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_tight_noniso_double.push_back(new TH1D(Form("h_tight_noniso_cluster_double_%d", ieta),
                Form("Tight NonIso Double %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_nontight_iso_double.push_back(new TH1D(Form("h_nontight_iso_cluster_double_%d", ieta),
                Form("NonTight Iso Double %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_nontight_noniso_double.push_back(new TH1D(Form("h_nontight_noniso_cluster_double_%d", ieta),
                Form("NonTight NonIso Double %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));

            h_tight_iso_double_signal.push_back(new TH1D(Form("h_tight_iso_cluster_double_signal_%d", ieta),
                Form("Tight Iso Double Signal %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_tight_noniso_double_signal.push_back(new TH1D(Form("h_tight_noniso_cluster_double_signal_%d", ieta),
                Form("Tight NonIso Double Signal %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_nontight_iso_double_signal.push_back(new TH1D(Form("h_nontight_iso_cluster_double_signal_%d", ieta),
                Form("NonTight Iso Double Signal %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_nontight_noniso_double_signal.push_back(new TH1D(Form("h_nontight_noniso_cluster_double_signal_%d", ieta),
                Form("NonTight NonIso Double Signal %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));

            // Inclusive signal (signal tagged for all MC including jet)
            h_tight_iso_single_incsig.push_back(new TH1D(Form("h_tight_iso_cluster_single_incsig_%d", ieta),
                Form("Tight Iso Single IncSig %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_tight_noniso_single_incsig.push_back(new TH1D(Form("h_tight_noniso_cluster_single_incsig_%d", ieta),
                Form("Tight NonIso Single IncSig %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_nontight_iso_single_incsig.push_back(new TH1D(Form("h_nontight_iso_cluster_single_incsig_%d", ieta),
                Form("NonTight Iso Single IncSig %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_nontight_noniso_single_incsig.push_back(new TH1D(Form("h_nontight_noniso_cluster_single_incsig_%d", ieta),
                Form("NonTight NonIso Single IncSig %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));

            h_tight_iso_double_incsig.push_back(new TH1D(Form("h_tight_iso_cluster_double_incsig_%d", ieta),
                Form("Tight Iso Double IncSig %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_tight_noniso_double_incsig.push_back(new TH1D(Form("h_tight_noniso_cluster_double_incsig_%d", ieta),
                Form("Tight NonIso Double IncSig %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_nontight_iso_double_incsig.push_back(new TH1D(Form("h_nontight_iso_cluster_double_incsig_%d", ieta),
                Form("NonTight Iso Double IncSig %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
            h_nontight_noniso_double_incsig.push_back(new TH1D(Form("h_nontight_noniso_cluster_double_incsig_%d", ieta),
                Form("NonTight NonIso Double IncSig %.1f<#eta<%.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                n_pT_bins, pT_bin_edges));
        }
    }

    // ---------------------------------------------------------------
    // BDT feature vector builders (sim only, from showershape_vertex_check.C)
    // ---------------------------------------------------------------
    auto buildNpbFeatureVector = [&](int icl, float vtxz) -> std::vector<float> {
        const float e11_over_e33 = (cluster_e33[icl] > 0) ? cluster_e11[icl] / cluster_e33[icl] : 0.0f;
        const float e32_over_e35 = (cluster_e35[icl] > 0) ? cluster_e32[icl] / cluster_e35[icl] : 0.0f;
        const float e11_over_e22 = (cluster_e22[icl] > 0) ? cluster_e11[icl] / cluster_e22[icl] : 0.0f;
        const float e11_over_e13 = (cluster_e13[icl] > 0) ? cluster_e11[icl] / cluster_e13[icl] : 0.0f;
        const float e11_over_e15 = (cluster_e15[icl] > 0) ? cluster_e11[icl] / cluster_e15[icl] : 0.0f;
        const float e11_over_e17 = (cluster_e17[icl] > 0) ? cluster_e11[icl] / cluster_e17[icl] : 0.0f;
        const float e11_over_e31 = (cluster_e31[icl] > 0) ? cluster_e11[icl] / cluster_e31[icl] : 0.0f;
        const float e11_over_e51 = (cluster_e51[icl] > 0) ? cluster_e11[icl] / cluster_e51[icl] : 0.0f;
        const float e11_over_e71 = (cluster_e71[icl] > 0) ? cluster_e11[icl] / cluster_e71[icl] : 0.0f;
        const float e22_over_e33 = (cluster_e33[icl] > 0) ? cluster_e22[icl] / cluster_e33[icl] : 0.0f;
        const float e22_over_e35 = (cluster_e35[icl] > 0) ? cluster_e22[icl] / cluster_e35[icl] : 0.0f;
        const float e22_over_e37 = (cluster_e37[icl] > 0) ? cluster_e22[icl] / cluster_e37[icl] : 0.0f;
        const float e22_over_e53 = (cluster_e53[icl] > 0) ? cluster_e22[icl] / cluster_e53[icl] : 0.0f;
        return {
            cluster_Et[icl], cluster_Eta[icl], vtxz,
            e11_over_e33, e32_over_e35, e11_over_e22,
            e11_over_e13, e11_over_e15, e11_over_e17,
            e11_over_e31, e11_over_e51, e11_over_e71,
            e22_over_e33, e22_over_e35, e22_over_e37, e22_over_e53,
            cluster_weta_cogx[icl], cluster_wphi_cogx[icl],
            cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl],
            cluster_w32[icl], cluster_w52[icl], cluster_w72[icl]
        };
    };

    auto buildSsBdtFeatureVector = [&](int icl, float vtxz) -> std::vector<float> {
        const float e11_over_e33 = (cluster_e33[icl] > 0) ? cluster_e11[icl] / cluster_e33[icl] : 0.0f;
        const float e32_over_e35 = (cluster_e35[icl] > 0) ? cluster_e32[icl] / cluster_e35[icl] : 0.0f;
        if (ss_bdt_model_name == "base_v3E")
        {
            return { cluster_Et[icl], cluster_weta_cogx[icl], cluster_wphi_cogx[icl], vtxz,
                     cluster_Eta[icl], e11_over_e33,
                     cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl], e32_over_e35 };
        }
        else if (ss_bdt_model_name == "base_v3")
        {
            return { cluster_weta_cogx[icl], cluster_wphi_cogx[icl], vtxz,
                     cluster_Eta[icl], e11_over_e33,
                     cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl], e32_over_e35 };
        }
        else if (ss_bdt_model_name == "base_v2E")
        {
            return { cluster_Et[icl], cluster_weta_cogx[icl], cluster_wphi_cogx[icl], vtxz,
                     cluster_Eta[icl], e11_over_e33,
                     cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl] };
        }
        else if (ss_bdt_model_name == "base_v2")
        {
            return { cluster_weta_cogx[icl], cluster_wphi_cogx[icl], vtxz,
                     cluster_Eta[icl], e11_over_e33,
                     cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl] };
        }
        else if (ss_bdt_model_name == "base_v1E")
        {
            return { cluster_Et[icl], cluster_weta_cogx[icl], vtxz,
                     cluster_Eta[icl], e11_over_e33,
                     cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl] };
        }
        else if (ss_bdt_model_name == "base_v1")
        {
            return { cluster_weta_cogx[icl], vtxz, cluster_Eta[icl], e11_over_e33,
                     cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl] };
        }
        else if (ss_bdt_model_name == "base_E")
        {
            return { cluster_Et[icl], vtxz, cluster_Eta[icl], e11_over_e33,
                     cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl] };
        }
        else if (ss_bdt_model_name == "base" || ss_bdt_model_name == "base_vr")
        {
            return { vtxz, cluster_Eta[icl], e11_over_e33,
                     cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl] };
        }
        else
        {
            // Default fallback to base_v3E
            return { cluster_Et[icl], cluster_weta_cogx[icl], cluster_wphi_cogx[icl], vtxz,
                     cluster_Eta[icl], e11_over_e33,
                     cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl], e32_over_e35 };
        }
    };

    // ---------------------------------------------------------------
    // Lambda: classify tight/non-tight given shower shape variables and a BDT score
    // Returns {tight, nontight}
    // ---------------------------------------------------------------
    auto classifyTightNonTight = [&](int icl, float clusterET, float e11_over_e33, float e32_over_e35,
                                     float bdt_score) -> std::pair<bool, bool>
    {
        float tw_max = tight_weta_cogx_max_b + tight_weta_cogx_max_s * clusterET;
        float twp_max = tight_wphi_cogx_max_b + tight_wphi_cogx_max_s * clusterET;
        float tet1_min = tight_et1_min_b + tight_et1_min_s * clusterET;

        bool is_weta_tight = (cluster_weta_cogx[icl] > tight_weta_cogx_min) && (cluster_weta_cogx[icl] < tw_max);
        bool is_wphi_tight = (cluster_wphi_cogx[icl] > tight_wphi_cogx_min) && (cluster_wphi_cogx[icl] < twp_max);
        bool is_et1_tight = (cluster_et1[icl] > tet1_min) && (cluster_et1[icl] < tight_et1_max);
        bool is_et2_tight = (cluster_et2[icl] > tight_et2_min) && (cluster_et2[icl] < tight_et2_max);
        bool is_et3_tight = (cluster_et3[icl] > tight_et3_min) && (cluster_et3[icl] < tight_et3_max);
        bool is_e11e33_tight = (e11_over_e33 > tight_e11_over_e33_min) && (e11_over_e33 < tight_e11_over_e33_max);
        bool is_e32e35_tight = (e32_over_e35 > tight_e32_over_e35_min) && (e32_over_e35 < tight_e32_over_e35_max);
        bool is_et4_tight = (cluster_et4[icl] > tight_et4_min) && (cluster_et4[icl] < tight_et4_max);
        bool is_prob_tight = (cluster_prob[icl] > tight_prob_min) && (cluster_prob[icl] < tight_prob_max);
        bool is_bdt_tight = (bdt_score > tight_bdt_min) && (bdt_score < tight_bdt_max);

        bool tight = is_weta_tight && is_wphi_tight && is_et1_tight && is_et2_tight && is_et3_tight &&
                     is_e11e33_tight && is_e32e35_tight && is_et4_tight && is_prob_tight && is_bdt_tight;

        bool nontight = false;
        if (cluster_weta_cogx[icl] > non_tight_weta_cogx_min && cluster_weta_cogx[icl] < non_tight_weta_cogx_max &&
            cluster_wphi_cogx[icl] > non_tight_wphi_cogx_min && cluster_wphi_cogx[icl] < non_tight_wphi_cogx_max &&
            cluster_prob[icl] > non_tight_prob_min && cluster_prob[icl] < non_tight_prob_max &&
            e11_over_e33 > non_tight_e11_over_e33_min && e11_over_e33 < non_tight_e11_over_e33_max &&
            e32_over_e35 > non_tight_e32_over_e35_min && e32_over_e35 < non_tight_e32_over_e35_max &&
            cluster_et1[icl] > non_tight_et1_min && cluster_et1[icl] < non_tight_et1_max &&
            cluster_et4[icl] > non_tight_et4_min && cluster_et4[icl] < non_tight_et4_max &&
            bdt_score > non_tight_bdt_min && bdt_score < non_tight_bdt_max)
        {
            int nfail = 0;
            if (!is_weta_tight)  nfail += weta_on;
            if (!is_wphi_tight)  nfail += wphi_on;
            if (!is_et1_tight)   nfail += et1_on;
            if (!is_et2_tight)   nfail += et2_on;
            if (!is_et3_tight)   nfail += et3_on;
            if (!is_e11e33_tight) nfail += e11_to_e33_on;
            if (!is_e32e35_tight) nfail += e32_to_e35_on;
            if (!is_et4_tight)   nfail += et4_on;
            if (!is_prob_tight)  nfail++;
            if (!is_bdt_tight)   nfail += bdt_on;

            if (nfail > n_nt_fail)
            {
                bool all_flags_fail = true;
                if (weta_fail && is_weta_tight) all_flags_fail = false;
                if (wphi_fail && is_wphi_tight) all_flags_fail = false;
                if (et1_fail && is_et1_tight)   all_flags_fail = false;
                if (e11_to_e33_fail && is_e11e33_tight) all_flags_fail = false;
                if (e32_to_e35_fail && is_e32e35_tight) all_flags_fail = false;
                if (bdt_fail && is_bdt_tight)   all_flags_fail = false;

                if (all_flags_fail) nontight = true;
            }
        }

        return {tight, nontight};
    };

    // ---------------------------------------------------------------
    // Event loop
    // ---------------------------------------------------------------
    int nentries = chain.GetEntries();
    int ientry = 0;

    while (reader.Next())
    {
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;

        // Trigger selection (data only)
        if (!issim)
        {
            bool any_trigger_fired = false;
            const auto ntrig = scaledtrigger.GetSize();
            for (int itrig : trigger_used)
            {
                if (itrig < 0 || (unsigned int)itrig >= ntrig) continue;
                if (scaledtrigger[itrig] != 0) { any_trigger_fired = true; break; }
            }
            if (!any_trigger_fired) { ientry++; continue; }
        }

        // Vertex reweighting (sim)
        weight = cross_weight;
        vertex_weight = 1.0;
        if (issim && vertex_reweight_on && h_vertex_reweight)
        {
            int bin = h_vertex_reweight->FindBin(*vertexz);
            if (bin < 1) bin = 1;
            if (bin > h_vertex_reweight->GetNbinsX()) bin = h_vertex_reweight->GetNbinsX();
            vertex_weight = h_vertex_reweight->GetBinContent(bin);
            if (!std::isfinite(vertex_weight) || vertex_weight <= 0.0) vertex_weight = 1.0;
            weight *= vertex_weight;
        }

        // MBD pileup metrics
        MbdPileupResult pileup_result = calculateMbdPileupMetrics(
            mbd_north_time, mbd_south_time, mbd_north_charge, mbd_south_charge);

        // Task 2: Fill MBD avg sigma vs run (data only, event-level)
        if (!issim && pileup_result.valid)
        {
            int run = *runnumber;
            h_mbd_avgsigma_vs_run_all->Fill(run, pileup_result.avgsigma);
            h2_mbd_avgsigma_vs_run->Fill(run, pileup_result.avgsigma);
            if (is_0mrad(run))
                h_mbd_avgsigma_vs_run_0mrad->Fill(run, pileup_result.avgsigma);
            if (is_1p5mrad(run))
                h_mbd_avgsigma_vs_run_1p5mrad->Fill(run, pileup_result.avgsigma);
        }

        // MC truth signal classification
        std::set<int> signal_set;
        std::set<int> background_set;
        std::map<int, int> particle_trkidmap;

        if (issim)
        {
            float maxphotonpT = 0;
            for (int ip = 0; ip < *nparticles; ip++)
            {
                particle_trkidmap[particle_trkid[ip]] = ip;
                float truthisoET = 0;
                if (conesize == 4)      truthisoET = particle_truth_iso_04[ip];
                else if (conesize == 3) truthisoET = particle_truth_iso_03[ip];
                else if (conesize == 2) truthisoET = particle_truth_iso_02[ip];

                if (particle_pid[ip] == 22)
                {
                    if (particle_Pt[ip] > maxphotonpT) maxphotonpT = particle_Pt[ip];
                    if (particle_photonclass[ip] < 3)
                    {
                        if (truthisoET < truthisocut) signal_set.insert(ip);
                    }
                    else
                    {
                        background_set.insert(ip);
                    }
                }
                else
                {
                    background_set.insert(ip);
                }
            }

            // pT-hat range cut for photon samples
            if (!isbackground)
            {
                if (maxphotonpT > max_photon_upper || maxphotonpT < max_photon_lower)
                { ientry++; continue; }
            }

            // pT-hat range cut for jet samples
            if (isbackground)
            {
                float maxjetpT = 0;
                for (int ij = 0; ij < *njet_truth; ij++)
                {
                    if (jet_truth_Pt[ij] > maxjetpT) maxjetpT = jet_truth_Pt[ij];
                }
                if (maxjetpT > max_jet_upper || maxjetpT < max_jet_lower)
                { ientry++; continue; }
            }
        }

        // Vertex cut: only applied per-task, not globally
        // Single interaction and Tasks 1/3 require original vertex in range
        // Double interaction requires double_vtx in range
        bool vtx_pass = (std::abs(*vertexz) <= vertexcut);

        // MBD hit requirement
        if (!(*mbdnorthhit >= 1 && *mbdsouthhit >= 1)) { ientry++; continue; }

        // For sim: draw random second vertex for double interaction
        float double_vtx = 0;
        if (issim && h_vertexz_data)
        {
            float random_vertex = h_vertexz_data->GetRandom();
            double_vtx = (*vertexz + random_vertex) / 2.0;
        }

        // MBD time correction (data)
        float mbd_mean_time = *mbd_time;
        if (!issim && mbd_t0_correction_on)
        {
            auto it = mbd_t0_correction.find(*runnumber);
            if (it != mbd_t0_correction.end()) mbd_mean_time -= it->second;
        }

        // ---------------------------------------------------------------
        // Cluster loop
        // ---------------------------------------------------------------
        for (int icl = 0; icl < *ncluster; icl++)
        {
            float clusterET = cluster_Et[icl];
            if (clusterET < reco_min_ET) continue;
            if (isbackground && clusterET > cluster_ET_upper) continue;

            // Sim: truth matching (for signal classification only - still fill all clusters in task 1/3)
            bool is_signal = false;
            bool is_signal_inclusive = false;
            if (issim)
            {
                if (particle_trkidmap.find(cluster_truthtrkID[icl]) != particle_trkidmap.end())
                {
                    int ip = particle_trkidmap[cluster_truthtrkID[icl]];
                    is_signal_inclusive = (signal_set.find(ip) != signal_set.end());
                    if (!isbackground)
                    {
                        is_signal = is_signal_inclusive;
                    }
                    else
                    {
                        // For jet background, none are signal
                        is_signal = false;
                    }
                }
            }

            float e11_over_e33 = (cluster_e33[icl] > 0) ? cluster_e11[icl] / cluster_e33[icl] : 0.0f;
            float e32_over_e35 = (cluster_e35[icl] > 0) ? cluster_e32[icl] / cluster_e35[icl] : 0.0f;

            // Isolation ET
            float recoisoET = -999;
            if (conesize == 4)      recoisoET = cluster_iso_04[icl];
            else if (conesize == 3) recoisoET = cluster_iso_03[icl];
            else if (conesize == 2) recoisoET = cluster_iso_02[icl];

            if (iso_threshold)
            {
                if (iso_hcalonly)
                {
                    recoisoET = cluster_iso_03_70_hcalin[icl] + cluster_iso_03_70_hcalout[icl];
                }
                else
                {
                    recoisoET = cluster_iso_03_70_emcal[icl] + cluster_iso_03_70_hcalin[icl] + cluster_iso_03_70_hcalout[icl];
                    if (iso_emcalinnerr > 0.04 && iso_emcalinnerr < 0.06)
                        recoisoET -= cluster_iso_005_70_emcal[icl];
                    else if (iso_emcalinnerr > 0.09 && iso_emcalinnerr < 0.11)
                        recoisoET -= cluster_iso_01_70_emcal[icl];
                    else if (iso_emcalinnerr > 0.19 && iso_emcalinnerr < 0.21)
                        recoisoET -= cluster_iso_02_70_emcal[icl];
                }
            }

            if (issim)
            {
                recoisoET = recoisoET * mc_iso_scale + mc_iso_shift;
            }

            // Iso/non-iso classification
            float recoiso_max = recoiso_max_b + recoiso_max_s * clusterET;
            float recononiso_min = recoiso_max + recononiso_min_shift;
            bool iso = (recoisoET > recoiso_min && recoisoET < recoiso_max);
            bool noniso = (recoisoET > recononiso_min && recoisoET < recononiso_max);

            // Eta bin
            int etabin = -1;
            for (int ieta = 0; ieta < n_eta_bins; ieta++)
            {
                if (cluster_Eta[icl] > eta_bins[ieta] && cluster_Eta[icl] < eta_bins[ieta + 1])
                { etabin = ieta; break; }
            }
            if (etabin == -1) continue;

            // NPB score
            float npb_score_val = (cluster_npb_score ? (*cluster_npb_score)[icl] : 1.0f);
            bool npb_pass = (npb_score_val >= npb_score_cut);

            // Common cuts (require NPB pass for standard ABCD, but for Task 1 we fill both pass/fail)
            bool common_pass_no_npb =
                cluster_prob[icl] > common_prob_min && cluster_prob[icl] < common_prob_max &&
                e11_over_e33 > common_e11_over_e33_min && e11_over_e33 < common_e11_over_e33_max &&
                cluster_weta_cogx[icl] < common_cluster_weta_cogx_bound;

            // Original BDT score for tight/non-tight classification
            float bdt_score_original = cluster_bdt[icl];

            // Tight/non-tight classification with original BDT
            auto [tight, nontight] = classifyTightNonTight(icl, clusterET, e11_over_e33, e32_over_e35, bdt_score_original);

            // ---------------------------------------------------------------
            // Task 1: ABCD yields as 2D (ET vs NPB score)
            // Requires common cuts (without NPB requirement) + original vertex in range
            // ---------------------------------------------------------------
            if (vtx_pass && common_pass_no_npb)
            {
                if (tight && iso)       h_tight_iso_npb[etabin]->Fill(clusterET, npb_score_val, weight);
                if (tight && noniso)    h_tight_noniso_npb[etabin]->Fill(clusterET, npb_score_val, weight);
                if (nontight && iso)    h_nontight_iso_npb[etabin]->Fill(clusterET, npb_score_val, weight);
                if (nontight && noniso) h_nontight_noniso_npb[etabin]->Fill(clusterET, npb_score_val, weight);
            }

            // For tasks 3 & 4, apply full common cut (including NPB)
            bool common_pass = common_pass_no_npb && npb_pass;
            if (!common_pass) continue;

            // ---------------------------------------------------------------
            // Task 3: Data ABCD as 2D (ET vs MBD avg sigma)
            // Requires original vertex in range
            // ---------------------------------------------------------------
            if (vtx_pass && !issim && pileup_result.valid)
            {
                if (tight && iso)       h_tight_iso_mbdsigma[etabin]->Fill(clusterET, pileup_result.avgsigma, weight);
                if (tight && noniso)    h_tight_noniso_mbdsigma[etabin]->Fill(clusterET, pileup_result.avgsigma, weight);
                if (nontight && iso)    h_nontight_iso_mbdsigma[etabin]->Fill(clusterET, pileup_result.avgsigma, weight);
                if (nontight && noniso) h_nontight_noniso_mbdsigma[etabin]->Fill(clusterET, pileup_result.avgsigma, weight);
            }

            // ---------------------------------------------------------------
            // Task 4: Single vs double interaction purity (sim only)
            // ---------------------------------------------------------------
            if (issim && npb_bdt_ptr && ss_bdt_ptr)
            {
                // Single interaction: requires original vertex in range
                if (vtx_pass)
                {
                    if (tight && iso)    { h_tight_iso_single[etabin]->Fill(clusterET, weight);
                                           if (is_signal) h_tight_iso_single_signal[etabin]->Fill(clusterET, weight);
                                           if (is_signal_inclusive) h_tight_iso_single_incsig[etabin]->Fill(clusterET, weight); }
                    if (tight && noniso) { h_tight_noniso_single[etabin]->Fill(clusterET, weight);
                                           if (is_signal) h_tight_noniso_single_signal[etabin]->Fill(clusterET, weight);
                                           if (is_signal_inclusive) h_tight_noniso_single_incsig[etabin]->Fill(clusterET, weight); }
                    if (nontight && iso) { h_nontight_iso_single[etabin]->Fill(clusterET, weight);
                                           if (is_signal) h_nontight_iso_single_signal[etabin]->Fill(clusterET, weight);
                                           if (is_signal_inclusive) h_nontight_iso_single_incsig[etabin]->Fill(clusterET, weight); }
                    if (nontight && noniso) { h_nontight_noniso_single[etabin]->Fill(clusterET, weight);
                                              if (is_signal) h_nontight_noniso_single_signal[etabin]->Fill(clusterET, weight);
                                              if (is_signal_inclusive) h_nontight_noniso_single_incsig[etabin]->Fill(clusterET, weight); }
                }

                // Double interaction: requires double_vtx in range (not original vertex)
                if (std::abs(double_vtx) <= vertexcut)
                {
                    std::vector<float> x_npb_shifted = buildNpbFeatureVector(icl, double_vtx);
                    float npb_score_shifted = npb_bdt_ptr->Compute(x_npb_shifted)[0];

                    if (npb_score_shifted >= npb_score_cut)
                    {
                        std::vector<float> x_bdt_shifted = buildSsBdtFeatureVector(icl, double_vtx);
                        float bdt_score_shifted = ss_bdt_ptr->Compute(x_bdt_shifted)[0];

                        auto [tight_d, nontight_d] = classifyTightNonTight(icl, clusterET, e11_over_e33, e32_over_e35, bdt_score_shifted);

                        if (tight_d && iso)    { h_tight_iso_double[etabin]->Fill(clusterET, weight);
                                                 if (is_signal) h_tight_iso_double_signal[etabin]->Fill(clusterET, weight);
                                                 if (is_signal_inclusive) h_tight_iso_double_incsig[etabin]->Fill(clusterET, weight); }
                        if (tight_d && noniso) { h_tight_noniso_double[etabin]->Fill(clusterET, weight);
                                                 if (is_signal) h_tight_noniso_double_signal[etabin]->Fill(clusterET, weight);
                                                 if (is_signal_inclusive) h_tight_noniso_double_incsig[etabin]->Fill(clusterET, weight); }
                        if (nontight_d && iso) { h_nontight_iso_double[etabin]->Fill(clusterET, weight);
                                                 if (is_signal) h_nontight_iso_double_signal[etabin]->Fill(clusterET, weight);
                                                 if (is_signal_inclusive) h_nontight_iso_double_incsig[etabin]->Fill(clusterET, weight); }
                        if (nontight_d && noniso) { h_nontight_noniso_double[etabin]->Fill(clusterET, weight);
                                                    if (is_signal) h_nontight_noniso_double_signal[etabin]->Fill(clusterET, weight);
                                                    if (is_signal_inclusive) h_nontight_noniso_double_incsig[etabin]->Fill(clusterET, weight); }
                    }
                }
            }
        } // end cluster loop

        ientry++;
    } // end event loop

    // ---------------------------------------------------------------
    // Post-loop: compute and print simplified purity
    // ---------------------------------------------------------------
    std::cout << "\n========================================" << std::endl;
    std::cout << "Simplified Purity Results" << std::endl;
    std::cout << "========================================" << std::endl;

    // Task 3: Data MBD sigma 2D integral summary
    if (!issim)
    {
        std::cout << "\n--- Task 3: Data MBD sigma 2D histograms (integral summary) ---" << std::endl;
        for (int ieta = 0; ieta < n_eta_bins; ieta++)
        {
            std::cout << Form("Eta bin %d: TightIso=%.0f  TightNonIso=%.0f  NonTightIso=%.0f  NonTightNonIso=%.0f",
                ieta,
                h_tight_iso_mbdsigma[ieta]->Integral(),
                h_tight_noniso_mbdsigma[ieta]->Integral(),
                h_nontight_iso_mbdsigma[ieta]->Integral(),
                h_nontight_noniso_mbdsigma[ieta]->Integral()) << std::endl;
        }
    }

    // Task 4: Sim single vs double purity
    if (issim)
    {
        std::cout << "\n--- Task 4: Sim Single vs Double Interaction Purity ---" << std::endl;
        for (int ieta = 0; ieta < n_eta_bins; ieta++)
        {
            std::cout << Form("\nEta bin %d (%.1f < eta < %.1f):", ieta, eta_bins[ieta], eta_bins[ieta + 1]) << std::endl;
            std::cout << "pT_bin | A_single A_signal purity_single | A_double A_d_signal purity_double" << std::endl;
            for (int ipt = 0; ipt < n_pT_bins; ipt++)
            {
                double A_single = h_tight_iso_single[ieta]->GetBinContent(ipt + 1);
                double A_single_sig = h_tight_iso_single_signal[ieta]->GetBinContent(ipt + 1);
                double purity_single = (A_single > 0) ? A_single_sig / A_single : -1;

                double A_double = h_tight_iso_double[ieta]->GetBinContent(ipt + 1);
                double A_double_sig = h_tight_iso_double_signal[ieta]->GetBinContent(ipt + 1);
                double purity_double = (A_double > 0) ? A_double_sig / A_double : -1;

                std::cout << Form("  %d   | %8.1f %8.1f   %6.3f      | %8.1f  %8.1f    %6.3f",
                    ipt, A_single, A_single_sig, purity_single,
                    A_double, A_double_sig, purity_double) << std::endl;
            }
        }
    }

    // Task 1: NPB 2D integral summary
    std::cout << "\n--- Task 1: NPB Score vs ET 2D histograms (integral summary) ---" << std::endl;
    for (int ieta = 0; ieta < n_eta_bins; ieta++)
    {
        std::cout << Form("Eta bin %d: TightIso=%.0f  TightNonIso=%.0f  NonTightIso=%.0f  NonTightNonIso=%.0f",
            ieta,
            h_tight_iso_npb[ieta]->Integral(),
            h_tight_noniso_npb[ieta]->Integral(),
            h_nontight_iso_npb[ieta]->Integral(),
            h_nontight_noniso_npb[ieta]->Integral()) << std::endl;
    }

    std::cout << "\n========================================" << std::endl;

    // Write and close
    fout->cd();
    fout->Write();
    fout->Close();

    // Cleanup
    if (npb_bdt_ptr) delete npb_bdt_ptr;
    if (ss_bdt_ptr) delete ss_bdt_ptr;

    std::cout << "Done. Output written to " << outfilename << std::endl;
}
