#include <yaml-cpp/yaml.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TFile.h>
#include <TH1.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <map>
#include <memory>

#include "MbdPileupHelper.h"

// EMCAL time sample period used in RecoEffCalculator_TTreeReader.C
const float TIME_SAMPLE_NS = 17.6;

void ShowerShapeCheck(const std::string &configname = "config_showershape.yaml", const std::string filetype = "jet10", bool doinclusive = true)
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

    // Non-physical background (NPB) selection:
    // - apply run-by-run MBD t0 correction (data)
    // - "bad time" is defined by cluster-MBD time: (t_cluster - t_MBD) < npb_delta_t_cut (default: -5 ns)
    int npb_cut_on = configYaml["analysis"]["npb_cut_on"].as<int>(1);
    float npb_weta_min = configYaml["analysis"]["npb_weta_min"].as<float>(0.4);
    float npb_delta_t_cut = configYaml["analysis"]["npb_delta_t_cut"].as<float>(-5.0);
    // Common cut selection on the per-cluster NPB score (default: > 0.5)
    float npb_score_cut = configYaml["analysis"]["npb_score_cut"].as<float>(0.5);

    int mbd_t0_correction_on = configYaml["analysis"]["mbd_t0_correction_on"].as<int>(1);
    std::string mbd_t0_correction_file =
        configYaml["analysis"]["mbd_t0_correction_file"].as<std::string>("/sphenix/user/shuhangli/ppg12/efficiencytool/MbdOut.corr");
    std::map<int, float> mbd_t0_correction;
    if (mbd_t0_correction_on)
    {
        std::ifstream file(mbd_t0_correction_file);
        if (!file.is_open())
        {
            std::cerr << "[MBD t0] WARNING: cannot open correction file: " << mbd_t0_correction_file
                      << " (will proceed with zero correction)" << std::endl;
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
            std::cout << "[MBD t0] Loaded " << mbd_t0_correction.size()
                      << " run-by-run corrections from " << mbd_t0_correction_file << std::endl;
        }
    }

    float max_photon_lower = 0;
    float max_photon_upper = 100;

    // unit in pb
    const float photon5cross = 146359.3;
    const float photon10cross = 6944.675;
    const float photon20cross = 130.4461;

    // Hanpu uses unit in b
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

    float energy_scale_lower = 0;
    float energy_scale_upper = 100;

    float cluster_ET_upper = 100;

    if (filetype == "photon5")
    {
        max_photon_lower = 0;
        max_photon_upper = 14;
        // max_photon_upper = 200;
        weight = photon5cross / photon20cross;
    }
    else if (filetype == "photon10")
    {
        max_photon_lower = 14;
        max_photon_upper = 30;

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
    // Vertex reweighting for simulation (required when enabled):
    //   results/vertex_reweight_bdt_none.root : h_vertexz_ratio_data_over_mccombined
    TH1* h_vertex_reweight = nullptr;
    int vertex_reweight_on = 1;
    std::string vertex_reweight_file = "results/vertex_reweight.root";
    std::cout << "loading vertex reweight file: " << vertex_reweight_file << std::endl;
    if (issim)
    {
        // optional config knobs (safe defaults)
        vertex_reweight_on = configYaml["analysis"]["vertex_reweight_on"].as<int>(1);
        vertex_reweight_file =
            configYaml["analysis"]["vertex_reweight_file"].as<std::string>("results/vertex_reweight.root");
        if (vertex_reweight_on)
        {
            TFile* fvtx = TFile::Open(vertex_reweight_file.c_str(), "READ");

            if (!fvtx || fvtx->IsZombie())
            {
                std::cerr << "[VertexReweight] ERROR: cannot open vertex reweight file: "
                          << vertex_reweight_file << std::endl;
                return;
            }

            TH1* htmp = dynamic_cast<TH1*>(fvtx->Get("h_vertexz_ratio_data_over_mccombined"));
            if (!htmp)
            {
                std::cerr << "[VertexReweight] ERROR: cannot find histogram 'h_vertexz_ratio_data_over_mccombined' in "
                          << vertex_reweight_file << std::endl;
                fvtx->Close();
                delete fvtx;
                return;
            }
            h_vertex_reweight = dynamic_cast<TH1*>(htmp->Clone("h_vertexz_ratio_data_over_mccombined_clone"));
            //fvtx->Close();
            //delete fvtx;



            if (!h_vertex_reweight)
            {
                std::cerr << "[VertexReweight] ERROR: failed to clone histogram from "
                          << vertex_reweight_file << std::endl;
                return;
            }

            h_vertex_reweight->SetDirectory(nullptr);
            h_vertex_reweight->Sumw2();
            std::cout << "[VertexReweight] Using histogram weights from "
                      << vertex_reweight_file << " : " << htmp->GetName() << std::endl;
        }
    }


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

    // Create TChain instead of TTree
    std::string treename = configYaml["input"]["tree"].as<std::string>();
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());

    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();

    int iso_threshold = configYaml["analysis"]["iso_threshold"].as<int>(0);
    int iso_hcalonly = configYaml["analysis"]["iso_hcalonly"].as<int>(0);
    // Inner exclusion radius options: 0, 0.05, 0.1, 0.2
    float iso_emcalinnerr = configYaml["analysis"]["iso_emcalinnerr"].as<float>(0.0);
    std::cout << "iso_emcalinnerr: " << iso_emcalinnerr << std::endl;

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

    // BDT bins with default values if not specified in config
    std::vector<float> bdt_bins_default = {0.0, 0.3, 0.7, 1.0};
    std::vector<float> bdt_bins = configYaml["analysis"]["bdt_bins"].as<std::vector<float>>(bdt_bins_default);
    int n_bdt_bins = bdt_bins.size() - 1;
    std::cout << "n_bdt_bins: " << n_bdt_bins << std::endl;
    for (int i = 0; i < n_bdt_bins + 1; i++)
    {
        std::cout << "bdt_bins[" << i << "]: " << bdt_bins[i] << std::endl;
    }

    float reco_min_ET = configYaml["analysis"]["reco_min_ET"].as<float>();

    float eff_dR = configYaml["analysis"]["eff_dR"].as<float>();

    // trigger_used can be either a scalar int or a YAML sequence of ints.
    // Example: trigger_used: [26, 29, 30, 31, 36, 37, 38]
    std::vector<int> trigger_used;
    {
        YAML::Node trigNode = configYaml["analysis"]["trigger_used"];
        if (trigNode && trigNode.IsSequence())
        {
            trigger_used = trigNode.as<std::vector<int>>();
        }
        else
        {
            // backward-compatible: allow single int
            trigger_used.push_back(configYaml["analysis"]["trigger_used"].as<int>());
        }
    }

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

    float tight_bdt_max = configYaml["analysis"]["tight"]["bdt_max"].as<float>(1.0);
    float tight_bdt_min = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0.0);

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

    float non_tight_bdt_max = configYaml["analysis"]["non_tight"]["bdt_max"].as<float>(1.0);
    float non_tight_bdt_min = configYaml["analysis"]["non_tight"]["bdt_min"].as<float>(0.0);

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
    float common_cluster_weta_cogx_bound = configYaml["analysis"]["common"]["cluster_weta_cogx_bound"].as<float>(0.8);

    int reweight = configYaml["analysis"]["unfold"]["reweight"].as<int>(0); // 0 for no reweighting, 1 for reweighting
    float clusterescale = configYaml["analysis"]["cluster_escale"].as<float>(1.0);

    // Optional per-cluster NPB score branch:
    // try `cluster_npb_score_<node>` first, then fallback to `cluster_npb_score`.
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
            std::cerr << "[NPBScore] WARNING: cannot find branch 'cluster_npb_score_<node>' or 'cluster_npb_score' in input tree. "<<std::endl;
            return;
        }
    }
    else
    {
        std::cout << "[NPBScore] Found branch 'cluster_npb_score_<node>' in input tree. " << std::endl;
        std::cout << "npb_score_branch: " << npb_score_branch << std::endl;
    }

    // TTreeReader setup
    TTreeReader reader(&chain);

    // Basic event variables
    TTreeReaderValue<int> mbdnorthhit(reader, "mbdnorthhit");
    TTreeReaderValue<int> mbdsouthhit(reader, "mbdsouthhit");
    TTreeReaderValue<int> runnumber(reader, "runnumber");
    TTreeReaderValue<float> mbd_time(reader, "mbd_time");
    
    // MBD PMT timing and charge arrays for pileup detection
    TTreeReaderArray<float> mbd_north_time(reader, "mbdnortht");
    TTreeReaderArray<float> mbd_south_time(reader, "mbdsoutht");
    TTreeReaderArray<float> mbd_north_charge(reader, "mbdnorthq");
    TTreeReaderArray<float> mbd_south_charge(reader, "mbdsouthq");
    
    TTreeReaderValue<int> pythiaid(reader, "pythiaid");
    TTreeReaderValue<int> nparticles(reader, "nparticles");
    TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
    TTreeReaderValue<float> vertexz(reader, "vertexz");
    TTreeReaderValue<float> vertexz_truth(reader, "vertexz_truth");
    TTreeReaderValue<float> energy_scale(reader, "energy_scale");
    TTreeReaderArray<Bool_t> scaledtrigger(reader, "scaledtrigger");
    TTreeReaderArray<Bool_t> livetrigger(reader, "livetrigger");

    // Particle arrays
    TTreeReaderArray<float> particle_E(reader, "particle_E");
    TTreeReaderArray<float> particle_Pt(reader, "particle_Pt");
    TTreeReaderArray<float> particle_Eta(reader, "particle_Eta");
    TTreeReaderArray<float> particle_Phi(reader, "particle_Phi");
    TTreeReaderArray<float> particle_truth_iso_02(reader, "particle_truth_iso_02");
    TTreeReaderArray<float> particle_truth_iso_03(reader, "particle_truth_iso_03");
    TTreeReaderArray<float> particle_truth_iso_04(reader, "particle_truth_iso_04");
    TTreeReaderArray<int> particle_pid(reader, "particle_pid");
    TTreeReaderArray<int> particle_trkid(reader, "particle_trkid");
    TTreeReaderArray<int> particle_photonclass(reader, "particle_photonclass");
    TTreeReaderArray<int> particle_converted(reader, "particle_converted");

    // Cluster arrays
    TTreeReaderArray<float> cluster_E(reader, Form("cluster_E_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_prob(reader, Form("cluster_prob_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_CNN_prob(reader, Form("cluster_CNN_prob_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_pid(reader, Form("cluster_pid_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_02(reader, Form("cluster_iso_02_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03(reader, Form("cluster_iso_03_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_04(reader, Form("cluster_iso_04_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e1(reader, Form("cluster_e1_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e2(reader, Form("cluster_e2_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e3(reader, Form("cluster_e3_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e4(reader, Form("cluster_e4_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et1(reader, Form("cluster_et1_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et2(reader, Form("cluster_et2_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et3(reader, Form("cluster_et3_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et4(reader, Form("cluster_et4_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_weta(reader, Form("cluster_weta_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_wphi(reader, Form("cluster_wphi_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_ietacent(reader, Form("cluster_ietacent_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iphicent(reader, Form("cluster_iphicent_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_detamax(reader, Form("cluster_detamax_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_dphimax(reader, Form("cluster_dphimax_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_weta_cogx(reader, Form("cluster_weta_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_wphi_cogx(reader, Form("cluster_wphi_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_weta_cog(reader, Form("cluster_weta_cog_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_wphi_cog(reader, Form("cluster_wphi_cog_%s", clusternodename.c_str()));

    // Cluster energy arrays
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
    TTreeReaderArray<float> cluster_e55(reader, Form("cluster_e55_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e57(reader, Form("cluster_e57_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e75(reader, Form("cluster_e75_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e77(reader, Form("cluster_e77_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w32(reader, Form("cluster_w32_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e32(reader, Form("cluster_e32_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w72(reader, Form("cluster_w72_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e72(reader, Form("cluster_e72_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w52(reader, Form("cluster_w52_%s", clusternodename.c_str()));

    // 2D arrays (flattened in storage)
    TTreeReaderArray<float> cluster_e_array(reader, Form("cluster_e_array_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_time_array(reader, Form("cluster_time_array_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ownership_array(reader, Form("cluster_ownership_array_%s", clusternodename.c_str()));

    // Cluster isolation arrays
    TTreeReaderArray<float> cluster_iso_03_emcal(reader, Form("cluster_iso_03_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_hcalin(reader, Form("cluster_iso_03_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_hcalout(reader, Form("cluster_iso_03_hcalout_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_emcal(reader, Form("cluster_iso_03_60_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_hcalin(reader, Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_hcalout(reader, Form("cluster_iso_03_60_hcalout_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_emcal(reader, Form("cluster_iso_03_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalin(reader, Form("cluster_iso_03_70_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalout(reader, Form("cluster_iso_03_70_hcalout_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_005_70_emcal(reader, Form("cluster_iso_005_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_01_70_emcal(reader, Form("cluster_iso_01_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_02_70_emcal(reader, Form("cluster_iso_02_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_120_emcal(reader, Form("cluster_iso_03_120_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_120_hcalin(reader, Form("cluster_iso_03_120_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_120_hcalout(reader, Form("cluster_iso_03_120_hcalout_%s", clusternodename.c_str()));

    // Cluster HCal arrays
    TTreeReaderArray<float> cluster_ihcal_et(reader, Form("cluster_ihcal_et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_ohcal_et(reader, Form("cluster_ohcal_et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_ihcal_et22(reader, Form("cluster_ihcal_et22_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_ohcal_et22(reader, Form("cluster_ohcal_et22_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_ihcal_et33(reader, Form("cluster_ihcal_et33_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_ohcal_et33(reader, Form("cluster_ohcal_et33_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ihcal_ieta(reader, Form("cluster_ihcal_ieta_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ihcal_iphi(reader, Form("cluster_ihcal_iphi_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ohcal_ieta(reader, Form("cluster_ohcal_ieta_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ohcal_iphi(reader, Form("cluster_ohcal_iphi_%s", clusternodename.c_str()));

    // BDT score - NEW
    std::string bdt_model_name = configYaml["input"]["bdt_model_name"].as<std::string>("base");
    TTreeReaderArray<float> cluster_bdt(reader, Form("cluster_bdt_%s_%s", clusternodename.c_str(), bdt_model_name.c_str()));

    // Optional per-cluster NPB score (may be absent in some trees)
    std::unique_ptr<TTreeReaderArray<float>> cluster_npb_score;
    if (!npb_score_branch.empty())
    {
        cluster_npb_score = std::make_unique<TTreeReaderArray<float>>(reader, npb_score_branch.c_str());
    }

    // Truth jet arrays
    TTreeReaderValue<int> njet_truth(reader, "njet_truth");
    TTreeReaderArray<float> jet_truth_E(reader, "jet_truth_E");
    TTreeReaderArray<float> jet_truth_Pt(reader, "jet_truth_Pt");
    TTreeReaderArray<float> jet_truth_Eta(reader, "jet_truth_Eta");
    TTreeReaderArray<float> jet_truth_Phi(reader, "jet_truth_Phi");

    // Reco jet arrays
    TTreeReaderValue<int> njet(reader, "njet");
    TTreeReaderArray<float> jet_E(reader, "jet_E");
    TTreeReaderArray<float> jet_Pt(reader, "jet_Pt");
    TTreeReaderArray<float> jet_Eta(reader, "jet_Eta");
    TTreeReaderArray<float> jet_Phi(reader, "jet_Phi");

    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    fout->cd();
    // TH2D *
    std::vector<std::vector<TH1D *>> h_tight_reco_iso;
    h_tight_reco_iso.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH1D *>> h_nt_reco_iso;
    h_nt_reco_iso.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH1D *>> h_npb_reco_iso;
    h_npb_reco_iso.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH1D *>> h_signal_truth_iso;
    h_signal_truth_iso.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH1D *>> h_background_truth_iso;
    h_background_truth_iso.resize(eta_bins.size() - 1);
    std::vector<TH2D *> h_ET_isoET;

    // NPB score vs cluster-MBD time histograms (2D)
    std::vector<std::vector<TH2D *>> h_npb_score_vs_time;
    h_npb_score_vs_time.resize(eta_bins.size() - 1);

    // NPB score vs cluster-MBD time histograms with strict pileup selection (MBD pileup detected)
    std::vector<std::vector<TH2D *>> h_npb_score_vs_time_strict;
    h_npb_score_vs_time_strict.resize(eta_bins.size() - 1);

    // NPB score vs cluster-MBD time histograms with clean event selection (avgsigma < 0.5, no pileup)
    std::vector<std::vector<TH2D *>> h_npb_score_vs_time_clean;
    h_npb_score_vs_time_clean.resize(eta_bins.size() - 1);

    // Delta t vs MBD pileup metric (avgsigma) histograms
    std::vector<std::vector<TH2D *>> h_delta_t_vs_mbd_pileup;
    h_delta_t_vs_mbd_pileup.resize(eta_bins.size() - 1);

    // MBD pileup metrics vs cluster-MBD time histograms
    std::vector<std::vector<TH2D *>> h_mbd_avgsigma_vs_time;
    h_mbd_avgsigma_vs_time.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH2D *>> h_mbd_prodsigma_vs_time;
    h_mbd_prodsigma_vs_time.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH2D *>> h_mbd_maxsigma_vs_time;
    h_mbd_maxsigma_vs_time.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH2D *>> h_mbd_avgdelta_vs_time;
    h_mbd_avgdelta_vs_time.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH2D *>> h_mbd_maxdelta_vs_time;
    h_mbd_maxdelta_vs_time.resize(eta_bins.size() - 1);

    // MBD avgsigma vs NPB score histograms
    std::vector<std::vector<TH2D *>> h_mbd_avgsigma_vs_npb_score;
    h_mbd_avgsigma_vs_npb_score.resize(eta_bins.size() - 1);

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

            h_npb_reco_iso[ieta].push_back(new TH1D(Form("h_npb_reco_iso_eta%d_pt%d", ieta, ipt),
                                                    Form("NPB Iso ET %.1f < eta < %.1f, %.1f < pT < %.1f", eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                    400, -10, 30));

            h_npb_score_vs_time[ieta].push_back(new TH2D(Form("h_npb_score_vs_time_eta%d_pt%d", ieta, ipt),
                                                         Form("NPB Score vs Cluster-MBD Time %.1f < eta < %.1f, %.1f < pT < %.1f;Cluster-MBD Time [ns];NPB Score",
                                                              eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                         40, -20, 20, 100, 0.0, 1.0));

            h_npb_score_vs_time_strict[ieta].push_back(new TH2D(Form("h_npb_score_vs_time_strict_eta%d_pt%d", ieta, ipt),
                                                         Form("NPB Score vs Cluster-MBD Time (MBD pileup) %.1f < eta < %.1f, %.1f < pT < %.1f;Cluster-MBD Time [ns];NPB Score",
                                                              eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                         40, -20, 20, 100, 0.0, 1.0));

            h_npb_score_vs_time_clean[ieta].push_back(new TH2D(Form("h_npb_score_vs_time_clean_eta%d_pt%d", ieta, ipt),
                                                         Form("NPB Score vs Cluster-MBD Time (clean, avgsigma<0.5) %.1f < eta < %.1f, %.1f < pT < %.1f;Cluster-MBD Time [ns];NPB Score",
                                                              eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                         40, -20, 20, 100, 0.0, 1.0));

            h_delta_t_vs_mbd_pileup[ieta].push_back(new TH2D(Form("h_delta_t_vs_mbd_pileup_eta%d_pt%d", ieta, ipt),
                                                         Form("Cluster-MBD Time vs MBD Pileup Score %.1f < eta < %.1f, %.1f < pT < %.1f;Cluster-MBD Time [ns];MBD Pileup avgsigma [ns]",
                                                              eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                         80, -20, 20, 100, 0.0, 10.0));

            // MBD pileup metrics vs cluster-MBD time histograms
            h_mbd_avgsigma_vs_time[ieta].push_back(new TH2D(Form("h_mbd_avgsigma_vs_time_eta%d_pt%d", ieta, ipt),
                                                           Form("MBD Avg RMS Time vs Cluster-MBD Time %.1f < eta < %.1f, %.1f < pT < %.1f;Cluster-MBD Time [ns];MBD Avg #sigma_{t} [ns]",
                                                                eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                           80, -20, 20, 100, 0.0, 10.0));

            h_mbd_prodsigma_vs_time[ieta].push_back(new TH2D(Form("h_mbd_prodsigma_vs_time_eta%d_pt%d", ieta, ipt),
                                                            Form("MBD Prod RMS Time vs Cluster-MBD Time %.1f < eta < %.1f, %.1f < pT < %.1f;Cluster-MBD Time [ns];MBD #sigma_{N} #times #sigma_{S} [ns^{2}]",
                                                                 eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                            80, -20, 20, 100, 0.0, 20.0));

            h_mbd_maxsigma_vs_time[ieta].push_back(new TH2D(Form("h_mbd_maxsigma_vs_time_eta%d_pt%d", ieta, ipt),
                                                           Form("MBD Max RMS Time vs Cluster-MBD Time %.1f < eta < %.1f, %.1f < pT < %.1f;Cluster-MBD Time [ns];MBD Max #sigma_{t} [ns]",
                                                                eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                           80, -20, 20, 100, 0.0, 10.0));

            h_mbd_avgdelta_vs_time[ieta].push_back(new TH2D(Form("h_mbd_avgdelta_vs_time_eta%d_pt%d", ieta, ipt),
                                                           Form("MBD Avg #Deltat vs Cluster-MBD Time %.1f < eta < %.1f, %.1f < pT < %.1f;Cluster-MBD Time [ns];MBD Avg #Deltat [ns]",
                                                                eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                           80, -20, 20, 100, 0.0, 20.0));

            h_mbd_maxdelta_vs_time[ieta].push_back(new TH2D(Form("h_mbd_maxdelta_vs_time_eta%d_pt%d", ieta, ipt),
                                                           Form("MBD Max #Deltat vs Cluster-MBD Time %.1f < eta < %.1f, %.1f < pT < %.1f;Cluster-MBD Time [ns];MBD Max #Deltat [ns]",
                                                                eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                           80, -20, 20, 100, 0.0, 20.0));

            h_mbd_avgsigma_vs_npb_score[ieta].push_back(new TH2D(Form("h_mbd_avgsigma_vs_npb_score_eta%d_pt%d", ieta, ipt),
                                                           Form("MBD Avg #sigma_{t} vs NPB Score %.1f < eta < %.1f, %.1f < pT < %.1f;NPB Score;MBD Avg #sigma_{t} [ns]",
                                                                eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                           100, 0.0, 1.0, 100, 0.0, 10.0));
        }
    }

    // Histograms for BDT bins (reco iso ET)
    // Merge tight/non-tight into one: BDT score already encodes shower-shape information.
    // We fill this histogram for clusters passing either tight or non-tight region.
    std::vector<std::vector<std::vector<TH1D *>>> h_reco_iso_bdt;
    h_reco_iso_bdt.resize(eta_bins.size() - 1);

    for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
    {
        h_reco_iso_bdt[ieta].resize(n_pT_bins);

        for (int ipt = 0; ipt < n_pT_bins; ipt++)
        {
            h_reco_iso_bdt[ieta][ipt].resize(n_bdt_bins);

            for (int ibdt = 0; ibdt < n_bdt_bins; ibdt++)
            {
                h_reco_iso_bdt[ieta][ipt][ibdt] = new TH1D(
                    Form("h_reco_iso_eta%d_pt%d_bdt%d", ieta, ipt, ibdt),
                    Form("Reco Iso ET %.1f < eta < %.1f, %.1f < pT < %.1f, %.2f < BDT < %.2f",
                         eta_bins[ieta], eta_bins[ieta + 1],
                         pT_bin_edges[ipt], pT_bin_edges[ipt + 1],
                         bdt_bins[ibdt], bdt_bins[ibdt + 1]),
                    400, -10, 30);
            }
        }
    }

    // BDT-binned shower shape histograms
    // Structure: [varname][eta][pt][bdt] - no cut dimension since BDT range is already the ID
    std::map<std::string, std::vector<std::vector<std::vector<TH1D *>>>> h_showershapes_bdt;

    // Define key shower shape variables to histogram
    std::vector<std::string> ss_vars = {
        "weta_cogx",
        "wphi_cogx",
        "et1",
        "et2",
        "et3",
        "et4",
        "prob",
        "e11_to_e33",
        "e32_to_e35",
        "npb_score"
    };

    // Binning for each variable [nbins, min, max]
    std::map<std::string, std::vector<double>> ss_binning = {
        {"weta_cogx", {100, 0, 2.0}},
        {"wphi_cogx", {100, 0, 2.0}},
        {"et1", {100, 0, 1.0}},
        {"et2", {100, 0, 1.0}},
        {"et3", {100, 0, 1.0}},
        {"et4", {100, 0, 0.5}},
        {"prob", {100, 0, 1.0}},
        {"e11_to_e33", {100, 0, 1.0}},
        {"e32_to_e35", {100, 0.5, 1.0}},
        {"npb_score", {100, 0.0, 1.0}}
    };

    for (const auto &varname : ss_vars)
    {
        h_showershapes_bdt[varname].resize(eta_bins.size() - 1);

        for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
        {
            h_showershapes_bdt[varname][ieta].resize(n_pT_bins);

            for (int ipt = 0; ipt < n_pT_bins; ipt++)
            {
                h_showershapes_bdt[varname][ieta][ipt].resize(n_bdt_bins);

                for (int ibdt = 0; ibdt < n_bdt_bins; ibdt++)
                {
                    auto &binning = ss_binning[varname];
                    h_showershapes_bdt[varname][ieta][ipt][ibdt] = new TH1D(
                        Form("h_%s_eta%d_pt%d_bdt%d", varname.c_str(), ieta, ipt, ibdt),
                        Form("%s %.1f < eta < %.1f, %.1f < pT < %.1f, %.2f < BDT < %.2f",
                             varname.c_str(),
                             eta_bins[ieta], eta_bins[ieta + 1],
                             pT_bin_edges[ipt], pT_bin_edges[ipt + 1],
                             bdt_bins[ibdt], bdt_bins[ibdt + 1]),
                        (int)binning[0], binning[1], binning[2]);
                }
            }
        }
    }

    // shower shape correlation with isoET before prelim cut
    // a master vector of 3D histograms
    std::map<std::string, std::vector<std::vector<std::vector<TH2F *>>>> h2d_all;

    // 2) Define your "histogram names" (the old map keys become string keys).
    //    For instance, you can directly populate them:

    static const int ncuts = 5;  // all, common, tight, non-tight, npb

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
    h2d_all["h2d_bdt"];
    h2d_all["h2d_npb_score"];

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
                    // Default axis ranges (most shower-shape vars)
                    int nx = 300;
                    double xmin = -1.0;
                    double xmax = 2.0;
                    int ny = 200;
                    double ymin = -10.0;
                    double ymax = 30.0;

                    // Special-case: NPB score is naturally in [0, 1]
                    if (basename == "h2d_npb_score")
                    {
                        nx = 100;
                        xmin = 0.0;
                        xmax = 1.0;
                    }

                    TH2F *h2D = new TH2F(
                        Form("%s_eta%d_pt%d_cut%d", basename.c_str(), ieta, ipt, icut),
                        Form("%s vs Iso ET; %s; Iso ET [GeV]",
                             basename.c_str(), basename.c_str()),
                        nx, xmin, xmax,
                        ny, ymin, ymax);
                    h2D->SetDirectory(fout);
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

    int nentries = chain.GetEntries();
    int ientry = 0;
    while (reader.Next())
    {
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        if (!issim)
        {
            // Accept event if ANY trigger in trigger_used fired.
            bool any_trigger_fired = false;
            const auto ntrig = scaledtrigger.GetSize();
            for (int itrig : trigger_used)
            {
                if (itrig < 0 || (unsigned int)itrig >= ntrig)
                {
                    continue;
                }
                if (scaledtrigger[itrig] != 0)
                {
                    any_trigger_fired = true;
                    break;
                }
            }

            if (!any_trigger_fired)
            {
                continue;
            }
        }

        // Apply vertex reweighting for MC
        weight = cross_weight;
        vertex_weight = 1.0;
        if (issim && vertex_reweight_on)
        {
            if (!h_vertex_reweight)
            {
                std::cerr << "[VertexReweight] ERROR: vertex reweighting is enabled but histogram is not loaded."
                          << std::endl;
                return;
            }

            int bin = h_vertex_reweight->FindBin(*vertexz);
            if (bin < 1) bin = 1;
            if (bin > h_vertex_reweight->GetNbinsX()) bin = h_vertex_reweight->GetNbinsX();
            vertex_weight = h_vertex_reweight->GetBinContent(bin);

            if (!std::isfinite(vertex_weight) || vertex_weight <= 0.0)
            {
                vertex_weight = 1.0;
            }
            weight *= vertex_weight;
        }

        // Calculate MBD pileup metrics for strict pileup selection
        MbdPileupResult pileup_result = calculateMbdPileupMetrics(
            mbd_north_time, mbd_south_time, mbd_north_charge, mbd_south_charge);
        bool is_mbd_pileup = isMbdPileup(pileup_result, PileupCutStrength::STRICT);

        std::set<int> signal_set;
        std::set<int> background_set;
        std::map<int, int> particle_trkidmap;
        if (issim)
        {
            float maxphotonpT = 0;
            int maxphotonclass = 0;

            for (int iparticle = 0; iparticle < *nparticles; iparticle++)
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
                for (int ijet = 0; ijet < *njet_truth; ijet++)
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
                 if ((*energy_scale > energy_scale_upper) || (*energy_scale < energy_scale_lower))
                 {
                     continue;
                 }
                 */
            }
        }

        // loop over clusters
        for (int icluster = 0; icluster < *ncluster; icluster++)
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
            float rhad33 = (cluster_ihcal_et33[icluster] + cluster_ohcal_et33[icluster]) / (cluster_Et[icluster] + (cluster_ihcal_et33[icluster] + cluster_ohcal_et33[icluster]));

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
                if (iso_hcalonly)
                {
                    recoisoET = cluster_iso_03_70_hcalin[icluster] + cluster_iso_03_70_hcalout[icluster];
                }
                else
                {
                    recoisoET = cluster_iso_03_70_emcal[icluster] + cluster_iso_03_70_hcalin[icluster] + cluster_iso_03_70_hcalout[icluster];

                    // Apply inner exclusion based on iso_emcalinnerr
                    if (iso_emcalinnerr > 0.04 && iso_emcalinnerr < 0.06)
                    {
                        recoisoET -= cluster_iso_005_70_emcal[icluster];
                    }
                    else if (iso_emcalinnerr > 0.09 && iso_emcalinnerr < 0.11)
                    {
                        recoisoET -= cluster_iso_01_70_emcal[icluster];
                    }
                    else if (iso_emcalinnerr > 0.19 && iso_emcalinnerr < 0.21)
                    {
                        recoisoET -= cluster_iso_02_70_emcal[icluster];
                    }
                }
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

            if (abs(*vertexz) > vertexcut)
                continue;

            if (!(*mbdnorthhit >= 1 && *mbdsouthhit >= 1))
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

            // Determine BDT bin
            int bdtbin = -1;
            float bdt_score = cluster_bdt[icluster];
            for (int ibdt = 0; ibdt < n_bdt_bins; ibdt++)
            {
                if (bdt_score > bdt_bins[ibdt] && bdt_score < bdt_bins[ibdt + 1])
                {
                    bdtbin = ibdt;
                    break;
                }
            }

            const float npb_score_val = (cluster_npb_score ? (*cluster_npb_score)[icluster] : 1.0f);

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
                //std::cout<<cluster_e1[icl]<<" "<<recoisoET<<" "<<weight<<std::endl;
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
                h2d_all["h2d_bdt"][idx][etabin][pTbin]->Fill(cluster_bdt[icl], recoisoET, weight);
                h2d_all["h2d_npb_score"][idx][etabin][pTbin]->Fill(npb_score_val, recoisoET, weight);

                h_ET_eta[idx][etabin]->Fill(cluster_Et[icl], cluster_Eta[icl], weight);
            };

            // Helper lambda to fill BDT-binned shower shape histograms
            auto fillShowerShapesBDT = [&](size_t icl)
            {
                if (bdtbin == -1) return;

                h_showershapes_bdt["weta_cogx"][etabin][pTbin][bdtbin]->Fill(cluster_weta_cogx[icl], weight);
                h_showershapes_bdt["wphi_cogx"][etabin][pTbin][bdtbin]->Fill(cluster_wphi_cogx[icl], weight);
                h_showershapes_bdt["et1"][etabin][pTbin][bdtbin]->Fill(cluster_et1[icl], weight);
                h_showershapes_bdt["et2"][etabin][pTbin][bdtbin]->Fill(cluster_et2[icl], weight);
                h_showershapes_bdt["et3"][etabin][pTbin][bdtbin]->Fill(cluster_et3[icl], weight);
                h_showershapes_bdt["et4"][etabin][pTbin][bdtbin]->Fill(cluster_et4[icl], weight);
                h_showershapes_bdt["prob"][etabin][pTbin][bdtbin]->Fill(cluster_prob[icl], weight);
                h_showershapes_bdt["e11_to_e33"][etabin][pTbin][bdtbin]->Fill(e11_over_e33, weight);
                h_showershapes_bdt["e32_to_e35"][etabin][pTbin][bdtbin]->Fill(e32_over_e35, weight);
                h_showershapes_bdt["npb_score"][etabin][pTbin][bdtbin]->Fill(npb_score_val, weight);
            };

            fillAllHists(0, icluster);
            fillShowerShapesBDT(icluster);

            // calculate cluster time
            float clusteravgtime = 0;
            float cluster_total_e = 0;
            for (int i = 0; i < 49; i++)
            {
                int idx = icluster * 49 + i;  // Flatten 2D array access
                if (cluster_ownership_array[idx] == 1)
                {
                    clusteravgtime += cluster_time_array[idx] * cluster_e_array[idx];
                    // std::cout<<"cluster_time_array[idx]: "<<cluster_time_array[idx]<<std::endl;
                    cluster_total_e += cluster_e_array[idx];
                }
            }
            clusteravgtime = cluster_total_e > 0 ? clusteravgtime / cluster_total_e : 0;
            // Match RecoEffCalculator_TTreeReader.C convention (convert to ns)
            clusteravgtime = clusteravgtime * TIME_SAMPLE_NS;

            // Run-by-run MBD time calibration (match RecoEffCalculator_TTreeReader.C)
            float mbd_mean_time = *mbd_time;
            float mbdoffset = 0.0;
            if (!issim && mbd_t0_correction_on)
            {
                auto it = mbd_t0_correction.find(*runnumber);
                if (it != mbd_t0_correction.end())
                {
                    mbdoffset = it->second;
                }
            }
            // For MC keep mbdoffset=0.0
            mbd_mean_time = mbd_mean_time - mbdoffset;

            // Requested NPB "bad time": cluster-MBD time < npb_delta_t_cut (default: -5 ns)
            const float delta_t_cluster_mbd = clusteravgtime - mbd_mean_time;
            bool badtime = (delta_t_cluster_mbd < npb_delta_t_cut);

            // Fill NPB score vs cluster-MBD time histogram
            h_npb_score_vs_time[etabin][pTbin]->Fill(delta_t_cluster_mbd, npb_score_val, weight);

            // Fill delta t vs MBD pileup metric (avgsigma) histogram
            if (pileup_result.valid)
            {
                h_delta_t_vs_mbd_pileup[etabin][pTbin]->Fill(delta_t_cluster_mbd, pileup_result.avgsigma, weight);

                // Fill all MBD pileup metrics vs cluster-MBD time histograms
                h_mbd_avgsigma_vs_time[etabin][pTbin]->Fill(delta_t_cluster_mbd, pileup_result.avgsigma, weight);
                h_mbd_prodsigma_vs_time[etabin][pTbin]->Fill(delta_t_cluster_mbd, pileup_result.prodsigma, weight);
                h_mbd_maxsigma_vs_time[etabin][pTbin]->Fill(delta_t_cluster_mbd, pileup_result.maxsigma, weight);
                h_mbd_avgdelta_vs_time[etabin][pTbin]->Fill(delta_t_cluster_mbd, pileup_result.avgdelta, weight);
                h_mbd_maxdelta_vs_time[etabin][pTbin]->Fill(delta_t_cluster_mbd, pileup_result.maxdelta, weight);

                // Fill MBD avgsigma vs NPB score histogram
                h_mbd_avgsigma_vs_npb_score[etabin][pTbin]->Fill(npb_score_val, pileup_result.avgsigma, weight);
            }

            // Check for back-to-back jets (veto if present)
            bool otherside_jet = false;
            for (int ijet = 0; ijet < *njet; ijet++)
            {
                if (jet_Pt[ijet] < 5) continue;  // Require jet pT > 5 GeV

                float dphi = cluster_Phi[icluster] - jet_Phi[ijet];
                while (dphi > M_PI) dphi -= 2 * M_PI;
                while (dphi < -M_PI) dphi += 2 * M_PI;

                if (std::abs(dphi) > (M_PI / 2))  // Back-to-back: Δφ > 90°
                {
                    otherside_jet = true;
                    break;
                }
            }

            // Fill NPB score vs cluster-MBD time histogram with strict pileup selection (MBD pileup detected)
            if (is_mbd_pileup)
            {
                h_npb_score_vs_time_strict[etabin][pTbin]->Fill(delta_t_cluster_mbd, npb_score_val, weight);
            }

            // Fill NPB score vs cluster-MBD time histogram with clean event selection (avgsigma < 0.5)
            if (pileup_result.valid && pileup_result.avgsigma < 0.5)
            {
                h_npb_score_vs_time_clean[etabin][pTbin]->Fill(delta_t_cluster_mbd, npb_score_val, weight);
            }

            // NPB tag
            bool isnpb = badtime && !otherside_jet && (cluster_weta_cogx[icluster] >= npb_weta_min);

            if (isnpb)
            {
                fillAllHists(4, icluster);
                h_npb_reco_iso[etabin][pTbin]->Fill(recoisoET, weight);
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
                // require NPB score for the common cut
                (npb_score_val > npb_score_cut) &&
                //(!(wr_cogx < common_wr_cogx_bound && cluster_weta_cogx[icluster] > common_cluster_weta_cogx_bound))
                (cluster_weta_cogx[icluster] < common_cluster_weta_cogx_bound))
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

            bool is_bdt_tight =
                (cluster_bdt[icluster] > tight_bdt_min) &&
                (cluster_bdt[icluster] < tight_bdt_max);

            // Combined condition
            if (is_cluster_weta_cogx_tight &&
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
                cluster_et4[icluster] < non_tight_et4_max &&
                cluster_bdt[icluster] > non_tight_bdt_min &&
                cluster_bdt[icluster] < non_tight_bdt_max)
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
                if (!is_bdt_tight)
                {
                    nfail += bdt_on;
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
                    if (bdt_fail)
                    {
                        if (is_bdt_tight)
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

            // Fill BDT-binned reco iso histograms (merged tight + non-tight)
            if (bdtbin != -1 && (tight || nontight))
            {
                h_reco_iso_bdt[etabin][pTbin][bdtbin]->Fill(recoisoET, weight);
            }
        }
        ientry++;
    }
    fout->cd();
    fout->Write();
    fout->Close();
}
