#include <yaml-cpp/yaml.h>
#include "CrossSectionWeights.h"
using namespace PPG12;
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TObjString.h>
#include <cmath>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <algorithm>

void SaveYamlToRoot(TFile *f, const char *yaml_filename)
{
    std::ifstream yaml_file(yaml_filename);
    std::string yaml_content((std::istreambuf_iterator<char>(yaml_file)),
                             std::istreambuf_iterator<char>());
    TObjString yaml_obj(yaml_content.c_str());
    f->cd();
    yaml_obj.Write("config");
}

void EtaMigrationStudy(
    const std::string &configname = "config_bdt_nom.yaml",
    const std::string &filetype = "photon10",
    bool do_double_interaction = false)
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    // This macro is sim-only (photon and jet MC samples)
    bool issim = true;
    const bool isbackground = filetype.find("jet") != std::string::npos;

    // ---------------------------------------------------------------
    // Input file setup
    // ---------------------------------------------------------------
    std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();

    std::string input_filetype = filetype;
    if (filetype == "photon10_nom") input_filetype = "photon10";
    if (filetype == "jet12_nom")    input_filetype = "jet12";

    std::string infilename = infilename_root_dir + input_filetype + infilename_branch_dir;
    std::cout << "infilename: " << infilename << std::endl;

    // Sample kinematic windows and cross-section weight
    PPG12::SampleConfig sc = PPG12::GetSampleConfig(filetype);
    if (!sc.valid)
    {
        std::cerr << "ERROR: unrecognized filetype '" << filetype << "'" << std::endl;
        return;
    }
    float max_photon_lower = sc.photon_pt_lower;
    float max_photon_upper = sc.photon_pt_upper;
    float max_jet_lower    = sc.jet_pt_lower;
    float max_jet_upper    = sc.jet_pt_upper;
    float cluster_ET_upper = sc.cluster_ET_upper;
    float cross_weight     = sc.weight;

    // ---------------------------------------------------------------
    // Config parsing
    // ---------------------------------------------------------------
    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();
    std::string bdt_model_name = configYaml["input"]["bdt_model_name"].as<std::string>("base");

    std::vector<float> bdt_et_bin_edges;
    std::vector<std::string> bdt_et_bin_models;
    bool use_et_binned_bdt = false;
    if (configYaml["input"]["bdt_et_bin_edges"] && configYaml["input"]["bdt_et_bin_models"])
    {
        for (auto v : configYaml["input"]["bdt_et_bin_edges"])
            bdt_et_bin_edges.push_back(v.as<float>());
        for (auto v : configYaml["input"]["bdt_et_bin_models"])
            bdt_et_bin_models.push_back(v.as<std::string>());
        use_et_binned_bdt = (bdt_et_bin_models.size() == bdt_et_bin_edges.size() - 1);
        if (!use_et_binned_bdt)
            std::cout << "WARNING: bdt_et_bin_edges/bdt_et_bin_models size mismatch; falling back to single model" << std::endl;
    }

    int use_topo_iso = configYaml["analysis"]["use_topo_iso"].as<int>(0);
    int conesize = configYaml["analysis"]["cone_size"].as<int>();

    float truthisocut = configYaml["analysis"]["truth_iso_max"].as<float>();

    float recoiso_min = configYaml["analysis"]["reco_iso_min"].as<float>();
    float recoiso_max_b = configYaml["analysis"]["reco_iso_max_b"].as<float>();
    float recoiso_max_s = configYaml["analysis"]["reco_iso_max_s"].as<float>();

    float recononiso_min_shift = configYaml["analysis"]["reco_noniso_min_shift"].as<float>();
    float recononiso_max = configYaml["analysis"]["reco_noniso_max"].as<float>();

    float vertexcut = configYaml["analysis"]["vertex_cut"].as<float>();
    std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();
    float eta_lo = eta_bins.front();
    float eta_hi = eta_bins.back();

    std::vector<float> pT_bins = configYaml["analysis"]["pT_bins"].as<std::vector<float>>();
    int n_pT_bins = pT_bins.size() - 1;
    double pTmin = pT_bins[0];
    double pTmax = pT_bins[n_pT_bins];

    std::vector<float> pT_bins_truth = configYaml["analysis"]["pT_bins_truth"].as<std::vector<float>>();
    int n_pT_bins_truth = pT_bins_truth.size() - 1;
    double pTmin_truth = pT_bins_truth[0];
    double pTmax_truth = pT_bins_truth[n_pT_bins_truth];

    float reco_min_ET = configYaml["analysis"]["reco_min_ET"].as<float>();
    float mc_iso_shift = configYaml["analysis"]["mc_iso_shift"].as<float>(0.0);
    float mc_iso_scale = configYaml["analysis"]["mc_iso_scale"].as<float>(1.0);
    float clusterescale = configYaml["analysis"]["cluster_escale"].as<float>(1.0);
    float clustereres = configYaml["analysis"]["cluster_eres"].as<float>(0.0);
    int nosat = configYaml["analysis"]["nosat"].as<int>(0);

    // Tight cuts
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

    float tight_bdt_max = configYaml["analysis"]["tight"]["bdt_max"].as<float>(1);
    float tight_bdt_min = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0);
    float tight_bdt_min_slope = configYaml["analysis"]["tight"]["bdt_min_slope"].as<float>(0);
    float tight_bdt_min_intercept = configYaml["analysis"]["tight"]["bdt_min_intercept"].as<float>(tight_bdt_min);

    // Non-tight cuts
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

    float non_tight_bdt_max = configYaml["analysis"]["non_tight"]["bdt_max"].as<float>(1);
    float non_tight_bdt_min = configYaml["analysis"]["non_tight"]["bdt_min"].as<float>(0);
    float non_tight_bdt_max_slope = configYaml["analysis"]["non_tight"]["bdt_max_slope"].as<float>(0);
    float non_tight_bdt_max_intercept = configYaml["analysis"]["non_tight"]["bdt_max_intercept"].as<float>(non_tight_bdt_max);

    int n_nt_fail = configYaml["analysis"]["n_nt_fail"].as<int>(1);
    int weta_fail = configYaml["analysis"]["weta_fail"].as<int>(0);
    int wphi_fail = configYaml["analysis"]["wphi_fail"].as<int>(0);
    int et1_fail = configYaml["analysis"]["et1_fail"].as<int>(0);
    int e11_to_e33_fail = configYaml["analysis"]["e11_to_e33_fail"].as<int>(0);
    int e32_to_e35_fail = configYaml["analysis"]["e32_to_e35_fail"].as<int>(0);
    int bdt_fail = configYaml["analysis"]["bdt_fail"].as<int>(0);

    int weta_on = configYaml["analysis"]["weta_on"].as<int>(1);
    int wphi_on = configYaml["analysis"]["wphi_on"].as<int>(1);
    int et1_on = configYaml["analysis"]["et1_on"].as<int>(1);
    int e11_to_e33_on = configYaml["analysis"]["e11_to_e33_on"].as<int>(1);
    int e32_to_e35_on = configYaml["analysis"]["e32_to_e35_on"].as<int>(1);
    int et2_on = configYaml["analysis"]["et2_on"].as<int>(1);
    int et3_on = configYaml["analysis"]["et3_on"].as<int>(1);
    int et4_on = configYaml["analysis"]["et4_on"].as<int>(1);
    int bdt_on = configYaml["analysis"]["bdt_on"].as<int>(1);

    // Common cuts
    float common_prob_max = configYaml["analysis"]["common"]["prob_max"].as<float>();
    float common_prob_min = configYaml["analysis"]["common"]["prob_min"].as<float>();
    float common_e11_over_e33_max = configYaml["analysis"]["common"]["e11_over_e33_max"].as<float>();
    float common_e11_over_e33_min = configYaml["analysis"]["common"]["e11_over_e33_min"].as<float>();
    float common_wr_cogx_bound = configYaml["analysis"]["common"]["wr_cogx_bound"].as<float>();
    float common_cluster_weta_cogx_bound = configYaml["analysis"]["common"]["cluster_weta_cogx_bound"].as<float>();
    int common_npb_cut_on = configYaml["analysis"]["common"]["npb_cut_on"].as<int>(0);
    float common_npb_score_cut = configYaml["analysis"]["common"]["npb_score_cut"].as<float>(0.5);

    // ---------------------------------------------------------------
    // Output file
    // ---------------------------------------------------------------
    std::string var_type = configYaml["output"]["var_type"].as<std::string>();
    std::string outfilename = "results/eta_migration_" + std::string(do_double_interaction ? "double_" : "") + filetype + "_" + var_type + ".root";
    std::cout << "outfilename: " << outfilename << std::endl;

    // ---------------------------------------------------------------
    // Vertex reweighting
    // ---------------------------------------------------------------
    TH1 *h_vertex_reweight = nullptr;
    int vertex_reweight_on = configYaml["analysis"]["vertex_reweight_on"].as<int>(1);
    std::string vertex_reweight_file =
        configYaml["analysis"]["vertex_reweight_file"].as<std::string>("results/vertex_reweight.root");

    if (vertex_reweight_on)
    {
        TFile *fvtx = TFile::Open(vertex_reweight_file.c_str(), "READ");
        if (fvtx && !fvtx->IsZombie())
        {
            TH1 *htmp = dynamic_cast<TH1 *>(fvtx->Get("h_vertexz_ratio_data_over_mccombined"));
            if (htmp)
            {
                h_vertex_reweight = dynamic_cast<TH1 *>(htmp->Clone("h_vertexz_ratio_clone"));
                h_vertex_reweight->SetDirectory(nullptr);
                h_vertex_reweight->Sumw2();
                std::cout << "[VertexReweight] Loaded from " << vertex_reweight_file << std::endl;
            }
            else
            {
                std::cerr << "[VertexReweight] WARNING: histogram not found in " << vertex_reweight_file << std::endl;
                vertex_reweight_on = 0;
            }
            fvtx->Close();
            delete fvtx;
        }
        else
        {
            std::cerr << "[VertexReweight] WARNING: cannot open " << vertex_reweight_file << ", disabling vertex reweighting" << std::endl;
            vertex_reweight_on = 0;
        }
    }

    // ---------------------------------------------------------------
    // Double interaction: load data vertex distribution
    // ---------------------------------------------------------------
    TH1 *h_vertexz_data_nocut = nullptr;
    if (do_double_interaction)
    {
        std::string vertex_scan_data_file =
            configYaml["analysis"]["vertex_scan_data_file"].as<std::string>("");

        if (!vertex_scan_data_file.empty())
        {
            TFile *fvtxscan = TFile::Open(vertex_scan_data_file.c_str(), "READ");
            if (fvtxscan && !fvtxscan->IsZombie())
            {
                TH1 *htmp = dynamic_cast<TH1 *>(fvtxscan->Get("h_vertexz"));
                if (htmp)
                {
                    h_vertexz_data_nocut = dynamic_cast<TH1 *>(htmp->Clone("h_vertexz_data_nocut_clone"));
                    h_vertexz_data_nocut->SetDirectory(nullptr);
                    std::cout << "[VtxScan] Loaded unrestricted vertex distribution ("
                              << h_vertexz_data_nocut->GetEntries() << " entries)" << std::endl;
                }
                else
                {
                    std::cerr << "[VtxScan] ERROR: h_vertexz not found in " << vertex_scan_data_file << std::endl;
                }
                fvtxscan->Close();
                delete fvtxscan;
            }
            else
            {
                std::cerr << "[VtxScan] ERROR: cannot open " << vertex_scan_data_file << std::endl;
            }
        }

        // Fallback: try vertex_reweight_file for h_vertexz_data
        if (!h_vertexz_data_nocut)
        {
            TFile *fvtx = TFile::Open(vertex_reweight_file.c_str(), "READ");
            if (fvtx && !fvtx->IsZombie())
            {
                TH1 *htmp = dynamic_cast<TH1 *>(fvtx->Get("h_vertexz_data"));
                if (htmp)
                {
                    h_vertexz_data_nocut = dynamic_cast<TH1 *>(htmp->Clone("h_vertexz_data_fallback_clone"));
                    h_vertexz_data_nocut->SetDirectory(nullptr);
                    std::cout << "[VtxScan] Fallback: loaded h_vertexz_data from " << vertex_reweight_file
                              << " (" << h_vertexz_data_nocut->GetEntries() << " entries)" << std::endl;
                }
                else
                {
                    std::cerr << "[VtxScan] ERROR: h_vertexz_data not found in " << vertex_reweight_file << std::endl;
                    std::cerr << "[VtxScan] Double interaction mode disabled (no vertex distribution)" << std::endl;
                    do_double_interaction = false;
                }
                fvtx->Close();
                delete fvtx;
            }
            else
            {
                std::cerr << "[VtxScan] Cannot open " << vertex_reweight_file << ", disabling double interaction" << std::endl;
                do_double_interaction = false;
            }
        }
    }

    // ---------------------------------------------------------------
    // TChain setup
    // ---------------------------------------------------------------
    std::string treename = configYaml["input"]["tree"].as<std::string>();
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());
    int nentries = chain.GetEntries();
    std::cout << "nentries: " << nentries << std::endl;
    if (nentries == 0)
    {
        std::cerr << "ERROR: no entries found in " << infilename << std::endl;
        return;
    }

    // ---------------------------------------------------------------
    // TTreeReader setup
    // ---------------------------------------------------------------
    TTreeReader reader(&chain);

    TTreeReaderValue<int> nparticles(reader, "nparticles");
    TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
    TTreeReaderValue<float> vertexz(reader, "vertexz");

    // Double-interaction MC uses the "_AntiKt_Truth_r04" jet branch naming.
    const bool use_r04_branches = (filetype == "jet12_double" || filetype == "photon10_double");
    const std::string njet_truth_bname   = use_r04_branches ? "njet_truth_AntiKt_Truth_r04"   : "njet_truth";
    const std::string jet_truth_Pt_bname = use_r04_branches ? "jet_truth_Pt_AntiKt_Truth_r04" : "jet_truth_Pt";

    // Truth jet arrays (used for jet pT window cut on background samples)
    std::unique_ptr<TTreeReaderValue<int>>      njet_truth_ptr;
    std::unique_ptr<TTreeReaderArray<float>>    jet_truth_Pt_ptr;
    if (isbackground)
    {
        njet_truth_ptr.reset(new TTreeReaderValue<int>(reader, njet_truth_bname.c_str()));
        jet_truth_Pt_ptr.reset(new TTreeReaderArray<float>(reader, jet_truth_Pt_bname.c_str()));
    }

    // Particle arrays
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
    TTreeReaderArray<int> cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_nsaturated(reader, Form("cluster_nsaturated_%s", clusternodename.c_str()));

    // Isolation
    TTreeReaderArray<float> cluster_iso_02(reader, Form("cluster_iso_02_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03(reader, Form("cluster_iso_03_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_04(reader, Form("cluster_iso_04_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_excl_04(reader, Form("cluster_iso_excl_04_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_topo_03(reader, Form("cluster_iso_topo_03_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_topo_04(reader, Form("cluster_iso_topo_04_%s", clusternodename.c_str()));

    // Shower shape variables
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

    // NPB score
    TTreeReaderArray<float> cluster_npb_score(reader, Form("cluster_npb_score_%s", clusternodename.c_str()));

    // BDT score arrays
    std::vector<std::string> all_bdt_models = {bdt_model_name};
    if (use_et_binned_bdt)
        all_bdt_models.insert(all_bdt_models.end(), bdt_et_bin_models.begin(), bdt_et_bin_models.end());
    std::sort(all_bdt_models.begin(), all_bdt_models.end());
    all_bdt_models.erase(std::unique(all_bdt_models.begin(), all_bdt_models.end()), all_bdt_models.end());

    std::map<std::string, TTreeReaderArray<float> *> bdt_arrays;
    for (auto &mname : all_bdt_models)
    {
        bdt_arrays[mname] = new TTreeReaderArray<float>(reader,
            Form("cluster_bdt_%s_%s", clusternodename.c_str(), mname.c_str()));
    }

    // ---------------------------------------------------------------
    // Histograms
    // ---------------------------------------------------------------
    const int nETbins_fine = 50;
    const float ETlo = 0, EThi = 50;
    const int nEtaBins = 220;
    const float etaLo = -1.1, etaHi = 1.1;

    // Helper lambda to create a set of histograms (nominal + optional _double)
    auto makeH1 = [&](const char *name, int nb, float lo, float hi) -> TH1D * {
        return new TH1D(name, "", nb, lo, hi);
    };
    auto makeH2 = [&](const char *name, int nbx, float lox, float hix, int nby, float loy, float hiy) -> TH2D * {
        return new TH2D(name, "", nbx, lox, hix, nby, loy, hiy);
    };

    // --- 2D Migration matrix ---
    TH2D *h2_eta_truth_vs_reco = makeH2("h2_eta_truth_vs_reco", nEtaBins, etaLo, etaHi, nEtaBins, etaLo, etaHi);
    h2_eta_truth_vs_reco->SetXTitle("#eta_{reco}");
    h2_eta_truth_vs_reco->SetYTitle("#eta_{truth}");

    TH2D *h2_eta_truth_vs_reco_tight_iso = makeH2("h2_eta_truth_vs_reco_tight_iso", nEtaBins, etaLo, etaHi, nEtaBins, etaLo, etaHi);

    // --- Fake rate (inward migration) ---
    TH1D *h_reco_ET_matched_all = makeH1("h_reco_ET_matched_all", nETbins_fine, ETlo, EThi);
    TH1D *h_reco_ET_truth_outside = makeH1("h_reco_ET_truth_outside", nETbins_fine, ETlo, EThi);
    TH1D *h_reco_ET_truth_outside_tight_iso = makeH1("h_reco_ET_truth_outside_tight_iso", nETbins_fine, ETlo, EThi);
    TH1D *h_reco_ET_truth_outside_tight_noniso = makeH1("h_reco_ET_truth_outside_tight_noniso", nETbins_fine, ETlo, EThi);
    TH1D *h_reco_ET_truth_outside_nontight_iso = makeH1("h_reco_ET_truth_outside_nontight_iso", nETbins_fine, ETlo, EThi);
    TH1D *h_reco_ET_truth_outside_nontight_noniso = makeH1("h_reco_ET_truth_outside_nontight_noniso", nETbins_fine, ETlo, EThi);

    // --- Loss rate (outward migration) ---
    TH1D *h_truth_pT_inside_all = makeH1("h_truth_pT_inside_all", nETbins_fine, ETlo, EThi);
    TH1D *h_truth_pT_reco_outside = makeH1("h_truth_pT_reco_outside", nETbins_fine, ETlo, EThi);

    // --- Response contamination ---
    TH2D *h2_response_clean = makeH2("h2_response_clean", nETbins_fine, ETlo, EThi, nETbins_fine, ETlo, EThi);
    h2_response_clean->SetXTitle("reco E_{T} [GeV]");
    h2_response_clean->SetYTitle("truth p_{T} [GeV]");

    TH2D *h2_response_contaminated = makeH2("h2_response_contaminated", nETbins_fine, ETlo, EThi, nETbins_fine, ETlo, EThi);
    h2_response_contaminated->SetXTitle("reco E_{T} [GeV]");
    h2_response_contaminated->SetYTitle("truth p_{T} [GeV]");

    // --- ABCD impact for good matches ---
    TH1D *h_good_tight_iso = makeH1("h_good_tight_iso", nETbins_fine, ETlo, EThi);
    TH1D *h_good_tight_noniso = makeH1("h_good_tight_noniso", nETbins_fine, ETlo, EThi);
    TH1D *h_good_nontight_iso = makeH1("h_good_nontight_iso", nETbins_fine, ETlo, EThi);
    TH1D *h_good_nontight_noniso = makeH1("h_good_nontight_noniso", nETbins_fine, ETlo, EThi);

    // --- Diagnostics ---
    TH1D *h_truth_eta_of_fakes = makeH1("h_truth_eta_of_fakes", nEtaBins, etaLo, etaHi);
    TH2D *h2_deta_vs_ET = makeH2("h2_deta_vs_ET", nETbins_fine, ETlo, EThi, 200, -0.5, 0.5);
    h2_deta_vs_ET->SetXTitle("reco E_{T} [GeV]");
    h2_deta_vs_ET->SetYTitle("#eta_{reco} - #eta_{truth}");

    TH1D *h_vtxz_all = makeH1("h_vtxz_all", 240, -120, 120);
    TH1D *h_vtxz_fakes = makeH1("h_vtxz_fakes", 240, -120, 120);

    // --- Background (jet MC) sideband flow histograms ---
    // Per-ABCD-region sideband flow histograms.
    //   native[r]  : truth inside, reco inside  → legitimate count
    //   inflow[r]  : truth outside, reco inside → contaminates measured region
    //   outflow[r] : truth inside, reco outside → would have been in r without migration
    //   inside[r]  : reco inside (= native + inflow)
    // Two parallel sets are booked: one with the nominal NPB preselection, one
    // with it skipped — the NPB score is -1 for outside-fiducial clusters, so
    // only the noNPB set yields non-zero outflow.
    struct BkgFlowSet {
        TH1D *inside[4]  = {nullptr, nullptr, nullptr, nullptr};
        TH1D *inflow[4]  = {nullptr, nullptr, nullptr, nullptr};
        TH1D *outflow[4] = {nullptr, nullptr, nullptr, nullptr};
        TH1D *native[4]  = {nullptr, nullptr, nullptr, nullptr};
    };
    auto bookBkgSet = [&](const char *suffix) -> BkgFlowSet {
        const char letters[4] = {'A', 'B', 'C', 'D'};
        BkgFlowSet s;
        for (int r = 0; r < 4; r++)
        {
            const char L = letters[r];
            s.inside[r]  = makeH1(Form("h_bkg_%c_inside%s",  L, suffix), nETbins_fine, ETlo, EThi);
            s.inflow[r]  = makeH1(Form("h_bkg_%c_inflow%s",  L, suffix), nETbins_fine, ETlo, EThi);
            s.outflow[r] = makeH1(Form("h_bkg_%c_outflow%s", L, suffix), nETbins_fine, ETlo, EThi);
            s.native[r]  = makeH1(Form("h_bkg_%c_native%s",  L, suffix), nETbins_fine, ETlo, EThi);
        }
        return s;
    };

    BkgFlowSet bkg_npb   = bookBkgSet("");
    BkgFlowSet bkg_noNPB = bookBkgSet("_noNPB");

    // Matching diagnostic histograms for background
    TH1D *h_bkg_all_reco_inside     = makeH1("h_bkg_all_reco_inside",     nETbins_fine, ETlo, EThi);
    TH1D *h_bkg_matched_reco_inside = makeH1("h_bkg_matched_reco_inside", nETbins_fine, ETlo, EThi);
    TH1D *h_bkg_truth_eta_inflow    = makeH1("h_bkg_truth_eta_inflow",    nEtaBins, etaLo, etaHi);
    TH2D *h2_bkg_eta_truth_vs_reco  = makeH2("h2_bkg_eta_truth_vs_reco",  nEtaBins, etaLo, etaHi, nEtaBins, etaLo, etaHi);

    // --- Double interaction histograms ---
    TH2D *h2_eta_truth_vs_reco_double = nullptr;
    TH2D *h2_eta_truth_vs_reco_tight_iso_double = nullptr;
    TH1D *h_reco_ET_matched_all_double = nullptr;
    TH1D *h_reco_ET_truth_outside_double = nullptr;
    TH1D *h_reco_ET_truth_outside_tight_iso_double = nullptr;
    TH1D *h_reco_ET_truth_outside_tight_noniso_double = nullptr;
    TH1D *h_reco_ET_truth_outside_nontight_iso_double = nullptr;
    TH1D *h_reco_ET_truth_outside_nontight_noniso_double = nullptr;
    TH1D *h_truth_pT_inside_all_double = nullptr;
    TH1D *h_truth_pT_reco_outside_double = nullptr;
    TH2D *h2_response_clean_double = nullptr;
    TH2D *h2_response_contaminated_double = nullptr;
    TH1D *h_good_tight_iso_double = nullptr;
    TH1D *h_good_tight_noniso_double = nullptr;
    TH1D *h_good_nontight_iso_double = nullptr;
    TH1D *h_good_nontight_noniso_double = nullptr;
    TH1D *h_truth_eta_of_fakes_double = nullptr;
    TH2D *h2_deta_vs_ET_double = nullptr;
    TH1D *h_vtxz_all_double = nullptr;
    TH1D *h_vtxz_fakes_double = nullptr;

    if (do_double_interaction)
    {
        h2_eta_truth_vs_reco_double = makeH2("h2_eta_truth_vs_reco_double", nEtaBins, etaLo, etaHi, nEtaBins, etaLo, etaHi);
        h2_eta_truth_vs_reco_tight_iso_double = makeH2("h2_eta_truth_vs_reco_tight_iso_double", nEtaBins, etaLo, etaHi, nEtaBins, etaLo, etaHi);
        h_reco_ET_matched_all_double = makeH1("h_reco_ET_matched_all_double", nETbins_fine, ETlo, EThi);
        h_reco_ET_truth_outside_double = makeH1("h_reco_ET_truth_outside_double", nETbins_fine, ETlo, EThi);
        h_reco_ET_truth_outside_tight_iso_double = makeH1("h_reco_ET_truth_outside_tight_iso_double", nETbins_fine, ETlo, EThi);
        h_reco_ET_truth_outside_tight_noniso_double = makeH1("h_reco_ET_truth_outside_tight_noniso_double", nETbins_fine, ETlo, EThi);
        h_reco_ET_truth_outside_nontight_iso_double = makeH1("h_reco_ET_truth_outside_nontight_iso_double", nETbins_fine, ETlo, EThi);
        h_reco_ET_truth_outside_nontight_noniso_double = makeH1("h_reco_ET_truth_outside_nontight_noniso_double", nETbins_fine, ETlo, EThi);
        h_truth_pT_inside_all_double = makeH1("h_truth_pT_inside_all_double", nETbins_fine, ETlo, EThi);
        h_truth_pT_reco_outside_double = makeH1("h_truth_pT_reco_outside_double", nETbins_fine, ETlo, EThi);
        h2_response_clean_double = makeH2("h2_response_clean_double", nETbins_fine, ETlo, EThi, nETbins_fine, ETlo, EThi);
        h2_response_contaminated_double = makeH2("h2_response_contaminated_double", nETbins_fine, ETlo, EThi, nETbins_fine, ETlo, EThi);
        h_good_tight_iso_double = makeH1("h_good_tight_iso_double", nETbins_fine, ETlo, EThi);
        h_good_tight_noniso_double = makeH1("h_good_tight_noniso_double", nETbins_fine, ETlo, EThi);
        h_good_nontight_iso_double = makeH1("h_good_nontight_iso_double", nETbins_fine, ETlo, EThi);
        h_good_nontight_noniso_double = makeH1("h_good_nontight_noniso_double", nETbins_fine, ETlo, EThi);
        h_truth_eta_of_fakes_double = makeH1("h_truth_eta_of_fakes_double", nEtaBins, etaLo, etaHi);
        h2_deta_vs_ET_double = makeH2("h2_deta_vs_ET_double", nETbins_fine, ETlo, EThi, 200, -0.5, 0.5);
        h_vtxz_all_double = makeH1("h_vtxz_all_double", 240, -120, 120);
        h_vtxz_fakes_double = makeH1("h_vtxz_fakes_double", 240, -120, 120);
    }

    // ---------------------------------------------------------------
    // Double interaction kinematics lambda
    // ---------------------------------------------------------------
    constexpr float cemc_radius_cm = 93.5f;
    auto recalculateClusterKinematics = [&](float old_et, float old_eta, float old_vtxz, float new_vtxz) -> std::pair<float, float> {
        const float old_z_over_r = std::sinh(old_eta);
        const float new_z_over_r = old_z_over_r + (old_vtxz - new_vtxz) / cemc_radius_cm;
        const float new_eta = std::asinh(new_z_over_r);
        const float cluster_energy = old_et * std::cosh(old_eta);
        const float new_et = cluster_energy / std::cosh(new_eta);
        return {new_et, new_eta};
    };

    // ---------------------------------------------------------------
    // ABCD classification lambda
    // ---------------------------------------------------------------
    // Takes all per-cluster shower shape variables and ET, returns (common_pass, tight, nontight, iso, noniso)
    struct ABCDResult {
        bool common_pass;
        bool tight;
        bool nontight;
        bool iso;
        bool noniso;
    };

    auto classifyABCD = [&](float clusterET, float recoisoET,
                            float cl_weta_cogx, float cl_wphi_cogx,
                            float cl_et1, float cl_et2, float cl_et3, float cl_et4,
                            float e11_over_e33, float e32_over_e35,
                            float cl_prob, float bdt_score, float cl_npb_score,
                            float wr_cogx, bool skip_npb = false) -> ABCDResult
    {
        ABCDResult r = {false, false, false, false, false};

        // Isolation
        float recoiso_max = recoiso_max_b + recoiso_max_s * clusterET;
        float recononiso_min = recoiso_max + recononiso_min_shift;
        if (recoisoET > recoiso_min && recoisoET < recoiso_max)
            r.iso = true;
        if (recoisoET > recononiso_min && recoisoET < recononiso_max)
            r.noniso = true;

        // Common cuts. The NPB common cut acts as a preselection that scores
        // clusters using a TMVA model trained only on in-fiducial clusters —
        // outside |eta|<0.7 the stored NPB score is the sentinel -1, so the cut
        // rejects every migrated cluster. `skip_npb` disables that cut so the
        // shower-shape/isolation sidebands remain populated for the outflow
        // measurement in the eta migration study.
        bool passes_common_shape =
            cl_prob > common_prob_min &&
            cl_prob < common_prob_max &&
            e11_over_e33 > common_e11_over_e33_min &&
            e11_over_e33 < common_e11_over_e33_max &&
            wr_cogx > common_wr_cogx_bound &&
            (cl_weta_cogx < common_cluster_weta_cogx_bound) &&
            (skip_npb || !common_npb_cut_on || cl_npb_score > common_npb_score_cut);

        if (passes_common_shape)
            r.common_pass = true;
        else
            return r;

        // ET-dependent tight cuts
        float t_weta_cogx_max = tight_weta_cogx_max_b + tight_weta_cogx_max_s * clusterET;
        float t_wphi_cogx_max = tight_wphi_cogx_max_b + tight_wphi_cogx_max_s * clusterET;
        float t_et1_min = tight_et1_min_b + tight_et1_min_s * clusterET;

        bool is_weta_tight = (cl_weta_cogx > tight_weta_cogx_min) && (cl_weta_cogx < t_weta_cogx_max);
        bool is_wphi_tight = (cl_wphi_cogx > tight_wphi_cogx_min) && (cl_wphi_cogx < t_wphi_cogx_max);
        bool is_et1_tight = (cl_et1 > t_et1_min) && (cl_et1 < tight_et1_max);
        bool is_et2_tight = (cl_et2 > tight_et2_min) && (cl_et2 < tight_et2_max);
        bool is_et3_tight = (cl_et3 > tight_et3_min) && (cl_et3 < tight_et3_max);
        bool is_e11_tight = (e11_over_e33 > tight_e11_over_e33_min) && (e11_over_e33 < tight_e11_over_e33_max);
        bool is_e32_tight = (e32_over_e35 > tight_e32_over_e35_min) && (e32_over_e35 < tight_e32_over_e35_max);
        bool is_et4_tight = (cl_et4 > tight_et4_min) && (cl_et4 < tight_et4_max);
        bool is_prob_tight = (cl_prob > tight_prob_min) && (cl_prob < tight_prob_max);

        float tight_bdt_min_et = tight_bdt_min_slope * clusterET + tight_bdt_min_intercept;
        bool is_bdt_tight = (bdt_score > tight_bdt_min_et) && (bdt_score < tight_bdt_max);

        if (is_weta_tight && is_wphi_tight && is_et1_tight && is_et2_tight && is_et3_tight &&
            is_e11_tight && is_e32_tight && is_et4_tight && is_prob_tight && is_bdt_tight)
        {
            r.tight = true;
        }

        // Non-tight classification
        if (cl_weta_cogx > non_tight_weta_cogx_min && cl_weta_cogx < non_tight_weta_cogx_max &&
            cl_wphi_cogx > non_tight_wphi_cogx_min && cl_wphi_cogx < non_tight_wphi_cogx_max &&
            cl_prob > non_tight_prob_min && cl_prob < non_tight_prob_max &&
            e11_over_e33 > non_tight_e11_over_e33_min && e11_over_e33 < non_tight_e11_over_e33_max &&
            e32_over_e35 > non_tight_e32_over_e35_min && e32_over_e35 < non_tight_e32_over_e35_max &&
            cl_et1 > non_tight_et1_min && cl_et1 < non_tight_et1_max &&
            cl_et4 > non_tight_et4_min && cl_et4 < non_tight_et4_max &&
            bdt_score > non_tight_bdt_min &&
            bdt_score < non_tight_bdt_max_slope * clusterET + non_tight_bdt_max_intercept)
        {
            int nfail = 0;
            if (!is_weta_tight) nfail += weta_on;
            if (!is_wphi_tight) nfail += wphi_on;
            if (!is_et1_tight)  nfail += et1_on;
            if (!is_et2_tight)  nfail += et2_on;
            if (!is_et3_tight)  nfail += et3_on;
            if (!is_e11_tight)  nfail += e11_to_e33_on;
            if (!is_e32_tight)  nfail += e32_to_e35_on;
            if (!is_et4_tight)  nfail += et4_on;
            if (!is_prob_tight) nfail++;
            if (!is_bdt_tight)  nfail += bdt_on;

            if (nfail > n_nt_fail)
            {
                bool all_flags_fail = true;
                if (weta_fail && is_weta_tight) all_flags_fail = false;
                if (wphi_fail && is_wphi_tight) all_flags_fail = false;
                if (et1_fail && is_et1_tight)   all_flags_fail = false;
                if (e11_to_e33_fail && is_e11_tight) all_flags_fail = false;
                if (e32_to_e35_fail && is_e32_tight) all_flags_fail = false;
                if (bdt_fail && is_bdt_tight)   all_flags_fail = false;

                if (all_flags_fail)
                    r.nontight = true;
            }
        }

        return r;
    };

    // ---------------------------------------------------------------
    // Event loop
    // ---------------------------------------------------------------
    TRandom3 rng(42);
    int ientry = 0;

    // Counters
    long long n_total_matched = 0;
    long long n_good = 0;
    long long n_fake_inward = 0;
    long long n_lost_outward = 0;
    long long n_both_outside = 0;
    long long n_fake_inward_double = 0;

    // Per-event truth-particle bookkeeping. Lifted outside the event loop so
    // the containers reuse capacity across events instead of allocating fresh
    // every time (jet samples can have hundreds of particles per event).
    // Vectors are keyed by iparticle (dense 0..nparticles-1).
    std::unordered_map<int, int> particle_trkidmap;
    std::vector<float> truth_eta_v;
    std::vector<float> truth_pt_v;
    std::vector<float> truth_iso_v;
    std::vector<unsigned char> truth_inside_fiducial;
    std::vector<unsigned char> is_registered_truth;

    while (reader.Next())
    {
        // Vertex reweighting
        float weight = cross_weight;
        float vertex_weight = 1.0;
        if (vertex_reweight_on && h_vertex_reweight)
        {
            int bin = h_vertex_reweight->FindBin(*vertexz);
            if (bin < 1) bin = 1;
            if (bin > h_vertex_reweight->GetNbinsX()) bin = h_vertex_reweight->GetNbinsX();
            vertex_weight = h_vertex_reweight->GetBinContent(bin);
        }
        if (!std::isfinite(vertex_weight) || vertex_weight <= 0.0)
            vertex_weight = 1.0;
        weight *= vertex_weight;

        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;

        // Vertex cut
        if (std::abs(*vertexz) > vertexcut)
        {
            ientry++;
            continue;
        }

        // Truth-particle bookkeeping.
        //   Signal mode: register direct/fragmentation isolated photons.
        //   Background mode: register every particle — cluster_truthtrkID picks
        //     the dominant contributor at match time.
        const int npart = *nparticles;
        particle_trkidmap.clear();
        particle_trkidmap.reserve(npart);
        truth_eta_v.assign(npart, 0.0f);
        truth_pt_v.assign(npart, 0.0f);
        truth_iso_v.assign(npart, 0.0f);
        truth_inside_fiducial.assign(npart, 0);
        is_registered_truth.assign(npart, 0);

        float maxphotonpT = 0;
        for (int iparticle = 0; iparticle < npart; iparticle++)
        {
            particle_trkidmap[particle_trkid[iparticle]] = iparticle;

            float truthisoET = 0;
            if (conesize == 4)
                truthisoET = particle_truth_iso_04[iparticle];
            else if (conesize == 3)
                truthisoET = particle_truth_iso_03[iparticle];
            else if (conesize == 2)
                truthisoET = particle_truth_iso_02[iparticle];

            if (isbackground)
            {
                float peta = particle_Eta[iparticle];
                truth_eta_v[iparticle] = peta;
                truth_pt_v[iparticle] = particle_Pt[iparticle];
                truth_iso_v[iparticle] = truthisoET;
                truth_inside_fiducial[iparticle] = (peta > eta_lo && peta < eta_hi) ? 1 : 0;
                is_registered_truth[iparticle] = 1;
            }
            else if (particle_pid[iparticle] == 22)
            {
                if (particle_Pt[iparticle] > maxphotonpT)
                    maxphotonpT = particle_Pt[iparticle];

                if (particle_photonclass[iparticle] < 3 && truthisoET < truthisocut)
                {
                    float peta = particle_Eta[iparticle];
                    truth_eta_v[iparticle] = peta;
                    truth_pt_v[iparticle] = particle_Pt[iparticle];
                    truth_iso_v[iparticle] = truthisoET;
                    truth_inside_fiducial[iparticle] = (peta > eta_lo && peta < eta_hi) ? 1 : 0;
                    is_registered_truth[iparticle] = 1;
                }
            }
        }

        // Sample kinematic window cut
        if (!isbackground)
        {
            if ((maxphotonpT > max_photon_upper) || (maxphotonpT < max_photon_lower))
            {
                ientry++;
                continue;
            }
        }
        else
        {
            // Jet pT window cut — prevents double-counting across overlapping jet samples
            float maxjetpT = 0;
            for (int ijet = 0; ijet < **njet_truth_ptr; ijet++)
            {
                if ((*jet_truth_Pt_ptr)[ijet] > maxjetpT)
                    maxjetpT = (*jet_truth_Pt_ptr)[ijet];
            }
            if (maxjetpT == 0 || maxjetpT > max_jet_upper || maxjetpT < max_jet_lower)
            {
                ientry++;
                continue;
            }
        }

        // Double interaction vertex
        float double_vtx = 0;
        if (do_double_interaction && h_vertexz_data_nocut)
        {
            float random_vertex = h_vertexz_data_nocut->GetRandom();
            double_vtx = (*vertexz + random_vertex) / 2.0;
        }

        h_vtxz_all->Fill(*vertexz, weight);
        if (do_double_interaction)
            h_vtxz_all_double->Fill(double_vtx, weight);

        // ---------------------------------------------------------------
        // Cluster loop
        // ---------------------------------------------------------------
        for (int icluster = 0; icluster < *ncluster; icluster++)
        {
            float clusterET = cluster_Et[icluster] * clusterescale;
            if (clustereres > 0)
                clusterET = clusterET * rng.Gaus(1, clustereres);

            if (clusterET < reco_min_ET)
                continue;
            if (isbackground && clusterET > cluster_ET_upper)
                continue;
            if (nosat && cluster_nsaturated[icluster] > 0)
                continue;

            float cluster_eta = cluster_Eta[icluster];

            // Shower shape variables
            float e11_over_e33 = cluster_e11[icluster] / cluster_e33[icluster];
            float e32_over_e35 = cluster_e32[icluster] / cluster_e35[icluster];
            float wr_cogx = cluster_wphi_cogx[icluster] / cluster_weta_cogx[icluster];

            // Reco isolation
            float recoisoET = -999;
            if (use_topo_iso == 2)
                recoisoET = cluster_iso_topo_04[icluster];
            else if (use_topo_iso == 1)
                recoisoET = cluster_iso_topo_03[icluster];
            else if (use_topo_iso == 4)
                recoisoET = cluster_iso_04[icluster];
            else if (use_topo_iso == 5)
                recoisoET = cluster_iso_excl_04[icluster];
            else if (conesize == 4)
                recoisoET = cluster_iso_04[icluster];
            else if (conesize == 3)
                recoisoET = cluster_iso_03[icluster];
            else if (conesize == 2)
                recoisoET = cluster_iso_02[icluster];

            // MC iso fudge
            recoisoET = recoisoET * mc_iso_scale;
            recoisoET += mc_iso_shift;

            // Select BDT model
            std::string selected_bdt_model = bdt_model_name;
            if (use_et_binned_bdt)
            {
                for (int ib = 0; ib < (int)bdt_et_bin_models.size(); ++ib)
                {
                    if (clusterET >= bdt_et_bin_edges[ib] && clusterET < bdt_et_bin_edges[ib + 1])
                    {
                        selected_bdt_model = bdt_et_bin_models[ib];
                        break;
                    }
                }
            }
            float bdt_score = (*bdt_arrays[selected_bdt_model])[icluster];

            // ABCD classification with nominal kinematics
            ABCDResult abcd = classifyABCD(clusterET, recoisoET,
                cluster_weta_cogx[icluster], cluster_wphi_cogx[icluster],
                cluster_et1[icluster], cluster_et2[icluster],
                cluster_et3[icluster], cluster_et4[icluster],
                e11_over_e33, e32_over_e35,
                cluster_prob[icluster], bdt_score, cluster_npb_score[icluster],
                wr_cogx);

            bool reco_inside = (cluster_eta > eta_lo && cluster_eta < eta_hi);

            // Background diagnostic: count every reco cluster inside fiducial, regardless of match
            if (isbackground && reco_inside)
                h_bkg_all_reco_inside->Fill(clusterET, weight);

            // Truth matching via cluster_truthtrkID
            auto trk_it = particle_trkidmap.find(cluster_truthtrkID[icluster]);
            if (trk_it == particle_trkidmap.end()) continue;
            int iparticle = trk_it->second;
            if (!is_registered_truth[iparticle]) continue;

            float truth_eta = truth_eta_v[iparticle];
            float truth_pt = truth_pt_v[iparticle];
            bool truth_inside = (truth_inside_fiducial[iparticle] != 0);

            n_total_matched++;

            h2_eta_truth_vs_reco->Fill(cluster_eta, truth_eta, weight);
            if (abcd.tight && abcd.iso)
                h2_eta_truth_vs_reco_tight_iso->Fill(cluster_eta, truth_eta, weight);
            h2_deta_vs_ET->Fill(clusterET, cluster_eta - truth_eta, weight);

            // Background-only diagnostics: matched-cluster fill, truth-vs-reco
            // eta matrix, and ABCD sideband inflow/outflow. The _noNPB pass is
            // needed because the NPB score is a sentinel for outside-fiducial
            // clusters (it would zero outflow under the standard cut).
            if (isbackground)
            {
                if (reco_inside)
                    h_bkg_matched_reco_inside->Fill(clusterET, weight);
                h2_bkg_eta_truth_vs_reco->Fill(cluster_eta, truth_eta, weight);

                ABCDResult abcd_noNPB = classifyABCD(clusterET, recoisoET,
                    cluster_weta_cogx[icluster], cluster_wphi_cogx[icluster],
                    cluster_et1[icluster], cluster_et2[icluster],
                    cluster_et3[icluster], cluster_et4[icluster],
                    e11_over_e33, e32_over_e35,
                    cluster_prob[icluster], bdt_score, cluster_npb_score[icluster],
                    wr_cogx, /*skip_npb=*/true);

                // Each ABCD result picks at most one region (A=0, B=1, C=2, D=3).
                auto fillSet = [&](const ABCDResult &res, BkgFlowSet &s)
                {
                    int region = -1;
                    if      (res.tight    && res.iso)    region = 0;
                    else if (res.tight    && res.noniso) region = 1;
                    else if (res.nontight && res.iso)    region = 2;
                    else if (res.nontight && res.noniso) region = 3;
                    if (region < 0) return;
                    if (reco_inside)                  s.inside[region]->Fill(clusterET, weight);
                    if (reco_inside &&  truth_inside) s.native[region]->Fill(clusterET, weight);
                    if (reco_inside && !truth_inside) s.inflow[region]->Fill(clusterET, weight);
                    if (!reco_inside && truth_inside) s.outflow[region]->Fill(clusterET, weight);
                };
                fillSet(abcd, bkg_npb);
                fillSet(abcd_noNPB, bkg_noNPB);

                if (reco_inside && !truth_inside)
                    h_bkg_truth_eta_inflow->Fill(truth_eta, weight);
            }

            // Classify migration category (signal-oriented accounting)
            if (truth_inside && reco_inside)
            {
                // Good: both inside fiducial
                n_good++;
                h2_response_clean->Fill(clusterET, truth_pt, weight);
                if (abcd.tight && abcd.iso)     h_good_tight_iso->Fill(clusterET, weight);
                if (abcd.tight && abcd.noniso)   h_good_tight_noniso->Fill(clusterET, weight);
                if (abcd.nontight && abcd.iso)   h_good_nontight_iso->Fill(clusterET, weight);
                if (abcd.nontight && abcd.noniso) h_good_nontight_noniso->Fill(clusterET, weight);

                // Loss rate denominator: truth inside + has matched cluster
                h_truth_pT_inside_all->Fill(truth_pt, weight);
            }
            else if (!truth_inside && reco_inside)
            {
                // Fake: inward migration (truth outside, reco inside)
                n_fake_inward++;
                h_reco_ET_truth_outside->Fill(clusterET, weight);
                h_truth_eta_of_fakes->Fill(truth_eta, weight);
                h_vtxz_fakes->Fill(*vertexz, weight);
                h2_response_contaminated->Fill(clusterET, truth_pt, weight);

                if (abcd.tight && abcd.iso)       h_reco_ET_truth_outside_tight_iso->Fill(clusterET, weight);
                if (abcd.tight && abcd.noniso)     h_reco_ET_truth_outside_tight_noniso->Fill(clusterET, weight);
                if (abcd.nontight && abcd.iso)     h_reco_ET_truth_outside_nontight_iso->Fill(clusterET, weight);
                if (abcd.nontight && abcd.noniso)  h_reco_ET_truth_outside_nontight_noniso->Fill(clusterET, weight);
            }
            else if (truth_inside && !reco_inside)
            {
                // Lost: outward migration (truth inside, reco outside)
                n_lost_outward++;
                h_truth_pT_inside_all->Fill(truth_pt, weight);
                h_truth_pT_reco_outside->Fill(truth_pt, weight);
            }
            else
            {
                // Both outside
                n_both_outside++;
            }

            // Denominator for inward fake rate: all matched with reco inside
            if (reco_inside)
                h_reco_ET_matched_all->Fill(clusterET, weight);

            // ---------------------------------------------------------------
            // Double interaction: re-classify with shifted vertex
            // ---------------------------------------------------------------
            if (do_double_interaction && h_vertexz_data_nocut)
            {
                auto [new_et, new_eta] = recalculateClusterKinematics(
                    clusterET, cluster_eta, *vertexz, double_vtx);

                if (new_et < reco_min_ET) continue;

                bool reco_inside_double = (new_eta > eta_lo && new_eta < eta_hi);

                // Re-evaluate iso threshold at shifted ET
                // Note: BDT score itself unchanged, only threshold re-evaluated
                float recoisoET_double = recoisoET; // iso energy doesn't change with vertex shift in this simplified study

                ABCDResult abcd_d = classifyABCD(new_et, recoisoET_double,
                    cluster_weta_cogx[icluster], cluster_wphi_cogx[icluster],
                    cluster_et1[icluster], cluster_et2[icluster],
                    cluster_et3[icluster], cluster_et4[icluster],
                    e11_over_e33, e32_over_e35,
                    cluster_prob[icluster], bdt_score, cluster_npb_score[icluster],
                    wr_cogx);

                h2_eta_truth_vs_reco_double->Fill(new_eta, truth_eta, weight);
                if (abcd_d.tight && abcd_d.iso)
                    h2_eta_truth_vs_reco_tight_iso_double->Fill(new_eta, truth_eta, weight);

                h2_deta_vs_ET_double->Fill(new_et, new_eta - truth_eta, weight);

                if (truth_inside && reco_inside_double)
                {
                    h2_response_clean_double->Fill(new_et, truth_pt, weight);
                    if (abcd_d.tight && abcd_d.iso)     h_good_tight_iso_double->Fill(new_et, weight);
                    if (abcd_d.tight && abcd_d.noniso)   h_good_tight_noniso_double->Fill(new_et, weight);
                    if (abcd_d.nontight && abcd_d.iso)   h_good_nontight_iso_double->Fill(new_et, weight);
                    if (abcd_d.nontight && abcd_d.noniso) h_good_nontight_noniso_double->Fill(new_et, weight);
                    h_truth_pT_inside_all_double->Fill(truth_pt, weight);
                }
                else if (!truth_inside && reco_inside_double)
                {
                    n_fake_inward_double++;
                    h_reco_ET_truth_outside_double->Fill(new_et, weight);
                    h_truth_eta_of_fakes_double->Fill(truth_eta, weight);
                    h_vtxz_fakes_double->Fill(double_vtx, weight);
                    h2_response_contaminated_double->Fill(new_et, truth_pt, weight);

                    if (abcd_d.tight && abcd_d.iso)       h_reco_ET_truth_outside_tight_iso_double->Fill(new_et, weight);
                    if (abcd_d.tight && abcd_d.noniso)     h_reco_ET_truth_outside_tight_noniso_double->Fill(new_et, weight);
                    if (abcd_d.nontight && abcd_d.iso)     h_reco_ET_truth_outside_nontight_iso_double->Fill(new_et, weight);
                    if (abcd_d.nontight && abcd_d.noniso)  h_reco_ET_truth_outside_nontight_noniso_double->Fill(new_et, weight);
                }
                else if (truth_inside && !reco_inside_double)
                {
                    h_truth_pT_inside_all_double->Fill(truth_pt, weight);
                    h_truth_pT_reco_outside_double->Fill(truth_pt, weight);
                }

                if (reco_inside_double)
                    h_reco_ET_matched_all_double->Fill(new_et, weight);

                h_vtxz_all_double->Fill(double_vtx, weight);
            }

        } // end cluster loop

        ientry++;
    } // end event loop

    // ---------------------------------------------------------------
    // Summary statistics
    // ---------------------------------------------------------------
    std::cout << "\n========== Eta Migration Study Summary ==========" << std::endl;
    std::cout << "Sample: " << filetype << std::endl;
    std::cout << "Total events processed: " << ientry << std::endl;
    std::cout << "Total matched truth-reco pairs: " << n_total_matched << std::endl;
    std::cout << "  Good (both inside |eta|<0.7):  " << n_good << std::endl;
    std::cout << "  Fakes (inward migration):      " << n_fake_inward << std::endl;
    std::cout << "  Losses (outward migration):    " << n_lost_outward << std::endl;
    std::cout << "  Both outside:                  " << n_both_outside << std::endl;
    if (n_good + n_fake_inward > 0)
    {
        float fake_rate = (float)n_fake_inward / (float)(n_good + n_fake_inward);
        std::cout << "Inward migration fake rate (unweighted): "
                  << fake_rate * 100.0 << "%" << std::endl;
    }
    if (n_good + n_lost_outward > 0)
    {
        float loss_rate = (float)n_lost_outward / (float)(n_good + n_lost_outward);
        std::cout << "Outward migration loss rate (unweighted): "
                  << loss_rate * 100.0 << "%" << std::endl;
    }
    if (do_double_interaction)
    {
        std::cout << "Double interaction inward fakes: " << n_fake_inward_double << std::endl;
    }
    std::cout << "=================================================" << std::endl;

    // ---------------------------------------------------------------
    // Write output
    // ---------------------------------------------------------------
    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    SaveYamlToRoot(fout, configname.c_str());

    // Write all nominal histograms
    h2_eta_truth_vs_reco->Write();
    h2_eta_truth_vs_reco_tight_iso->Write();
    h_reco_ET_matched_all->Write();
    h_reco_ET_truth_outside->Write();
    h_reco_ET_truth_outside_tight_iso->Write();
    h_reco_ET_truth_outside_tight_noniso->Write();
    h_reco_ET_truth_outside_nontight_iso->Write();
    h_reco_ET_truth_outside_nontight_noniso->Write();
    h_truth_pT_inside_all->Write();
    h_truth_pT_reco_outside->Write();
    h2_response_clean->Write();
    h2_response_contaminated->Write();
    h_good_tight_iso->Write();
    h_good_tight_noniso->Write();
    h_good_nontight_iso->Write();
    h_good_nontight_noniso->Write();
    h_truth_eta_of_fakes->Write();
    h2_deta_vs_ET->Write();
    h_vtxz_all->Write();
    h_vtxz_fakes->Write();

    auto writeBkgSet = [](const BkgFlowSet &s) {
        for (int r = 0; r < 4; r++) s.inside[r]->Write();
        for (int r = 0; r < 4; r++) s.inflow[r]->Write();
        for (int r = 0; r < 4; r++) s.outflow[r]->Write();
        for (int r = 0; r < 4; r++) s.native[r]->Write();
    };
    writeBkgSet(bkg_npb);
    writeBkgSet(bkg_noNPB);
    h_bkg_all_reco_inside->Write();
    h_bkg_matched_reco_inside->Write();
    h_bkg_truth_eta_inflow->Write();
    h2_bkg_eta_truth_vs_reco->Write();

    // Write double interaction histograms
    if (do_double_interaction)
    {
        h2_eta_truth_vs_reco_double->Write();
        h2_eta_truth_vs_reco_tight_iso_double->Write();
        h_reco_ET_matched_all_double->Write();
        h_reco_ET_truth_outside_double->Write();
        h_reco_ET_truth_outside_tight_iso_double->Write();
        h_reco_ET_truth_outside_tight_noniso_double->Write();
        h_reco_ET_truth_outside_nontight_iso_double->Write();
        h_reco_ET_truth_outside_nontight_noniso_double->Write();
        h_truth_pT_inside_all_double->Write();
        h_truth_pT_reco_outside_double->Write();
        h2_response_clean_double->Write();
        h2_response_contaminated_double->Write();
        h_good_tight_iso_double->Write();
        h_good_tight_noniso_double->Write();
        h_good_nontight_iso_double->Write();
        h_good_nontight_noniso_double->Write();
        h_truth_eta_of_fakes_double->Write();
        h2_deta_vs_ET_double->Write();
        h_vtxz_all_double->Write();
        h_vtxz_fakes_double->Write();
    }

    fout->Close();
    std::cout << "Output written to: " << outfilename << std::endl;

    // Cleanup BDT arrays
    for (auto &kv : bdt_arrays)
        delete kv.second;
}
