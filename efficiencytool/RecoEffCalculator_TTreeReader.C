#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <vector>
#include <sstream>
#include <cmath>
#include <limits>
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
#include <TEfficiency.h>
#include <TObjString.h>
#include <TRandom3.h>
#include <yaml-cpp/yaml.h>
#include "CrossSectionWeights.h"
#include "TruthVertexReweightLoader.h"
using namespace PPG12;
// unfolding
#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>

// R__LOAD_LIBRARY(/sphenix/user/egm2153/calib_study/JetValidation/analysis/roounfold/libRooUnfold.so)
const float TIME_SAMPLE_NS = 17.6;
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

void RecoEffCalculator_TTreeReader(const std::string &configname = "config_bdt_nom.yaml", const std::string filetype = "jet40", bool do_vertex_scan = false, float mix_weight = 1.0, const std::string vtxscan_sim_override = "")
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

    // _nom aliases read from the base single-interaction (SI) sample input path
    // (same physics, output label only). Pairs the *_nom label with the existing
    // nominal sample trees so every SI sample has a matching *_double DI partner
    // for the cross-section DI blending pipeline (oneforall_tree_double.sh).
    std::string input_filetype = filetype;
    if (filetype == "photon5_nom")  input_filetype = "photon5";
    if (filetype == "photon10_nom") input_filetype = "photon10";
    if (filetype == "photon20_nom") input_filetype = "photon20";
    if (filetype == "jet8_nom")     input_filetype = "jet8";
    if (filetype == "jet12_nom")    input_filetype = "jet12";
    if (filetype == "jet20_nom")    input_filetype = "jet20";
    if (filetype == "jet30_nom")    input_filetype = "jet30";
    if (filetype == "jet40_nom")    input_filetype = "jet40";

    std::string infilename = infilename_root_dir + input_filetype + infilename_branch_dir;
    //infilename = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/auau_test/caloana_with_bdt.root";
    if (!issim)
    {
        infilename = configYaml["input"]["data_file"].as<std::string>();
    }

    std::cout << "infilename: " << infilename << std::endl;

    // Sample kinematic windows and cross-section weight from shared lookup
    PPG12::SampleConfig sc = PPG12::GetSampleConfig(filetype);
    float max_photon_lower = sc.photon_pt_lower;
    float max_photon_upper = sc.photon_pt_upper;
    float max_jet_lower    = sc.jet_pt_lower;
    float max_jet_upper    = sc.jet_pt_upper;
    float cluster_ET_upper = sc.cluster_ET_upper;
    float weight           = sc.weight;

    float energy_scale_lower = 0;
    float energy_scale_upper = 100;

    float vertex_weight = 1.0;
    float cross_weight = weight;
    cross_weight *= mix_weight;   // single/double interaction blending fraction (1.0 for data)

    /*
 1  Constant     1.07368e+00   1.10770e-03  -6.32308e-06  -1.15696e-02
2  Mean        -1.73712e+00   8.31735e-02  -1.32021e-04  -1.25750e-03
3  Sigma        4.47289e+01   2.37499e-01   2.37499e-01  -5.45865e-02
*/
    // Vertex reweighting for simulation:
    //   On-the-fly only — built per-config from the first-pass vtxscan files
    //   (data and sim h_vertexz). The static fallback file is deprecated as of
    //   the 8ac4d1b commit (it was outdated and produced biased h_truth_pT).
    TH1* h_vertex_reweight = nullptr;
    int vertex_reweight_on = 1;
    std::string vtx_scan_data_file = "";
    // Data-driven iterative truth-vertex reweight (mutually exclusive with
    // the f(z_r) reco-level reweight above).
    TH1* h_truth_vtx_reweight = nullptr;
    int truth_vertex_reweight_on = 0;
    std::string truth_vertex_reweight_file;
    if (issim)
    {
        vertex_reweight_on = configYaml["analysis"]["vertex_reweight_on"].as<int>(1);
        vtx_scan_data_file = configYaml["analysis"]["vertex_scan_data_file"].as<std::string>("");
        truth_vertex_reweight_on = configYaml["analysis"]["truth_vertex_reweight_on"].as<int>(0);
        truth_vertex_reweight_file =
            configYaml["analysis"]["truth_vertex_reweight_file"].as<std::string>("");
        if (truth_vertex_reweight_on && vertex_reweight_on)
        {
            std::cout << "[TruthVertexReweight] forcing vertex_reweight_on=0 "
                      << "(mutually exclusive)." << std::endl;
            vertex_reweight_on = 0;
        }
    }

    // std::string infilename = "/sphenix/tg/tg01/commissioning/CaloCalibWG/sli/ppg12/ana450/condorout/combine.root";
    // Build input chain and connect reader
    std::string treename = configYaml["input"]["tree"].as<std::string>();
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());
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

    // Derive the vertex-scan output filename from outfilename
    std::string vtxscan_outfilename = outfilename;
    {
        size_t pos = vtxscan_outfilename.rfind(".root");
        if (pos != std::string::npos) vtxscan_outfilename.replace(pos, 5, "_vtxscan.root");
    }

    // In Pass 2, allow the caller to supply a pre-merged combined vtxscan file
    // so vertex reweighting uses the blended MC distribution.
    if (!vtxscan_sim_override.empty() && !do_vertex_scan && issim)
    {
        std::cout << "[VertexReweight] Using vtxscan_sim_override: "
                  << vtxscan_sim_override << std::endl;
        vtxscan_outfilename = vtxscan_sim_override;
    }

    // Vertex reweighting in second pass:
    // Build the data/sim ratio on-the-fly from the first-pass scan files
    // (h_vertexz in both data and sim vtxscan ROOT files). This is now the
    // ONLY supported mode — the static fallback file was deprecated because
    // it became outdated relative to the current data vertex distribution.
    if (!do_vertex_scan && issim && vertex_reweight_on)
    {
        std::string resolved_vtx_scan_data_file = vtx_scan_data_file;
        if (resolved_vtx_scan_data_file.empty())
        {
            resolved_vtx_scan_data_file =
                configYaml["output"]["data_outfile"].as<std::string>() + "_" + var_type + "_vtxscan.root";
            std::cout << "[VertexReweight] vertex_scan_data_file is empty. "
                      << "Trying first-pass scan files:\n"
                      << "  data: " << resolved_vtx_scan_data_file << "\n"
                      << "  sim : " << vtxscan_outfilename << std::endl;
        }

        bool built_from_scan = false;
        TFile* fvtx_data = TFile::Open(resolved_vtx_scan_data_file.c_str(), "READ");
        TFile* fvtx_sim = TFile::Open(vtxscan_outfilename.c_str(), "READ");
        if (fvtx_data && !fvtx_data->IsZombie() && fvtx_sim && !fvtx_sim->IsZombie())
        {
            TH1* hdata_vtx = dynamic_cast<TH1*>(fvtx_data->Get("h_vertexz"));
            TH1* hsim_vtx  = dynamic_cast<TH1*>(fvtx_sim->Get("h_vertexz"));
            if (hdata_vtx && hsim_vtx)
            {
                TH1* hdata_clone = dynamic_cast<TH1*>(hdata_vtx->Clone("h_vtx_data_clone"));
                TH1* hsim_clone  = dynamic_cast<TH1*>(hsim_vtx->Clone("h_vtx_sim_clone"));
                hdata_clone->SetDirectory(nullptr);
                hsim_clone->SetDirectory(nullptr);

                // Rebin to wider bins to improve statistics in the tails.
                // Source histogram is 200 bins of 1 cm from -100 to +100 cm.
                // Rebin factor 5 → 40 bins of 5 cm width. This makes the data/sim ratio
                // well-defined further out into the tails where 1 cm bins would be empty.
                const int rebin_factor = configYaml["analysis"]["vertex_reweight_rebin"].as<int>(5);
                if (rebin_factor > 1)
                {
                    hdata_clone->Rebin(rebin_factor);
                    hsim_clone->Rebin(rebin_factor);
                }

                if (hdata_clone->Integral() > 0) hdata_clone->Scale(1.0 / hdata_clone->Integral());
                if (hsim_clone->Integral() > 0)  hsim_clone->Scale(1.0 / hsim_clone->Integral());
                hdata_clone->Divide(hsim_clone);

                // Smooth the ratio to interpolate over residual low-statistics fluctuations.
                // TH1::Smooth(n) applies the 353QH twice-smoothing algorithm n times.
                const int smooth_passes = configYaml["analysis"]["vertex_reweight_smooth"].as<int>(1);
                if (smooth_passes > 0)
                {
                    hdata_clone->Smooth(smooth_passes);
                }

                h_vertex_reweight = hdata_clone;
                h_vertex_reweight->SetDirectory(nullptr);
                built_from_scan = true;
                std::cout << "[VertexReweight] On-the-fly weights derived from "
                          << resolved_vtx_scan_data_file << " / " << vtxscan_outfilename << std::endl;
                std::cout << "[VertexReweight] Rebinned by factor " << rebin_factor
                          << " (now " << h_vertex_reweight->GetNbinsX() << " bins, "
                          << h_vertex_reweight->GetBinWidth(1) << " cm wide); smoothed "
                          << smooth_passes << " pass(es)." << std::endl;
                delete hsim_clone;
            }
        }

        if (fvtx_data) { fvtx_data->Close(); delete fvtx_data; }
        if (fvtx_sim)  { fvtx_sim->Close();  delete fvtx_sim;  }

        if (!built_from_scan)
        {
            // The static fallback file (vertex_reweight_bdt_none.root) was found to
            // be outdated — built from an older, wider data vertex distribution that
            // does not match the current 1.5 mrad / 0 mrad data. Using it gives a
            // bloated h_truth_pT integral (h_vertexcut/h_truth ~ 0.36 instead of
            // ~0.99 expected for the current 1.5 mrad data). On-the-fly is now the
            // ONLY supported path (for the reco-vertex reweight branch). The
            // nominal pipeline now uses truth-vertex reweight (truth_vertex_reweight_on=1),
            // which bypasses this branch entirely. This message fires only if someone
            // explicitly sets vertex_reweight_on=1 without generating vtxscan files.
            std::cerr << "[VertexReweight] FATAL: first-pass vtxscan files are missing or invalid:\n"
                      << "  data: " << resolved_vtx_scan_data_file << "\n"
                      << "  sim : " << vtxscan_outfilename << "\n"
                      << "Either (a) set truth_vertex_reweight_on=1 to use the single-pass "
                      << "truth-vertex reweight (recommended), or (b) run oneforall_tree.sh / "
                      << "oneforall_tree_double.sh first to produce the vtxscan files." << std::endl;
            return;
        }
    }

    if (issim && truth_vertex_reweight_on)
    {
        h_truth_vtx_reweight = LoadTruthVertexReweight(truth_vertex_reweight_file);
        if (!h_truth_vtx_reweight) return;
    }

    // TChain is used instead of a single TTree

    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();

    std::string bdt_model_name = configYaml["input"]["bdt_model_name"].as<std::string>("base");
    //std::string bdt_model_name = "base";

    std::vector<float> bdt_et_bin_edges;
    std::vector<std::string> bdt_et_bin_models;
    bool use_et_binned_bdt = false;

    if (configYaml["input"]["bdt_et_bin_edges"] && configYaml["input"]["bdt_et_bin_models"]) {
        for (auto v : configYaml["input"]["bdt_et_bin_edges"])
            bdt_et_bin_edges.push_back(v.as<float>());
        for (auto v : configYaml["input"]["bdt_et_bin_models"])
            bdt_et_bin_models.push_back(v.as<std::string>());
        use_et_binned_bdt = (bdt_et_bin_models.size() == bdt_et_bin_edges.size() - 1);
        if (!use_et_binned_bdt)
            std::cout << "WARNING: bdt_et_bin_edges/bdt_et_bin_models size mismatch; falling back to single model" << std::endl;
    }
    int iso_threshold = configYaml["analysis"]["iso_threshold"].as<int>(0);
    int iso_hcalonly = configYaml["analysis"]["iso_hcalonly"].as<int>(0);
    int use_topo_iso = configYaml["analysis"]["use_topo_iso"].as<int>(0);
    // we have 0, 0.05, 0.1, 0.2 options
    float iso_emcalinnerr = configYaml["analysis"]["iso_emcalinnerr"].as<float>(0.0);
    std::cout<<"iso_emcalinnerr: "<<iso_emcalinnerr<<std::endl;
    float iso_topo_innerr = configYaml["analysis"]["iso_topo_innerr"].as<float>(0.0);
    std::cout<<"iso_topo_innerr: "<<iso_topo_innerr<<std::endl;

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
    int nosat = configYaml["analysis"]["nosat"].as<int>(0);
    int common_b2bjet_cut = configYaml["analysis"]["common_b2bjet_cut"].as<int>(0);
    float common_b2bjet_pt_min = configYaml["analysis"]["common_b2bjet_pt_min"].as<float>(7.0);

    float cluster_mbd_time_min = configYaml["analysis"]["cluster_mbd_time_min"].as<float>(-999.0f);
    float cluster_mbd_time_max = configYaml["analysis"]["cluster_mbd_time_max"].as<float>(999.0f);
    std::cout << "cluster_mbd_time_min: " << cluster_mbd_time_min
              << "  cluster_mbd_time_max: " << cluster_mbd_time_max << std::endl;

    float truthisocut = configYaml["analysis"]["truth_iso_max"].as<float>();

    float recoiso_min = configYaml["analysis"]["reco_iso_min"].as<float>();
    float recoiso_max_b = configYaml["analysis"]["reco_iso_max_b"].as<float>();
    float recoiso_max_s = configYaml["analysis"]["reco_iso_max_s"].as<float>();

    float recononiso_min_shift = configYaml["analysis"]["reco_noniso_min_shift"].as<float>();
    float recononiso_max = configYaml["analysis"]["reco_noniso_max"].as<float>();

    float vertexcut = configYaml["analysis"]["vertex_cut"].as<float>();
    // Truth-vertex denominator cut for the MBD-eff plumbing (line ~1785).
    // Defaults to the reco vertex_cut so legacy configs are unchanged.
    // Setting vertex_cut_truth >> vertexcut (e.g. 9999) makes the MBD-eff
    // denominator span all truth events while the analysis-fiducial
    // |z_reco|<vertexcut stays in place — the standard MBD-eff
    // parameterization paired with beam-delivered (allz) lumi.
    float vertexcut_truth = configYaml["analysis"]["vertex_cut_truth"].as<float>(vertexcut);

    // Per-event lumi weight: pre-scales MC fills to a merge-target lumi so
    // plain hadd across periods reproduces the all-range expectation. When
    // lumi_target == lumi (default), lumi_weight = 1 → standalone per-period
    // analysis unchanged. Merge-feeder configs set lumi_target = sum(L_periods)
    // and the resulting per-period MC files are pre-scaled contributions to
    // the all-range MC. Data is normalized downstream in CalculatePhotonYield
    // via 1/binwidth/lumi, so the lumi_weight is gated to MC only (issim).
    float lumi = configYaml["analysis"]["lumi"].as<float>(1.0);
    float lumi_target = configYaml["analysis"]["lumi_target"].as<float>(lumi);
    if (issim) {
        cross_weight *= lumi / lumi_target;
    }
    std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();
    int n_eta_bins = eta_bins.size() - 1;
    double eta_bin_edges[n_eta_bins + 1];
    std::copy(eta_bins.begin(), eta_bins.end(), eta_bin_edges);

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

    float mc_iso_shift = configYaml["analysis"]["mc_iso_shift"].as<float>(0.0);
    float mc_iso_scale = configYaml["analysis"]["mc_iso_scale"].as<float>(1.0);
    float mbd_avg_sigma_min = configYaml["analysis"]["mbd_avg_sigma_min"].as<float>(
        -std::numeric_limits<float>::infinity());
    float mbd_avg_sigma_max = configYaml["analysis"]["mbd_avg_sigma_max"].as<float>(
        std::numeric_limits<float>::infinity());

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

    float tight_bdt_max = configYaml["analysis"]["tight"]["bdt_max"].as<float>(1);
    float tight_bdt_min = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0);
    float tight_bdt_min_slope = configYaml["analysis"]["tight"]["bdt_min_slope"].as<float>(0);
    float tight_bdt_min_intercept = configYaml["analysis"]["tight"]["bdt_min_intercept"].as<float>(tight_bdt_min);

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

    float non_tight_bdt_max = configYaml["analysis"]["non_tight"]["bdt_max"].as<float>(1);
    float non_tight_bdt_min = configYaml["analysis"]["non_tight"]["bdt_min"].as<float>(0);
    float non_tight_bdt_max_slope = configYaml["analysis"]["non_tight"]["bdt_max_slope"].as<float>(0);
    float non_tight_bdt_max_intercept = configYaml["analysis"]["non_tight"]["bdt_max_intercept"].as<float>(non_tight_bdt_max);
    float non_tight_bdt_min_slope = configYaml["analysis"]["non_tight"]["bdt_min_slope"].as<float>(0);
    float non_tight_bdt_min_intercept = configYaml["analysis"]["non_tight"]["bdt_min_intercept"].as<float>(non_tight_bdt_min);

    // common cuts for both tight and non tight

    float common_prob_max = configYaml["analysis"]["common"]["prob_max"].as<float>();
    float common_prob_min = configYaml["analysis"]["common"]["prob_min"].as<float>();

    float common_e11_over_e33_max = configYaml["analysis"]["common"]["e11_over_e33_max"].as<float>();
    float common_e11_over_e33_min = configYaml["analysis"]["common"]["e11_over_e33_min"].as<float>();

    float common_wr_cogx_bound = configYaml["analysis"]["common"]["wr_cogx_bound"].as<float>();
    float common_cluster_weta_cogx_bound = configYaml["analysis"]["common"]["cluster_weta_cogx_bound"].as<float>();

    int common_npb_cut_on = configYaml["analysis"]["common"]["npb_cut_on"].as<int>(0);
    float common_npb_score_cut = configYaml["analysis"]["common"]["npb_score_cut"].as<float>(0.5);

    int reweight = configYaml["analysis"]["unfold"]["reweight"].as<int>(); // 0 for no reweighting, 1 for reweighting

    float clusterescale = configYaml["analysis"]["cluster_escale"].as<float>(1.0);
    float clustereres = configYaml["analysis"]["cluster_eres"].as<float>(0.0);

    // ------------------------------------------------------------
    // Tower-mask config: veto clusters whose center tower is flagged
    // in an external binary TH2I mask. Applied symmetrically to data
    // and MC so efficiency corrections remain self-consistent. Used
    // for the acceptance systematic variants that emulate the real
    // detector's dead-tower list.
    // ------------------------------------------------------------
    int tower_mask_on = configYaml["analysis"]["tower_mask_on"].as<int>(0);
    std::string tower_mask_file = configYaml["analysis"]["tower_mask_file"].as<std::string>("");
    std::string tower_mask_name = configYaml["analysis"]["tower_mask_name"].as<std::string>("");
    TH2I *h_tower_mask = nullptr;
    TFile *f_tower_mask = nullptr;
    if (tower_mask_on && !tower_mask_file.empty() && !tower_mask_name.empty()) {
        f_tower_mask = TFile::Open(tower_mask_file.c_str(), "READ");
        if (!f_tower_mask || f_tower_mask->IsZombie()) {
            std::cerr << "[tower-mask] FATAL: cannot open " << tower_mask_file << std::endl;
            return;
        }
        h_tower_mask = dynamic_cast<TH2I *>(f_tower_mask->Get(tower_mask_name.c_str()));
        if (!h_tower_mask) {
            std::cerr << "[tower-mask] FATAL: histogram '" << tower_mask_name
                      << "' not found in " << tower_mask_file << std::endl;
            return;
        }
        h_tower_mask->SetDirectory(nullptr);  // detach from file so close doesn't delete
        int n_masked = (int) h_tower_mask->Integral();
        std::cout << "[tower-mask] loaded " << tower_mask_name
                  << " from " << tower_mask_file
                  << ": " << n_masked << " masked towers" << std::endl;
    }

    // polynomial 3 for the reweighting
    // TF1 *f_reweight = new TF1("f_reweight", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0, 100);
    TF1 *f_reweight = new TF1("f_reweight", "([0] + [1]*x + [3]*x*x) / (1 + [2]*x + [4]*x*x)", 0, 100);
    // need to make this into the config file in the future!!!
    // f_reweight->SetParameters(1.04713, 0.00623875, -0.00106856, 2.64199e-06);
    // f_reweight->SetParameters(1.04713, 0.00623875, -0.00106856, 2.64199e-06);
    // f_reweight->SetParameters(0.714962, -0.0856443, -0.125383, 0.00345831, 0.00462972);
    f_reweight->SetParameters(0.787183, -0.109761, -0.137895, 0.00449023, 0.00533778);

    //load MBD t0 correction
    std::cout << "loading MBD t0 correction" << std::endl;
    std::map<int, float> mbd_t0_correction;
    ifstream file = ifstream("/sphenix/user/shuhangli/ppg12/efficiencytool/MbdOut.corr");
    string line;
    while (getline(file, line))
    {
        std::istringstream iss(line);
        int runnumber;
        float t0;
        iss >> runnumber >> t0;
        mbd_t0_correction[runnumber] = t0;
        std::cout << "runnumber: " << runnumber << " t0: " << t0 << std::endl;
    }
    std::cout << "loaded MBD t0 correction" << std::endl;

    // TTreeReader setup
    TTreeReader reader(&chain);
    
    // Basic event variables
    TTreeReaderValue<int> mbdnorthhit(reader, "mbdnorthhit");
    TTreeReaderValue<int> mbdsouthhit(reader, "mbdsouthhit");
    // Per-channel MBD info (64 channels per side)
    TTreeReaderArray<float> mbdnorthq(reader, "mbdnorthq");
    TTreeReaderArray<float> mbdsouthq(reader, "mbdsouthq");
    TTreeReaderArray<float> mbdnortht(reader, "mbdnortht");
    TTreeReaderArray<float> mbdsoutht(reader, "mbdsoutht");
    TTreeReaderValue<int> pythiaid(reader, "pythiaid");
    TTreeReaderValue<int> nparticles(reader, "nparticles");
    TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
    TTreeReaderValue<int> runnumber(reader, "runnumber");
    TTreeReaderArray<Bool_t> scaledtrigger(reader, "scaledtrigger");
    TTreeReaderArray<Bool_t> livetrigger(reader, "livetrigger");
    TTreeReaderValue<float> energy_scale(reader, "energy_scale");
    TTreeReaderValue<float> vertexz(reader, "vertexz");
    TTreeReaderValue<float> vertexz_truth(reader, "vertexz_truth");
    // Optional; only present in reprocessed double-interaction slimtrees.
    std::unique_ptr<TTreeReaderValue<float>> vertexz_truth_mb_ptr;
    if (chain.GetBranch("vertexz_truth_mb"))
        vertexz_truth_mb_ptr.reset(new TTreeReaderValue<float>(reader, "vertexz_truth_mb"));
    TTreeReaderValue<float> mbdnorthtmean(reader, "mbdnorthtmean");
    TTreeReaderValue<float> mbdsouthtmean(reader, "mbdsouthtmean");
    TTreeReaderValue<float> mbd_time(reader, "mbd_time");
    TTreeReaderArray<float> trigger_prescale(reader, "trigger_prescale");

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
    TTreeReaderArray<float> cluster_iso_005(reader, Form("cluster_iso_005_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_0075(reader, Form("cluster_iso_0075_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_excl_04(reader, Form("cluster_iso_excl_04_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_excl_005(reader, Form("cluster_iso_excl_005_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_excl_0075(reader, Form("cluster_iso_excl_0075_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_excl_01(reader, Form("cluster_iso_excl_01_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_excl_02(reader, Form("cluster_iso_excl_02_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_topo_03(reader, Form("cluster_iso_topo_03_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_topo_04(reader, Form("cluster_iso_topo_04_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_topo_005 (reader, Form("cluster_iso_topo_005_%s",  clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_topo_0075(reader, Form("cluster_iso_topo_0075_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_topo_01  (reader, Form("cluster_iso_topo_01_%s",   clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_topo_02  (reader, Form("cluster_iso_topo_02_%s",   clusternodename.c_str()));
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
    TTreeReaderArray<int> cluster_nsaturated(reader, Form("cluster_nsaturated_%s", clusternodename.c_str()));

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

    // Cluster arrays for 2D data - using regular TTreeReaderArray for 2D arrays
    TTreeReaderArray<float> cluster_e_array(reader, Form("cluster_e_array_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_adc_array(reader, Form("cluster_adc_array_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_time_array(reader, Form("cluster_time_array_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_status_array(reader, Form("cluster_status_array_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ownership_array(reader, Form("cluster_ownership_array_%s", clusternodename.c_str()));

    // Cluster isolation arrays
    TTreeReaderArray<float> cluster_iso_03_emcal(reader, Form("cluster_iso_03_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_hcalin(reader, Form("cluster_iso_03_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_hcalout(reader, Form("cluster_iso_03_hcalout_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_emcal(reader, Form("cluster_iso_03_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_02_70_emcal(reader, Form("cluster_iso_02_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_005_70_emcal(reader, Form("cluster_iso_005_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_0075_70_emcal(reader, Form("cluster_iso_0075_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_01_70_emcal(reader, Form("cluster_iso_01_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalin(reader, Form("cluster_iso_03_70_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalout(reader, Form("cluster_iso_03_70_hcalout_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_120_emcal(reader, Form("cluster_iso_03_120_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_120_hcalin(reader, Form("cluster_iso_03_120_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_120_hcalout(reader, Form("cluster_iso_03_120_hcalout_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_emcal(reader, Form("cluster_iso_03_60_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_hcalin(reader, Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_hcalout(reader, Form("cluster_iso_03_60_hcalout_%s", clusternodename.c_str()));

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

    // BDT score — build one TTreeReaderArray per unique model name
    std::vector<std::string> all_bdt_models = {bdt_model_name};
    if (use_et_binned_bdt)
        all_bdt_models.insert(all_bdt_models.end(), bdt_et_bin_models.begin(), bdt_et_bin_models.end());
    std::sort(all_bdt_models.begin(), all_bdt_models.end());
    all_bdt_models.erase(std::unique(all_bdt_models.begin(), all_bdt_models.end()), all_bdt_models.end());

    std::map<std::string, TTreeReaderArray<float>*> bdt_arrays;
    for (auto& mname : all_bdt_models) {
        bdt_arrays[mname] = new TTreeReaderArray<float>(reader,
            Form("cluster_bdt_%s_%s", clusternodename.c_str(), mname.c_str()));
    }

    // NPB score
    TTreeReaderArray<float> cluster_npb_score(reader, Form("cluster_npb_score_%s", clusternodename.c_str()));

    // Double-interaction MC uses r04 jet branch names
    // use_r04_branches: fresh slimtrees (run28 production) use AntiKt_Truth_r04 naming for
    // ALL samples. Default true matches current production state; legacy slimtrees with
    // the older plain "njet_truth" naming can override to false via config.
    bool use_r04_branches = configYaml["analysis"]["use_r04_branches"].as<bool>(true);
    std::string njet_truth_bname    = use_r04_branches ? "njet_truth_AntiKt_Truth_r04"          : "njet_truth";
    std::string jet_truth_E_bname   = use_r04_branches ? "jet_truth_E_AntiKt_Truth_r04"         : "jet_truth_E";
    std::string jet_truth_Pt_bname  = use_r04_branches ? "jet_truth_Pt_AntiKt_Truth_r04"        : "jet_truth_Pt";
    std::string jet_truth_Eta_bname = use_r04_branches ? "jet_truth_Eta_AntiKt_Truth_r04"       : "jet_truth_Eta";
    std::string jet_truth_Phi_bname = use_r04_branches ? "jet_truth_Phi_AntiKt_Truth_r04"       : "jet_truth_Phi";

    // Truth jet arrays
    TTreeReaderValue<int> njet_truth(reader, njet_truth_bname.c_str());
    TTreeReaderArray<float> jet_truth_E(reader, jet_truth_E_bname.c_str());
    TTreeReaderArray<float> jet_truth_Pt(reader, jet_truth_Pt_bname.c_str());
    TTreeReaderArray<float> jet_truth_Eta(reader, jet_truth_Eta_bname.c_str());
    TTreeReaderArray<float> jet_truth_Phi(reader, jet_truth_Phi_bname.c_str());

    // MC slimtrees store JetCalib-calibrated jets (_AntiKt_unsubtracted_r04_calib);
    // data slimtrees (ana521 production) use plain _AntiKt_unsubtracted_r04 without _calib.
    // Keep use_r04_branches=false as the legacy escape for old trees.
    const std::string jet_suffix = use_r04_branches
        ? (issim ? "_AntiKt_unsubtracted_r04_calib" : "_AntiKt_unsubtracted_r04")
        : std::string("");
    std::string njet_bname    = use_r04_branches ? ("njet"    + jet_suffix) : std::string("njet");
    std::string jet_E_bname   = use_r04_branches ? ("jet_E"   + jet_suffix) : std::string("jet_E");
    std::string jet_Pt_bname  = use_r04_branches ? ("jet_Pt"  + jet_suffix) : std::string("jet_Pt");
    std::string jet_Eta_bname = use_r04_branches ? ("jet_Eta" + jet_suffix) : std::string("jet_Eta");
    std::string jet_Phi_bname = use_r04_branches ? ("jet_Phi" + jet_suffix) : std::string("jet_Phi");

    // Reco jet arrays
    TTreeReaderValue<int> njet(reader, njet_bname.c_str());
    TTreeReaderArray<float> jet_E(reader, jet_E_bname.c_str());
    TTreeReaderArray<float> jet_Pt(reader, jet_Pt_bname.c_str());
    TTreeReaderArray<float> jet_Eta(reader, jet_Eta_bname.c_str());
    TTreeReaderArray<float> jet_Phi(reader, jet_Phi_bname.c_str());


    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    TH1::SetDefaultSumw2(kTRUE);
    TH1F *h_max_photon_pT = new TH1F("h_max_photon_pT", "Max Photon pT", 1000, 0, 100);
    TH1F *h_photon_pT = new TH1F("h_photon_pT", "Photon pT", 1000, 0, 100);
    TH1F *h_max_direct_pT = new TH1F("h_max_direct_pT", "Max Direct Photon pT", 1000, 0, 100);
    TH1F *h_direct_pT = new TH1F("h_direct_pT", "Direct Photon pT", 1000, 0, 100);
    TH1F *h_max_frag_pT = new TH1F("h_max_frag_pT", "Max Fragmentation Photon pT", 1000, 0, 100);
    TH1F *h_frag_pT = new TH1F("h_frag_pT", "Fragmentation Photon pT", 1000, 0, 100);
    TH1F *h_max_decay_pT = new TH1F("h_max_decay_pT", "Max Decay Photon pT", 1000, 0, 100);
    TH1F *h_decay_photon_pT = new TH1F("h_decay_photon_pT", "Decay Photon pT", 1000, 0, 100);
    TH1F *h_vertexz = new TH1F("h_vertexz", "Vertex z", 200, -100, 100);
    TH1F *h_cluster_common_Et = new TH1F("h_cluster_common_E", "Cluster Common E", 1000, 0, 100);
    TH1F *h_cluster_common_leading_Et = new TH1F("h_cluster_common_leading_E", "Cluster Common Leading E", 1000, 0, 100);

    TH1F *h_max_truth_jet_pT = new TH1F("h_max_truth_jet_pT", "Max Truth Jet pT", 1000, 0, 100);

    TH1F *h_reco_jet_pT     = new TH1F("h_reco_jet_pT",     "Reco Jet pT",      1000, 0, 100);
    TH1F *h_max_reco_jet_pT = new TH1F("h_max_reco_jet_pT", "Max Reco Jet pT",  1000, 0, 100);

    TH1F *h_max_photon_pT_vertexcut = new TH1F("h_max_photon_pT_vertexcut", "Max Photon pT Vertex Cut", 1000, 0, 100);

    TH1F *h_max_photon_pT_vertexcut_mbd_cut = new TH1F("h_max_photon_pT_vertexcut_mbd_cut", "Max Photon pT Vertex Cut MBD Cut", 1000, 0, 100);

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
    std::vector<TEfficiency *> eff_all_eta;

    // truth pythia
    std::vector<TH1D *> h_truth_pT;
    std::vector<TH1D *> h_truth_pT_novtx; // truth without vertex_weight (pure generator prediction)

    std::vector<TH1D *> h_truth_pT_vertexcut;

    std::vector<TH1D *> h_truth_pT_vertexcut_mbd_cut;
    std::vector<TH1D *> h_truth_pT_vertexcut_mbd_north_cut;
    std::vector<TH1D *> h_truth_pT_vertexcut_mbd_south_cut;
    std::vector<TH1D *> h_truth_pT_vertexcut_mbd_only_north;
    std::vector<TH1D *> h_truth_pT_vertexcut_mbd_only_south;
    std::vector<TH1D *> h_truth_pT_vertexcut_mbd_neither;

    std::vector<TH1D *> h_tight_iso_cluster_signal;
    std::vector<TH1D *> h_tight_noniso_cluster_signal;
    std::vector<TH1D *> h_nontight_iso_cluster_signal;
    std::vector<TH1D *> h_nontight_noniso_cluster_signal;
    std::vector<TH1D *> h_tight_iso_cluster_notmatch;
    std::vector<TH1D *> h_tight_noniso_cluster_notmatch;
    std::vector<TH1D *> h_nontight_iso_cluster_notmatch;
    std::vector<TH1D *> h_nontight_noniso_cluster_notmatch;
    std::vector<TH1D *> h_all_cluster_signal;
    std::vector<TH2D *> h_all_cluster_Et_max_b2bjet; // max cluster Et vs max backtobjets Et
    std::vector<TH1D *> h_tight_cluster_signal;

    std::vector<TH1D *> h_tight_iso_cluster_background;
    std::vector<TH1D *> h_tight_noniso_cluster_background;
    std::vector<TH1D *> h_nontight_iso_cluster_background;
    std::vector<TH1D *> h_nontight_noniso_cluster_background;
    std::vector<TH2D *> h_tight_recoisoET_background;
    std::vector<TH2D *> h_nontight_recoisoET_background;

    std::vector<TH2D *> h_singal_reco_isoET;
    std::vector<TH2D *> h_singal_truth_isoET;
    std::vector<TH2D *> h_background_truth_isoET;

    // here are for the plots we gonna make for both data and simulation
    std::vector<TH1D *> h_tight_iso_cluster;
    std::vector<TH1D *> h_tight_noniso_cluster;
    std::vector<TH1D *> h_nontight_iso_cluster;
    std::vector<TH1D *> h_nontight_noniso_cluster;
    std::vector<TH1D *> h_common_cluster;
    std::vector<TH1D *> h_all_cluster;
    std::vector<TH1D *> h_tight_cluster;
    // delta t between cluster and mbd
    std::vector<TH2D *> h_delta_t_all_cluster;
    std::vector<TH2D *> h_delta_t_tight_iso_cluster;
    std::vector<TH2D *> h_delta_t_tight_noniso_cluster;
    std::vector<TH2D *> h_delta_t_nontight_iso_cluster;
    std::vector<TH2D *> h_delta_t_nontight_noniso_cluster;
    std::vector<TH2D *> h_delta_t_common_cluster;
    std::vector<TH2D *> h_delta_t_tight_cluster;
    std::vector<TH2D *> h_delta_t_tight_cluster_b2bjet;
    std::vector<TH2D *> h_delta_t_nontight_iso_cluster_b2bjet;
    std::vector<TH2D *> h_delta_t_nontight_noniso_cluster_b2bjet;
    std::vector<TH2D *> h_delta_t_npb_cluster;
    //mbd t vs cluster t
    TH2D* h_mbd_t_vs_cluster_t = new TH2D("h_mbd_t_vs_cluster_t", "MBD T vs Cluster T", 160, -40, 40, 160, -40, 40);

    // Per-event std dev of MBD per-channel times (q>0.1) north vs south
    TH2D* h_mbd_time_std_north_vs_south = new TH2D(
        "h_mbd_time_std_north_vs_south",
        "MBD per-event time std dev (channels with q>0.1);#sigma(t_{north}) [ns];#sigma(t_{south}) [ns]",
        200, 0, 20,
        200, 0, 20);

    // Tower-index (ietacent x iphicent) acceptance maps at 4 selection levels,
    // filled for clusters with ET >= tower_map_et_min (default 8 GeV = analysis
    // cross-section lower edge). Uses the same per-event `weight` as the
    // adjacent fills (MC: cross_section x lumi/lumi_target x vertex_weight x
    // truth_vertex_reweight; data: prescale). MergeSim plain-hadds these so the
    // merged MC output is the inclusive (signal + jet background) tower map --
    // apples-to-apples with data for dead-tower diagnosis.
    const float tower_map_et_min = configYaml["analysis"]["tower_map_et_min"].as<float>(8.0f);
    TH2F *h_etaphi_tower_preselect = new TH2F(
        "h_etaphi_tower_preselect",
        Form("Cluster tower map (ET>=%.1f GeV, preselect);cluster i#eta;cluster i#phi", tower_map_et_min),
        96, 0, 96, 256, 0, 256);
    TH2F *h_etaphi_tower_common = new TH2F(
        "h_etaphi_tower_common",
        Form("Cluster tower map (ET>=%.1f GeV, common);cluster i#eta;cluster i#phi", tower_map_et_min),
        96, 0, 96, 256, 0, 256);
    TH2F *h_etaphi_tower_tight = new TH2F(
        "h_etaphi_tower_tight",
        Form("Cluster tower map (ET>=%.1f GeV, tight);cluster i#eta;cluster i#phi", tower_map_et_min),
        96, 0, 96, 256, 0, 256);
    TH2F *h_etaphi_tower_tight_iso = new TH2F(
        "h_etaphi_tower_tight_iso",
        Form("Cluster tower map (ET>=%.1f GeV, tight + iso);cluster i#eta;cluster i#phi", tower_map_et_min),
        96, 0, 96, 256, 0, 256);

    // unfold response matrix
    std::vector<RooUnfoldResponse *> responses_full;
    std::vector<RooUnfoldResponse *> responses_half;
    // vector for the response matrix th2
    std::vector<TH2D *> h_response_full_list;
    std::vector<TH2D *> h_response_half_list;
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

    // xjgamma, why not?
    // this is for sim only
    std::vector<TH2D *> h_tight_iso_xjgamma_signal;
    std::vector<TH2D *> h_tight_noniso_xjgamma_signal;
    std::vector<TH2D *> h_nontight_iso_xjgamma_signal;
    std::vector<TH2D *> h_nontight_noniso_xjgamma_signal;
    std::vector<TH2D *> h_all_xjgamma_signal;
    std::vector<TH2D *> h_tight_xjgamma_signal;

    // using truth jet
    std::vector<TH2D *> h_tight_iso_truthjet_xjgamma_signal;
    std::vector<TH2D *> h_tight_noniso_truthjet_xjgamma_signal;
    std::vector<TH2D *> h_nontight_iso_truthjet_xjgamma_signal;
    std::vector<TH2D *> h_nontight_noniso_truthjet_xjgamma_signal;
    std::vector<TH2D *> h_all_truthjet_xjgamma_signal;
    std::vector<TH2D *> h_tight_truthjet_xjgamma_signal;

    std::vector<TH2D *> h_tight_iso_xjgamma_background;
    std::vector<TH2D *> h_tight_noniso_xjgamma_background;
    std::vector<TH2D *> h_nontight_iso_xjgamma_background;
    std::vector<TH2D *> h_nontight_noniso_xjgamma_background;

    //
    std::vector<TH2D *> h_tight_iso_pid_pt;
    std::vector<TH2D *> h_tight_noniso_pid_pt;
    std::vector<TH2D *> h_nontight_iso_pid_pt;
    std::vector<TH2D *> h_nontight_noniso_pid_pt;
    std::vector<TH2D *> h_common_pid_pt;

    // for both data and sim
    std::vector<TH2D *> h_tight_iso_xjgamma;
    std::vector<TH2D *> h_tight_noniso_xjgamma;
    std::vector<TH2D *> h_nontight_iso_xjgamma;
    std::vector<TH2D *> h_nontight_noniso_xjgamma;
    std::vector<TH2D *> h_common_xjgamma;
    std::vector<TH2D *> h_all_xjgamma;
    std::vector<TH2D *> h_tight_xjgamma;

    // isolation profile for debugging reasons
    std::vector<std::vector<TH1D *>> h_tight_cluster_pT;
    h_tight_cluster_pT.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH1D *>> h_nontight_cluster_pT;
    h_nontight_cluster_pT.resize(eta_bins.size() - 1);
    // truth iso vs reco iso for different pT bins and eta bins
    std::vector<std::vector<TH2D *>> h_iso_truth_reco;
    h_iso_truth_reco.resize(eta_bins.size() - 1);
    std::vector<std::vector<TH2D *>> h_background_iso_truth_reco;
    h_background_iso_truth_reco.resize(eta_bins.size() - 1);
    // response vs. isoET
    std::vector<std::vector<TH2D *>> h_response_isoET;
    h_response_isoET.resize(eta_bins.size() - 1);

    //vertex bins for testing efficiency
    std::vector<double> vertex_bins = {0, 30, 60, 100};
    std::vector<TH2D*> h_vertex_efficiency_denominator;
    std::vector<TH2D*> h_vertex_efficiency_reco;
    std::vector<TH2D*> h_vertex_efficiency_id;
    std::vector<TH2D*> h_vertex_efficiency_iso;
    //2d histogram using et and eta bins
    for (int ivtx = 0; ivtx < (int)vertex_bins.size() - 1; ivtx++)
    {
        h_vertex_efficiency_denominator.push_back(new TH2D(Form("h_vertex_efficiency_denominator_%d", ivtx),
                                                           Form("Vertex Efficiency Denominator %.1f < vertex < %.1f", vertex_bins[ivtx], vertex_bins[ivtx + 1]),
                                                           n_pT_bins_truth, pT_bin_edges_truth, n_eta_bins, eta_bin_edges));
        h_vertex_efficiency_reco.push_back(new TH2D(Form("h_vertex_efficiency_reco_%d", ivtx),
                                                     Form("Vertex Efficiency Reco %.1f < vertex < %.1f", vertex_bins[ivtx], vertex_bins[ivtx + 1]),
                                                     n_pT_bins_truth, pT_bin_edges_truth, n_eta_bins, eta_bin_edges));
        h_vertex_efficiency_id.push_back(new TH2D(Form("h_vertex_efficiency_id_%d", ivtx),
                                                   Form("Vertex Efficiency ID %.1f < vertex < %.1f", vertex_bins[ivtx], vertex_bins[ivtx + 1]),
                                                   n_pT_bins_truth, pT_bin_edges_truth, n_eta_bins, eta_bin_edges));
        h_vertex_efficiency_iso.push_back(new TH2D(Form("h_vertex_efficiency_iso_%d", ivtx),
                                                    Form("Vertex Efficiency Iso %.1f < vertex < %.1f", vertex_bins[ivtx], vertex_bins[ivtx + 1]),
                                                    n_pT_bins_truth, pT_bin_edges_truth, n_eta_bins, eta_bin_edges));
    }

    std::vector<TH2D*> h_vertex_tight_iso_cluster_signal_pt_eta;
    std::vector<TH2D*> h_vertex_tight_noniso_cluster_signal_pt_eta;
    std::vector<TH2D*> h_vertex_nontight_iso_cluster_signal_pt_eta;
    std::vector<TH2D*> h_vertex_nontight_noniso_cluster_signal_pt_eta;
    std::vector<TH2D*> h_vertex_tight_iso_cluster_pt_eta;
    std::vector<TH2D*> h_vertex_tight_noniso_cluster_pt_eta;
    std::vector<TH2D*> h_vertex_nontight_iso_cluster_pt_eta;
    std::vector<TH2D*> h_vertex_nontight_noniso_cluster_pt_eta;

    for (int ivtx = 0; ivtx < (int)vertex_bins.size() - 1; ivtx++)
    {
        h_vertex_tight_iso_cluster_signal_pt_eta.push_back(new TH2D(Form("h_vertex_tight_iso_cluster_signal_pt_eta_%d", ivtx),
                                                                   Form("Vertex Tight Iso Cluster Signal %.1f < vertex < %.1f", vertex_bins[ivtx], vertex_bins[ivtx + 1]),
                                                                   n_pT_bins, pT_bin_edges, n_eta_bins, eta_bin_edges));
        h_vertex_tight_noniso_cluster_signal_pt_eta.push_back(new TH2D(Form("h_vertex_tight_noniso_cluster_signal_pt_eta_%d", ivtx),
                                                                       Form("Vertex Tight Non-Iso Cluster Signal %.1f < vertex < %.1f", vertex_bins[ivtx], vertex_bins[ivtx + 1]),
                                                                       n_pT_bins, pT_bin_edges, n_eta_bins, eta_bin_edges));
        h_vertex_nontight_iso_cluster_signal_pt_eta.push_back(new TH2D(Form("h_vertex_nontight_iso_cluster_signal_pt_eta_%d", ivtx),
                                                                       Form("Vertex Non-Tight Iso Cluster Signal %.1f < vertex < %.1f", vertex_bins[ivtx], vertex_bins[ivtx + 1]),
                                                                       n_pT_bins, pT_bin_edges, n_eta_bins, eta_bin_edges));
        h_vertex_nontight_noniso_cluster_signal_pt_eta.push_back(new TH2D(Form("h_vertex_nontight_noniso_cluster_signal_pt_eta_%d", ivtx),
                                                                         Form("Vertex Non-Tight Non-Iso Cluster Signal %.1f < vertex < %.1f", vertex_bins[ivtx], vertex_bins[ivtx + 1]),
                                                                         n_pT_bins, pT_bin_edges, n_eta_bins, eta_bin_edges));
        h_vertex_tight_iso_cluster_pt_eta.push_back(new TH2D(Form("h_vertex_tight_iso_cluster_pt_eta_%d", ivtx),
                                                             Form("Vertex Tight Iso Cluster %.1f < vertex < %.1f", vertex_bins[ivtx], vertex_bins[ivtx + 1]),
                                                             n_pT_bins, pT_bin_edges, n_eta_bins, eta_bin_edges));
        h_vertex_tight_noniso_cluster_pt_eta.push_back(new TH2D(Form("h_vertex_tight_noniso_cluster_pt_eta_%d", ivtx),
                                                                 Form("Vertex Tight Non-Iso Cluster %.1f < vertex < %.1f", vertex_bins[ivtx], vertex_bins[ivtx + 1]),
                                                                 n_pT_bins, pT_bin_edges, n_eta_bins, eta_bin_edges));
        h_vertex_nontight_iso_cluster_pt_eta.push_back(new TH2D(Form("h_vertex_nontight_iso_cluster_pt_eta_%d", ivtx),
                                                                 Form("Vertex Non-Tight Iso Cluster %.1f < vertex < %.1f", vertex_bins[ivtx], vertex_bins[ivtx + 1]),
                                                                 n_pT_bins, pT_bin_edges, n_eta_bins, eta_bin_edges));
        h_vertex_nontight_noniso_cluster_pt_eta.push_back(new TH2D(Form("h_vertex_nontight_noniso_cluster_pt_eta_%d", ivtx),
                                                                   Form("Vertex Non-Tight Non-Iso Cluster %.1f < vertex < %.1f", vertex_bins[ivtx], vertex_bins[ivtx + 1]),
                                                                   n_pT_bins, pT_bin_edges, n_eta_bins, eta_bin_edges));
    }

    

    for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
    {
        eff_reco_eta.push_back(new TEfficiency(Form("eff_reco_eta_%d", ieta),
                                               Form("Reco Efficiency %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                               n_pT_bins_truth, pT_bin_edges_truth));
        eff_reco_eta[ieta]->SetStatisticOption(effopt);
        eff_reco_eta[ieta]->SetUseWeightedEvents();

        eff_id_eta.push_back(new TEfficiency(Form("eff_id_eta_%d", ieta),
                                             Form("ID Efficiency %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                             n_pT_bins_truth, pT_bin_edges_truth));
        eff_id_eta[ieta]->SetStatisticOption(effopt);
        eff_id_eta[ieta]->SetUseWeightedEvents();

        eff_converts_eta.push_back(new TEfficiency(Form("eff_converts_eta_%d", ieta),
                                                   Form("Conversion Prob %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                   n_pT_bins_truth, pT_bin_edges_truth));
        eff_converts_eta[ieta]->SetStatisticOption(effopt);
        eff_converts_eta[ieta]->SetUseWeightedEvents();

        eff_iso_eta.push_back(new TEfficiency(Form("eff_iso_eta_%d", ieta),
                                              Form("Iso Efficiency %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                              n_pT_bins_truth, pT_bin_edges_truth));
        eff_iso_eta[ieta]->SetStatisticOption(effopt);
        eff_iso_eta[ieta]->SetUseWeightedEvents();

        eff_all_eta.push_back(new TEfficiency(Form("eff_all_eta_%d", ieta),
                                              Form("All Efficiency %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                              n_pT_bins_truth, pT_bin_edges_truth));
        eff_all_eta[ieta]->SetStatisticOption(effopt);
        eff_all_eta[ieta]->SetUseWeightedEvents();

        h_truth_pT.push_back(new TH1D(Form("h_truth_pT_%d", ieta),
                                      Form("Truth pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                      n_pT_bins_truth, pT_bin_edges_truth));
        h_truth_pT_novtx.push_back(new TH1D(Form("h_truth_pT_novtx_%d", ieta),
                                             Form("Truth pT (no vtx wt) %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                             n_pT_bins_truth, pT_bin_edges_truth));

        h_truth_pT_vertexcut.push_back(new TH1D(Form("h_truth_pT_vertexcut_%d", ieta),
                                                Form("Truth pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                n_pT_bins_truth, pT_bin_edges_truth));

        h_truth_pT_vertexcut_mbd_cut.push_back(new TH1D(Form("h_truth_pT_vertexcut_mbd_cut_%d", ieta),
                                                        Form("Truth pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                        n_pT_bins_truth, pT_bin_edges_truth));
        h_truth_pT_vertexcut_mbd_north_cut.push_back(new TH1D(Form("h_truth_pT_vertexcut_mbd_north_cut_%d", ieta),
                                                              Form("Truth pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                              n_pT_bins_truth, pT_bin_edges_truth));
        h_truth_pT_vertexcut_mbd_south_cut.push_back(new TH1D(Form("h_truth_pT_vertexcut_mbd_south_cut_%d", ieta),
                                                              Form("Truth pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                              n_pT_bins_truth, pT_bin_edges_truth));
        h_truth_pT_vertexcut_mbd_only_north.push_back(new TH1D(Form("h_truth_pT_vertexcut_mbd_only_north_%d", ieta),
                                                               Form("Truth pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                               n_pT_bins_truth, pT_bin_edges_truth));
        h_truth_pT_vertexcut_mbd_only_south.push_back(new TH1D(Form("h_truth_pT_vertexcut_mbd_only_south_%d", ieta),
                                                               Form("Truth pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                               n_pT_bins_truth, pT_bin_edges_truth));
        h_truth_pT_vertexcut_mbd_neither.push_back(new TH1D(Form("h_truth_pT_vertexcut_mbd_neither_%d", ieta),
                                                            Form("Truth pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                            n_pT_bins_truth, pT_bin_edges_truth));

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
        h_tight_iso_cluster_notmatch.push_back(new TH1D(Form("h_tight_iso_cluster_notmatch_%d", ieta),
                                                        Form("Tight Iso Not Matched %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                        n_pT_bins, pT_bin_edges));
        h_tight_noniso_cluster_notmatch.push_back(new TH1D(Form("h_tight_noniso_cluster_notmatch_%d", ieta),
                                                           Form("Tight Non-Iso Not Matched %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                           n_pT_bins, pT_bin_edges));
        h_nontight_iso_cluster_notmatch.push_back(new TH1D(Form("h_nontight_iso_cluster_notmatch_%d", ieta),
                                                           Form("Non-Tight Iso Not Matched %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                           n_pT_bins, pT_bin_edges));
        h_nontight_noniso_cluster_notmatch.push_back(new TH1D(Form("h_nontight_noniso_cluster_notmatch_%d", ieta),
                                                              Form("Non-Tight Non-Iso Not Matched %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                              n_pT_bins, pT_bin_edges));
        h_all_cluster_signal.push_back(new TH1D(Form("h_all_cluster_signal_%d", ieta),
                                                Form("All Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                n_pT_bins, pT_bin_edges));
        h_all_cluster_Et_max_b2bjet.push_back(new TH2D(Form("h_all_cluster_Et_max_b2bjet_%d", ieta),
                                                       Form("All Cluster Et Max Backtobjets %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                       n_pT_bins, pT_bin_edges, 100, 0, 100));

        h_tight_cluster_signal.push_back(new TH1D(Form("h_tight_cluster_signal_%d", ieta),
                                                  Form("Tight Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                  n_pT_bins, pT_bin_edges));

        h_tight_iso_cluster_background.push_back(new TH1D(Form("h_tight_iso_cluster_background_%d", ieta),
                                                          Form("Tight Iso Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                          n_pT_bins, pT_bin_edges));
        h_tight_noniso_cluster_background.push_back(new TH1D(Form("h_tight_noniso_cluster_background_%d", ieta),
                                                             Form("Tight Non-Iso Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                             n_pT_bins, pT_bin_edges));
        h_nontight_iso_cluster_background.push_back(new TH1D(Form("h_nontight_iso_cluster_background_%d", ieta),
                                                             Form("Non-Tight Iso Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                             n_pT_bins, pT_bin_edges));
        h_nontight_noniso_cluster_background.push_back(new TH1D(Form("h_nontight_noniso_cluster_background_%d", ieta),
                                                                Form("Non-Tight Non-Iso Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                                n_pT_bins, pT_bin_edges));
        h_tight_recoisoET_background.push_back(new TH2D(Form("h_tight_recoisoET_background_%d", ieta),
                                                        Form("Background Tight Cluster ET vs Reco Iso ET %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                        n_pT_bins, pT_bin_edges, 4400, -5, 50));
        h_nontight_recoisoET_background.push_back(new TH2D(Form("h_nontight_recoisoET_background_%d", ieta),
                                                           Form("Background Non-Tight Cluster ET vs Reco Iso ET %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                           n_pT_bins, pT_bin_edges, 4400, -5, 50));

        h_singal_reco_isoET.push_back(new TH2D(Form("h_singal_reco_isoET_%d", ieta),
                                               Form("Signal Reco Iso ET %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                               400, 0, 50, 4400, -5, 50));

        h_singal_truth_isoET.push_back(new TH2D(Form("h_singal_truth_isoET_%d", ieta),
                                                Form("Signal Truth Iso ET %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                400, 0, 50, 4400, -5, 50));

        h_background_truth_isoET.push_back(new TH2D(Form("h_background_truth_isoET_%d", ieta),
                                                    Form("Background Truth Iso ET %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
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
        h_common_cluster.push_back(new TH1D(Form("h_common_cluster_%d", ieta),
                                            Form("Common Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                            n_pT_bins, pT_bin_edges));
        h_all_cluster.push_back(new TH1D(Form("h_all_cluster_%d", ieta),
                                         Form("All Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                         n_pT_bins, pT_bin_edges));
        h_tight_cluster.push_back(new TH1D(Form("h_tight_cluster_%d", ieta),
                                           Form("Tight Cluster %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                           n_pT_bins, pT_bin_edges));

        h_delta_t_all_cluster.push_back(new TH2D(Form("h_delta_t_all_cluster_%d", ieta),
                                                  Form("Delta T Between All Cluster and MBD %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                n_pT_bins, pT_bin_edges, 160, -40, 40));
        h_delta_t_tight_iso_cluster.push_back(new TH2D(Form("h_delta_t_tight_iso_cluster_%d", ieta),
                                                  Form("Delta T Between Tight Iso Cluster and MBD %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                n_pT_bins, pT_bin_edges, 160, -40, 40));
        h_delta_t_tight_noniso_cluster.push_back(new TH2D(Form("h_delta_t_tight_noniso_cluster_%d", ieta),
                                                  Form("Delta T Between Tight Non-Iso Cluster and MBD %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                n_pT_bins, pT_bin_edges, 160, -40, 40));
        h_delta_t_nontight_iso_cluster.push_back(new TH2D(Form("h_delta_t_nontight_iso_cluster_%d", ieta),
                                                  Form("Delta T Between Non-Tight Iso Cluster and MBD %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                n_pT_bins, pT_bin_edges, 160, -40, 40));
        h_delta_t_nontight_noniso_cluster.push_back(new TH2D(Form("h_delta_t_nontight_noniso_cluster_%d", ieta),
                                                  Form("Delta T Between Non-Tight Non-Iso Cluster and MBD %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                n_pT_bins, pT_bin_edges, 160, -40, 40));
        h_delta_t_common_cluster.push_back(new TH2D(Form("h_delta_t_common_cluster_%d", ieta),
                                                  Form("Delta T Between Common Cluster and MBD %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                n_pT_bins, pT_bin_edges, 160, -40, 40));
        h_delta_t_tight_cluster.push_back(new TH2D(Form("h_delta_t_tight_cluster_%d", ieta),
                                                  Form("Delta T Between Tight Cluster and MBD %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                n_pT_bins, pT_bin_edges, 160, -40, 40));
        h_delta_t_tight_cluster_b2bjet.push_back(new TH2D(Form("h_delta_t_tight_cluster_b2bjet_%d", ieta),
                                                  Form("Delta T Between Tight Cluster and B2BJet %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                n_pT_bins, pT_bin_edges, 160, -40, 40));
        h_delta_t_nontight_iso_cluster_b2bjet.push_back(new TH2D(Form("h_delta_t_nontight_iso_cluster_b2bjet_%d", ieta),
                                                  Form("Delta T Between Non-Tight Iso Cluster and MBD %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                n_pT_bins, pT_bin_edges, 160, -40, 40));
        h_delta_t_nontight_noniso_cluster_b2bjet.push_back(new TH2D(Form("h_delta_t_nontight_noniso_cluster_b2bjet_%d", ieta),
                                                  Form("Delta T Between Non-Tight Non-Iso Cluster and MBD %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                n_pT_bins, pT_bin_edges, 160, -40, 40));
        h_delta_t_npb_cluster.push_back(new TH2D(Form("h_delta_t_npb_cluster_%d", ieta),
                                                  Form("Delta T Between NPB Cluster and MBD %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                n_pT_bins, pT_bin_edges, 160, -40, 40));
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

        h_response_full_list.push_back(h_response_full);

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

        h_response_half_list.push_back(h_response_half);

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

        // xjgamma
        h_tight_iso_xjgamma_signal.push_back(new TH2D(Form("h_tight_iso_xjgamma_signal_%d", ieta),
                                                      Form("Tight Iso XJGamma Signal %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                      300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_tight_noniso_xjgamma_signal.push_back(new TH2D(Form("h_tight_noniso_xjgamma_signal_%d", ieta),
                                                         Form("Tight Non-Iso XJGamma Signal %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                         300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_nontight_iso_xjgamma_signal.push_back(new TH2D(Form("h_nontight_iso_xjgamma_signal_%d", ieta),
                                                         Form("Non-Tight Iso XJGamma Signal %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                         300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_nontight_noniso_xjgamma_signal.push_back(new TH2D(Form("h_nontight_noniso_xjgamma_signal_%d", ieta),
                                                            Form("Non-Tight Non-Iso XJGamma Signal %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                            150, 0, 1.5, n_pT_bins, pT_bin_edges));
        h_all_xjgamma_signal.push_back(new TH2D(Form("h_all_xjgamma_signal_%d", ieta),
                                                Form("All XJGamma Signal %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_tight_xjgamma_signal.push_back(new TH2D(Form("h_tight_xjgamma_signal_%d", ieta),
                                                  Form("Tight XJGamma Signal %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                  300, 0, 3.0, n_pT_bins, pT_bin_edges));
        // using truth jet
        h_tight_iso_truthjet_xjgamma_signal.push_back(new TH2D(Form("h_tight_iso_truthjet_xjgamma_signal_%d", ieta),
                                                               Form("Tight Iso TruthJet XJGamma Signal %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                               300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_tight_noniso_truthjet_xjgamma_signal.push_back(new TH2D(Form("h_tight_noniso_truthjet_xjgamma_signal_%d", ieta),
                                                                  Form("Tight Non-Iso TruthJet XJGamma Signal %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                                  300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_nontight_iso_truthjet_xjgamma_signal.push_back(new TH2D(Form("h_nontight_iso_truthjet_xjgamma_signal_%d", ieta),
                                                                  Form("Non-Tight Iso TruthJet XJGamma Signal %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                                  300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_nontight_noniso_truthjet_xjgamma_signal.push_back(new TH2D(Form("h_nontight_noniso_truthjet_xjgamma_signal_%d", ieta),
                                                                     Form("Non-Tight Non-Iso TruthJet XJGamma Signal %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                                     300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_all_truthjet_xjgamma_signal.push_back(new TH2D(Form("h_all_truthjet_xjgamma_signal_%d", ieta),
                                                         Form("All TruthJet XJGamma Signal %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                         300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_tight_truthjet_xjgamma_signal.push_back(new TH2D(Form("h_tight_truthjet_xjgamma_signal_%d", ieta),
                                                           Form("Tight TruthJet XJGamma Signal %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                           300, 0, 3.0, n_pT_bins, pT_bin_edges));
        // background
        h_tight_iso_xjgamma_background.push_back(new TH2D(Form("h_tight_iso_xjgamma_background_%d", ieta),
                                                          Form("Tight Iso XJGamma Background %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                          300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_tight_noniso_xjgamma_background.push_back(new TH2D(Form("h_tight_noniso_xjgamma_background_%d", ieta),
                                                             Form("Tight Non-Iso XJGamma Background %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                             300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_nontight_iso_xjgamma_background.push_back(new TH2D(Form("h_nontight_iso_xjgamma_background_%d", ieta),
                                                             Form("Non-Tight Iso XJGamma Background %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                             300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_nontight_noniso_xjgamma_background.push_back(new TH2D(Form("h_nontight_noniso_xjgamma_background_%d", ieta),
                                                                Form("Non-Tight Non-Iso XJGamma Background %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                                300, 0, 3.0, n_pT_bins, pT_bin_edges));
        // pid vs pT
        h_tight_iso_pid_pt.push_back(new TH2D(Form("h_tight_iso_pid_pt_%d", ieta),
                                              Form("Tight Iso PID pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                              8000, -4000, 4000, n_pT_bins, pT_bin_edges));
        h_tight_noniso_pid_pt.push_back(new TH2D(Form("h_tight_noniso_pid_pt_%d", ieta),
                                                 Form("Tight Non-Iso PID pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                 8000, -4000, 4000, n_pT_bins, pT_bin_edges));
        h_nontight_iso_pid_pt.push_back(new TH2D(Form("h_nontight_iso_pid_pt_%d", ieta),
                                                    Form("Non-Tight Iso PID pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                    8000, -4000, 4000, n_pT_bins, pT_bin_edges));       
        h_nontight_noniso_pid_pt.push_back(new TH2D(Form("h_nontight_noniso_pid_pt_%d", ieta),
                                                        Form("Non-Tight Non-Iso PID pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                        8000, -4000, 4000, n_pT_bins, pT_bin_edges));                   
        h_common_pid_pt.push_back(new TH2D(Form("h_common_pid_pt_%d", ieta),
                                            Form("Common PID pT %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                            8000, -4000, 4000, n_pT_bins, pT_bin_edges));           


        h_tight_iso_xjgamma.push_back(new TH2D(Form("h_tight_iso_xjgamma_%d", ieta),
                                               Form("Tight Iso XJGamma %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                               300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_tight_noniso_xjgamma.push_back(new TH2D(Form("h_tight_noniso_xjgamma_%d", ieta),
                                                  Form("Tight Non-Iso XJGamma %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                  300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_nontight_iso_xjgamma.push_back(new TH2D(Form("h_nontight_iso_xjgamma_%d", ieta),
                                                  Form("Non-Tight Iso XJGamma %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                  300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_nontight_noniso_xjgamma.push_back(new TH2D(Form("h_nontight_noniso_xjgamma_%d", ieta),
                                                     Form("Non-Tight Non-Iso XJGamma %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                                     300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_common_xjgamma.push_back(new TH2D(Form("h_common_xjgamma_%d", ieta),
                                            Form("Common XJGamma %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                            300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_all_xjgamma.push_back(new TH2D(Form("h_all_xjgamma_%d", ieta),
                                         Form("All XJGamma %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                         300, 0, 3.0, n_pT_bins, pT_bin_edges));
        h_tight_xjgamma.push_back(new TH2D(Form("h_tight_xjgamma_%d", ieta),
                                           Form("Tight XJGamma %.1f < eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]),
                                           300, 0, 3.0, n_pT_bins, pT_bin_edges));

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

            h_background_iso_truth_reco[ieta].push_back(new TH2D(Form("h_background_iso_truth_reco_%d_%d", ieta, ipt),
                                                              Form("Background Truth Iso ET %.1f < eta < %.1f %.1f < pT < %.1f", eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                              300, 0, 30, 400, -10, 30));
            
            h_response_isoET[ieta].push_back(new TH2D(Form("h_response_isoET_%d_%d", ieta, ipt),
                                                      Form("Response Iso ET %.1f < eta < %.1f %.1f < pT < %.1f", eta_bins[ieta], eta_bins[ieta + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                                                      150, 0, 1.5, 400, -10, 30));
        }
    }

    int run_min = configYaml["analysis"]["run_min"].as<int>(-1);
    int run_max = configYaml["analysis"]["run_max"].as<int>(-1);
    if (run_min > 0 || run_max > 0)
        std::cout << "Run range filter: [" << run_min << ", " << run_max << "]" << std::endl;

    // Optional run list: if set, only events from runs in the list are processed.
    std::set<int> allowed_runnumbers;
    std::string run_list_file = configYaml["analysis"]["run_list_file"].as<std::string>("");
    if (!run_list_file.empty())
    {
        std::ifstream rlf(run_list_file);
        if (!rlf.is_open())
        {
            std::cerr << "[RunList] ERROR: cannot open run list file: " << run_list_file << std::endl;
            return;
        }
        int rn;
        while (rlf >> rn)
            allowed_runnumbers.insert(rn);
        std::cout << "[RunList] Loaded " << allowed_runnumbers.size()
                  << " runs from " << run_list_file << std::endl;
    }

    TRandom3 *rand = new TRandom3(0);
    std::set<int> skiprunnumbers = {0};
    int nentries = chain.GetEntries();
    int ientry = 0;

    while (reader.Next())
    {
        // Event-by-event MC vertex reweighting.
        // We directly update `weight` so all existing Fill(..., weight) calls use the per-event weight.
        weight = cross_weight;
        vertex_weight = 1.0;
        if (issim)
        {
            if (vertex_reweight_on && !do_vertex_scan)
            {
                if (!h_vertex_reweight)
                {
                    std::cerr << "[VertexReweight] ERROR: vertex reweighting is enabled but histogram is not loaded."
                              << std::endl;
                    return;
                }
                // Linearly interpolate between bin centers for smooth reweight values.
                // For events outside the histogram range (|vertexz| > histogram extent),
                // set weight to 0 — these events are outside the data vertex acceptance.
                const float vtx_lo = h_vertex_reweight->GetXaxis()->GetXmin();
                const float vtx_hi = h_vertex_reweight->GetXaxis()->GetXmax();
                if (*vertexz < vtx_lo || *vertexz > vtx_hi)
                {
                    vertex_weight = 0.0;
                }
                else
                {
                    // Interpolate handles values between bin centers; for values inside
                    // the first/last bin (between edge and center), it falls back to
                    // GetBinContent. Safe for our use case.
                    vertex_weight = h_vertex_reweight->Interpolate(*vertexz);
                }
            }

            // Fallback: drop events with non-finite or negative weight (data has
            // essentially zero events at this vertex position, so MC events here
            // should NOT contribute). Previously used 1.0 which kept raw MC events
            // at large |z| and biased h_truth_pT.
            if (!std::isfinite(vertex_weight) || vertex_weight < 0.0)
            {
                vertex_weight = 0.0;
            }

            weight *= vertex_weight;

            // Truth-vertex reweight: single MC → w(z_h); double MC → w(z_h)·w(z_mb).
            if (truth_vertex_reweight_on && !do_vertex_scan)
            {
                weight *= vertexz_truth_mb_ptr
                    ? TruthVertexWeight(h_truth_vtx_reweight, *vertexz_truth, **vertexz_truth_mb_ptr)
                    : TruthVertexWeight(h_truth_vtx_reweight, *vertexz_truth);
            }
        }

        if (ientry < 0)
        {
            ientry++;
            continue;
        }
        //std::cout<<"ientry: "<<ientry<<std::endl;
        // new calib entry
        // if (ientry == 16095884)
        //     continue;

        // if (ientry == 22626876)
        //     continue;

        // std::cout<<ientry<<std::endl;
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        
        //std::cout<<"vertexz: "<<*vertexz<<std::endl;


        


        if (!issim)
        {
            if (skiprunnumbers.find(*runnumber) != skiprunnumbers.end())
            {
                ientry++;
                continue;
            }

            bool in_run_range = (run_min <= 0 || *runnumber >= run_min) &&
                                (run_max <= 0 || *runnumber <= run_max);
            if (!in_run_range)
            {
                ientry++;
                continue;
            }

            if (!allowed_runnumbers.empty() &&
                allowed_runnumbers.find(*runnumber) == allowed_runnumbers.end())
            {
                ientry++;
                continue;
            }

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
                ientry++;
                continue;
            }
            /*
            float trigger_weight = 1.0;
            // If you later want a prescale-based weight, decide how to combine
            // prescales when multiple triggers are accepted (min/max/first-fired/etc).
            // For now we keep the old single-trigger logic as an example.
            if (!trigger_used.empty())
            {
                const int itrig = trigger_used.front();
                if (itrig >= 0 && (unsigned int)itrig < trigger_prescale.GetSize() && trigger_prescale[itrig] != 0)
                {
                    trigger_weight = trigger_prescale[itrig];
                }
                else
                {
                    std::cout << "Warning: trigger prescale is 0 or trigger index out of range" << std::endl;
                }
            }
            weight = trigger_weight;
            // check for nan and inf
            if (std::isnan(weight) || std::isinf(weight))
            {
                std::cout << "Warning: weight is nan or inf" << std::endl;
                if (!trigger_used.empty())
                {
                    const int itrig = trigger_used.front();
                    if (itrig >= 0 && (unsigned int)itrig < trigger_prescale.GetSize())
                    {
                        std::cout << trigger_prescale[itrig] << std::endl;
                    }
                }
                continue;
            }
            */
        }

        // Scan mode: fill vertex histogram and skip all downstream processing
        if (do_vertex_scan)
        {
            h_vertexz->Fill(*vertexz, mix_weight);
            ientry++;
            continue;
        }

        std::map<int, int> particle_trkidmap;
        // map for photon reco eff
        std::map<int, bool> photon_converts;
        // map for photon reco eff
        std::map<int, bool> photon_reco;
        // map for iso eff
        std::map<int, bool> photon_iso;
        std::map<int, float> photon_iso_ET;
        std::map<int, float> all_iso_ET;
        // map for id eff
        std::map<int, bool> photon_id;
        // map for n cluster per truth photon
        std::map<int, int> photon_ncluster;
        if (issim)
        {
            float maxphotonpT = 0;
            int maxphotonclass = 0;
            for (int iparticle = 0; iparticle < *nparticles; iparticle++)
            {
                particle_trkidmap[particle_trkid[iparticle]] = iparticle;
                // if (!isbackground)
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
                            // we check for conversion here
                            photon_converts[iparticle] = (particle_converted[iparticle] > 0);
                            // initialize reco to false
                            photon_reco[iparticle] = false;
                            // initialize id to false
                            photon_id[iparticle] = false;
                            // initialize iso to false
                            photon_iso[iparticle] = false;

                            photon_ncluster[iparticle] = 0;

                            // fill the truth histogram
                            float photonpT = particle_Pt[iparticle];
                            float photon_eta = particle_Eta[iparticle];

                            int etabin = -1;
                            for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
                            {
                                if ((photon_eta > eta_bins[ieta]) && (photon_eta < eta_bins[ieta + 1]))
                                {
                                    etabin = ieta;
                                    break;
                                }
                            }
                            if (etabin == -1)
                            {
                                continue;
                            }
                        }
                        photon_iso_ET[iparticle] = truthisoET;
                    }
                }
                all_iso_ET[iparticle] = truthisoET;
            }
            if (!isbackground)
            {
                if ((maxphotonpT > max_photon_upper) || (maxphotonpT < max_photon_lower))
                {
                    ientry++;
                    continue;
                }
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
                if (maxjetpT == 0) continue;

                if ((maxjetpT > max_jet_upper) )
                {
                    ientry++;
                    continue;
                }
                if (maxjetpT < max_jet_lower){
                    ientry++;
                    continue;
                }

                // energy scale cut for now
                /*
                if ((*energy_scale > energy_scale_upper) || (*energy_scale < energy_scale_lower))
                {
                    continue;
                }
                */
                //std::cout<<"accepted maxjetpT: "<<maxjetpT<<std::endl;
                h_max_truth_jet_pT->Fill(maxjetpT, weight);
            }
            //std::cout<<"at line 1350"<<std::endl;
            if (!isbackground)
            {
                for (int iparticle = 0; iparticle < *nparticles; iparticle++)
                {
                    particle_trkidmap[particle_trkid[iparticle]] = iparticle;

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

                        if (particle_photonclass[iparticle] < 3) // direct or fragmentation
                        {
                            if (truthisoET < truthisocut)
                            {
                                float photonpT = particle_Pt[iparticle];
                                float photon_eta = particle_Eta[iparticle];

                                int etabin = -1;
                                for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
                                {
                                    if ((photon_eta > eta_bins[ieta]) && (photon_eta < eta_bins[ieta + 1]))
                                    {
                                        etabin = ieta;
                                        break;
                                    }
                                }
                                if (etabin == -1)
                                {
                                    continue;
                                }

                                h_truth_pT[etabin]->Fill(photonpT, weight);
                                h_truth_pT_novtx[etabin]->Fill(photonpT, cross_weight);

                                // check if truth vertex is within the vertex cut
                                if (abs(*vertexz_truth) < vertexcut_truth)
                                {
                                    h_truth_pT_vertexcut[etabin]->Fill(photonpT, weight);
                                    if (*mbdnorthhit >= 1 && abs(*vertexz) < vertexcut)
                                    {
                                        h_truth_pT_vertexcut_mbd_north_cut[etabin]->Fill(photonpT, weight);
                                    }
                                    if (*mbdsouthhit >= 1 && abs(*vertexz) < vertexcut)
                                    {
                                        h_truth_pT_vertexcut_mbd_south_cut[etabin]->Fill(photonpT, weight);
                                    }
                                    if (*mbdnorthhit >= 1 && *mbdsouthhit >= 1 && abs(*vertexz) < vertexcut)
                                    {
                                        h_truth_pT_vertexcut_mbd_cut[etabin]->Fill(photonpT, weight);
                                    }
                                    if (*mbdnorthhit >= 1 && *mbdsouthhit < 1)
                                    {
                                        h_truth_pT_vertexcut_mbd_only_north[etabin]->Fill(photonpT, weight);
                                    }
                                    if (*mbdsouthhit >= 1 && *mbdnorthhit < 1)
                                    {
                                        h_truth_pT_vertexcut_mbd_only_south[etabin]->Fill(photonpT, weight);
                                    }
                                    if (*mbdnorthhit < 1 && *mbdsouthhit < 1)
                                    {
                                        h_truth_pT_vertexcut_mbd_neither[etabin]->Fill(photonpT, weight);
                                    }

                                    //if (abs(*vertexz) < vertexcut)
                                    //{
                                    //    h_truth_pT_vertexcut_mbd_cut[etabin]->Fill(photonpT, weight);
                                    //}
                                }
                            }
                        }
                    }
                }
            }
        }
        //std::cout<<"at line 1442"<<std::endl;
        // check vertex cut
        //std::cout<<"vertexz: "<<*vertexz<<std::endl;
        if (abs(*vertexz) > vertexcut)
        {
            ientry++;
            continue;
        }
        int vertex_bin = -1;
        for (int ivtx = 0; ivtx < (int)vertex_bins.size() - 1; ivtx++)
        {
            if ((abs(*vertexz) >= vertex_bins[ivtx]) && (abs(*vertexz) < vertex_bins[ivtx + 1]))
            {
                vertex_bin = ivtx;
                break;
            }
        }
        // Guard: skip events whose |vertexz| falls outside all vertex_bins ranges.
        // This can happen when vertex_cut > vertex_bins.back() or at the exact boundary.
        if (vertex_bin < 0)
        {
            ientry++;
            continue;
        }

        //std::cout<<"accepted vertexz: "<<*vertexz<<std::endl;
        //std::cout<<"mbdnorthhit: "<<*mbdnorthhit<<std::endl;
        //std::cout<<"mbdsouthhit: "<<*mbdsouthhit<<std::endl;
        //if (!(*mbdnorthhit >= 1 && *mbdsouthhit >= 1))
        //{
        //    ientry++;
        //    continue;
        //}

        // ------------------------------------------------------------------
        // Event-level MBD channel-time std dev, using channels with q > 0.1
        // ------------------------------------------------------------------
        auto compute_std_time = [](const TTreeReaderArray<float>& t,
                                   const TTreeReaderArray<float>& q,
                                   const float qmin) -> double {
            const int n = std::min((int)t.GetSize(), (int)q.GetSize());
            if (n <= 0) return 0.0;
            double sum = 0.0;
            double sumsq = 0.0;
            int nsel = 0;
            for (int i = 0; i < n; ++i)
            {
                if (q[i] <= qmin) continue;
                const double ti = t[i];
                if (!std::isfinite(ti)) continue;
                sum += ti;
                sumsq += ti * ti;
                ++nsel;
            }
            if (nsel < 2) return 0.0;
            const double mean = sum / nsel;
            // population variance: E[x^2] - mean^2
            const double var = std::max(0.0, (sumsq / nsel) - mean * mean);
            return std::sqrt(var);
        };

        const double mbd_std_north = compute_std_time(mbdnortht, mbdnorthq, 0.1f);
        const double mbd_std_south = compute_std_time(mbdsoutht, mbdsouthq, 0.1f);
        double mbd_avg_sigma = 0.0;
        if (mbd_std_north > 0.0 || mbd_std_south > 0.0)
        {
            h_mbd_time_std_north_vs_south->Fill(mbd_std_north, mbd_std_south, weight);
            if (mbd_std_north > 0.0 && mbd_std_south > 0.0)
            {
                mbd_avg_sigma = 0.5 * (mbd_std_north + mbd_std_south);
            }
            else
            {
                mbd_avg_sigma = std::max(mbd_std_north, mbd_std_south);
            }
        }
        if (std::isfinite(mbd_avg_sigma))
        {
            if (mbd_avg_sigma < mbd_avg_sigma_min || mbd_avg_sigma > mbd_avg_sigma_max)
            {
                ientry++;
                continue;
            }
        }

        float mbd_mean_time = -999;
        {
            //mbd_mean_time = (*mbdnorthtmean * *mbdnorthhit + *mbdsouthtmean * *mbdsouthhit) / (*mbdnorthhit + *mbdsouthhit);

            //// now check if MBD time is available, and if so, apply the stricter cut
            //double mbdoffset = -2.07;
            //if (*runnumber > 49375 && *runnumber < 49700) mbdoffset += +2.5;
            //if (*runnumber > 48600 && *runnumber < 49375) mbdoffset += -7.5;
            //if (*runnumber > 48180 && *runnumber < 48270) mbdoffset += -7.5;
            //if (*runnumber > 48050 && *runnumber < 48090) mbdoffset += -7.5;

            mbd_mean_time = *mbd_time;
            float mbdoffset =0;
            if(mbd_t0_correction.find(*runnumber) != mbd_t0_correction.end())
            {
                mbdoffset = mbd_t0_correction[*runnumber];
            }

            if(issim)
            {
                mbdoffset = 0.0;
            }

            mbd_mean_time =  mbd_mean_time - mbdoffset;
            
        }
        //std::cout<<"at line 1459"<<std::endl;
        // vertex_weight is already computed above for MC (and equals 1 for data).

        std::vector<float> jetphi;

        float maxrecojetpT = 0;
        for (int ijet = 0; ijet < *njet; ijet++)
        {
            h_reco_jet_pT->Fill(jet_Pt[ijet], weight);
            if (jet_Pt[ijet] > maxrecojetpT) maxrecojetpT = jet_Pt[ijet];

            if (jet_Pt[ijet] < 10)
                continue;
            jetphi.push_back(jet_Phi[ijet]);
        }
        if (maxrecojetpT > 0) h_max_reco_jet_pT->Fill(maxrecojetpT, weight);

        h_vertexz->Fill(*vertexz, weight);
        float leading_common_cluster_ET = 0;
        int leading_cluster_ET_index = -1;
        // Tower-mask veto helper: true if the cluster center lands on a
        // tower flagged in the mask (applied symmetrically to data/MC).
        auto is_tower_masked = [&](int icluster) -> bool {
            if (!h_tower_mask) return false;
            int ieta = (int) cluster_ietacent[icluster];
            int iphi = (int) cluster_iphicent[icluster];
            if (ieta < 0 || ieta >= h_tower_mask->GetNbinsX()) return false;
            if (iphi < 0 || iphi >= h_tower_mask->GetNbinsY()) return false;
            return h_tower_mask->GetBinContent(ieta + 1, iphi + 1) > 0;
        };

        float leading_cluster_ET = 0;
        for (int icluster = 0; icluster < *ncluster; icluster++)
        {
            if (is_tower_masked(icluster)) continue;
            if (issim)
            {
                cluster_Et[icluster] = cluster_Et[icluster] * clusterescale;
                if (clustereres > 0)
                {
                    cluster_Et[icluster] = cluster_Et[icluster] * rand->Gaus(1, clustereres);
                }
            }
            // need ET > 10 GeV
            if (cluster_Et[icluster] < reco_min_ET)
                continue;
            if (nosat)
            {
                if (cluster_nsaturated[icluster] > 0)
                    continue;
            }
            if (cluster_Et[icluster] > leading_common_cluster_ET)
            {
                leading_cluster_ET = cluster_Et[icluster];
                leading_cluster_ET_index = icluster;
            }
        }

        for (int icluster = 0; icluster < *ncluster; icluster++)
        {
            if (is_tower_masked(icluster)) continue;
            bool is_leading_cluster = (icluster == leading_cluster_ET_index);
            // need ET > 10 GeV
            if (cluster_Et[icluster] < reco_min_ET)
                continue;
            
            // for jet 10 event we want to remove some high ET clusters to reduce the fluctuation
            if (isbackground)
            {
                //std::cout<<"cluster_Et[icluster]: "<<cluster_Et[icluster]<<" cluster_ET_upper: "<<cluster_ET_upper<<std::endl;
                if (cluster_Et[icluster] > cluster_ET_upper)
                {
                    continue;
                }
            }
            if (nosat)
            {
                if (cluster_nsaturated[icluster] > 0)
                    continue;
            }

            // Tower acceptance map -- preselect level (ET gate + isbackground
            // ET-upper veto + saturation already applied above; no shape cuts).
            if (cluster_Et[icluster] >= tower_map_et_min)
            {
                h_etaphi_tower_preselect->Fill(
                    cluster_ietacent[icluster], cluster_iphicent[icluster], weight);
            }

            //std::cout<<"cluster_Et[icluster]: "<<cluster_Et[icluster]<<std::endl;
            float rhad22 = (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]) / (cluster_Et[icluster] + (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]));
            float rhad33 = (cluster_ihcal_et33[icluster] + cluster_ohcal_et33[icluster]) / (cluster_Et[icluster] + (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]));

            float reta77 = cluster_e37[icluster] / cluster_e77[icluster];
            float rphi77 = cluster_e73[icluster] / cluster_e77[icluster];

            float reta55 = cluster_e35[icluster] / cluster_e55[icluster];
            float rphi55 = cluster_e53[icluster] / cluster_e55[icluster];

            float reta = cluster_e33[icluster] / cluster_e73[icluster];

            float rphi = cluster_e33[icluster] / cluster_e37[icluster];

            float re11_E = cluster_e11[icluster] / cluster_E[icluster];

            float e11_over_e33 = cluster_e11[icluster] / cluster_e33[icluster];

            float e32_over_e35 = cluster_e32[icluster] / cluster_e35[icluster];

            float wr_cogx = cluster_wphi_cogx[icluster] / cluster_weta_cogx[icluster];

            // reco cut
            float recoisoET = -999;
            if (use_topo_iso == 2)
            {
                recoisoET = cluster_iso_topo_04[icluster];
                // Inner-R subtraction on topo path (iso_topo_innerr == 0.0 => no subtraction)
                if (iso_topo_innerr > 0.04 && iso_topo_innerr < 0.06)
                    recoisoET -= cluster_iso_topo_005[icluster];
                else if (iso_topo_innerr > 0.07 && iso_topo_innerr < 0.08)
                    recoisoET -= cluster_iso_topo_0075[icluster];
                else if (iso_topo_innerr > 0.09 && iso_topo_innerr < 0.11)
                    recoisoET -= cluster_iso_topo_01[icluster];
                else if (iso_topo_innerr > 0.19 && iso_topo_innerr < 0.21)
                    recoisoET -= cluster_iso_topo_02[icluster];
            }
            else if (use_topo_iso == 1)
            {
                recoisoET = cluster_iso_topo_03[icluster];
            }
            else if (use_topo_iso == 3)
            {
                // Tower-based donut: 70 MeV threshold outer cone R=0.3
                // (emcal + hcalin + hcalout) minus EMCal-only inner ring at 70 MeV.
                // iso_emcalinnerr must be one of {0.0, 0.05, 0.075, 0.10, 0.20}.
                recoisoET = cluster_iso_03_70_emcal[icluster]
                          + cluster_iso_03_70_hcalin[icluster]
                          + cluster_iso_03_70_hcalout[icluster];
                if (iso_emcalinnerr > 0.04 && iso_emcalinnerr < 0.06)
                    recoisoET -= cluster_iso_005_70_emcal[icluster];
                else if (iso_emcalinnerr > 0.07 && iso_emcalinnerr < 0.08)
                    recoisoET -= cluster_iso_0075_70_emcal[icluster];
                else if (iso_emcalinnerr > 0.09 && iso_emcalinnerr < 0.11)
                    recoisoET -= cluster_iso_01_70_emcal[icluster];
                else if (iso_emcalinnerr > 0.19 && iso_emcalinnerr < 0.21)
                    recoisoET -= cluster_iso_02_70_emcal[icluster];
            }
            else if (use_topo_iso == 4)
            {
                // All-calo 120 MeV donut: outer R=0.4 minus all-calo inner {0.05, 0.075, 0.20}.
                recoisoET = cluster_iso_04[icluster];
                if (iso_emcalinnerr > 0.04 && iso_emcalinnerr < 0.06)
                    recoisoET -= cluster_iso_005[icluster];
                else if (iso_emcalinnerr > 0.07 && iso_emcalinnerr < 0.08)
                    recoisoET -= cluster_iso_0075[icluster];
                else if (iso_emcalinnerr > 0.19 && iso_emcalinnerr < 0.21)
                    recoisoET -= cluster_iso_02[icluster];
            }
            else if (use_topo_iso == 5)
            {
                // All-calo 120 MeV donut with EMCal tower-ownership exclusion.
                recoisoET = cluster_iso_excl_04[icluster];
                if (iso_emcalinnerr > 0.04 && iso_emcalinnerr < 0.06)
                    recoisoET -= cluster_iso_excl_005[icluster];
                else if (iso_emcalinnerr > 0.07 && iso_emcalinnerr < 0.08)
                    recoisoET -= cluster_iso_excl_0075[icluster];
                else if (iso_emcalinnerr > 0.09 && iso_emcalinnerr < 0.11)
                    recoisoET -= cluster_iso_excl_01[icluster];
                else if (iso_emcalinnerr > 0.19 && iso_emcalinnerr < 0.21)
                    recoisoET -= cluster_iso_excl_02[icluster];
            }
            else if (conesize == 4)
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

            if (!use_topo_iso && iso_threshold)
            {
                if (iso_hcalonly)
                {
                    recoisoET = cluster_iso_03_70_hcalin[icluster] + cluster_iso_03_70_hcalout[icluster];
                }
                else
                {
                    //recoisoET = cluster_iso_03_70_emcal[icluster] + cluster_iso_03_70_hcalin[icluster] + cluster_iso_03_70_hcalout[icluster];
                    recoisoET = cluster_iso_03_120_emcal[icluster] + cluster_iso_03_120_hcalin[icluster] + cluster_iso_03_120_hcalout[icluster];
                    //recoisoET = cluster_iso_03_60_emcal[icluster] + cluster_iso_03_60_hcalin[icluster] + cluster_iso_03_60_hcalout[icluster];
                    if (iso_emcalinnerr > 0.04 && iso_emcalinnerr < 0.06)
                    {
                        //std::cout<<"using 0.05 emcal inner"<<std::endl;
                        recoisoET -= cluster_iso_005_70_emcal[icluster];
                    }
                    else if (iso_emcalinnerr > 0.09 && iso_emcalinnerr < 0.11)
                    {
                        //std::cout<<"using 0.1 emcal inner"<<std::endl;
                        recoisoET -= cluster_iso_01_70_emcal[icluster];
                    }
                    else if (iso_emcalinnerr > 0.19 && iso_emcalinnerr < 0.21)
                    {
                        //std::cout<<"using 0.2 emcal inner"<<std::endl;
                        recoisoET -= cluster_iso_02_70_emcal[icluster];
                    }
                    //if(recoisoET < 0)
                    //{
                    //    std::cout<<"recoisoET < 0, icluster: "<<icluster<<" recoisoET: "<<recoisoET<<std::endl;
                    //}

                }
            }
            // fudge the MC isoET
            if (issim)
            {
                recoisoET = recoisoET * mc_iso_scale;
                recoisoET += mc_iso_shift;
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

            // calculate cluster time
            float clusteravgtime = 0;
            float cluster_total_e = 0;
            for (int i = 0; i < 49; i++)
            {

                if (cluster_ownership_array[icluster * 49 + i] == 1)
                {
                    int status = cluster_status_array[icluster * 49 + i];
                    //check bit 5 is set then it is ZS tower
                    if (status & (1 << 5))
                    {
                        continue;
                    }
                    clusteravgtime += cluster_time_array[icluster * 49 + i] * cluster_e_array[icluster * 49 + i];
                    // std::cout<<"cluster_time_array[icluster][i]: "<<cluster_time_array[icluster * 49 + i]<<std::endl;
                    cluster_total_e += cluster_e_array[icluster * 49 + i];
                }
            }
            clusteravgtime = cluster_total_e > 0 ? clusteravgtime / cluster_total_e : 0;
            clusteravgtime = clusteravgtime * TIME_SAMPLE_NS;

            bool badtime = clusteravgtime < bg_timing_cut;

            // Cluster-MBD time window cut (data only)
            if (!issim && mbd_mean_time > -990)
            {
                float delta_t = clusteravgtime - mbd_mean_time;
                if (delta_t < cluster_mbd_time_min || delta_t > cluster_mbd_time_max)
                    continue;
            }

            bool otherside_jet = false;

            float max_b2bjet_pT = -1;

            float max_b2btruthjet_pT = -1;

            float b2bjet_dphi = 3 * M_PI / 4;

            float jet_eta = 0.6;

            for (int ijet = 0; ijet < *njet; ijet++)
            {
                float dphi = cluster_Phi[icluster] - jet_Phi[ijet];

                while (dphi > M_PI)
                    dphi = dphi - 2 * M_PI;
                while (dphi < -M_PI)
                    dphi = dphi + 2 * M_PI;

                if (abs(dphi) > (M_PI / 2))
                {
                    otherside_jet = true;
                    // break;
                }
                float calibrated_jet_pT = jet_Pt[ijet];
                if (abs(jet_Eta[ijet]) < jet_eta)
                {
                    if (abs(dphi) > b2bjet_dphi)
                    {
                        if (calibrated_jet_pT > max_b2bjet_pT && calibrated_jet_pT > 5)
                        {
                            max_b2bjet_pT = calibrated_jet_pT;
                        }
                    }
                }
            }
            if (issim)
            {
                for (int ijet = 0; ijet < *njet_truth; ijet++)
                {
                    float dphi = cluster_Phi[icluster] - jet_truth_Phi[ijet];

                    while (dphi > M_PI)
                        dphi = dphi - 2 * M_PI;
                    while (dphi < -M_PI)
                        dphi = dphi + 2 * M_PI;

                    if (abs(jet_truth_Eta[ijet]) < jet_eta)
                    {
                        if (abs(dphi) > b2bjet_dphi)
                        {
                            if (jet_truth_Pt[ijet] > max_b2btruthjet_pT && jet_truth_Pt[ijet] > 5)
                            {
                                max_b2btruthjet_pT = jet_truth_Pt[ijet];
                            }
                        }
                    }
                }
            }
            bool passes_common_b2bjet = (!common_b2bjet_cut) || (max_b2bjet_pT >= common_b2bjet_pt_min);
            float xjgamma = -1;
            float xjgamma_truthjet = -1;

            if (is_leading_cluster)
            {
                xjgamma = max_b2bjet_pT / leading_cluster_ET;
                xjgamma_truthjet = max_b2btruthjet_pT / leading_cluster_ET;
            }

            bool isnpb = badtime && !otherside_jet;
            if (cluster_weta_cogx[icluster] < npb_weta_min)
            {
                isnpb = false;
            }

            bool common_pass = false;
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

            // common cuts
            bool passes_common_shape =
                cluster_prob[icluster] > common_prob_min &&
                cluster_prob[icluster] < common_prob_max &&
                e11_over_e33 > common_e11_over_e33_min &&
                e11_over_e33 < common_e11_over_e33_max &&
                wr_cogx > common_wr_cogx_bound &&
                (cluster_weta_cogx[icluster] < common_cluster_weta_cogx_bound) &&
                (!common_npb_cut_on || cluster_npb_score[icluster] > common_npb_score_cut);

            if (passes_common_shape && passes_common_b2bjet)
            {
                common_pass = true;
            }

            // Tower acceptance map -- common level
            if (common_pass && cluster_Et[icluster] >= tower_map_et_min)
            {
                h_etaphi_tower_common->Fill(
                    cluster_ietacent[icluster], cluster_iphicent[icluster], weight);
            }
            //std::cout<<"cluster_prob[icluster]: "<<cluster_prob[icluster]<<std::endl;
            //std::cout<<"e11_over_e33: "<<e11_over_e33<<std::endl;
            //std::cout<<"cluster_weta_cogx[icluster]: "<<cluster_weta_cogx[icluster]<<std::endl;
            //std::cout<<"common_pass: "<<common_pass<<std::endl;

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
                
                // Select BDT model based on cluster ET (ET-binned or fallback)
                std::string selected_bdt_model = bdt_model_name;
                if (use_et_binned_bdt) {
                    float et = cluster_Et[icluster];
                    for (int ib = 0; ib < (int)bdt_et_bin_models.size(); ++ib) {
                        if (et >= bdt_et_bin_edges[ib] && et < bdt_et_bin_edges[ib + 1]) {
                            selected_bdt_model = bdt_et_bin_models[ib];
                            break;
                        }
                    }
                }
                float bdt_score = (*bdt_arrays[selected_bdt_model])[icluster];

                //std::cout<<"bdt_score: "<<bdt_score<<std::endl;
                float tight_bdt_min_et = tight_bdt_min_slope * cluster_Et[icluster] + tight_bdt_min_intercept;
                bool is_bdt_tight =
                    (bdt_score > tight_bdt_min_et) &&
                    (bdt_score < tight_bdt_max);

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

                // Tower acceptance maps -- tight and tight+iso levels
                if (tight && cluster_Et[icluster] >= tower_map_et_min)
                {
                    h_etaphi_tower_tight->Fill(
                        cluster_ietacent[icluster], cluster_iphicent[icluster], weight);
                    if (iso)
                    {
                        h_etaphi_tower_tight_iso->Fill(
                            cluster_ietacent[icluster], cluster_iphicent[icluster], weight);
                    }
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
                    bdt_score > non_tight_bdt_min_slope * cluster_Et[icluster] + non_tight_bdt_min_intercept &&
                    bdt_score < non_tight_bdt_max_slope * cluster_Et[icluster] + non_tight_bdt_max_intercept)
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
            }

            if (tight && iso)
            {
                h_tight_iso_cluster[etabin]->Fill(cluster_Et[icluster], weight);
                h_tight_iso_xjgamma[etabin]->Fill(xjgamma, cluster_Et[icluster], weight);
                h_delta_t_tight_iso_cluster[etabin]->Fill(cluster_Et[icluster], clusteravgtime - mbd_mean_time, weight);
                // h_pT_reco_response[etabin]->Fill(cluster_Et[icluster]);
                //std::cout<<"vertex_bin: "<<vertex_bin<<std::endl;
                h_vertex_tight_iso_cluster_pt_eta[vertex_bin]->Fill(cluster_Et[icluster], cluster_Eta[icluster], weight);
                if (isnpb)
                    h_tight_iso_cluster_background[etabin]->Fill(cluster_Et[icluster], weight);
            }
            if (tight)
            {
                h_tight_cluster[etabin]->Fill(cluster_Et[icluster], weight);
                h_tight_xjgamma[etabin]->Fill(xjgamma, cluster_Et[icluster], weight);
                h_delta_t_tight_cluster[etabin]->Fill(cluster_Et[icluster], clusteravgtime - mbd_mean_time, weight);
                if(max_b2bjet_pT > 10){
                    h_delta_t_tight_cluster_b2bjet[etabin]->Fill(cluster_Et[icluster], clusteravgtime - mbd_mean_time, weight);
                }
            }
            if (tight && noniso)
            {
                h_tight_noniso_cluster[etabin]->Fill(cluster_Et[icluster], weight);
                h_tight_noniso_xjgamma[etabin]->Fill(xjgamma, cluster_Et[icluster], weight);
                h_delta_t_tight_noniso_cluster[etabin]->Fill(cluster_Et[icluster], clusteravgtime - mbd_mean_time, weight);
                h_vertex_tight_noniso_cluster_pt_eta[vertex_bin]->Fill(cluster_Et[icluster], cluster_Eta[icluster], weight);
                if (isnpb)
                    h_tight_noniso_cluster_background[etabin]->Fill(cluster_Et[icluster], weight);
            }
            if (nontight && iso)
            {
                h_nontight_iso_cluster[etabin]->Fill(cluster_Et[icluster], weight);
                h_nontight_iso_xjgamma[etabin]->Fill(xjgamma, cluster_Et[icluster], weight);
                h_delta_t_nontight_iso_cluster[etabin]->Fill(cluster_Et[icluster], clusteravgtime - mbd_mean_time, weight);
                if(max_b2bjet_pT > 10){
                    h_delta_t_nontight_iso_cluster_b2bjet[etabin]->Fill(cluster_Et[icluster], clusteravgtime - mbd_mean_time, weight);
                }
                h_vertex_nontight_iso_cluster_pt_eta[vertex_bin]->Fill(cluster_Et[icluster], cluster_Eta[icluster], weight);
                if (isnpb)
                    h_nontight_iso_cluster_background[etabin]->Fill(cluster_Et[icluster], weight);
            }
            if (nontight && noniso)
            {
                h_nontight_noniso_cluster[etabin]->Fill(cluster_Et[icluster], weight);
                h_nontight_noniso_xjgamma[etabin]->Fill(xjgamma, cluster_Et[icluster], weight);
                h_delta_t_nontight_noniso_cluster[etabin]->Fill(cluster_Et[icluster], clusteravgtime - mbd_mean_time, weight);
                if(max_b2bjet_pT > 10){
                    h_delta_t_nontight_noniso_cluster_b2bjet[etabin]->Fill(cluster_Et[icluster], clusteravgtime - mbd_mean_time, weight);
                }
                h_vertex_nontight_noniso_cluster_pt_eta[vertex_bin]->Fill(cluster_Et[icluster], cluster_Eta[icluster], weight);
                if (isnpb)
                    h_nontight_noniso_cluster_background[etabin]->Fill(cluster_Et[icluster], weight);
            }
            if (common_pass)
            {
                h_common_cluster[etabin]->Fill(cluster_Et[icluster], weight);
                h_common_xjgamma[etabin]->Fill(xjgamma, cluster_Et[icluster], weight);
                h_delta_t_common_cluster[etabin]->Fill(cluster_Et[icluster], clusteravgtime - mbd_mean_time, weight);
                h_cluster_common_Et->Fill(cluster_Et[icluster], weight);
                if (cluster_Et[icluster] > leading_common_cluster_ET)
                {
                    leading_common_cluster_ET = cluster_Et[icluster];
                }
            }
            h_all_cluster[etabin]->Fill(cluster_Et[icluster], weight);
            h_all_xjgamma[etabin]->Fill(xjgamma, cluster_Et[icluster], weight);
            h_delta_t_all_cluster[etabin]->Fill(cluster_Et[icluster], clusteravgtime - mbd_mean_time, weight);
            h_mbd_t_vs_cluster_t->Fill(mbd_mean_time, clusteravgtime, weight);
            if (isnpb)
            {
                h_delta_t_npb_cluster[etabin]->Fill(cluster_Et[icluster], clusteravgtime - mbd_mean_time, weight);
            }
            if (pTbin != -1)
            {
                if (tight)
                {
                    h_tight_cluster_pT[etabin][pTbin]->Fill(recoisoET, weight);
                }
                if (nontight)
                {
                    h_nontight_cluster_pT[etabin][pTbin]->Fill(recoisoET, weight);
                }
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
                    dphi -= 2 * M_PI;
                }
                else if (dphi < -M_PI)
                {
                    dphi += 2 * M_PI;
                }
                float dR = sqrt(deta * deta + dphi * dphi);
                //check the pid in each sideband region
                int pid = particle_pid[iparticle];
                if(common_pass)
                {
                    h_common_pid_pt[etabin]->Fill(pid, cluster_Et[icluster], weight);
                }
                if(tight && iso)
                {
                    h_tight_iso_pid_pt[etabin]->Fill(pid, cluster_Et[icluster], weight);
                }
                if(tight && noniso)
                {
                    h_tight_noniso_pid_pt[etabin]->Fill(pid, cluster_Et[icluster], weight);
                }
                if(nontight && iso)
                {
                    h_nontight_iso_pid_pt[etabin]->Fill(pid, cluster_Et[icluster], weight);
                }
                if(nontight && noniso)
                {
                    h_nontight_noniso_pid_pt[etabin]->Fill(pid, cluster_Et[icluster], weight);
                }


                // truth iso vs. reco iso

                if (photon_iso_ET.find(iparticle) != photon_iso_ET.end())
                {
                    if (pTbin != -1)
                    {
                        h_iso_truth_reco[etabin][pTbin]->Fill(photon_iso_ET[iparticle], recoisoET, weight);
                    }
                    if (photon_reco.find(iparticle) == photon_reco.end())
                    {
                        // then it is non truth signal, if it pass the reco, iso, and tight cuts, then it is a fake
                        // Note: deltaR cut removed — trkID match is sufficient for truth association
                        // (see efficiencytool/reports/deltaR_matching_study.tex)
                        if (iso && tight)
                        {
                            h_pT_reco_fake[etabin]->Fill(cluster_Et[icluster], weight);
                        }
                    }
                }

                // iparticle has to be in the photon map
                if (photon_reco.find(iparticle) == photon_reco.end())
                {
                    // this is non-photon background
                    if (pTbin != -1)
                    {
                        h_background_iso_truth_reco[etabin][pTbin]->Fill(all_iso_ET[iparticle], recoisoET, weight);
                    }
                    h_background_truth_isoET[etabin]->Fill(particle_Pt[iparticle], recoisoET, weight);
                    if (tight)    h_tight_recoisoET_background[etabin]->Fill(cluster_Et[icluster], recoisoET, weight);
                    if (nontight) h_nontight_recoisoET_background[etabin]->Fill(cluster_Et[icluster], recoisoET, weight);
                    // not matched to a direct or fragmentation photon — fill ABCD background yield
                    if (tight && iso)      h_tight_iso_cluster_notmatch[etabin]->Fill(cluster_Et[icluster], weight);
                    if (tight && noniso)   h_tight_noniso_cluster_notmatch[etabin]->Fill(cluster_Et[icluster], weight);
                    if (nontight && iso)   h_nontight_iso_cluster_notmatch[etabin]->Fill(cluster_Et[icluster], weight);
                    if (nontight && noniso) h_nontight_noniso_cluster_notmatch[etabin]->Fill(cluster_Et[icluster], weight);
                    continue;
                }

                photon_ncluster[iparticle]++;

                // deltaR cut removed — trkID match alone is sufficient for truth association.
                // The dR < 0.1 cut was rejecting ~3.5% of well-reconstructed photons due to
                // vertex resolution (dz = vtx_reco - vtx_truth), not genuine reconstruction failure.
                // See efficiencytool/reports/deltaR_matching_study.tex for the full study.
                // if (dR < eff_dR)  // OLD: geometric cut on top of trkID match
                {

                    // if(photon_converts[iparticle]) continue;

                    photon_reco[iparticle] = true;

                    h_pT_truth_reco[etabin]->Fill(particle_Pt[iparticle], cluster_Et[icluster] / particle_Pt[iparticle], weight);
                    if (pTbin != -1)
                    {
                        h_response_isoET[etabin][pTbin]->Fill(cluster_Et[icluster] / particle_Pt[iparticle], recoisoET, weight);
                    }
                    if (iso)
                    {
                        photon_iso[iparticle] = true;
                    }

                    if (tight)
                    {
                        photon_id[iparticle] = true;
                    }

                    if (particle_Pt[iparticle] > pTmin_truth && particle_Pt[iparticle] < pTmax_truth && cluster_Et[icluster] > pTmin && cluster_Et[icluster] < pTmax)
                    {
                        h_all_cluster_signal[etabin]->Fill(cluster_Et[icluster], weight);
                        // fill the max backtobjets
                        h_all_cluster_Et_max_b2bjet[etabin]->Fill(cluster_Et[icluster], max_b2bjet_pT, weight);
                    }

                    if (tight)
                    {
                        if (particle_Pt[iparticle] > pTmin_truth && particle_Pt[iparticle] < pTmax_truth && cluster_Et[icluster] > pTmin && cluster_Et[icluster] < pTmax)
                        {
                            h_tight_cluster_signal[etabin]->Fill(cluster_Et[icluster], weight);
                            h_tight_xjgamma_signal[etabin]->Fill(xjgamma, cluster_Et[icluster], weight);
                            if (iso)
                            {
                                h_tight_iso_cluster_signal[etabin]->Fill(cluster_Et[icluster], weight);
                                h_tight_iso_xjgamma_signal[etabin]->Fill(xjgamma, cluster_Et[icluster], weight);
                                h_tight_iso_truthjet_xjgamma_signal[etabin]->Fill(xjgamma_truthjet, cluster_Et[icluster], weight);
                                h_vertex_tight_iso_cluster_signal_pt_eta[vertex_bin]->Fill(cluster_Et[icluster], cluster_Eta[icluster], weight);
                                // fill the response matrix

                                float response_reweight = 1.0;
                                if (reweight)
                                {
                                    response_reweight = particle_Pt[iparticle] > 30 ? f_reweight->Eval(30) : f_reweight->Eval(particle_Pt[iparticle]);
                                }
                                h_pT_truth_response[etabin]->Fill(particle_Pt[iparticle], weight*response_reweight);
                                h_pT_reco_response[etabin]->Fill(cluster_Et[icluster], weight*response_reweight);
                                responses_full[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
                                h_response_full_list[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
                                if (ientry < (nentries / 2))
                                {
                                    h_pT_truth_half_response[etabin]->Fill(particle_Pt[iparticle], weight);
                                    h_pT_reco_half_response[etabin]->Fill(cluster_Et[icluster], weight);
                                    responses_half[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
                                    h_response_half_list[etabin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
                                }
                                else
                                {
                                    h_pT_truth_secondhalf_response[etabin]->Fill(particle_Pt[iparticle], weight);
                                    h_pT_reco_secondhalf_response[etabin]->Fill(cluster_Et[icluster], weight);
                                }
                            }
                        }
                    }
                    if (tight && noniso)
                    {
                        h_tight_noniso_cluster_signal[etabin]->Fill(cluster_Et[icluster], weight);
                        h_tight_noniso_xjgamma_signal[etabin]->Fill(xjgamma, cluster_Et[icluster], weight);
                        h_tight_noniso_truthjet_xjgamma_signal[etabin]->Fill(xjgamma_truthjet, cluster_Et[icluster], weight);
                        h_vertex_tight_noniso_cluster_signal_pt_eta[vertex_bin]->Fill(cluster_Et[icluster], cluster_Eta[icluster], weight);
                    }
                    if (nontight && iso)
                    {
                        h_nontight_iso_cluster_signal[etabin]->Fill(cluster_Et[icluster], weight);
                        h_nontight_iso_xjgamma_signal[etabin]->Fill(xjgamma, cluster_Et[icluster], weight);
                        h_nontight_iso_truthjet_xjgamma_signal[etabin]->Fill(xjgamma_truthjet, cluster_Et[icluster], weight);
                        h_vertex_nontight_iso_cluster_signal_pt_eta[vertex_bin]->Fill(cluster_Et[icluster], cluster_Eta[icluster], weight);
                    }
                    if (nontight && noniso)
                    {
                        h_nontight_noniso_cluster_signal[etabin]->Fill(cluster_Et[icluster], weight);
                        h_nontight_noniso_xjgamma_signal[etabin]->Fill(xjgamma, cluster_Et[icluster], weight);
                        h_nontight_noniso_truthjet_xjgamma_signal[etabin]->Fill(xjgamma_truthjet, cluster_Et[icluster], weight);
                        h_vertex_nontight_noniso_cluster_signal_pt_eta[vertex_bin]->Fill(cluster_Et[icluster], cluster_Eta[icluster], weight);
                    }

                    h_singal_reco_isoET[etabin]->Fill(cluster_Et[icluster], recoisoET, weight);
                    h_singal_truth_isoET[etabin]->Fill(particle_Pt[iparticle], recoisoET, weight);
                }
            }
        } // end of cluster loop
        if (leading_common_cluster_ET > 0)
        {
            h_cluster_common_leading_Et->Fill(leading_common_cluster_ET, weight);
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
            eff_converts_eta[etabin]->FillWeighted(photon_converts[it->first], weight, photon_pT);

            h_ncluster_truth[etabin]->Fill(photon_pT, photon_ncluster[it->first], weight);

            // std::cout<<"photon_pT: "<<photon_pT<<" photon_eta: "<<photon_eta<<" photon_converts: "<<photon_converts[it->first]<<std::endl;
            // std::cout<<"effieicncy: "<<eff_converts->GetEfficiency(0)<<std::endl;
            eff_reco_eta[etabin]->FillWeighted(photon_reco[it->first], weight, photon_pT);
            h_vertex_efficiency_denominator[vertex_bin]->Fill(photon_pT, photon_eta, weight);
            if(photon_reco[it->first]) h_vertex_efficiency_reco[vertex_bin]->Fill(photon_pT, photon_eta, weight);


            bool totalpass = photon_reco[it->first] && photon_iso[it->first] && photon_id[it->first];
            eff_all_eta[etabin]->FillWeighted(totalpass, weight, photon_pT);

            if (photon_reco[it->first])
            {
               if(photon_id[it->first]) h_vertex_efficiency_id[vertex_bin]->Fill(photon_pT, photon_eta, weight);
               if(photon_iso[it->first]) h_vertex_efficiency_iso[vertex_bin]->Fill(photon_pT, photon_eta, weight);
                // if(!photon_converts[it->first])
                {
                    eff_iso_eta[etabin]->FillWeighted(photon_iso[it->first], weight, photon_pT);
                }

                if (photon_iso[it->first])
                {
                    eff_id_eta[etabin]->FillWeighted(photon_id[it->first], weight, photon_pT);
                }
            }
        }
        ientry++;
    }

    // Scan mode: write only the vertex histogram and return early
    if (do_vertex_scan)
    {
        TFile* fscan = new TFile(vtxscan_outfilename.c_str(), "RECREATE");
        h_vertexz->Write();
        fscan->Close();
        std::cout << "[VertexScan] Written vertex histogram to " << vtxscan_outfilename << std::endl;
        return;
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
    if (!issim)
    {
        SaveYamlToRoot(fout, configname.c_str());
    }

    fresponse->Write();
    fresponse->Close();
}
