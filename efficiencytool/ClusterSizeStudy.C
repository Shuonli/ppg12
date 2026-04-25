#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <yaml-cpp/yaml.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "CrossSectionWeights.h"
#include "TruthVertexReweightLoader.h"

namespace {

constexpr int   kNtow       = 49;
constexpr int   kGridN      = 7;
constexpr float kEtaFidMax  = 0.7f;
constexpr int   kNptBins    = 40;
constexpr float kPtMin      = 0.0f;
constexpr float kPtMax      = 40.0f;

// Selection categories produced per cluster node.
// - "all"     : baseline (|eta|<0.7 only), always filled
// - "common"  : passing the common cut (shower-shape + NPB + weta bound)
// - "tight"   : passing the common cut plus the tight (ET-dependent BDT + ET-fraction) cut
//
// cut-application is driven by a YAML config file.  when configname is empty,
// only "all" is filled (back-compat with the no-cut baseline report).
struct NodeHists {
  TProfile* prof_n_owned[3]   = {nullptr, nullptr, nullptr};
  TProfile* prof_width_eta[3] = {nullptr, nullptr, nullptr};
  TProfile* prof_width_phi[3] = {nullptr, nullptr, nullptr};
  // 2D histograms: cluster-size metric vs ET. Y-axis binning is integer-centered
  // (0.5..49.5 for n_owned, 0.5..7.5 for width_eta/width_phi). X-axis matches the
  // TProfile binning so ProjectionY over ptRanges[] boundaries is exact.
  TH2F* h2_n_owned[3]   = {nullptr, nullptr, nullptr};
  TH2F* h2_width_eta[3] = {nullptr, nullptr, nullptr};
  TH2F* h2_width_phi[3] = {nullptr, nullptr, nullptr};
};

const char* kSelSuffix[3] = {"", "_common", "_tight"};

NodeHists MakeNodeHists(const std::string& key, bool have_cuts) {
  NodeHists nh;
  const int nsel = have_cuts ? 3 : 1;
  for (int s = 0; s < nsel; ++s) {
    const std::string suffix = kSelSuffix[s];
    nh.prof_n_owned[s] = new TProfile(
        Form("prof_n_owned_%s%s", key.c_str(), suffix.c_str()),
        ";E_{T} [GeV];#LTn_{owned}#GT",
        kNptBins, kPtMin, kPtMax);
    nh.prof_width_eta[s] = new TProfile(
        Form("prof_width_eta_%s%s", key.c_str(), suffix.c_str()),
        ";E_{T} [GeV];#LT#Deltai_{#eta}#GT",
        kNptBins, kPtMin, kPtMax);
    nh.prof_width_phi[s] = new TProfile(
        Form("prof_width_phi_%s%s", key.c_str(), suffix.c_str()),
        ";E_{T} [GeV];#LT#Deltai_{#phi}#GT",
        kNptBins, kPtMin, kPtMax);
    nh.h2_n_owned[s] = new TH2F(
        Form("h2_n_owned_%s%s", key.c_str(), suffix.c_str()),
        ";E_{T} [GeV];n_{owned}",
        kNptBins, kPtMin, kPtMax,
        kNtow, 0.5f, kNtow + 0.5f);
    nh.h2_width_eta[s] = new TH2F(
        Form("h2_width_eta_%s%s", key.c_str(), suffix.c_str()),
        ";E_{T} [GeV];#Deltai_{#eta}",
        kNptBins, kPtMin, kPtMax,
        kGridN, 0.5f, kGridN + 0.5f);
    nh.h2_width_phi[s] = new TH2F(
        Form("h2_width_phi_%s%s", key.c_str(), suffix.c_str()),
        ";E_{T} [GeV];#Deltai_{#phi}",
        kNptBins, kPtMin, kPtMax,
        kGridN, 0.5f, kGridN + 0.5f);
    nh.prof_n_owned[s]->Sumw2();
    nh.prof_width_eta[s]->Sumw2();
    nh.prof_width_phi[s]->Sumw2();
    nh.h2_n_owned[s]->Sumw2();
    nh.h2_width_eta[s]->Sumw2();
    nh.h2_width_phi[s]->Sumw2();
  }
  return nh;
}

void WriteNodeHists(const NodeHists& nh) {
  for (int s = 0; s < 3; ++s) {
    if (nh.prof_n_owned[s])   nh.prof_n_owned[s]->Write();
    if (nh.prof_width_eta[s]) nh.prof_width_eta[s]->Write();
    if (nh.prof_width_phi[s]) nh.prof_width_phi[s]->Write();
    if (nh.h2_n_owned[s])     nh.h2_n_owned[s]->Write();
    if (nh.h2_width_eta[s])   nh.h2_width_eta[s]->Write();
    if (nh.h2_width_phi[s])   nh.h2_width_phi[s]->Write();
  }
}

std::string MapInputFiletype(const std::string& filetype) {
  if (filetype == "photon10_nom") return "photon10";
  if (filetype == "jet12_nom")    return "jet12";
  return filetype;
}

std::string SampleFilePath(const std::string& input_filetype, bool issim) {
  if (!issim) {
    return "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/anatreemaker/"
           "macro_maketree/data/ana521/condorout/part_*_with_bdt_split.root";
  }
  return "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/FunWithxgboost/" +
         input_filetype + "/bdt_split.root";
}

// Configuration of per-node cut parameters loaded from the showershape config.
struct CutConfig {
  // common cut
  float common_prob_min            = 0.f;
  float common_prob_max            = 1.f;
  float common_e11_over_e33_min    = 0.f;
  float common_e11_over_e33_max    = 1.f;
  float common_weta_cogx_bound     = 2.f;
  int   npb_cut_on                 = 0;
  float npb_score_cut              = 0.5f;
  // tight cut (common-superset)
  float tight_weta_cogx_min        = 0.f;
  float tight_weta_cogx_max_b      = 1.f;
  float tight_weta_cogx_max_s      = 0.f;
  float tight_wphi_cogx_min        = 0.f;
  float tight_wphi_cogx_max_b      = 1.f;
  float tight_wphi_cogx_max_s      = 0.f;
  float tight_et1_min_b            = 0.f;
  float tight_et1_min_s            = 0.f;
  float tight_et1_max              = 1.f;
  float tight_et2_min              = 0.f;
  float tight_et2_max              = 1.f;
  float tight_et3_min              = 0.f;
  float tight_et3_max              = 1.f;
  float tight_et4_min              = 0.f;
  float tight_et4_max              = 1.f;
  float tight_e11_over_e33_min     = 0.f;
  float tight_e11_over_e33_max     = 1.f;
  float tight_e32_over_e35_min     = 0.f;
  float tight_e32_over_e35_max     = 1.f;
  float tight_prob_min             = 0.f;
  float tight_prob_max             = 1.f;
  float tight_bdt_min_intercept    = 0.86f;
  float tight_bdt_min_slope        = -0.02f;
  float tight_bdt_max              = 1.f;
  // BDT model name used for the tight BDT cut (nominal: base_v1E)
  std::string bdt_model_for_tight = "base_v1E";
};

CutConfig LoadCutConfig(const YAML::Node& cfg) {
  CutConfig c;
  auto A = cfg["analysis"];
  auto C = A["common"];
  c.common_prob_min         = C["prob_min"].as<float>(0.f);
  c.common_prob_max         = C["prob_max"].as<float>(1.f);
  c.common_e11_over_e33_min = C["e11_over_e33_min"].as<float>(0.f);
  c.common_e11_over_e33_max = C["e11_over_e33_max"].as<float>(1.f);
  c.common_weta_cogx_bound  = C["cluster_weta_cogx_bound"].as<float>(2.f);
  c.npb_cut_on              = C["npb_cut_on"].as<int>(0);
  c.npb_score_cut           = C["npb_score_cut"].as<float>(0.5f);

  auto T = A["tight"];
  c.tight_weta_cogx_min     = T["weta_cogx_min"].as<float>(0.f);
  c.tight_weta_cogx_max_b   = T["weta_cogx_max_b"].as<float>(1.f);
  c.tight_weta_cogx_max_s   = T["weta_cogx_max_s"].as<float>(0.f);
  c.tight_wphi_cogx_min     = T["wphi_cogx_min"].as<float>(0.f);
  c.tight_wphi_cogx_max_b   = T["wphi_cogx_max_b"].as<float>(1.f);
  c.tight_wphi_cogx_max_s   = T["wphi_cogx_max_s"].as<float>(0.f);
  c.tight_et1_min_b         = T["et1_min_b"].as<float>(0.f);
  c.tight_et1_min_s         = T["et1_min_s"].as<float>(0.f);
  c.tight_et1_max           = T["et1_max"].as<float>(1.f);
  c.tight_et2_min           = T["et2_min"].as<float>(0.f);
  c.tight_et2_max           = T["et2_max"].as<float>(1.f);
  c.tight_et3_min           = T["et3_min"].as<float>(0.f);
  c.tight_et3_max           = T["et3_max"].as<float>(1.f);
  c.tight_et4_min           = T["et4_min"].as<float>(0.f);
  c.tight_et4_max           = T["et4_max"].as<float>(1.f);
  c.tight_e11_over_e33_min  = T["e11_over_e33_min"].as<float>(0.f);
  c.tight_e11_over_e33_max  = T["e11_over_e33_max"].as<float>(1.f);
  c.tight_e32_over_e35_min  = T["e32_over_e35_min"].as<float>(0.f);
  c.tight_e32_over_e35_max  = T["e32_over_e35_max"].as<float>(1.f);
  c.tight_prob_min          = T["prob_min"].as<float>(0.f);
  c.tight_prob_max          = T["prob_max"].as<float>(1.f);
  c.tight_bdt_min_intercept = T["bdt_min_intercept"].as<float>(0.86f);
  c.tight_bdt_min_slope     = T["bdt_min_slope"].as<float>(-0.02f);
  c.tight_bdt_max           = T["bdt_max"].as<float>(1.f);

  // BDT model selection: the showershape config lists ET-binned models;
  // for the cluster-size cut study we use a single model (nominal base_v1E)
  // across the full ET range, matching the dominant PPG12 pipeline choice.
  if (auto m = cfg["input"]["bdt_et_bin_models"]) {
    auto vm = m.as<std::vector<std::string>>();
    if (!vm.empty()) c.bdt_model_for_tight = vm.front();
  }
  return c;
}

// Load a run-list file (one runnumber per line) into a set. Empty file or
// empty set means "no restriction".
std::set<int> LoadRunList(const std::string& path) {
  std::set<int> out;
  if (path.empty()) return out;
  std::ifstream f(path);
  if (!f.is_open()) {
    std::cerr << "[ClusterSizeStudy] WARN: run_list_file not found: " << path << std::endl;
    return out;
  }
  int r;
  while (f >> r) out.insert(r);
  std::cout << "[ClusterSizeStudy] loaded " << out.size() << " runs from " << path << std::endl;
  return out;
}

}  // namespace

void ClusterSizeStudy(const std::string& filetype,
                      const std::string& outpath_override = "",
                      float mix_weight = 1.0f,
                      const std::string& configname = "")
{
  const bool issim = (filetype != "data");

  // ------------------------------------------------------------------
  // Load YAML config (optional). When empty, behave as the baseline macro.
  // ------------------------------------------------------------------
  const bool have_cfg = !configname.empty();
  YAML::Node cfg;
  CutConfig cuts;
  float  vertex_cut            = 60.0f;
  int    run_min               = -1;
  int    run_max               = -1;
  std::set<int> allowed_runs;
  int    truth_vertex_reweight_on = 0;
  std::string truth_vertex_reweight_file;
  TH1*   h_truth_vtx_reweight  = nullptr;
  bool   need_vertexz_truth_mb = false;
  if (have_cfg) {
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    cfg = YAML::LoadFile(configname);
    cuts = LoadCutConfig(cfg);
    vertex_cut = cfg["analysis"]["vertex_cut"].as<float>(60.f);
    run_min = cfg["analysis"]["run_min"].as<int>(-1);
    run_max = cfg["analysis"]["run_max"].as<int>(-1);
    std::string run_list_file =
        cfg["analysis"]["run_list_file"].as<std::string>("");
    allowed_runs = LoadRunList(run_list_file);
    truth_vertex_reweight_on =
        cfg["analysis"]["truth_vertex_reweight_on"].as<int>(0);
    truth_vertex_reweight_file =
        cfg["analysis"]["truth_vertex_reweight_file"].as<std::string>("");
    if (issim && truth_vertex_reweight_on) {
      h_truth_vtx_reweight = LoadTruthVertexReweight(truth_vertex_reweight_file);
      if (!h_truth_vtx_reweight) {
        std::cerr << "[ClusterSizeStudy] FATAL: truth-vertex reweight requested but load failed"
                  << std::endl;
        return;
      }
    }
    std::cout << "[ClusterSizeStudy] cuts: bdt_model=" << cuts.bdt_model_for_tight
              << " npb_cut_on=" << cuts.npb_cut_on
              << " tight_bdt_slope=" << cuts.tight_bdt_min_slope
              << " tight_bdt_intercept=" << cuts.tight_bdt_min_intercept
              << std::endl;
    std::cout << "[ClusterSizeStudy] period: run_min=" << run_min
              << " run_max=" << run_max
              << " vertex_cut=" << vertex_cut
              << " truth_vertex_reweight_on=" << truth_vertex_reweight_on
              << std::endl;
  }

  // ------------------------------------------------------------------
  // Sample kinematics and cross-section weight
  // ------------------------------------------------------------------
  PPG12::SampleConfig sc;
  std::string input_filetype = MapInputFiletype(filetype);

  if (issim) {
    sc = PPG12::GetSampleConfig(filetype);
    if (!sc.valid) {
      std::cerr << "[ClusterSizeStudy] Unknown filetype: " << filetype << std::endl;
      return;
    }
  }

  const float weight_base  = issim ? sc.weight : 1.0f;
  const float cross_weight = weight_base * mix_weight;

  TChain chain("slimtree");
  const std::string infile = SampleFilePath(input_filetype, issim);
  int nadded = chain.Add(infile.c_str());
  Long64_t nentries = chain.GetEntries();
  std::cout << "[ClusterSizeStudy] filetype=" << filetype
            << " (input=" << input_filetype << ")"
            << " files=" << nadded
            << " entries=" << nentries
            << " weight=" << weight_base
            << " mix_weight=" << mix_weight
            << std::endl;

  if (nentries == 0) {
    std::cerr << "[ClusterSizeStudy] No entries loaded for " << filetype << std::endl;
    return;
  }

  // Detect whether the secondary truth-vertex branch is present. Falls back
  // to w(z_h) only when absent (matches ShowerShapeCheck.C's optional-pointer
  // guard behavior on pre-reprocess bdt_split.root files).
  need_vertexz_truth_mb =
      (have_cfg && issim && truth_vertex_reweight_on &&
       chain.GetBranch("vertexz_truth_mb") != nullptr);

  // ------------------------------------------------------------------
  // Output path
  // ------------------------------------------------------------------
  std::string outpath = outpath_override;
  if (outpath.empty()) {
    outpath = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/"
              "results/cluster_size_" + filetype + ".root";
  }
  TFile outfile(outpath.c_str(), "RECREATE");
  if (outfile.IsZombie()) {
    std::cerr << "[ClusterSizeStudy] Cannot open " << outpath << std::endl;
    return;
  }

  const std::vector<std::pair<std::string, std::string>> nodes = {
      {"split",   "CLUSTERINFO_CEMC"},
      {"nosplit", "CLUSTERINFO_CEMC_NO_SPLIT"},
  };

  std::map<std::string, NodeHists> hists;
  for (const auto& kv : nodes) {
    hists[kv.first] = MakeNodeHists(kv.first, have_cfg);
  }

  TTreeReader reader(&chain);

  // ------------------------------------------------------------------
  // Event-level readers
  // ------------------------------------------------------------------
  std::unique_ptr<TTreeReaderValue<Int_t>>    r_runnumber;
  std::unique_ptr<TTreeReaderValue<Float_t>>  r_vertexz;
  std::unique_ptr<TTreeReaderValue<Float_t>>  r_vertexz_truth;
  std::unique_ptr<TTreeReaderValue<Float_t>>  r_vertexz_truth_mb;
  if (have_cfg) {
    r_runnumber.reset(new TTreeReaderValue<Int_t>(reader, "runnumber"));
    r_vertexz.reset(new TTreeReaderValue<Float_t>(reader, "vertexz"));
    if (issim && truth_vertex_reweight_on) {
      r_vertexz_truth.reset(new TTreeReaderValue<Float_t>(reader, "vertexz_truth"));
      if (need_vertexz_truth_mb) {
        r_vertexz_truth_mb.reset(new TTreeReaderValue<Float_t>(reader, "vertexz_truth_mb"));
      }
    }
  }

  // ------------------------------------------------------------------
  // Per-node cluster-branch readers
  // ------------------------------------------------------------------
  std::map<std::string, std::unique_ptr<TTreeReaderValue<Int_t>>>    r_ncluster;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cet;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_ceta;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Int_t>>>    r_ownership;
  // cut-related arrays (loaded only if have_cfg)
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cluster_prob;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cluster_weta_cogx;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cluster_wphi_cogx;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cluster_e11;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cluster_e33;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cluster_e32;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cluster_e35;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cluster_et1;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cluster_et2;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cluster_et3;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cluster_et4;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cluster_npb_score;
  std::map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>>  r_cluster_bdt;
  // Track whether the cut-related branches exist (drives per-node cut skip)
  std::map<std::string, bool> can_apply_cuts;

  for (const auto& kv : nodes) {
    const std::string& key  = kv.first;
    const std::string& node = kv.second;
    r_ncluster[key].reset(new TTreeReaderValue<Int_t>(
        reader, Form("ncluster_%s", node.c_str())));
    r_cet[key].reset(new TTreeReaderArray<Float_t>(
        reader, Form("cluster_Et_%s", node.c_str())));
    r_ceta[key].reset(new TTreeReaderArray<Float_t>(
        reader, Form("cluster_Eta_%s", node.c_str())));
    r_ownership[key].reset(new TTreeReaderArray<Int_t>(
        reader, Form("cluster_ownership_array_%s", node.c_str())));

    can_apply_cuts[key] = false;
    if (!have_cfg) continue;

    const std::string prob_br  = std::string("cluster_prob_")          + node;
    const std::string weta_br  = std::string("cluster_weta_cogx_")     + node;
    const std::string wphi_br  = std::string("cluster_wphi_cogx_")     + node;
    const std::string e11_br   = std::string("cluster_e11_")           + node;
    const std::string e33_br   = std::string("cluster_e33_")           + node;
    const std::string e32_br   = std::string("cluster_e32_")           + node;
    const std::string e35_br   = std::string("cluster_e35_")           + node;
    const std::string et1_br   = std::string("cluster_et1_")           + node;
    const std::string et2_br   = std::string("cluster_et2_")           + node;
    const std::string et3_br   = std::string("cluster_et3_")           + node;
    const std::string et4_br   = std::string("cluster_et4_")           + node;
    const std::string npb_br   = std::string("cluster_npb_score_")     + node;
    const std::string bdt_br   = std::string("cluster_bdt_")           + node +
                                 "_" + cuts.bdt_model_for_tight;

    const bool has_all_cut_br =
        chain.GetBranch(prob_br.c_str()) &&
        chain.GetBranch(weta_br.c_str()) &&
        chain.GetBranch(wphi_br.c_str()) &&
        chain.GetBranch(e11_br.c_str())  &&
        chain.GetBranch(e33_br.c_str())  &&
        chain.GetBranch(e32_br.c_str())  &&
        chain.GetBranch(e35_br.c_str())  &&
        chain.GetBranch(et1_br.c_str())  &&
        chain.GetBranch(et2_br.c_str())  &&
        chain.GetBranch(et3_br.c_str())  &&
        chain.GetBranch(et4_br.c_str())  &&
        chain.GetBranch(npb_br.c_str())  &&
        chain.GetBranch(bdt_br.c_str());

    can_apply_cuts[key] = has_all_cut_br;
    if (!has_all_cut_br) {
      std::cout << "[ClusterSizeStudy] " << key << " missing cut branches; "
                << "cut TProfiles will not be filled for this node." << std::endl;
      continue;
    }

    r_cluster_prob[key].reset(new TTreeReaderArray<Float_t>(reader, prob_br.c_str()));
    r_cluster_weta_cogx[key].reset(new TTreeReaderArray<Float_t>(reader, weta_br.c_str()));
    r_cluster_wphi_cogx[key].reset(new TTreeReaderArray<Float_t>(reader, wphi_br.c_str()));
    r_cluster_e11[key].reset(new TTreeReaderArray<Float_t>(reader, e11_br.c_str()));
    r_cluster_e33[key].reset(new TTreeReaderArray<Float_t>(reader, e33_br.c_str()));
    r_cluster_e32[key].reset(new TTreeReaderArray<Float_t>(reader, e32_br.c_str()));
    r_cluster_e35[key].reset(new TTreeReaderArray<Float_t>(reader, e35_br.c_str()));
    r_cluster_et1[key].reset(new TTreeReaderArray<Float_t>(reader, et1_br.c_str()));
    r_cluster_et2[key].reset(new TTreeReaderArray<Float_t>(reader, et2_br.c_str()));
    r_cluster_et3[key].reset(new TTreeReaderArray<Float_t>(reader, et3_br.c_str()));
    r_cluster_et4[key].reset(new TTreeReaderArray<Float_t>(reader, et4_br.c_str()));
    r_cluster_npb_score[key].reset(new TTreeReaderArray<Float_t>(reader, npb_br.c_str()));
    r_cluster_bdt[key].reset(new TTreeReaderArray<Float_t>(reader, bdt_br.c_str()));
  }

  // ------------------------------------------------------------------
  // Truth-matching readers (sim only)
  // ------------------------------------------------------------------
  std::unique_ptr<TTreeReaderValue<Int_t>>   r_nparticles;
  std::unique_ptr<TTreeReaderArray<Int_t>>   r_particle_pid;
  std::unique_ptr<TTreeReaderArray<Float_t>> r_particle_Pt;
  std::unique_ptr<TTreeReaderArray<Float_t>> r_particle_Eta;
  std::unique_ptr<TTreeReaderValue<Int_t>>   r_njet_truth;
  std::unique_ptr<TTreeReaderArray<Float_t>> r_jet_truth_Pt;
  std::unique_ptr<TTreeReaderArray<Float_t>> r_jet_truth_Eta;

  if (issim) {
    r_nparticles.reset(new TTreeReaderValue<Int_t>(reader, "nparticles"));
    r_particle_pid.reset(new TTreeReaderArray<Int_t>(reader, "particle_pid"));
    r_particle_Pt.reset(new TTreeReaderArray<Float_t>(reader, "particle_Pt"));
    r_particle_Eta.reset(new TTreeReaderArray<Float_t>(reader, "particle_Eta"));
    const bool has_r04 =
        chain.GetBranch("njet_truth_AntiKt_Truth_r04") != nullptr;
    const char* njet_name = has_r04 ? "njet_truth_AntiKt_Truth_r04"
                                    : "njet_truth";
    const char* jetpt_name = has_r04 ? "jet_truth_Pt_AntiKt_Truth_r04"
                                     : "jet_truth_Pt";
    const char* jeteta_name = has_r04 ? "jet_truth_Eta_AntiKt_Truth_r04"
                                      : "jet_truth_Eta";
    std::cout << "[ClusterSizeStudy] jet branch: " << njet_name << std::endl;
    r_njet_truth.reset(new TTreeReaderValue<Int_t>(reader, njet_name));
    r_jet_truth_Pt.reset(new TTreeReaderArray<Float_t>(reader, jetpt_name));
    r_jet_truth_Eta.reset(new TTreeReaderArray<Float_t>(reader, jeteta_name));
  }

  // ------------------------------------------------------------------
  // Event loop
  // ------------------------------------------------------------------
  Long64_t ievent  = 0;
  Long64_t nfilled[3] = {0, 0, 0};
  Long64_t n_pass_cfg = 0;
  while (reader.Next()) {
    ++ievent;
    if (ievent % 200000 == 0) {
      std::cout << "  event " << ievent << "/" << nentries
                << " filled(all/common/tight)=" << nfilled[0]
                << "/" << nfilled[1]
                << "/" << nfilled[2] << std::endl;
    }

    // run range filter (data only)
    if (have_cfg && !issim) {
      const int run = **r_runnumber;
      if (run_min >= 0 && run < run_min) continue;
      if (run_max >= 0 && run > run_max) continue;
      if (!allowed_runs.empty() && allowed_runs.find(run) == allowed_runs.end()) continue;
    }

    // vertex cut
    float vtxz = have_cfg ? (**r_vertexz) : 0.f;
    if (have_cfg && std::fabs(vtxz) > vertex_cut) continue;

    // event-level weight
    double w_event = cross_weight;

    // Truth-vertex reweight for MC
    if (have_cfg && issim && truth_vertex_reweight_on) {
      const float tw = (need_vertexz_truth_mb && r_vertexz_truth_mb)
        ? TruthVertexWeight(h_truth_vtx_reweight, **r_vertexz_truth, **r_vertexz_truth_mb)
        : TruthVertexWeight(h_truth_vtx_reweight, **r_vertexz_truth);
      w_event *= tw;
      if (tw <= 0.0) continue;
    }

    // MC sample kinematic window (photon_pt or jet_pt)
    if (issim) {
      float max_photon_pt = 0.0f;
      const int npart = **r_nparticles;
      for (int ip = 0; ip < npart; ++ip) {
        if ((*r_particle_pid)[ip] == 22 &&
            std::abs((*r_particle_Eta)[ip]) < kEtaFidMax &&
            (*r_particle_Pt)[ip] > max_photon_pt) {
          max_photon_pt = (*r_particle_Pt)[ip];
        }
      }
      float max_jet_pt = 0.0f;
      const int njet = **r_njet_truth;
      for (int ij = 0; ij < njet; ++ij) {
        if (std::abs((*r_jet_truth_Eta)[ij]) >= kEtaFidMax) continue;
        if ((*r_jet_truth_Pt)[ij] > max_jet_pt) {
          max_jet_pt = (*r_jet_truth_Pt)[ij];
        }
      }

      if (sc.isbackground) {
        if (max_jet_pt < sc.jet_pt_lower || max_jet_pt >= sc.jet_pt_upper) continue;
      } else {
        if (max_photon_pt < sc.photon_pt_lower ||
            max_photon_pt >= sc.photon_pt_upper) continue;
      }
    }

    ++n_pass_cfg;

    for (const auto& kv : nodes) {
      const std::string& key = kv.first;
      const int ncl = **r_ncluster[key];

      for (int icl = 0; icl < ncl; ++icl) {
        const float et  = (*r_cet[key])[icl];
        const float eta = (*r_ceta[key])[icl];
        if (et <= 0) continue;
        if (std::abs(eta) >= kEtaFidMax) continue;
        if (issim && et > sc.cluster_ET_upper) continue;

        int n_owned     = 0;
        int min_eta_idx = 999, max_eta_idx = -1;
        int min_phi_idx = 999, max_phi_idx = -1;
        const int base = icl * kNtow;
        for (int i = 0; i < kNtow; ++i) {
          if ((*r_ownership[key])[base + i] != 1) continue;
          ++n_owned;
          const int ei = i / kGridN;
          const int pi = i % kGridN;
          if (ei < min_eta_idx) min_eta_idx = ei;
          if (ei > max_eta_idx) max_eta_idx = ei;
          if (pi < min_phi_idx) min_phi_idx = pi;
          if (pi > max_phi_idx) max_phi_idx = pi;
        }
        if (n_owned == 0) continue;

        const int width_eta = max_eta_idx - min_eta_idx + 1;
        const int width_phi = max_phi_idx - min_phi_idx + 1;

        // Baseline fill (all clusters passing the fiducial)
        auto& H = hists[key];
        const double w = w_event;
        H.prof_n_owned[0]  ->Fill(et, n_owned,   w);
        H.prof_width_eta[0]->Fill(et, width_eta, w);
        H.prof_width_phi[0]->Fill(et, width_phi, w);
        H.h2_n_owned[0]    ->Fill(et, n_owned,   w);
        H.h2_width_eta[0]  ->Fill(et, width_eta, w);
        H.h2_width_phi[0]  ->Fill(et, width_phi, w);
        ++nfilled[0];

        // Cut-variant fills (only if config available and this node has all required branches)
        if (!have_cfg || !can_apply_cuts[key]) continue;

        const float prob_v = (*r_cluster_prob[key])[icl];
        const float weta_v = (*r_cluster_weta_cogx[key])[icl];
        const float wphi_v = (*r_cluster_wphi_cogx[key])[icl];
        const float e11_v  = (*r_cluster_e11[key])[icl];
        const float e33_v  = (*r_cluster_e33[key])[icl];
        const float e32_v  = (*r_cluster_e32[key])[icl];
        const float e35_v  = (*r_cluster_e35[key])[icl];
        const float et1_v  = (*r_cluster_et1[key])[icl];
        const float et2_v  = (*r_cluster_et2[key])[icl];
        const float et3_v  = (*r_cluster_et3[key])[icl];
        const float et4_v  = (*r_cluster_et4[key])[icl];
        const float npb_v  = (*r_cluster_npb_score[key])[icl];
        const float bdt_v  = (*r_cluster_bdt[key])[icl];

        const float e11_o_e33 = (e33_v > 0) ? (e11_v / e33_v) : 0.f;
        const float e32_o_e35 = (e35_v > 0) ? (e32_v / e35_v) : 0.f;

        // common cut (mirrors ShowerShapeCheck.C common block)
        const bool common_pass =
            (prob_v > cuts.common_prob_min) && (prob_v < cuts.common_prob_max) &&
            (e11_o_e33 > cuts.common_e11_over_e33_min) &&
            (e11_o_e33 < cuts.common_e11_over_e33_max) &&
            ((!cuts.npb_cut_on) || (npb_v > cuts.npb_score_cut)) &&
            (weta_v < cuts.common_weta_cogx_bound);

        if (!common_pass) continue;
        H.prof_n_owned[1]  ->Fill(et, n_owned,   w);
        H.prof_width_eta[1]->Fill(et, width_eta, w);
        H.prof_width_phi[1]->Fill(et, width_phi, w);
        H.h2_n_owned[1]    ->Fill(et, n_owned,   w);
        H.h2_width_eta[1]  ->Fill(et, width_eta, w);
        H.h2_width_phi[1]  ->Fill(et, width_phi, w);
        ++nfilled[1];

        // tight cut (ET-dependent parametric thresholds)
        const float tight_weta_cogx_max = cuts.tight_weta_cogx_max_b + cuts.tight_weta_cogx_max_s * et;
        const float tight_wphi_cogx_max = cuts.tight_wphi_cogx_max_b + cuts.tight_wphi_cogx_max_s * et;
        const float tight_et1_min       = cuts.tight_et1_min_b       + cuts.tight_et1_min_s       * et;
        const float tight_bdt_min_et    = cuts.tight_bdt_min_intercept + cuts.tight_bdt_min_slope * et;

        const bool tight_pass =
            (weta_v > cuts.tight_weta_cogx_min) && (weta_v < tight_weta_cogx_max) &&
            (wphi_v > cuts.tight_wphi_cogx_min) && (wphi_v < tight_wphi_cogx_max) &&
            (et1_v  > tight_et1_min)            && (et1_v  < cuts.tight_et1_max) &&
            (et2_v  > cuts.tight_et2_min)       && (et2_v  < cuts.tight_et2_max) &&
            (et3_v  > cuts.tight_et3_min)       && (et3_v  < cuts.tight_et3_max) &&
            (et4_v  > cuts.tight_et4_min)       && (et4_v  < cuts.tight_et4_max) &&
            (e11_o_e33 > cuts.tight_e11_over_e33_min) && (e11_o_e33 < cuts.tight_e11_over_e33_max) &&
            (e32_o_e35 > cuts.tight_e32_over_e35_min) && (e32_o_e35 < cuts.tight_e32_over_e35_max) &&
            (prob_v > cuts.tight_prob_min) && (prob_v < cuts.tight_prob_max) &&
            (bdt_v  > tight_bdt_min_et)    && (bdt_v  < cuts.tight_bdt_max);

        if (tight_pass) {
          H.prof_n_owned[2]  ->Fill(et, n_owned,   w);
          H.prof_width_eta[2]->Fill(et, width_eta, w);
          H.prof_width_phi[2]->Fill(et, width_phi, w);
          H.h2_n_owned[2]    ->Fill(et, n_owned,   w);
          H.h2_width_eta[2]  ->Fill(et, width_eta, w);
          H.h2_width_phi[2]  ->Fill(et, width_phi, w);
          ++nfilled[2];
        }
      }
    }
  }

  // ------------------------------------------------------------------
  // Persist
  // ------------------------------------------------------------------
  outfile.cd();

  TNamed meta_filetype("filetype", filetype.c_str());
  meta_filetype.Write();
  TNamed meta_input("input_filetype", input_filetype.c_str());
  meta_input.Write();
  TParameter<float> meta_weight("weight_rel", weight_base);
  meta_weight.Write();
  TParameter<float> meta_mix("mix_weight", mix_weight);
  meta_mix.Write();
  const float ref_xsec = issim
      ? (sc.isbackground ? PPG12::jet50cross : PPG12::photon20cross)
      : 1.0f;
  TParameter<float> meta_ref_xsec("ref_xsec_pb", ref_xsec);
  meta_ref_xsec.Write();
  TParameter<int> meta_isbg("isbackground", (issim && sc.isbackground) ? 1 : 0);
  meta_isbg.Write();
  TParameter<Long64_t> meta_nevents("nevents_processed", ievent);
  meta_nevents.Write();
  TParameter<Long64_t> meta_nevents_pass("nevents_pass_cfg", n_pass_cfg);
  meta_nevents_pass.Write();
  TParameter<Long64_t> meta_nfilled_all("ncluster_filled_all",    nfilled[0]);
  meta_nfilled_all.Write();
  TParameter<Long64_t> meta_nfilled_common("ncluster_filled_common", nfilled[1]);
  meta_nfilled_common.Write();
  TParameter<Long64_t> meta_nfilled_tight("ncluster_filled_tight",  nfilled[2]);
  meta_nfilled_tight.Write();
  TNamed meta_config("config", configname.c_str());
  meta_config.Write();
  TParameter<int> meta_run_min("run_min", run_min);
  meta_run_min.Write();
  TParameter<int> meta_run_max("run_max", run_max);
  meta_run_max.Write();
  TParameter<float> meta_vertex_cut("vertex_cut", vertex_cut);
  meta_vertex_cut.Write();
  TParameter<int> meta_tv_on("truth_vertex_reweight_on", truth_vertex_reweight_on);
  meta_tv_on.Write();
  TParameter<int> meta_mb_used("vertexz_truth_mb_used", need_vertexz_truth_mb ? 1 : 0);
  meta_mb_used.Write();
  TNamed meta_tv_file("truth_vertex_reweight_file", truth_vertex_reweight_file.c_str());
  meta_tv_file.Write();

  for (const auto& kv : hists) WriteNodeHists(kv.second);

  outfile.Close();
  std::cout << "[ClusterSizeStudy] Wrote " << outpath
            << " events=" << ievent
            << " pass_cfg=" << n_pass_cfg
            << " fills(all/common/tight)=" << nfilled[0]
            << "/" << nfilled[1]
            << "/" << nfilled[2]
            << std::endl;
}
