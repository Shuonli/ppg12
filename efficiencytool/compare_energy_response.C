// compare_energy_response.C
// Photon energy response R = ET_reco / pT_truth (and R_E = E_reco/E_truth)
// across single, double, and mixed (0 mrad / 1.5 mrad, with/without vertex
// reweighting) MC, using the same selection as compare_efficiency_deltaR.C.
//
// Matching: trkID-only (noDR). The delta-R cut is dropped, consistent with
// the companion efficiency study.
//
// Usage:
//   root -l -b -q 'compare_energy_response.C("config_bdt_nom.yaml", -1)'
//   root -l -b -q 'compare_energy_response.C("config_bdt_nom.yaml", 1000000)'
//
// Output: /sphenix/user/shuhangli/ppg12/efficiencytool/results/energy_response_comparison.root

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TSystem.h>
#include <TMath.h>
#include <TROOT.h>

#include <yaml-cpp/yaml.h>

// ------------- Double-sided Crystal-Ball PDF --------------------------
// par[0] = N
// par[1] = mu
// par[2] = sigma
// par[3] = alpha_L (>0)    par[4] = n_L (>1)
// par[5] = alpha_H (>0)    par[6] = n_H (>1)
// Gaussian core on  -alpha_L < (x-mu)/sigma < alpha_H
// Power-law tail on both sides outside that window.
static Double_t double_cb(Double_t *x, Double_t *par)
{
    Double_t N  = par[0];
    Double_t mu = par[1];
    Double_t sg = par[2];
    Double_t aL = par[3], nL = par[4];
    Double_t aH = par[5], nH = par[6];
    if (sg <= 0) return 0;
    Double_t t = (x[0] - mu) / sg;
    if (t > -aL && t < aH)
        return N * TMath::Exp(-0.5 * t * t);
    if (t <= -aL) {
        Double_t A = TMath::Power(nL / TMath::Abs(aL), nL) * TMath::Exp(-0.5 * aL * aL);
        Double_t B = nL / TMath::Abs(aL) - TMath::Abs(aL);
        return N * A * TMath::Power(B - t, -nL);
    }
    // t >= aH
    Double_t A = TMath::Power(nH / TMath::Abs(aH), nH) * TMath::Exp(-0.5 * aH * aH);
    Double_t B = nH / TMath::Abs(aH) - TMath::Abs(aH);
    return N * A * TMath::Power(B + t, -nH);
}

// ------------- Cluster-weighted pileup fractions (from calc_pileup_range.C)
// 0 mrad crossing angle : 22.4% double, 77.6% single
// 1.5 mrad crossing angle: 7.9% double, 92.1% single
static const float FRAC_DOUBLE_0MRAD   = 0.224f;
static const float FRAC_DOUBLE_1P5MRAD = 0.079f;

// ------------- Input MC sample paths
static const char *PATH_PHOTON10_NOM    = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root";
static const char *PATH_PHOTON10_DOUBLE = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_double/bdt_split.root";

// ------------- Data vertex-distribution references (Pass-1 products from the
// main pipeline; see ShowerShapeCheck.C vertex-reweight block).
static const char *PATH_DATA_VTX_0MRAD   = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_0rad_vtxscan.root";
static const char *PATH_DATA_VTX_1P5MRAD = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_1p5rad_vtxscan.root";

// =====================================================================
// Sample set: a linear combination of the two MC inputs with optional
// vertex reweighting.  During the event loop, a cluster is filled into
// each set it belongs to, with weight = alpha(sample) * vtx_reweight(set).
// =====================================================================
struct SetSpec {
    TString     name;
    float       alpha_nom;      // weight on photon10 nominal clusters
    float       alpha_double;   // weight on photon10_double clusters
    bool        do_vtx_reweight;
    TString     data_vtx_path;  // used iff do_vtx_reweight
    TH1        *h_vtx_weight;   // data/sim ratio; filled after Pass 1
};

// =====================================================================
// Helper: build a data/sim vertex-reweight histogram following the exact
// recipe in ShowerShapeCheck.C (rebin by 5, normalize, divide, smooth 1x).
// =====================================================================
static TH1 *BuildVertexReweight(TH1 *h_sim_combined, const char *data_path,
                                 int rebin_factor, int smooth_passes,
                                 const char *label)
{
    TFile *fvtx_data = TFile::Open(data_path, "READ");
    if (!fvtx_data || fvtx_data->IsZombie()) {
        std::cerr << "[VertexReweight] FATAL: cannot open " << data_path << std::endl;
        return nullptr;
    }
    TH1 *hdata_vtx = dynamic_cast<TH1*>(fvtx_data->Get("h_vertexz"));
    if (!hdata_vtx) {
        std::cerr << "[VertexReweight] FATAL: h_vertexz missing in " << data_path << std::endl;
        fvtx_data->Close();
        return nullptr;
    }

    TH1 *hdata_clone = dynamic_cast<TH1*>(hdata_vtx->Clone(Form("h_vtx_data_%s", label)));
    TH1 *hsim_clone  = dynamic_cast<TH1*>(h_sim_combined->Clone(Form("h_vtx_sim_%s", label)));
    hdata_clone->SetDirectory(nullptr);
    hsim_clone->SetDirectory(nullptr);

    if (rebin_factor > 1) {
        hdata_clone->Rebin(rebin_factor);
        hsim_clone->Rebin(rebin_factor);
    }

    if (hdata_clone->Integral() > 0) hdata_clone->Scale(1.0 / hdata_clone->Integral());
    if (hsim_clone->Integral() > 0)  hsim_clone->Scale(1.0 / hsim_clone->Integral());
    hdata_clone->Divide(hsim_clone);

    if (smooth_passes > 0) hdata_clone->Smooth(smooth_passes);

    fvtx_data->Close();
    delete fvtx_data;
    delete hsim_clone;

    std::cout << "[VertexReweight] " << label << ": ratio built from " << data_path
              << "  nbins=" << hdata_clone->GetNbinsX()
              << "  bin width=" << hdata_clone->GetBinWidth(1) << " cm"
              << "  mean weight=" << hdata_clone->Integral() / hdata_clone->GetNbinsX() << std::endl;
    return hdata_clone;
}

// =====================================================================
// Main function
// =====================================================================
void compare_energy_response(TString configname = "config_bdt_nom.yaml",
                             int maxEvents = -1)
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");

    // -----------------------------------------------------------------
    // Load configuration (same as compare_efficiency_deltaR.C)
    // -----------------------------------------------------------------
    YAML::Node configYaml = YAML::LoadFile(configname.Data());

    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();
    std::string bdt_model_name  = configYaml["input"]["bdt_model_name"].as<std::string>("base");

    std::vector<float>       bdt_et_bin_edges;
    std::vector<std::string> bdt_et_bin_models;
    bool use_et_binned_bdt = false;
    if (configYaml["input"]["bdt_et_bin_edges"] && configYaml["input"]["bdt_et_bin_models"]) {
        for (auto v : configYaml["input"]["bdt_et_bin_edges"])
            bdt_et_bin_edges.push_back(v.as<float>());
        for (auto v : configYaml["input"]["bdt_et_bin_models"])
            bdt_et_bin_models.push_back(v.as<std::string>());
        use_et_binned_bdt = (bdt_et_bin_models.size() == bdt_et_bin_edges.size() - 1);
    }

    // Analysis cuts (subset needed here)
    int use_topo_iso       = configYaml["analysis"]["use_topo_iso"].as<int>(0);
    int conesize           = configYaml["analysis"]["cone_size"].as<int>();
    float truthisocut      = configYaml["analysis"]["truth_iso_max"].as<float>();
    float recoiso_min      = configYaml["analysis"]["reco_iso_min"].as<float>();
    float recoiso_max_b    = configYaml["analysis"]["reco_iso_max_b"].as<float>();
    float recoiso_max_s    = configYaml["analysis"]["reco_iso_max_s"].as<float>();
    float vertexcut        = configYaml["analysis"]["vertex_cut"].as<float>();
    float reco_min_ET      = configYaml["analysis"]["reco_min_ET"].as<float>();
    float mc_iso_shift     = configYaml["analysis"]["mc_iso_shift"].as<float>(0.0);
    float mc_iso_scale     = configYaml["analysis"]["mc_iso_scale"].as<float>(1.0);

    std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();

    std::vector<float> pT_bins_truth;
    if (configYaml["analysis"]["pT_bins_truth"]) {
        pT_bins_truth = configYaml["analysis"]["pT_bins_truth"].as<std::vector<float>>();
    } else {
        pT_bins_truth = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45};
    }
    int NptBins = (int)pT_bins_truth.size() - 1;
    double ptEdges[100];
    for (int i = 0; i <= NptBins; i++) ptEdges[i] = pT_bins_truth[i];
    float pTmin_truth = pT_bins_truth.front();
    float pTmax_truth = pT_bins_truth.back();

    // Common cuts
    float common_prob_min             = configYaml["analysis"]["common"]["prob_min"].as<float>();
    float common_prob_max             = configYaml["analysis"]["common"]["prob_max"].as<float>();
    float common_e11_over_e33_min     = configYaml["analysis"]["common"]["e11_over_e33_min"].as<float>();
    float common_e11_over_e33_max     = configYaml["analysis"]["common"]["e11_over_e33_max"].as<float>();
    float common_wr_cogx_bound        = configYaml["analysis"]["common"]["wr_cogx_bound"].as<float>();
    float common_cluster_weta_cogx_bound = configYaml["analysis"]["common"]["cluster_weta_cogx_bound"].as<float>();
    int   common_npb_cut_on           = configYaml["analysis"]["common"]["npb_cut_on"].as<int>(0);
    float common_npb_score_cut        = configYaml["analysis"]["common"]["npb_score_cut"].as<float>(0.5);

    // Tight cuts (BDT + shower-shape) — same parametric forms as nominal
    float tight_bdt_max           = configYaml["analysis"]["tight"]["bdt_max"].as<float>(1);
    float tight_bdt_min           = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0);
    float tight_bdt_min_slope     = configYaml["analysis"]["tight"]["bdt_min_slope"].as<float>(0);
    float tight_bdt_min_intercept = configYaml["analysis"]["tight"]["bdt_min_intercept"].as<float>(tight_bdt_min);
    float tight_weta_cogx_min     = configYaml["analysis"]["tight"]["weta_cogx_min"].as<float>();
    float tight_weta_cogx_max_b   = configYaml["analysis"]["tight"]["weta_cogx_max_b"].as<float>();
    float tight_weta_cogx_max_s   = configYaml["analysis"]["tight"]["weta_cogx_max_s"].as<float>();
    float tight_wphi_cogx_min     = configYaml["analysis"]["tight"]["wphi_cogx_min"].as<float>();
    float tight_wphi_cogx_max_b   = configYaml["analysis"]["tight"]["wphi_cogx_max_b"].as<float>();
    float tight_wphi_cogx_max_s   = configYaml["analysis"]["tight"]["wphi_cogx_max_s"].as<float>();
    float tight_e11_over_e33_min  = configYaml["analysis"]["tight"]["e11_over_e33_min"].as<float>();
    float tight_e11_over_e33_max  = configYaml["analysis"]["tight"]["e11_over_e33_max"].as<float>();
    float tight_e32_over_e35_min  = configYaml["analysis"]["tight"]["e32_over_e35_min"].as<float>();
    float tight_e32_over_e35_max  = configYaml["analysis"]["tight"]["e32_over_e35_max"].as<float>();
    float tight_et1_max           = configYaml["analysis"]["tight"]["et1_max"].as<float>();
    float tight_et1_min_b         = configYaml["analysis"]["tight"]["et1_min_b"].as<float>();
    float tight_et1_min_s         = configYaml["analysis"]["tight"]["et1_min_s"].as<float>();
    float tight_et2_min           = configYaml["analysis"]["tight"]["et2_min"].as<float>(0.0);
    float tight_et2_max           = configYaml["analysis"]["tight"]["et2_max"].as<float>(1.0);
    float tight_et3_min           = configYaml["analysis"]["tight"]["et3_min"].as<float>(0.0);
    float tight_et3_max           = configYaml["analysis"]["tight"]["et3_max"].as<float>(1.0);
    float tight_et4_min           = configYaml["analysis"]["tight"]["et4_min"].as<float>();
    float tight_et4_max           = configYaml["analysis"]["tight"]["et4_max"].as<float>();
    float tight_prob_min          = configYaml["analysis"]["tight"]["prob_min"].as<float>();
    float tight_prob_max          = configYaml["analysis"]["tight"]["prob_max"].as<float>();

    // Vertex-reweight knobs (same defaults as ShowerShapeCheck.C)
    int vtx_rebin = configYaml["analysis"]["vertex_reweight_rebin"].as<int>(5);
    int vtx_smooth = configYaml["analysis"]["vertex_reweight_smooth"].as<int>(1);

    std::cout << "[compare_energy_response] config: " << configname << std::endl;
    std::cout << "  vertexcut = " << vertexcut << " cm" << std::endl;
    std::cout << "  pT_bins_truth: "; for (auto v : pT_bins_truth) std::cout << v << " "; std::cout << std::endl;
    std::cout << "  0 mrad DI fraction: " << FRAC_DOUBLE_0MRAD << std::endl;
    std::cout << "  1.5 mrad DI fraction: " << FRAC_DOUBLE_1P5MRAD << std::endl;

    // -----------------------------------------------------------------
    // Define sample sets (each = linear combo of nom + double, optional vtx rw)
    // -----------------------------------------------------------------
    std::vector<SetSpec> sets;
    sets.push_back({"single",             1.0f, 0.0f, false, "", nullptr});
    sets.push_back({"double",             0.0f, 1.0f, false, "", nullptr});
    sets.push_back({"mixed_0mrad",        1.0f - FRAC_DOUBLE_0MRAD,   FRAC_DOUBLE_0MRAD,   false, "",                       nullptr});
    sets.push_back({"mixed_0mrad_vtxrw",  1.0f - FRAC_DOUBLE_0MRAD,   FRAC_DOUBLE_0MRAD,   true,  PATH_DATA_VTX_0MRAD,     nullptr});
    sets.push_back({"mixed_1p5mrad",      1.0f - FRAC_DOUBLE_1P5MRAD, FRAC_DOUBLE_1P5MRAD, false, "",                       nullptr});
    sets.push_back({"mixed_1p5mrad_vtxrw",1.0f - FRAC_DOUBLE_1P5MRAD, FRAC_DOUBLE_1P5MRAD, true,  PATH_DATA_VTX_1P5MRAD,   nullptr});

    // -----------------------------------------------------------------
    // Response binning
    //   axis 0: truth pT, using pT_bins_truth
    //   axis 1: response r = ET_reco / pT_truth (or E_reco/E_truth), 0..2.5 in 0.01 steps
    // -----------------------------------------------------------------
    const int   NrespBins = 250;
    const float respMin   = 0.0f;
    const float respMax   = 2.5f;

    // Selection levels (renamed 2026-04-16):
    //   reco       : trkID-matched only (no ΔR, no common)   [was "matched"]
    //   common_cut : trkID + common cuts                     [was "reco"]
    //   tight_iso  : trkID + common + iso + tight BDT        [was "selected"]
    const std::vector<TString> levels = {"reco", "common_cut", "tight_iso"};

    // Histogram maps
    std::map<TString, TH2F*> hRespET;   // ET_reco / pT_truth
    std::map<TString, TH2F*> hRespE;    // E_reco   / E_truth

    auto make2D = [&](const TString &name, const TString &ytitle) -> TH2F*
    {
        TH2F *h = new TH2F(name.Data(),
            Form(";truth p_{T} [GeV];%s", ytitle.Data()),
            NptBins, ptEdges,
            NrespBins, respMin, respMax);
        h->Sumw2();
        return h;
    };

    for (auto &s : sets) {
        for (auto &lvl : levels) {
            TString k1 = Form("h_respET_%s_%s", s.name.Data(), lvl.Data());
            TString k2 = Form("h_respE_%s_%s",  s.name.Data(), lvl.Data());
            hRespET[k1] = make2D(k1, "E_{T}^{reco} / p_{T}^{truth}");
            hRespE [k2] = make2D(k2, "E^{reco} / E^{truth}");
        }
    }

    // 1D vertex histograms for the two MC inputs (unweighted — we will
    // scale by (alpha_nom, alpha_double) in memory to build the per-set
    // combined sim distributions).
    TH1D *h_vtxz_single  = new TH1D("h_vtxz_single",  ";z_{vtx} [cm];entries", 200, -100, 100);
    TH1D *h_vtxz_double  = new TH1D("h_vtxz_double",  ";z_{vtx} [cm];entries", 200, -100, 100);
    h_vtxz_single->Sumw2();
    h_vtxz_double->Sumw2();

    // -----------------------------------------------------------------
    // Pass 1: vertex scan — loop over both input files, fill h_vtxz_*.
    // This is intentionally lightweight (only vertex+event gate, no
    // cluster or truth-particle processing).
    // -----------------------------------------------------------------
    auto fillVtxScan = [&](const char *path, TH1D *hout) {
        std::cout << "[Pass1] vertex scan: " << path << std::endl;
        TFile *fin = TFile::Open(path);
        if (!fin || fin->IsZombie()) {
            std::cerr << "ERROR: Cannot open " << path << std::endl;
            return;
        }
        TTreeReader reader("slimtree", fin);
        TTreeReaderValue<float> vertexz(reader, "vertexz");
        long long n = 0;
        while (reader.Next()) {
            if (maxEvents > 0 && n >= maxEvents) break;
            n++;
            if (fabs(*vertexz) > vertexcut) continue;
            hout->Fill(*vertexz);
        }
        std::cout << "  filled " << hout->GetEntries() << " entries (from " << n << " events)" << std::endl;
        fin->Close();
        delete fin;
    };
    fillVtxScan(PATH_PHOTON10_NOM,    h_vtxz_single);
    fillVtxScan(PATH_PHOTON10_DOUBLE, h_vtxz_double);

    // -----------------------------------------------------------------
    // Build combined sim vertex distributions per vtxrw set, then
    // construct data/sim ratio for reweighting.
    // -----------------------------------------------------------------
    for (auto &s : sets) {
        if (!s.do_vtx_reweight) continue;
        TH1D *h_combined = (TH1D*)h_vtxz_single->Clone(Form("h_vtxz_simcomb_%s", s.name.Data()));
        h_combined->SetDirectory(nullptr);
        h_combined->Scale(s.alpha_nom);
        h_combined->Add(h_vtxz_double, s.alpha_double);
        s.h_vtx_weight = BuildVertexReweight(h_combined, s.data_vtx_path.Data(),
                                             vtx_rebin, vtx_smooth, s.name.Data());
        // Preserve the combined sim for diagnostics
        h_combined->SetName(Form("h_vtxz_sim_combined_%s", s.name.Data()));
    }

    // -----------------------------------------------------------------
    // Pass 2: main event loop, processing each input sample once and
    // filling EVERY set that gives it non-zero weight.
    // -----------------------------------------------------------------
    struct InputSample {
        const char *path;
        int         which;   // 0 = nom, 1 = double
    };
    std::vector<InputSample> inputs = {
        {PATH_PHOTON10_NOM,    0},
        {PATH_PHOTON10_DOUBLE, 1},
    };

    for (auto &isamp : inputs)
    {
        std::cout << "\n[Pass2] processing " << isamp.path << " (which=" << isamp.which << ")" << std::endl;
        TFile *fin = TFile::Open(isamp.path);
        if (!fin || fin->IsZombie()) {
            std::cerr << "ERROR: Cannot open " << isamp.path << std::endl;
            continue;
        }

        TTreeReader reader("slimtree", fin);

        // ---- Event-level ----
        TTreeReaderValue<int>   nparticles(reader, "nparticles");
        TTreeReaderValue<int>   nclusters(reader, Form("ncluster_%s", clusternodename.c_str()));
        TTreeReaderValue<float> vertexz(reader, "vertexz");

        // ---- Truth particles ----
        TTreeReaderArray<float> particle_Pt(reader, "particle_Pt");
        TTreeReaderArray<float> particle_Eta(reader, "particle_Eta");
        TTreeReaderArray<float> particle_Phi(reader, "particle_Phi");
        TTreeReaderArray<int>   particle_pid(reader, "particle_pid");
        TTreeReaderArray<int>   particle_trkid(reader, "particle_trkid");
        TTreeReaderArray<int>   particle_photonclass(reader, "particle_photonclass");
        TTreeReaderArray<float> particle_truth_iso_03(reader, "particle_truth_iso_03");
        TTreeReaderArray<float> particle_truth_iso_04(reader, "particle_truth_iso_04");

        // ---- Clusters ----
        TTreeReaderArray<float> cluster_Et(reader,  Form("cluster_Et_%s",  clusternodename.c_str()));
        TTreeReaderArray<float> cluster_E(reader,   Form("cluster_E_%s",   clusternodename.c_str()));
        TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusternodename.c_str()));
        TTreeReaderArray<int>   cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_prob(reader, Form("cluster_prob_%s", clusternodename.c_str()));

        TTreeReaderArray<float> cluster_e11(reader, Form("cluster_e11_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_e33(reader, Form("cluster_e33_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_e32(reader, Form("cluster_e32_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_e35(reader, Form("cluster_e35_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_et1(reader, Form("cluster_et1_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_et2(reader, Form("cluster_et2_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_et3(reader, Form("cluster_et3_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_et4(reader, Form("cluster_et4_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_weta_cogx(reader, Form("cluster_weta_cogx_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_wphi_cogx(reader, Form("cluster_wphi_cogx_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_npb_score(reader, Form("cluster_npb_score_%s", clusternodename.c_str()));

        // Isolation
        TTreeReaderArray<float> cluster_iso_02(reader, Form("cluster_iso_02_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_iso_03(reader, Form("cluster_iso_03_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_iso_04(reader, Form("cluster_iso_04_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_iso_excl_04(reader, Form("cluster_iso_excl_04_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_iso_topo_03(reader, Form("cluster_iso_topo_03_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_iso_topo_04(reader, Form("cluster_iso_topo_04_%s", clusternodename.c_str()));

        // BDT score branches — one per unique model name
        std::vector<std::string> all_bdt_models = {bdt_model_name};
        if (use_et_binned_bdt)
            all_bdt_models.insert(all_bdt_models.end(), bdt_et_bin_models.begin(), bdt_et_bin_models.end());
        std::sort(all_bdt_models.begin(), all_bdt_models.end());
        all_bdt_models.erase(std::unique(all_bdt_models.begin(), all_bdt_models.end()), all_bdt_models.end());
        std::map<std::string, TTreeReaderArray<float>*> bdt_arrays;
        for (auto &mname : all_bdt_models) {
            bdt_arrays[mname] = new TTreeReaderArray<float>(reader,
                Form("cluster_bdt_%s_%s", clusternodename.c_str(), mname.c_str()));
        }

        long long nEvents = 0;
        while (reader.Next())
        {
            if (maxEvents > 0 && nEvents >= maxEvents) break;
            nEvents++;
            if (nEvents % 500000 == 0)
                std::cout << "  event " << nEvents << std::endl;

            float vtx = *vertexz;
            if (fabs(vtx) > vertexcut) continue;

            int ncl   = *nclusters;
            int npart = *nparticles;

            // Build trkID -> list of cluster indices
            std::map<int, std::vector<int>> trkid_to_clusters;
            for (int icl = 0; icl < ncl; icl++) {
                if (cluster_Et[icl] < reco_min_ET) continue;
                trkid_to_clusters[cluster_truthtrkID[icl]].push_back(icl);
            }

            for (int ip = 0; ip < npart; ip++)
            {
                if (particle_pid[ip] != 22) continue;
                if (particle_photonclass[ip] >= 3) continue;

                // Truth isolation (fiducial cross-section definition)
                float truthisoET = (conesize == 4) ? particle_truth_iso_04[ip]
                                                   : particle_truth_iso_03[ip];
                if (truthisoET >= truthisocut) continue;

                // Eta acceptance
                float truth_eta = particle_Eta[ip];
                bool in_eta = false;
                for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++) {
                    if (truth_eta > eta_bins[ieta] && truth_eta < eta_bins[ieta+1]) {
                        in_eta = true;
                        break;
                    }
                }
                if (!in_eta) continue;

                float truth_pT = particle_Pt[ip];
                if (truth_pT < pTmin_truth || truth_pT >= pTmax_truth) continue;

                int   truth_trkid = particle_trkid[ip];
                float truth_phi   = particle_Phi[ip];
                float truth_E     = truth_pT * cosh(truth_eta);

                auto it = trkid_to_clusters.find(truth_trkid);
                if (it == trkid_to_clusters.end()) continue;  // no trkID match

                // Pick the closest ΔR cluster (diagnostic — not a selection cut)
                float best_dR = 1e9;
                int   best_icl = -1;
                for (int icl : it->second) {
                    float deta = cluster_Eta[icl] - truth_eta;
                    float dphi = cluster_Phi[icl] - truth_phi;
                    if (dphi >  TMath::Pi()) dphi -= 2*TMath::Pi();
                    if (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
                    float dR = sqrt(deta*deta + dphi*dphi);
                    if (dR < best_dR) { best_dR = dR; best_icl = icl; }
                }
                if (best_icl < 0) continue;

                int icl = best_icl;
                float clusterET = cluster_Et[icl];
                float clusterE  = cluster_E[icl];
                if (truth_pT <= 0 || truth_E <= 0) continue;
                float rET = clusterET / truth_pT;
                float rE  = clusterE  / truth_E;

                // Common / iso / tight flags — computed once, shared across sets
                float e11_over_e33 = cluster_e11[icl] / cluster_e33[icl];
                float e32_over_e35 = cluster_e32[icl] / cluster_e35[icl];
                float wr_cogx      = cluster_wphi_cogx[icl] / cluster_weta_cogx[icl];

                bool passes_common =
                    cluster_prob[icl] > common_prob_min &&
                    cluster_prob[icl] < common_prob_max &&
                    e11_over_e33 > common_e11_over_e33_min &&
                    e11_over_e33 < common_e11_over_e33_max &&
                    wr_cogx > common_wr_cogx_bound &&
                    cluster_weta_cogx[icl] < common_cluster_weta_cogx_bound &&
                    (!common_npb_cut_on || cluster_npb_score[icl] > common_npb_score_cut);

                bool passes_iso = false;
                bool passes_tight = false;
                if (passes_common) {
                    float recoisoET = -999;
                    if      (use_topo_iso == 2)      recoisoET = cluster_iso_topo_04[icl];
                    else if (use_topo_iso == 1)      recoisoET = cluster_iso_topo_03[icl];
                    else if (use_topo_iso == 4)      recoisoET = cluster_iso_04[icl];
                    else if (use_topo_iso == 5)      recoisoET = cluster_iso_excl_04[icl];
                    else if (conesize == 4)          recoisoET = cluster_iso_04[icl];
                    else if (conesize == 3)          recoisoET = cluster_iso_03[icl];
                    else if (conesize == 2)          recoisoET = cluster_iso_02[icl];
                    recoisoET = recoisoET * mc_iso_scale + mc_iso_shift;
                    float recoiso_max = recoiso_max_b + recoiso_max_s * clusterET;
                    passes_iso = (recoisoET > recoiso_min && recoisoET < recoiso_max);

                    float tight_weta_cogx_max = tight_weta_cogx_max_b + tight_weta_cogx_max_s * clusterET;
                    float tight_wphi_cogx_max = tight_wphi_cogx_max_b + tight_wphi_cogx_max_s * clusterET;
                    float tight_et1_min       = tight_et1_min_b + tight_et1_min_s * clusterET;

                    std::string selected_bdt_model = bdt_model_name;
                    if (use_et_binned_bdt) {
                        for (int ib = 0; ib < (int)bdt_et_bin_models.size(); ++ib) {
                            if (clusterET >= bdt_et_bin_edges[ib] && clusterET < bdt_et_bin_edges[ib+1]) {
                                selected_bdt_model = bdt_et_bin_models[ib];
                                break;
                            }
                        }
                    }
                    float bdt_score = (*bdt_arrays[selected_bdt_model])[icl];
                    float tight_bdt_min_et = tight_bdt_min_slope * clusterET + tight_bdt_min_intercept;

                    passes_tight =
                        (cluster_weta_cogx[icl] > tight_weta_cogx_min) &&
                        (cluster_weta_cogx[icl] < tight_weta_cogx_max) &&
                        (cluster_wphi_cogx[icl] > tight_wphi_cogx_min) &&
                        (cluster_wphi_cogx[icl] < tight_wphi_cogx_max) &&
                        (cluster_et1[icl] > tight_et1_min) &&
                        (cluster_et1[icl] < tight_et1_max) &&
                        (cluster_et2[icl] > tight_et2_min) &&
                        (cluster_et2[icl] < tight_et2_max) &&
                        (cluster_et3[icl] > tight_et3_min) &&
                        (cluster_et3[icl] < tight_et3_max) &&
                        (e11_over_e33 > tight_e11_over_e33_min) &&
                        (e11_over_e33 < tight_e11_over_e33_max) &&
                        (e32_over_e35 > tight_e32_over_e35_min) &&
                        (e32_over_e35 < tight_e32_over_e35_max) &&
                        (cluster_et4[icl] > tight_et4_min) &&
                        (cluster_et4[icl] < tight_et4_max) &&
                        (cluster_prob[icl] > tight_prob_min) &&
                        (cluster_prob[icl] < tight_prob_max) &&
                        (bdt_score > tight_bdt_min_et) &&
                        (bdt_score < tight_bdt_max);
                }

                // Fill every set with non-zero alpha
                for (auto &s : sets) {
                    float alpha = (isamp.which == 0) ? s.alpha_nom : s.alpha_double;
                    if (alpha <= 0) continue;
                    float w = alpha;
                    if (s.do_vtx_reweight && s.h_vtx_weight) {
                        const float vlo = s.h_vtx_weight->GetXaxis()->GetXmin();
                        const float vhi = s.h_vtx_weight->GetXaxis()->GetXmax();
                        if (vtx > vlo && vtx < vhi) {
                            w *= s.h_vtx_weight->Interpolate(vtx);
                        }  // else weight unchanged — out of range, trust alpha
                    }

                    // reco level: trkID match only, no additional cuts
                    hRespET[Form("h_respET_%s_reco", s.name.Data())]->Fill(truth_pT, rET, w);
                    hRespE [Form("h_respE_%s_reco",  s.name.Data())]->Fill(truth_pT, rE,  w);

                    if (passes_common) {
                        hRespET[Form("h_respET_%s_common_cut", s.name.Data())]->Fill(truth_pT, rET, w);
                        hRespE [Form("h_respE_%s_common_cut",  s.name.Data())]->Fill(truth_pT, rE,  w);

                        if (passes_iso && passes_tight) {
                            hRespET[Form("h_respET_%s_tight_iso", s.name.Data())]->Fill(truth_pT, rET, w);
                            hRespE [Form("h_respE_%s_tight_iso",  s.name.Data())]->Fill(truth_pT, rE,  w);
                        }
                    }
                }
            }  // truth particles
        }  // events

        for (auto &kv : bdt_arrays) delete kv.second;
        fin->Close();
        delete fin;
    }

    // -----------------------------------------------------------------
    // Summary stats: for each (set, level), extract mean/sigma per pT bin
    // by two independent fits:
    //   (a) truncated Gaussian on [0.3, 1.7] — core resolution, legacy default
    //       stored in  h_{resp}_{mean,sigma}_{set}_{level}
    //   (b) Crystal-Ball  on [0.3, 2.0]      — captures low-side tail from
    //       brem / shower leakage / double-interaction kinematic shift
    //       stored in  h_{resp}_cb_{mean,sigma}_{set}_{level}
    //       (plus h_{resp}_cb_{alpha,n}_{set}_{level} for diagnostics)
    // Fallbacks: if either fit fails, the histogram mean/RMS is stored.
    // -----------------------------------------------------------------
    std::map<TString, TH1F*> hMeanET, hSigmaET, hMeanE, hSigmaE;
    std::map<TString, TH1F*> hMeanET_cb, hSigmaET_cb, hMeanE_cb, hSigmaE_cb;
    std::map<TString, TH1F*> hAlphaET_cb, hNET_cb;  // diagnostics (ET only)
    // Double-sided CB (for tight_iso level where the one-sided CB misses
    // a residual high-side tail from isolation leakage; populated for all
    // levels for consistency)
    std::map<TString, TH1F*> hMeanET_dscb, hSigmaET_dscb;
    std::map<TString, TH1F*> hAlphaLET_dscb, hAlphaHET_dscb;
    auto makeSummary = [&](const char *prefix, const TString &setname, const TString &lvl)
    {
        TString n_m = Form("%s_mean_%s_%s",  prefix, setname.Data(), lvl.Data());
        TString n_s = Form("%s_sigma_%s_%s", prefix, setname.Data(), lvl.Data());
        TH1F *hm = new TH1F(n_m.Data(), Form(";truth p_{T} [GeV];<r>"),      NptBins, ptEdges);
        TH1F *hs = new TH1F(n_s.Data(), Form(";truth p_{T} [GeV];#sigma(r)"), NptBins, ptEdges);
        hm->Sumw2(); hs->Sumw2();
        return std::make_pair(hm, hs);
    };
    auto makeDiag = [&](const char *prefix, const TString &setname, const TString &lvl,
                        const char *suffix, const char *ytit)
    {
        TString n = Form("%s_%s_%s_%s", prefix, suffix, setname.Data(), lvl.Data());
        TH1F *h = new TH1F(n.Data(), Form(";truth p_{T} [GeV];%s", ytit), NptBins, ptEdges);
        h->Sumw2();
        return h;
    };

    auto fitSliceStats = [&](TH2F *h2,
                             TH1F *hm,    TH1F *hs,
                             TH1F *hm_cb, TH1F *hs_cb,
                             TH1F *ha_cb = nullptr, TH1F *hn_cb = nullptr,
                             TH1F *hm_dscb = nullptr, TH1F *hs_dscb = nullptr,
                             TH1F *haL_dscb = nullptr, TH1F *haH_dscb = nullptr)
    {
        for (int ipt = 1; ipt <= NptBins; ipt++) {
            TH1D *proj = h2->ProjectionY(Form("_tmp_%s_pt%d", h2->GetName(), ipt), ipt, ipt);
            if (proj->GetEntries() < 20 || proj->Integral() <= 0) { delete proj; continue; }

            // --- (a) Truncated Gaussian fit in [0.3, 1.7] ---
            TF1 f(Form("fresp_%s_pt%d", h2->GetName(), ipt), "gaus", 0.3, 1.7);
            f.SetParameter(0, proj->GetMaximum());
            f.SetParameter(1, proj->GetMean());
            f.SetParameter(2, std::max(0.05, (double)proj->GetStdDev()));
            int rc = proj->Fit(&f, "RQN");
            double mu_g, sg_g, dmu_g, dsg_g;
            if (rc == 0) {
                mu_g = f.GetParameter(1); dmu_g = f.GetParError(1);
                sg_g = f.GetParameter(2); dsg_g = f.GetParError(2);
            } else {
                mu_g = proj->GetMean();   dmu_g = proj->GetMeanError();
                sg_g = proj->GetStdDev(); dsg_g = proj->GetStdDevError();
            }
            hm->SetBinContent(ipt, mu_g); hm->SetBinError(ipt, dmu_g);
            hs->SetBinContent(ipt, sg_g); hs->SetBinError(ipt, dsg_g);

            // --- (b) Crystal-Ball fit in [0.3, 2.0] ---
            // ROOT "crystalball": [0]=const, [1]=mean, [2]=sigma, [3]=alpha (>0 -> low-side tail), [4]=n
            TF1 fcb(Form("fcb_%s_pt%d", h2->GetName(), ipt), "crystalball", 0.3, 2.0);
            fcb.SetParameter(0, proj->GetMaximum());
            fcb.SetParameter(1, mu_g);
            fcb.SetParameter(2, std::max(0.03, sg_g));
            fcb.SetParameter(3, 1.0);
            fcb.SetParameter(4, 5.0);
            fcb.SetParLimits(2, 0.005, 0.4);
            fcb.SetParLimits(3, 0.1, 10.0);
            fcb.SetParLimits(4, 1.05, 100.0);
            int rc_cb = proj->Fit(&fcb, "RQN");
            double mu_cb = mu_g, sg_cb = sg_g;
            if (rc_cb == 0) {
                mu_cb = fcb.GetParameter(1);
                sg_cb = fcb.GetParameter(2);
                hm_cb->SetBinContent(ipt, mu_cb);
                hm_cb->SetBinError  (ipt, fcb.GetParError(1));
                hs_cb->SetBinContent(ipt, sg_cb);
                hs_cb->SetBinError  (ipt, fcb.GetParError(2));
                if (ha_cb) { ha_cb->SetBinContent(ipt, fcb.GetParameter(3)); ha_cb->SetBinError(ipt, fcb.GetParError(3)); }
                if (hn_cb) { hn_cb->SetBinContent(ipt, fcb.GetParameter(4)); hn_cb->SetBinError(ipt, fcb.GetParError(4)); }
            } else {
                hm_cb->SetBinContent(ipt, mu_g); hm_cb->SetBinError(ipt, dmu_g);
                hs_cb->SetBinContent(ipt, sg_g); hs_cb->SetBinError(ipt, dsg_g);
            }

            // --- (c) Double-sided Crystal-Ball fit in [0.3, 2.0] ---
            if (hm_dscb && hs_dscb) {
                TF1 fdscb(Form("fdscb_%s_pt%d", h2->GetName(), ipt), double_cb, 0.3, 2.0, 7);
                fdscb.SetParNames("N", "mu", "sigma", "alphaL", "nL", "alphaH", "nH");
                fdscb.SetParameter(0, proj->GetMaximum());
                fdscb.SetParameter(1, mu_cb);
                fdscb.SetParameter(2, std::max(0.03, sg_cb));
                fdscb.SetParameter(3, 1.0);
                fdscb.SetParameter(4, 5.0);
                fdscb.SetParameter(5, 1.5);
                fdscb.SetParameter(6, 5.0);
                fdscb.SetParLimits(2, 0.005, 0.4);
                fdscb.SetParLimits(3, 0.1, 10.0);
                fdscb.SetParLimits(4, 1.05, 100.0);
                fdscb.SetParLimits(5, 0.1, 10.0);
                fdscb.SetParLimits(6, 1.05, 100.0);
                int rc_dscb = proj->Fit(&fdscb, "RQN");
                if (rc_dscb == 0) {
                    hm_dscb->SetBinContent(ipt, fdscb.GetParameter(1));
                    hm_dscb->SetBinError  (ipt, fdscb.GetParError(1));
                    hs_dscb->SetBinContent(ipt, fdscb.GetParameter(2));
                    hs_dscb->SetBinError  (ipt, fdscb.GetParError(2));
                    if (haL_dscb) {
                        haL_dscb->SetBinContent(ipt, fdscb.GetParameter(3));
                        haL_dscb->SetBinError  (ipt, fdscb.GetParError(3));
                    }
                    if (haH_dscb) {
                        haH_dscb->SetBinContent(ipt, fdscb.GetParameter(5));
                        haH_dscb->SetBinError  (ipt, fdscb.GetParError(5));
                    }
                } else {
                    // Fallback: inherit CB (or Gauss if CB also failed)
                    hm_dscb->SetBinContent(ipt, mu_cb); hm_dscb->SetBinError(ipt, dmu_g);
                    hs_dscb->SetBinContent(ipt, sg_cb); hs_dscb->SetBinError(ipt, dsg_g);
                }
            }

            delete proj;
        }
    };

    for (auto &s : sets) {
        for (auto &lvl : levels) {
            auto smET      = makeSummary("h_respET",      s.name, lvl);
            auto smE       = makeSummary("h_respE",       s.name, lvl);
            auto smET_cb   = makeSummary("h_respET_cb",   s.name, lvl);
            auto smE_cb    = makeSummary("h_respE_cb",    s.name, lvl);
            auto smET_dscb = makeSummary("h_respET_dscb", s.name, lvl);
            TH1F *ha   = makeDiag("h_respET_cb",   s.name, lvl, "alpha",  "#alpha_{CB}");
            TH1F *hn   = makeDiag("h_respET_cb",   s.name, lvl, "n",      "n_{CB}");
            TH1F *haL  = makeDiag("h_respET_dscb", s.name, lvl, "alphaL", "#alpha_{L}");
            TH1F *haH  = makeDiag("h_respET_dscb", s.name, lvl, "alphaH", "#alpha_{H}");
            TString k1 = Form("h_respET_%s_%s", s.name.Data(), lvl.Data());
            TString k2 = Form("h_respE_%s_%s",  s.name.Data(), lvl.Data());
            fitSliceStats(hRespET[k1], smET.first, smET.second,
                          smET_cb.first, smET_cb.second, ha, hn,
                          smET_dscb.first, smET_dscb.second, haL, haH);
            fitSliceStats(hRespE [k2], smE.first,  smE.second,
                          smE_cb.first,  smE_cb.second);
            hMeanET       [smET.first      ->GetName()] = smET.first;
            hSigmaET      [smET.second     ->GetName()] = smET.second;
            hMeanE        [smE.first       ->GetName()] = smE.first;
            hSigmaE       [smE.second      ->GetName()] = smE.second;
            hMeanET_cb    [smET_cb.first   ->GetName()] = smET_cb.first;
            hSigmaET_cb   [smET_cb.second  ->GetName()] = smET_cb.second;
            hMeanE_cb     [smE_cb.first    ->GetName()] = smE_cb.first;
            hSigmaE_cb    [smE_cb.second   ->GetName()] = smE_cb.second;
            hAlphaET_cb   [ha->GetName()] = ha;
            hNET_cb       [hn->GetName()] = hn;
            hMeanET_dscb  [smET_dscb.first ->GetName()] = smET_dscb.first;
            hSigmaET_dscb [smET_dscb.second->GetName()] = smET_dscb.second;
            hAlphaLET_dscb[haL->GetName()] = haL;
            hAlphaHET_dscb[haH->GetName()] = haH;
        }
    }

    // -----------------------------------------------------------------
    // Write output
    // -----------------------------------------------------------------
    TString outFilename = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/energy_response_comparison.root";
    TFile *fout = new TFile(outFilename, "RECREATE");

    h_vtxz_single->Write();
    h_vtxz_double->Write();
    for (auto &s : sets) {
        if (s.h_vtx_weight) {
            s.h_vtx_weight->SetName(Form("h_vtx_weight_%s", s.name.Data()));
            s.h_vtx_weight->Write();
        }
    }
    for (auto &kv : hRespET)     kv.second->Write();
    for (auto &kv : hRespE)      kv.second->Write();
    for (auto &kv : hMeanET)     kv.second->Write();
    for (auto &kv : hSigmaET)    kv.second->Write();
    for (auto &kv : hMeanE)      kv.second->Write();
    for (auto &kv : hSigmaE)     kv.second->Write();
    for (auto &kv : hMeanET_cb)     kv.second->Write();
    for (auto &kv : hSigmaET_cb)    kv.second->Write();
    for (auto &kv : hMeanE_cb)      kv.second->Write();
    for (auto &kv : hSigmaE_cb)     kv.second->Write();
    for (auto &kv : hAlphaET_cb)    kv.second->Write();
    for (auto &kv : hNET_cb)        kv.second->Write();
    for (auto &kv : hMeanET_dscb)   kv.second->Write();
    for (auto &kv : hSigmaET_dscb)  kv.second->Write();
    for (auto &kv : hAlphaLET_dscb) kv.second->Write();
    for (auto &kv : hAlphaHET_dscb) kv.second->Write();
    fout->Close();
    delete fout;

    std::cout << "\n[done] output: " << outFilename << std::endl;

    // -----------------------------------------------------------------
    // Integrated summary (pT-integrated mean/RMS of rET at common_cut level,
    // i.e. post common cuts; the working point that enters ABCD)
    // -----------------------------------------------------------------
    std::cout << "\nIntegrated response (rET = ET_reco/pT_truth, level=common_cut):" << std::endl;
    std::cout << std::string(78, '-') << std::endl;
    std::cout << Form("%-22s  %8s  %8s  %12s",
                      "set", "<rET>", "sigma", "entries") << std::endl;
    std::cout << std::string(78, '-') << std::endl;
    for (auto &s : sets) {
        TH2F *h = hRespET[Form("h_respET_%s_common_cut", s.name.Data())];
        TH1D *proj = h->ProjectionY("_proj_common_cut");
        double m  = proj->GetMean();
        double sd = proj->GetStdDev();
        double n  = proj->Integral();
        std::cout << Form("%-22s  %8.4f  %8.4f  %12.1f",
                          s.name.Data(), m, sd, n) << std::endl;
        delete proj;
    }
}
