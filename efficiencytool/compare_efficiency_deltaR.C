// compare_efficiency_deltaR.C
// Compare per-stage photon efficiency with and without the deltaR < 0.1
// truth-matching cut, for single-interaction MC and mixed double-interaction
// MC (0 mrad).  Special attention to the photon ID (BDT) efficiency.
//
// Usage:
//   root -l -b -q 'compare_efficiency_deltaR.C("config_bdt_nom.yaml", -1)'
//   root -l -b -q 'compare_efficiency_deltaR.C("config_bdt_nom.yaml", 500000)'

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
#include <TEfficiency.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TMath.h>
#include <TROOT.h>

#include <yaml-cpp/yaml.h>

// =====================================================================
// Sample descriptor
// =====================================================================
struct SampleInfo {
    TString name;
    TString path;
    float mix_weight;
};

// =====================================================================
// Helper: compute integrated efficiency from a TEfficiency
// =====================================================================
static double IntegratedEfficiency(TEfficiency *eff)
{
    if (!eff) return -1;
    TH1 *hPass = (TH1*)eff->GetPassedHistogram();
    TH1 *hTot  = (TH1*)eff->GetTotalHistogram();
    double nPass = 0, nTot = 0;
    for (int i = 1; i <= hTot->GetNbinsX(); i++) {
        nPass += hPass->GetBinContent(i);
        nTot  += hTot->GetBinContent(i);
    }
    return nTot > 0 ? nPass / nTot : 0;
}

// =====================================================================
// Main function
// =====================================================================
void compare_efficiency_deltaR(TString configname = "config_bdt_nom.yaml",
                               int maxEvents = -1)
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");

    // =================================================================
    // Load configuration
    // =================================================================
    YAML::Node configYaml = YAML::LoadFile(configname.Data());

    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();
    std::string bdt_model_name  = configYaml["input"]["bdt_model_name"].as<std::string>("base");

    // ET-binned BDT models
    std::vector<float>       bdt_et_bin_edges;
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

    // Analysis cuts
    int use_topo_iso       = configYaml["analysis"]["use_topo_iso"].as<int>(0);
    int iso_threshold      = configYaml["analysis"]["iso_threshold"].as<int>(0);
    int iso_hcalonly       = configYaml["analysis"]["iso_hcalonly"].as<int>(0);
    float iso_emcalinnerr  = configYaml["analysis"]["iso_emcalinnerr"].as<float>(0.0);
    int conesize           = configYaml["analysis"]["cone_size"].as<int>();
    float truthisocut      = configYaml["analysis"]["truth_iso_max"].as<float>();
    float recoiso_min      = configYaml["analysis"]["reco_iso_min"].as<float>();
    float recoiso_max_b    = configYaml["analysis"]["reco_iso_max_b"].as<float>();
    float recoiso_max_s    = configYaml["analysis"]["reco_iso_max_s"].as<float>();
    float vertexcut        = configYaml["analysis"]["vertex_cut"].as<float>();
    float reco_min_ET      = configYaml["analysis"]["reco_min_ET"].as<float>();
    float eff_dR           = configYaml["analysis"]["eff_dR"].as<float>();
    float mc_iso_shift     = configYaml["analysis"]["mc_iso_shift"].as<float>(0.0);
    float mc_iso_scale     = configYaml["analysis"]["mc_iso_scale"].as<float>(1.0);

    std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();

    // pT bins for truth efficiency
    std::vector<float> pT_bins_truth;
    if (configYaml["analysis"]["pT_bins_truth"]) {
        pT_bins_truth = configYaml["analysis"]["pT_bins_truth"].as<std::vector<float>>();
    } else {
        pT_bins_truth = {7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45};
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

    // Tight BDT cuts (parametric)
    float tight_bdt_max           = configYaml["analysis"]["tight"]["bdt_max"].as<float>(1);
    float tight_bdt_min           = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0);
    float tight_bdt_min_slope     = configYaml["analysis"]["tight"]["bdt_min_slope"].as<float>(0);
    float tight_bdt_min_intercept = configYaml["analysis"]["tight"]["bdt_min_intercept"].as<float>(tight_bdt_min);

    // Tight shower-shape cuts
    float tight_weta_cogx_min   = configYaml["analysis"]["tight"]["weta_cogx_min"].as<float>();
    float tight_weta_cogx_max_b = configYaml["analysis"]["tight"]["weta_cogx_max_b"].as<float>();
    float tight_weta_cogx_max_s = configYaml["analysis"]["tight"]["weta_cogx_max_s"].as<float>();
    float tight_wphi_cogx_min   = configYaml["analysis"]["tight"]["wphi_cogx_min"].as<float>();
    float tight_wphi_cogx_max_b = configYaml["analysis"]["tight"]["wphi_cogx_max_b"].as<float>();
    float tight_wphi_cogx_max_s = configYaml["analysis"]["tight"]["wphi_cogx_max_s"].as<float>();
    float tight_e11_over_e33_min = configYaml["analysis"]["tight"]["e11_over_e33_min"].as<float>();
    float tight_e11_over_e33_max = configYaml["analysis"]["tight"]["e11_over_e33_max"].as<float>();
    float tight_e32_over_e35_min = configYaml["analysis"]["tight"]["e32_over_e35_min"].as<float>();
    float tight_e32_over_e35_max = configYaml["analysis"]["tight"]["e32_over_e35_max"].as<float>();
    float tight_et1_max          = configYaml["analysis"]["tight"]["et1_max"].as<float>();
    float tight_et1_min_b        = configYaml["analysis"]["tight"]["et1_min_b"].as<float>();
    float tight_et1_min_s        = configYaml["analysis"]["tight"]["et1_min_s"].as<float>();
    float tight_et2_min          = configYaml["analysis"]["tight"]["et2_min"].as<float>(0.0);
    float tight_et2_max          = configYaml["analysis"]["tight"]["et2_max"].as<float>(1.0);
    float tight_et3_min          = configYaml["analysis"]["tight"]["et3_min"].as<float>(0.0);
    float tight_et3_max          = configYaml["analysis"]["tight"]["et3_max"].as<float>(1.0);
    float tight_et4_min          = configYaml["analysis"]["tight"]["et4_min"].as<float>();
    float tight_et4_max          = configYaml["analysis"]["tight"]["et4_max"].as<float>();
    float tight_prob_min         = configYaml["analysis"]["tight"]["prob_min"].as<float>();
    float tight_prob_max         = configYaml["analysis"]["tight"]["prob_max"].as<float>();
    float tight_w32_min          = configYaml["analysis"]["tight"]["w32_min"].as<float>();
    float tight_w32_max          = configYaml["analysis"]["tight"]["w32_max"].as<float>();

    std::cout << "eff_dR threshold: " << eff_dR << std::endl;
    std::cout << "tight_bdt_min_intercept: " << tight_bdt_min_intercept
              << "  tight_bdt_min_slope: " << tight_bdt_min_slope << std::endl;

    // =================================================================
    // Define samples
    // =================================================================
    float DOUBLE_FRAC = 0.187f;
    float SINGLE_FRAC = 1.0f - DOUBLE_FRAC; // 0.813

    std::vector<SampleInfo> single_samples = {
        {"photon10", "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root", 1.0f},
    };

    std::vector<SampleInfo> double_samples = {
        {"photon10_double", "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_double/bdt_split.root", 1.0f},
    };

    // Mixed: same files but with pileup-fraction weights
    std::vector<SampleInfo> mixed_samples = {
        {"photon10",        "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root",        SINGLE_FRAC},
        {"photon10_double", "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_double/bdt_split.root", DOUBLE_FRAC},
    };

    // =================================================================
    // Efficiency bookkeeping
    // =================================================================
    // We track 3 sample sets x 2 matching modes x 4 stages = 24 TEfficiency,
    // plus 3 "newclusters" TEfficiency objects.
    const std::vector<TString> stages   = {"reco", "iso", "id", "all"};
    const std::vector<TString> matchings = {"withDR", "noDR"};
    const std::vector<TString> setnames  = {"single", "double", "mixed"};

    std::map<TString, TEfficiency*> effMap;

    auto makeEff = [&](const TString &name) -> TEfficiency* {
        TEfficiency *e = new TEfficiency(name.Data(),
            Form("%s;truth p_{T} [GeV];Efficiency", name.Data()),
            NptBins, ptEdges);
        e->SetUseWeightedEvents();
        effMap[name] = e;
        return e;
    };

    // Create all 24 standard efficiencies
    for (auto &sset : setnames)
        for (auto &stage : stages)
            for (auto &match : matchings)
                makeEff(Form("eff_%s_%s_%s", stage.Data(), match.Data(), sset.Data()));

    // Create 3 "newclusters" efficiencies
    for (auto &sset : setnames)
        makeEff(Form("eff_id_newclusters_%s", sset.Data()));

    // Also create a deltaR distribution histogram for diagnostics
    TH1F *h_deltaR_single = new TH1F("h_deltaR_single", "min #DeltaR (trkID matched), single;#DeltaR;counts", 200, 0, 0.4);
    TH1F *h_deltaR_double = new TH1F("h_deltaR_double", "min #DeltaR (trkID matched), double;#DeltaR;counts", 200, 0, 0.4);
    h_deltaR_single->Sumw2();
    h_deltaR_double->Sumw2();

    // =================================================================
    // Process samples
    // =================================================================
    // We process 3 sets.  For "mixed", we process two files with different
    // weights and accumulate into the same TEfficiency objects.

    struct SetSpec {
        TString setname;
        std::vector<SampleInfo> samples;
    };

    std::vector<SetSpec> allSets = {
        {"single", single_samples},
        {"double", double_samples},
        {"mixed",  mixed_samples},
    };

    for (auto &spec : allSets)
    {
        TString setname = spec.setname;
        std::cout << "\n========================================"  << std::endl;
        std::cout << "Processing sample set: " << setname         << std::endl;
        std::cout << "========================================"  << std::endl;

        for (auto &si : spec.samples)
        {
            std::cout << "\n  Sample: " << si.name << "  path: " << si.path
                      << "  mix_weight: " << si.mix_weight << std::endl;

            TFile *fin = TFile::Open(si.path);
            if (!fin || fin->IsZombie()) {
                std::cerr << "ERROR: Cannot open " << si.path << std::endl;
                continue;
            }

            TTreeReader reader("slimtree", fin);

            // --- Event-level ---
            TTreeReaderValue<int>   nparticles(reader, "nparticles");
            TTreeReaderValue<int>   nclusters(reader, Form("ncluster_%s", clusternodename.c_str()));
            TTreeReaderValue<float> vertexz(reader, "vertexz");
            TTreeReaderValue<float> vertexz_truth(reader, "vertexz_truth");

            // --- Truth particles ---
            TTreeReaderArray<float> particle_Pt(reader, "particle_Pt");
            TTreeReaderArray<float> particle_Eta(reader, "particle_Eta");
            TTreeReaderArray<float> particle_Phi(reader, "particle_Phi");
            TTreeReaderArray<int>   particle_pid(reader, "particle_pid");
            TTreeReaderArray<int>   particle_trkid(reader, "particle_trkid");
            TTreeReaderArray<int>   particle_photonclass(reader, "particle_photonclass");
            TTreeReaderArray<float> particle_truth_iso_03(reader, "particle_truth_iso_03");
            TTreeReaderArray<float> particle_truth_iso_04(reader, "particle_truth_iso_04");

            // --- Clusters ---
            TTreeReaderArray<float> cluster_Et(reader,  Form("cluster_Et_%s",  clusternodename.c_str()));
            TTreeReaderArray<float> cluster_E(reader,   Form("cluster_E_%s",   clusternodename.c_str()));
            TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
            TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusternodename.c_str()));
            TTreeReaderArray<int>   cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusternodename.c_str()));
            TTreeReaderArray<float> cluster_prob(reader, Form("cluster_prob_%s", clusternodename.c_str()));

            // Shower-shape variables
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
            TTreeReaderArray<float> cluster_w32(reader, Form("cluster_w32_%s", clusternodename.c_str()));
            TTreeReaderArray<float> cluster_npb_score(reader, Form("cluster_npb_score_%s", clusternodename.c_str()));

            // Isolation branches
            TTreeReaderArray<float> cluster_iso_02(reader, Form("cluster_iso_02_%s", clusternodename.c_str()));
            TTreeReaderArray<float> cluster_iso_03(reader, Form("cluster_iso_03_%s", clusternodename.c_str()));
            TTreeReaderArray<float> cluster_iso_04(reader, Form("cluster_iso_04_%s", clusternodename.c_str()));
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

            // --- Get the TEfficiency pointers for this set ---
            TEfficiency *eff_reco_withDR = effMap[Form("eff_reco_withDR_%s", setname.Data())];
            TEfficiency *eff_reco_noDR   = effMap[Form("eff_reco_noDR_%s",   setname.Data())];
            TEfficiency *eff_iso_withDR  = effMap[Form("eff_iso_withDR_%s",  setname.Data())];
            TEfficiency *eff_iso_noDR    = effMap[Form("eff_iso_noDR_%s",    setname.Data())];
            TEfficiency *eff_id_withDR   = effMap[Form("eff_id_withDR_%s",   setname.Data())];
            TEfficiency *eff_id_noDR     = effMap[Form("eff_id_noDR_%s",     setname.Data())];
            TEfficiency *eff_all_withDR  = effMap[Form("eff_all_withDR_%s",  setname.Data())];
            TEfficiency *eff_all_noDR    = effMap[Form("eff_all_noDR_%s",    setname.Data())];
            TEfficiency *eff_id_newclusters = effMap[Form("eff_id_newclusters_%s", setname.Data())];

            TH1F *h_deltaR_diag = (setname == "single") ? h_deltaR_single : h_deltaR_double;

            // --- Event loop ---
            long long nEvents_processed = 0;
            while (reader.Next())
            {
                if (maxEvents > 0 && nEvents_processed >= maxEvents) break;
                nEvents_processed++;

                if (nEvents_processed % 200000 == 0)
                    std::cout << "    Event " << nEvents_processed << std::endl;

                float vtx_reco = *vertexz;
                if (fabs(vtx_reco) > vertexcut) continue;

                float weight = si.mix_weight;
                int ncl   = *nclusters;
                int npart = *nparticles;

                // Build trkID -> list of cluster indices
                std::map<int, std::vector<int>> trkid_to_clusters;
                for (int icl = 0; icl < ncl; icl++) {
                    if (cluster_Et[icl] < reco_min_ET) continue; // skip low-ET clusters
                    trkid_to_clusters[cluster_truthtrkID[icl]].push_back(icl);
                }

                // Loop over truth particles
                for (int ip = 0; ip < npart; ip++)
                {
                    if (particle_pid[ip] != 22) continue;
                    if (particle_photonclass[ip] >= 3) continue;

                    // Truth isolation
                    float truthisoET = (conesize == 4) ? particle_truth_iso_04[ip]
                                     : particle_truth_iso_03[ip];
                    if (truthisoET >= truthisocut) continue;

                    // Eta acceptance
                    float truth_eta = particle_Eta[ip];
                    bool in_eta = false;
                    for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++) {
                        if (truth_eta > eta_bins[ieta] && truth_eta < eta_bins[ieta + 1]) {
                            in_eta = true;
                            break;
                        }
                    }
                    if (!in_eta) continue;

                    // pT range
                    float truth_pT = particle_Pt[ip];
                    if (truth_pT < pTmin_truth || truth_pT >= pTmax_truth) continue;

                    int truth_trkid = particle_trkid[ip];
                    float truth_phi = particle_Phi[ip];

                    // -----------------------------------------------
                    // Find best trkID-matched cluster (closest deltaR)
                    // -----------------------------------------------
                    auto itMap = trkid_to_clusters.find(truth_trkid);
                    bool has_trkid = (itMap != trkid_to_clusters.end());

                    float best_dR = 999;
                    int best_icl = -1;

                    if (has_trkid) {
                        for (int icl : itMap->second) {
                            float deta = cluster_Eta[icl] - truth_eta;
                            float dphi = cluster_Phi[icl] - truth_phi;
                            if (dphi > TMath::Pi())  dphi -= 2 * TMath::Pi();
                            if (dphi < -TMath::Pi()) dphi += 2 * TMath::Pi();
                            float dR = sqrt(deta * deta + dphi * dphi);
                            if (dR < best_dR) {
                                best_dR = dR;
                                best_icl = icl;
                            }
                        }
                    }

                    bool passes_dR = has_trkid && (best_dR < eff_dR);

                    // Fill deltaR diagnostic histogram (only for single/double pure sets)
                    if (has_trkid && setname != "mixed") {
                        h_deltaR_diag->Fill(best_dR, weight);
                    }

                    // -----------------------------------------------
                    // Evaluate cuts on the best-matched cluster
                    // -----------------------------------------------
                    bool passes_common = false;
                    bool passes_iso    = false;
                    bool passes_tight  = false;

                    if (has_trkid && best_icl >= 0)
                    {
                        int icl = best_icl;
                        float clusterET = cluster_Et[icl];

                        // Derived shower-shape variables
                        float e11_over_e33 = cluster_e11[icl] / cluster_e33[icl];
                        float e32_over_e35 = cluster_e32[icl] / cluster_e35[icl];
                        float wr_cogx      = cluster_wphi_cogx[icl] / cluster_weta_cogx[icl];

                        // Common cuts (same logic as RecoEffCalculator_TTreeReader.C)
                        passes_common =
                            cluster_prob[icl] > common_prob_min &&
                            cluster_prob[icl] < common_prob_max &&
                            e11_over_e33 > common_e11_over_e33_min &&
                            e11_over_e33 < common_e11_over_e33_max &&
                            wr_cogx > common_wr_cogx_bound &&
                            cluster_weta_cogx[icl] < common_cluster_weta_cogx_bound &&
                            (!common_npb_cut_on || cluster_npb_score[icl] > common_npb_score_cut);

                        if (passes_common)
                        {
                            // --- Isolation ---
                            float recoisoET = -999;
                            if (use_topo_iso == 2) {
                                recoisoET = cluster_iso_topo_04[icl];
                            } else if (use_topo_iso == 1) {
                                recoisoET = cluster_iso_topo_03[icl];
                            } else if (conesize == 4) {
                                recoisoET = cluster_iso_04[icl];
                            } else if (conesize == 3) {
                                recoisoET = cluster_iso_03[icl];
                            } else if (conesize == 2) {
                                recoisoET = cluster_iso_02[icl];
                            }

                            // MC isolation fudge
                            recoisoET = recoisoET * mc_iso_scale + mc_iso_shift;

                            float recoiso_max = recoiso_max_b + recoiso_max_s * clusterET;
                            passes_iso = (recoisoET > recoiso_min && recoisoET < recoiso_max);

                            // --- Tight (BDT + shower shape) ---
                            // ET-dependent thresholds
                            float tight_weta_cogx_max = tight_weta_cogx_max_b + tight_weta_cogx_max_s * clusterET;
                            float tight_wphi_cogx_max = tight_wphi_cogx_max_b + tight_wphi_cogx_max_s * clusterET;
                            float tight_et1_min       = tight_et1_min_b + tight_et1_min_s * clusterET;

                            bool is_weta_tight = (cluster_weta_cogx[icl] > tight_weta_cogx_min) &&
                                                 (cluster_weta_cogx[icl] < tight_weta_cogx_max);
                            bool is_wphi_tight = (cluster_wphi_cogx[icl] > tight_wphi_cogx_min) &&
                                                 (cluster_wphi_cogx[icl] < tight_wphi_cogx_max);
                            bool is_et1_tight  = (cluster_et1[icl] > tight_et1_min) &&
                                                 (cluster_et1[icl] < tight_et1_max);
                            bool is_et2_tight  = (cluster_et2[icl] > tight_et2_min) &&
                                                 (cluster_et2[icl] < tight_et2_max);
                            bool is_et3_tight  = (cluster_et3[icl] > tight_et3_min) &&
                                                 (cluster_et3[icl] < tight_et3_max);
                            bool is_e11_tight  = (e11_over_e33 > tight_e11_over_e33_min) &&
                                                 (e11_over_e33 < tight_e11_over_e33_max);
                            bool is_e32_tight  = (e32_over_e35 > tight_e32_over_e35_min) &&
                                                 (e32_over_e35 < tight_e32_over_e35_max);
                            bool is_et4_tight  = (cluster_et4[icl] > tight_et4_min) &&
                                                 (cluster_et4[icl] < tight_et4_max);
                            bool is_prob_tight = (cluster_prob[icl] > tight_prob_min) &&
                                                 (cluster_prob[icl] < tight_prob_max);

                            // Select BDT model based on cluster ET
                            std::string selected_bdt_model = bdt_model_name;
                            if (use_et_binned_bdt) {
                                for (int ib = 0; ib < (int)bdt_et_bin_models.size(); ++ib) {
                                    if (clusterET >= bdt_et_bin_edges[ib] && clusterET < bdt_et_bin_edges[ib + 1]) {
                                        selected_bdt_model = bdt_et_bin_models[ib];
                                        break;
                                    }
                                }
                            }
                            float bdt_score = (*bdt_arrays[selected_bdt_model])[icl];
                            float tight_bdt_min_et = tight_bdt_min_slope * clusterET + tight_bdt_min_intercept;
                            bool is_bdt_tight = (bdt_score > tight_bdt_min_et) && (bdt_score < tight_bdt_max);

                            passes_tight = is_weta_tight && is_wphi_tight && is_et1_tight &&
                                           is_et2_tight && is_et3_tight && is_e11_tight &&
                                           is_e32_tight && is_et4_tight && is_prob_tight &&
                                           is_bdt_tight;
                        }
                    }

                    // -----------------------------------------------
                    // Build factorized efficiency booleans
                    // -----------------------------------------------
                    // WITH deltaR
                    bool reco_withDR = passes_dR && passes_common;
                    bool iso_withDR  = reco_withDR && passes_iso;
                    bool id_withDR   = iso_withDR && passes_tight;
                    bool all_withDR  = id_withDR;

                    // WITHOUT deltaR
                    bool reco_noDR = has_trkid && passes_common;
                    bool iso_noDR  = reco_noDR && passes_iso;
                    bool id_noDR   = iso_noDR && passes_tight;
                    bool all_noDR  = id_noDR;

                    // -----------------------------------------------
                    // Fill TEfficiency objects
                    // -----------------------------------------------
                    // eff_reco and eff_all: unconditional denominator (all truth photons)
                    eff_reco_withDR->FillWeighted(reco_withDR, weight, truth_pT);
                    eff_reco_noDR->FillWeighted(reco_noDR,     weight, truth_pT);
                    eff_all_withDR->FillWeighted(all_withDR,   weight, truth_pT);
                    eff_all_noDR->FillWeighted(all_noDR,       weight, truth_pT);

                    // eff_iso: conditional on reco=true
                    if (reco_withDR)
                        eff_iso_withDR->FillWeighted(passes_iso, weight, truth_pT);
                    if (reco_noDR)
                        eff_iso_noDR->FillWeighted(passes_iso, weight, truth_pT);

                    // eff_id: conditional on reco+iso=true
                    if (iso_withDR)
                        eff_id_withDR->FillWeighted(passes_tight, weight, truth_pT);
                    if (iso_noDR)
                        eff_id_noDR->FillWeighted(passes_tight, weight, truth_pT);

                    // "new clusters" = those in noDR but NOT in withDR
                    bool is_new = reco_noDR && !passes_dR;
                    if (is_new && passes_iso) {
                        eff_id_newclusters->FillWeighted(passes_tight, weight, truth_pT);
                    }

                } // end truth particle loop
            } // end event loop

            std::cout << "  Processed " << nEvents_processed << " events for "
                      << si.name << std::endl;

            // Clean up BDT reader arrays
            for (auto &kv : bdt_arrays) delete kv.second;
            fin->Close();
            delete fin;

        } // end sample loop
    } // end set loop

    // =================================================================
    // Write output
    // =================================================================
    TString outFilename = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/efficiency_comparison_deltaR.root";
    TFile *fout = new TFile(outFilename, "RECREATE");

    for (auto &kv : effMap) {
        kv.second->Write();
    }
    h_deltaR_single->Write();
    h_deltaR_double->Write();
    fout->Close();
    delete fout;

    std::cout << "\nOutput written to " << outFilename << std::endl;

    // =================================================================
    // Summary printout
    // =================================================================
    std::cout << "\n================================================================" << std::endl;
    std::cout << "  Efficiency comparison: with deltaR < " << eff_dR << "  vs  no deltaR cut" << std::endl;
    std::cout << "================================================================" << std::endl;

    for (auto &sset : setnames)
    {
        std::cout << "\n--- Sample set: " << sset << " ---" << std::endl;
        for (auto &stage : stages)
        {
            TEfficiency *e_withDR = effMap[Form("eff_%s_withDR_%s", stage.Data(), sset.Data())];
            TEfficiency *e_noDR  = effMap[Form("eff_%s_noDR_%s",  stage.Data(), sset.Data())];
            double eps_withDR = IntegratedEfficiency(e_withDR);
            double eps_noDR   = IntegratedEfficiency(e_noDR);
            double ratio = (eps_withDR > 0) ? eps_noDR / eps_withDR : 0;
            std::cout << Form("  %-6s  withDR = %.4f   noDR = %.4f   ratio(noDR/withDR) = %.4f",
                              stage.Data(), eps_withDR, eps_noDR, ratio)
                      << std::endl;
        }

        TEfficiency *e_new = effMap[Form("eff_id_newclusters_%s", sset.Data())];
        double eps_new = IntegratedEfficiency(e_new);
        std::cout << Form("  id(newclusters) = %.4f", eps_new) << std::endl;
    }

    std::cout << "\nDone." << std::endl;
}
