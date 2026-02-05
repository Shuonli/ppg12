#include <iostream>
#include <string>
#include <fstream>
#include <set>
#include <sstream>
#include <map>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TSystem.h>
#include <TProfile.h>
#include <yaml-cpp/yaml.h>
#include <cmath>

const float TIME_SAMPLE_NS = 17.6;

namespace
{
    TProfile *load_profile_clone(const std::string &filename, const std::string &profilename)
    {
        if (filename.empty())
            return nullptr;

        TFile *f = TFile::Open(filename.c_str(), "READ");
        if (!f || f->IsZombie())
        {
            std::cerr << "WARNING: Cannot open leading-time correction file: " << filename << std::endl;
            return nullptr;
        }

        TProfile *p = (TProfile *)f->Get(profilename.c_str());
        if (!p)
        {
            std::cerr << "WARNING: Cannot find TProfile '" << profilename << "' in file: " << filename << std::endl;
            f->Close();
            delete f;
            return nullptr;
        }

        // Clone so we can close the file safely
        TProfile *pclone = (TProfile *)p->Clone(Form("%s_clone", profilename.c_str()));
        pclone->SetDirectory(nullptr);

        f->Close();
        delete f;
        return pclone;
    }

    double correct_leading_time_ns(double raw_time_ns,
                                   double leading_adc,
                                   const TProfile *prof_time_vs_adc,
                                   double prof_global_mean_y)
    {
        if (!prof_time_vs_adc)
            return raw_time_ns;

        const int bin = prof_time_vs_adc->GetXaxis()->FindBin(leading_adc);
        if (bin < 1 || bin > prof_time_vs_adc->GetNbinsX())
            return raw_time_ns;

        // Only apply correction when that ADC bin actually has entries
        if (prof_time_vs_adc->GetBinEntries(bin) <= 0)
            return raw_time_ns;

        const double mean_time_at_adc = prof_time_vs_adc->GetBinContent(bin);

        // Remove ADC-dependent shift while preserving overall mean timing
        return raw_time_ns - (mean_time_at_adc - prof_global_mean_y);
    }
} // namespace

// Cross-section values (pb for photon, b for jet)
const float photon5cross = 146359.3;
const float photon10cross = 6944.675;
const float photon20cross = 130.4461;
const float jet10cross = 3.997e+06;
const float jet15cross = 4.073e+05;
const float jet20cross = 6.218e+04;
const float jet30cross = 2.502e+03;
const float jet50cross = 7.2695;

void plot_cluster_time(const std::string &configname = "config_showershape.yaml",
                       const std::string filetype = "data",
                       const bool use_leading_tower_time = false,
                       const bool apply_leading_time_corr = false,
                       const std::string &leading_time_corr_profile_file = "/sphenix/user/shuhangli/ppg12/plotting/figures/leading_tower_time_vs_adc_profile.root",
                       const std::string &leading_time_corr_profile_name = "prof_leading_time_vs_adc")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    bool issim = (filetype != "data");
    std::cout << "Cluster time definition: "
              << (use_leading_tower_time ? "LEADING TOWER time" : "ENERGY-WEIGHTED average tower time")
              << std::endl;
    if (apply_leading_time_corr && !use_leading_tower_time)
    {
        std::cout << "WARNING: apply_leading_time_corr=true but use_leading_tower_time=false; correction will be ignored" << std::endl;
    }
    if (use_leading_tower_time && apply_leading_time_corr && !leading_time_corr_profile_file.empty())
    {
        std::cout << "Leading tower time correction profile: " << leading_time_corr_profile_file
                  << " (object: " << leading_time_corr_profile_name << ")" << std::endl;
    }

    // Truth-level pT range cuts (to avoid double counting between samples)
    float max_photon_lower = 0;
    float max_photon_upper = 100;
    float max_jet_lower = 0;
    float max_jet_upper = 100;
    float weight = 1.0;
    float cross_weight = 1.0;
    float vertex_weight = 1.0;

    // Set ranges and weights based on file type
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
        weight = jet10cross / jet50cross;
    }
    else if (filetype == "jet15")
    {
        max_jet_lower = 15;
        max_jet_upper = 20;
        weight = jet15cross / jet50cross;
    }
    else if (filetype == "jet20")
    {
        max_jet_lower = 20;
        max_jet_upper = 30;
        weight = jet20cross / jet50cross;
    }
    else if (filetype == "jet30")
    {
        max_jet_lower = 30;
        max_jet_upper = 50;
        weight = jet30cross / jet50cross;
    }
    else if (filetype == "jet50")
    {
        max_jet_lower = 50;
        max_jet_upper = 100;
        weight = jet50cross / jet50cross;
    }

    cross_weight = weight;

    std::string infilename;
    if (!issim)
    {
        infilename = configYaml["input"]["data_file"].as<std::string>();
    }
    else
    {
        std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
        std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
        infilename = infilename_root_dir + filetype + infilename_branch_dir;
    }

    std::cout << "Input file: " << infilename << std::endl;

    // Build input chain
    std::string treename = configYaml["input"]["tree"].as<std::string>();
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());

    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();

    // Read analysis parameters from config
    float vertexcut = configYaml["analysis"]["vertex_cut"].as<float>();
    float reco_min_ET = configYaml["analysis"]["reco_min_ET"].as<float>();

    // Get pT bins
    std::vector<float> pT_bins = configYaml["analysis"]["pT_bins"].as<std::vector<float>>();
    int n_pT_bins = pT_bins.size() - 1;
    double pT_bin_edges[n_pT_bins + 1];
    std::copy(pT_bins.begin(), pT_bins.end(), pT_bin_edges);

    // Common selection parameters
    float common_prob_max = configYaml["analysis"]["common"]["prob_max"].as<float>();
    float common_prob_min = configYaml["analysis"]["common"]["prob_min"].as<float>();
    float common_e11_over_e33_max = configYaml["analysis"]["common"]["e11_over_e33_max"].as<float>();
    float common_e11_over_e33_min = configYaml["analysis"]["common"]["e11_over_e33_min"].as<float>();
    float common_cluster_weta_cogx_bound = configYaml["analysis"]["common"]["cluster_weta_cogx_bound"].as<float>();
    int common_b2bjet_cut = configYaml["analysis"]["common_b2bjet_cut"].as<int>(0);
    float common_b2bjet_pt_min = configYaml["analysis"]["common_b2bjet_pt_min"].as<float>(7.0);

    // Tight selection parameters
    float tight_weta_cogx_min = configYaml["analysis"]["tight"]["weta_cogx_min"].as<float>();
    float tight_weta_cogx_max_b = configYaml["analysis"]["tight"]["weta_cogx_max_b"].as<float>();
    float tight_weta_cogx_max_s = configYaml["analysis"]["tight"]["weta_cogx_max_s"].as<float>();
    float tight_wphi_cogx_min = configYaml["analysis"]["tight"]["wphi_cogx_min"].as<float>();
    float tight_wphi_cogx_max_b = configYaml["analysis"]["tight"]["wphi_cogx_max_b"].as<float>();
    float tight_wphi_cogx_max_s = configYaml["analysis"]["tight"]["wphi_cogx_max_s"].as<float>();
    float tight_et1_min_b = configYaml["analysis"]["tight"]["et1_min_b"].as<float>();
    float tight_et1_min_s = configYaml["analysis"]["tight"]["et1_min_s"].as<float>();
    float tight_et1_max = configYaml["analysis"]["tight"]["et1_max"].as<float>();
    float tight_et2_min = configYaml["analysis"]["tight"]["et2_min"].as<float>(0.0);
    float tight_et2_max = configYaml["analysis"]["tight"]["et2_max"].as<float>(1.0);
    float tight_et3_min = configYaml["analysis"]["tight"]["et3_min"].as<float>(0.0);
    float tight_et3_max = configYaml["analysis"]["tight"]["et3_max"].as<float>(1.0);
    float tight_et4_min = configYaml["analysis"]["tight"]["et4_min"].as<float>();
    float tight_et4_max = configYaml["analysis"]["tight"]["et4_max"].as<float>();
    float tight_e11_over_e33_min = configYaml["analysis"]["tight"]["e11_over_e33_min"].as<float>();
    float tight_e11_over_e33_max = configYaml["analysis"]["tight"]["e11_over_e33_max"].as<float>();
    float tight_e32_over_e35_min = configYaml["analysis"]["tight"]["e32_over_e35_min"].as<float>();
    float tight_e32_over_e35_max = configYaml["analysis"]["tight"]["e32_over_e35_max"].as<float>();
    float tight_prob_min = configYaml["analysis"]["tight"]["prob_min"].as<float>();
    float tight_prob_max = configYaml["analysis"]["tight"]["prob_max"].as<float>();
    float tight_bdt_min = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0);
    float tight_bdt_max = configYaml["analysis"]["tight"]["bdt_max"].as<float>(1);

    // Isolation parameters
    float recoiso_min = configYaml["analysis"]["reco_iso_min"].as<float>();
    float recoiso_max_b = configYaml["analysis"]["reco_iso_max_b"].as<float>();
    float recoiso_max_s = configYaml["analysis"]["reco_iso_max_s"].as<float>();
    float recononiso_min_shift = configYaml["analysis"]["reco_noniso_min_shift"].as<float>();
    float recononiso_max = configYaml["analysis"]["reco_noniso_max"].as<float>();

    int conesize = configYaml["analysis"]["cone_size"].as<int>();
    std::string bdt_model_name = configYaml["input"]["bdt_model_name"].as<std::string>("base");

    // Keep isoET definition consistent with RecoEffCalculator_TTreeReader.C
    const int iso_threshold = configYaml["analysis"]["iso_threshold"].as<int>(0);
    const int iso_hcalonly = configYaml["analysis"]["iso_hcalonly"].as<int>(0);
    // options: 0, 0.05, 0.1, 0.2
    const float iso_emcalinnerr = configYaml["analysis"]["iso_emcalinnerr"].as<float>(0.0);
    const float mc_iso_shift = configYaml["analysis"]["mc_iso_shift"].as<float>(0.0);
    const float mc_iso_scale = configYaml["analysis"]["mc_iso_scale"].as<float>(1.2);

    // Vertex reweighting for simulation (optional, consistent with RecoEffCalculator_TTreeReader.C):
    //   results/vertex_reweight_bdt_none.root : h_vertexz_ratio_data_over_mccombined
    TH1 *h_vertex_reweight = nullptr;
    int vertex_reweight_on = 1;
    std::string vertex_reweight_file = "results/vertex_reweight_bdt_none.root";
    if (issim)
    {
        vertex_reweight_on = configYaml["analysis"]["vertex_reweight_on"].as<int>(1);
        vertex_reweight_file =
            configYaml["analysis"]["vertex_reweight_file"].as<std::string>("results/vertex_reweight.root");

        if (vertex_reweight_on)
        {
            TFile *fvtx = TFile::Open(vertex_reweight_file.c_str(), "READ");
            if (!fvtx || fvtx->IsZombie())
            {
                std::cerr << "[VertexReweight] ERROR: cannot open vertex reweight file: "
                          << vertex_reweight_file << std::endl;
                return;
            }

            TH1 *htmp = dynamic_cast<TH1 *>(fvtx->Get("h_vertexz_ratio_data_over_mccombined"));
            if (!htmp)
            {
                std::cerr << "[VertexReweight] ERROR: cannot find histogram 'h_vertexz_ratio_data_over_mccombined' in "
                          << vertex_reweight_file << std::endl;
                fvtx->Close();
                delete fvtx;
                return;
            }

            h_vertex_reweight = dynamic_cast<TH1 *>(htmp->Clone("h_vertexz_ratio_data_over_mccombined_clone"));
            fvtx->Close();
            delete fvtx;

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

    // First, check which branches exist in ALL files in the chain
    std::cout << "\n========================================" << std::endl;
    std::cout << "Checking branch availability in all files..." << std::endl;
    std::cout << "========================================" << std::endl;

    TObjArray *fileElements = chain.GetListOfFiles();
    int nFiles = fileElements->GetEntries();
    std::cout << "Total files in chain: " << nFiles << std::endl;

    bool all_have_mbd_time = true;
    bool all_have_jet_time = true;

    for (int iFile = 0; iFile < nFiles; iFile++)
    {
        TChainElement *element = (TChainElement *)fileElements->At(iFile);
        std::string filename = element->GetTitle();
        TFile *f = TFile::Open(filename.c_str());
        if (!f || f->IsZombie())
        {
            std::cerr << "ERROR: Cannot open file: " << filename << std::endl;
            continue;
        }

        TTree *tree = (TTree *)f->Get(treename.c_str());
        if (!tree)
        {
            std::cerr << "ERROR: Cannot find tree " << treename << " in file: " << filename << std::endl;
            f->Close();
            continue;
        }

        bool has_mbd = (tree->GetBranch("mbd_time") != nullptr);
        bool has_jet = (tree->GetBranch("jet_time") != nullptr);

        if (!has_mbd || !has_jet)
        {
            std::cout << "File " << iFile << ": " << filename << std::endl;
            std::cout << "  mbd_time: " << (has_mbd ? "FOUND" : "NOT FOUND") << std::endl;
            std::cout << "  jet_time: " << (has_jet ? "FOUND" : "NOT FOUND") << std::endl;
        }

        all_have_mbd_time = all_have_mbd_time && has_mbd;
        all_have_jet_time = all_have_jet_time && has_jet;

        f->Close();
    }

    std::cout << "========================================" << std::endl;
    std::cout << "All files have mbd_time: " << (all_have_mbd_time ? "YES" : "NO") << std::endl;
    std::cout << "All files have jet_time: " << (all_have_jet_time ? "YES" : "NO") << std::endl;
    std::cout << "========================================\n" << std::endl;

    if (!all_have_mbd_time)
    {
        std::cerr << "ERROR: mbd_time branch not found in all files! Cannot proceed." << std::endl;
        return;
    }

    if (!all_have_jet_time)
    {
        std::cerr << "ERROR: jet_time branch not found in all files! Cannot proceed." << std::endl;
        return;
    }

    // TTreeReader setup
    TTreeReader reader(&chain);
    TTreeReaderValue<int> mbdnorthhit(reader, "mbdnorthhit");
    TTreeReaderValue<int> mbdsouthhit(reader, "mbdsouthhit");
    TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
    TTreeReaderValue<int> runnumber(reader, "runnumber");
    TTreeReaderValue<float> vertexz(reader, "vertexz");
    TTreeReaderValue<float> mbd_time(reader, "mbd_time");

    // Jet arrays including jet_time
    TTreeReaderValue<int> njet(reader, "njet");
    TTreeReaderArray<float> jet_Pt(reader, "jet_Pt");
    TTreeReaderArray<float> jet_Eta(reader, "jet_Eta");
    TTreeReaderArray<float> jet_Phi(reader, "jet_Phi");
    TTreeReaderArray<float> jet_time(reader, "jet_time");  // NEW BRANCH

    // Cluster arrays
    TTreeReaderArray<float> cluster_E(reader, Form("cluster_E_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_prob(reader, Form("cluster_prob_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e11(reader, Form("cluster_e11_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e33(reader, Form("cluster_e33_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e35(reader, Form("cluster_e35_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e32(reader, Form("cluster_e32_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_weta_cogx(reader, Form("cluster_weta_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_wphi_cogx(reader, Form("cluster_wphi_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et1(reader, Form("cluster_et1_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et2(reader, Form("cluster_et2_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et3(reader, Form("cluster_et3_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et4(reader, Form("cluster_et4_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_nsaturated(reader, Form("cluster_nsaturated_%s", clusternodename.c_str()));
    // Isolation branches (keep consistent with RecoEffCalculator_TTreeReader.C)
    TTreeReaderArray<float> cluster_iso_02(reader, Form("cluster_iso_02_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03(reader, Form("cluster_iso_03_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_04(reader, Form("cluster_iso_04_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_emcal(reader, Form("cluster_iso_03_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalin(reader, Form("cluster_iso_03_70_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalout(reader, Form("cluster_iso_03_70_hcalout_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_005_70_emcal(reader, Form("cluster_iso_005_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_01_70_emcal(reader, Form("cluster_iso_01_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_02_70_emcal(reader, Form("cluster_iso_02_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_bdt(reader, Form("cluster_bdt_%s_%s", clusternodename.c_str(), bdt_model_name.c_str()));

    // Cluster time arrays
    TTreeReaderArray<float> cluster_e_array(reader, Form("cluster_e_array_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_adc_array(reader, Form("cluster_adc_array_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_time_array(reader, Form("cluster_time_array_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_status_array(reader, Form("cluster_status_array_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ownership_array(reader, Form("cluster_ownership_array_%s", clusternodename.c_str()));

    // Truth particle information (only for simulation)
    TTreeReaderValue<int> *nparticle = nullptr;
    TTreeReaderArray<float> *particle_Pt = nullptr;
    TTreeReaderArray<float> *particle_Eta = nullptr;
    TTreeReaderArray<int> *particle_pid = nullptr;
    TTreeReaderValue<int> *njet_truth = nullptr;
    TTreeReaderArray<float> *jet_truth_Pt = nullptr;

    if (issim)
    {
        nparticle = new TTreeReaderValue<int>(reader, "nparticles");
        particle_Pt = new TTreeReaderArray<float>(reader, "particle_Pt");
        particle_Eta = new TTreeReaderArray<float>(reader, "particle_Eta");
        particle_pid = new TTreeReaderArray<int>(reader, "particle_pid");
        njet_truth = new TTreeReaderValue<int>(reader, "njet_truth");
        jet_truth_Pt = new TTreeReaderArray<float>(reader, "jet_truth_Pt");
    }

    // Load MBD t0 correction
    std::map<int, float> mbd_t0_correction;
    std::ifstream mbd_file("/sphenix/user/shuhangli/ppg12/efficiencytool/MbdOut.corr");
    std::string mbd_line;
    while (std::getline(mbd_file, mbd_line))
    {
        std::istringstream iss(mbd_line);
        int runnumber_corr;
        float t0;
        iss >> runnumber_corr >> t0;
        mbd_t0_correction[runnumber_corr] = t0;
    }

    // Create output file
    std::string suffix = "";
    if (use_leading_tower_time)
    {
        suffix = apply_leading_time_corr ? "_leadingTowerTimeCorr" : "_leadingTowerTime";
    }
    std::string outfilename = "results/cluster_time_analysis_" + filetype + suffix + ".root";
    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");

    // Optional: load leading-tower ADC slewing correction profile (produced by plotting/plot_tower_timing.C)
    TProfile *prof_leading_time_vs_adc = nullptr;
    double prof_leading_time_vs_adc_global_mean = 0.0;
    if (use_leading_tower_time && apply_leading_time_corr && !leading_time_corr_profile_file.empty())
    {
        prof_leading_time_vs_adc = load_profile_clone(leading_time_corr_profile_file, leading_time_corr_profile_name);
        if (prof_leading_time_vs_adc)
        {
            // axis=2 is mean of Y-values (time)
            prof_leading_time_vs_adc_global_mean = prof_leading_time_vs_adc->GetMean(2);
            std::cout << "Loaded leading-tower time correction profile (global mean time = "
                      << prof_leading_time_vs_adc_global_mean << " ns)" << std::endl;
        }
        else
        {
            std::cout << "WARNING: Leading-tower time correction requested but profile could not be loaded; continuing without correction"
                      << std::endl;
        }
    }

    // Create inclusive histograms across all pTs for jet-time analysis
    // For each selection: all (no selection), common, tight_iso, tight_noniso, NPB
    TH2D *h_all_delta_t_jet_all = new TH2D("h_all_delta_t_jet_all",
                                            "All Clusters: Cluster-Jet Time vs Cluster pT;Cluster pT (GeV);Cluster-Jet Time (ns)",
                                            n_pT_bins, pT_bin_edges, 160, -40, 40);
    TH2D *h_common_delta_t_jet_all = new TH2D("h_common_delta_t_jet_all",
                                               "Common Selection: Cluster-Jet Time vs Cluster pT;Cluster pT (GeV);Cluster-Jet Time (ns)",
                                               n_pT_bins, pT_bin_edges, 160, -40, 40);
    TH2D *h_tight_iso_delta_t_jet_all = new TH2D("h_tight_iso_delta_t_jet_all",
                                                  "Tight+Iso Selection: Cluster-Jet Time vs Cluster pT;Cluster pT (GeV);Cluster-Jet Time (ns)",
                                                  n_pT_bins, pT_bin_edges, 160, -40, 40);
    TH2D *h_tight_noniso_delta_t_jet_all = new TH2D("h_tight_noniso_delta_t_jet_all",
                                                     "Tight+NonIso Selection: Cluster-Jet Time vs Cluster pT;Cluster pT (GeV);Cluster-Jet Time (ns)",
                                                     n_pT_bins, pT_bin_edges, 160, -40, 40);
    TH2D *h_tight_delta_t_jet_all = new TH2D("h_tight_delta_t_jet_all",
                                              "Tight Selection: Cluster-Jet Time vs Cluster pT;Cluster pT (GeV);Cluster-Jet Time (ns)",
                                              n_pT_bins, pT_bin_edges, 160, -40, 40);
    TH2D *h_nontight_delta_t_jet_all = new TH2D("h_nontight_delta_t_jet_all",
                                                 "Non-Tight Selection: Cluster-Jet Time vs Cluster pT;Cluster pT (GeV);Cluster-Jet Time (ns)",
                                                 n_pT_bins, pT_bin_edges, 160, -40, 40);
    TH2D *h_npb_delta_t_jet_all = new TH2D("h_npb_delta_t_jet_all",
                                            "NPB Selection: Cluster-Jet Time vs Cluster pT;Cluster pT (GeV);Cluster-Jet Time (ns)",
                                            n_pT_bins, pT_bin_edges, 160, -40, 40);

    // MBD-time eta-binned histograms (pT>10 GeV requirement)
    // For 4 selections: all (no selection), common, weta>1.0, NPB

    // MBD-time vs cluster pT (ET) inclusive histograms
    TH2D *h_all_delta_t_mbd_all = new TH2D("h_all_delta_t_mbd_all",
                                          "All Clusters: Cluster-MBD Time vs Cluster pT;Cluster pT (GeV);Cluster-MBD Time (ns)",
                                          n_pT_bins, pT_bin_edges, 160, -40, 40);

    // MBD-time vs cluster pT (ET) inclusive histogram for NPB selection (pT>10 GeV requirement, see fill section)
    TH2D *h_npb_delta_t_mbd_all = new TH2D("h_npb_delta_t_mbd_all",
                                          "NPB Selection: Cluster-MBD Time vs Cluster pT;Cluster pT (GeV);Cluster-MBD Time (ns)",
                                          n_pT_bins, pT_bin_edges, 160, -40, 40);

    // Cluster time vs cluster pT (ET) inclusive histogram (pT>10 GeV requirement, see fill section)
    TH2D *h_all_cluster_t_vs_pt_all = new TH2D("h_all_cluster_t_vs_pt_all",
                                              "All Clusters: Cluster Time vs Cluster pT;Cluster pT (GeV);Cluster Time (ns)",
                                              n_pT_bins, pT_bin_edges, 160, -40, 40);

    // All clusters MBD histograms
    TH2D *h_all_delta_t_mbd_vs_eta = new TH2D("h_all_delta_t_mbd_vs_eta",
                                               "Cluster-MBD Time vs Eta (All Clusters, pT>10 GeV);Cluster #eta;Cluster-MBD Time (ns)",
                                               50, -1.0, 1.0, 160, -40, 40);
    TH2D *h_all_cluster_t_vs_eta = new TH2D("h_all_cluster_t_vs_eta",
                                             "Cluster Time vs Eta (All Clusters, pT>10 GeV);Cluster #eta;Cluster Time (ns)",
                                             50, -1.0, 1.0, 160, -40, 40);

    // Common selection MBD histograms
    TH2D *h_common_delta_t_mbd_vs_eta = new TH2D("h_common_delta_t_mbd_vs_eta",
                                                  "Cluster-MBD Time vs Eta (Common Selection, pT>10 GeV);Cluster #eta;Cluster-MBD Time (ns)",
                                                  50, -1.0, 1.0, 160, -40, 40);
    TH2D *h_common_cluster_t_vs_eta = new TH2D("h_common_cluster_t_vs_eta",
                                                "Cluster Time vs Eta (Common Selection, pT>10 GeV);Cluster #eta;Cluster Time (ns)",
                                                50, -1.0, 1.0, 160, -40, 40);

    // weta>1.0 selection MBD histograms
    TH2D *h_weta1p0_delta_t_mbd_vs_eta = new TH2D("h_weta1p0_delta_t_mbd_vs_eta",
                                                   "Cluster-MBD Time vs Eta (weta_cogx>1.0, pT>10 GeV);Cluster #eta;Cluster-MBD Time (ns)",
                                                   50, -1.0, 1.0, 160, -40, 40);
    TH2D *h_weta1p0_cluster_t_vs_eta = new TH2D("h_weta1p0_cluster_t_vs_eta",
                                                 "Cluster Time vs Eta (weta_cogx>1.0, pT>10 GeV);Cluster #eta;Cluster Time (ns)",
                                                 50, -1.0, 1.0, 160, -40, 40);

    // NPB selection MBD histograms
    TH2D *h_npb_delta_t_mbd_vs_eta = new TH2D("h_npb_delta_t_mbd_vs_eta",
                                               "Cluster-MBD Time vs Eta (NPB Selection, pT>10 GeV);Cluster #eta;Cluster-MBD Time (ns)",
                                               50, -1.0, 1.0, 160, -40, 40);
    TH2D *h_npb_cluster_t_vs_eta = new TH2D("h_npb_cluster_t_vs_eta",
                                             "Cluster Time vs Eta (NPB Selection, pT>10 GeV);Cluster #eta;Cluster Time (ns)",
                                             50, -1.0, 1.0, 160, -40, 40);

    // Skip bad runs
    std::set<int> skiprunnumbers = {};

    int nentries = chain.GetEntries();
    int ientry = 0;

    std::cout << "Processing " << nentries << " entries..." << std::endl;

    while (reader.Next())
    {
        // Event-by-event MC vertex reweighting.
        // Update `weight` so all existing Fill(..., weight) calls use the per-event weight.
        weight = cross_weight;
        vertex_weight = 1.0;
        if (issim)
        {
            if (vertex_reweight_on)
            {
                if (!h_vertex_reweight)
                {
                    std::cerr << "[VertexReweight] ERROR: vertex reweighting is enabled but histogram is not loaded."
                              << std::endl;
                    return;
                }
                int bin = h_vertex_reweight->FindBin(*vertexz);
                if (bin < 1)
                    bin = 1;
                if (bin > h_vertex_reweight->GetNbinsX())
                    bin = h_vertex_reweight->GetNbinsX();
                vertex_weight = h_vertex_reweight->GetBinContent(bin);
            }

            if (!std::isfinite(vertex_weight) || vertex_weight <= 0.0)
            {
                std::cout << "Warning: vertex weight is nan or inf" << std::endl;
                std::cout << "vertexz: " << *vertexz << std::endl;
                vertex_weight = 1.0;
            }

            weight *= vertex_weight;
        }

        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;

        // Skip bad runs for data
        if (!issim && skiprunnumbers.find(*runnumber) != skiprunnumbers.end())
        {
            ientry++;
            continue;
        }

        // Vertex cut
        if (std::abs(*vertexz) > vertexcut)
        {
            ientry++;
            continue;
        }

        // MBD hit requirement
        if (!(*mbdnorthhit >= 1 && *mbdsouthhit >= 1))
        {
            ientry++;
            continue;
        }

        // Calculate MBD mean time with correction
        float mbd_mean_time = *mbd_time;
        float mbdoffset = 0;
        if (mbd_t0_correction.find(*runnumber) != mbd_t0_correction.end())
        {
            mbdoffset = mbd_t0_correction[*runnumber];
        }
        if (issim)
        {
            mbdoffset = 0.0;
        }
        mbd_mean_time = mbd_mean_time - mbdoffset;

        // Apply truth-level pT range cuts for simulation
        if (issim)
        {
            // Find max truth photon pT
            float maxphotonpT = -1;
            for (int iparticle = 0; iparticle < **nparticle; iparticle++)
            {
                int pid = (*particle_pid)[iparticle];
                if (pid == 22 && std::abs((*particle_Eta)[iparticle]) < 0.7)
                {
                    if ((*particle_Pt)[iparticle] > maxphotonpT)
                    {
                        maxphotonpT = (*particle_Pt)[iparticle];
                    }
                }
            }

            // Find max truth jet pT
            float maxjetpT = -1;
            for (int ijet = 0; ijet < **njet_truth; ijet++)
            {
                if ((*jet_truth_Pt)[ijet] > maxjetpT)
                {
                    maxjetpT = (*jet_truth_Pt)[ijet];
                }
            }

            // Check if event passes truth pT range cuts
            bool passes_photon_range = (maxphotonpT >= max_photon_lower && maxphotonpT < max_photon_upper);
            bool passes_jet_range = (maxjetpT >= max_jet_lower && maxjetpT < max_jet_upper);

            // For photon samples, require photon in range
            // For jet samples, require jet in range
            if (filetype.find("photon") != std::string::npos)
            {
                if (!passes_photon_range)
                {
                    ientry++;
                    continue;
                }
            }
            else if (filetype.find("jet") != std::string::npos)
            {
                if (!passes_jet_range)
                {
                    ientry++;
                    continue;
                }
            }
        }

        // Loop over clusters
        for (int icluster = 0; icluster < *ncluster; icluster++)
        {
            // Minimum ET cut
            if (cluster_Et[icluster] < reco_min_ET)
                continue;

            // Calculate cluster time:
            // - default: energy-weighted average tower time over towers in cluster (excluding zero-suppressed)
            // - optional: leading tower time (tower with max energy in cluster, excluding zero-suppressed)
            float clustertime_ns = 0.0;
            if (!use_leading_tower_time)
            {
                float clusteravgtime = 0;
                float cluster_total_e = 0;
                for (int i = 0; i < 49; i++)
                {
                    const int idx = icluster * 49 + i;
                    if (cluster_ownership_array[idx] == 1)
                    {
                        const int status = cluster_status_array[idx];
                        if (status & (1 << 5))
                        {
                            continue;
                        }
                        clusteravgtime += cluster_time_array[idx] * cluster_e_array[idx];
                        cluster_total_e += cluster_e_array[idx];
                    }
                }
                clusteravgtime = cluster_total_e > 0 ? clusteravgtime / cluster_total_e : 0;
                clustertime_ns = clusteravgtime * TIME_SAMPLE_NS;
            }
            else
            {
                float max_tower_e = -1.0;
                float leading_time_samples = 0.0;
                float leading_adc = 0.0;
                for (int i = 0; i < 49; i++)
                {
                    const int idx = icluster * 49 + i;
                    if (cluster_ownership_array[idx] != 1)
                        continue;

                    const int status = cluster_status_array[idx];
                    if (status & (1 << 5))
                        continue;

                    const float e = cluster_e_array[idx];
                    if (e > max_tower_e)
                    {
                        max_tower_e = e;
                        leading_time_samples = cluster_time_array[idx];
                        leading_adc = cluster_adc_array[idx];
                    }
                }
                const double raw_time_ns = (max_tower_e > 0) ? (leading_time_samples * TIME_SAMPLE_NS) : 0.0;
                if (apply_leading_time_corr)
                {
                    clustertime_ns = correct_leading_time_ns(raw_time_ns,
                                                            leading_adc,
                                                            prof_leading_time_vs_adc,
                                                            prof_leading_time_vs_adc_global_mean);
                }
                else
                {
                    clustertime_ns = raw_time_ns;
                }
            }

            // Calculate delta_t between cluster and MBD time
            float delta_t_mbd = clustertime_ns - mbd_mean_time;

            // Find leading b2b jet for cluster-jet timing
            // Only use jets with |dphi| > 3pi/4 and pT > threshold
            float b2b_jet_time = 0;
            float b2b_jet_pt = -1;
            float b2bjet_dphi_cut = 3 * M_PI / 4;
            float b2bjet_pt_threshold = 5.0; // GeV
            bool has_b2b_jet = false;

            for (int ijet = 0; ijet < *njet; ijet++)
            {
                float dphi = cluster_Phi[icluster] - jet_Phi[ijet];
                while (dphi > M_PI)
                    dphi = dphi - 2 * M_PI;
                while (dphi < -M_PI)
                    dphi = dphi + 2 * M_PI;

                // Check if jet is in b2b region and above threshold
                if (std::abs(dphi) > b2bjet_dphi_cut && jet_Pt[ijet] > b2bjet_pt_threshold)
                {
                    if (jet_Pt[ijet] > b2b_jet_pt)
                    {
                        b2b_jet_pt = jet_Pt[ijet];
                        b2b_jet_time = jet_time[ijet];
                        has_b2b_jet = true;
                    }
                }
            }

            // Calculate delta_t between cluster and b2b jet (only if b2b jet exists)
            float delta_t_jet = 0;
            if (has_b2b_jet)
            {
                delta_t_jet = clustertime_ns - b2b_jet_time;
            }

            float cluster_eta = cluster_Eta[icluster];
            float clusterET = cluster_Et[icluster];

            // Calculate e11_over_e33 and e32_over_e35
            float e11_over_e33 = cluster_e11[icluster] / cluster_e33[icluster];
            float e32_over_e35 = cluster_e32[icluster] / cluster_e35[icluster];

            // Calculate isolation
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

            // Optional "threshold/70%" iso definition used in RecoEffCalculator_TTreeReader.C
            if (iso_threshold)
            {
                if (iso_hcalonly)
                {
                    recoisoET = cluster_iso_03_70_hcalin[icluster] + cluster_iso_03_70_hcalout[icluster];
                }
                else
                {
                    recoisoET = cluster_iso_03_70_emcal[icluster] + cluster_iso_03_70_hcalin[icluster] + cluster_iso_03_70_hcalout[icluster];
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

            // Optional MC isoET fudge (consistent with RecoEffCalculator_TTreeReader.C)
            if (issim)
            {
                recoisoET = recoisoET * mc_iso_scale;
                recoisoET += mc_iso_shift;
            }

            // Check for b2b jet and otherside_jet
            bool otherside_jet = false;
            float max_b2bjet_pT = -1;
            float b2bjet_dphi = 3 * M_PI / 4;
            float jet_eta_cut = 0.6;

            for (int ijet = 0; ijet < *njet; ijet++)
            {
                float dphi = cluster_Phi[icluster] - jet_Phi[ijet];
                while (dphi > M_PI)
                    dphi = dphi - 2 * M_PI;
                while (dphi < -M_PI)
                    dphi = dphi + 2 * M_PI;

                if (std::abs(dphi) > (M_PI / 2))
                {
                    otherside_jet = true;
                }

                if (std::abs(jet_Eta[ijet]) < jet_eta_cut)
                {
                    if (std::abs(dphi) > b2bjet_dphi)
                    {
                        if (jet_Pt[ijet] > max_b2bjet_pT && jet_Pt[ijet] > 5)
                        {
                            max_b2bjet_pT = jet_Pt[ijet];
                        }
                    }
                }
            }

            bool passes_common_b2bjet = true;
            if (common_b2bjet_cut)
            {
                passes_common_b2bjet = (max_b2bjet_pT >= common_b2bjet_pt_min);
            }

            // Common selection
            bool passes_common_shape =
                cluster_prob[icluster] > common_prob_min &&
                cluster_prob[icluster] < common_prob_max &&
                e11_over_e33 > common_e11_over_e33_min &&
                e11_over_e33 < common_e11_over_e33_max &&
                (cluster_weta_cogx[icluster] < common_cluster_weta_cogx_bound);

            bool common_pass = passes_common_shape && passes_common_b2bjet;

            // Tight selection
            float tight_weta_cogx_max = tight_weta_cogx_max_b + tight_weta_cogx_max_s * clusterET;
            float tight_wphi_cogx_max = tight_wphi_cogx_max_b + tight_wphi_cogx_max_s * clusterET;
            float tight_et1_min = tight_et1_min_b + tight_et1_min_s * clusterET;

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

            bool tight = false;
            if (common_pass &&
                is_cluster_weta_cogx_tight &&
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

            // Isolation cuts
            float recoiso_max = recoiso_max_b + recoiso_max_s * clusterET;
            float recononiso_min = recoiso_max + recononiso_min_shift;

            bool iso = false;
            bool noniso = false;

            if (recoisoET > recoiso_min && recoisoET < recoiso_max)
            {
                iso = true;
            }
            if (recoisoET > recononiso_min && recoisoET < recononiso_max)
            {
                noniso = true;
            }

            // NPB selection
            // For cluster-jet timing: only use weta cut
            // For cluster-MBD timing: use original NPB selection (!otherside_jet && weta >= 0.6)
            float npb_weta_min = 0.6;
            bool isnpb_jet = (cluster_weta_cogx[icluster] >= npb_weta_min);
            bool isnpb_mbd = !otherside_jet && (cluster_weta_cogx[icluster] >= npb_weta_min);

            // Fill cluster-jet histograms (only if has valid b2b jet)
            if (has_b2b_jet)
            {
                // All clusters (no selection)
                h_all_delta_t_jet_all->Fill(clusterET, delta_t_jet, weight);

                if (common_pass)
                {
                    h_common_delta_t_jet_all->Fill(clusterET, delta_t_jet, weight);
                }

                if (tight && iso)
                {
                    h_tight_iso_delta_t_jet_all->Fill(clusterET, delta_t_jet, weight);
                }

                if (tight && noniso)
                {
                    h_tight_noniso_delta_t_jet_all->Fill(clusterET, delta_t_jet, weight);
                }

                // NPB for cluster-jet: only weta cut
                if (isnpb_jet)
                {
                    h_npb_delta_t_jet_all->Fill(clusterET, delta_t_jet, weight);
                }
            }

            // Fill eta-binned MBD histograms (only for pT>10)
            if (clusterET > 10.0)
            {
                // Inclusive cluster pT vs MBD delta-t (all clusters)
                h_all_delta_t_mbd_all->Fill(clusterET, delta_t_mbd, weight);
                h_all_cluster_t_vs_pt_all->Fill(clusterET, clustertime_ns, weight);

                // Inclusive cluster pT vs MBD delta-t (NPB selection)
                if (isnpb_mbd)
                {
                    h_npb_delta_t_mbd_all->Fill(clusterET, delta_t_mbd, weight);
                }

                // All clusters (no selection)
                h_all_delta_t_mbd_vs_eta->Fill(cluster_eta, delta_t_mbd, weight);
                h_all_cluster_t_vs_eta->Fill(cluster_eta, clustertime_ns, weight);

                // Common selection
                if (common_pass)
                {
                    h_common_delta_t_mbd_vs_eta->Fill(cluster_eta, delta_t_mbd, weight);
                    h_common_cluster_t_vs_eta->Fill(cluster_eta, clustertime_ns, weight);
                }

                // weta>1.0 selection
                if (cluster_weta_cogx[icluster] > 1.0)
                {
                    h_weta1p0_delta_t_mbd_vs_eta->Fill(cluster_eta, delta_t_mbd, weight);
                    h_weta1p0_cluster_t_vs_eta->Fill(cluster_eta, clustertime_ns, weight);
                }

                // NPB selection for MBD: use original NPB (!otherside_jet && weta >= 0.6)
                if (isnpb_mbd)
                {
                    h_npb_delta_t_mbd_vs_eta->Fill(cluster_eta, delta_t_mbd, weight);
                    h_npb_cluster_t_vs_eta->Fill(cluster_eta, clustertime_ns, weight);
                }
            }
        }

        ientry++;
    }

    // Write histograms
    fout->cd();

    // Jet-time inclusive histograms
    h_all_delta_t_jet_all->Write();
    h_common_delta_t_jet_all->Write();
    h_tight_iso_delta_t_jet_all->Write();
    h_tight_noniso_delta_t_jet_all->Write();
    h_npb_delta_t_jet_all->Write();

    // MBD-time vs cluster pT (inclusive)
    h_all_delta_t_mbd_all->Write();
    h_all_cluster_t_vs_pt_all->Write();
    h_npb_delta_t_mbd_all->Write();

    // MBD-time eta-binned histograms
    h_all_delta_t_mbd_vs_eta->Write();
    h_all_cluster_t_vs_eta->Write();
    h_common_delta_t_mbd_vs_eta->Write();
    h_common_cluster_t_vs_eta->Write();
    h_weta1p0_delta_t_mbd_vs_eta->Write();
    h_weta1p0_cluster_t_vs_eta->Write();
    h_npb_delta_t_mbd_vs_eta->Write();
    h_npb_cluster_t_vs_eta->Write();

    fout->Close();
    delete fout;

    // Cleanup dynamically allocated truth readers
    if (issim)
    {
        delete nparticle;
        delete particle_Pt;
        delete particle_Eta;
        delete particle_pid;
        delete njet_truth;
        delete jet_truth_Pt;
    }

    delete prof_leading_time_vs_adc;
    delete h_vertex_reweight;

    std::cout << "========================================" << std::endl;
    std::cout << "Output written to: " << outfilename << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Jet-time inclusive histograms (all pTs):" << std::endl;
    std::cout << "  - 5 selections (all, common, tight+iso, tight+noniso, NPB)" << std::endl;
    std::cout << "  - Total: 5 histograms" << std::endl;
    std::cout << "MBD-time eta-binned histograms (pT>10):" << std::endl;
    std::cout << "  - 4 selections (all, common, weta>1.0, NPB) × 2 types (delta_t, cluster_t)" << std::endl;
    std::cout << "  - Total: 8 histograms" << std::endl;
    std::cout << "MBD-time vs cluster pT (pT>10):" << std::endl;
    std::cout << "  - 2 histograms: h_all_delta_t_mbd_all, h_npb_delta_t_mbd_all" << std::endl;
    std::cout << "Cluster time vs cluster pT (pT>10):" << std::endl;
    std::cout << "  - 1 histogram: h_all_cluster_t_vs_pt_all" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Grand total: 16 histograms" << std::endl;
    std::cout << "========================================" << std::endl;
}
