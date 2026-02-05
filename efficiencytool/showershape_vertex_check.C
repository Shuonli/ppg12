#include <yaml-cpp/yaml.h>
#include <TMVA/RBDT.hxx>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TRandom3.h>
#include <cmath>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <iostream>

// NPB phase-space (matches config_npb_training.yaml)
const float NPB_ET_MIN  = 6.0;
const float NPB_ET_MAX  = 40.0;
const float NPB_ETA_MAX = 0.7;

void showershape_vertex_check(const std::string &configname = "config_showershape.yaml",
                              const std::string filetype = "jet10",
                              bool doinclusive = true)
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    bool issim = true;
    bool isbackground = false;

    if (filetype == "data")
    {
        std::cerr << "ERROR: this macro is MC-only. Use a simulation filetype." << std::endl;
        return;
    }

    std::string inclusive_str = doinclusive ? "_inclusive" : "";

    std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
    std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;

    std::cout << "infilename: " << infilename << std::endl;

    // ---------------------------------------------------------------
    // Cross-section weights (same as ShowerShapeCheck.C)
    // ---------------------------------------------------------------
    float max_photon_lower = 0;
    float max_photon_upper = 100;

    const float photon5cross  = 146359.3;
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
    // Vertex reweighting histogram (MC event weight)
    // ---------------------------------------------------------------
    TH1 *h_vertex_reweight = nullptr;
    int vertex_reweight_on = configYaml["analysis"]["vertex_reweight_on"].as<int>(1);
    std::string vertex_reweight_file =
        configYaml["analysis"]["vertex_reweight_file"].as<std::string>("results/vertex_reweight.root");

    std::cout << "loading vertex reweight file: " << vertex_reweight_file << std::endl;

    TFile *fvtx = TFile::Open(vertex_reweight_file.c_str(), "READ");
    if (!fvtx || fvtx->IsZombie())
    {
        std::cerr << "[VertexReweight] ERROR: cannot open vertex reweight file: "
                  << vertex_reweight_file << std::endl;
        return;
    }

    if (vertex_reweight_on)
    {
        TH1 *htmp = dynamic_cast<TH1 *>(fvtx->Get("h_vertexz_ratio_data_over_mccombined"));
        if (!htmp)
        {
            std::cerr << "[VertexReweight] ERROR: cannot find histogram "
                         "'h_vertexz_ratio_data_over_mccombined'" << std::endl;
            fvtx->Close();
            return;
        }
        h_vertex_reweight = dynamic_cast<TH1 *>(htmp->Clone("h_vertexz_ratio_clone"));
        h_vertex_reweight->SetDirectory(nullptr);
        h_vertex_reweight->Sumw2();
        std::cout << "[VertexReweight] Loaded reweight histogram from "
                  << vertex_reweight_file << std::endl;
    }

    // ---------------------------------------------------------------
    // Data vertex distribution for random sampling
    // ---------------------------------------------------------------
    TH1D *h_vertexz_data = nullptr;
    {
        TH1D *htmp_data = dynamic_cast<TH1D *>(fvtx->Get("h_vertexz_data"));
        if (!htmp_data)
        {
            std::cerr << "[VertexData] ERROR: cannot find histogram 'h_vertexz_data' in "
                      << vertex_reweight_file << std::endl;
            fvtx->Close();
            return;
        }
        h_vertexz_data = dynamic_cast<TH1D *>(htmp_data->Clone("h_vertexz_data_clone"));
        h_vertexz_data->SetDirectory(nullptr);
        std::cout << "[VertexData] Loaded data vertex distribution ("
                  << h_vertexz_data->GetEntries() << " entries)" << std::endl;
    }

    TRandom3 rng(42);

    // ---------------------------------------------------------------
    // NPB TMVA model
    // ---------------------------------------------------------------
    std::string npb_tmva_file = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/npb_models/npb_score_tmva.root";
    if (configYaml["analysis"] && configYaml["analysis"]["npb_tmva_file"])
    {
        npb_tmva_file = configYaml["analysis"]["npb_tmva_file"].as<std::string>();
    }
    if (gSystem->AccessPathName(npb_tmva_file.c_str()))
    {
        std::cerr << "ERROR: NPB TMVA model not found at '" << npb_tmva_file
                  << "'. Cannot proceed." << std::endl;
        return;
    }
    std::cout << "[NPB] Loading TMVA model from " << npb_tmva_file << std::endl;
    TMVA::Experimental::RBDT npb_bdt("myBDT", npb_tmva_file);

    // ---------------------------------------------------------------
    // Showershape BDT TMVA model (for re-scoring with shifted vertex)
    // ---------------------------------------------------------------
    std::string ss_bdt_model_name = "base_v3E";  // default, will be updated from config
    if (configYaml["input"] && configYaml["input"]["bdt_model_name"])
    {
        ss_bdt_model_name = configYaml["input"]["bdt_model_name"].as<std::string>();
    }
    std::string ss_bdt_file = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_"
                              + ss_bdt_model_name + "_single_tmva.root";
    if (gSystem->AccessPathName(ss_bdt_file.c_str()))
    {
        std::cerr << "ERROR: Showershape BDT model not found at '" << ss_bdt_file
                  << "'. Cannot proceed." << std::endl;
        return;
    }
    std::cout << "[SS-BDT] Loading showershape BDT model (" << ss_bdt_model_name
              << ") from " << ss_bdt_file << std::endl;
    TMVA::Experimental::RBDT ss_bdt("myBDT", ss_bdt_file);

    // ---------------------------------------------------------------
    // Config parsing (simplified)
    // ---------------------------------------------------------------
    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();
    std::string bdt_model_name = configYaml["input"]["bdt_model_name"].as<std::string>("base");
    std::string bdt_branch_name = Form("cluster_bdt_%s_%s", clusternodename.c_str(), bdt_model_name.c_str());
    std::cout << "[BDT] Using showershape BDT branch: " << bdt_branch_name << std::endl;

    float vertexcut = configYaml["analysis"]["vertex_cut"].as<float>();
    std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();
    std::vector<float> pT_bins = configYaml["analysis"]["pT_bins"].as<std::vector<float>>();
    int n_pT_bins = pT_bins.size() - 1;
    int conesize = configYaml["analysis"]["cone_size"].as<int>();
    float reco_min_ET = configYaml["analysis"]["reco_min_ET"].as<float>();
    float truthisocut = configYaml["analysis"]["truth_iso_max"].as<float>();

    // ---------------------------------------------------------------
    // TChain + TTreeReader
    // ---------------------------------------------------------------
    std::string treename = configYaml["input"]["tree"].as<std::string>();
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());

    TTreeReader reader(&chain);

    // Event-level
    TTreeReaderValue<int> mbdnorthhit(reader, "mbdnorthhit");
    TTreeReaderValue<int> mbdsouthhit(reader, "mbdsouthhit");
    TTreeReaderValue<int> nparticles(reader, "nparticles");
    TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
    TTreeReaderValue<float> vertexz(reader, "vertexz");

    // Truth particles
    TTreeReaderArray<float> particle_Pt(reader, "particle_Pt");
    TTreeReaderArray<float> particle_Eta(reader, "particle_Eta");
    TTreeReaderArray<float> particle_truth_iso_02(reader, "particle_truth_iso_02");
    TTreeReaderArray<float> particle_truth_iso_03(reader, "particle_truth_iso_03");
    TTreeReaderArray<float> particle_truth_iso_04(reader, "particle_truth_iso_04");
    TTreeReaderArray<int> particle_pid(reader, "particle_pid");
    TTreeReaderArray<int> particle_trkid(reader, "particle_trkid");
    TTreeReaderArray<int> particle_photonclass(reader, "particle_photonclass");

    // Truth jets
    TTreeReaderValue<int> njet_truth(reader, "njet_truth");
    TTreeReaderArray<float> jet_truth_Pt(reader, "jet_truth_Pt");

    // Cluster arrays (25 NPB features + truth matching)
    TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e11(reader, Form("cluster_e11_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e22(reader, Form("cluster_e22_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e13(reader, Form("cluster_e13_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e15(reader, Form("cluster_e15_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e17(reader, Form("cluster_e17_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e31(reader, Form("cluster_e31_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e51(reader, Form("cluster_e51_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e71(reader, Form("cluster_e71_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e33(reader, Form("cluster_e33_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e32(reader, Form("cluster_e32_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e35(reader, Form("cluster_e35_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e37(reader, Form("cluster_e37_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e53(reader, Form("cluster_e53_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_weta_cogx(reader, Form("cluster_weta_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_wphi_cogx(reader, Form("cluster_wphi_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et1(reader, Form("cluster_et1_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et2(reader, Form("cluster_et2_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et3(reader, Form("cluster_et3_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et4(reader, Form("cluster_et4_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w32(reader, Form("cluster_w32_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w52(reader, Form("cluster_w52_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w72(reader, Form("cluster_w72_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusternodename.c_str()));

    // Showershape BDT score (from tree, already computed)
    TTreeReaderArray<float> cluster_bdt(reader, bdt_branch_name.c_str());

    // ---------------------------------------------------------------
    // Output file
    // ---------------------------------------------------------------
    std::string outfilename = configYaml["output"]["eff_outfile"].as<std::string>()
                              + "_vertex_check_" + filetype + inclusive_str + ".root";
    std::cout << "outfilename: " << outfilename << std::endl;

    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    fout->cd();

    // ---------------------------------------------------------------
    // Histogram booking
    // ---------------------------------------------------------------
    int n_eta_bins = (int)eta_bins.size() - 1;

    // Per [eta][pT]: new NPB score vs old NPB score
    std::vector<std::vector<TH2D *>> h_npb_2d(n_eta_bins);
    // Per [eta][pT]: new BDT score vs old BDT score
    std::vector<std::vector<TH2D *>> h_bdt_2d(n_eta_bins);

    // pT-inclusive (per eta)
    std::vector<TH2D *> h_npb_2d_inc(n_eta_bins);
    std::vector<TH2D *> h_bdt_2d_inc(n_eta_bins);

    for (int ieta = 0; ieta < n_eta_bins; ieta++)
    {
        h_npb_2d_inc[ieta] = new TH2D(
            Form("h_npb_2d_eta%d_inclusive", ieta),
            Form("NPB Score: Original vs Shifted %.1f < #eta < %.1f;Original NPB Score;Shifted NPB Score",
                 eta_bins[ieta], eta_bins[ieta + 1]),
            100, 0.0, 1.0, 100, 0.0, 1.0);

        h_bdt_2d_inc[ieta] = new TH2D(
            Form("h_bdt_2d_eta%d_inclusive", ieta),
            Form("BDT Score (%s): Original vs Shifted %.1f < #eta < %.1f;Original BDT Score;Shifted BDT Score",
                 bdt_model_name.c_str(), eta_bins[ieta], eta_bins[ieta + 1]),
            100, 0.0, 1.0, 100, 0.0, 1.0);

        for (int ipt = 0; ipt < n_pT_bins; ipt++)
        {
            h_npb_2d[ieta].push_back(new TH2D(
                Form("h_npb_2d_eta%d_pt%d", ieta, ipt),
                Form("NPB Score: Original vs Shifted %.1f < #eta < %.1f, %.1f < p_{T} < %.1f GeV;Original NPB Score;Shifted NPB Score",
                     eta_bins[ieta], eta_bins[ieta + 1], pT_bins[ipt], pT_bins[ipt + 1]),
                100, 0.0, 1.0, 100, 0.0, 1.0));

            h_bdt_2d[ieta].push_back(new TH2D(
                Form("h_bdt_2d_eta%d_pt%d", ieta, ipt),
                Form("BDT Score (%s): Original vs Shifted %.1f < #eta < %.1f, %.1f < p_{T} < %.1f GeV;Original BDT Score;Shifted BDT Score",
                     bdt_model_name.c_str(), eta_bins[ieta], eta_bins[ieta + 1], pT_bins[ipt], pT_bins[ipt + 1]),
                100, 0.0, 1.0, 100, 0.0, 1.0));
        }
    }

    // ---------------------------------------------------------------
    // NPB feature vector builder
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
            cluster_Et[icl],          // index 0
            cluster_Eta[icl],         // index 1
            vtxz,                     // index 2 -- modified for double-interaction
            e11_over_e33,             // index 3
            e32_over_e35,             // index 4
            e11_over_e22,             // index 5
            e11_over_e13,             // index 6
            e11_over_e15,             // index 7
            e11_over_e17,             // index 8
            e11_over_e31,             // index 9
            e11_over_e51,             // index 10
            e11_over_e71,             // index 11
            e22_over_e33,             // index 12
            e22_over_e35,             // index 13
            e22_over_e37,             // index 14
            e22_over_e53,             // index 15
            cluster_weta_cogx[icl],   // index 16
            cluster_wphi_cogx[icl],   // index 17
            cluster_et1[icl],         // index 18
            cluster_et2[icl],         // index 19
            cluster_et3[icl],         // index 20
            cluster_et4[icl],         // index 21
            cluster_w32[icl],         // index 22
            cluster_w52[icl],         // index 23
            cluster_w72[icl],         // index 24
        };
    };

    // ---------------------------------------------------------------
    // Showershape BDT feature vector builder
    // Feature order depends on model variant (matches apply_BDT.C)
    // ---------------------------------------------------------------
    auto buildSsBdtFeatureVector = [&](int icl, float vtxz) -> std::vector<float> {
        const float e11_over_e33 = (cluster_e33[icl] > 0) ? cluster_e11[icl] / cluster_e33[icl] : 0.0f;
        const float e32_over_e35 = (cluster_e35[icl] > 0) ? cluster_e32[icl] / cluster_e35[icl] : 0.0f;

        // Feature order for base_v3E: ET, weta, wphi, vertexz, eta, e11/e33, et1-4, e32/e35
        if (ss_bdt_model_name == "base_v3E")
        {
            return {
                cluster_Et[icl],
                cluster_weta_cogx[icl],
                cluster_wphi_cogx[icl],
                vtxz,
                cluster_Eta[icl],
                e11_over_e33,
                cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl],
                e32_over_e35
            };
        }
        // Feature order for base_v3: weta, wphi, vertexz, eta, e11/e33, et1-4, e32/e35
        else if (ss_bdt_model_name == "base_v3")
        {
            return {
                cluster_weta_cogx[icl],
                cluster_wphi_cogx[icl],
                vtxz,
                cluster_Eta[icl],
                e11_over_e33,
                cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl],
                e32_over_e35
            };
        }
        // Feature order for base_v2E: ET, weta, wphi, vertexz, eta, e11/e33, et1-4
        else if (ss_bdt_model_name == "base_v2E")
        {
            return {
                cluster_Et[icl],
                cluster_weta_cogx[icl],
                cluster_wphi_cogx[icl],
                vtxz,
                cluster_Eta[icl],
                e11_over_e33,
                cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl]
            };
        }
        // Feature order for base_v2: weta, wphi, vertexz, eta, e11/e33, et1-4
        else if (ss_bdt_model_name == "base_v2")
        {
            return {
                cluster_weta_cogx[icl],
                cluster_wphi_cogx[icl],
                vtxz,
                cluster_Eta[icl],
                e11_over_e33,
                cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl]
            };
        }
        // Feature order for base_v1E: ET, weta, vertexz, eta, e11/e33, et1-4
        else if (ss_bdt_model_name == "base_v1E")
        {
            return {
                cluster_Et[icl],
                cluster_weta_cogx[icl],
                vtxz,
                cluster_Eta[icl],
                e11_over_e33,
                cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl]
            };
        }
        // Feature order for base_v1: weta, vertexz, eta, e11/e33, et1-4
        else if (ss_bdt_model_name == "base_v1")
        {
            return {
                cluster_weta_cogx[icl],
                vtxz,
                cluster_Eta[icl],
                e11_over_e33,
                cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl]
            };
        }
        // Feature order for base_E: ET, vertexz, eta, e11/e33, et1-4
        else if (ss_bdt_model_name == "base_E")
        {
            return {
                cluster_Et[icl],
                vtxz,
                cluster_Eta[icl],
                e11_over_e33,
                cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl]
            };
        }
        // Feature order for base: vertexz, eta, e11/e33, et1-4
        else if (ss_bdt_model_name == "base" || ss_bdt_model_name == "base_vr")
        {
            return {
                vtxz,
                cluster_Eta[icl],
                e11_over_e33,
                cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl]
            };
        }
        // Feature order for base_v0E: ET, vertexz, eta, e11/e33, et2-4
        else if (ss_bdt_model_name == "base_v0E")
        {
            return {
                cluster_Et[icl],
                vtxz,
                cluster_Eta[icl],
                e11_over_e33,
                cluster_et2[icl], cluster_et3[icl], cluster_et4[icl]
            };
        }
        // Feature order for base_v0: vertexz, eta, e11/e33, et2-4
        else if (ss_bdt_model_name == "base_v0")
        {
            return {
                vtxz,
                cluster_Eta[icl],
                e11_over_e33,
                cluster_et2[icl], cluster_et3[icl], cluster_et4[icl]
            };
        }
        // Default fallback to base_v3E
        else
        {
            std::cerr << "WARNING: Unknown BDT model '" << ss_bdt_model_name
                      << "', using base_v3E feature order" << std::endl;
            return {
                cluster_Et[icl],
                cluster_weta_cogx[icl],
                cluster_wphi_cogx[icl],
                vtxz,
                cluster_Eta[icl],
                e11_over_e33,
                cluster_et1[icl], cluster_et2[icl], cluster_et3[icl], cluster_et4[icl],
                e32_over_e35
            };
        }
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

        // Cross-section weight + vertex reweighting
        weight = cross_weight;
        vertex_weight = 1.0;
        if (vertex_reweight_on && h_vertex_reweight)
        {
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

        // MC truth: build signal set and particle track ID map
        std::set<int> signal_set;
        std::set<int> background_set;
        std::map<int, int> particle_trkidmap;

        float maxphotonpT = 0;
        for (int iparticle = 0; iparticle < *nparticles; iparticle++)
        {
            particle_trkidmap[particle_trkid[iparticle]] = iparticle;

            float truthisoET = 0;
            if (conesize == 4)
                truthisoET = particle_truth_iso_04[iparticle];
            else if (conesize == 3)
                truthisoET = particle_truth_iso_03[iparticle];
            else if (conesize == 2)
                truthisoET = particle_truth_iso_02[iparticle];

            if (particle_pid[iparticle] == 22)
            {
                if (particle_Pt[iparticle] > maxphotonpT)
                    maxphotonpT = particle_Pt[iparticle];

                if (particle_photonclass[iparticle] < 3) // direct or fragmentation
                {
                    if (truthisoET < truthisocut)
                        signal_set.insert(iparticle);
                }
                else
                {
                    background_set.insert(iparticle);
                }
            }
            else
            {
                background_set.insert(iparticle);
            }
        }

        // pT-hat range cut for photon samples
        if (!isbackground)
        {
            if (maxphotonpT > max_photon_upper || maxphotonpT < max_photon_lower)
            {
                ientry++;
                continue;
            }
        }

        // pT-hat range cut for jet samples
        if (isbackground)
        {
            float maxjetpT = 0;
            for (int ijet = 0; ijet < *njet_truth; ijet++)
            {
                if (jet_truth_Pt[ijet] > maxjetpT)
                    maxjetpT = jet_truth_Pt[ijet];
            }
            if (maxjetpT > max_jet_upper || maxjetpT < max_jet_lower)
            {
                ientry++;
                continue;
            }
        }

        // Draw one random vertex per event from data distribution
        float random_vertex = h_vertexz_data->GetRandom();
        float double_vtx = (*vertexz + random_vertex) / 2.0;

        // Vertex cut on original vertex
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

        // Cluster loop
        for (int icl = 0; icl < *ncluster; icl++)
        {
            // ET cut
            if (cluster_Et[icl] < reco_min_ET)
                continue;

            // Background ET ceiling to reduce fluctuations
            if (isbackground && cluster_Et[icl] > cluster_ET_upper)
                continue;

            // Truth matching
            if (particle_trkidmap.find(cluster_truthtrkID[icl]) == particle_trkidmap.end())
                continue;

            int iparticle = particle_trkidmap[cluster_truthtrkID[icl]];
            if (!isbackground)
            {
                if (signal_set.find(iparticle) == signal_set.end())
                    continue;
            }
            else
            {
                if (!doinclusive)
                {
                    if (background_set.find(iparticle) == background_set.end())
                        continue;
                }
            }

            // NPB phase-space cut
            if (cluster_Et[icl] < NPB_ET_MIN || cluster_Et[icl] > NPB_ET_MAX)
                continue;
            if (std::fabs(cluster_Eta[icl]) > NPB_ETA_MAX)
                continue;

            // Find eta bin
            int etabin = -1;
            for (int ieta = 0; ieta < n_eta_bins; ieta++)
            {
                if (cluster_Eta[icl] > eta_bins[ieta] && cluster_Eta[icl] < eta_bins[ieta + 1])
                {
                    etabin = ieta;
                    break;
                }
            }
            if (etabin == -1) continue;

            // Find pT bin
            int ptbin = -1;
            for (int ipt = 0; ipt < n_pT_bins; ipt++)
            {
                if (cluster_Et[icl] > pT_bins[ipt] && cluster_Et[icl] < pT_bins[ipt + 1])
                {
                    ptbin = ipt;
                    break;
                }
            }

            // Compute original NPB score
            std::vector<float> x_npb_original = buildNpbFeatureVector(icl, *vertexz);
            float npb_score_original = npb_bdt.Compute(x_npb_original)[0];

            // Compute shifted NPB score with double-interaction vertex
            std::vector<float> x_npb_shifted = buildNpbFeatureVector(icl, double_vtx);
            float npb_score_shifted = npb_bdt.Compute(x_npb_shifted)[0];

            // Compute original showershape BDT score
            std::vector<float> x_bdt_original = buildSsBdtFeatureVector(icl, *vertexz);
            float bdt_score_original = ss_bdt.Compute(x_bdt_original)[0];

            // Compute shifted showershape BDT score with double-interaction vertex
            std::vector<float> x_bdt_shifted = buildSsBdtFeatureVector(icl, double_vtx);
            float bdt_score_shifted = ss_bdt.Compute(x_bdt_shifted)[0];

            // Fill pT-inclusive histograms
            h_npb_2d_inc[etabin]->Fill(npb_score_original, npb_score_shifted, weight);
            h_bdt_2d_inc[etabin]->Fill(bdt_score_original, bdt_score_shifted, weight);

            // Fill per-pT-bin histograms
            if (ptbin != -1)
            {
                h_npb_2d[etabin][ptbin]->Fill(npb_score_original, npb_score_shifted, weight);
                h_bdt_2d[etabin][ptbin]->Fill(bdt_score_original, bdt_score_shifted, weight);
            }
        }
        ientry++;
    }

    fout->cd();
    fout->Write();
    fout->Close();
    std::cout << "Done. Output written to " << outfilename << std::endl;
}
