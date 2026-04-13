// apply_BDT.C
#include <TMVA/RBDT.hxx> // fast interpreter
#include <iostream>
#include <yaml-cpp/yaml.h>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <utility>

// ---------------------------------------------------------------------
// Container size for cluster-level arrays used by this macro.
// NOTE: This intentionally matches the pre-refactor value. CaloAna24.h
// uses nclustermax = 10000 for the producer, but this macro historically
// used 4096 for the reader side (slimtree ncluster is typically << 100).
// ---------------------------------------------------------------------
static const int nclustercontainermx = 4096;

// Model names must match trained file stems: model_{name}{suffix}_single_tmva.root
static const std::vector<std::string> kModelNames = {
    "base",
    "base_vr",
    "base_v0",
    "base_v1",
    "base_v2",
    "base_v3",
    "base_E",
    "base_v0E",
    "base_v1E",
    "base_v2E",
    "base_v3E",
};

// ---------------------------------------------------------------------
// VariantState
//
// Per-cluster-node state: input branch backing arrays, RBDT model lists,
// and output branch backing arrays. Allocated on the heap (large struct;
// avoid stack overflow when multiple variants are active).
// ---------------------------------------------------------------------
struct VariantState
{
    std::string node;         // e.g. "CLUSTERINFO_CEMC" or "CLUSTERINFO_CEMC_NO_SPLIT"
    std::string model_suffix; // e.g. "_split" or "_nosplit" or "" (legacy)

    // --- input branch arrays (owned, heap) ---
    int   ncluster = 0;

    std::vector<float> cluster_Et;
    std::vector<float> cluster_Eta;
    std::vector<float> cluster_Phi;
    std::vector<float> cluster_et1;
    std::vector<float> cluster_et2;
    std::vector<float> cluster_et3;
    std::vector<float> cluster_et4;
    std::vector<float> cluster_weta_cogx;
    std::vector<float> cluster_wphi_cogx;
    std::vector<float> cluster_prob;
    std::vector<float> cluster_e11;
    std::vector<float> cluster_e22;
    std::vector<float> cluster_e13;
    std::vector<float> cluster_e15;
    std::vector<float> cluster_e17;
    std::vector<float> cluster_e31;
    std::vector<float> cluster_e51;
    std::vector<float> cluster_e71;
    std::vector<float> cluster_e33;
    std::vector<float> cluster_e35;
    std::vector<float> cluster_e37;
    std::vector<float> cluster_e53;
    std::vector<float> cluster_e32;
    std::vector<float> cluster_w32;
    std::vector<float> cluster_w52;
    std::vector<float> cluster_w72;
    std::vector<float> cluster_iso_02;
    std::vector<float> cluster_iso_03;
    std::vector<float> cluster_iso_04;
    std::vector<float> cluster_iso_03_60_emcal;
    std::vector<float> cluster_iso_03_60_hcalin;
    std::vector<float> cluster_iso_03_60_hcalout;

    // --- BDT models ---
    std::vector<TMVA::Experimental::RBDT> bdt_list;
    TMVA::Experimental::RBDT *npb_bdt = nullptr; // owned
    bool do_npb_score = false;

    // --- output branch arrays (one per photon-ID model + one NPB) ---
    // cluster_bdt[imodel][icluster] — contiguous, branch-bound
    std::vector<std::vector<float>> cluster_bdt;
    std::vector<float> cluster_npb_score;

    VariantState()
        : cluster_Et(nclustercontainermx, 0.f),
          cluster_Eta(nclustercontainermx, 0.f),
          cluster_Phi(nclustercontainermx, 0.f),
          cluster_et1(nclustercontainermx, 0.f),
          cluster_et2(nclustercontainermx, 0.f),
          cluster_et3(nclustercontainermx, 0.f),
          cluster_et4(nclustercontainermx, 0.f),
          cluster_weta_cogx(nclustercontainermx, 0.f),
          cluster_wphi_cogx(nclustercontainermx, 0.f),
          cluster_prob(nclustercontainermx, 0.f),
          cluster_e11(nclustercontainermx, 0.f),
          cluster_e22(nclustercontainermx, 0.f),
          cluster_e13(nclustercontainermx, 0.f),
          cluster_e15(nclustercontainermx, 0.f),
          cluster_e17(nclustercontainermx, 0.f),
          cluster_e31(nclustercontainermx, 0.f),
          cluster_e51(nclustercontainermx, 0.f),
          cluster_e71(nclustercontainermx, 0.f),
          cluster_e33(nclustercontainermx, 0.f),
          cluster_e35(nclustercontainermx, 0.f),
          cluster_e37(nclustercontainermx, 0.f),
          cluster_e53(nclustercontainermx, 0.f),
          cluster_e32(nclustercontainermx, 0.f),
          cluster_w32(nclustercontainermx, 0.f),
          cluster_w52(nclustercontainermx, 0.f),
          cluster_w72(nclustercontainermx, 0.f),
          cluster_iso_02(nclustercontainermx, 0.f),
          cluster_iso_03(nclustercontainermx, 0.f),
          cluster_iso_04(nclustercontainermx, 0.f),
          cluster_iso_03_60_emcal(nclustercontainermx, 0.f),
          cluster_iso_03_60_hcalin(nclustercontainermx, 0.f),
          cluster_iso_03_60_hcalout(nclustercontainermx, 0.f),
          cluster_npb_score(nclustercontainermx, 0.f)
    {
        cluster_bdt.assign(kModelNames.size(), std::vector<float>(nclustercontainermx, 0.f));
    }

    ~VariantState()
    {
        delete npb_bdt;
        npb_bdt = nullptr;
    }
};

// ---------------------------------------------------------------------
// parse_variants
//
// Returns list of (cluster_node_name, model_suffix) pairs.
//
// If `input.cluster_variants` is present in the YAML, use it verbatim.
// Otherwise synthesize a single-element list from the legacy fields
// `input.cluster_node_name` + `analysis.use_split_bdt` so existing
// configs keep working unchanged.
// ---------------------------------------------------------------------
static std::vector<std::pair<std::string, std::string>>
parse_variants(const YAML::Node &cfg)
{
    std::vector<std::pair<std::string, std::string>> out;
    if (cfg["input"] && cfg["input"]["cluster_variants"] &&
        cfg["input"]["cluster_variants"].IsSequence() &&
        cfg["input"]["cluster_variants"].size() > 0)
    {
        for (const auto &v : cfg["input"]["cluster_variants"])
        {
            std::string node = v["node"].as<std::string>();
            std::string suffix = v["model_suffix"] ? v["model_suffix"].as<std::string>() : std::string("");
            out.emplace_back(node, suffix);
        }
    }
    else
    {
        std::string node = cfg["input"]["cluster_node_name"].as<std::string>();
        const bool use_split = cfg["analysis"] && cfg["analysis"]["use_split_bdt"] &&
                               cfg["analysis"]["use_split_bdt"].as<int>(0) != 0;
        std::string suffix = use_split ? "_split" : "";
        out.emplace_back(node, suffix);
    }
    return out;
}

void apply_BDT(const std::string &configname = "config_nom.yaml", const std::string filetype = "jet12_double", const std::string inputfilename = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/auau_test/caloana.root")
{
    using namespace TMVA::Experimental;
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    bool issim = true;
    if (filetype == "data")
    {
        issim = false;
    }

    std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
    std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;
    if (!issim)
    {
        infilename = inputfilename;
    }

    // recoisoET definition (match `BDTinput.C`)
    const int conesize = configYaml["analysis"]["cone_size"].as<int>(3);
    const int iso_threshold = configYaml["analysis"]["iso_threshold"].as<int>(0);

    // NPB training phase-space (match `config_npb_training.yaml`)
    const float npb_et_min = 6.0;
    const float npb_et_max = 40.0;
    const float npb_eta_max = 0.7;

    // ------------------------------------------------------------
    // Build the list of cluster-node variants to score.
    // ------------------------------------------------------------
    std::vector<std::pair<std::string, std::string>> variant_specs = parse_variants(configYaml);
    std::cout << "apply_BDT: " << variant_specs.size() << " variant(s):" << std::endl;
    for (const auto &vs : variant_specs)
    {
        std::cout << "  - node='" << vs.first << "' model_suffix='" << vs.second << "'" << std::endl;
    }

    // ------------------------------------------------------------
    // Optional per-variant NPB TMVA override; otherwise derive from suffix.
    // The legacy single-path override `analysis.npb_tmva_file` is still
    // honored when only one variant is active (back-compat).
    // ------------------------------------------------------------
    std::string legacy_npb_override;
    if (configYaml["analysis"] && configYaml["analysis"]["npb_tmva_file"])
    {
        legacy_npb_override = configYaml["analysis"]["npb_tmva_file"].as<std::string>();
    }

    // Open input + output files
    TFile *ftreein = new TFile(infilename.c_str(), "READ");
    TTree *slimtree = (TTree *)ftreein->Get(configYaml["input"]["tree"].as<std::string>().c_str());

    // TODO: filename is kept as "bdt_split.root" for backward compatibility
    // with downstream configs (efficiencytool, plotting). After refactor it
    // contains BOTH split and no-split BDT branches. Rename once downstream
    // configs migrate.
    std::string outfile_name = filetype + "/bdt_split.root";
    if (!issim)
    {
        // TODO: filename kept as "_with_bdt_split.root" for backward compat.
        std::string namebase = inputfilename.substr(0, inputfilename.find_last_of("."));
        outfile_name = namebase + "_with_bdt_split.root";
    }
    TFile *fout = new TFile(outfile_name.c_str(), "RECREATE");

    // Clone full structure, zero entries so far
    TTree *outtree = slimtree->CloneTree(0);

    // ------------------------------------------------------------
    // Common (non-per-variant) branches
    // ------------------------------------------------------------
    float vertexz = 0;
    slimtree->SetBranchStatus("vertexz", 1);
    slimtree->SetBranchAddress("vertexz", &vertexz);

    // ------------------------------------------------------------
    // Allocate VariantState per variant and wire up branches
    // ------------------------------------------------------------
    std::vector<std::unique_ptr<VariantState>> variants;
    variants.reserve(variant_specs.size());

    for (size_t iv = 0; iv < variant_specs.size(); ++iv)
    {
        auto vs = std::unique_ptr<VariantState>(new VariantState());
        vs->node = variant_specs[iv].first;
        vs->model_suffix = variant_specs[iv].second;

        // --- Load photon-ID RBDT models ---
        vs->bdt_list.clear();
        vs->bdt_list.reserve(kModelNames.size());
        std::cout << "[variant " << iv << " / node=" << vs->node
                  << "] loading photon-ID models (suffix='" << vs->model_suffix << "'):" << std::endl;
        for (size_t im = 0; im < kModelNames.size(); ++im)
        {
            std::string mf = "binned_models/model_" + kModelNames[im] + vs->model_suffix + "_single_tmva.root";
            std::cout << "    " << mf << std::endl;
            vs->bdt_list.push_back(TMVA::Experimental::RBDT("myBDT", mf));
        }

        // --- Load NPB TMVA model ---
        // Prefer: legacy single-variant override -> per-suffix default.
        std::string npb_tmva_file;
        if (variant_specs.size() == 1 && !legacy_npb_override.empty())
        {
            npb_tmva_file = legacy_npb_override;
        }
        else if (vs->model_suffix == "_split")
        {
            npb_tmva_file = "npb_models/npb_score_split_tmva.root";
        }
        else if (vs->model_suffix == "_nosplit")
        {
            npb_tmva_file = "npb_models/npb_score_nosplit_tmva.root";
        }
        else
        {
            // legacy "" suffix
            npb_tmva_file = "npb_models/npb_score_tmva.root";
        }

        vs->do_npb_score = !gSystem->AccessPathName(npb_tmva_file.c_str());
        if (!vs->do_npb_score)
        {
            std::cout << "WARNING: NPB TMVA model not found at '" << npb_tmva_file
                      << "' for variant " << vs->node
                      << ". Skipping NPB score application." << std::endl;
        }
        else
        {
            std::cout << "    NPB: " << npb_tmva_file << std::endl;
            vs->npb_bdt = new TMVA::Experimental::RBDT("myBDT", npb_tmva_file);
        }

        // --- Wire input branches for this variant ---
        const std::string &nn = vs->node;
        auto bind = [&](const char *prefix, void *addr) {
            std::string bn = std::string(prefix) + "_" + nn;
            slimtree->SetBranchStatus(bn.c_str(), 1);
            slimtree->SetBranchAddress(bn.c_str(), addr);
        };

        {
            std::string bn = std::string("ncluster_") + nn;
            slimtree->SetBranchStatus(bn.c_str(), 1);
            slimtree->SetBranchAddress(bn.c_str(), &vs->ncluster);
        }
        bind("cluster_Et",        vs->cluster_Et.data());
        bind("cluster_Eta",       vs->cluster_Eta.data());
        bind("cluster_Phi",       vs->cluster_Phi.data());
        bind("cluster_et1",       vs->cluster_et1.data());
        bind("cluster_et2",       vs->cluster_et2.data());
        bind("cluster_et3",       vs->cluster_et3.data());
        bind("cluster_et4",       vs->cluster_et4.data());
        bind("cluster_weta_cogx", vs->cluster_weta_cogx.data());
        bind("cluster_wphi_cogx", vs->cluster_wphi_cogx.data());
        bind("cluster_prob",      vs->cluster_prob.data());
        bind("cluster_e11",       vs->cluster_e11.data());
        bind("cluster_e22",       vs->cluster_e22.data());
        bind("cluster_e13",       vs->cluster_e13.data());
        bind("cluster_e15",       vs->cluster_e15.data());
        bind("cluster_e17",       vs->cluster_e17.data());
        bind("cluster_e31",       vs->cluster_e31.data());
        bind("cluster_e51",       vs->cluster_e51.data());
        bind("cluster_e71",       vs->cluster_e71.data());
        bind("cluster_e33",       vs->cluster_e33.data());
        bind("cluster_e32",       vs->cluster_e32.data());
        bind("cluster_e35",       vs->cluster_e35.data());
        bind("cluster_e37",       vs->cluster_e37.data());
        bind("cluster_e53",       vs->cluster_e53.data());
        bind("cluster_w32",       vs->cluster_w32.data());
        bind("cluster_w52",       vs->cluster_w52.data());
        bind("cluster_w72",       vs->cluster_w72.data());
        bind("cluster_iso_02",    vs->cluster_iso_02.data());
        bind("cluster_iso_03",    vs->cluster_iso_03.data());
        bind("cluster_iso_04",    vs->cluster_iso_04.data());
        bind("cluster_iso_03_60_emcal",  vs->cluster_iso_03_60_emcal.data());
        bind("cluster_iso_03_60_hcalin", vs->cluster_iso_03_60_hcalin.data());
        bind("cluster_iso_03_60_hcalout",vs->cluster_iso_03_60_hcalout.data());

        // --- Create output branches for this variant ---
        for (size_t im = 0; im < kModelNames.size(); ++im)
        {
            std::string bname = Form("cluster_bdt_%s_%s", nn.c_str(), kModelNames[im].c_str());
            std::string leaf  = Form("%s[ncluster_%s]/F", bname.c_str(), nn.c_str());
            outtree->Branch(bname.c_str(), vs->cluster_bdt[im].data(), leaf.c_str());
        }
        {
            std::string bname = Form("cluster_npb_score_%s", nn.c_str());
            std::string leaf  = Form("%s[ncluster_%s]/F", bname.c_str(), nn.c_str());
            outtree->Branch(bname.c_str(), vs->cluster_npb_score.data(), leaf.c_str());
        }

        variants.push_back(std::move(vs));
    }

    // ------------------------------------------------------------
    // Event loop — single pass over slimtree entries.
    // ------------------------------------------------------------
    int nentries = slimtree->GetEntries();
    for (int ientry = 0; ientry < nentries; ientry++)
    {
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        slimtree->GetEntry(ientry);

        // Score each variant in turn from its own input arrays
        for (size_t iv = 0; iv < variants.size(); ++iv)
        {
            VariantState &var = *variants[iv];

            for (int icluster = 0; icluster < var.ncluster; icluster++)
            {
                const float cluster_Et_BDT         = var.cluster_Et[icluster];
                const float cluster_Eta_BDT        = var.cluster_Eta[icluster];
                const float vertexz_BDT            = vertexz;
                const float e11_over_e33_BDT      = (var.cluster_e11[icluster] > 0) ? var.cluster_e11[icluster] / var.cluster_e33[icluster] : 0;
                const float e32_over_e35_BDT      = (var.cluster_e32[icluster] > 0) ? var.cluster_e32[icluster] / var.cluster_e35[icluster] : 0;
                const float cluster_weta_cogx_BDT = var.cluster_weta_cogx[icluster];
                const float cluster_wphi_cogx_BDT = var.cluster_wphi_cogx[icluster];
                const float cluster_et1_BDT       = var.cluster_et1[icluster];
                const float cluster_et2_BDT       = var.cluster_et2[icluster];
                const float cluster_et3_BDT       = var.cluster_et3[icluster];
                const float cluster_et4_BDT       = var.cluster_et4[icluster];

                // default: invalid (outside training phase-space or model missing)
                var.cluster_npb_score[icluster] = -1;

                // Photon-ID feature vectors — identical to pre-refactor version.
                std::vector<std::vector<float>> x_list = {
                    //base
                    {
                        vertexz_BDT,
                        cluster_Eta_BDT,
                        e11_over_e33_BDT,
                        cluster_et1_BDT, cluster_et2_BDT, cluster_et3_BDT, cluster_et4_BDT
                    },

                    //base_vr
                    {
                        vertexz_BDT,
                        cluster_Eta_BDT,
                        e11_over_e33_BDT,
                        cluster_et1_BDT, cluster_et2_BDT, cluster_et3_BDT, cluster_et4_BDT,
                    },

                    //base_v0
                    {
                        vertexz_BDT,
                        cluster_Eta_BDT,
                        e11_over_e33_BDT,
                        cluster_et2_BDT, cluster_et3_BDT, cluster_et4_BDT,
                    },

                    //base_v1:
                    {
                        cluster_weta_cogx_BDT,
                        vertexz_BDT,
                        cluster_Eta_BDT,
                        e11_over_e33_BDT,
                        cluster_et1_BDT, cluster_et2_BDT, cluster_et3_BDT, cluster_et4_BDT
                    },

                    //base_v2:
                    {
                        cluster_weta_cogx_BDT,
                        cluster_wphi_cogx_BDT,
                        vertexz_BDT,
                        cluster_Eta_BDT,
                        e11_over_e33_BDT,
                        cluster_et1_BDT, cluster_et2_BDT, cluster_et3_BDT, cluster_et4_BDT,
                    },

                    //base_v3
                    {
                        cluster_weta_cogx_BDT,
                        cluster_wphi_cogx_BDT,
                        vertexz_BDT,
                        cluster_Eta_BDT,
                        e11_over_e33_BDT,
                        cluster_et1_BDT, cluster_et2_BDT, cluster_et3_BDT, cluster_et4_BDT,
                        e32_over_e35_BDT
                    },

                    //base_E
                    {
                        cluster_Et_BDT,
                        vertexz_BDT,
                        cluster_Eta_BDT,
                        e11_over_e33_BDT,
                        cluster_et1_BDT, cluster_et2_BDT, cluster_et3_BDT, cluster_et4_BDT
                    },
                    //base_v0E
                    {
                        cluster_Et_BDT,
                        vertexz_BDT,
                        cluster_Eta_BDT,
                        e11_over_e33_BDT,
                        cluster_et2_BDT, cluster_et3_BDT, cluster_et4_BDT,
                    },

                    //base_v1E
                    {
                        cluster_Et_BDT,
                        cluster_weta_cogx_BDT,
                        vertexz_BDT,
                        cluster_Eta_BDT,
                        e11_over_e33_BDT,
                        cluster_et1_BDT, cluster_et2_BDT, cluster_et3_BDT, cluster_et4_BDT
                    },

                    //base_v2E
                    {
                        cluster_Et_BDT,
                        cluster_weta_cogx_BDT,
                        cluster_wphi_cogx_BDT,
                        vertexz_BDT,
                        cluster_Eta_BDT,
                        e11_over_e33_BDT,
                        cluster_et1_BDT, cluster_et2_BDT, cluster_et3_BDT, cluster_et4_BDT
                    },

                    //base_v3E
                    {
                        cluster_Et_BDT,
                        cluster_weta_cogx_BDT,
                        cluster_wphi_cogx_BDT,
                        vertexz_BDT,
                        cluster_Eta_BDT,
                        e11_over_e33_BDT,
                        cluster_et1_BDT, cluster_et2_BDT, cluster_et3_BDT, cluster_et4_BDT,
                        e32_over_e35_BDT
                    },
                };

                for (size_t i = 0; i < kModelNames.size(); i++)
                {
                    if (cluster_Et_BDT > 7)
                    {
                        var.cluster_bdt[i][icluster] = var.bdt_list[i].Compute(x_list[i])[0];
                    }
                    else
                    {
                        var.cluster_bdt[i][icluster] = -1;
                    }

                    // Debug output for first few clusters
                    if (icluster < 3 && ientry < 3) {
                        std::cout << "Entry " << ientry << ", Variant " << var.node
                                 << ", Cluster " << icluster
                                 << ", ET=" << cluster_Et_BDT
                                 << ", Model " << kModelNames[i]
                                 << ": score=" << var.cluster_bdt[i][icluster] << std::endl;
                    }
                }

                // NPB score application (single TMVA model)
                if (var.do_npb_score && var.npb_bdt &&
                    cluster_Et_BDT >= npb_et_min && cluster_Et_BDT <= npb_et_max &&
                    std::fabs(cluster_Eta_BDT) <= npb_eta_max)
                {
                    // recoisoET (match `BDTinput.C`)
                    float recoisoET = -999;
                    if (conesize == 4)
                        recoisoET = var.cluster_iso_04[icluster];
                    else if (conesize == 3)
                        recoisoET = var.cluster_iso_03[icluster];
                    else if (conesize == 2)
                        recoisoET = var.cluster_iso_02[icluster];

                    if (iso_threshold)
                    {
                        recoisoET = var.cluster_iso_03_60_emcal[icluster] +
                                    var.cluster_iso_03_60_hcalin[icluster] +
                                    var.cluster_iso_03_60_hcalout[icluster];
                    }

                    // ratios (protect against divide-by-zero)
                    const float e11_over_e33 = (var.cluster_e33[icluster] > 0) ? var.cluster_e11[icluster] / var.cluster_e33[icluster] : 0.0f;
                    const float e32_over_e35 = (var.cluster_e35[icluster] > 0) ? var.cluster_e32[icluster] / var.cluster_e35[icluster] : 0.0f;
                    const float e11_over_e22 = (var.cluster_e22[icluster] > 0) ? var.cluster_e11[icluster] / var.cluster_e22[icluster] : 0.0f;
                    const float e11_over_e13 = (var.cluster_e13[icluster] > 0) ? var.cluster_e11[icluster] / var.cluster_e13[icluster] : 0.0f;
                    const float e11_over_e15 = (var.cluster_e15[icluster] > 0) ? var.cluster_e11[icluster] / var.cluster_e15[icluster] : 0.0f;
                    const float e11_over_e17 = (var.cluster_e17[icluster] > 0) ? var.cluster_e11[icluster] / var.cluster_e17[icluster] : 0.0f;
                    const float e11_over_e31 = (var.cluster_e31[icluster] > 0) ? var.cluster_e11[icluster] / var.cluster_e31[icluster] : 0.0f;
                    const float e11_over_e51 = (var.cluster_e51[icluster] > 0) ? var.cluster_e11[icluster] / var.cluster_e51[icluster] : 0.0f;
                    const float e11_over_e71 = (var.cluster_e71[icluster] > 0) ? var.cluster_e11[icluster] / var.cluster_e71[icluster] : 0.0f;
                    const float e22_over_e33 = (var.cluster_e33[icluster] > 0) ? var.cluster_e22[icluster] / var.cluster_e33[icluster] : 0.0f;
                    const float e22_over_e35 = (var.cluster_e35[icluster] > 0) ? var.cluster_e22[icluster] / var.cluster_e35[icluster] : 0.0f;
                    const float e22_over_e37 = (var.cluster_e37[icluster] > 0) ? var.cluster_e22[icluster] / var.cluster_e37[icluster] : 0.0f;
                    const float e22_over_e53 = (var.cluster_e53[icluster] > 0) ? var.cluster_e22[icluster] / var.cluster_e53[icluster] : 0.0f;

                    // Feature order must match `config_npb_training.yaml` feature_list
                    std::vector<float> x_npb = {
                        cluster_Et_BDT,
                        cluster_Eta_BDT,
                        vertexz_BDT,
                        e11_over_e33,
                        e32_over_e35,
                        e11_over_e22,
                        e11_over_e13,
                        e11_over_e15,
                        e11_over_e17,
                        e11_over_e31,
                        e11_over_e51,
                        e11_over_e71,
                        e22_over_e33,
                        e22_over_e35,
                        e22_over_e37,
                        e22_over_e53,
                        cluster_weta_cogx_BDT,
                        cluster_wphi_cogx_BDT,
                        cluster_et1_BDT,
                        cluster_et2_BDT,
                        cluster_et3_BDT,
                        cluster_et4_BDT,
                        var.cluster_w32[icluster],
                        var.cluster_w52[icluster],
                        var.cluster_w72[icluster],
                    };

                    var.cluster_npb_score[icluster] = var.npb_bdt->Compute(x_npb)[0];
                }
            } // end cluster loop
        } // end variant loop

        // Fill ONCE per event after all variants have been scored
        outtree->Fill();
    }

    // Write the output tree to the file
    fout->cd();
    fout->Write();
    fout->Close();
}
