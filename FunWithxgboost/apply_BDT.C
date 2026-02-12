// apply_BDT.C
#include <TMVA/RBDT.hxx> // fast interpreter
#include <iostream>
#include <yaml-cpp/yaml.h>
#include <vector>
#include <string>
#include <cmath>

void apply_BDT(const std::string &configname = "config_nom.yaml", const std::string filetype = "data", const std::string inputfilename = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/condorout/part_10.root")
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
    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();

    // ------------------------------------------------------------
    // NPB score model (trained in `train_npb_score.py`)
    // ------------------------------------------------------------
    std::string npb_tmva_file = "npb_models/npb_score_tmva.root";
    if (configYaml["analysis"] && configYaml["analysis"]["npb_tmva_file"])
    {
        npb_tmva_file = configYaml["analysis"]["npb_tmva_file"].as<std::string>();
    }
    const bool do_npb_score = !gSystem->AccessPathName(npb_tmva_file.c_str());
    if (!do_npb_score)
    {
        std::cout << "WARNING: NPB TMVA model not found at '" << npb_tmva_file
                  << "'. Skipping NPB score application." << std::endl;
    }
    TMVA::Experimental::RBDT *npb_bdt = nullptr;
    if (do_npb_score)
    {
        npb_bdt = new TMVA::Experimental::RBDT("myBDT", npb_tmva_file);
    }

    // recoisoET definition (match `BDTinput.C`)
    const int conesize = configYaml["analysis"]["cone_size"].as<int>(3);
    const int iso_threshold = configYaml["analysis"]["iso_threshold"].as<int>(0);

    // NPB training phase-space (match `config_npb_training.yaml`)
    const float npb_et_min = 6.0;
    const float npb_et_max = 40.0;
    const float npb_eta_max = 0.7;

    // 1) Load the list of models
    std::vector<std::string> model_names = {
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
    std::vector<std::string> model_files = {
        "binned_models/model_base_single_tmva.root",
        "binned_models/model_base_vr_single_tmva.root",
        "binned_models/model_base_v0_single_tmva.root",
        "binned_models/model_base_v1_single_tmva.root",
        "binned_models/model_base_v2_single_tmva.root",
        "binned_models/model_base_v3_single_tmva.root",
        "binned_models/model_base_E_single_tmva.root",
        "binned_models/model_base_v0E_single_tmva.root",
        "binned_models/model_base_v1E_single_tmva.root",
        "binned_models/model_base_v2E_single_tmva.root",
        "binned_models/model_base_v3E_single_tmva.root",
    };

    std::vector<TMVA::Experimental::RBDT> bdt_list;

    for(int i = 0; i < model_names.size(); i++)
    {
        bdt_list.push_back(TMVA::Experimental::RBDT("myBDT", model_files[i]));
    }

    TFile *ftreein = new TFile(infilename.c_str(), "READ");
    TTree *slimtree = (TTree *)ftreein->Get(configYaml["input"]["tree"].as<std::string>().c_str());


    std::string outfile_name = filetype + "/bdt_1214.root";
    if (!issim)
    {
        //remove the .root
        std::string namebase = inputfilename.substr(0, inputfilename.find_last_of("."));
        outfile_name = namebase + "_with_bdt.root";
    }
    TFile *fout = new TFile(outfile_name.c_str(), "RECREATE");

    // clone the full structure, zero entries so far
    TTree *outtree = slimtree->CloneTree(0);


    static const int nclustercontainermx = 4096;


    //std::string leaf = Form("%s[ncluster_%s]/F",
    //                        bname.c_str(), clusternodename.c_str());
    //outtree->Branch(bname.c_str(), cluster_bdt, leaf.c_str());

    float cluster_Et_BDT = 0;
    float cluster_Eta_BDT = 0;
    float vertexz_BDT = 0;
    float e11_over_e33_BDT = 0;
    float e32_over_e35_BDT = 0;
    float cluster_weta_cogx_BDT = 0;
    float cluster_wphi_cogx_BDT = 0;
    float cluster_et1_BDT = 0;
    float cluster_et2_BDT = 0;
    float cluster_et3_BDT = 0;
    float cluster_et4_BDT = 0;

    // read from tree
    float vertexz;
    int ncluster;

    // Base inputs (also used by the photon-ID BDTs)
    float cluster_Et[nclustercontainermx], cluster_Eta[nclustercontainermx], cluster_Phi[nclustercontainermx];
    float cluster_et1[nclustercontainermx], cluster_et2[nclustercontainermx], cluster_et3[nclustercontainermx], cluster_et4[nclustercontainermx];
    float cluster_weta_cogx[nclustercontainermx], cluster_wphi_cogx[nclustercontainermx];

    // Extra inputs for NPB score
    float cluster_prob[nclustercontainermx];
    float cluster_e11[nclustercontainermx], cluster_e22[nclustercontainermx], cluster_e13[nclustercontainermx], cluster_e15[nclustercontainermx], cluster_e17[nclustercontainermx];
    float cluster_e31[nclustercontainermx], cluster_e51[nclustercontainermx], cluster_e71[nclustercontainermx];
    float cluster_e33[nclustercontainermx], cluster_e35[nclustercontainermx], cluster_e37[nclustercontainermx], cluster_e53[nclustercontainermx];
    float cluster_e32[nclustercontainermx];
    float cluster_w32[nclustercontainermx], cluster_w52[nclustercontainermx], cluster_w72[nclustercontainermx];

    // Isolation inputs for recoisoET (match `BDTinput.C`)
    float cluster_iso_02[nclustercontainermx], cluster_iso_03[nclustercontainermx], cluster_iso_04[nclustercontainermx];
    float cluster_iso_03_60_emcal[nclustercontainermx], cluster_iso_03_60_hcalin[nclustercontainermx], cluster_iso_03_60_hcalout[nclustercontainermx];

    //slimtree->SetBranchStatus("*", 0);
    slimtree->SetBranchStatus(Form("ncluster_%s", clusternodename.c_str()));
    slimtree->SetBranchAddress(Form("ncluster_%s", clusternodename.c_str()), &ncluster);

    slimtree->SetBranchStatus("vertexz", 1);
    slimtree->SetBranchAddress("vertexz", &vertexz);

    slimtree->SetBranchStatus(Form("cluster_Et_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_Et_%s", clusternodename.c_str()), &cluster_Et);
    slimtree->SetBranchStatus(Form("cluster_Eta_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_Eta_%s", clusternodename.c_str()), &cluster_Eta);
    slimtree->SetBranchStatus(Form("cluster_Phi_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_Phi_%s", clusternodename.c_str()), &cluster_Phi);
    slimtree->SetBranchStatus(Form("cluster_et1_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_et1_%s", clusternodename.c_str()), &cluster_et1);
    slimtree->SetBranchStatus(Form("cluster_et2_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_et2_%s", clusternodename.c_str()), &cluster_et2);
    slimtree->SetBranchStatus(Form("cluster_et3_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_et3_%s", clusternodename.c_str()), &cluster_et3);
    slimtree->SetBranchStatus(Form("cluster_et4_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_et4_%s", clusternodename.c_str()), &cluster_et4);
    slimtree->SetBranchStatus(Form("cluster_weta_cogx_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_weta_cogx_%s", clusternodename.c_str()), &cluster_weta_cogx);
    slimtree->SetBranchStatus(Form("cluster_wphi_cogx_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_wphi_cogx_%s", clusternodename.c_str()), &cluster_wphi_cogx);
    slimtree->SetBranchStatus(Form("cluster_prob_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_prob_%s", clusternodename.c_str()), &cluster_prob);
    slimtree->SetBranchStatus(Form("cluster_e11_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e11_%s", clusternodename.c_str()), &cluster_e11);
    slimtree->SetBranchStatus(Form("cluster_e22_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e22_%s", clusternodename.c_str()), &cluster_e22);
    slimtree->SetBranchStatus(Form("cluster_e13_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e13_%s", clusternodename.c_str()), &cluster_e13);
    slimtree->SetBranchStatus(Form("cluster_e15_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e15_%s", clusternodename.c_str()), &cluster_e15);
    slimtree->SetBranchStatus(Form("cluster_e17_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e17_%s", clusternodename.c_str()), &cluster_e17);
    slimtree->SetBranchStatus(Form("cluster_e31_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e31_%s", clusternodename.c_str()), &cluster_e31);
    slimtree->SetBranchStatus(Form("cluster_e51_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e51_%s", clusternodename.c_str()), &cluster_e51);
    slimtree->SetBranchStatus(Form("cluster_e71_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e71_%s", clusternodename.c_str()), &cluster_e71);
    slimtree->SetBranchStatus(Form("cluster_e33_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e33_%s", clusternodename.c_str()), &cluster_e33);
    slimtree->SetBranchStatus(Form("cluster_e32_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e32_%s", clusternodename.c_str()), &cluster_e32);
    slimtree->SetBranchStatus(Form("cluster_e35_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e35_%s", clusternodename.c_str()), &cluster_e35);
    slimtree->SetBranchStatus(Form("cluster_e37_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e37_%s", clusternodename.c_str()), &cluster_e37);
    slimtree->SetBranchStatus(Form("cluster_e53_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e53_%s", clusternodename.c_str()), &cluster_e53);

    slimtree->SetBranchStatus(Form("cluster_w32_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_w32_%s", clusternodename.c_str()), &cluster_w32);
    slimtree->SetBranchStatus(Form("cluster_w52_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_w52_%s", clusternodename.c_str()), &cluster_w52);
    slimtree->SetBranchStatus(Form("cluster_w72_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_w72_%s", clusternodename.c_str()), &cluster_w72);

    slimtree->SetBranchStatus(Form("cluster_iso_02_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_02_%s", clusternodename.c_str()), &cluster_iso_02);
    slimtree->SetBranchStatus(Form("cluster_iso_03_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_%s", clusternodename.c_str()), &cluster_iso_03);
    slimtree->SetBranchStatus(Form("cluster_iso_04_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_04_%s", clusternodename.c_str()), &cluster_iso_04);
    slimtree->SetBranchStatus(Form("cluster_iso_03_60_emcal_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_60_emcal_%s", clusternodename.c_str()), &cluster_iso_03_60_emcal);
    slimtree->SetBranchStatus(Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()), &cluster_iso_03_60_hcalin);
    slimtree->SetBranchStatus(Form("cluster_iso_03_60_hcalout_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_60_hcalout_%s", clusternodename.c_str()), &cluster_iso_03_60_hcalout);


    static const int NMODEL = model_names.size();
    float cluster_bdt[NMODEL][nclustercontainermx];
    float cluster_npb_score[nclustercontainermx];

    for (int i = 0; i < NMODEL; i++)
    {
        std::string bname = Form("cluster_bdt_%s_%s",
                                 clusternodename.c_str(), model_names[i].c_str());
        std::string leaf  = Form("%s[ncluster_%s]/F",
                                 bname.c_str(), clusternodename.c_str());
    
        outtree->Branch(bname.c_str(), cluster_bdt[i], leaf.c_str());
    }

    // NPB score output branch
    {
        std::string bname = Form("cluster_npb_score_%s", clusternodename.c_str());
        std::string leaf  = Form("%s[ncluster_%s]/F", bname.c_str(), clusternodename.c_str());
        outtree->Branch(bname.c_str(), cluster_npb_score, leaf.c_str());
    }

    int nentries = slimtree->GetEntries();
    //nentries = 100000;
    //float bdt_value = 0;
    for (int ientry = 0; ientry <nentries; ientry++)
    {

        if (ientry == 16152886)
            continue;
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        slimtree->GetEntry(ientry);


        //resize the cluster_bdt vector
        //cluster_bdt.resize(model_names.size(), std::vector<float>(ncluster));
        
        // loop over clusters   
        for (int icluster = 0; icluster < ncluster; icluster++)
        {
            cluster_Et_BDT = cluster_Et[icluster];
            cluster_Eta_BDT = cluster_Eta[icluster];
            vertexz_BDT = vertexz;
            e11_over_e33_BDT = (cluster_e11[icluster] > 0) ? cluster_e11[icluster] / cluster_e33[icluster] : 0;
            e32_over_e35_BDT = (cluster_e32[icluster] > 0) ? cluster_e32[icluster] / cluster_e35[icluster] : 0;
            cluster_weta_cogx_BDT = cluster_weta_cogx[icluster];
            cluster_wphi_cogx_BDT = cluster_wphi_cogx[icluster];
            cluster_et1_BDT = cluster_et1[icluster];
            cluster_et2_BDT = cluster_et2[icluster];
            cluster_et3_BDT = cluster_et3[icluster];
            cluster_et4_BDT = cluster_et4[icluster];

            // default: invalid (outside training phase-space or model missing)
            cluster_npb_score[icluster] = -1;


            std::vector<vector<float>> x_list = {
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

                
            for (int i = 0; i < model_names.size(); i++)
            {
                if (cluster_Et_BDT > 7)
                {
                    cluster_bdt[i][icluster] = bdt_list[i].Compute(x_list[i])[0];
                }
                else
                {
                    cluster_bdt[i][icluster] = -1;
                }
                
                // Debug output for first few clusters to understand the pattern
                if (icluster < 3 && ientry < 3) {
                    std::cout << "Entry " << ientry << ", Cluster " << icluster 
                             << ", ET=" << cluster_Et_BDT 
                             << ", Model " << model_names[i] 
                             << ": score=" << cluster_bdt[i][icluster] << std::endl;
                }
            }

            // NPB score application (single TMVA model)
            if (do_npb_score && npb_bdt &&
                cluster_Et_BDT >= npb_et_min && cluster_Et_BDT <= npb_et_max &&
                std::fabs(cluster_Eta_BDT) <= npb_eta_max)
            {
                // recoisoET (match `BDTinput.C`)
                float recoisoET = -999;
                if (conesize == 4)
                    recoisoET = cluster_iso_04[icluster];
                else if (conesize == 3)
                    recoisoET = cluster_iso_03[icluster];
                else if (conesize == 2)
                    recoisoET = cluster_iso_02[icluster];

                if (iso_threshold)
                {
                    recoisoET = cluster_iso_03_60_emcal[icluster] +
                               cluster_iso_03_60_hcalin[icluster] +
                               cluster_iso_03_60_hcalout[icluster];
                }

                // ratios (protect against divide-by-zero)
                const float e11_over_e33 = (cluster_e33[icluster] > 0) ? cluster_e11[icluster] / cluster_e33[icluster] : 0.0f;
                const float e32_over_e35 = (cluster_e35[icluster] > 0) ? cluster_e32[icluster] / cluster_e35[icluster] : 0.0f;
                const float e11_over_e22 = (cluster_e22[icluster] > 0) ? cluster_e11[icluster] / cluster_e22[icluster] : 0.0f;
                const float e11_over_e13 = (cluster_e13[icluster] > 0) ? cluster_e11[icluster] / cluster_e13[icluster] : 0.0f;
                const float e11_over_e15 = (cluster_e15[icluster] > 0) ? cluster_e11[icluster] / cluster_e15[icluster] : 0.0f;
                const float e11_over_e17 = (cluster_e17[icluster] > 0) ? cluster_e11[icluster] / cluster_e17[icluster] : 0.0f;
                const float e11_over_e31 = (cluster_e31[icluster] > 0) ? cluster_e11[icluster] / cluster_e31[icluster] : 0.0f;
                const float e11_over_e51 = (cluster_e51[icluster] > 0) ? cluster_e11[icluster] / cluster_e51[icluster] : 0.0f;
                const float e11_over_e71 = (cluster_e71[icluster] > 0) ? cluster_e11[icluster] / cluster_e71[icluster] : 0.0f;
                const float e22_over_e33 = (cluster_e33[icluster] > 0) ? cluster_e22[icluster] / cluster_e33[icluster] : 0.0f;
                const float e22_over_e35 = (cluster_e35[icluster] > 0) ? cluster_e22[icluster] / cluster_e35[icluster] : 0.0f;
                const float e22_over_e37 = (cluster_e37[icluster] > 0) ? cluster_e22[icluster] / cluster_e37[icluster] : 0.0f;
                const float e22_over_e53 = (cluster_e53[icluster] > 0) ? cluster_e22[icluster] / cluster_e53[icluster] : 0.0f;

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
                    cluster_w32[icluster],
                    cluster_w52[icluster],
                    cluster_w72[icluster],
                };

                cluster_npb_score[icluster] = npb_bdt->Compute(x_npb)[0];
            }
            
            //cluster_bdt[icluster] = bdt_value;
            // Fill the output tree
            //std::cout << "BDT response for entry " << ientry << ", icluster " << icluster << " = " << bdt_value << std::endl;
        }
        outtree->Fill();
    }
    // Write the output tree to the file
    fout->cd();
    fout->Write();
    fout->Close();
}
