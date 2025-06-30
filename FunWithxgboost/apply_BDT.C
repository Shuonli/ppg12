// apply_BDT.C
#include <TMVA/RBDT.hxx> // fast interpreter
#include <iostream>
#include <yaml-cpp/yaml.h>

void apply_BDT(const std::string &configname = "config_nom.yaml", const std::string filetype = "jet30")
{
    using namespace TMVA::Experimental;
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    bool issim = true;

    if (filetype == "data")
    {
        issim = false;
    }
    bool isbackground = false;
    float max_photon_lower = 0;
    float max_photon_upper = 100;
    // unit in pb
    const float photon5cross = 2.017e+08 * 0.000442571;
    const float photon10cross = 3.688e+07 * 0.000181474;
    const float photon20cross = 1.571e+05 * 0.000673448;

    // Hanpu uses unit in b
    const float jet10cross = 3.646e-6;
    const float jet15cross = 36864930.0 * 0.011059973 * 1e-12;
    const float jet20cross = 1392140.9 * 0.042 * 1e-12;
    const float jet30cross = 2.505e-9;

    float max_jet_lower = 0;
    float max_jet_upper = 100;

    float energy_scale_lower = 0;
    float energy_scale_upper = 100;

    float cluster_ET_upper = 100;

    float weight = 1.0;
    float vertex_weight = 1.0;
    float cross_weight = 1.0;

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
        max_jet_upper = 19;
        energy_scale_lower = 10;
        energy_scale_upper = 16;
        cluster_ET_upper = 25;
        weight = jet10cross / jet30cross;
        isbackground = true;
    }
    else if (filetype == "jet15")
    {
        max_jet_lower = 19;
        max_jet_upper = 23;
        energy_scale_lower = 16;
        energy_scale_upper = 20;
        cluster_ET_upper = 25;
        weight = jet15cross / jet30cross;
        isbackground = true;
    }
    else if (filetype == "jet20")
    {
        max_jet_lower = 23;
        max_jet_upper = 30;
        energy_scale_lower = 20;
        energy_scale_upper = 30;
        weight = jet20cross / jet30cross;
        isbackground = true;
    }
    else if (filetype == "jet30")
    {
        max_jet_lower = 30;
        max_jet_upper = 100;
        energy_scale_lower = 30;
        energy_scale_upper = 100;
        weight = 1.0;
        isbackground = true;
    }

    std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();

    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();

    std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;

    if (!issim)
    {
        infilename = configYaml["input"]["data_file"].as<std::string>();
    }
    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();

    // 1) Load the model once
    TMVA::Experimental::RBDT bdt("myBDT", "myBDT.root"); // name, file created above

    TFile *ftreein = new TFile(infilename.c_str(), "READ");
    TTree *slimtree = (TTree *)ftreein->Get(configYaml["input"]["tree"].as<std::string>().c_str());

    TFile *fout = new TFile("outfile_with_bdt.root", "RECREATE");

    // clone the full structure, zero entries so far
    TTree *outtree = slimtree->CloneTree(0);

    static const int nclustercontainermx = 50;

    float cluster_bdt[nclustercontainermx];

    std::string bname = Form("cluster_bdt_%s", clusternodename.c_str());
    std::string leaf = Form("%s[ncluster_%s]/F",
                            bname.c_str(), clusternodename.c_str());
    outtree->Branch(bname.c_str(), cluster_bdt, leaf.c_str());

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

    float cluster_E[nclustercontainermx], cluster_Et[nclustercontainermx], cluster_Eta[nclustercontainermx], cluster_Phi[nclustercontainermx], cluster_prob[nclustercontainermx], cluster_iso_02[nclustercontainermx], cluster_iso_03[nclustercontainermx], cluster_iso_04[nclustercontainermx], cluster_e1[nclustercontainermx], cluster_e2[nclustercontainermx], cluster_e3[nclustercontainermx], cluster_e4[nclustercontainermx], cluster_et1[nclustercontainermx], cluster_et2[nclustercontainermx], cluster_et3[nclustercontainermx], cluster_et4[nclustercontainermx], cluster_weta[nclustercontainermx], cluster_wphi[nclustercontainermx], cluster_ietacent[nclustercontainermx], cluster_iphicent[nclustercontainermx], cluster_e11[nclustercontainermx], cluster_e22[nclustercontainermx], cluster_e13[nclustercontainermx], cluster_e15[nclustercontainermx], cluster_e17[nclustercontainermx], cluster_e31[nclustercontainermx], cluster_e51[nclustercontainermx], cluster_e71[nclustercontainermx], cluster_e33[nclustercontainermx], cluster_e35[nclustercontainermx], cluster_e37[nclustercontainermx], cluster_e53[nclustercontainermx], cluster_e73[nclustercontainermx], cluster_e55[nclustercontainermx], cluster_e57[nclustercontainermx], cluster_e75[nclustercontainermx], cluster_e77[nclustercontainermx], cluster_w32[nclustercontainermx], cluster_e32[nclustercontainermx], cluster_w72[nclustercontainermx], cluster_e72[nclustercontainermx], cluster_ihcal_et[nclustercontainermx], cluster_ohcal_et[nclustercontainermx], cluster_ihcal_et22[nclustercontainermx], cluster_ohcal_et22[nclustercontainermx], cluster_ihcal_et33[nclustercontainermx], cluster_ohcal_et33[nclustercontainermx];
    float cluster_w52[nclustercontainermx];
    int cluster_truthtrkID[nclustercontainermx], cluster_pid[nclustercontainermx];

    int cluster_detamax[nclustercontainermx], cluster_dphimax[nclustercontainermx], cluster_ihcal_ieta[nclustercontainermx], cluster_ihcal_iphi[nclustercontainermx], cluster_ohcal_ieta[nclustercontainermx], cluster_ohcal_iphi[nclustercontainermx];
    float cluster_CNN_prob[nclustercontainermx];

    float cluster_weta_cogx[nclustercontainermx], cluster_wphi_cogx[nclustercontainermx];

    slimtree->SetBranchStatus("*", 0);
    slimtree->SetBranchStatus(Form("ncluster_%s", clusternodename.c_str()));
    slimtree->SetBranchAddress(Form("ncluster_%s", clusternodename.c_str()), &ncluster);

    slimtree->SetBranchStatus("vertexz", 1);
    slimtree->SetBranchAddress("vertexz", &vertexz);

    slimtree->SetBranchStatus(Form("cluster_E_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_E_%s", clusternodename.c_str()), &cluster_E);
    slimtree->SetBranchStatus(Form("cluster_Et_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_Et_%s", clusternodename.c_str()), &cluster_Et);
    slimtree->SetBranchStatus(Form("cluster_Eta_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_Eta_%s", clusternodename.c_str()), &cluster_Eta);
    slimtree->SetBranchStatus(Form("cluster_Phi_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_Phi_%s", clusternodename.c_str()), &cluster_Phi);
    slimtree->SetBranchStatus(Form("cluster_prob_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_prob_%s", clusternodename.c_str()), &cluster_prob);
    slimtree->SetBranchStatus(Form("cluster_CNN_prob_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_CNN_prob_%s", clusternodename.c_str()), &cluster_CNN_prob);
    slimtree->SetBranchStatus(Form("cluster_truthtrkID_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_truthtrkID_%s", clusternodename.c_str()), &cluster_truthtrkID);
    slimtree->SetBranchStatus(Form("cluster_pid_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_pid_%s", clusternodename.c_str()), &cluster_pid);
    slimtree->SetBranchStatus(Form("cluster_iso_02_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_02_%s", clusternodename.c_str()), &cluster_iso_02);
    slimtree->SetBranchStatus(Form("cluster_iso_03_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_03_%s", clusternodename.c_str()), &cluster_iso_03);
    slimtree->SetBranchStatus(Form("cluster_iso_04_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_iso_04_%s", clusternodename.c_str()), &cluster_iso_04);
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
    slimtree->SetBranchStatus(Form("cluster_e35_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e35_%s", clusternodename.c_str()), &cluster_e35);
    slimtree->SetBranchStatus(Form("cluster_e37_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e37_%s", clusternodename.c_str()), &cluster_e37);
    slimtree->SetBranchStatus(Form("cluster_e53_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e53_%s", clusternodename.c_str()), &cluster_e53);
    slimtree->SetBranchStatus(Form("cluster_e73_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e73_%s", clusternodename.c_str()), &cluster_e73);
    slimtree->SetBranchStatus(Form("cluster_e55_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e55_%s", clusternodename.c_str()), &cluster_e55);
    slimtree->SetBranchStatus(Form("cluster_e57_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e57_%s", clusternodename.c_str()), &cluster_e57);
    slimtree->SetBranchStatus(Form("cluster_e75_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e75_%s", clusternodename.c_str()), &cluster_e75);
    slimtree->SetBranchStatus(Form("cluster_e77_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e77_%s", clusternodename.c_str()), &cluster_e77);
    slimtree->SetBranchStatus(Form("cluster_w32_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_w32_%s", clusternodename.c_str()), &cluster_w32);
    slimtree->SetBranchStatus(Form("cluster_e32_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e32_%s", clusternodename.c_str()), &cluster_e32);
    slimtree->SetBranchStatus(Form("cluster_w72_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_w72_%s", clusternodename.c_str()), &cluster_w72);
    slimtree->SetBranchStatus(Form("cluster_e72_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_e72_%s", clusternodename.c_str()), &cluster_e72);
    slimtree->SetBranchStatus(Form("cluster_w52_%s", clusternodename.c_str()), 1);
    slimtree->SetBranchAddress(Form("cluster_w52_%s", clusternodename.c_str()), &cluster_w52);

    int nentries = slimtree->GetEntries();
    float bdt_value = 0;
    for (int ientry = 0; ientry <nentries; ientry++)
    {

        if (ientry == 16152886)
            continue;
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        slimtree->GetEntry(ientry);

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

            std::vector<float> x = {
                cluster_Et_BDT, cluster_Eta_BDT, vertexz_BDT, e11_over_e33_BDT,
                e32_over_e35_BDT, cluster_weta_cogx_BDT, cluster_wphi_cogx_BDT,
                cluster_et1_BDT, cluster_et2_BDT, cluster_et3_BDT, cluster_et4_BDT};

            // compute BDT score
            if (cluster_Et_BDT > 10)
            {
                bdt_value = bdt.Compute(x)[0];
            }
            else
            {
                bdt_value = -1; // or some other value indicating low Et
            }
            cluster_bdt[icluster] = bdt_value;
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
