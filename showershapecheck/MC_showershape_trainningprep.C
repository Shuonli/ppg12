#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>

void MC_showershape_trainningprep()
{

    std::string infilename = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run22/jet10/condorout/caloana.root";
    // std::string infilename = "/sphenix/tg/tg01/commissioning/CaloCalibWG/sli/ppg12/ana450/condorout/combine.root";
    TFile *ftreein = new TFile(infilename.c_str(), "READ");
    std::string outfilename = "MC_showershape_treeana.root";

    TTree *slimtree = (TTree *)ftreein->Get("slimtree");

    std::string clusternodename = "CLUSTERINFO_CEMC";

    int mbdnorthhit, mbdsouthhit;
    int pythiaid, nparticles;
    int ncluster;

    float vertexz;
    float vertexz_truth;

    static const int nparticlesmax = 10000;
    static const int nclustercontainermx = 10000;

    float particle_E[nparticlesmax], particle_Pt[nparticlesmax], particle_Eta[nparticlesmax], particle_Phi[nparticlesmax], particle_truth_iso_02[nparticlesmax], particle_truth_iso_03[nparticlesmax], particle_truth_iso_04[nparticlesmax];
    int particle_pid[nparticlesmax], particle_trkid[nparticlesmax], particle_photonclass[nparticlesmax], particle_converted[nparticlesmax];
    float cluster_E[nclustercontainermx], cluster_Et[nclustercontainermx], cluster_Eta[nclustercontainermx], cluster_Phi[nclustercontainermx], cluster_prob[nclustercontainermx], cluster_iso_02[nclustercontainermx], cluster_iso_03[nclustercontainermx], cluster_iso_04[nclustercontainermx], cluster_e1[nclustercontainermx], cluster_e2[nclustercontainermx], cluster_e3[nclustercontainermx], cluster_e4[nclustercontainermx], cluster_et1[nclustercontainermx], cluster_et2[nclustercontainermx], cluster_et3[nclustercontainermx], cluster_et4[nclustercontainermx], cluster_weta[nclustercontainermx], cluster_wphi[nclustercontainermx], cluster_ietacent[nclustercontainermx], cluster_iphicent[nclustercontainermx], cluster_e11[nclustercontainermx], cluster_e22[nclustercontainermx], cluster_e13[nclustercontainermx], cluster_e15[nclustercontainermx], cluster_e17[nclustercontainermx], cluster_e31[nclustercontainermx], cluster_e51[nclustercontainermx], cluster_e71[nclustercontainermx], cluster_e33[nclustercontainermx], cluster_e35[nclustercontainermx], cluster_e37[nclustercontainermx], cluster_e53[nclustercontainermx], cluster_e73[nclustercontainermx], cluster_e55[nclustercontainermx], cluster_e57[nclustercontainermx], cluster_e75[nclustercontainermx], cluster_e77[nclustercontainermx], cluster_w32[nclustercontainermx], cluster_e32[nclustercontainermx], cluster_w72[nclustercontainermx], cluster_e72[nclustercontainermx], cluster_ihcal_et[nclustercontainermx], cluster_ohcal_et[nclustercontainermx], cluster_ihcal_et22[nclustercontainermx], cluster_ohcal_et22[nclustercontainermx], cluster_ihcal_et33[nclustercontainermx], cluster_ohcal_et33[nclustercontainermx];
    float cluster_w52[nclustercontainermx];
    int cluster_truthtrkID[nclustercontainermx], cluster_pid[nclustercontainermx];

    int cluster_detamax[nclustercontainermx], cluster_dphimax[nclustercontainermx], cluster_ihcal_ieta[nclustercontainermx], cluster_ihcal_iphi[nclustercontainermx], cluster_ohcal_ieta[nclustercontainermx], cluster_ohcal_iphi[nclustercontainermx];

    slimtree->SetBranchAddress("mbdnorthhit", &mbdnorthhit);
    slimtree->SetBranchAddress("mbdsouthhit", &mbdsouthhit);
    slimtree->SetBranchAddress("pythiaid", &pythiaid);
    slimtree->SetBranchAddress("nparticles", &nparticles);
    slimtree->SetBranchAddress(Form("ncluster_%s", clusternodename.c_str()), &ncluster);
    slimtree->SetBranchAddress("vertexz", &vertexz);
    slimtree->SetBranchAddress("vertexz_truth", &vertexz_truth);

    slimtree->SetBranchAddress("particle_E", &particle_E);
    slimtree->SetBranchAddress("particle_Pt", &particle_Pt);
    slimtree->SetBranchAddress("particle_Eta", &particle_Eta);
    slimtree->SetBranchAddress("particle_Phi", &particle_Phi);
    slimtree->SetBranchAddress("particle_truth_iso_02", &particle_truth_iso_02);
    slimtree->SetBranchAddress("particle_truth_iso_03", &particle_truth_iso_03);
    slimtree->SetBranchAddress("particle_truth_iso_04", &particle_truth_iso_04);
    slimtree->SetBranchAddress("particle_pid", &particle_pid);
    slimtree->SetBranchAddress("particle_trkid", &particle_trkid);
    slimtree->SetBranchAddress("particle_photonclass", &particle_photonclass);
    slimtree->SetBranchAddress("particle_converted", &particle_converted);

    slimtree->SetBranchAddress(Form("cluster_E_%s", clusternodename.c_str()), &cluster_E);
    slimtree->SetBranchAddress(Form("cluster_Et_%s", clusternodename.c_str()), &cluster_Et);
    slimtree->SetBranchAddress(Form("cluster_Eta_%s", clusternodename.c_str()), &cluster_Eta);
    slimtree->SetBranchAddress(Form("cluster_Phi_%s", clusternodename.c_str()), &cluster_Phi);
    slimtree->SetBranchAddress(Form("cluster_prob_%s", clusternodename.c_str()), &cluster_prob);
    slimtree->SetBranchAddress(Form("cluster_truthtrkID_%s", clusternodename.c_str()), &cluster_truthtrkID);
    slimtree->SetBranchAddress(Form("cluster_pid_%s", clusternodename.c_str()), &cluster_pid);
    slimtree->SetBranchAddress(Form("cluster_iso_02_%s", clusternodename.c_str()), &cluster_iso_02);
    slimtree->SetBranchAddress(Form("cluster_iso_03_%s", clusternodename.c_str()), &cluster_iso_03);
    slimtree->SetBranchAddress(Form("cluster_iso_04_%s", clusternodename.c_str()), &cluster_iso_04);
    slimtree->SetBranchAddress(Form("cluster_e1_%s", clusternodename.c_str()), &cluster_e1);
    slimtree->SetBranchAddress(Form("cluster_e2_%s", clusternodename.c_str()), &cluster_e2);
    slimtree->SetBranchAddress(Form("cluster_e3_%s", clusternodename.c_str()), &cluster_e3);
    slimtree->SetBranchAddress(Form("cluster_e4_%s", clusternodename.c_str()), &cluster_e4);
    slimtree->SetBranchAddress(Form("cluster_et1_%s", clusternodename.c_str()), &cluster_et1);
    slimtree->SetBranchAddress(Form("cluster_et2_%s", clusternodename.c_str()), &cluster_et2);
    slimtree->SetBranchAddress(Form("cluster_et3_%s", clusternodename.c_str()), &cluster_et3);
    slimtree->SetBranchAddress(Form("cluster_et4_%s", clusternodename.c_str()), &cluster_et4);
    slimtree->SetBranchAddress(Form("cluster_weta_%s", clusternodename.c_str()), &cluster_weta);
    slimtree->SetBranchAddress(Form("cluster_wphi_%s", clusternodename.c_str()), &cluster_wphi);
    slimtree->SetBranchAddress(Form("cluster_ietacent_%s", clusternodename.c_str()), &cluster_ietacent);
    slimtree->SetBranchAddress(Form("cluster_iphicent_%s", clusternodename.c_str()), &cluster_iphicent);
    slimtree->SetBranchAddress(Form("cluster_detamax_%s", clusternodename.c_str()), &cluster_detamax);
    slimtree->SetBranchAddress(Form("cluster_dphimax_%s", clusternodename.c_str()), &cluster_dphimax);

    slimtree->SetBranchAddress(Form("cluster_e11_%s", clusternodename.c_str()), &cluster_e11);
    slimtree->SetBranchAddress(Form("cluster_e22_%s", clusternodename.c_str()), &cluster_e22);
    slimtree->SetBranchAddress(Form("cluster_e13_%s", clusternodename.c_str()), &cluster_e13);
    slimtree->SetBranchAddress(Form("cluster_e15_%s", clusternodename.c_str()), &cluster_e15);
    slimtree->SetBranchAddress(Form("cluster_e17_%s", clusternodename.c_str()), &cluster_e17);
    slimtree->SetBranchAddress(Form("cluster_e31_%s", clusternodename.c_str()), &cluster_e31);
    slimtree->SetBranchAddress(Form("cluster_e51_%s", clusternodename.c_str()), &cluster_e51);
    slimtree->SetBranchAddress(Form("cluster_e71_%s", clusternodename.c_str()), &cluster_e71);
    slimtree->SetBranchAddress(Form("cluster_e33_%s", clusternodename.c_str()), &cluster_e33);
    slimtree->SetBranchAddress(Form("cluster_e35_%s", clusternodename.c_str()), &cluster_e35);
    slimtree->SetBranchAddress(Form("cluster_e37_%s", clusternodename.c_str()), &cluster_e37);
    slimtree->SetBranchAddress(Form("cluster_e53_%s", clusternodename.c_str()), &cluster_e53);
    slimtree->SetBranchAddress(Form("cluster_e73_%s", clusternodename.c_str()), &cluster_e73);
    slimtree->SetBranchAddress(Form("cluster_e55_%s", clusternodename.c_str()), &cluster_e55);
    slimtree->SetBranchAddress(Form("cluster_e57_%s", clusternodename.c_str()), &cluster_e57);
    slimtree->SetBranchAddress(Form("cluster_e75_%s", clusternodename.c_str()), &cluster_e75);
    slimtree->SetBranchAddress(Form("cluster_e77_%s", clusternodename.c_str()), &cluster_e77);
    slimtree->SetBranchAddress(Form("cluster_w32_%s", clusternodename.c_str()), &cluster_w32);
    slimtree->SetBranchAddress(Form("cluster_e32_%s", clusternodename.c_str()), &cluster_e32);
    slimtree->SetBranchAddress(Form("cluster_w72_%s", clusternodename.c_str()), &cluster_w72);
    slimtree->SetBranchAddress(Form("cluster_e72_%s", clusternodename.c_str()), &cluster_e72);
    slimtree->SetBranchAddress(Form("cluster_w52_%s", clusternodename.c_str()), &cluster_w52);

    slimtree->SetBranchAddress(Form("cluster_ihcal_et_%s", clusternodename.c_str()), &cluster_ihcal_et);
    slimtree->SetBranchAddress(Form("cluster_ohcal_et_%s", clusternodename.c_str()), &cluster_ohcal_et);
    slimtree->SetBranchAddress(Form("cluster_ihcal_et22_%s", clusternodename.c_str()), &cluster_ihcal_et22);
    slimtree->SetBranchAddress(Form("cluster_ohcal_et22_%s", clusternodename.c_str()), &cluster_ohcal_et22);
    slimtree->SetBranchAddress(Form("cluster_ihcal_et33_%s", clusternodename.c_str()), &cluster_ihcal_et33);
    slimtree->SetBranchAddress(Form("cluster_ohcal_et33_%s", clusternodename.c_str()), &cluster_ohcal_et33);
    slimtree->SetBranchAddress(Form("cluster_ihcal_ieta_%s", clusternodename.c_str()), &cluster_ihcal_ieta);
    slimtree->SetBranchAddress(Form("cluster_ihcal_iphi_%s", clusternodename.c_str()), &cluster_ihcal_iphi);
    slimtree->SetBranchAddress(Form("cluster_ohcal_ieta_%s", clusternodename.c_str()), &cluster_ohcal_ieta);
    slimtree->SetBranchAddress(Form("cluster_ohcal_iphi_%s", clusternodename.c_str()), &cluster_ohcal_iphi);

    // a csv file to store the training data
    std::string outcsvfilename = "MC_showershape_training_jet10.csv";
    std::ofstream outcsvfile;
    outcsvfile.open(outcsvfilename.c_str());
    //make sure the file is empty
    outcsvfile << "cluster_ET,rhad22,rhad33,reta77,rphi77,reta55,rphi55,reta,rphi,re11_E,w72,w52,weta,wphi,et1,et2,et3,et4,prob,is_signal_class" << std::endl;
    int nentries = slimtree->GetEntries();
    for (int ientry = 0; ientry < nentries; ientry++)
    {
        if (ientry % 1000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        slimtree->GetEntry(ientry);
        std::map<int, int> particle_trkidmap;
        for (int iparticle = 0; iparticle < nparticles; iparticle++)
        {
            particle_trkidmap[particle_trkid[iparticle]] = iparticle;
        }
        for (int icluster = 0; icluster < ncluster; icluster++)
        {
            // need ET > 10 GeV
            if (cluster_Et[icluster] < 10)
                continue;

            float rhad22 = (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]) / (cluster_Et[icluster]+ (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]) );
            float rhad33 = (cluster_ihcal_et33[icluster] + cluster_ohcal_et33[icluster]) / (cluster_Et[icluster] + (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]) );

            float reta77 = cluster_e37[icluster] / cluster_e77[icluster];
            float rphi77 = cluster_e73[icluster] / cluster_e77[icluster];

            float reta55 = cluster_e35[icluster] / cluster_e55[icluster];
            float rphi55 = cluster_e53[icluster] / cluster_e55[icluster];

            float reta = cluster_e33[icluster] / cluster_e73[icluster];

            float rphi = cluster_e33[icluster] / cluster_e37[icluster];

            float re11_E = cluster_e11[icluster] / cluster_E[icluster];

            bool is_signal_class = true;

            if (particle_trkidmap.find(cluster_truthtrkID[icluster]) == particle_trkidmap.end())
            {
                // std::cout<<"trackid: "<<cluster_truthtrkID[icluster]<<std::endl;
                // std::cout << "Error: cluster_truthtrkID not found in particle_trkidmap" << std::endl;
                continue;
            }
            int iparticle = particle_trkidmap[cluster_truthtrkID[icluster]];
            if (iparticle < 0 || iparticle >= nparticles)
            {
                std::cout << "Error: iparticle = " << iparticle << " out of range" << std::endl;
                continue;
            }
            if (particle_Pt[iparticle] < 8)
                continue;
            //if converted
            if (particle_converted[iparticle] == 1)
                continue;
            if (particle_pid[iparticle] != 22)
                is_signal_class = false;
            if (particle_photonclass[iparticle] == 3)
                is_signal_class = false;

            // cluster_ET, rhad22, rhad33, reta77, rphi77, reta55, rphi55, reta, rphi, re11_E, w72, w52, weta, wphi, et1, et2, et3, et4, prob, is_signal_class
            outcsvfile << cluster_Et[icluster] << ","
                       << rhad22 << ","
                       << rhad33 << ","
                       << reta77 << ","
                       << rphi77 << ","
                       << reta55 << ","
                       << rphi55 << ","
                       << reta << ","
                       << rphi << ","
                       << re11_E << ","
                       << cluster_w72[icluster] << ","
                       << cluster_w52[icluster] << ","
                       << cluster_weta[icluster] << ","
                       << cluster_wphi[icluster] << ","
                       << cluster_et1[icluster] << ","
                       << cluster_et2[icluster] << ","
                       << cluster_et3[icluster] << ","
                       << cluster_et4[icluster] << ","
                       << cluster_prob[icluster] << ","
                       << is_signal_class << std::endl;
        }
    }
    outcsvfile.close();
}