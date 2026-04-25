#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

void Cluster_timing()
{
    std::cout << "Cluster_timing" << std::endl;

    std::string infilename = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana462/condorout/combine0128.root";
    // std::string infilename = "/sphenix/tg/tg01/commissioning/CaloCalibWG/sli/ppg12/ana462/ana462_newlist/OUTTREE_run2pp_ana462_2024p010-00047810.root";
    //  std::string infilename = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run22/photon10/condorout_waveform/caloana.root";
    //   std::string infilename = "/sphenix/tg/tg01/commissioning/CaloCalibWG/sli/ppg12/ana450/condorout/combine.root";
    std::cout << "Opening File" << std::endl;

    TFile *ftreein = new TFile(infilename.c_str(), "READ");
    std::string outfilename = "clustertime.root";
    TTree *slimtree = (TTree *)ftreein->Get("slimtree");

    std::string clusternodename = "CLUSTERINFO_CEMC_NO_SPLIT";

    int mbdnorthhit, mbdsouthhit;
    int pythiaid, nparticles;
    int ncluster;

    float vertexz;
    float vertexz_truth;

    static const int nparticlesmax = 10000;
    static const int nclustercontainermx = 100;

    static const int arraysize = 49;

    float particle_E[nparticlesmax], particle_Pt[nparticlesmax], particle_Eta[nparticlesmax], particle_Phi[nparticlesmax], particle_truth_iso_02[nparticlesmax], particle_truth_iso_03[nparticlesmax], particle_truth_iso_04[nparticlesmax];
    int particle_pid[nparticlesmax], particle_trkid[nparticlesmax], particle_photonclass[nparticlesmax], particle_converted[nparticlesmax];
    float cluster_E[nclustercontainermx], cluster_Et[nclustercontainermx], cluster_Eta[nclustercontainermx], cluster_Phi[nclustercontainermx], cluster_prob[nclustercontainermx], cluster_iso_02[nclustercontainermx], cluster_iso_03[nclustercontainermx], cluster_iso_04[nclustercontainermx], cluster_e1[nclustercontainermx], cluster_e2[nclustercontainermx], cluster_e3[nclustercontainermx], cluster_e4[nclustercontainermx], cluster_et1[nclustercontainermx], cluster_et2[nclustercontainermx], cluster_et3[nclustercontainermx], cluster_et4[nclustercontainermx], cluster_weta[nclustercontainermx], cluster_wphi[nclustercontainermx], cluster_ietacent[nclustercontainermx], cluster_iphicent[nclustercontainermx], cluster_e11[nclustercontainermx], cluster_e22[nclustercontainermx], cluster_e13[nclustercontainermx], cluster_e15[nclustercontainermx], cluster_e17[nclustercontainermx], cluster_e31[nclustercontainermx], cluster_e51[nclustercontainermx], cluster_e71[nclustercontainermx], cluster_e33[nclustercontainermx], cluster_e35[nclustercontainermx], cluster_e37[nclustercontainermx], cluster_e53[nclustercontainermx], cluster_e73[nclustercontainermx], cluster_e55[nclustercontainermx], cluster_e57[nclustercontainermx], cluster_e75[nclustercontainermx], cluster_e77[nclustercontainermx], cluster_w32[nclustercontainermx], cluster_e32[nclustercontainermx], cluster_w72[nclustercontainermx], cluster_e72[nclustercontainermx], cluster_ihcal_et[nclustercontainermx], cluster_ohcal_et[nclustercontainermx], cluster_ihcal_et22[nclustercontainermx], cluster_ohcal_et22[nclustercontainermx], cluster_ihcal_et33[nclustercontainermx], cluster_ohcal_et33[nclustercontainermx];

    int cluster_truthtrkID[nclustercontainermx], cluster_pid[nclustercontainermx];

    int cluster_detamax[nclustercontainermx], cluster_dphimax[nclustercontainermx], cluster_ihcal_ieta[nclustercontainermx], cluster_ihcal_iphi[nclustercontainermx], cluster_ohcal_ieta[nclustercontainermx], cluster_ohcal_iphi[nclustercontainermx];

    int cluster_ownership_array[nclustercontainermx][arraysize] = {0};

    float cluster_time_array[nclustercontainermx][arraysize] = {0};

    float cluster_e_array[nclustercontainermx][arraysize] = {0};

    float cluster_adc_array[nclustercontainermx][arraysize] = {0};

    int cluster_e_array_idx[nclustercontainermx][arraysize] = {0};

    int cluster_status_array[nclustercontainermx][arraysize] = {0};

    float cluster_weta_cogx[nclustercontainermx], cluster_wphi_cogx[nclustercontainermx];

    // jet stuff

    int njet;

    static const int njetmax = 100;

    float jet_E[njetmax], jet_Pt[njetmax], jet_Eta[njetmax], jet_Phi[njetmax];

    std::cout << "Setting Branch Addresses" << std::endl;

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
    slimtree->SetBranchAddress(Form("cluster_weta_cogx_%s", clusternodename.c_str()), &cluster_weta_cogx);
    slimtree->SetBranchAddress(Form("cluster_wphi_cogx_%s", clusternodename.c_str()), &cluster_wphi_cogx);

    slimtree->SetBranchAddress(Form("cluster_e_array_%s", clusternodename.c_str()), cluster_e_array);
    // slimtree->SetBranchAddress(Form("cluster_adc_array_%s", clusternodename.c_str()), &cluster_adc_array);
    // slimtree->SetBranchAddress(Form("cluster_e_array_idx_%s", clusternodename.c_str()), &cluster_e_array_idx);
    // slimtree->SetBranchAddress(Form("cluster_status_array_%s", clusternodename.c_str()), &cluster_status_array);
    slimtree->SetBranchAddress(Form("cluster_time_array_%s", clusternodename.c_str()), &cluster_time_array);
    slimtree->SetBranchAddress(Form("cluster_ownership_array_%s", clusternodename.c_str()), &cluster_ownership_array);

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

    slimtree->SetBranchAddress("njet", &njet);
    slimtree->SetBranchAddress("jet_E", &jet_E);
    slimtree->SetBranchAddress("jet_Pt", &jet_Pt);
    slimtree->SetBranchAddress("jet_Eta", &jet_Eta);
    slimtree->SetBranchAddress("jet_Phi", &jet_Phi);

    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    std::cout << "Creating output file: " << outfilename << std::endl;

    TH2F *h_cluster_time_withjet = new TH2F("h_cluster_time_withjet", "h_cluster_time_withjet", 1000, -6, 6, 100, 0, 100);

    TH2F *h_cluster_time_large_wetacogx = new TH2F("h_cluster_time_large_wetacogx", "h_cluster_time_large_wetacogx", 1000, -6, 6, 100, 0, 100);

    TH2F *h_cluster_time_small_wetacogx = new TH2F("h_cluster_time_small_wetacogx", "h_cluster_time_small_wetacogx", 1000, -6, 6, 100, 0, 100);

    TH2F *h_cluster_time_signal = new TH2F("h_cluster_time_signal", "h_cluster_time_signal", 1000, -6, 6, 100, 0, 100);

    TH2F *h_cluster_time = new TH2F("h_cluster_time", "h_cluster_time", 1000, -6, 6, 100, 0, 100);

    TH3F *h_cluster_time_large_wetacogx_eta = new TH3F("h_cluster_time_large_wetacogx_eta", "h_cluster_time_large_wetacogx_eta", 1000, -6, 6, 100, 0, 100, 100, -1, 1);
    TH3F *h_cluster_time_eta = new TH3F("h_cluster_time_eta", "h_cluster_time_eta", 1000, -6, 6, 100, 0, 100, 100, -1, 1);
    TH3F *h_cluster_time_withjet_eta = new TH3F("h_cluster_time_withjet_eta", "h_cluster_time_withjet_eta", 1000, -6, 6, 100, 0, 100, 100, -1, 1);

    int nentries = slimtree->GetEntries();
    std::cout << "nentries: " << nentries << std::endl;
    for (int ientry = 0; ientry < nentries; ientry++)
    {
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        slimtree->GetEntry(ientry);

        // check vertex cut
        if (abs(vertexz) > 30)
            continue;

        if (!(mbdnorthhit >= 1 && mbdsouthhit >= 1))
            continue;
        // loop over jets
        std::vector<float> jetphi;

        for (int ijet = 0; ijet < njet; ijet++)
        {
            if (jet_Pt[ijet] < 10)
                continue;
            jetphi.push_back(jet_Phi[ijet]);
        }

        for (int icluster = 0; icluster < ncluster; icluster++)
        {
            // need ET > 10 GeV
            if (cluster_Et[icluster] < 10)
                continue;

            bool otherside_jet = false;

            for (int ijet = 0; ijet < (int)jetphi.size(); ijet++)
            {
                float dphi = cluster_Phi[icluster] - jetphi[ijet];

                while (dphi > M_PI)
                    dphi = dphi - 2 * M_PI;
                while (dphi < -M_PI)
                    dphi = dphi + 2 * M_PI;

                if (abs(dphi) > (7 * M_PI / 8))
                {
                    otherside_jet = true;
                    break;
                }
            }
            // check cluster eta
            if (abs(cluster_Eta[icluster]) > 0.7)
                continue;
            // find cluster average time

            float clusteravgtime = 0;
            float cluster_total_e = 0;
            for (int i = 0; i < 49; i++)
            {

                if (cluster_ownership_array[icluster][i] == 1)
                {
                    clusteravgtime += cluster_time_array[icluster][i] * cluster_e_array[icluster][i];
                    // std::cout<<"cluster_time_array[icluster][i]: "<<cluster_time_array[icluster][i]<<std::endl;
                    cluster_total_e += cluster_e_array[icluster][i];
                }
            }
            clusteravgtime = cluster_total_e > 0 ? clusteravgtime / cluster_total_e : 0;
            // std::cout<<"clusteravgtime: "<<clusteravgtime<< " cluster_total_e: "<<cluster_total_e<<std::endl;

            // std::cout<<"clusteravgtime: "<<clusteravgtime<<std::endl;

            h_cluster_time->Fill(clusteravgtime, cluster_Et[icluster]);
            h_cluster_time_eta->Fill(clusteravgtime, cluster_Et[icluster], cluster_Eta[icluster]);

            if (otherside_jet)
            {
                h_cluster_time_withjet->Fill(clusteravgtime, cluster_Et[icluster]);
                h_cluster_time_withjet_eta->Fill(clusteravgtime, cluster_Et[icluster], cluster_Eta[icluster]);
            }

            if (cluster_weta_cogx[icluster] > 1.0)
            {
                h_cluster_time_large_wetacogx->Fill(clusteravgtime, cluster_Et[icluster]);
                h_cluster_time_large_wetacogx_eta->Fill(clusteravgtime, cluster_Et[icluster], cluster_Eta[icluster]);
            }
            else if (cluster_weta_cogx[icluster] < 0.4)
            {
                h_cluster_time_small_wetacogx->Fill(clusteravgtime, cluster_Et[icluster]);
                if (otherside_jet)
                {
                    h_cluster_time_signal->Fill(clusteravgtime, cluster_Et[icluster]);
                }
            }
        }
    }

    fout->cd();
    fout->Write();
    fout->Close();
}
