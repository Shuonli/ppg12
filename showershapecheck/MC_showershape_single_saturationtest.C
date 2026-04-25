#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>

bool get_status_bit(int bit, uint8_t _status)
{
    if (bit < 0 || bit > 7)
    {
        return false; // default behavior
    }
    return (_status & ((uint8_t)1 << bit)) != 0;
}

void MC_showershape_single_saturationtest()
{

    std::string infilename = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run24/photon/condor_sat/caloana.root";
    // std::string infilename = "/sphenix/tg/tg01/commissioning/CaloCalibWG/sli/ppg12/ana450/condorout/combine.root";
    TFile *ftreein = new TFile(infilename.c_str(), "READ");
    std::string outfilename = "MC_single_sat_new.root";

    TTree *slimtree = (TTree *)ftreein->Get("slimtree");

    std::vector<float> eta_bins = {-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9};

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

    int cluster_truthtrkID[nclustercontainermx], cluster_pid[nclustercontainermx];

    int cluster_detamax[nclustercontainermx], cluster_dphimax[nclustercontainermx], cluster_ihcal_ieta[nclustercontainermx], cluster_ihcal_iphi[nclustercontainermx], cluster_ohcal_ieta[nclustercontainermx], cluster_ohcal_iphi[nclustercontainermx];

    int cluster_status_array[nclustercontainermx][49] = {0};
    float cluster_e_array[nclustercontainermx][49] = {0};

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

    slimtree->SetBranchAddress(Form("cluster_e_array_%s", clusternodename.c_str()), &cluster_e_array);
    slimtree->SetBranchAddress(Form("cluster_status_array_%s", clusternodename.c_str()), &cluster_status_array);
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

    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    // TEfficiency for conversion and reco
    TEfficiency::EStatOption effopt = TEfficiency::kBUniform;

    std::vector<TEfficiency *> eff_sat_prob;
    std::vector<TEfficiency *> eff_sat_badchi2_prob;
    std::vector<TEfficiency *> eff_reco_eff;
    std::vector<TEfficiency *> eff_iso_eta;

    for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
    {
        eff_sat_prob.push_back(new TEfficiency(Form("eff_sat_prob_eta_%d", ieta), ";p_{T} [GeV];Saturation Prob", 50, 10, 60));
        eff_sat_prob[ieta]->SetStatisticOption(effopt);
        eff_sat_badchi2_prob.push_back(new TEfficiency(Form("eff_sat_badchi2_prob_eta_%d", ieta), ";p_{T} [GeV];Saturation with bad chi2 Prob", 50, 10, 60));
        eff_sat_badchi2_prob[ieta]->SetStatisticOption(effopt);
        eff_reco_eff.push_back(new TEfficiency(Form("eff_reco_eff_eta_%d", ieta), ";p_{T} [GeV];Reco Efficiency", 50, 10, 60));
        eff_reco_eff[ieta]->SetStatisticOption(effopt);
    }
    TEfficiency *eff_tower_sat_prob = new TEfficiency("eff_tower_sat_prob", ";E [GeV];Saturation Prob", 20, 10, 50);
    eff_tower_sat_prob->SetStatisticOption(effopt);
    TEfficiency *eff_tower_sat_badchi2_prob = new TEfficiency("eff_tower_sat_badchi2_prob", ";E [GeV];Saturation with bad chi2 Prob", 20, 10, 50);
    eff_tower_sat_badchi2_prob->SetStatisticOption(effopt);

    int nentries = slimtree->GetEntries();
    for (int ientry = 0; ientry < nentries; ientry++)
    {
        if (ientry % 10000 == 0)
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;

        
        slimtree->GetEntry(ientry);
        //vertex cut
        if (abs(vertexz) > 30)
            continue;
        std::map<int, int> particle_trkidmap;
        // map for photon reco eff
        std::map<int, bool> photon_reco;
        std::map<int, bool> photon_sat;
        std::map<int, bool> photon_sat_badchi2;
        for (int iparticle = 0; iparticle < nparticles; iparticle++)
        {
            particle_trkidmap[particle_trkid[iparticle]] = iparticle;
            photon_reco[iparticle] = false;
            photon_sat[iparticle] = false;
            photon_sat_badchi2[iparticle] = false;
        }
        for (int icluster = 0; icluster < ncluster; icluster++)
        {
            // need ET > 10 GeV
            if (cluster_Et[icluster] < 3)
                continue;

            if (particle_trkidmap.find(cluster_truthtrkID[icluster]) == particle_trkidmap.end())
            {
                // std::cout<<"trackid: "<<cluster_truthtrkID[icluster]<<std::endl;
                // std::cout << "Error: cluster_truthtrkID not found in particle_trkidmap" << std::endl;
                continue;
            }
            // find the particle with the same trkid
            int iparticle = particle_trkidmap[cluster_truthtrkID[icluster]];
            if (iparticle < 0 || iparticle >= nparticles)
            {
                std::cout << "Error: iparticle = " << iparticle << " out of range" << std::endl;
                continue;
            }
            float deta = cluster_Eta[icluster] - particle_Eta[iparticle];
            float dphi = cluster_Phi[icluster] - particle_Phi[iparticle];
            if (dphi > M_PI)
            {
                dphi = 2 * M_PI - dphi;
            }
            float dR = sqrt(deta * deta + dphi * dphi);
            if (dR < 0.1)
            {
                //if (cluster_Et[icluster] > 0.5 * particle_Pt[iparticle])
                if (cluster_Et[icluster] > 8.0)
                {
                    photon_reco[iparticle] = true;
                }
            }
            //check the status bit
            bool has_sat = false;
            bool has_badchi2 = false;
            for(int itower = 0; itower < 49; itower++)
            {   
                bool is_sat = get_status_bit(7, (uint8_t)cluster_status_array[icluster][itower]);
                bool is_badchi2 = get_status_bit(2, (uint8_t)cluster_status_array[icluster][itower]);

                float tower_E = cluster_e_array[icluster][itower];

                
                    eff_tower_sat_prob->Fill(is_sat, tower_E);
                    eff_tower_sat_badchi2_prob->Fill(is_badchi2&&is_sat, tower_E);
                if(is_badchi2&& tower_E>0)
                {
                    std::cout<<"bad chi2 tower E: "<<tower_E<<std::endl;
                }
                
                if(get_status_bit(7, (uint8_t)cluster_status_array[icluster][itower]))
                {
                    has_sat = true;
                    //std::cout<<"Tower "<<itower<<" has saturation , status: "<< cluster_status_array[icluster][itower]<<std::endl;
                }
                if(get_status_bit(2, (uint8_t)cluster_status_array[icluster][itower]))
                {
                    has_badchi2 = true;
                }
            }
            if(has_sat)
            {
                photon_sat[iparticle] = true;
            }
            if(has_badchi2 && has_sat)
            {
                photon_sat_badchi2[iparticle] = true;
            }


            
        }
         for (auto it = photon_reco.begin(); it != photon_reco.end(); ++it)
        {
            float photon_pT = particle_Pt[it->first];
            float photon_eta = particle_Eta[it->first];
            int etabin = -1;
            for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
            {
                if (photon_eta > eta_bins[ieta] && photon_eta < eta_bins[ieta + 1])
                {
                    etabin = ieta;
                    break;
                }
            }
            if (etabin == -1)
            {
                continue;
            }
            eff_sat_prob[etabin]->Fill(photon_sat[it->first], photon_pT);
            eff_sat_badchi2_prob[etabin]->Fill(photon_sat_badchi2[it->first], photon_pT);

            eff_reco_eff[etabin]->Fill(photon_reco[it->first], photon_pT);

        }
    }

    fout->cd();
    fout->Write();
    fout->Close();
}
