#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <iostream>
#include <vector>

void slimtreeana() {
    bool issim = false;
    // Open the ROOT file
    ///TFile *file = TFile::Open("/sphenix/user/shuhangli/pi0pythiasim/macro/condorphotonjet10cnn/caloana.root", "READ");
    //TFile *file = TFile::Open("/sphenix/user/shuhangli/pi0pythiasim/macro/condorincludivejet10cnn/caloananewmodel.root", "READ");
    TFile *file = TFile::Open("/sphenix/user/shuhangli/2024calofastana/macro/condor47-48/condorphotontest/condorout/treeall.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Get the tree from the file
    TTree *tree = (TTree*)file->Get("slimtree");
    if (!tree) {
        std::cerr << "Tree not found!" << std::endl;
        return;
    }

    std::vector<std::string> clusternodenames = 
    {
    "CLUSTERINFO_CEMC", 
    "CLUSTERINFO_CEMC_CNN", 
    "CLUSTERINFO_CEMC_CNN_single", 
    "CLUSTERINFO_CEMC_CNN_noise",
    "CLUSTERINFO_CEMC_CNN_single_noise"
    };
    //use 4GeV for the truth iso cut now
    const int truthisocut = 4;
    const int nnodes = (int)clusternodenames.size();
    const int maxclusters = 10000;
    float vertexz;
    int scaledtrigger[32];
    int ncluster[nnodes];
    float cluster_E[nnodes][maxclusters];
    float cluster_Et[nnodes][maxclusters];
    float cluster_Eta[nnodes][maxclusters];
    float cluster_Phi[nnodes][maxclusters];
    float cluster_prob[nnodes][maxclusters];
    float cluster_iso_04[nnodes][maxclusters];
    int cluster_truthtrkID[nnodes][maxclusters];

    int nparticles;
    const int maxparticles = 10000;
    float particle_E[maxparticles];
    float particle_Pt[maxparticles];
    float particle_Eta[maxparticles];
    float particle_Phi[maxparticles];
    int particle_pid[maxparticles];
    int particle_trkid[maxparticles];
    int particle_isprompt_photon[maxparticles];
    float particle_truth_iso_04[maxparticles];
    int particle_converted[maxparticles];
    
    
    tree->SetBranchAddress("vertexz", &vertexz);
    if(!issim){
        tree->SetBranchAddress("scaledtrigger", scaledtrigger);
    }
    
    if(issim)
    {
        tree->SetBranchAddress("nparticles", &nparticles);
        tree->SetBranchAddress("particle_E", particle_E);
        tree->SetBranchAddress("particle_Pt", particle_Pt);
        tree->SetBranchAddress("particle_Eta", particle_Eta);
        tree->SetBranchAddress("particle_Phi", particle_Phi);
        tree->SetBranchAddress("particle_pid", particle_pid);
        tree->SetBranchAddress("particle_trkid", particle_trkid);
        tree->SetBranchAddress("particle_isprompt_photon", particle_isprompt_photon);
        tree->SetBranchAddress("particle_truth_iso_04", particle_truth_iso_04);
        tree->SetBranchAddress("particle_converted", particle_converted);

    }
    //struct for particle
    struct particle
    {
        float E;
        float Pt;
        float Eta;
        float Phi;
        int pid;
        int trkid;
        int isprompt_photon;
        float truth_iso_04;
        int converted;
    };
    

    for(int i = 0; i < nnodes; i++)
    {
        tree->SetBranchAddress(Form("ncluster_%s", clusternodenames[i].c_str()), &ncluster[i]);
        tree->SetBranchAddress(Form("cluster_E_%s", clusternodenames[i].c_str()), cluster_E[i]);
        tree->SetBranchAddress(Form("cluster_Et_%s", clusternodenames[i].c_str()), cluster_Et[i]);
        tree->SetBranchAddress(Form("cluster_Eta_%s", clusternodenames[i].c_str()), cluster_Eta[i]);
        tree->SetBranchAddress(Form("cluster_Phi_%s", clusternodenames[i].c_str()), cluster_Phi[i]);
        tree->SetBranchAddress(Form("cluster_prob_%s", clusternodenames[i].c_str()), cluster_prob[i]);
        tree->SetBranchAddress(Form("cluster_iso_04_%s", clusternodenames[i].c_str()), cluster_iso_04[i]);
        if(issim)tree->SetBranchAddress(Form("cluster_truthtrkID_%s", clusternodenames[i].c_str()), cluster_truthtrkID[i]);
    }
    
    TFile* fhistoout[nnodes];
    const int nprobbins = 10;
    TH2D* h_cluster_iso_04_Et[nnodes][nprobbins];
    TH2D* h_prob_ET[nnodes];
    TH2D* h_signal_cluster_iso_04_Et[nnodes][nprobbins];
    TH2D* h_nonphoton_cluster_iso_04_Et[nnodes][nprobbins];
    TH2D* h_nonisophoton_cluster_iso_04_E[nnodes][nprobbins];
    TH2D* h_sim_prob_ET[nnodes];
    for (int i = 0; i < nnodes; i++) {
      fhistoout[i] = new TFile(
          //Form("/sphenix/user/shuhangli/pi0pythiasim/macro/isophotonsimhisto/inclusive_photon_histo_%s.root", clusternodenames[i].c_str()), "RECREATE");
          Form("/sphenix/user/shuhangli/pi0pythiasim/macro/isophotondatahisto/photon_histo_%s.root", clusternodenames[i].c_str()), "RECREATE");

        h_prob_ET[i] = new TH2D(
                "h_prob_ET",
                //Form("h_prob_ET_%s_%d", clusternodenames[i].c_str(), j),
                "",
                100, 0, 1, 100, 0, 50);

      for (int j = 0; j < nprobbins; j++) {

          h_cluster_iso_04_Et[i][j] = new TH2D(
              Form("h_cluster_iso_04_Et_%d", j),
              //Form("h_cluster_iso_04_Et_%s_%d", clusternodenames[i].c_str(), j),
              "",
              55, -15, 40, 100, 0, 50);
            
        
        if (issim) 
        {
            h_signal_cluster_iso_04_Et[i][j] = new TH2D(
                Form("h_signal_cluster_iso_04_Et_%d", j),
                //Form("h_signal_cluster_iso_04_Et_%s_%d", clusternodenames[i].c_str(), j),
                "",
                55, -15, 40, 100, 0, 50);
            h_nonphoton_cluster_iso_04_Et[i][j] = new TH2D(
                Form("h_nonphoton_cluster_iso_04_Et_%d", j),
                //Form("h_nonphoton_cluster_iso_04_Et_%s_%d", clusternodenames[i].c_str(), j),
                "",
                55, -15, 40, 100, 0, 50);
            h_nonisophoton_cluster_iso_04_E[i][j] = new TH2D(
                Form("h_nonisophoton_cluster_iso_04_E_%d", j),
                //Form("h_nonisophoton_cluster_iso_04_E_%s_%d", clusternodenames[i].c_str(), j),
                "",
                55, -15, 40, 100, 0, 50);

        }

      }
    }

    // Loop over all events in the tree
    Long64_t nentries = tree->GetEntries();
    for (Long64_t ie = 0; ie < nentries; ie++) {
        if(ie%1000==0)std::cout<<ie<<"/"<<nentries<<std::endl;
        if(ie>nentries/2) continue;
        tree->GetEntry(ie);
        if(!issim)
        {
            if(scaledtrigger[26] == 0) continue;
        }
        //particle loop
        //map from trkid to particle
        std::map<int, particle> truth_particles;
        if(issim)
        {
            for (int j = 0; j < nparticles; j++) {
                particle p;
                p.E = particle_E[j];
                p.Pt = particle_Pt[j];
                p.Eta = particle_Eta[j];
                p.Phi = particle_Phi[j];
                p.pid = particle_pid[j];
                p.trkid = particle_trkid[j];
                p.isprompt_photon = particle_isprompt_photon[j];
                p.truth_iso_04 = particle_truth_iso_04[j];
                p.converted = particle_converted[j];
                truth_particles[p.trkid] = p;
            }
        }
        //cluster loop
        for(int i = 0; i < nnodes; i++)
        {
            //check vertex <30
            if (abs(vertexz) > 60) continue;
            for (int j = 0; j < ncluster[i]; j++) {
                //check cluster Eta
                if (abs(cluster_Eta[i][j]) >0.6) continue;
                float clusterprob = cluster_prob[i][j];
                int clusterprobbin = (int)(clusterprob * nprobbins);
                if (clusterprobbin >= nprobbins) clusterprobbin = nprobbins - 1;
                h_cluster_iso_04_Et[i][clusterprobbin]->Fill(cluster_iso_04[i][j], cluster_Et[i][j]);
                h_prob_ET[i]->Fill(clusterprob, cluster_Et[i][j]);
                if(issim)
                {
                    //check if not found
                    if (truth_particles.find(cluster_truthtrkID[i][j]) == truth_particles.end()) {
                        //std::cerr << "Truth particle not found!" << std::endl;
                        continue;
                    }
                    //find the truth particle in the map
                    particle p = truth_particles[cluster_truthtrkID[i][j]];
                    float particleiso = p.truth_iso_04;
                    //check if the particle is a prompt photon
                    if (p.pid == 22 && p.isprompt_photon == 1) 
                    {
                        if(p.converted == 1)continue;
                        if(particleiso < truthisocut)
                        {
                            h_signal_cluster_iso_04_Et[i][clusterprobbin]->Fill(cluster_iso_04[i][j], cluster_Et[i][j]);
                        }
                        else
                        {
                            h_nonisophoton_cluster_iso_04_E[i][clusterprobbin]->Fill(cluster_iso_04[i][j], cluster_E[i][j]);
                        }
                    }
                    else
                    {
                        h_nonphoton_cluster_iso_04_Et[i][clusterprobbin]->Fill(cluster_iso_04[i][j], cluster_Et[i][j]);
                    }
                }
            }
        }


        
    }

    for (int i = 0; i < nnodes; i++) {
      fhistoout[i]->Write();
      fhistoout[i]->Close();
    }
}
