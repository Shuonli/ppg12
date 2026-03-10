// BDTScoreVsET.C
// Creates 2D histograms of BDT score vs. reco cluster ET:
//   h_incl_bdt_<model>_<ieta>   -- inclusive (jet MC, NPB-passed clusters)
//   h_signal_bdt_<model>_<ieta> -- signal (photon MC, truth-matched direct/frag photons)
//
// Usage: root -l -q 'BDTScoreVsET.C("config_bdt_nom.yaml")'

#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <vector>
#include <map>
#include <cmath>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TSystem.h>
#include <TObjString.h>
#include <yaml-cpp/yaml.h>

// Hardcoded BDT model names — edit this list to add more models.
static const std::vector<std::string> BDT_MODEL_NAMES = {"base_v3E", "base", "base_E", "base_v1E", "base_v2E"};

void SaveYamlToRoot_BDT(TFile *f, const char *yaml_filename)
{
    std::ifstream yaml_file(yaml_filename);
    std::string yaml_content((std::istreambuf_iterator<char>(yaml_file)),
                             std::istreambuf_iterator<char>());
    TObjString yaml_obj(yaml_content.c_str());
    f->cd();
    yaml_obj.Write("config");
}

void BDTScoreVsET(const std::string &configname = "config_bdt_nom.yaml")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    // -----------------------------------------------------------------------
    // Read config
    // -----------------------------------------------------------------------
    std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();
    std::string treename = configYaml["input"]["tree"].as<std::string>();
    std::string var_type = configYaml["output"]["var_type"].as<std::string>();

    float vertexcut = configYaml["analysis"]["vertex_cut"].as<float>();
    std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();
    int n_eta_bins = (int)eta_bins.size() - 1;
    float reco_min_ET = configYaml["analysis"]["reco_min_ET"].as<float>();
    float eff_dR = configYaml["analysis"]["eff_dR"].as<float>();
    float truthisocut = configYaml["analysis"]["truth_iso_max"].as<float>();
    int conesize = configYaml["analysis"]["cone_size"].as<int>();

    int common_npb_cut_on = configYaml["analysis"]["common"]["npb_cut_on"].as<int>(0);
    float common_npb_score_cut = configYaml["analysis"]["common"]["npb_score_cut"].as<float>(0.5);

    // -----------------------------------------------------------------------
    // Cross sections
    // -----------------------------------------------------------------------
    // photon samples normalised to photon20; jet samples normalised to jet50
    const float photon5cross  = 146359.3;
    const float photon10cross = 6944.675;
    const float photon20cross = 130.4461;
    const float jet10cross = 3.997e+06;
    const float jet15cross = 4.073e+05;
    const float jet20cross = 6.218e+04;
    const float jet30cross = 2.502e+03;
    const float jet50cross = 7.2695;

    // -----------------------------------------------------------------------
    // Vertex reweighting (applied to all sim samples)
    // -----------------------------------------------------------------------
    TH1 *h_vertex_reweight = nullptr;
    int vertex_reweight_on = configYaml["analysis"]["vertex_reweight_on"].as<int>(1);
    std::string vertex_reweight_file =
        configYaml["analysis"]["vertex_reweight_file"].as<std::string>("results/vertex_reweight.root");

    if (vertex_reweight_on)
    {
        TFile *fvtx = TFile::Open(vertex_reweight_file.c_str(), "READ");
        if (!fvtx || fvtx->IsZombie())
        {
            std::cerr << "[VertexReweight] ERROR: cannot open " << vertex_reweight_file << std::endl;
            return;
        }
        TH1 *htmp = dynamic_cast<TH1*>(fvtx->Get("h_vertexz_ratio_data_over_mccombined"));
        if (!htmp)
        {
            std::cerr << "[VertexReweight] ERROR: histogram not found in " << vertex_reweight_file << std::endl;
            fvtx->Close();
            return;
        }
        h_vertex_reweight = dynamic_cast<TH1*>(htmp->Clone("h_vtx_reweight_clone"));
        h_vertex_reweight->SetDirectory(nullptr);
        std::cout << "[VertexReweight] Loaded from " << vertex_reweight_file << std::endl;
    }

    // -----------------------------------------------------------------------
    // Create output histograms
    //   X axis: cluster ET  (60 bins, 5–35 GeV)
    //   Y axis: BDT score   (100 bins, 0–1)
    // -----------------------------------------------------------------------
    std::map<std::string, std::vector<TH2D*>> h_incl_bdt;
    std::map<std::string, std::vector<TH2D*>> h_signal_bdt;

    for (const auto &model : BDT_MODEL_NAMES)
    {
        for (int ieta = 0; ieta < n_eta_bins; ieta++)
        {
            h_incl_bdt[model].push_back(new TH2D(
                Form("h_incl_bdt_%s_%d", model.c_str(), ieta),
                Form("Inclusive BDT vs E_{T} [%s], %.1f < #eta < %.1f;"
                     " Cluster E_{T} [GeV]; BDT Score",
                     model.c_str(), eta_bins[ieta], eta_bins[ieta + 1]),
                60, 5, 35, 1000, 0, 1));

            h_signal_bdt[model].push_back(new TH2D(
                Form("h_signal_bdt_%s_%d", model.c_str(), ieta),
                Form("Signal BDT vs E_{T} [%s], %.1f < #eta < %.1f;"
                     " Cluster E_{T} [GeV]; BDT Score",
                     model.c_str(), eta_bins[ieta], eta_bins[ieta + 1]),
                60, 5, 35, 1000, 0, 1));
        }
    }

    // -----------------------------------------------------------------------
    // Process all sample types
    // -----------------------------------------------------------------------
    std::vector<std::string> all_filetypes = {
        "photon5", "photon10", "photon20",
        "jet10",   "jet15",    "jet20",   "jet30",   "jet50"
    };

    for (const auto &filetype : all_filetypes)
    {
        std::cout << "\n=== Processing sample: " << filetype << " ===" << std::endl;

        bool isbackground = (filetype.find("jet") != std::string::npos);

        // Per-sample weights and pT-hat window
        float cross_weight = 1.0;
        float max_photon_lower = 0, max_photon_upper = 100;
        float max_jet_lower    = 0, max_jet_upper    = 100;
        float cluster_ET_upper = 100;

        if      (filetype == "photon5")  { max_photon_lower =  0; max_photon_upper =  14; cross_weight = photon5cross  / photon20cross; }
        else if (filetype == "photon10") { max_photon_lower = 14; max_photon_upper =  30; cross_weight = photon10cross / photon20cross; }
        else if (filetype == "photon20") { max_photon_lower = 30; max_photon_upper = 200; cross_weight = 1.0; }
        else if (filetype == "jet10")    { max_jet_lower = 10; max_jet_upper =  15; cluster_ET_upper = 18; cross_weight = jet10cross / jet50cross; }
        else if (filetype == "jet15")    { max_jet_lower = 15; max_jet_upper =  20; cluster_ET_upper = 23; cross_weight = jet15cross / jet50cross; }
        else if (filetype == "jet20")    { max_jet_lower = 20; max_jet_upper =  30; cluster_ET_upper = 33; cross_weight = jet20cross / jet50cross; }
        else if (filetype == "jet30")    { max_jet_lower = 30; max_jet_upper =  50; cluster_ET_upper = 45; cross_weight = jet30cross / jet50cross; }
        else if (filetype == "jet50")    { max_jet_lower = 50; max_jet_upper = 100;                        cross_weight = 1.0; }

        // Build input chain
        std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;
        TChain chain(treename.c_str());
        chain.Add(infilename.c_str());
        std::cout << "Input: " << infilename << std::endl;
        std::cout << "Entries: " << chain.GetEntries() << std::endl;

        // -------------------------------------------------------------------
        // TTreeReader setup
        // -------------------------------------------------------------------
        TTreeReader reader(&chain);

        // Event-level branches
        TTreeReaderValue<int>   mbdnorthhit(reader, "mbdnorthhit");
        TTreeReaderValue<int>   mbdsouthhit(reader, "mbdsouthhit");
        TTreeReaderValue<int>   nparticles (reader, "nparticles");
        TTreeReaderValue<int>   ncluster   (reader, Form("ncluster_%s", clusternodename.c_str()));
        TTreeReaderValue<float> vertexz    (reader, "vertexz");

        // Truth-jet branches (for jet pT-hat cut)
        TTreeReaderValue<int>      njet_truth   (reader, "njet_truth");
        TTreeReaderArray<float>    jet_truth_Pt (reader, "jet_truth_Pt");

        // Particle branches (for photon truth matching)
        TTreeReaderArray<float> particle_Pt           (reader, "particle_Pt");
        TTreeReaderArray<float> particle_Eta          (reader, "particle_Eta");
        TTreeReaderArray<float> particle_Phi          (reader, "particle_Phi");
        TTreeReaderArray<float> particle_truth_iso_02 (reader, "particle_truth_iso_02");
        TTreeReaderArray<float> particle_truth_iso_03 (reader, "particle_truth_iso_03");
        TTreeReaderArray<float> particle_truth_iso_04 (reader, "particle_truth_iso_04");
        TTreeReaderArray<int>   particle_pid          (reader, "particle_pid");
        TTreeReaderArray<int>   particle_trkid        (reader, "particle_trkid");
        TTreeReaderArray<int>   particle_photonclass  (reader, "particle_photonclass");

        // Cluster branches
        TTreeReaderArray<float> cluster_Et       (reader, Form("cluster_Et_%s",       clusternodename.c_str()));
        TTreeReaderArray<float> cluster_Eta      (reader, Form("cluster_Eta_%s",      clusternodename.c_str()));
        TTreeReaderArray<float> cluster_Phi      (reader, Form("cluster_Phi_%s",      clusternodename.c_str()));
        TTreeReaderArray<int>   cluster_truthtrkID (reader, Form("cluster_truthtrkID_%s", clusternodename.c_str()));
        TTreeReaderArray<float> cluster_npb_score (reader, Form("cluster_npb_score_%s", clusternodename.c_str()));

        // BDT score branch for each model (allocated on heap so they stay in scope)
        std::vector<TTreeReaderArray<float>*> bdt_readers;
        for (const auto &model : BDT_MODEL_NAMES)
        {
            bdt_readers.push_back(new TTreeReaderArray<float>(
                reader,
                Form("cluster_bdt_%s_%s", clusternodename.c_str(), model.c_str())));
        }

        // -------------------------------------------------------------------
        // Event loop
        // -------------------------------------------------------------------
        int nentries = (int)chain.GetEntries();
        int ientry   = 0;

        while (reader.Next())
        {
            if (ientry % 10000 == 0)
                std::cout << "  Entry " << ientry << " / " << nentries << std::endl;

            // Per-event weight starts with cross-section weight
            float weight = cross_weight;

            // Vertex reweighting
            if (vertex_reweight_on && h_vertex_reweight)
            {
                int bin = h_vertex_reweight->FindBin(*vertexz);
                if (bin < 1) bin = 1;
                if (bin > h_vertex_reweight->GetNbinsX()) bin = h_vertex_reweight->GetNbinsX();
                float vw = h_vertex_reweight->GetBinContent(bin);
                if (std::isfinite(vw) && vw > 0.0f) weight *= vw;
            }

            // ----------------------------------------------------------
            // pT-hat window cut + build truth-matching map (photon samples)
            // ----------------------------------------------------------
            std::map<int, int>  particle_trkidmap; // trkid -> particle index
            std::map<int, bool> photon_signal;      // particle index -> true if signal photon

            if (!isbackground)
            {
                // Photon sample: find max photon pT and build signal map
                float maxphotonpT = 0;
                for (int ip = 0; ip < *nparticles; ip++)
                {
                    particle_trkidmap[particle_trkid[ip]] = ip;

                    if (particle_pid[ip] == 22 && particle_Pt[ip] > maxphotonpT)
                        maxphotonpT = particle_Pt[ip];

                    // Signal: direct or fragmentation photon with truth iso < cut
                    if (particle_pid[ip] == 22 && particle_photonclass[ip] < 3)
                    {
                        float truthisoET = 0;
                        if      (conesize == 4) truthisoET = particle_truth_iso_04[ip];
                        else if (conesize == 3) truthisoET = particle_truth_iso_03[ip];
                        else if (conesize == 2) truthisoET = particle_truth_iso_02[ip];
                        if (truthisoET < truthisocut)
                            photon_signal[ip] = true;
                    }
                }

                // pT-hat cut
                if (maxphotonpT > max_photon_upper || maxphotonpT < max_photon_lower)
                { ientry++; continue; }
            }
            else
            {
                // Jet sample: find max truth jet pT for pT-hat cut
                float maxjetpT = 0;
                for (int ijet = 0; ijet < *njet_truth; ijet++)
                {
                    if (jet_truth_Pt[ijet] > maxjetpT) maxjetpT = jet_truth_Pt[ijet];
                }
                if (maxjetpT == 0)               { ientry++; continue; }
                if (maxjetpT > max_jet_upper)    { ientry++; continue; }
                if (maxjetpT < max_jet_lower)    { ientry++; continue; }
            }

            // Vertex cut
            if (std::abs(*vertexz) > vertexcut) { ientry++; continue; }

            // MBD cut
            if (!(*mbdnorthhit >= 1 && *mbdsouthhit >= 1)) { ientry++; continue; }

            // ----------------------------------------------------------
            // Cluster loop
            // ----------------------------------------------------------
            for (int ic = 0; ic < *ncluster; ic++)
            {
                float clusterET  = cluster_Et[ic];
                float cluster_eta = cluster_Eta[ic];

                if (clusterET < reco_min_ET) continue;
                if (isbackground && clusterET > cluster_ET_upper) continue;

                // Eta bin
                int etabin = -1;
                for (int ieta = 0; ieta < n_eta_bins; ieta++)
                {
                    if (cluster_eta > eta_bins[ieta] && cluster_eta < eta_bins[ieta + 1])
                    { etabin = ieta; break; }
                }
                if (etabin == -1) continue;

                // NPB cut (applied to both inclusive and signal)
                if (common_npb_cut_on && cluster_npb_score[ic] <= common_npb_score_cut)
                    continue;

                if (isbackground)
                {
                    // Inclusive fill (jet MC)
                    for (int im = 0; im < (int)BDT_MODEL_NAMES.size(); im++)
                    {
                        float bdt_score = (*bdt_readers[im])[ic];
                        h_incl_bdt[BDT_MODEL_NAMES[im]][etabin]->Fill(clusterET, bdt_score, weight);
                    }
                }
                else
                {
                    // Signal fill (photon MC, truth-matched)
                    auto it_trk = particle_trkidmap.find(cluster_truthtrkID[ic]);
                    if (it_trk == particle_trkidmap.end()) continue;
                    int iparticle = it_trk->second;
                    if (photon_signal.find(iparticle) == photon_signal.end()) continue;

                    // dR cut between cluster and truth photon
                    float deta = cluster_Eta[ic]  - particle_Eta[iparticle];
                    float dphi = cluster_Phi[ic]  - particle_Phi[iparticle];
                    while (dphi >  M_PI) dphi -= 2 * M_PI;
                    while (dphi < -M_PI) dphi += 2 * M_PI;
                    float dR = std::sqrt(deta * deta + dphi * dphi);
                    if (dR >= eff_dR) continue;

                    for (int im = 0; im < (int)BDT_MODEL_NAMES.size(); im++)
                    {
                        float bdt_score = (*bdt_readers[im])[ic];
                        h_signal_bdt[BDT_MODEL_NAMES[im]][etabin]->Fill(clusterET, bdt_score, weight);
                    }
                }
            } // end cluster loop

            ientry++;
        } // end event loop

        // Clean up per-sample BDT reader pointers
        for (auto *ptr : bdt_readers) delete ptr;
        bdt_readers.clear();

    } // end filetype loop

    // -----------------------------------------------------------------------
    // Write output
    // -----------------------------------------------------------------------
    std::string outfilename = "results/BDTScoreVsET_" + var_type + ".root";
    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");

    for (const auto &model : BDT_MODEL_NAMES)
    {
        for (int ieta = 0; ieta < n_eta_bins; ieta++)
        {
            h_incl_bdt[model][ieta]->Write();
            h_signal_bdt[model][ieta]->Write();
        }
    }

    SaveYamlToRoot_BDT(fout, configname.c_str());
    fout->Close();

    std::cout << "\nOutput written to: " << outfilename << std::endl;
}
