// SaturationStudy.C
//
// Produces:
//   1. Cluster ET distributions split by saturation status (nsaturated == 0 vs > 0)
//      for data, signal MC (photon5/10/20), and inclusive MC (jet10-50).
//   2. For signal MC only: 2D histogram of truth photon energy vs reco cluster energy
//      for truth-matched photon clusters that contain at least one saturated tower,
//      and separately for unsaturated clusters.
//
// Usage (ROOT interpreter):
//   root -l -q 'SaturationStudy.C("config_showershape.yaml","data")'
//   root -l -q 'SaturationStudy.C("config_showershape.yaml","photon20")'
//   root -l -q 'SaturationStudy.C("config_showershape.yaml","jet10")'

#include <yaml-cpp/yaml.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <cmath>
#include <map>
#include <limits>

void SaturationStudy(const std::string &configname = "config_showershape.yaml",
                     const std::string filetype = "data")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node cfg = YAML::LoadFile(configname);

    const bool issim   = (filetype != "data");
    const bool issignal = (filetype.find("photon") != std::string::npos);
    const bool isbackground = (filetype.find("jet") != std::string::npos);

    // -----------------------------------------------------------------------
    // Input file
    // -----------------------------------------------------------------------
    std::string infilename;
    if (!issim)
    {
        infilename = cfg["input"]["data_file"].as<std::string>();
    }
    else
    {
        infilename = cfg["input"]["photon_jet_file_root_dir"].as<std::string>()
                   + filetype
                   + cfg["input"]["photon_jet_file_branch_dir"].as<std::string>();
    }
    std::cout << "Input: " << infilename << std::endl;

    // -----------------------------------------------------------------------
    // Cross-section weights (all normalised to photon20 / jet50 respectively)
    // -----------------------------------------------------------------------
    const float photon5cross  = 146359.3f;
    const float photon10cross = 6944.675f;
    const float photon20cross = 130.4461f;
    const float jet10cross    = 3.997e+06f;
    const float jet15cross    = 4.073e+05f;
    const float jet20cross    = 6.218e+04f;
    const float jet30cross    = 2.502e+03f;
    const float jet50cross    = 7.2695f;

    float cross_weight = 1.0f;
    float max_photon_lower = 0, max_photon_upper = 200;
    float max_jet_lower    = 0, max_jet_upper    = 200;
    float cluster_ET_upper = 200;

    if      (filetype == "photon5")  { max_photon_lower =  0; max_photon_upper =  14; cross_weight = photon5cross  / photon20cross; }
    else if (filetype == "photon10") { max_photon_lower = 14; max_photon_upper =  30; cross_weight = photon10cross / photon20cross; }
    else if (filetype == "photon20") { max_photon_lower = 30; max_photon_upper = 200; cross_weight = 1.0f; }
    else if (filetype == "jet10")    { max_jet_lower = 10; max_jet_upper = 15;  cluster_ET_upper = 18;  cross_weight = jet10cross / jet50cross; }
    else if (filetype == "jet15")    { max_jet_lower = 15; max_jet_upper = 20;  cluster_ET_upper = 23;  cross_weight = jet15cross / jet50cross; }
    else if (filetype == "jet20")    { max_jet_lower = 20; max_jet_upper = 30;  cluster_ET_upper = 33;  cross_weight = jet20cross / jet50cross; }
    else if (filetype == "jet30")    { max_jet_lower = 30; max_jet_upper = 50;  cluster_ET_upper = 45;  cross_weight = jet30cross / jet50cross; }
    else if (filetype == "jet50")    { max_jet_lower = 50; max_jet_upper = 100;                         cross_weight = 1.0f; }

    // -----------------------------------------------------------------------
    // Vertex reweighting (simulation only)
    // -----------------------------------------------------------------------
    TH1 *h_vertex_reweight = nullptr;
    if (issim)
    {
        int vtx_on = cfg["analysis"]["vertex_reweight_on"].as<int>(1);
        if (vtx_on)
        {
            std::string vtxfile = cfg["analysis"]["vertex_reweight_file"].as<std::string>("results/vertex_reweight.root");
            TFile *fvtx = TFile::Open(vtxfile.c_str(), "READ");
            if (!fvtx || fvtx->IsZombie())
            {
                std::cerr << "[VertexReweight] Cannot open " << vtxfile << std::endl;
                return;
            }
            TH1 *htmp = dynamic_cast<TH1*>(fvtx->Get("h_vertexz_ratio_data_over_mccombined"));
            if (!htmp)
            {
                std::cerr << "[VertexReweight] Histogram not found in " << vtxfile << std::endl;
                return;
            }
            h_vertex_reweight = dynamic_cast<TH1*>(htmp->Clone("h_vtx_rw_clone"));
            h_vertex_reweight->SetDirectory(nullptr);
            fvtx->Close();
            std::cout << "[VertexReweight] Loaded from " << vtxfile << std::endl;
        }
    }

    // -----------------------------------------------------------------------
    // Analysis cuts from config
    // -----------------------------------------------------------------------
    float vertexcut   = cfg["analysis"]["vertex_cut"].as<float>();
    float reco_min_ET = cfg["analysis"]["reco_min_ET"].as<float>();
    float eff_dR      = cfg["analysis"]["eff_dR"].as<float>();

    std::vector<int> trigger_used;
    {
        YAML::Node trigNode = cfg["analysis"]["trigger_used"];
        if (trigNode && trigNode.IsSequence())
            trigger_used = trigNode.as<std::vector<int>>();
        else
            trigger_used.push_back(cfg["analysis"]["trigger_used"].as<int>());
    }

    std::string clusternodename = cfg["input"]["cluster_node_name"].as<std::string>();
    std::string var_type        = cfg["output"]["var_type"].as<std::string>();

    // -----------------------------------------------------------------------
    // Output file
    // -----------------------------------------------------------------------
    std::string outfilename;
    if (issim)
        outfilename = cfg["output"]["eff_outfile"].as<std::string>() + "_saturation_" + filetype + "_" + var_type + ".root";
    else
        outfilename = cfg["output"]["data_outfile"].as<std::string>() + "_saturation_" + var_type + ".root";
    std::cout << "Output: " << outfilename << std::endl;

    // -----------------------------------------------------------------------
    // TChain + TTreeReader
    // -----------------------------------------------------------------------
    std::string treename = cfg["input"]["tree"].as<std::string>();
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());

    TTreeReader reader(&chain);

    // Event variables
    TTreeReaderValue<int>   mbdnorthhit(reader, "mbdnorthhit");
    TTreeReaderValue<int>   mbdsouthhit(reader, "mbdsouthhit");
    TTreeReaderValue<int>   runnumber  (reader, "runnumber");
    TTreeReaderValue<float> vertexz    (reader, "vertexz");
    TTreeReaderArray<Bool_t> scaledtrigger(reader, "scaledtrigger");

    // Sim-only branches (truth particles + jet truth)
    TTreeReaderValue<int>   nparticles     (reader, "nparticles");
    TTreeReaderArray<float> particle_E     (reader, "particle_E");
    TTreeReaderArray<float> particle_Pt    (reader, "particle_Pt");
    TTreeReaderArray<float> particle_Eta   (reader, "particle_Eta");
    TTreeReaderArray<float> particle_Phi   (reader, "particle_Phi");
    TTreeReaderArray<int>   particle_pid   (reader, "particle_pid");
    TTreeReaderArray<int>   particle_trkid (reader, "particle_trkid");
    TTreeReaderArray<int>   particle_photonclass(reader, "particle_photonclass");
    TTreeReaderValue<int>   njet_truth     (reader, "njet_truth");
    TTreeReaderArray<float> jet_truth_Pt   (reader, "jet_truth_Pt");

    // Cluster branches
    TTreeReaderValue<int>   ncluster   (reader, Form("ncluster_%s",   clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Et (reader, Form("cluster_Et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_E  (reader, Form("cluster_E_%s",  clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s",clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s",clusternodename.c_str()));
    TTreeReaderArray<int>   cluster_nsaturated(reader, Form("cluster_nsaturated_%s", clusternodename.c_str()));
    TTreeReaderArray<int>   cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusternodename.c_str()));

    // -----------------------------------------------------------------------
    // Histograms
    // -----------------------------------------------------------------------
    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");

    // Cluster ET distributions (all clusters vs saturated clusters)
    const int  nET  = 90;
    const float ETlo = 0, EThi = 45;

    TH1D *h_ET_all       = new TH1D("h_ET_all",       "All clusters;Cluster E_{T} [GeV];Entries",       nET, ETlo, EThi);
    TH1D *h_ET_sat       = new TH1D("h_ET_sat",        "Saturated clusters (n_{sat}>0);Cluster E_{T} [GeV];Entries", nET, ETlo, EThi);
    TH1D *h_ET_nosat     = new TH1D("h_ET_nosat",      "Unsaturated clusters;Cluster E_{T} [GeV];Entries", nET, ETlo, EThi);
    TH1D *h_nsat_per_cluster = new TH1D("h_nsat_per_cluster", "Saturated towers per cluster;n_{sat};Entries", 30, 0, 30);

    // Fraction of clusters that are saturated vs ET (use TProfile or fill numerator+denominator)
    TH1D *h_sat_fraction_num   = new TH1D("h_sat_fraction_num",   "Sat fraction numerator;E_{T} [GeV]",   nET, ETlo, EThi);
    TH1D *h_sat_fraction_denom = new TH1D("h_sat_fraction_denom", "Sat fraction denominator;E_{T} [GeV]", nET, ETlo, EThi);

    // Signal MC only: 2D energy resolution histograms for truth-matched photon clusters
    // axes: truth photon E (x) vs reco cluster E (y)
    const int  nEbins = 100;
    const float Elo = 0, Ehi = 50;
    TH2D *h2_Etruth_Ereco_sat   = nullptr;
    TH2D *h2_Etruth_Ereco_nosat = nullptr;
    TH2D *h2_Etruth_Ereco_all   = nullptr;   // all truth-matched, regardless of saturation
    if (issignal)
    {
        h2_Etruth_Ereco_sat   = new TH2D("h2_Etruth_Ereco_sat",
            "Truth-matched photon, saturated cluster;E_{truth} [GeV];E_{reco} [GeV]",
            nEbins, Elo, Ehi, nEbins, Elo, Ehi);
        h2_Etruth_Ereco_nosat = new TH2D("h2_Etruth_Ereco_nosat",
            "Truth-matched photon, unsaturated cluster;E_{truth} [GeV];E_{reco} [GeV]",
            nEbins, Elo, Ehi, nEbins, Elo, Ehi);
        h2_Etruth_Ereco_all   = new TH2D("h2_Etruth_Ereco_all",
            "Truth-matched photon, all clusters;E_{truth} [GeV];E_{reco} [GeV]",
            nEbins, Elo, Ehi, nEbins, Elo, Ehi);
        // Sumw2 required so GetBinError returns sqrt(sum_w_i^2),
        // enabling correct N_eff = sum_w^2 / sum_w2 in plotting code.
        h2_Etruth_Ereco_sat->Sumw2();
        h2_Etruth_Ereco_nosat->Sumw2();
        h2_Etruth_Ereco_all->Sumw2();
    }

    // Vertex z for monitoring
    TH1D *h_vertexz = new TH1D("h_vertexz", "Vertex z;z [cm];Entries", 200, -100, 100);

    h_ET_all->Sumw2();
    h_ET_sat->Sumw2();
    h_ET_nosat->Sumw2();
    h_sat_fraction_num->Sumw2();
    h_sat_fraction_denom->Sumw2();

    // -----------------------------------------------------------------------
    // Event loop
    // -----------------------------------------------------------------------
    const long long nentries = chain.GetEntries();
    long long ientry = 0;
    std::cout << "Total entries: " << nentries << std::endl;

    while (reader.Next())
    {
        if (ientry % 100000 == 0)
            std::cout << "Processing " << ientry << " / " << nentries << std::endl;
        ientry++;

        // Per-event weight (will be updated by vertex reweighting for MC)
        float weight = cross_weight;

        // -- Data: trigger selection --
        if (!issim)
        {
            bool any_fired = false;
            for (int itrig : trigger_used)
            {
                if (itrig >= 0 && (unsigned int)itrig < scaledtrigger.GetSize() && scaledtrigger[itrig])
                {
                    any_fired = true;
                    break;
                }
            }
            if (!any_fired) continue;
        }

        // -- Sim: pT-hat range cut + vertex reweighting --
        if (issim)
        {
            // pT-hat range cut
            if (issignal)
            {
                float maxphotonpT = 0;
                for (int ip = 0; ip < *nparticles; ip++)
                {
                    if (particle_pid[ip] == 22 && particle_Pt[ip] > maxphotonpT)
                        maxphotonpT = particle_Pt[ip];
                }
                if (maxphotonpT < max_photon_lower || maxphotonpT > max_photon_upper) continue;
            }
            else if (isbackground)
            {
                float maxjetpT = 0;
                for (int ij = 0; ij < *njet_truth; ij++)
                {
                    if (jet_truth_Pt[ij] > maxjetpT) maxjetpT = jet_truth_Pt[ij];
                }
                if (maxjetpT < max_jet_lower || maxjetpT > max_jet_upper) continue;
            }

            // vertex reweighting
            if (h_vertex_reweight)
            {
                int bin = h_vertex_reweight->FindBin(*vertexz);
                if (bin < 1) bin = 1;
                if (bin > h_vertex_reweight->GetNbinsX()) bin = h_vertex_reweight->GetNbinsX();
                float vw = h_vertex_reweight->GetBinContent(bin);
                if (std::isfinite(vw) && vw > 0) weight *= vw;
            }
        }

        // -- Common event selection --
        if (std::fabs(*vertexz) > vertexcut) continue;
        if (!(*mbdnorthhit >= 1 && *mbdsouthhit >= 1)) continue;

        h_vertexz->Fill(*vertexz, weight);

        // -- Build truth-trkid map for signal MC --
        std::map<int, int> trkidmap;   // trkid -> particle index
        if (issignal)
        {
            for (int ip = 0; ip < *nparticles; ip++)
                trkidmap[particle_trkid[ip]] = ip;
        }

        // -- Cluster loop --
        for (int ic = 0; ic < *ncluster; ic++)
        {
            float clET = cluster_Et[ic];
            float clE  = cluster_E[ic];
            if (clET < reco_min_ET) continue;
            if (isbackground && clET > cluster_ET_upper) continue;

            int nsat = cluster_nsaturated[ic];
            bool is_saturated = (nsat > 0);

            h_ET_all->Fill(clET, weight);
            h_nsat_per_cluster->Fill(nsat, weight);
            h_sat_fraction_denom->Fill(clET, weight);

            if (is_saturated)
            {
                h_ET_sat->Fill(clET, weight);
                h_sat_fraction_num->Fill(clET, weight);
            }
            else
            {
                h_ET_nosat->Fill(clET, weight);
            }

            // -- Signal MC: truth matching for energy resolution --
            if (issignal)
            {
                int trkid = cluster_truthtrkID[ic];
                auto it = trkidmap.find(trkid);
                if (it == trkidmap.end()) continue;
                int ip = it->second;

                // pid must be photon
                if (particle_pid[ip] != 22) continue;

                // dR matching cut
                float deta = cluster_Eta[ic] - particle_Eta[ip];
                float dphi = cluster_Phi[ic] - particle_Phi[ip];
                while (dphi >  M_PI) dphi -= 2*M_PI;
                while (dphi < -M_PI) dphi += 2*M_PI;
                float dR = std::sqrt(deta*deta + dphi*dphi);
                if (dR > eff_dR) continue;

                float truthE = particle_E[ip];

                h2_Etruth_Ereco_all->Fill(truthE, clE, weight);
                if (is_saturated)
                    h2_Etruth_Ereco_sat->Fill(truthE, clE, weight);
                else
                    h2_Etruth_Ereco_nosat->Fill(truthE, clE, weight);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Write and close
    // -----------------------------------------------------------------------
    fout->cd();
    h_ET_all->Write();
    h_ET_sat->Write();
    h_ET_nosat->Write();
    h_nsat_per_cluster->Write();
    h_sat_fraction_num->Write();
    h_sat_fraction_denom->Write();
    h_vertexz->Write();

    if (issignal)
    {
        h2_Etruth_Ereco_all->Write();
        h2_Etruth_Ereco_sat->Write();
        h2_Etruth_Ereco_nosat->Write();
    }

    fout->Close();
    std::cout << "Written to " << outfilename << std::endl;
}
