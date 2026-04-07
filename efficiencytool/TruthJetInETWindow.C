// TruthJetInETWindow.C
// For clusters in each of a list of reco ET windows, fills the distribution of
// truth jet pT in the same event, split into signal (truth-matched direct+frag
// photons with isoETtruth < 4 GeV) and background categories.
// All ET windows are filled in a single pass over the data.
// Processes one MC sample at a time (use hadd to combine).
//
// Output ROOT file: results/truthjet_ETwindows_<filetype>.root
//   h_maxtruthjet_pt_{sig,bkg,all}_win{i}  - max truth jet pT per event with a qualifying cluster
//   h_alltruthjet_pt_{sig,bkg,all}_win{i}  - pT of every truth jet (|eta|<1.1) in those events
//   h2_clusterET_maxjet                    - 2D: cluster ET vs max truth jet pT (full ET range)
//   h_ET_edges                             - stores window edge values for the plotting macro
//
// Plotting: hadd all per-sample files, then run PlotTruthJetWindows.C.
//
// Usage:
//   root -l -b -q 'TruthJetInETWindow.C("config_bdt_nom.yaml","photon5")'

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TSystem.h>
#include <yaml-cpp/yaml.h>
#include "CrossSectionWeights.h"
using namespace PPG12;

// ── per-sample parameters ─────────────────────────────────────────────────────
struct SampleParams
{
    std::string filetype;
    bool   isbackground;
    float  max_photon_lower;
    float  max_photon_upper;
    float  max_jet_lower;
    float  max_jet_upper;
    float  cluster_ET_upper;
    float  cross_weight;   // relative to jet50 (bkg) or photon20 (sig)
};

static std::vector<SampleParams> BuildSampleList()
{
    return {
        // photon samples (signal)
        {"photon5",  false,  0,  14,  0, 100, 100, photon5cross  / photon20cross},
        {"photon10", false, 14,  30,  0, 100, 100, photon10cross / photon20cross},
        {"photon20", false, 30, 200,  0, 100, 100, 1.0f},
        // jet samples (background)
        {"jet5",     true,   0, 200,  7,  14,  14, jet5cross  / jet50cross},
        {"jet12",    true,   0, 200, 14,  21,  23, jet12cross / jet50cross},
        {"jet20",    true,   0, 200, 21,  32,  34, jet20cross / jet50cross},
        {"jet30",    true,   0, 200, 32,  42,  44, jet30cross / jet50cross},
        {"jet40",    true,   0, 200, 42, 100, 100, jet40cross / jet50cross},
    };
}

// ── main macro ────────────────────────────────────────────────────────────────
void TruthJetInETWindow(const std::string &configname = "config_isoroc.yaml",
                        const std::string &filetype   = "photon5")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node cfg = YAML::LoadFile(configname);

    // ── config ────────────────────────────────────────────────────────────────
    const std::string root_dir      = cfg["input"]["photon_jet_file_root_dir"].as<std::string>();
    const std::string branch_dir    = cfg["input"]["photon_jet_file_branch_dir"].as<std::string>();
    const std::string treename      = cfg["input"]["tree"].as<std::string>();
    const std::string clusternodename = cfg["input"]["cluster_node_name"].as<std::string>();
    const float vertexcut           = cfg["analysis"]["vertex_cut"].as<float>(60.f);
    const float truthisocut         = cfg["analysis"]["truth_iso_max"].as<float>(4.f);
    const float eff_dR              = cfg["analysis"]["eff_dR"].as<float>(0.1f);
    const float eta_max_cluster     = cfg["analysis"]["eta_bins"].as<std::vector<float>>().back();
    const float eta_max_jet         = 1.1f;   // truth jet eta acceptance for filling

    // ── tight + iso selection config (only used if bdt_model_name is set) ────
    const std::string bdt_model_name = cfg["input"]["bdt_model_name"].as<std::string>("");
    const float tight_bdt_min  = cfg["analysis"]["tight"]["bdt_min"].as<float>(0.6f);
    const float tight_bdt_max  = cfg["analysis"]["tight"]["bdt_max"].as<float>(1.0f);
    const float reco_iso_max_b = cfg["analysis"]["reco_iso_max_b"].as<float>(0.f);
    const float reco_iso_max_s = cfg["analysis"]["reco_iso_max_s"].as<float>(0.f);
    const int   npb_cut_on     = cfg["analysis"]["common"]["npb_cut_on"].as<int>(0);
    const float npb_score_cut  = cfg["analysis"]["common"]["npb_score_cut"].as<float>(0.5f);
    const bool  tight_on       = !bdt_model_name.empty();

    const int vtx_reweight_on = cfg["analysis"]["vertex_reweight_on"].as<int>(1);
    const std::string vtx_reweight_file =
        cfg["analysis"]["vertex_reweight_file"].as<std::string>("results/vertex_reweight.root");

    // vertex reweight histogram
    TH1 *h_vtxrw = nullptr;
    if (vtx_reweight_on)
    {
        TFile *fvtx = TFile::Open(vtx_reweight_file.c_str(), "READ");
        if (fvtx && !fvtx->IsZombie())
        {
            h_vtxrw = dynamic_cast<TH1 *>(fvtx->Get("h_vertexz_ratio_data_over_mccombined"));
            if (h_vtxrw) h_vtxrw->SetDirectory(nullptr);
            fvtx->Close();
        }
        if (!h_vtxrw)
            std::cerr << "[WARNING] vertex reweight histogram not found, running without\n";
    }

    // ── ET windows ────────────────────────────────────────────────────────────
    // Edges define N_win = N_edges-1 windows: [edges[i], edges[i+1])
    const std::vector<float> ET_edges = {10, 15, 20, 25, 30};
    const int N_win = (int)ET_edges.size() - 1;

    // ── find sample params ────────────────────────────────────────────────────
    const auto all_samples = BuildSampleList();
    SampleParams sp;
    bool found = false;
    for (const auto &s : all_samples)
    {
        if (s.filetype == filetype) { sp = s; found = true; break; }
    }
    if (!found)
    {
        std::cerr << "ERROR: unknown filetype '" << filetype << "'\n";
        std::cerr << "Valid: photon5 photon10 photon20 jet5 jet12 jet20 jet30 jet40\n";
        return;
    }

    // ── output ────────────────────────────────────────────────────────────────
    const std::string outfilename = "results/truthjet_ETwindows_" + filetype + ".root";
    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    TH1::SetDefaultSumw2(kTRUE);

    // Per-window histograms
    std::vector<TH1D*> h_maxjet_sig(N_win), h_maxjet_bkg(N_win), h_maxjet_all(N_win);
    std::vector<TH1D*> h_alljet_sig(N_win), h_alljet_bkg(N_win), h_alljet_all(N_win);
    std::vector<TH1D*> h_clusterET_sig(N_win), h_clusterET_bkg(N_win);
    std::vector<TH1D*> h_maxjet_tight(N_win), h_alljet_tight(N_win);
    std::vector<TH1D*> h_maxjet_tightonly(N_win), h_alljet_tightonly(N_win);

    for (int iw = 0; iw < N_win; iw++)
    {
        float lo = ET_edges[iw], hi = ET_edges[iw + 1];
        const char *tag = Form("win%d", iw);

        h_maxjet_sig[iw] = new TH1D(Form("h_maxtruthjet_pt_sig_%s", tag),
            Form("Max truth jet p_{T} (signal) %.0f-%.0f GeV;p_{T}^{jet,max} [GeV];Weighted counts", lo, hi), 100, 0, 100);
        h_maxjet_bkg[iw] = new TH1D(Form("h_maxtruthjet_pt_bkg_%s", tag),
            Form("Max truth jet p_{T} (bkg) %.0f-%.0f GeV;p_{T}^{jet,max} [GeV];Weighted counts",    lo, hi), 100, 0, 100);
        h_maxjet_all[iw] = new TH1D(Form("h_maxtruthjet_pt_all_%s", tag),
            Form("Max truth jet p_{T} (all) %.0f-%.0f GeV;p_{T}^{jet,max} [GeV];Weighted counts",    lo, hi), 100, 0, 100);

        h_alljet_sig[iw] = new TH1D(Form("h_alltruthjet_pt_sig_%s", tag),
            Form("All truth jet p_{T} (signal) %.0f-%.0f GeV;p_{T}^{jet} [GeV];Weighted counts", lo, hi), 100, 0, 100);
        h_alljet_bkg[iw] = new TH1D(Form("h_alltruthjet_pt_bkg_%s", tag),
            Form("All truth jet p_{T} (bkg) %.0f-%.0f GeV;p_{T}^{jet} [GeV];Weighted counts",    lo, hi), 100, 0, 100);
        h_alljet_all[iw] = new TH1D(Form("h_alltruthjet_pt_all_%s", tag),
            Form("All truth jet p_{T} (all) %.0f-%.0f GeV;p_{T}^{jet} [GeV];Weighted counts",    lo, hi), 100, 0, 100);

        h_clusterET_sig[iw] = new TH1D(Form("h_clusterET_sig_%s", tag),
            Form("Cluster E_{T} (signal) %.0f-%.0f GeV;E_{T} [GeV];Weighted counts", lo, hi), 20, lo, hi);
        h_clusterET_bkg[iw] = new TH1D(Form("h_clusterET_bkg_%s", tag),
            Form("Cluster E_{T} (bkg) %.0f-%.0f GeV;E_{T} [GeV];Weighted counts",    lo, hi), 20, lo, hi);

        h_maxjet_tight[iw] = new TH1D(Form("h_maxtruthjet_pt_tight_%s", tag),
            Form("Max truth jet p_{T} (tight+iso) %.0f-%.0f GeV;p_{T}^{jet,max} [GeV];Weighted counts", lo, hi), 100, 0, 100);
        h_alljet_tight[iw] = new TH1D(Form("h_alltruthjet_pt_tight_%s", tag),
            Form("All truth jet p_{T} (tight+iso) %.0f-%.0f GeV;p_{T}^{jet} [GeV];Weighted counts",     lo, hi), 100, 0, 100);

        h_maxjet_tightonly[iw] = new TH1D(Form("h_maxtruthjet_pt_tightonly_%s", tag),
            Form("Max truth jet p_{T} (tight only) %.0f-%.0f GeV;p_{T}^{jet,max} [GeV];Weighted counts", lo, hi), 100, 0, 100);
        h_alljet_tightonly[iw] = new TH1D(Form("h_alltruthjet_pt_tightonly_%s", tag),
            Form("All truth jet p_{T} (tight only) %.0f-%.0f GeV;p_{T}^{jet} [GeV];Weighted counts",     lo, hi), 100, 0, 100);
    }

    // 2D: cluster ET vs max truth jet pT (full ET range, filled for every cluster)
    TH2D *h2 = new TH2D("h2_clusterET_maxjet",
                         "Cluster E_{T} vs max truth jet p_{T};Cluster E_{T} [GeV];Max truth jet p_{T} [GeV]",
                         80, 5, 45, 100, 0, 100);

    std::cout << "\n=== TruthJetInETWindow  (" << N_win << " windows, sample=" << filetype << ") ===\n";

    // ── process single sample ─────────────────────────────────────────────────
    {
        std::string infile = root_dir + sp.filetype + branch_dir;
        std::cout << "  Processing " << sp.filetype << "  (" << infile << ")\n";

        TChain chain(treename.c_str());
        chain.Add(infile.c_str());
        if (chain.GetEntries() == 0)
        {
            std::cerr << "  [ERROR] no entries found in " << infile << "\n";
            return;
        }

        TTreeReader reader(&chain);

        TTreeReaderValue<int>   nparticles (reader, "nparticles");
        TTreeReaderValue<int>   ncluster   (reader, Form("ncluster_%s", clusternodename.c_str()));
        TTreeReaderValue<float> vertexz    (reader, "vertexz");
        TTreeReaderValue<int>   njet_truth (reader, "njet_truth");

        TTreeReaderArray<float> jet_truth_Pt  (reader, "jet_truth_Pt");
        TTreeReaderArray<float> jet_truth_Eta (reader, "jet_truth_Eta");

        TTreeReaderArray<float> particle_Pt          (reader, "particle_Pt");
        TTreeReaderArray<float> particle_Eta         (reader, "particle_Eta");
        TTreeReaderArray<float> particle_Phi         (reader, "particle_Phi");
        TTreeReaderArray<float> particle_truth_iso_03(reader, "particle_truth_iso_03");
        TTreeReaderArray<int>   particle_pid         (reader, "particle_pid");
        TTreeReaderArray<int>   particle_trkid       (reader, "particle_trkid");
        TTreeReaderArray<int>   particle_photonclass  (reader, "particle_photonclass");

        TTreeReaderArray<float> cluster_Et  (reader, Form("cluster_Et_%s",       clusternodename.c_str()));
        TTreeReaderArray<float> cluster_Eta (reader, Form("cluster_Eta_%s",      clusternodename.c_str()));
        TTreeReaderArray<float> cluster_Phi (reader, Form("cluster_Phi_%s",      clusternodename.c_str()));
        TTreeReaderArray<int>   cluster_trkID(reader, Form("cluster_truthtrkID_%s", clusternodename.c_str()));

        // Optional tight-selection branches (only present in bdt_split.root files)
        std::string bdt_branch = "cluster_bdt_" + clusternodename + "_" + bdt_model_name;
        std::string iso_branch = "cluster_iso_03_" + clusternodename;
        std::string npb_branch = "cluster_npb_score_" + clusternodename;
        TTreeReaderArray<float> *cluster_bdt_arr = nullptr;
        TTreeReaderArray<float> *cluster_iso_arr = nullptr;
        TTreeReaderArray<float> *cluster_npb_arr = nullptr;
        if (tight_on)
        {
            cluster_bdt_arr = new TTreeReaderArray<float>(reader, bdt_branch.c_str());
            cluster_iso_arr = new TTreeReaderArray<float>(reader, iso_branch.c_str());
            if (npb_cut_on)
                cluster_npb_arr = new TTreeReaderArray<float>(reader, npb_branch.c_str());
        }

        const long nentries = chain.GetEntries();
        long ientry = 0;

        while (reader.Next())
        {
            if (ientry % 100000 == 0)
                std::cout << "    " << ientry << " / " << nentries << "\r" << std::flush;

            // ── per-event weight ─────────────────────────────────────────────
            float weight = sp.cross_weight;
            if (h_vtxrw)
            {
                int bin = h_vtxrw->FindBin(*vertexz);
                bin = std::max(1, std::min(bin, h_vtxrw->GetNbinsX()));
                float vw = h_vtxrw->GetBinContent(bin);
                if (std::isfinite(vw) && vw > 0.f) weight *= vw;
            }

            // ── max photon pT filter (photon samples) ────────────────────────
            if (!sp.isbackground)
            {
                float maxphotonpT = 0.f;
                for (int ip = 0; ip < *nparticles; ip++)
                    if (particle_pid[ip] == 22 && particle_Pt[ip] > maxphotonpT)
                        maxphotonpT = particle_Pt[ip];
                if (maxphotonpT < sp.max_photon_lower || maxphotonpT >= sp.max_photon_upper)
                { ientry++; continue; }
            }

            // ── max jet pT filter (jet samples) ─────────────────────────────
            float maxjetpT = 0.f;
            for (int ij = 0; ij < *njet_truth; ij++)
                if (jet_truth_Pt[ij] > maxjetpT) maxjetpT = jet_truth_Pt[ij];

            if (sp.isbackground)
            {
                if (maxjetpT == 0.f || maxjetpT < sp.max_jet_lower || maxjetpT >= sp.max_jet_upper)
                { ientry++; continue; }
            }

            // ── vertex cut ───────────────────────────────────────────────────
            if (std::abs(*vertexz) > vertexcut) { ientry++; continue; }

            // ── build truth photon map ───────────────────────────────────────
            std::map<int, int> trkid_to_ipart;  // trkid → particle index
            std::map<int, int> photon_map;       // particle index → class (1=direct,2=frag)
            for (int ip = 0; ip < *nparticles; ip++)
            {
                trkid_to_ipart[particle_trkid[ip]] = ip;
                if (particle_pid[ip] == 22 &&
                    (particle_photonclass[ip] == 1 || particle_photonclass[ip] == 2) &&
                    particle_truth_iso_03[ip] < truthisocut)
                {
                    photon_map[ip] = particle_photonclass[ip];
                }
            }

            // ── collect truth jet info for this event ────────────────────────
            // (used below for every qualifying cluster in the event)
            float max_truthjet_in_eta = 0.f;
            for (int ij = 0; ij < *njet_truth; ij++)
                if (std::abs(jet_truth_Eta[ij]) < eta_max_jet && jet_truth_Pt[ij] > max_truthjet_in_eta)
                    max_truthjet_in_eta = jet_truth_Pt[ij];

            // ── cluster loop ─────────────────────────────────────────────────
            for (int ic = 0; ic < *ncluster; ic++)
            {
                float clET  = cluster_Et[ic];
                float clEta = cluster_Eta[ic];

                // basic fiducial
                if (std::abs(clEta) >= eta_max_cluster) continue;
                if (sp.isbackground && clET > sp.cluster_ET_upper) continue;

                // fill 2D for full ET range (no window requirement)
                h2->Fill(clET, max_truthjet_in_eta, weight);

                // find which window this cluster falls into (-1 if none)
                int iwin = -1;
                for (int iw = 0; iw < N_win; iw++)
                {
                    if (clET >= ET_edges[iw] && clET < ET_edges[iw + 1])
                    { iwin = iw; break; }
                }
                if (iwin == -1) continue;

                // ── truth classification ─────────────────────────────────────
                int cluster_class = 0; // 0=bkg, 1=direct, 2=frag
                auto it = trkid_to_ipart.find(cluster_trkID[ic]);
                if (it != trkid_to_ipart.end())
                {
                    int ip = it->second;
                    float deta = clEta - particle_Eta[ip];
                    float dphi = cluster_Phi[ic] - particle_Phi[ip];
                    while (dphi >  M_PI) dphi -= 2 * M_PI;
                    while (dphi < -M_PI) dphi += 2 * M_PI;
                    float dR = std::sqrt(deta * deta + dphi * dphi);
                    if (dR < eff_dR)
                    {
                        auto jt = photon_map.find(ip);
                        if (jt != photon_map.end()) cluster_class = jt->second;
                    }
                }

                const bool is_sig = (cluster_class == 1 || cluster_class == 2);

                // ── tight + iso selection ────────────────────────────────────
                bool is_tightonly = false;
                bool is_tight     = false;
                if (tight_on && cluster_bdt_arr && cluster_iso_arr)
                {
                    bool bdt_pass = ((*cluster_bdt_arr)[ic] > tight_bdt_min &&
                                     (*cluster_bdt_arr)[ic] < tight_bdt_max) &&
                                    (!npb_cut_on || !cluster_npb_arr ||
                                     (*cluster_npb_arr)[ic] > npb_score_cut);
                    float iso_max = reco_iso_max_b + reco_iso_max_s * clET;
                    bool  iso_pass = ((*cluster_iso_arr)[ic] < iso_max);
                    is_tightonly = bdt_pass;
                    is_tight     = bdt_pass && iso_pass;
                }

                // ── fill 1D histograms for this window ───────────────────────
                h_maxjet_all[iwin]->Fill(max_truthjet_in_eta, weight);
                if (is_sig)   h_maxjet_sig[iwin]->Fill(max_truthjet_in_eta, weight);
                else          h_maxjet_bkg[iwin]->Fill(max_truthjet_in_eta, weight);
                if (is_tightonly) h_maxjet_tightonly[iwin]->Fill(max_truthjet_in_eta, weight);
                if (is_tight)     h_maxjet_tight[iwin]->Fill(max_truthjet_in_eta, weight);

                for (int ij = 0; ij < *njet_truth; ij++)
                {
                    if (std::abs(jet_truth_Eta[ij]) >= eta_max_jet) continue;
                    h_alljet_all[iwin]->Fill(jet_truth_Pt[ij], weight);
                    if (is_sig)       h_alljet_sig[iwin]->Fill(jet_truth_Pt[ij], weight);
                    else              h_alljet_bkg[iwin]->Fill(jet_truth_Pt[ij], weight);
                    if (is_tightonly) h_alljet_tightonly[iwin]->Fill(jet_truth_Pt[ij], weight);
                    if (is_tight)     h_alljet_tight[iwin]->Fill(jet_truth_Pt[ij], weight);
                }

                if (is_sig) h_clusterET_sig[iwin]->Fill(clET, weight);
                else        h_clusterET_bkg[iwin]->Fill(clET, weight);

            } // end cluster loop
            ientry++;
        } // end event loop
        std::cout << "    done (" << nentries << " entries)\n";

        delete cluster_bdt_arr;
        delete cluster_iso_arr;
        delete cluster_npb_arr;
    } // end sample block

    // ── write histograms ──────────────────────────────────────────────────────
    fout->cd();

    // Store window edges as a histogram so the plotting macro can read them back
    TH1D *h_ET_edges_hist = new TH1D("h_ET_edges", "ET window edges", N_win,
                                      ET_edges.data());
    h_ET_edges_hist->Write();
    h2->Write();

    for (int iw = 0; iw < N_win; iw++)
    {
        h_maxjet_sig[iw]->Write();
        h_maxjet_bkg[iw]->Write();
        h_maxjet_all[iw]->Write();
        h_alljet_sig[iw]->Write();
        h_alljet_bkg[iw]->Write();
        h_alljet_all[iw]->Write();
        h_clusterET_sig[iw]->Write();
        h_clusterET_bkg[iw]->Write();
        h_maxjet_tightonly[iw]->Write();
        h_alljet_tightonly[iw]->Write();
        h_maxjet_tight[iw]->Write();
        h_alljet_tight[iw]->Write();
    }

    fout->Close();
    std::cout << "\nOutput written to: " << outfilename << "\n";
    std::cout << "hadd all per-sample files, then run PlotTruthJetWindows.C.\n";
}
