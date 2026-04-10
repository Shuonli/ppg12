// study_vertex_double_interaction.C
// Compare vertex distributions and deta-dz correlations between
// single-interaction and full GEANT double-interaction MC (photon10).
//
// Produces:
// 1. Event-level vertex distributions (reco, truth, dz) for both samples
// 2. Per-cluster deta vs dz correlation for truth-matched clusters
// 3. Supporting diagnostics (deltaR, ET ratio, deta vs eta)
//
// Usage: root -l -b -q 'study_vertex_double_interaction.C(-1)'
//        root -l -b -q 'study_vertex_double_interaction.C(2000000)'

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>

#include "CrossSectionWeights.h"

void study_vertex_double_interaction(int maxEvents = -1)
{
    gStyle->SetOptStat(0);

    // =====================================================================
    // Constants
    // =====================================================================
    const TString treeName = "slimtree";
    const TString clusterNode = "CLUSTERINFO_CEMC";
    const float vertex_cut = 60.0;   // |vtx_reco| < 60 cm
    const float eta_max    = 0.7;
    const float pT_min     = 8.0;
    const float pT_max     = 36.0;
    const float truth_iso_max = 4.0;
    const float R_CEMC     = 93.5;   // EMCal barrel radius [cm]

    // =====================================================================
    // Samples
    // =====================================================================
    struct SampleInfo { TString name; TString path; };
    std::vector<SampleInfo> samples = {
        {"single", "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root"},
        {"double", "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_double/bdt_split.root"},
    };

    const TString outFile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/vertex_double_interaction_study.root";
    TFile *fout = new TFile(outFile, "RECREATE");

    const int nSamples = (int)samples.size();

    // =====================================================================
    // Histograms — event-level (filled BEFORE vertex cut)
    // =====================================================================
    std::vector<TH1F*> h_vtxz_reco(nSamples);
    std::vector<TH1F*> h_vtxz_truth(nSamples);
    std::vector<TH1F*> h_dz_event(nSamples);
    std::vector<TH1F*> h_dz_abs_event(nSamples);
    std::vector<TH2F*> h_vtx_reco_vs_truth(nSamples);

    // =====================================================================
    // Histograms — per truth-matched cluster (after vertex cut + selection)
    // =====================================================================
    std::vector<TH1F*> h_deta(nSamples);
    std::vector<TH1F*> h_deltaR(nSamples);
    std::vector<TH1F*> h_dz_matched(nSamples);
    std::vector<TH2F*> h_deta_vs_dz(nSamples);
    std::vector<TH2F*> h_deltaR_vs_dz(nSamples);
    std::vector<TH1F*> h_ET_ratio(nSamples);
    std::vector<TH2F*> h_deta_vs_eta(nSamples);
    std::vector<TH2F*> h_dET_vs_dz(nSamples);

    for (int is = 0; is < nSamples; is++)
    {
        TString s = samples[is].name;

        // -- Event-level --
        h_vtxz_reco[is] = new TH1F(Form("h_vtxz_reco_%s", s.Data()),
            Form("Reco vertex [%s];z_{vtx}^{reco} [cm];Events", s.Data()),
            240, -120, 120);
        h_vtxz_reco[is]->Sumw2();

        h_vtxz_truth[is] = new TH1F(Form("h_vtxz_truth_%s", s.Data()),
            Form("Truth vertex [%s];z_{vtx}^{truth} [cm];Events", s.Data()),
            240, -120, 120);
        h_vtxz_truth[is]->Sumw2();

        h_dz_event[is] = new TH1F(Form("h_dz_event_%s", s.Data()),
            Form("dz = vtx_{reco} - vtx_{truth} [%s];dz [cm];Events", s.Data()),
            400, -200, 200);
        h_dz_event[is]->Sumw2();

        h_dz_abs_event[is] = new TH1F(Form("h_dz_abs_event_%s", s.Data()),
            Form("|dz| [%s];|dz| [cm];Events", s.Data()),
            200, 0, 200);
        h_dz_abs_event[is]->Sumw2();

        h_vtx_reco_vs_truth[is] = new TH2F(Form("h_vtx_reco_vs_truth_%s", s.Data()),
            Form("Reco vs truth vertex [%s];z_{vtx}^{truth} [cm];z_{vtx}^{reco} [cm]", s.Data()),
            120, -120, 120, 120, -120, 120);
        h_vtx_reco_vs_truth[is]->Sumw2();

        // -- Per-cluster (wide ranges to accommodate double interaction) --
        h_deta[is] = new TH1F(Form("h_deta_%s", s.Data()),
            Form("#Delta#eta [%s];#Delta#eta (cluster - truth);Counts", s.Data()),
            400, -1.5, 1.5);
        h_deta[is]->Sumw2();

        h_deltaR[is] = new TH1F(Form("h_deltaR_%s", s.Data()),
            Form("#DeltaR [%s];#DeltaR;Counts", s.Data()),
            500, 0, 2.0);
        h_deltaR[is]->Sumw2();

        h_dz_matched[is] = new TH1F(Form("h_dz_matched_%s", s.Data()),
            Form("dz (matched clusters) [%s];dz [cm];Counts", s.Data()),
            400, -200, 200);
        h_dz_matched[is]->Sumw2();

        h_deta_vs_dz[is] = new TH2F(Form("h_deta_vs_dz_%s", s.Data()),
            Form("#Delta#eta vs dz [%s];dz [cm];#Delta#eta", s.Data()),
            200, -150, 150, 200, -1.5, 1.5);
        h_deta_vs_dz[is]->Sumw2();

        h_deltaR_vs_dz[is] = new TH2F(Form("h_deltaR_vs_dz_%s", s.Data()),
            Form("#DeltaR vs dz [%s];dz [cm];#DeltaR", s.Data()),
            200, -150, 150, 200, 0, 2.0);
        h_deltaR_vs_dz[is]->Sumw2();

        h_ET_ratio[is] = new TH1F(Form("h_ET_ratio_%s", s.Data()),
            Form("E_{T}^{cluster}/p_{T}^{truth} [%s];E_{T}/p_{T}^{truth};Counts", s.Data()),
            200, 0, 2);
        h_ET_ratio[is]->Sumw2();

        h_deta_vs_eta[is] = new TH2F(Form("h_deta_vs_eta_%s", s.Data()),
            Form("#Delta#eta vs #eta_{truth} [%s];#eta_{truth};#Delta#eta", s.Data()),
            28, -0.7, 0.7, 200, -1.5, 1.5);
        h_deta_vs_eta[is]->Sumw2();

        h_dET_vs_dz[is] = new TH2F(Form("h_dET_vs_dz_%s", s.Data()),
            Form("#DeltaE_{T}/E_{T} vs dz [%s];dz [cm];(E_{T}^{cluster} - p_{T}^{truth})/p_{T}^{truth}", s.Data()),
            200, -150, 150, 200, -1.0, 1.0);
        h_dET_vs_dz[is]->Sumw2();
    }

    // =====================================================================
    // Process samples
    // =====================================================================
    for (int is = 0; is < nSamples; is++)
    {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Processing " << samples[is].name << ": " << samples[is].path << std::endl;
        std::cout << "========================================" << std::endl;

        TFile *fin = TFile::Open(samples[is].path);
        if (!fin || fin->IsZombie())
        {
            std::cerr << "ERROR: Cannot open " << samples[is].path << std::endl;
            continue;
        }

        TTreeReader reader(treeName, fin);
        TTreeReaderValue<int>   nparticles(reader, "nparticles");
        TTreeReaderValue<int>   nclusters(reader, Form("ncluster_%s", clusterNode.Data()));
        TTreeReaderValue<float> vertexz(reader, "vertexz");
        TTreeReaderValue<float> vertexz_truth(reader, "vertexz_truth");

        TTreeReaderArray<float> particle_Pt(reader, "particle_Pt");
        TTreeReaderArray<float> particle_Eta(reader, "particle_Eta");
        TTreeReaderArray<float> particle_Phi(reader, "particle_Phi");
        TTreeReaderArray<int>   particle_pid(reader, "particle_pid");
        TTreeReaderArray<int>   particle_trkid(reader, "particle_trkid");
        TTreeReaderArray<int>   particle_photonclass(reader, "particle_photonclass");
        TTreeReaderArray<float> particle_truth_iso_03(reader, "particle_truth_iso_03");

        TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusterNode.Data()));
        TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusterNode.Data()));
        TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusterNode.Data()));
        TTreeReaderArray<int>   cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusterNode.Data()));

        long long nEvents  = 0;
        long long nTruth   = 0;
        long long nMatched = 0;
        long long nNoTrkID = 0;
        long long nFailDR  = 0;

        while (reader.Next())
        {
            if (maxEvents > 0 && nEvents >= maxEvents) break;
            nEvents++;

            if (nEvents % 500000 == 0)
                std::cout << "  [" << samples[is].name << "] Event " << nEvents << std::endl;

            float vtx_reco  = *vertexz;
            float vtx_truth = *vertexz_truth;
            float dz = vtx_reco - vtx_truth;

            // Event-level histograms (before vertex cut)
            h_vtxz_reco[is]->Fill(vtx_reco);
            h_vtxz_truth[is]->Fill(vtx_truth);
            h_dz_event[is]->Fill(dz);
            h_dz_abs_event[is]->Fill(fabs(dz));
            h_vtx_reco_vs_truth[is]->Fill(vtx_truth, vtx_reco);

            // Vertex cut
            if (fabs(vtx_reco) > vertex_cut) continue;

            int ncl   = *nclusters;
            int npart = *nparticles;

            // Build trkID -> cluster index map
            std::map<int, std::vector<int>> trkid_to_clusters;
            for (int icl = 0; icl < ncl; icl++)
                trkid_to_clusters[cluster_truthtrkID[icl]].push_back(icl);

            // Loop over truth particles
            for (int ip = 0; ip < npart; ip++)
            {
                if (particle_pid[ip] != 22) continue;
                if (particle_photonclass[ip] >= 3) continue;
                if (particle_truth_iso_03[ip] >= truth_iso_max) continue;
                if (fabs(particle_Eta[ip]) >= eta_max) continue;
                if (particle_Pt[ip] < pT_min || particle_Pt[ip] >= pT_max) continue;

                float truth_pT  = particle_Pt[ip];
                float truth_eta = particle_Eta[ip];
                float truth_phi = particle_Phi[ip];
                int   truth_trkid = particle_trkid[ip];

                nTruth++;

                auto itMap = trkid_to_clusters.find(truth_trkid);
                if (itMap == trkid_to_clusters.end())
                {
                    nNoTrkID++;
                    continue;
                }

                // Find closest trkID-matched cluster by deltaR
                float best_dR   = 999;
                int   best_icl  = -1;
                float best_deta = 999;

                for (int icl : itMap->second)
                {
                    float deta_cl = cluster_Eta[icl] - truth_eta;
                    float dphi_cl = cluster_Phi[icl] - truth_phi;
                    if (dphi_cl >  TMath::Pi()) dphi_cl -= 2 * TMath::Pi();
                    if (dphi_cl < -TMath::Pi()) dphi_cl += 2 * TMath::Pi();
                    float dR = sqrt(deta_cl * deta_cl + dphi_cl * dphi_cl);
                    if (dR < best_dR)
                    {
                        best_dR   = dR;
                        best_icl  = icl;
                        best_deta = deta_cl;
                    }
                }

                nMatched++;
                if (best_dR >= 0.1) nFailDR++;

                float ET_ratio = cluster_Et[best_icl] / truth_pT;

                // Fill per-cluster histograms
                h_deta[is]->Fill(best_deta);
                h_deltaR[is]->Fill(best_dR);
                h_dz_matched[is]->Fill(dz);
                h_deta_vs_dz[is]->Fill(dz, best_deta);
                h_deltaR_vs_dz[is]->Fill(dz, best_dR);
                h_ET_ratio[is]->Fill(ET_ratio);
                h_deta_vs_eta[is]->Fill(truth_eta, best_deta);
                h_dET_vs_dz[is]->Fill(dz, ET_ratio - 1.0);
            }
        }

        fin->Close();

        // Summary
        std::cout << "\n--- Summary for " << samples[is].name << " ---" << std::endl;
        std::cout << "Events processed: " << nEvents << std::endl;
        std::cout << "Truth photons selected: " << nTruth << std::endl;
        std::cout << "  No trkID match:    " << nNoTrkID
                  << " (" << (nTruth > 0 ? 100.0 * nNoTrkID / nTruth : 0) << "%)" << std::endl;
        std::cout << "  trkID matched:     " << nMatched << std::endl;
        std::cout << "  Fail dR < 0.1:     " << nFailDR
                  << " (" << (nMatched > 0 ? 100.0 * nFailDR / nMatched : 0) << "%)" << std::endl;
        std::cout << "  Vertex dz RMS:     " << h_dz_event[is]->GetRMS() << " cm" << std::endl;
        std::cout << "  Vertex dz mean:    " << h_dz_event[is]->GetMean() << " cm" << std::endl;
        std::cout << "  Frac |dz| > 9.35:  "
                  << (h_dz_abs_event[is]->Integral() > 0
                      ? 100.0 * h_dz_abs_event[is]->Integral(
                            h_dz_abs_event[is]->FindBin(9.35),
                            h_dz_abs_event[is]->GetNbinsX() + 1)
                        / h_dz_abs_event[is]->Integral()
                      : 0)
                  << "%" << std::endl;
        std::cout << "  Cluster deta RMS:  " << h_deta[is]->GetRMS() << std::endl;
    }

    // =====================================================================
    // Write output
    // =====================================================================
    fout->cd();
    for (int is = 0; is < nSamples; is++)
    {
        h_vtxz_reco[is]->Write();
        h_vtxz_truth[is]->Write();
        h_dz_event[is]->Write();
        h_dz_abs_event[is]->Write();
        h_vtx_reco_vs_truth[is]->Write();
        h_deta[is]->Write();
        h_deltaR[is]->Write();
        h_dz_matched[is]->Write();
        h_deta_vs_dz[is]->Write();
        h_deltaR_vs_dz[is]->Write();
        h_ET_ratio[is]->Write();
        h_deta_vs_eta[is]->Write();
        h_dET_vs_dz[is]->Write();
    }
    fout->Close();
    std::cout << "\nOutput saved to: " << outFile << std::endl;
}
