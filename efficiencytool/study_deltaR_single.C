// study_deltaR_single.C
// Study the deltaR truth-matching requirement in single-interaction MC.
// Characterize the deltaR distribution, its dependence on truth-reco vertex
// mismatch, and the efficiency impact of varying the deltaR threshold.
//
// Usage: root -l -b -q 'study_deltaR_single.C(-1)'
//        root -l -b -q 'study_deltaR_single.C(500000)'

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TEfficiency.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>

#include "CrossSectionWeights.h"

void study_deltaR_single(int maxEvents = -1)
{
    gStyle->SetOptStat(0);

    // =====================================================================
    // Constants
    // =====================================================================
    const TString treeName = "slimtree";
    const TString clusterNode = "CLUSTERINFO_CEMC";
    const float vertex_cut = 60.0;   // |vtx_reco| < 60 cm
    const float eta_max = 0.7;
    const float pT_min = 8.0;
    const float pT_max = 36.0;
    const float truth_iso_max = 4.0; // truth isolation cone R=0.3
    const float R_CEMC = 93.5;       // EMCal barrel radius in cm

    // Samples
    struct SampleInfo { TString name; TString path; float xsec; };
    std::vector<SampleInfo> samples = {
        {"photon5",  "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon5/bdt_split.root",  PPG12::photon5cross},
        {"photon10", "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root", PPG12::photon10cross},
        {"photon20", "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon20/bdt_split.root", PPG12::photon20cross},
    };

    // deltaR thresholds to scan
    const std::vector<float> dR_thresholds = {0.05f, 0.08f, 0.1f, 0.15f, 0.2f, 0.3f};

    // pT bin edges (from plotcommon.h)
    const int NptBins = 12;
    const double ptEdges[13] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36};

    // Output
    const TString outFile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/deltaR_single_study.root";

    const int nSamples = (int)samples.size();

    // =====================================================================
    // Helper: make threshold label string (0.1 -> "0p1")
    // =====================================================================
    auto makeLabel = [](float t) -> TString {
        TString s = Form("%.2f", t);
        s.ReplaceAll(".", "p");
        // strip trailing zeros after 'p'
        while (s.EndsWith("0") && !s.EndsWith("p0")) s.Remove(s.Length() - 1);
        return s;
    };

    // =====================================================================
    // Open output file
    // =====================================================================
    TFile *fout = new TFile(outFile, "RECREATE");

    // =====================================================================
    // Histogram bookkeeping
    // We create per-sample and "combined" histograms.
    // Index: 0..nSamples-1 = per-sample, nSamples = combined
    // =====================================================================
    const int nSets = nSamples + 1; // 3 samples + combined
    std::vector<TString> setLabels;
    for (auto &s : samples) setLabels.push_back(s.name);
    setLabels.push_back("combined");

    // --- 1D distributions ---
    std::vector<TH1F*> h_deltaR(nSets);
    std::vector<TH1F*> h_deltaR_fine(nSets);
    std::vector<TH1F*> h_dz(nSets);
    std::vector<TH1F*> h_deta(nSets);
    std::vector<TH1F*> h_dphi(nSets);
    std::vector<TH1F*> h_ET_ratio_pass(nSets);
    std::vector<TH1F*> h_ET_ratio_fail(nSets);
    std::vector<TH1F*> h_pT_truth_all(nSets);
    std::vector<TH1F*> h_pT_truth_notrkid(nSets);

    // --- 2D distributions ---
    std::vector<TH2F*> h_deltaR_vs_pT(nSets);
    std::vector<TH2F*> h_deltaR_vs_dz(nSets);
    std::vector<TH2F*> h_deltaR_vs_eta(nSets);
    std::vector<TH2F*> h_deta_vs_dz(nSets);

    for (int iset = 0; iset < nSets; iset++)
    {
        TString s = setLabels[iset];

        h_deltaR[iset] = new TH1F(Form("h_deltaR_%s", s.Data()),
            Form("Min #DeltaR (trkID matched) [%s];min #DeltaR;Weighted counts", s.Data()),
            300, 0, 0.3);
        h_deltaR[iset]->Sumw2();

        h_deltaR_fine[iset] = new TH1F(Form("h_deltaR_fine_%s", s.Data()),
            Form("Min #DeltaR fine [%s];min #DeltaR;Weighted counts", s.Data()),
            1000, 0, 0.5);
        h_deltaR_fine[iset]->Sumw2();

        h_dz[iset] = new TH1F(Form("h_dz_%s", s.Data()),
            Form("dz = vtx_{reco} - vtx_{truth} [%s];dz [cm];Weighted counts", s.Data()),
            200, -30, 30);
        h_dz[iset]->Sumw2();

        h_deta[iset] = new TH1F(Form("h_deta_%s", s.Data()),
            Form("#Delta#eta (cluster - truth) [%s];#Delta#eta;Weighted counts", s.Data()),
            200, -0.15, 0.15);
        h_deta[iset]->Sumw2();

        h_dphi[iset] = new TH1F(Form("h_dphi_%s", s.Data()),
            Form("#Delta#phi (cluster - truth) [%s];#Delta#phi;Weighted counts", s.Data()),
            200, -0.15, 0.15);
        h_dphi[iset]->Sumw2();

        h_ET_ratio_pass[iset] = new TH1F(Form("h_ET_ratio_pass_%s", s.Data()),
            Form("E_{T}^{cluster}/p_{T}^{truth} (dR < 0.1) [%s];E_{T}/p_{T};Weighted counts", s.Data()),
            200, 0, 2);
        h_ET_ratio_pass[iset]->Sumw2();

        h_ET_ratio_fail[iset] = new TH1F(Form("h_ET_ratio_fail_%s", s.Data()),
            Form("E_{T}^{cluster}/p_{T}^{truth} (dR #geq 0.1) [%s];E_{T}/p_{T};Weighted counts", s.Data()),
            200, 0, 2);
        h_ET_ratio_fail[iset]->Sumw2();

        h_pT_truth_all[iset] = new TH1F(Form("h_pT_truth_all_%s", s.Data()),
            Form("Truth p_{T}, all selected [%s];p_{T} [GeV];Weighted counts", s.Data()),
            NptBins, ptEdges);
        h_pT_truth_all[iset]->Sumw2();

        h_pT_truth_notrkid[iset] = new TH1F(Form("h_pT_truth_notrkid_%s", s.Data()),
            Form("Truth p_{T}, no trkID match [%s];p_{T} [GeV];Weighted counts", s.Data()),
            NptBins, ptEdges);
        h_pT_truth_notrkid[iset]->Sumw2();

        // 2D
        h_deltaR_vs_pT[iset] = new TH2F(Form("h_deltaR_vs_pT_%s", s.Data()),
            Form("#DeltaR vs truth p_{T} [%s];truth p_{T} [GeV];#DeltaR", s.Data()),
            28, 8, 36, 300, 0, 0.3);
        h_deltaR_vs_pT[iset]->Sumw2();

        h_deltaR_vs_dz[iset] = new TH2F(Form("h_deltaR_vs_dz_%s", s.Data()),
            Form("#DeltaR vs dz [%s];dz [cm];#DeltaR", s.Data()),
            60, -30, 30, 300, 0, 0.3);
        h_deltaR_vs_dz[iset]->Sumw2();

        h_deltaR_vs_eta[iset] = new TH2F(Form("h_deltaR_vs_eta_%s", s.Data()),
            Form("#DeltaR vs truth #eta [%s];truth #eta;#DeltaR", s.Data()),
            28, -0.7, 0.7, 300, 0, 0.3);
        h_deltaR_vs_eta[iset]->Sumw2();

        h_deta_vs_dz[iset] = new TH2F(Form("h_deta_vs_dz_%s", s.Data()),
            Form("#Delta#eta vs dz [%s];dz [cm];#Delta#eta", s.Data()),
            60, -30, 30, 200, -0.15, 0.15);
        h_deta_vs_dz[iset]->Sumw2();
    }

    // --- Threshold scan TEfficiency (combined only) ---
    const int nThresh = (int)dR_thresholds.size();
    std::vector<TEfficiency*> eff_reco_dR(nThresh);
    for (int it = 0; it < nThresh; it++)
    {
        TString tLabel = makeLabel(dR_thresholds[it]);
        TString ename = Form("eff_reco_dR%s_combined", tLabel.Data());
        eff_reco_dR[it] = new TEfficiency(ename.Data(),
            Form("Reco eff. (dR < %s) [combined];truth p_{T} [GeV];Efficiency", tLabel.Data()),
            NptBins, ptEdges);
        eff_reco_dR[it]->SetUseWeightedEvents();
    }

    TEfficiency *eff_reco_noDR = new TEfficiency("eff_reco_noDR_combined",
        "Reco eff. (trkID match, no dR cut) [combined];truth p_{T} [GeV];Efficiency",
        NptBins, ptEdges);
    eff_reco_noDR->SetUseWeightedEvents();

    TEfficiency *eff_reco_trkid = new TEfficiency("eff_reco_trkid_combined",
        "Reco eff. (any trkID-matched cluster) [combined];truth p_{T} [GeV];Efficiency",
        NptBins, ptEdges);
    eff_reco_trkid->SetUseWeightedEvents();

    // --- Fail-fraction vs |dz| (combined only) ---
    TH1F *h_dz_abs_all_combined = new TH1F("h_dz_abs_all_combined",
        "|dz|, all trkID-matched [combined];|dz| [cm];Weighted counts",
        30, 0, 30);
    h_dz_abs_all_combined->Sumw2();

    TH1F *h_dz_abs_fail_combined = new TH1F("h_dz_abs_fail_combined",
        "|dz|, dR #geq 0.1 [combined];|dz| [cm];Weighted counts",
        30, 0, 30);
    h_dz_abs_fail_combined->Sumw2();

    // =====================================================================
    // Global counters (unweighted, for summary printout)
    // =====================================================================
    long long total_truth_all = 0;
    long long total_no_trkid = 0;
    long long total_fail_dR01 = 0;
    long long total_fail_dR02 = 0;

    // =====================================================================
    // Process each sample
    // =====================================================================
    for (int isample = 0; isample < nSamples; isample++)
    {
        const SampleInfo &si = samples[isample];
        float weight = si.xsec;
        int iCombined = nSamples; // index for combined set

        std::cout << "\n========================================" << std::endl;
        std::cout << "Processing " << si.name << ": " << si.path << std::endl;
        std::cout << "Cross-section weight: " << weight << std::endl;
        std::cout << "========================================" << std::endl;

        TFile *fin = TFile::Open(si.path);
        if (!fin || fin->IsZombie())
        {
            std::cerr << "ERROR: Cannot open " << si.path << std::endl;
            continue;
        }

        TTreeReader reader(treeName, fin);
        TTreeReaderValue<int> nparticles(reader, "nparticles");
        TTreeReaderValue<int> nclusters(reader, Form("ncluster_%s", clusterNode.Data()));
        TTreeReaderValue<float> vertexz(reader, "vertexz");
        TTreeReaderValue<float> vertexz_truth(reader, "vertexz_truth");

        TTreeReaderArray<float> particle_Pt(reader, "particle_Pt");
        TTreeReaderArray<float> particle_Eta(reader, "particle_Eta");
        TTreeReaderArray<float> particle_Phi(reader, "particle_Phi");
        TTreeReaderArray<int> particle_pid(reader, "particle_pid");
        TTreeReaderArray<int> particle_trkid(reader, "particle_trkid");
        TTreeReaderArray<int> particle_photonclass(reader, "particle_photonclass");
        TTreeReaderArray<float> particle_truth_iso_03(reader, "particle_truth_iso_03");

        TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusterNode.Data()));
        TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusterNode.Data()));
        TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusterNode.Data()));
        TTreeReaderArray<int> cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusterNode.Data()));

        // Per-sample counters (unweighted)
        long long nEvents = 0;
        long long nTruth = 0;
        long long nNoTrkID = 0;
        long long nFail01 = 0;
        long long nFail02 = 0;

        while (reader.Next())
        {
            if (maxEvents > 0 && nEvents >= maxEvents) break;
            nEvents++;

            if (nEvents % 100000 == 0)
            {
                std::cout << "  [" << si.name << "] Event " << nEvents;
                if (maxEvents > 0) std::cout << " / " << maxEvents;
                std::cout << std::endl;
            }

            float vtx_reco = *vertexz;
            float vtx_truth = *vertexz_truth;

            if (fabs(vtx_reco) > vertex_cut) continue;

            float dz = vtx_reco - vtx_truth;
            int ncl = *nclusters;
            int npart = *nparticles;

            // Build trkID -> cluster index map
            std::map<int, std::vector<int>> trkid_to_clusters;
            for (int icl = 0; icl < ncl; icl++)
            {
                trkid_to_clusters[cluster_truthtrkID[icl]].push_back(icl);
            }

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
                total_truth_all++;

                // Fill denominator histograms (per-sample + combined)
                h_pT_truth_all[isample]->Fill(truth_pT, weight);
                h_pT_truth_all[iCombined]->Fill(truth_pT, weight);

                // Find best trkID-matched cluster
                auto itMap = trkid_to_clusters.find(truth_trkid);
                if (itMap == trkid_to_clusters.end())
                {
                    // No trkID-matched cluster
                    nNoTrkID++;
                    total_no_trkid++;

                    h_pT_truth_notrkid[isample]->Fill(truth_pT, weight);
                    h_pT_truth_notrkid[iCombined]->Fill(truth_pT, weight);

                    // TEfficiency: all fail
                    for (int it = 0; it < nThresh; it++)
                        eff_reco_dR[it]->FillWeighted(false, weight, truth_pT);
                    eff_reco_noDR->FillWeighted(false, weight, truth_pT);
                    eff_reco_trkid->FillWeighted(false, weight, truth_pT);
                    continue;
                }

                // trkID match exists
                eff_reco_trkid->FillWeighted(true, weight, truth_pT);

                // Find closest cluster by deltaR among trkID matches
                float best_dR = 999;
                int best_icl = -1;
                float best_deta = 999;
                float best_dphi = 999;

                for (int icl : itMap->second)
                {
                    float deta_cl = cluster_Eta[icl] - truth_eta;
                    float dphi_cl = cluster_Phi[icl] - truth_phi;
                    if (dphi_cl > TMath::Pi()) dphi_cl -= 2 * TMath::Pi();
                    else if (dphi_cl < -TMath::Pi()) dphi_cl += 2 * TMath::Pi();
                    float dR = sqrt(deta_cl * deta_cl + dphi_cl * dphi_cl);
                    if (dR < best_dR)
                    {
                        best_dR = dR;
                        best_icl = icl;
                        best_deta = deta_cl;
                        best_dphi = dphi_cl;
                    }
                }

                // --- Fill 1D histograms (per-sample + combined) ---
                h_deltaR[isample]->Fill(best_dR, weight);
                h_deltaR[iCombined]->Fill(best_dR, weight);

                h_deltaR_fine[isample]->Fill(best_dR, weight);
                h_deltaR_fine[iCombined]->Fill(best_dR, weight);

                h_dz[isample]->Fill(dz, weight);
                h_dz[iCombined]->Fill(dz, weight);

                h_deta[isample]->Fill(best_deta, weight);
                h_deta[iCombined]->Fill(best_deta, weight);

                h_dphi[isample]->Fill(best_dphi, weight);
                h_dphi[iCombined]->Fill(best_dphi, weight);

                // ET ratio
                float ET_ratio = cluster_Et[best_icl] / truth_pT;
                if (best_dR < 0.1)
                {
                    h_ET_ratio_pass[isample]->Fill(ET_ratio, weight);
                    h_ET_ratio_pass[iCombined]->Fill(ET_ratio, weight);
                }
                else
                {
                    h_ET_ratio_fail[isample]->Fill(ET_ratio, weight);
                    h_ET_ratio_fail[iCombined]->Fill(ET_ratio, weight);
                }

                // --- Fill 2D histograms (per-sample + combined) ---
                h_deltaR_vs_pT[isample]->Fill(truth_pT, best_dR, weight);
                h_deltaR_vs_pT[iCombined]->Fill(truth_pT, best_dR, weight);

                h_deltaR_vs_dz[isample]->Fill(dz, best_dR, weight);
                h_deltaR_vs_dz[iCombined]->Fill(dz, best_dR, weight);

                h_deltaR_vs_eta[isample]->Fill(truth_eta, best_dR, weight);
                h_deltaR_vs_eta[iCombined]->Fill(truth_eta, best_dR, weight);

                h_deta_vs_dz[isample]->Fill(dz, best_deta, weight);
                h_deta_vs_dz[iCombined]->Fill(dz, best_deta, weight);

                // --- |dz| fail-fraction (combined only) ---
                h_dz_abs_all_combined->Fill(fabs(dz), weight);
                if (best_dR >= 0.1)
                    h_dz_abs_fail_combined->Fill(fabs(dz), weight);

                // --- Threshold scan TEfficiency (combined only) ---
                for (int it = 0; it < nThresh; it++)
                    eff_reco_dR[it]->FillWeighted(best_dR < dR_thresholds[it], weight, truth_pT);
                eff_reco_noDR->FillWeighted(true, weight, truth_pT);

                // Unweighted counters for summary
                if (best_dR >= 0.1) { nFail01++; total_fail_dR01++; }
                if (best_dR >= 0.2) { nFail02++; total_fail_dR02++; }

            } // end truth particle loop

        } // end event loop

        fin->Close();

        // Per-sample summary
        std::cout << "\n--- Summary for " << si.name << " ---" << std::endl;
        std::cout << "Events processed: " << nEvents << std::endl;
        std::cout << "Truth photons selected: " << nTruth << std::endl;
        if (nTruth > 0)
        {
            std::cout << "  No trkID match:  " << nNoTrkID
                      << " (" << 100.0 * nNoTrkID / nTruth << "%)" << std::endl;
            std::cout << "  Fail dR < 0.1:   " << nFail01
                      << " (" << 100.0 * nFail01 / nTruth << "%)" << std::endl;
            std::cout << "  Fail dR < 0.2:   " << nFail02
                      << " (" << 100.0 * nFail02 / nTruth << "%)" << std::endl;
        }

    } // end sample loop

    // =====================================================================
    // Global summary
    // =====================================================================
    std::cout << "\n========================================" << std::endl;
    std::cout << "COMBINED SUMMARY (unweighted counts)" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Total truth photons: " << total_truth_all << std::endl;
    if (total_truth_all > 0)
    {
        std::cout << "  No trkID match:  " << total_no_trkid
                  << " (" << 100.0 * total_no_trkid / total_truth_all << "%)" << std::endl;
        std::cout << "  Fail dR < 0.1:   " << total_fail_dR01
                  << " (" << 100.0 * total_fail_dR01 / total_truth_all << "%)" << std::endl;
        std::cout << "  Fail dR < 0.2:   " << total_fail_dR02
                  << " (" << 100.0 * total_fail_dR02 / total_truth_all << "%)" << std::endl;
    }

    // =====================================================================
    // Write everything to output
    // =====================================================================
    fout->cd();

    for (int iset = 0; iset < nSets; iset++)
    {
        h_deltaR[iset]->Write();
        h_deltaR_fine[iset]->Write();
        h_dz[iset]->Write();
        h_deta[iset]->Write();
        h_dphi[iset]->Write();
        h_ET_ratio_pass[iset]->Write();
        h_ET_ratio_fail[iset]->Write();
        h_pT_truth_all[iset]->Write();
        h_pT_truth_notrkid[iset]->Write();
        h_deltaR_vs_pT[iset]->Write();
        h_deltaR_vs_dz[iset]->Write();
        h_deltaR_vs_eta[iset]->Write();
        h_deta_vs_dz[iset]->Write();
    }

    for (int it = 0; it < nThresh; it++)
        eff_reco_dR[it]->Write();
    eff_reco_noDR->Write();
    eff_reco_trkid->Write();

    h_dz_abs_all_combined->Write();
    h_dz_abs_fail_combined->Write();

    fout->Close();
    std::cout << "\nOutput saved to: " << outFile << std::endl;
}
