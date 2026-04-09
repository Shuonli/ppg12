// investigate_deltaR.C
// Investigate truth-reco cluster matching in single vs double-interaction MC
// Goal: understand whether efficiency loss in double MC is from vertex-shift dR failure
// or from extra truth photons (from second collision) that have no corresponding cluster.
//
// Usage: root -l -b -q 'investigate_deltaR.C(100000)'

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TStyle.h>
#include <TROOT.h>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>

void investigate_deltaR(int maxEvents = 100000)
{
    gStyle->SetOptStat(0);

    const TString singleFile = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root";
    const TString doubleFile = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_double/bdt_split.root";
    const TString outFile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/deltaR_study.root";
    const TString treeName = "slimtree";
    const TString clusterNode = "CLUSTERINFO_CEMC";

    const float eff_dR = 0.1;
    const float truth_iso_max = 4.0;
    const float eta_max = 0.7;
    const float pT_min = 8.0;
    const float pT_max = 36.0;
    const float vertex_cut = 60.0;

    TFile *fout = new TFile(outFile, "RECREATE");

    // ========================================================================
    // Histogram definitions
    // ========================================================================
    // Labels: "single" vs "double"
    TString labels[2] = {"single", "double"};

    // 1. min_dR distributions for all truth photons
    TH1F *h_mindR[2];
    TH1F *h_mindR_matched[2];
    TH1F *h_mindR_unmatched[2];
    // Fine binning near threshold
    TH1F *h_mindR_zoom[2];
    TH1F *h_mindR_unmatched_zoom[2];

    // 2. 2D: min_dR vs truth_pT
    TH2F *h_mindR_vs_pT[2];
    TH2F *h_mindR_vs_pT_unmatched[2];

    // 3. Signed delta_eta (truth - nearest cluster)
    TH1F *h_deta_all[2];
    TH1F *h_deta_matched[2];
    TH1F *h_deta_unmatched[2];

    // 4. delta_phi
    TH1F *h_dphi_all[2];

    // 5. Vertex distributions
    TH1F *h_vertexz_reco[2];
    TH1F *h_vertexz_truth[2];
    TH1F *h_vertex_diff[2];  // reco - truth
    TH2F *h_vertex_2d[2];    // reco vs truth

    // 6. Unmatched classification: "near miss" (0.1 < dR < 0.2) vs "far" (dR > 0.5) vs "no cluster"
    TH1F *h_unmatched_category[2];  // 4 bins

    // 7. Number of truth photons per event (with cuts)
    TH1F *h_nphotons[2];

    // 8. min_dR for unmatched, binned by |vertex_diff|
    TH1F *h_mindR_unmatched_smallvtx[2];  // |vtx_diff| < 5 cm
    TH1F *h_mindR_unmatched_largevtx[2];  // |vtx_diff| > 5 cm

    // 9. deta vs vertex_diff for unmatched photons
    TH2F *h_deta_vs_vtxdiff_unmatched[2];

    // 10. Has-cluster flag: does the truth photon have ANY cluster with matching trkID?
    TH1F *h_has_trkid_match[2]; // 0=no cluster with matching trkid, 1=has cluster with matching trkid

    // 11. For truth photons that DO have a trkid-matched cluster, what is dR?
    TH1F *h_trkid_dR[2];
    TH1F *h_trkid_dR_unmatched[2]; // subset where dR > eff_dR

    // 12. Efficiency vs pT
    TH1F *h_pT_all[2];
    TH1F *h_pT_matched[2];

    for (int i = 0; i < 2; i++)
    {
        h_mindR[i] = new TH1F(Form("h_mindR_%s", labels[i].Data()),
                               Form("Min #DeltaR (truth to nearest cluster) [%s];min #DeltaR;Truth photons", labels[i].Data()),
                               200, 0, 2.0);
        h_mindR_matched[i] = new TH1F(Form("h_mindR_matched_%s", labels[i].Data()),
                                       Form("Min #DeltaR (matched, dR<%.2f) [%s];min #DeltaR;Truth photons", eff_dR, labels[i].Data()),
                                       100, 0, eff_dR);
        h_mindR_unmatched[i] = new TH1F(Form("h_mindR_unmatched_%s", labels[i].Data()),
                                         Form("Min #DeltaR (unmatched, dR>%.2f) [%s];min #DeltaR;Truth photons", eff_dR, labels[i].Data()),
                                         200, 0, 2.0);

        h_mindR_zoom[i] = new TH1F(Form("h_mindR_zoom_%s", labels[i].Data()),
                                     Form("Min #DeltaR zoom [%s];min #DeltaR;Truth photons", labels[i].Data()),
                                     100, 0, 0.5);
        h_mindR_unmatched_zoom[i] = new TH1F(Form("h_mindR_unmatched_zoom_%s", labels[i].Data()),
                                              Form("Min #DeltaR (unmatched) zoom [%s];min #DeltaR;Truth photons", labels[i].Data()),
                                              100, 0, 0.5);

        h_mindR_vs_pT[i] = new TH2F(Form("h_mindR_vs_pT_%s", labels[i].Data()),
                                      Form("Min #DeltaR vs truth p_{T} [%s];truth p_{T} [GeV];min #DeltaR", labels[i].Data()),
                                      28, 8, 36, 200, 0, 2.0);
        h_mindR_vs_pT_unmatched[i] = new TH2F(Form("h_mindR_vs_pT_unmatched_%s", labels[i].Data()),
                                                Form("Min #DeltaR vs truth p_{T} (unmatched) [%s];truth p_{T} [GeV];min #DeltaR", labels[i].Data()),
                                                28, 8, 36, 200, 0, 2.0);

        h_deta_all[i] = new TH1F(Form("h_deta_all_%s", labels[i].Data()),
                                   Form("#Delta#eta (cluster - truth) [%s];#Delta#eta;Truth photons", labels[i].Data()),
                                   200, -0.5, 0.5);
        h_deta_matched[i] = new TH1F(Form("h_deta_matched_%s", labels[i].Data()),
                                       Form("#Delta#eta (cluster - truth, matched) [%s];#Delta#eta;Truth photons", labels[i].Data()),
                                       200, -0.5, 0.5);
        h_deta_unmatched[i] = new TH1F(Form("h_deta_unmatched_%s", labels[i].Data()),
                                         Form("#Delta#eta (cluster - truth, unmatched) [%s];#Delta#eta;Truth photons", labels[i].Data()),
                                         200, -0.5, 0.5);

        h_dphi_all[i] = new TH1F(Form("h_dphi_all_%s", labels[i].Data()),
                                   Form("#Delta#phi (cluster - truth) [%s];#Delta#phi;Truth photons", labels[i].Data()),
                                   200, -0.5, 0.5);

        h_vertexz_reco[i] = new TH1F(Form("h_vertexz_reco_%s", labels[i].Data()),
                                       Form("Reco vertex z [%s];vertex z [cm];Events", labels[i].Data()),
                                       200, -100, 100);
        h_vertexz_truth[i] = new TH1F(Form("h_vertexz_truth_%s", labels[i].Data()),
                                        Form("Truth vertex z [%s];vertex z [cm];Events", labels[i].Data()),
                                        200, -100, 100);
        h_vertex_diff[i] = new TH1F(Form("h_vertex_diff_%s", labels[i].Data()),
                                      Form("Reco - Truth vertex z [%s];#Deltavtx_{z} [cm];Events", labels[i].Data()),
                                      200, -50, 50);
        h_vertex_2d[i] = new TH2F(Form("h_vertex_2d_%s", labels[i].Data()),
                                    Form("Reco vs Truth vertex z [%s];Truth vtx_{z} [cm];Reco vtx_{z} [cm]", labels[i].Data()),
                                    100, -80, 80, 100, -80, 80);

        h_unmatched_category[i] = new TH1F(Form("h_unmatched_category_%s", labels[i].Data()),
                                            Form("Unmatched photon category [%s];;Count", labels[i].Data()),
                                            4, 0, 4);
        h_unmatched_category[i]->GetXaxis()->SetBinLabel(1, "Near miss (0.1<dR<0.2)");
        h_unmatched_category[i]->GetXaxis()->SetBinLabel(2, "Medium (0.2<dR<0.5)");
        h_unmatched_category[i]->GetXaxis()->SetBinLabel(3, "Far (dR>0.5)");
        h_unmatched_category[i]->GetXaxis()->SetBinLabel(4, "No cluster (trkID)");

        h_nphotons[i] = new TH1F(Form("h_nphotons_%s", labels[i].Data()),
                                   Form("Truth photons per event [%s];N truth photons;Events", labels[i].Data()),
                                   10, 0, 10);

        h_mindR_unmatched_smallvtx[i] = new TH1F(Form("h_mindR_unmatched_smallvtx_%s", labels[i].Data()),
                                                   Form("Min #DeltaR (unmatched, |#Deltavtx|<5cm) [%s];min #DeltaR;Truth photons", labels[i].Data()),
                                                   200, 0, 2.0);
        h_mindR_unmatched_largevtx[i] = new TH1F(Form("h_mindR_unmatched_largevtx_%s", labels[i].Data()),
                                                   Form("Min #DeltaR (unmatched, |#Deltavtx|>5cm) [%s];min #DeltaR;Truth photons", labels[i].Data()),
                                                   200, 0, 2.0);

        h_deta_vs_vtxdiff_unmatched[i] = new TH2F(Form("h_deta_vs_vtxdiff_unmatched_%s", labels[i].Data()),
                                                    Form("#Delta#eta vs #Deltavtx_{z} (unmatched) [%s];Reco-Truth vtx_{z} [cm];#Delta#eta (cluster-truth)", labels[i].Data()),
                                                    100, -50, 50, 100, -0.5, 0.5);

        h_has_trkid_match[i] = new TH1F(Form("h_has_trkid_match_%s", labels[i].Data()),
                                          Form("Truth photon has trkID-matched cluster [%s];;Count", labels[i].Data()),
                                          2, 0, 2);
        h_has_trkid_match[i]->GetXaxis()->SetBinLabel(1, "No trkID match");
        h_has_trkid_match[i]->GetXaxis()->SetBinLabel(2, "Has trkID match");

        h_trkid_dR[i] = new TH1F(Form("h_trkid_dR_%s", labels[i].Data()),
                                   Form("#DeltaR for trkID-matched cluster [%s];#DeltaR;Truth photons", labels[i].Data()),
                                   200, 0, 2.0);
        h_trkid_dR_unmatched[i] = new TH1F(Form("h_trkid_dR_unmatched_%s", labels[i].Data()),
                                             Form("#DeltaR for trkID-matched cluster (dR>%.2f) [%s];#DeltaR;Truth photons", eff_dR, labels[i].Data()),
                                             200, 0, 2.0);

        h_pT_all[i] = new TH1F(Form("h_pT_all_%s", labels[i].Data()),
                                 Form("Truth p_{T} (all) [%s];p_{T} [GeV];Truth photons", labels[i].Data()),
                                 28, 8, 36);
        h_pT_matched[i] = new TH1F(Form("h_pT_matched_%s", labels[i].Data()),
                                     Form("Truth p_{T} (matched) [%s];p_{T} [GeV];Truth photons", labels[i].Data()),
                                     28, 8, 36);
    }

    // ========================================================================
    // Process each sample
    // ========================================================================
    TString files[2] = {singleFile, doubleFile};

    for (int isample = 0; isample < 2; isample++)
    {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Processing " << labels[isample] << " MC: " << files[isample] << std::endl;
        std::cout << "========================================" << std::endl;

        TFile *fin = TFile::Open(files[isample]);
        if (!fin || fin->IsZombie())
        {
            std::cerr << "Cannot open " << files[isample] << std::endl;
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
        TTreeReaderArray<int> particle_photonclass(reader, "particle_photonclass");
        TTreeReaderArray<int> particle_trkid(reader, "particle_trkid");
        TTreeReaderArray<float> particle_truth_iso_03(reader, "particle_truth_iso_03");

        TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusterNode.Data()));
        TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusterNode.Data()));
        TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusterNode.Data()));
        TTreeReaderArray<int> cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusterNode.Data()));

        // Counters
        long long nEvents = 0;
        long long nTruthPhotons = 0;
        long long nMatched_dR = 0;
        long long nUnmatched_dR = 0;
        long long nNoTrkIDCluster = 0;
        long long nHasTrkIDCluster = 0;
        long long nHasTrkID_passdR = 0;
        long long nHasTrkID_faildR = 0;
        long long nNearMiss = 0;   // 0.1 < dR < 0.2
        long long nMedium = 0;     // 0.2 < dR < 0.5
        long long nFar = 0;        // dR > 0.5

        while (reader.Next())
        {
            if (nEvents >= maxEvents) break;
            nEvents++;

            if (nEvents % 20000 == 0)
            {
                std::cout << "  Event " << nEvents << " / " << maxEvents << std::endl;
            }

            float vtx_reco = *vertexz;
            float vtx_truth = *vertexz_truth;
            float vtx_diff = vtx_reco - vtx_truth;

            // Vertex cut (on reco vertex, same as in RecoEffCalculator)
            if (fabs(vtx_reco) > vertex_cut) continue;

            h_vertexz_reco[isample]->Fill(vtx_reco);
            h_vertexz_truth[isample]->Fill(vtx_truth);
            h_vertex_diff[isample]->Fill(vtx_diff);
            h_vertex_2d[isample]->Fill(vtx_truth, vtx_reco);

            int ncl = *nclusters;
            int npart = *nparticles;

            // Build trkid -> cluster index map (for fast lookup)
            std::map<int, std::vector<int>> trkid_to_clusters;
            for (int icl = 0; icl < ncl; icl++)
            {
                trkid_to_clusters[cluster_truthtrkID[icl]].push_back(icl);
            }

            // Count and process truth photons
            int nPhotonsThisEvent = 0;

            for (int ipart = 0; ipart < npart; ipart++)
            {
                // Same selection as RecoEffCalculator
                if (particle_pid[ipart] != 22) continue;
                if (particle_photonclass[ipart] >= 3) continue;  // direct or fragmentation only
                if (particle_truth_iso_03[ipart] >= truth_iso_max) continue;
                if (fabs(particle_Eta[ipart]) >= eta_max) continue;
                if (particle_Pt[ipart] < pT_min || particle_Pt[ipart] >= pT_max) continue;

                nTruthPhotons++;
                nPhotonsThisEvent++;

                float truth_eta = particle_Eta[ipart];
                float truth_phi = particle_Phi[ipart];
                float truth_pT = particle_Pt[ipart];
                int truth_trkid = particle_trkid[ipart];

                // --- Method 1: Find ANY nearest cluster (minimum dR to all clusters) ---
                float min_dR_any = 999;
                float best_deta_any = 999;
                float best_dphi_any = 999;
                int best_cl_any = -1;

                for (int icl = 0; icl < ncl; icl++)
                {
                    if (cluster_Et[icl] < 5.0) continue;  // reco_min_ET cut
                    float deta = cluster_Eta[icl] - truth_eta;
                    float dphi = cluster_Phi[icl] - truth_phi;
                    if (dphi > TMath::Pi()) dphi -= 2 * TMath::Pi();
                    else if (dphi < -TMath::Pi()) dphi += 2 * TMath::Pi();
                    float dR = sqrt(deta * deta + dphi * dphi);
                    if (dR < min_dR_any)
                    {
                        min_dR_any = dR;
                        best_deta_any = deta;
                        best_dphi_any = dphi;
                        best_cl_any = icl;
                    }
                }

                // --- Method 2: trkID-based matching (same as RecoEffCalculator) ---
                bool has_trkid_cluster = false;
                float trkid_dR = 999;
                float trkid_deta = 999;
                float trkid_dphi = 999;

                if (trkid_to_clusters.find(truth_trkid) != trkid_to_clusters.end())
                {
                    has_trkid_cluster = true;
                    // Find the closest cluster with this trkid
                    for (int icl_idx : trkid_to_clusters[truth_trkid])
                    {
                        float deta = cluster_Eta[icl_idx] - truth_eta;
                        float dphi = cluster_Phi[icl_idx] - truth_phi;
                        if (dphi > TMath::Pi()) dphi -= 2 * TMath::Pi();
                        else if (dphi < -TMath::Pi()) dphi += 2 * TMath::Pi();
                        float dR = sqrt(deta * deta + dphi * dphi);
                        if (dR < trkid_dR)
                        {
                            trkid_dR = dR;
                            trkid_deta = deta;
                            trkid_dphi = dphi;
                        }
                    }
                }

                // Fill histograms
                h_mindR[isample]->Fill(min_dR_any);
                h_mindR_zoom[isample]->Fill(min_dR_any);
                h_mindR_vs_pT[isample]->Fill(truth_pT, min_dR_any);
                h_pT_all[isample]->Fill(truth_pT);

                if (best_cl_any >= 0)
                {
                    h_deta_all[isample]->Fill(best_deta_any);
                    h_dphi_all[isample]->Fill(best_dphi_any);
                }

                // trkID matching diagnostics
                if (has_trkid_cluster)
                {
                    nHasTrkIDCluster++;
                    h_has_trkid_match[isample]->Fill(1);
                    h_trkid_dR[isample]->Fill(trkid_dR);

                    if (trkid_dR < eff_dR)
                    {
                        // Matched! (same logic as RecoEffCalculator)
                        nMatched_dR++;
                        nHasTrkID_passdR++;
                        h_mindR_matched[isample]->Fill(min_dR_any);
                        h_deta_matched[isample]->Fill(trkid_deta);
                        h_pT_matched[isample]->Fill(truth_pT);
                    }
                    else
                    {
                        // Has trkID cluster but dR too large (vertex shift effect!)
                        nUnmatched_dR++;
                        nHasTrkID_faildR++;
                        h_mindR_unmatched[isample]->Fill(min_dR_any);
                        h_mindR_unmatched_zoom[isample]->Fill(min_dR_any);
                        h_mindR_vs_pT_unmatched[isample]->Fill(truth_pT, min_dR_any);
                        h_deta_unmatched[isample]->Fill(trkid_deta);
                        h_trkid_dR_unmatched[isample]->Fill(trkid_dR);
                        h_deta_vs_vtxdiff_unmatched[isample]->Fill(vtx_diff, trkid_deta);

                        // Categorize
                        if (trkid_dR < 0.2) { nNearMiss++; h_unmatched_category[isample]->Fill(0); }
                        else if (trkid_dR < 0.5) { nMedium++; h_unmatched_category[isample]->Fill(1); }
                        else { nFar++; h_unmatched_category[isample]->Fill(2); }

                        if (fabs(vtx_diff) < 5.0)
                            h_mindR_unmatched_smallvtx[isample]->Fill(trkid_dR);
                        else
                            h_mindR_unmatched_largevtx[isample]->Fill(trkid_dR);
                    }
                }
                else
                {
                    // No cluster with matching trkID at all
                    nNoTrkIDCluster++;
                    nUnmatched_dR++;
                    h_has_trkid_match[isample]->Fill(0);
                    h_mindR_unmatched[isample]->Fill(min_dR_any);
                    h_mindR_unmatched_zoom[isample]->Fill(min_dR_any);
                    h_mindR_vs_pT_unmatched[isample]->Fill(truth_pT, min_dR_any);
                    h_unmatched_category[isample]->Fill(3);  // "No cluster (trkID)"
                    if (best_cl_any >= 0)
                    {
                        h_deta_unmatched[isample]->Fill(best_deta_any);
                    }

                    if (fabs(vtx_diff) < 5.0)
                        h_mindR_unmatched_smallvtx[isample]->Fill(min_dR_any);
                    else
                        h_mindR_unmatched_largevtx[isample]->Fill(min_dR_any);
                }

            } // end particle loop

            h_nphotons[isample]->Fill(nPhotonsThisEvent);

        } // end event loop

        fin->Close();

        // Print summary
        std::cout << "\n--- Summary for " << labels[isample] << " MC ---" << std::endl;
        std::cout << "Events processed: " << nEvents << std::endl;
        std::cout << "Truth photons (pid==22, class<3, iso<4, |eta|<0.7, 8<pT<36): " << nTruthPhotons << std::endl;
        std::cout << "  Matched (trkID exists + dR < " << eff_dR << "): " << nMatched_dR
                  << " (" << 100.0 * nMatched_dR / nTruthPhotons << "%)" << std::endl;
        std::cout << "  Unmatched total: " << nUnmatched_dR
                  << " (" << 100.0 * nUnmatched_dR / nTruthPhotons << "%)" << std::endl;
        std::cout << std::endl;
        std::cout << "  Breakdown of unmatched:" << std::endl;
        std::cout << "    No trkID-matched cluster at all: " << nNoTrkIDCluster
                  << " (" << 100.0 * nNoTrkIDCluster / nTruthPhotons << "% of all, "
                  << 100.0 * nNoTrkIDCluster / std::max((long long)1, nUnmatched_dR) << "% of unmatched)" << std::endl;
        std::cout << "    Has trkID cluster but dR > " << eff_dR << ": " << nHasTrkID_faildR
                  << " (" << 100.0 * nHasTrkID_faildR / nTruthPhotons << "% of all, "
                  << 100.0 * nHasTrkID_faildR / std::max((long long)1, nUnmatched_dR) << "% of unmatched)" << std::endl;
        std::cout << std::endl;
        std::cout << "  Of those with trkID but dR too large:" << std::endl;
        std::cout << "    Near miss (0.1 < dR < 0.2): " << nNearMiss << std::endl;
        std::cout << "    Medium   (0.2 < dR < 0.5): " << nMedium << std::endl;
        std::cout << "    Far      (dR > 0.5):        " << nFar << std::endl;
        std::cout << std::endl;
        std::cout << "  Vertex info:" << std::endl;
        std::cout << "    Mean reco vtx_z: " << h_vertexz_reco[isample]->GetMean() << " cm" << std::endl;
        std::cout << "    Mean truth vtx_z: " << h_vertexz_truth[isample]->GetMean() << " cm" << std::endl;
        std::cout << "    Mean (reco - truth): " << h_vertex_diff[isample]->GetMean()
                  << " +/- " << h_vertex_diff[isample]->GetRMS() << " cm" << std::endl;
        std::cout << "  Reco efficiency: " << 100.0 * nMatched_dR / nTruthPhotons << "%" << std::endl;
    }

    // ========================================================================
    // Save all histograms
    // ========================================================================
    fout->cd();
    for (int i = 0; i < 2; i++)
    {
        h_mindR[i]->Write();
        h_mindR_matched[i]->Write();
        h_mindR_unmatched[i]->Write();
        h_mindR_zoom[i]->Write();
        h_mindR_unmatched_zoom[i]->Write();
        h_mindR_vs_pT[i]->Write();
        h_mindR_vs_pT_unmatched[i]->Write();
        h_deta_all[i]->Write();
        h_deta_matched[i]->Write();
        h_deta_unmatched[i]->Write();
        h_dphi_all[i]->Write();
        h_vertexz_reco[i]->Write();
        h_vertexz_truth[i]->Write();
        h_vertex_diff[i]->Write();
        h_vertex_2d[i]->Write();
        h_unmatched_category[i]->Write();
        h_nphotons[i]->Write();
        h_mindR_unmatched_smallvtx[i]->Write();
        h_mindR_unmatched_largevtx[i]->Write();
        h_deta_vs_vtxdiff_unmatched[i]->Write();
        h_has_trkid_match[i]->Write();
        h_trkid_dR[i]->Write();
        h_trkid_dR_unmatched[i]->Write();
        h_pT_all[i]->Write();
        h_pT_matched[i]->Write();
    }

    fout->Close();
    std::cout << "\nOutput saved to: " << outFile << std::endl;
}
