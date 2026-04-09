// investigate_vertex_shift.C
// Quantifies the vertex shift effect in double-interaction MC on cluster kinematics
// and truth-reco matching efficiency.
//
// Usage: root -l -b -q 'investigate_vertex_shift.C()'

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLine.h>
#include <TProfile.h>
#include <iostream>
#include <vector>
#include <cmath>

void investigate_vertex_shift() {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // --- Config ---
    const char* singleFile = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root";
    const char* doubleFile = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_double/bdt_split.root";
    const char* outFile    = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/vertex_shift_study.root";

    const double R_CEMC    = 93.5;  // cm, EMCal radius
    const double eff_dR    = 0.1;   // matching threshold from config
    const double vertex_cut = 60.0;
    const double eta_min   = -0.7;
    const double eta_max   =  0.7;
    const double ET_min    =  8.0;
    const int    maxEvents = 100000;
    const TString node     = "CLUSTERINFO_CEMC";

    // --- Output file ---
    TFile *fout = new TFile(outFile, "RECREATE");

    // ============================================================
    // PART 1: Vertex distributions
    // ============================================================
    cout << "=== PART 1: Vertex distributions ===" << endl;

    // Histograms
    TH1F *h_vtxz_single      = new TH1F("h_vtxz_single",      ";vertex z [cm];Events", 200, -100, 100);
    TH1F *h_vtxz_double      = new TH1F("h_vtxz_double",      ";vertex z [cm];Events", 200, -100, 100);
    TH1F *h_vtxz_truth_single = new TH1F("h_vtxz_truth_single", ";truth vertex z [cm];Events", 200, -100, 100);
    TH1F *h_vtxz_truth_double = new TH1F("h_vtxz_truth_double", ";truth vertex z [cm];Events", 200, -100, 100);
    TH1F *h_dz_single        = new TH1F("h_dz_single",        ";#Deltaz = vtxz - vtxz_{truth} [cm];Events", 200, -80, 80);
    TH1F *h_dz_double        = new TH1F("h_dz_double",        ";#Deltaz = vtxz - vtxz_{truth} [cm];Events", 200, -80, 80);
    TH1F *h_abs_dz_double    = new TH1F("h_abs_dz_double",    ";|#Deltaz| [cm];Events", 200, 0, 80);

    h_vtxz_single->SetLineColor(kBlue);
    h_vtxz_double->SetLineColor(kRed);
    h_vtxz_truth_single->SetLineColor(kBlue);
    h_vtxz_truth_double->SetLineColor(kRed);
    h_dz_single->SetLineColor(kBlue);
    h_dz_double->SetLineColor(kRed);

    // Fill single MC vertex
    {
        TFile *f = TFile::Open(singleFile);
        TTree *t = (TTree*)f->Get("slimtree");
        Float_t vertexz, vertexz_truth;
        t->SetBranchAddress("vertexz", &vertexz);
        t->SetBranchAddress("vertexz_truth", &vertexz_truth);
        Long64_t nEntries = TMath::Min((Long64_t)maxEvents, t->GetEntries());
        for (Long64_t i = 0; i < nEntries; i++) {
            t->GetEntry(i);
            h_vtxz_single->Fill(vertexz);
            h_vtxz_truth_single->Fill(vertexz_truth);
            h_dz_single->Fill(vertexz - vertexz_truth);
        }
        f->Close();
    }

    // Fill double MC vertex
    {
        TFile *f = TFile::Open(doubleFile);
        TTree *t = (TTree*)f->Get("slimtree");
        Float_t vertexz, vertexz_truth;
        t->SetBranchAddress("vertexz", &vertexz);
        t->SetBranchAddress("vertexz_truth", &vertexz_truth);
        Long64_t nEntries = TMath::Min((Long64_t)maxEvents, t->GetEntries());
        for (Long64_t i = 0; i < nEntries; i++) {
            t->GetEntry(i);
            h_vtxz_double->Fill(vertexz);
            h_vtxz_truth_double->Fill(vertexz_truth);
            h_dz_double->Fill(vertexz - vertexz_truth);
            h_abs_dz_double->Fill(fabs(vertexz - vertexz_truth));
        }
        f->Close();
    }

    // Print vertex shift statistics
    cout << "Single MC:  <vtxz> = " << h_vtxz_single->GetMean()
         << " +/- " << h_vtxz_single->GetRMS() << " cm" << endl;
    cout << "Double MC:  <vtxz> = " << h_vtxz_double->GetMean()
         << " +/- " << h_vtxz_double->GetRMS() << " cm" << endl;
    cout << "Single MC:  <vtxz_truth> = " << h_vtxz_truth_single->GetMean()
         << " +/- " << h_vtxz_truth_single->GetRMS() << " cm" << endl;
    cout << "Double MC:  <vtxz_truth> = " << h_vtxz_truth_double->GetMean()
         << " +/- " << h_vtxz_truth_double->GetRMS() << " cm" << endl;
    cout << "Single MC:  <dz> = " << h_dz_single->GetMean()
         << " +/- " << h_dz_single->GetRMS() << " cm" << endl;
    cout << "Double MC:  <dz> = " << h_dz_double->GetMean()
         << " +/- " << h_dz_double->GetRMS() << " cm" << endl;
    cout << "Double MC:  <|dz|> = " << h_abs_dz_double->GetMean()
         << " +/- " << h_abs_dz_double->GetRMS() << " cm" << endl;
    cout << "Double MC:  median |dz| = " << h_abs_dz_double->GetQuantiles(1, new double[1], new double[1]{0.5})
         << endl;

    // Get proper median
    double median_val;
    double prob = 0.5;
    h_abs_dz_double->GetQuantiles(1, &median_val, &prob);
    cout << "Double MC:  median |dz| = " << median_val << " cm" << endl;

    // ============================================================
    // PART 2: Expected eta shift (analytical)
    // ============================================================
    cout << endl << "=== PART 2: Expected eta shift ===" << endl;

    // For a cluster at position z_cluster on CEMC (radius R),
    // the true eta = asinh((z_cluster - vtxz_truth) / R)
    // the reco eta = asinh((z_cluster - vtxz_reco) / R)
    // The shift dz = vtxz_reco - vtxz_truth shifts z_rel by -dz
    //
    // Exact: delta_eta = asinh((z_rel + dz)/R) - asinh(z_rel/R)
    //   where z_rel = R * sinh(eta_true), so:
    //   delta_eta = asinh(sinh(eta_true) + dz/R) - eta_true
    //
    // Small dz approximation: delta_eta ~ -dz / (R * cosh(eta))
    //   (Note: the cluster z is fixed, but vtx moves, so z_rel decreases by dz,
    //    hence delta_eta ~ -dz / (R * cosh(eta)))

    // Analytical curves: delta_eta vs dz for several eta values
    const int nEtaVals = 5;
    double etaVals[nEtaVals] = {0.0, 0.2, 0.4, 0.6, 0.7};

    TH2F *h_frame_deta = new TH2F("h_frame_deta",
        ";#Deltaz = vtxz_{reco} - vtxz_{truth} [cm];#Delta#eta (reco - truth)",
        100, -60, 60, 100, -0.08, 0.08);

    // Create TGraphs for exact formula
    const int nDz = 121;
    double dzArr[nDz];
    for (int i = 0; i < nDz; i++) dzArr[i] = -60.0 + i * 1.0;

    int colors[nEtaVals] = {kBlack, kBlue, kGreen+2, kRed, kMagenta};
    TGraph *g_deta_exact[nEtaVals];
    TGraph *g_deta_approx[nEtaVals];

    for (int ie = 0; ie < nEtaVals; ie++) {
        double eta0 = etaVals[ie];
        double z_rel_true = R_CEMC * sinh(eta0);

        g_deta_exact[ie]  = new TGraph(nDz);
        g_deta_approx[ie] = new TGraph(nDz);

        for (int id = 0; id < nDz; id++) {
            double dz = dzArr[id];
            // Vertex shifts by +dz, so z_rel = z_cluster - vtxz_reco = (z_cluster - vtxz_truth) - dz
            double z_rel_reco = z_rel_true - dz;
            double eta_reco = asinh(z_rel_reco / R_CEMC);
            double delta_eta_exact = eta_reco - eta0;
            double delta_eta_approx = -dz / (R_CEMC * cosh(eta0));

            g_deta_exact[ie]->SetPoint(id, dz, delta_eta_exact);
            g_deta_approx[ie]->SetPoint(id, dz, delta_eta_approx);
        }

        g_deta_exact[ie]->SetLineColor(colors[ie]);
        g_deta_exact[ie]->SetLineWidth(2);
        g_deta_exact[ie]->SetLineStyle(1);

        g_deta_approx[ie]->SetLineColor(colors[ie]);
        g_deta_approx[ie]->SetLineWidth(2);
        g_deta_approx[ie]->SetLineStyle(2);

        char name_exact[64], name_approx[64];
        sprintf(name_exact, "g_deta_exact_eta%.1f", eta0);
        sprintf(name_approx, "g_deta_approx_eta%.1f", eta0);
        g_deta_exact[ie]->SetName(name_exact);
        g_deta_approx[ie]->SetName(name_approx);
    }

    // Print characteristic delta_eta for typical dz values
    cout << "delta_eta for typical vertex shifts (exact formula):" << endl;
    cout << "  eta=0.0: dz=10cm -> delta_eta = " << asinh((R_CEMC*sinh(0.0)-10)/R_CEMC) - 0.0 << endl;
    cout << "  eta=0.0: dz=20cm -> delta_eta = " << asinh((R_CEMC*sinh(0.0)-20)/R_CEMC) - 0.0 << endl;
    cout << "  eta=0.0: dz=30cm -> delta_eta = " << asinh((R_CEMC*sinh(0.0)-30)/R_CEMC) - 0.0 << endl;
    cout << "  eta=0.4: dz=10cm -> delta_eta = " << asinh((R_CEMC*sinh(0.4)-10)/R_CEMC) - 0.4 << endl;
    cout << "  eta=0.4: dz=20cm -> delta_eta = " << asinh((R_CEMC*sinh(0.4)-20)/R_CEMC) - 0.4 << endl;
    cout << "  eta=0.4: dz=30cm -> delta_eta = " << asinh((R_CEMC*sinh(0.4)-30)/R_CEMC) - 0.4 << endl;
    cout << "  eta=0.7: dz=10cm -> delta_eta = " << asinh((R_CEMC*sinh(0.7)-10)/R_CEMC) - 0.7 << endl;
    cout << "  eta=0.7: dz=20cm -> delta_eta = " << asinh((R_CEMC*sinh(0.7)-20)/R_CEMC) - 0.7 << endl;

    // What dz gives delta_eta = dR_threshold at eta=0?
    // delta_eta = -dz/R at eta=0 => dz_max = R * dR_threshold
    double dz_for_dR_at_eta0 = R_CEMC * eff_dR;
    cout << endl << "dz that produces |delta_eta| = " << eff_dR
         << " at eta=0: " << dz_for_dR_at_eta0 << " cm" << endl;


    // ============================================================
    // PART 3: Actual cluster eta and ET comparisons
    // ============================================================
    cout << endl << "=== PART 3: Cluster eta/ET distributions ===" << endl;

    TH1F *h_clusEta_single = new TH1F("h_clusEta_single", ";cluster #eta;Clusters (norm)", 140, -0.8, 0.8);
    TH1F *h_clusEta_double = new TH1F("h_clusEta_double", ";cluster #eta;Clusters (norm)", 140, -0.8, 0.8);
    TH1F *h_clusET_single  = new TH1F("h_clusET_single",  ";cluster E_{T} [GeV];Clusters (norm)", 100, 5, 40);
    TH1F *h_clusET_double  = new TH1F("h_clusET_double",  ";cluster E_{T} [GeV];Clusters (norm)", 100, 5, 40);

    // Truth-reco delta_eta in data
    TH1F *h_deta_single  = new TH1F("h_deta_single", ";#Delta#eta (reco cluster - truth particle);Matched pairs (norm)", 200, -0.15, 0.15);
    TH1F *h_deta_double  = new TH1F("h_deta_double", ";#Delta#eta (reco cluster - truth particle);Matched pairs (norm)", 200, -0.15, 0.15);
    TH1F *h_dphi_single  = new TH1F("h_dphi_single", ";#Delta#phi (reco cluster - truth particle);Matched pairs (norm)", 200, -0.15, 0.15);
    TH1F *h_dphi_double  = new TH1F("h_dphi_double", ";#Delta#phi (reco cluster - truth particle);Matched pairs (norm)", 200, -0.15, 0.15);
    TH1F *h_dR_single    = new TH1F("h_dR_single",   ";#DeltaR (reco cluster - truth particle);Matched pairs (norm)", 200, 0, 0.3);
    TH1F *h_dR_double    = new TH1F("h_dR_double",   ";#DeltaR (reco cluster - truth particle);Matched pairs (norm)", 200, 0, 0.3);

    // delta_eta vs |dz|
    TH2F *h2_deta_vs_dz_double = new TH2F("h2_deta_vs_dz_double",
        ";|#Deltaz| [cm];#Delta#eta (reco - truth)", 60, 0, 60, 100, -0.15, 0.15);
    TProfile *hp_deta_vs_dz_double = new TProfile("hp_deta_vs_dz_double",
        ";|#Deltaz| [cm];<|#Delta#eta|>", 60, 0, 60, 0, 0.15);

    // Matching efficiency vs dz
    TH1F *h_dz_all_double    = new TH1F("h_dz_all_double",    ";|#Deltaz| [cm];Signal photons", 60, 0, 60);
    TH1F *h_dz_matched_double = new TH1F("h_dz_matched_double", ";|#Deltaz| [cm];Matched", 60, 0, 60);
    TH1F *h_dz_failed_double = new TH1F("h_dz_failed_double",  ";|#Deltaz| [cm];Failed match", 60, 0, 60);

    // Matching eff vs ET for single and double
    const int nPtBins = 12;
    double ptBins[nPtBins+1] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36};
    TH1F *h_truth_pt_single   = new TH1F("h_truth_pt_single",   ";truth p_{T} [GeV];", nPtBins, ptBins);
    TH1F *h_matched_pt_single = new TH1F("h_matched_pt_single", ";truth p_{T} [GeV];", nPtBins, ptBins);
    TH1F *h_truth_pt_double   = new TH1F("h_truth_pt_double",   ";truth p_{T} [GeV];", nPtBins, ptBins);
    TH1F *h_matched_pt_double = new TH1F("h_matched_pt_double", ";truth p_{T} [GeV];", nPtBins, ptBins);

    h_clusEta_single->SetLineColor(kBlue);
    h_clusEta_double->SetLineColor(kRed);
    h_clusET_single->SetLineColor(kBlue);
    h_clusET_double->SetLineColor(kRed);
    h_deta_single->SetLineColor(kBlue);
    h_deta_double->SetLineColor(kRed);
    h_dphi_single->SetLineColor(kBlue);
    h_dphi_double->SetLineColor(kRed);
    h_dR_single->SetLineColor(kBlue);
    h_dR_double->SetLineColor(kRed);

    // ============================================================
    // Process SINGLE MC
    // ============================================================
    cout << "Processing single MC..." << endl;
    {
        TFile *f = TFile::Open(singleFile);
        TTree *t = (TTree*)f->Get("slimtree");

        Float_t vertexz, vertexz_truth;
        Int_t ncluster, nparticles;
        Float_t cluster_Et[500], cluster_Eta[500], cluster_Phi[500];
        Int_t   cluster_truthtrkID[500], cluster_pid[500];
        Float_t particle_Pt[200], particle_Eta[200], particle_Phi[200];
        Int_t   particle_pid[200], particle_trkid[200], particle_photonclass[200];
        Float_t particle_truth_iso_03[200];

        t->SetBranchAddress("vertexz", &vertexz);
        t->SetBranchAddress("vertexz_truth", &vertexz_truth);
        TString nclusName = Form("ncluster_%s", node.Data());
        t->SetBranchAddress(nclusName, &ncluster);
        t->SetBranchAddress(Form("cluster_Et_%s", node.Data()), cluster_Et);
        t->SetBranchAddress(Form("cluster_Eta_%s", node.Data()), cluster_Eta);
        t->SetBranchAddress(Form("cluster_Phi_%s", node.Data()), cluster_Phi);
        t->SetBranchAddress(Form("cluster_truthtrkID_%s", node.Data()), cluster_truthtrkID);
        t->SetBranchAddress(Form("cluster_pid_%s", node.Data()), cluster_pid);
        t->SetBranchAddress("nparticles", &nparticles);
        t->SetBranchAddress("particle_Pt", particle_Pt);
        t->SetBranchAddress("particle_Eta", particle_Eta);
        t->SetBranchAddress("particle_Phi", particle_Phi);
        t->SetBranchAddress("particle_pid", particle_pid);
        t->SetBranchAddress("particle_trkid", particle_trkid);
        t->SetBranchAddress("particle_photonclass", particle_photonclass);
        t->SetBranchAddress("particle_truth_iso_03", particle_truth_iso_03);

        Long64_t nEntries = TMath::Min((Long64_t)maxEvents, t->GetEntries());
        int nMatched = 0, nSignalPhotons = 0;

        for (Long64_t i = 0; i < nEntries; i++) {
            t->GetEntry(i);
            if (fabs(vertexz) > vertex_cut) continue;

            // Fill cluster distributions (all clusters with ET > ET_min in eta range)
            for (int ic = 0; ic < ncluster && ic < 500; ic++) {
                if (cluster_Et[ic] > ET_min && fabs(cluster_Eta[ic]) < eta_max) {
                    h_clusEta_single->Fill(cluster_Eta[ic]);
                    h_clusET_single->Fill(cluster_Et[ic]);
                }
            }

            // Find signal photons (pid=1 or 2, photonclass=1 = direct photon, isolated)
            for (int ip = 0; ip < nparticles && ip < 200; ip++) {
                if (!(particle_pid[ip] == 22 || abs(particle_pid[ip]) == 11)) continue;
                // Use photonclass: 1 = direct photon
                if (particle_photonclass[ip] != 1) continue;
                if (fabs(particle_Eta[ip]) > eta_max) continue;
                if (particle_Pt[ip] < ET_min) continue;
                if (particle_truth_iso_03[ip] > 4.0) continue;

                nSignalPhotons++;
                h_truth_pt_single->Fill(particle_Pt[ip]);

                // Find best matching cluster by truthtrkID
                int best_ic = -1;
                double best_dR = 999;
                for (int ic = 0; ic < ncluster && ic < 500; ic++) {
                    if (cluster_truthtrkID[ic] != particle_trkid[ip]) continue;
                    if (cluster_Et[ic] < ET_min) continue;

                    double deta = cluster_Eta[ic] - particle_Eta[ip];
                    double dphi = cluster_Phi[ic] - particle_Phi[ip];
                    while (dphi >  TMath::Pi()) dphi -= 2*TMath::Pi();
                    while (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
                    double dR = sqrt(deta*deta + dphi*dphi);

                    if (dR < best_dR) {
                        best_dR = dR;
                        best_ic = ic;
                    }
                }

                if (best_ic >= 0) {
                    double deta = cluster_Eta[best_ic] - particle_Eta[ip];
                    double dphi = cluster_Phi[best_ic] - particle_Phi[ip];
                    while (dphi >  TMath::Pi()) dphi -= 2*TMath::Pi();
                    while (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();

                    h_deta_single->Fill(deta);
                    h_dphi_single->Fill(dphi);
                    h_dR_single->Fill(best_dR);

                    if (best_dR < eff_dR) {
                        nMatched++;
                        h_matched_pt_single->Fill(particle_Pt[ip]);
                    }
                }
            }
        }
        cout << "Single MC: signal photons = " << nSignalPhotons
             << ", matched (dR<" << eff_dR << ") = " << nMatched
             << ", match rate = " << (nSignalPhotons > 0 ? (double)nMatched/nSignalPhotons : 0) << endl;
        f->Close();
    }

    // ============================================================
    // Process DOUBLE MC
    // ============================================================
    cout << "Processing double MC..." << endl;
    {
        TFile *f = TFile::Open(doubleFile);
        TTree *t = (TTree*)f->Get("slimtree");

        Float_t vertexz, vertexz_truth;
        Int_t ncluster, nparticles;
        Float_t cluster_Et[500], cluster_Eta[500], cluster_Phi[500];
        Int_t   cluster_truthtrkID[500], cluster_pid[500], cluster_embed_id[500];
        Float_t particle_Pt[200], particle_Eta[200], particle_Phi[200];
        Int_t   particle_pid[200], particle_trkid[200], particle_photonclass[200];
        Float_t particle_truth_iso_03[200];

        t->SetBranchAddress("vertexz", &vertexz);
        t->SetBranchAddress("vertexz_truth", &vertexz_truth);
        TString nclusName = Form("ncluster_%s", node.Data());
        t->SetBranchAddress(nclusName, &ncluster);
        t->SetBranchAddress(Form("cluster_Et_%s", node.Data()), cluster_Et);
        t->SetBranchAddress(Form("cluster_Eta_%s", node.Data()), cluster_Eta);
        t->SetBranchAddress(Form("cluster_Phi_%s", node.Data()), cluster_Phi);
        t->SetBranchAddress(Form("cluster_truthtrkID_%s", node.Data()), cluster_truthtrkID);
        t->SetBranchAddress(Form("cluster_pid_%s", node.Data()), cluster_pid);
        t->SetBranchAddress(Form("cluster_embed_id_%s", node.Data()), cluster_embed_id);
        t->SetBranchAddress("nparticles", &nparticles);
        t->SetBranchAddress("particle_Pt", particle_Pt);
        t->SetBranchAddress("particle_Eta", particle_Eta);
        t->SetBranchAddress("particle_Phi", particle_Phi);
        t->SetBranchAddress("particle_pid", particle_pid);
        t->SetBranchAddress("particle_trkid", particle_trkid);
        t->SetBranchAddress("particle_photonclass", particle_photonclass);
        t->SetBranchAddress("particle_truth_iso_03", particle_truth_iso_03);

        Long64_t nEntries = TMath::Min((Long64_t)maxEvents, t->GetEntries());
        int nMatched = 0, nSignalPhotons = 0, nFailed = 0;

        for (Long64_t i = 0; i < nEntries; i++) {
            t->GetEntry(i);
            if (fabs(vertexz) > vertex_cut) continue;

            double dz = vertexz - vertexz_truth;

            // Fill cluster distributions
            for (int ic = 0; ic < ncluster && ic < 500; ic++) {
                if (cluster_Et[ic] > ET_min && fabs(cluster_Eta[ic]) < eta_max) {
                    h_clusEta_double->Fill(cluster_Eta[ic]);
                    h_clusET_double->Fill(cluster_Et[ic]);
                }
            }

            // Find signal photons
            for (int ip = 0; ip < nparticles && ip < 200; ip++) {
                if (!(particle_pid[ip] == 22 || abs(particle_pid[ip]) == 11)) continue;
                if (particle_photonclass[ip] != 1) continue;
                if (fabs(particle_Eta[ip]) > eta_max) continue;
                if (particle_Pt[ip] < ET_min) continue;
                if (particle_truth_iso_03[ip] > 4.0) continue;

                nSignalPhotons++;
                h_truth_pt_double->Fill(particle_Pt[ip]);
                h_dz_all_double->Fill(fabs(dz));

                // Find best matching cluster
                int best_ic = -1;
                double best_dR = 999;
                for (int ic = 0; ic < ncluster && ic < 500; ic++) {
                    if (cluster_truthtrkID[ic] != particle_trkid[ip]) continue;
                    if (cluster_Et[ic] < ET_min) continue;

                    double deta = cluster_Eta[ic] - particle_Eta[ip];
                    double dphi = cluster_Phi[ic] - particle_Phi[ip];
                    while (dphi >  TMath::Pi()) dphi -= 2*TMath::Pi();
                    while (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
                    double dRval = sqrt(deta*deta + dphi*dphi);

                    if (dRval < best_dR) {
                        best_dR = dRval;
                        best_ic = ic;
                    }
                }

                if (best_ic >= 0) {
                    double deta = cluster_Eta[best_ic] - particle_Eta[ip];
                    double dphi = cluster_Phi[best_ic] - particle_Phi[ip];
                    while (dphi >  TMath::Pi()) dphi -= 2*TMath::Pi();
                    while (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();

                    h_deta_double->Fill(deta);
                    h_dphi_double->Fill(dphi);
                    h_dR_double->Fill(best_dR);
                    h2_deta_vs_dz_double->Fill(fabs(dz), deta);
                    hp_deta_vs_dz_double->Fill(fabs(dz), fabs(deta));

                    if (best_dR < eff_dR) {
                        nMatched++;
                        h_matched_pt_double->Fill(particle_Pt[ip]);
                        h_dz_matched_double->Fill(fabs(dz));
                    } else {
                        nFailed++;
                        h_dz_failed_double->Fill(fabs(dz));
                    }
                } else {
                    // No truthtrkID match found at all
                    nFailed++;
                    h_dz_failed_double->Fill(fabs(dz));
                }
            }
        }
        cout << "Double MC: signal photons = " << nSignalPhotons
             << ", matched (dR<" << eff_dR << ") = " << nMatched
             << ", failed = " << nFailed
             << ", match rate = " << (nSignalPhotons > 0 ? (double)nMatched/nSignalPhotons : 0) << endl;
        f->Close();
    }

    // ============================================================
    // PART 4: Compute efficiencies and failure rates
    // ============================================================
    cout << endl << "=== PART 4: Matching efficiency vs pT ===" << endl;

    TH1F *h_eff_single = (TH1F*)h_matched_pt_single->Clone("h_eff_single");
    h_eff_single->Divide(h_truth_pt_single);
    TH1F *h_eff_double = (TH1F*)h_matched_pt_double->Clone("h_eff_double");
    h_eff_double->Divide(h_truth_pt_double);

    cout << "pT bin | single eff | double eff | ratio" << endl;
    for (int ib = 1; ib <= nPtBins; ib++) {
        double eff_s = h_eff_single->GetBinContent(ib);
        double eff_d = h_eff_double->GetBinContent(ib);
        double ratio = (eff_s > 0) ? eff_d / eff_s : 0;
        cout << Form("  [%4.0f, %4.0f] GeV: %.4f  %.4f  %.4f",
            ptBins[ib-1], ptBins[ib], eff_s, eff_d, ratio) << endl;
    }

    // Matching failure rate vs |dz| bins
    cout << endl << "=== PART 5: Failure rate vs |dz| ===" << endl;
    TH1F *h_fail_rate_vs_dz = (TH1F*)h_dz_all_double->Clone("h_fail_rate_vs_dz");
    h_fail_rate_vs_dz->SetTitle(";|#Deltaz| [cm];Matching failure rate");
    for (int ib = 1; ib <= h_fail_rate_vs_dz->GetNbinsX(); ib++) {
        double all = h_dz_all_double->GetBinContent(ib);
        double matched = h_dz_matched_double->GetBinContent(ib);
        double fail = (all > 0) ? 1.0 - matched/all : 0;
        double err = (all > 0) ? sqrt(fail*(1-fail)/all) : 0;
        h_fail_rate_vs_dz->SetBinContent(ib, fail);
        h_fail_rate_vs_dz->SetBinError(ib, err);
    }

    // Print in bins
    cout << "|dz| range | all | matched | fail rate" << endl;
    double dz_edges[] = {0, 5, 10, 15, 20, 30, 40, 60};
    for (int ie = 0; ie < 7; ie++) {
        int bin1 = h_dz_all_double->FindBin(dz_edges[ie] + 0.01);
        int bin2 = h_dz_all_double->FindBin(dz_edges[ie+1] - 0.01);
        double all = h_dz_all_double->Integral(bin1, bin2);
        double matched = h_dz_matched_double->Integral(bin1, bin2);
        double fail = (all > 0) ? 1.0 - matched/all : 0;
        cout << Form("  [%3.0f, %3.0f] cm: %6.0f  %6.0f  %.4f",
            dz_edges[ie], dz_edges[ie+1], all, matched, fail) << endl;
    }

    // ============================================================
    // PART 6: Fraction of double MC events with |dz| > threshold
    // ============================================================
    cout << endl << "=== PART 6: Fraction of events with large vertex shift ===" << endl;
    double total_dz = h_abs_dz_double->Integral();
    for (double thresh : {5.0, 9.35, 10.0, 15.0, 20.0, 30.0}) {
        int bin = h_abs_dz_double->FindBin(thresh);
        double above = h_abs_dz_double->Integral(bin, h_abs_dz_double->GetNbinsX()+1);
        cout << Form("  |dz| > %5.1f cm: %.1f%% of events", thresh, 100*above/total_dz) << endl;
    }

    // ============================================================
    // Write all histograms
    // ============================================================
    fout->cd();

    h_vtxz_single->Write();
    h_vtxz_double->Write();
    h_vtxz_truth_single->Write();
    h_vtxz_truth_double->Write();
    h_dz_single->Write();
    h_dz_double->Write();
    h_abs_dz_double->Write();

    h_frame_deta->Write();
    for (int ie = 0; ie < nEtaVals; ie++) {
        g_deta_exact[ie]->Write();
        g_deta_approx[ie]->Write();
    }

    h_clusEta_single->Write();
    h_clusEta_double->Write();
    h_clusET_single->Write();
    h_clusET_double->Write();

    h_deta_single->Write();
    h_deta_double->Write();
    h_dphi_single->Write();
    h_dphi_double->Write();
    h_dR_single->Write();
    h_dR_double->Write();

    h2_deta_vs_dz_double->Write();
    hp_deta_vs_dz_double->Write();

    h_dz_all_double->Write();
    h_dz_matched_double->Write();
    h_dz_failed_double->Write();
    h_fail_rate_vs_dz->Write();

    h_truth_pt_single->Write();
    h_matched_pt_single->Write();
    h_truth_pt_double->Write();
    h_matched_pt_double->Write();
    h_eff_single->Write();
    h_eff_double->Write();

    fout->Close();
    cout << endl << "Output written to " << outFile << endl;
}
