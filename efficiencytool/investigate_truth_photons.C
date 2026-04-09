// investigate_truth_photons.C
// Compare truth photon content between single (photon10) and double-interaction
// (photon10_double) MC to understand why the efficiency denominator is inflated
// in double-interaction samples.
//
// Usage: root -l -b -q 'investigate_truth_photons.C(100000)'

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <iostream>
#include <cmath>

void investigate_truth_photons(int maxEvents = 100000)
{
    gStyle->SetOptStat(0);

    // Input files
    const char *singleFile = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root";
    const char *doubleFile = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10_double/bdt_split.root";

    // Output
    TFile *fout = new TFile("results/investigate_truth_photons.root", "RECREATE");

    // Truth selection criteria (from RecoEffCalculator_TTreeReader.C):
    //   pid == 22
    //   photonclass < 3  (1=direct, 2=fragmentation; exclude 3=decay, 0=unknown)
    //   truth_iso_03 < 4.0 GeV
    //   |eta| < 0.7
    //   For photon10: max photon pT in [14, 30] (event-level cut)
    //   For photon10_double: max photon pT in [10, 100] (wider window)
    // We will apply the per-particle selection but NOT the event-level max photon pT
    // window, since we want to compare the raw truth content.
    // We will also show results WITH the max photon pT cut for reference.

    const float eta_min = -0.7;
    const float eta_max = 0.7;
    const float truthisocut = 4.0;   // GeV
    const int   conesize = 3;        // R=0.3

    // pT window for photon10 (single) event filter
    const float single_pt_lower = 14.0;
    const float single_pt_upper = 30.0;
    // pT window for photon10_double event filter
    const float double_pt_lower = 10.0;
    const float double_pt_upper = 100.0;

    // Histograms
    // --- Per-event counts ---
    TH1F *h_ntruth_single    = new TH1F("h_ntruth_single",    "N truth photons (denom) per event; N; Events", 10, -0.5, 9.5);
    TH1F *h_ntruth_double    = new TH1F("h_ntruth_double",    "N truth photons (denom) per event; N; Events", 10, -0.5, 9.5);
    TH1F *h_ntruth_double_ptcut = new TH1F("h_ntruth_double_ptcut", "N truth photons (denom, with pT window) per event; N; Events", 10, -0.5, 9.5);

    // --- Total nparticles ---
    TH1F *h_npart_single = new TH1F("h_npart_single", "nparticles per event; nparticles; Events", 100, 0, 500);
    TH1F *h_npart_double = new TH1F("h_npart_double", "nparticles per event; nparticles; Events", 100, 0, 500);

    // --- Truth photon pT ---
    TH1F *h_pt_single       = new TH1F("h_pt_single",       "Truth photon pT (denom); pT [GeV]; Counts", 60, 0, 60);
    TH1F *h_pt_double       = new TH1F("h_pt_double",       "Truth photon pT (denom); pT [GeV]; Counts", 60, 0, 60);
    TH1F *h_pt_double_extra = new TH1F("h_pt_double_extra", "Extra truth photon pT (2nd+ photon); pT [GeV]; Counts", 60, 0, 60);

    // --- Photonclass distribution ---
    TH1F *h_class_single = new TH1F("h_class_single", "photonclass (pid==22); class; Counts", 6, -0.5, 5.5);
    TH1F *h_class_double = new TH1F("h_class_double", "photonclass (pid==22); class; Counts", 6, -0.5, 5.5);

    // --- All photons (pid==22) pT, regardless of class ---
    TH1F *h_allphoton_pt_single = new TH1F("h_allphoton_pt_single", "All photon pT (pid==22); pT [GeV]; Counts", 60, 0, 60);
    TH1F *h_allphoton_pt_double = new TH1F("h_allphoton_pt_double", "All photon pT (pid==22); pT [GeV]; Counts", 60, 0, 60);

    // --- Truth iso distribution for denom photons ---
    TH1F *h_iso_single = new TH1F("h_iso_single", "Truth iso E_{T} (R=0.3, denom photons); iso E_{T} [GeV]; Counts", 50, 0, 20);
    TH1F *h_iso_double = new TH1F("h_iso_double", "Truth iso E_{T} (R=0.3, denom photons); iso E_{T} [GeV]; Counts", 50, 0, 20);

    // --- Track ID distribution for denom photons ---
    TH1F *h_trkid_single = new TH1F("h_trkid_single", "particle_trkid (denom photons); trkid; Counts", 200, 0, 200);
    TH1F *h_trkid_double = new TH1F("h_trkid_double", "particle_trkid (denom photons); trkid; Counts", 200, 0, 200);

    // --- N denom photons vs max photon pT ---
    TH2F *h_ntruth_vs_maxpt_double = new TH2F("h_ntruth_vs_maxpt_double",
        "N denom photons vs max photon pT (double); max photon pT [GeV]; N denom photons",
        40, 0, 40, 6, -0.5, 5.5);

    // --- vertexz_truth distribution ---
    TH1F *h_vtxz_single = new TH1F("h_vtxz_single", "vertexz_truth; z [cm]; Events", 100, -50, 50);
    TH1F *h_vtxz_double = new TH1F("h_vtxz_double", "vertexz_truth; z [cm]; Events", 100, -50, 50);

    // Colors
    h_ntruth_single->SetLineColor(kBlue);   h_ntruth_single->SetLineWidth(2);
    h_ntruth_double->SetLineColor(kRed);    h_ntruth_double->SetLineWidth(2);
    h_ntruth_double_ptcut->SetLineColor(kRed); h_ntruth_double_ptcut->SetLineStyle(2); h_ntruth_double_ptcut->SetLineWidth(2);

    // =========================================================================
    // Process function — runs on one file
    // =========================================================================
    auto processFile = [&](const char *fname, bool isDouble,
                           TH1F *h_ntruth, TH1F *h_npart, TH1F *h_pt,
                           TH1F *h_class, TH1F *h_allpt, TH1F *h_iso,
                           TH1F *h_trkid, TH1F *h_vtxz,
                           TH1F *h_ntruth_ptcut_arg = nullptr,
                           TH1F *h_pt_extra_arg = nullptr,
                           TH2F *h_ntruth_vs_maxpt_arg = nullptr)
    {
        TFile *f = TFile::Open(fname);
        if (!f || f->IsZombie()) { std::cerr << "Cannot open " << fname << std::endl; return; }
        TTree *tree = (TTree *)f->Get("slimtree");
        if (!tree) { std::cerr << "Cannot find slimtree in " << fname << std::endl; return; }

        int nparticles;
        float vertexz_truth;
        const int MAXPART = 5000;
        float particle_Pt[MAXPART], particle_Eta[MAXPART];
        float particle_truth_iso_03[MAXPART];
        int particle_pid[MAXPART], particle_trkid[MAXPART];
        int particle_photonclass[MAXPART];

        tree->SetBranchAddress("nparticles", &nparticles);
        tree->SetBranchAddress("vertexz_truth", &vertexz_truth);
        tree->SetBranchAddress("particle_Pt", particle_Pt);
        tree->SetBranchAddress("particle_Eta", particle_Eta);
        tree->SetBranchAddress("particle_truth_iso_03", particle_truth_iso_03);
        tree->SetBranchAddress("particle_pid", particle_pid);
        tree->SetBranchAddress("particle_trkid", particle_trkid);
        tree->SetBranchAddress("particle_photonclass", particle_photonclass);

        Long64_t nEntries = tree->GetEntries();
        if (maxEvents > 0 && nEntries > maxEvents) nEntries = maxEvents;

        std::cout << "\n=== Processing " << fname << " ===" << std::endl;
        std::cout << "Entries to process: " << nEntries << std::endl;

        // Counters
        long long total_events = 0;
        long long total_denom_photons = 0;
        long long events_with_0_denom = 0;
        long long events_with_1_denom = 0;
        long long events_with_2plus_denom = 0;
        long long total_pid22 = 0;
        long long total_pid22_class12 = 0;
        long long total_pid22_class12_eta = 0;
        long long total_pid22_class12_eta_iso = 0;

        // For event-level pT window
        float pt_lower = isDouble ? double_pt_lower : single_pt_lower;
        float pt_upper = isDouble ? double_pt_upper : single_pt_upper;
        long long events_pass_ptwindow = 0;
        long long denom_photons_pass_ptwindow = 0;

        for (Long64_t i = 0; i < nEntries; i++)
        {
            tree->GetEntry(i);
            total_events++;
            h_npart->Fill(nparticles);
            h_vtxz->Fill(vertexz_truth);

            // Find max photon pT (for event-level cut)
            float maxphotonpT = 0;
            for (int ip = 0; ip < nparticles && ip < MAXPART; ip++)
            {
                if (particle_pid[ip] == 22)
                {
                    if (particle_Pt[ip] > maxphotonpT)
                        maxphotonpT = particle_Pt[ip];
                }
            }

            bool pass_ptwindow = (maxphotonpT >= pt_lower && maxphotonpT <= pt_upper);

            // Count denominator photons
            int n_denom = 0;
            int n_denom_ptcut = 0;
            bool first_denom = true;

            for (int ip = 0; ip < nparticles && ip < MAXPART; ip++)
            {
                // All photons
                if (particle_pid[ip] == 22)
                {
                    total_pid22++;
                    h_allpt->Fill(particle_Pt[ip]);
                    h_class->Fill(particle_photonclass[ip]);

                    // Direct or fragmentation
                    if (particle_photonclass[ip] >= 1 && particle_photonclass[ip] <= 2)
                    {
                        total_pid22_class12++;

                        // Eta cut
                        if (particle_Eta[ip] > eta_min && particle_Eta[ip] < eta_max)
                        {
                            total_pid22_class12_eta++;

                            // Isolation cut
                            if (particle_truth_iso_03[ip] < truthisocut)
                            {
                                total_pid22_class12_eta_iso++;
                                n_denom++;

                                h_pt->Fill(particle_Pt[ip]);
                                h_iso->Fill(particle_truth_iso_03[ip]);
                                h_trkid->Fill(particle_trkid[ip]);

                                if (pass_ptwindow) n_denom_ptcut++;

                                // Track extra photons in double
                                if (isDouble && h_pt_extra_arg && !first_denom)
                                {
                                    h_pt_extra_arg->Fill(particle_Pt[ip]);
                                }
                                first_denom = false;
                            }
                        }
                    }
                }
            }

            h_ntruth->Fill(n_denom);
            if (n_denom == 0) events_with_0_denom++;
            else if (n_denom == 1) events_with_1_denom++;
            else events_with_2plus_denom++;

            total_denom_photons += n_denom;

            if (pass_ptwindow)
            {
                events_pass_ptwindow++;
                denom_photons_pass_ptwindow += n_denom_ptcut;
                if (h_ntruth_ptcut_arg) h_ntruth_ptcut_arg->Fill(n_denom_ptcut);
            }

            if (isDouble && h_ntruth_vs_maxpt_arg)
            {
                h_ntruth_vs_maxpt_arg->Fill(maxphotonpT, n_denom);
            }
        }

        // Print summary
        std::cout << "\n--- Summary (" << (isDouble ? "DOUBLE" : "SINGLE") << ") ---" << std::endl;
        std::cout << "Total events processed:        " << total_events << std::endl;
        std::cout << "Total pid==22 particles:       " << total_pid22 << std::endl;
        std::cout << "  class 1 or 2:                " << total_pid22_class12 << std::endl;
        std::cout << "  + |eta| < 0.7:               " << total_pid22_class12_eta << std::endl;
        std::cout << "  + iso < 4 GeV (DENOM):       " << total_pid22_class12_eta_iso << std::endl;
        std::cout << "Total denom photons:           " << total_denom_photons << std::endl;
        std::cout << "Avg denom photons per event:   " << (float)total_denom_photons / total_events << std::endl;
        std::cout << "Events with 0 denom photons:   " << events_with_0_denom
                  << " (" << 100.0 * events_with_0_denom / total_events << "%)" << std::endl;
        std::cout << "Events with 1 denom photon:    " << events_with_1_denom
                  << " (" << 100.0 * events_with_1_denom / total_events << "%)" << std::endl;
        std::cout << "Events with 2+ denom photons:  " << events_with_2plus_denom
                  << " (" << 100.0 * events_with_2plus_denom / total_events << "%)" << std::endl;
        std::cout << "\nWith event-level max photon pT window [" << pt_lower << ", " << pt_upper << "]:" << std::endl;
        std::cout << "Events passing pT window:      " << events_pass_ptwindow
                  << " (" << 100.0 * events_pass_ptwindow / total_events << "%)" << std::endl;
        std::cout << "Denom photons (pT window):     " << denom_photons_pass_ptwindow << std::endl;
        if (events_pass_ptwindow > 0)
            std::cout << "Avg denom per event (pT win):  " << (float)denom_photons_pass_ptwindow / events_pass_ptwindow << std::endl;

        f->Close();
    };

    // =========================================================================
    // Run both samples
    // =========================================================================
    processFile(singleFile, false,
                h_ntruth_single, h_npart_single, h_pt_single,
                h_class_single, h_allphoton_pt_single, h_iso_single,
                h_trkid_single, h_vtxz_single);

    processFile(doubleFile, true,
                h_ntruth_double, h_npart_double, h_pt_double,
                h_class_double, h_allphoton_pt_double, h_iso_double,
                h_trkid_double, h_vtxz_double,
                h_ntruth_double_ptcut, h_pt_double_extra, h_ntruth_vs_maxpt_double);

    // =========================================================================
    // Make comparison plots
    // =========================================================================

    // --- Plot 1: N denom photons per event ---
    TCanvas *c1 = new TCanvas("c1", "N denom photons", 800, 600);
    h_ntruth_double->Draw("hist");
    h_ntruth_single->Draw("hist same");
    TLegend *leg1 = new TLegend(0.55, 0.65, 0.88, 0.85);
    leg1->AddEntry(h_ntruth_single, Form("Single (mean=%.3f)", h_ntruth_single->GetMean()), "l");
    leg1->AddEntry(h_ntruth_double, Form("Double (mean=%.3f)", h_ntruth_double->GetMean()), "l");
    leg1->Draw();
    c1->SaveAs("results/investigate_ndenom.pdf");

    // --- Plot 2: nparticles per event ---
    TCanvas *c2 = new TCanvas("c2", "nparticles", 800, 600);
    h_npart_double->SetLineColor(kRed); h_npart_double->SetLineWidth(2);
    h_npart_single->SetLineColor(kBlue); h_npart_single->SetLineWidth(2);
    h_npart_double->Draw("hist");
    h_npart_single->Draw("hist same");
    TLegend *leg2 = new TLegend(0.55, 0.65, 0.88, 0.85);
    leg2->AddEntry(h_npart_single, Form("Single (mean=%.1f)", h_npart_single->GetMean()), "l");
    leg2->AddEntry(h_npart_double, Form("Double (mean=%.1f)", h_npart_double->GetMean()), "l");
    leg2->Draw();
    c2->SaveAs("results/investigate_nparticles.pdf");

    // --- Plot 3: denom photon pT ---
    TCanvas *c3 = new TCanvas("c3", "Denom photon pT", 800, 600);
    c3->SetLogy();
    h_pt_single->SetLineColor(kBlue); h_pt_single->SetLineWidth(2);
    h_pt_double->SetLineColor(kRed);  h_pt_double->SetLineWidth(2);
    h_pt_double_extra->SetLineColor(kMagenta); h_pt_double_extra->SetLineWidth(2); h_pt_double_extra->SetLineStyle(2);
    h_pt_double->Draw("hist");
    h_pt_single->Draw("hist same");
    h_pt_double_extra->Draw("hist same");
    TLegend *leg3 = new TLegend(0.45, 0.65, 0.88, 0.88);
    leg3->AddEntry(h_pt_single, Form("Single (N=%.0f)", h_pt_single->GetEntries()), "l");
    leg3->AddEntry(h_pt_double, Form("Double (N=%.0f)", h_pt_double->GetEntries()), "l");
    leg3->AddEntry(h_pt_double_extra, Form("Double extra (2nd+ photon, N=%.0f)", h_pt_double_extra->GetEntries()), "l");
    leg3->Draw();
    c3->SaveAs("results/investigate_denom_pt.pdf");

    // --- Plot 4: photonclass ---
    TCanvas *c4 = new TCanvas("c4", "Photonclass", 800, 600);
    c4->SetLogy();
    h_class_single->SetLineColor(kBlue); h_class_single->SetLineWidth(2);
    h_class_double->SetLineColor(kRed);  h_class_double->SetLineWidth(2);
    // Normalize to same number of events for shape comparison
    float norm_s = 1.0 / h_class_single->Integral();
    float norm_d = 1.0 / h_class_double->Integral();
    TH1F *h_class_single_n = (TH1F *)h_class_single->Clone("h_class_single_n");
    TH1F *h_class_double_n = (TH1F *)h_class_double->Clone("h_class_double_n");
    h_class_single_n->Scale(norm_s);
    h_class_double_n->Scale(norm_d);
    h_class_double_n->Draw("hist");
    h_class_single_n->Draw("hist same");
    TLegend *leg4 = new TLegend(0.55, 0.65, 0.88, 0.85);
    leg4->AddEntry(h_class_single_n, "Single (normalized)", "l");
    leg4->AddEntry(h_class_double_n, "Double (normalized)", "l");
    leg4->Draw();
    c4->SaveAs("results/investigate_photonclass.pdf");

    // --- Plot 5: truth iso ---
    TCanvas *c5 = new TCanvas("c5", "Truth iso", 800, 600);
    c5->SetLogy();
    h_iso_single->SetLineColor(kBlue); h_iso_single->SetLineWidth(2);
    h_iso_double->SetLineColor(kRed);  h_iso_double->SetLineWidth(2);
    h_iso_double->Draw("hist");
    h_iso_single->Draw("hist same");
    TLegend *leg5 = new TLegend(0.55, 0.65, 0.88, 0.85);
    leg5->AddEntry(h_iso_single, "Single", "l");
    leg5->AddEntry(h_iso_double, "Double", "l");
    leg5->Draw();
    c5->SaveAs("results/investigate_truthiso.pdf");

    // --- Plot 6: trkid ---
    TCanvas *c6 = new TCanvas("c6", "Track ID", 800, 600);
    c6->SetLogy();
    h_trkid_single->SetLineColor(kBlue); h_trkid_single->SetLineWidth(2);
    h_trkid_double->SetLineColor(kRed);  h_trkid_double->SetLineWidth(2);
    h_trkid_double->Draw("hist");
    h_trkid_single->Draw("hist same");
    TLegend *leg6 = new TLegend(0.55, 0.65, 0.88, 0.85);
    leg6->AddEntry(h_trkid_single, "Single", "l");
    leg6->AddEntry(h_trkid_double, "Double", "l");
    leg6->Draw();
    c6->SaveAs("results/investigate_trkid.pdf");

    // --- Plot 7: vertexz_truth ---
    TCanvas *c7 = new TCanvas("c7", "vertexz_truth", 800, 600);
    h_vtxz_single->SetLineColor(kBlue); h_vtxz_single->SetLineWidth(2);
    h_vtxz_double->SetLineColor(kRed);  h_vtxz_double->SetLineWidth(2);
    h_vtxz_double->Draw("hist");
    h_vtxz_single->Draw("hist same");
    TLegend *leg7 = new TLegend(0.55, 0.65, 0.88, 0.85);
    leg7->AddEntry(h_vtxz_single, "Single", "l");
    leg7->AddEntry(h_vtxz_double, "Double", "l");
    leg7->Draw();
    c7->SaveAs("results/investigate_vtxz.pdf");

    // --- Plot 8: all photon pT ---
    TCanvas *c8 = new TCanvas("c8", "All photon pT", 800, 600);
    c8->SetLogy();
    h_allphoton_pt_single->SetLineColor(kBlue); h_allphoton_pt_single->SetLineWidth(2);
    h_allphoton_pt_double->SetLineColor(kRed);  h_allphoton_pt_double->SetLineWidth(2);
    h_allphoton_pt_double->Draw("hist");
    h_allphoton_pt_single->Draw("hist same");
    TLegend *leg8 = new TLegend(0.45, 0.65, 0.88, 0.85);
    leg8->AddEntry(h_allphoton_pt_single, Form("Single (all pid==22, N=%.0f)", h_allphoton_pt_single->GetEntries()), "l");
    leg8->AddEntry(h_allphoton_pt_double, Form("Double (all pid==22, N=%.0f)", h_allphoton_pt_double->GetEntries()), "l");
    leg8->Draw();
    c8->SaveAs("results/investigate_allphoton_pt.pdf");

    // Write and close
    fout->Write();
    fout->Close();

    std::cout << "\n=== Output saved to results/investigate_truth_photons.root ===" << std::endl;
    std::cout << "=== PDF plots saved to results/investigate_*.pdf ===" << std::endl;
}
