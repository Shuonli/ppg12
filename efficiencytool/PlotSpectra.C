#include <iostream>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TNamed.h>
#include <TMath.h>

using namespace std;

void PlotSpectra(string filetype,
                 double min_photon_pt = -1,
                 double max_photon_pt = -1,
                 double min_jet_pt = -1,
                 double max_jet_pt = -1)
{
    cout << "\n========== PlotSpectra Macro ==========" << endl;
    cout << "Sample type: " << filetype << endl;
    cout << "Photon pT filter: ";
    if (min_photon_pt > 0 || max_photon_pt > 0) {
        cout << (min_photon_pt > 0 ? min_photon_pt : 0) << " - "
             << (max_photon_pt > 0 ? max_photon_pt : 999) << " GeV" << endl;
    } else {
        cout << "None" << endl;
    }
    cout << "Jet pT filter: ";
    if (min_jet_pt > 0 || max_jet_pt > 0) {
        cout << (min_jet_pt > 0 ? min_jet_pt : 0) << " - "
             << (max_jet_pt > 0 ? max_jet_pt : 999) << " GeV" << endl;
    } else {
        cout << "None" << endl;
    }
    cout << "======================================\n" << endl;

    // ===== Cross-section constants =====
    // Photon samples (pb)
    const float photon5cross = 146359.3;
    const float photon10cross = 6944.675;
    const float photon20cross = 130.4461;

    // Hanpu uses unit in b
    const float jet10cross = 3.997e+06;
    const float jet15cross = 4.073e+05;
    const float jet20cross = 6.218e+04;
    const float jet30cross = 2.502e+03;
    const float jet50cross = 7.2695;

    // ===== Calculate cross-section weight =====
    float weight = 1.0;

    if (filetype == "photon5") {
        weight = photon5cross / photon20cross;
    } else if (filetype == "photon10") {
        weight = photon10cross / photon20cross;
    } else if (filetype == "photon20") {
        weight = 1.0;
    } else if (filetype == "jet10") {
        weight = jet10cross / jet50cross;
    } else if (filetype == "jet15") {
        weight = jet15cross / jet50cross;
    } else if (filetype == "jet20") {
        weight = jet20cross / jet50cross;
    } else if (filetype == "jet30") {
        weight = jet30cross / jet50cross;
    } else if (filetype == "jet50") {
        weight = 1.0;
    } else {
        cout << "ERROR: Unknown filetype: " << filetype << endl;
        cout << "Valid types: photon5, photon10, photon20, jet10, jet15, jet20, jet30, jet50" << endl;
        return;
    }

    cout << "Cross-section weight: " << weight << endl;

    // ===== Build input file path =====
    string base_dir = "/sphenix/user/shuhangli/ppg12/FunWithxgboost/";
    string branch_dir = "/bdt_1214.root";
    string input_file = base_dir + filetype + branch_dir;

    cout << "Input file: " << input_file << endl;

    // ===== Create TChain and TTreeReader =====
    TChain chain("slimtree");
    chain.Add(input_file.c_str());

    if (chain.GetEntries() == 0) {
        cout << "ERROR: No entries found in tree!" << endl;
        cout << "Check if file exists: " << input_file << endl;
        return;
    }

    cout << "Total entries in tree: " << chain.GetEntries() << endl;

    TTreeReader reader(&chain);

    // ===== Set up TTreeReader branches =====
    // Event-level variables
    TTreeReaderValue<float> vertexz_truth(reader, "vertexz_truth");
    TTreeReaderValue<int> nparticles(reader, "nparticles");
    TTreeReaderValue<int> njet_truth(reader, "njet_truth");
    TTreeReaderValue<int> ncluster(reader, "ncluster_CLUSTERINFO_CEMC_NO_SPLIT");

    // Truth particle arrays
    TTreeReaderArray<float> particle_Pt(reader, "particle_Pt");
    TTreeReaderArray<float> particle_Eta(reader, "particle_Eta");
    TTreeReaderArray<int> particle_pid(reader, "particle_pid");

    // Truth jet arrays
    TTreeReaderArray<float> jet_truth_Pt(reader, "jet_truth_Pt");

    // Reco cluster arrays
    TTreeReaderArray<float> cluster_Et(reader, "cluster_Et_CLUSTERINFO_CEMC_NO_SPLIT");

    // ===== Create histograms =====
    // Truth photon spectra
    TH1F *h_truth_photon_pt = new TH1F("h_truth_photon_pt",
        "Truth Photon p_{T};p_{T} [GeV];Weighted Entries",
        100, 0, 100);
    h_truth_photon_pt->Sumw2();

    TH1F *h_max_truth_photon_pt = new TH1F("h_max_truth_photon_pt",
        "Max Truth Photon p_{T} per Event;Max p_{T} [GeV];Weighted Events",
        50, 0, 100);
    h_max_truth_photon_pt->Sumw2();

    // Truth jet spectra
    TH1F *h_truth_jet_pt = new TH1F("h_truth_jet_pt",
        "Truth Jet p_{T};p_{T} [GeV];Weighted Entries",
        100, 0, 100);
    h_truth_jet_pt->Sumw2();

    TH1F *h_max_truth_jet_pt = new TH1F("h_max_truth_jet_pt",
        "Max Truth Jet p_{T} per Event;Max p_{T} [GeV];Weighted Events",
        50, 0, 100);
    h_max_truth_jet_pt->Sumw2();

    // Reco cluster spectrum
    TH1F *h_reco_cluster_et = new TH1F("h_reco_cluster_et",
        "Reco Cluster E_{T};E_{T} [GeV];Weighted Entries",
        100, 0, 100);
    h_reco_cluster_et->Sumw2();

    // QA histograms
    TH1F *h_vertex_z = new TH1F("h_vertex_z",
        "Vertex z;z [cm];Weighted Events",
        100, -50, 50);
    h_vertex_z->Sumw2();

    TH1F *h_nparticles = new TH1F("h_nparticles",
        "Number of Truth Particles;N_{particles};Weighted Events",
        100, 0, 100);
    h_nparticles->Sumw2();

    TH1F *h_nclusters = new TH1F("h_nclusters",
        "Number of Reco Clusters;N_{clusters};Weighted Events",
        100, 0, 100);
    h_nclusters->Sumw2();

    // Event counting histogram
    TH1F *h_event_counts = new TH1F("h_event_counts",
        "Event Counts;Category;Events",
        5, 0, 5);
    h_event_counts->GetXaxis()->SetBinLabel(1, "Total");
    h_event_counts->GetXaxis()->SetBinLabel(2, "Pass Vertex");
    h_event_counts->GetXaxis()->SetBinLabel(3, "Pass Photon Filter");
    h_event_counts->GetXaxis()->SetBinLabel(4, "Pass Jet Filter");
    h_event_counts->GetXaxis()->SetBinLabel(5, "Pass All");

    // ===== Event loop counters =====
    int total_events = 0;
    int pass_vertex = 0;
    int pass_photon_filter = 0;
    int pass_jet_filter = 0;
    int pass_all = 0;

    const float vertex_cut = 30.0;  // cm
    const float eta_cut = 0.7;      // for photon counting

    cout << "\nProcessing events..." << endl;

    // ===== Event loop =====
    while (reader.Next()) {
        total_events++;

        // Progress indicator
        if (total_events % 10000 == 0) {
            cout << "Processing event " << total_events << "..." << endl;
        }

        // 1. Vertex cut
        if (abs(*vertexz_truth) > vertex_cut) continue;
        pass_vertex++;

        // 2. Find max photon pT in event (only photons with |eta| < 0.7)
        float max_photon_pt_event = 0.0;
        for (int ip = 0; ip < *nparticles; ip++) {
            if (particle_pid[ip] == 22 && abs(particle_Eta[ip]) < eta_cut) {
                if (particle_Pt[ip] > max_photon_pt_event) {
                    max_photon_pt_event = particle_Pt[ip];
                }
            }
        }

        // 3. Find max jet pT in event
        float max_jet_pt_event = 0.0;
        for (int ij = 0; ij < *njet_truth; ij++) {
            if (jet_truth_Pt[ij] > max_jet_pt_event) {
                max_jet_pt_event = jet_truth_Pt[ij];
            }
        }

        // 4. Apply photon pT filter
        if (min_photon_pt > 0 && max_photon_pt_event < min_photon_pt) continue;
        if (max_photon_pt > 0 && max_photon_pt_event > max_photon_pt) continue;
        pass_photon_filter++;

        // 5. Apply jet pT filter
        if (min_jet_pt > 0 && max_jet_pt_event < min_jet_pt) continue;
        if (max_jet_pt > 0 && max_jet_pt_event > max_jet_pt) continue;
        pass_jet_filter++;
        pass_all++;

        // ===== Fill event-level histograms =====
        h_max_truth_photon_pt->Fill(max_photon_pt_event, weight);
        h_max_truth_jet_pt->Fill(max_jet_pt_event, weight);
        h_vertex_z->Fill(*vertexz_truth, weight);
        h_nparticles->Fill(*nparticles, weight);
        h_nclusters->Fill(*ncluster, weight);

        // ===== Fill object-level histograms (all objects in passing events) =====
        // Fill all truth photons
        for (int ip = 0; ip < *nparticles; ip++) {
            if (particle_pid[ip] == 22) {
                h_truth_photon_pt->Fill(particle_Pt[ip], weight);
            }
        }

        // Fill all truth jets
        for (int ij = 0; ij < *njet_truth; ij++) {
            h_truth_jet_pt->Fill(jet_truth_Pt[ij], weight);
        }

        // Fill all reco clusters
        for (int ic = 0; ic < *ncluster; ic++) {
            h_reco_cluster_et->Fill(cluster_Et[ic], weight);
        }
    }

    // Fill event counting histogram
    h_event_counts->SetBinContent(1, total_events);
    h_event_counts->SetBinContent(2, pass_vertex);
    h_event_counts->SetBinContent(3, pass_photon_filter);
    h_event_counts->SetBinContent(4, pass_jet_filter);
    h_event_counts->SetBinContent(5, pass_all);

    // ===== Generate output filename =====
    string output_name = "results/PlotSpectra_" + filetype;

    if (min_photon_pt > 0 || max_photon_pt > 0) {
        char photon_suffix[100];
        sprintf(photon_suffix, "_photon%.0fto%s",
                min_photon_pt > 0 ? min_photon_pt : 0,
                max_photon_pt > 0 ? Form("%.0f", max_photon_pt) : "inf");
        output_name += photon_suffix;
    }

    if (min_jet_pt > 0 || max_jet_pt > 0) {
        char jet_suffix[100];
        sprintf(jet_suffix, "_jet%.0fto%s",
                min_jet_pt > 0 ? min_jet_pt : 0,
                max_jet_pt > 0 ? Form("%.0f", max_jet_pt) : "inf");
        output_name += jet_suffix;
    }

    output_name += ".root";

    // ===== Save output =====
    TFile *fout = new TFile(output_name.c_str(), "RECREATE");

    // Write all histograms
    h_truth_photon_pt->Write();
    h_max_truth_photon_pt->Write();
    h_truth_jet_pt->Write();
    h_max_truth_jet_pt->Write();
    h_reco_cluster_et->Write();
    h_vertex_z->Write();
    h_nparticles->Write();
    h_nclusters->Write();
    h_event_counts->Write();

    // Save metadata as TNamed objects
    TNamed("filetype", filetype.c_str()).Write();
    TNamed("cross_section_weight", Form("%.6e", weight)).Write();
    TNamed("min_photon_pt", Form("%.2f", min_photon_pt)).Write();
    TNamed("max_photon_pt", Form("%.2f", max_photon_pt)).Write();
    TNamed("min_jet_pt", Form("%.2f", min_jet_pt)).Write();
    TNamed("max_jet_pt", Form("%.2f", max_jet_pt)).Write();

    fout->Close();

    // ===== Print summary =====
    cout << "\n========== Event Summary ==========" << endl;
    cout << "Sample: " << filetype << endl;
    cout << "Cross-section weight: " << weight << endl;
    cout << "Total events: " << total_events << endl;
    cout << "Pass vertex cut: " << pass_vertex << " ("
         << 100.0*pass_vertex/total_events << "%)" << endl;
    cout << "Pass photon filter: " << pass_photon_filter << " ("
         << 100.0*pass_photon_filter/total_events << "%)" << endl;
    cout << "Pass jet filter: " << pass_jet_filter << " ("
         << 100.0*pass_jet_filter/total_events << "%)" << endl;
    cout << "Pass all filters: " << pass_all << " ("
         << 100.0*pass_all/total_events << "%)" << endl;
    cout << "Output file: " << output_name << endl;
    cout << "==================================\n" << endl;

    cout << "Done! Histograms saved to: " << output_name << endl;
}
