#include <iostream>
#include <string>
#include <map>
#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TSystem.h>
#include <yaml-cpp/yaml.h>
#include <cmath>

const float TIME_SAMPLE_NS = 17.6;

// Cross-section constants (in pb for photons, b for jets)
const float photon5cross = 146359.3;
const float photon10cross = 6944.675;
const float photon20cross = 130.4461;
const float jet10cross = 3.997e+06;
const float jet15cross = 4.073e+05;
const float jet20cross = 6.218e+04;
const float jet30cross = 2.502e+03;
const float jet50cross = 7.2695;

namespace
{
    inline float deltaPhi(float a, float b)
    {
        float d = a - b;
        while (d > (float)M_PI)
            d -= (float)(2.0 * M_PI);
        while (d <= (float)-M_PI)
            d += (float)(2.0 * M_PI);
        return d;
    }

    inline float deltaR(float eta1, float phi1, float eta2, float phi2)
    {
        const float deta = eta1 - eta2;
        const float dphi = deltaPhi(phi1, phi2);
        return std::sqrt(deta * deta + dphi * dphi);
    }
} // namespace

void AnalyzeTruthPhotonTowers(
    const std::string &configname = "config_bdt_none.yaml",
    const std::string &filetype = "photon20",
    const std::string &outfilename = "results/truth_photon_tower_analysis.root")
{
    // Load yaml-cpp library
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");

    std::cout << "==================================================" << std::endl;
    std::cout << "  Tower Analysis for Truth-Matched Photons" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "Config file: " << configname << std::endl;
    std::cout << "File type: " << filetype << std::endl;
    std::cout << "Output file: " << outfilename << std::endl;
    std::cout << std::endl;

    // Load YAML configuration
    YAML::Node configYaml = YAML::LoadFile(configname);

    // Get input file path from YAML
    std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
    std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;

    // Get cluster node name from YAML
    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();

    // Get tree name from YAML
    std::string treename = configYaml["input"]["tree"].as<std::string>();

    std::cout << "Input file: " << infilename << std::endl;
    std::cout << "Tree name: " << treename << std::endl;
    std::cout << "Cluster node: " << clusternodename << std::endl;
    std::cout << std::endl;

    // Determine cross-section weight based on file type
    float weight = 1.0;
    // Truth-level leading-photon pT range cuts (to avoid double counting between samples),
    // matching the convention used in RecoEffCalculator_TTreeReader.C
    float max_photon_lower = 0.0;
    float max_photon_upper = 1.0e9;
    if (filetype == "photon5")
    {
        max_photon_lower = 0.0;
        max_photon_upper = 10.0;
        weight = photon5cross / photon20cross;
    }
    else if (filetype == "photon10")
    {
        max_photon_lower = 10.0;
        max_photon_upper = 20.0;
        weight = photon10cross / photon20cross;
    }
    else if (filetype == "photon20")
    {
        max_photon_lower = 20.0;
        max_photon_upper = 30.0;
        weight = 1.0;
    }
    else if (filetype == "jet10")
    {
        weight = jet10cross / jet50cross;
    }
    else if (filetype == "jet15")
    {
        weight = jet15cross / jet50cross;
    }
    else if (filetype == "jet20")
    {
        weight = jet20cross / jet50cross;
    }
    else if (filetype == "jet30")
    {
        weight = jet30cross / jet50cross;
    }
    else if (filetype == "jet50")
    {
        weight = 1.0;
    }
    else
    {
        std::cout << "WARNING: Unknown filetype '" << filetype << "', using weight = 1.0" << std::endl;
    }

    std::cout << "Cross-section weight: " << weight << std::endl;
    if (filetype.find("photon") != std::string::npos)
    {
        std::cout << "Leading-truth-photon pT window: [" << max_photon_lower << ", " << max_photon_upper << "] GeV" << std::endl;
    }
    std::cout << std::endl;

    // Setup TChain and TTreeReader
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());

    if (chain.GetEntries() == 0)
    {
        std::cerr << "ERROR: No entries found in tree!" << std::endl;
        return;
    }

    std::cout << "Total entries in tree: " << chain.GetEntries() << std::endl;
    std::cout << std::endl;

    TTreeReader reader(&chain);

    // Event-level branches
    TTreeReaderValue<int> nparticles(reader, "nparticles");
    TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));

    // Truth particle branches
    TTreeReaderArray<int> particle_pid(reader, "particle_pid");
    TTreeReaderArray<int> particle_trkid(reader, "particle_trkid");
    TTreeReaderArray<int> particle_photonclass(reader, "particle_photonclass");
    TTreeReaderArray<int> particle_converted(reader, "particle_converted");
    TTreeReaderArray<float> particle_Pt(reader, "particle_Pt");
    TTreeReaderArray<float> particle_Eta(reader, "particle_Eta");
    TTreeReaderArray<float> particle_Phi(reader, "particle_Phi");

    // Cluster branches
    TTreeReaderArray<int> cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusternodename.c_str()));

    // Tower array branches (flattened 2D arrays, 49 towers per cluster)
    TTreeReaderArray<float> cluster_e_array(reader, Form("cluster_e_array_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_adc_array(reader, Form("cluster_adc_array_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_time_array(reader, Form("cluster_time_array_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_status_array(reader, Form("cluster_status_array_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ownership_array(reader, Form("cluster_ownership_array_%s", clusternodename.c_str()));

    // Create output histograms

    // Common eta axis binning (used for cluster-eta-tagged tower/cluster timing plots)
    const int n_eta_bins = 60;
    const double eta_min = -1.2;
    const double eta_max = 1.2;

    TH2D *h_tower_adc_vs_time = new TH2D(
        "h_tower_adc_vs_time",
        "Truth-Matched Photon Towers: ADC vs Time;Tower Time (ns);Tower ADC",
        200, -20, 20,  // Time: -20 to 20 ns
        200, 0, 20000   // ADC: 0 to 20000
    );

    TH2D *h_tower_e_vs_time = new TH2D(
        "h_tower_e_vs_time",
        "Truth-Matched Photon Towers: Energy vs Time;Tower Time (ns);Tower Energy (GeV)",
        200, -20, 20,  // Time: -20 to 20 ns
        200, 0, 20      // Energy: 0 to 20 GeV
    );

    // Tower-time plots with an additional dimension: cluster eta (used as the eta tag for each tower)
    // Axis convention (x,y,z) = (time, ADC/E, cluster eta)
    TH3D *h3_tower_adc_vs_time_vs_eta = new TH3D(
        "h3_tower_adc_vs_time_vs_eta",
        "Truth-Matched Photon Towers: ADC vs Time vs Cluster #eta;Tower Time (ns);Tower ADC;Cluster #eta",
        200, -20, 20,     // Time [ns]
        200, 0, 20000,    // ADC
        n_eta_bins, eta_min, eta_max
    );

    TH3D *h3_tower_e_vs_time_vs_eta = new TH3D(
        "h3_tower_e_vs_time_vs_eta",
        "Truth-Matched Photon Towers: Energy vs Time vs Cluster #eta;Tower Time (ns);Tower Energy (GeV);Cluster #eta",
        200, -20, 20,     // Time [ns]
        200, 0, 20,       // Energy [GeV]
        n_eta_bins, eta_min, eta_max
    );

    // Leading tower histograms
    TH2D *h_leading_tower_adc_vs_time = new TH2D(
        "h_leading_tower_adc_vs_time",
        "Leading Tower: ADC vs Time;Tower Time (ns);Tower ADC",
        200, -20, 20,  // Time: -20 to 20 ns
        200, 0, 20000   // ADC: 0 to 20000
    );

    TH2D *h_leading_tower_e_vs_time = new TH2D(
        "h_leading_tower_e_vs_time",
        "Leading Tower: Energy vs Time;Tower Time (ns);Tower Energy (GeV)",
        200, -20, 20,  // Time: -20 to 20 ns
        200, 0, 20      // Energy: 0 to 20 GeV
    );

    TH3D *h3_leading_tower_adc_vs_time_vs_eta = new TH3D(
        "h3_leading_tower_adc_vs_time_vs_eta",
        "Leading Tower: ADC vs Time vs Cluster #eta;Tower Time (ns);Tower ADC;Cluster #eta",
        200, -20, 20,     // Time [ns]
        200, 0, 20000,    // ADC
        n_eta_bins, eta_min, eta_max
    );

    TH3D *h3_leading_tower_e_vs_time_vs_eta = new TH3D(
        "h3_leading_tower_e_vs_time_vs_eta",
        "Leading Tower: Energy vs Time vs Cluster #eta;Tower Time (ns);Tower Energy (GeV);Cluster #eta",
        200, -20, 20,     // Time [ns]
        200, 0, 20,       // Energy [GeV]
        n_eta_bins, eta_min, eta_max
    );

    // Saturated leading tower histograms (status bit 7 set)
    TH2D *h_leading_tower_saturated_e_vs_time = new TH2D(
        "h_leading_tower_saturated_e_vs_time",
        "Leading Tower (Saturated): Energy vs Time;Tower Time (ns);Tower Energy (GeV)",
        200, -20, 20,  // Time: -20 to 20 ns
        200, 0, 20     // Energy: 0 to 20 GeV
    );

    TH3D *h3_leading_tower_saturated_e_vs_time_vs_eta = new TH3D(
        "h3_leading_tower_saturated_e_vs_time_vs_eta",
        "Leading Tower (Saturated): Energy vs Time vs Cluster #eta;Tower Time (ns);Tower Energy (GeV);Cluster #eta",
        200, -20, 20,     // Time [ns]
        200, 0, 20,       // Energy [GeV]
        n_eta_bins, eta_min, eta_max
    );

    // Non-saturated leading tower histograms (status bit 7 NOT set)
    TH2D *h_leading_tower_nonsaturated_e_vs_time = new TH2D(
        "h_leading_tower_nonsaturated_e_vs_time",
        "Leading Tower (Non-Saturated): Energy vs Time;Tower Time (ns);Tower Energy (GeV)",
        200, -20, 20,  // Time: -20 to 20 ns
        200, 0, 20     // Energy: 0 to 20 GeV
    );

    TH3D *h3_leading_tower_nonsaturated_e_vs_time_vs_eta = new TH3D(
        "h3_leading_tower_nonsaturated_e_vs_time_vs_eta",
        "Leading Tower (Non-Saturated): Energy vs Time vs Cluster #eta;Tower Time (ns);Tower Energy (GeV);Cluster #eta",
        200, -20, 20,     // Time [ns]
        200, 0, 20,       // Energy [GeV]
        n_eta_bins, eta_min, eta_max
    );

    // Subleading tower histograms (tower with 2nd-largest energy in the cluster)
    TH2D *h_subleading_tower_adc_vs_time = new TH2D(
        "h_subleading_tower_adc_vs_time",
        "Subleading Tower: ADC vs Time;Tower Time (ns);Tower ADC",
        200, -20, 20,  // Time: -20 to 20 ns
        200, 0, 20000  // ADC: 0 to 20000
    );

    TH2D *h_subleading_tower_e_vs_time = new TH2D(
        "h_subleading_tower_e_vs_time",
        "Subleading Tower: Energy vs Time;Tower Time (ns);Tower Energy (GeV)",
        200, -20, 20,  // Time: -20 to 20 ns
        200, 0, 20     // Energy: 0 to 20 GeV
    );

    TH3D *h3_subleading_tower_adc_vs_time_vs_eta = new TH3D(
        "h3_subleading_tower_adc_vs_time_vs_eta",
        "Subleading Tower: ADC vs Time vs Cluster #eta;Tower Time (ns);Tower ADC;Cluster #eta",
        200, -20, 20,   // Time [ns]
        200, 0, 20000,  // ADC
        n_eta_bins, eta_min, eta_max
    );

    TH3D *h3_subleading_tower_e_vs_time_vs_eta = new TH3D(
        "h3_subleading_tower_e_vs_time_vs_eta",
        "Subleading Tower: Energy vs Time vs Cluster #eta;Tower Time (ns);Tower Energy (GeV);Cluster #eta",
        200, -20, 20,  // Time [ns]
        200, 0, 20,    // Energy [GeV]
        n_eta_bins, eta_min, eta_max
    );

    // Non-leading tower histograms (all towers except the leading tower)
    TH2D *h_nonleading_tower_adc_vs_time = new TH2D(
        "h_nonleading_tower_adc_vs_time",
        "Non-Leading Towers: ADC vs Time;Tower Time (ns);Tower ADC",
        200, -20, 20,  // Time: -20 to 20 ns
        200, 0, 20000   // ADC: 0 to 20000
    );

    TH2D *h_nonleading_tower_e_vs_time = new TH2D(
        "h_nonleading_tower_e_vs_time",
        "Non-Leading Towers: Energy vs Time;Tower Time (ns);Tower Energy (GeV)",
        200, -20, 20,  // Time: -20 to 20 ns
        200, 0, 20      // Energy: 0 to 20 GeV
    );

    TH3D *h3_nonleading_tower_adc_vs_time_vs_eta = new TH3D(
        "h3_nonleading_tower_adc_vs_time_vs_eta",
        "Non-Leading Towers: ADC vs Time vs Cluster #eta;Tower Time (ns);Tower ADC;Cluster #eta",
        200, -20, 20,     // Time [ns]
        200, 0, 20000,    // ADC
        n_eta_bins, eta_min, eta_max
    );

    TH3D *h3_nonleading_tower_e_vs_time_vs_eta = new TH3D(
        "h3_nonleading_tower_e_vs_time_vs_eta",
        "Non-Leading Towers: Energy vs Time vs Cluster #eta;Tower Time (ns);Tower Energy (GeV);Cluster #eta",
        200, -20, 20,     // Time [ns]
        200, 0, 20,       // Energy [GeV]
        n_eta_bins, eta_min, eta_max
    );

    // Delta-t vs non-leading tower energy (non-leading = all towers except the leading tower)
    TH2D *h2_delta_t_leading_nonleading_vs_e = new TH2D(
        "h2_delta_t_leading_nonleading_vs_e",
        "Delta-t: Leading - Non-Leading vs Non-Leading Energy;#Deltat = t_{leading} - t_{non-leading} [ns];Non-Leading Tower Energy [GeV]",
        200, -10, 10, // Delta-t [ns]
        200, 0, 20    // Non-leading energy [GeV]
    );

    // 3D histogram: number of towers in the cluster vs cluster energy and cluster eta
    // Axis convention (x,y,z) = (cluster energy, cluster eta, N towers in cluster)
    TH3D *h3_n_towers_in_cluster_vs_cluster_e_vs_eta = new TH3D(
        "h3_n_towers_in_cluster_vs_cluster_e_vs_eta",
        "N towers in cluster vs Cluster Energy vs Cluster #eta;Cluster Energy (GeV);Cluster #eta;N towers in cluster",
        200, 0, 50,                  // Cluster energy [GeV]
        n_eta_bins, eta_min, eta_max,// Cluster eta
        50, 0, 50                    // N towers (0..49 possible; last bin includes overflow up to 50)
    );

    // Cluster-time histograms with an additional dimension: cluster eta
    // Axis convention (x,y,z) = (energy, time, eta)
    // Cluster-time vs cluster-energy vs cluster-eta
    TH3D *h3_cluster_time_ew_all_vs_cluster_e_vs_eta = new TH3D(
        "h3_cluster_time_ew_all_vs_cluster_e_vs_eta",
        "Cluster Time (E-weighted, all towers) vs Cluster Energy vs Cluster #eta;Cluster Energy (GeV);Cluster Time (ns);Cluster #eta",
        200, 0, 50,   // Cluster energy [GeV]
        200, -20, 20, // Cluster time [ns]
        n_eta_bins, eta_min, eta_max
    );

    TH3D *h3_cluster_time_leading_vs_cluster_e_vs_eta = new TH3D(
        "h3_cluster_time_leading_vs_cluster_e_vs_eta",
        "Cluster Time (leading tower) vs Cluster Energy vs Cluster #eta;Cluster Energy (GeV);Leading Tower Time (ns);Cluster #eta",
        200, 0, 50,   // Cluster energy [GeV]
        200, -20, 20, // Time [ns]
        n_eta_bins, eta_min, eta_max
    );

    TH3D *h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta = new TH3D(
        "h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta",
        "Cluster Time (E-weighted, non-leading towers) vs Cluster Energy vs Cluster #eta;Cluster Energy (GeV);Non-Leading E-weighted Time (ns);Cluster #eta",
        200, 0, 50,   // Cluster energy [GeV]
        200, -20, 20, // Cluster time [ns]
        n_eta_bins, eta_min, eta_max
    );

    // Cluster-time vs matched truth photon energy vs cluster-eta
    // Truth photon energy is computed as E = pT*cosh(eta) (massless approximation) to avoid reliance on a particle_E branch.
    TH3D *h3_cluster_time_ew_all_vs_truthphoton_e_vs_eta = new TH3D(
        "h3_cluster_time_ew_all_vs_truthphoton_e_vs_eta",
        "Cluster Time (E-weighted, all towers) vs Matched Truth Photon Energy vs Cluster #eta;Truth Photon Energy (GeV);Cluster Time (ns);Cluster #eta",
        200, 0, 50,   // Truth photon energy [GeV]
        200, -20, 20, // Cluster time [ns]
        n_eta_bins, eta_min, eta_max
    );

    TH3D *h3_cluster_time_leading_vs_truthphoton_e_vs_eta = new TH3D(
        "h3_cluster_time_leading_vs_truthphoton_e_vs_eta",
        "Cluster Time (leading tower) vs Matched Truth Photon Energy vs Cluster #eta;Truth Photon Energy (GeV);Leading Tower Time (ns);Cluster #eta",
        200, 0, 50,   // Truth photon energy [GeV]
        200, -20, 20, // Time [ns]
        n_eta_bins, eta_min, eta_max
    );

    TH3D *h3_cluster_time_ew_nonleading_vs_truthphoton_e_vs_eta = new TH3D(
        "h3_cluster_time_ew_nonleading_vs_truthphoton_e_vs_eta",
        "Cluster Time (E-weighted, non-leading towers) vs Matched Truth Photon Energy vs Cluster #eta;Truth Photon Energy (GeV);Non-Leading E-weighted Time (ns);Cluster #eta",
        200, 0, 50,   // Truth photon energy [GeV]
        200, -20, 20, // Cluster time [ns]
        n_eta_bins, eta_min, eta_max
    );

    // Statistics counters
    int n_events = 0;
    int n_photon_clusters = 0;
    int n_towers_filled = 0;
    double n_towers_weighted = 0.0;

    // Event loop
    std::cout << "Starting event loop..." << std::endl;
    while (reader.Next())
    {
        if (n_events % 10000 == 0)
        {
            std::cout << "Processing event " << n_events << " / " << chain.GetEntries() << std::endl;
        }
        n_events++;

        // Build particle track ID map for this event
        std::map<int, int> particle_trkidmap;
        for (int iparticle = 0; iparticle < *nparticles; iparticle++)
        {
            particle_trkidmap[particle_trkid[iparticle]] = iparticle;
        }

        // Find the leading truth photon in this event (highest pT),
        // using the same truth-photon definition as the cluster selection.
        int leading_truth_photon_idx = -1;
        float leading_truth_photon_pt = -1.0;
        for (int iparticle = 0; iparticle < *nparticles; iparticle++)
        {
            if (particle_pid[iparticle] != 22)
            {
                continue;
            }

            // Keep consistent with the photon class requirement used below (direct/frag)
            if (particle_photonclass[iparticle] != 1 && particle_photonclass[iparticle] != 2)
            {
                continue;
            }

            // Require unconverted photon
            if (particle_converted[iparticle] > 0)
            {
                continue;
            }

            const float pt = particle_Pt[iparticle];
            if (pt > leading_truth_photon_pt)
            {
                leading_truth_photon_pt = pt;
                leading_truth_photon_idx = iparticle;
            }
        }

        // If no eligible leading photon exists, skip this event
        if (leading_truth_photon_idx < 0)
        {
            continue;
        }
        //else
        //{
        //    std::cout << "Leading truth photon: " << leading_truth_photon_idx << " with pT: " << leading_truth_photon_pt << std::endl;
        //}

        // Apply truth-level sample window cut using the leading truth photon pT
        // (only for photon samples; jets are unaffected)
        if (filetype.find("photon") != std::string::npos)
        {
            if (leading_truth_photon_pt > max_photon_upper || leading_truth_photon_pt < max_photon_lower)
            {
                continue;
            }
        }

        // Loop over clusters
        for (int icluster = 0; icluster < *ncluster; icluster++)
        {
            // Check if cluster has valid truth match
            if (particle_trkidmap.find(cluster_truthtrkID[icluster]) == particle_trkidmap.end())
            {
                continue; // No truth match
            }

            int iparticle = particle_trkidmap[cluster_truthtrkID[icluster]];

            // Only keep clusters matched to the event's leading truth photon
            if (iparticle != leading_truth_photon_idx)
            {
                continue;
            }

            // Check if matched particle is a photon (direct or fragmentation, not decay)
            if (particle_pid[iparticle] != 22)
            {
                continue; // Not a photon
            }

            if (particle_photonclass[iparticle] != 1 && particle_photonclass[iparticle] != 2)
            {
                continue; // Decay photon, skip (only accept direct or fragmentation)
            }

            if (particle_converted[iparticle] > 0)
            {
                continue; // Converted photon, skip
            }
            //std::cout << "Particle: " << iparticle << " with pT: " << particle_Pt[iparticle] << " and eta: " << particle_Eta[iparticle] << " and phi: " << particle_Phi[iparticle] << std::endl;
            //std::cout << "Cluster: " << icluster << " with Et: " << cluster_Et[icluster] << " and eta: " << cluster_Eta[icluster] << " and phi: " << cluster_Phi[icluster] << std::endl;
            // Truth-reco matching requirement in (eta,phi)
            // Require dR(truth photon, reco cluster) < 0.1
            const float dR = deltaR(particle_Eta[iparticle], particle_Phi[iparticle],
                                    cluster_Eta[icluster], cluster_Phi[icluster]);
            if (dR > 0.1)
            {
                continue;
            }

            // This cluster is truth-matched to a direct or fragmentation, unconverted photon
            n_photon_clusters++;

            // Matched truth photon energy (massless approx)
            const float truth_photon_e = particle_Pt[iparticle] * (float)std::cosh(particle_Eta[iparticle]);

            // First pass: find leading tower by energy
            int leading_tower_idx = -1;
            float max_tower_e = -1.0;
            float leading_tower_time = 0.0;

            // Also track the subleading (2nd-largest energy) tower
            int subleading_tower_idx = -1;
            float subleading_tower_e = -1.0;
            float subleading_tower_time = 0.0;

            // Accumulators for cluster-level quantities (using the same tower filters as histogramming)
            float cluster_e_sum = 0.0;
            double time_e_sum_all = 0.0;
            double e_sum_all = 0.0;
            int n_towers_in_cluster = 0;

            for (int i = 0; i < 49; i++)
            {
                int tower_index = icluster * 49 + i;

                // Filter 1: Check if tower belongs to cluster
                if (cluster_ownership_array[tower_index] != 1)
                {
                    continue;
                }

                // Filter 2: Skip zero-suppressed towers (bit 5 set)
                int status = cluster_status_array[tower_index];
                if (status & (1 << 5))
                {
                    continue;
                }

                // Find tower with maximum energy
                float tower_e = cluster_e_array[tower_index];
                float tower_time_ns = cluster_time_array[tower_index] * TIME_SAMPLE_NS;

                // Cluster energy and energy-weighted time (all towers)
                cluster_e_sum += tower_e;
                time_e_sum_all += (double)tower_e * (double)tower_time_ns;
                e_sum_all += (double)tower_e;
                n_towers_in_cluster++;

                // Track leading and subleading by energy
                if (tower_e > max_tower_e)
                {
                    // previous leading becomes subleading
                    subleading_tower_e = max_tower_e;
                    subleading_tower_idx = leading_tower_idx;
                    subleading_tower_time = leading_tower_time;

                    max_tower_e = tower_e;
                    leading_tower_idx = i;
                    leading_tower_time = tower_time_ns;
                }
                else if (tower_e > subleading_tower_e)
                {
                    subleading_tower_e = tower_e;
                    subleading_tower_idx = i;
                    subleading_tower_time = tower_time_ns;
                }
            }

            // Towers-per-cluster vs (cluster E, cluster eta): fill full multiplicity distribution
            h3_n_towers_in_cluster_vs_cluster_e_vs_eta->Fill(cluster_e_sum, cluster_Eta[icluster], (double)n_towers_in_cluster, weight);

            // Fill leading tower histograms if found
            if (leading_tower_idx >= 0)
            {
                int leading_index = icluster * 49 + leading_tower_idx;
                float leading_adc = cluster_adc_array[leading_index];
                int leading_status = cluster_status_array[leading_index];

                h_leading_tower_adc_vs_time->Fill(leading_tower_time, leading_adc, weight);
                h_leading_tower_e_vs_time->Fill(leading_tower_time, max_tower_e, weight);
                h3_leading_tower_adc_vs_time_vs_eta->Fill(leading_tower_time, leading_adc, cluster_Eta[icluster], weight);
                h3_leading_tower_e_vs_time_vs_eta->Fill(leading_tower_time, max_tower_e, cluster_Eta[icluster], weight);

                // Fill saturated leading tower histograms (bit 7 set)
                if (leading_status & (1 << 7))
                {
                    h_leading_tower_saturated_e_vs_time->Fill(leading_tower_time, max_tower_e, weight);
                    h3_leading_tower_saturated_e_vs_time_vs_eta->Fill(leading_tower_time, max_tower_e, cluster_Eta[icluster], weight);
                }
                else
                {
                    // Fill non-saturated leading tower histograms (bit 7 NOT set)
                    h_leading_tower_nonsaturated_e_vs_time->Fill(leading_tower_time, max_tower_e, weight);
                    h3_leading_tower_nonsaturated_e_vs_time_vs_eta->Fill(leading_tower_time, max_tower_e, cluster_Eta[icluster], weight);
                }
            }

            // Fill subleading tower histograms if found (requires at least 2 eligible towers)
            if (subleading_tower_idx >= 0)
            {
                int subleading_index = icluster * 49 + subleading_tower_idx;
                float subleading_adc = cluster_adc_array[subleading_index];

                h_subleading_tower_adc_vs_time->Fill(subleading_tower_time, subleading_adc, weight);
                h_subleading_tower_e_vs_time->Fill(subleading_tower_time, subleading_tower_e, weight);
                h3_subleading_tower_adc_vs_time_vs_eta->Fill(subleading_tower_time, subleading_adc, cluster_Eta[icluster], weight);
                h3_subleading_tower_e_vs_time_vs_eta->Fill(subleading_tower_time, subleading_tower_e, cluster_Eta[icluster], weight);
            }

            // Cluster time definitions vs cluster energy
            // - All towers: energy-weighted time over all towers passing the filters
            if (e_sum_all > 0.0)
            {
                float cluster_time_ew_all = (float)(time_e_sum_all / e_sum_all);
                h3_cluster_time_ew_all_vs_cluster_e_vs_eta->Fill(cluster_e_sum, cluster_time_ew_all, cluster_Eta[icluster], weight);
                h3_cluster_time_ew_all_vs_truthphoton_e_vs_eta->Fill(truth_photon_e, cluster_time_ew_all, cluster_Eta[icluster], weight);
            }

            // - Leading tower only
            if (leading_tower_idx >= 0)
            {
                h3_cluster_time_leading_vs_cluster_e_vs_eta->Fill(cluster_e_sum, leading_tower_time, cluster_Eta[icluster], weight);
                h3_cluster_time_leading_vs_truthphoton_e_vs_eta->Fill(truth_photon_e, leading_tower_time, cluster_Eta[icluster], weight);
            }

            // - Non-leading towers: energy-weighted time excluding the leading tower
            if (leading_tower_idx >= 0)
            {
                double time_e_sum_nonlead = 0.0;
                double e_sum_nonlead = 0.0;
                for (int i = 0; i < 49; i++)
                {
                    if (i == leading_tower_idx)
                    {
                        continue;
                    }

                    int tower_index = icluster * 49 + i;

                    // Same filters
                    if (cluster_ownership_array[tower_index] != 1)
                    {
                        continue;
                    }
                    int status = cluster_status_array[tower_index];
                    if (status & (1 << 5))
                    {
                        continue;
                    }

                    float tower_e = cluster_e_array[tower_index];
                    float tower_time_ns = cluster_time_array[tower_index] * TIME_SAMPLE_NS;
                    time_e_sum_nonlead += (double)tower_e * (double)tower_time_ns;
                    e_sum_nonlead += (double)tower_e;
                }

                if (e_sum_nonlead > 0.0)
                {
                    float cluster_time_ew_nonlead = (float)(time_e_sum_nonlead / e_sum_nonlead);
                    h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta->Fill(cluster_e_sum, cluster_time_ew_nonlead, cluster_Eta[icluster], weight);
                    h3_cluster_time_ew_nonleading_vs_truthphoton_e_vs_eta->Fill(truth_photon_e, cluster_time_ew_nonlead, cluster_Eta[icluster], weight);
                }
            }

            // Second pass: fill all towers
            for (int i = 0; i < 49; i++)
            {
                int tower_index = icluster * 49 + i;

                // Filter 1: Check if tower belongs to cluster
                if (cluster_ownership_array[tower_index] != 1)
                {
                    continue;
                }

                // Filter 2: Skip zero-suppressed towers (bit 5 set)
                int status = cluster_status_array[tower_index];
                if (status & (1 << 5))
                {
                    continue;
                }

                // Extract tower quantities
                float tower_e = cluster_e_array[tower_index];
                float tower_adc = cluster_adc_array[tower_index];
                float tower_time_ns = cluster_time_array[tower_index] * TIME_SAMPLE_NS;

                // Fill all-tower histograms with cross-section weight
                h_tower_adc_vs_time->Fill(tower_time_ns, tower_adc, weight);
                h_tower_e_vs_time->Fill(tower_time_ns, tower_e, weight);
                h3_tower_adc_vs_time_vs_eta->Fill(tower_time_ns, tower_adc, cluster_Eta[icluster], weight);
                h3_tower_e_vs_time_vs_eta->Fill(tower_time_ns, tower_e, cluster_Eta[icluster], weight);

                // Fill non-leading tower histograms (exclude leading tower)
                if (leading_tower_idx >= 0 && i != leading_tower_idx)
                {
                    h_nonleading_tower_adc_vs_time->Fill(tower_time_ns, tower_adc, weight);
                    h_nonleading_tower_e_vs_time->Fill(tower_time_ns, tower_e, weight);
                    h3_nonleading_tower_adc_vs_time_vs_eta->Fill(tower_time_ns, tower_adc, cluster_Eta[icluster], weight);
                    h3_nonleading_tower_e_vs_time_vs_eta->Fill(tower_time_ns, tower_e, cluster_Eta[icluster], weight);
                }

                // Fill delta-t vs non-leading tower energy (all towers except leading)
                if (leading_tower_idx >= 0 && i != leading_tower_idx)
                {
                    float delta_t = leading_tower_time - tower_time_ns;
                    h2_delta_t_leading_nonleading_vs_e->Fill(delta_t, tower_e, weight);
                }

                n_towers_filled++;
                n_towers_weighted += weight;
            }
        }
    }

    std::cout << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "  Analysis Complete" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "Events processed: " << n_events << std::endl;
    std::cout << "Truth-matched photon clusters (direct/frag, unconverted): " << n_photon_clusters << std::endl;
    std::cout << "Towers filled (unweighted): " << n_towers_filled << std::endl;
    std::cout << "Towers filled (weighted): " << n_towers_weighted << std::endl;
    if (n_photon_clusters > 0)
    {
        std::cout << "Average towers per photon cluster: "
                  << (float)n_towers_filled / n_photon_clusters << std::endl;
    }
    std::cout << std::endl;

    // Write output file
    std::cout << "Writing output to: " << outfilename << std::endl;
    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    h_tower_adc_vs_time->Write();
    h_tower_e_vs_time->Write();
    h_leading_tower_adc_vs_time->Write();
    h_leading_tower_e_vs_time->Write();
    h_subleading_tower_adc_vs_time->Write();
    h_subleading_tower_e_vs_time->Write();
    h_nonleading_tower_adc_vs_time->Write();
    h_nonleading_tower_e_vs_time->Write();
    h3_tower_adc_vs_time_vs_eta->Write();
    h3_tower_e_vs_time_vs_eta->Write();
    h3_leading_tower_adc_vs_time_vs_eta->Write();
    h3_leading_tower_e_vs_time_vs_eta->Write();
    h_leading_tower_saturated_e_vs_time->Write();
    h3_leading_tower_saturated_e_vs_time_vs_eta->Write();
    h_leading_tower_nonsaturated_e_vs_time->Write();
    h3_leading_tower_nonsaturated_e_vs_time_vs_eta->Write();
    h3_subleading_tower_adc_vs_time_vs_eta->Write();
    h3_subleading_tower_e_vs_time_vs_eta->Write();
    h3_nonleading_tower_adc_vs_time_vs_eta->Write();
    h3_nonleading_tower_e_vs_time_vs_eta->Write();
    h2_delta_t_leading_nonleading_vs_e->Write();
    h3_n_towers_in_cluster_vs_cluster_e_vs_eta->Write();
    h3_cluster_time_ew_all_vs_cluster_e_vs_eta->Write();
    h3_cluster_time_leading_vs_cluster_e_vs_eta->Write();
    h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta->Write();
    h3_cluster_time_ew_all_vs_truthphoton_e_vs_eta->Write();
    h3_cluster_time_leading_vs_truthphoton_e_vs_eta->Write();
    h3_cluster_time_ew_nonleading_vs_truthphoton_e_vs_eta->Write();
    fout->Close();

    std::cout << "Done!" << std::endl;
    std::cout << "==================================================" << std::endl;

    delete h_tower_adc_vs_time;
    delete h_tower_e_vs_time;
    delete h_leading_tower_adc_vs_time;
    delete h_leading_tower_e_vs_time;
    delete h_subleading_tower_adc_vs_time;
    delete h_subleading_tower_e_vs_time;
    delete h_nonleading_tower_adc_vs_time;
    delete h_nonleading_tower_e_vs_time;
    delete h3_tower_adc_vs_time_vs_eta;
    delete h3_tower_e_vs_time_vs_eta;
    delete h3_leading_tower_adc_vs_time_vs_eta;
    delete h3_leading_tower_e_vs_time_vs_eta;
    delete h_leading_tower_saturated_e_vs_time;
    delete h3_leading_tower_saturated_e_vs_time_vs_eta;
    delete h_leading_tower_nonsaturated_e_vs_time;
    delete h3_leading_tower_nonsaturated_e_vs_time_vs_eta;
    delete h3_subleading_tower_adc_vs_time_vs_eta;
    delete h3_subleading_tower_e_vs_time_vs_eta;
    delete h3_nonleading_tower_adc_vs_time_vs_eta;
    delete h3_nonleading_tower_e_vs_time_vs_eta;
    delete h2_delta_t_leading_nonleading_vs_e;
    delete h3_n_towers_in_cluster_vs_cluster_e_vs_eta;
    delete h3_cluster_time_ew_all_vs_cluster_e_vs_eta;
    delete h3_cluster_time_leading_vs_cluster_e_vs_eta;
    delete h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta;
    delete h3_cluster_time_ew_all_vs_truthphoton_e_vs_eta;
    delete h3_cluster_time_leading_vs_truthphoton_e_vs_eta;
    delete h3_cluster_time_ew_nonleading_vs_truthphoton_e_vs_eta;
    delete fout;
}
