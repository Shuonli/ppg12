// gen_bins.C
// Pass 2: Separate pThat-bin generation (ROOT macro)
// Bins: [5-10, 10-15, 15-20, 20-30, 30-50] GeV
// Processes: HardQCD:all + PromptPhoton:all, no bias
//
// Usage: root -l -b -q 'gen_bins.C+(10000)'
//        root -l -b -q gen_bins.C+    (default 10000 events per bin)

#ifdef __CLING__
R__LOAD_LIBRARY(libpythia8)
R__LOAD_LIBRARY(libfastjet)
#endif

#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <iostream>

// Walk the mother chain; return true if any ancestor is a hadron
bool hasHadronAncestor_bins(const Pythia8::Event& event, int idx)
{
    int mother = idx;
    while (mother > 0) {
        int nextMother = event[mother].mother1();
        if (nextMother <= 0) break;
        if (event[nextMother].isHadron()) return true;
        mother = nextMother;
    }
    return false;
}

void gen_bins(int nEvents = 10000)
{
    // ---------------------------------------------------------------
    // Configuration
    // ---------------------------------------------------------------
    const double sqrtS        = 200.;
    const double jetR         = 0.4;
    const double jetPtMin     = 5.;
    const double jetEtaMax    = 1.1;
    const double photonPtMin  = 5.;
    const double photonEtaMax = 1.1;

    const std::vector<double> binEdges = {5., 10., 15., 20., 30., 50.};
    const int nBins = static_cast<int>(binEdges.size()) - 1;  // 5

    const char* pythia8share = gSystem->Getenv("PYTHIA8");
    std::string xmldoc = pythia8share ? std::string(pythia8share) + "/xmldoc" : "";

    // ---------------------------------------------------------------
    // Output ROOT file
    // ---------------------------------------------------------------
    TFile* outFile = TFile::Open("output_bins.root", "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "ERROR: cannot open output_bins.root" << std::endl;
        return;
    }

    // ---------------------------------------------------------------
    // Meta tree
    // ---------------------------------------------------------------
    TTree* metaTree = new TTree("meta", "Per-bin metadata");
    Double_t m_bin_low, m_bin_high, m_sigma_gen, m_weight_sum, m_lumi_eff;
    Long64_t m_n_accepted;

    metaTree->Branch("bin_low",    &m_bin_low,    "bin_low/D");
    metaTree->Branch("bin_high",   &m_bin_high,   "bin_high/D");
    metaTree->Branch("sigma_gen",  &m_sigma_gen,  "sigma_gen/D");
    metaTree->Branch("weight_sum", &m_weight_sum, "weight_sum/D");
    metaTree->Branch("n_accepted", &m_n_accepted, "n_accepted/L");
    metaTree->Branch("lumi_eff",   &m_lumi_eff,   "lumi_eff/D");

    // ---------------------------------------------------------------
    // FastJet
    // ---------------------------------------------------------------
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);

    // ---------------------------------------------------------------
    // Loop over pThat bins
    // ---------------------------------------------------------------
    for (int iBin = 0; iBin < nBins; ++iBin) {
        double ptLow  = binEdges[iBin];
        double ptHigh = binEdges[iBin + 1];

        // Per-bin event tree
        std::string treeName  = "events_bin" + std::to_string(iBin);
        std::string treeTitle = Form("pThat [%.0f-%.0f] GeV", ptLow, ptHigh);
        TTree* evTree = new TTree(treeName.c_str(), treeTitle.c_str());

        Double_t t_event_weight, t_pthat;
        std::vector<float> t_jet_pt, t_jet_eta, t_jet_phi;
        std::vector<float> t_photon_pt, t_photon_eta, t_photon_phi;
        std::vector<int>   t_photon_is_prompt;

        evTree->Branch("event_weight",     &t_event_weight,    "event_weight/D");
        evTree->Branch("pthat",            &t_pthat,           "pthat/D");
        evTree->Branch("jet_pt",           &t_jet_pt);
        evTree->Branch("jet_eta",          &t_jet_eta);
        evTree->Branch("jet_phi",          &t_jet_phi);
        evTree->Branch("photon_pt",        &t_photon_pt);
        evTree->Branch("photon_eta",       &t_photon_eta);
        evTree->Branch("photon_phi",       &t_photon_phi);
        evTree->Branch("photon_is_prompt", &t_photon_is_prompt);

        // Init Pythia for this bin
        Pythia8::Pythia pythia(xmldoc);

        pythia.settings.parm("Beams:eCM", sqrtS);
        pythia.settings.mode("Beams:idA", 2212);
        pythia.settings.mode("Beams:idB", 2212);
        pythia.readString("HardQCD:all = on");
        pythia.readString("PromptPhoton:all = on");
        pythia.settings.parm("PhaseSpace:pTHatMin", ptLow);
        pythia.settings.parm("PhaseSpace:pTHatMax", ptHigh);
        pythia.readString("Print:quiet = on");

        if (!pythia.init()) {
            std::cerr << "ERROR: Pythia init failed for bin " << iBin << std::endl;
            return;
        }

        Long64_t nAccepted = 0;

        for (int iEv = 0; iEv < nEvents; ++iEv) {
            if (!pythia.next()) continue;

            t_event_weight = 1.0;
            t_pthat        = pythia.info.pTHat();

            t_jet_pt.clear();  t_jet_eta.clear();  t_jet_phi.clear();
            t_photon_pt.clear(); t_photon_eta.clear(); t_photon_phi.clear();
            t_photon_is_prompt.clear();

            std::vector<fastjet::PseudoJet> particles;
            const Pythia8::Event& event = pythia.event;

            for (int i = 0; i < event.size(); ++i) {
                const Pythia8::Particle& p = event[i];
                if (!p.isFinal()) continue;

                particles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));

                if (p.id() == 22) {
                    float pt  = p.pT();
                    float eta = p.eta();
                    if (pt > photonPtMin && std::abs(eta) < photonEtaMax) {
                        t_photon_pt.push_back(pt);
                        t_photon_eta.push_back(eta);
                        t_photon_phi.push_back(p.phi());
                        t_photon_is_prompt.push_back(hasHadronAncestor_bins(event, i) ? 0 : 1);
                    }
                }
            }

            fastjet::ClusterSequence cs(particles, jetDef);
            for (const auto& jet : fastjet::sorted_by_pt(cs.inclusive_jets(jetPtMin))) {
                if (std::abs(jet.rap()) < jetEtaMax) {
                    t_jet_pt.push_back(jet.pt());
                    t_jet_eta.push_back(jet.rap());
                    t_jet_phi.push_back(jet.phi_std());
                }
            }

            evTree->Fill();
            ++nAccepted;
        }

        // Fill meta row
        m_bin_low    = ptLow;
        m_bin_high   = ptHigh;
        m_sigma_gen  = pythia.info.sigmaGen();
        m_weight_sum = pythia.info.weightSum();
        m_n_accepted = nAccepted;
        m_lumi_eff   = (m_sigma_gen > 0.) ? m_weight_sum / m_sigma_gen : 0.;
        metaTree->Fill();

        std::cout << std::scientific << std::setprecision(4)
                  << "Bin " << iBin
                  << " [" << ptLow << "-" << ptHigh << "] GeV"
                  << "  sigma=" << m_sigma_gen << " mb"
                  << "  nAccepted=" << nAccepted
                  << "  lumi_eff=" << m_lumi_eff << " mb^-1"
                  << std::defaultfloat << std::endl;

        pythia.stat();

        outFile->cd();
        evTree->Write("", TObject::kOverwrite);
    }

    outFile->cd();
    metaTree->Write("", TObject::kOverwrite);
    outFile->Close();

    std::cout << "Done. Output: output_bins.root" << std::endl;
}
