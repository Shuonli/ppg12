// gen_reweighted.C
// Pass 1: Mode 5-style reweighted Pythia8 generation (ROOT macro)
// Sub-run 0: SoftQCD:nonDiffractive, pThat 0-20 GeV, no bias
// Sub-run 1: HardQCD:all + PromptPhoton:all, pThat >= 20 GeV, pT^4 bias
//
// Usage: root -l -b -q 'gen_reweighted.C+(10000)'
//        root -l -b -q gen_reweighted.C+    (default 10000 events)

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
bool hasHadronAncestor_rw(const Pythia8::Event& event, int idx)
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

void gen_reweighted(int nEvents = 10000)
{
    // ---------------------------------------------------------------
    // Configuration
    // ---------------------------------------------------------------
    const double sqrtS        = 200.;
    const double pThatSplit   = 20.;
    const double jetR         = 0.4;
    const double jetPtMin     = 5.;
    const double jetEtaMax    = 1.1;
    const double photonPtMin  = 5.;
    const double photonEtaMax = 1.1;

    // Pythia8 xmldoc path (PYTHIA8 env var points to share/Pythia8)
    const char* pythia8share = gSystem->Getenv("PYTHIA8");
    std::string xmldoc = pythia8share ? std::string(pythia8share) + "/xmldoc" : "";

    // ---------------------------------------------------------------
    // Output ROOT file
    // ---------------------------------------------------------------
    TFile* outFile = TFile::Open("output_reweighted.root", "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "ERROR: cannot open output_reweighted.root" << std::endl;
        return;
    }

    // ---------------------------------------------------------------
    // TTree branches (shared across both sub-runs)
    // ---------------------------------------------------------------
    TTree* evTree = new TTree("events", "Reweighted Pythia8 events");

    Int_t    t_bin_id;
    Double_t t_event_weight;
    Double_t t_pthat;
    std::vector<float> t_jet_pt, t_jet_eta, t_jet_phi;
    std::vector<float> t_photon_pt, t_photon_eta, t_photon_phi;
    std::vector<int>   t_photon_is_prompt;

    evTree->Branch("bin_id",           &t_bin_id,          "bin_id/I");
    evTree->Branch("event_weight",     &t_event_weight,    "event_weight/D");
    evTree->Branch("pthat",            &t_pthat,           "pthat/D");
    evTree->Branch("jet_pt",           &t_jet_pt);
    evTree->Branch("jet_eta",          &t_jet_eta);
    evTree->Branch("jet_phi",          &t_jet_phi);
    evTree->Branch("photon_pt",        &t_photon_pt);
    evTree->Branch("photon_eta",       &t_photon_eta);
    evTree->Branch("photon_phi",       &t_photon_phi);
    evTree->Branch("photon_is_prompt", &t_photon_is_prompt);

    // Meta tree: one row per sub-run, filled after event loop
    TTree* metaTree = new TTree("meta", "Per-sub-run metadata");
    Int_t    m_bin_id;
    Double_t m_sigma_gen, m_weight_sum, m_lumi_eff;
    Long64_t m_n_accepted;

    metaTree->Branch("bin_id",     &m_bin_id,     "bin_id/I");
    metaTree->Branch("sigma_gen",  &m_sigma_gen,  "sigma_gen/D");
    metaTree->Branch("weight_sum", &m_weight_sum, "weight_sum/D");
    metaTree->Branch("n_accepted", &m_n_accepted, "n_accepted/L");
    metaTree->Branch("lumi_eff",   &m_lumi_eff,   "lumi_eff/D");

    // ---------------------------------------------------------------
    // FastJet
    // ---------------------------------------------------------------
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);

    // ---------------------------------------------------------------
    // Loop over sub-runs
    // ---------------------------------------------------------------
    for (int iBin = 0; iBin < 2; ++iBin) {
        Pythia8::Pythia pythia(xmldoc);

        pythia.settings.parm("Beams:eCM", sqrtS);
        pythia.settings.mode("Beams:idA", 2212);
        pythia.settings.mode("Beams:idB", 2212);

        if (iBin == 0) {
            // --- Soft sub-run ---
            pythia.readString("SoftQCD:nonDiffractive = on");
        } else {
            // --- Hard sub-run ---
            pythia.readString("HardQCD:all = on");
            pythia.readString("PromptPhoton:all = on");
            pythia.settings.parm("PhaseSpace:pTHatMin", pThatSplit);
            pythia.readString("PhaseSpace:bias2Selection = on");
            pythia.readString("PhaseSpace:bias2SelectionPow = 4");
            pythia.settings.parm("PhaseSpace:bias2SelectionRef", pThatSplit);
        }

        pythia.readString("Print:quiet = on");
        if (!pythia.init()) {
            std::cerr << "ERROR: Pythia init failed for bin " << iBin << std::endl;
            return;
        }

        Long64_t nAccepted = 0;

        for (int iEv = 0; iEv < nEvents; ++iEv) {
            if (!pythia.next()) continue;

            double pThat = pythia.info.pTHat();
            // Soft-bin overlap removal
            if (iBin == 0 && pThat > pThatSplit) continue;

            t_bin_id       = iBin;
            t_event_weight = pythia.info.weight();
            t_pthat        = pThat;

            t_jet_pt.clear();  t_jet_eta.clear();  t_jet_phi.clear();
            t_photon_pt.clear(); t_photon_eta.clear(); t_photon_phi.clear();
            t_photon_is_prompt.clear();

            // Collect final-state particles
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
                        t_photon_is_prompt.push_back(hasHadronAncestor_rw(event, i) ? 0 : 1);
                    }
                }
            }

            // Jet finding
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
        m_bin_id     = iBin;
        m_sigma_gen  = pythia.info.sigmaGen();
        m_weight_sum = pythia.info.weightSum();
        m_n_accepted = nAccepted;
        m_lumi_eff   = (m_sigma_gen > 0.) ? m_weight_sum / m_sigma_gen : 0.;
        metaTree->Fill();

        std::cout << std::scientific << std::setprecision(4)
                  << "Sub-run " << iBin
                  << "  sigma=" << m_sigma_gen << " mb"
                  << "  weightSum=" << m_weight_sum
                  << "  nAccepted=" << nAccepted
                  << "  lumi_eff=" << m_lumi_eff << " mb^-1"
                  << std::defaultfloat << std::endl;

        pythia.stat();
    }

    outFile->cd();
    evTree->Write("", TObject::kOverwrite);
    metaTree->Write("", TObject::kOverwrite);
    outFile->Close();

    std::cout << "Done. Output: output_reweighted.root" << std::endl;
}
