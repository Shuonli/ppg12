// =============================================================================
// TruthSpectrumCheck.C
// -----------------------------------------------------------------------------
// Sanity-check macro for the PPG12 cross-section pipeline.
// Fills truth photon / jet pT spectra for a single MC sample, weighted by
// weight = xsec_pb / N_evt_read, so the per-sample luminosity cancels and
// overlays of stitched samples are meaningful (d\sigma/dp_T in pb/GeV).
//
// Output:
//   results/truth_spectrum_{filetype}.root containing
//     h_truth_photon_pt               (every truth photon, |eta|<0.7)
//     h_max_truth_photon_pt           (per-event max truth photon, |eta|<0.7)
//     h_max_truth_photon_pt_filtered  (same as above, restricted to declared
//                                       photon_pt_lower/upper window -- mirrors
//                                       the cut applied by the efficiency tool
//                                       so cross-sample stitching is visible)
//     h_truth_jet_pt                  (every truth jet)
//     h_max_truth_jet_pt              (per-event max truth jet)
//     h_max_truth_jet_pt_filtered     (same as above, restricted to declared
//                                       jet_pt_lower/upper window)
//   + TNamed metadata: filetype, xsec_pb, nevt_read, weight, jet_branch_scheme,
//     photon_pt_lower/upper, jet_pt_lower/upper
//
// Branch auto-detection:
//   Current "combined.root" production uses *_AntiKt_Truth_r04 names.
//   Legacy files (e.g. stale jet10/jet15 combined_1214.root) use njet_truth /
//   jet_truth_Pt. The macro probes for the r04 branch first and falls back.
//
// Usage:
//   root -l -b -q 'TruthSpectrumCheck.C("photon5","results")'
// =============================================================================

#include <iostream>
#include <string>
#include <sys/stat.h>

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TNamed.h>
#include <TMath.h>

#include "CrossSectionWeights.h"

using namespace std;
using namespace PPG12;

// Resolve absolute xsec (pb) for a given filetype. The existing
// PPG12::GetSampleConfig weights are normalized to photon20cross / jet50cross;
// for a d\sigma/dp_T overlay we need the absolute cross section.
static float xsec_pb_for(const string &ft)
{
    if (ft == "photon5"  || ft == "photon5_double"  || ft == "photon5_nom")  return photon5cross;
    if (ft == "photon10" || ft == "photon10_double" || ft == "photon10_nom") return photon10cross;
    if (ft == "photon20" || ft == "photon20_double" || ft == "photon20_nom") return photon20cross;
    if (ft == "jet5")  return jet5cross;
    if (ft == "jet8"  || ft == "jet8_double"  || ft == "jet8_nom")  return jet8cross;
    if (ft == "jet10") return jet10cross;
    if (ft == "jet12" || ft == "jet12_double" || ft == "jet12_nom") return jet12cross;
    if (ft == "jet15") return jet15cross;
    if (ft == "jet20" || ft == "jet20_double" || ft == "jet20_nom") return jet20cross;
    if (ft == "jet30" || ft == "jet30_double" || ft == "jet30_nom") return jet30cross;
    if (ft == "jet40" || ft == "jet40_double" || ft == "jet40_nom") return jet40cross;
    if (ft == "jet50") return jet50cross;
    return -1.f;
}

void TruthSpectrumCheck(const string &filetype  = "photon5",
                        const string &outdir    = "results",
                        const string &input_override = "")
{
    cout << "\n========== TruthSpectrumCheck ==========\n";
    cout << "filetype: " << filetype << "\n";

    const float xsec_pb = xsec_pb_for(filetype);
    if (xsec_pb <= 0) {
        cout << "ERROR: unknown filetype '" << filetype << "'\n";
        return;
    }
    cout << "xsec_pb = " << xsec_pb << " pb\n";

    // -- Input file: canonical combined.root under run28, with optional
    //    per-sample override to support stale-stub fallbacks (jet10/jet15).
    string base_dir =
        "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/";
    string input_file = input_override.empty()
        ? (base_dir + filetype + "/condorout/combined.root")
        : input_override;
    cout << "input: " << input_file << "\n";

    TFile *fin = TFile::Open(input_file.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        cout << "ERROR: cannot open input file\n";
        return;
    }
    TTree *tree = (TTree *)fin->Get("slimtree");
    if (!tree) {
        cout << "ERROR: tree 'slimtree' missing in " << input_file << "\n";
        fin->Close();
        return;
    }
    const Long64_t nent = tree->GetEntries();
    if (nent <= 0) {
        cout << "ERROR: empty tree (" << nent << " entries) in " << input_file << "\n";
        fin->Close();
        return;
    }
    cout << "entries: " << nent << "\n";

    // -- Decide truth-jet branch scheme by probing the tree.
    //    Prefer r04 (current production); fall back to legacy names.
    const bool has_r04 = (tree->GetBranch("njet_truth_AntiKt_Truth_r04") != nullptr);
    const bool has_leg = (tree->GetBranch("njet_truth") != nullptr);
    string njet_bname, jetpt_bname, jeteta_bname, jet_scheme;
    if (has_r04) {
        njet_bname   = "njet_truth_AntiKt_Truth_r04";
        jetpt_bname  = "jet_truth_Pt_AntiKt_Truth_r04";
        jeteta_bname = "jet_truth_Eta_AntiKt_Truth_r04";
        jet_scheme   = "r04";
    } else if (has_leg) {
        njet_bname   = "njet_truth";
        jetpt_bname  = "jet_truth_Pt";
        jeteta_bname = "jet_truth_Eta";
        jet_scheme   = "legacy";
    } else {
        cout << "ERROR: no truth-jet branches found\n";
        fin->Close();
        return;
    }
    cout << "truth-jet branch scheme: " << jet_scheme << "\n";

    // Set up reader. We give the reader a raw pointer; its lifetime is tied to
    // this function scope.
    TTreeReader reader(tree);
    TTreeReaderValue<int>   nparticles(reader, "nparticles");
    TTreeReaderArray<float> particle_Pt (reader, "particle_Pt");
    TTreeReaderArray<float> particle_Eta(reader, "particle_Eta");
    TTreeReaderArray<int>   particle_pid(reader, "particle_pid");
    TTreeReaderValue<int>   njet_truth(reader, njet_bname.c_str());
    TTreeReaderArray<float> jet_truth_Pt (reader, jetpt_bname.c_str());
    TTreeReaderArray<float> jet_truth_Eta(reader, jeteta_bname.c_str());

    // -- Histograms: 100 bins, 0-60 GeV for all four.
    const int    nb   = 100;
    const double xlo  = 0.0;
    const double xhi  = 60.0;
    const string yttl = "d#sigma/dp_{T} [pb/GeV]";

    TH1D *h_truth_photon_pt = new TH1D("h_truth_photon_pt",
        Form("Truth photon p_{T} (%s);truth p_{T}^{#gamma} [GeV];%s", filetype.c_str(), yttl.c_str()),
        nb, xlo, xhi);
    h_truth_photon_pt->Sumw2();

    TH1D *h_max_truth_photon_pt = new TH1D("h_max_truth_photon_pt",
        Form("Max truth photon p_{T} per event (%s);max truth p_{T}^{#gamma} [GeV];%s", filetype.c_str(), yttl.c_str()),
        nb, xlo, xhi);
    h_max_truth_photon_pt->Sumw2();

    TH1D *h_truth_jet_pt = new TH1D("h_truth_jet_pt",
        Form("Truth jet p_{T} (%s);truth p_{T}^{jet} [GeV];%s", filetype.c_str(), yttl.c_str()),
        nb, xlo, xhi);
    h_truth_jet_pt->Sumw2();

    TH1D *h_max_truth_jet_pt = new TH1D("h_max_truth_jet_pt",
        Form("Max truth jet p_{T} per event (%s);max truth p_{T}^{jet} [GeV];%s", filetype.c_str(), yttl.c_str()),
        nb, xlo, xhi);
    h_max_truth_jet_pt->Sumw2();

    // Same histograms but ONLY filled when the event falls inside the sample's
    // declared photon / jet pT windows from GetSampleConfig. These mirror the
    // window cuts applied by RecoEffCalculator_TTreeReader.C (see lines 1649
    // and 1686 there), so the combined sum over samples in the overlay should
    // be smooth *if and only if* the declared windows partition truth pT
    // space without overlap. The known CrossSectionWeights.h bugs
    // (photon10_double: [10,100]; jet12_double: [10,100]) will show up as a
    // bump in the combined sum built from these filtered histograms.
    TH1D *h_max_truth_photon_pt_filt = new TH1D("h_max_truth_photon_pt_filtered",
        Form("Max truth #gamma p_{T} (%s) filtered to [photon_pt_lower, photon_pt_upper];max truth p_{T}^{#gamma} [GeV];%s",
             filetype.c_str(), yttl.c_str()),
        nb, xlo, xhi);
    h_max_truth_photon_pt_filt->Sumw2();

    TH1D *h_max_truth_jet_pt_filt = new TH1D("h_max_truth_jet_pt_filtered",
        Form("Max truth jet p_{T} (%s) filtered to [jet_pt_lower, jet_pt_upper];max truth p_{T}^{jet} [GeV];%s",
             filetype.c_str(), yttl.c_str()),
        nb, xlo, xhi);
    h_max_truth_jet_pt_filt->Sumw2();

    // Weight: d\sigma/dp_T in pb/GeV per entry counted.
    //   bin_content = (xsec_pb / N_evt) * (# entries in bin) / bin_width
    // We set weight = (xsec_pb / N_evt) / bin_width so that the overlay is a
    // differential cross section directly.
    const double bin_width = (xhi - xlo) / nb;
    const double weight    = static_cast<double>(xsec_pb) / static_cast<double>(nent) / bin_width;
    cout << "per-entry weight = xsec/N_evt/binw = " << weight << " pb/GeV\n";

    const float eta_cut = 0.7f;

    // Declared windows from GetSampleConfig -- applied to max truth photon /
    // max truth jet to produce the *_filtered histograms.
    PPG12::SampleConfig sc = PPG12::GetSampleConfig(filetype);
    const float ph_lo = sc.valid ? sc.photon_pt_lower : 0.f;
    const float ph_hi = sc.valid ? sc.photon_pt_upper : 999.f;
    const float jt_lo = sc.valid ? sc.jet_pt_lower    : 0.f;
    const float jt_hi = sc.valid ? sc.jet_pt_upper    : 999.f;
    cout << "declared windows: photon[" << ph_lo << "," << ph_hi
         << "]  jet[" << jt_lo << "," << jt_hi << "]\n";

    cout << "looping over " << nent << " entries...\n";
    Long64_t ie = 0;
    while (reader.Next()) {
        ++ie;
        if (ie % 1000000 == 0) cout << "  " << ie << " / " << nent << "\n";

        // Truth photons (|eta| < 0.7)
        float max_ph_pt = -1.f;
        for (int ip = 0; ip < *nparticles; ++ip) {
            if (particle_pid[ip] != 22) continue;
            if (TMath::Abs(particle_Eta[ip]) > eta_cut) continue;
            const float pt = particle_Pt[ip];
            h_truth_photon_pt->Fill(pt, weight);
            if (pt > max_ph_pt) max_ph_pt = pt;
        }
        if (max_ph_pt > 0) {
            h_max_truth_photon_pt->Fill(max_ph_pt, weight);
            // Pipeline-style window cut on max photon pT (matches
            // RecoEffCalculator_TTreeReader.C line 1649).
            if (max_ph_pt >= ph_lo && max_ph_pt <= ph_hi)
                h_max_truth_photon_pt_filt->Fill(max_ph_pt, weight);
        }

        // Truth jets (no extra eta cut -- anatreemaker already restricts)
        float max_jet_pt = -1.f;
        for (int ij = 0; ij < *njet_truth; ++ij) {
            const float jpt = jet_truth_Pt[ij];
            h_truth_jet_pt->Fill(jpt, weight);
            if (jpt > max_jet_pt) max_jet_pt = jpt;
        }
        if (max_jet_pt > 0) {
            h_max_truth_jet_pt->Fill(max_jet_pt, weight);
            // Pipeline-style window cut on max jet pT (matches
            // RecoEffCalculator_TTreeReader.C line 1686).
            if (max_jet_pt >= jt_lo && max_jet_pt <= jt_hi)
                h_max_truth_jet_pt_filt->Fill(max_jet_pt, weight);
        }
    }

    // -- Ensure outdir exists.
    ::mkdir(outdir.c_str(), 0755);

    const string outname = outdir + "/truth_spectrum_" + filetype + ".root";
    TFile *fout = new TFile(outname.c_str(), "RECREATE");
    h_truth_photon_pt        ->Write();
    h_max_truth_photon_pt    ->Write();
    h_max_truth_photon_pt_filt->Write();
    h_truth_jet_pt           ->Write();
    h_max_truth_jet_pt       ->Write();
    h_max_truth_jet_pt_filt  ->Write();

    TNamed("filetype",      filetype.c_str()).Write();
    TNamed("input_file",    input_file.c_str()).Write();
    TNamed("xsec_pb",       Form("%.6e", xsec_pb)).Write();
    TNamed("nevt_read",     Form("%lld",  nent)).Write();
    TNamed("per_entry_weight_pb_per_GeV", Form("%.6e", weight)).Write();
    TNamed("bin_width_GeV", Form("%.4f", bin_width)).Write();
    TNamed("jet_branch_scheme", jet_scheme.c_str()).Write();
    TNamed("photon_pt_lower", Form("%.3f", ph_lo)).Write();
    TNamed("photon_pt_upper", Form("%.3f", ph_hi)).Write();
    TNamed("jet_pt_lower",    Form("%.3f", jt_lo)).Write();
    TNamed("jet_pt_upper",    Form("%.3f", jt_hi)).Write();

    fout->Close();
    fin->Close();

    cout << "wrote " << outname << "\n";
    cout << "========================================\n";
}
