// StudyMBDVertexEff.C
//
// Standalone study of MBD-vertex reconstruction efficiency vs leading
// truth pT, separately for the photon channel (signal MC: photon5/10/20)
// and the jet channel (background MC: jet5/8/12/20/30/40).
//
// Numerator (per event):
//   loose   : vertexz != -9999       (the MBD finder produced *any* value)
//   bundled : mbdnorthhit>=1 && mbdsouthhit>=1 && |vertexz|<60 cm
//             (matches the eps_MBD numerator in RecoEffCalculator)
//
// Denominator: all events with a leading truth particle inside the
// per-sample truth-pT window (from PPG12::GetSampleConfig). The window
// cut is applied identically in RecoEffCalculator_TTreeReader.C:1738
// (photon) and :1775/:1780 (jet) to prevent double-counting between
// neighbouring samples.
//   - photon: max particle_Pt over pid==22  (any photonclass, any eta)
//   - jet:    max jet_truth_Pt_AntiKt_Truth_r04 (any eta)
//
// Weighting:
//   per-event base weight = xsec_sample / N_generated (xsec from
//   CrossSectionWeights.h; N_generated = default_nsimevents = 1e7)
//   Optionally multiplied by TruthVertexWeight(h_w_iterative, vertexz_truth)
//   from efficiencytool/truth_vertex_reweight/output/{0mrad,1p5mrad}/reweight.root.
//
// Outputs (efficiencytool/results/mbd_vertex_eff_study.root):
//   h_den_{photon,jet}_{noreweight,0mrad,1p5mrad}
//   h_num_{photon,jet}_{loose,bundled}_{noreweight,0mrad,1p5mrad}
//   g_eff_{photon,jet}_{loose,bundled}_{noreweight,0mrad,1p5mrad}
//
// Run:
//   cd efficiencytool
//   root -l -b -q StudyMBDVertexEff.C+

#include "CrossSectionWeights.h"
#include "TruthVertexReweightLoader.h"

#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TNamed.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>


namespace {

constexpr float kVtxRecoSentinel = -1000.0f;  // anything > this means "MBD vertex found"
constexpr float kVtxRecoFiducial = 60.0f;

const char* kInputDir = "/sphenix/user/shuhangli/ppg12/FunWithxgboost";
const char* kTreeName = "slimtree";
const char* kOutFile  =
    "/sphenix/user/shuhangli/ppg12/efficiencytool/results/mbd_vertex_eff_study.root";

const std::vector<std::string> kPhotonSamples = {
    "photon5", "photon10", "photon20"};
const std::vector<std::string> kJetSamples = {
    "jet5", "jet8", "jet12", "jet20", "jet30", "jet40"};

const std::vector<std::string> kReweightLabels = {
    "noreweight", "0mrad", "1p5mrad"};
const std::map<std::string, std::string> kReweightFiles = {
    {"0mrad",
     "/sphenix/user/shuhangli/ppg12/efficiencytool/truth_vertex_reweight/output/0mrad/reweight.root"},
    {"1p5mrad",
     "/sphenix/user/shuhangli/ppg12/efficiencytool/truth_vertex_reweight/output/1p5mrad/reweight.root"},
};

// Leading-pT binning. 8-36 GeV is the analysis range; extended to 3 GeV
// on the low side so photon5 / jet5 contribute to dedicated bins.
const std::vector<double> kPtEdges = {
    3, 5, 7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36};


float SampleXsec(const std::string& s)
{
    using namespace PPG12;
    if (s == "photon5")  return photon5cross;
    if (s == "photon10") return photon10cross;
    if (s == "photon20") return photon20cross;
    if (s == "jet5")     return jet5cross;
    if (s == "jet8")     return jet8cross;
    if (s == "jet12")    return jet12cross;
    if (s == "jet20")    return jet20cross;
    if (s == "jet30")    return jet30cross;
    if (s == "jet40")    return jet40cross;
    return 0.0f;
}

struct ChannelHists {
    std::map<std::string, TH1F*> den;
    std::map<std::string, TH1F*> num_loose;
    std::map<std::string, TH1F*> num_bundled;
};

ChannelHists MakeChannelHists(const std::string& channel)
{
    const int nbins = static_cast<int>(kPtEdges.size()) - 1;
    ChannelHists ch;
    for (const auto& rw : kReweightLabels) {
        TH1F* h_d = new TH1F(
            Form("h_den_%s_%s", channel.c_str(), rw.c_str()),
            ";leading truth p_{T} [GeV];weighted events",
            nbins, kPtEdges.data());
        TH1F* h_l = new TH1F(
            Form("h_num_%s_loose_%s", channel.c_str(), rw.c_str()),
            ";leading truth p_{T} [GeV];weighted events",
            nbins, kPtEdges.data());
        TH1F* h_b = new TH1F(
            Form("h_num_%s_bundled_%s", channel.c_str(), rw.c_str()),
            ";leading truth p_{T} [GeV];weighted events",
            nbins, kPtEdges.data());
        h_d->Sumw2();  h_l->Sumw2();  h_b->Sumw2();
        ch.den[rw]         = h_d;
        ch.num_loose[rw]   = h_l;
        ch.num_bundled[rw] = h_b;
    }
    return ch;
}

void FillPhotonSample(const std::string& sample, ChannelHists& ch,
                      const std::map<std::string, TH1*>& reweights)
{
    TString path = Form("%s/%s/bdt_split.root", kInputDir, sample.c_str());
    TFile* f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "  ERROR: cannot open " << path << std::endl;
        return;
    }
    TTree* tree = dynamic_cast<TTree*>(f->Get(kTreeName));
    if (!tree) {
        std::cerr << "  ERROR: tree '" << kTreeName << "' missing in "
                  << path << std::endl;
        f->Close();
        return;
    }

    Long64_t nEntries = tree->GetEntries();
    std::cout << "  [" << sample << "] " << path
              << "  nEntries = " << nEntries << std::endl;

    TTreeReader reader(tree);
    TTreeReaderValue<float>     vertexz(reader, "vertexz");
    TTreeReaderValue<float>     vertexz_truth(reader, "vertexz_truth");
    TTreeReaderValue<int>       mbdnorthhit(reader, "mbdnorthhit");
    TTreeReaderValue<int>       mbdsouthhit(reader, "mbdsouthhit");
    TTreeReaderArray<float>     particle_Pt(reader, "particle_Pt");
    TTreeReaderArray<int>       particle_pid(reader, "particle_pid");

    const PPG12::SampleConfig sc = PPG12::GetSampleConfig(sample);
    const float xsec = SampleXsec(sample);
    if (xsec <= 0 || !sc.valid) {
        std::cerr << "  ERROR: unknown sample " << sample << std::endl;
        f->Close();
        return;
    }
    const float base_w  = xsec / PPG12::default_nsimevents;
    const float win_lo  = sc.photon_pt_lower;
    const float win_hi  = sc.photon_pt_upper;

    // Map reweight label -> histogram pointer (or null for "noreweight")
    std::map<std::string, TH1*> rw_map = {{"noreweight", nullptr}};
    for (const auto& [k, h] : reweights) rw_map[k] = h;

    Long64_t nKept = 0, nLoose = 0, nBundled = 0;
    while (reader.Next()) {
        float lead_pt = 0.0f;
        const size_t np = particle_Pt.GetSize();
        for (size_t i = 0; i < np; ++i) {
            if (particle_pid[i] != 22) continue;
            if (particle_Pt[i] > lead_pt) lead_pt = particle_Pt[i];
        }
        if (lead_pt <= 0.0f) continue;
        // Sample truth-pT window (matches RecoEffCalculator_TTreeReader.C:1738)
        if (lead_pt < win_lo || lead_pt > win_hi) continue;

        const bool m_loose   = (*vertexz > kVtxRecoSentinel);
        const bool m_bundled = (*mbdnorthhit >= 1) && (*mbdsouthhit >= 1)
                            && (std::abs(*vertexz) < kVtxRecoFiducial);
        const float vzt = *vertexz_truth;

        for (const auto& rw_label : kReweightLabels) {
            TH1* h_rw = rw_map[rw_label];
            const float rw_factor = (h_rw ? TruthVertexWeight(h_rw, vzt) : 1.0f);
            const float w = base_w * rw_factor;
            ch.den[rw_label]->Fill(lead_pt, w);
            if (m_loose)   ch.num_loose[rw_label]->Fill(lead_pt, w);
            if (m_bundled) ch.num_bundled[rw_label]->Fill(lead_pt, w);
        }

        ++nKept;
        if (m_loose)   ++nLoose;
        if (m_bundled) ++nBundled;
    }

    const double keep_frac = double(nKept) / std::max<Long64_t>(nEntries, 1);
    const double p_loose   = double(nLoose)   / std::max<Long64_t>(nKept, 1);
    const double p_bundled = double(nBundled) / std::max<Long64_t>(nKept, 1);
    printf("    kept %lld / %lld (%.3f); P(loose)=%.4f  P(bundled)=%.4f\n",
           (long long) nKept, (long long) nEntries,
           keep_frac, p_loose, p_bundled);

    f->Close();
}

void FillJetSample(const std::string& sample, ChannelHists& ch,
                   const std::map<std::string, TH1*>& reweights)
{
    TString path = Form("%s/%s/bdt_split.root", kInputDir, sample.c_str());
    TFile* f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "  ERROR: cannot open " << path << std::endl;
        return;
    }
    TTree* tree = dynamic_cast<TTree*>(f->Get(kTreeName));
    if (!tree) {
        std::cerr << "  ERROR: tree '" << kTreeName << "' missing in "
                  << path << std::endl;
        f->Close();
        return;
    }

    Long64_t nEntries = tree->GetEntries();
    std::cout << "  [" << sample << "] " << path
              << "  nEntries = " << nEntries << std::endl;

    TTreeReader reader(tree);
    TTreeReaderValue<float>     vertexz(reader, "vertexz");
    TTreeReaderValue<float>     vertexz_truth(reader, "vertexz_truth");
    TTreeReaderValue<int>       mbdnorthhit(reader, "mbdnorthhit");
    TTreeReaderValue<int>       mbdsouthhit(reader, "mbdsouthhit");
    TTreeReaderArray<float>     jet_pt(reader, "jet_truth_Pt_AntiKt_Truth_r04");

    const PPG12::SampleConfig sc = PPG12::GetSampleConfig(sample);
    const float xsec = SampleXsec(sample);
    if (xsec <= 0 || !sc.valid) {
        std::cerr << "  ERROR: unknown sample " << sample << std::endl;
        f->Close();
        return;
    }
    const float base_w  = xsec / PPG12::default_nsimevents;
    const float win_lo  = sc.jet_pt_lower;
    const float win_hi  = sc.jet_pt_upper;

    std::map<std::string, TH1*> rw_map = {{"noreweight", nullptr}};
    for (const auto& [k, h] : reweights) rw_map[k] = h;

    Long64_t nKept = 0, nLoose = 0, nBundled = 0;
    while (reader.Next()) {
        float lead_pt = 0.0f;
        const size_t nj = jet_pt.GetSize();
        for (size_t i = 0; i < nj; ++i) {
            if (jet_pt[i] > lead_pt) lead_pt = jet_pt[i];
        }
        if (lead_pt <= 0.0f) continue;
        // Sample truth-pT window (matches RecoEffCalculator_TTreeReader.C:1775,1780)
        if (lead_pt < win_lo || lead_pt > win_hi) continue;

        const bool m_loose   = (*vertexz > kVtxRecoSentinel);
        const bool m_bundled = (*mbdnorthhit >= 1) && (*mbdsouthhit >= 1)
                            && (std::abs(*vertexz) < kVtxRecoFiducial);
        const float vzt = *vertexz_truth;

        for (const auto& rw_label : kReweightLabels) {
            TH1* h_rw = rw_map[rw_label];
            const float rw_factor = (h_rw ? TruthVertexWeight(h_rw, vzt) : 1.0f);
            const float w = base_w * rw_factor;
            ch.den[rw_label]->Fill(lead_pt, w);
            if (m_loose)   ch.num_loose[rw_label]->Fill(lead_pt, w);
            if (m_bundled) ch.num_bundled[rw_label]->Fill(lead_pt, w);
        }

        ++nKept;
        if (m_loose)   ++nLoose;
        if (m_bundled) ++nBundled;
    }

    const double keep_frac = double(nKept) / std::max<Long64_t>(nEntries, 1);
    const double p_loose   = double(nLoose)   / std::max<Long64_t>(nKept, 1);
    const double p_bundled = double(nBundled) / std::max<Long64_t>(nKept, 1);
    printf("    kept %lld / %lld (%.3f); P(loose)=%.4f  P(bundled)=%.4f\n",
           (long long) nKept, (long long) nEntries,
           keep_frac, p_loose, p_bundled);

    f->Close();
}

}  // namespace


void StudyMBDVertexEff()
{
    TStopwatch sw; sw.Start();

    // Load reweights
    std::cout << "Loading truth-vertex reweights ..." << std::endl;
    std::map<std::string, TH1*> reweights;
    std::vector<TFile*> rw_files;
    for (const auto& [k, p] : kReweightFiles) {
        TH1* h = LoadTruthVertexReweight(p);
        if (!h) {
            std::cerr << "  WARN: failed to load " << p << std::endl;
            continue;
        }
        printf("  %s: nbins=%d, range [%.1f, %.1f]\n",
               k.c_str(), h->GetNbinsX(),
               h->GetXaxis()->GetXmin(),
               h->GetXaxis()->GetXmax());
        reweights[k] = h;
    }

    ChannelHists ch_photon = MakeChannelHists("photon");
    ChannelHists ch_jet    = MakeChannelHists("jet");

    std::cout << "\n=== Photon channel ===" << std::endl;
    for (const auto& s : kPhotonSamples) FillPhotonSample(s, ch_photon, reweights);

    std::cout << "\n=== Jet channel ===" << std::endl;
    for (const auto& s : kJetSamples) FillJetSample(s, ch_jet, reweights);

    std::cout << "\nWriting " << kOutFile << " ..." << std::endl;
    TFile* fout = TFile::Open(kOutFile, "RECREATE");

    auto write_channel = [&](const std::string& channel,
                             const ChannelHists& ch) {
        for (const auto& rw : kReweightLabels) {
            ch.den.at(rw)->Write();
            ch.num_loose.at(rw)->Write();
            ch.num_bundled.at(rw)->Write();
            for (const auto& def : {"loose", "bundled"}) {
                TH1F* num = (std::string(def) == "loose")
                                ? ch.num_loose.at(rw)
                                : ch.num_bundled.at(rw);
                TGraphAsymmErrors* g = new TGraphAsymmErrors();
                g->BayesDivide(num, ch.den.at(rw));
                g->SetName(Form("g_eff_%s_%s_%s",
                                channel.c_str(), def, rw.c_str()));
                g->Write();
            }
        }
    };
    write_channel("photon", ch_photon);
    write_channel("jet",    ch_jet);

    TNamed prov(
        "provenance",
        "StudyMBDVertexEff.C | photon5/10/20 + jet5/8/12/20/30/40 | "
        "loose: vz_reco != -9999 | "
        "bundled: |vz_reco|<60 cm && MBD N&S>=1 | "
        "denom: events with leading truth particle (gamma_prompt or R=0.4 jet) in |eta|<0.7");
    prov.Write();
    fout->Close();

    std::cout << "\n=== Summary (no-reweight, integrated) ===" << std::endl;
    auto print_sum = [&](const std::string& channel,
                         const ChannelHists& ch) {
        const double den = ch.den.at("noreweight")->Integral();
        const double n_l = ch.num_loose.at("noreweight")->Integral();
        const double n_b = ch.num_bundled.at("noreweight")->Integral();
        if (den > 0) {
            printf("  %-7s loose = %.4f   bundled = %.4f\n",
                   channel.c_str(), n_l / den, n_b / den);
        }
    };
    print_sum("photon", ch_photon);
    print_sum("jet",    ch_jet);

    sw.Stop();
    printf("\nReal time = %.1f s, CPU time = %.1f s\n",
           sw.RealTime(), sw.CpuTime());
    std::cout << "Done." << std::endl;
}
