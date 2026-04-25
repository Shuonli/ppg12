// ============================================================================
// VertexReweightAlt.C
//
// Alternative MC vertex reweighting for the PPG12 isolated-photon cross-
// section analysis using the *truth* vertex instead of the reconstructed
// vertex. The motivation is documented in
//     reports/truth_vertex_reco_check.tex
//
// The nominal reweight uses reco vertex z and is applied in
// efficiencytool/RecoEffCalculator_TTreeReader.C (~L1331-1372). Any MC event
// whose reco vertex falls outside the histogram range gets weight 0, dropping
// ~23% of the sample. In particular MC events where the MBD vertex finder
// failed are stored as vertexz = -9999 and are always dropped.
//
// This macro derives an alternative weight keyed on the truth vertex so that
// events with -9999 reco vertex can be recovered. Two methods are computed
// for each (crossing-angle) x (single/double/mixed) variant:
//
//   (1) one-step direct:      w1(z_t) = Ddata(z_r) / T_truth_mc(z_t)
//       (smoothed, area-normalized ratio bin-by-bin)
//
//   (2) two-step:
//       (a) f(z_r) = Ddata(z_r) / Mreco_mc(z_r)    (reco-space weight)
//       (b) pass over MC, for each event with a valid reco vertex look up
//           f(vertexz), then fill a truth-vertex histogram weighted by that
//           f-value. This gives Ttilde(z_t) - the predicted truth
//           distribution assuming the reco-space reweight is correct.
//       (c) w2(z_t) = Ttilde(z_t) / T_truth_mc(z_t)
//
// The two should agree in the limit where the reco-truth vertex mapping is
// faithful. Differences quantify reco-vertex-finding bias, and crucially w2
// stays defined for events with vertexz = -9999 because the weight is
// indexed on the truth vertex.
//
// Inputs:
//   * data pool: anatreemaker/macro_maketree/data/ana521/condorout/
//                  part_*_with_bdt_split.root  (77 files, slimtree)
//   * MC single: anatreemaker/macro_maketree/sim/run28/jet12/condorout/combined.root
//   * MC double: anatreemaker/macro_maketree/sim/run28/jet12_double/condorout/combined.root
//
// Crossing-angle run ranges:
//   0 mrad   : 47289..51274
//   1.5 mrad : 51274..54000
//
// Cluster-weighted double-interaction fractions:
//   0 mrad   : single 0.776, double 0.224
//   1.5 mrad : single 0.921, double 0.079
//
// Preselection: slimtree cluster_Et_CLUSTERINFO_CEMC contains at least one
// entry with Et >= 7 GeV. Same cut for data and MC. NOTE: MC events with
// reco vertex = -9999 have cluster Et recomputed with the wrong z and thus
// collapse to ~0 GeV; those events therefore automatically fail the
// cluster-Et>=7 preselection. This is expected. We still fill the hM_all /
// hT_all sanity histograms separately (those use the same Et>=7 cut) so the
// "all" flavor only catches events with a valid reco vertex in practice;
// truly-missing reco vertices show up as entries with cluster Et collapsed
// near 0 and do not pass preselection.
//
// Output: results/vertex_reweight_alt_jet.root
//
// Compile / run:
//   cd /sphenix/user/shuhangli/ppg12
//   root -l -b -q 'efficiencytool/VertexReweightAlt.C+'
// ============================================================================

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

// ---------- configuration ----------
constexpr double kEtMin         = 7.0;    // preselection cluster Et threshold (GeV)
constexpr double kVtxLo         = -200.0; // histogram range (cm)
constexpr double kVtxHi         =  200.0;
constexpr int    kVtxNBins      = 80;     // 5 cm bins
constexpr double kRecoVtxSentinel = -999.0; // consider vertexz > kRecoVtxSentinel "valid"

// Per-crossing-angle run ranges
constexpr int kRun0mradMin   = 47289;
constexpr int kRun0mradMax   = 51274;
constexpr int kRun1p5mradMin = 51274;
constexpr int kRun1p5mradMax = 54000;

// Cluster-weighted double-interaction fractions
constexpr double kFracSingle0mrad   = 0.776;
constexpr double kFracDouble0mrad   = 0.224;
constexpr double kFracSingle1p5mrad = 0.921;
constexpr double kFracDouble1p5mrad = 0.079;

// Variants: 6 total = {0mrad,1p5mrad} x {single,double,mixed}
const std::vector<std::string> kVariants = {
    "single_0mrad",  "double_0mrad",  "mixed_0mrad",
    "single_1p5mrad","double_1p5mrad","mixed_1p5mrad"
};

// Return the mix weight for a given MC file ("jet12" = single, "jet12_double"
// = double) under a given variant. Return negative value => skip this event
// for this variant.
double mix_weight_for(const std::string& variant, bool is_double)
{
    // singles
    if (variant == "single_0mrad" || variant == "single_1p5mrad") {
        return is_double ? -1.0 : 1.0;
    }
    // doubles
    if (variant == "double_0mrad" || variant == "double_1p5mrad") {
        return is_double ? 1.0 : -1.0;
    }
    // mixed
    if (variant == "mixed_0mrad") {
        return is_double ? kFracDouble0mrad : kFracSingle0mrad;
    }
    if (variant == "mixed_1p5mrad") {
        return is_double ? kFracDouble1p5mrad : kFracSingle1p5mrad;
    }
    return -1.0;
}

// Return true if a given data run falls in the crossing-angle bucket
// implied by variant suffix.
bool run_in_variant(const std::string& variant, int runnumber)
{
    const bool is_0mrad   = variant.find("0mrad")   != std::string::npos
                         && variant.find("1p5mrad") == std::string::npos;
    const bool is_1p5mrad = variant.find("1p5mrad") != std::string::npos;
    if (is_0mrad) {
        return (runnumber >= kRun0mradMin && runnumber < kRun0mradMax);
    }
    if (is_1p5mrad) {
        return (runnumber >= kRun1p5mradMin && runnumber < kRun1p5mradMax);
    }
    return false;
}

// Area-normalize to unit integral (if nonzero). Options "width" NOT used so
// that downstream Interpolate returns comparable values across bin widths
// when all histograms share the same binning (they do).
void area_normalize(TH1D* h)
{
    if (!h) return;
    const double integ = h->Integral();
    if (integ > 0.0) h->Scale(1.0 / integ);
}

// Simple one-pass 3-bin running-average smoothing.
// Preserves total integral approximately (edges use 2-bin average). We create
// a copy of bin contents, then overwrite. Errors are downweighted similarly.
void smooth3(TH1D* h)
{
    if (!h) return;
    const int nb = h->GetNbinsX();
    std::vector<double> c(nb + 2, 0.0);
    std::vector<double> e(nb + 2, 0.0);
    for (int i = 1; i <= nb; ++i) {
        c[i] = h->GetBinContent(i);
        e[i] = h->GetBinError(i);
    }
    for (int i = 1; i <= nb; ++i) {
        if (i == 1) {
            h->SetBinContent(i, 0.5 * (c[1] + c[2]));
            h->SetBinError  (i, 0.5 * std::hypot(e[1], e[2]));
        } else if (i == nb) {
            h->SetBinContent(i, 0.5 * (c[nb - 1] + c[nb]));
            h->SetBinError  (i, 0.5 * std::hypot(e[nb - 1], e[nb]));
        } else {
            h->SetBinContent(i, (c[i - 1] + c[i] + c[i + 1]) / 3.0);
            h->SetBinError  (i, std::sqrt(e[i - 1]*e[i - 1] + e[i]*e[i]
                                          + e[i + 1]*e[i + 1]) / 3.0);
        }
    }
}

// Bin-by-bin ratio num/den with den=0 protection.
// Assumes same binning; result lives in a clone of num.
TH1D* safe_ratio(const TH1D* num, const TH1D* den, const std::string& newname)
{
    TH1D* r = static_cast<TH1D*>(num->Clone(newname.c_str()));
    r->SetDirectory(nullptr);
    r->Reset();
    const int nb = num->GetNbinsX();
    for (int i = 1; i <= nb; ++i) {
        const double n = num->GetBinContent(i);
        const double d = den->GetBinContent(i);
        if (d > 0.0) {
            r->SetBinContent(i, n / d);
            // propagate relative errors in quadrature
            const double en = num->GetBinError(i);
            const double ed = den->GetBinError(i);
            double rel2 = 0.0;
            if (n > 0.0) rel2 += (en / n) * (en / n);
            rel2 += (ed / d) * (ed / d);
            r->SetBinError(i, (n / d) * std::sqrt(rel2));
        } else {
            r->SetBinContent(i, 0.0);
            r->SetBinError  (i, 0.0);
        }
    }
    return r;
}

// Interpolate, clipped to histogram range; outside => 0.
double safe_interp(const TH1D* h, double x)
{
    if (!h) return 0.0;
    const double xmin = h->GetXaxis()->GetXmin();
    const double xmax = h->GetXaxis()->GetXmax();
    // Interpolate uses bin-center anchors, so clip a hair inside the edges
    const double eps = 1e-6;
    if (x <= xmin + eps || x >= xmax - eps) return 0.0;
    return const_cast<TH1D*>(h)->Interpolate(x);
}

TH1D* make_hist(const std::string& name)
{
    TH1D* h = new TH1D(name.c_str(), name.c_str(), kVtxNBins, kVtxLo, kVtxHi);
    h->Sumw2();
    h->SetDirectory(nullptr);
    return h;
}

TH2D* make_hist2(const std::string& name)
{
    // x = reco vertex, y = truth vertex, same binning as 1D histograms.
    TH2D* h = new TH2D(name.c_str(), name.c_str(),
                       kVtxNBins, kVtxLo, kVtxHi,
                       kVtxNBins, kVtxLo, kVtxHi);
    h->Sumw2();
    h->SetDirectory(nullptr);
    return h;
}

// D'Agostini iterative Bayesian unfolding.
// Inputs (all same 1D binning kVtxNBins):
//   hD      : data reco distribution (area-normalized)
//   hT0     : MC truth prior (area-normalized)
//   hR      : MC 2D joint (x=z_r, y=z_t), raw counts; rows get column-normalized
//             internally to give P(z_r | z_t).
//   n_iters : number of iterations.
//   out_name: name for the returned weight histogram.
// Returns: w_iter(z_t) = T_final(z_t) / T0(z_t), ready to be applied as a
// truth-vertex reweight. Output is NOT re-normalized (only the shape ratio
// matters for a weight).
TH1D* iterate_dagostini(const TH1D* hD, const TH1D* hT0, const TH2D* hR,
                        int n_iters, const std::string& out_name)
{
    if (!hD || !hT0 || !hR) return nullptr;
    const int NZ = kVtxNBins;

    // Extract response matrix P(zr_i | zt_j) = R[i][j] / sum_k R[k][j]
    // We treat the TH2D as (zr, zt) → (x, y).
    // Note: TH2D binning is 1-indexed.
    std::vector<std::vector<double>> P(NZ, std::vector<double>(NZ, 0.0));
    for (int j = 1; j <= NZ; ++j) {
        double col_sum = 0.0;
        for (int i = 1; i <= NZ; ++i) col_sum += hR->GetBinContent(i, j);
        if (col_sum <= 0.0) continue;
        for (int i = 1; i <= NZ; ++i) {
            P[i-1][j-1] = hR->GetBinContent(i, j) / col_sum;
        }
    }

    // Initial truth prior (copy, normalized)
    std::vector<double> T(NZ, 0.0);
    std::vector<double> D(NZ, 0.0);
    double Dsum = 0.0, Tsum = 0.0;
    for (int j = 1; j <= NZ; ++j) T[j-1] = hT0->GetBinContent(j);
    for (int i = 1; i <= NZ; ++i) D[i-1] = hD->GetBinContent(i);
    for (auto x : T) Tsum += x;
    for (auto x : D) Dsum += x;
    if (Tsum <= 0.0 || Dsum <= 0.0) return nullptr;
    for (auto& x : T) x /= Tsum;  // normalize prior
    for (auto& x : D) x /= Dsum;  // normalize target

    // Keep a copy of the initial prior for the weight construction
    std::vector<double> T0 = T;

    for (int iter = 0; iter < n_iters; ++iter) {
        // Forward: M_i = sum_j P[i][j] * T[j]
        std::vector<double> M(NZ, 0.0);
        for (int i = 0; i < NZ; ++i) {
            double s = 0.0;
            for (int j = 0; j < NZ; ++j) s += P[i][j] * T[j];
            M[i] = s;
        }
        // Update: T_new[j] = T[j] * sum_i P[i][j] * D[i] / M[i]
        std::vector<double> Tnew(NZ, 0.0);
        for (int j = 0; j < NZ; ++j) {
            double s = 0.0;
            for (int i = 0; i < NZ; ++i) {
                if (M[i] > 0.0) s += P[i][j] * D[i] / M[i];
            }
            Tnew[j] = T[j] * s;
        }
        T.swap(Tnew);
    }

    // Build the weight w_iter(z_t) = T_final(z_t) / T0(z_t)
    TH1D* w = static_cast<TH1D*>(hT0->Clone(out_name.c_str()));
    w->SetDirectory(nullptr);
    w->Reset();
    for (int j = 1; j <= NZ; ++j) {
        if (T0[j-1] > 0.0) {
            w->SetBinContent(j, T[j-1] / T0[j-1]);
            w->SetBinError(j, 0.0); // statistical error from unfolding not propagated here
        }
    }
    return w;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
void VertexReweightAlt()
{
    TRandom3 rng(42); // unused but honoring the convention
    TStopwatch watch;
    watch.Start();

    const std::string kRepoRoot = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12";
    const std::string kDataGlob = kRepoRoot
        + "/anatreemaker/macro_maketree/data/ana521/condorout/part_*_with_bdt_split.root";
    const std::string kMCSingle = kRepoRoot
        + "/anatreemaker/macro_maketree/sim/run28/jet12/condorout/combined.root";
    const std::string kMCDouble = kRepoRoot
        + "/anatreemaker/macro_maketree/sim/run28/jet12_double/condorout/combined.root";
    const std::string kOutPath  = kRepoRoot
        + "/results/vertex_reweight_alt_jet.root";

    // ---- book histograms ----
    std::map<std::string, TH1D*> hD;        // data reco vertex
    std::map<std::string, TH1D*> hM;        // MC reco vertex, MC events with valid reco vtx
    std::map<std::string, TH1D*> hT;        // MC truth vertex, same selection as hM
    std::map<std::string, TH1D*> hM_all;    // MC reco vertex, all MC events (preselection only)
    std::map<std::string, TH1D*> hT_all;    // MC truth vertex, all MC events (preselection only)
    std::map<std::string, TH1D*> hTtilde;   // MC truth vertex, f(z_r)-weighted (2-step intermediate)
    std::map<std::string, TH1D*> hM_closure_w1;   // MC reco vertex after applying w1(z_t)
    std::map<std::string, TH1D*> hM_closure_w2;   // MC reco vertex after applying w2(z_t)
    std::map<std::string, TH1D*> hM_closure_iter; // MC reco vertex after applying w_iter(z_t)
    std::map<std::string, TH2D*> hR;        // MC joint (z_r, z_t), mix-weighted, valid reco

    for (const auto& v : kVariants) {
        hD             [v] = make_hist ("hD_"              + v);
        hM             [v] = make_hist ("hM_"              + v);
        hT             [v] = make_hist ("hT_"              + v);
        hM_all         [v] = make_hist ("hM_all_"          + v);
        hT_all         [v] = make_hist ("hT_all_"          + v);
        hTtilde        [v] = make_hist ("hTtilde_"         + v);
        hM_closure_w1  [v] = make_hist ("hM_closure_w1_"   + v);
        hM_closure_w2  [v] = make_hist ("hM_closure_w2_"   + v);
        hM_closure_iter[v] = make_hist ("hM_closure_iter_" + v);
        hR             [v] = make_hist2("hR_"              + v);
    }

    // ===================================================================
    // Pass 0: data
    // ===================================================================
    std::cout << "[VertexReweightAlt] === reading data ===" << std::endl;
    TChain dataChain("slimtree");
    const int nAdded = dataChain.Add(kDataGlob.c_str());
    std::cout << "[VertexReweightAlt] data: added " << nAdded
              << " files from " << kDataGlob << std::endl;

    // Limit branch reads for speed. We cannot use SetBranchStatus on a chain
    // before the first GetEntry; TTreeReader + a helper tree does this for
    // us. Here we do explicit SetBranchStatus on the chain.
    dataChain.SetBranchStatus("*", 0);
    dataChain.SetBranchStatus("vertexz", 1);
    dataChain.SetBranchStatus("runnumber", 1);
    dataChain.SetBranchStatus("cluster_Et_CLUSTERINFO_CEMC", 1);

    TTreeReader dReader(&dataChain);
    TTreeReaderValue<float>      dVtxZ(dReader, "vertexz");
    TTreeReaderValue<int>        dRun (dReader, "runnumber");
    TTreeReaderArray<float>      dEt  (dReader, "cluster_Et_CLUSTERINFO_CEMC");

    Long64_t nDataTotal = dataChain.GetEntries();
    std::cout << "[VertexReweightAlt] data total entries: " << nDataTotal << std::endl;
    const Long64_t dPrintEvery = std::max<Long64_t>(1, nDataTotal / 20);

    Long64_t iD = 0;
    Long64_t nDataPass = 0;
    while (dReader.Next()) {
        ++iD;
        if (iD % dPrintEvery == 0) {
            std::cout << "  data " << iD << " / " << nDataTotal << std::endl;
        }
        // Et preselection
        bool pass_et = false;
        for (size_t k = 0; k < dEt.GetSize(); ++k) {
            if (dEt[k] >= kEtMin) { pass_et = true; break; }
        }
        if (!pass_et) continue;
        const double vz = *dVtxZ;
        const int    rn = *dRun;
        for (const auto& v : kVariants) {
            if (run_in_variant(v, rn)) {
                hD[v]->Fill(vz, 1.0);
            }
        }
        ++nDataPass;
    }
    std::cout << "[VertexReweightAlt] data events passing Et>=7: "
              << nDataPass << std::endl;

    // ===================================================================
    // Pass 1a/b: MC single (jet12) and MC double (jet12_double)
    // Fill hM, hT, hM_all, hT_all for every variant the file contributes to.
    // ===================================================================
    auto fillMCPass1 = [&](const std::string& path, bool is_double)
    {
        const std::string tag = is_double ? "jet12_double" : "jet12";
        std::cout << "[VertexReweightAlt] === MC pass1: " << tag << " ===" << std::endl;
        TFile* fin = TFile::Open(path.c_str(), "READ");
        if (!fin || fin->IsZombie()) {
            std::cerr << "  ERROR: cannot open " << path << std::endl;
            return;
        }
        TTree* tr = dynamic_cast<TTree*>(fin->Get("slimtree"));
        if (!tr) {
            std::cerr << "  ERROR: no slimtree in " << path << std::endl;
            fin->Close();
            return;
        }
        tr->SetBranchStatus("*", 0);
        tr->SetBranchStatus("vertexz", 1);
        tr->SetBranchStatus("vertexz_truth", 1);
        tr->SetBranchStatus("cluster_Et_CLUSTERINFO_CEMC", 1);

        TTreeReader            reader(tr);
        TTreeReaderValue<float> vz   (reader, "vertexz");
        TTreeReaderValue<float> vztru(reader, "vertexz_truth");
        TTreeReaderArray<float> et   (reader, "cluster_Et_CLUSTERINFO_CEMC");

        const Long64_t nTot = tr->GetEntries();
        const Long64_t printEvery = std::max<Long64_t>(1, nTot / 20);
        std::cout << "  entries: " << nTot << std::endl;

        Long64_t ie = 0, nPass = 0, nPassValid = 0;
        while (reader.Next()) {
            ++ie;
            if (ie % printEvery == 0) {
                std::cout << "  " << tag << " pass1 " << ie << " / " << nTot << std::endl;
            }
            bool pass_et = false;
            for (size_t k = 0; k < et.GetSize(); ++k) {
                if (et[k] >= kEtMin) { pass_et = true; break; }
            }
            if (!pass_et) continue;
            ++nPass;
            const double zr = *vz;
            const double zt = *vztru;
            const bool valid_reco = (zr > kRecoVtxSentinel);
            if (valid_reco) ++nPassValid;
            for (const auto& v : kVariants) {
                const double w = mix_weight_for(v, is_double);
                if (w < 0.0) continue;
                hM_all[v]->Fill(zr, w);
                hT_all[v]->Fill(zt, w);
                if (valid_reco) {
                    hM[v]->Fill(zr, w);
                    hT[v]->Fill(zt, w);
                    hR[v]->Fill(zr, zt, w);
                }
            }
        }
        std::cout << "  " << tag << ": pass Et=" << nPass
                  << " of which valid reco=" << nPassValid << std::endl;
        fin->Close();
        delete fin;
    };

    fillMCPass1(kMCSingle, /*is_double=*/false);
    fillMCPass1(kMCDouble, /*is_double=*/true);

    // ===================================================================
    // Validation
    // ===================================================================
    for (const auto& v : kVariants) {
        if (hD[v]->Integral() <= 0.0) {
            std::cout << "[VertexReweightAlt] WARN: hD_" << v
                      << " has zero integral" << std::endl;
        }
        if (hT[v]->Integral() <= 0.0) {
            std::cout << "[VertexReweightAlt] WARN: hT_" << v
                      << " has zero integral" << std::endl;
        }
    }

    // ===================================================================
    // Build f(z_r) = Ddata / Mreco(normalized) per variant, and
    // w1(z_t) = Ddata / Ttruth(normalized) per variant.
    // Both sides are smoothed and area-normalized BEFORE the ratio.
    // ===================================================================
    std::map<std::string, TH1D*> h_f;    // reco-space weight
    std::map<std::string, TH1D*> h_w1;   // one-step direct truth weight
    std::map<std::string, TH1D*> h_w2;   // two-step truth weight (built later)

    // Keep smoothed+normalized copies around for closure pass.
    std::map<std::string, TH1D*> hD_sn;
    std::map<std::string, TH1D*> hM_sn;
    std::map<std::string, TH1D*> hT_sn;

    for (const auto& v : kVariants) {
        TH1D* hDs = static_cast<TH1D*>(hD[v]->Clone(("hD_sn_" + v).c_str()));
        TH1D* hMs = static_cast<TH1D*>(hM[v]->Clone(("hM_sn_" + v).c_str()));
        TH1D* hTs = static_cast<TH1D*>(hT[v]->Clone(("hT_sn_" + v).c_str()));
        hDs->SetDirectory(nullptr);
        hMs->SetDirectory(nullptr);
        hTs->SetDirectory(nullptr);
        smooth3(hDs);
        smooth3(hMs);
        smooth3(hTs);
        area_normalize(hDs);
        area_normalize(hMs);
        area_normalize(hTs);
        hD_sn[v] = hDs;
        hM_sn[v] = hMs;
        hT_sn[v] = hTs;

        h_f [v] = safe_ratio(hDs, hMs, "h_f_"  + v);
        h_w1[v] = safe_ratio(hDs, hTs, "h_w1_" + v);
    }

    // ===================================================================
    // Pass 2: re-read MC to build Ttilde and closure histograms.
    // For each MC event with valid reco vertex + Et>=7:
    //   fval = f_<variant>(vertexz)  (clipped to range)
    //   hTtilde[v]->Fill(vertexz_truth, fval * mix_weight)
    // Then build w2 = hTtilde / hT (smoothed, normalized).
    //
    // For closure, we also need to apply w1 and w2 to each MC event (keyed
    // on vertexz_truth) and fill the *reco* vertex distribution so we can
    // compare to data later. We do that in-line here to save another pass.
    // Note that w2 is not yet available; we therefore fill closure_w1 now
    // and closure_w2 in pass 3 after w2 is built.
    // ===================================================================
    auto fillMCPass2 = [&](const std::string& path, bool is_double)
    {
        const std::string tag = is_double ? "jet12_double" : "jet12";
        std::cout << "[VertexReweightAlt] === MC pass2: " << tag << " ===" << std::endl;
        TFile* fin = TFile::Open(path.c_str(), "READ");
        if (!fin || fin->IsZombie()) {
            std::cerr << "  ERROR: cannot open " << path << std::endl;
            return;
        }
        TTree* tr = dynamic_cast<TTree*>(fin->Get("slimtree"));
        if (!tr) {
            std::cerr << "  ERROR: no slimtree in " << path << std::endl;
            fin->Close();
            return;
        }
        tr->SetBranchStatus("*", 0);
        tr->SetBranchStatus("vertexz", 1);
        tr->SetBranchStatus("vertexz_truth", 1);
        tr->SetBranchStatus("cluster_Et_CLUSTERINFO_CEMC", 1);

        TTreeReader            reader(tr);
        TTreeReaderValue<float> vz   (reader, "vertexz");
        TTreeReaderValue<float> vztru(reader, "vertexz_truth");
        TTreeReaderArray<float> et   (reader, "cluster_Et_CLUSTERINFO_CEMC");

        const Long64_t nTot = tr->GetEntries();
        const Long64_t printEvery = std::max<Long64_t>(1, nTot / 20);

        Long64_t ie = 0;
        while (reader.Next()) {
            ++ie;
            if (ie % printEvery == 0) {
                std::cout << "  " << tag << " pass2 " << ie << " / " << nTot << std::endl;
            }
            bool pass_et = false;
            for (size_t k = 0; k < et.GetSize(); ++k) {
                if (et[k] >= kEtMin) { pass_et = true; break; }
            }
            if (!pass_et) continue;
            const double zr = *vz;
            const double zt = *vztru;
            const bool valid_reco = (zr > kRecoVtxSentinel);
            for (const auto& v : kVariants) {
                const double mw = mix_weight_for(v, is_double);
                if (mw < 0.0) continue;
                // Ttilde: uses f(z_r), requires a valid reco vertex
                if (valid_reco) {
                    const double fval = safe_interp(h_f[v], zr);
                    if (fval > 0.0) {
                        hTtilde[v]->Fill(zt, fval * mw);
                    }
                }
                // closure_w1: apply w1(z_t) to all events, fill reco vertex
                const double w1val = safe_interp(h_w1[v], zt);
                if (w1val > 0.0 && valid_reco) {
                    hM_closure_w1[v]->Fill(zr, w1val * mw);
                }
            }
        }
        fin->Close();
        delete fin;
    };

    fillMCPass2(kMCSingle, /*is_double=*/false);
    fillMCPass2(kMCDouble, /*is_double=*/true);

    // Build w2 from hTtilde and hT (both smoothed+normalized before ratio).
    for (const auto& v : kVariants) {
        TH1D* hTtS = static_cast<TH1D*>(hTtilde[v]->Clone(("hTtilde_sn_" + v).c_str()));
        hTtS->SetDirectory(nullptr);
        smooth3(hTtS);
        area_normalize(hTtS);
        // hT_sn already smoothed+normalized from above
        h_w2[v] = safe_ratio(hTtS, hT_sn[v], "h_w2_" + v);
        delete hTtS;
    }

    // ===================================================================
    // Build w_iter (D'Agostini iterative unfolding) per variant.
    // Uses the 2D response matrix hR (z_r, z_t) and the smoothed+normalized
    // 1D hD_sn, hT_sn. 5 iterations is plenty for a 24-bin problem.
    // ===================================================================
    const int kNIter = 5;
    std::map<std::string, TH1D*> h_w_iter;
    for (const auto& v : kVariants) {
        h_w_iter[v] = iterate_dagostini(hD_sn[v], hT_sn[v], hR[v],
                                        kNIter, "h_w_iter_" + v);
        if (h_w_iter[v]) {
            // Smooth the final weight once for display consistency; the shape
            // is dominated by the unfolded T, not by this pass.
            smooth3(h_w_iter[v]);
        } else {
            std::cerr << "[VertexReweightAlt] WARN: unfolding failed for " << v << std::endl;
            h_w_iter[v] = make_hist("h_w_iter_" + v); // empty
        }
    }

    // ===================================================================
    // Pass 3: closure_w2 — apply w2(z_t) to MC, fill reco vertex.
    // ===================================================================
    auto fillMCPass3 = [&](const std::string& path, bool is_double)
    {
        const std::string tag = is_double ? "jet12_double" : "jet12";
        std::cout << "[VertexReweightAlt] === MC pass3: " << tag << " ===" << std::endl;
        TFile* fin = TFile::Open(path.c_str(), "READ");
        if (!fin || fin->IsZombie()) { std::cerr << "  ERROR: open\n"; return; }
        TTree* tr = dynamic_cast<TTree*>(fin->Get("slimtree"));
        if (!tr) { std::cerr << "  ERROR: no slimtree\n"; fin->Close(); return; }
        tr->SetBranchStatus("*", 0);
        tr->SetBranchStatus("vertexz", 1);
        tr->SetBranchStatus("vertexz_truth", 1);
        tr->SetBranchStatus("cluster_Et_CLUSTERINFO_CEMC", 1);

        TTreeReader            reader(tr);
        TTreeReaderValue<float> vz   (reader, "vertexz");
        TTreeReaderValue<float> vztru(reader, "vertexz_truth");
        TTreeReaderArray<float> et   (reader, "cluster_Et_CLUSTERINFO_CEMC");

        const Long64_t nTot = tr->GetEntries();
        const Long64_t printEvery = std::max<Long64_t>(1, nTot / 20);
        Long64_t ie = 0;
        while (reader.Next()) {
            ++ie;
            if (ie % printEvery == 0) {
                std::cout << "  " << tag << " pass3 " << ie << " / " << nTot << std::endl;
            }
            bool pass_et = false;
            for (size_t k = 0; k < et.GetSize(); ++k) {
                if (et[k] >= kEtMin) { pass_et = true; break; }
            }
            if (!pass_et) continue;
            const double zr = *vz;
            const double zt = *vztru;
            const bool valid_reco = (zr > kRecoVtxSentinel);
            if (!valid_reco) continue;
            for (const auto& v : kVariants) {
                const double mw = mix_weight_for(v, is_double);
                if (mw < 0.0) continue;
                const double w2val = safe_interp(h_w2[v], zt);
                if (w2val > 0.0) {
                    hM_closure_w2[v]->Fill(zr, w2val * mw);
                }
                const double witerval = safe_interp(h_w_iter[v], zt);
                if (witerval > 0.0) {
                    hM_closure_iter[v]->Fill(zr, witerval * mw);
                }
            }
        }
        fin->Close();
        delete fin;
    };

    fillMCPass3(kMCSingle, /*is_double=*/false);
    fillMCPass3(kMCDouble, /*is_double=*/true);

    // Normalize closure histograms for easy comparison with data
    for (const auto& v : kVariants) {
        area_normalize(hM_closure_w1[v]);
        area_normalize(hM_closure_w2[v]);
        area_normalize(hM_closure_iter[v]);
    }

    // ===================================================================
    // Print summary table
    // ===================================================================
    auto f_stats = [&](const TH1D* h) {
        double mean = 0.0, rms = 0.0, sumw = 0.0;
        const int nb = h->GetNbinsX();
        for (int i = 1; i <= nb; ++i) {
            const double c = h->GetBinContent(i);
            sumw += c;
            mean += c * h->GetBinCenter(i);
        }
        if (sumw > 0.0) mean /= sumw;
        for (int i = 1; i <= nb; ++i) {
            const double c = h->GetBinContent(i);
            const double d = h->GetBinCenter(i) - mean;
            rms += c * d * d;
        }
        if (sumw > 0.0) rms = std::sqrt(rms / sumw);
        return std::make_pair(mean, rms);
    };

    std::printf("\n==================== VertexReweightAlt summary ====================\n");
    std::printf("%-16s %10s %10s %10s %10s %10s %10s %10s\n",
                "variant", "nD", "nM_reco", "nM_all", "f_mean", "f_rms",
                "w1(z=0)", "w2(z=0)");
    for (const auto& v : kVariants) {
        const double nD     = hD[v]->Integral();
        const double nMreco = hM[v]->Integral();
        const double nMall  = hM_all[v]->Integral();
        auto fs = f_stats(h_f[v]);
        const double w1at0 = safe_interp(h_w1[v], 0.0);
        const double w2at0 = safe_interp(h_w2[v], 0.0);
        std::printf("%-16s %10.3g %10.3g %10.3g %10.4f %10.4f %10.4f %10.4f\n",
                    v.c_str(), nD, nMreco, nMall, fs.first, fs.second,
                    w1at0, w2at0);
    }
    std::printf("===================================================================\n\n");

    // ===================================================================
    // Write all histograms
    // ===================================================================
    // Ensure results directory exists (it normally does, but be defensive)
    gSystem->mkdir((kRepoRoot + "/results").c_str(), true);

    TFile fout(kOutPath.c_str(), "RECREATE");
    if (fout.IsZombie()) {
        std::cerr << "[VertexReweightAlt] ERROR: cannot open output "
                  << kOutPath << std::endl;
        return;
    }

    // Provenance
    TNamed prov("provenance",
        "VertexReweightAlt.C: truth-vertex-based MC reweighting for PPG12. "
        "Inputs: data part_*_with_bdt_split.root (chain), jet12/combined.root, "
        "jet12_double/combined.root. Variants: {0mrad,1p5mrad}x{single,double,mixed}. "
        "Two methods: w1=D/T (direct), w2=Ttilde/T with Ttilde built from f=D/Mreco. "
        "See reports/truth_vertex_reco_check.tex for motivation.");
    prov.Write();

    auto write_map = [&](std::map<std::string, TH1D*>& mp) {
        for (auto& kv : mp) {
            if (kv.second) kv.second->Write(kv.second->GetName(), TObject::kOverwrite);
        }
    };
    // Store raw + intermediates + weights
    write_map(hD);
    write_map(hM);
    write_map(hT);
    write_map(hM_all);
    write_map(hT_all);
    write_map(hTtilde);
    write_map(hD_sn);
    write_map(hM_sn);
    write_map(hT_sn);
    write_map(h_f);
    write_map(h_w1);
    write_map(h_w2);
    write_map(h_w_iter);
    write_map(hM_closure_w1);
    write_map(hM_closure_w2);
    write_map(hM_closure_iter);

    // Also write the response matrices (2D) for provenance
    for (auto& kv : hR) {
        if (kv.second) kv.second->Write(kv.second->GetName(), TObject::kOverwrite);
    }

    fout.Close();
    std::cout << "[VertexReweightAlt] wrote " << kOutPath << std::endl;

    watch.Stop();
    std::cout << "[VertexReweightAlt] elapsed: real=" << watch.RealTime()
              << "s cpu=" << watch.CpuTime() << "s" << std::endl;
}
