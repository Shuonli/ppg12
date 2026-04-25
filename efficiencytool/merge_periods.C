// merge_periods.C
//
// Combine two per-crossing-angle MergeSim outputs (e.g. bdt_0rad + bdt_nom)
// into a single full-run-range MC. Per-period MC is assumed to already be
// pre-scaled to the merge-target lumi via the per-event lumi_weight =
// lumi/lumi_target applied in RecoEffCalculator_TTreeReader.C. With the
// scaling baked into per-event Fill weights, plain TFileMerger (= hadd)
// correctly combines TH1/TH2/TH3, TEfficiency, and RooUnfoldResponse via
// each class's native Merge() method.
//
// Validation: hadd verified to sum integrals exactly for TH1, TEfficiency
// (passed/total), and RooUnfoldResponse (Hresponse/Hmeasured/Htruth) in
// session 2026-04-20, with Sumw2 propagated through the upstream weighted
// fills.
//
// File-path derivation matches MergeSim.C output naming:
//   {eff_outfile}_{var_type}.root        (signal photon eff)
//   {eff_outfile}_jet_{var_type}.root    (inclusive jet eff)
//   {response_outfile}_{var_type}.root   (response matrix)
//
// Usage:
//   root -l -b -q 'merge_periods.C("config_bdt_0rad.yaml","config_bdt_nom.yaml","config_bdt_all.yaml")'

#include <TFile.h>
#include <TFileMerger.h>
#include <TH1.h>
#include <TNamed.h>
#include <TSystem.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace {

// Post-merge integral validator: read back a canonical histogram from the
// output and compare to the sum across inputs. Catches the silent partial-
// merge mode where TFileMerger swallows a transient I/O failure on one
// input and writes a single-period output that still RECREATEs cleanly.
// Returns true if integrals match within tolerance (or the validation
// histogram doesn't exist in this file type — e.g. response files).
bool validate_merge_integrals(const std::string &in0, const std::string &in1,
                              const std::string &out, const char *probe_name,
                              double rel_tol = 1e-3)
{
    auto integral = [](const std::string &path, const char *hname) -> double {
        TFile *f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); return -1.0; }
        if (f->TestBit(TFile::kRecovered)) { f->Close(); return -2.0; }
        TH1 *h = (TH1 *) f->Get(hname);
        if (!h) { f->Close(); return -3.0; }
        double i = h->Integral(0, h->GetNbinsX() + 1);  // include under/over
        f->Close();
        return i;
    };
    const double i0 = integral(in0, probe_name);
    const double i1 = integral(in1, probe_name);
    const double io = integral(out, probe_name);
    if (i0 < 0 || i1 < 0) {
        std::cerr << "[VALIDATE] cannot read probe " << probe_name
                  << " from per-period inputs (rc " << i0 << "/" << i1
                  << "); skipping integral check\n";
        return true;  // probe not present => not applicable
    }
    if (io < 0) {
        std::cerr << "[FATAL] merge output unreadable for probe " << probe_name
                  << " in " << out << " (rc=" << io << ")\n";
        return false;
    }
    const double exp = i0 + i1;
    const double rel = std::abs(io - exp) / (exp > 0 ? exp : 1.0);
    if (rel > rel_tol) {
        std::cerr << "[FATAL] merge integral mismatch for " << out << ": "
                  << "merged=" << io << "  expected " << i0 << "+" << i1 << "=" << exp
                  << "  rel=" << rel << "  (>" << rel_tol << ")\n"
                  << "  TFileMerger likely silently dropped one period from a "
                  << "transient I/O failure. Re-run merge_periods.sh.\n";
        return false;
    }
    std::cout << "  [VALIDATE] " << probe_name << ": merged=" << io
              << "  exp=" << exp << "  rel=" << rel << "  OK\n";
    return true;
}

bool merge_one_pair(const std::string &in0, const std::string &in1,
                    const std::string &out, const char *probe_name)
{
    std::cout << "  in0 = " << in0 << "\n  in1 = " << in1 << "\n  out = " << out << std::endl;
    // Pre-merge zombie check on inputs (would otherwise be silently skipped
    // by TFileMerger's default error handling, producing a single-period
    // output that passes "Merge() returned kTRUE").
    for (const std::string &in : {in0, in1}) {
        TFile *f = TFile::Open(in.c_str(), "READ");
        bool bad = (!f || f->IsZombie() || f->TestBit(TFile::kRecovered));
        if (f) f->Close();
        if (bad) {
            std::cerr << "[FATAL] input " << in << " is zombie/recovered. Aborting merge.\n";
            return false;
        }
    }
    TFileMerger m;
    m.OutputFile(out.c_str(), "RECREATE");
    if (!m.AddFile(in0.c_str())) { std::cerr << "[ERROR] failed to add " << in0 << "\n"; return false; }
    if (!m.AddFile(in1.c_str())) { std::cerr << "[ERROR] failed to add " << in1 << "\n"; return false; }
    if (!m.Merge())              { std::cerr << "[ERROR] TFileMerger::Merge failed for " << out << "\n"; return false; }
    // Post-merge integral validation (TFileMerger returns kTRUE even when
    // it silently skipped an input on a transient gpfs I/O failure).
    if (probe_name && !validate_merge_integrals(in0, in1, out, probe_name)) {
        return false;
    }
    return true;
}

std::string make_path(const std::string &base, const std::string &mid,
                      const std::string &var_type)
{
    if (mid.empty()) return base + "_" + var_type + ".root";
    return base + "_" + mid + "_" + var_type + ".root";
}

void write_provenance(const std::string &out, double L0, double L1,
                      double Lt, const std::string &vt0, const std::string &vt1)
{
    TFile f(out.c_str(), "UPDATE");
    if (f.IsZombie()) return;
    TNamed(("merge_lumi_period0_" + vt0).c_str(), Form("%.6f", L0)).Write();
    TNamed(("merge_lumi_period1_" + vt1).c_str(), Form("%.6f", L1)).Write();
    TNamed("merge_lumi_target",                   Form("%.6f", Lt)).Write();
    f.Close();
}

} // namespace

int merge_periods(const std::string &cfg0    = "config_bdt_0rad.yaml",
                  const std::string &cfg1    = "config_bdt_nom.yaml",
                  const std::string &cfg_out = "config_bdt_all.yaml")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node y0 = YAML::LoadFile(cfg0);
    YAML::Node y1 = YAML::LoadFile(cfg1);
    YAML::Node yo = YAML::LoadFile(cfg_out);

    const double L0  = y0["analysis"]["lumi"].as<double>();
    const double L1  = y1["analysis"]["lumi"].as<double>();
    const double Lt0 = y0["analysis"]["lumi_target"].as<double>(L0);
    const double Lt1 = y1["analysis"]["lumi_target"].as<double>(L1);

    if (Lt0 != Lt1) {
        std::cerr << "[WARN] per-period configs disagree on lumi_target: "
                  << "cfg0=" << Lt0 << "  cfg1=" << Lt1
                  << " — per-event scaling will not align cleanly.\n";
    }
    if (std::abs((L0 + L1) - Lt0) > 1e-3) {
        std::cerr << "[WARN] L_0+L_1=" << (L0+L1) << " differs from lumi_target=" << Lt0
                  << ". Per-event lumi_weight assumed sum(L) = lumi_target — reconcile if intended.\n";
    }

    const std::string vt0 = y0["output"]["var_type"].as<std::string>();
    const std::string vt1 = y1["output"]["var_type"].as<std::string>();
    const std::string vto = yo["output"]["var_type"].as<std::string>();

    const std::string eff_b0 = y0["output"]["eff_outfile"].as<std::string>();
    const std::string eff_b1 = y1["output"]["eff_outfile"].as<std::string>();
    const std::string eff_bo = yo["output"]["eff_outfile"].as<std::string>();
    const std::string rsp_b0 = y0["output"]["response_outfile"].as<std::string>();
    const std::string rsp_b1 = y1["output"]["response_outfile"].as<std::string>();
    const std::string rsp_bo = yo["output"]["response_outfile"].as<std::string>();

    std::cout << "=== merge_periods (TFileMerger / hadd-style) ===\n"
              << "  L(" << vt0 << ") = " << L0 << " pb-1   lumi_target = " << Lt0 << "\n"
              << "  L(" << vt1 << ") = " << L1 << " pb-1   lumi_target = " << Lt1 << "\n"
              << "  combined var_type = " << vto << std::endl;

    struct Pair { std::string in0, in1, out; const char *label; const char *probe; };
    // Probe histogram per merge job — chosen because they are NPB-/iso-/cut-
    // independent and straightforward to integrate. Set to nullptr to skip
    // the integral validator (response files don't have a single canonical
    // TH1 to integrate; rely on AddFile/Merge return codes for those).
    std::vector<Pair> jobs = {
        { make_path(eff_b0, "",    vt0), make_path(eff_b1, "",    vt1), make_path(eff_bo, "",    vto), "eff_signal", "h_all_cluster_signal_0" },
        { make_path(eff_b0, "jet", vt0), make_path(eff_b1, "jet", vt1), make_path(eff_bo, "jet", vto), "eff_jet",    "h_all_cluster_signal_0" },
        { make_path(rsp_b0, "",    vt0), make_path(rsp_b1, "",    vt1), make_path(rsp_bo, "",    vto), "response",   nullptr                  },
    };

    int rc = 0;
    for (const auto &j : jobs) {
        std::cout << "\n--- " << j.label << " ---" << std::endl;
        if (!merge_one_pair(j.in0, j.in1, j.out, j.probe)) { rc = 1; continue; }
        write_provenance(j.out, L0, L1, Lt0, vt0, vt1);
    }
    return rc;
}
