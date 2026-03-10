// test_vertex_reweight.C
//
// Unit tests for the vertex-z reweighting logic in RecoEffCalculator_TTreeReader.C.
// All inputs are synthetic — no real data files are required.
//
// Behaviour under test (exact code from the main macro):
//   weight = cross_weight;
//   vertex_weight = 1.0;
//   if (vertex_reweight_on) {
//       int bin = h_vtx_rw->FindBin(vertexz);
//       if (bin < 1)                      bin = 1;
//       if (bin > h_vtx_rw->GetNbinsX()) bin = h_vtx_rw->GetNbinsX();
//       vertex_weight = h_vtx_rw->GetBinContent(bin);
//   }
//   if (!std::isfinite(vertex_weight) || vertex_weight <= 0.0)
//       vertex_weight = 1.0;
//   weight *= vertex_weight;
//
// Tests covered:
//   Group 1 – reweight OFF  : weight always equals cross_weight
//   Group 2 – normal in-range: weight = cross_weight * bin_content
//   Group 3 – out-of-range   : vz << lo → clamp to bin 1; vz >> hi → clamp to last bin
//   Group 4 – guard cases    : bin content = 0 or < 0 → vertex_weight forced to 1.0
//   Group 5 – loading errors : missing file, wrong histogram name
//
// Usage:
//   root -l -b -q test_vertex_reweight.C

#include <TFile.h>
#include <TH1D.h>
#include <cmath>
#include <iostream>
#include <string>

// ── assertion helpers ─────────────────────────────────────────────────────
static int n_pass = 0, n_fail = 0;

void check(bool cond, const char *label)
{
    const char *tag = cond ? "  PASS" : "  FAIL";
    std::cout << tag << "  " << label << "\n";
    cond ? ++n_pass : ++n_fail;
}

void checkNear(double got, double expected, double tol, const char *label)
{
    bool ok = std::fabs(got - expected) < tol;
    std::cout << (ok ? "  PASS" : "  FAIL")
              << "  " << label
              << "  (got " << got << ", expected " << expected << ")\n";
    ok ? ++n_pass : ++n_fail;
}

// ── Exact replication of per-event vertex weight logic ────────────────────
// Lines 1204-1229 of RecoEffCalculator_TTreeReader.C
double computeWeight(TH1 *h_vtx_rw, float vertexz,
                     bool vertex_reweight_on, double cross_weight)
{
    double weight        = cross_weight;
    double vertex_weight = 1.0;

    if (vertex_reweight_on)
    {
        int bin = h_vtx_rw->FindBin(vertexz);
        if (bin < 1)                      bin = 1;
        if (bin > h_vtx_rw->GetNbinsX()) bin = h_vtx_rw->GetNbinsX();
        vertex_weight = h_vtx_rw->GetBinContent(bin);
    }

    if (!std::isfinite(vertex_weight) || vertex_weight <= 0.0)
        vertex_weight = 1.0;

    weight *= vertex_weight;
    return weight;
}

// ─────────────────────────────────────────────────────────────────────────
void test_vertex_reweight()
{
    std::cout << "\n=== test_vertex_reweight ===\n\n";

    // ── Build synthetic vertex reweight histogram ─────────────────────────
    // 20 bins, -10 to +10 cm (1 cm / bin).
    // Bin i content = 0.5 + 0.1*i  →  0.6, 0.7, …, 2.5
    // Bin layout:
    //   bin 1 : [-10,-9], centre -9.5, content 0.6
    //   bin 2 : [ -9,-8], centre -8.5, content 0.7
    //   bin 3 : [ -8,-7], centre -7.5, content → overridden to 0.0  (guard test)
    //   bin 4 : [ -7,-6], centre -6.5, content → overridden to -0.5 (guard test)
    //   bin 5 : [ -6,-5], centre -5.5, content 1.0
    //   …
    //   bin 11: [  0, 1], centre  0.5, content 1.6
    //   …
    //   bin 20: [  9,10], centre  9.5, content 2.5

    const char *vtxfile  = "/tmp/test_vtxrw_tmp.root";
    const char *histname = "h_vertexz_ratio_data_over_mccombined";

    TH1D *h_ref = new TH1D(histname, histname, 20, -10.0, 10.0);
    for (int i = 1; i <= 20; i++)
        h_ref->SetBinContent(i, 0.5 + 0.1 * i);

    h_ref->SetBinContent(3,  0.0);   // zero     → vertex_weight must be guarded to 1.0
    h_ref->SetBinContent(4, -0.5);   // negative → vertex_weight must be guarded to 1.0

    {
        TFile *fw = TFile::Open(vtxfile, "RECREATE");
        h_ref->Write();
        fw->Close();
    }

    // ── Load histogram (mirrors loading code in RecoEffCalculator) ────────
    TH1 *h_vtx_rw = nullptr;
    {
        TFile *fvtx = TFile::Open(vtxfile, "READ");
        check(fvtx && !fvtx->IsZombie(),
              "Synthetic vertex reweight file opens successfully");

        TH1 *htmp = dynamic_cast<TH1*>(fvtx->Get(histname));
        check(htmp != nullptr,
              "h_vertexz_ratio_data_over_mccombined found in file");

        if (htmp)
        {
            h_vtx_rw = dynamic_cast<TH1*>(htmp->Clone("h_vtxrw_clone"));
            h_vtx_rw->SetDirectory(nullptr);
            h_vtx_rw->Sumw2();   // matches Sumw2() call in main code
        }
        fvtx->Close();
    }
    if (!h_vtx_rw)
    {
        std::cerr << "Cannot continue: histogram not loaded.\n";
        return;
    }

    // Representative cross-section weight (photon5cross / photon20cross)
    const double cross = 1122.7;

    // =====================================================================
    // Group 1: vertex_reweight_on = 0  →  weight == cross_weight always
    // =====================================================================
    std::cout << "\n--- Group 1: vertex_reweight_on = 0 ---\n";

    for (float vz : {-15.f, -9.5f, 0.5f, 7.5f, 999.f})
    {
        double w = computeWeight(h_vtx_rw, vz, /*on=*/false, cross);
        checkNear(w, cross, 1e-9,
                  Form("OFF, vz=%.1f  →  weight == cross_weight", (double)vz));
    }

    // =====================================================================
    // Group 2: normal in-range vertexz  →  weight = cross * bin_content
    // =====================================================================
    std::cout << "\n--- Group 2: vertex_reweight_on = 1, in-range vz ---\n";

    // Use bin centres of non-special bins (skip bins 3 & 4)
    struct { float vz; int expected_bin; } cases[] = {
        {-9.5f,  1},   // bin 1, content = 0.6
        {-8.5f,  2},   // bin 2, content = 0.7
        {-5.5f,  5},   // bin 5, content = 1.0
        { 0.5f, 11},   // bin 11, content = 1.6
        { 7.5f, 18},   // bin 18, content = 2.3
        { 9.5f, 20},   // bin 20, content = 2.5
    };

    for (auto &c : cases)
    {
        int   bin    = h_vtx_rw->FindBin(c.vz);
        double vw    = h_vtx_rw->GetBinContent(bin);
        double w     = computeWeight(h_vtx_rw, c.vz, /*on=*/true, cross);
        check(bin == c.expected_bin,
              Form("FindBin(%.1f) == %d", (double)c.vz, c.expected_bin));
        checkNear(w, cross * vw, 1e-6,
                  Form("ON, vz=%.1f → bin %d, weight = cross * %.3f",
                       (double)c.vz, bin, vw));
    }

    // =====================================================================
    // Group 3: out-of-range clamping
    // =====================================================================
    std::cout << "\n--- Group 3: out-of-range vz → clamped to boundary bin ---\n";

    // Below range → bin 1
    {
        float vz = -999.0f;
        double expected_vw = h_vtx_rw->GetBinContent(1);
        double w = computeWeight(h_vtx_rw, vz, /*on=*/true, cross);
        checkNear(w, cross * expected_vw, 1e-6,
                  Form("ON, vz=-999 (below lo) → clamped to bin 1, content=%.3f",
                       expected_vw));
    }

    // Above range → last bin
    {
        float vz   = +999.0f;
        int   nb   = h_vtx_rw->GetNbinsX();
        double expected_vw = h_vtx_rw->GetBinContent(nb);
        double w = computeWeight(h_vtx_rw, vz, /*on=*/true, cross);
        checkNear(w, cross * expected_vw, 1e-6,
                  Form("ON, vz=+999 (above hi) → clamped to bin %d, content=%.3f",
                       nb, expected_vw));
    }

    // Underflow bin (FindBin returns 0) → clamped to 1
    {
        float vz = -10.5f;   // falls in ROOT underflow bin (bin 0)
        int raw_bin = h_vtx_rw->FindBin(vz);
        check(raw_bin == 0,
              Form("Pre-condition: FindBin(%.1f) == 0 (underflow)", (double)vz));
        double expected_vw = h_vtx_rw->GetBinContent(1);
        double w = computeWeight(h_vtx_rw, vz, /*on=*/true, cross);
        checkNear(w, cross * expected_vw, 1e-6,
                  "ON, underflow bin → clamped to bin 1");
    }

    // =====================================================================
    // Group 4: bad bin content → vertex_weight guarded to 1.0
    // =====================================================================
    std::cout << "\n--- Group 4: bad bin content → guarded to 1.0 ---\n";

    // Bin 3 content = 0.0
    {
        float vz  = -7.5f;   // centre of bin 3
        int   bin = h_vtx_rw->FindBin(vz);
        check(bin == 3,
              Form("Pre-condition: FindBin(%.1f) == 3", (double)vz));
        check(h_vtx_rw->GetBinContent(bin) == 0.0,
              "Pre-condition: bin 3 content is 0");
        double w = computeWeight(h_vtx_rw, vz, /*on=*/true, cross);
        checkNear(w, cross * 1.0, 1e-6,
                  "ON, bin content=0 → vertex_weight guarded to 1.0");
    }

    // Bin 4 content = -0.5
    {
        float vz  = -6.5f;   // centre of bin 4
        int   bin = h_vtx_rw->FindBin(vz);
        check(bin == 4,
              Form("Pre-condition: FindBin(%.1f) == 4", (double)vz));
        check(h_vtx_rw->GetBinContent(bin) < 0.0,
              Form("Pre-condition: bin 4 content is negative (%.2f)",
                   h_vtx_rw->GetBinContent(bin)));
        double w = computeWeight(h_vtx_rw, vz, /*on=*/true, cross);
        checkNear(w, cross * 1.0, 1e-6,
                  "ON, bin content<0 → vertex_weight guarded to 1.0");
    }

    // =====================================================================
    // Group 5: file / histogram loading errors
    // =====================================================================
    std::cout << "\n--- Group 5: loading errors ---\n";

    // Non-existent file
    {
        int saved = gErrorIgnoreLevel;
        gErrorIgnoreLevel = kFatal;   // suppress ROOT's error printout
        TFile *fbad = TFile::Open("/tmp/nonexistent_vtxrw_xyzabc.root", "READ");
        gErrorIgnoreLevel = saved;
        bool is_bad = (!fbad || fbad->IsZombie());
        check(is_bad,
              "Non-existent file → TFile::Open returns null/zombie "
              "(code returns early)");
        if (fbad) { fbad->Close(); delete fbad; }
    }

    // Wrong histogram name
    {
        TFile *fvtx = TFile::Open(vtxfile, "READ");
        TH1 *htmp = dynamic_cast<TH1*>(fvtx->Get("h_wrong_histogram_name"));
        check(htmp == nullptr,
              "Wrong histogram name → Get() returns nullptr "
              "(code returns early before null deref)");
        fvtx->Close();
    }

    // Correct name, correct file → should not be nullptr
    {
        TFile *fvtx = TFile::Open(vtxfile, "READ");
        TH1 *htmp = dynamic_cast<TH1*>(fvtx->Get(histname));
        check(htmp != nullptr,
              "Correct name → Get() returns valid pointer");
        fvtx->Close();
    }

    // =====================================================================
    // Summary
    // =====================================================================
    std::cout << "\n─────────────────────────────────────────\n";
    std::cout << "Results: " << n_pass << " passed,  " << n_fail << " failed\n";
    if (n_fail == 0)
        std::cout << "All tests PASSED.\n";
    else
        std::cout << "*** " << n_fail << " test(s) FAILED ***\n";
    std::cout << "─────────────────────────────────────────\n\n";

    // Remove temp file
    gSystem->Unlink(vtxfile);
    delete h_ref;
}
