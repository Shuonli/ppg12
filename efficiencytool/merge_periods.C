// merge_periods.C
//
// Combine two per-crossing-angle MergeSim outputs (e.g. bdt_0rad + bdt_nom)
// into a single full-run-range MC by lumi-weighted scale + add. Per-period
// MC is left untouched upstream, so DI mix fractions and truth-vertex
// reweights stay correct per period; the combined output reproduces the
// data composition for the union of the two run ranges.
//
// Weights are computed from the per-period configs:
//   w_i = L_i / (L_0 + L_1)
// The combined config's `analysis.lumi` field is NOT used for weighting --
// keeps the script self-consistent if that field is stale.
//
// Handles TH1/TH2/TH3 (incl. TProfile descendants), TEfficiency (rebuilt
// from scaled passed/total), and RooUnfoldResponse (rebuilt from scaled
// measured/truth/response). Other classes are skipped with a warning.
//
// Usage:
//   root -l -b -q 'merge_periods.C("config_bdt_0rad.yaml","config_bdt_nom.yaml","config_bdt_all.yaml")'

#include <RooUnfoldResponse.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TKey.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TSystem.h>
#include <iostream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace {

void merge_one_pair(const std::string &in0_path, const std::string &in1_path,
                    const std::string &out_path, double w0, double w1)
{
    std::cout << "  in0 = " << in0_path << "  (weight " << w0 << ")\n"
              << "  in1 = " << in1_path << "  (weight " << w1 << ")\n"
              << "  out = " << out_path << std::endl;

    TFile *f0 = TFile::Open(in0_path.c_str(), "READ");
    if (!f0 || f0->IsZombie()) {
        std::cerr << "[ERROR] cannot open " << in0_path << std::endl;
        return;
    }
    TFile *f1 = TFile::Open(in1_path.c_str(), "READ");
    if (!f1 || f1->IsZombie()) {
        std::cerr << "[ERROR] cannot open " << in1_path << std::endl;
        f0->Close();
        return;
    }
    TFile *fout = TFile::Open(out_path.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "[ERROR] cannot open output " << out_path << std::endl;
        f0->Close();
        f1->Close();
        return;
    }

    int n_h = 0, n_eff = 0, n_resp = 0, n_skip = 0, n_missing = 0;

    TIter next(f0->GetListOfKeys());
    while (TKey *k0 = (TKey *)next()) {
        const char *name = k0->GetName();
        const std::string cn = k0->GetClassName();

        if (cn == "TH1D" || cn == "TH1F" || cn == "TH2D" || cn == "TH2F" ||
            cn == "TH3D" || cn == "TH3F" || cn == "TProfile" ||
            cn == "TProfile2D")
        {
            TH1 *h0 = (TH1 *)k0->ReadObj();
            h0->SetDirectory(nullptr);
            h0->Scale(w0);
            TH1 *h1 = (TH1 *)f1->Get(name);
            if (h1) h0->Add(h1, w1);
            else { ++n_missing; std::cerr << "[WARN] " << name << " missing in input1\n"; }
            fout->cd();
            h0->Write(name, TObject::kOverwrite);
            delete h0;
            ++n_h;
        }
        else if (cn == "TEfficiency")
        {
            TEfficiency *e0 = (TEfficiency *)k0->ReadObj();
            TEfficiency *e1 = (TEfficiency *)f1->Get(name);
            TH1 *p0 = (TH1 *)e0->GetCopyPassedHisto();
            TH1 *t0 = (TH1 *)e0->GetCopyTotalHisto();
            p0->SetDirectory(nullptr); t0->SetDirectory(nullptr);
            p0->Scale(w0); t0->Scale(w0);
            if (e1) {
                TH1 *p1 = (TH1 *)e1->GetCopyPassedHisto();
                TH1 *t1 = (TH1 *)e1->GetCopyTotalHisto();
                p0->Add(p1, w1); t0->Add(t1, w1);
                delete p1; delete t1;
            } else {
                ++n_missing; std::cerr << "[WARN] " << name << " missing in input1\n";
            }
            // Bins where total is zero would make TEfficiency complain; only
            // happens if both periods are empty, which is fine to skip.
            TEfficiency *eff = new TEfficiency(*p0, *t0);
            eff->SetName(name);
            eff->SetTitle(e0->GetTitle());
            eff->SetUseWeightedEvents();
            eff->SetStatisticOption(TEfficiency::kFNormal);
            fout->cd();
            eff->Write(name, TObject::kOverwrite);
            delete eff;
            delete p0; delete t0;
            delete e0; if (e1) delete e1;
            ++n_eff;
        }
        else if (cn == "RooUnfoldResponse")
        {
            RooUnfoldResponse *r0 = (RooUnfoldResponse *)k0->ReadObj();
            RooUnfoldResponse *r1 = (RooUnfoldResponse *)f1->Get(name);
            // Hresponse/Hmeasured/Htruth return non-owning pointers; clone.
            TH2 *h_resp = (TH2 *)r0->Hresponse()->Clone();
            TH1 *h_meas = (TH1 *)r0->Hmeasured()->Clone();
            TH1 *h_truth = (TH1 *)r0->Htruth()->Clone();
            h_resp->SetDirectory(nullptr);
            h_meas->SetDirectory(nullptr);
            h_truth->SetDirectory(nullptr);
            h_resp->Scale(w0); h_meas->Scale(w0); h_truth->Scale(w0);
            if (r1) {
                h_resp->Add((TH2 *)r1->Hresponse(), w1);
                h_meas->Add(r1->Hmeasured(), w1);
                h_truth->Add(r1->Htruth(), w1);
            } else {
                ++n_missing; std::cerr << "[WARN] " << name << " missing in input1\n";
            }
            RooUnfoldResponse *r = new RooUnfoldResponse(
                h_meas, h_truth, h_resp, name, r0->GetTitle(), false);
            fout->cd();
            r->Write(name, TObject::kOverwrite);
            delete r;
            delete h_resp; delete h_meas; delete h_truth;
            // RooUnfoldResponse from Get/ReadObj is owned by the file in
            // this TStreamer setup; do not delete to avoid double-free.
            ++n_resp;
        }
        else
        {
            std::cerr << "[SKIP] " << cn << " '" << name << "'\n";
            ++n_skip;
        }
    }

    fout->cd();
    TNamed("merge_periods_w0", Form("%.10f", w0)).Write();
    TNamed("merge_periods_w1", Form("%.10f", w1)).Write();
    TNamed("merge_periods_in0", in0_path.c_str()).Write();
    TNamed("merge_periods_in1", in1_path.c_str()).Write();

    fout->Close();
    f0->Close();
    f1->Close();

    std::cout << "  wrote: TH=" << n_h << "  TEfficiency=" << n_eff
              << "  RooUnfoldResponse=" << n_resp
              << "  missing-in-1=" << n_missing
              << "  skipped=" << n_skip << "\n";
}

std::string make_path(const std::string &base, const std::string &mid,
                      const std::string &var_type)
{
    if (mid.empty()) return base + "_" + var_type + ".root";
    return base + "_" + mid + "_" + var_type + ".root";
}

} // namespace

void merge_periods(const std::string &cfg0    = "config_bdt_0rad.yaml",
                   const std::string &cfg1    = "config_bdt_nom.yaml",
                   const std::string &cfg_out = "config_bdt_all.yaml")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node y0 = YAML::LoadFile(cfg0);
    YAML::Node y1 = YAML::LoadFile(cfg1);
    YAML::Node yo = YAML::LoadFile(cfg_out);

    const double L0 = y0["analysis"]["lumi"].as<double>();
    const double L1 = y1["analysis"]["lumi"].as<double>();
    const double Lsum = L0 + L1;
    const double w0 = L0 / Lsum;
    const double w1 = L1 / Lsum;

    const std::string vt0 = y0["output"]["var_type"].as<std::string>();
    const std::string vt1 = y1["output"]["var_type"].as<std::string>();
    const std::string vto = yo["output"]["var_type"].as<std::string>();

    const std::string eff_base_0 = y0["output"]["eff_outfile"].as<std::string>();
    const std::string eff_base_1 = y1["output"]["eff_outfile"].as<std::string>();
    const std::string eff_base_o = yo["output"]["eff_outfile"].as<std::string>();
    const std::string rsp_base_0 = y0["output"]["response_outfile"].as<std::string>();
    const std::string rsp_base_1 = y1["output"]["response_outfile"].as<std::string>();
    const std::string rsp_base_o = yo["output"]["response_outfile"].as<std::string>();

    std::cout << "=== merge_periods ===\n"
              << "  L(" << vt0 << ") = " << L0 << " pb-1\n"
              << "  L(" << vt1 << ") = " << L1 << " pb-1\n"
              << "  L(combined)    = " << Lsum << " pb-1\n"
              << "  w0 = " << w0 << "   w1 = " << w1 << "\n"
              << "  combined var_type = " << vto << std::endl;

    const double Lo_cfg = yo["analysis"]["lumi"].as<double>();
    if (std::abs(Lo_cfg - Lsum) > 1e-3) {
        std::cout << "  [NOTE] " << cfg_out << " lumi=" << Lo_cfg
                  << " differs from sum=" << Lsum
                  << " (script uses sum for MC weights; data side will "
                  << "use the cfg_out lumi -- reconcile if intended)\n";
    }

    struct Pair { std::string in0, in1, out; const char *label; };
    std::vector<Pair> jobs = {
        // signal eff (photon merged)
        { make_path(eff_base_0, "",     vt0),
          make_path(eff_base_1, "",     vt1),
          make_path(eff_base_o, "",     vto), "eff_signal" },
        // inclusive jet eff (jet8/12/20/30/40 merged)
        { make_path(eff_base_0, "jet",  vt0),
          make_path(eff_base_1, "jet",  vt1),
          make_path(eff_base_o, "jet",  vto), "eff_jet"    },
        // response matrix
        { make_path(rsp_base_0, "",     vt0),
          make_path(rsp_base_1, "",     vt1),
          make_path(rsp_base_o, "",     vto), "response"   },
    };

    for (const auto &j : jobs) {
        std::cout << "\n--- " << j.label << " ---" << std::endl;
        merge_one_pair(j.in0, j.in1, j.out, w0, w1);
    }
}
