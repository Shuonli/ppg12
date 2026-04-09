// test_apply_npb_tmva_from_txt.C
//
// Read the shower-shape text file produced by `BDTinput.C` (e.g. shapes_data_npb.txt),
// apply the TMVA NPB score model (same TMVA::Experimental::RBDT usage as `apply_BDT.C`),
// and plot the resulting NPB score distribution(s).
//
// Usage (example):
//   root -l -b -q 'FunWithxgboost/test_apply_npb_tmva_from_txt.C("shapes_data_npb.txt","npb_models/npb_score_tmva.root","npb_score_from_txt.root")'
//
#include <TMVA/RBDT.hxx>
// STL
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
// ROOT
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TSystem.h>

namespace
{
  bool try_parse_float(const std::string& s, float& out)
  {
    char* end = nullptr;
    out = std::strtof(s.c_str(), &end);
    return end != s.c_str() && *end == '\0';
  }
}  // namespace

void test_apply_npb_tmva_from_txt(const std::string& txtfile = "shapes_data_npb.txt",
                                  const std::string& npb_tmva_file = "npb_models/npb_score_tmva.root",
                                  const std::string& outroot = "npb_score_from_txt.root",
                                  const bool apply_training_phase_space = true,
                                  const std::string& model_name = "myBDT",
                                  const int require_pid = 999999,
                                  const bool skip_non_finite_rows = true)
{
  using namespace TMVA::Experimental;

  if (gSystem->AccessPathName(npb_tmva_file.c_str()))
  {
    std::cerr << "ERROR: TMVA model not found: " << npb_tmva_file << std::endl;
    return;
  }

  std::ifstream fin(txtfile);
  if (!fin.is_open())
  {
    std::cerr << "ERROR: cannot open input txt file: " << txtfile << std::endl;
    return;
  }

  // Same training phase-space used in `apply_BDT.C` / `config_npb_training.yaml`
  const float npb_et_min = 6.0;
  const float npb_et_max = 40.0;
  const float npb_eta_max = 0.7;

  // Load TMVA model (same pattern as `apply_BDT.C`)
  // IMPORTANT: the first argument must match the model name stored inside the TMVA file.
  // `apply_BDT.C` uses "myBDT", so that is the default here too.
  RBDT npb_bdt(model_name, npb_tmva_file);

  // Histograms
  TH1D* h_all = new TH1D("h_npb_score_all", "NPB score (from txt);NPB score;Clusters", 100, 0.0, 1.0);
  h_all->Sumw2();

  // Optional: split by pid if present (MC photonclass) or pid=-1 (data NPB)
  std::map<int, TH1D*> h_by_pid;

  std::string line;
  long long n_lines = 0;
  long long n_parsed = 0;
  long long n_used = 0;
  long long n_pid_filtered = 0;
  long long n_non_finite = 0;

  while (std::getline(fin, line))
  {
    ++n_lines;
    if (line.empty()) continue;

    // Skip header lines that start with non-numeric token (BDTinput writes a header)
    {
      std::istringstream iss(line);
      std::string tok0;
      if (!(iss >> tok0)) continue;
      float tmp = 0;
      if (!try_parse_float(tok0, tmp)) continue;
    }

    // The BDTinput.C header (column order) is:
    // cluster_Et cluster_Eta cluster_Phi vertexz
    // e11_over_e33 e32_over_e35 e11_over_e22 e11_over_e13
    // e11_over_e15 e11_over_e17 e11_over_e31
    // e11_over_e51 e11_over_e71 e22_over_e33
    // e22_over_e35 e22_over_e37 e22_over_e53
    // cluster_prob cluster_weta_cogx cluster_wphi_cogx
    // cluster_et1 cluster_et2 cluster_et3 cluster_et4
    // cluster_w32 cluster_w52 cluster_w72
    // recoisoET is_tight pid
    float cluster_Et = 0, cluster_Eta = 0, cluster_Phi = 0, vertexz = 0;
    float e11_over_e33 = 0, e32_over_e35 = 0, e11_over_e22 = 0, e11_over_e13 = 0;
    float e11_over_e15 = 0, e11_over_e17 = 0, e11_over_e31 = 0, e11_over_e51 = 0, e11_over_e71 = 0;
    float e22_over_e33 = 0, e22_over_e35 = 0, e22_over_e37 = 0, e22_over_e53 = 0;
    float cluster_prob = 0, cluster_weta_cogx = 0, cluster_wphi_cogx = 0;
    float cluster_et1 = 0, cluster_et2 = 0, cluster_et3 = 0, cluster_et4 = 0;
    float cluster_w32 = 0, cluster_w52 = 0, cluster_w72 = 0;
    float recoisoET = 0;
    int is_tight = 0;
    int pid = 0;

    std::istringstream iss(line);
    if (!(iss >> cluster_Et >> cluster_Eta >> cluster_Phi >> vertexz
              >> e11_over_e33 >> e32_over_e35 >> e11_over_e22 >> e11_over_e13
              >> e11_over_e15 >> e11_over_e17 >> e11_over_e31
              >> e11_over_e51 >> e11_over_e71 >> e22_over_e33
              >> e22_over_e35 >> e22_over_e37 >> e22_over_e53
              >> cluster_prob >> cluster_weta_cogx >> cluster_wphi_cogx
              >> cluster_et1 >> cluster_et2 >> cluster_et3 >> cluster_et4
              >> cluster_w32 >> cluster_w52 >> cluster_w72
              >> recoisoET >> is_tight >> pid))
    {
      continue;
    }
    ++n_parsed;

    if (require_pid != 999999 && pid != require_pid)
    {
      ++n_pid_filtered;
      continue;
    }

    if (apply_training_phase_space)
    {
      if (cluster_Et < npb_et_min || cluster_Et > npb_et_max) continue;
      if (std::fabs(cluster_Eta) > npb_eta_max) continue;
    }

    // Feature order MUST match `apply_BDT.C` x_npb vector (and config_npb_training.yaml feature_list):
    // [Et, Eta, vertexz, e11/e33, e32/e35, e11/e22, e11/e13, e11/e15, e11/e17, e11/e31,
    //  e11/e51, e11/e71, e22/e33, e22/e35, e22/e37, e22/e53, weta_cogx, wphi_cogx,
    //  et1, et2, et3, et4, w32, w52, w72]
    std::vector<float> x_npb = {
      cluster_Et,
      cluster_Eta,
      vertexz,
      e11_over_e33,
      e32_over_e35,
      e11_over_e22,
      e11_over_e13,
      e11_over_e15,
      e11_over_e17,
      e11_over_e31,
      e11_over_e51,
      e11_over_e71,
      e22_over_e33,
      e22_over_e35,
      e22_over_e37,
      e22_over_e53,
      cluster_weta_cogx,
      cluster_wphi_cogx,
      cluster_et1,
      cluster_et2,
      cluster_et3,
      cluster_et4,
      cluster_w32,
      cluster_w52,
      cluster_w72
    };

    bool ok = true;
    for (const float v : x_npb)
    {
      if (!std::isfinite(v))
      {
        ok = false;
        break;
      }
    }
    if (!ok)
    {
      ++n_non_finite;
      if (skip_non_finite_rows)
      {
        continue;
      }
      // If not skipping, coerce non-finite values to 0 (closer to apply_BDT.C's divide-by-zero protections)
      for (float& v : x_npb)
      {
        if (!std::isfinite(v)) v = 0.0f;
      }
    }

    const float score = npb_bdt.Compute(x_npb)[0];
    ++n_used;

    h_all->Fill(score);

    if (!h_by_pid.count(pid))
    {
      h_by_pid[pid] = new TH1D(Form("h_npb_score_pid%d", pid),
                               Form("NPB score pid=%d;NPB score;Clusters", pid),
                               100, 0.0, 1.0);
      h_by_pid[pid]->Sumw2();
    }
    h_by_pid[pid]->Fill(score);
  }

  std::cout << "Read lines:   " << n_lines << "\n"
            << "Parsed rows:  " << n_parsed << "\n"
            << "PID filter:   " << require_pid << " (999999 means no filter)\n"
            << "PID skipped:  " << n_pid_filtered << "\n"
            << "Non-finite:   " << n_non_finite
            << (skip_non_finite_rows ? " (skipped)" : " (coerced to 0)") << "\n"
            << "Used rows:    " << n_used << (apply_training_phase_space ? " (phase-space cuts ON)" : " (phase-space cuts OFF)") << "\n"
            << "Model:        " << npb_tmva_file << "\n"
            << "Input:        " << txtfile << std::endl;

  // Output ROOT file + plot
  TFile* fout = TFile::Open(outroot.c_str(), "RECREATE");
  if (!fout || fout->IsZombie())
  {
    std::cerr << "ERROR: cannot create output ROOT file: " << outroot << std::endl;
    return;
  }
  fout->cd();
  h_all->Write();
  for (auto& kv : h_by_pid) kv.second->Write();

  TCanvas* c = new TCanvas("c_npb_score", "c_npb_score", 900, 700);
  c->SetTicks(1, 1);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);

  h_all->SetLineColor(kBlack);
  h_all->SetLineWidth(2);
  h_all->Draw("hist");

  // Overlay per-pid histos if there are only a few; otherwise keep it simple.
  if (h_by_pid.size() > 1 && h_by_pid.size() <= 6)
  {
    TLegend* leg = new TLegend(0.62, 0.62, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_all, "all", "l");
    int color = 2;
    for (auto& kv : h_by_pid)
    {
      if (kv.first == 0 && h_by_pid.size() > 1) { /* keep */ }
      kv.second->SetLineColor(color++);
      kv.second->SetLineWidth(2);
      kv.second->SetLineStyle(2);
      kv.second->Draw("hist same");
      leg->AddEntry(kv.second, Form("pid=%d", kv.first), "l");
      if (color == 5) color = 6; // skip yellow
    }
    leg->Draw();
  }

  c->Write();
  c->SaveAs("npb_score_from_txt.pdf");

  fout->Write();
  fout->Close();

  std::cout << "Wrote: " << outroot << " and npb_score_from_txt.pdf" << std::endl;
}


