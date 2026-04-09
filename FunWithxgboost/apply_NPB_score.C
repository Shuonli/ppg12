// apply_NPB_score.C
//
// Apply the NPB classification score (trained in `train_npb_score.py`) to a ROOT TTree,
// in the same spirit as `apply_BDT.C`, using TMVA's fast `RBDT` inference.
//
// Usage (example):
//   root -l -b -q 'FunWithxgboost/apply_NPB_score.C("FunWithxgboost/config_nom.yaml","FunWithxgboost/npb_models/npb_score_tmva.root","FunWithxgboost/npb_models/npb_score_metadata.yaml","photon10","")'
//
#include <TMVA/RBDT.hxx> // fast interpreter
#include <yaml-cpp/yaml.h>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{
  constexpr int nclustercontainermx = 4096;

  inline float safe_div(float num, float den, float default_val = 0.0f)
  {
    return (std::isfinite(den) && std::fabs(den) > 0.0f) ? (num / den) : default_val;
  }

  inline bool hasBranch(TTree *t, const std::string &bname)
  {
    return (t && t->GetBranch(bname.c_str()) != nullptr);
  }

  inline void enableAndBind(TTree *t, const std::string &bname, void *addr, bool required)
  {
    if (!hasBranch(t, bname))
    {
      if (required)
      {
        std::cerr << "ERROR: required branch missing: " << bname << std::endl;
      }
      return;
    }
    t->SetBranchStatus(bname.c_str(), 1);
    t->SetBranchAddress(bname.c_str(), addr);
  }
} // namespace

void apply_NPB_score(const std::string &analysis_configname = "config_nom.yaml",
                     const std::string &tmva_model_file = "npb_models/npb_score_tmva.root",
                     const std::string &npb_metadata_yaml = "npb_models/npb_score_metadata.yaml",
                     const std::string &filetype = "photon10",
                     const std::string &inputfilename = "")
{
  using namespace TMVA::Experimental;

  // Keep the same yaml-cpp load convention as `apply_BDT.C`
  gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");

  YAML::Node analysisCfg = YAML::LoadFile(analysis_configname);
  YAML::Node npbMeta = YAML::LoadFile(npb_metadata_yaml);

  const bool issim = (filetype != "data");

  const std::string infilename_root_dir = analysisCfg["input"]["photon_jet_file_root_dir"].as<std::string>();
  const std::string infilename_branch_dir = analysisCfg["input"]["photon_jet_file_branch_dir"].as<std::string>();
  std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;
  if (!issim)
  {
    if (inputfilename.empty())
    {
      std::cerr << "ERROR: filetype=data requires a non-empty inputfilename" << std::endl;
      return;
    }
    infilename = inputfilename;
  }

  const std::string treename = analysisCfg["input"]["tree"].as<std::string>();
  const std::string clusternodename = analysisCfg["input"]["cluster_node_name"].as<std::string>();

  // recoisoET computation settings (match `BDTinput.C` logic)
  const int conesize = analysisCfg["analysis"]["cone_size"].as<int>(3);
  const int iso_threshold = analysisCfg["analysis"]["iso_threshold"].as<int>(0);

  // Apply the same kinematic range as training by default
  float et_min = 6.0;
  float et_max = 1.0e9;
  float eta_max = 1.0e9;
  if (npbMeta["config"] && npbMeta["config"]["data"])
  {
    et_min = npbMeta["config"]["data"]["et_min"].as<float>(et_min);
    et_max = npbMeta["config"]["data"]["et_max"].as<float>(et_max);
    eta_max = npbMeta["config"]["data"]["eta_max"].as<float>(eta_max);
  }

  // Load feature list (order matters!)
  std::vector<std::string> feature_list;
  if (npbMeta["features"])
  {
    feature_list = npbMeta["features"].as<std::vector<std::string>>();
  }
  else if (npbMeta["config"] && npbMeta["config"]["features"] && npbMeta["config"]["features"]["feature_list"])
  {
    feature_list = npbMeta["config"]["features"]["feature_list"].as<std::vector<std::string>>();
  }
  else
  {
    std::cerr << "ERROR: could not find feature list in " << npb_metadata_yaml << std::endl;
    return;
  }

  // TMVA method name used when exporting the model (default matches training code in this repo)
  const std::string method_name = "myBDT";
  TMVA::Experimental::RBDT bdt(method_name, tmva_model_file);

  // Open input
  TFile *fin = new TFile(infilename.c_str(), "READ");
  if (!fin || fin->IsZombie())
  {
    std::cerr << "ERROR: could not open input file: " << infilename << std::endl;
    return;
  }

  TTree *slimtree = (TTree *)fin->Get(treename.c_str());
  if (!slimtree)
  {
    std::cerr << "ERROR: could not find tree '" << treename << "' in " << infilename << std::endl;
    return;
  }

  // Output name (same convention as `apply_BDT.C`)
  std::string outfile_name = filetype + "/bdt_npb_score.root";
  if (!issim)
  {
    std::string namebase = inputfilename.substr(0, inputfilename.find_last_of("."));
    outfile_name = namebase + "_with_npb_score.root";
  }
  TFile *fout = new TFile(outfile_name.c_str(), "RECREATE");
  TTree *outtree = slimtree->CloneTree(0);

  // Required branches
  int ncluster = 0;
  float vertexz = 0.0;

  float cluster_Et[nclustercontainermx];
  float cluster_Eta[nclustercontainermx];
  float cluster_Phi[nclustercontainermx];
  float cluster_prob[nclustercontainermx];
  float cluster_weta_cogx[nclustercontainermx];
  float cluster_wphi_cogx[nclustercontainermx];
  float cluster_et1[nclustercontainermx];
  float cluster_et2[nclustercontainermx];
  float cluster_et3[nclustercontainermx];
  float cluster_et4[nclustercontainermx];

  float cluster_iso_02[nclustercontainermx];
  float cluster_iso_03[nclustercontainermx];
  float cluster_iso_04[nclustercontainermx];
  float cluster_iso_03_60_emcal[nclustercontainermx];
  float cluster_iso_03_60_hcalin[nclustercontainermx];
  float cluster_iso_03_60_hcalout[nclustercontainermx];

  float cluster_e11[nclustercontainermx];
  float cluster_e22[nclustercontainermx];
  float cluster_e13[nclustercontainermx];
  float cluster_e15[nclustercontainermx];
  float cluster_e17[nclustercontainermx];
  float cluster_e31[nclustercontainermx];
  float cluster_e51[nclustercontainermx];
  float cluster_e71[nclustercontainermx];
  float cluster_e33[nclustercontainermx];
  float cluster_e35[nclustercontainermx];
  float cluster_e37[nclustercontainermx];
  float cluster_e53[nclustercontainermx];
  float cluster_e32[nclustercontainermx];

  float cluster_w32[nclustercontainermx];
  float cluster_w52[nclustercontainermx];
  float cluster_w72[nclustercontainermx];

  slimtree->SetBranchStatus("*", 1);
  enableAndBind(slimtree, Form("ncluster_%s", clusternodename.c_str()), &ncluster, true);
  enableAndBind(slimtree, "vertexz", &vertexz, true);

  enableAndBind(slimtree, Form("cluster_Et_%s", clusternodename.c_str()), &cluster_Et, true);
  enableAndBind(slimtree, Form("cluster_Eta_%s", clusternodename.c_str()), &cluster_Eta, true);
  enableAndBind(slimtree, Form("cluster_Phi_%s", clusternodename.c_str()), &cluster_Phi, true);
  enableAndBind(slimtree, Form("cluster_prob_%s", clusternodename.c_str()), &cluster_prob, true);
  enableAndBind(slimtree, Form("cluster_weta_cogx_%s", clusternodename.c_str()), &cluster_weta_cogx, true);
  enableAndBind(slimtree, Form("cluster_wphi_cogx_%s", clusternodename.c_str()), &cluster_wphi_cogx, true);
  enableAndBind(slimtree, Form("cluster_et1_%s", clusternodename.c_str()), &cluster_et1, true);
  enableAndBind(slimtree, Form("cluster_et2_%s", clusternodename.c_str()), &cluster_et2, true);
  enableAndBind(slimtree, Form("cluster_et3_%s", clusternodename.c_str()), &cluster_et3, true);
  enableAndBind(slimtree, Form("cluster_et4_%s", clusternodename.c_str()), &cluster_et4, true);

  enableAndBind(slimtree, Form("cluster_iso_02_%s", clusternodename.c_str()), &cluster_iso_02, true);
  enableAndBind(slimtree, Form("cluster_iso_03_%s", clusternodename.c_str()), &cluster_iso_03, true);
  enableAndBind(slimtree, Form("cluster_iso_04_%s", clusternodename.c_str()), &cluster_iso_04, true);
  enableAndBind(slimtree, Form("cluster_iso_03_60_emcal_%s", clusternodename.c_str()), &cluster_iso_03_60_emcal, true);
  enableAndBind(slimtree, Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()), &cluster_iso_03_60_hcalin, true);
  enableAndBind(slimtree, Form("cluster_iso_03_60_hcalout_%s", clusternodename.c_str()), &cluster_iso_03_60_hcalout, true);

  enableAndBind(slimtree, Form("cluster_e11_%s", clusternodename.c_str()), &cluster_e11, true);
  enableAndBind(slimtree, Form("cluster_e22_%s", clusternodename.c_str()), &cluster_e22, true);
  enableAndBind(slimtree, Form("cluster_e13_%s", clusternodename.c_str()), &cluster_e13, true);
  enableAndBind(slimtree, Form("cluster_e15_%s", clusternodename.c_str()), &cluster_e15, true);
  enableAndBind(slimtree, Form("cluster_e17_%s", clusternodename.c_str()), &cluster_e17, true);
  enableAndBind(slimtree, Form("cluster_e31_%s", clusternodename.c_str()), &cluster_e31, true);
  enableAndBind(slimtree, Form("cluster_e51_%s", clusternodename.c_str()), &cluster_e51, true);
  enableAndBind(slimtree, Form("cluster_e71_%s", clusternodename.c_str()), &cluster_e71, true);
  enableAndBind(slimtree, Form("cluster_e33_%s", clusternodename.c_str()), &cluster_e33, true);
  enableAndBind(slimtree, Form("cluster_e35_%s", clusternodename.c_str()), &cluster_e35, true);
  enableAndBind(slimtree, Form("cluster_e37_%s", clusternodename.c_str()), &cluster_e37, true);
  enableAndBind(slimtree, Form("cluster_e53_%s", clusternodename.c_str()), &cluster_e53, true);
  enableAndBind(slimtree, Form("cluster_e32_%s", clusternodename.c_str()), &cluster_e32, true);

  enableAndBind(slimtree, Form("cluster_w32_%s", clusternodename.c_str()), &cluster_w32, true);
  enableAndBind(slimtree, Form("cluster_w52_%s", clusternodename.c_str()), &cluster_w52, true);
  enableAndBind(slimtree, Form("cluster_w72_%s", clusternodename.c_str()), &cluster_w72, true);

  // Output branch
  float cluster_npb_score[nclustercontainermx];
  std::string bname = Form("cluster_npb_score_%s", clusternodename.c_str());
  std::string leaf = Form("%s[ncluster_%s]/F", bname.c_str(), clusternodename.c_str());
  outtree->Branch(bname.c_str(), cluster_npb_score, leaf.c_str());

  const int nentries = slimtree->GetEntries();
  std::cout << "Applying NPB score to " << nentries << " entries, writing: " << outfile_name << std::endl;
  std::cout << "Model: " << tmva_model_file << " (method=" << method_name << ")" << std::endl;

  for (int ientry = 0; ientry < nentries; ++ientry)
  {
    if (ientry % 10000 == 0)
      std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;

    slimtree->GetEntry(ientry);

    const int ncl = std::min(ncluster, nclustercontainermx);
    for (int icluster = 0; icluster < ncl; ++icluster)
    {
      // Default: invalid score (outside training phase space)
      cluster_npb_score[icluster] = -1.0f;

      const float et = cluster_Et[icluster];
      const float eta = cluster_Eta[icluster];
      if (!(et >= et_min && et <= et_max && std::fabs(eta) <= eta_max))
      {
        continue;
      }

      // recoisoET (match `BDTinput.C`)
      float recoisoET = -999.0f;
      if (conesize == 4)
        recoisoET = cluster_iso_04[icluster];
      else if (conesize == 3)
        recoisoET = cluster_iso_03[icluster];
      else if (conesize == 2)
        recoisoET = cluster_iso_02[icluster];

      if (iso_threshold)
      {
        recoisoET = cluster_iso_03_60_emcal[icluster] + cluster_iso_03_60_hcalin[icluster] +
                    cluster_iso_03_60_hcalout[icluster];
      }

      // Build feature vector in the exact order used during training/export.
      std::vector<float> x;
      x.reserve(feature_list.size());
      for (const auto &fname : feature_list)
      {
        if (fname == "cluster_Et")
          x.push_back(cluster_Et[icluster]);
        else if (fname == "cluster_Eta")
          x.push_back(cluster_Eta[icluster]);
        else if (fname == "cluster_Phi")
          x.push_back(cluster_Phi[icluster]);
        else if (fname == "vertexz")
          x.push_back(vertexz);
        else if (fname == "e11_over_e33")
          x.push_back(safe_div(cluster_e11[icluster], cluster_e33[icluster]));
        else if (fname == "e32_over_e35")
          x.push_back(safe_div(cluster_e32[icluster], cluster_e35[icluster]));
        else if (fname == "e11_over_e22")
          x.push_back(safe_div(cluster_e11[icluster], cluster_e22[icluster]));
        else if (fname == "e11_over_e13")
          x.push_back(safe_div(cluster_e11[icluster], cluster_e13[icluster]));
        else if (fname == "e11_over_e15")
          x.push_back(safe_div(cluster_e11[icluster], cluster_e15[icluster]));
        else if (fname == "e11_over_e17")
          x.push_back(safe_div(cluster_e11[icluster], cluster_e17[icluster]));
        else if (fname == "e11_over_e31")
          x.push_back(safe_div(cluster_e11[icluster], cluster_e31[icluster]));
        else if (fname == "e11_over_e51")
          x.push_back(safe_div(cluster_e11[icluster], cluster_e51[icluster]));
        else if (fname == "e11_over_e71")
          x.push_back(safe_div(cluster_e11[icluster], cluster_e71[icluster]));
        else if (fname == "e22_over_e33")
          x.push_back(safe_div(cluster_e22[icluster], cluster_e33[icluster]));
        else if (fname == "e22_over_e35")
          x.push_back(safe_div(cluster_e22[icluster], cluster_e35[icluster]));
        else if (fname == "e22_over_e37")
          x.push_back(safe_div(cluster_e22[icluster], cluster_e37[icluster]));
        else if (fname == "e22_over_e53")
          x.push_back(safe_div(cluster_e22[icluster], cluster_e53[icluster]));
        else if (fname == "cluster_prob")
          x.push_back(cluster_prob[icluster]);
        else if (fname == "cluster_weta_cogx")
          x.push_back(cluster_weta_cogx[icluster]);
        else if (fname == "cluster_wphi_cogx")
          x.push_back(cluster_wphi_cogx[icluster]);
        else if (fname == "cluster_et1")
          x.push_back(cluster_et1[icluster]);
        else if (fname == "cluster_et2")
          x.push_back(cluster_et2[icluster]);
        else if (fname == "cluster_et3")
          x.push_back(cluster_et3[icluster]);
        else if (fname == "cluster_et4")
          x.push_back(cluster_et4[icluster]);
        else if (fname == "cluster_w32")
          x.push_back(cluster_w32[icluster]);
        else if (fname == "cluster_w52")
          x.push_back(cluster_w52[icluster]);
        else if (fname == "cluster_w72")
          x.push_back(cluster_w72[icluster]);
        else if (fname == "recoisoET")
          x.push_back(recoisoET);
        else
        {
          // Unknown feature name: push a placeholder to keep alignment, but warn once per job.
          static bool warned = false;
          if (!warned)
          {
            std::cerr << "WARNING: unknown feature name '" << fname
                      << "' encountered in feature_list; filling 0 for this and future occurrences." << std::endl;
            warned = true;
          }
          x.push_back(0.0f);
        }
      }

      // Compute score
      const auto out = bdt.Compute(x);
      cluster_npb_score[icluster] = out.empty() ? -1.0f : out[0];
    }

    outtree->Fill();
  }

  fout->cd();
  outtree->Write();
  fout->Close();
  fin->Close();

  std::cout << "Done. Output: " << outfile_name << std::endl;
}







