#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

#include <yaml-cpp/yaml.h>

namespace
{
  void write_csv_header(std::ofstream& out, const int arraysize)
  {
    out << "event,run,cluster_idx,cluster_Et,cluster_Eta,vertexz,particle_photonclass";
    for (int i = 0; i < arraysize; ++i)
    {
      out << ",e" << i;
    }
    for (int i = 0; i < arraysize; ++i)
    {
      out << ",own" << i;
    }
    out << "\n";
  }
}  // namespace

// Export one row per cluster into CSV, including:
// - cluster_Et, cluster_Eta, vertexz
// - cluster_e_array (49 floats) and cluster_ownership_array (49 ints)
// - ground-truth label: particle_photonclass (matched via cluster_truthtrkID -> particle_trkid)
//
// Output:
// - <outprefix>_<filetype>.csv   (single file; use particle_photonclass as the label column)
//
// Notes:
// - The 2D branches cluster_e_array_* and cluster_ownership_array_* are stored as
//   [ncluster][49] in the TTree. ROOT exposes these via TTreeReaderArray as a
//   flattened 1D array, so index = cluster_idx * arraysize + cell_idx.
void export_cluster_arrays_csv(const std::string& configname = "config_nom.yaml",
                               const std::string& filetype = "photon5",
                               const std::string& outprefix = "clusters",
                               const Long64_t maxEvents = -1,
                               const int arraysize = 49)
{
  // yaml-cpp load (try default first, then sPHENIX user path as fallback)
  if (gSystem->Load("libyaml-cpp.so") < 0)
  {
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
  }

  YAML::Node configYaml = YAML::LoadFile(configname);
  const std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
  const std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
  const std::string tree_name = configYaml["input"]["tree"].as<std::string>();
  const std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();

  std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;
  if (filetype == "data")
  {
    std::cerr << "Error: filetype=='data' does not have truth labels (particle_photonclass)." << std::endl;
    return;
  }

  TFile* fin = TFile::Open(infilename.c_str(), "READ");
  if (!fin || fin->IsZombie())
  {
    std::cerr << "Error: failed to open input file: " << infilename << std::endl;
    return;
  }

  TTree* t = dynamic_cast<TTree*>(fin->Get(tree_name.c_str()));
  if (!t)
  {
    std::cerr << "Error: failed to get TTree '" << tree_name << "' from file: " << infilename << std::endl;
    fin->Close();
    return;
  }

  const std::string out_name = outprefix + "_" + filetype + ".csv";
  std::ofstream out(out_name);
  if (!out.is_open())
  {
    std::cerr << "Error: failed to open output CSV: " << out_name << std::endl;
    fin->Close();
    return;
  }

  write_csv_header(out, arraysize);

  TTreeReader reader(t);

  // Event-level
  TTreeReaderValue<int> runnumber(reader, "runnumber");
  TTreeReaderValue<float> vertexz(reader, "vertexz");

  // Truth particle arrays
  TTreeReaderArray<int> particle_trkid(reader, "particle_trkid");
  TTreeReaderArray<int> particle_photonclass(reader, "particle_photonclass");

  // Cluster arrays (dynamic names)
  TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
  TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusternodename.c_str()));
  TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
  TTreeReaderArray<int> cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusternodename.c_str()));

  // 2D arrays flattened to 1D: [ncluster][arraysize] -> [ncluster*arraysize]
  TTreeReaderArray<float> cluster_e_array(reader, Form("cluster_e_array_%s", clusternodename.c_str()));
  TTreeReaderArray<int> cluster_ownership_array(reader, Form("cluster_ownership_array_%s", clusternodename.c_str()));

  Long64_t event_idx = 0;
  while (reader.Next())
  {
    if (maxEvents >= 0 && event_idx >= maxEvents)
    {
      break;
    }

    if (event_idx % 10000 == 0)
    {
      std::cout << "Processing event " << event_idx << std::endl;
    }

    // Build map: truth track id -> particle index
    std::unordered_map<int, int> trkid_to_idx;
    trkid_to_idx.reserve(particle_trkid.GetSize());
    for (int ip = 0; ip < (int)particle_trkid.GetSize(); ++ip)
    {
      trkid_to_idx[particle_trkid[ip]] = ip;
    }

    const int ncl = *ncluster;
    for (int icl = 0; icl < ncl; ++icl)
    {
      int label = -1;
      const int truth_trkid = cluster_truthtrkID[icl];
      auto it = trkid_to_idx.find(truth_trkid);
      if (it != trkid_to_idx.end())
      {
        const int ip = it->second;
        if (ip >= 0 && ip < (int)particle_photonclass.GetSize())
        {
          label = particle_photonclass[ip];
        }
      }

      out << event_idx << "," << *runnumber << "," << icl << ","
          << cluster_Et[icl] << "," << cluster_Eta[icl] << "," << *vertexz << ","
          << label;

      // Write 49 energies
      for (int itower = 0; itower < arraysize; ++itower)
      {
        const Long64_t idx = (Long64_t)icl * arraysize + itower;
        float v = 0.0f;
        if (idx >= 0 && idx < (Long64_t)cluster_e_array.GetSize())
        {
          v = cluster_e_array[(int)idx];
        }
        out << "," << v;
      }

      // Write 49 ownership flags
      for (int itower = 0; itower < arraysize; ++itower)
      {
        const Long64_t idx = (Long64_t)icl * arraysize + itower;
        int v = 0;
        if (idx >= 0 && idx < (Long64_t)cluster_ownership_array.GetSize())
        {
          v = cluster_ownership_array[(int)idx];
        }
        out << "," << v;
      }

      out << "\n";
    }

    ++event_idx;
  }

  out.close();
  fin->Close();

  std::cout << "Wrote:\n"
            << "  " << out_name << "\n";
}


