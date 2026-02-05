// merge_cluster_time.C
#include "TFileMerger.h"
#include "TError.h"
#include <iostream>

void merge_cluster_time(const bool use_leading_tower_time = false,
                        const bool apply_leading_time_corr = false)
{
    std::string suffix = "";
    if (use_leading_tower_time)
    {
        suffix = apply_leading_time_corr ? "_leadingTowerTimeCorr" : "_leadingTowerTime";
    }

    // Merge photon samples
    TFileMerger merger_photon;
    merger_photon.OutputFile(Form("results/cluster_time_analysis_photon%s.root", suffix.c_str()), "RECREATE");
    merger_photon.AddFile(Form("results/cluster_time_analysis_photon5%s.root", suffix.c_str()));
    merger_photon.AddFile(Form("results/cluster_time_analysis_photon10%s.root", suffix.c_str()));
    merger_photon.AddFile(Form("results/cluster_time_analysis_photon20%s.root", suffix.c_str()));

    if (!merger_photon.Merge())
    {
        ::Error("merge_cluster_time", "Photon merge failed!");
        return;
    }
    std::cout << "Photon files merged successfully into "
              << Form("results/cluster_time_analysis_photon%s.root", suffix.c_str()) << std::endl;

    // Merge jet samples
    TFileMerger merger_jet;
    merger_jet.OutputFile(Form("results/cluster_time_analysis_jet%s.root", suffix.c_str()), "RECREATE");
    merger_jet.AddFile(Form("results/cluster_time_analysis_jet10%s.root", suffix.c_str()));
    merger_jet.AddFile(Form("results/cluster_time_analysis_jet15%s.root", suffix.c_str()));
    merger_jet.AddFile(Form("results/cluster_time_analysis_jet20%s.root", suffix.c_str()));
    merger_jet.AddFile(Form("results/cluster_time_analysis_jet30%s.root", suffix.c_str()));
    merger_jet.AddFile(Form("results/cluster_time_analysis_jet50%s.root", suffix.c_str()));

    if (!merger_jet.Merge())
    {
        ::Error("merge_cluster_time", "Jet merge failed!");
        return;
    }
    std::cout << "Jet files merged successfully into "
              << Form("results/cluster_time_analysis_jet%s.root", suffix.c_str()) << std::endl;

    // Optionally merge both photon and jet into a single inclusive file
    //TFileMerger merger_all;
    //merger_all.OutputFile("results/cluster_time_analysis_all.root", "RECREATE");
    //merger_all.AddFile("results/cluster_time_analysis_photon.root");
    //merger_all.AddFile("results/cluster_time_analysis_jet.root");

    //if (!merger_all.Merge())
   // {
    //    ::Error("merge_cluster_time", "All merge failed!");
    //    return;
    //}
    //std::cout << "All files merged successfully into results/cluster_time_analysis_all.root" << std::endl;
}
