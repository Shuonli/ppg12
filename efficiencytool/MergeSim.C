// mergeFiles.C
#include "TFileMerger.h"
#include "TError.h"
#include <iostream>
#include <yaml-cpp/yaml.h>

void MergeSim(const std::string &configname = "config.yaml") {
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);
    // Create a TFileMerger instance
    TFileMerger merger;

    std::string var_type = configYaml["output"]["var_type"].as<std::string>();

    std::string infilename1 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "photon5" + "_" + var_type + ".root";
    std::string infilename2 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "photon10" + "_" + var_type + ".root";
    std::string infilename3 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "photon20" + "_" + var_type + ".root";

    std::string outfilename = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + var_type + ".root";

    std::string infilenamejet10 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "jet10" + "_" + var_type + ".root";
    std::string infilenamejet15 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "jet15" + "_" + var_type + ".root";
    std::string infilenamejet20 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "jet20" + "_" + var_type + ".root";
    std::string infilenamejet30 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "jet30" + "_" + var_type + ".root";

    std::string outfilenamejet = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "jet" + "_" + var_type + ".root";

    std::string infilenameresponse1 = configYaml["output"]["response_outfile"].as<std::string>() + "_" + "photon5" + "_" + var_type + ".root";
    std::string infilenameresponse2 = configYaml["output"]["response_outfile"].as<std::string>() + "_" + "photon10" + "_" + var_type + ".root";
    std::string infilenameresponse3 = configYaml["output"]["response_outfile"].as<std::string>() + "_" + "photon20" + "_" + var_type + ".root";

    std::string outfilenameresponse = configYaml["output"]["response_outfile"].as<std::string>() + "_" + var_type + ".root";
    
    // Set the name of the output file.
    merger.OutputFile(outfilename.c_str(), "RECREATE");

    // Add the files you want to merge
    merger.AddFile(infilename1.c_str());
    merger.AddFile(infilename2.c_str());
    merger.AddFile(infilename3.c_str());


    // Perform the merge
    if (!merger.Merge()) {
        ::Error("mergeFiles", "Merge failed!");
        return;
    }
    
    std::cout << "Files merged successfully into " << outfilename << std::endl;

    //merge the inclusive samples

    TFileMerger merger_jet;
    
    merger_jet.OutputFile(outfilenamejet.c_str(), "RECREATE");

    merger_jet.AddFile(infilenamejet10.c_str());
    merger_jet.AddFile(infilenamejet15.c_str());
    merger_jet.AddFile(infilenamejet20.c_str());
    merger_jet.AddFile(infilenamejet30.c_str());

    // Perform the merge
    if (!merger_jet.Merge()) {
        ::Error("mergeFiles", "Merge failed!");
        return;
    }

    std::cout << "Files merged successfully into " << outfilenamejet << std::endl;


    // Create a TFileMerger instance
    TFileMerger merger_response;

    // Set the name of the output file.
    merger_response.OutputFile(outfilenameresponse.c_str(), "RECREATE");

    // Add the files you want to merge
    merger_response.AddFile(infilenameresponse1.c_str());
    merger_response.AddFile(infilenameresponse2.c_str());
    merger_response.AddFile(infilenameresponse3.c_str());

    // Perform the merge
    if (!merger_response.Merge()) {
        ::Error("mergeFiles", "Merge failed!");
        return;
    }

    std::cout << "Files merged successfully into " << outfilenameresponse << std::endl;

}
