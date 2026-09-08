// mergeFiles.C
#include "TFileMerger.h"
#include "TError.h"
#include <iostream>
#include <yaml-cpp/yaml.h>

void MergeSim(const std::string &configname = "config.yaml") {
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    std::string var_type = configYaml["output"]["var_type"].as<std::string>();
    std::string eff_out  = configYaml["output"]["eff_outfile"].as<std::string>();
    std::string resp_out = configYaml["output"]["response_outfile"].as<std::string>();

    // [PPG12] Herwig signal cross-check: replace the Pythia photon5/10/20
    // signal trio with Herwig photon10_herwig + photon20_herwig (no Herwig
    // photon5 analog — shape mismatch + xsec issues, see reports/herwig_*).
    // Background jet samples remain Pythia (no Herwig jet samples processed).
    // Pythia path below this if-block is unchanged.
    if (var_type.find("herwig") != std::string::npos) {
        std::cout << "[MergeSim] HERWIG cross-check mode (var_type=" << var_type << ")\n";

        // --- Photon-jet signal: 2 Herwig samples ---
        std::string out_photon = eff_out + "_" + var_type + ".root";
        TFileMerger m_photon;
        m_photon.OutputFile(out_photon.c_str(), "RECREATE");
        m_photon.AddFile((eff_out + "_photon10_herwig_" + var_type + ".root").c_str());
        m_photon.AddFile((eff_out + "_photon20_herwig_" + var_type + ".root").c_str());
        if (!m_photon.Merge()) { ::Error("MergeSim", "HERWIG photon merge failed"); return; }
        std::cout << "  -> " << out_photon << std::endl;

        // --- Inclusive jets: Pythia jets (same set as nominal) ---
        std::string out_jet = eff_out + "_jet_" + var_type + ".root";
        TFileMerger m_jet;
        m_jet.OutputFile(out_jet.c_str(), "RECREATE");
        for (auto j : {"jet8", "jet12", "jet20", "jet30", "jet40"}) {
            m_jet.AddFile((eff_out + "_" + j + "_" + var_type + ".root").c_str());
        }
        if (!m_jet.Merge()) { ::Error("MergeSim", "HERWIG jet merge failed"); return; }
        std::cout << "  -> " << out_jet << std::endl;

        // --- Response matrix (Herwig signal only) ---
        std::string out_resp = resp_out + "_" + var_type + ".root";
        TFileMerger m_resp;
        m_resp.OutputFile(out_resp.c_str(), "RECREATE");
        m_resp.AddFile((resp_out + "_photon10_herwig_" + var_type + ".root").c_str());
        m_resp.AddFile((resp_out + "_photon20_herwig_" + var_type + ".root").c_str());
        if (!m_resp.Merge()) { ::Error("MergeSim", "HERWIG response merge failed"); return; }
        std::cout << "  -> " << out_resp << std::endl;

        return;
    }

    // ============================================================
    // Pythia path (unchanged)
    // ============================================================
    TFileMerger merger;

    std::string infilename1 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "photon5" + "_" + var_type + ".root";
    std::string infilename2 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "photon10" + "_" + var_type + ".root";
    std::string infilename3 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "photon20" + "_" + var_type + ".root";

    std::string outfilename = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + var_type + ".root";

    // jet5 is dropped from the nominal DI-blending pipeline (no jet5_double
    // MC partner), so the merged inclusive jet sample starts at jet8.
    std::string infilenamejet8 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "jet8" + "_" + var_type + ".root";
    std::string infilenamejet12 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "jet12" + "_" + var_type + ".root";
    std::string infilenamejet20 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "jet20" + "_" + var_type + ".root";
    std::string infilenamejet30 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "jet30" + "_" + var_type + ".root";
    std::string infilenamejet40 = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + "jet40" + "_" + var_type + ".root";

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

    merger_jet.AddFile(infilenamejet8.c_str());
    merger_jet.AddFile(infilenamejet12.c_str());
    merger_jet.AddFile(infilenamejet20.c_str());
    merger_jet.AddFile(infilenamejet30.c_str());
    merger_jet.AddFile(infilenamejet40.c_str());

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
