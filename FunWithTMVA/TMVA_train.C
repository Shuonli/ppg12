#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TString.h>
#include <iostream>

// Include TMVA headers
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>

void TMVA_train(int ivtx, int ieta, int iet)
{

    // Open the input ROOT file containing your trees
    TFile *inputFileSignal5 = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_tree_photon5_0.root", "READ");
    TTree *signalTree5 = (TTree *)inputFileSignal5->Get("tree");
    TFile *inputFileSignal10 = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_tree_photon10_0.root", "READ");
    TTree *signalTree10 = (TTree *)inputFileSignal10->Get("tree");
    TFile *inputFileSignal20 = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_tree_photon20_0.root", "READ");
    TTree *signalTree20 = (TTree *)inputFileSignal20->Get("tree");
    TFile *inputFileBackground10 = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_tree_jet10_0.root", "READ");
    TTree *backgroundTree10 = (TTree *)inputFileBackground10->Get("tree");
    TFile *inputFileBackground15 = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_tree_jet15_0.root", "READ");
    TTree *backgroundTree15 = (TTree *)inputFileBackground15->Get("tree");
    TFile *inputFileBackground20 = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_tree_jet20_0.root", "READ");
    TTree *backgroundTree20 = (TTree *)inputFileBackground20->Get("tree");
    TFile *inputFileBackground30 = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_tree_jet30_0.root", "READ");
    TTree *backgroundTree30 = (TTree *)inputFileBackground30->Get("tree");

    // unit in pb
    const float photon5cross = 2.017e+08 * 0.000442571;
    const float photon10cross = 3.690e+07 * 0.000181474;
    const float photon20cross = 1.571e+05 * 0.000673448;

    // Hanpu uses unit in pb
    const float jet10cross = 3.646e-6 * 1e12;
    const float jet15cross = 36864930.0 * 0.011059973;
    const float jet20cross = 1392140.9 * 0.042;
    const float jet30cross = 2.505e-9 * 1e12;

    float extra_bg_weight = 1.0;

    float photon5weight = photon5cross / photon20cross;
    float photon10weight = photon10cross / photon20cross;
    float photon20weight = 1.0;
    float jet10weight = jet10cross / photon20cross * extra_bg_weight;
    float jet15weight = jet15cross / photon20cross * extra_bg_weight;
    float jet20weight = jet20cross / photon20cross * extra_bg_weight;
    float jet30weight = jet30cross / photon20cross * extra_bg_weight;

    // Open the output file for TMVA results
    TFile *outputFile = TFile::Open(Form("results/cor_TMVA_BinOptimization_v%d_eta%d_ET%d.root", ivtx, ieta, iet), "RECREATE");
    outputFile->cd();

    // Define binning for vertex, eta, and ET (example bin edges)
    /*
        const int nVertexBins = 8;
        double vertexBins[nVertexBins + 1] = {-40.0, -30, -20.0, 10 ,0.0, 10.0 ,20.0,30.0 ,40.0};

        const int nEtaBins = 6;
        double etaBins[nEtaBins + 1] = {-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9};

        const int nETBins = 3;
        double etBins[nETBins + 1] = {5, 10, 20, 50};
    */

    const int nVertexBins = 1;
    double vertexBins[nVertexBins + 1] = {-30.0, 30.0};

    const int nEtaBins = 1;
    double etaBins[nEtaBins + 1] = {-0.7, 0.7};

    const int nETBins = 5;
    double etBins[nETBins + 1] = {8, 10, 14, 20, 30, 35};

    std::ofstream outFile(Form("results/TMVA_BinOptimization_v%d_eta%d_ET%d.txt", ivtx, ieta, iet));
    if (!outFile)
    {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    std::streambuf *originalCoutBuffer = std::cout.rdbuf();

    // Loop over all bins in vertex, eta, and ET
    for (int iVertex = 0; iVertex < nVertexBins; iVertex++)
    {
        for (int iEta = 0; iEta < nEtaBins; iEta++)
        {
            for (int iET = 0; iET < nETBins; iET++)
            {
                if (!(iVertex == ivtx && iEta == ieta && iET == iet))
                    continue;
                // Build a selection cut string for the current bin.
                TString binCut = Form("vertex_z > %f && vertex_z <= %f && cluster_eta > %f && cluster_eta <= %f && cluster_ET > %f && cluster_ET <= %f",
                                      vertexBins[iVertex], vertexBins[iVertex + 1],
                                      etaBins[iEta], etaBins[iEta + 1],
                                      etBins[iET], etBins[iET + 1]);
                std::cout << "Processing bin with cuts: " << binCut.Data() << std::endl;

                // Create a unique name for this bin’s analysis
                TString dataloaderName = Form("dataset_bin_v%d_eta%d_ET%d", iVertex, iEta, iET);
                TString factoryName = Form("TMVA_Factory_bin_v%d_eta%d_ET%d", iVertex, iEta, iET);

                // Create a TMVA factory for this bin
                TMVA::Factory *factory = new TMVA::Factory(factoryName, outputFile,
                                                           "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");

                // Create a DataLoader for this bin
                TMVA::DataLoader *dataloader = new TMVA::DataLoader(dataloaderName);

                // Add the variable on which you want to optimize a cut.
                // In this example, "x" is the variable to cut on.
                // dataloader->AddVariable("e17_to_e77", 'F');
                // dataloader->AddVariable("e37_to_e77", 'F');
                // dataloader->AddVariable("e32_to_e35", 'F');

                //dataloader->AddVariable("prob", 'F');
                //dataloader->AddVariable("CNN_prob", 'F');
                dataloader->AddVariable("e17_to_e77", 'F');
                dataloader->AddVariable("e37_to_e77", 'F');
                dataloader->AddVariable("e32_to_e35", 'F');
                //dataloader->AddVariable("e33_to_e35", 'F');
                dataloader->AddVariable("e11_to_e33", 'F');
                //dataloader->AddVariable("e11_to_E", 'F');
                //dataloader->AddVariable("e33_to_E", 'F');
                //dataloader->AddVariable("hcalet33_to_ettot", 'F');
                // dataloader->AddVariable("ihcalet33_to_ettot", 'F');
                // dataloader->AddVariable("ohcalet33_to_ettot", 'F');
                // dataloader->AddVariable("hcalet22_to_ettot", 'F');
                // dataloader->AddVariable("ihcalet22_to_ettot", 'F');
                // dataloader->AddVariable("ohcalet22_to_ettot", 'F');
                // dataloader->AddVariable("detamax", 'F');
                // dataloader->AddVariable("dphimax", 'F');
                dataloader->AddVariable("et1", 'F');
                dataloader->AddVariable("et2", 'F');
                dataloader->AddVariable("et3", 'F');
                dataloader->AddVariable("et4", 'F');
                // dataloader->AddVariable("weta", 'F');
                // dataloader->AddVariable("wphi", 'F');
                dataloader->AddVariable("w32", 'F');
                //dataloader->AddVariable("w52", 'F');
                //dataloader->AddVariable("w72", 'F');
                // dataloader->AddVariable("wr", 'F');
                // dataloader->AddVariable("wrr", 'F');
                // dataloader->AddVariable("weta_cog", 'F');
                // dataloader->AddVariable("wphi_cog", 'F');
                dataloader->AddVariable("weta_cogx", 'F');
                dataloader->AddVariable("wphi_cogx", 'F');

                // dataloader->AddSpectator("vertex_z", 'F');

                //dataloader->AddVariable("cluster_iso_03");
                // dataloader->AddVariable("cluster_iso_03_emcal");
                //dataloader->AddVariable("cluster_iso_03_hcalin");
                //dataloader->AddVariable("cluster_iso_03_hcalout");
                //dataloader->AddVariable("cluster_iso_03_60_emcal");
                //dataloader->AddVariable("cluster_iso_03_60_hcalin");
                //dataloader->AddVariable("cluster_iso_03_60_hcalout");
                //dataloader->AddVariable("cluster_iso_03_120_emcal");
                //dataloader->AddVariable("cluster_iso_03_120_hcalin");
                //dataloader->AddVariable("cluster_iso_03_120_hcalout");
                dataloader->AddVariable("cluster_iso_03_60_emcal + cluster_iso_03_60_hcalin + cluster_iso_03_60_hcalout", "cluster_iso_03_60_total","");
                //dataloader->AddVariable("cluster_iso_03_120_emcal + cluster_iso_03_120_hcalin + cluster_iso_03_120_hcalout", "cluster_iso_03_120_total","");
                /*
                dataloader->AddSpectator("cluster_iso_02");
                dataloader->AddSpectator("cluster_iso_03");
                dataloader->AddSpectator("cluster_iso_04");
                dataloader->AddSpectator("cluster_iso_03_emcal");
                dataloader->AddSpectator("cluster_iso_03_hcalin");
                dataloader->AddSpectator("cluster_iso_03_hcalout");
                dataloader->AddSpectator("cluster_iso_03_60_emcal");
                dataloader->AddSpectator("cluster_iso_03_60_hcalin");
                dataloader->AddSpectator("cluster_iso_03_60_hcalout");
                dataloader->AddSpectator("cluster_iso_03_120_emcal");
                dataloader->AddSpectator("cluster_iso_03_120_hcalin");
                dataloader->AddSpectator("cluster_iso_03_120_hcalout");
                dataloader->AddSpectator("cluster_iso_03_60_emcal + cluster_iso_03_60_hcalin + cluster_iso_03_60_hcalout", "cluster_iso_03_60_total");
                dataloader->AddSpectator("cluster_iso_03_120_emcal + cluster_iso_03_120_hcalin + cluster_iso_03_120_hcalout", "cluster_iso_03_120_total");
                */
                TCut selection(binCut);

                // Add trees for signal and background using the bin selection.
                dataloader->AddTree(signalTree5, "Signal", photon5weight, selection);
                dataloader->AddTree(signalTree10, "Signal", photon10weight, selection);
                // dataloader->AddTree(signalTree10, "Signal", 1.0);
                dataloader->AddTree(signalTree20, "Signal", photon20weight, selection);
                dataloader->AddTree(backgroundTree10, "Background", jet10weight, selection);
                dataloader->AddTree(backgroundTree15, "Background", jet15weight, selection);
                // dataloader->AddTree(backgroundTree15, "Background", 1.0);
                dataloader->AddTree(backgroundTree20, "Background", jet20weight, selection);
                dataloader->AddTree(backgroundTree30, "Background", jet30weight, selection);
                // Prepare the training and test trees.
                // The options here set the number of training events to all available events,
                // with a random split.
                dataloader->PrepareTrainingAndTestTree("", "",
                                                       "nTrain_Signal=20000:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents");

                // Book the CutsGA method (Genetic Algorithm–optimized cuts)
                // Parameters:
                // - PopSize: number of candidate solutions per generation,
                // - Cycles: number of generations,
                // - CutRangeMin/Max: allowed range for cut values.
                TString methodName = Form("CutsGA_bin_v%d_eta%d_ET%d", iVertex, iEta, iET);
                factory->BookMethod(dataloader, TMVA::Types::kCuts, methodName,
                                    "!H:!V:VarProp=FSmart:PopSize=2:Cycles=1:CutRangeMin=-1:CutRangeMax=2");

                // Train, test, and evaluate the method for this bin.
                factory->TrainAllMethods();
                factory->TestAllMethods();
                factory->EvaluateAllMethods();

                TMVA::IMethod *methodCutsBase = factory->GetMethod(dataloaderName, methodName);
                TMVA::MethodCuts *methodCuts = dynamic_cast<TMVA::MethodCuts *>(methodCutsBase);
                std::cout.rdbuf(outFile.rdbuf());

                // print the eta vertex and ET range
                std::cout << "Vertex range: " << vertexBins[iVertex] << " to " << vertexBins[iVertex + 1] << std::endl;
                std::cout << "Eta range: " << etaBins[iEta] << " to " << etaBins[iEta + 1] << std::endl;
                std::cout << "ET range: " << etBins[iET] << " to " << etBins[iET + 1] << std::endl;
                std::cout << "Method: " << methodName.Data() << std::endl;

                methodCuts->PrintCuts(0.9);
                methodCuts->PrintCuts(0.8);
                methodCuts->PrintCuts(0.7);
                std::cout.rdbuf(originalCoutBuffer);
                std::vector<Double_t> cutsMin;
                std::vector<Double_t> cutsMax;
                Double_t true_effS = methodCuts->GetCuts(0.9, cutsMin, cutsMax);

                std::cout << "true_effS: " << true_effS << std::endl;
                for (size_t i = 0; i < cutsMin.size(); ++i)
                {
                    std::cout << "Cut " << i << ": " << cutsMin[i] << " to " << cutsMax[i] << std::endl;
                }

                // Clean up for this bin
                delete factory;
                delete dataloader;
            }
        }
    }

    // Close the output file
    outputFile->Close();
}
