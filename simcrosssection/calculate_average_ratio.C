void calculate_average_ratio() {
    const TString baseDir = "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run22/photon20/condorout_waveform/OutDir";
    const int startNum = 0;
    const int endNum = 1999;

    double totalRatio = 0.0;
    int validFiles = 0;

    for (int i = startNum; i <= endNum; ++i) {
        TString filePath = Form("%s%d/caloana.root", baseDir.Data(), i);

        // Check if the file exists
        if (gSystem->AccessPathName(filePath)) {
            printf("Skipping missing file: %s\n", filePath.Data());
            continue;
        }

        // Open the ROOT file
        TFile *file = TFile::Open(filePath, "READ");
        if (!file || file->IsZombie()) {
            printf("Error opening file: %s\n", filePath.Data());
            continue;
        }

        // Retrieve the histogram
        TH1D *hist = (TH1D*)file->Get("sim_cross_counting");
        if (!hist) {
            printf("Histogram 'sim_cross_counting' not found in %s\n", filePath.Data());
            file->Close();
            continue;
        }

        // Check histogram has enough bins
        if (hist->GetNbinsX() < 2) {
            printf("Histogram in %s has fewer than 2 bins\n", filePath.Data());
            file->Close();
            continue;
        }

        const double bin1 = hist->GetBinContent(1);
        const double bin2 = hist->GetBinContent(2);

        // Avoid division by zero
        if (bin1 == 0) {
            printf("Bin 1 is zero in %s, skipping...\n", filePath.Data());
            file->Close();
            continue;
        }

        // Calculate ratio and accumulate
        const double ratio = bin2 / bin1;
        totalRatio += ratio;
        validFiles++;

        file->Close();
    }

    if (validFiles == 0) {
        printf("No valid files processed. Check paths and histograms.\n");
        return;
    }


    std::cout<<"Total ratio: "<<totalRatio<<std::endl;
    const double average = totalRatio / validFiles;
    std::cout<<"Average ratio: "<<average<<std::endl;
    printf("Computed from %d valid files out of %d attempted.\n", validFiles, endNum - startNum + 1);
}