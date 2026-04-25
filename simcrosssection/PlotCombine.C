#include "draw.C"

void PlotCombine(){

    SetAtlasStyle();
    std::string photon5energyfile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_photon5.root";
    std::string photon10energyfile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_photon10.root";
    std::string photon20energyfile = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_photon20.root";
    //float photon5cross = 2.017e+08 * 0.000442571;
    //float photon10cross = 3.690e+07 * 0.000181474;
    //float photon20cross = 1.571e+05 * 0.000673448;

    float photon5cross = 1.0;
    float photon10cross = 1.0;
    float photon20cross = 1.0;

    TFile* fphoton5 = TFile::Open(photon5energyfile.c_str(), "READ");
    TFile* fphoton10 = TFile::Open(photon10energyfile.c_str(), "READ");
    TFile* fphoton20 = TFile::Open(photon20energyfile.c_str(), "READ");

    TH1D* hphoton5 = (TH1D*)fphoton5->Get("h_max_photon_pT")->Clone("hphoton5");
    TH1D* hphoton10 = (TH1D*)fphoton10->Get("h_max_photon_pT")->Clone("hphoton10");
    TH1D* hphoton20 = (TH1D*)fphoton20->Get("h_max_photon_pT")->Clone("hphoton20");

    //set axis limits
    hphoton5->GetXaxis()->SetRangeUser(0, 100);
    hphoton10->GetXaxis()->SetRangeUser(0, 100);
    hphoton20->GetXaxis()->SetRangeUser(0, 100);

    hphoton5->Sumw2();
    hphoton10->Sumw2();
    hphoton20->Sumw2();

    //normalized by entries
    //hphoton5->Scale(1.0/hphoton5->GetEntries());
    //hphoton10->Scale(1.0/hphoton10->GetEntries());
    //hphoton20->Scale(1.0/hphoton20->GetEntries());

    //use photon20 as base
    hphoton5->Scale(photon5cross / photon20cross);
    hphoton10->Scale(photon10cross / photon20cross);
    hphoton20->Scale(1.0);

    TH1D* h20over10 = (TH1D*)hphoton20->Clone("h20over10");
    h20over10->Divide(hphoton10);
    TH1D* h20over5 = (TH1D*)hphoton20->Clone("h20over5");
    h20over5->Divide(hphoton5);



    std::vector<TH1F*> h_input;

    h_input.push_back((TH1F*)hphoton5);
    h_input.push_back((TH1F*)hphoton10);
    h_input.push_back((TH1F*)hphoton20);



    std::vector<int> colors = {kPink+8, kSpring-7, kAzure-3, kViolet+3, kOrange +10};

    std::vector<int> markerstyle = {20, 20, 20, 20, 20};

    std::vector<std::string> legend = {"photon 5 GeV", "photon 10 GeV", "photon 20 GeV"};

    std::vector<std::string> text = {"Pythia, #sqrt{s} = 200 GeV"};

    
    
        draw_1D_multiple_plot(h_input, colors, markerstyle,
                           false, 1 , false,
                           true, 0, 50, false,
                           true, 100, 2E10, true,
                           true, "Max photon pT in event [GeV]", "",
                           false, "",
                           true, text, 0.2, 0.83, 0.04,
                           true, legend, 0.6, 0.8, 0.04,
                           "test");
    

    std::vector<TH1F*> h_input_ratio;


    h_input_ratio.push_back((TH1F*)h20over5);  
    h_input_ratio.push_back((TH1F*)h20over10);

    //h20over10->Draw();
    
    /*
    std::vector<std::string> legend_ratio = {"photon20/photon5", "photon20/photon10"};
    
    draw_1D_multiple_plot(h_input_ratio, colors, markerstyle,
                           false, 1 , true ,
                           false, 0.5, 1.5, false,
                           false, 0, 1, false,
                           true, "event energy scale", "",
                           true, "",
                           true, text, 0.2, 0.8, 0.03,
                           true, legend_ratio, 0.6, 0.8, 0.03,
                           "test_ratio");
    */

}
