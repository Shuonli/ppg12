#include <iostream>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TSystem.h>
#include <TEfficiency.h>
#include <yaml-cpp/yaml.h>

void MakeJetPHOXhisto(const std::string &configname = "/sphenix/user/shuhangli/ppg12/efficiencytool/config.yaml")
{

    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    std::vector<float> pT_bins = configYaml["analysis"]["pT_bins_truth"].as<std::vector<float>>();
    int n_pT_bins = pT_bins.size() - 1;
    double pT_bin_edges[n_pT_bins + 1];

    for (int i = 0; i < n_pT_bins + 1; i++)
    {
        pT_bin_edges[i] = pT_bins[i];
    }

    TFile *fjetphox = new TFile("/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/ggdrhic.root", "READ");

    Int_t iprov, ntrack;
    // Float_t e[3],px[3],py[3],pz[3];
    Double_t e[3], px[3], py[3], pz[3];
    Double_t x3;
    // Float_t pt[3],y[3];
    Double_t pt[3], y[3];
    // Float_t x1,x2;
    Double_t x1, x2;
    Float_t pdf_weight[1000];
    Float_t weight;

    TTree *t2 = (TTree *)fjetphox->Get("t2");

    // we get the value stored into the header
    TList *list = t2->GetUserInfo();
    list->Print();
    TVectorT<float> *v = static_cast<TVectorT<float> *>(list->At(0));
    float nb_evt = (*v)[0];
    float xsec = (*v)[1];
    float sqrt_s = (*v)[2];
    float norma = xsec / nb_evt;

    //
    t2->SetBranchAddress("iprov", &iprov);
    t2->SetBranchAddress("ntrack", &ntrack);
    t2->SetBranchAddress("x3", &x3);
    t2->SetBranchAddress("energy", e);
    t2->SetBranchAddress("px", px);
    t2->SetBranchAddress("py", py);
    t2->SetBranchAddress("pz", pz);
    t2->SetBranchAddress("pdf_weight", pdf_weight);
    //

    TH1F *h_truth_pT = new TH1F("h_truth_pT", "truth pT", n_pT_bins, pT_bin_edges);
    Int_t entries = (Int_t)t2->GetEntries();
    //
    for (Int_t i = 0; i < entries; i++)
    {
        t2->GetEntry(i);
        for (Int_t j = 0; j < ntrack; j++)
        {
            pt[j] = sqrt(px[j] * px[j] + py[j] * py[j]);
            y[j] = log((e[j] + pz[j]) / (e[j] - pz[j])) * 0.5;
        }
        weight = pdf_weight[0];
        if (abs(y[0]) < 0.7)
        {
            h_truth_pT->Fill(pt[0], weight);
        }
    }
    h_truth_pT->Scale(norma);

    // loop over bins and scale by bin width
    for (int ibin = 1; ibin <= h_truth_pT->GetNbinsX(); ibin++)
    {
        float binwidth = h_truth_pT->GetBinWidth(ibin);
        h_truth_pT->SetBinContent(ibin, h_truth_pT->GetBinContent(ibin) / binwidth);
        h_truth_pT->SetBinError(ibin, h_truth_pT->GetBinError(ibin) / binwidth);
    }
    h_truth_pT->Draw();

    TFile *systOut = new TFile(Form("rootFiles/jetPHOX.root"), "RECREATE");
    h_truth_pT->Write();
    systOut->Write();
    systOut->Close();

}
