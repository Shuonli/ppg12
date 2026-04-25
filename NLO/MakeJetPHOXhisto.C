#include <iostream>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TSystem.h>
#include <yaml-cpp/yaml.h>

// scale: "05" (mu=0.5*ET), "10" (nominal), "20" (mu=2*ET)
// pdf_tag: "" (default, CT14lo, reads ggd/orhic_{scale}.root) or "_nlo" (CT14nlo rerun)
void MakeJetPHOXhisto(const std::string &scale = "10",
                      const std::string &pdf_tag = "",
                      const std::string &configname = "/sphenix/user/shuhangli/ppg12/efficiencytool/config_bdt_nom.yaml")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    std::vector<float> pT_bins = configYaml["analysis"]["pT_bins_truth"].as<std::vector<float>>();
    int n_pT_bins = pT_bins.size() - 1;
    double pT_bin_edges[n_pT_bins + 1];
    for (int i = 0; i <= n_pT_bins; i++)
        pT_bin_edges[i] = pT_bins[i];

    const float nseg = 100;
    const std::string pawres = "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/";

    TH1F *h_truth_pT    = new TH1F("h_truth_pT",    "truth pT (dir+frag)", n_pT_bins, pT_bin_edges);
    TH2F *h_truth_eta_pT = new TH2F("h_truth_eta_pT", "truth eta vs pT",   100, -1, 1, 100, 0, 100);
    h_truth_pT->Sumw2();

    // Process one JETPHOX file and add its contribution to h and h2
    auto fillFromFile = [&](const std::string &fname)
    {
        TFile *f = new TFile(fname.c_str(), "READ");
        if (!f || f->IsZombie()) { std::cerr << "Cannot open " << fname << std::endl; return; }

        TTree *t2 = (TTree *)f->Get("t2");
        TList *list = t2->GetUserInfo();
        list->Print();
        TVectorT<float> *v = static_cast<TVectorT<float> *>(list->At(0));
        float norma = (*v)[1] / (*v)[0] / nseg;  // xsec / nb_evt / nseg

        Int_t iprov, ntrack;
        Double_t e[3], px[3], py[3], pz[3], x3, x1, x2;
        Double_t pt[3], y[3];
        Float_t pdf_weight[1000], weight;

        t2->SetBranchAddress("iprov",      &iprov);
        t2->SetBranchAddress("ntrack",     &ntrack);
        t2->SetBranchAddress("x3",         &x3);
        t2->SetBranchAddress("energy",     e);
        t2->SetBranchAddress("px",         px);
        t2->SetBranchAddress("py",         py);
        t2->SetBranchAddress("pz",         pz);
        t2->SetBranchAddress("pdf_weight", pdf_weight);

        TH1F *htmp = new TH1F("htmp", "", n_pT_bins, pT_bin_edges);
        htmp->Sumw2();

        Int_t entries = (Int_t)t2->GetEntries();
        for (Int_t i = 0; i < entries; i++)
        {
            if (i % 100000 == 0)
                std::cout << fname << ": " << i << " / " << entries << std::endl;
            t2->GetEntry(i);
            for (Int_t j = 0; j < ntrack; j++)
            {
                pt[j] = sqrt(px[j]*px[j] + py[j]*py[j]);
                y[j]  = log((e[j]+pz[j])/(e[j]-pz[j])) * 0.5;
            }
            weight = pdf_weight[0];
            h_truth_eta_pT->Fill(y[0], pt[0], weight);
            if (abs(y[0]) < 0.7)
                htmp->Fill(pt[0], weight);
        }

        htmp->Scale(norma);
        for (int ibin = 1; ibin <= htmp->GetNbinsX(); ibin++)
        {
            float bw = htmp->GetBinWidth(ibin);
            htmp->SetBinContent(ibin, htmp->GetBinContent(ibin) / bw);
            htmp->SetBinError  (ibin, htmp->GetBinError  (ibin) / bw);
        }
        std::cout << fname << ": integral = " << htmp->Integral() << std::endl;

        h_truth_pT->Add(htmp);
        delete htmp;
        f->Close();
        delete f;
    };

    fillFromFile(pawres + "ggdrhic" + pdf_tag + "_" + scale + ".root");  // direct photon
    fillFromFile(pawres + "ggorhic" + pdf_tag + "_" + scale + ".root");  // fragmentation

    std::cout << "Combined integral: " << h_truth_pT->Integral() << std::endl;

    std::string outputname = "rootFiles/jetPHOX" + pdf_tag + "_" + scale + ".root";
    TFile *fout = new TFile(outputname.c_str(), "RECREATE");
    h_truth_pT->Write();
    h_truth_eta_pT->Write();
    fout->Write();
    fout->Close();
    std::cout << "Written to " << outputname << std::endl;
}
