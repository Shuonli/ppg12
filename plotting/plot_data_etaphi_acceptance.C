#include "plotcommon.h"

void plot_data_etaphi_acceptance()
{
    init_plot();

    const char *glob_pattern = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/condorout/part_*_with_bdt_split.root";

    TChain *t = new TChain("slimtree");
    int n = t->Add(glob_pattern);
    std::cout << "added " << n << " files" << std::endl;

    TH2F *h = new TH2F("h_etaphi", "", 96, -1.1, 1.1, 256, -TMath::Pi(), TMath::Pi());
    h->SetXTitle("#it{#eta^{cluster}}");
    h->SetYTitle("#it{#phi^{cluster}} [rad]");
    h->SetZTitle("clusters / (#Delta#eta #times #Delta#phi)");

    TTreeReader R(t);
    TTreeReaderValue<int>   ncluster(R, "ncluster_CLUSTERINFO_CEMC");
    TTreeReaderArray<float> cEt (R, "cluster_Et_CLUSTERINFO_CEMC");
    TTreeReaderArray<float> cEta(R, "cluster_Eta_CLUSTERINFO_CEMC");
    TTreeReaderArray<float> cPhi(R, "cluster_Phi_CLUSTERINFO_CEMC");
    TTreeReaderValue<float> vtxz(R, "vertexz");
    TTreeReaderArray<bool>  scaledtrig(R, "scaledtrigger");

    long long nev = 0, nev_vz = 0, nev_trig = 0;
    long long nclus_fill = 0;

    while (R.Next())
    {
        nev++;
        if (std::fabs(*vtxz) > 10.0) continue;
        nev_vz++;
        // require Photon_4_GeV (bit 30, Map2) -- matches analysis selection
        if (scaledtrig.GetSize() <= 30 || scaledtrig[30] == 0) continue;
        nev_trig++;
        for (int i = 0; i < *ncluster; ++i)
        {
            if (cEt[i] < 10.0) continue;
            h->Fill(cEta[i], cPhi[i]);
            nclus_fill++;
        }
    }

    std::cout << "events total        = " << nev       << "\n"
              << "events |vz|<10      = " << nev_vz    << "\n"
              << "events + bit30 fired= " << nev_trig  << "\n"
              << "clusters ET>10 kept = " << nclus_fill << std::endl;

    TCanvas *c = new TCanvas("c", "", 900, 700);
    c->SetRightMargin(0.16);
    c->SetLogz();

    h->Draw("COLZ");

    TLatex lx;
    lx.SetNDC();
    lx.SetTextSize(0.035);
    lx.DrawLatex(0.14, 0.945, strleg1.c_str());
    lx.DrawLatex(0.14, 0.905, "Data (Photon 4 GeV trig), |z_{vtx}| < 10 cm, #it{E}_{T}^{cluster} > 10 GeV");

    const char *outdir = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures";
    c->SaveAs(Form("%s/data_etaphi_acceptance.pdf", outdir));
    c->SaveAs(Form("%s/data_etaphi_acceptance.png", outdir));

    TFile *fout = TFile::Open(Form("%s/data_etaphi_acceptance.root", outdir), "RECREATE");
    h->Write();
    fout->Close();
}
