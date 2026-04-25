#include "plotcommon.h"

void plot_photon10_etaphi_acceptance()
{
    init_plot();

    const char *fname = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/photon10/condorout/combined.root";

    TFile *f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) { std::cerr << "cannot open " << fname << std::endl; return; }

    TTree *t = (TTree *) f->Get("slimtree");
    if (!t) { std::cerr << "no slimtree" << std::endl; return; }

    // 96 eta towers cover [-1.1, 1.1]; 256 phi towers cover [-pi, pi]
    TH2F *h = new TH2F("h_etaphi", "", 96, -1.1, 1.1, 256, -TMath::Pi(), TMath::Pi());
    h->SetXTitle("#it{#eta^{cluster}}");
    h->SetYTitle("#it{#phi^{cluster}} [rad]");
    h->SetZTitle("clusters / (#Delta#eta #times #Delta#phi)");

    TTreeReader R(t);
    TTreeReaderValue<int>  ncluster(R, "ncluster_CLUSTERINFO_CEMC");
    TTreeReaderArray<float> cEt (R, "cluster_Et_CLUSTERINFO_CEMC");
    TTreeReaderArray<float> cEta(R, "cluster_Eta_CLUSTERINFO_CEMC");
    TTreeReaderArray<float> cPhi(R, "cluster_Phi_CLUSTERINFO_CEMC");
    TTreeReaderValue<float> vtxz(R, "vertexz");

    long long nev = 0;
    long long nev_pass = 0;
    long long nclus_fill = 0;

    while (R.Next())
    {
        nev++;
        if (std::fabs(*vtxz) > 10.0) continue;
        nev_pass++;
        for (int i = 0; i < *ncluster; ++i)
        {
            if (cEt[i] < 10.0) continue;
            h->Fill(cEta[i], cPhi[i]);
            nclus_fill++;
        }
    }

    std::cout << "events total = " << nev
              << ", pass |vz|<10 = " << nev_pass
              << ", clusters filled (ET>10) = " << nclus_fill << std::endl;

    TCanvas *c = new TCanvas("c", "", 900, 700);
    c->SetRightMargin(0.16);
    c->SetLogz();

    h->Draw("COLZ");

    TLatex lx;
    lx.SetNDC();
    lx.SetTextSize(0.035);
    lx.DrawLatex(0.14, 0.945, strleg1.c_str());
    lx.DrawLatex(0.14, 0.905, "PYTHIA photon10, |z_{vtx}| < 10 cm, #it{E}_{T}^{cluster} > 10 GeV");

    const char *outdir = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures";
    c->SaveAs(Form("%s/photon10_etaphi_acceptance.pdf", outdir));
    c->SaveAs(Form("%s/photon10_etaphi_acceptance.png", outdir));

    TFile *fout = TFile::Open(Form("%s/photon10_etaphi_acceptance.root", outdir), "RECREATE");
    h->Write();
    fout->Close();

    f->Close();
}
