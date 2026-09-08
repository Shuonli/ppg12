#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TSystem.h>
#include <yaml-cpp/yaml.h>

// Iso-scan variant of MakeJetPHOXhisto_chunked.C. IDENTICAL normalisation
// (per chunk: norma = xsec / nb_evt / nseg, |y|<0.7, pdf_weight[0], dir+frag
// summed, divided by bin width) -- the ONLY change is the per-chunk segment
// counts, which must match run_condor_SL7_ct18iso_chunked.sh /
// hadd_ct18iso_chunked.sh.
//
// IMPORTANT: nseg per chunk MUST equal the ACTUAL number of segments that
// hadd_ct18iso_chunked.sh merged (it prints them). If condor jobs failed and
// fewer segments were merged, EDIT the `chunks` vector below to the real
// counts, otherwise the per-chunk normalisation is wrong by actual/nseg.
//
// pdf_tag: "_ct18iso200" (inclusive, ET=200) or "_ct18iso4" (isolated, ET=4).
void MakeJetPHOXhisto_chunked_iso(const std::string &scale = "10",
                                  const std::string &pdf_tag = "_ct18iso200",
                                  const std::string &configname = "/sphenix/user/shuhangli/ppg12/efficiencytool/config_bdt_nom.yaml")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    std::vector<float> pT_bins = configYaml["analysis"]["pT_bins_truth"].as<std::vector<float>>();
    int n_pT_bins = pT_bins.size() - 1;
    double pT_bin_edges[n_pT_bins + 1];
    for (int i = 0; i <= n_pT_bins; i++)
        pT_bin_edges[i] = pT_bins[i];

    const std::string pawres = "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/";

    // (chunk label, ACTUAL number of segments merged into that chunk file).
    // MUST mirror the hadd_ct18iso_chunked.sh report.
    // chunk A is nseg=1: its OutDir0 segment was clobbered by an interactive
    // debug run and set aside, leaving only the full 1M-event OutDir1 segment.
    std::vector<std::pair<std::string, int>> chunks = {
        {"A", 1},
        {"B", 4},
        {"C", 6},
        {"D", 6},
        {"E", 6},
    };

    TH1F *h_truth_pT     = new TH1F("h_truth_pT",     "truth pT (dir+frag)", n_pT_bins, pT_bin_edges);
    TH2F *h_truth_eta_pT = new TH2F("h_truth_eta_pT", "truth eta vs pT",     100, -1, 1, 100, 0, 100);
    h_truth_pT->Sumw2();

    auto fillFromFile = [&](const std::string &fname, float nseg)
    {
        TFile *f = new TFile(fname.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Cannot open " << fname << std::endl;
            return;
        }

        TTree *t2 = (TTree *)f->Get("t2");
        TList *list = t2->GetUserInfo();
        TVectorT<float> *v = static_cast<TVectorT<float> *>(list->At(0));
        float norma = (*v)[1] / (*v)[0] / nseg;  // xsec / nb_evt / nseg, per-chunk

        std::cout << "  -> " << fname
                  << " : xsec=" << (*v)[1] << " nb_evt=" << (*v)[0]
                  << " nseg=" << nseg << " norma=" << norma << std::endl;

        Int_t iprov, ntrack;
        Double_t e[3], px[3], py[3], pz[3], x3, x1, x2;
        Double_t pt[3], y[3];
        Float_t pdf_weight[1000], weight;

        t2->SetBranchAddress("iprov", &iprov);
        t2->SetBranchAddress("ntrack", &ntrack);
        t2->SetBranchAddress("x3", &x3);
        t2->SetBranchAddress("energy", e);
        t2->SetBranchAddress("px", px);
        t2->SetBranchAddress("py", py);
        t2->SetBranchAddress("pz", pz);
        t2->SetBranchAddress("pdf_weight", pdf_weight);

        TH1F *htmp = new TH1F("htmp", "", n_pT_bins, pT_bin_edges);
        htmp->Sumw2();

        Int_t entries = (Int_t)t2->GetEntries();
        for (Int_t i = 0; i < entries; i++)
        {
            if (i % 500000 == 0)
                std::cout << "      " << fname << ": " << i << " / " << entries << std::endl;
            t2->GetEntry(i);
            for (Int_t j = 0; j < ntrack; j++)
            {
                pt[j] = sqrt(px[j]*px[j] + py[j]*py[j]);
                y[j]  = log((e[j]+pz[j]) / (e[j]-pz[j])) * 0.5;
            }
            weight = pdf_weight[0];
            h_truth_eta_pT->Fill(y[0], pt[0], weight);
            if (std::abs(y[0]) < 0.7)
                htmp->Fill(pt[0], weight);
        }

        htmp->Scale(norma);
        for (int ibin = 1; ibin <= htmp->GetNbinsX(); ibin++)
        {
            float bw = htmp->GetBinWidth(ibin);
            htmp->SetBinContent(ibin, htmp->GetBinContent(ibin) / bw);
            htmp->SetBinError  (ibin, htmp->GetBinError  (ibin) / bw);
        }
        std::cout << "  integral(this chunk file) = " << htmp->Integral("width") << std::endl;

        h_truth_pT->Add(htmp);
        delete htmp;
        f->Close();
        delete f;
    };

    for (const auto &cn : chunks) {
        const std::string &label = cn.first;
        const int          nseg  = cn.second;
        std::cout << "=== chunk " << label << " (nseg=" << nseg << ") ===" << std::endl;
        fillFromFile(pawres + "ggdrhic" + pdf_tag + "_chunk" + label + "_" + scale + ".root", nseg);
        fillFromFile(pawres + "ggorhic" + pdf_tag + "_chunk" + label + "_" + scale + ".root", nseg);
    }

    std::cout << "Combined integral over all chunks: " << h_truth_pT->Integral("width") << std::endl;

    std::string outputname = "rootFiles/jetPHOX" + pdf_tag + "_" + scale + "_chunked.root";
    TFile *fout = new TFile(outputname.c_str(), "RECREATE");
    h_truth_pT->Write();
    h_truth_eta_pT->Write();
    fout->Write();
    fout->Close();
    std::cout << "Written to " << outputname << std::endl;
}
