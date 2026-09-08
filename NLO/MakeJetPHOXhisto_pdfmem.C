#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TSystem.h>
#include <yaml-cpp/yaml.h>

// Per-PDF-member histogramming version of MakeJetPHOXhisto_chunked.C.
// Mirrors the existing chunked pipeline but instead of using only
// pdf_weight[0] (the central member) it loops over m = 0 .. nb_member-1
// and fills one histogram per LHAPDF member. The downstream PDF
// uncertainty band (Hessian or replica recipe) is then computed from
// these histograms by compute_pdf_unc.py.
//
// scale: usually "10" (PDF uncertainty is essentially scale-invariant
//        to LO so only nominal scale is needed).
// pdf_tag: "_ct18", "_nlo", "_nnpdf", "_nnpdf4", "_cteq", "_msht"
//
// Output: rootFiles/jetPHOX<pdf_tag>_<scale>_pdfmem.root, containing
// h_truth_pT_m0, h_truth_pT_m1, ... (one TH1F per LHAPDF member).
void MakeJetPHOXhisto_pdfmem(const std::string &scale = "10",
                             const std::string &pdf_tag = "_ct18",
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

    // Same chunk layout as MakeJetPHOXhisto_chunked.C
    std::vector<std::pair<std::string, int>> chunks = {
        {"A", 1},
        {"B", 2},
        {"C", 3},
        {"D", 3},
        {"E", 3},
    };

    // Determine number of PDF members from the first input file.
    int n_members = 0;
    {
        std::string fname = pawres + "ggdrhic" + pdf_tag + "_chunkA_" + scale + ".root";
        TFile *fprobe = new TFile(fname.c_str(), "READ");
        if (!fprobe || fprobe->IsZombie()) {
            std::cerr << "Cannot open " << fname << " to determine n_members" << std::endl;
            return;
        }
        TTree *tprobe = (TTree *)fprobe->Get("t2");
        Int_t nb_member;
        tprobe->SetBranchAddress("nb_member", &nb_member);
        tprobe->GetEntry(0);
        n_members = nb_member;
        fprobe->Close();
        delete fprobe;
    }
    std::cout << "PDF set " << pdf_tag << " has " << n_members << " members" << std::endl;

    // Allocate one histogram per member.
    std::vector<TH1F*> h_members;
    h_members.reserve(n_members);
    for (int m = 0; m < n_members; ++m) {
        TH1F *h = new TH1F(Form("h_truth_pT_m%d", m),
                           Form("truth pT (dir+frag), PDF member %d", m),
                           n_pT_bins, pT_bin_edges);
        h->Sumw2();
        h_members.push_back(h);
    }

    auto fillFromFile = [&](const std::string &fname, float /*nseg_hint*/)
    {
        TFile *f = new TFile(fname.c_str(), "READ");
        if (!f || f->IsZombie()) { std::cerr << "Cannot open " << fname << std::endl; return; }

        TTree *t2 = (TTree *)f->Get("t2");
        TList *list = t2->GetUserInfo();
        TVectorT<float> *v = static_cast<TVectorT<float> *>(list->At(0));
        // Auto-detect nseg from tree entries / per-segment nb_evt. hadd keeps
        // the first file's UserInfo (per-segment), so this is robust to
        // post-hoc additions of segments to a chunk's condorout dir.
        float nseg = (float)t2->GetEntries() / (*v)[0];
        float norma = (*v)[1] / (*v)[0] / nseg;

        std::cout << "  -> " << fname << " : entries=" << t2->GetEntries() << " nb_evt=" << (*v)[0] << " nseg=" << nseg << " norma=" << norma << std::endl;

        Int_t iprov, ntrack, nb_member;
        Double_t e[3], px[3], py[3], pz[3], x3;
        Float_t pdf_weight[1000];

        t2->SetBranchAddress("iprov",      &iprov);
        t2->SetBranchAddress("ntrack",     &ntrack);
        t2->SetBranchAddress("nb_member",  &nb_member);
        t2->SetBranchAddress("x3",         &x3);
        t2->SetBranchAddress("energy",     e);
        t2->SetBranchAddress("px",         px);
        t2->SetBranchAddress("py",         py);
        t2->SetBranchAddress("pz",         pz);
        t2->SetBranchAddress("pdf_weight", pdf_weight);

        // Per-file accumulators (one per member) so we can apply the
        // file-specific norma at the end.
        std::vector<TH1F*> htmp_members;
        htmp_members.reserve(n_members);
        for (int m = 0; m < n_members; ++m) {
            TH1F *htmp = new TH1F(Form("htmp_m%d_%p", m, (void*)f), "",
                                  n_pT_bins, pT_bin_edges);
            htmp->Sumw2();
            htmp_members.push_back(htmp);
        }

        Int_t entries = (Int_t)t2->GetEntries();
        for (Int_t i = 0; i < entries; i++)
        {
            if (i % 500000 == 0)
                std::cout << "      " << i << " / " << entries << std::endl;
            t2->GetEntry(i);

            double pt0 = sqrt(px[0]*px[0] + py[0]*py[0]);
            double y0  = log((e[0]+pz[0]) / (e[0]-pz[0])) * 0.5;
            if (std::abs(y0) >= 0.7) continue;

            // Loop over all PDF members in this event.
            const int nm_evt = (nb_member < n_members) ? nb_member : n_members;
            for (int m = 0; m < nm_evt; ++m) {
                htmp_members[m]->Fill(pt0, pdf_weight[m]);
            }
        }

        for (int m = 0; m < n_members; ++m) {
            htmp_members[m]->Scale(norma);
            for (int ibin = 1; ibin <= htmp_members[m]->GetNbinsX(); ibin++) {
                float bw = htmp_members[m]->GetBinWidth(ibin);
                htmp_members[m]->SetBinContent(ibin, htmp_members[m]->GetBinContent(ibin) / bw);
                htmp_members[m]->SetBinError  (ibin, htmp_members[m]->GetBinError  (ibin) / bw);
            }
            h_members[m]->Add(htmp_members[m]);
            delete htmp_members[m];
        }

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

    std::cout << "Combined integral (member 0) over all chunks: "
              << h_members[0]->Integral("width") << std::endl;

    std::string outputname = "rootFiles/jetPHOX" + pdf_tag + "_" + scale + "_pdfmem.root";
    TFile *fout = new TFile(outputname.c_str(), "RECREATE");
    for (int m = 0; m < n_members; ++m) {
        h_members[m]->Write();
    }
    fout->Close();
    std::cout << "Written to " << outputname << " with " << n_members << " members" << std::endl;
}
