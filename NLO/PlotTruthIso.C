#include <iostream>
#include <string>
#include <array>
#include <cmath>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>
#include <TVectorT.h>
#include <TList.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <TAxis.h>
#include <TStyle.h>
#include "sPhenixStyle.C"

void PlotTruthIso()
{
    SetsPhenixStyle();
    /* -------------------------------------------------- configuration */
    const double isoVals[] = {0.5, 1, 2, 3, 4, 200};
    const char *isoFilesDir[] = {
        "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/ggdrhic_10_0_5iso.root",
        "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/ggdrhic_10_1iso.root",
        "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/ggdrhic_10_2iso.root",
        "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/ggdrhic_10_3iso.root",
        "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/ggdrhic_10_4iso.root",
        "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/ggdrhic_10_200iso.root"};

    // Matching fragmentation files (ggorhic*.root)
    const char *isoFilesFrag[] = {
        "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/ggorhic_10_0_5iso.root",
        "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/ggorhic_10_1iso.root",
        "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/ggorhic_10_2iso.root",
        "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/ggorhic_10_3iso.root",
        "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/ggorhic_10_4iso.root",
        "/sphenix/user/shuhangli/jetphox/newbuild/jetphox_1.3.1_4/pawres/ggorhic_10_200iso.root"};

    static_assert(std::size(isoVals) == std::size(isoFilesDir),  "isoVals/Dir mismatch");
    static_assert(std::size(isoVals) == std::size(isoFilesFrag), "isoVals/Frag mismatch");
    constexpr int nIso = std::size(isoVals);

    const double ptLow[]  = {10, 15, 20};
    const double ptHigh[] = {15, 20, 30};
    constexpr int nPtBins = std::size(ptLow);

    const double etaMax = 0.7;
    const double nseg   = 10.; // JetPhox segmentation

    /* -------------- helper: read one ROOT file -> yields per pT bin -- */
    auto processFileToYields = [&](const char *fname)
        -> std::array<double, nPtBins>
    {
        TFile f(fname, "READ");
        if (f.IsZombie()) throw std::runtime_error(std::string("Cannot open ")+fname);

        TTree *t2 = (TTree*)f.Get("t2");
        if (!t2) throw std::runtime_error(std::string("Missing TTree t2 in ")+fname);

        TList *lst = t2->GetUserInfo();
        if (!lst || lst->GetEntries()<1) throw std::runtime_error("No UserInfo / xsec block");
        auto *v = dynamic_cast<TVectorT<float>*>(lst->At(0));
        if (!v || v->GetNrows()<2) throw std::runtime_error("Bad xsec vector");
        const double nb_evt = (*v)[0];
        const double xsec   = (*v)[1];
        const double norma  = xsec / nb_evt / nseg;

        Int_t    iprov=0, ntrack=0;
        Double_t e[3]{}, px[3]{}, py[3]{}, pz[3]{}, x3=0;
        Double_t pt[3]{}, y[3]{};
        Double_t x1=0, x2=0;
        Float_t  pdf_weight[1000]{};
        Float_t  weight=1.0;

        t2->SetBranchAddress("iprov",      &iprov);
        t2->SetBranchAddress("ntrack",     &ntrack);
        t2->SetBranchAddress("x3",         &x3);
        t2->SetBranchAddress("energy",     e);
        t2->SetBranchAddress("px",         px);
        t2->SetBranchAddress("py",         py);
        t2->SetBranchAddress("pz",         pz);
        t2->SetBranchAddress("pdf_weight", pdf_weight);

        TH1F h("tmp", ";p_{T};d#sigma/dp_{T}", 33, 7, 40);
        Int_t entries = (Int_t)t2->GetEntries();
        for (Int_t i = 0; i < entries; i++)
        {
            t2->GetEntry(i);
            for (Int_t j = 0; j < ntrack; j++) {
                pt[j] = std::sqrt(px[j]*px[j] + py[j]*py[j]);
                y[j]  = 0.5*std::log((e[j]+pz[j])/(e[j]-pz[j]));
            }
            weight = pdf_weight[0];
            if (std::abs(y[0]) < etaMax) {
                h.Fill(pt[0], weight);
            }
        }

        // normalize to per-GeV
        for (int b = 1; b <= h.GetNbinsX(); ++b) {
            const double bw = h.GetBinWidth(b);
            h.SetBinContent(b, h.GetBinContent(b)/bw);
            h.SetBinError  (b, h.GetBinError  (b)/bw);
        }
        h.Scale(norma);

        std::array<double, nPtBins> yld{};
        for (int ip = 0; ip < nPtBins; ++ip) {
            const int binLo = h.FindBin(ptLow[ip]  + 1e-4);
            const int binHi = h.FindBin(ptHigh[ip] - 1e-4);
            yld[ip] = h.Integral(binLo, binHi, "width");
        }
        return yld;
    };

    /* ---------------------------------------------- compute yields --- */
    std::array<std::array<double, nPtBins>, nIso> yDir{};
    std::array<std::array<double, nPtBins>, nIso> yFrag{};

    for (int i = 0; i < nIso; ++i) {
        yDir [i] = processFileToYields(isoFilesDir [i]);
        yFrag[i] = processFileToYields(isoFilesFrag[i]);
    }

    /* ---------------------------------------------- compute fractions */
    std::array<std::array<double, nPtBins>, nIso> fDir{};
    std::array<std::array<double, nPtBins>, nIso> fFrag{};
    for (int ii = 0; ii < nIso; ++ii) {
        for (int ip = 0; ip < nPtBins; ++ip) {
            const double sum = yDir[ii][ip] + yFrag[ii][ip];
            if (sum > 0) {
                fDir [ii][ip] = yDir [ii][ip] / sum;
                fFrag[ii][ip] = yFrag[ii][ip] / sum;
            } else {
                fDir [ii][ip] = 0.0;
                fFrag[ii][ip] = 0.0;
            }
        }
    }

    /* ------------------------------------------- plot vs iso_Et (ET bins) */
    std::array<TGraph*, nPtBins> grDirIso{};
    std::array<TGraph*, nPtBins> grFragIso{};

    for (int ip = 0; ip < nPtBins; ++ip) {
        grDirIso [ip] = new TGraph(nIso);
        grFragIso[ip] = new TGraph(nIso);
        for (int ii = 0; ii < nIso; ++ii) {
            grDirIso [ip]->SetPoint(ii, isoVals[ii], fDir [ii][ip]);
            grFragIso[ip]->SetPoint(ii, isoVals[ii], fFrag[ii][ip]);
        }
        const int color = 1 + ip;
        grDirIso [ip]->SetMarkerStyle(20 + ip);
        grDirIso [ip]->SetMarkerColor(color);
        grDirIso [ip]->SetLineColor  (color);
        grDirIso [ip]->SetLineWidth  (2);

        grFragIso[ip]->SetMarkerStyle(24 + ip);
        grFragIso[ip]->SetMarkerColor(color);
        grFragIso[ip]->SetLineColor  (color);
        grFragIso[ip]->SetLineStyle  (kDashed);
        grFragIso[ip]->SetLineWidth  (2);
    }

    TCanvas *c1 = new TCanvas("cFracVsIso", "Direct/Frag Fractions vs iso", 900, 650);
    c1->SetGridy();
    TH1F *frame1 = c1->DrawFrame(
        0.0, 0.0,   // xmin, ymin
        5.0, 1.05,  // xmax (zoom on 0–5), ymax
        ";Isolation cut (GeV);Fraction"
    );
    frame1->GetYaxis()->SetNdivisions(505);

    TLegend *leg1 = new TLegend(0.48, 0.15, 0.88, 0.42);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetHeader("#mu=E_{T}^{#gamma}", "C");

    for (int ip = 0; ip < nPtBins; ++ip) {
        grDirIso [ip]->Draw("PL SAME");
        grFragIso[ip]->Draw("PL SAME");
        leg1->AddEntry(grDirIso[ip],  Form("%.0f < E_{T}^{#gamma} < %.0f GeV (dir)",  ptLow[ip], ptHigh[ip]), "pl");
        leg1->AddEntry(grFragIso[ip], Form("%.0f < E_{T}^{#gamma} < %.0f GeV (frag)", ptLow[ip], ptHigh[ip]), "pl");
    }
    leg1->Draw();
    // c1->SaveAs("frac_vs_iso.pdf");

    /* -------------------------------------- plot vs photon ET (iso bins) */
    // x = ET bin center, y = f_dir (direct fraction)
    std::array<TGraph*, nIso> grDirVsPt{};
    double ptCenter[nPtBins];
    for (int ip = 0; ip < nPtBins; ++ip) ptCenter[ip] = 0.5*(ptLow[ip] + ptHigh[ip]);

    for (int ii = 0; ii < nIso; ++ii) {
        grDirVsPt[ii] = new TGraph(nPtBins);
        for (int ip = 0; ip < nPtBins; ++ip) {
            grDirVsPt[ii]->SetPoint(ip, ptCenter[ip], fDir[ii][ip]);
        }
        const int color = 1 + (ii % 6);
        grDirVsPt[ii]->SetMarkerStyle(20 + (ii % 6));
        grDirVsPt[ii]->SetMarkerColor(color);
        grDirVsPt[ii]->SetLineColor  (color);
        grDirVsPt[ii]->SetLineWidth  (2);
    }

    TCanvas *c2 = new TCanvas("cFracVsPt", "Direct Fraction vs photon ET", 900, 650);
    c2->SetGridy();
    TH1F *frame2 = c2->DrawFrame(
        ptLow[0]*0.9, 0.0,
        ptHigh[nPtBins-1]*1.1, 1.05,
        ";E_{T}^{#gamma} (GeV);Direct fraction"
    );
    frame2->GetYaxis()->SetNdivisions(505);

    TLegend *leg2 = new TLegend(0.52, 0.15, 0.88, 0.42);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetHeader("Curves by iso cut", "C");

    for (int ii = 0; ii < nIso; ++ii) {
        grDirVsPt[ii]->Draw("PL SAME");
        leg2->AddEntry(grDirVsPt[ii], Form("iso = %.3g", isoVals[ii]), "pl");
    }
    leg2->Draw();
    // c2->SaveAs("direct_frac_vs_pt.pdf");

    /* ------------------------------ print table to stdout (optional) -- */
    std::cout << "\nFractions f_dir (rows=iso, cols=ET bins):\n";
    for (int ii = 0; ii < nIso; ++ii) {
        std::cout << "iso=" << isoVals[ii] << " : ";
        for (int ip = 0; ip < nPtBins; ++ip)
            std::cout << fDir[ii][ip] << "  ";
        std::cout << "\n";
    }
}
