#include "plotcommon.h"

// Overlay 70/80/90% iso-efficiency cuts vs cluster pT and fit each
// with a linear function. Reads h_singal_reco_isoET_0 from
// MC_efficiency_bdt_<var>.root and inverts the truth-iso cumulative
// distribution per pT bin.
//
// Output: figures/ETCut_FitResults.pdf (the figure referenced by
// reconstruction.tex Figure~\ref{fig:isoET_cut} caption).

namespace {
float findEffCutoff(TH2D *h, float eff, float pTlow, float pThigh)
{
    int xlo = std::max(1, h->GetXaxis()->FindBin(pTlow));
    int xhi = std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(pThigh));
    double total = 0;
    for (int i = xlo; i <= xhi; ++i)
        for (int j = 1; j <= h->GetNbinsY(); ++j)
            total += h->GetBinContent(i, j);
    if (total <= 0) return 0;
    double tgt = total * eff, cum = 0;
    for (int j = 1; j <= h->GetNbinsY(); ++j) {
        for (int i = xlo; i <= xhi; ++i) cum += h->GetBinContent(i, j);
        if (cum >= tgt) return h->GetYaxis()->GetBinUpEdge(j);
    }
    return h->GetYaxis()->GetBinUpEdge(h->GetNbinsY());
}
}

void plot_etcut_overlay(const std::string &var = "nom")
{
    init_plot();

    TFile *fin = TFile::Open(Form(
        "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_%s.root",
        var.c_str()), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open MC_efficiency_bdt_" << var << ".root\n";
        return;
    }
    TH2D *h2 = (TH2D *)fin->Get("h_singal_reco_isoET_0");
    if (!h2) { std::cerr << "ERROR: h_singal_reco_isoET_0 missing\n"; return; }

    const int nbins = 13;
    const double xlo = 10, xhi = 36;
    std::vector<double> effs   = {0.70, 0.80, 0.90};
    std::vector<int>    cols   = {kAzure + 7, kGreen + 2, kPink + 8};
    std::vector<int>    marks  = {21, 22, 20};
    std::vector<TH1F*>  hists;
    std::vector<TF1*>   fits;
    std::vector<std::pair<double,double>> pars;  // (intercept, slope)

    double ywid = h2->GetYaxis()->GetBinWidth(1);
    for (size_t k = 0; k < effs.size(); ++k) {
        TH1F *h = new TH1F(Form("h_eff_%zu", k), "", nbins, xlo, xhi);
        for (int i = 1; i <= nbins; ++i) {
            double pTlo = h->GetXaxis()->GetBinLowEdge(i);
            double pThi = h->GetXaxis()->GetBinUpEdge(i);
            h->SetBinContent(i, findEffCutoff(h2, effs[k], pTlo, pThi));
            h->SetBinError(i, ywid / 2.0);
        }
        TF1 *f = new TF1(Form("fit_%zu", k), "pol1", xlo, xhi);
        f->SetLineColor(cols[k]);
        f->SetLineStyle(2);
        h->Fit(f, "RQN");
        hists.push_back(h);
        fits.push_back(f);
        pars.push_back({f->GetParameter(0), f->GetParameter(1)});
    }

    TCanvas *c = new TCanvas("c_etcut", "", 700, 600);
    frame_et_rec->SetTitle(";Cluster #it{p}_{T} [GeV];#it{E}_{T}^{iso} cutoff [GeV]");
    frame_et_rec->GetXaxis()->SetRangeUser(xlo, xhi);
    frame_et_rec->GetYaxis()->SetRangeUser(0, 4.5);
    frame_et_rec->Draw("axis");

    for (size_t k = 0; k < effs.size(); ++k) {
        hists[k]->SetMarkerStyle(marks[k]);
        hists[k]->SetMarkerColor(cols[k]);
        hists[k]->SetMarkerSize(1.4);
        hists[k]->SetLineColor(cols[k]);
        hists[k]->Draw("P same");
        fits[k]->Draw("same");
    }

    myText(0.20, 0.88, 1, strleg1.c_str(), 0.040, 0);
    myText(0.20, 0.83, 1, strleg2.c_str(), 0.040, 0);
    myText(0.20, 0.78, 1, strleg3.c_str(), 0.040, 0);
    myText(0.88, 0.88, 1, strSigMC.c_str(), 0.040, 1);

    TLegend *leg = new TLegend(0.20, 0.55, 0.55, 0.74);
    legStyle(leg, 0.20, 0.040);
    for (size_t k = 0; k < effs.size(); ++k) {
        leg->AddEntry(hists[k],
            Form("%d%%: %.3f + %.4f #it{p}_{T}",
                 (int)(100*effs[k]), pars[k].first, pars[k].second), "pl");
    }
    leg->Draw("same");

    c->SaveAs("figures/ETCut_FitResults.pdf");
    fin->Close();
}
