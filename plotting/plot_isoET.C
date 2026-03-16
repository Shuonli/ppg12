#include "plotcommon.h"

// plot_isoET.C
//
// Overlays isolation ET distributions for data vs inclusive jet MC,
// one canvas per pT bin.  Two output PDFs per pT bin:
//
//   iso_ET_tight_pt<ipt>.pdf   – tight selection,    data vs jet MC + Data/MC ratio
//   iso_ET_nontight_pt<ipt>.pdf – non-tight selection, data vs jet MC + Data/MC ratio
//
// Input files (produced by RecoEffCalculator_TTreeReader.C + hadd):
//   results/data_histo_<suffix>.root        – data
//   results/MC_efficiency_jet_<suffix>.root – inclusive jet MC (hadd of jet10–50)
//
// Histograms read:
//   h_tight_isoET_<etabin>_<ipt>     (TH1D, isoET in the tight shower-shape region)
//   h_nontight_isoET_<etabin>_<ipt>  (TH1D, isoET in the non-tight region)
//
// Usage:
//   root -l -q 'plot_isoET.C("bdt_base_v3E_2")'

void plot_isoET(const std::string suffix  = "bdt_emcalinnerr0_isoscale08",
                const int         etabin  = 0,
                const int         rebin   = 2)
{
    init_plot();

    const std::string resdir = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    const std::string figdir = "figures";

    // -----------------------------------------------------------------------
    // Open files
    // -----------------------------------------------------------------------
    TFile *fdata = TFile::Open(Form("%s/data_histo_%s.root",            resdir.c_str(), suffix.c_str()));
    TFile *fjet  = TFile::Open(Form("%s/MC_efficiency_jet_%s.root",     resdir.c_str(), suffix.c_str()));

    if (!fdata || fdata->IsZombie()) { std::cerr << "Cannot open data file\n";    return; }
    if (!fjet  || fjet->IsZombie())  { std::cerr << "Cannot open jet MC file\n";  return; }

    // -----------------------------------------------------------------------
    // Colours
    // -----------------------------------------------------------------------
    const int col_data  = kBlack;
    const int col_mc    = kAzure + 2;

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    // Normalise clone to probability density (unit area / bin width).
    auto normToWidth = [](TH1D *h, const char *name) -> TH1D* {
        TH1D *hc = dynamic_cast<TH1D*>(h->Clone(name));
        double integ = hc->Integral("width");
        if (integ > 0) hc->Scale(1.0 / integ);
        return hc;
    };

    // Two-pad split canvas: top 70%, bottom 30%.
    auto makeSplitCanvas = [](const char *name) -> TCanvas* {
        TCanvas *c = new TCanvas(name, name, 600, 700);
        TPad *p1 = new TPad("p1", "top",    0, 0.30, 1, 1.00);
        p1->SetBottomMargin(0.02); p1->SetLeftMargin(0.15); p1->Draw();
        TPad *p2 = new TPad("p2", "bottom", 0, 0.00, 1, 0.30);
        p2->SetTopMargin(0.04);   p2->SetBottomMargin(0.35);
        p2->SetLeftMargin(0.15);  p2->Draw();
        return c;
    };

    // Style: data (points), MC (line).
    auto styleData = [&](TH1D *h) {
        h->SetMarkerStyle(20); h->SetMarkerSize(0.9);
        h->SetMarkerColor(col_data); h->SetLineColor(col_data);
    };
    auto styleMC = [&](TH1D *h) {
        h->SetLineColor(col_mc); h->SetLineWidth(2); h->SetMarkerSize(0);
    };

    // Fit exponential [0]*exp(-x/[1]) to h from isoET>0; draws result on current pad.
    // Returns the fitted TF1 (tau = parameter [1]).
    auto fitExp = [](TH1D *h, int col, const char *fname) -> TF1* {
        TF1 *f = new TF1(fname, "[0]*exp(-x/[1])", 1.0, 15.0);
        f->SetParameter(0, h->GetMaximum());
        f->SetParameter(1, 2.0);
        f->SetLineColor(col); f->SetLineWidth(2); f->SetLineStyle(2);
        h->Fit(f, "REMQ");   // fit range, quiet, no auto-draw
        f->Draw("same");
        return f;
    };

    // Build the ratio top-pad frame from frame_isoET.
    auto makeIsoFrame = [](const char *name, float ylo, float yhi,
                           bool label_x) -> TH1F* {
        TH1F *fr = dynamic_cast<TH1F*>(frame_isoET->Clone(name));
        fr->GetXaxis()->SetRangeUser(-1, 6);
        fr->GetYaxis()->SetRangeUser(ylo, yhi);
        if (!label_x) fr->GetXaxis()->SetLabelSize(0);
        return fr;
    };

    // -----------------------------------------------------------------------
    // One set of plots per pT bin
    // -----------------------------------------------------------------------
    for (int ipt = 0; ipt < NptBins; ++ipt)
    {
        // --- Load histograms ------------------------------------------------
        TH1D *htight_data = dynamic_cast<TH1D*>(fdata->Get(Form("h_tight_isoET_%d_%d",    etabin, ipt)));
        TH1D *htight_mc   = dynamic_cast<TH1D*>(fjet ->Get(Form("h_tight_isoET_%d_%d",    etabin, ipt)));
        TH1D *hnt_data    = dynamic_cast<TH1D*>(fdata->Get(Form("h_nontight_isoET_%d_%d", etabin, ipt)));
        TH1D *hnt_mc      = dynamic_cast<TH1D*>(fjet ->Get(Form("h_nontight_isoET_%d_%d", etabin, ipt)));

        if (!htight_data || !htight_mc || !hnt_data || !hnt_mc)
        {
            std::cerr << "Missing histogram for pT bin " << ipt
                      << " (eta bin " << etabin << "), skipping\n";
            continue;
        }

        // --- Rebin + normalise clones ---------------------------------------
        htight_data->Rebin(rebin); htight_mc->Rebin(rebin);
        hnt_data   ->Rebin(rebin); hnt_mc   ->Rebin(rebin);

        TH1D *hd_tight  = normToWidth(htight_data, Form("hd_tight_%d",  ipt));
        TH1D *hmc_tight = normToWidth(htight_mc,   Form("hmc_tight_%d", ipt));
        TH1D *hd_nt     = normToWidth(hnt_data,    Form("hd_nt_%d",     ipt));
        TH1D *hmc_nt    = normToWidth(hnt_mc,      Form("hmc_nt_%d",    ipt));

        styleData(hd_tight);  styleMC(hmc_tight);
        styleData(hd_nt);     styleMC(hmc_nt);

        // Ratios (data / MC)
        TH1D *hratio_tight = dynamic_cast<TH1D*>(hd_tight->Clone(Form("hratio_tight_%d", ipt)));
        hratio_tight->Divide(hmc_tight);
        TH1D *hratio_nt    = dynamic_cast<TH1D*>(hd_nt->Clone(Form("hratio_nt_%d", ipt)));
        hratio_nt->Divide(hmc_nt);

        hratio_tight->SetMarkerStyle(20); hratio_tight->SetMarkerSize(0.7);
        hratio_tight->SetMarkerColor(col_data); hratio_tight->SetLineColor(col_data);
        hratio_nt->SetMarkerStyle(20); hratio_nt->SetMarkerSize(0.7);
        hratio_nt->SetMarkerColor(col_data); hratio_nt->SetLineColor(col_data);

        // pT bin label
        const char *ptlabel = Form("%.0f < #it{E}_{T}^{#gamma,rec} < %.0f GeV",
                                    ptRanges[ipt], ptRanges[ipt + 1]);

        // ===================================================================
        // Plot A: tight selection
        // ===================================================================
        {
            TCanvas *c = makeSplitCanvas(Form("c_tight_isoET_%d", ipt));

            // -- top pad --
            TPad *p1 = dynamic_cast<TPad*>(c->FindObject("p1"));
            p1->cd(); 
            //p1->SetLogy();

            double ymax = std::max(hd_tight->GetMaximum(), hmc_tight->GetMaximum());
            TH1F *fr = makeIsoFrame(Form("fr_tight_top_%d", ipt), 1e-4, ymax * 1.2, false);
            fr->SetYTitle("Probability density / GeV");
            fr->GetYaxis()->SetTitleSize(0.048);
            fr->GetYaxis()->SetTitleOffset(1.5);
            fr->Draw("axis");

            hmc_tight->Draw("same hist");
            hd_tight ->Draw("same ex0");

            //TF1 *fexp_t = fitExp(hd_tight, kRed+1, Form("fexp_tight_%d", ipt));
            //double tau_t  = fexp_t->GetParameter(1);
            //double tauE_t = fexp_t->GetParError(1);

            myText(0.18, 0.88, 1, strleg1.c_str(), 0.042);
            myText(0.18, 0.83, 1, strleg2.c_str(), 0.042);
            myText(0.18, 0.78, 1, strleg3.c_str(), 0.042);
            myText(0.18, 0.73, 1, ptlabel,          0.040);

            TLegend *leg = new TLegend(0.50, 0.52, 0.90, 0.73);
            legStyle(leg, 0.17, 0.040);
            leg->SetHeader("Tight #it{ID}");
            leg->AddEntry(hd_tight,  "Data",               "pl");
            leg->AddEntry(hmc_tight, strIncMC.c_str(),     "l");
            //leg->AddEntry(fexp_t,    Form("Exp. fit: #tau = %.2f #pm %.2f GeV", tau_t, tauE_t), "l");
            leg->Draw();

            // -- bottom pad (ratio) --
            c->cd();
            TPad *p2 = dynamic_cast<TPad*>(c->FindObject("p2"));
            p2->cd();
            TH1F *fr2 = makeIsoFrame(Form("fr_tight_bot_%d", ipt), 0.0, 2.5, true);
            fr2->SetYTitle("Data / MC");
            fr2->GetYaxis()->SetNdivisions(504);
            fr2->GetYaxis()->SetTitleSize(0.12);  fr2->GetYaxis()->SetTitleOffset(0.55);
            fr2->GetYaxis()->SetLabelSize(0.10);
            fr2->GetXaxis()->SetTitleSize(0.13);  fr2->GetXaxis()->SetLabelSize(0.11);
            fr2->Draw("axis");
            hratio_tight->Draw("same ex0");
            lineone->Draw("same");

            c->SaveAs(Form("%s/iso_ET_tight_pt%d_%s.pdf", figdir.c_str(), ipt, suffix.c_str()));
            delete c;
        }

        // ===================================================================
        // Plot B: non-tight selection
        // ===================================================================
        {
            TCanvas *c = makeSplitCanvas(Form("c_nontight_isoET_%d", ipt));

            // -- top pad --
            TPad *p1 = dynamic_cast<TPad*>(c->FindObject("p1"));
            p1->cd(); 
            //p1->SetLogy();

            double ymax = std::max(hd_nt->GetMaximum(), hmc_nt->GetMaximum());
            TH1F *fr = makeIsoFrame(Form("fr_nt_top_%d", ipt), 1e-4, ymax * 1.2, false);
            fr->SetYTitle("Probability density / GeV");
            fr->GetYaxis()->SetTitleSize(0.048);
            fr->GetYaxis()->SetTitleOffset(1.5);
            fr->Draw("axis");

            hmc_nt->Draw("same hist");
            hd_nt ->Draw("same ex0");

            //TF1 *fexp_nt = fitExp(hd_nt, kRed+1, Form("fexp_nt_%d", ipt));
            //double tau_nt  = fexp_nt->GetParameter(1);
            //double tauE_nt = fexp_nt->GetParError(1);

            myText(0.18, 0.88, 1, strleg1.c_str(), 0.042);
            myText(0.18, 0.83, 1, strleg2.c_str(), 0.042);
            myText(0.18, 0.78, 1, strleg3.c_str(), 0.042);
            myText(0.18, 0.73, 1, ptlabel,          0.040);

            TLegend *leg = new TLegend(0.50, 0.52, 0.90, 0.73);
            legStyle(leg, 0.17, 0.040);
            leg->SetHeader("Non-tight #it{ID}");
            leg->AddEntry(hd_nt,  "Data",               "pl");
            leg->AddEntry(hmc_nt, strIncMC.c_str(),     "l");
            //leg->AddEntry(fexp_nt, Form("Exp. fit: #tau = %.2f #pm %.2f GeV", tau_nt, tauE_nt), "l");
            leg->Draw();

            // -- bottom pad (ratio) --
            c->cd();
            TPad *p2 = dynamic_cast<TPad*>(c->FindObject("p2"));
            p2->cd();
            TH1F *fr2 = makeIsoFrame(Form("fr_nt_bot_%d", ipt), 0.0, 2.5, true);
            fr2->SetYTitle("Data / MC");
            fr2->GetYaxis()->SetNdivisions(504);
            fr2->GetYaxis()->SetTitleSize(0.12);  fr2->GetYaxis()->SetTitleOffset(0.55);
            fr2->GetYaxis()->SetLabelSize(0.10);
            fr2->GetXaxis()->SetTitleSize(0.13);  fr2->GetXaxis()->SetLabelSize(0.11);
            fr2->Draw("axis");
            hratio_nt->Draw("same ex0");
            lineone->Draw("same");

            c->SaveAs(Form("%s/iso_ET_nontight_pt%d_%s.pdf", figdir.c_str(), ipt, suffix.c_str()));
            delete c;
        }

        // ===================================================================
        // Plot C: tight + non-tight overlaid (data vs MC, no ratio)
        // ===================================================================
        {
            const int col_nt_data = kRed - 4;
            const int col_nt_mc   = kOrange + 2;

            hd_nt ->SetMarkerColor(col_nt_data); hd_nt ->SetLineColor(col_nt_data);
            hmc_nt->SetLineColor(col_nt_mc);

            TCanvas *c = new TCanvas(Form("c_comb_isoET_%d", ipt),
                                     Form("c_comb_isoET_%d", ipt), 600, 600);
            c->SetLeftMargin(0.15); c->SetLogy();

            double ymax = std::max({hd_tight->GetMaximum(), hmc_tight->GetMaximum(),
                                    hd_nt->GetMaximum(),    hmc_nt->GetMaximum()});
            TH1F *fr = makeIsoFrame(Form("fr_comb_%d", ipt), 1e-4, ymax * 8, true);
            fr->SetYTitle("Probability density / GeV");
            fr->GetYaxis()->SetTitleSize(0.048);
            fr->GetYaxis()->SetTitleOffset(1.5);
            fr->Draw("axis");

            hmc_tight->Draw("same hist");
            hd_tight ->Draw("same ex0");
            hmc_nt   ->Draw("same hist");
            hd_nt    ->Draw("same ex0");

            TF1 *fexp_ct = fitExp(hd_tight, kRed+1,    Form("fexp_ct_%d",  ipt));
            TF1 *fexp_cn = fitExp(hd_nt,    kOrange+3,  Form("fexp_cn_%d", ipt));
            double tau_ct  = fexp_ct->GetParameter(1), tauE_ct  = fexp_ct->GetParError(1);
            double tau_cn  = fexp_cn->GetParameter(1), tauE_cn  = fexp_cn->GetParError(1);

            myText(0.18, 0.88, 1, strleg1.c_str(), 0.042);
            myText(0.18, 0.83, 1, strleg2.c_str(), 0.042);
            myText(0.18, 0.78, 1, strleg3.c_str(), 0.042);
            myText(0.18, 0.73, 1, ptlabel,          0.038);

            TLegend *leg = new TLegend(0.36, 0.40, 0.90, 0.72);
            legStyle(leg, 0.17, 0.036);
            leg->AddEntry(hd_tight,  "Data (tight)",                               "pl");
            leg->AddEntry(hmc_tight, Form("%s (tight)", strIncMC.c_str()),          "l");
            leg->AddEntry(fexp_ct,   Form("Exp. fit (tight): #tau = %.2f #pm %.2f GeV", tau_ct, tauE_ct), "l");
            leg->AddEntry(hd_nt,     "Data (non-tight)",                            "pl");
            leg->AddEntry(hmc_nt,    Form("%s (non-tight)", strIncMC.c_str()),      "l");
            leg->AddEntry(fexp_cn,   Form("Exp. fit (non-tight): #tau = %.2f #pm %.2f GeV", tau_cn, tauE_cn), "l");
            leg->Draw();

            c->SaveAs(Form("%s/iso_ET_comb_pt%d_%s.pdf", figdir.c_str(), ipt, suffix.c_str()));
            delete c;
        }
    }

    std::cout << "Saved plots to " << figdir << "/iso_ET_*.pdf\n";
}
