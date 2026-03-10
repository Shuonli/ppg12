// PlotTruthJetWindows.C
// Reads the output of TruthJetInETWindow.C and produces:
//   1. Multi-page PDF: max truth jet pT per ET window (sig vs bkg vs all)
//   2. Multi-page PDF: all truth jet pT per ET window (sig vs bkg vs all)
//   3. Single PDF: 2D cluster ET vs max truth jet pT
//
// Usage:
//   root -l -b -q 'PlotTruthJetWindows.C("results/truthjet_ETwindows_all.root","truthjet_plots")'

#include "../plotting/plotcommon.h"

void PlotTruthJetWindows(
    const std::string &infile  = "results/truthjet_ETwindows_all.root",
    const std::string &outdir  = "truthjet_plots")
{
    init_plot();
    gSystem->mkdir(outdir.c_str(), kTRUE);

    TFile *fin = TFile::Open(infile.c_str(), "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << infile << std::endl;
        return;
    }

    // Read window edges back from stored histogram
    TH1D *h_edges = dynamic_cast<TH1D *>(fin->Get("h_ET_edges"));
    if (!h_edges)
    {
        std::cerr << "ERROR: h_ET_edges not found in " << infile << std::endl;
        return;
    }
    const int N_win = h_edges->GetNbinsX();
    std::vector<float> ET_edges(N_win + 1);
    for (int i = 0; i <= N_win; i++)
        ET_edges[i] = h_edges->GetXaxis()->GetBinLowEdge(i + 1);

    // ── canvas style helpers ─────────────────────────────────────────────────
    const float tsize     = 0.038;
    const float tsize_leg = 0.034;

    auto MakeCanvas = [&](const char *name) -> TCanvas * {
        TCanvas *c = new TCanvas(name, name, 600, 600);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.05);
        c->SetTopMargin(0.07);
        c->SetBottomMargin(0.14);
        c->SetGrid();
        return c;
    };

    auto StyleHist = [&](TH1D *h, int color, int style) {
        h->SetLineColor(color);
        h->SetLineStyle(style);
        h->SetLineWidth(2);
        h->GetXaxis()->SetTitleSize(tsize);
        h->GetYaxis()->SetTitleSize(tsize);
        h->GetXaxis()->SetLabelSize(tsize);
        h->GetYaxis()->SetLabelSize(tsize);
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetYaxis()->SetTitleOffset(1.6);
        h->GetXaxis()->CenterTitle();
        h->GetYaxis()->CenterTitle();
    };

    // ── 1 & 2: per-window 1D plots ───────────────────────────────────────────
    for (int iw = 0; iw < N_win; iw++)
    {
        float lo = ET_edges[iw], hi = ET_edges[iw + 1];
        const char *tag    = Form("win%d", iw);
        const char *strWin = Form("%.0f < #it{E}_{T}^{cluster} < %.0f GeV", lo, hi);
        const char *wintag = Form("%.0f_%.0f", lo, hi);

        TH1D *hmax_sig   = dynamic_cast<TH1D *>(fin->Get(Form("h_maxtruthjet_pt_sig_%s",   tag)));
        TH1D *hmax_bkg   = dynamic_cast<TH1D *>(fin->Get(Form("h_maxtruthjet_pt_bkg_%s",   tag)));
        TH1D *hmax_all   = dynamic_cast<TH1D *>(fin->Get(Form("h_maxtruthjet_pt_all_%s",   tag)));
        TH1D *hmax_tight     = dynamic_cast<TH1D *>(fin->Get(Form("h_maxtruthjet_pt_tight_%s",     tag)));
        TH1D *hmax_tightonly = dynamic_cast<TH1D *>(fin->Get(Form("h_maxtruthjet_pt_tightonly_%s", tag)));
        TH1D *hall_sig       = dynamic_cast<TH1D *>(fin->Get(Form("h_alltruthjet_pt_sig_%s",       tag)));
        TH1D *hall_bkg       = dynamic_cast<TH1D *>(fin->Get(Form("h_alltruthjet_pt_bkg_%s",       tag)));
        TH1D *hall_all       = dynamic_cast<TH1D *>(fin->Get(Form("h_alltruthjet_pt_all_%s",       tag)));
        TH1D *hall_tight     = dynamic_cast<TH1D *>(fin->Get(Form("h_alltruthjet_pt_tight_%s",     tag)));
        TH1D *hall_tightonly = dynamic_cast<TH1D *>(fin->Get(Form("h_alltruthjet_pt_tightonly_%s", tag)));

        if (!hmax_all || !hmax_sig || !hmax_bkg || !hall_all || !hall_sig || !hall_bkg)
        {
            std::cerr << "WARNING: missing histograms for window " << iw << ", skipping\n";
            continue;
        }

        // ── max truth jet pT ─────────────────────────────────────────────────
        {
            TCanvas *c = MakeCanvas(Form("c_maxjet_%s", wintag));
            c->SetLogy();

            StyleHist(hmax_all, kBlack,     1);
            StyleHist(hmax_sig, kBlue,      1);
            StyleHist(hmax_bkg, kRed,       2);
            if (hmax_tightonly) StyleHist(hmax_tightonly, kOrange + 7, 7);
            if (hmax_tight)     StyleHist(hmax_tight,     kGreen + 2,  9);

            double ymax = std::max({hmax_sig->GetMaximum(),
                                    hmax_bkg->GetMaximum(),
                                    hmax_all->GetMaximum()}) * 8.0;
            if (hmax_tightonly) ymax = std::max(ymax, hmax_tightonly->GetMaximum() * 8.0);
            if (hmax_tight)     ymax = std::max(ymax, hmax_tight->GetMaximum() * 8.0);
            double ymin = std::max(1e-4 * ymax, 1e-3);
            hmax_all->GetYaxis()->SetRangeUser(ymin, ymax);
            hmax_all->Draw("hist");
            hmax_sig->Draw("hist same");
            hmax_bkg->Draw("hist same");
            if (hmax_tightonly) hmax_tightonly->Draw("hist same");
            if (hmax_tight)     hmax_tight->Draw("hist same");

            TLegend *leg = new TLegend(0.38, 0.60, 0.93, 0.90);
            legStyle(leg, 0.12, tsize_leg, 42, strWin);
            leg->AddEntry(hmax_all,       "All clusters",         "L");
            leg->AddEntry(hmax_sig,       "Signal clusters",      "L");
            leg->AddEntry(hmax_bkg,       "Bkg clusters",         "L");
            if (hmax_tightonly) leg->AddEntry(hmax_tightonly, "Tight only clusters",    "L");
            if (hmax_tight)     leg->AddEntry(hmax_tight,     "Tight+iso clusters",     "L");
            leg->Draw();

            myText(0.93, 0.93, 1, strleg1.c_str(), tsize, true);
            myText(0.93, 0.88, 1, strleg2.c_str(), tsize, true);
            myText(0.93, 0.83, 1, strleg3.c_str(), tsize, true);

            c->SaveAs(Form("%s/truthjet_maxpt_%s.pdf", outdir.c_str(), wintag));
            delete leg;
            delete c;
        }

        // ── all truth jets pT ────────────────────────────────────────────────
        {
            TCanvas *c = MakeCanvas(Form("c_alljet_%s", wintag));
            c->SetLogy();

            StyleHist(hall_all, kBlack,     1);
            StyleHist(hall_sig, kBlue,      1);
            StyleHist(hall_bkg, kRed,       2);
            if (hall_tightonly) StyleHist(hall_tightonly, kOrange + 7, 7);
            if (hall_tight)     StyleHist(hall_tight,     kGreen + 2,  9);

            double ymax = std::max({hall_sig->GetMaximum(),
                                    hall_bkg->GetMaximum(),
                                    hall_all->GetMaximum()}) * 8.0;
            if (hall_tightonly) ymax = std::max(ymax, hall_tightonly->GetMaximum() * 8.0);
            if (hall_tight)     ymax = std::max(ymax, hall_tight->GetMaximum() * 8.0);
            double ymin = std::max(1e-4 * ymax, 1e-3);
            hall_all->GetYaxis()->SetRangeUser(ymin, ymax);
            hall_all->Draw("hist");
            hall_sig->Draw("hist same");
            hall_bkg->Draw("hist same");
            if (hall_tightonly) hall_tightonly->Draw("hist same");
            if (hall_tight)     hall_tight->Draw("hist same");

            TLegend *leg = new TLegend(0.38, 0.60, 0.93, 0.90);
            legStyle(leg, 0.12, tsize_leg, 42, strWin);
            leg->AddEntry(hall_all,       "All clusters",         "L");
            leg->AddEntry(hall_sig,       "Signal clusters",      "L");
            leg->AddEntry(hall_bkg,       "Bkg clusters",         "L");
            if (hall_tightonly) leg->AddEntry(hall_tightonly, "Tight only clusters",    "L");
            if (hall_tight)     leg->AddEntry(hall_tight,     "Tight+iso clusters",     "L");
            leg->Draw();

            myText(0.93, 0.93, 1, strleg1.c_str(), tsize, true);
            myText(0.93, 0.88, 1, strleg2.c_str(), tsize, true);
            myText(0.93, 0.83, 1, strleg3.c_str(), tsize, true);

            c->SaveAs(Form("%s/truthjet_allpt_%s.pdf", outdir.c_str(), wintag));
            delete leg;
            delete c;
        }
    }

    // ── 3: 2D cluster ET vs max truth jet pT ─────────────────────────────────
    TH2D *h2 = dynamic_cast<TH2D *>(fin->Get("h2_clusterET_maxjet"));
    if (h2)
    {
        TCanvas *c2d = MakeCanvas("c_2d");
        c2d->SetRightMargin(0.13);
        gStyle->SetPalette(kBird);
        h2->GetXaxis()->SetTitleSize(tsize);
        h2->GetYaxis()->SetTitleSize(tsize);
        h2->GetXaxis()->SetLabelSize(tsize);
        h2->GetYaxis()->SetLabelSize(tsize);
        h2->GetXaxis()->SetTitleOffset(1.1);
        h2->GetYaxis()->SetTitleOffset(1.6);
        h2->GetXaxis()->CenterTitle();
        h2->GetYaxis()->CenterTitle();
        h2->Draw("colz");

        myText(0.18, 0.93, 1, strleg1.c_str(), tsize, false);
        myText(0.18, 0.88, 1, strleg2.c_str(), tsize, false);
        myText(0.18, 0.83, 1, strleg3.c_str(), tsize, false);

        c2d->SaveAs((outdir + "/truthjet_2d_clET_maxjet.pdf").c_str());
        delete c2d;
    }

    fin->Close();
    std::cout << "Plots saved to: " << outdir << "/" << std::endl;
}
