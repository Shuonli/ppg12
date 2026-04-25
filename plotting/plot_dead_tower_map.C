#include "plotcommon.h"

// Categorical tower map per selection level. Z-score uses the empirical
// R = p_data / p_MC distribution (core-Gaussian fit to the bulk), not
// per-bin Poisson. Identical definition to make_tower_masks.C.
//
//   0 = normal (|z| <= 2 or MC empty, drawn blank/white)
//   1 = -5 < z < -2  (moderate data deficit)
//   2 =        z < -5 (severe data deficit; often empty under new method
//                      because R>=0 bounds the lower tail)
//   3 = hard dead    (n_data=0 AND n_MC>0; always categorised here)

static std::pair<double, double> fit_logR(TH2F *h_mc, TH2F *h_da, const char *lbl)
{
    // log(R) is roughly Gaussian (count-ratio statistics are log-normal).
    // Use 3-sigma iterative clip directly on log(R).
    const double Nm = h_mc->Integral(), Nd = h_da->Integral();
    std::vector<double> logR;
    for (int ix = 1; ix <= h_mc->GetNbinsX(); ++ix)
      for (int iy = 1; iy <= h_mc->GetNbinsY(); ++iy) {
        double nm = h_mc->GetBinContent(ix,iy), nd = h_da->GetBinContent(ix,iy);
        if (nm <= 0 || nd <= 0) continue;
        double r = (nd/Nd) / (nm/Nm);
        if (r <= 0) continue;
        logR.push_back(std::log(r));
      }
    if (logR.empty()) return {0.0, 1.0};
    double mean = 0, sigma = 0;
    for (double x : logR) mean += x; mean /= logR.size();
    for (double x : logR) sigma += (x - mean) * (x - mean);
    sigma = std::sqrt(sigma / logR.size());
    for (int it = 0; it < 3; ++it) {
        double m2 = 0, s2 = 0; int nk = 0;
        for (double x : logR) if (std::fabs(x - mean) < 3.0 * sigma) { m2 += x; nk++; }
        if (nk == 0) break;
        mean = m2 / nk;
        for (double x : logR) if (std::fabs(x - mean) < 3.0 * sigma) s2 += (x - mean) * (x - mean);
        sigma = std::sqrt(s2 / nk);
    }
    std::cout << "[" << lbl << "] log(R) fit: mean=" << mean << " sigma=" << sigma
              << " (<R>=exp(mean)=" << std::exp(mean) << ")" << std::endl;
    return {mean, sigma};
}

static TH2F *build_category(TH2F *h_mc, TH2F *h_da, const std::string &lvl)
{
    auto stats = fit_logR(h_mc, h_da, lvl.c_str());
    double mean_lR = stats.first, sigma_lR = stats.second;

    double N_mc = h_mc->Integral();
    double N_da = h_da->Integral();
    const int nx = h_mc->GetNbinsX();
    const int ny = h_mc->GetNbinsY();

    TH2F *h_cat = (TH2F *) h_mc->Clone(Form("h_dead_cat_%s", lvl.c_str()));
    h_cat->Reset();
    h_cat->SetTitle("");
    h_cat->SetXTitle("cluster i#it{#eta} (tower)");
    h_cat->SetYTitle("cluster i#it{#phi} (tower)");
    h_cat->SetZTitle("dead-tower category");
    h_cat->GetZaxis()->SetRangeUser(0, 3);

    int n_hard = 0, n_zlt5 = 0, n_zlt2 = 0;
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            double nm = h_mc->GetBinContent(ix, iy);
            double nd = h_da->GetBinContent(ix, iy);
            if (nm <= 0) continue;
            if (nd <= 0) {
                h_cat->SetBinContent(ix, iy, 3.0);
                n_hard++;
                continue;
            }
            double p_m = nm / N_mc;
            double p_d = nd / N_da;
            double r   = p_d / p_m;
            if (r <= 0) continue;
            double z   = (std::log(r) - mean_lR) / sigma_lR;
            if (z < -5) {
                h_cat->SetBinContent(ix, iy, 2.0);
                n_zlt5++;
            } else if (z < -2) {
                h_cat->SetBinContent(ix, iy, 1.0);
                n_zlt2++;
            }
        }
    }
    std::cout << "[" << lvl << "] hard-dead=" << n_hard
              << "  z<-5=" << n_zlt5
              << "  -5<z<-2=" << n_zlt2 << std::endl;
    return h_cat;
}

void plot_dead_tower_map()
{
    init_plot();
    const char *infile = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root";
    TFile *f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) { std::cerr << "cannot open " << infile << std::endl; return; }

    const std::vector<std::string> levels = {"preselect", "common", "tight", "tight_iso"};
    const char *outdir = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures";

    const int NCOL = 4;
    int palette[NCOL];
    palette[0] = TColor::GetColor("#F4F4F4");
    palette[1] = TColor::GetColor("#FFCC33");
    palette[2] = TColor::GetColor("#FF6600");
    palette[3] = TColor::GetColor("#B00000");

    TFile *fo = TFile::Open(Form("%s/dead_tower_map.root", outdir), "RECREATE");

    for (const auto &lvl : levels)
    {
        TH2F *h_mc = (TH2F *) f->Get(Form("h_etaphi_tower_%s_mc_inclusive", lvl.c_str()));
        TH2F *h_da = (TH2F *) f->Get(Form("h_etaphi_tower_%s_data", lvl.c_str()));
        if (!h_mc || !h_da) { std::cerr << "missing " << lvl << ", skip" << std::endl; continue; }
        TH2F *h_cat = build_category(h_mc, h_da, lvl);
        h_cat->Write();

        gStyle->SetPalette(NCOL, palette);
        TCanvas *cc = new TCanvas(Form("c_cat_%s", lvl.c_str()), "", 1100, 820);
        cc->SetRightMargin(0.22);
        cc->SetLeftMargin(0.12);
        cc->SetTopMargin(0.11);
        h_cat->GetXaxis()->SetRangeUser(0, 96);
        h_cat->GetYaxis()->SetRangeUser(0, 256);
        h_cat->GetZaxis()->SetRangeUser(-0.5, 3.5);
        h_cat->GetZaxis()->SetNdivisions(4);
        h_cat->Draw("COL");

        int n_cat[NCOL] = {0, 0, 0, 0};
        for (int ix = 1; ix <= h_cat->GetNbinsX(); ++ix)
            for (int iy = 1; iy <= h_cat->GetNbinsY(); ++iy)
                n_cat[(int)h_cat->GetBinContent(ix, iy)]++;

        TLatex lx; lx.SetNDC();
        lx.SetTextSize(0.034);
        lx.DrawLatex(0.13, 0.955,
            Form("#bf{#it{sPHENIX}} Internal -- tower category map (%s, R-stat z)", lvl.c_str()));
        lx.SetTextSize(0.025);
        lx.DrawLatex(0.13, 0.915,
            Form("hard dead = %d,  z<-5 = %d,  -5<z<-2 = %d", n_cat[3], n_cat[2], n_cat[1]));

        auto drawbox = [&](double x, double y, int col, const char *txt) {
            TPave *b = new TPave(x, y, x + 0.018, y + 0.025, 1, "NDC");
            b->SetFillColor(col); b->SetLineColor(kBlack); b->SetLineWidth(1); b->SetBorderSize(1);
            b->Draw();
            TLatex t; t.SetNDC(); t.SetTextSize(0.022);
            t.DrawLatex(x + 0.022, y + 0.005, txt);
        };
        drawbox(0.80, 0.84, palette[3], "hard dead (n_{data}=0)");
        drawbox(0.80, 0.80, palette[2], "z < -5 (R-stat)");
        drawbox(0.80, 0.76, palette[1], "-5 < z < -2");
        drawbox(0.80, 0.72, palette[0], "|z| #leq 2 or MC=0");

        TBox *fid = new TBox(17, 0, 79, 256);
        fid->SetFillStyle(0);
        fid->SetLineColor(kBlue + 1);
        fid->SetLineStyle(2);
        fid->SetLineWidth(2);
        fid->Draw();
        TLatex fl; fl.SetNDC(); fl.SetTextSize(0.02); fl.SetTextColor(kBlue + 1);
        fl.DrawLatex(0.80, 0.68, "blue dashed: |#it{#eta}|<0.7");

        cc->SaveAs(Form("%s/dead_tower_map_%s.pdf", outdir, lvl.c_str()));
        cc->SaveAs(Form("%s/dead_tower_map_%s.png", outdir, lvl.c_str()));
    }
    fo->Close();
    f->Close();
}
