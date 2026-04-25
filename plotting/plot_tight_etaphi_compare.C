#include "plotcommon.h"

// Fill a tower-index (ieta, iphi) TH2 of "tight" clusters matching the nominal
// config_bdt_nom.yaml selection. Using tower-center integer indices removes the
// vertex-z dependence that smears the physical cluster_Eta/Phi values.
//
// Selections:
//   event:    |vertexz| < 10 cm
//             (if is_data) scaledtrigger[30] != 0  (Photon_4_GeV bit 30)
//   cluster:  cluster_Et > 10 GeV
//             common:   e11/e33 < 0.98, wr_cogx > 0, weta_cogx < 2.0, NPB > 0.5
//             tight:    weta_cogx < 1.0, wphi_cogx < 1.0, e32/e35 > 0.8,
//                       et1 > 0.5, BDT > (0.8666... - 0.006666...*ET)
//
// BDT branch: base_v1E for ET <= 35, base_E for ET > 35.
// Axes: X = cluster_ietacent (0..95), Y = cluster_iphicent (0..255).
static TH2F *fill_tight_etaphi(const char *tree_files, bool is_data,
                               const char *hname, long &n_events, long &n_tight)
{
    TChain *t = new TChain("slimtree");
    int nf = t->Add(tree_files);
    std::cout << "[" << hname << "] added " << nf << " files" << std::endl;

    const char *nd = "_CLUSTERINFO_CEMC";

    TTreeReader R(t);
    TTreeReaderValue<int>   ncl(R, Form("ncluster%s", nd));
    TTreeReaderArray<float> cEt (R, Form("cluster_Et%s", nd));
    TTreeReaderArray<float> cIeta(R, Form("cluster_ietacent%s", nd));
    TTreeReaderArray<float> cIphi(R, Form("cluster_iphicent%s", nd));
    TTreeReaderArray<float> cNPB(R, Form("cluster_npb_score%s", nd));
    TTreeReaderArray<float> cProb(R, Form("cluster_prob%s", nd));
    TTreeReaderArray<float> cE11(R, Form("cluster_e11%s", nd));
    TTreeReaderArray<float> cE33(R, Form("cluster_e33%s", nd));
    TTreeReaderArray<float> cE32(R, Form("cluster_e32%s", nd));
    TTreeReaderArray<float> cE35(R, Form("cluster_e35%s", nd));
    TTreeReaderArray<float> cEt1(R, Form("cluster_et1%s", nd));
    TTreeReaderArray<float> cWetaCogx(R, Form("cluster_weta_cogx%s", nd));
    TTreeReaderArray<float> cWphiCogx(R, Form("cluster_wphi_cogx%s", nd));
    TTreeReaderArray<float> cBDTv1E(R, Form("cluster_bdt%s_base_v1E", nd));
    TTreeReaderArray<float> cBDTE  (R, Form("cluster_bdt%s_base_E",   nd));
    TTreeReaderValue<float> vtxz(R, "vertexz");
    std::unique_ptr<TTreeReaderArray<bool>> trig;
    if (is_data) trig = std::make_unique<TTreeReaderArray<bool>>(R, "scaledtrigger");

    // 96 x 256 tower grid; bin i at i+0.5
    TH2F *h = new TH2F(hname, "", 96, 0, 96, 256, 0, 256);
    h->SetXTitle("cluster i#it{#eta} (tower)");
    h->SetYTitle("cluster i#it{#phi} (tower)");
    h->SetZTitle("tight clusters / tower");

    n_events = 0;
    n_tight  = 0;

    while (R.Next())
    {
        n_events++;
        if (std::fabs(*vtxz) > 10.0) continue;
        if (is_data) {
            if (trig->GetSize() <= 30 || !(*trig)[30]) continue;
        }
        for (int i = 0; i < *ncl; ++i)
        {
            float ET = cEt[i];
            if (ET < 10.0) continue;

            float e11e33 = (cE33[i] > 0) ? cE11[i] / cE33[i] : 0;
            float e32e35 = (cE35[i] > 0) ? cE32[i] / cE35[i] : 0;
            float wetax = cWetaCogx[i];
            float wphix = cWphiCogx[i];
            float wr    = (wetax > 0) ? wphix / wetax : 0;
            if (!(e11e33 > 0.0 && e11e33 < 0.98)) continue;
            if (!(wr > 0.0))                     continue;
            if (!(wetax < 2.0))                  continue;
            if (!(cNPB[i] > 0.5))                continue;
            if (!(wetax < 1.0))                  continue;
            if (!(wphix < 1.0))                  continue;
            if (!(e32e35 > 0.8))                 continue;
            if (!(cEt1[i] > 0.5))                continue;
            double bdt_min_et = 0.8666666666666668 - 0.006666666666666672 * ET;
            double bdt = (ET <= 35.0) ? cBDTv1E[i] : cBDTE[i];
            if (!(bdt > bdt_min_et && bdt <= 1.0)) continue;

            // ietacent / iphicent are stored as float in the slimtree; truncate
            // to integer tower index then add 0.5 for bin-centre.
            int ieta_tower = (int) cIeta[i];
            int iphi_tower = (int) cIphi[i];
            h->Fill(ieta_tower + 0.5, iphi_tower + 0.5);
            n_tight++;
        }
    }

    return h;
}

void plot_tight_etaphi_compare()
{
    init_plot();

    const char *mc_path   = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root";
    const char *data_glob = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/condorout/part_*_with_bdt_split.root";

    long n_ev_mc, n_tight_mc;
    long n_ev_da, n_tight_da;

    TH2F *h_mc = fill_tight_etaphi(mc_path,   false, "h_tight_mc",   n_ev_mc, n_tight_mc);
    TH2F *h_da = fill_tight_etaphi(data_glob, true,  "h_tight_data", n_ev_da, n_tight_da);

    std::cout << "-------------------------------------------\n"
              << "MC   photon10: events = " << n_ev_mc
              << "  tight clusters (ET>10) = " << n_tight_mc << "\n"
              << "Data          : events = " << n_ev_da
              << "  tight clusters (ET>10) = " << n_tight_da << "\n"
              << "-------------------------------------------" << std::endl;

    const char *outdir = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/plotting/figures";

    auto draw_single = [&](TH2F *h, const char *title, const char *subtitle,
                           const char *savename, bool logz,
                           double xmin, double xmax, double zmin, double zmax)
    {
        TCanvas *cc = new TCanvas(Form("c_%s", savename), "", 1000, 750);
        cc->SetRightMargin(0.17);
        cc->SetLeftMargin(0.12);
        cc->SetTopMargin(0.09);
        if (logz) cc->SetLogz();
        h->GetXaxis()->SetRangeUser(xmin, xmax);
        if (zmin != zmax) h->GetZaxis()->SetRangeUser(zmin, zmax);
        h->Draw("COLZ");
        TLatex lx; lx.SetNDC();
        lx.SetTextSize(0.04);
        lx.DrawLatex(0.13, 0.935, title);
        lx.SetTextSize(0.03);
        lx.DrawLatex(0.13, 0.895, subtitle);
        cc->SaveAs(Form("%s/tight_%s.pdf", outdir, savename));
        cc->SaveAs(Form("%s/tight_%s.png", outdir, savename));
    };

    // Full tower range: ieta in [0, 95]
    draw_single(h_mc, "#bf{#it{sPHENIX}} Internal -- MC photon10 (tight)",
                "|z_{vtx}|<10 cm, #it{E}_{T}^{cluster}>10 GeV, common+tight+BDT",
                "mc_full", true, 0, 96, 0, 0);
    draw_single(h_da, "#bf{#it{sPHENIX}} Internal -- Data (Photon 4 GeV trig, tight)",
                "|z_{vtx}|<10 cm, #it{E}_{T}^{cluster}>10 GeV, common+tight+BDT",
                "data_full", true, 0, 96, 0, 0);

    // -----------------------------------------------------------
    // Ratio and z-score (same machinery as before), fiducial |eta|<0.7
    // -----------------------------------------------------------
    double N_mc = 0, N_da = 0;
    const int nx = h_mc->GetNbinsX();
    const int ny = h_mc->GetNbinsY();
    // Full range: all 96 x 256 = 24576 bins.
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            N_mc += h_mc->GetBinContent(ix, iy);
            N_da += h_da->GetBinContent(ix, iy);
        }
    }
    std::cout << "full-range totals: N_mc = " << N_mc << "  N_da = " << N_da << std::endl;

    TH2F *h_ratio = (TH2F *) h_mc->Clone("h_ratio_tight");
    TH2F *h_z     = (TH2F *) h_mc->Clone("h_z_tight");
    h_ratio->Reset(); h_z->Reset();
    h_ratio->SetZTitle("p_{data}/p_{MC}");
    h_z->SetZTitle("z-score");

    TH1F *h_zdist = new TH1F("h_zdist_tight", ";z-score;bins", 120, -12, 12);

    int n_fid = 0, n_zlo = 0, n_zhi = 0, n_zvlo = 0, n_zvhi = 0;
    int n_empty_mc = 0, n_empty_da = 0;
    double def_zlt_m2 = 0, def_zlt_m5 = 0;

    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            double nm = h_mc->GetBinContent(ix, iy);
            double nd = h_da->GetBinContent(ix, iy);
            n_fid++;
            if (nm <= 0) { n_empty_mc++; continue; }
            if (nd <= 0)   n_empty_da++;
            double p_m = nm / N_mc;
            double p_d = nd / N_da;
            double R   = p_d / p_m;
            double relerr = std::sqrt(1.0 / std::max(nd, 1.0) + 1.0 / nm);
            double z = (R - 1.0) / (R * relerr);
            h_ratio->SetBinContent(ix, iy, R);
            h_z->SetBinContent(ix, iy, z);
            h_zdist->Fill(z);
            if (z < -2) {
                n_zlo++;
                double mu = nm * (N_da / N_mc);
                def_zlt_m2 += std::max(0.0, mu - nd);
            }
            if (z >  2) n_zhi++;
            if (z < -5) {
                n_zvlo++;
                double mu = nm * (N_da / N_mc);
                def_zlt_m5 += std::max(0.0, mu - nd);
            }
            if (z >  5) n_zvhi++;
        }
    }

    double f_anom = (double)(n_zlo + n_zhi) / n_fid;
    std::cout << "-------------------------------------------\n"
              << "TIGHT-cluster full-range summary (|eta|<1.1, ieta 0..95)\n"
              << "  total bins          : " << n_fid << "\n"
              << "  empty MC bins       : " << n_empty_mc << " (skipped)\n"
              << "  empty data bins     : " << n_empty_da << "\n"
              << "  |z|>2 total         : " << (n_zlo + n_zhi)
              << " (" << 100.0 * f_anom << " %; Gauss 4.55%)\n"
              << "     z<-2 (data<MC)   : " << n_zlo
              << " (" << 100.0 * n_zlo / (double)n_fid << " %)\n"
              << "     z>+2 (MC<data)   : " << n_zhi
              << " (" << 100.0 * n_zhi / (double)n_fid << " %)\n"
              << "  |z|>5 total         : " << (n_zvlo + n_zvhi)
              << " (" << 100.0 * (n_zvlo + n_zvhi) / (double)n_fid << " %)\n"
              << "     z<-5             : " << n_zvlo << "\n"
              << "     z>+5             : " << n_zvhi << "\n"
              << "Data deficit z<-2 : " << def_zlt_m2
              << " clusters (" << 100.0 * def_zlt_m2 / N_da << " % of full-range data)\n"
              << "Data deficit z<-5 : " << def_zlt_m5
              << " clusters (" << 100.0 * def_zlt_m5 / N_da << " % of full-range data)\n"
              << "-------------------------------------------" << std::endl;

    draw_single(h_ratio, "Tight clusters: normalized ratio  (data / MC)",
                Form("z<-2: %.2f%%  z>+2: %.2f%%  def: %.1f%% of data",
                     100.0 * n_zlo / (double)n_fid,
                     100.0 * n_zhi / (double)n_fid,
                     100.0 * def_zlt_m2 / N_da),
                "ratio", false, 0, 96, 0.0, 2.0);
    draw_single(h_z, "Tight clusters: z-score (data - MC, Poisson)",
                Form("def(z<-2)= %.1f%%  def(z<-5)= %.1f%% of data",
                     100.0 * def_zlt_m2 / N_da,
                     100.0 * def_zlt_m5 / N_da),
                "zscore", false, 0, 96, -6.0, 6.0);

    // z-score 1D
    TCanvas *cz = new TCanvas("cz", "", 900, 700);
    cz->SetLogy();
    h_zdist->SetXTitle("z-score (data vs MC tight, Poisson)");
    h_zdist->SetYTitle("bins / 0.2");
    h_zdist->SetLineColor(kBlack);
    h_zdist->SetLineWidth(2);
    h_zdist->Draw("HIST");
    double norm = h_zdist->Integral() * 0.2;
    TF1 *gaus = new TF1("gaus_ref_t", "[0]*TMath::Gaus(x,0,1,true)", -12, 12);
    gaus->SetParameter(0, norm);
    gaus->SetLineColor(kRed);
    gaus->SetLineWidth(2);
    gaus->SetLineStyle(7);
    gaus->Draw("same");
    TLatex lx; lx.SetNDC(); lx.SetTextSize(0.035);
    lx.DrawLatex(0.15, 0.93, strleg1.c_str());
    lx.DrawLatex(0.15, 0.89, "Tight clusters, |#eta|<0.7, ET>10 GeV");
    cz->SaveAs(Form("%s/tight_zdist.pdf", outdir));
    cz->SaveAs(Form("%s/tight_zdist.png", outdir));

    TFile *fo = TFile::Open(Form("%s/tight_etaphi_compare.root", outdir), "RECREATE");
    h_mc  ->Write("h_tight_mc");
    h_da  ->Write("h_tight_data");
    h_ratio->Write("h_ratio_tight");
    h_z    ->Write("h_z_tight");
    h_zdist->Write("h_zdist_tight");
    fo->Close();
}
