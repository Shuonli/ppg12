#include "plotcommon.h"

// plot_saturation.C
//
// Produces two sets of plots from SaturationStudy.C output:
//
//   1. Saturated cluster E_T distributions: data vs inclusive jet MC (shape comparison + ratio)
//      - All clusters
//      - Saturated clusters (n_sat > 0)
//      - Fraction of saturated clusters vs E_T
//
//   2. Photon energy resolution from signal MC (truth-matched photon clusters):
//      - 2D colz: truth E vs reco E  (all / saturated / unsaturated)
//      - Profile overlay: mean E_reco / E_truth vs E_truth, comparing sat vs nosat
//
// Usage:
//   root -l -q 'plot_saturation.C("showershape")'

void plot_saturation(const std::string suffix = "showershape")
{
    init_plot();

    const std::string resdir  = "/sphenix/user/shuhangli/ppg12/efficiencytool/results";
    const std::string figdir  = "figures";

    const int col_data   = kBlack;
    const int col_jetmc  = kAzure + 2;
    const int col_sigmc  = kPink + 5;
    const int col_sat    = kRed - 4;
    const int col_nosat  = kGreen - 2;

    // -----------------------------------------------------------------------
    // Open files
    // -----------------------------------------------------------------------
    TFile *fsig  = TFile::Open(Form("%s/MC_efficiency_saturation_photon_%s.root", resdir.c_str(), suffix.c_str()));
    TFile *fjet  = TFile::Open(Form("%s/MC_efficiency_saturation_jet_%s.root",    resdir.c_str(), suffix.c_str()));
    TFile *fdata = TFile::Open(Form("%s/data_histo_saturation_%s.root",           resdir.c_str(), suffix.c_str()));

    if (!fsig  || fsig->IsZombie())  { std::cerr << "Cannot open signal MC file"  << std::endl; return; }
    if (!fjet  || fjet->IsZombie())  { std::cerr << "Cannot open jet MC file"     << std::endl; return; }
    if (!fdata || fdata->IsZombie()) { std::cerr << "Cannot open data file"       << std::endl; return; }

    // -----------------------------------------------------------------------
    // Load ET distribution histograms
    // -----------------------------------------------------------------------
    auto load1D = [](TFile *f, const char *name) -> TH1D* {
        TH1D *h = dynamic_cast<TH1D*>(f->Get(name));
        if (!h) std::cerr << "Cannot find " << name << " in " << f->GetName() << std::endl;
        return h;
    };

    TH1D *h_all_data    = load1D(fdata, "h_ET_all");
    TH1D *h_sat_data    = load1D(fdata, "h_ET_sat");
    TH1D *h_nosat_data  = load1D(fdata, "h_ET_nosat");
    TH1D *h_fnum_data   = load1D(fdata, "h_sat_fraction_num");
    TH1D *h_fden_data   = load1D(fdata, "h_sat_fraction_denom");

    TH1D *h_all_jet     = load1D(fjet,  "h_ET_all");
    TH1D *h_sat_jet     = load1D(fjet,  "h_ET_sat");
    TH1D *h_nosat_jet   = load1D(fjet,  "h_ET_nosat");
    TH1D *h_fnum_jet    = load1D(fjet,  "h_sat_fraction_num");
    TH1D *h_fden_jet    = load1D(fjet,  "h_sat_fraction_denom");

    // Load 2D energy resolution histograms from signal MC
    TH2D *h2_all  = dynamic_cast<TH2D*>(fsig->Get("h2_Etruth_Ereco_all"));
    TH2D *h2_sat  = dynamic_cast<TH2D*>(fsig->Get("h2_Etruth_Ereco_sat"));
    TH2D *h2_nosat = dynamic_cast<TH2D*>(fsig->Get("h2_Etruth_Ereco_nosat"));

    if (!h_all_data || !h_sat_data || !h_all_jet || !h_sat_jet) {
        std::cerr << "Missing ET distribution histograms." << std::endl;
        return;
    }
    if (!h2_all || !h2_sat || !h2_nosat) {
        std::cerr << "Missing 2D energy resolution histograms from signal MC." << std::endl;
        return;
    }

    // -----------------------------------------------------------------------
    // Common annotation helpers
    // -----------------------------------------------------------------------
    const float xL = 0.18, yT = 0.88, dy = 0.055, fs = 0.042;

    auto drawLabels = [&](const char *sample_str) {
        myText(xL, yT,        1, strleg1.c_str(),    fs);
        myText(xL, yT - dy,   1, strleg2.c_str(),    fs);
        myText(xL, yT - 2*dy, 1, sample_str,         fs);
    };

    // -----------------------------------------------------------------------
    // Helper: normalise a histogram clone to unit area
    // -----------------------------------------------------------------------
    auto norm = [](TH1D *h, const char *newname) -> TH1D* {
        TH1D *hc = dynamic_cast<TH1D*>(h->Clone(newname));
        if (hc->Integral() > 0) hc->Scale(1.0 / hc->Integral());
        return hc;
    };

    // -----------------------------------------------------------------------
    // Helper: build a two-pad canvas (top 70%, bottom 30%)
    // -----------------------------------------------------------------------
    auto makeSplitCanvas = [](const char *name) -> TCanvas* {
        TCanvas *c = new TCanvas(name, name, 600, 700);
        c->cd();
        TPad *p1 = new TPad("p1", "top",    0, 0.30, 1, 1.00); p1->SetBottomMargin(0.02); p1->Draw();
        TPad *p2 = new TPad("p2", "bottom", 0, 0.00, 1, 0.30); p2->SetTopMargin(0.04); p2->SetBottomMargin(0.32); p2->Draw();
        return c;
    };

    // -----------------------------------------------------------------------
    // Style: set histogram appearance
    // -----------------------------------------------------------------------
    auto styleData = [&](TH1D *h) {
        h->SetMarkerStyle(20); h->SetMarkerSize(0.9);
        h->SetMarkerColor(col_data); h->SetLineColor(col_data);
    };
    auto styleJet = [&](TH1D *h) {
        h->SetLineColor(col_jetmc); h->SetLineWidth(2);
        h->SetMarkerSize(0);
    };

    // =====================================================================
    // SECTION 1: Cluster ET distributions – data vs inclusive jet MC
    // =====================================================================

    // --- 1a. All clusters (normalised shape comparison + ratio) -----------
    {
        TH1D *hd = norm(h_all_data, "h_all_data_norm");
        TH1D *hj = norm(h_all_jet,  "h_all_jet_norm");
        styleData(hd);
        styleJet(hj);

        TH1D *hratio = dynamic_cast<TH1D*>(hd->Clone("h_all_ratio"));
        hratio->Divide(hj);
        hratio->SetYTitle("Data / MC");

        TCanvas *c = makeSplitCanvas("c_ET_all");

        // top pad
        TPad *p1 = dynamic_cast<TPad*>(c->FindObject("p1"));
        p1->cd(); p1->SetLogy();
        TH1F *fr = dynamic_cast<TH1F*>(frame_et_rec->Clone("fr_all"));
        fr->GetXaxis()->SetRangeUser(8, 40);
        fr->GetYaxis()->SetRangeUser(2e-4, 1.0);
        fr->SetYTitle("Normalised entries / GeV");
        fr->GetXaxis()->SetLabelSize(0);
        fr->Draw("axis");
        hj->Draw("same hist");
        hd->Draw("same ex0");

        TLegend *leg = new TLegend(0.55, 0.55, 0.90, 0.75);
        legStyle(leg, 0.17, 0.042);
        leg->AddEntry(hd, "Data",             "pl");
        leg->AddEntry(hj, strIncMC.c_str(),   "l");
        leg->Draw();
        drawLabels("All clusters");

        // bottom pad: ratio
        p1->cd(); // back to top for label placement
        c->cd();
        TPad *p2 = dynamic_cast<TPad*>(c->FindObject("p2"));
        p2->cd();
        TH1F *fr2 = dynamic_cast<TH1F*>(frame_et_rec->Clone("fr2_all"));
        fr2->GetXaxis()->SetRangeUser(8, 40);
        fr2->GetYaxis()->SetRangeUser(0.0, 2.5);
        fr2->SetYTitle("Data / MC");
        fr2->GetYaxis()->SetTitleSize(0.12);  fr2->GetYaxis()->SetTitleOffset(0.45);
        fr2->GetYaxis()->SetLabelSize(0.10);
        fr2->GetXaxis()->SetTitleSize(0.13);  fr2->GetXaxis()->SetLabelSize(0.12);
        fr2->Draw("axis");
        hratio->SetMarkerStyle(20); hratio->SetMarkerSize(0.7);
        hratio->SetMarkerColor(col_data); hratio->SetLineColor(col_data);
        hratio->Draw("same ex0");
        lineone->Draw("same");

        c->SaveAs(Form("%s/sat_ET_all.pdf", figdir.c_str()));
    }

    // --- 1b. Saturated clusters (normalised shape comparison + ratio) -----
    {
        TH1D *hd = norm(h_sat_data, "h_sat_data_norm");
        TH1D *hj = norm(h_sat_jet,  "h_sat_jet_norm");
        styleData(hd);
        styleJet(hj);

        TH1D *hratio = dynamic_cast<TH1D*>(hd->Clone("h_sat_ratio"));
        hratio->Divide(hj);

        TCanvas *c = makeSplitCanvas("c_ET_sat");

        TPad *p1 = dynamic_cast<TPad*>(c->FindObject("p1"));
        p1->cd(); p1->SetLogy();
        TH1F *fr = dynamic_cast<TH1F*>(frame_et_rec->Clone("fr_sat"));
        fr->GetXaxis()->SetRangeUser(8, 40);
        fr->GetYaxis()->SetRangeUser(2e-4, 1.0);
        fr->SetYTitle("Normalised entries / GeV");
        fr->GetXaxis()->SetLabelSize(0);
        fr->Draw("axis");
        hj->Draw("same hist");
        hd->Draw("same ex0");

        TLegend *leg = new TLegend(0.55, 0.55, 0.90, 0.75);
        legStyle(leg, 0.17, 0.042);
        leg->AddEntry(hd, "Data",           "pl");
        leg->AddEntry(hj, strIncMC.c_str(), "l");
        leg->Draw();
        drawLabels("Saturated clusters (n_{sat} > 0)");

        c->cd();
        TPad *p2 = dynamic_cast<TPad*>(c->FindObject("p2"));
        p2->cd();
        TH1F *fr2 = dynamic_cast<TH1F*>(frame_et_rec->Clone("fr2_sat"));
        fr2->GetXaxis()->SetRangeUser(8, 40);
        fr2->GetYaxis()->SetRangeUser(0.0, 2.5);
        fr2->SetYTitle("Data / MC");
        fr2->GetYaxis()->SetTitleSize(0.12);  fr2->GetYaxis()->SetTitleOffset(0.45);
        fr2->GetYaxis()->SetLabelSize(0.10);
        fr2->GetXaxis()->SetTitleSize(0.13);  fr2->GetXaxis()->SetLabelSize(0.12);
        fr2->Draw("axis");
        hratio->SetMarkerStyle(20); hratio->SetMarkerSize(0.7);
        hratio->SetMarkerColor(col_data); hratio->SetLineColor(col_data);
        hratio->Draw("same ex0");
        lineone->Draw("same");

        c->SaveAs(Form("%s/sat_ET_sat.pdf", figdir.c_str()));
    }

    // --- 1c. Saturation fraction vs ET: data vs jet MC -------------------
    {
        TH1D *hf_data = dynamic_cast<TH1D*>(h_fnum_data->Clone("hf_data"));
        hf_data->Divide(h_fden_data);
        TH1D *hf_jet  = dynamic_cast<TH1D*>(h_fnum_jet->Clone("hf_jet"));
        hf_jet->Divide(h_fden_jet);

        styleData(hf_data);
        styleJet(hf_jet);

        TCanvas *c = new TCanvas("c_sat_frac", "c_sat_frac", 600, 600);
        TH1F *fr = dynamic_cast<TH1F*>(frame_et_rec->Clone("fr_frac"));
        fr->GetXaxis()->SetRangeUser(8, 40);
        fr->GetYaxis()->SetRangeUser(0.0, 0.6);
        fr->SetYTitle("Fraction of clusters with n_{sat} > 0");
        fr->Draw("axis");
        hf_jet->Draw("same hist");
        hf_data->Draw("same ex0");

        TLegend *leg = new TLegend(0.55, 0.62, 0.90, 0.78);
        legStyle(leg, 0.17, 0.042);
        leg->AddEntry(hf_data, "Data",           "pl");
        leg->AddEntry(hf_jet,  strIncMC.c_str(), "l");
        leg->Draw();
        drawLabels(strleg3.c_str());

        c->SaveAs(Form("%s/sat_fraction_vs_ET.pdf", figdir.c_str()));
    }

    // --- 1d. All / saturated / unsaturated overlaid (data, log scale) -----
    {
        TH1D *ha = dynamic_cast<TH1D*>(h_all_data->Clone("ha_d"));
        TH1D *hs = dynamic_cast<TH1D*>(h_sat_data->Clone("hs_d"));
        TH1D *hn = dynamic_cast<TH1D*>(h_nosat_data->Clone("hn_d"));
        if (ha->Integral() > 0) { hs->Scale(1./ha->Integral()); hn->Scale(1./ha->Integral()); ha->Scale(1./ha->Integral()); }

        ha->SetLineColor(kBlack);    ha->SetLineWidth(2);  ha->SetMarkerSize(0);
        hs->SetLineColor(col_sat);   hs->SetLineWidth(2);  hs->SetMarkerSize(0);
        hn->SetLineColor(col_nosat); hn->SetLineWidth(2);  hn->SetMarkerSize(0);

        TCanvas *c = new TCanvas("c_sat_breakdown_data", "c_sat_breakdown_data", 600, 600);
        c->SetLogy();
        TH1F *fr = dynamic_cast<TH1F*>(frame_et_rec->Clone("fr_breakdown"));
        fr->GetXaxis()->SetRangeUser(8, 40);
        fr->GetYaxis()->SetRangeUser(1e-5, 1.0);
        fr->SetYTitle("Fraction of clusters / GeV");
        fr->Draw("axis");
        ha->Draw("same hist");
        hs->Draw("same hist");
        hn->Draw("same hist");

        TLegend *leg = new TLegend(0.55, 0.60, 0.90, 0.78);
        legStyle(leg, 0.17, 0.042);
        leg->AddEntry(ha, "All",               "l");
        leg->AddEntry(hn, "n_{sat} = 0",       "l");
        leg->AddEntry(hs, "n_{sat} > 0",       "l");
        leg->Draw();
        drawLabels("Data");

        c->SaveAs(Form("%s/sat_ET_breakdown_data.pdf", figdir.c_str()));
    }

    // =====================================================================
    // SECTION 2: Photon energy resolution from signal MC
    // =====================================================================

    gStyle->SetPalette(kBird);

    // --- 2a. 2D colz: truth E vs reco E, all truth-matched ---------------
    {
        TCanvas *c = new TCanvas("c_res2d_all", "c_res2d_all", 640, 600);
        c->SetRightMargin(0.15);
        c->SetLogz();
        TH2D *h = dynamic_cast<TH2D*>(h2_all->Clone("h2_all_plot"));
        h->GetXaxis()->SetRangeUser(0, 40);
        h->GetYaxis()->SetRangeUser(0, 40);
        h->SetMinimum(1);
        h->SetXTitle("Truth photon #it{E} [GeV]");
        h->SetYTitle("Reco cluster #it{E} [GeV]");
        h->Draw("colz");

        // diagonal line (perfect resolution)
        TLine *diag = new TLine(0, 0, 40, 40);
        diag->SetLineColor(kRed); diag->SetLineStyle(7); diag->SetLineWidth(2);
        diag->Draw("same");

        myText(xL, yT,        1, strleg1.c_str(),    fs);
        myText(xL, yT - dy,   1, strleg2.c_str(),    fs);
        myText(xL, yT - 2*dy, 1, strSigMC.c_str(),   fs);
        myText(xL, yT - 3*dy, 1, "All truth-matched #gamma", fs);

        c->SaveAs(Form("%s/photon_res2d_all.pdf", figdir.c_str()));
    }

    // --- 2b. 2D colz: saturated clusters ----------------------------------
    {
        TCanvas *c = new TCanvas("c_res2d_sat", "c_res2d_sat", 640, 600);
        c->SetRightMargin(0.15);
        c->SetLogz();
        TH2D *h = dynamic_cast<TH2D*>(h2_sat->Clone("h2_sat_plot"));
        h->GetXaxis()->SetRangeUser(0, 40);
        h->GetYaxis()->SetRangeUser(0, 40);
        h->SetMinimum(1);
        h->SetXTitle("Truth photon #it{E} [GeV]");
        h->SetYTitle("Reco cluster #it{E} [GeV]");
        h->Draw("colz");

        TLine *diag = new TLine(0, 0, 40, 40);
        diag->SetLineColor(kRed); diag->SetLineStyle(7); diag->SetLineWidth(2);
        diag->Draw("same");

        myText(xL, yT,        1, strleg1.c_str(),         fs);
        myText(xL, yT - dy,   1, strleg2.c_str(),         fs);
        myText(xL, yT - 2*dy, 1, strSigMC.c_str(),        fs);
        myText(xL, yT - 3*dy, 1, "Saturated (n_{sat} > 0)", fs);

        c->SaveAs(Form("%s/photon_res2d_sat.pdf", figdir.c_str()));
    }

    // --- 2c. 2D colz: unsaturated clusters --------------------------------
    {
        TCanvas *c = new TCanvas("c_res2d_nosat", "c_res2d_nosat", 640, 600);
        c->SetRightMargin(0.15);
        c->SetLogz();
        TH2D *h = dynamic_cast<TH2D*>(h2_nosat->Clone("h2_nosat_plot"));
        h->GetXaxis()->SetRangeUser(0, 40);
        h->GetYaxis()->SetRangeUser(0, 40);
        h->SetMinimum(1);
        h->SetXTitle("Truth photon #it{E} [GeV]");
        h->SetYTitle("Reco cluster #it{E} [GeV]");
        h->Draw("colz");

        TLine *diag = new TLine(0, 0, 40, 40);
        diag->SetLineColor(kRed); diag->SetLineStyle(7); diag->SetLineWidth(2);
        diag->Draw("same");

        myText(xL, yT,        1, strleg1.c_str(),      fs);
        myText(xL, yT - dy,   1, strleg2.c_str(),      fs);
        myText(xL, yT - 2*dy, 1, strSigMC.c_str(),     fs);
        myText(xL, yT - 3*dy, 1, "Unsaturated (n_{sat} = 0)", fs);

        c->SaveAs(Form("%s/photon_res2d_nosat.pdf", figdir.c_str()));
    }

    // Shared helper: for each truth-E column in h2, compute the weighted
    // mean and sigma of response r = E_reco / E_truth directly from bin
    // contents. Returns {mean_graph, sigma_graph}.
    auto makeResponseGraphs = [](TH2D *h2, const char *tag)
        -> std::pair<TGraphErrors*, TGraphErrors*>
    {
        TGraphErrors *gr_mean  = new TGraphErrors();
        TGraphErrors *gr_sigma = new TGraphErrors();
        int npt = 0;
        const int nxb = h2->GetNbinsX();
        const int nyb = h2->GetNbinsY();
        for (int ix = 1; ix <= nxb; ix++) {
            double Etruth = h2->GetXaxis()->GetBinCenter(ix);
            if (Etruth <= 0) continue;
            double sum_w = 0, sum_w2 = 0, sum_r = 0, sum_r2 = 0;
            for (int iy = 1; iy <= nyb; iy++) {
                double Ereco = h2->GetYaxis()->GetBinCenter(iy);
                double w     = h2->GetBinContent(ix, iy);
                double we    = h2->GetBinError(ix, iy);  // sqrt(sum_w_i^2), requires Sumw2()
                if (w <= 0) continue;
                double r  = Ereco / Etruth;
                sum_w  += w;
                sum_w2 += we * we;  // accumulate sum of w_i^2 across y-bins
                sum_r  += w * r;
                sum_r2 += w * r * r;
            }
            // N_eff = (sum_w)^2 / sum(w_i^2) — effective number of events
            // Use N_eff (not sum_w) for error estimates and minimum threshold.
            if (sum_w2 <= 0) continue;
            double Neff = sum_w * sum_w / sum_w2;
            if (Neff < 5) continue;
            double mean  = sum_r / sum_w;
            double var   = std::max(0.0, sum_r2 / sum_w - mean * mean);
            double sigma = std::sqrt(var);
            // statistical errors on mean and sigma using effective event count
            double err_mean  = sigma / std::sqrt(Neff);
            double err_sigma = sigma / std::sqrt(2.0 * Neff);
            gr_mean->SetPoint(npt,  Etruth, mean);
            gr_mean->SetPointError(npt, 0, err_mean);
            gr_sigma->SetPoint(npt, Etruth, sigma);
            gr_sigma->SetPointError(npt, 0, err_sigma);
            npt++;
        }
        return {gr_mean, gr_sigma};
    };

    auto [gr_mean_all,   gr_sigma_all]   = makeResponseGraphs(h2_all,   "all");
    auto [gr_mean_sat,   gr_sigma_sat]   = makeResponseGraphs(h2_sat,   "sat");
    auto [gr_mean_nosat, gr_sigma_nosat] = makeResponseGraphs(h2_nosat, "nosat");

    // common graph styling
    auto styleGraph = [](TGraphErrors *gr, int marker, int color) {
        gr->SetMarkerStyle(marker); gr->SetMarkerSize(0.9);
        gr->SetMarkerColor(color);  gr->SetLineColor(color);
    };
    styleGraph(gr_mean_all,   20, kBlack);
    styleGraph(gr_mean_sat,   21, col_sat);
    styleGraph(gr_mean_nosat, 24, col_nosat);
    styleGraph(gr_sigma_all,   20, kBlack);
    styleGraph(gr_sigma_sat,   21, col_sat);
    styleGraph(gr_sigma_nosat, 24, col_nosat);

    // --- 2d. Mean response <E_reco / E_truth> vs truth E_T ---------------
    {
        TCanvas *c = new TCanvas("c_res_mean", "c_res_mean", 600, 600);
        TH1F *fr = dynamic_cast<TH1F*>(frame_et_truth->Clone("fr_resmean"));
        fr->GetXaxis()->SetRangeUser(5, 40);
        fr->GetYaxis()->SetRangeUser(0.60, 1.35);
        fr->SetXTitle("Truth photon #it{E} [GeV]");
        fr->SetYTitle("#LT#it{E}_{reco} / #it{E}_{truth}#GT");
        fr->Draw("axis");

        gr_mean_all->Draw("same p");
        gr_mean_nosat->Draw("same p");
        gr_mean_sat->Draw("same p");
        lineone->Draw("same");

        TLegend *leg = new TLegend(0.55, 0.22, 0.90, 0.42);
        legStyle(leg, 0.17, 0.042);
        leg->AddEntry(gr_mean_all,   "All",         "pl");
        leg->AddEntry(gr_mean_nosat, "n_{sat} = 0", "pl");
        leg->AddEntry(gr_mean_sat,   "n_{sat} > 0", "pl");
        leg->Draw();

        myText(xL, yT,        1, strleg1.c_str(),  fs);
        myText(xL, yT - dy,   1, strleg2.c_str(),  fs);
        myText(xL, yT - 2*dy, 1, strSigMC.c_str(), fs);

        c->SaveAs(Form("%s/photon_res_mean.pdf", figdir.c_str()));
    }

    // --- 2e. Response width sigma(E_reco/E_truth) vs truth E_T -----------
    {
        TCanvas *c = new TCanvas("c_res_width", "c_res_width", 600, 600);
        TH1F *fr = dynamic_cast<TH1F*>(frame_et_truth->Clone("fr_reswidth"));
        fr->GetXaxis()->SetRangeUser(5, 40);
        fr->GetYaxis()->SetRangeUser(0.0, 0.55);
        fr->SetXTitle("Truth photon #it{E} [GeV]");
        fr->SetYTitle("#sigma(#it{E}_{reco} / #it{E}_{truth})");
        fr->Draw("axis");

        gr_sigma_all->Draw("same p");
        gr_sigma_nosat->Draw("same p");
        gr_sigma_sat->Draw("same p");

        TLegend *leg = new TLegend(0.55, 0.62, 0.90, 0.80);
        legStyle(leg, 0.17, 0.042);
        leg->AddEntry(gr_sigma_all,   "All",         "pl");
        leg->AddEntry(gr_sigma_nosat, "n_{sat} = 0", "pl");
        leg->AddEntry(gr_sigma_sat,   "n_{sat} > 0", "pl");
        leg->Draw();

        myText(xL, yT,        1, strleg1.c_str(),  fs);
        myText(xL, yT - dy,   1, strleg2.c_str(),  fs);
        myText(xL, yT - 2*dy, 1, strSigMC.c_str(), fs);

        c->SaveAs(Form("%s/photon_res_width.pdf", figdir.c_str()));
    }

    std::cout << "Plots saved to " << figdir << "/" << std::endl;
}
