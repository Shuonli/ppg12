#include "plotcommon.h"
#include <yaml-cpp/yaml.h>

void plot_R_decomposition(string configname = "config_bdt_nom.yaml")
{
    init_plot();

    // ---- load config ----
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    std::string configpath = "/sphenix/user/shuhangli/ppg12/efficiencytool/" + configname;
    YAML::Node cfg = YAML::LoadFile(configpath);
    std::string var_type = cfg["output"]["var_type"].as<std::string>();
    std::cout << "[plot_R_decomposition] var_type = " << var_type << std::endl;

    // ---- open files ----
    // Merged MC (signal+background): source of leakage fractions and notmatch histograms
    std::string mcEffPath = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_" + var_type + ".root";
    // Jet-only MC: used as "data" surrogate in MC closure test
    std::string jetEffPath = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_" + var_type + ".root";
    // CalculatePhotonYield output: the official R
    std::string finalMcPath = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_" + var_type + "_mc.root";

    TFile *fMCeff = TFile::Open(mcEffPath.c_str(), "READ");
    TFile *fJet   = TFile::Open(jetEffPath.c_str(), "READ");
    TFile *fFinal = TFile::Open(finalMcPath.c_str(), "READ");

    if (!fMCeff || fMCeff->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << mcEffPath << std::endl;
        return;
    }
    if (!fJet || fJet->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << jetEffPath << std::endl;
        return;
    }
    if (!fFinal || fFinal->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << finalMcPath << std::endl;
        return;
    }

    // ---- Merged MC: signal ABCD for leakage fractions ----
    TH1D *hAs_mc = (TH1D *)fMCeff->Get("h_tight_iso_cluster_signal_0");
    TH1D *hBs_mc = (TH1D *)fMCeff->Get("h_tight_noniso_cluster_signal_0");
    TH1D *hCs_mc = (TH1D *)fMCeff->Get("h_nontight_iso_cluster_signal_0");
    TH1D *hDs_mc = (TH1D *)fMCeff->Get("h_nontight_noniso_cluster_signal_0");

    // ---- Jet file: total ABCD and signal (fragmentation photons in jet samples) ----
    TH1D *hA_jet = (TH1D *)fJet->Get("h_tight_iso_cluster_0");
    TH1D *hB_jet = (TH1D *)fJet->Get("h_tight_noniso_cluster_0");
    TH1D *hC_jet = (TH1D *)fJet->Get("h_nontight_iso_cluster_0");
    TH1D *hD_jet = (TH1D *)fJet->Get("h_nontight_noniso_cluster_0");

    TH1D *hAs_jet = (TH1D *)fJet->Get("h_tight_iso_cluster_signal_0");
    TH1D *hBs_jet = (TH1D *)fJet->Get("h_tight_noniso_cluster_signal_0");
    TH1D *hCs_jet = (TH1D *)fJet->Get("h_nontight_iso_cluster_signal_0");
    TH1D *hDs_jet = (TH1D *)fJet->Get("h_nontight_noniso_cluster_signal_0");

    // ---- Official h_R from CalculatePhotonYield ----
    TH1D *h_R_official = (TH1D *)fFinal->Get("h_R");

    // Validate
    if (!hAs_mc || !hBs_mc || !hCs_mc || !hDs_mc)
    {
        std::cerr << "ERROR: missing signal histograms in " << mcEffPath << std::endl;
        return;
    }
    if (!hA_jet || !hB_jet || !hC_jet || !hD_jet ||
        !hAs_jet || !hBs_jet || !hCs_jet || !hDs_jet)
    {
        std::cerr << "ERROR: missing ABCD histograms in " << jetEffPath << std::endl;
        return;
    }
    if (!h_R_official)
    {
        std::cerr << "ERROR: h_R not found in " << finalMcPath << std::endl;
        return;
    }

    int nbins = hA_jet->GetNbinsX();

    // ---- build TGraphErrors for each R decomposition ----
    //
    // Following CalculatePhotonYield.C (isMC branch):
    //   datainput  = MC_efficiency_jet_{var_type}.root  (jet-only, used as "data")
    //   siminput   = MC_efficiency_{var_type}.root      (merged, for leakage fractions)
    //
    // h_R_official: background_data = total_jet - signal_jet,
    //               R = (A_bkg*D_bkg) / (B_bkg*C_bkg)
    //
    // R_raw:  naive R from jet-file total yields (signal-contaminated)
    //         = (A_jet * D_jet) / (B_jet * C_jet)
    //
    // R_bkg:  truth background R from jet file (tot - sig from jet file)
    //         = ((A-As_jet) * (D-Ds_jet)) / ((B-Bs_jet) * (C-Cs_jet))
    //         This IS what h_R_official computes
    //
    // R_corrected: leakage-subtracted using merged-MC fractions
    //         cB = Bs_mc/As_mc,  cC = Cs_mc/As_mc,  cD = Ds_mc/As_mc
    //         X_corr = X_jet - cX * As_jet  (approximate signal subtraction)
    //         This tests whether the pipeline's leakage fractions from
    //         signal MC correctly subtract signal from the jet "data"
    //

    TGraphErrors *g_R_bkg = new TGraphErrors();        // pure bkg (truth, jet file)
    TGraphErrors *g_R_raw = new TGraphErrors();         // total, no corrections
    TGraphErrors *g_R_corrected = new TGraphErrors();   // leakage-subtracted via MC fractions
    TGraphErrors *g_R_official = new TGraphErrors();    // from CalculatePhotonYield
    TGraphErrors *g_ratio = new TGraphErrors();         // R_corrected / R_bkg

    int npt = 0;
    std::cout << Form("%8s %10s %10s %10s %10s %10s",
                      "pT", "R_bkg", "R_raw", "R_corr", "R_offic", "corr/bkg")
              << std::endl;

    for (int i = 1; i <= nbins; i++)
    {
        double ptCenter = hA_jet->GetXaxis()->GetBinCenter(i);

        // Jet file totals
        double A = hA_jet->GetBinContent(i);
        double B = hB_jet->GetBinContent(i);
        double C = hC_jet->GetBinContent(i);
        double D = hD_jet->GetBinContent(i);

        // Jet file signal (fragmentation photons in jet samples)
        double As_j = hAs_jet->GetBinContent(i);
        double Bs_j = hBs_jet->GetBinContent(i);
        double Cs_j = hCs_jet->GetBinContent(i);
        double Ds_j = hDs_jet->GetBinContent(i);

        // Merged MC signal (for leakage fractions)
        double As_mc = hAs_mc->GetBinContent(i);
        double Bs_mc = hBs_mc->GetBinContent(i);
        double Cs_mc = hCs_mc->GetBinContent(i);
        double Ds_mc = hDs_mc->GetBinContent(i);

        // Background = jet total - jet signal
        double Ab = A - As_j;
        double Bb = B - Bs_j;
        double Cb = C - Cs_j;
        double Db = D - Ds_j;

        // Guard against zero denominators
        if (B * C < 1e-6 || Bb < 1e-6 || Cb < 1e-6)
            continue;

        // R_raw = (A * D) / (B * C) -- signal-contaminated
        double R_raw = (A * D) / (B * C);

        // R_bkg = (Ab * Db) / (Bb * Cb) -- truth background
        double R_bkg = (Ab * Db) / (Bb * Cb);

        // R_corrected using pipeline leakage fractions
        // cX = X_sig_mc / A_sig_mc
        double cB = (As_mc > 0) ? Bs_mc / As_mc : 0;
        double cC = (As_mc > 0) ? Cs_mc / As_mc : 0;
        double cD = (As_mc > 0) ? Ds_mc / As_mc : 0;

        // Estimate signal in region A of jet file using pipeline logic:
        // For the R decomposition, we use As_j (the truth signal in jet A) as
        // the signal yield to subtract. The pipeline's ABCD solver finds NsigA
        // self-consistently, but for this plot we use the truth value.
        double NsigA = As_j;

        double A_corr = A - NsigA;
        double B_corr = B - cB * NsigA;
        double C_corr = C - cC * NsigA;
        double D_corr = D - cD * NsigA;

        double R_corr = -1;
        if (B_corr > 1e-6 && C_corr > 1e-6)
            R_corr = (A_corr * D_corr) / (B_corr * C_corr);

        double R_off = h_R_official->GetBinContent(i);

        double ratio = (R_bkg > 1e-6) ? R_corr / R_bkg : 0;

        std::cout << Form("[%2.0f-%2.0f] %10.4f %10.4f %10.4f %10.4f %10.4f",
                          hA_jet->GetXaxis()->GetBinLowEdge(i), hA_jet->GetXaxis()->GetBinUpEdge(i),
                          R_bkg, R_raw, R_corr, R_off, ratio)
                  << std::endl;

        g_R_bkg->SetPoint(npt, ptCenter, R_bkg);
        g_R_bkg->SetPointError(npt, 0, 0);

        g_R_raw->SetPoint(npt, ptCenter, R_raw);
        g_R_raw->SetPointError(npt, 0, 0);

        g_R_corrected->SetPoint(npt, ptCenter, R_corr);
        g_R_corrected->SetPointError(npt, 0, 0);

        g_R_official->SetPoint(npt, ptCenter, R_off);
        g_R_official->SetPointError(npt, 0, h_R_official->GetBinError(i));

        g_ratio->SetPoint(npt, ptCenter, ratio);
        g_ratio->SetPointError(npt, 0, 0);

        npt++;
    }

    // ---- styling ----
    // R_bkg: pure background, truth -- filled circles, blue
    g_R_bkg->SetMarkerStyle(20);
    g_R_bkg->SetMarkerSize(1.4);
    g_R_bkg->SetMarkerColor(kAzure + 2);
    g_R_bkg->SetLineColor(kAzure + 2);
    g_R_bkg->SetLineWidth(2);

    // R_raw: total, no corrections -- open squares, gray
    g_R_raw->SetMarkerStyle(25);
    g_R_raw->SetMarkerSize(1.3);
    g_R_raw->SetMarkerColor(kGray + 2);
    g_R_raw->SetLineColor(kGray + 2);
    g_R_raw->SetLineWidth(2);

    // R_corrected: leakage-subtracted -- filled triangles-up, red
    g_R_corrected->SetMarkerStyle(22);
    g_R_corrected->SetMarkerSize(1.5);
    g_R_corrected->SetMarkerColor(kRed + 1);
    g_R_corrected->SetLineColor(kRed + 1);
    g_R_corrected->SetLineWidth(2);

    // R_official: from CalculatePhotonYield -- open circles, black
    g_R_official->SetMarkerStyle(24);
    g_R_official->SetMarkerSize(1.4);
    g_R_official->SetMarkerColor(kBlack);
    g_R_official->SetLineColor(kBlack);
    g_R_official->SetLineWidth(2);

    // Ratio: R_corrected / R_bkg
    g_ratio->SetMarkerStyle(20);
    g_ratio->SetMarkerSize(1.3);
    g_ratio->SetMarkerColor(kViolet + 1);
    g_ratio->SetLineColor(kViolet + 1);
    g_ratio->SetLineWidth(2);

    // ---- canvas with upper + lower pads ----
    TCanvas *c = new TCanvas("c_R_decomp", "", 650, 750);

    // Upper pad (70%)
    TPad *pUp = new TPad("pUp", "", 0, 0.30, 1, 1.0);
    pUp->SetBottomMargin(0.02);
    pUp->SetTopMargin(0.06);
    pUp->SetLeftMargin(0.15);
    pUp->SetRightMargin(0.04);
    pUp->Draw();

    // Lower pad (30%)
    TPad *pLow = new TPad("pLow", "", 0, 0.0, 1, 0.30);
    pLow->SetTopMargin(0.02);
    pLow->SetBottomMargin(0.32);
    pLow->SetLeftMargin(0.15);
    pLow->SetRightMargin(0.04);
    pLow->Draw();

    // ---- upper pad: R vs pT ----
    pUp->cd();

    // Determine y-axis range from data
    double ymax_up = 2.2;
    for (int ip = 0; ip < g_R_raw->GetN(); ip++)
    {
        double x, y;
        g_R_raw->GetPoint(ip, x, y);
        if (y * 1.2 > ymax_up)
            ymax_up = y * 1.2;
    }

    TH1F *frameUp = new TH1F("frameUp", "", 100, 8, 36);
    frameUp->SetMinimum(0.0);
    frameUp->SetMaximum(ymax_up);
    frameUp->GetYaxis()->SetTitle("#it{R} = #frac{#it{A} #times #it{D}}{#it{B} #times #it{C}}");
    frameUp->GetYaxis()->SetTitleSize(0.055);
    frameUp->GetYaxis()->SetLabelSize(0.048);
    frameUp->GetYaxis()->SetTitleOffset(1.15);
    frameUp->GetXaxis()->SetLabelSize(0);
    frameUp->GetXaxis()->SetTickLength(0.03);
    frameUp->Draw("AXIS");

    // Dashed line at R = 1
    TLine *line1 = new TLine(8, 1, 36, 1);
    line1->SetLineColor(kGray + 2);
    line1->SetLineStyle(2);
    line1->SetLineWidth(2);
    line1->Draw("same");

    g_R_raw->Draw("P same");
    g_R_bkg->Draw("P same");
    g_R_corrected->Draw("P same");
    g_R_official->Draw("P same");

    // Labels
    float xpos = 0.20, ypos = 0.92, dy = 0.060, fontsize = 0.048;
    myText(xpos, ypos, 1, strleg1.c_str(), fontsize, 0);
    myText(xpos, ypos - dy, 1, strMC.c_str(), fontsize, 0);
    myText(0.93, ypos, 1, strleg3.c_str(), fontsize, 1);

    // Legend
    TLegend *leg = new TLegend(0.42, 0.55, 0.93, 0.88);
    legStyle(leg, 0.14, 0.044);
    leg->AddEntry(g_R_raw, "#it{R}_{raw} (total, no correction)", "pl");
    leg->AddEntry(g_R_bkg, "#it{R}_{bkg} (pure bkg, truth)", "pl");
    leg->AddEntry(g_R_corrected, "#it{R}_{corr} (leakage subtracted)", "pl");
    leg->AddEntry(g_R_official, "#it{R} (CalculatePhotonYield)", "pl");
    leg->Draw("same");

    // ---- lower pad: R_corrected / R_bkg ratio ----
    pLow->cd();

    TH1F *frameLow = new TH1F("frameLow", "", 100, 8, 36);
    frameLow->SetMinimum(0.0);
    frameLow->SetMaximum(2.0);
    frameLow->GetYaxis()->SetTitle("#it{R}_{corr} / #it{R}_{bkg}");
    frameLow->GetYaxis()->SetTitleSize(0.11);
    frameLow->GetYaxis()->SetLabelSize(0.10);
    frameLow->GetYaxis()->SetTitleOffset(0.50);
    frameLow->GetYaxis()->SetNdivisions(505);
    frameLow->GetXaxis()->SetTitle("#it{E}_{T}^{#gamma} [GeV]");
    frameLow->GetXaxis()->SetTitleSize(0.12);
    frameLow->GetXaxis()->SetLabelSize(0.10);
    frameLow->GetXaxis()->SetTitleOffset(0.95);
    frameLow->GetXaxis()->SetTickLength(0.08);
    frameLow->Draw("AXIS");

    TLine *line1_low = new TLine(8, 1, 36, 1);
    line1_low->SetLineColor(kGray + 2);
    line1_low->SetLineStyle(2);
    line1_low->SetLineWidth(2);
    line1_low->Draw("same");

    g_ratio->Draw("P same");

    // ---- save ----
    c->SaveAs(Form("figures/R_decomposition_%s.pdf", var_type.c_str()));

    std::cout << "[plot_R_decomposition] Saved figures/R_decomposition_" << var_type << ".pdf" << std::endl;

    fMCeff->Close();
    fJet->Close();
    fFinal->Close();
}
