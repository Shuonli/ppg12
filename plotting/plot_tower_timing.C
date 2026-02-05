#include "plotcommon.h"
#include <vector>
#include <algorithm>
#include <memory>
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TSystem.h"

const int col_samples[] = {kBlack, kRed, kBlue, kGreen + 2, kMagenta + 2, kOrange + 1};

namespace
{
    struct MeanRmsGraphs
    {
        TGraphErrors *g_mean = nullptr;
        TGraphErrors *g_rms = nullptr;
    };

    // Forward decl (used by helpers below)
    TH2D *project_xy_in_eta(const TH3D *h3,
                            const std::string &name,
                            double eta_lo,
                            double eta_hi);

    MeanRmsGraphs make_mean_rms_vs_x(const TH2D *h2,
                                    const std::string &name_prefix,
                                    int min_entries_per_xbin = 5)
    {
        MeanRmsGraphs out;
        if (!h2)
            return out;

        const int nbx = h2->GetNbinsX();
        std::vector<double> xs, exs, means, emeans, rmss, ermss;
        xs.reserve(nbx);
        exs.reserve(nbx);
        means.reserve(nbx);
        emeans.reserve(nbx);
        rmss.reserve(nbx);
        ermss.reserve(nbx);

        // For each X bin, project Y and compute mean/RMS of Y distribution
        for (int ibx = 1; ibx <= nbx; ++ibx)
        {
            std::unique_ptr<TH1D> proj((TH1D *)h2->ProjectionY(Form("%s_py_%d", name_prefix.c_str(), ibx), ibx, ibx));
            if (!proj)
                continue;

            const double n = proj->GetEntries();
            if (n < min_entries_per_xbin)
                continue;

            const double x = h2->GetXaxis()->GetBinCenter(ibx);
            const double ex = 0.5 * h2->GetXaxis()->GetBinWidth(ibx);
            const double mean = proj->GetMean();
            const double emean = proj->GetMeanError();
            const double rms = proj->GetRMS();
            const double erms = proj->GetRMSError();

            xs.push_back(x);
            exs.push_back(ex);
            means.push_back(mean);
            emeans.push_back(emean);
            rmss.push_back(rms);
            ermss.push_back(erms);
        }

        out.g_mean = new TGraphErrors(xs.size());
        out.g_rms = new TGraphErrors(xs.size());
        out.g_mean->SetName(Form("%s_g_mean", name_prefix.c_str()));
        out.g_rms->SetName(Form("%s_g_rms", name_prefix.c_str()));

        for (int i = 0; i < (int)xs.size(); ++i)
        {
            out.g_mean->SetPoint(i, xs[i], means[i]);
            out.g_mean->SetPointError(i, exs[i], emeans[i]);
            out.g_rms->SetPoint(i, xs[i], rmss[i]);
            out.g_rms->SetPointError(i, exs[i], ermss[i]);
        }
        return out;
    }

    struct EtaBin
    {
        double lo = 0.0;
        double hi = 0.0;
        std::string label;
    };

    std::vector<EtaBin> default_eta_bins()
    {
        return {
            {-0.7, -0.35, "-0.7 < #eta < -0.35"},
            {-0.1, 0.1, "-0.1 < #eta < 0.1"},
            {0.35, 0.7, "0.35 < #eta < 0.7"}};
    }

    // Ground-truth method: restrict eta, project TH3 -> TH2(X,Y), then use TH2::ProfileY to get <X>(Y).
    // This matches the original 2D logic used elsewhere in this macro.
    TProfile *make_mean_x_vs_y_in_eta_via_projection(const TH3D *h3,
                                                     const std::string &name,
                                                     double eta_lo,
                                                     double eta_hi)
    {
        TH2D *h2 = project_xy_in_eta(h3, Form("%s__h2tmp", name.c_str()), eta_lo, eta_hi);
        if (!h2)
            return nullptr;
        TProfile *p = h2->ProfileY(name.c_str(), 1, -1, "");
        if (p)
            p->SetDirectory(nullptr);
        delete h2;
        return p;
    }

    // Convenience: project TH3 (X=time, Y=value, Z=eta) to TH2 (X,Y) in an eta window.
    TH2D *project_xy_in_eta(const TH3D *h3,
                            const std::string &name,
                            double eta_lo,
                            double eta_hi)
    {
        if (!h3)
            return nullptr;

        TH3D *h3tmp = (TH3D *)h3->Clone(Form("%s__tmp_clone", name.c_str()));
        if (!h3tmp)
            return nullptr;
        h3tmp->SetDirectory(nullptr);

        const double eps = 1e-6;
        TAxis *zax = h3tmp->GetZaxis();
        const int nbz = zax->GetNbins();

        // IMPORTANT:
        // - If eta_lo/hi exceed axis bounds, include under/overflow bins (0 and nbins+1).
        //   This matters for "inclusive eta" comparisons vs standalone TH2 histograms.
        int zlo = 1;
        int zhi = nbz;
        if (eta_lo <= zax->GetXmin())
            zlo = 0; // include underflow
        else
            zlo = zax->FindBin(eta_lo + eps);

        if (eta_hi >= zax->GetXmax())
            zhi = nbz + 1; // include overflow
        else
            zhi = zax->FindBin(eta_hi - eps);

        if (zlo < 0)
            zlo = 0;
        if (zhi > nbz + 1)
            zhi = nbz + 1;
        if (zhi < zlo)
            std::swap(zlo, zhi);

        zax->SetRange(zlo, zhi);
        // ROOT quirk: Project3D("xy") swaps axes (X becomes original Y, Y becomes original X).
        // We want TH2 with X=original X and Y=original Y, so we must use "yx".
        TH2D *h2 = (TH2D *)h3tmp->Project3D("yx");
        if (h2)
        {
            h2->SetName(name.c_str());
            h2->SetDirectory(nullptr);
        }
        delete h3tmp;
        return h2;
    }
} // namespace

void plot_tower_timing(const std::string &inputfile = "../efficiencytool/results/truth_photon_tower_analysis.root",
                       const std::string &outputdir = "figures/")
{
    init_plot();
    gSystem->mkdir(outputdir.c_str(), true);

    std::cout << "==================================================" << std::endl;
    std::cout << "  Tower Timing Analysis Plotting" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "Input file: " << inputfile << std::endl;
    std::cout << "Output directory: " << outputdir << std::endl;
    std::cout << std::endl;

    // Open input file
    TFile *fin = TFile::Open(inputfile.c_str());
    if (!fin || fin->IsZombie())
    {
        std::cerr << "ERROR: Cannot open input file: " << inputfile << std::endl;
        return;
    }

    // Get histograms
    TH2D *h_tower_adc_vs_time = (TH2D *)fin->Get("h_tower_adc_vs_time");
    TH2D *h_tower_e_vs_time = (TH2D *)fin->Get("h_tower_e_vs_time");
    TH2D *h_leading_tower_adc_vs_time = (TH2D *)fin->Get("h_leading_tower_adc_vs_time");
    TH2D *h_leading_tower_e_vs_time = (TH2D *)fin->Get("h_leading_tower_e_vs_time");
    TH2D *h_subleading_tower_adc_vs_time = (TH2D *)fin->Get("h_subleading_tower_adc_vs_time");
    TH2D *h_subleading_tower_e_vs_time = (TH2D *)fin->Get("h_subleading_tower_e_vs_time");
    TH2D *h_nonleading_tower_adc_vs_time = (TH2D *)fin->Get("h_nonleading_tower_adc_vs_time");
    TH2D *h_nonleading_tower_e_vs_time = (TH2D *)fin->Get("h_nonleading_tower_e_vs_time");

    // New 3D histograms (time vs ADC/E vs cluster-eta tag)
    TH3D *h3_tower_adc_vs_time_vs_eta = (TH3D *)fin->Get("h3_tower_adc_vs_time_vs_eta");
    TH3D *h3_tower_e_vs_time_vs_eta = (TH3D *)fin->Get("h3_tower_e_vs_time_vs_eta");
    TH3D *h3_leading_tower_adc_vs_time_vs_eta = (TH3D *)fin->Get("h3_leading_tower_adc_vs_time_vs_eta");
    TH3D *h3_leading_tower_e_vs_time_vs_eta = (TH3D *)fin->Get("h3_leading_tower_e_vs_time_vs_eta");
    TH3D *h3_subleading_tower_adc_vs_time_vs_eta = (TH3D *)fin->Get("h3_subleading_tower_adc_vs_time_vs_eta");
    TH3D *h3_subleading_tower_e_vs_time_vs_eta = (TH3D *)fin->Get("h3_subleading_tower_e_vs_time_vs_eta");
    TH3D *h3_nonleading_tower_adc_vs_time_vs_eta = (TH3D *)fin->Get("h3_nonleading_tower_adc_vs_time_vs_eta");
    TH3D *h3_nonleading_tower_e_vs_time_vs_eta = (TH3D *)fin->Get("h3_nonleading_tower_e_vs_time_vs_eta");

    TH2D *h2_delta_t_leading_nonleading_vs_e = (TH2D *)fin->Get("h2_delta_t_leading_nonleading_vs_e");

    // Cluster-time histograms are now 3D: (cluster E, cluster time, eta)
    TH3D *h3_cluster_time_ew_all_vs_cluster_e_vs_eta = (TH3D *)fin->Get("h3_cluster_time_ew_all_vs_cluster_e_vs_eta");
    TH3D *h3_cluster_time_leading_vs_cluster_e_vs_eta = (TH3D *)fin->Get("h3_cluster_time_leading_vs_cluster_e_vs_eta");
    TH3D *h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta = (TH3D *)fin->Get("h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta");

    if (!h_tower_adc_vs_time || !h_tower_e_vs_time)
    {
        std::cerr << "ERROR: Cannot find required histograms" << std::endl;
        if (!h_tower_adc_vs_time)
            std::cerr << "  - Missing: h_tower_adc_vs_time" << std::endl;
        if (!h_tower_e_vs_time)
            std::cerr << "  - Missing: h_tower_e_vs_time" << std::endl;
        fin->Close();
        return;
    }

    std::cout << "Found histograms:" << std::endl;
    std::cout << "  - h_tower_adc_vs_time: " << h_tower_adc_vs_time->GetEntries() << " entries" << std::endl;
    std::cout << "  - h_tower_e_vs_time: " << h_tower_e_vs_time->GetEntries() << " entries" << std::endl;
    if (h3_tower_adc_vs_time_vs_eta)
        std::cout << "  - h3_tower_adc_vs_time_vs_eta: " << h3_tower_adc_vs_time_vs_eta->GetEntries() << " entries" << std::endl;
    else
        std::cout << "  - h3_tower_adc_vs_time_vs_eta: (missing)" << std::endl;
    if (h3_tower_e_vs_time_vs_eta)
        std::cout << "  - h3_tower_e_vs_time_vs_eta: " << h3_tower_e_vs_time_vs_eta->GetEntries() << " entries" << std::endl;
    else
        std::cout << "  - h3_tower_e_vs_time_vs_eta: (missing)" << std::endl;
    if (h_leading_tower_adc_vs_time)
        std::cout << "  - h_leading_tower_adc_vs_time: " << h_leading_tower_adc_vs_time->GetEntries() << " entries" << std::endl;
    else
        std::cout << "  - h_leading_tower_adc_vs_time: (missing)" << std::endl;
    if (h_leading_tower_e_vs_time)
        std::cout << "  - h_leading_tower_e_vs_time: " << h_leading_tower_e_vs_time->GetEntries() << " entries" << std::endl;
    else
        std::cout << "  - h_leading_tower_e_vs_time: (missing)" << std::endl;
    if (h_subleading_tower_adc_vs_time)
        std::cout << "  - h_subleading_tower_adc_vs_time: " << h_subleading_tower_adc_vs_time->GetEntries() << " entries" << std::endl;
    else
        std::cout << "  - h_subleading_tower_adc_vs_time: (missing)" << std::endl;
    if (h_subleading_tower_e_vs_time)
        std::cout << "  - h_subleading_tower_e_vs_time: " << h_subleading_tower_e_vs_time->GetEntries() << " entries" << std::endl;
    else
        std::cout << "  - h_subleading_tower_e_vs_time: (missing)" << std::endl;
    if (h_nonleading_tower_adc_vs_time)
        std::cout << "  - h_nonleading_tower_adc_vs_time: " << h_nonleading_tower_adc_vs_time->GetEntries() << " entries" << std::endl;
    else
        std::cout << "  - h_nonleading_tower_adc_vs_time: (missing)" << std::endl;
    if (h_nonleading_tower_e_vs_time)
        std::cout << "  - h_nonleading_tower_e_vs_time: " << h_nonleading_tower_e_vs_time->GetEntries() << " entries" << std::endl;
    else
        std::cout << "  - h_nonleading_tower_e_vs_time: (missing)" << std::endl;
    if (h2_delta_t_leading_nonleading_vs_e)
        std::cout << "  - h2_delta_t_leading_nonleading_vs_e: " << h2_delta_t_leading_nonleading_vs_e->GetEntries() << " entries" << std::endl;
    else
        std::cout << "  - h2_delta_t_leading_nonleading_vs_e: (missing)" << std::endl;
    if (h3_cluster_time_ew_all_vs_cluster_e_vs_eta)
        std::cout << "  - h3_cluster_time_ew_all_vs_cluster_e_vs_eta: " << h3_cluster_time_ew_all_vs_cluster_e_vs_eta->GetEntries() << " entries" << std::endl;
    else
        std::cout << "  - h3_cluster_time_ew_all_vs_cluster_e_vs_eta: (missing)" << std::endl;
    if (h3_cluster_time_leading_vs_cluster_e_vs_eta)
        std::cout << "  - h3_cluster_time_leading_vs_cluster_e_vs_eta: " << h3_cluster_time_leading_vs_cluster_e_vs_eta->GetEntries() << " entries" << std::endl;
    else
        std::cout << "  - h3_cluster_time_leading_vs_cluster_e_vs_eta: (missing)" << std::endl;
    if (h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta)
        std::cout << "  - h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta: " << h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta->GetEntries() << " entries" << std::endl;
    else
        std::cout << "  - h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta: (missing)" << std::endl;
    std::cout << std::endl;

    // ========================================================================
    // Part 1: Average Tower Time vs Tower Energy (Profile)
    // ========================================================================

    std::cout << "Creating average tower time vs energy profile..." << std::endl;

    TCanvas *c_time_vs_energy = new TCanvas("c_time_vs_energy", "", 800, 600);
    c_time_vs_energy->SetLeftMargin(0.14);
    c_time_vs_energy->SetRightMargin(0.05);
    c_time_vs_energy->SetBottomMargin(0.12);

    // Create profile: average time (X) as function of energy (Y)
    // Our histogram has time on X-axis, energy on Y-axis
    // ProfileY gives us average X (time) for each Y (energy) bin
    TProfile *prof_time_vs_e = h_tower_e_vs_time->ProfileY("prof_time_vs_e", 1, -1, "");
    TProfile *prof_leading_time_vs_e = nullptr;
    TProfile *prof_subleading_time_vs_e = nullptr;
    TProfile *prof_nonleading_time_vs_e = nullptr;
    if (h_leading_tower_e_vs_time)
    {
        prof_leading_time_vs_e = h_leading_tower_e_vs_time->ProfileY("prof_leading_time_vs_e", 1, -1, "");
    }
    if (h_subleading_tower_e_vs_time)
    {
        prof_subleading_time_vs_e = h_subleading_tower_e_vs_time->ProfileY("prof_subleading_time_vs_e", 1, -1, "");
    }
    if (h_nonleading_tower_e_vs_time)
    {
        prof_nonleading_time_vs_e = h_nonleading_tower_e_vs_time->ProfileY("prof_nonleading_time_vs_e", 1, -1, "");
    }

    prof_time_vs_e->SetLineColor(kBlue);
    prof_time_vs_e->SetLineWidth(3);
    prof_time_vs_e->SetMarkerColor(kBlue);
    prof_time_vs_e->SetMarkerStyle(20);
    prof_time_vs_e->SetMarkerSize(1.0);

    // Set axis titles and ranges
    prof_time_vs_e->SetTitle(";Tower Energy [GeV];Average Tower Time [ns]");
    prof_time_vs_e->GetXaxis()->SetTitleSize(0.05);
    prof_time_vs_e->GetXaxis()->SetLabelSize(0.045);
    prof_time_vs_e->GetXaxis()->SetTitleOffset(1.1);
    prof_time_vs_e->GetYaxis()->SetTitleSize(0.05);
    prof_time_vs_e->GetYaxis()->SetLabelSize(0.045);
    prof_time_vs_e->GetYaxis()->SetTitleOffset(1.3);
    prof_time_vs_e->GetXaxis()->SetRangeUser(0, 15);

    prof_time_vs_e->Draw();
    if (prof_leading_time_vs_e)
    {
        prof_leading_time_vs_e->SetLineColor(kGreen + 2);
        prof_leading_time_vs_e->SetLineWidth(3);
        prof_leading_time_vs_e->SetMarkerColor(kGreen + 2);
        prof_leading_time_vs_e->SetMarkerStyle(24);
        prof_leading_time_vs_e->SetMarkerSize(1.0);
        prof_leading_time_vs_e->Draw("same");
    }
    if (prof_subleading_time_vs_e)
    {
        prof_subleading_time_vs_e->SetLineColor(kMagenta + 2);
        prof_subleading_time_vs_e->SetLineWidth(3);
        prof_subleading_time_vs_e->SetMarkerColor(kMagenta + 2);
        prof_subleading_time_vs_e->SetMarkerStyle(26);
        prof_subleading_time_vs_e->SetMarkerSize(1.0);
        prof_subleading_time_vs_e->Draw("same");
    }
    if (prof_nonleading_time_vs_e)
    {
        prof_nonleading_time_vs_e->SetLineColor(kOrange + 1);
        prof_nonleading_time_vs_e->SetLineWidth(3);
        prof_nonleading_time_vs_e->SetMarkerColor(kOrange + 1);
        prof_nonleading_time_vs_e->SetMarkerStyle(25);
        prof_nonleading_time_vs_e->SetMarkerSize(1.0);
        prof_nonleading_time_vs_e->Draw("same");
    }

    // Add labels
    myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
    myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);
    myText(0.20, 0.78, 1, "(Direct + Frag)", 0.045);
    myMarkerLineText(0.55, 0.25, 1, kBlue, 20, kBlue, 1, "All Towers", 0.04, true);
    if (prof_leading_time_vs_e)
        myMarkerLineText(0.55, 0.20, 1, kGreen + 2, 24, kGreen + 2, 1, "Leading Tower", 0.04, true);
    if (prof_subleading_time_vs_e)
        myMarkerLineText(0.55, 0.15, 1, kMagenta + 2, 26, kMagenta + 2, 1, "Subleading Tower", 0.04, true);
    if (prof_nonleading_time_vs_e)
        myMarkerLineText(0.55, 0.10, 1, kOrange + 1, 25, kOrange + 1, 1, "Non-Leading Towers", 0.04, true);

    c_time_vs_energy->SaveAs(Form("%s/tower_energy_vs_time_profile.pdf", outputdir.c_str()));

    std::cout << "Created tower energy vs time profile plot" << std::endl;

    // ========================================================================
    // Part 1a: Average Tower Time vs Tower Energy (inclusive eta from 3D)
    // ========================================================================

    if (h3_tower_e_vs_time_vs_eta)
    {
        std::cout << "Creating <time>(E) profiles from 3D histograms (inclusive eta via Z-projection)..." << std::endl;

        // IMPORTANT: TH3::Project3D respects axis ranges. To guarantee "inclusive eta",
        // use the same safe pattern as eta-sliced profiles: clone -> set Z range -> Project3D("yx") -> ProfileY.
        const double eta_lo = -1.0e6;
        const double eta_hi = 1.0e6;
        TProfile *prof_all_from3d = make_mean_x_vs_y_in_eta_via_projection(h3_tower_e_vs_time_vs_eta,
                                                                           "prof_time_vs_e_from3d_inclusive_eta",
                                                                           eta_lo, eta_hi);
        TProfile *prof_lead_from3d = nullptr;
        TProfile *prof_sublead_from3d = nullptr;
        TProfile *prof_nonlead_from3d = nullptr;
        if (h3_leading_tower_e_vs_time_vs_eta)
        {
            prof_lead_from3d = make_mean_x_vs_y_in_eta_via_projection(h3_leading_tower_e_vs_time_vs_eta,
                                                                      "prof_leading_time_vs_e_from3d_inclusive_eta",
                                                                      eta_lo, eta_hi);
        }
        if (h3_subleading_tower_e_vs_time_vs_eta)
        {
            prof_sublead_from3d = make_mean_x_vs_y_in_eta_via_projection(h3_subleading_tower_e_vs_time_vs_eta,
                                                                         "prof_subleading_time_vs_e_from3d_inclusive_eta",
                                                                         eta_lo, eta_hi);
        }
        if (h3_nonleading_tower_e_vs_time_vs_eta)
        {
            prof_nonlead_from3d = make_mean_x_vs_y_in_eta_via_projection(h3_nonleading_tower_e_vs_time_vs_eta,
                                                                         "prof_nonleading_time_vs_e_from3d_inclusive_eta",
                                                                         eta_lo, eta_hi);
        }

        if (!prof_all_from3d)
        {
            std::cout << "WARNING: Could not build inclusive-eta <time>(E) profile from h3_tower_e_vs_time_vs_eta" << std::endl;
            delete prof_lead_from3d;
            delete prof_sublead_from3d;
            delete prof_nonlead_from3d;
        }
        else
        {
            // Create comparison canvas
            TCanvas *c_compare_3d = new TCanvas("c_compare_3d", "", 800, 600);
            c_compare_3d->SetLeftMargin(0.14);
            c_compare_3d->SetRightMargin(0.05);
            c_compare_3d->SetBottomMargin(0.12);

            prof_all_from3d->SetLineColor(kBlue);
            prof_all_from3d->SetLineWidth(3);
            prof_all_from3d->SetMarkerColor(kBlue);
            prof_all_from3d->SetMarkerStyle(20);
            prof_all_from3d->SetMarkerSize(1.0);
            prof_all_from3d->SetTitle(";Tower Energy [GeV];Average Tower Time [ns]");
            prof_all_from3d->GetXaxis()->SetTitleSize(0.05);
            prof_all_from3d->GetXaxis()->SetLabelSize(0.045);
            prof_all_from3d->GetXaxis()->SetTitleOffset(1.1);
            prof_all_from3d->GetYaxis()->SetTitleSize(0.05);
            prof_all_from3d->GetYaxis()->SetLabelSize(0.045);
            prof_all_from3d->GetYaxis()->SetTitleOffset(1.3);
            prof_all_from3d->GetXaxis()->SetRangeUser(0, 20);
            prof_all_from3d->GetYaxis()->SetRangeUser(0, 3);
            prof_all_from3d->Draw();

            if (prof_lead_from3d)
            {
                prof_lead_from3d->SetLineColor(kGreen + 2);
                prof_lead_from3d->SetLineWidth(3);
                prof_lead_from3d->SetMarkerColor(kGreen + 2);
                prof_lead_from3d->SetMarkerStyle(24);
                prof_lead_from3d->SetMarkerSize(1.0);
                prof_lead_from3d->Draw("same");
            }

            if (prof_sublead_from3d)
            {
                prof_sublead_from3d->SetLineColor(kMagenta + 2);
                prof_sublead_from3d->SetLineWidth(3);
                prof_sublead_from3d->SetMarkerColor(kMagenta + 2);
                prof_sublead_from3d->SetMarkerStyle(26);
                prof_sublead_from3d->SetMarkerSize(1.0);
                prof_sublead_from3d->Draw("same");
            }

            if (prof_nonlead_from3d)
            {
                prof_nonlead_from3d->SetLineColor(kOrange + 1);
                prof_nonlead_from3d->SetLineWidth(3);
                prof_nonlead_from3d->SetMarkerColor(kOrange + 1);
                prof_nonlead_from3d->SetMarkerStyle(25);
                prof_nonlead_from3d->SetMarkerSize(1.0);
                prof_nonlead_from3d->Draw("same");
            }

            myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
            myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);
            myText(0.20, 0.78, 1, "(Direct + Frag)", 0.045);
            myText(0.20, 0.73, 1, "Inclusive #eta (from 3D)", 0.045);
            myMarkerLineText(0.55, 0.25, 1, kBlue, 20, kBlue, 1, "All Towers", 0.04, true);
            if (prof_lead_from3d)
                myMarkerLineText(0.55, 0.20, 1, kGreen + 2, 24, kGreen + 2, 1, "Leading Tower", 0.04, true);
            if (prof_sublead_from3d)
                myMarkerLineText(0.55, 0.15, 1, kMagenta + 2, 26, kMagenta + 2, 1, "Subleading Tower", 0.04, true);
            if (prof_nonlead_from3d)
                myMarkerLineText(0.55, 0.10, 1, kOrange + 1, 25, kOrange + 1, 1, "Non-Leading Towers", 0.04, true);

            c_compare_3d->SaveAs(Form("%s/tower_energy_vs_time_profile_inclusive_eta_from3d.pdf", outputdir.c_str()));
            std::cout << "Created inclusive eta profile from 3D histogram" << std::endl;

            delete prof_all_from3d;
            delete prof_lead_from3d;
            delete prof_sublead_from3d;
            delete prof_nonlead_from3d;
            delete c_compare_3d;
        }
    }

    // ========================================================================
    // Part 1b: Average Tower Time vs Tower Energy in eta slices (uses 3D profile)
    // ========================================================================

    if (h3_tower_e_vs_time_vs_eta)
    {
        std::cout << "Creating <time>(E) profiles in eta slices using 3D profiling..." << std::endl;
        const auto etaBins = default_eta_bins();

        for (size_t ib = 0; ib < etaBins.size(); ++ib)
        {
            const auto &b = etaBins[ib];

            TCanvas *c = new TCanvas(Form("c_time_vs_energy_eta_%zu", ib), "", 800, 600);
            c->SetLeftMargin(0.14);
            c->SetRightMargin(0.05);
            c->SetBottomMargin(0.12);

            TH1D *h_frame = new TH1D(Form("h_frame_e_eta_%zu", ib), "", 100, 0, 15);
            h_frame->SetTitle(";Tower Energy [GeV];Average Tower Time [ns]");
            h_frame->GetXaxis()->SetTitleSize(0.05);
            h_frame->GetXaxis()->SetLabelSize(0.045);
            h_frame->GetXaxis()->SetTitleOffset(1.1);
            h_frame->GetYaxis()->SetTitleSize(0.05);
            h_frame->GetYaxis()->SetLabelSize(0.045);
            h_frame->GetYaxis()->SetTitleOffset(1.3);
            h_frame->GetYaxis()->SetRangeUser(0, 2);
            h_frame->Draw("axis");

            // Recommended ROOT method: select eta on Z, project to TH2(time,value), then ProfileY => <time>(value)
            TProfile *p_all = make_mean_x_vs_y_in_eta_via_projection(h3_tower_e_vs_time_vs_eta,
                                                                     Form("prof_time_vs_e_all_eta_%zu", ib),
                                                                     b.lo, b.hi);
            TProfile *p_lead = nullptr;
            TProfile *p_sublead = nullptr;
            TProfile *p_nonlead = nullptr;
            if (h3_leading_tower_e_vs_time_vs_eta)
            {
                p_lead = make_mean_x_vs_y_in_eta_via_projection(h3_leading_tower_e_vs_time_vs_eta,
                                                               Form("prof_time_vs_e_lead_eta_%zu", ib),
                                                               b.lo, b.hi);
            }
            if (h3_subleading_tower_e_vs_time_vs_eta)
            {
                p_sublead = make_mean_x_vs_y_in_eta_via_projection(h3_subleading_tower_e_vs_time_vs_eta,
                                                                  Form("prof_time_vs_e_sublead_eta_%zu", ib),
                                                                  b.lo, b.hi);
            }
            if (h3_nonleading_tower_e_vs_time_vs_eta)
            {
                p_nonlead = make_mean_x_vs_y_in_eta_via_projection(h3_nonleading_tower_e_vs_time_vs_eta,
                                                                  Form("prof_time_vs_e_nonlead_eta_%zu", ib),
                                                                  b.lo, b.hi);
            }

            if (p_all)
            {
                p_all->SetLineColor(kBlue);
                p_all->SetLineWidth(3);
                p_all->SetMarkerColor(kBlue);
                p_all->SetMarkerStyle(20);
                p_all->SetMarkerSize(1.0);
                p_all->Draw("same");
            }
            if (p_lead)
            {
                p_lead->SetLineColor(kGreen + 2);
                p_lead->SetLineWidth(3);
                p_lead->SetMarkerColor(kGreen + 2);
                p_lead->SetMarkerStyle(24);
                p_lead->SetMarkerSize(1.0);
                p_lead->Draw("same");
            }
            if (p_sublead)
            {
                p_sublead->SetLineColor(kMagenta + 2);
                p_sublead->SetLineWidth(3);
                p_sublead->SetMarkerColor(kMagenta + 2);
                p_sublead->SetMarkerStyle(26);
                p_sublead->SetMarkerSize(1.0);
                p_sublead->Draw("same");
            }
            if (p_nonlead)
            {
                p_nonlead->SetLineColor(kOrange + 1);
                p_nonlead->SetLineWidth(3);
                p_nonlead->SetMarkerColor(kOrange + 1);
                p_nonlead->SetMarkerStyle(25);
                p_nonlead->SetMarkerSize(1.0);
                p_nonlead->Draw("same");
            }

            myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
            myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);
            myText(0.20, 0.78, 1, "(Direct + Frag)", 0.045);
            myText(0.20, 0.73, 1, b.label.c_str(), 0.045);
            myMarkerLineText(0.55, 0.25, 1, kBlue, 20, kBlue, 1, "All Towers", 0.04, true);
            if (p_lead)
                myMarkerLineText(0.55, 0.20, 1, kGreen + 2, 24, kGreen + 2, 1, "Leading Tower", 0.04, true);
            if (p_sublead)
                myMarkerLineText(0.55, 0.15, 1, kMagenta + 2, 26, kMagenta + 2, 1, "Subleading Tower", 0.04, true);
            if (p_nonlead)
                myMarkerLineText(0.55, 0.10, 1, kOrange + 1, 25, kOrange + 1, 1, "Non-Leading Towers", 0.04, true);

            c->SaveAs(Form("%s/tower_energy_vs_time_profile_eta_%zu.pdf", outputdir.c_str(), ib));

            delete h_frame;
            delete c;
            delete p_all;
            delete p_lead;
            delete p_sublead;
            delete p_nonlead;
        }
    }
    else
    {
        std::cout << "WARNING: Missing h3_tower_e_vs_time_vs_eta; skipping eta-sliced <time>(E) profiles" << std::endl;
    }

    // ========================================================================
    // Part 1b: Cluster time (mean/RMS) vs cluster energy overlays
    // ========================================================================

    if (h3_cluster_time_ew_all_vs_cluster_e_vs_eta &&
        h3_cluster_time_leading_vs_cluster_e_vs_eta &&
        h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta)
    {
        std::cout << "Creating cluster time vs cluster energy (mean/RMS) overlay plots (inclusive eta)..." << std::endl;

        const int min_entries = 10;
        // Inclusive eta projection
        TH2D *h2_all = (TH2D *)h3_cluster_time_ew_all_vs_cluster_e_vs_eta->Project3D("yx");
        TH2D *h2_lead = (TH2D *)h3_cluster_time_leading_vs_cluster_e_vs_eta->Project3D("yx");
        TH2D *h2_nonlead = (TH2D *)h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta->Project3D("yx");
        if (h2_all)
            h2_all->SetName("h2_cluster_time_ew_all_vs_cluster_e__inclusive_eta");
        if (h2_lead)
            h2_lead->SetName("h2_cluster_time_leading_vs_cluster_e__inclusive_eta");
        if (h2_nonlead)
            h2_nonlead->SetName("h2_cluster_time_ew_nonleading_vs_cluster_e__inclusive_eta");

        MeanRmsGraphs g_all = make_mean_rms_vs_x(h2_all, "cluster_time_all", min_entries);
        MeanRmsGraphs g_lead = make_mean_rms_vs_x(h2_lead, "cluster_time_leading", min_entries);
        MeanRmsGraphs g_nonlead = make_mean_rms_vs_x(h2_nonlead, "cluster_time_nonleading", min_entries);

        // Mean plot
        TCanvas *c_cluster_time_mean = new TCanvas("c_cluster_time_mean", "", 800, 600);
        c_cluster_time_mean->SetLeftMargin(0.14);
        c_cluster_time_mean->SetRightMargin(0.05);
        c_cluster_time_mean->SetBottomMargin(0.12);

        TH1D *h_frame_cluster_mean = new TH1D("h_frame_cluster_mean", "", 100, 0, 50);
        h_frame_cluster_mean->SetTitle(";Cluster Energy [GeV];Mean Cluster Time [ns]");
        h_frame_cluster_mean->GetXaxis()->SetTitleSize(0.05);
        h_frame_cluster_mean->GetXaxis()->SetLabelSize(0.045);
        h_frame_cluster_mean->GetXaxis()->SetTitleOffset(1.1);
        h_frame_cluster_mean->GetYaxis()->SetTitleSize(0.05);
        h_frame_cluster_mean->GetYaxis()->SetLabelSize(0.045);
        h_frame_cluster_mean->GetYaxis()->SetTitleOffset(1.3);
        h_frame_cluster_mean->GetYaxis()->SetRangeUser(0, 3);
        h_frame_cluster_mean->Draw("axis");

        if (g_all.g_mean)
        {
            g_all.g_mean->SetLineColor(kBlue);
            g_all.g_mean->SetMarkerColor(kBlue);
            g_all.g_mean->SetMarkerStyle(20);
            g_all.g_mean->SetLineWidth(2);
            g_all.g_mean->Draw("P same");
        }
        if (g_lead.g_mean)
        {
            g_lead.g_mean->SetLineColor(kGreen + 2);
            g_lead.g_mean->SetMarkerColor(kGreen + 2);
            g_lead.g_mean->SetMarkerStyle(24);
            g_lead.g_mean->SetLineWidth(2);
            g_lead.g_mean->Draw("P same");
        }
        if (g_nonlead.g_mean)
        {
            g_nonlead.g_mean->SetLineColor(kOrange + 1);
            g_nonlead.g_mean->SetMarkerColor(kOrange + 1);
            g_nonlead.g_mean->SetMarkerStyle(25);
            g_nonlead.g_mean->SetLineWidth(2);
            g_nonlead.g_mean->Draw("P same");
        }

        myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
        myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);
        myText(0.20, 0.78, 1, "(Direct + Frag)", 0.045);
        myMarkerLineText(0.55, 0.25, 1, kBlue, 20, kBlue, 1, "E-weighted (all towers)", 0.04, true);
        myMarkerLineText(0.55, 0.20, 1, kGreen + 2, 24, kGreen + 2, 1, "Leading tower time", 0.04, true);
        myMarkerLineText(0.55, 0.15, 1, kOrange + 1, 25, kOrange + 1, 1, "E-weighted (non-leading)", 0.04, true);

        c_cluster_time_mean->SaveAs(Form("%s/cluster_time_mean_vs_cluster_energy.pdf", outputdir.c_str()));

        // RMS plot
        TCanvas *c_cluster_time_rms = new TCanvas("c_cluster_time_rms", "", 800, 600);
        c_cluster_time_rms->SetLeftMargin(0.14);
        c_cluster_time_rms->SetRightMargin(0.05);
        c_cluster_time_rms->SetBottomMargin(0.12);

        TH1D *h_frame_cluster_rms = new TH1D("h_frame_cluster_rms", "", 100, 0, 50);
        h_frame_cluster_rms->SetTitle(";Cluster Energy [GeV];RMS of Cluster Time [ns]");
        h_frame_cluster_rms->GetXaxis()->SetTitleSize(0.05);
        h_frame_cluster_rms->GetXaxis()->SetLabelSize(0.045);
        h_frame_cluster_rms->GetXaxis()->SetTitleOffset(1.1);
        h_frame_cluster_rms->GetYaxis()->SetTitleSize(0.05);
        h_frame_cluster_rms->GetYaxis()->SetLabelSize(0.045);
        h_frame_cluster_rms->GetYaxis()->SetTitleOffset(1.3);
        h_frame_cluster_rms->GetYaxis()->SetRangeUser(0, 3);
        h_frame_cluster_rms->Draw("axis");

        if (g_all.g_rms)
        {
            g_all.g_rms->SetLineColor(kBlue);
            g_all.g_rms->SetMarkerColor(kBlue);
            g_all.g_rms->SetMarkerStyle(20);
            g_all.g_rms->SetLineWidth(2);
            g_all.g_rms->Draw("P same");
        }
        if (g_lead.g_rms)
        {
            g_lead.g_rms->SetLineColor(kGreen + 2);
            g_lead.g_rms->SetMarkerColor(kGreen + 2);
            g_lead.g_rms->SetMarkerStyle(24);
            g_lead.g_rms->SetLineWidth(2);
            g_lead.g_rms->Draw("P same");
        }
        if (g_nonlead.g_rms)
        {
            g_nonlead.g_rms->SetLineColor(kOrange + 1);
            g_nonlead.g_rms->SetMarkerColor(kOrange + 1);
            g_nonlead.g_rms->SetMarkerStyle(25);
            g_nonlead.g_rms->SetLineWidth(2);
            g_nonlead.g_rms->Draw("P same");
        }

        myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
        myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);
        myText(0.20, 0.78, 1, "(Direct + Frag)", 0.045);
        myMarkerLineText(0.55, 0.25, 1, kBlue, 20, kBlue, 1, "E-weighted (all towers)", 0.04, true);
        myMarkerLineText(0.55, 0.20, 1, kGreen + 2, 24, kGreen + 2, 1, "Leading tower time", 0.04, true);
        myMarkerLineText(0.55, 0.15, 1, kOrange + 1, 25, kOrange + 1, 1, "E-weighted (non-leading)", 0.04, true);

        c_cluster_time_rms->SaveAs(Form("%s/cluster_time_rms_vs_cluster_energy.pdf", outputdir.c_str()));

        // Cleanup for this block
        delete h_frame_cluster_mean;
        delete h_frame_cluster_rms;
        delete c_cluster_time_mean;
        delete c_cluster_time_rms;
        delete g_all.g_mean;
        delete g_all.g_rms;
        delete g_lead.g_mean;
        delete g_lead.g_rms;
        delete g_nonlead.g_mean;
        delete g_nonlead.g_rms;

        std::cout << "Created cluster time mean/RMS vs cluster energy overlay plots" << std::endl;

        delete h2_all;
        delete h2_lead;
        delete h2_nonlead;
    }
    else
    {
        std::cout << "WARNING: Missing one or more h3_cluster_time_* histograms; skipping cluster-time plots" << std::endl;
    }

    // ========================================================================
    // Part 1c: Cluster time (mean/RMS) vs cluster energy in eta slices
    // ========================================================================

    if (h3_cluster_time_ew_all_vs_cluster_e_vs_eta &&
        h3_cluster_time_leading_vs_cluster_e_vs_eta &&
        h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta)
    {
        std::cout << "Creating cluster time (mean/RMS) vs cluster energy in eta slices..." << std::endl;
        const auto etaBins = default_eta_bins();
        const int min_entries = 10;

        for (size_t ib = 0; ib < etaBins.size(); ++ib)
        {
            const auto &b = etaBins[ib];
            TH2D *h2_all = project_xy_in_eta(h3_cluster_time_ew_all_vs_cluster_e_vs_eta,
                                             Form("h2_cluster_time_all_eta_%zu", ib),
                                             b.lo, b.hi);
            TH2D *h2_lead = project_xy_in_eta(h3_cluster_time_leading_vs_cluster_e_vs_eta,
                                              Form("h2_cluster_time_lead_eta_%zu", ib),
                                              b.lo, b.hi);
            TH2D *h2_nonlead = project_xy_in_eta(h3_cluster_time_ew_nonleading_vs_cluster_e_vs_eta,
                                                 Form("h2_cluster_time_nonlead_eta_%zu", ib),
                                                 b.lo, b.hi);

            MeanRmsGraphs g_all = make_mean_rms_vs_x(h2_all, Form("cluster_time_all_eta_%zu", ib), min_entries);
            MeanRmsGraphs g_lead = make_mean_rms_vs_x(h2_lead, Form("cluster_time_leading_eta_%zu", ib), min_entries);
            MeanRmsGraphs g_nonlead = make_mean_rms_vs_x(h2_nonlead, Form("cluster_time_nonleading_eta_%zu", ib), min_entries);

            // Mean
            TCanvas *c_mean = new TCanvas(Form("c_cluster_time_mean_eta_%zu", ib), "", 800, 600);
            c_mean->SetLeftMargin(0.14);
            c_mean->SetRightMargin(0.05);
            c_mean->SetBottomMargin(0.12);

            TH1D *h_frame_mean = new TH1D(Form("h_frame_cluster_mean_eta_%zu", ib), "", 100, 0, 50);
            h_frame_mean->SetTitle(";Cluster Energy [GeV];Mean Cluster Time [ns]");
            h_frame_mean->GetXaxis()->SetTitleSize(0.05);
            h_frame_mean->GetXaxis()->SetLabelSize(0.045);
            h_frame_mean->GetXaxis()->SetTitleOffset(1.1);
            h_frame_mean->GetYaxis()->SetTitleSize(0.05);
            h_frame_mean->GetYaxis()->SetLabelSize(0.045);
            h_frame_mean->GetYaxis()->SetTitleOffset(1.3);
            h_frame_mean->GetYaxis()->SetRangeUser(0, 3);
            h_frame_mean->Draw("axis");

            if (g_all.g_mean)
            {
                g_all.g_mean->SetLineColor(kBlue);
                g_all.g_mean->SetMarkerColor(kBlue);
                g_all.g_mean->SetMarkerStyle(20);
                g_all.g_mean->SetLineWidth(2);
                g_all.g_mean->Draw("P same");
            }
            if (g_lead.g_mean)
            {
                g_lead.g_mean->SetLineColor(kGreen + 2);
                g_lead.g_mean->SetMarkerColor(kGreen + 2);
                g_lead.g_mean->SetMarkerStyle(24);
                g_lead.g_mean->SetLineWidth(2);
                g_lead.g_mean->Draw("P same");
            }
            if (g_nonlead.g_mean)
            {
                g_nonlead.g_mean->SetLineColor(kOrange + 1);
                g_nonlead.g_mean->SetMarkerColor(kOrange + 1);
                g_nonlead.g_mean->SetMarkerStyle(25);
                g_nonlead.g_mean->SetLineWidth(2);
                g_nonlead.g_mean->Draw("P same");
            }

            myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
            myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);
            myText(0.20, 0.78, 1, "(Direct + Frag)", 0.045);
            myText(0.20, 0.73, 1, b.label.c_str(), 0.045);
            myMarkerLineText(0.55, 0.25, 1, kBlue, 20, kBlue, 1, "E-weighted (all towers)", 0.04, true);
            myMarkerLineText(0.55, 0.20, 1, kGreen + 2, 24, kGreen + 2, 1, "Leading tower time", 0.04, true);
            myMarkerLineText(0.55, 0.15, 1, kOrange + 1, 25, kOrange + 1, 1, "E-weighted (non-leading)", 0.04, true);

            c_mean->SaveAs(Form("%s/cluster_time_mean_vs_cluster_energy_eta_%zu.pdf", outputdir.c_str(), ib));

            // RMS
            TCanvas *c_rms = new TCanvas(Form("c_cluster_time_rms_eta_%zu", ib), "", 800, 600);
            c_rms->SetLeftMargin(0.14);
            c_rms->SetRightMargin(0.05);
            c_rms->SetBottomMargin(0.12);

            TH1D *h_frame_rms = new TH1D(Form("h_frame_cluster_rms_eta_%zu", ib), "", 100, 0, 50);
            h_frame_rms->SetTitle(";Cluster Energy [GeV];RMS of Cluster Time [ns]");
            h_frame_rms->GetXaxis()->SetTitleSize(0.05);
            h_frame_rms->GetXaxis()->SetLabelSize(0.045);
            h_frame_rms->GetXaxis()->SetTitleOffset(1.1);
            h_frame_rms->GetYaxis()->SetTitleSize(0.05);
            h_frame_rms->GetYaxis()->SetLabelSize(0.045);
            h_frame_rms->GetYaxis()->SetTitleOffset(1.3);
            h_frame_rms->GetYaxis()->SetRangeUser(0, 3);
            h_frame_rms->Draw("axis");

            if (g_all.g_rms)
            {
                g_all.g_rms->SetLineColor(kBlue);
                g_all.g_rms->SetMarkerColor(kBlue);
                g_all.g_rms->SetMarkerStyle(20);
                g_all.g_rms->SetLineWidth(2);
                g_all.g_rms->Draw("P same");
            }
            if (g_lead.g_rms)
            {
                g_lead.g_rms->SetLineColor(kGreen + 2);
                g_lead.g_rms->SetMarkerColor(kGreen + 2);
                g_lead.g_rms->SetMarkerStyle(24);
                g_lead.g_rms->SetLineWidth(2);
                g_lead.g_rms->Draw("P same");
            }
            if (g_nonlead.g_rms)
            {
                g_nonlead.g_rms->SetLineColor(kOrange + 1);
                g_nonlead.g_rms->SetMarkerColor(kOrange + 1);
                g_nonlead.g_rms->SetMarkerStyle(25);
                g_nonlead.g_rms->SetLineWidth(2);
                g_nonlead.g_rms->Draw("P same");
            }

            myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
            myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);
            myText(0.20, 0.78, 1, "(Direct + Frag)", 0.045);
            myText(0.20, 0.73, 1, b.label.c_str(), 0.045);
            myMarkerLineText(0.55, 0.25, 1, kBlue, 20, kBlue, 1, "E-weighted (all towers)", 0.04, true);
            myMarkerLineText(0.55, 0.20, 1, kGreen + 2, 24, kGreen + 2, 1, "Leading tower time", 0.04, true);
            myMarkerLineText(0.55, 0.15, 1, kOrange + 1, 25, kOrange + 1, 1, "E-weighted (non-leading)", 0.04, true);

            c_rms->SaveAs(Form("%s/cluster_time_rms_vs_cluster_energy_eta_%zu.pdf", outputdir.c_str(), ib));

            // Cleanup
            delete h2_all;
            delete h2_lead;
            delete h2_nonlead;
            delete h_frame_mean;
            delete h_frame_rms;
            delete c_mean;
            delete c_rms;
            delete g_all.g_mean;
            delete g_all.g_rms;
            delete g_lead.g_mean;
            delete g_lead.g_rms;
            delete g_nonlead.g_mean;
            delete g_nonlead.g_rms;
        }
    }

    // ========================================================================
    // Part 2: Average Tower Time vs Tower ADC (Profile)
    // ========================================================================

    std::cout << "Creating average tower time vs ADC profile..." << std::endl;

    TCanvas *c_time_vs_adc = new TCanvas("c_time_vs_adc", "", 800, 600);
    c_time_vs_adc->SetLeftMargin(0.14);
    c_time_vs_adc->SetRightMargin(0.05);
    c_time_vs_adc->SetBottomMargin(0.12);

    // Create profile: average time (X) as function of ADC (Y)
    // Our histogram has time on X-axis, ADC on Y-axis
    // ProfileY gives us average X (time) for each Y (ADC) bin
    TProfile *prof_time_vs_adc = h_tower_adc_vs_time->ProfileY("prof_time_vs_adc", 1, -1, "");
    TProfile *prof_leading_time_vs_adc = nullptr;
    TProfile *prof_subleading_time_vs_adc = nullptr;
    TProfile *prof_nonleading_time_vs_adc = nullptr;
    if (h_leading_tower_adc_vs_time)
    {
        prof_leading_time_vs_adc = h_leading_tower_adc_vs_time->ProfileY("prof_leading_time_vs_adc", 1, -1, "");
    }
    if (h_subleading_tower_adc_vs_time)
    {
        prof_subleading_time_vs_adc = h_subleading_tower_adc_vs_time->ProfileY("prof_subleading_time_vs_adc", 1, -1, "");
    }
    if (h_nonleading_tower_adc_vs_time)
    {
        prof_nonleading_time_vs_adc = h_nonleading_tower_adc_vs_time->ProfileY("prof_nonleading_time_vs_adc", 1, -1, "");
    }

    prof_time_vs_adc->SetLineColor(kRed);
    prof_time_vs_adc->SetLineWidth(3);
    prof_time_vs_adc->SetMarkerColor(kRed);
    prof_time_vs_adc->SetMarkerStyle(21);
    prof_time_vs_adc->SetMarkerSize(1.0);

    // Set axis titles and ranges
    prof_time_vs_adc->SetTitle(";Tower ADC;Average Tower Time [ns]");
    prof_time_vs_adc->GetXaxis()->SetTitleSize(0.05);
    prof_time_vs_adc->GetXaxis()->SetLabelSize(0.045);
    prof_time_vs_adc->GetXaxis()->SetTitleOffset(1.1);
    prof_time_vs_adc->GetYaxis()->SetTitleSize(0.05);
    prof_time_vs_adc->GetYaxis()->SetLabelSize(0.045);
    prof_time_vs_adc->GetYaxis()->SetTitleOffset(1.3);
    prof_time_vs_adc->GetXaxis()->SetRangeUser(0, 15000);
    prof_time_vs_adc->GetYaxis()->SetRangeUser(0, 3);

    prof_time_vs_adc->Draw();
    if (prof_leading_time_vs_adc)
    {
        prof_leading_time_vs_adc->SetLineColor(kMagenta + 2);
        prof_leading_time_vs_adc->SetLineWidth(3);
        prof_leading_time_vs_adc->SetMarkerColor(kMagenta + 2);
        prof_leading_time_vs_adc->SetMarkerStyle(25);
        prof_leading_time_vs_adc->SetMarkerSize(1.0);
        prof_leading_time_vs_adc->Draw("same");

        // Write leading-tower <time>(ADC) profile to a small ROOT file for later corrections
        // (used by efficiencytool/plot_cluster_time.C when using leading tower time)
        const std::string prof_outfile = Form("%s/leading_tower_time_vs_adc_profile.root", outputdir.c_str());
        TFile *fprof = new TFile(prof_outfile.c_str(), "RECREATE");
        if (fprof && !fprof->IsZombie())
        {
            prof_leading_time_vs_adc->Write("prof_leading_time_vs_adc");
            fprof->Close();
            std::cout << "Wrote leading tower time vs ADC profile to: " << prof_outfile << std::endl;
        }
        else
        {
            std::cout << "WARNING: Could not create profile output file: " << prof_outfile << std::endl;
        }
        delete fprof;
    }
    else
    {
        std::cout << "WARNING: Missing h_leading_tower_adc_vs_time; skipping leading tower time vs ADC profile" << std::endl;
    }
    if (prof_subleading_time_vs_adc)
    {
        prof_subleading_time_vs_adc->SetLineColor(kGreen + 2);
        prof_subleading_time_vs_adc->SetLineWidth(3);
        prof_subleading_time_vs_adc->SetMarkerColor(kGreen + 2);
        prof_subleading_time_vs_adc->SetMarkerStyle(24);
        prof_subleading_time_vs_adc->SetMarkerSize(1.0);
        prof_subleading_time_vs_adc->Draw("same");
    }
    if (prof_nonleading_time_vs_adc)
    {
        prof_nonleading_time_vs_adc->SetLineColor(kOrange + 1);
        prof_nonleading_time_vs_adc->SetLineWidth(3);
        prof_nonleading_time_vs_adc->SetMarkerColor(kOrange + 1);
        prof_nonleading_time_vs_adc->SetMarkerStyle(26);
        prof_nonleading_time_vs_adc->SetMarkerSize(1.0);
        prof_nonleading_time_vs_adc->Draw("same");
    }

    // Add labels
    myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
    myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);
    myText(0.20, 0.78, 1, "(Direct + Frag)", 0.045);
    myMarkerLineText(0.55, 0.25, 1, kRed, 21, kRed, 1, "All Towers", 0.04, true);
    if (prof_leading_time_vs_adc)
        myMarkerLineText(0.55, 0.20, 1, kMagenta + 2, 25, kMagenta + 2, 1, "Leading Tower", 0.04, true);
    if (prof_subleading_time_vs_adc)
        myMarkerLineText(0.55, 0.15, 1, kGreen + 2, 24, kGreen + 2, 1, "Subleading Tower", 0.04, true);
    if (prof_nonleading_time_vs_adc)
        myMarkerLineText(0.55, 0.10, 1, kOrange + 1, 26, kOrange + 1, 1, "Non-Leading Towers", 0.04, true);

    c_time_vs_adc->SaveAs(Form("%s/tower_adc_vs_time_profile.pdf", outputdir.c_str()));

    std::cout << "Created tower ADC vs time profile plot" << std::endl;

    // ========================================================================
    // Part 2a: Average Tower Time vs Tower ADC (inclusive eta from 3D)
    // ========================================================================

    if (h3_tower_adc_vs_time_vs_eta)
    {
        std::cout << "Creating <time>(ADC) profiles from 3D histograms (inclusive eta)..." << std::endl;

        const double eta_lo = -1.0e6;
        const double eta_hi = 1.0e6;
        TProfile *prof_all_adc_from3d = make_mean_x_vs_y_in_eta_via_projection(h3_tower_adc_vs_time_vs_eta,
                                                                               "prof_time_vs_adc_from3d_inclusive_eta",
                                                                               eta_lo, eta_hi);
        TProfile *prof_lead_adc_from3d = nullptr;
        TProfile *prof_sublead_adc_from3d = nullptr;
        TProfile *prof_nonlead_adc_from3d = nullptr;
        if (h3_leading_tower_adc_vs_time_vs_eta)
        {
            prof_lead_adc_from3d = make_mean_x_vs_y_in_eta_via_projection(h3_leading_tower_adc_vs_time_vs_eta,
                                                                          "prof_leading_time_vs_adc_from3d_inclusive_eta",
                                                                          eta_lo, eta_hi);
        }
        if (h3_subleading_tower_adc_vs_time_vs_eta)
        {
            prof_sublead_adc_from3d = make_mean_x_vs_y_in_eta_via_projection(h3_subleading_tower_adc_vs_time_vs_eta,
                                                                             "prof_subleading_time_vs_adc_from3d_inclusive_eta",
                                                                             eta_lo, eta_hi);
        }
        if (h3_nonleading_tower_adc_vs_time_vs_eta)
        {
            prof_nonlead_adc_from3d = make_mean_x_vs_y_in_eta_via_projection(h3_nonleading_tower_adc_vs_time_vs_eta,
                                                                             "prof_nonleading_time_vs_adc_from3d_inclusive_eta",
                                                                             eta_lo, eta_hi);
        }

        if (!prof_all_adc_from3d)
        {
            std::cout << "WARNING: Could not build inclusive-eta <time>(ADC) profile from h3_tower_adc_vs_time_vs_eta" << std::endl;
            delete prof_lead_adc_from3d;
            delete prof_sublead_adc_from3d;
            delete prof_nonlead_adc_from3d;
        }
        else
        {
            // Create comparison canvas
            TCanvas *c_compare_adc_3d = new TCanvas("c_compare_adc_3d", "", 800, 600);
            c_compare_adc_3d->SetLeftMargin(0.14);
            c_compare_adc_3d->SetRightMargin(0.05);
            c_compare_adc_3d->SetBottomMargin(0.12);

            prof_all_adc_from3d->SetLineColor(kRed);
            prof_all_adc_from3d->SetLineWidth(3);
            prof_all_adc_from3d->SetMarkerColor(kRed);
            prof_all_adc_from3d->SetMarkerStyle(21);
            prof_all_adc_from3d->SetMarkerSize(1.0);
            prof_all_adc_from3d->SetTitle(";Tower ADC;Average Tower Time [ns]");
            prof_all_adc_from3d->GetXaxis()->SetTitleSize(0.05);
            prof_all_adc_from3d->GetXaxis()->SetLabelSize(0.045);
            prof_all_adc_from3d->GetXaxis()->SetTitleOffset(1.1);
            prof_all_adc_from3d->GetYaxis()->SetTitleSize(0.05);
            prof_all_adc_from3d->GetYaxis()->SetLabelSize(0.045);
            prof_all_adc_from3d->GetYaxis()->SetTitleOffset(1.3);
            prof_all_adc_from3d->GetXaxis()->SetRangeUser(0, 15000);
            prof_all_adc_from3d->GetYaxis()->SetRangeUser(0, 3);
            prof_all_adc_from3d->Draw();

            if (prof_lead_adc_from3d)
            {
                prof_lead_adc_from3d->SetLineColor(kMagenta + 2);
                prof_lead_adc_from3d->SetLineWidth(3);
                prof_lead_adc_from3d->SetMarkerColor(kMagenta + 2);
                prof_lead_adc_from3d->SetMarkerStyle(25);
                prof_lead_adc_from3d->SetMarkerSize(1.0);
                prof_lead_adc_from3d->Draw("same");
            }

            if (prof_sublead_adc_from3d)
            {
                prof_sublead_adc_from3d->SetLineColor(kGreen + 2);
                prof_sublead_adc_from3d->SetLineWidth(3);
                prof_sublead_adc_from3d->SetMarkerColor(kGreen + 2);
                prof_sublead_adc_from3d->SetMarkerStyle(24);
                prof_sublead_adc_from3d->SetMarkerSize(1.0);
                prof_sublead_adc_from3d->Draw("same");
            }

            if (prof_nonlead_adc_from3d)
            {
                prof_nonlead_adc_from3d->SetLineColor(kOrange + 1);
                prof_nonlead_adc_from3d->SetLineWidth(3);
                prof_nonlead_adc_from3d->SetMarkerColor(kOrange + 1);
                prof_nonlead_adc_from3d->SetMarkerStyle(26);
                prof_nonlead_adc_from3d->SetMarkerSize(1.0);
                prof_nonlead_adc_from3d->Draw("same");
            }

            myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
            myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);
            myText(0.20, 0.78, 1, "(Direct + Frag)", 0.045);
            myText(0.20, 0.73, 1, "Inclusive #eta (from 3D)", 0.045);
            myMarkerLineText(0.55, 0.25, 1, kRed, 21, kRed, 1, "All Towers", 0.04, true);
            if (prof_lead_adc_from3d)
                myMarkerLineText(0.55, 0.20, 1, kMagenta + 2, 25, kMagenta + 2, 1, "Leading Tower", 0.04, true);
            if (prof_sublead_adc_from3d)
                myMarkerLineText(0.55, 0.15, 1, kGreen + 2, 24, kGreen + 2, 1, "Subleading Tower", 0.04, true);
            if (prof_nonlead_adc_from3d)
                myMarkerLineText(0.55, 0.10, 1, kOrange + 1, 26, kOrange + 1, 1, "Non-Leading Towers", 0.04, true);

            c_compare_adc_3d->SaveAs(Form("%s/tower_adc_vs_time_profile_inclusive_eta_from3d.pdf", outputdir.c_str()));
            std::cout << "Created inclusive eta ADC profile from 3D histogram" << std::endl;

            delete prof_all_adc_from3d;
            delete prof_lead_adc_from3d;
            delete prof_sublead_adc_from3d;
            delete prof_nonlead_adc_from3d;
            delete c_compare_adc_3d;
        }
    }

    // ========================================================================
    // Part 2b: Average Tower Time vs Tower ADC in eta slices (uses 3D profile)
    // ========================================================================

    if (h3_tower_adc_vs_time_vs_eta)
    {
        std::cout << "Creating <time>(ADC) profiles in eta slices using 3D profiling..." << std::endl;
        const auto etaBins = default_eta_bins();

        for (size_t ib = 0; ib < etaBins.size(); ++ib)
        {
            const auto &b = etaBins[ib];

            TCanvas *c = new TCanvas(Form("c_time_vs_adc_eta_%zu", ib), "", 800, 600);
            c->SetLeftMargin(0.14);
            c->SetRightMargin(0.05);
            c->SetBottomMargin(0.12);

            TH1D *h_frame = new TH1D(Form("h_frame_adc_eta_%zu", ib), "", 100, 0, 15000);
            h_frame->SetTitle(";Tower ADC;Average Tower Time [ns]");
            h_frame->GetXaxis()->SetTitleSize(0.05);
            h_frame->GetXaxis()->SetLabelSize(0.045);
            h_frame->GetXaxis()->SetTitleOffset(1.1);
            h_frame->GetYaxis()->SetTitleSize(0.05);
            h_frame->GetYaxis()->SetLabelSize(0.045);
            h_frame->GetYaxis()->SetTitleOffset(1.3);
            h_frame->GetYaxis()->SetRangeUser(0, 3);
            h_frame->Draw("axis");

            // Recommended ROOT method: select eta on Z, project to TH2(time,value), then ProfileY => <time>(value)
            TProfile *p_all = make_mean_x_vs_y_in_eta_via_projection(h3_tower_adc_vs_time_vs_eta,
                                                                     Form("prof_time_vs_adc_all_eta_%zu", ib),
                                                                     b.lo, b.hi);
            TProfile *p_lead = nullptr;
            TProfile *p_sublead = nullptr;
            TProfile *p_nonlead = nullptr;
            if (h3_leading_tower_adc_vs_time_vs_eta)
            {
                p_lead = make_mean_x_vs_y_in_eta_via_projection(h3_leading_tower_adc_vs_time_vs_eta,
                                                               Form("prof_time_vs_adc_lead_eta_%zu", ib),
                                                               b.lo, b.hi);
            }
            if (h3_subleading_tower_adc_vs_time_vs_eta)
            {
                p_sublead = make_mean_x_vs_y_in_eta_via_projection(h3_subleading_tower_adc_vs_time_vs_eta,
                                                                  Form("prof_time_vs_adc_sublead_eta_%zu", ib),
                                                                  b.lo, b.hi);
            }
            if (h3_nonleading_tower_adc_vs_time_vs_eta)
            {
                p_nonlead = make_mean_x_vs_y_in_eta_via_projection(h3_nonleading_tower_adc_vs_time_vs_eta,
                                                                  Form("prof_time_vs_adc_nonlead_eta_%zu", ib),
                                                                  b.lo, b.hi);
            }

            if (p_all)
            {
                p_all->SetLineColor(kRed);
                p_all->SetLineWidth(3);
                p_all->SetMarkerColor(kRed);
                p_all->SetMarkerStyle(21);
                p_all->SetMarkerSize(1.0);
                p_all->Draw("same");
            }
            if (p_lead)
            {
                p_lead->SetLineColor(kMagenta + 2);
                p_lead->SetLineWidth(3);
                p_lead->SetMarkerColor(kMagenta + 2);
                p_lead->SetMarkerStyle(25);
                p_lead->SetMarkerSize(1.0);
                p_lead->Draw("same");
            }
            if (p_sublead)
            {
                p_sublead->SetLineColor(kGreen + 2);
                p_sublead->SetLineWidth(3);
                p_sublead->SetMarkerColor(kGreen + 2);
                p_sublead->SetMarkerStyle(24);
                p_sublead->SetMarkerSize(1.0);
                p_sublead->Draw("same");
            }
            if (p_nonlead)
            {
                p_nonlead->SetLineColor(kOrange + 1);
                p_nonlead->SetLineWidth(3);
                p_nonlead->SetMarkerColor(kOrange + 1);
                p_nonlead->SetMarkerStyle(26);
                p_nonlead->SetMarkerSize(1.0);
                p_nonlead->Draw("same");
            }

            myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
            myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);
            myText(0.20, 0.78, 1, "(Direct + Frag)", 0.045);
            myText(0.20, 0.73, 1, b.label.c_str(), 0.045);
            myMarkerLineText(0.55, 0.25, 1, kRed, 21, kRed, 1, "All Towers", 0.04, true);
            if (p_lead)
                myMarkerLineText(0.55, 0.20, 1, kMagenta + 2, 25, kMagenta + 2, 1, "Leading Tower", 0.04, true);
            if (p_sublead)
                myMarkerLineText(0.55, 0.15, 1, kGreen + 2, 24, kGreen + 2, 1, "Subleading Tower", 0.04, true);
            if (p_nonlead)
                myMarkerLineText(0.55, 0.10, 1, kOrange + 1, 26, kOrange + 1, 1, "Non-Leading Towers", 0.04, true);

            c->SaveAs(Form("%s/tower_adc_vs_time_profile_eta_%zu.pdf", outputdir.c_str(), ib));

            delete h_frame;
            delete c;
            delete p_all;
            delete p_lead;
            delete p_sublead;
            delete p_nonlead;
        }
    }
    else
    {
        std::cout << "WARNING: Missing h3_tower_adc_vs_time_vs_eta; skipping eta-sliced <time>(ADC) profiles" << std::endl;
    }

    // ========================================================================
    // Part 3: 2D Histograms with profiles overlaid
    // ========================================================================

    std::cout << "Creating 2D histogram plots..." << std::endl;

    // Energy vs Time 2D
    TCanvas *c_e_vs_time_2d = new TCanvas("c_e_vs_time_2d", "", 900, 700);
    c_e_vs_time_2d->SetLeftMargin(0.12);
    c_e_vs_time_2d->SetRightMargin(0.15);
    c_e_vs_time_2d->SetBottomMargin(0.12);
    c_e_vs_time_2d->SetLogz();

    h_tower_e_vs_time->SetTitle(";Tower Time [ns];Tower Energy [GeV]");
    h_tower_e_vs_time->GetXaxis()->SetTitleSize(0.045);
    h_tower_e_vs_time->GetYaxis()->SetTitleSize(0.045);
    h_tower_e_vs_time->GetXaxis()->SetTitleOffset(1.1);
    h_tower_e_vs_time->GetYaxis()->SetTitleOffset(1.3);
    h_tower_e_vs_time->GetXaxis()->SetRangeUser(-10, 10);
    h_tower_e_vs_time->GetYaxis()->SetRangeUser(0, 3);
    h_tower_e_vs_time->Draw("colz");

    myText(0.18, 0.88, 1, strleg1.c_str(), 0.04);
    myText(0.18, 0.83, 1, "Truth-Matched Photons", 0.04);
    myText(0.18, 0.78, 1, "(Direct + Frag)", 0.04);

    c_e_vs_time_2d->SaveAs(Form("%s/tower_energy_vs_time_2d.pdf", outputdir.c_str()));

    std::cout << "Created tower energy vs time 2D plot" << std::endl;

    // ADC vs Time 2D
    TCanvas *c_adc_vs_time_2d = new TCanvas("c_adc_vs_time_2d", "", 900, 700);
    c_adc_vs_time_2d->SetLeftMargin(0.12);
    c_adc_vs_time_2d->SetRightMargin(0.15);
    c_adc_vs_time_2d->SetBottomMargin(0.12);
    c_adc_vs_time_2d->SetLogz();

    h_tower_adc_vs_time->SetTitle(";Tower Time [ns];Tower ADC");
    h_tower_adc_vs_time->GetXaxis()->SetTitleSize(0.045);
    h_tower_adc_vs_time->GetYaxis()->SetTitleSize(0.045);
    h_tower_adc_vs_time->GetXaxis()->SetTitleOffset(1.1);
    h_tower_adc_vs_time->GetYaxis()->SetTitleOffset(1.3);
    h_tower_adc_vs_time->GetXaxis()->SetRangeUser(-10, 10);
    h_tower_adc_vs_time->GetYaxis()->SetRangeUser(0, 5000);
    h_tower_adc_vs_time->Draw("colz");

    myText(0.18, 0.88, 1, strleg1.c_str(), 0.04);
    myText(0.18, 0.83, 1, "Truth-Matched Photons", 0.04);
    myText(0.18, 0.78, 1, "(Direct + Frag)", 0.04);

    c_adc_vs_time_2d->SaveAs(Form("%s/tower_adc_vs_time_2d.pdf", outputdir.c_str()));

    std::cout << "Created tower ADC vs time 2D plot" << std::endl;

    // ========================================================================
    // Part 3a: 2D histograms in eta slices (time vs energy / time vs ADC)
    // ========================================================================

    if (h3_tower_e_vs_time_vs_eta && h3_tower_adc_vs_time_vs_eta)
    {
        std::cout << "Creating 2D time-vs-(E/ADC) plots in eta slices..." << std::endl;
        const auto etaBins = default_eta_bins();

        for (size_t ib = 0; ib < etaBins.size(); ++ib)
        {
            const auto &b = etaBins[ib];

            // Energy vs Time 2D in eta bin
            TH2D *h2_e = project_xy_in_eta(h3_tower_e_vs_time_vs_eta,
                                           Form("h2_tower_e_vs_time_eta_%zu", ib),
                                           b.lo, b.hi);
            if (h2_e)
            {
                TCanvas *c = new TCanvas(Form("c_e_vs_time_2d_eta_%zu", ib), "", 900, 700);
                c->SetLeftMargin(0.12);
                c->SetRightMargin(0.15);
                c->SetBottomMargin(0.12);
                c->SetLogz();

                h2_e->SetTitle(";Tower Time [ns];Tower Energy [GeV]");
                h2_e->GetXaxis()->SetTitleSize(0.045);
                h2_e->GetYaxis()->SetTitleSize(0.045);
                h2_e->GetXaxis()->SetTitleOffset(1.1);
                h2_e->GetYaxis()->SetTitleOffset(1.3);
                h2_e->GetXaxis()->SetRangeUser(-10, 10);
                h2_e->GetYaxis()->SetRangeUser(0, 3);
                h2_e->Draw("colz");

                myText(0.18, 0.88, 1, strleg1.c_str(), 0.04);
                myText(0.18, 0.83, 1, "Truth-Matched Photons", 0.04);
                myText(0.18, 0.78, 1, "(Direct + Frag)", 0.04);
                myText(0.18, 0.73, 1, b.label.c_str(), 0.04);

                c->SaveAs(Form("%s/tower_energy_vs_time_2d_eta_%zu.pdf", outputdir.c_str(), ib));
                delete c;
                delete h2_e;
            }

            // ADC vs Time 2D in eta bin
            TH2D *h2_adc = project_xy_in_eta(h3_tower_adc_vs_time_vs_eta,
                                             Form("h2_tower_adc_vs_time_eta_%zu", ib),
                                             b.lo, b.hi);
            if (h2_adc)
            {
                TCanvas *c = new TCanvas(Form("c_adc_vs_time_2d_eta_%zu", ib), "", 900, 700);
                c->SetLeftMargin(0.12);
                c->SetRightMargin(0.15);
                c->SetBottomMargin(0.12);
                c->SetLogz();

                h2_adc->SetTitle(";Tower Time [ns];Tower ADC");
                h2_adc->GetXaxis()->SetTitleSize(0.045);
                h2_adc->GetYaxis()->SetTitleSize(0.045);
                h2_adc->GetXaxis()->SetTitleOffset(1.1);
                h2_adc->GetYaxis()->SetTitleOffset(1.3);
                h2_adc->GetXaxis()->SetRangeUser(-10, 10);
                h2_adc->GetYaxis()->SetRangeUser(0, 5000);
                h2_adc->Draw("colz");

                myText(0.18, 0.88, 1, strleg1.c_str(), 0.04);
                myText(0.18, 0.83, 1, "Truth-Matched Photons", 0.04);
                myText(0.18, 0.78, 1, "(Direct + Frag)", 0.04);
                myText(0.18, 0.73, 1, b.label.c_str(), 0.04);

                c->SaveAs(Form("%s/tower_adc_vs_time_2d_eta_%zu.pdf", outputdir.c_str(), ib));
                delete c;
                delete h2_adc;
            }
        }
    }

    // ========================================================================
    // Part 4: Time distribution projections for different energy/ADC ranges
    // ========================================================================

    std::cout << "Creating time distribution plots for different energy ranges..." << std::endl;

    // Define energy ranges
    const int n_e_ranges = 4;
    float e_mins[n_e_ranges] = {0.1, 0.5, 1.0, 2.0};
    float e_maxs[n_e_ranges] = {0.5, 1.0, 2.0, 5.0};
    std::string e_labels[n_e_ranges] = {
        "0.1 < E < 0.5 GeV",
        "0.5 < E < 1.0 GeV",
        "1.0 < E < 2.0 GeV",
        "2.0 < E < 5.0 GeV"
    };

    TCanvas *c_time_e_ranges = new TCanvas("c_time_e_ranges", "", 700, 600);
    c_time_e_ranges->SetLogy();

    std::vector<TH1D *> h_time_projs;
    double ymax_time = 0.0;

    for (int irange = 0; irange < n_e_ranges; irange++)
    {
        h_tower_e_vs_time->GetYaxis()->SetRangeUser(e_mins[irange], e_maxs[irange]);
        TH1D *h_proj = (TH1D *)h_tower_e_vs_time->ProjectionX(Form("h_time_e_proj_%d", irange));
        h_proj->Rebin(4);
        double integral = h_proj->Integral();
        if (integral > 0)
            h_proj->Scale(1.0 / integral);

        h_proj->SetLineColor(col_samples[irange]);
        h_proj->SetMarkerColor(col_samples[irange]);
        h_proj->SetLineWidth(2);
        h_proj->SetMarkerStyle(20);

        ymax_time = std::max(ymax_time, h_proj->GetMaximum());
        h_time_projs.push_back(h_proj);
    }

    TH1D *h_frame_time = new TH1D("h_frame_time", "", 100, -10, 10);
    h_frame_time->SetTitle(";Tower Time [ns];Normalized Counts");
    h_frame_time->GetXaxis()->SetTitleSize(0.045);
    h_frame_time->GetYaxis()->SetTitleSize(0.045);
    h_frame_time->GetYaxis()->SetTitleOffset(1.3);
    h_frame_time->GetYaxis()->SetRangeUser(1e-4, 2.0 * ymax_time);
    h_frame_time->Draw("axis");

    for (size_t i = 0; i < h_time_projs.size(); i++)
    {
        h_time_projs[i]->Draw("same hist");
    }

    myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
    myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);

    for (int i = 0; i < n_e_ranges; i++)
    {
        myMarkerLineText(0.65, 0.88 - i * 0.05, 1, col_samples[i], 20, col_samples[i], 1,
                         e_labels[i].c_str(), 0.04, true);
    }

    c_time_e_ranges->SaveAs(Form("%s/tower_time_energy_ranges.pdf", outputdir.c_str()));

    std::cout << "Created time distribution for energy ranges" << std::endl;

    // Reset ranges
    h_tower_e_vs_time->GetYaxis()->SetRange(0, 0);

    // ========================================================================
    // Part 5: Delta-t (leading - non-leading towers)
    // ========================================================================

    TCanvas *c_delta_t = nullptr;
    TH1D *h_delta_t_all = nullptr;
    if (h2_delta_t_leading_nonleading_vs_e)
    {
        std::cout << "Creating delta-t (leading - non-leading) distribution plot..." << std::endl;

        // Focus x-range around the peak (adjust here if you want a different window)
        const double delta_t_xmin = -5.0;
        const double delta_t_xmax = 5.0;

        c_delta_t = new TCanvas("c_delta_t", "", 800, 600);
        c_delta_t->SetLeftMargin(0.14);
        c_delta_t->SetRightMargin(0.05);
        c_delta_t->SetBottomMargin(0.12);

        h_delta_t_all = (TH1D *)h2_delta_t_leading_nonleading_vs_e->ProjectionX("h_delta_t_all");
        h_delta_t_all->SetLineColor(kBlack);
        h_delta_t_all->SetLineWidth(2);
        h_delta_t_all->SetMarkerColor(kBlack);
        h_delta_t_all->SetMarkerStyle(20);

        h_delta_t_all->SetTitle(";#Deltat = t_{leading} - t_{non-leading} [ns];Counts (weighted)");
        h_delta_t_all->GetXaxis()->SetTitleSize(0.05);
        h_delta_t_all->GetXaxis()->SetLabelSize(0.045);
        h_delta_t_all->GetXaxis()->SetTitleOffset(1.1);
        h_delta_t_all->GetYaxis()->SetTitleSize(0.05);
        h_delta_t_all->GetYaxis()->SetLabelSize(0.045);
        h_delta_t_all->GetYaxis()->SetTitleOffset(1.3);
        h_delta_t_all->GetXaxis()->SetRangeUser(delta_t_xmin, delta_t_xmax);

        h_delta_t_all->Draw("hist");

        myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
        myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);
        myText(0.20, 0.78, 1, "(Direct + Frag)", 0.045);
        myText(0.20, 0.73, 1, Form("Mean = %.3f ns", h_delta_t_all->GetMean()), 0.045);

        c_delta_t->SaveAs(Form("%s/delta_t_leading_minus_nonleading.pdf", outputdir.c_str()));
        std::cout << "Created delta-t distribution plot" << std::endl;
    }
    else
    {
        std::cout << "WARNING: Missing h2_delta_t_leading_nonleading_vs_e; skipping delta-t plots" << std::endl;
    }

    // ========================================================================
    // Part 6: Overlay delta-t for non-leading tower energy ranges
    // ========================================================================

    TCanvas *c_delta_t_e_ranges = nullptr;
    std::vector<TH1D *> h_delta_t_projs;
    TH1D *h_frame_delta_t = nullptr;
    if (h2_delta_t_leading_nonleading_vs_e)
    {
        std::cout << "Creating delta-t overlays for non-leading tower energy ranges..." << std::endl;

        const int n_sub_e_ranges = 4;
        float sube_mins[n_sub_e_ranges] = {0.1, 0.5, 1.0, 2.0};
        float sube_maxs[n_sub_e_ranges] = {0.5, 1.0, 2.0, 5.0};
        std::string sube_labels[n_sub_e_ranges] = {
            "0.1 < E_{non-leading} < 0.5 GeV",
            "0.5 < E_{non-leading} < 1.0 GeV",
            "1.0 < E_{non-leading} < 2.0 GeV",
            "2.0 < E_{non-leading} < 5.0 GeV"};

        // Match the window used in the inclusive delta-t plot
        const double delta_t_xmin = -5.0;
        const double delta_t_xmax = 5.0;

        c_delta_t_e_ranges = new TCanvas("c_delta_t_e_ranges", "", 800, 600);
        c_delta_t_e_ranges->SetLeftMargin(0.14);
        c_delta_t_e_ranges->SetRightMargin(0.05);
        c_delta_t_e_ranges->SetBottomMargin(0.12);
        c_delta_t_e_ranges->SetLogy();

        double ymax_delta_t = 0.0;
        for (int irange = 0; irange < n_sub_e_ranges; irange++)
        {
            h2_delta_t_leading_nonleading_vs_e->GetYaxis()->SetRangeUser(sube_mins[irange], sube_maxs[irange]);
            TH1D *h_proj = (TH1D *)h2_delta_t_leading_nonleading_vs_e->ProjectionX(Form("h_delta_t_subE_%d", irange));
            h_proj->Rebin(2);
            double integral = h_proj->Integral();
            if (integral > 0)
                h_proj->Scale(1.0 / integral);

            h_proj->GetXaxis()->SetRangeUser(delta_t_xmin, delta_t_xmax);
            h_proj->SetLineColor(col_samples[irange]);
            h_proj->SetMarkerColor(col_samples[irange]);
            h_proj->SetLineWidth(2);
            h_proj->SetMarkerStyle(20);

            ymax_delta_t = std::max(ymax_delta_t, h_proj->GetMaximum());
            h_delta_t_projs.push_back(h_proj);
        }

        // Reset Y range
        h2_delta_t_leading_nonleading_vs_e->GetYaxis()->SetRange(0, 0);

        h_frame_delta_t = new TH1D("h_frame_delta_t", "", 100, delta_t_xmin, delta_t_xmax);
        h_frame_delta_t->SetTitle(";#Deltat = t_{leading} - t_{non-leading} [ns];Normalized Counts");
        h_frame_delta_t->GetXaxis()->SetTitleSize(0.05);
        h_frame_delta_t->GetXaxis()->SetLabelSize(0.045);
        h_frame_delta_t->GetXaxis()->SetTitleOffset(1.1);
        h_frame_delta_t->GetYaxis()->SetTitleSize(0.05);
        h_frame_delta_t->GetYaxis()->SetLabelSize(0.045);
        h_frame_delta_t->GetYaxis()->SetTitleOffset(1.3);
        h_frame_delta_t->GetYaxis()->SetRangeUser(1e-4, 2.0 * ymax_delta_t);
        h_frame_delta_t->Draw("axis");

        for (auto h : h_delta_t_projs)
            h->Draw("same hist");

        myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
        myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);

        for (int i = 0; i < n_sub_e_ranges; i++)
        {
            myMarkerLineText(0.65, 0.88 - i * 0.05, 1, col_samples[i], 20, col_samples[i], 1,
                             sube_labels[i].c_str(), 0.04, true);
        }

        c_delta_t_e_ranges->SaveAs(Form("%s/delta_t_nonleading_energy_ranges.pdf", outputdir.c_str()));
        std::cout << "Created delta-t overlays for non-leading energy ranges" << std::endl;
    }

    // ========================================================================
    // Part 7: Mean delta-t vs non-leading tower energy
    // ========================================================================

    TCanvas *c_mean_dt_vs_sube = nullptr;
    TProfile *prof_mean_dt_vs_sube = nullptr;
    if (h2_delta_t_leading_nonleading_vs_e)
    {
        std::cout << "Creating mean delta-t vs non-leading tower energy profile..." << std::endl;

        // Mean delta-t (X) as a function of non-leading energy (Y) -> ProfileY
        prof_mean_dt_vs_sube = h2_delta_t_leading_nonleading_vs_e->ProfileY("prof_mean_dt_vs_sube", 1, -1, "");

        c_mean_dt_vs_sube = new TCanvas("c_mean_dt_vs_sube", "", 800, 600);
        c_mean_dt_vs_sube->SetLeftMargin(0.14);
        c_mean_dt_vs_sube->SetRightMargin(0.05);
        c_mean_dt_vs_sube->SetBottomMargin(0.12);

        prof_mean_dt_vs_sube->SetLineColor(kBlue);
        prof_mean_dt_vs_sube->SetLineWidth(3);
        prof_mean_dt_vs_sube->SetMarkerColor(kBlue);
        prof_mean_dt_vs_sube->SetMarkerStyle(20);
        prof_mean_dt_vs_sube->SetMarkerSize(1.0);

        prof_mean_dt_vs_sube->SetTitle(";Non-Leading Tower Energy [GeV];Mean #Deltat = t_{leading} - t_{non-leading} [ns]");
        prof_mean_dt_vs_sube->GetXaxis()->SetTitleSize(0.05);
        prof_mean_dt_vs_sube->GetXaxis()->SetLabelSize(0.045);
        prof_mean_dt_vs_sube->GetXaxis()->SetTitleOffset(1.1);
        prof_mean_dt_vs_sube->GetYaxis()->SetTitleSize(0.05);
        prof_mean_dt_vs_sube->GetYaxis()->SetLabelSize(0.045);
        prof_mean_dt_vs_sube->GetYaxis()->SetTitleOffset(1.3);

        // Plot window (adjust if needed)
        prof_mean_dt_vs_sube->GetXaxis()->SetRangeUser(0.0, 20.0);
        prof_mean_dt_vs_sube->GetYaxis()->SetRangeUser(-2.0, 2.0);

        prof_mean_dt_vs_sube->Draw();

        myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
        myText(0.20, 0.83, 1, "Truth-Matched Photons", 0.045);
        myText(0.20, 0.78, 1, "(Direct + Frag)", 0.045);

        c_mean_dt_vs_sube->SaveAs(Form("%s/mean_delta_t_vs_nonleading_energy.pdf", outputdir.c_str()));
        std::cout << "Created mean delta-t vs non-leading energy plot" << std::endl;
    }

    // ========================================================================
    // Summary
    // ========================================================================

    std::cout << "\n========================================" << std::endl;
    std::cout << "Summary:" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Created plots:" << std::endl;
    std::cout << "  - Tower energy vs time profile (2D)" << std::endl;
    std::cout << "  - Tower energy vs time profile (inclusive eta from 3D)" << std::endl;
    std::cout << "  - Tower energy vs time profile in eta slices (3D profile)" << std::endl;
    std::cout << "  - Cluster time mean vs cluster energy (3 methods)" << std::endl;
    std::cout << "  - Cluster time RMS vs cluster energy (3 methods)" << std::endl;
    std::cout << "  - Tower ADC vs time profile (2D)" << std::endl;
    std::cout << "  - Tower ADC vs time profile (inclusive eta from 3D)" << std::endl;
    std::cout << "  - Tower ADC vs time profile in eta slices (3D profile)" << std::endl;
    std::cout << "  - Tower energy vs time 2D with profile overlay" << std::endl;
    std::cout << "  - Tower ADC vs time 2D with profile overlay" << std::endl;
    std::cout << "  - Tower energy/ADC vs time 2D in eta slices" << std::endl;
    std::cout << "  - Time distributions for different energy ranges" << std::endl;
    std::cout << "  - Delta-t (leading - non-leading) distribution" << std::endl;
    std::cout << "  - Delta-t overlays in non-leading energy ranges" << std::endl;
    std::cout << "  - Mean delta-t vs non-leading energy" << std::endl;
    std::cout << "========================================" << std::endl;

    // Cleanup
    delete h_frame_time;
    for (auto h : h_time_projs)
        delete h;
    delete h_frame_delta_t;
    for (auto h : h_delta_t_projs)
        delete h;
    delete c_delta_t_e_ranges;
    delete c_delta_t;
    delete c_mean_dt_vs_sube;
    delete c_time_e_ranges;
    delete c_adc_vs_time_2d;
    delete c_e_vs_time_2d;
    delete c_time_vs_adc;
    delete c_time_vs_energy;
    delete prof_time_vs_adc;
    delete prof_time_vs_e;
    delete prof_leading_time_vs_adc;
    delete prof_leading_time_vs_e;
    delete prof_subleading_time_vs_adc;
    delete prof_subleading_time_vs_e;
    delete prof_nonleading_time_vs_adc;
    delete prof_nonleading_time_vs_e;
    delete prof_mean_dt_vs_sube;
    delete h_delta_t_all;

    fin->Close();

    std::cout << "\nAll plots saved to: " << outputdir << std::endl;
}
