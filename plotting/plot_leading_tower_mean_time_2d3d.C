#include "plotcommon.h"

#include <iostream>
#include <string>

#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TSystem.h"

namespace
{
    double sumw2_2d_including_overflow(const TH2D *h)
    {
        if (!h)
            return 0.0;
        const int nbx = h->GetNbinsX();
        const int nby = h->GetNbinsY();
        double s = 0.0;
        for (int ix = 0; ix <= nbx + 1; ++ix)
        {
            for (int iy = 0; iy <= nby + 1; ++iy)
            {
                s += h->GetBinContent(ix, iy);
            }
        }
        return s;
    }

    void compare_h2(const TH2D *a, const TH2D *b, const std::string &label_a, const std::string &label_b)
    {
        if (!a || !b)
            return;

        std::cout << label_a << " axes: X[" << a->GetXaxis()->GetXmin() << "," << a->GetXaxis()->GetXmax()
                  << "] nbins=" << a->GetNbinsX()
                  << " ; Y[" << a->GetYaxis()->GetXmin() << "," << a->GetYaxis()->GetXmax()
                  << "] nbins=" << a->GetNbinsY() << std::endl;
        std::cout << label_b << " axes: X[" << b->GetXaxis()->GetXmin() << "," << b->GetXaxis()->GetXmax()
                  << "] nbins=" << b->GetNbinsX()
                  << " ; Y[" << b->GetYaxis()->GetXmin() << "," << b->GetYaxis()->GetXmax()
                  << "] nbins=" << b->GetNbinsY() << std::endl;

        const int nbx = a->GetNbinsX();
        const int nby = a->GetNbinsY();
        if (nbx != b->GetNbinsX() || nby != b->GetNbinsY())
        {
            std::cout << "WARNING: TH2 binning mismatch: " << label_a << " is ("
                      << nbx << "," << nby << ") while " << label_b << " is ("
                      << b->GetNbinsX() << "," << b->GetNbinsY() << ")" << std::endl;
        }

        const int nbx_cmp = std::min(nbx, b->GetNbinsX());
        const int nby_cmp = std::min(nby, b->GetNbinsY());

        double max_abs = 0.0;
        double max_rel = 0.0;
        int max_abs_ix = -1, max_abs_iy = -1;
        int max_rel_ix = -1, max_rel_iy = -1;

        // include under/overflow bins (0..nb+1)
        for (int ix = 0; ix <= nbx_cmp + 1; ++ix)
        {
            for (int iy = 0; iy <= nby_cmp + 1; ++iy)
            {
                const double va = a->GetBinContent(ix, iy);
                const double vb = b->GetBinContent(ix, iy);
                const double d = std::abs(va - vb);
                if (d > max_abs)
                {
                    max_abs = d;
                    max_abs_ix = ix;
                    max_abs_iy = iy;
                }
                const double denom = std::max(1e-12, std::abs(va));
                const double r = d / denom;
                if (r > max_rel)
                {
                    max_rel = r;
                    max_rel_ix = ix;
                    max_rel_iy = iy;
                }
            }
        }

        std::cout << "TH2 compare (" << label_a << " vs " << label_b << "):\n"
                  << "  sumw(a) = " << sumw2_2d_including_overflow(a) << "\n"
                  << "  sumw(b) = " << sumw2_2d_including_overflow(b) << "\n"
                  << "  max |a-b| = " << max_abs << " at bin(ix,iy)=(" << max_abs_ix << "," << max_abs_iy << ")\n"
                  << "  max rel |a-b|/|a| = " << max_rel << " at bin(ix,iy)=(" << max_rel_ix << "," << max_rel_iy << ")\n";
    }

    // Compare by coordinates rather than assuming identical binning.
    // If swap_xy=true, compare a(x,y) with b(y,x) using axis FindBin on b.
    void compare_h2_by_coordinates(const TH2D *a, const TH2D *b, const std::string &label, bool swap_xy)
    {
        if (!a || !b)
            return;

        double max_abs = 0.0;
        double max_rel = 0.0;
        double max_x = 0.0, max_y = 0.0;

        // Compare only in-bin (exclude under/overflow for stability)
        for (int iax = 1; iax <= a->GetNbinsX(); ++iax)
        {
            const double x = a->GetXaxis()->GetBinCenter(iax);
            const int ibx = swap_xy ? b->GetXaxis()->FindBin(0.0) : b->GetXaxis()->FindBin(x);
            for (int iay = 1; iay <= a->GetNbinsY(); ++iay)
            {
                const double y = a->GetYaxis()->GetBinCenter(iay);
                const int iby = swap_xy ? b->GetYaxis()->FindBin(x) : b->GetYaxis()->FindBin(y);
                const int ibx2 = swap_xy ? b->GetXaxis()->FindBin(y) : ibx;
                const double va = a->GetBinContent(iax, iay);
                const double vb = b->GetBinContent(ibx2, iby);
                const double d = std::abs(va - vb);
                if (d > max_abs)
                {
                    max_abs = d;
                    max_x = x;
                    max_y = y;
                }
                const double denom = std::max(1e-12, std::abs(va));
                const double r = d / denom;
                if (r > max_rel)
                    max_rel = r;
            }
        }

        std::cout << "TH2 coord-compare (" << label << (swap_xy ? ", swapped XY" : ", same XY")
                  << "): max |diff|=" << max_abs << " at (x,y)=(" << max_x << "," << max_y << ")"
                  << " ; max rel=" << max_rel << std::endl;
    }

    void compare_profiles(const TProfile *a, const TProfile *b, const std::string &label_a, const std::string &label_b)
    {
        if (!a || !b)
            return;
        const int nb = std::min(a->GetNbinsX(), b->GetNbinsX());
        double max_abs = 0.0;
        int max_ib = -1;
        for (int ib = 1; ib <= nb; ++ib)
        {
            const double ea = a->GetBinEntries(ib);
            const double eb = b->GetBinEntries(ib);
            if (ea <= 0.0 && eb <= 0.0)
                continue;
            const double d = std::abs(a->GetBinContent(ib) - b->GetBinContent(ib));
            if (d > max_abs)
            {
                max_abs = d;
                max_ib = ib;
            }
        }
        std::cout << "TProfile compare (" << label_a << " vs " << label_b << "): max |diff| = " << max_abs;
        if (max_ib > 0)
            std::cout << " ns at x=" << a->GetXaxis()->GetBinCenter(max_ib);
        std::cout << std::endl;
    }
} // namespace

// Plot ONLY: mean tower time vs leading-tower energy
// - From 2D histogram: h_leading_tower_e_vs_time (time on X, energy on Y)
// - From 3D histogram: h3_leading_tower_e_vs_time_vs_eta (time on X, energy on Y, eta on Z)
//
// Produces two separate PDFs:
// - leading_tower_mean_time_vs_energy_from2d.pdf
// - leading_tower_mean_time_vs_energy_from3d_inclusive_eta.pdf
void plot_leading_tower_mean_time_2d3d(const std::string &inputfile = "../efficiencytool/results/truth_photon_tower_analysis.root",
                                       const std::string &outputdir = "figures/")
{
    init_plot();
    gSystem->mkdir(outputdir.c_str(), true);

    std::cout << "==================================================" << std::endl;
    std::cout << "  Leading Tower: <time> vs Energy (2D vs 3D)" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "Input file: " << inputfile << std::endl;
    std::cout << "Output directory: " << outputdir << std::endl;

    TFile *fin = TFile::Open(inputfile.c_str());
    if (!fin || fin->IsZombie())
    {
        std::cerr << "ERROR: Cannot open input file: " << inputfile << std::endl;
        return;
    }

    TH2D *h2_lead = (TH2D *)fin->Get("h_leading_tower_e_vs_time");
    TH3D *h3_lead = (TH3D *)fin->Get("h3_leading_tower_e_vs_time_vs_eta");

    if (!h2_lead)
        std::cout << "WARNING: Missing h_leading_tower_e_vs_time" << std::endl;
    if (!h3_lead)
        std::cout << "WARNING: Missing h3_leading_tower_e_vs_time_vs_eta" << std::endl;
    if (!h2_lead && !h3_lead)
    {
        std::cerr << "ERROR: No input histograms found; aborting." << std::endl;
        fin->Close();
        return;
    }

    // ----------------------------------------------------------------------
    // 1) From 2D
    // ----------------------------------------------------------------------
    TProfile *p2_keep = nullptr;
    if (h2_lead)
    {
        std::cout << "2D entries: " << h2_lead->GetEntries() << std::endl;
        std::cout << "2D sumw: " << sumw2_2d_including_overflow(h2_lead) << std::endl;

        // mean time (X) vs energy (Y) => ProfileY
        TProfile *p2 = h2_lead->ProfileY("p_leading_time_vs_e_from2d");
        p2->SetDirectory(nullptr);
        p2_keep = p2; // keep for later comparison

        TCanvas *c2 = new TCanvas("c_leading_time_vs_e_2d", "", 800, 600);
        c2->SetLeftMargin(0.14);
        c2->SetRightMargin(0.05);
        c2->SetBottomMargin(0.12);

        p2->SetTitle(";Leading Tower Energy [GeV];Mean Leading Tower Time [ns]");
        p2->SetLineColor(kBlue);
        p2->SetLineWidth(3);
        p2->SetMarkerColor(kBlue);
        p2->SetMarkerStyle(20);
        p2->SetMarkerSize(1.0);

        // Keep ranges loose; user can zoom later
        p2->GetXaxis()->SetTitleSize(0.05);
        p2->GetXaxis()->SetLabelSize(0.045);
        p2->GetXaxis()->SetTitleOffset(1.1);
        p2->GetYaxis()->SetTitleSize(0.05);
        p2->GetYaxis()->SetLabelSize(0.045);
        p2->GetYaxis()->SetTitleOffset(1.3);
        p2->GetXaxis()->SetRangeUser(0, 20);
        p2->GetYaxis()->SetRangeUser(-5, 10);

        p2->Draw();
        myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
        myText(0.20, 0.83, 1, "Leading tower (from 2D)", 0.045);

        c2->SaveAs(Form("%s/leading_tower_mean_time_vs_energy_from2d.pdf", outputdir.c_str()));

        delete c2;
    }

    // ----------------------------------------------------------------------
    // 2) From 3D (inclusive eta)
    // ----------------------------------------------------------------------
    if (h3_lead)
    {
        std::cout << "3D entries: " << h3_lead->GetEntries() << std::endl;
        std::cout << "3D eta axis: [" << h3_lead->GetZaxis()->GetXmin() << ", " << h3_lead->GetZaxis()->GetXmax() << "]"
                  << " with nbins=" << h3_lead->GetZaxis()->GetNbins() << std::endl;

        // IMPORTANT:
        // - We want "inclusive eta", including under/overflow in Z.
        // - TH3::Project3D respects axis ranges, so we work on a clone and set the Z range explicitly.
        TH3D *h3tmp = (TH3D *)h3_lead->Clone("h3_leading_e_vs_time_vs_eta__tmp");
        h3tmp->SetDirectory(nullptr);
        const int nbz = h3tmp->GetZaxis()->GetNbins();
        h3tmp->GetZaxis()->SetRange(0, nbz + 1); // include under/overflow

        // ROOT quirk: Project3D("xy") swaps axes (X becomes original Y, Y becomes original X).
        // We want X=time (original X) and Y=energy (original Y), so we must use "yx".
        TH2D *h2_from3d = (TH2D *)h3tmp->Project3D("yx"); // (time, energy)
        if (!h2_from3d)
        {
            std::cerr << "ERROR: Project3D(\"yx\") returned null for 3D leading tower hist." << std::endl;
            delete h3tmp;
        }
        else
        {
            h2_from3d->SetName("h2_leading_tower_e_vs_time__from3d_inclusive_eta");
            h2_from3d->SetDirectory(nullptr);

            TProfile *p3 = h2_from3d->ProfileY("p_leading_time_vs_e_from3d_inclusive_eta");
            p3->SetDirectory(nullptr);

            // Diagnostics: compare raw TH2 contents and the resulting profiles
            if (h2_lead)
            {
                compare_h2(h2_lead, h2_from3d, "h_leading_tower_e_vs_time(2D)", "Project3D(xy) of h3_leading...(3D)");
                compare_h2_by_coordinates(h2_lead, h2_from3d, "h2_lead vs h2_from3d", false);
                compare_h2_by_coordinates(h2_lead, h2_from3d, "h2_lead vs h2_from3d", true);
            }
            if (p2_keep)
                compare_profiles(p2_keep, p3, "ProfileY(2D)", "ProfileY(Project3D(xy)(3D))");

            TCanvas *c3 = new TCanvas("c_leading_time_vs_e_3d", "", 800, 600);
            c3->SetLeftMargin(0.14);
            c3->SetRightMargin(0.05);
            c3->SetBottomMargin(0.12);

            p3->SetTitle(";Leading Tower Energy [GeV];Mean Leading Tower Time [ns]");
            p3->SetLineColor(kRed);
            p3->SetLineWidth(3);
            p3->SetMarkerColor(kRed);
            p3->SetMarkerStyle(21);
            p3->SetMarkerSize(1.0);

            p3->GetXaxis()->SetTitleSize(0.05);
            p3->GetXaxis()->SetLabelSize(0.045);
            p3->GetXaxis()->SetTitleOffset(1.1);
            p3->GetYaxis()->SetTitleSize(0.05);
            p3->GetYaxis()->SetLabelSize(0.045);
            p3->GetYaxis()->SetTitleOffset(1.3);
            p3->GetXaxis()->SetRangeUser(0, 20);
            p3->GetYaxis()->SetRangeUser(-5, 10);

            p3->Draw();
            myText(0.20, 0.88, 1, strleg1.c_str(), 0.045);
            myText(0.20, 0.83, 1, "Leading tower (from 3D, inclusive #eta)", 0.045);

            c3->SaveAs(Form("%s/leading_tower_mean_time_vs_energy_from3d_inclusive_eta.pdf", outputdir.c_str()));

            delete c3;
            delete p3;
            delete h2_from3d;
            delete h3tmp;
        }
    }

    delete p2_keep;
    fin->Close();
    std::cout << "Done. Plots saved to: " << outputdir << std::endl;
}


