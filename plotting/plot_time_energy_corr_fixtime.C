#include "plotcommon.h"
#include <vector>
#include <algorithm>

const int col_ratio[] = {kAzure + 2, kPink + 5, kSpring - 6, kRed - 4, kBlue - 4};

void plot_time_energy_corr_fixtime(const char* infile = "/sphenix/user/shuhangli/ppg12/efficiencytool/time_energy_corr.root")
{
    init_plot();

    std::string base_dir = "/sphenix/user/shuhangli/ppg12/efficiencytool/";

    std::string m8_name = "time_energy_corr_m8.root";
    std::string p8_name = "time_energy_corr_p8.root";
    std::string no_shift_name = "time_energy_corr_0.root";

    std::string f_m8_path = base_dir + m8_name.c_str();
    std::string f_p8_path = base_dir + p8_name.c_str();
    std::string f_no_shift_path = base_dir + no_shift_name.c_str();


    TFile *f_m8 = new TFile(f_m8_path.c_str());    
    TFile *f_p8 = new TFile(f_p8_path.c_str());
    TFile *f_no_shift = new TFile(f_no_shift_path.c_str());

    if (!f_m8 || f_m8->IsZombie() || !f_p8 || f_p8->IsZombie() || !f_no_shift || f_no_shift->IsZombie())
    {
        std::cerr << "Failed to open files: " << m8_name << ", " << p8_name << ", " << no_shift_name << std::endl;
        return;
    }

    TH3D *h_time_ratio_truth_m8 = dynamic_cast<TH3D*>(f_m8->Get("h_time_ratio_truth"));
    TH3D *h_time_ratio_truth_p8 = dynamic_cast<TH3D*>(f_p8->Get("h_time_ratio_truth"));
    TH3D *h_time_ratio_truth_no_shift = dynamic_cast<TH3D*>(f_no_shift->Get("h_time_ratio_truth"));

    if (!h_time_ratio_truth_m8 || h_time_ratio_truth_m8->IsZombie() || !h_time_ratio_truth_p8 || h_time_ratio_truth_p8->IsZombie() || !h_time_ratio_truth_no_shift || h_time_ratio_truth_no_shift->IsZombie())
    {
        std::cerr << "Failed to get h_time_ratio_truth from files: " << m8_name << ", " << p8_name << ", " << no_shift_name << std::endl;
        return;
    }

    float time_min = -3;
    float time_max = 3;

    h_time_ratio_truth_m8->GetXaxis()->SetRangeUser(time_min, time_max);
    h_time_ratio_truth_p8->GetXaxis()->SetRangeUser(time_min, time_max);
    h_time_ratio_truth_no_shift->GetXaxis()->SetRangeUser(time_min, time_max);

    TH1D *h_ratio_m8 = h_time_ratio_truth_m8->ProjectionY("h_ratio_m8");
    TH1D *h_ratio_p8 = h_time_ratio_truth_p8->ProjectionY("h_ratio_p8");
    TH1D *h_ratio_no_shift = h_time_ratio_truth_no_shift->ProjectionY("h_ratio_no_shift");

    h_ratio_m8->SetLineColor(col_ratio[0]);
    h_ratio_p8->SetLineColor(col_ratio[1]);
    h_ratio_no_shift->SetLineColor(col_ratio[2]);

    h_ratio_m8->SetMarkerColor(col_ratio[0]);
    h_ratio_p8->SetMarkerColor(col_ratio[1]);
    h_ratio_no_shift->SetMarkerColor(col_ratio[2]);

    h_ratio_m8->SetMarkerStyle(20);
    h_ratio_p8->SetMarkerStyle(20);
    h_ratio_no_shift->SetMarkerStyle(20);

    h_ratio_m8->SetLineWidth(2);
    h_ratio_p8->SetLineWidth(2);
    h_ratio_no_shift->SetLineWidth(2);

    float scale_min = 0.8;
    float scale_max = 1.2;
    //normalize to unit area
    h_ratio_m8->Scale(1.0 / h_ratio_m8->Integral(h_ratio_m8->GetXaxis()->FindBin(scale_min), h_ratio_m8->GetXaxis()->FindBin(scale_max)));
    h_ratio_p8->Scale(1.0 / h_ratio_p8->Integral(h_ratio_p8->GetXaxis()->FindBin(scale_min), h_ratio_p8->GetXaxis()->FindBin(scale_max)));
    h_ratio_no_shift->Scale(1.0 / h_ratio_no_shift->Integral(h_ratio_no_shift->GetXaxis()->FindBin(scale_min), h_ratio_no_shift->GetXaxis()->FindBin(scale_max)));

    double ymax = std::max(h_ratio_m8->GetMaximum(), std::max(h_ratio_p8->GetMaximum(), h_ratio_no_shift->GetMaximum()));

    TCanvas *c = new TCanvas("c", "c", 800, 600);

    TH1F *frame = new TH1F("frame_time_ratio", "", 20, 0.8, 1.2);
    frame->GetXaxis()->SetTitle("Reco E_{T} / Truth p_{T}");
    frame->GetYaxis()->SetTitle("Normalized counts");
    frame->GetYaxis()->SetRangeUser(0.0, ymax * 1.7);
    frame->Draw("axis");

    h_ratio_m8->Draw("hist same");
    h_ratio_p8->Draw("hist same");
    h_ratio_no_shift->Draw("hist same");

    float xpos(0.2), xpos2(0.915), ypos(0.885), ypos2(0.65), dy(0.054), dy1(0.055), fontsize(0.046), fontsize1(0.048);
    myText(xpos2, ypos - 0 * dy, 1, strMC.c_str(), fontsize, 1);

    TLegend *leg = new TLegend(0.45, ypos2, xpos2, ypos2 + 3 * dy1);
    legStyle(leg, 0.20, 0.06);
    leg->AddEntry(h_ratio_m8, "+8 ns", "l");
    leg->AddEntry(h_ratio_p8, "-8 ns", "l");
    leg->AddEntry(h_ratio_no_shift, "No shift", "l");
    leg->Draw("same");

    c->Update();

    c->SaveAs("figures/time_energy_corr_fixtime.pdf");

}
