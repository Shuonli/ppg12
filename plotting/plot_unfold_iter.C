#include "plotcommon.h"

void plot_unfold_iter()
{
    init_plot();

    string savePath = "figures/";

    TFile *fin = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_nom.root");
    static const int ntotal_iterations = 10;

    string leg1 = "quadratic sum of rellative statistical uncertainty";
    string leg2 = "quadratic sum of relative deviation";

    std::string h_result_name_base = "h_unfold_sub_leak_";

    TH1F *h_stat = new TH1F("h_stat", "Statistical Uncertainty", ntotal_iterations, 0.5, ntotal_iterations + 0.5);

    TH1F *h_iter_delta = new TH1F("h_iter_delta", "Iteration Delta", ntotal_iterations, 0.5, ntotal_iterations + 0.5);

    for (int i = 1; i < ntotal_iterations; i++)
    {
        std::string h_prev_name = h_result_name_base + std::to_string(i);
        std::string h_this_name = h_result_name_base + std::to_string(i + 1);

        TH1F *h_prev = (TH1F *)fin->Get(h_prev_name.c_str());
        TH1F *h_this = (TH1F *)fin->Get(h_this_name.c_str());

        TH1F *h_dev_rel = (TH1F *)calcDelta(h_this, h_prev, "h_dev_rel_" + std::to_string(i)).second;

        // loop over bin find the quad sum
        double rel_quad_sum = 0;
        double rel_stat_sum = 0;

        //for (int ibin = 1; ibin <= h_dev_rel->GetNbinsX(); ibin++)
        for (int ibin = 2; ibin <= 6; ibin++)
        {
            double rel = h_dev_rel->GetBinContent(ibin);
            double stat = h_this->GetBinError(ibin) / h_this->GetBinContent(ibin);
            //std::cout << "ibin: " << ibin << " rel: " << rel << " stat: " << stat << std::endl;
            if (rel == rel)
            {
                rel_quad_sum += rel * rel;
            }
            if (stat == stat)
            {
                rel_stat_sum += stat * stat;
            }
        }
        rel_quad_sum = sqrt(rel_quad_sum);
        rel_stat_sum = sqrt(rel_stat_sum);

        std::cout << "rel_quad_sum: " << rel_quad_sum << " rel_stat_sum: " << rel_stat_sum << std::endl;

        h_stat->SetBinContent(i + 1, rel_stat_sum);
        h_stat->SetBinError(i + 1, 0);
        h_iter_delta->SetBinContent(i + 1, rel_quad_sum);
        h_iter_delta->SetBinError(i + 1, 0);
    }

    TCanvas *c1 = new TCanvas("can", "", 600, 600);

    frame_et_truth->SetYTitle("#Delta(what name should I use?)");
    frame_et_truth->GetYaxis()->SetRangeUser(0, 0.5);
    frame_et_truth->GetXaxis()->SetRangeUser(0, 10);
    frame_et_truth->Draw("axis");

    h_stat->SetMarkerStyle(20);
    h_stat->SetMarkerColor(kBlack);
    h_stat->SetLineColor(kBlack);
    h_stat->Draw("same p");

    h_iter_delta->SetMarkerStyle(20);
    h_iter_delta->SetMarkerColor(kBlue);
    h_iter_delta->SetLineColor(kBlue);
    h_iter_delta->Draw("same p");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.5, 0.80, 1, strleg3.c_str(), 0.04);

    myMarkerLineText(0.25, 0.25, 1, kBlack, 20, kBlack, 1, leg1.c_str(), 0.04, true);
    myMarkerLineText(0.25, 0.20, 1, kBlue, 20, kBlue, 1, leg2.c_str(), 0.04, true);

    c1->SaveAs(Form("%s/unfold_iter.pdf", savePath.c_str()));
}
