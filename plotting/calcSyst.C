#include "plotcommon.h"

void calcSyst()
{

    init_plot();

    std::string file_path = "rootFiles/";
    std::string savePath = "figures/";

    std::vector<std::string> sys_names =
        {
            "iso",
            "nt",
            "tight",
            "nor"
            };

    std::vector<std::string> sys_names_leg =
        {
            "iso variation",
            "non-tight variation",
            "tight variation",
            "unfolding variation"
            };

    std::vector<int> colors = {kPink + 8, kSpring - 7, kAzure - 3,  kOrange + 10, kViolet + 3, kYellow + 4};

    std::vector<TH1F *> h_dev_low_list;
    std::vector<TH1F *> h_dev_high_list;

    std::vector<TH1F *> h_dev_rel_low_list;
    std::vector<TH1F *> h_dev_rel_high_list;

    for (int isys = 0; isys < (int)sys_names.size(); isys++)
    {
        std::string filename = file_path + "syst_" + sys_names[isys] + ".root";
        std::cout << "reading " << filename << std::endl;
        TFile *f = new TFile(filename.c_str(), "READ");

        TH1F *h_dev_low = (TH1F *)f->Get("h_dev_low");
        TH1F *h_dev_high = (TH1F *)f->Get("h_dev_high");

        TH1F *h_dev_rel_low = (TH1F *)f->Get("h_dev_rel_low");
        TH1F *h_dev_rel_high = (TH1F *)f->Get("h_dev_rel_high");

        h_dev_low_list.push_back(h_dev_low);
        h_dev_high_list.push_back(h_dev_high);

        h_dev_rel_low_list.push_back(h_dev_rel_low);
        h_dev_rel_high_list.push_back(h_dev_rel_high);
    }

    TH1F *h_sum_low = (TH1F *)h_dev_low_list.at(0)->Clone("h_sum_low");
    h_sum_low->Reset();
    TH1F *h_sum_high = (TH1F *)h_dev_high_list.at(0)->Clone("h_sum_high");
    h_sum_high->Reset();

    TH1F *h_sum_rel_low = (TH1F *)h_dev_rel_low_list.at(0)->Clone("h_sum_rel_low");
    h_sum_rel_low->Reset();
    TH1F *h_sum_rel_high = (TH1F *)h_dev_rel_high_list.at(0)->Clone("h_sum_rel_high");
    h_sum_rel_high->Reset();

    // calculate sum syst
    for (int i = 1; i < h_sum_low->GetNbinsX() + 1; i++)
    {
        float error_low2 = 0;
        float error_high2 = 0;

        float error_low_rel2 = 0;
        float error_high_rel2 = 0;

        for (int isys = 0; isys < (int)sys_names.size(); isys++)
        {
            float error_low_sys = h_dev_low_list.at(isys)->GetBinContent(i);
            float error_high_sys = h_dev_high_list.at(isys)->GetBinContent(i);

            error_low2 += error_low_sys * error_low_sys;
            error_high2 += error_high_sys * error_high_sys;

            float error_low_sys_rel = h_dev_rel_low_list.at(isys)->GetBinContent(i);
            float error_high_sys_rel = h_dev_rel_high_list.at(isys)->GetBinContent(i);

            error_low_rel2 += error_low_sys_rel * error_low_sys_rel;
            error_high_rel2 += error_high_sys_rel * error_high_sys_rel;
        }

        error_low2 = sqrt(error_low2);
        error_high2 = sqrt(error_high2);

        error_low_rel2 = sqrt(error_low_rel2);
        error_high_rel2 = sqrt(error_high_rel2);

        h_sum_low->SetBinContent(i, error_low2);
        h_sum_high->SetBinContent(i, error_high2);

        h_sum_rel_low->SetBinContent(i, error_low_rel2);
        h_sum_rel_high->SetBinContent(i, error_high_rel2);
    }
    // place holder for plotting
    // low
    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);

    frame_et_truth->GetYaxis()->SetRangeUser(0, 0.7);
    frame_et_truth->GetXaxis()->SetRangeUser(10, 30);
    frame_et_truth->SetYTitle("Relative difference");
    frame_et_truth->Draw("axis");

    h_sum_rel_low->SetMarkerStyle(20);
    h_sum_rel_low->SetMarkerColor(kBlack);
    h_sum_rel_low->SetLineColor(kBlack);
    h_sum_rel_low->Draw("same ][ HIST");

    myMarkerLineText(0.25, 0.80, 1, kBlack, 20, kBlack, 0, "sum", 0.05, true);

    for (int isys = 0; isys < (int)sys_names.size(); isys++)
    {
        h_dev_rel_low_list.at(isys)->SetMarkerStyle(20);
        h_dev_rel_low_list.at(isys)->SetMarkerColor(colors.at(isys));
        h_dev_rel_low_list.at(isys)->SetLineColor(colors.at(isys));
        h_dev_rel_low_list.at(isys)->Draw("same ][ HIST");

        myMarkerLineText(0.25, 0.75 - 0.05 * isys, 1, colors.at(isys), 20, colors.at(isys), 0, sys_names_leg.at(isys).c_str(), 0.05, true);
    }
    myText(0.25, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.25, 0.85, 1, strleg2.c_str(), 0.05);

    c1->SaveAs(Form("%ssyst_sum_rel.pdf", savePath.c_str()));

    // save to file
    TFile *fout = new TFile(Form("%ssyst_sum.root",file_path.c_str()), "RECREATE");
    h_sum_low->Write();
    h_sum_high->Write();
    h_sum_rel_low->Write();
    h_sum_rel_high->Write();
}