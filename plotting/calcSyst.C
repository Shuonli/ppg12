#include "plotcommon.h"

void calcSyst()
{

    init_plot();

    std::string file_path = "rootFiles/";
    std::string savePath = "figures/";

    std::vector<std::string> sys_names =
        {
            "purity",
            "escale",
            "eres",
            "tight",
            "nor",
            "mbd"
            };

    std::vector<std::string> sys_names_leg =
        {
            "Purity",
            "Energy scale",
            "Energy resolution",
            "Efficiency",
            "Unfolding",
            "MBD efficiency"
        };

    std::vector<int> colors = {kPink+5,kGreen-2,kAzure+7,kRed-4,kBlue-3,kYellow+2,kPink-5,kGreen+3,kBlue-3,kBlack};

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
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

    frame_et_truth->GetYaxis()->SetRangeUser(-0.58, 0.5);
    frame_et_truth->GetXaxis()->SetRangeUser(10.1, 25.9);
    frame_et_truth->SetXTitle("#it{E}_{T}^{#gamma} [GeV]");
    frame_et_truth->SetYTitle("Relative difference");
    frame_et_truth->Draw("axis");

    h_sum_rel_low->Scale(-1);

    h_sum_rel_low->SetMarkerStyle(20);
    h_sum_rel_low->SetLineWidth(4);
    h_sum_rel_low->SetMarkerColor(kBlack);
    h_sum_rel_low->SetLineColor(kBlack);
    h_sum_rel_low->Draw("same ][ HIST");

    h_sum_rel_high->SetMarkerStyle(20);
    h_sum_rel_high->SetLineWidth(4);
    h_sum_rel_high->SetMarkerColor(kBlack);
    h_sum_rel_high->SetLineColor(kBlack);
    h_sum_rel_high->Draw("same ][ HIST");

    myMarkerLineText(0.25, 0.42, 0, kBlack, 20, kBlack, 0, "Total", 0.05, true);

    int switchover = 3;

    for (int isys = 0; isys < (int)sys_names.size(); isys++)
    {
        h_dev_rel_low_list.at(isys)->SetMarkerStyle(20);
        h_dev_rel_low_list.at(isys)->Scale(-1);
        h_dev_rel_low_list.at(isys)->SetMarkerColor(colors.at(isys));
        h_dev_rel_low_list.at(isys)->SetLineColor(colors.at(isys));
        h_dev_rel_low_list.at(isys)->Draw("same ][ HIST");

        h_dev_rel_high_list.at(isys)->SetMarkerStyle(20);
        h_dev_rel_high_list.at(isys)->SetMarkerColor(colors.at(isys));
        h_dev_rel_high_list.at(isys)->SetLineColor(colors.at(isys));
        h_dev_rel_high_list.at(isys)->Draw("same ][ HIST");

        float xshift = 0;
        float yshift = 0;

        if(isys >= switchover)
        {
            xshift = 0.33;
            yshift = 0.05 *switchover;
        }

        myMarkerLineText(0.25 + xshift, 0.37 - 0.05 * isys + yshift, 0, colors.at(isys), 20, colors.at(isys), 0, sys_names_leg.at(isys).c_str(), 0.05, true);
    }
    myText(0.20, 0.88, 1, strleg1.c_str(), 0.05);
    myText(0.20, 0.83, 1, strleg2.c_str(), 0.05);

    c1->SaveAs(Form("%ssyst_sum_rel.pdf", savePath.c_str()));

    // add lumi error
    float lumierror_up = 4.3/26.1;
    float lumierror_down = 1.1/26.1;
    
    for (int i = 1; i < h_sum_low->GetNbinsX() + 1; i++)
    {
        float error_low2 = h_sum_low->GetBinContent(i);
        float error_high2 = h_sum_high->GetBinContent(i);

        float error_low_rel2 = h_sum_rel_low->GetBinContent(i);
        float error_high_rel2 = h_sum_rel_high->GetBinContent(i);

        error_low2 = sqrt(error_low2 * error_low2 + lumierror_down * lumierror_down);
        error_high2 = sqrt(error_high2 * error_high2 + lumierror_up * lumierror_up);

        error_low_rel2 = sqrt(error_low_rel2 * error_low_rel2 + lumierror_down * lumierror_down);
        error_high_rel2 = sqrt(error_high_rel2 * error_high_rel2 + lumierror_up * lumierror_up);

        h_sum_low->SetBinContent(i, error_low2);
        h_sum_high->SetBinContent(i, error_high2);

        h_sum_rel_low->SetBinContent(i, error_low_rel2);
        h_sum_rel_high->SetBinContent(i, error_high_rel2);
    }
    

    // save to file
    TFile *fout = new TFile(Form("%ssyst_sum.root",file_path.c_str()), "RECREATE");
    h_sum_low->Write();
    h_sum_high->Write();
    h_sum_rel_low->Write();
    h_sum_rel_high->Write();
}
