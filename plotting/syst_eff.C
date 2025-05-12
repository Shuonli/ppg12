#include "plotcommon.h"

void syst_eff()
{

    init_plot();

    string varStr = "eff";
    string savePath = "figures";
    std::string file_path = "rootFiles/";

    std::vector<std::string> sys_names =
        {
            "tight",
            "fudge"
            };

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

    TH1F *h_dev_low = (TH1F *)h_dev_low_list.at(0)->Clone("h_dev_low");
    TH1F *h_dev_high = (TH1F *)h_dev_high_list.at(0)->Clone("h_dev_high");

    TH1F *h_dev_rel_low = (TH1F *)h_dev_rel_low_list.at(0)->Clone("h_dev_rel_low");
    TH1F *h_dev_rel_high = (TH1F *)h_dev_rel_high_list.at(0)->Clone("h_dev_rel_high");

    for (int i = 1; i < h_dev_low->GetNbinsX() + 1; i++)
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

        h_dev_low->SetBinContent(i, error_low2);
        h_dev_high->SetBinContent(i, error_high2);

        h_dev_rel_low->SetBinContent(i, error_low_rel2);
        h_dev_rel_high->SetBinContent(i, error_high_rel2);
    }

    TFile *systOut = new TFile(Form("rootFiles/syst_%s.root", varStr.c_str()), "RECREATE");
    h_dev_low->Write();
    h_dev_high->Write();
    h_dev_rel_low->Write();
    h_dev_rel_high->Write();
    systOut->Close();
}