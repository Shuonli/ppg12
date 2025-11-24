#include "plotcommon.h"
#include <algorithm>

void analyze_and_plot(TGraphErrors* graph, const string& name, const string& savePath, float minrun, float maxrun)
{
    vector<double> yvals;
    for (int i = 0; i < graph->GetN(); ++i)
    {
        double x, y;
        graph->GetPoint(i, x, y);
        yvals.push_back(y);
    }

    // Sort and compute 95th percentile
    sort(yvals.begin(), yvals.end());
    double percentile95 = yvals[int(0.95 * yvals.size())];
    double percentile97_5 = yvals[int(0.975 * yvals.size())];
    double percentile2_5 = yvals[int(0.025 * yvals.size())];
    
    double mean = std::accumulate(yvals.begin(), yvals.end(), 0.0) / yvals.size();
    double sq_sum = std::inner_product(yvals.begin(), yvals.end(), yvals.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / yvals.size() - mean * mean);

    cout << name << " mean: " << mean << ", spread (RMS/mean): " << stdev / mean
         << ", 95th percentile: " << percentile95 << endl;
    cout << name << " 97.5th percentile: " << percentile97_5 << ", 2.5th percentile: " << percentile2_5 << endl;

    // Create plot
    TCanvas *c = new TCanvas(("c_" + name).c_str(), name.c_str(), 1200, 600);
    TH1D *frame = new TH1D(("frame_" + name).c_str(), (";" + string("Run Number;") + name).c_str(), 1, minrun, maxrun);
    frame->GetYaxis()->SetRangeUser(*min_element(yvals.begin(), yvals.end()) * 0.9, *max_element(yvals.begin(), yvals.end()) * 1.1);
    frame->GetXaxis()->SetNdivisions(505);
    frame->Draw("axis");

    graph->SetMarkerStyle(20);
    graph->Draw("P same");

    // Draw 95 percentile line
    //TLine *line = new TLine(minrun, percentile95, maxrun, percentile95);
    //line->SetLineStyle(2);
    //line->SetLineColor(kBlue + 2);
    //line->SetLineWidth(2);
    //line->Draw();
//
    //TLatex *latex = new TLatex(minrun + 10, percentile95 * 1.01, "95^{th} percentile");
    //latex->SetTextColor(kBlue + 2);
    //latex->SetTextSize(0.04);
    //latex->Draw();
    //draw 97.5th and 2.5th percentile lines
    TLine *line97_5 = new TLine(minrun, percentile97_5, maxrun, percentile97_5);
    line97_5->SetLineStyle(2);
    line97_5->SetLineColor(kRed);
    line97_5->SetLineWidth(2);
    line97_5->Draw();
    TLatex *latex97_5 = new TLatex(minrun + 10, percentile97_5 * 1.01, "97.5^{th} percentile");
    latex97_5->SetTextColor(kRed);
    latex97_5->SetTextSize(0.04);
    latex97_5->Draw();
    TLine *line2_5 = new TLine(minrun, percentile2_5, maxrun, percentile2_5);
    line2_5->SetLineStyle(2);
    line2_5->SetLineColor(kGreen);
    line2_5->SetLineWidth(2);
    line2_5->Draw();
    TLatex *latex2_5 = new TLatex(minrun + 10, percentile2_5 * 0.99, "2.5^{th} percentile");
    latex2_5->SetTextColor(kGreen);
    latex2_5->SetTextSize(0.04);
    latex2_5->Draw();
    // print the name on plot
    TLatex *title = new TLatex(0.25, 0.9, name.c_str());
    title->SetNDC();
    title->SetTextSize(0.05);
    title->SetTextFont(42);
    title->Draw();
    


    c->SaveAs((savePath + name + ".pdf").c_str());
}

void plot_rbrQA()
{
    init_plot();
    string savePath = "figures/";
    float minrun = 48200;
    float maxrun = 48300;

    TFile *fdata = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/rbrQA.root");
    TFile *fdata_nc = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/rbrQA_nc.root");

    TGraphErrors *gr_common = (TGraphErrors *)fdata->Get("gr_common");
    TGraphErrors *gr_tight_iso = (TGraphErrors *)fdata->Get("gr_tight_iso");
    TGraphErrors *gr_tight_noniso = (TGraphErrors *)fdata->Get("gr_tight_noniso");
    TGraphErrors *gr_nontight_iso = (TGraphErrors *)fdata->Get("gr_nontight_iso");
    TGraphErrors *gr_nontight_noniso = (TGraphErrors *)fdata->Get("gr_nontight_noniso");

    TGraphErrors *gr_common_nc = (TGraphErrors *)fdata_nc->Get("gr_common")->Clone("gr_common_nc");

    // Original plot for gr_common and gr_common_nc
    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    TH1D *frame = new TH1D("frame", ";Run Number;Common Clusters", 1, minrun, maxrun);
    frame->GetYaxis()->SetRangeUser(60000, 140000);
    frame->GetXaxis()->SetNdivisions(505);
    frame->Draw("axis");

    gr_common->SetMarkerColor(kBlack);
    gr_common->SetMarkerStyle(20);
    gr_common->Draw("P same");

    gr_common_nc->SetMarkerColor(kRed);
    gr_common_nc->SetMarkerStyle(20);
    gr_common_nc->Draw("P same");

    myMarkerLineText(0.25, 0.9, 1, kBlack, 20, kBlack, 1, "w/ pileup correction", 0.05, true);
    myMarkerLineText(0.25, 0.85, 1, kRed, 20, kRed, 1, "w/o pileup correction", 0.05, true);

    c1->SaveAs((savePath + "common_clusters_vs_run.pdf").c_str());

    // Call analyze_and_plot for each additional graph
    analyze_and_plot(gr_tight_iso, "tight_iso", savePath, minrun, maxrun);
    analyze_and_plot(gr_tight_noniso, "tight_noniso", savePath, minrun, maxrun);
    analyze_and_plot(gr_nontight_iso, "nontight_iso", savePath, minrun, maxrun);
    analyze_and_plot(gr_nontight_noniso, "nontight_noniso", savePath, minrun, maxrun);

    delete fdata;
    delete fdata_nc;
}
