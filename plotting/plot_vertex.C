#include "plotcommon.h"

void plot_vertex()
{
    init_plot();
    TFile *f_mc = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_nom.root", "READ");
    TFile *f_data = new TFile("/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_nom.root", "READ");


    TH1F * h_vertex_z_mc = (TH1F *)f_mc->Get("h_vertexz");
    TH1F * h_vertex_z_data = (TH1F *)f_data->Get("h_vertexz");

    h_vertex_z_mc->Scale(1.0 / h_vertex_z_mc->Integral());
    h_vertex_z_data->Scale(1.0 / h_vertex_z_data->Integral());



    h_vertex_z_mc->SetLineColor(kRed);
    h_vertex_z_data->SetLineColor(kBlue);

    TH1F* h_vertex_z_data_mc = (TH1F*)h_vertex_z_data->Clone("h_vertex_z_data_mc");
    h_vertex_z_data_mc->Divide(h_vertex_z_mc);
    //fit with a gaussian
    TF1 *f1 = new TF1("f1", "gaus", -30, 30);
    h_vertex_z_data_mc->Fit("f1", "REM", "", -30, 30);

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    h_vertex_z_data->Draw("hist");
    h_vertex_z_mc->Draw("hist same");

    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
    
    
    h_vertex_z_data_mc->Draw("hist");
    f1->Draw("same");
    



}