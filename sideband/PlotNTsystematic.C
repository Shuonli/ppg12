#include "/sphenix/u/bseidlitz/plotstyle/AtlasUtils.C"
#include "/sphenix/u/shuhang98/AtlasStyle.C"

void PlotNTsystematic() {
    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    // Load the ROOT files
    TFile *file1 = TFile::Open("photon_yield_CLUSTERINFO_CEMC_CNN_single_2_5_8_10.root");
    TFile *file2 = TFile::Open("photon_yield_CLUSTERINFO_CEMC_CNN_single_4_7_8_10.root");

    if (!file1 || !file2 || file1->IsZombie() || file2->IsZombie()) {
        std::cerr << "Error opening one of the files!" << std::endl;
        return;
    }

    // Extract TGraphErrors from both files
    TGraphErrors *graph1 = (TGraphErrors*)file1->Get("photonyield");
    TGraphErrors *graph2 = (TGraphErrors*)file2->Get("photonyield");

    if (!graph1 || !graph2) {
        std::cerr << "Error: TGraphErrors 'photonyield' not found in one of the files!" << std::endl;
        std::cerr << "graph1: " << graph1 << ", graph2: " << graph2 << std::endl;
        return;
    }

    // Create canvas and pads
    TCanvas *c1 = new TCanvas("c1", "Photon Yield Overlay and Ratio", 800, 1200);

    // Upper pad for photon yield overlay
    TPad *pad1 = new TPad("pad1", "Photon Yield", 0, 0.4, 1, 1.0);
    pad1->SetBottomMargin(0.02); // Upper and lower margin for separation
    pad1->Draw();
    pad1->cd();
    
    // Draw the first graph with axis settings
    graph1->SetLineColor(kRed);
    graph1->SetMarkerStyle(20);
    graph1->SetMarkerColor(kRed);
    //set x label offset
    graph1->GetXaxis()->SetLabelOffset(999);
    graph1->GetYaxis()->SetTitle("N_{A}^{sig}");
    graph1->GetYaxis()->SetRangeUser(15, 2800);
    graph1->Draw("AP");

    // Draw the second graph on the same pad
    graph2->SetLineColor(kBlue);
    graph2->SetMarkerStyle(21);
    graph2->SetMarkerColor(kBlue);
    graph2->Draw("P SAME");

    // Add legend
    //TLegend *legend = new TLegend(0.5, 0.65, 0.9, 0.9);
    //legend->AddEntry(graph1, "Non-tight prob:[0.2,0.5]", "lp");
    //legend->AddEntry(graph2, "Non-tight prob:[0.4,0.7]", "lp");
    //legend->SetFillStyle(0);
    //legend->SetBorderSize(0);
    //legend->Draw();

    myMarkerLineText(0.5, 0.85, 1, kBlue, 21, kBlue, 1, "Non-tight prob:[0.4,0.7]");
    myMarkerLineText(0.5, 0.80, 1, kRed, 20, kRed, 1, "Non-tight prob:[0.2,0.5]");

    // Go back to the main canvas and create the lower pad for the ratio plot
    c1->cd();
    TPad *pad2 = new TPad("pad2", "Ratio", 0, 0.0, 1, 0.4);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();

    // Calculate the ratio of the two graphs
    int nPoints1 = graph1->GetN();
    int nPoints2 = graph2->GetN();

    if (nPoints1 != nPoints2) {
        std::cerr << "Error: Graphs have different number of points!" << std::endl;
        return;
    }

    TGraphErrors *ratioGraph = new TGraphErrors(nPoints1);
    
    for (int i = 0; i < nPoints1; i++) {
        double x1, y1, ex1, ey1;
        double x2, y2, ex2, ey2;

        graph1->GetPoint(i, x1, y1);
        graph2->GetPoint(i, x2, y2);

        ex1 = graph1->GetErrorX(i);
        ey1 = graph1->GetErrorY(i);
        ex2 = graph2->GetErrorX(i);
        ey2 = graph2->GetErrorY(i);

        if (y2 != 0) {
            double ratio = y1 / y2;
            double ratioError = ratio * sqrt((ey1 / y1) * (ey1 / y1) + (ey2 / y2) * (ey2 / y2));

            ratioGraph->SetPoint(i, x1, ratio);
            ratioGraph->SetPointError(i, ex1, ratioError);
        } else {
            ratioGraph->SetPoint(i, x1, 0);
            ratioGraph->SetPointError(i, ex1, 0);
        }
    }

    // Draw the ratio plot
    ratioGraph->SetLineColor(kBlack);
    ratioGraph->SetMarkerStyle(20);
    ratioGraph->SetMarkerColor(kBlack);
    ratioGraph->GetYaxis()->SetRangeUser(0., 1.5);
    ratioGraph->Draw("AP");
    //set title size larger
    ratioGraph->GetXaxis()->SetTitleSize(0.07);
    ratioGraph->GetXaxis()->SetLabelSize(0.07);
    ratioGraph->GetYaxis()->SetTitleOffset(1.);
    ratioGraph->GetYaxis()->SetLabelSize(0.07);
    ratioGraph->GetYaxis()->SetTitleSize(0.07);

    ratioGraph->GetXaxis()->SetTitle("Raw E_{T} [GeV]");
    ratioGraph->GetYaxis()->SetTitle("Ratio");

    //a line at y = 1
    TLine *line = new TLine(ratioGraph->GetXaxis()->GetXmin(), 1, ratioGraph->GetXaxis()->GetXmax(), 1);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->Draw();

    // Update and save the canvas
    c1->Update();
    //c1->SaveAs("photon_yield_overlay_and_ratio.png");

    // Clean up
    //delete file1;
    //delete file2;
    //delete c1;
}
