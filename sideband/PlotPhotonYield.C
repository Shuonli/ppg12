#include "/sphenix/u/shuhang98/AtlasUtils.C"
#include "/sphenix/u/shuhang98/AtlasStyle.C"

double myfunc(double *x, double *params) {
    double NsigA = x[0];  // Independent variable (N_A^sig)

    // Parameters passed to the function
    double NA = params[0];
    double NB = params[1];
    double NC = params[2];
    double ND = params[3];
    double cB = params[4];
    double cC = params[5];
    double cD = params[6];
    double R = params[7];

    // The equation:
    double numerator = NB - cB * NsigA;
    double denominator = ND - cD * NsigA;

    // Avoid division by zero
    if (denominator == 0) {
        denominator = 1e-12;
    }

    double fraction = (NC - cC * NsigA) / denominator;

    // The function to find the root of
    double result = NsigA - (NA - R * numerator * fraction);
    return result;
}

void addhisto(TH1D* &h1, TH1D* h2){
    if(!h1){
         h1 = (TH1D*)h2->Clone();
         //std::cout << "h1: " << h1 << std::endl;
    }
    else h1->Add(h2);
}

void PlotPhotonYield(){
    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    std::string nodename = "CLUSTERINFO_CEMC_CNN_single";

    std::string datainfile = "/sphenix/user/shuhangli/pi0pythiasim/macro/isophotondatahisto/photon_histo_" + nodename + ".root";
    std::string simsignalinfile = "/sphenix/user/shuhangli/pi0pythiasim/macro/isophotonsimhisto/sim_photon_histo_" + nodename + ".root";
    std::string siminclusiveinfile = "/sphenix/user/shuhangli/pi0pythiasim/macro/isophotonsimhisto/inclusive_photon_histo_" + nodename + ".root";

    TFile* fdata = new TFile(datainfile.c_str());
    TFile* fsimsignal = new TFile(simsignalinfile.c_str());
    TFile* fsiminclusive = new TFile(siminclusiveinfile.c_str());

    std::vector<float> clusterpTbins = {10, 12, 15, 20, 25};

    std::pair<int, int> nontightrange = {2, 5};
    std::pair<int, int> tightrange = {8, 10};

    std::pair<float, float> ETisorange = {-5, 3};
    std::pair<float, float> ETnoisorange = {5, 30};
    
    //const int npTbins = clusterpTbins.size() - 1;
    const int npTbins = 4;

    TH1D* hdataiso[npTbins] = {nullptr};
    TH1D* hdatatightiso[npTbins] = {nullptr};
    TH1D* hdatanontightiso[npTbins] = {nullptr};
    TH1D* hsimsignaliso[npTbins] = {nullptr};
    TH1D* hsimsignaltightiso[npTbins] = {nullptr};
    TH1D* hsimsignalnontightiso[npTbins] = {nullptr};
    TH1D* hsimbackgroundiso[npTbins] = {nullptr};
    TH1D* hsimbackgroundtightiso[npTbins] = {nullptr};
    TH1D* hsimbackgroundnontightiso[npTbins] = {nullptr};

    TH1D* hsiminclusivealliso[npTbins] = {nullptr};
    TH1D* hsiminclusivealltightiso[npTbins] = {nullptr};
    TH1D* hsiminclusiveallnontightiso[npTbins] = {nullptr};
    TH1D* hsiminclusivesignaltightiso[npTbins] = {nullptr};
    TH1D* hsiminclusivesignalnontightiso[npTbins] = {nullptr};
    TH1D* hsiminclusivebackgroundtightiso[npTbins] = {nullptr};
    TH1D* hsiminclusivebackgroundnontightiso[npTbins] = {nullptr};

    TCanvas* canvastight_cluster_iso[npTbins];
    TCanvas* canvasnontight_cluster_iso[npTbins];

    float CB[npTbins], CB_err[npTbins], CC[npTbins], CC_err[npTbins], CD[npTbins], CD_err[npTbins], R[npTbins], R_err[npTbins];
    float NA[npTbins], NA_err[npTbins], NB[npTbins], NB_err[npTbins], NC[npTbins], NC_err[npTbins], ND[npTbins], ND_err[npTbins];

    float yield[npTbins], yield_err[npTbins];
    

    std::string outputfile = "demo_photon_yield_" + nodename + "_" + std::to_string(nontightrange.first) + "_" + std::to_string(nontightrange.second) + "_" + std::to_string(tightrange.first) + "_" + std::to_string(tightrange.second) + ".root";
    TFile* foutput = new TFile(outputfile.c_str(), "RECREATE");
    TGraphErrors* gyield = new TGraphErrors();
    gyield->SetName("photonyield");
    gyield->SetTitle("#Photon Yield");
    gyield->GetXaxis()->SetTitle("E_{T} [GeV]");

    TGraphErrors* grawyield = new TGraphErrors();
    grawyield->SetName("rawphotonyield");
    grawyield->SetTitle("#Photon Yield");
    grawyield->GetXaxis()->SetTitle("E_{T} [GeV]");

    for(int i = 0; i< npTbins; i++){
        //std::cout << "pT bin: " << clusterpTbins[i] << " - " << clusterpTbins[i+1] << std::endl;
        for(int iprob = 0; iprob<10; iprob++){
            //std::cout << "prob bin: " << iprob << std::endl;
            //data all
            TH2D* htemp2D = (TH2D*)fdata->Get(Form("h_cluster_iso_04_Et_%d", iprob));
            htemp2D->GetYaxis()->SetRangeUser(clusterpTbins[i], clusterpTbins[i+1]);
            addhisto(hdataiso[i], (TH1D*)htemp2D->ProjectionX()->Clone("hdataiso"));
            //std::cout << "data all: " << hdataiso[i]->GetEntries() << std::endl;
            if(iprob >= nontightrange.first && iprob < nontightrange.second){
                addhisto(hdatanontightiso[i], (TH1D*)htemp2D->ProjectionX()->Clone("hdatanontightiso"));
            }
            if(iprob >= tightrange.first && iprob < tightrange.second){
                addhisto(hdatatightiso[i], (TH1D*)htemp2D->ProjectionX()->Clone("hdatatightiso"));
            }
            //sim all
            TH2D* hsimtemp2D = (TH2D*)fsimsignal->Get(Form("h_cluster_iso_04_Et_%d", iprob));
            hsimtemp2D->GetYaxis()->SetRangeUser(clusterpTbins[i], clusterpTbins[i+1]);
            addhisto(hsimsignaliso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsimsignaliso"));
            if(iprob >= nontightrange.first && iprob < nontightrange.second){
                addhisto(hsimsignalnontightiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsimsignalnontightiso"));
            }
            if(iprob >= tightrange.first && iprob < tightrange.second){
                addhisto(hsimsignaltightiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsimsignaltightiso"));
            }
            //sim inclusive
            hsimtemp2D = (TH2D*)fsiminclusive->Get(Form("h_cluster_iso_04_Et_%d", iprob));
            hsimtemp2D->GetYaxis()->SetRangeUser(clusterpTbins[i], clusterpTbins[i+1]);
            addhisto(hsiminclusivealliso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsiminclusivealliso"));
            if(iprob >= nontightrange.first && iprob < nontightrange.second){
                addhisto(hsiminclusiveallnontightiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsiminclusiveallnontightiso"));
            }
            if(iprob >= tightrange.first && iprob < tightrange.second){
                addhisto(hsiminclusivealltightiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsiminclusivealltightiso"));
            }
            //sim inclusive signal
            hsimtemp2D = (TH2D*)fsiminclusive->Get(Form("h_signal_cluster_iso_04_Et_%d", iprob));
            hsimtemp2D->GetYaxis()->SetRangeUser(clusterpTbins[i], clusterpTbins[i+1]);
            //addhisto(hsiminclusivebackgroundiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsiminclusivesignaliso"));
            if(iprob >= nontightrange.first && iprob < nontightrange.second){
                addhisto(hsiminclusivesignalnontightiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsiminclusivesignalnontightiso"));
            }
            if(iprob >= tightrange.first && iprob < tightrange.second){
                addhisto(hsiminclusivesignaltightiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsiminclusivesignaltightiso"));
            }
            //sim inclusive background
            hsimtemp2D = (TH2D*)fsiminclusive->Get(Form("h_nonphoton_cluster_iso_04_Et_%d", iprob));
            hsimtemp2D->GetYaxis()->SetRangeUser(clusterpTbins[i], clusterpTbins[i+1]);
            //addhisto(hsimbackgroundiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsimbackgroundiso"));
            if(iprob >= nontightrange.first && iprob < nontightrange.second){
                addhisto(hsiminclusivebackgroundnontightiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hinclusivebackgroundnontightiso"));
            }
            if(iprob >= tightrange.first && iprob < tightrange.second){
                addhisto(hsiminclusivebackgroundtightiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hinclusivebackgroundtightiso"));
            }
            //non iso photon is also the background
            hsimtemp2D = (TH2D*)fsiminclusive->Get(Form("h_nonisophoton_cluster_iso_04_E_%d", iprob));
            hsimtemp2D->GetYaxis()->SetRangeUser(clusterpTbins[i], clusterpTbins[i+1]);
            addhisto(hsimbackgroundiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsimbackgroundiso"));
            if(iprob >= nontightrange.first && iprob < nontightrange.second){
                addhisto(hsiminclusivebackgroundnontightiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hinclusivebackgroundnontightiso"));
            }
            if(iprob >= tightrange.first && iprob < tightrange.second){
                addhisto(hsiminclusivebackgroundtightiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hinclusivebackgroundtightiso"));
            }
            //sim signal
            hsimtemp2D = (TH2D*)fsimsignal->Get(Form("h_signal_cluster_iso_04_Et_%d", iprob));
            hsimtemp2D->GetYaxis()->SetRangeUser(clusterpTbins[i], clusterpTbins[i+1]);
            addhisto(hsimsignaliso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsimsignaliso"));
            if(iprob >= nontightrange.first && iprob < nontightrange.second){
                addhisto(hsimsignalnontightiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsimsignalnontightiso"));
            }
            if(iprob >= tightrange.first && iprob < tightrange.second){
                addhisto(hsimsignaltightiso[i], (TH1D*)hsimtemp2D->ProjectionX()->Clone("hsimsignaltightiso"));
            }

        }
        float N_sig = hsimsignaliso[i]->Integral();
        //now we can calculate the parameters
        float NA_sig = hsimsignaltightiso[i]->Integral(hsimsignaltightiso[i]->FindBin(ETisorange.first), hsimsignaltightiso[i]->FindBin(ETisorange.second));
        //B is tight but not iso, ETnoisorange
        float NB_sig = hsimsignaltightiso[i]->Integral(hsimsignaltightiso[i]->FindBin(ETnoisorange.first), hsimsignaltightiso[i]->FindBin(ETnoisorange.second));
        //C is non-tight and iso, ETisorange
        float NC_sig = hsimsignalnontightiso[i]->Integral(hsimsignalnontightiso[i]->FindBin(ETisorange.first), hsimsignalnontightiso[i]->FindBin(ETisorange.second));
        //D is non-tight and not iso, ETnoisorange
        float ND_sig = hsimsignalnontightiso[i]->Integral(hsimsignalnontightiso[i]->FindBin(ETnoisorange.first), hsimsignalnontightiso[i]->FindBin(ETnoisorange.second));
        //CK is NK_sig / NA_sig
        CB[i] = NB_sig / NA_sig;
        CC[i] = NC_sig / NA_sig;
        CD[i] = ND_sig / NA_sig;
        //photon id efficiency is NA_sig + NB_sig / (sum of all)
        float IDeff = (NA_sig + NB_sig) / N_sig;
        //float IDeff = NA_sig /(NC_sig + NA_sig);
        CB_err[i] = CB[i] * sqrt(1/NB_sig + 1/NA_sig);
        CC_err[i] = CC[i] * sqrt(1/NC_sig + 1/NA_sig);
        CD_err[i] = CD[i] * sqrt(1/ND_sig + 1/NA_sig);
        //from data
        NA[i] = hdatatightiso[i]->Integral(hdatatightiso[i]->FindBin(ETisorange.first), hdatatightiso[i]->FindBin(ETisorange.second));
        NB[i] = hdatatightiso[i]->Integral(hdatatightiso[i]->FindBin(ETnoisorange.first), hdatatightiso[i]->FindBin(ETnoisorange.second));
        NC[i] = hdatanontightiso[i]->Integral(hdatanontightiso[i]->FindBin(ETisorange.first), hdatanontightiso[i]->FindBin(ETisorange.second));
        ND[i] = hdatanontightiso[i]->Integral(hdatanontightiso[i]->FindBin(ETnoisorange.first), hdatanontightiso[i]->FindBin(ETnoisorange.second));

        NA_err[i] = sqrt(NA[i]);
        NB_err[i] = sqrt(NB[i]);
        NC_err[i] = sqrt(NC[i]);
        ND_err[i] = sqrt(ND[i]);

        double xmin = 0.; 
        double xmax = 1000000;

        //std::cout << "NA: " << NA[i] << " +/- " << NA_err[i] << std::endl;
        //std::cout << "NB: " << NB[i] << " +/- " << NB_err[i] << std::endl;
        //std::cout << "NC: " << NC[i] << " +/- " << NC_err[i] << std::endl;
        //std::cout << "ND: " << ND[i] << " +/- " << ND_err[i] << std::endl;
        //std::cout << "CB: " << CB[i] << " +/- " << CB_err[i] << std::endl;
        //std::cout << "CC: " << CC[i] << " +/- " << CC_err[i] << std::endl;
        //std::cout << "CD: " << CD[i] << " +/- " << CD_err[i] << std::endl;
        //std::cout << "IDeff: " << IDeff << std::endl;

        TF1 *f = new TF1("myfunc", myfunc, 0.0, xmax, 8); 
        f->SetParameter(0, NA[i]);
        f->SetParameter(1, NB[i]);
        f->SetParameter(2, NC[i]);
        f->SetParameter(3, ND[i]);
        f->SetParameter(4, CB[i]);
        f->SetParameter(5, CC[i]);
        f->SetParameter(6, CD[i]);
        //f->SetParameter(7, R_cluster[i]);
        f->SetParameter(7, 1);

        double root = f->GetX(0, 0, NA[i]);
        root /= IDeff;
        std::cout << "root: " << root <<" with efficiency: " << IDeff << std::endl;
        yield[i] = root;
        //calculate uncertainty
        int nsamples = 1000;
        //calculate the error from resampling
        TH1D* h_N = new TH1D("h_N", "h_N", 10, -10000, 10000);
        TRandom3* r = new TRandom3(0);
         for (int j = 0; j < nsamples; j++){
            float NA_sampled = r->Gaus(NA[i], NA_err[i]);
            float NB_sampled = r->Gaus(NB[i], NB_err[i]);
            float NC_sampled = r->Gaus(NC[i], NC_err[i]);
            float ND_sampled = r->Gaus(ND[i], ND_err[i]);
            float CC_sampled = r->Gaus(CC[i], CC_err[i]);
            float CB_sampled = r->Gaus(CB[i], CB_err[i]);
            float CD_sampled = r->Gaus(CD[i], CD_err[i]);
            //float R_sampled = r->Gaus(R[i], R_err[i]);

            TF1 *f_sampled = new TF1("myfunc", myfunc, 0.0, xmax, 8);
            f_sampled->SetParameter(0, NA_sampled);
            f_sampled->SetParameter(1, NB_sampled);
            f_sampled->SetParameter(2, NC_sampled);
            f_sampled->SetParameter(3, ND_sampled);
            f_sampled->SetParameter(4, CB_sampled);
            f_sampled->SetParameter(5, CC_sampled);
            f_sampled->SetParameter(6, CD_sampled);
            //f_sampled->SetParameter(7, R_sampled);
            f_sampled->SetParameter(7, 1);

            double root_sampled = f_sampled->GetX(0, -20*NA_sampled, 20*NA_sampled);

            root_sampled /= IDeff;

            //if root is not nan
            if(root_sampled == root_sampled){
                h_N->Fill(root_sampled);
            }
         }


        yield_err[i] = h_N->GetRMS();
        std::cout << "yield: " << yield[i] << " +/- " << yield_err[i] << std::endl;

        grawyield->SetPoint(i, (clusterpTbins[i] + clusterpTbins[i+1])/2, NA[i]);
        grawyield->SetPointError(i, (clusterpTbins[i+1] - clusterpTbins[i])/2, NA_err[i]);

    }

    //plotting
    for(int i = 0; i< npTbins; i++){
        gyield->SetPoint(i, (clusterpTbins[i] + clusterpTbins[i+1])/2, yield[i]);
        gyield->SetPointError(i, (clusterpTbins[i+1] - clusterpTbins[i])/2, yield_err[i]);

        
        gStyle->SetCanvasPreferGL(true);
        {
        canvastight_cluster_iso[i] = new TCanvas(Form("canvastight_cluster_iso_%d", i), Form("canvastight_cluster_iso_%d", i), 900, 800);
        canvastight_cluster_iso[i]->SetSupportGL(true);
        hsiminclusivealltightiso[i]->GetXaxis()->SetRangeUser(-5, 20);
        hsiminclusivealltightiso[i]->GetYaxis()->SetRangeUser(0, hsiminclusivealltightiso[i]->GetMaximum()*1.2);
        
        hsiminclusivealltightiso[i]->SetXTitle("E_{T}^{iso} [GeV]");
        hsiminclusivealltightiso[i]->SetYTitle("Counts");
        hsiminclusivealltightiso[i]->SetLineColor(kBlack);
        hsiminclusivealltightiso[i]->SetLineWidth(5);
        //setfill color and transparency
        hsiminclusivesignaltightiso[i]->SetLineColor(kPink+8);
        hsiminclusivesignaltightiso[i]->SetLineWidth(5);
        //hsiminclusivesignaltightiso[i]->SetFillColor(kPink+8);
        hsiminclusivesignaltightiso[i]->SetFillColorAlpha(kPink+8,0.35);
        //hsiminclusivesignaltightiso[i]->SetFillStyle(3001);
        hsiminclusivebackgroundtightiso[i]->SetLineColor(kSpring+3);
        hsiminclusivebackgroundtightiso[i]->SetLineWidth(5);
        //hsiminclusivebackgroundtightiso[i]->SetFillColor(kSpring+3);
        hsiminclusivebackgroundtightiso[i]->SetFillColorAlpha(kSpring+3,0.35);
        //hsiminclusivebackgroundtightiso[i]->SetFillStyle(3001);

        hsiminclusivealltightiso[i]->Draw("hist");
        hsiminclusivesignaltightiso[i]->Draw("SAME hist F");
        hsiminclusivebackgroundtightiso[i]->Draw("SAME hist F");

        //myMarkerLineText(0.55, 0.85, 0, kBlack, 20, kBlack, 1, "Total");
        //myMarkerLineText(0.55, 0.80, 0, kBlack, 20, kPink+8, 1, "Truth Isolated Photon");
        //myMarkerLineText(0.55, 0.75, 0, kBlack, 20, kSpring+3, 1, "Background");
        myBoxTextAlpha(0.57, 0.85, 0.05, kBlack, 0, kBlack, 1, "Total");
        myBoxTextAlpha(0.57, 0.80, 0.05, kPink+8, 0.35, kPink+8, 1, "Signal Photons");
        myBoxTextAlpha(0.57, 0.75, 0.05, kSpring+3, 0.35, kSpring+3, 1, "Background");
        //a blue line at x = 3
     
        myText(0.5, 0.90, kBlack, Form("Cluster E_{T} %.0f - %.0f GeV", clusterpTbins[i], clusterpTbins[i + 1]));
        
        }
        //non-tight
        {
        canvasnontight_cluster_iso[i] = new TCanvas(Form("canvasnontight_cluster_iso_%d", i), Form("canvasnontight_cluster_iso_%d", i), 900, 800);
        canvasnontight_cluster_iso[i]->SetSupportGL(true);
        hsiminclusiveallnontightiso[i]->GetXaxis()->SetRangeUser(-5, 20);
        hsiminclusiveallnontightiso[i]->GetYaxis()->SetRangeUser(0, hsiminclusiveallnontightiso[i]->GetMaximum()*1.2);

        hsiminclusiveallnontightiso[i]->SetXTitle("E_{T}^{iso} [GeV]");
        hsiminclusiveallnontightiso[i]->SetYTitle("Counts");
        hsiminclusiveallnontightiso[i]->SetLineColor(kBlack);
        hsiminclusiveallnontightiso[i]->SetLineWidth(5);
        hsiminclusivesignalnontightiso[i]->SetLineColor(kPink + 8);
        //setfill color and transparency
        hsiminclusivesignalnontightiso[i]->SetLineColor(kPink+8);
        hsiminclusivesignalnontightiso[i]->SetLineWidth(5);        
        //hsiminclusivesignalnontightiso[i]->SetFillColor(kPink+8);
        hsiminclusivesignalnontightiso[i]->SetFillColorAlpha(kPink+8,0.35);
        //hsiminclusivesignalnontightiso[i]->SetFillStyle(3001);
        hsiminclusivebackgroundnontightiso[i]->SetLineColor(kSpring+3);
        hsiminclusivebackgroundnontightiso[i]->SetLineWidth(5);
        //hsiminclusivebackgroundnontightiso[i]->SetFillColor(kSpring+3);
        hsiminclusivebackgroundnontightiso[i]->SetFillColorAlpha(kSpring+3,0.35);
        //hsiminclusivebackgroundnontightiso[i]->SetFillStyle(3001);

        hsiminclusiveallnontightiso[i]->Draw("hist");
        hsiminclusivesignalnontightiso[i]->Draw("SAME hist F");
        hsiminclusivebackgroundnontightiso[i]->Draw("SAME hist F");

        //myMarkerLineText(0.55, 0.85, 0, kBlack, 20, kBlack, 1, "Total");
        //myMarkerLineText(0.55, 0.80, 0, kBlack, 20, kPink+8, 1, "Truth  Isolated Photon");
        //myMarkerLineText(0.55, 0.75, 0, kBlack, 20, kSpring+3, 1, "Background");

        //void myBoxText(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor,Int_t bstyle,Int_t lcolor,Int_t lstyle, const char *text)
        myBoxTextAlpha(0.57, 0.85, 0.05, kBlack, 0, kBlack, 1, "Total");
        myBoxTextAlpha(0.57, 0.80, 0.05, kPink+8, 0.35, kPink+8, 1, "Signal Photons");
        myBoxTextAlpha(0.57, 0.75, 0.05, kSpring+3, 0.35, kSpring+3, 1, "Background");

        myText(0.5, 0.90, kBlack, Form("Cluster E_{T} %.0f - %.0f GeV", clusterpTbins[i], clusterpTbins[i + 1]));

        }
        
    }
    foutput->cd();
    gyield->Write();
    grawyield->Write();
    foutput->Close();
}