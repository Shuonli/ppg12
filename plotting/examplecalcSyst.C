
#include "plotCommon.h"
#include "systemDept.h"

void calcSyst(){
    initStyle();

    string figPath = "../figures/13TeV";

    TFile* fin_nch = new TFile("../../Correlation/track_cluster/rootFiles/13TeV2/corr_nominal.root");
    TH1F* h_Nch_raw = (TH1F*) fin_nch->Get("h_Nch_raw");
    float nch_avg[Nbins_out_Nch];
    calcAvgNch(h_Nch_raw,Nbins_out_Nch,dNch_out_range,nch_avg);

    bool sym_mix    = 1;
    bool sym_NoExcl = 0;
    bool sym_bin    = 0;
    bool sym_Mu    = 0;
    bool sym_ref    = 1;
    bool sym_sum    = 0;
    bool sym_eff    = 0;

    TFile* fin  = new TFile("../output/13TeV/ANA_Results_nominal.root","READ");
    //TFile* fin3 = new TFile("./RootFiles/_result_nominalv3.root","READ");

    // Files with systematic histograms
    TFile* fin_ref = new TFile("../output/13TeV/Syst_Ref_SS.root", "READ");
    TFile* fin_Mu  = new TFile("../output/13TeV/Syst_Mu.root", "READ");
    TFile* fin_FSR = new TFile("../output/13TeV/Syst_D1Fstar.root", "READ");
    TFile* fin_D1S = new TFile("../output/13TeV/Syst_D1scale.root", "READ");
    TFile* fin_LMG = new TFile("../output/13TeV/Syst_Ref_SS.root","READ");
    TFile* fin_mix = new TFile("../output/13TeV/Syst_Mix.root", "READ");
    TFile* fin_bin = new TFile("../output/13TeV/Syst_Bin.root", "READ");
    TFile* fin_eff = new TFile("../output/13TeV/Syst_NoEff.root", "READ");


    TFile* fout = new TFile("../output/13TeV/Results_withSyst.root","RECREATE");

    TGraphAsymmErrors* g_ob_syst;

    TH1F* hsys_D1S_low;
    TH1F* hsys_D1S_high;
    TH1F* hsys_LMG_low;
    TH1F* hsys_LMG_high;
    TH1F* hsys_Mu_low;
    TH1F* hsys_Mu_high;
    TH1F* hsys_ref_low;
    TH1F* hsys_ref_high;
    TH1F* hsys_FSR_low;
    TH1F* hsys_FSR_high;
    TH1F* hsys_mix_low;
    TH1F* hsys_mix_high;
    TH1F* hsys_bin_low;
    TH1F* hsys_bin_high;
    TH1F* hsys_eff_low;
    TH1F* hsys_eff_high;

    TH1F* hsys_sum_low;
    TH1F* hsys_sum_high;


    // v11Raw_Fn
    {
        string obs = "v11Raw_Fn";

        TH1F* h_ob = (TH1F*) fin->Get(Form("h_%s",obs.c_str()));
        TGraphAsymmErrors* g_ob_stat =  new TGraphAsymmErrors(h_ob);

        TGraphAsymmErrors* g_ob      = (TGraphAsymmErrors*) g_ob_stat->Clone( Form("g_%s_StatError",obs.c_str()) );
        TGraphAsymmErrors* g_ob_syst = (TGraphAsymmErrors*) g_ob->Clone( Form("g_%s_SystError",obs.c_str()) );

        // asymmetric errors
        hsys_Mu_low  = (TH1F*)fin_Mu->Get(Form("h_%s_sysLow", obs.c_str()) );
        hsys_Mu_high = (TH1F*)fin_Mu->Get(Form("h_%s_sysHigh",obs.c_str()) );
        hsys_Mu_low ->SetMarkerColor(2);
        hsys_Mu_low ->SetMarkerStyle(22);
        hsys_Mu_high ->SetMarkerColor(2);
        hsys_Mu_high ->SetMarkerStyle(22);

        hsys_ref_low  = (TH1F*)fin_ref->Get(Form("h_%s_sysLow", obs.c_str()) );
        hsys_ref_high = (TH1F*)fin_ref->Get(Form("h_%s_sysHigh",obs.c_str()) );
        hsys_ref_low ->SetLineColor(4);
        hsys_ref_low ->SetLineStyle(2);
        hsys_ref_high ->SetLineColor(4);
        hsys_ref_high ->SetLineStyle(2);

        hsys_D1S_low  = (TH1F*)fin_D1S->Get(Form("h_%s_sysLow", obs.c_str()) );
        hsys_D1S_high = (TH1F*)fin_D1S->Get(Form("h_%s_sysHigh",obs.c_str()) );
        hsys_D1S_low ->SetLineColor(kMagenta+1);
        hsys_D1S_low ->SetLineStyle(5);
        hsys_D1S_high ->SetLineColor(kMagenta+1);
        hsys_D1S_high ->SetLineStyle(5);

        hsys_FSR_low  = (TH1F*)fin_FSR->Get( Form("h_%s_sysLow", obs.c_str()));
        hsys_FSR_high = (TH1F*)fin_FSR->Get( Form("h_%s_sysHigh",obs.c_str()));
        hsys_FSR_low ->SetLineColor(kSpring+4);
        hsys_FSR_low ->SetLineStyle(3);
        hsys_FSR_high->SetLineColor(kSpring+4);
        hsys_FSR_high->SetLineStyle(3);

        hsys_mix_low   = (TH1F*)fin_mix->Get( Form("h_%s_sysLow", obs.c_str()));
        hsys_mix_high  = (TH1F*)fin_mix->Get( Form("h_%s_sysHigh",obs.c_str()));
        hsys_mix_low->SetLineColor(kRed);
        hsys_mix_low->SetLineStyle(2);
        hsys_mix_high->SetLineColor(kRed);
        hsys_mix_high->SetLineStyle(2);

        hsys_eff_low   = (TH1F*)fin_eff->Get( Form("h_%s_sysLow", obs.c_str()));
        hsys_eff_high  = (TH1F*)fin_eff->Get( Form("h_%s_sysHigh",obs.c_str()));
        hsys_eff_low->SetLineColor(kSpring+4);
        hsys_eff_low->SetLineStyle(2);
        hsys_eff_high->SetLineColor(kSpring+4);
        hsys_eff_high->SetLineStyle(2);

        // symmetric errors
        hsys_LMG_low  = (TH1F*)fin_LMG->Get( Form("h_%s_sysSym", obs.c_str()));
        hsys_LMG_low ->SetLineColor(kOrange-2);
        hsys_LMG_low ->SetLineStyle(1);
        hsys_LMG_high = (TH1F*)hsys_LMG_low->Clone( "hsys_LMG_high");


        hsys_bin_low  = (TH1F*)fin_bin->Get( Form("h_%s_sysSym", obs.c_str()));
        hsys_bin_low->SetLineColor(kMagenta);
        hsys_bin_low->SetLineStyle(1);
        hsys_bin_high = (TH1F*)hsys_bin_low->Clone( "hsys_bin_high");

        /// combined systematics
        hsys_sum_low = (TH1F*)hsys_Mu_low->Clone( Form("hsys_sum_%s_low", obs.c_str()) );
        hsys_sum_low->Reset();
        hsys_sum_low->SetLineColor(1);
        hsys_sum_low->SetLineStyle(1);
        hsys_sum_low->SetLineWidth(2);

        hsys_sum_high = (TH1F*)hsys_Mu_high->Clone( Form("hsys_sum_%s_high", obs.c_str()) );
        hsys_sum_high->Reset();
        hsys_sum_high->SetLineColor(1);
        hsys_sum_high->SetLineStyle(1);
        hsys_sum_high->SetLineWidth(2);


        for (int i=1; i< hsys_sum_low->GetNbinsX()+1; i++) {
            float error_D1S_low = hsys_D1S_low->GetBinContent(i);
            float error_D1S_high = hsys_D1S_high->GetBinContent(i);

            float error_LMG_low = hsys_LMG_low->GetBinContent(i);
            float error_LMG_high = hsys_LMG_high->GetBinContent(i);


            float error_Mu_low = hsys_Mu_low->GetBinContent(i);
            float error_Mu_high = hsys_Mu_high->GetBinContent(i);
            if (sym_Mu){
              error_Mu_low  = TMath::Max(error_Mu_low,error_Mu_high);
              error_Mu_high = TMath::Max(error_Mu_low,error_Mu_high);
            }

            float error_FSR_low = hsys_FSR_low->GetBinContent(i);
            float error_FSR_high = hsys_FSR_high->GetBinContent(i);

            float error_ref_low = hsys_ref_low->GetBinContent(i);
            float error_ref_high = hsys_ref_high->GetBinContent(i);
            if (sym_ref){
              error_ref_low  = TMath::Max(error_ref_low,error_ref_high);
              error_ref_high = TMath::Max(error_ref_low,error_ref_high);
            }

            float error_eff_low = hsys_eff_low->GetBinContent(i);
            float error_eff_high = hsys_eff_high->GetBinContent(i);
            if (sym_eff){
              error_eff_low  = TMath::Max(error_eff_low,error_eff_high);
              error_eff_high = TMath::Max(error_eff_low,error_eff_high);
            }

            float error_mix_low = hsys_mix_low->GetBinContent(i);
            float error_mix_high = hsys_mix_high->GetBinContent(i);
            if (sym_mix){
              error_mix_low  = TMath::Max(error_mix_low,error_mix_high);
              error_mix_high = TMath::Max(error_mix_low,error_mix_high);
            }

            float error_bin_low = hsys_bin_low->GetBinContent(i);
            float error_bin_high = hsys_bin_high->GetBinContent(i);
            if (sym_bin){
              error_bin_low  = TMath::Max(error_bin_low,error_bin_high);
              error_bin_high = TMath::Max(error_bin_low,error_bin_high);
            }


            // absolute errors
            float error_low = sqrt( pow(error_bin_low,2) +
                                    pow(error_mix_low,2) +
                                    //pow(error_LMG_low,2) +
                                    pow(error_eff_low,2) +
                                    //pow(error_ref_low,2)
                                    pow(error_Mu_low,2)
                                    //pow(error_FSR_low,2) +
                                    //pow(error_D1S_low,2) +
                                  );

            float error_high = sqrt( pow(error_bin_high,2) +
                                     pow(error_mix_high,2) +
                                    // pow(error_LMG_high,2) +
                                     pow(error_eff_high,2) +
                                     //pow(error_ref_high,2)
                                     pow(error_Mu_high,2)
                                     //pow(error_FSR_high,2) +
                                     //pow(error_D1S_high,2) +
                                  );


            if (sym_sum){
              error_low  = TMath::Max(error_low,error_high);
              error_high = TMath::Max(error_low,error_high);
            }

            float v2_nominal = g_ob->GetY()[i-1];
            float v2_StatErr = g_ob->GetEYhigh()[i-1];

            g_ob_syst->GetEYlow()[i-1]  = error_low;
            g_ob_syst->GetEYhigh()[i-1] = error_high;

            g_ob_syst->GetEXlow()[i-1]  = 0.08;
            g_ob_syst->GetEXhigh()[i-1] = 0.08;


            hsys_D1S_low->SetBinContent(i, -1*fabs(error_D1S_low));
            hsys_LMG_low->SetBinContent(i, -1*fabs(error_LMG_low));
            hsys_Mu_low->SetBinContent(i, -1*fabs(error_Mu_low));
            hsys_FSR_low->SetBinContent(i, -1*fabs(error_FSR_low));
            hsys_ref_low->SetBinContent(i, -1*fabs(error_ref_low));
            hsys_eff_low->SetBinContent(i, -1*fabs(error_eff_low));
            hsys_bin_low->SetBinContent(i, -1*fabs(error_bin_low));
            hsys_mix_low->SetBinContent(i, -1*fabs(error_mix_low) );


            hsys_D1S_high->SetBinContent(i, fabs(error_D1S_high ));
            hsys_LMG_high->SetBinContent(i, fabs(error_LMG_high ));
            hsys_Mu_high->SetBinContent(i, fabs(error_Mu_high ));
            hsys_FSR_high->SetBinContent(i, fabs(error_FSR_high ));
            hsys_ref_high->SetBinContent(i, fabs(error_ref_high ));
            hsys_eff_high->SetBinContent(i, fabs(error_eff_high ));
            hsys_mix_high->SetBinContent(i, fabs(error_mix_high ));
            hsys_bin_high->SetBinContent(i, fabs(error_bin_high ));


            hsys_sum_low->SetBinContent(i, -1*fabs(error_low ) );
            hsys_sum_high->SetBinContent(i,fabs(error_high   ) );

        }


        TCanvas* c1 = new TCanvas("c1","New Canvas",50,50,600,600);
        h_frame_err_Fn->Draw("AXIS");
        h_frame_err_Fn->GetYaxis()->SetRangeUser(-0.01,0.02);

        hsys_sum_low->SetLineWidth(8);

        hsys_sum_low->Draw("HIST ][ SAME");
        hsys_bin_low->Draw("HIST ][ SAME");
        hsys_eff_low->Draw("HIST ][ SAME");
        //hsys_ref_low->Draw("HIST ][ SAME");
        hsys_Mu_low->Draw("PHIST ][  SAME");
        //  hsys_LMG_low->Draw("HIST ][ SAME");
        hsys_mix_low->Draw("HIST ][ SAME");
        //hsys_FSR_low->Draw("HIST ][ SAME");
        //hsys_D1S_low->Draw("HIST ][ SAME");

        hsys_sum_high->SetLineWidth(8);
        hsys_sum_high->Draw("HIST ][ SAME");
        hsys_bin_high->Draw("HIST ][ SAME");
        hsys_eff_high->Draw("HIST ][ SAME");
        //hsys_ref_high->Draw("HIST ][ SAME");
        hsys_Mu_high->Draw("PHIST ][  SAME");
        //  hsys_LMG_high->Draw("HISTSAME");
        hsys_mix_high->Draw("HIST ][ SAME");
        //hsys_FSR_high->Draw("HISTSAME");
        //hsys_D1S_high->Draw("HISTSAME");


        myText(          0.6,0.88+0.02, 1,"#font[72]{ATLAS} Internal",0.04);
        myText(          0.6,0.84+0.02, 1,sysTag.c_str(),0.04);
        myText(          0.6,0.80+0.02, 1, "#it{F}_{1} systematics",0.04);
        myMarkerLineText(0.25,0.88+0.02,0.0, 1, 1, 1, 1, "Total Uncertainty", 0.04, true);
        myMarkerLineText(0.25,0.84+0.02,1.5, kRed   , 22, kRed, 1,"#it{#mu}", 0.04, true);
        //myMarkerLineText(0.25,0.84+0.02,0.0, 4, 2, 4, 2, "LM Reference sym.", 0.04);
        myMarkerLineText(0.25,0.80+0.02,0.0, kMagenta, 1, kMagenta, 1,"#Delta#phi binning sym.", 0.04, true);
        myMarkerLineText(0.25,0.76+0.02,0.0, kSpring+4, 1, kSpring+4, 2,"tracking efficiency", 0.04);
      //  myMarkerLineText(0.25,0.72+0.02,0.0, kOrange-2, 1,kOrange-2, 1,"#it{F}_{2}^{LM} input",  0.04, true);
      myMarkerLineText(0.25,0.72+0.02,0.0, kRed     , 1,kRed      , 2,"Mixed Event",  0.04, true);



        c1->SaveAs(Form("%s/systBreak_%s.pdf",figPath.c_str(),obs.c_str()));

        fout->cd();
        g_ob     ->Write();
        g_ob_syst->Write();

    }







    // v22Raw_Fn
    {
        string obs = "v22Raw_Fn";

        TH1F* h_ob = (TH1F*) fin->Get(Form("h_%s",obs.c_str()));
        TGraphAsymmErrors* g_ob_stat =  new TGraphAsymmErrors(h_ob);

        TGraphAsymmErrors* g_ob      = (TGraphAsymmErrors*) g_ob_stat->Clone( Form("g_%s_StatError",obs.c_str()) );
        TGraphAsymmErrors* g_ob_syst = (TGraphAsymmErrors*) g_ob->Clone( Form("g_%s_SystError",obs.c_str()) );

        // asymmetric errors
        hsys_Mu_low  = (TH1F*)fin_Mu->Get(Form("h_%s_sysLow", obs.c_str()) );
        hsys_Mu_high = (TH1F*)fin_Mu->Get(Form("h_%s_sysHigh",obs.c_str()) );
        hsys_Mu_low ->SetMarkerColor(2);
        hsys_Mu_low ->SetMarkerStyle(22);
        hsys_Mu_high ->SetMarkerColor(2);
        hsys_Mu_high ->SetMarkerStyle(22);

        hsys_ref_low  = (TH1F*)fin_ref->Get(Form("h_%s_sysLow", obs.c_str()) );
        hsys_ref_high = (TH1F*)fin_ref->Get(Form("h_%s_sysHigh",obs.c_str()) );
        hsys_ref_low ->SetLineColor(4);
        hsys_ref_low ->SetLineStyle(2);
        hsys_ref_high ->SetLineColor(4);
        hsys_ref_high ->SetLineStyle(2);

        hsys_D1S_low  = (TH1F*)fin_D1S->Get(Form("h_%s_sysLow", obs.c_str()) );
        hsys_D1S_high = (TH1F*)fin_D1S->Get(Form("h_%s_sysHigh",obs.c_str()) );
        hsys_D1S_low ->SetLineColor(kMagenta+1);
        hsys_D1S_low ->SetLineStyle(5);
        hsys_D1S_high ->SetLineColor(kMagenta+1);
        hsys_D1S_high ->SetLineStyle(5);

        hsys_FSR_low  = (TH1F*)fin_FSR->Get( Form("h_%s_sysLow", obs.c_str()));
        hsys_FSR_high = (TH1F*)fin_FSR->Get( Form("h_%s_sysHigh",obs.c_str()));
        hsys_FSR_low ->SetLineColor(kSpring+4);
        hsys_FSR_low ->SetLineStyle(3);
        hsys_FSR_high->SetLineColor(kSpring+4);
        hsys_FSR_high->SetLineStyle(3);

        hsys_mix_low   = (TH1F*)fin_mix->Get( Form("h_%s_sysLow", obs.c_str()));
        hsys_mix_high  = (TH1F*)fin_mix->Get( Form("h_%s_sysHigh",obs.c_str()));
        hsys_mix_low->SetLineColor(kRed);
        hsys_mix_low->SetLineStyle(2);
        hsys_mix_high->SetLineColor(kRed);
        hsys_mix_high->SetLineStyle(2);

        hsys_eff_low   = (TH1F*)fin_eff->Get( Form("h_%s_sysLow", obs.c_str()));
        hsys_eff_high  = (TH1F*)fin_eff->Get( Form("h_%s_sysHigh",obs.c_str()));
        hsys_eff_low->SetLineColor(kSpring+4);
        hsys_eff_low->SetLineStyle(2);
        hsys_eff_high->SetLineColor(kSpring+4);
        hsys_eff_high->SetLineStyle(2);

        // symmetric errors
        hsys_LMG_low  = (TH1F*)fin_LMG->Get( Form("h_%s_sysSym", obs.c_str()));
        hsys_LMG_low ->SetLineColor(kOrange-2);
        hsys_LMG_low ->SetLineStyle(1);
        hsys_LMG_high = (TH1F*)hsys_LMG_low->Clone( "hsys_LMG_high");


        hsys_bin_low  = (TH1F*)fin_bin->Get( Form("h_%s_sysSym", obs.c_str()));
        hsys_bin_low->SetLineColor(kMagenta);
        hsys_bin_low->SetLineStyle(1);
        hsys_bin_high = (TH1F*)hsys_bin_low->Clone( "hsys_bin_high");

        /// combined systematics
        hsys_sum_low = (TH1F*)hsys_Mu_low->Clone( Form("hsys_sum_%s_low", obs.c_str()) );
        hsys_sum_low->Reset();
        hsys_sum_low->SetLineColor(1);
        hsys_sum_low->SetLineStyle(1);
        hsys_sum_low->SetLineWidth(2);

        hsys_sum_high = (TH1F*)hsys_Mu_high->Clone( Form("hsys_sum_%s_high", obs.c_str()) );
        hsys_sum_high->Reset();
        hsys_sum_high->SetLineColor(1);
        hsys_sum_high->SetLineStyle(1);
        hsys_sum_high->SetLineWidth(2);


        for (int i=1; i< hsys_sum_low->GetNbinsX()+1; i++) {
            float error_D1S_low = hsys_D1S_low->GetBinContent(i);
            float error_D1S_high = hsys_D1S_high->GetBinContent(i);

            float error_LMG_low = hsys_LMG_low->GetBinContent(i);
            float error_LMG_high = hsys_LMG_high->GetBinContent(i);


            float error_Mu_low = hsys_Mu_low->GetBinContent(i);
            float error_Mu_high = hsys_Mu_high->GetBinContent(i);
            if (sym_Mu){
              error_Mu_low  = TMath::Max(error_Mu_low,error_Mu_high);
              error_Mu_high = TMath::Max(error_Mu_low,error_Mu_high);
            }

            float error_FSR_low = hsys_FSR_low->GetBinContent(i);
            float error_FSR_high = hsys_FSR_high->GetBinContent(i);

            float error_ref_low = hsys_ref_low->GetBinContent(i);
            float error_ref_high = hsys_ref_high->GetBinContent(i);
            if (sym_ref){
              error_ref_low  = TMath::Max(error_ref_low,error_ref_high);
              error_ref_high = TMath::Max(error_ref_low,error_ref_high);
            }

            float error_eff_low = hsys_eff_low->GetBinContent(i);
            float error_eff_high = hsys_eff_high->GetBinContent(i);
            if (sym_eff){
              error_eff_low  = TMath::Max(error_eff_low,error_eff_high);
              error_eff_high = TMath::Max(error_eff_low,error_eff_high);
            }

            float error_mix_low = hsys_mix_low->GetBinContent(i);
            float error_mix_high = hsys_mix_high->GetBinContent(i);
            if (sym_mix){
              error_mix_low  = TMath::Max(error_mix_low,error_mix_high);
              error_mix_high = TMath::Max(error_mix_low,error_mix_high);
            }

            float error_bin_low = hsys_bin_low->GetBinContent(i);
            float error_bin_high = hsys_bin_high->GetBinContent(i);
            if (sym_bin){
              error_bin_low  = TMath::Max(error_bin_low,error_bin_high);
              error_bin_high = TMath::Max(error_bin_low,error_bin_high);
            }


            // absolute errors
            float error_low = sqrt( pow(error_bin_low,2) +
                                    pow(error_mix_low,2) +
                                    //pow(error_LMG_low,2) +
                                    pow(error_eff_low,2) +
                                    pow(error_ref_low,2)
                                    //pow(error_Mu_low,2)
                                    //pow(error_FSR_low,2) +
                                    //pow(error_D1S_low,2) +
                                  );

            float error_high = sqrt( pow(error_bin_high,2) +
                                     pow(error_mix_high,2) +
                                    // pow(error_LMG_high,2) +
                                     pow(error_eff_high,2) +
                                     pow(error_ref_high,2)
                                     //pow(error_Mu_high,2)
                                     //pow(error_FSR_high,2) +
                                     //pow(error_D1S_high,2) +
                                  );


            if (sym_sum){
              error_low  = TMath::Max(error_low,error_high);
              error_high = TMath::Max(error_low,error_high);
            }

            float v2_nominal = g_ob->GetY()[i-1];
            float v2_StatErr = g_ob->GetEYhigh()[i-1];

            g_ob_syst->GetEYlow()[i-1]  = error_low;
            g_ob_syst->GetEYhigh()[i-1] = error_high;

            g_ob_syst->GetEXlow()[i-1]  = 0.08;
            g_ob_syst->GetEXhigh()[i-1] = 0.08;


            hsys_D1S_low->SetBinContent(i, -1*fabs(error_D1S_low));
            hsys_LMG_low->SetBinContent(i, -1*fabs(error_LMG_low));
            hsys_Mu_low->SetBinContent(i, -1*fabs(error_Mu_low));
            hsys_FSR_low->SetBinContent(i, -1*fabs(error_FSR_low));
            hsys_ref_low->SetBinContent(i, -1*fabs(error_ref_low));
            hsys_eff_low->SetBinContent(i, -1*fabs(error_eff_low));
            hsys_bin_low->SetBinContent(i, -1*fabs(error_bin_low));
            hsys_mix_low->SetBinContent(i, -1*fabs(error_mix_low) );


            hsys_D1S_high->SetBinContent(i, fabs(error_D1S_high ));
            hsys_LMG_high->SetBinContent(i, fabs(error_LMG_high ));
            hsys_Mu_high->SetBinContent(i, fabs(error_Mu_high ));
            hsys_FSR_high->SetBinContent(i, fabs(error_FSR_high ));
            hsys_ref_high->SetBinContent(i, fabs(error_ref_high ));
            hsys_eff_high->SetBinContent(i, fabs(error_eff_high ));
            hsys_mix_high->SetBinContent(i, fabs(error_mix_high ));
            hsys_bin_high->SetBinContent(i, fabs(error_bin_high ));


            hsys_sum_low->SetBinContent(i, -1*fabs(error_low ) );
            hsys_sum_high->SetBinContent(i,fabs(error_high   ) );

        }


        TCanvas* c1 = new TCanvas("c1","New Canvas",50,50,600,600);
        h_frame_err_Fn->Draw("AXIS");
        h_frame_err_Fn->GetYaxis()->SetRangeUser(-0.005,0.02);

        hsys_sum_low->SetLineWidth(8);

        hsys_sum_low->Draw("HIST ][ SAME");
        hsys_bin_low->Draw("HIST ][ SAME");
        hsys_eff_low->Draw("HIST ][ SAME");
        //hsys_ref_low->Draw("HIST ][ SAME");
        //hsys_Mu_low->Draw("PHIST ][  SAME");
        //  hsys_LMG_low->Draw("HIST ][ SAME");
        hsys_mix_low->Draw("HIST ][ SAME");
        //hsys_FSR_low->Draw("HIST ][ SAME");
        //hsys_D1S_low->Draw("HIST ][ SAME");

        hsys_sum_high->SetLineWidth(8);
        hsys_sum_high->Draw("HIST ][ SAME");
        hsys_bin_high->Draw("HIST ][ SAME");
        hsys_eff_high->Draw("HIST ][ SAME");
        //hsys_ref_high->Draw("HIST ][ SAME");
        //hsys_Mu_high->Draw("PHIST ][  SAME");
        //  hsys_LMG_high->Draw("HISTSAME");
        hsys_mix_high->Draw("HIST ][ SAME");
        //hsys_FSR_high->Draw("HISTSAME");
        //hsys_D1S_high->Draw("HISTSAME");


        myText(          0.6,0.88+0.02, 1,"#font[72]{ATLAS} Internal",0.04);
        myText(          0.6,0.84+0.02, 1,sysTag.c_str(),0.04);
        myText(          0.6,0.80+0.02, 1, "raw #it{F}_{2} systematics",0.04);
        myMarkerLineText(0.25,0.88+0.02,0.0, 1, 1, 1, 1, "Total Uncertainty", 0.04, true);
        //myMarkerLineText(0.25,0.84+0.02,1.5, kRed   , 22, kRed, 1,"#it{#mu}", 0.04, true);
        //myMarkerLineText(0.25,0.84+0.02,0.0, 4, 2, 4, 2, "LM Reference sym.", 0.04);
        myMarkerLineText(0.25,0.80+0.02,0.0, kMagenta, 1, kMagenta, 1,"#Delta#phi binning sym.", 0.04, true);
        myMarkerLineText(0.25,0.76+0.02,0.0, kSpring+4, 1, kSpring+4, 2,"tracking efficiency", 0.04);
      //  myMarkerLineText(0.25,0.72+0.02,0.0, kOrange-2, 1,kOrange-2, 1,"#it{F}_{2}^{LM} input",  0.04, true);
      myMarkerLineText(0.25,0.72+0.02,0.0, kRed     , 1,kRed      , 2,"Mixed Event",  0.04, true);



        c1->SaveAs(Form("%s/systBreak_%s.pdf",figPath.c_str(),obs.c_str()));


        g_ob     ->Write();
        g_ob_syst->Write();

    }







        // v22Sub_Fn
        {
            string obs = "v22Sub_Fn";

            TH1F* h_ob = (TH1F*) fin->Get(Form("h_%s",obs.c_str()));
            TGraphAsymmErrors* g_ob_stat =  new TGraphAsymmErrors(h_ob);

            TGraphAsymmErrors* g_ob      = (TGraphAsymmErrors*) g_ob_stat->Clone( Form("g_%s_StatError",obs.c_str()) );
            TGraphAsymmErrors* g_ob_syst = (TGraphAsymmErrors*) g_ob->Clone( Form("g_%s_SystError",obs.c_str()) );

            // asymmetric errors
            hsys_Mu_low  = (TH1F*)fin_Mu->Get(Form("h_%s_sysLow", obs.c_str()) );
            hsys_Mu_high = (TH1F*)fin_Mu->Get(Form("h_%s_sysHigh",obs.c_str()) );
            hsys_Mu_low ->SetMarkerColor(2);
            hsys_Mu_low ->SetMarkerStyle(22);
            hsys_Mu_high ->SetMarkerColor(2);
            hsys_Mu_high ->SetMarkerStyle(22);

            hsys_ref_low  = (TH1F*)fin_ref->Get(Form("h_%s_sysLow", obs.c_str()) );
            hsys_ref_high = (TH1F*)fin_ref->Get(Form("h_%s_sysHigh",obs.c_str()) );
            hsys_ref_low ->SetLineColor(4);
            hsys_ref_low ->SetLineStyle(2);
            hsys_ref_high ->SetLineColor(4);
            hsys_ref_high ->SetLineStyle(2);

            hsys_D1S_low  = (TH1F*)fin_D1S->Get(Form("h_%s_sysLow", obs.c_str()) );
            hsys_D1S_high = (TH1F*)fin_D1S->Get(Form("h_%s_sysHigh",obs.c_str()) );
            hsys_D1S_low ->SetLineColor(kMagenta+1);
            hsys_D1S_low ->SetLineStyle(5);
            hsys_D1S_high ->SetLineColor(kMagenta+1);
            hsys_D1S_high ->SetLineStyle(5);

            hsys_FSR_low  = (TH1F*)fin_FSR->Get( Form("h_%s_sysLow", obs.c_str()));
            hsys_FSR_high = (TH1F*)fin_FSR->Get( Form("h_%s_sysHigh",obs.c_str()));
            hsys_FSR_low ->SetLineColor(kSpring+4);
            hsys_FSR_low ->SetLineStyle(3);
            hsys_FSR_high->SetLineColor(kSpring+4);
            hsys_FSR_high->SetLineStyle(3);

            hsys_mix_low   = (TH1F*)fin_mix->Get( Form("h_%s_sysLow", obs.c_str()));
            hsys_mix_high  = (TH1F*)fin_mix->Get( Form("h_%s_sysHigh",obs.c_str()));
            hsys_mix_low->SetLineColor(kRed);
            hsys_mix_low->SetLineStyle(2);
            hsys_mix_high->SetLineColor(kRed);
            hsys_mix_high->SetLineStyle(2);

            hsys_eff_low   = (TH1F*)fin_eff->Get( Form("h_%s_sysLow", obs.c_str()));
            hsys_eff_high  = (TH1F*)fin_eff->Get( Form("h_%s_sysHigh",obs.c_str()));
            hsys_eff_low->SetLineColor(kSpring+4);
            hsys_eff_low->SetLineStyle(2);
            hsys_eff_high->SetLineColor(kSpring+4);
            hsys_eff_high->SetLineStyle(2);

            // symmetric errors
            hsys_LMG_low  = (TH1F*)fin_LMG->Get( Form("h_%s_sysSym", obs.c_str()));
            hsys_LMG_low ->SetLineColor(kOrange-2);
            hsys_LMG_low ->SetLineStyle(1);
            hsys_LMG_high = (TH1F*)hsys_LMG_low->Clone( "hsys_LMG_high");


            hsys_bin_low  = (TH1F*)fin_bin->Get( Form("h_%s_sysSym", obs.c_str()));
            hsys_bin_low->SetLineColor(kMagenta);
            hsys_bin_low->SetLineStyle(1);
            hsys_bin_high = (TH1F*)hsys_bin_low->Clone( "hsys_bin_high");

            /// combined systematics
            hsys_sum_low = (TH1F*)hsys_Mu_low->Clone( Form("hsys_sum_%s_low", obs.c_str()) );
            hsys_sum_low->Reset();
            hsys_sum_low->SetLineColor(1);
            hsys_sum_low->SetLineStyle(1);
            hsys_sum_low->SetLineWidth(2);

            hsys_sum_high = (TH1F*)hsys_Mu_high->Clone( Form("hsys_sum_%s_high", obs.c_str()) );
            hsys_sum_high->Reset();
            hsys_sum_high->SetLineColor(1);
            hsys_sum_high->SetLineStyle(1);
            hsys_sum_high->SetLineWidth(2);


            for (int i=1; i< hsys_sum_low->GetNbinsX()+1; i++) {
                float error_D1S_low = hsys_D1S_low->GetBinContent(i);
                float error_D1S_high = hsys_D1S_high->GetBinContent(i);

                float error_LMG_low = hsys_LMG_low->GetBinContent(i);
                float error_LMG_high = hsys_LMG_high->GetBinContent(i);


                float error_Mu_low = hsys_Mu_low->GetBinContent(i);
                float error_Mu_high = hsys_Mu_high->GetBinContent(i);
                if (sym_Mu){
                  error_Mu_low  = TMath::Max(error_Mu_low,error_Mu_high);
                  error_Mu_high = TMath::Max(error_Mu_low,error_Mu_high);
                }

                float error_FSR_low = hsys_FSR_low->GetBinContent(i);
                float error_FSR_high = hsys_FSR_high->GetBinContent(i);

                float error_ref_low = hsys_ref_low->GetBinContent(i);
                float error_ref_high = hsys_ref_high->GetBinContent(i);
                if (sym_ref){
                  error_ref_low  = TMath::Max(error_ref_low,error_ref_high);
                  error_ref_high = TMath::Max(error_ref_low,error_ref_high);
                }

                float error_eff_low = hsys_eff_low->GetBinContent(i);
                float error_eff_high = hsys_eff_high->GetBinContent(i);
                if (sym_eff){
                  error_eff_low  = TMath::Max(error_eff_low,error_eff_high);
                  error_eff_high = TMath::Max(error_eff_low,error_eff_high);
                }

                float error_mix_low = hsys_mix_low->GetBinContent(i);
                float error_mix_high = hsys_mix_high->GetBinContent(i);
                if (sym_mix){
                  error_mix_low  = TMath::Max(error_mix_low,error_mix_high);
                  error_mix_high = TMath::Max(error_mix_low,error_mix_high);
                }

                float error_bin_low = hsys_bin_low->GetBinContent(i);
                float error_bin_high = hsys_bin_high->GetBinContent(i);
                if (sym_bin){
                  error_bin_low  = TMath::Max(error_bin_low,error_bin_high);
                  error_bin_high = TMath::Max(error_bin_low,error_bin_high);
                }


                // absolute errors
                float error_low = sqrt( pow(error_bin_low,2) +
                                        pow(error_mix_low,2) +
                                        //pow(error_LMG_low,2) +
                                        pow(error_eff_low,2) +
                                        pow(error_ref_low,2)
                                        //pow(error_Mu_low,2)
                                        //pow(error_FSR_low,2) +
                                        //pow(error_D1S_low,2) +
                                      );

                float error_high = sqrt( pow(error_bin_high,2) +
                                        pow(error_mix_high,2) +
                                        // pow(error_LMG_high,2) +
                                         pow(error_eff_high,2) +
                                         pow(error_ref_high,2)
                                         //pow(error_Mu_high,2)
                                         //pow(error_FSR_high,2) +
                                         //pow(error_D1S_high,2) +
                                      );


                if (sym_sum){
                  error_low  = TMath::Max(error_low,error_high);
                  error_high = TMath::Max(error_low,error_high);
                }

                float v2_nominal = g_ob->GetY()[i-1];
                float v2_StatErr = g_ob->GetEYhigh()[i-1];

                g_ob_syst->GetEYlow()[i-1]  = error_low;
                g_ob_syst->GetEYhigh()[i-1] = error_high;

                g_ob_syst->GetEXlow()[i-1]  = 0.08;
                g_ob_syst->GetEXhigh()[i-1] = 0.08;


                hsys_D1S_low->SetBinContent(i, -1*fabs(error_D1S_low));
                hsys_LMG_low->SetBinContent(i, -1*fabs(error_LMG_low));
                hsys_Mu_low->SetBinContent(i, -1*fabs(error_Mu_low));
                hsys_FSR_low->SetBinContent(i, -1*fabs(error_FSR_low));
                hsys_ref_low->SetBinContent(i, -1*fabs(error_ref_low));
                hsys_eff_low->SetBinContent(i, -1*fabs(error_eff_low));
                hsys_bin_low->SetBinContent(i, -1*fabs(error_bin_low));
                hsys_mix_low->SetBinContent(i, -1*fabs(error_mix_low) );


                hsys_D1S_high->SetBinContent(i, fabs(error_D1S_high ));
                hsys_LMG_high->SetBinContent(i, fabs(error_LMG_high ));
                hsys_Mu_high->SetBinContent(i, fabs(error_Mu_high ));
                hsys_FSR_high->SetBinContent(i, fabs(error_FSR_high ));
                hsys_ref_high->SetBinContent(i, fabs(error_ref_high ));
                hsys_eff_high->SetBinContent(i, fabs(error_eff_high ));
                hsys_mix_high->SetBinContent(i, fabs(error_mix_high ));
                hsys_bin_high->SetBinContent(i, fabs(error_bin_high ));


                hsys_sum_low->SetBinContent(i, -1*fabs(error_low ) );
                hsys_sum_high->SetBinContent(i,fabs(error_high   ) );

            }


            TCanvas* c1 = new TCanvas("c1","New Canvas",50,50,600,600);
            h_frame_err_Fn->Draw("AXIS");
            h_frame_err_Fn->GetYaxis()->SetRangeUser(-0.05,0.05);
            h_frame_err_Fn->GetXaxis()->SetRangeUser(60,160);

            hsys_sum_low->SetLineWidth(8);

            hsys_sum_low->Draw("HIST ][ SAME");
            hsys_bin_low->Draw("HIST ][ SAME");
            hsys_eff_low->Draw("HIST ][ SAME");
            hsys_ref_low->Draw("HIST ][ SAME");
            //hsys_Mu_low->Draw("PHIST ][  SAME");
            //  hsys_LMG_low->Draw("HIST ][ SAME");
            hsys_mix_low->Draw("HIST ][ SAME");
            //hsys_FSR_low->Draw("HIST ][ SAME");
            //hsys_D1S_low->Draw("HIST ][ SAME");

            hsys_sum_high->SetLineWidth(8);
            hsys_sum_high->Draw("HIST ][ SAME");
            hsys_bin_high->Draw("HIST ][ SAME");
            hsys_eff_high->Draw("HIST ][ SAME");
            hsys_ref_high->Draw("HIST ][ SAME");
          //  hsys_Mu_high->Draw("PHIST ][  SAME");
            //  hsys_LMG_high->Draw("HISTSAME");
            hsys_mix_high->Draw("HIST ][ SAME");
            //hsys_FSR_high->Draw("HISTSAME");
            //hsys_D1S_high->Draw("HISTSAME");


            myText(          0.6,0.88+0.02, 1,"#font[72]{ATLAS} Internal",0.04);
            myText(          0.6,0.84+0.02, 1,sysTag.c_str(),0.04);
            myText(          0.6,0.80+0.02, 1, "temp. fit #it{F}_{2} syst.",0.04);
            myMarkerLineText(0.25,0.88+0.02,0.0, 1, 1, 1, 1, "Total Uncertainty", 0.04, true);
          //  myMarkerLineText(0.25,0.84+0.02,1.5, kRed   , 22, kRed, 1,"#it{#mu}", 0.04, true);
            myMarkerLineText(0.25,0.80+0.02,0.0, kMagenta, 1, kMagenta, 1,"#Delta#phi binning sym.", 0.04, true);
            myMarkerLineText(0.25,0.76+0.02,0.0, kSpring+4, 1, kSpring+4, 2,"tracking efficiency", 0.04);
            myMarkerLineText(0.25,0.72+0.02,0.0, 4, 2, 4, 2, "LM Reference sym.", 0.04);
          //  myMarkerLineText(0.25,0.72+0.02,0.0, kOrange-2, 1,kOrange-2, 1,"#it{F}_{2}^{LM} input",  0.04, true);



            c1->SaveAs(Form("%s/systBreak_%s.pdf",figPath.c_str(),obs.c_str()));


            g_ob     ->Write();
            g_ob_syst->Write();

        }








      // v22SubD1_Fn
      {
          string obs = "v22SubD1_Fn";

          TH1F* h_ob = (TH1F*) fin->Get(Form("h_%s",obs.c_str()));
          TGraphAsymmErrors* g_ob_stat =  new TGraphAsymmErrors(h_ob);

          TGraphAsymmErrors* g_ob      = (TGraphAsymmErrors*) g_ob_stat->Clone( Form("g_%s_StatError",obs.c_str()) );
          TGraphAsymmErrors* g_ob_syst = (TGraphAsymmErrors*) g_ob->Clone( Form("g_%s_SystError",obs.c_str()) );

          // asymmetric errors
          hsys_Mu_low  = (TH1F*)fin_Mu->Get(Form("h_%s_sysLow", obs.c_str()) );
          hsys_Mu_high = (TH1F*)fin_Mu->Get(Form("h_%s_sysHigh",obs.c_str()) );
          hsys_Mu_low ->SetMarkerColor(2);
          hsys_Mu_low ->SetMarkerStyle(22);
          hsys_Mu_high ->SetMarkerColor(2);
          hsys_Mu_high ->SetMarkerStyle(22);

          hsys_ref_low  = (TH1F*)fin_ref->Get(Form("h_%s_sysLow", obs.c_str()) );
          hsys_ref_high = (TH1F*)fin_ref->Get(Form("h_%s_sysHigh",obs.c_str()) );
          hsys_ref_low ->SetLineColor(4);
          hsys_ref_low ->SetLineStyle(2);
          hsys_ref_high ->SetLineColor(4);
          hsys_ref_high ->SetLineStyle(2);

          hsys_D1S_low  = (TH1F*)fin_D1S->Get(Form("h_%s_sysLow", obs.c_str()) );
          hsys_D1S_high = (TH1F*)fin_D1S->Get(Form("h_%s_sysHigh",obs.c_str()) );
          hsys_D1S_low ->SetLineColor(kMagenta+1);
          hsys_D1S_low ->SetMarkerColor(kMagenta+1);
          hsys_D1S_low ->SetMarkerStyle(20);
          hsys_D1S_low ->SetLineStyle(5);
          hsys_D1S_high ->SetLineColor(kMagenta+1);
          hsys_D1S_high ->SetMarkerColor(kMagenta+1);
          hsys_D1S_high ->SetMarkerStyle(20);
          hsys_D1S_high ->SetLineStyle(5);

          hsys_FSR_low  = (TH1F*)fin_FSR->Get( Form("h_%s_sysLow", obs.c_str()));
          hsys_FSR_high = (TH1F*)fin_FSR->Get( Form("h_%s_sysHigh",obs.c_str()));
          hsys_FSR_low ->SetLineColor(kSpring+4);
          hsys_FSR_low ->SetLineStyle(3);
          hsys_FSR_high->SetLineColor(kSpring+4);
          hsys_FSR_high->SetLineStyle(3);

          hsys_mix_low   = (TH1F*)fin_mix->Get( Form("h_%s_sysLow", obs.c_str()));
          hsys_mix_high  = (TH1F*)fin_mix->Get( Form("h_%s_sysHigh",obs.c_str()));
          hsys_mix_low->SetLineColor(kRed);
          hsys_mix_low->SetLineStyle(2);
          hsys_mix_high->SetLineColor(kRed);
          hsys_mix_high->SetLineStyle(2);

          hsys_eff_low   = (TH1F*)fin_eff->Get( Form("h_%s_sysLow", obs.c_str()));
          hsys_eff_high  = (TH1F*)fin_eff->Get( Form("h_%s_sysHigh",obs.c_str()));
          hsys_eff_low->SetLineColor(kSpring+4);
          hsys_eff_low->SetLineStyle(2);
          hsys_eff_high->SetLineColor(kSpring+4);
          hsys_eff_high->SetLineStyle(2);

          // symmetric errors
          hsys_LMG_low  = (TH1F*)fin_LMG->Get( Form("h_%s_sysSym", obs.c_str()));
          hsys_LMG_low ->SetLineColor(kOrange-2);
          hsys_LMG_low ->SetLineStyle(1);
          hsys_LMG_high = (TH1F*)hsys_LMG_low->Clone( "hsys_LMG_high");


          hsys_bin_low  = (TH1F*)fin_bin->Get( Form("h_%s_sysSym", obs.c_str()));
          hsys_bin_low->SetLineColor(kMagenta);
          hsys_bin_low->SetLineStyle(1);
          hsys_bin_high = (TH1F*)hsys_bin_low->Clone( "hsys_bin_high");

          /// combined systematics
          hsys_sum_low = (TH1F*)hsys_Mu_low->Clone( Form("hsys_sum_%s_low", obs.c_str()) );
          hsys_sum_low->Reset();
          hsys_sum_low->SetLineColor(1);
          hsys_sum_low->SetLineStyle(1);
          hsys_sum_low->SetLineWidth(2);

          hsys_sum_high = (TH1F*)hsys_Mu_high->Clone( Form("hsys_sum_%s_high", obs.c_str()) );
          hsys_sum_high->Reset();
          hsys_sum_high->SetLineColor(1);
          hsys_sum_high->SetLineStyle(1);
          hsys_sum_high->SetLineWidth(2);


          for (int i=1; i< hsys_sum_low->GetNbinsX()+1; i++) {
              float error_D1S_low = hsys_D1S_low->GetBinContent(i);
              float error_D1S_high = hsys_D1S_high->GetBinContent(i);

              float error_LMG_low = hsys_LMG_low->GetBinContent(i);
              float error_LMG_high = hsys_LMG_high->GetBinContent(i);


              float error_Mu_low = hsys_Mu_low->GetBinContent(i);
              float error_Mu_high = hsys_Mu_high->GetBinContent(i);
              if (sym_Mu){
                error_Mu_low  = TMath::Max(error_Mu_low,error_Mu_high);
                error_Mu_high = TMath::Max(error_Mu_low,error_Mu_high);
              }

              float error_FSR_low = hsys_FSR_low->GetBinContent(i);
              float error_FSR_high = hsys_FSR_high->GetBinContent(i);

              float error_ref_low = hsys_ref_low->GetBinContent(i);
              float error_ref_high = hsys_ref_high->GetBinContent(i);
              if (sym_ref){
                error_ref_low  = TMath::Max(error_ref_low,error_ref_high);
                error_ref_high = TMath::Max(error_ref_low,error_ref_high);
              }

              float error_eff_low = hsys_eff_low->GetBinContent(i);
              float error_eff_high = hsys_eff_high->GetBinContent(i);
              if (sym_eff){
                error_eff_low  = TMath::Max(error_eff_low,error_eff_high);
                error_eff_high = TMath::Max(error_eff_low,error_eff_high);
              }

              float error_mix_low = hsys_mix_low->GetBinContent(i);
              float error_mix_high = hsys_mix_high->GetBinContent(i);
              if (sym_mix){
                error_mix_low  = TMath::Max(error_mix_low,error_mix_high);
                error_mix_high = TMath::Max(error_mix_low,error_mix_high);
              }

              float error_bin_low = hsys_bin_low->GetBinContent(i);
              float error_bin_high = hsys_bin_high->GetBinContent(i);
              if (sym_bin){
                error_bin_low  = TMath::Max(error_bin_low,error_bin_high);
                error_bin_high = TMath::Max(error_bin_low,error_bin_high);
              }


              // absolute errors
              float error_low = sqrt( pow(error_bin_low,2) +
                                      pow(error_mix_low,2) +
                                      //pow(error_LMG_low,2) +
                                      pow(error_eff_low,2) +
                                      //pow(error_ref_low,2) +
                                      //pow(error_Mu_low,2) +
                                      pow(error_FSR_low,2) +
                                      pow(error_D1S_low,2)
                                    );

              float error_high = sqrt( pow(error_bin_high,2) +
                                       pow(error_mix_high,2) +
                                      // pow(error_LMG_high,2) +
                                       pow(error_eff_high,2) +
                                      // pow(error_ref_high,2) +
                                       //pow(error_Mu_high,2) +
                                       pow(error_FSR_high,2) +
                                       pow(error_D1S_high,2)
                                    );


              if (sym_sum){
                error_low  = TMath::Max(error_low,error_high);
                error_high = TMath::Max(error_low,error_high);
              }

              float v2_nominal = g_ob->GetY()[i-1];
              float v2_StatErr = g_ob->GetEYhigh()[i-1];

              g_ob_syst->GetEYlow()[i-1]  = error_low;
              g_ob_syst->GetEYhigh()[i-1] = error_high;

              g_ob_syst->GetEXlow()[i-1]  = 0.08;
              g_ob_syst->GetEXhigh()[i-1] = 0.08;


              hsys_D1S_low->SetBinContent(i, -1*fabs(error_D1S_low));
              hsys_LMG_low->SetBinContent(i, -1*fabs(error_LMG_low));
              hsys_Mu_low->SetBinContent(i, -1*fabs(error_Mu_low));
              hsys_FSR_low->SetBinContent(i, -1*fabs(error_FSR_low));
              hsys_ref_low->SetBinContent(i, -1*fabs(error_ref_low));
              hsys_eff_low->SetBinContent(i, -1*fabs(error_eff_low));
              hsys_bin_low->SetBinContent(i, -1*fabs(error_bin_low));
              hsys_mix_low->SetBinContent(i, -1*fabs(error_mix_low) );


              hsys_D1S_high->SetBinContent(i, fabs(error_D1S_high ));
              hsys_LMG_high->SetBinContent(i, fabs(error_LMG_high ));
              hsys_Mu_high->SetBinContent(i, fabs(error_Mu_high ));
              hsys_FSR_high->SetBinContent(i, fabs(error_FSR_high ));
              hsys_ref_high->SetBinContent(i, fabs(error_ref_high ));
              hsys_eff_high->SetBinContent(i, fabs(error_eff_high ));
              hsys_mix_high->SetBinContent(i, fabs(error_mix_high ));
              hsys_bin_high->SetBinContent(i, fabs(error_bin_high ));


              hsys_sum_low->SetBinContent(i, -1*fabs(error_low ) );
              hsys_sum_high->SetBinContent(i,fabs(error_high   ) );

          }


          TCanvas* c1 = new TCanvas("c1","New Canvas",50,50,600,600);
          h_frame_err_Fn->Draw("AXIS");
          h_frame_err_Fn->GetYaxis()->SetRangeUser(-0.03,0.05);
          h_frame_err_Fn->GetXaxis()->SetRangeUser(0,160);

          hsys_sum_low->SetLineWidth(8);

          hsys_sum_low->Draw("HIST ][ SAME");
          hsys_bin_low->Draw("HIST ][ SAME");
          hsys_eff_low->Draw("HIST ][ SAME");
          //hsys_ref_low->Draw("HIST ][ SAME");
          //hsys_Mu_low->Draw("PHIST ][  SAME");
          //  hsys_LMG_low->Draw("HIST ][ SAME");
          hsys_mix_low->Draw("HIST ][ SAME");
          hsys_FSR_low->Draw("HIST ][ SAME");
          hsys_D1S_low->Draw("PHIST ][ SAME");

          hsys_sum_high->SetLineWidth(8);
          hsys_sum_high->Draw("HIST ][ SAME");
          hsys_bin_high->Draw("HIST ][ SAME");
          hsys_eff_high->Draw("HIST ][ SAME");
          //hsys_ref_high->Draw("HIST ][ SAME");
          //hsys_Mu_high->Draw("PHIST ][  SAME");
          //  hsys_LMG_high->Draw("HISTSAME");
          hsys_mix_high->Draw("HIST ][ SAME");
          hsys_FSR_high->Draw("HISTSAME");
          hsys_D1S_high->Draw("PHISTSAME");


          myText(          0.6,0.88+0.02, 1,"#font[72]{ATLAS} Internal",0.04);
          myText(          0.6,0.84+0.02, 1,sysTag.c_str(),0.04);
          myText(          0.6,0.80+0.02, 1, "#it{d}_{1} sub #it{F}_{2} syst.",0.04);
          myMarkerLineText(0.25,0.88+0.02,0.0, 1, 1, 1, 1, "Total Uncertainty", 0.04, true);
          //myMarkerLineText(0.25,0.84+0.02,1.5, kRed   , 22, kRed, 1,"#it{#mu}", 0.04, true);
          myMarkerLineText(0.25,0.80+0.02,0.0, kMagenta, 1, kMagenta, 1,"#Delta#phi binning sym.", 0.04, true);
          myMarkerLineText(0.25,0.76+0.02,0.0, kSpring+4, 1, kSpring+4, 2,"tracking efficiency", 0.04);
          myMarkerLineText(0.25,0.72+0.02,1.5, kMagenta, 20, kMagenta, 1, "#it{d}_{2}/#it{d}_{1} ", 0.04);
          myMarkerLineText(0.25,0.68+0.02,0 , kSpring+4, 1,kSpring+4, 3,"#it{F}_{2}-#it{F}_{1}",  0.04, true);
          myMarkerLineText(0.25,0.64+0.02,0.0, kRed     , 1,kRed      , 2,"Mixed Event",  0.04, true);



          c1->SaveAs(Form("%s/systBreak_%s.pdf",figPath.c_str(),obs.c_str()));


          g_ob     ->Write();
          g_ob_syst->Write();

      }








    // v22Crt_Fn
    {
        string obs = "v22Crt_Fn";

        TH1F* h_ob = (TH1F*) fin->Get(Form("h_%s",obs.c_str()));
        TGraphAsymmErrors* g_ob_stat =  new TGraphAsymmErrors(h_ob);

        TGraphAsymmErrors* g_ob      = (TGraphAsymmErrors*) g_ob_stat->Clone( Form("g_%s_StatError",obs.c_str()) );
        TGraphAsymmErrors* g_ob_syst = (TGraphAsymmErrors*) g_ob->Clone( Form("g_%s_SystError",obs.c_str()) );

        // asymmetric errors
        hsys_Mu_low  = (TH1F*)fin_Mu->Get(Form("h_%s_sysLow", obs.c_str()) );
        hsys_Mu_high = (TH1F*)fin_Mu->Get(Form("h_%s_sysHigh",obs.c_str()) );
        hsys_Mu_low ->SetMarkerColor(2);
        hsys_Mu_low ->SetMarkerStyle(22);
        hsys_Mu_high ->SetMarkerColor(2);
        hsys_Mu_high ->SetMarkerStyle(22);

        hsys_ref_low  = (TH1F*)fin_ref->Get(Form("h_%s_sysLow", obs.c_str()) );
        hsys_ref_high = (TH1F*)fin_ref->Get(Form("h_%s_sysHigh",obs.c_str()) );
        hsys_ref_low ->SetLineColor(4);
        hsys_ref_low ->SetLineStyle(2);
        hsys_ref_high ->SetLineColor(4);
        hsys_ref_high ->SetLineStyle(2);

        hsys_D1S_low  = (TH1F*)fin_D1S->Get(Form("h_%s_sysLow", obs.c_str()) );
        hsys_D1S_high = (TH1F*)fin_D1S->Get(Form("h_%s_sysHigh",obs.c_str()) );
        hsys_D1S_low ->SetLineColor(kMagenta+1);
        hsys_D1S_low ->SetMarkerColor(kMagenta+1);
        hsys_D1S_low ->SetMarkerStyle(20);
        hsys_D1S_low ->SetLineStyle(5);
        hsys_D1S_high ->SetLineColor(kMagenta+1);
        hsys_D1S_high ->SetMarkerColor(kMagenta+1);
        hsys_D1S_high ->SetMarkerStyle(20);
        hsys_D1S_high ->SetLineStyle(5);

        hsys_FSR_low  = (TH1F*)fin_FSR->Get( Form("h_%s_sysLow", obs.c_str()));
        hsys_FSR_high = (TH1F*)fin_FSR->Get( Form("h_%s_sysHigh",obs.c_str()));
        hsys_FSR_low ->SetLineColor(kSpring+4);
        hsys_FSR_low ->SetLineStyle(3);
        hsys_FSR_high->SetLineColor(kSpring+4);
        hsys_FSR_high->SetLineStyle(3);

        hsys_mix_low   = (TH1F*)fin_mix->Get( Form("h_%s_sysLow", obs.c_str()));
        hsys_mix_high  = (TH1F*)fin_mix->Get( Form("h_%s_sysHigh",obs.c_str()));
        hsys_mix_low->SetLineColor(kRed);
        hsys_mix_low->SetLineStyle(2);
        hsys_mix_high->SetLineColor(kRed);
        hsys_mix_high->SetLineStyle(2);

        hsys_eff_low   = (TH1F*)fin_eff->Get( Form("h_%s_sysLow", obs.c_str()));
        hsys_eff_high  = (TH1F*)fin_eff->Get( Form("h_%s_sysHigh",obs.c_str()));
        hsys_eff_low->SetLineColor(kSpring+4);
        hsys_eff_low->SetLineStyle(2);
        hsys_eff_high->SetLineColor(kSpring+4);
        hsys_eff_high->SetLineStyle(2);

        // symmetric errors
        hsys_LMG_low  = (TH1F*)fin_LMG->Get( Form("h_%s_sysSym", obs.c_str()));
        hsys_LMG_low ->SetLineColor(kOrange-2);
        hsys_LMG_low ->SetLineStyle(1);
        hsys_LMG_high = (TH1F*)hsys_LMG_low->Clone( "hsys_LMG_high");


        hsys_bin_low  = (TH1F*)fin_bin->Get( Form("h_%s_sysSym", obs.c_str()));
        hsys_bin_low->SetLineColor(kMagenta);
        hsys_bin_low->SetLineStyle(1);
        hsys_bin_high = (TH1F*)hsys_bin_low->Clone( "hsys_bin_high");

        /// combined systematics
        hsys_sum_low = (TH1F*)hsys_Mu_low->Clone( Form("hsys_sum_%s_low", obs.c_str()) );
        hsys_sum_low->Reset();
        hsys_sum_low->SetLineColor(1);
        hsys_sum_low->SetLineStyle(1);
        hsys_sum_low->SetLineWidth(2);

        hsys_sum_high = (TH1F*)hsys_Mu_high->Clone( Form("hsys_sum_%s_high", obs.c_str()) );
        hsys_sum_high->Reset();
        hsys_sum_high->SetLineColor(1);
        hsys_sum_high->SetLineStyle(1);
        hsys_sum_high->SetLineWidth(2);


        for (int i=1; i< hsys_sum_low->GetNbinsX()+1; i++) {
            float error_D1S_low = hsys_D1S_low->GetBinContent(i);
            float error_D1S_high = hsys_D1S_high->GetBinContent(i);

            float error_LMG_low = hsys_LMG_low->GetBinContent(i);
            float error_LMG_high = hsys_LMG_high->GetBinContent(i);


            float error_Mu_low = hsys_Mu_low->GetBinContent(i);
            float error_Mu_high = hsys_Mu_high->GetBinContent(i);
            if (sym_Mu){
              error_Mu_low  = TMath::Max(error_Mu_low,error_Mu_high);
              error_Mu_high = TMath::Max(error_Mu_low,error_Mu_high);
            }

            float error_FSR_low = hsys_FSR_low->GetBinContent(i);
            float error_FSR_high = hsys_FSR_high->GetBinContent(i);

            float error_ref_low = hsys_ref_low->GetBinContent(i);
            float error_ref_high = hsys_ref_high->GetBinContent(i);
            if (sym_ref){
              error_ref_low  = TMath::Max(error_ref_low,error_ref_high);
              error_ref_high = TMath::Max(error_ref_low,error_ref_high);
            }

            float error_eff_low = hsys_eff_low->GetBinContent(i);
            float error_eff_high = hsys_eff_high->GetBinContent(i);
            if (sym_eff){
              error_eff_low  = TMath::Max(error_eff_low,error_eff_high);
              error_eff_high = TMath::Max(error_eff_low,error_eff_high);
            }

            float error_mix_low = hsys_mix_low->GetBinContent(i);
            float error_mix_high = hsys_mix_high->GetBinContent(i);
            if (sym_mix){
              error_mix_low  = TMath::Max(error_mix_low,error_mix_high);
              error_mix_high = TMath::Max(error_mix_low,error_mix_high);
            }

            float error_bin_low = hsys_bin_low->GetBinContent(i);
            float error_bin_high = hsys_bin_high->GetBinContent(i);
            if (sym_bin){
              error_bin_low  = TMath::Max(error_bin_low,error_bin_high);
              error_bin_high = TMath::Max(error_bin_low,error_bin_high);
            }


            // absolute errors
            float error_low = sqrt( pow(error_bin_low,2) +
                                    pow(error_mix_low,2) +
                                    //pow(error_LMG_low,2) +
                                    pow(error_eff_low,2) +
                                    pow(error_ref_low,2) +
                                    //pow(error_Mu_low,2)
                                    pow(error_FSR_low,2) +
                                    pow(error_D1S_low,2)
                                  );

            float error_high = sqrt( pow(error_bin_high,2) +
                                     pow(error_mix_high,2) +
                                     //pow(error_LMG_high,2) +
                                     pow(error_eff_high,2) +
                                     pow(error_ref_high,2) +
                                     //pow(error_Mu_high,2)
                                     pow(error_FSR_high,2) +
                                     pow(error_D1S_high,2)
                                  );


            if (sym_sum){
              error_low  = TMath::Max(error_low,error_high);
              error_high = TMath::Max(error_low,error_high);
            }

            float v2_nominal = g_ob->GetY()[i-1];
            float v2_StatErr = g_ob->GetEYhigh()[i-1];

            g_ob_syst->GetEYlow()[i-1]  = error_low;
            g_ob_syst->GetEYhigh()[i-1] = error_high;

            g_ob_syst->GetEXlow()[i-1]  = 0.08;
            g_ob_syst->GetEXhigh()[i-1] = 0.08;


            hsys_D1S_low->SetBinContent(i, -1*fabs(error_D1S_low));
            hsys_LMG_low->SetBinContent(i, -1*fabs(error_LMG_low));
            hsys_Mu_low->SetBinContent(i, -1*fabs(error_Mu_low));
            hsys_FSR_low->SetBinContent(i, -1*fabs(error_FSR_low));
            hsys_ref_low->SetBinContent(i, -1*fabs(error_ref_low));
            hsys_eff_low->SetBinContent(i, -1*fabs(error_eff_low));
            hsys_bin_low->SetBinContent(i, -1*fabs(error_bin_low));
            hsys_mix_low->SetBinContent(i, -1*fabs(error_mix_low) );


            hsys_D1S_high->SetBinContent(i, fabs(error_D1S_high ));
            hsys_LMG_high->SetBinContent(i, fabs(error_LMG_high ));
            hsys_Mu_high->SetBinContent(i, fabs(error_Mu_high ));
            hsys_FSR_high->SetBinContent(i, fabs(error_FSR_high ));
            hsys_ref_high->SetBinContent(i, fabs(error_ref_high ));
            hsys_eff_high->SetBinContent(i, fabs(error_eff_high ));
            hsys_mix_high->SetBinContent(i, fabs(error_mix_high ));
            hsys_bin_high->SetBinContent(i, fabs(error_bin_high ));


            hsys_sum_low->SetBinContent(i, -1*fabs(error_low ) );
            hsys_sum_high->SetBinContent(i,fabs(error_high   ) );

        }


        TCanvas* c1 = new TCanvas("c1","New Canvas",50,50,600,600);
        h_frame_err_Fn->Draw("AXIS");
        h_frame_err_Fn->GetYaxis()->SetRangeUser(-0.015,0.05);
        h_frame_err_Fn->GetXaxis()->SetRangeUser(60,160);

        hsys_sum_low->SetLineWidth(8);

        hsys_sum_low->Draw("HIST ][ SAME");
        hsys_bin_low->Draw("HIST ][ SAME");
        hsys_eff_low->Draw("HIST ][ SAME");
        hsys_ref_low->Draw("HIST ][ SAME");
        //hsys_Mu_low->Draw("PHIST ][  SAME");
        //hsys_LMG_low->Draw("HIST ][ SAME");
        hsys_mix_low->Draw("HIST ][ SAME");
        hsys_FSR_low->Draw("HIST ][ SAME");
        hsys_D1S_low->Draw("PHIST ][ SAME");

        hsys_sum_high->SetLineWidth(8);
        hsys_sum_high->Draw("HIST ][ SAME");
        hsys_bin_high->Draw("HIST ][ SAME");
        hsys_eff_high->Draw("HIST ][ SAME");
        hsys_ref_high->Draw("HIST ][ SAME");
        //hsys_Mu_high->Draw("PHIST ][  SAME");
        //hsys_LMG_high->Draw("HISTSAME");
        hsys_mix_high->Draw("HIST ][ SAME");
        hsys_FSR_high->Draw("HISTSAME");
        hsys_D1S_high->Draw("PHISTSAME");


        myText(          0.6,0.88+0.02, 1,"#font[72]{ATLAS} Internal",0.04);
        myText(          0.6,0.84+0.02, 1,sysTag.c_str(),0.04);
        myText(          0.6,0.80+0.02, 1, "Crtd. temp. fit #it{F}_{2} syst.",0.04);
        myMarkerLineText(0.25,0.88+0.02,0.0, 1, 1, 1, 1, "Total Uncertainty", 0.04, true);
        myMarkerLineText(0.25,0.84+0.02,1.5, kRed   , 22, kRed, 1,"#it{#mu}", 0.04, true);
        myMarkerLineText(0.25,0.80+0.02,0.0, kMagenta, 1, kMagenta, 1,"#Delta#phi binning sym.", 0.04, true);
        myMarkerLineText(0.25,0.76+0.02,0.0, kSpring+4, 1, kSpring+4, 2,"tracking efficiency", 0.04);
        //myMarkerLineText(0.25,0.72+0.02,0.0, kOrange-2, 1,kOrange-2, 1,"#it{F}_{2}^{LM} input",  0.04, true);
        myMarkerLineText(0.25,0.72+0.02,1.5, kMagenta, 20, kMagenta, 1, "#it{d}_{2}/#it{d}_{1} ", 0.04);
        myMarkerLineText(0.25,0.68+0.02,0 , kSpring+4, 1,kSpring+4, 3,"#it{F}_{2}-#it{F}_{1}",  0.04, true);
        myMarkerLineText(0.25,0.64+0.02,0.0, kRed     , 1,kRed      , 2,"Mixed Event",  0.04, true);

        myMarkerLineText(0.25,0.60+0.02,0.0, 4, 2, 4, 2, "LM Reference sym.", 0.04);


        c1->SaveAs(Form("%s/systBreak_%s.pdf",figPath.c_str(),obs.c_str()));


        g_ob     ->Write();
        g_ob_syst->Write();

    }






}
