#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>
#include <yaml-cpp/yaml.h>


double myfunc(double *x, double *params)
{
    double NsigA = x[0]; // Independent variable (N_A^sig)

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
    if (denominator == 0)
    {
        denominator = 1e-12;
    }

    double fraction = (NC - cC * NsigA) / denominator;

    // The function to find the root of
    double result = NsigA - (NA - R * numerator * fraction);
    return result;
}

void CalculatePhotonYield(const std::string &configname = "config.yaml")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    std::string var_type = configYaml["output"]["var_type"].as<std::string>();

    std::string outfilename = configYaml["output"]["final_outfile"].as<std::string>() + "_" + var_type + ".root";

    // calculate and unfold the isolated photon spectrum given result from RecoEffCalculator.C

    std::string datainput = configYaml["output"]["data_outfile"].as<std::string>() + "_" + var_type + ".root";

    std::string siminput = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + var_type + ".root";

    std::string unfoldinput = configYaml["output"]["response_outfile"].as<std::string>() + "_" + var_type + ".root";

    std::string histogram_postfix = "_0";

    std::string tight_iso_cluster_name = "h_tight_iso_cluster"; // A

    std::string tight_noniso_cluster_name = "h_tight_noniso_cluster"; // B

    std::string nontight_iso_cluster_name = "h_nontight_iso_cluster"; // C

    std::string nontight_noniso_cluster_name = "h_nontight_noniso_cluster"; // D

    std::string tight_iso_cluster_signal_name = tight_iso_cluster_name + "_signal" + histogram_postfix; // A

    std::string tight_noniso_cluster_signal_name = tight_noniso_cluster_name + "_signal" + histogram_postfix; // B

    std::string nontight_iso_cluster_signal_name = nontight_iso_cluster_name + "_signal" + histogram_postfix; // C

    std::string nontight_noniso_cluster_signal_name = nontight_noniso_cluster_name + "_signal" + histogram_postfix; // D


    std::string reco_efficiency_name = "eff_reco_eta" + histogram_postfix;

    std::string reco_efficiency_iso_name = "eff_iso_eta" + histogram_postfix;

    std::string reco_efficiency_id_name = "eff_id_eta" + histogram_postfix;

    // calculate the signal leakage from simulation

    TFile *fsimin = new TFile(siminput.c_str(), "READ");

    TH1F *h_tight_iso_cluster_signal = (TH1F *)fsimin->Get(tight_iso_cluster_signal_name.c_str());

    TH1F *h_tight_noniso_cluster_signal = (TH1F *)fsimin->Get(tight_noniso_cluster_signal_name.c_str());

    TH1F *h_nontight_iso_cluster_signal = (TH1F *)fsimin->Get(nontight_iso_cluster_signal_name.c_str());

    TH1F *h_nontight_noniso_cluster_signal = (TH1F *)fsimin->Get(nontight_noniso_cluster_signal_name.c_str());

    TH1F *h_leak_B = (TH1F *)h_tight_noniso_cluster_signal->Clone("h_leak_B");
    TH1F *h_leak_C = (TH1F *)h_nontight_iso_cluster_signal->Clone("h_leak_C");
    TH1F *h_leak_D = (TH1F *)h_nontight_noniso_cluster_signal->Clone("h_leak_D");

    TEfficiency *eff_reco = (TEfficiency *)fsimin->Get(reco_efficiency_name.c_str());
    TEfficiency *eff_iso = (TEfficiency *)fsimin->Get(reco_efficiency_iso_name.c_str());
    TEfficiency *eff_id = (TEfficiency *)fsimin->Get(reco_efficiency_id_name.c_str());

    // h_tight_iso_cluster_signal->Sumw2();
    // h_leak_B->Sumw2();
    // h_leak_C->Sumw2();
    // h_leak_D->Sumw2();

    h_leak_B->Divide(h_tight_iso_cluster_signal);
    h_leak_C->Divide(h_tight_iso_cluster_signal);
    h_leak_D->Divide(h_tight_iso_cluster_signal);
    // calculate error manually
    /*
    for (int ibin = 1; ibin <= h_leak_B->GetNbinsX(); ibin++)
    {
        double A = h_tight_iso_cluster_signal->GetBinContent(ibin);
        double B = h_tight_noniso_cluster_signal->GetBinContent(ibin);
        double C = h_nontight_iso_cluster_signal->GetBinContent(ibin);
        double D = h_nontight_noniso_cluster_signal->GetBinContent(ibin);


        double errA = h_tight_iso_cluster_signal->GetBinError(ibin);
        double errB = h_tight_noniso_cluster_signal->GetBinError(ibin);
        double errC = h_nontight_iso_cluster_signal->GetBinError(ibin);
        double errD = h_nontight_noniso_cluster_signal->GetBinError(ibin);

        double errBA = h_leak_B->GetBinContent(ibin) * sqrt(pow(errB / B, 2) + pow(errA / A, 2));
        double errCA = h_leak_C->GetBinContent(ibin) * sqrt(pow(errC / C, 2) + pow(errA / A, 2));
        double errDA = h_leak_D->GetBinContent(ibin) * sqrt(pow(errD / D, 2) + pow(errA / A, 2));
        std::cout << "errBA: " << errBA << " errCA: " << errCA << " errDA: " << errDA << std::endl;
        //nan gard
        if (isnan(errBA))
        {
            errBA = 0;
        }
        if (isnan(errCA))
        {
            errCA = 0;
        }
        if (isnan(errDA))
        {
            errDA = 0;
        }
        h_leak_B->SetBinError(ibin, errBA);
        h_leak_C->SetBinError(ibin, errCA);
        h_leak_D->SetBinError(ibin, errDA);
    }
    */

    // get the raw spectrum from data

    std::string tight_iso_cluster_signal_name_data = tight_iso_cluster_name + histogram_postfix; // A
    // std::cout<<"tight_iso_cluster_signal_name_data: "<<tight_iso_cluster_signal_name_data<<std::endl;

    std::string tight_noniso_cluster_signal_name_data = tight_noniso_cluster_name + histogram_postfix; // B

    std::string nontight_iso_cluster_signal_name_data = nontight_iso_cluster_name + histogram_postfix; // C

    std::string nontight_noniso_cluster_signal_name_data = nontight_noniso_cluster_name + histogram_postfix; // D

    TFile *fdatain = new TFile(datainput.c_str(), "READ");

    TH1F *h_tight_iso_cluster_signal_data = (TH1F *)fdatain->Get(tight_iso_cluster_signal_name_data.c_str());

    TH1F *h_tight_noniso_cluster_signal_data = (TH1F *)fdatain->Get(tight_noniso_cluster_signal_name_data.c_str());

    TH1F *h_nontight_iso_cluster_signal_data = (TH1F *)fdatain->Get(nontight_iso_cluster_signal_name_data.c_str());

    TH1F *h_nontight_noniso_cluster_signal_data = (TH1F *)fdatain->Get(nontight_noniso_cluster_signal_name_data.c_str());

    // a histogram with same binning as h_tight_iso_cluster_signal_data
    //TH1F *h_purity = (TH1F *)h_tight_iso_cluster_signal_data->Clone("h_purity");
    //h_purity->Reset();
    TH1F *h_data_sub = (TH1F *)h_tight_iso_cluster_signal_data->Clone("h_data_sub");
    h_data_sub->Reset();

    // purity histogram with leakage correction
    //TH1F *h_purity_leak = (TH1F *)h_tight_iso_cluster_signal_data->Clone("h_purity_leak");
    //h_purity_leak->Reset();
    TH1F *h_data_sub_leak = (TH1F *)h_tight_iso_cluster_signal_data->Clone("h_data_sub_leak");
    h_data_sub_leak->Reset();

    // calculate the purity
    TRandom3 *rand = new TRandom3(0);
    for (int ibin = 1; ibin <= h_tight_iso_cluster_signal_data->GetNbinsX(); ibin++)
    {
        std::cout << "ibin: " << ibin << "bin center: " << h_tight_iso_cluster_signal_data->GetBinCenter(ibin) << std::endl;
        int nsamples = 10000;

        double A = h_tight_iso_cluster_signal_data->GetBinContent(ibin);
        double B = h_tight_noniso_cluster_signal_data->GetBinContent(ibin);
        double C = h_nontight_iso_cluster_signal_data->GetBinContent(ibin);
        double D = h_nontight_noniso_cluster_signal_data->GetBinContent(ibin);

        double errA = h_tight_iso_cluster_signal_data->GetBinError(ibin);
        double errB = h_tight_noniso_cluster_signal_data->GetBinError(ibin);
        double errC = h_nontight_iso_cluster_signal_data->GetBinError(ibin);
        double errD = h_nontight_noniso_cluster_signal_data->GetBinError(ibin);

        double cB = h_leak_B->GetBinContent(ibin);
        double cC = h_leak_C->GetBinContent(ibin);
        double cD = h_leak_D->GetBinContent(ibin);

        double errcb = h_leak_B->GetBinError(ibin);
        double errcc = h_leak_C->GetBinError(ibin);
        double errcd = h_leak_D->GetBinError(ibin);

        TH1D *h_result = new TH1D("h_result", "h_result", 100, 0, A*10);

        TH1D *h_result_leak = new TH1D("h_result_leak", "h_result_leak", 100, 0, A*10);

        for (int isample = 0; isample < nsamples; isample++)
        {
            
            double NA = rand->Gaus(A, errA);
            double NB = rand->Gaus(B, errB);
            double NC = rand->Gaus(C, errC);
            double ND = rand->Gaus(D, errD);
            
            /*
            double NA = rand->Gaus(A, 0);
            double NB = rand->Gaus(B, 0);
            double NC = rand->Gaus(C, 0);
            double ND = rand->Gaus(D, 0);
            */
            double CB = 0;
            double CC = 0;
            double CD = 0;

            TF1 *f = new TF1("myfunc", myfunc, 0.0, A, 8);
            f->SetParameter(0, NA);
            f->SetParameter(1, NB);
            f->SetParameter(2, NC);
            f->SetParameter(3, ND);
            f->SetParameter(4, CB);
            f->SetParameter(5, CC);
            f->SetParameter(6, CD);
            f->SetParameter(7, 1);

            double root = f->GetX(0, 0, NA);

            double root_manual = NA - NB * NC / ND;

            //fill if the root is not nan nor inf
            if (root == root && abs(root) != std::numeric_limits<double>::infinity())
            {
                h_result->Fill(root);
            }

            

            //std::cout<< "NA: " << NA << " NB: " << NB << " NC: " << NC << " ND: " << ND << std::endl;
            //std::cout << "root: " << root << " NA: " << NA << std::endl;
            //std::cout << "root_manual: " << root_manual << " NA: " << NA << std::endl;

            CB = rand->Gaus(cB, errcb);
            CC = rand->Gaus(cC, errcc);
            CD = rand->Gaus(cD, errcd);

            f->SetParameter(4, CB);
            f->SetParameter(5, CC);
            f->SetParameter(6, CD);

            double root_leak = f->GetX(0, 0, NA);

            if (root_leak == root_leak && abs(root_leak) != std::numeric_limits<double>::infinity())
            {
                h_result_leak->Fill(root_leak);
            }


            //std::cout << "root_leak: " << root_leak << " NA: " << NA << std::endl;
        }
        std::cout<<"NA: "<<A<<" NB: "<<B<<" NC: "<<C<<" ND: "<<D<<std::endl;
        std::cout<<"result mean: "<<h_result->GetMean()<<" result std: "<<h_result->GetRMS()<<std::endl;
        std::cout<<"result leak mean: "<<h_result_leak->GetMean()<<" result leak std: "<<h_result_leak->GetRMS()<<std::endl;


        h_data_sub->SetBinContent(ibin, h_result->GetMean());
        h_data_sub->SetBinError(ibin, h_result->GetRMS());

        h_data_sub_leak->SetBinContent(ibin, h_result_leak->GetMean());
        h_data_sub_leak->SetBinError(ibin, h_result_leak->GetRMS());


        
    }
    // purity correction with and without leakage correction
    //TGraph for the purity
    TGraphAsymmErrors *gpurity = new TGraphAsymmErrors(h_data_sub, h_tight_iso_cluster_signal_data);
    gpurity->SetName("gpurity");
    gpurity->SetTitle("Purity");
    gpurity->GetXaxis()->SetTitle("E_{T} [GeV]");

    TGraphAsymmErrors *gpurity_leak = new TGraphAsymmErrors(h_data_sub_leak, h_tight_iso_cluster_signal_data);
    gpurity_leak->SetName("gpurity_leak");
    gpurity_leak->SetTitle("Purity with leakage correction");
    gpurity_leak->GetXaxis()->SetTitle("E_{T} [GeV]");

 
    //place holder for reweighting before unfolding
    TFile *funfold = new TFile(unfoldinput.c_str(), "READ");

    std::string responsematrixname = "response_matrix_full" + histogram_postfix;
    RooUnfoldResponse *response = (RooUnfoldResponse *)funfold->Get(responsematrixname.c_str());

    int niterations_total = 10;
    int resultit = 2;
    int resultleak = 1; //0 for no leakage correction
    std::vector<TH1D *> h_unfold_sub_list;
    std::vector<TH1D *> h_unfold_sub_leak_list;

    for (int i = 0; i < niterations_total; i++)
    {
        RooUnfoldBayes *unfold = new RooUnfoldBayes(response, h_data_sub, i + 1);
        TH1D *h_unfold_sub = (TH1D *)unfold->Hunfold()->Clone(Form("h_unfold_sub_%d", i+1));
        h_unfold_sub->SetName(Form("h_unfold_sub_%d", i+1));
        h_unfold_sub_list.push_back(h_unfold_sub);



        RooUnfoldBayes *unfold_leak = new RooUnfoldBayes(response, h_data_sub_leak, i + 1);
        TH1D *h_unfold_sub_leak = (TH1D *)unfold_leak->Hunfold()->Clone(Form("h_unfold_sub_leak_%d", i+1));
        h_unfold_sub_leak->SetName(Form("h_unfold_sub_leak_%d", i+1));
        h_unfold_sub_leak_list.push_back(h_unfold_sub_leak);

    }

    //efficiency correction
    for (int i = 0; i < niterations_total; i++)
    {
        TH1D *h_unfold_sub = h_unfold_sub_list[i];
        TH1D *h_unfold_sub_leak = h_unfold_sub_leak_list[i];

        

        for (int ibin = 1; ibin <= h_unfold_sub->GetNbinsX(); ibin++)
        {
            float eff_reco_val = eff_reco->GetEfficiency(ibin);
            float eff_iso_val = eff_iso->GetEfficiency(ibin);
            float eff_id_val = eff_id->GetEfficiency(ibin);

            float total_eff = eff_reco_val * eff_iso_val * eff_id_val;

            h_unfold_sub->SetBinContent(ibin, h_unfold_sub->GetBinContent(ibin) / total_eff);
            h_unfold_sub->SetBinError(ibin, h_unfold_sub->GetBinError(ibin) / total_eff);

            h_unfold_sub_leak->SetBinContent(ibin, h_unfold_sub_leak->GetBinContent(ibin) / total_eff);
            h_unfold_sub_leak->SetBinError(ibin, h_unfold_sub_leak->GetBinError(ibin) / total_eff);
            
        }
    }

    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");

    for (int i = 0; i < niterations_total; i++)
    {
        h_unfold_sub_list[i]->Write();
        h_unfold_sub_leak_list[i]->Write();

        if(i+1 == resultit)
        {
            if(resultleak == 0)
            {
                TH1F *h_unfold_sub_result = (TH1F *)h_unfold_sub_list[i]->Clone("h_unfold_sub_result");
                h_unfold_sub_result->Write("h_unfold_sub_result");
            }
            else
            {
                TH1F *h_unfold_sub_result = (TH1F *)h_unfold_sub_leak_list[i]->Clone("h_unfold_sub_result");
                h_unfold_sub_result->Write("h_unfold_sub_result");
            }
        }
    }
    gpurity->Write();
    gpurity_leak->Write();
    h_leak_B->Write();
    h_leak_C->Write();
    h_leak_D->Write();

    h_data_sub->Write();
    h_data_sub_leak->Write();

    fout->Write();
    fout->Close();
 

    // unfold the spectrum based on data's response matrix
}