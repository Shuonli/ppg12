#include "/sphenix/u/bseidlitz/plotstyle/AtlasUtils.C"
#include "/sphenix/u/shuhang98/AtlasStyle.C"


void truthisoeff(){

    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    /*
    TFile *finboth10 = new TFile("/sphenix/user/shuhangli/pi0pythiasim/condorboth10/photoniso_both.root", "READ");
    TFile *finboth20 = new TFile("/sphenix/user/shuhangli/pi0pythiasim/condorboth20/photoniso_both.root", "READ");

    TFile *finprompt10 = new TFile("/sphenix/user/shuhangli/pi0pythiasim/condorprompt10/photoniso_prompt.root", "READ");
    TFile *finprompt20 = new TFile("/sphenix/user/shuhangli/pi0pythiasim/condorprompt20/photoniso_prompt.root", "READ");

    TFile *finfrag10 = new TFile("/sphenix/user/shuhangli/pi0pythiasim/condorfrag10/photoniso_frag.root", "READ");
    TFile *finfrag20 = new TFile("/sphenix/user/shuhangli/pi0pythiasim/condorfrag20/photoniso_frag.root", "READ");
    */

    TFile *finboth10 = new TFile("/sphenix/user/shuhangli/ppg12/showershapecheck/MC_showershape_treeana.root", "READ");

    TFile *finprompt10 = new TFile("/sphenix/user/shuhangli/ppg12/showershapecheck/MC_showershape_treeana.root", "READ");

    TFile *finfrag10 = new TFile("/sphenix/user/shuhangli/ppg12/showershapecheck/MC_showershape_treeana.root", "READ");

    //photonpT_isoet02_direct

    std::string r = "4";

    std::string bothname = "photonpT_isoet0" + r + "_direct";
    std::string promptname = "photonpT_isoet0" + r + "_direct";
    std::string fragname = "photonpT_isoet0" + r + "_frag";
    TH2D *hboth = (TH2D*)finboth10->Get(bothname.c_str())->Clone("hphotonboth");


    TH2D *hprompt = (TH2D*)finprompt10->Get(promptname.c_str())->Clone("hphotonprompt");


    TH2D *hfrag = (TH2D*)finfrag10->Get(fragname.c_str())->Clone("hphotonfrag");


    std::vector<std::pair<float, float>> pTrange = {{10,15}, {15,20}, {20,25}, {25,30}};

    std::vector<int> colors = {kPink+8, kSpring-7, kAzure-3, kViolet+3, kOrange +10};
    
    const int nptbins = (int)pTrange.size();

    TH1D *hbothisoeff[nptbins];
    TH1D *hpromptisoeff[nptbins];
    TH1D *hfragisoeff[nptbins];

    TCanvas *cbothisoeff = new TCanvas("cbothisoeff", "", 800, 750);
    cbothisoeff->SetTopMargin(0.1);
    TCanvas *cpromptisoeff = new TCanvas("cpromptisoeff", "", 800, 750);
    cpromptisoeff->SetTopMargin(0.1);
    TCanvas *cfragisoeff = new TCanvas("cfragisoeff", "", 800, 750);
    cfragisoeff->SetTopMargin(0.1);

    TH2D *hbothframe = new TH2D("hbothframe", "", 10, 0, 20, 1, 0.87, 1.02);
    hbothframe->GetXaxis()->SetTitle("Truth iso_{E_{T}} cut [GeV]");
    hbothframe->GetYaxis()->SetTitle("Efficiency");
    cbothisoeff->cd();
    hbothframe->Draw();
    myText(0.2, 0.8, kBlack, "Inclusive Photons", 0.04);
    myText(0.1, 0.95, kBlack, "#bf{#it{sPHENIX}} Internal");
    myText(0.4, 0.95, kBlack, "Pythia photon+jet samples", 0.04);

    TH2D *hpromptframe = new TH2D("hpromptframe", "", 10, 0, 20, 1, 0.93, 1.02);
    hpromptframe->GetXaxis()->SetTitle("Truth iso_{E_{T}} cut [GeV]");
    hpromptframe->GetYaxis()->SetTitle("Efficiency");
    cpromptisoeff->cd();
    hpromptframe->Draw();
    myText(0.2, 0.8, kBlack, "Direct Photons", 0.04);
    myText(0.1, 0.95, kBlack, "#bf{#it{sPHENIX}} Internal");
    myText(0.4, 0.95, kBlack, "Pythia photon+jet samples", 0.04);
    myText(0.6, 0.8, kBlack, Form("Isolation R = 0.%s", r.c_str()), 0.04);


    TH2D *hfragframe = new TH2D("hfragframe", "", 10, 0, 20, 1, 0.75, 1.1);
    hfragframe->GetXaxis()->SetTitle("Truth iso_{E_{T}} cut [GeV]");
    hfragframe->GetYaxis()->SetTitle("Efficiency");
    cfragisoeff->cd();
    hfragframe->Draw();
    myText(0.2, 0.8, kBlack, "Fragmentation Photons", 0.04);
    myText(0.1, 0.95, kBlack, "#bf{#it{sPHENIX}} Internal");
    myText(0.4, 0.95, kBlack, "Pythia photon+jet samples", 0.04);
    myText(0.6, 0.8, kBlack, Form("Isolation R = 0.%s", r.c_str()), 0.04);


    for(int ipt = 0; ipt < nptbins; ipt++){
        float ptlow = pTrange[ipt].first;
        float ptup = pTrange[ipt].second;
        int binlow = hboth->GetXaxis()->FindBin(ptlow);
        int binup = hboth->GetXaxis()->FindBin(ptup);
        TH1D *hbothiso = hboth->ProjectionY(Form("hbothiso_%d", ipt), binlow, binup);
        TH1D *hpromptiso = hprompt->ProjectionY(Form("hpromptiso_%d", ipt), binlow, binup);
        TH1D *hfragiso = hfrag->ProjectionY(Form("hfragiso_%d", ipt), binlow, binup);

        hbothisoeff[ipt] = new TH1D(Form("hbothisoeff_%d", ipt), "", 20, 0.5, 20.5);
        hpromptisoeff[ipt] = new TH1D(Form("hpromptisoeff_%d", ipt), "", 20, 0.5, 20.5);
        hfragisoeff[ipt] = new TH1D(Form("hfragisoeff_%d", ipt), "", 20, 0.5, 20.5);

        int totalboth = hbothiso->GetEntries();
        int totalprompt = hpromptiso->GetEntries();
        int totalfrag = hfragiso->GetEntries();

        //set each bin
        for(int ibin = 1; ibin <= hbothisoeff[ipt]->GetNbinsX(); ibin++){
            float bincenter = hbothisoeff[ipt]->GetBinCenter(ibin);
            float bothiso = hbothiso->Integral(1, hbothiso->FindBin(bincenter));
            float promptiso = hpromptiso->Integral(1, hpromptiso->FindBin(bincenter));
            float fragiso = hfragiso->Integral(1, hfragiso->FindBin(bincenter));
            //binomial error
            float bothisoeff = bothiso/totalboth;
            float promptisoeff = promptiso/totalprompt;
            float fragisoeff = fragiso/totalfrag;

            float bothisoefferr = sqrt(bothisoeff*(1-bothisoeff)/totalboth);
            float promptisoefferr = sqrt(promptisoeff*(1-promptisoeff)/totalprompt);
            float fragisoefferr = sqrt(fragisoeff*(1-fragisoeff)/totalfrag);

            hbothisoeff[ipt]->SetBinContent(ibin, bothisoeff);
            hbothisoeff[ipt]->SetBinError(ibin, bothisoefferr);

            hpromptisoeff[ipt]->SetBinContent(ibin, promptisoeff);
            hpromptisoeff[ipt]->SetBinError(ibin, promptisoefferr);

            hfragisoeff[ipt]->SetBinContent(ibin, fragisoeff);
            hfragisoeff[ipt]->SetBinError(ibin, fragisoefferr);
        }
        float ypos = 0.55 - 0.05*ipt;
        
        cbothisoeff->cd();
        hbothisoeff[ipt]->SetMarkerColor(colors[ipt]);
        hbothisoeff[ipt]->SetLineColor(colors[ipt]);
        hbothisoeff[ipt]->Draw("same");
        myMarkerText(0.6, ypos, colors[ipt],20,Form("%.0f < p_{T} < %.0f GeV", ptlow, ptup), 1, 0.05);

        cpromptisoeff->cd();
        hpromptisoeff[ipt]->SetMarkerColor(colors[ipt]);
        hpromptisoeff[ipt]->SetLineColor(colors[ipt]);
        hpromptisoeff[ipt]->Draw("same");
        myMarkerText(0.6, ypos, colors[ipt],20,  Form("%.0f < p_{T} < %.0f GeV", ptlow, ptup), 1, 0.05);

        cfragisoeff->cd();
        hfragisoeff[ipt]->SetMarkerColor(colors[ipt]);
        hfragisoeff[ipt]->SetLineColor(colors[ipt]);
        hfragisoeff[ipt]->Draw("same");
        myMarkerText(0.6, ypos, colors[ipt], 20, Form("%.0f < p_{T} < %.0f GeV", ptlow, ptup), 1, 0.05);
        



    }


}
