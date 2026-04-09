#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TBox.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TPolyLine.h"
#include "TLegend.h"

// author:  J.Nagle - email jamie.nagle@colorado.edu for more information

#include <iostream>
#include <fstream>

void nlo_plot_count (bool AuAucase = true, bool ptgtoption = true, bool plotheavy = false, bool onlyheavy = false) {

  // read ascii files for pp@200 GeV NLO pQCD calculations from Werner V.

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1",10,10,700,900);
  c1->SetGrid();
  c1->SetLogy(1);

  // template histogram for overlays  
  TH1F *h1 = new TH1F("h1","h1",100,0.0,60.0);
  if (AuAucase) {
    // if we consider 50 billion AuAu MB, or 10 billion of these 20% central
    h1->SetMaximum(1.e+8);
    h1->SetMinimum(1.0/10.e1); // that sets the bottom scale
    h1->SetYTitle("Counts/Event with p_{T} > p_{T}(cut) [AuAu L=7/nb]");
    if (!ptgtoption)     h1->SetYTitle("Counts");
  } else {
    h1->SetMaximum(1.e+8);
    h1->SetMinimum(1.0/10.e1);
    h1->SetYTitle("Counts/Event with p_{T} > p_{T}(cut) [pp L=100/pb]");
    if (!ptgtoption)    h1->SetYTitle("Counts/Event [p-p @ 200 GeV]");
  }
  h1->GetYaxis()->SetTitleOffset(1.5);
  h1->SetXTitle("Transverse Momentum (GeV/c)");
  h1->Draw();

  // Werner's numbers are in pb/GeV^2 which gets multiplied up x28 
  // to correspond to 24 pb^-1 or 1e12 pp@200 GeV events
  // Also need to multiply by 2*pi*pt*dpt*deta to get counts per bin with dpt = 2.5 and deta = 2.0
  // To get counts per event scale down by 1 / 1e12

  double ppscalefactor = 2.0*3.14159 * 2.0 * 100 *1 * 2.5;;  // need to include pt factor point-by-point

  // if we have 50 billion AuAu minimum bias events, then 10 billion 20% central AuAu events
  // the Ncoll value = 800
  // thus the yield is 800 x higher per event, but with pizero suppression by 0.2

  //assume 7 nb^-1 for AuAu
  double AuAucentralscalefactor = 2.0*3.14159 * 2.0 * 7 *1E-3 * 197 * 197 * 2.5;

  double scalefactor = ppscalefactor;
  double pizerosupp  = 1.0;
  if (AuAucase) {
    scalefactor = AuAucentralscalefactor;
    pizerosupp = 0.2;
  }

  // for heavy flavor - the numbers are just dsigma/dpt (in 1 GeV bins) in units pb/GeV
  // thus multipy by x28 to get 1e12 pp events
  // then muliply by dpt = 1 GeV (so just one)
  double heavyscalefactor = 28.0 / 1.e12;
  if (AuAucase) heavyscalefactor = heavyscalefactor * 800.0; // no suppression for now (?)...

  ifstream myfile;

  TGraph *tjet = new TGraph();
  TGraph *tphoton = new TGraph();
  TGraph *tfragphoton = new TGraph();
  TGraph *tpi0 = new TGraph();

  if (! onlyheavy) {

  // =====================================================
  // NLO jets file - default scale
  myfile.open("jets_newphenix_sc1.dat");
  if (! myfile) {
    cout << "No file found!" << endl;
    return;
  }
  for (int index=0; index<38; index++) {
    double pt; double yield;
    myfile >> pt >> yield;
    // rescale points 
    yield = yield * scalefactor*pt;
    tjet->SetPoint(index,pt,yield);
  }
  myfile.close();

  // option for counts with pt > ptvalue
  if (ptgtoption) {
    for (int index=0;index<38;index++) {
      double newsum = 0;
      double x; double y;
      for (int ii=37; ii>=index; ii--) {
	// do the sum above
	tjet->GetPoint(ii,x,y);
	newsum += y;
      }
      tjet->SetPoint(index,x,newsum);
    }
  }

  tjet->SetMarkerStyle(21);
  tjet->SetMarkerColor(2);
  tjet->SetLineColor(2);
  //tjet->Draw("p,l,same");

  // NLO jets file - scale 1/2
  myfile.open("jets_newphenix_sc05.dat");
  if (! myfile) {
    cout << "No file found!" << endl;
    return;
  }
  TGraph *tjet05 = new TGraph();
  for (int index=0; index<38; index++) {
    double pt; double yield;
    myfile >> pt >> yield;
    // rescale points 
    yield = yield * scalefactor*pt;
    tjet05->SetPoint(index,pt,yield);
  }
  myfile.close();

  // option for counts with pt > ptvalue
  if (ptgtoption) {
    for (int index=0;index<38;index++) {
      double newsum = 0;
      double x; double y;
      for (int ii=37; ii>=index; ii--) {
	// do the sum above
	tjet05->GetPoint(ii,x,y);
	newsum += y;
      }
      tjet05->SetPoint(index,x,newsum);
    }
  }

  tjet05->SetMarkerStyle(21);
  tjet05->SetMarkerColor(2);
  tjet05->SetLineColor(2);
  //tjet05->Draw("l,same");

  // NLO jets file - scale 2
  myfile.open("jets_newphenix_sc2.dat");
  if (! myfile) {
    cout << "No file found!" << endl;
    return;
  }
  TGraph *tjet2 = new TGraph();
  for (int index=0; index<38; index++) {
    double pt; double yield;
    myfile >> pt >> yield;
    // rescale points 
    yield = yield * scalefactor*pt;
    tjet2->SetPoint(index,pt,yield);
  }
  myfile.close();

  // option for counts with pt > ptvalue
  if (ptgtoption) {
    for (int index=0;index<38;index++) {
      double newsum = 0;
      double x; double y;
      for (int ii=37; ii>=index; ii--) {
	// do the sum above
	tjet2->GetPoint(ii,x,y);
	newsum += y;
      }
      tjet2->SetPoint(index,x,newsum);
    }
  }

  tjet2->SetMarkerStyle(21);
  tjet2->SetMarkerColor(2);
  tjet2->SetLineColor(2);
  //tjet2->Draw("l,same");


  // =====================================================
  // NLO pizero file - default scale
  myfile.open("pi0_newphenix_sc1.dat");
  if (! myfile) {
    cout << "No file found!" << endl;
    return;
  }

  for (int index=0; index<38; index++) {
    double pt; double yield;
    myfile >> pt >> yield;
    // rescale points 
    yield = yield * scalefactor*pt * pizerosupp;
    tpi0->SetPoint(index,pt,yield);
  }
  myfile.close();

  // option for counts with pt > ptvalue
  if (ptgtoption) {
    for (int index=0;index<38;index++) {
      double newsum = 0;
      double x; double y;
      for (int ii=37; ii>=index; ii--) {
	// do the sum above
	tpi0->GetPoint(ii,x,y);
	newsum += y;
      }
      tpi0->SetPoint(index,x,newsum);
    }
  }

  tpi0->SetMarkerStyle(21);
  tpi0->SetMarkerColor(1);
  tpi0->SetLineColor(1);
  //tpi0->Draw("p,l,same");

  // NLO pi0 file - scale 1/2
  myfile.open("pi0_newphenix_sc05.dat");
  if (! myfile) {
    cout << "No file found!" << endl;
    return;
  }
  TGraph *tpi005 = new TGraph();
  for (int index=0; index<38; index++) {
    double pt; double yield;
    myfile >> pt >> yield;
    // rescale points 
    yield = yield * scalefactor*pt * pizerosupp;
    tpi005->SetPoint(index,pt,yield);
  }
  myfile.close();

  // option for counts with pt > ptvalue
  if (ptgtoption) {
    for (int index=0;index<38;index++) {
      double newsum = 0;
      double x; double y;
      for (int ii=37; ii>=index; ii--) {
	// do the sum above
	tpi005->GetPoint(ii,x,y);
	newsum += y;
      }
      tpi005->SetPoint(index,x,newsum);
    }
  }

  tpi005->SetMarkerStyle(21);
  tpi005->SetMarkerColor(1);
  tpi005->SetLineColor(1);
  //tpi005->Draw("l,same");

  // NLO pi0 file - scale 2
  myfile.open("pi0_newphenix_sc2.dat");
  if (! myfile) {
    cout << "No file found!" << endl;
    return;
  }
  TGraph *tpi02 = new TGraph();
  for (int index=0; index<38; index++) {
    double pt; double yield;
    myfile >> pt >> yield;
    // rescale points 
    yield = yield * scalefactor*pt * pizerosupp;
    tpi02->SetPoint(index,pt,yield);
  }
  myfile.close();

  // option for counts with pt > ptvalue
  if (ptgtoption) {
    for (int index=0;index<38;index++) {
      double newsum = 0;
      double x; double y;
      for (int ii=37; ii>=index; ii--) {
	// do the sum above
	tpi02->GetPoint(ii,x,y);
	newsum += y;
      }
      tpi02->SetPoint(index,x,newsum);
    }
  }

  tpi02->SetMarkerStyle(21);
  tpi02->SetMarkerColor(1);
  tpi02->SetLineColor(1);
  //tpi02->Draw("l,same");


  // =====================================================
  // NLO direct photons file - default scale
  myfile.open("photons_newphenix_sc1.dat");
  if (! myfile) {
    cout << "No file found!" << endl;
    return;
  }

  for (int index=0; index<38; index++) {
    double pt; double yield; double foo;
    myfile >> pt >> yield >> foo >> foo;
    // rescale points 
    yield = yield * scalefactor*pt;
    tphoton->SetPoint(index,pt,yield);
  }
  myfile.close();

  // option for counts with pt > ptvalue
  if (ptgtoption) {
    for (int index=0;index<38;index++) {
      double newsum = 0;
      double x; double y;
      for (int ii=37; ii>=index; ii--) {
	// do the sum above
	tphoton->GetPoint(ii,x,y);
	newsum += y;
      }
      tphoton->SetPoint(index,x,newsum);
    }
  }

  tphoton->SetMarkerStyle(21);
  tphoton->SetMarkerColor(4);
  tphoton->SetLineColor(4);
  tphoton->Draw("p,l,same");

  // NLO photons file - scale 1/2
  myfile.open("photons_newphenix_sc05.dat");
  if (! myfile) {
    cout << "No file found!" << endl;
    return;
  }
  TGraph *tphoton05 = new TGraph();
  for (int index=0; index<38; index++) {
    double pt; double yield; double foo;
    myfile >> pt >> yield >> foo >> foo;
    // rescale points 
    yield = yield * scalefactor*pt;
    tphoton05->SetPoint(index,pt,yield);
  }
  myfile.close();

  // option for counts with pt > ptvalue
  if (ptgtoption) {
    for (int index=0;index<38;index++) {
      double newsum = 0;
      double x; double y;
      for (int ii=37; ii>=index; ii--) {
	// do the sum above
	tphoton05->GetPoint(ii,x,y);
	newsum += y;
      }
      tphoton05->SetPoint(index,x,newsum);
    }
  }

  tphoton05->SetMarkerStyle(21);
  tphoton05->SetMarkerColor(4);
  tphoton05->SetLineColor(4);
  tphoton05->Draw("l,same");

  // NLO photons file - scale 2
  myfile.open("photons_newphenix_sc2.dat");
  if (! myfile) {
    cout << "No file found!" << endl;
    return;
  }
  TGraph *tphoton2 = new TGraph();
  for (int index=0; index<38; index++) {
    double pt; double yield; double foo;
    myfile >> pt >> yield >> foo >> foo;
    // rescale points 
    yield = yield * scalefactor*pt;
    tphoton2->SetPoint(index,pt,yield);
  }
  myfile.close();

  // option for counts with pt > ptvalue
  if (ptgtoption) {
    for (int index=0;index<38;index++) {
      double newsum = 0;
      double x; double y;
      for (int ii=37; ii>=index; ii--) {
	// do the sum above
	tphoton2->GetPoint(ii,x,y);
	newsum += y;
      }
      tphoton2->SetPoint(index,x,newsum);
    }
  }

  tphoton2->SetMarkerStyle(21);
  tphoton2->SetMarkerColor(4);
  tphoton2->SetLineColor(4);
  tphoton2->Draw("l,same");

  // =====================================================
  // NLO fragmentation photons file - default scale
  myfile.open("photons_newphenix_sc1.dat");
  if (! myfile) {
    cout << "No file found!" << endl;
    return;
  }

  for (int index=0; index<38; index++) {
    double pt; double yield; double foo;
    myfile >> pt >> foo >> yield >> foo;
    // rescale points 
    yield = yield * scalefactor*pt;
    tfragphoton->SetPoint(index,pt,yield);
  }
  myfile.close();

  // option for counts with pt > ptvalue
  if (ptgtoption) {
    for (int index=0;index<38;index++) {
      double newsum = 0;
      double x; double y;
      for (int ii=37; ii>=index; ii--) {
	// do the sum above
	tfragphoton->GetPoint(ii,x,y);
	newsum += y;
      }
      tfragphoton->SetPoint(index,x,newsum);
    }
  }

  tfragphoton->SetMarkerStyle(25);
  tfragphoton->SetMarkerColor(4);
  tfragphoton->SetLineColor(4);
  tfragphoton->Draw("p,l,same");

  // NLO photons file - scale 1/2
  myfile.open("photons_newphenix_sc05.dat");
  if (! myfile) {
    cout << "No file found!" << endl;
    return;
  }
  TGraph *tfragphoton05 = new TGraph();
  for (int index=0; index<38; index++) {
    double pt; double yield; double foo;
    myfile >> pt >> foo >> yield >> foo;
    // rescale points 
    yield = yield * scalefactor*pt;
    tfragphoton05->SetPoint(index,pt,yield);
  }
  myfile.close();

  // option for counts with pt > ptvalue
  if (ptgtoption) {
    for (int index=0;index<38;index++) {
      double newsum = 0;
      double x; double y;
      for (int ii=37; ii>=index; ii--) {
	// do the sum above
	tfragphoton05->GetPoint(ii,x,y);
	newsum += y;
      }
      tfragphoton05->SetPoint(index,x,newsum);
    }
  }

  tfragphoton05->SetMarkerStyle(25);
  tfragphoton05->SetMarkerColor(4);
  tfragphoton05->SetLineColor(4);
  tfragphoton05->Draw("l,same");

  // NLO photons file - scale 2
  myfile.open("photons_newphenix_sc2.dat");
  if (! myfile) {
    cout << "No file found!" << endl;
    return;
  }
  TGraph *tfragphoton2 = new TGraph();
  for (int index=0; index<38; index++) {
    double pt; double yield; double foo;
    myfile >> pt >> foo>> yield >> foo;
    // rescale points 
    yield = yield * scalefactor*pt;
    tfragphoton2->SetPoint(index,pt,yield);
  }
  myfile.close();

  // option for counts with pt > ptvalue
  if (ptgtoption) {
    for (int index=0;index<38;index++) {
      double newsum = 0;
      double x; double y;
      for (int ii=37; ii>=index; ii--) {
	// do the sum above
	tfragphoton2->GetPoint(ii,x,y);
	newsum += y;
      }
      tfragphoton2->SetPoint(index,x,newsum);
    }
  }

  tfragphoton2->SetMarkerStyle(25);
  tfragphoton2->SetMarkerColor(4);
  tfragphoton2->SetLineColor(4);
  tfragphoton2->Draw("l,same");

  } // end loop for only heavy

  //---------------------------------------
  TGraph *tcharm = new TGraph();
  TGraph *tbeauty = new TGraph();

  TGraph *tcharmh = new TGraph();
  TGraph *tbeautyh = new TGraph();

  TGraph *tcharme = new TGraph();
  TGraph *tbeautye = new TGraph();

  if (plotheavy) {
    // then also plot heavy quarks from Cacciari
    myfile.open("Heavy/c-quark-eta_less_1.dat");
    if (! myfile) {
      cout << "No file found!" << endl;
      return;
    }
    for (int index=0; index<100; index++) {
      double pt; double yield; 
      myfile >> pt >> yield;
      // rescale points 
      yield = yield * heavyscalefactor;
      tcharm->SetPoint(index,pt,yield);
    }
    myfile.close();
    
    // option for counts with pt > ptvalue
    if (ptgtoption) {
      for (int index=0;index<100;index++) {
	double newsum = 0;
	double x; double y;
	for (int ii=99; ii>=index; ii--) {
	  // do the sum above
	  tcharm->GetPoint(ii,x,y);
	  newsum += y;
	}
	tcharm->SetPoint(index,x,newsum);
      }
    }

    tcharm->SetMarkerStyle(20);
    tcharm->SetMarkerColor(1);
    tcharm->SetLineColor(1);
    tcharm->SetLineWidth(5);
    tcharm->SetLineStyle(1);
    tcharm->Draw("l,,same");

    // now beauty
    myfile.open("Heavy/b-quark-eta_less_1.dat");
    if (! myfile) {
      cout << "No file found!" << endl;
      return;
    }
    for (int index=0; index<100; index++) {
      double pt; double yield; 
      myfile >> pt >> yield;
      // rescale points 
      yield = yield * heavyscalefactor;
      tbeauty->SetPoint(index,pt,yield);
    }
    myfile.close();
    
    // option for counts with pt > ptvalue
    if (ptgtoption) {
      for (int index=0;index<100;index++) {
	double newsum = 0;
	double x; double y;
	for (int ii=99; ii>=index; ii--) {
	  // do the sum above
	  tbeauty->GetPoint(ii,x,y);
	  newsum += y;
	}
	tbeauty->SetPoint(index,x,newsum);
      }
    }

    tbeauty->SetMarkerStyle(20);
    tbeauty->SetMarkerColor(2);
    tbeauty->SetLineColor(2);
    tbeauty->SetLineWidth(5);
    tbeauty->SetLineStyle(1);
    tbeauty->Draw("l,,same");

    // now heavy flavor electrons
    //------------------------------------------

    double heavyraa = 0.2;
    double heavyraaB = 0.5;

    // then also plot heavy quarks->e from Cacciari
    myfile.open("Heavy/c-hadrons-eta_less_1-nppar0.083-electron.dat");
    if (! myfile) {
      cout << "No file found!" << endl;
      return;
    }
    for (int index=0; index<100; index++) {
      double pt; double yield; 
      myfile >> pt >> yield;
      // rescale points 
      yield = yield * heavyscalefactor;
      if (AuAucase) yield = yield * heavyraa;
      tcharme->SetPoint(index,pt,yield);
    }
    myfile.close();
    
    // option for counts with pt > ptvalue
    if (ptgtoption) {
      for (int index=0;index<100;index++) {
	double newsum = 0;
	double x; double y;
	for (int ii=99; ii>=index; ii--) {
	  // do the sum above
	  tcharme->GetPoint(ii,x,y);
	  newsum += y;
	}
	tcharme->SetPoint(index,x,newsum);
      }
    }

    tcharme->SetMarkerStyle(29);
    tcharme->SetMarkerColor(1);
    tcharme->SetLineColor(1);
    tcharme->SetLineWidth(5);
    tcharme->SetLineStyle(3);
    tcharme->Draw("l,,same");

    // now beauty
    myfile.open("Heavy/b-hadrons-eta_less_1-nppar23.5-electron.dat");
    if (! myfile) {
      cout << "No file found!" << endl;
      return;
    }
    for (int index=0; index<100; index++) {
      double pt; double yield; 
      myfile >> pt >> yield;
      // rescale points 
      yield = yield * heavyscalefactor;
      if (AuAucase) yield = yield * heavyraaB;
      tbeautye->SetPoint(index,pt,yield);
    }
    myfile.close();
    
    // option for counts with pt > ptvalue
    if (ptgtoption) {
      for (int index=0;index<100;index++) {
	double newsum = 0;
	double x; double y;
	for (int ii=99; ii>=index; ii--) {
	  // do the sum above
	  tbeautye->GetPoint(ii,x,y);
	  newsum += y;
	}
	tbeautye->SetPoint(index,x,newsum);
      }
    }

    tbeautye->SetMarkerStyle(29);
    tbeautye->SetMarkerColor(2);
    tbeautye->SetLineColor(2);
    tbeautye->SetLineWidth(5);
    tbeautye->SetLineStyle(3);
    tbeautye->Draw("l,,same");


    // now heavy flavor hadrons !!!!!!
    //------------------------------------------
    // then also plot heavy quarks->e from Cacciari
    myfile.open("Heavy/c-hadrons-eta_less_1-nppar0.083.dat");
    if (! myfile) {
      cout << "No file found!" << endl;
      return;
    }
    for (int index=0; index<100; index++) {
      double pt; double yield; 
      myfile >> pt >> yield;
      // rescale points 
      yield = yield * heavyscalefactor;
      if (AuAucase) yield = yield * heavyraa;
      tcharmh->SetPoint(index,pt,yield);
    }
    myfile.close();
    
    // option for counts with pt > ptvalue
    if (ptgtoption) {
      for (int index=0;index<100;index++) {
	double newsum = 0;
	double x; double y;
	for (int ii=99; ii>=index; ii--) {
	  // do the sum above
	  tcharmh->GetPoint(ii,x,y);
	  newsum += y;
	}
	tcharmh->SetPoint(index,x,newsum);
      }
    }

    tcharmh->SetMarkerStyle(24);
    tcharmh->SetMarkerColor(1);
    tcharmh->SetLineColor(1);
    tcharmh->SetLineWidth(5);
    tcharmh->SetLineStyle(2);
    tcharmh->Draw("l,,same");

    // now beauty
    myfile.open("Heavy/b-hadrons-eta_less_1-nppar23.5.dat");
    if (! myfile) {
      cout << "No file found!" << endl;
      return;
    }
    for (int index=0; index<100; index++) {
      double pt; double yield; 
      myfile >> pt >> yield;
      // rescale points 
      yield = yield * heavyscalefactor;
      if (AuAucase) yield = yield * heavyraaB;
      tbeautyh->SetPoint(index,pt,yield);
    }
    myfile.close();
    
    // option for counts with pt > ptvalue
    if (ptgtoption) {
      for (int index=0;index<100;index++) {
	double newsum = 0;
	double x; double y;
	for (int ii=99; ii>=index; ii--) {
	  // do the sum above
	  tbeautyh->GetPoint(ii,x,y);
	  newsum += y;
	}
	tbeautyh->SetPoint(index,x,newsum);
      }
    }

    tbeautyh->SetMarkerStyle(24);
    tbeautyh->SetMarkerColor(2);
    tbeautyh->SetLineColor(2);
    tbeautyh->SetLineWidth(5);
    tbeautyh->SetLineStyle(2);
    tbeautyh->Draw("l,,same");

  } // end of plotheavy if loop

  TLegend *tleg = new TLegend(0.5,0.7,0.9,0.9,"Hard Processes pQCD @ 200 GeV","brNDC");
  if (! onlyheavy) {
    tleg->AddEntry(tjet,"NLO pQCD W. Vogelsang","");
    //tleg->AddEntry(tjet,"Light q + g jets","p,l");
    tleg->AddEntry(tphoton,"Direct #gamma","p,l");
    tleg->AddEntry(tfragphoton,"Fragmentation #gamma","p,l");
    //if (AuAucase) tleg->AddEntry(tpi0,"#pi^{0} (R_{AA}=0.2)","p,l");
    if (!AuAucase) tleg->AddEntry(tpi0,"#pi^{0}","p,l");
  }
  if (onlyheavy) tleg->AddEntry(tcharm,"FONLL pQCD - M. Carriari","");
  if (plotheavy && AuAucase) tleg->AddEntry(tcharm,"Charm Quark Jets (R_{AA}=1.0)","p,l");
  if (plotheavy && !AuAucase) tleg->AddEntry(tcharm,"Charm Quark Jets","p,l");
  if (plotheavy && AuAucase) tleg->AddEntry(tcharmh,"Charm Hadrons (R_{AA}=0.2)","p,l");
  if (plotheavy && !AuAucase) tleg->AddEntry(tcharmh,"Charm Hadrons","p,l");
  if (plotheavy && AuAucase) tleg->AddEntry(tcharme,"Charm #rightarrow Electron (R_{AA}=0.2)","p,l");
  if (plotheavy && !AuAucase) tleg->AddEntry(tcharme,"Charm #rightarrow Electron","p,l");

  if (plotheavy && AuAucase) tleg->AddEntry(tbeauty,"Beauty Quark Jets (R_{AA}=1.0)","p,l");
  if (plotheavy && !AuAucase) tleg->AddEntry(tbeauty,"Beauty Quark Jets","p,l");
  if (plotheavy && AuAucase) tleg->AddEntry(tbeautyh,"Beauty Hadrons (R_{AA}=0.5)","p,l");
  if (plotheavy && !AuAucase) tleg->AddEntry(tbeautyh,"Beauty Hadrons","p,l");
  if (plotheavy && AuAucase) tleg->AddEntry(tbeautye,"Beauty #rightarrow Electron (R_{AA}=0.5)","p,l");
  if (plotheavy && !AuAucase) tleg->AddEntry(tbeautye,"Beauty #rightarrow Electron","p,l");

  tleg->Draw("same");

  //=================================================================================  
  //=================================================================================  
  // now plot the sub-process fraction that contribute to the inclusive jets
  //=================================================================================  
  //=================================================================================  

  TCanvas *c2 = new TCanvas("c2","c2",0,0,700,700);
  c2->cd();

  TNtuple *nt3 = new TNtuple("nt3","nt3","pt:sub1:sub2:sub3");  // subprocesses are qq, qg, gg
  nt3->Reset();
  nt3->ReadFile("jets_newphenix_sc1_ratios.dat");

  TH1F *htemp = new TH1F("htemp","htemp",100,0.0,100.0);
  htemp->SetMaximum(1.0);
  htemp->SetXTitle("Transverse Momentum (GeV/c)");
  htemp->SetYTitle("SubProcess Fraction for Jets in pp@200GeV");
  htemp->DrawCopy();
  nt3->SetMarkerStyle(20);
  nt3->SetMarkerColor(1);
  nt3->SetLineColor(1);
  nt3->Draw("sub1:pt","","p,l,same");
  nt3->SetMarkerColor(4);
  nt3->SetLineColor(4);
  nt3->Draw("sub2:pt","","p,l,same");
  nt3->SetMarkerColor(2);
  nt3->SetLineColor(2);
  nt3->Draw("sub3:pt","","p,l,same");
  TLegend *tl = new TLegend(0.1,0.6,0.5,0.9,"NLO pQCD pp @ 200 GeV inclusive jets","brNDC");
  tl->AddEntry("","Black = Quark-Quark","");
  tl->AddEntry("","Blue  = Quark-Gluon","");
  tl->AddEntry("","Red   = Gluon-Gluon","");
  tl->Draw("same");
  

}
