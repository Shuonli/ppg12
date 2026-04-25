#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>
#include <fun4all/Fun4AllUtils.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <stdlib.h>

#include <g4main/PHG4Reco.h>
#include <g4main/PHG4TruthSubsystem.h>

#include <GlobalVariables.C>

#include <G4_Input.C>
#include <phpythia8/PHPy8JetTrigger.h>
#include <phpythia8/PHPy8ParticleTrigger.h>

// #include <jethistogrammer/jetHistogrammer.h>

#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <g4jets/TruthJetInput.h>
#include <time.h>

// R__LOAD_LIBRARY(libjetHistogrammer.so)

void Fun4AllPythia(int nEvents = 10,
                   const string &photontrigger = "PhotonJet20",
                   int doCrossSection = 1,
                   const string &outname = "jetHistogrammer")
{

  clock_t tStart = clock();

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  Input::PYTHIA8 = true;
  Input::VERBOSITY = 0;

  float jetTrig = 0;

  string pythia8_config_file = string(getenv("CALIBRATIONROOT")) + "/Generators/JetStructure_TG/";
  if (photontrigger == "PhotonJet5")
  {
    pythia8_config_file += "phpythia8_5GeV_JS_MDC2.cfg";
  }
  else if (photontrigger == "PhotonJet10")
  {
    pythia8_config_file += "phpythia8_15GeV_JS_MDC2.cfg";
  }
  else if (photontrigger == "PhotonJet20")
  {
    pythia8_config_file += "phpythia8_30GeV_JS_MDC2.cfg";
  }
  else
  {
    std::cout << "Invalid photon trigger " << photontrigger << std::endl;
    gSystem->Exit(1);
  }
  PYTHIA8::config_file = pythia8_config_file;

  InputInit();

  PHPy8ParticleTrigger *p8_photon_jet_trigger = new PHPy8ParticleTrigger();
  p8_photon_jet_trigger->AddParticles(22);
  p8_photon_jet_trigger->SetEtaHighLow(1.5, -1.5);     // sample a rapidity range higher than the sPHENIX tracking pseudorapidity
  p8_photon_jet_trigger->SetStableParticleOnly(false); // process unstable particles that include quarks

  std::vector<int> partentsId{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1};
  p8_photon_jet_trigger->AddParents(partentsId);

  if (photontrigger == "PhotonJet5")
    p8_photon_jet_trigger->SetPtLow(5.);
  else if (photontrigger == "PhotonJet10")
    p8_photon_jet_trigger->SetPtLow(10.);
  else if (photontrigger == "PhotonJet20")
    p8_photon_jet_trigger->SetPtLow(20.);
  else
  {
    std::cout << "Invalid photon trigger " << photontrigger << std::endl;
    gSystem->Exit(1);
  }

  p8_photon_jet_trigger->PrintConfig();
  // Input::ApplysPHENIXBeamParameter(INPUTGENERATOR::Pythia8);
  INPUTGENERATOR::Pythia8->register_trigger(p8_photon_jet_trigger);
  INPUTGENERATOR::Pythia8->Verbosity(doCrossSection * 1);

  InputRegister();
  InputManagers();
  PHG4Reco *reco = new PHG4Reco();
  reco->set_field(0);
  reco->SetWorldMaterial("G4_Galactic");
  reco->SetWorldSizeX(100);
  reco->SetWorldSizeY(100);
  reco->SetWorldSizeZ(100);
  reco->save_DST_geometry(false);

  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  reco->registerSubsystem(truth);
  se->registerSubsystem(reco);

  // JetReco *truthJetReco = new JetReco();
  // TruthJetInput *jetInput = new TruthJetInput(Jet::PARTICLE);
  // truthJetReco -> add_input(jetInput);
  // truthJetReco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Truth_r04");
  // truthJetReco->set_algo_node("ANTIKT");
  // truthJetReco->set_input_node("TRUTH");
  // truthJetReco->Verbosity(0);

  // se -> registerSubsystem(truthJetReco);

  // string out = outname + Jet_Trigger + ".root";
  // jetHistogrammer *jetHistos = new jetHistogrammer("jetHistogrammer",out.c_str());
  // jetHistos -> setJetTrig(jetTrig);
  // se -> registerSubsystem(jetHistos);

  se->run(nEvents);

  se->End();
  std::cout << "All done" << std::endl;

  delete se;
  std::cout << "Total runtime: " << double(clock() - tStart) / (double)CLOCKS_PER_SEC << std::endl;
  ;

  gSystem->Exit(0);
  return 0;
}
