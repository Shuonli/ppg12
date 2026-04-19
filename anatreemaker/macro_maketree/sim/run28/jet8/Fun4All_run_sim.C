#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 00, 0)

#include <GlobalVariables.C>
#include <G4_Input.C>
#include <Calo_Calib.C>
#include <G4_CEmc_Spacal.C>

// #include <caloreco/CaloGeomMapping.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/RawClusterCNNClassifier.h>
#include <caloreco/RawClusterBuilderTopo.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/CaloGeomMapping.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/RawClusterLikelihoodProfile.h>

#include <calowaveformsim/CaloWaveformSim.h>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>

#include <fun4allraw/Fun4AllPrdfInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <globalvertex/GlobalVertexReco.h>
#include <mbd/MbdReco.h>

#include <phool/recoConsts.h>

#include <clusteriso/ClusterIso.h>

#include <fun4allutils/TimerStats.h>

#include <phool/PHRandomSeed.h>

#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <jetbase/JetCalib.h>

#include <g4jets/TruthJetInput.h>

#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>

#include </sphenix/user/shuhangli/ppg12/anatreemaker/source/CaloAna24.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libtriggervalid.so)
R__LOAD_LIBRARY(libCaloAna24.so)
R__LOAD_LIBRARY(libclusteriso.so)
R__LOAD_LIBRARY(libfun4allutils.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libjetbackground.so)
#endif

void Fun4All_run_sim(
    const int nEvents = 0,
    const string &inputFile0 = "test.list",
    //const string &inputFile1 = "dst_calo_cluster.list",
    //const string &inputFile3 = "dst_mbd_epd.list",
    //const string &inputFile4 = "dst_truth_jet.list",
    //const string &inputFile2 = "dst_truth.list",

    const string &outputFile = "output_sim.root",

    const string &outputDSTFile = "DST_CALO_WAVEFORM_pp-0000000011-00000.root",
    const string &outDSTdir = ".",
    const string &cdbtag = "MDC2")
{

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);
  recoConsts *rc = recoConsts::instance();

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", cdbtag);
  rc->set_uint64Flag("TIMESTAMP", 28);
  rc->set_IntFlag("RUNNUMBER", 28);
  CDBInterface::instance()->Verbosity(1);

  Input::VERBOSITY = 1;
  Input::READHITS = true;
  INPUTREADHITS::listfile[0] = inputFile0;
  //INPUTREADHITS::listfile[1] = inputFile1;
  //INPUTREADHITS::listfile[2] = inputFile2;
  //INPUTREADHITS::listfile[3] = inputFile3;
  //INPUTREADHITS::listfile[4] = inputFile4;

  InputInit();
  InputRegister();

  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  Enable::DSTOUT = false;
  Enable::DSTOUT_COMPRESS = false;
  DstOut::OutputDir = outDSTdir;
  DstOut::OutputFile = outputDSTFile;

  MbdReco *mbdreco = new MbdReco();
  se->registerSubsystem(mbdreco);

  GlobalVertexReco *gblvertex = new GlobalVertexReco();
  gblvertex->Verbosity(0);
  se->registerSubsystem(gblvertex);
  //======================
  // What to run
  //======================
  // the only reason that we are including the calo_cluster node is we want to use the CEMC geom node from it,
  // but we have calib node name confliting(it has a towerinfov1 calib node with the same name we want to usebut we want to make it v2) if we do that...
  // so I will call the cemc tower reco here just to have the geom node.
  // by doing this it also remove the dependncy for running the calo_cluster before this pass so we can process it independently from G4Hits ;) and all of our output node name is exactly same with real data
  if(false){
  CEMC_Cells();

  CaloWaveformSim *caloWaveformSim = new CaloWaveformSim("HCALOUTWaveformSim");
  caloWaveformSim->set_detector_type(CaloTowerDefs::HCALOUT);
  caloWaveformSim->set_detector("HCALOUT");
  caloWaveformSim->set_nsamples(12);
  caloWaveformSim->set_pedestalsamples(12);
  caloWaveformSim->set_timewidth(0.2);
  caloWaveformSim->set_peakpos(6);
  //caloWaveformSim->set_pedestal_scale(0.69);
  // caloWaveformSim->Verbosity(2);
  // caloWaveformSim->set_noise_type(CaloWaveformSim::NOISE_NONE);
  se->registerSubsystem(caloWaveformSim);

  caloWaveformSim = new CaloWaveformSim("HCALINWaveformSim");
  caloWaveformSim->set_detector_type(CaloTowerDefs::HCALIN);
  caloWaveformSim->set_detector("HCALIN");
  caloWaveformSim->set_nsamples(12);
  caloWaveformSim->set_pedestalsamples(12);
  caloWaveformSim->set_timewidth(0.2);
  caloWaveformSim->set_peakpos(6);
  //  caloWaveformSim->set_noise_type(CaloWaveformSim::NOISE_NONE);
  se->registerSubsystem(caloWaveformSim);

  caloWaveformSim = new CaloWaveformSim("CEMCWaveformSim");
  caloWaveformSim->set_detector_type(CaloTowerDefs::CEMC);
  caloWaveformSim->set_detector("CEMC");
  caloWaveformSim->set_nsamples(12);
  caloWaveformSim->set_pedestalsamples(12);
  caloWaveformSim->set_timewidth(0.2);
  caloWaveformSim->set_peakpos(6);
  caloWaveformSim->set_pedestal_scale(0.69);

  //  caloWaveformSim->set_noise_type(CaloWaveformSim::NOISE_NONE);

  caloWaveformSim->get_light_collection_model().load_data_file(
      string(getenv("CALIBRATIONROOT")) +
          string("/CEMC/LightCollection/Prototype3Module.xml"),
      "data_grid_light_guide_efficiency", "data_grid_fiber_trans");

  se->registerSubsystem(caloWaveformSim);

  CaloTowerBuilder *ca2 = new CaloTowerBuilder("HCALOUTTowerBuilder");
  ca2->set_detector_type(CaloTowerDefs::HCALOUT);
  ca2->set_nsamples(12);
  ca2->set_dataflag(false);
  ca2->set_processing_type(CaloWaveformProcessing::TEMPLATE);
  ca2->set_builder_type(CaloTowerDefs::kWaveformTowerSimv1);
  // 30 ADC SZS
  ca2->set_softwarezerosuppression(true, 30);
  se->registerSubsystem(ca2);

  ca2 = new CaloTowerBuilder("HCALINTowerBuilder");
  ca2->set_detector_type(CaloTowerDefs::HCALIN);
  ca2->set_nsamples(12);
  ca2->set_dataflag(false);
  ca2->set_processing_type(CaloWaveformProcessing::TEMPLATE);
  ca2->set_builder_type(CaloTowerDefs::kWaveformTowerSimv1);
  ca2->set_softwarezerosuppression(true, 30);
  se->registerSubsystem(ca2);

  ca2 = new CaloTowerBuilder("CEMCTowerBuilder");
  ca2->set_detector_type(CaloTowerDefs::CEMC);
  ca2->set_nsamples(12);
  ca2->set_dataflag(false);
  ca2->set_processing_type(CaloWaveformProcessing::TEMPLATE);
  ca2->set_builder_type(CaloTowerDefs::kWaveformTowerSimv1);
  // a large uniform ZS threshold for CEMC, 60 ADC now
  ca2->set_softwarezerosuppression(true, 60);
  se->registerSubsystem(ca2);

  /////////////////////////////////////////////////////
  // Set status of towers, Calibrate towers,  Cluster
  /////////////////////////////////////////////////////
  std::cout << "status setters" << std::endl;
  CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
  statusEMC->set_detector_type(CaloTowerDefs::CEMC);
  //statusEMC->set_time_cut(1);
  se->registerSubsystem(statusEMC);

  CaloTowerStatus *statusHCalIn = new CaloTowerStatus("HCALINSTATUS");
  statusHCalIn->set_detector_type(CaloTowerDefs::HCALIN);
  //statusHCalIn->set_time_cut(2);
  se->registerSubsystem(statusHCalIn);

  CaloTowerStatus *statusHCALOUT = new CaloTowerStatus("HCALOUTSTATUS");
  statusHCALOUT->set_detector_type(CaloTowerDefs::HCALOUT);
  //statusHCALOUT->set_time_cut(2);
  se->registerSubsystem(statusHCALOUT);

  ////////////////////
  // Calibrate towers
  std::cout << "Calibrating EMCal" << std::endl;
  CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
  calibEMC->set_detector_type(CaloTowerDefs::CEMC);
  calibEMC->set_outputNodePrefix("TOWERINFO_CALIB_");
  se->registerSubsystem(calibEMC);

  std::cout << "Calibrating OHcal" << std::endl;
  CaloTowerCalib *calibOHCal = new CaloTowerCalib("HCALOUTCALIB");
  calibOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
  calibOHCal->set_outputNodePrefix("TOWERINFO_CALIB_");
  se->registerSubsystem(calibOHCal);

  std::cout << "Calibrating IHcal" << std::endl;
  CaloTowerCalib *calibIHCal = new CaloTowerCalib("HCALINCALIB");
  calibIHCal->set_detector_type(CaloTowerDefs::HCALIN);
  calibIHCal->set_outputNodePrefix("TOWERINFO_CALIB_");
  se->registerSubsystem(calibIHCal);

  std::cout<<"runnumber is:" << rc->get_IntFlag("RUNNUMBER")<<std::endl;
  ////////////////
  // MC Calibration
  std::string MC_Calib = CDBInterface::instance()->getUrl("CEMC_MC_RECALIB");
  if (MC_Calib.empty())
  {
    std::cout << "No MC calibration found :( )" << std::endl;
    gSystem->Exit(0);
  }
  //CaloTowerCalib *calibEMC_MC = new CaloTowerCalib("CEMCCALIB_MC");
  //calibEMC_MC->set_detector_type(CaloTowerDefs::CEMC);
  //calibEMC_MC->set_inputNodePrefix("TOWERINFO_CALIB_");
  //calibEMC_MC->set_outputNodePrefix("TOWERINFO_CALIB_");
  //calibEMC_MC->set_directURL(MC_Calib);
  //calibEMC_MC->set_doZScrosscalib(false);
  }
  //--------------
  // Timing module is last to register
  //--------------
  TimerStats *ts = new TimerStats();
  ts->OutFileName("jobtime.root");
  se->registerSubsystem(ts);

  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();
  TRandom3 randGen;
  // get seed
  unsigned int seed = PHRandomSeed();
  randGen.SetSeed(seed);
  // a int from 0 to 3259
  int sequence = randGen.Integer(3260);
  // pad the name
  std::ostringstream opedname;
  opedname << "pedestal-54256-0" << std::setw(4) << std::setfill('0') << sequence << ".root";

  std::string pedestalname = opedname.str();

  Fun4AllInputManager *hitsin = new Fun4AllNoSyncDstInputManager("DST2");
  hitsin->AddFile(pedestalname);
  hitsin->Repeat();
  se->registerInputManager(hitsin);

  /*
  RawClusterCNNClassifier *cnn = new RawClusterCNNClassifier();
  se->registerSubsystem(cnn);
  */

  Process_Calo_Calib();

  //////////////////
  // Clusters
  std::string clusternodename = "CLUSTERINFO_CEMC_NO_SPLIT";
  std::cout << "Building clusters" << std::endl;
  RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
  ClusterBuilder->Detector("CEMC");
  ClusterBuilder->set_threshold_energy(0.070); // for when using basic calibration
  std::string emc_prof = getenv("CALIBRATIONROOT");
  emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  ClusterBuilder->LoadProfile(emc_prof);
  ClusterBuilder->setSubclusterSplitting(false);
  ClusterBuilder->setOutputClusterNodeName(clusternodename);
  ClusterBuilder->set_UseTowerInfo(1); // to use towerinfo objects rather than old RawTower
  se->registerSubsystem(ClusterBuilder);
  if(false)
  {
    // cluster new prob
    RawClusterLikelihoodProfile *ClusterProleBuilder = new RawClusterLikelihoodProfile("RawClusterLikelihoodProleGamma");
    ClusterProleBuilder->set_profile_filepath("/sphenix/user/jpark4/CDBfiles/EMCalProb/EMCal_ProfileLikelihoodD2_Thres70MeV.root"); // default set to single gamma
    ClusterProleBuilder->set_outputNodeName("CLUSTERINFO_CEMC_NO_SPLIT");                                                           // could keep the same name to overwrite with new prob. values
    ClusterProleBuilder->set_profile_dimension(5);                                                                                  // 5x5
    se->registerSubsystem(ClusterProleBuilder);
  }
  {
    ClusterIso *cliso2 = new ClusterIso("ClusterIso2", 5, 2, false, true);
    cliso2->Verbosity(0);
    cliso2->set_cluster_node_name(clusternodename);
    cliso2->set_use_towerinfo(true);
    cliso2->setMinTowerEnergy(0.070);
    se->registerSubsystem(cliso2);

    // 3
    ClusterIso *cliso3 = new ClusterIso("ClusterIso3", 5, 3, false, true);
    cliso3->Verbosity(0);
    cliso3->set_use_towerinfo(true);
    cliso3->set_cluster_node_name(clusternodename);
    cliso3->setMinTowerEnergy(0.070);
    se->registerSubsystem(cliso3);
    // 4
    ClusterIso *cliso4 = new ClusterIso("ClusterIso4", 5, 4, false, true);
    // set verbosity to max int to see all the printouts
    cliso4->Verbosity(0);
    cliso4->set_use_towerinfo(true);
    cliso4->set_cluster_node_name(clusternodename);
    cliso4->setMinTowerEnergy(0.070);
    se->registerSubsystem(cliso4);
  }
  clusternodename = "CLUSTERINFO_CEMC";
  if(false)
  {
    // cluster new prob
    RawClusterLikelihoodProfile *ClusterProleBuilder = new RawClusterLikelihoodProfile("RawClusterLikelihoodProleGamma");
    ClusterProleBuilder->set_profile_filepath("/sphenix/user/jpark4/CDBfiles/EMCalProb/EMCal_ProfileLikelihoodD2_Thres70MeV.root"); // default set to single gamma
    ClusterProleBuilder->set_outputNodeName("CLUSTERINFO_CEMC");                                                                    // could keep the same name to overwrite with new prob. values
    ClusterProleBuilder->set_profile_dimension(5);                                                                                  // 5x5
    se->registerSubsystem(ClusterProleBuilder);
  }

  {
    ClusterIso *cliso2 = new ClusterIso("ClusterIso2", 5, 2, false, true);
    cliso2->Verbosity(0);
    cliso2->set_cluster_node_name(clusternodename);
    cliso2->set_use_towerinfo(true);
    cliso2->setMinTowerEnergy(0.070);
    se->registerSubsystem(cliso2);

    // 3
    ClusterIso *cliso3 = new ClusterIso("ClusterIso3", 5, 3, false, true);
    cliso3->Verbosity(0);
    cliso3->set_use_towerinfo(true);
    cliso3->set_cluster_node_name(clusternodename);
    cliso3->setMinTowerEnergy(0.070);
    se->registerSubsystem(cliso3);
    // 4
    ClusterIso *cliso4 = new ClusterIso("ClusterIso4", 5, 4, false, true);
    // set verbosity to max int to see all the printouts
    cliso4->Verbosity(0);
    cliso4->set_use_towerinfo(true);
    cliso4->set_cluster_node_name(clusternodename);
    cliso4->setMinTowerEnergy(0.070);
    se->registerSubsystem(cliso4);
  }

  // Topo clusters (all calo, with splitting)
  {
    int verbosity = 0;
    RawClusterBuilderTopo *ClusterBuilder = new RawClusterBuilderTopo("HcalRawClusterBuilderTopo");
    ClusterBuilder->Verbosity(verbosity);
    ClusterBuilder->set_nodename("TOPOCLUSTER_ALLCALO");
    ClusterBuilder->set_enable_HCal(true);
    ClusterBuilder->set_enable_EMCal(true);
    ClusterBuilder->set_noise(0.0053, 0.0351, 0.0684); // 3sigma of pedestal noise
    ClusterBuilder->set_significance(4.0, 2.0, 1.0);
    ClusterBuilder->allow_corner_neighbor(true);
    ClusterBuilder->set_do_split(true);
    ClusterBuilder->set_minE_local_max(1.0, 2.0, 0.5);
    ClusterBuilder->set_R_shower(0.025);
    ClusterBuilder->set_use_only_good_towers(true);
    ClusterBuilder->set_absE(true);
    se->registerSubsystem(ClusterBuilder);
  }

  // Topo clusters (all calo, soft thresholds)
  {
    int verbosity = 0;
    RawClusterBuilderTopo *ClusterBuilder2 = new RawClusterBuilderTopo("HcalRawClusterBuilderTopo_Soft");
    ClusterBuilder2->Verbosity(verbosity);
    ClusterBuilder2->set_nodename("TOPOCLUSTER_ALLCALO_SOFT");
    ClusterBuilder2->set_enable_HCal(true);
    ClusterBuilder2->set_enable_EMCal(true);
    ClusterBuilder2->set_noise(0.0053, 0.0351, 0.0684); // 3sigma of pedestal noise
    ClusterBuilder2->set_significance(3.0, 2.0, 0.0);
    ClusterBuilder2->allow_corner_neighbor(true);
    ClusterBuilder2->set_do_split(true);
    ClusterBuilder2->set_minE_local_max(1.0, 2.0, 0.5);
    ClusterBuilder2->set_R_shower(0.025);
    ClusterBuilder2->set_use_only_good_towers(true);
    ClusterBuilder2->set_absE(true);
    se->registerSubsystem(ClusterBuilder2);
  }

  // jet stuff
  std::string jetreco_input_prefix = "TOWERINFO_CALIB";
  RetowerCEMC *_retowerCEMC = new RetowerCEMC();
  _retowerCEMC->Verbosity(0);
  _retowerCEMC->set_towerinfo(true);
  _retowerCEMC->set_frac_cut(0.5); // fraction of retower that must be masked to mask the full retower
  _retowerCEMC->set_towerNodePrefix(jetreco_input_prefix);
  se->registerSubsystem(_retowerCEMC);

  JetReco *_jetRecoUnsub = new JetReco();
  std::vector<float> doRecoJet_radius = {0.4};
  _jetRecoUnsub->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER, jetreco_input_prefix));
  _jetRecoUnsub->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO, jetreco_input_prefix));
  _jetRecoUnsub->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO, jetreco_input_prefix));
  for (int ir = 0; ir < doRecoJet_radius.size(); ++ir)
  {
    _jetRecoUnsub->add_algo(new FastJetAlgoSub(Jet::ANTIKT, doRecoJet_radius[ir]), "AntiKt_unsubtracted_r0" + std::to_string((int)(10 * doRecoJet_radius[ir])));
  }
  _jetRecoUnsub->set_algo_node("ANTIKT");
  _jetRecoUnsub->set_input_node("TOWER");
  _jetRecoUnsub->Verbosity(0);
  se->registerSubsystem(_jetRecoUnsub);

  JetCalib *jetCalib04 = new JetCalib("JetCalib04");
  jetCalib04->set_InputNode("AntiKt_unsubtracted_r04");
  jetCalib04->set_OutputNode("AntiKt_unsubtracted_r04_calib");
  jetCalib04->set_JetRadius(0.4);
  jetCalib04->set_ZvrtxNode("GlobalVertexMap");
  jetCalib04->set_ApplyZvrtxDependentCalib(true);
  jetCalib04->set_ApplyEtaDependentCalib(true);
  se->registerSubsystem(jetCalib04);


  int verbosity = 0;
  JetReco *truthjets4 = new JetReco("TRUTHJETRECO4");
  truthjets4->add_input(new TruthJetInput(Jet::PARTICLE));
  truthjets4->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Truth_r04");
  truthjets4->set_algo_node("ANTIKT");
  truthjets4->set_input_node("TRUTH");
  truthjets4->Verbosity(verbosity);
  se->registerSubsystem(truthjets4);

  CaloAna24 *caloana24 = new CaloAna24();
  // caloana24->set_isSingleParticle(true);
  // caloana24->set_clusterpTmin(1);
  // caloana24->set_particlepTmin(1);
  caloana24->set_jet_node_name("AntiKt_unsubtracted_r04_calib");
  se->registerSubsystem(caloana24);

  se->run(nEvents);
  CDBInterface::instance()->Print(); // print used DB files
  se->End();
  cout << "JOB COMPLETE." << endl;
  se->PrintTimer();
  gSystem->Exit(0);
}
