#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 00, 0)

#include <GlobalVariables.C>
#include <G4_Input.C>
#include <Calo_Calib.C>

// #include <caloreco/CaloGeomMapping.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/RawClusterCNNClassifier.h>
#include <caloreco/RawClusterBuilderTopo.h>
#include <caloreco/RawClusterBuilderTemplate.h>
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
#include <fun4all/Fun4AllUtils.h>

#include <globalvertex/GlobalVertexReco.h>
#include <mbd/MbdReco.h>
#include <phool/recoConsts.h>

#include <clusteriso/ClusterIso.h>

#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>

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
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libtriggervalid.so)
R__LOAD_LIBRARY(libCaloAna24.so)
R__LOAD_LIBRARY(libclusteriso.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libjetbackground.so)
#endif

void Fun4All_runDST(
    const std::string &fname0 = "dstLists/dst_jet-00047289.list",
    const std::string &fname1 = "dstLists/dst_jetcalo-00047289.list",
    // 48342
    // 48244 MBD only
    int nEvents = 0,
    const std::string &dbtag = "ProdA_2024")
{

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);

  recoConsts *rc = recoConsts::instance();

  ifstream file(fname0);
  string first_file;
  getline(file, first_file);

  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(first_file);
  int runnumber = runseg.first;

  // conditions DB flags and timestamp
  rc->set_StringFlag("CDB_GLOBALTAG", dbtag);
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  CDBInterface::instance()->Verbosity(1);

  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  Fun4AllInputManager *In = new Fun4AllDstInputManager("in0");
  In->Verbosity(1);
  In->AddListFile(fname0);
  se->registerInputManager(In);

  In = new Fun4AllDstInputManager("in1");
  In->Verbosity(1);
  In->AddListFile(fname1);
  se->registerInputManager(In);

  // global reco
  /*
  GlobalVertexReco* gblvertex = new GlobalVertexReco();
  gblvertex->Verbosity(0);
  se->registerSubsystem(gblvertex);
  */
  /*
  MbdReco *mbdreco = new MbdReco();
  se->registerSubsystem(mbdreco);


  GlobalVertexReco* gblvertex = new GlobalVertexReco();
  gblvertex->Verbosity(0);
  se->registerSubsystem(gblvertex);
  */

  Process_Calo_Calib();

  /*
  ClusterIso *cliso2 = new ClusterIso("ClusterIso2", 5, 2, false, true);
  cliso2->Verbosity(0);
  cliso2->set_use_towerinfo(true);
  se->registerSubsystem(cliso2);

  // 3
  ClusterIso *cliso3 = new ClusterIso("ClusterIso3", 5, 3, false, true);
  cliso3->Verbosity(0);
  cliso3->set_use_towerinfo(true);
  se->registerSubsystem(cliso3);
  // 4
  ClusterIso *cliso4 = new ClusterIso("ClusterIso4", 5, 4, false, true);
  cliso4->Verbosity(0);
  cliso4->set_use_towerinfo(true);
  se->registerSubsystem(cliso4);
  */
  /*
   std::string clusternodename = "TOPOCLUSTER_EMCAL_SPLIT";
   RawClusterBuilderTopo *ClusterBuilder1 = new RawClusterBuilderTopo("HcalRawClusterBuilderTopo1");
   ClusterBuilder1->Verbosity(0);
   ClusterBuilder1->set_nodename(clusternodename);
   ClusterBuilder1->set_enable_HCal(false);
   ClusterBuilder1->set_enable_EMCal(true);
   ClusterBuilder1->set_noise(0.01, 0.03, 0.03);
   ClusterBuilder1->set_significance(4.0, 2.0, 1.0);
   ClusterBuilder1->allow_corner_neighbor(true);
   ClusterBuilder1->set_do_split(true);
   ClusterBuilder1->set_minE_local_max(1.0, 2.0, 0.5);
   ClusterBuilder1->set_R_shower(0.025);
   ClusterBuilder1->set_use_only_good_towers(true);
   ClusterBuilder1->set_absE(false);
   se->registerSubsystem(ClusterBuilder1);

   cliso2 = new ClusterIso("ClusterIso2", 5, 2, false, true);
   cliso2->Verbosity(0);
   cliso2->set_cluster_node_name(clusternodename);
   cliso2->set_use_towerinfo(true);
   se->registerSubsystem(cliso2);

   // 3
   cliso3 = new ClusterIso("ClusterIso3", 5, 3, false, true);
   cliso3->Verbosity(0);
   cliso3->set_use_towerinfo(true);
   cliso3->set_cluster_node_name(clusternodename);
   se->registerSubsystem(cliso3);
   // 4
   cliso4 = new ClusterIso("ClusterIso4", 5, 4, false, true);
   cliso4->Verbosity(0);
   cliso4->set_use_towerinfo(true);
   cliso4->set_cluster_node_name(clusternodename);
   se->registerSubsystem(cliso4);



   clusternodename = "TOPOCLUSTER_EMCAL";

   ClusterBuilder1 = new RawClusterBuilderTopo("HcalRawClusterBuilderTopo1");
   ClusterBuilder1->Verbosity(0);
   ClusterBuilder1->set_nodename(clusternodename);
   ClusterBuilder1->set_enable_HCal(false);
   ClusterBuilder1->set_enable_EMCal(true);
   ClusterBuilder1->set_noise(0.01, 0.03, 0.03);
   ClusterBuilder1->set_significance(4.0, 2.0, 1.0);
   ClusterBuilder1->allow_corner_neighbor(true);
   ClusterBuilder1->set_do_split(false);
   ClusterBuilder1->set_minE_local_max(1.0, 2.0, 0.5);
   ClusterBuilder1->set_R_shower(0.025);
   ClusterBuilder1->set_use_only_good_towers(true);
   ClusterBuilder1->set_absE(false);
   se->registerSubsystem(ClusterBuilder1);

   cliso2 = new ClusterIso("ClusterIso2", 5, 2, false, true);
   cliso2->Verbosity(0);
   cliso2->set_cluster_node_name(clusternodename);
   cliso2->set_use_towerinfo(true);
   se->registerSubsystem(cliso2);

   // 3
   cliso3 = new ClusterIso("ClusterIso3", 5, 3, false, true);
   cliso3->Verbosity(0);
   cliso3->set_use_towerinfo(true);
   cliso3->set_cluster_node_name(clusternodename);
   se->registerSubsystem(cliso3);
   // 4
   cliso4 = new ClusterIso("ClusterIso4", 5, 4, false, true);
   cliso4->Verbosity(0);
   cliso4->set_use_towerinfo(true);
   cliso4->set_cluster_node_name(clusternodename);
   se->registerSubsystem(cliso4);
   */
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

  if(false){
    // cluster new prob
    RawClusterLikelihoodProfile *ClusterProleBuilder = new RawClusterLikelihoodProfile("RawClusterLikelihoodProleGamma");
    ClusterProleBuilder->set_profile_filepath("/sphenix/user/jpark4/CDBfiles/EMCalProb/EMCal_ProfileLikelihoodD2_Thres70MeV.root"); // default set to single gamma
    ClusterProleBuilder->set_inputNodeName("CLUSTERINFO_CEMC_NO_SPLIT");
    ClusterProleBuilder->set_outputNodeName("CLUSTERINFO_CEMC_NO_SPLIT");                                                           // could keep the same name to overwrite with new prob. values
    ClusterProleBuilder->set_profile_dimension(5); // 5x5
    ClusterProleBuilder->set_min_cluster_e(4.0);                                                                               
    se->registerSubsystem(ClusterProleBuilder);
  }
  {
    ClusterIso *cliso2 = new ClusterIso("ClusterIso2", 5, 2, false, true);
    cliso2->Verbosity(0);
    cliso2->set_cluster_node_name(clusternodename);
    cliso2->set_use_towerinfo(true);
    se->registerSubsystem(cliso2);

    // 3
    ClusterIso *cliso3 = new ClusterIso("ClusterIso3", 5, 3, false, true);
    cliso3->Verbosity(0);
    cliso3->set_use_towerinfo(true);
    cliso3->set_cluster_node_name(clusternodename);
    se->registerSubsystem(cliso3);
    // 4
    ClusterIso *cliso4 = new ClusterIso("ClusterIso4", 5, 4, false, true);
    // set verbosity to max int to see all the printouts
    cliso4->Verbosity(0);
    cliso4->set_use_towerinfo(true);
    cliso4->set_cluster_node_name(clusternodename);
    se->registerSubsystem(cliso4);
  }

  clusternodename = "CLUSTERINFO_CEMC";
  if(false){
    // cluster new prob
    RawClusterLikelihoodProfile *ClusterProleBuilder = new RawClusterLikelihoodProfile("RawClusterLikelihoodProleGamma");
    ClusterProleBuilder->set_profile_filepath("/sphenix/user/jpark4/CDBfiles/EMCalProb/EMCal_ProfileLikelihoodD2_Thres70MeV.root"); // default set to single gamma
    ClusterProleBuilder->set_outputNodeName("CLUSTERINFO_CEMC");                                                                    // could keep the same name to overwrite with new prob. values
    ClusterProleBuilder->set_profile_dimension(5);
    ClusterProleBuilder->set_min_cluster_e(4.0);                                                                    // 5x5
    se->registerSubsystem(ClusterProleBuilder);
  }

  {
    ClusterIso *cliso2 = new ClusterIso("ClusterIso2", 5, 2, false, true);
    cliso2->Verbosity(0);
    cliso2->set_cluster_node_name(clusternodename);
    cliso2->set_use_towerinfo(true);
    se->registerSubsystem(cliso2);

    // 3
    ClusterIso *cliso3 = new ClusterIso("ClusterIso3", 5, 3, false, true);
    cliso3->Verbosity(0);
    cliso3->set_use_towerinfo(true);
    cliso3->set_cluster_node_name(clusternodename);
    se->registerSubsystem(cliso3);
    // 4
    ClusterIso *cliso4 = new ClusterIso("ClusterIso4", 5, 4, false, true);
    // set verbosity to max int to see all the printouts
    cliso4->Verbosity(0);
    cliso4->set_use_towerinfo(true);
    cliso4->set_cluster_node_name(clusternodename);
    se->registerSubsystem(cliso4);
  }

  /*
    ClusterProleBuilder = new RawClusterLikelihoodProfile("RawClusterLikelihoodProleGamma");
    ClusterProleBuilder->set_profile_filepath("/sphenix/user/jpark4/CDBfiles/EMCalProb/ProfileLikelihoodD2_single_gamma.root"); // default set to single gamma
    ClusterProleBuilder->set_outputNodeName("CLUSTERINFO_CEMC"); // could keep the same name to overwrite with new prob. values
    ClusterProleBuilder->set_profile_dimension(5); // 5x5
    se->registerSubsystem(ClusterProleBuilder);

    ClusterProleBuilder = new RawClusterLikelihoodProfile("RawClusterLikelihoodProleGamma");
    ClusterProleBuilder->set_profile_filepath("/sphenix/user/jpark4/CDBfiles/EMCalProb/ProfileLikelihoodD2_single_gamma.root"); // default set to single gamma
    ClusterProleBuilder->set_outputNodeName("TOPOCLUSTER_EMCAL_SPLIT"); // could keep the same name to overwrite with new prob. values
    ClusterProleBuilder->set_profile_dimension(5); // 5x5
    se->registerSubsystem(ClusterProleBuilder);

    ClusterProleBuilder = new RawClusterLikelihoodProfile("RawClusterLikelihoodProleGamma");
    ClusterProleBuilder->set_profile_filepath("/sphenix/user/jpark4/CDBfiles/EMCalProb/ProfileLikelihoodD2_single_gamma.root"); // default set to single gamma
    ClusterProleBuilder->set_outputNodeName("TOPOCLUSTER_EMCAL"); // could keep the same name to overwrite with new prob. values
    ClusterProleBuilder->set_profile_dimension(5); // 5x5
    se->registerSubsystem(ClusterProleBuilder);
  */
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

  std::string filename = first_file.substr(first_file.find_last_of("/\\") + 1);
  std::string OutFile = Form("OUTTREE_%s", filename.c_str());
  CaloAna24 *caloana24 = new CaloAna24("CaloAna24", OutFile);
  caloana24->set_isMC(false);
  se->registerSubsystem(caloana24);

  InputManagers();
  Fun4AllInputManager *intrue2 = new Fun4AllRunNodeInputManager("DST_GEO");
  std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
  intrue2->AddFile(geoLocation);
  se->registerInputManager(intrue2);

  se->run(nEvents);
  //se->run(20);
  CDBInterface::instance()->Print(); // print used DB files
  se->End();
  cout << "JOB COMPLETE." << endl;
  se->PrintTimer();
  gSystem->Exit(0);
}
