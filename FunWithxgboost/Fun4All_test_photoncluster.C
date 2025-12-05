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
#include <caloreco/PhotonClusterBuilder.h>

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

#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>


void Fun4All_test_photoncluster(
    const int nEvents = 10,
    const string &inputFile0 = "g4hits.list",
    const string &inputFile1 = "dst_calo_cluster.list",
    const string &inputFile3 = "dst_mbd_epd.list",
    const string &inputFile4 = "dst_truth_jet.list",
    const string &inputFile2 = "dst_truth.list",

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
    INPUTREADHITS::listfile[1] = inputFile1;
    //INPUTREADHITS::listfile[2] = inputFile2;
    INPUTREADHITS::listfile[3] = inputFile3;
    INPUTREADHITS::listfile[4] = inputFile4;
  
    InputInit();
    InputRegister();
  
    FlagHandler *flag = new FlagHandler();
    se->registerSubsystem(flag);
  
    Enable::DSTOUT = false;
    Enable::DSTOUT_COMPRESS = false;
    DstOut::OutputDir = outDSTdir;
    DstOut::OutputFile = outputDSTFile;

    InputManagers();
  
    MbdReco *mbdreco = new MbdReco();
    se->registerSubsystem(mbdreco);
  
    GlobalVertexReco *gblvertex = new GlobalVertexReco();
    gblvertex->Verbosity(0);
    se->registerSubsystem(gblvertex);


    Process_Calo_Calib();

    //photon cluster
    PhotonClusterBuilder *photonclusterbuilder = new PhotonClusterBuilder();
    photonclusterbuilder->Verbosity(1);
    photonclusterbuilder->set_do_bdt(true);
    photonclusterbuilder->set_bdt_model_file("/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_single_tmva.root");
    std::vector<std::string> bdt_feature_list = {"vertex_z", "cluster_eta", "e11_over_e33", "et1", "et2", "et3", "et4"};
    photonclusterbuilder->set_bdt_feature_list(bdt_feature_list);
    photonclusterbuilder->set_ET_threshold(0.00);
    photonclusterbuilder->set_shower_shape_min_tower_energy(0.00);
    se->registerSubsystem(photonclusterbuilder);

    se->run(nEvents);
    se->End();

    delete se;
    std::cout << "All done!" << std::endl;
    gSystem->Exit(0);
}