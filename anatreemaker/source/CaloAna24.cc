

#include "CaloAna24.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

// Fun4All
#include <ffaobjects/EventHeader.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// ROOT stuff
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLorentzVector.h>
#include <TTree.h>

// For clusters and geometry
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

// Tower stuff
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

// GL1 Information
#include <ffarawobjects/Gl1Packet.h>

// for cluster vertex correction
#include <CLHEP/Geometry/Point3D.h>

// for the vertex
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Shower.h>
// caloEvalStack for cluster to truth matching
#include <g4eval/CaloEvalStack.h>
#include <g4eval/CaloRawClusterEval.h>

#include <mbd/MbdGeom.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <TLorentzVector.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h> // for GenVertex, GenVertex::part...
#pragma GCC diagnostic pop
//____________________________________________________________________________..
CaloAna24::CaloAna24(const std::string &name) : SubsysReco(name)
{
  std::cout << "CaloAna24::CaloAna24(const std::string &name) Calling ctor"
            << std::endl;
}

//____________________________________________________________________________..
CaloAna24::~CaloAna24()
{
  std::cout << "CaloAna24::~CaloAna24() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int CaloAna24::Init(PHCompositeNode *topNode)
{
  onnxmodule = onnxSession(m_modelPath);
  fout = new TFile("caloana.root", "RECREATE");

  slimtree = new TTree("slimtree", "slimtree");
  slimtree->Branch("mbdnorthhit", &mbdnorthhit, "mbdnorthhit/I");
  slimtree->Branch("mbdsouthhit", &mbdsouthhit, "mbdsouthhit/I");
  slimtree->Branch("vertexz", &vertexz, "vertexz/F");
  slimtree->Branch("vertexz_truth", &vertexz_truth, "vertexz_truth/F");
  slimtree->Branch("pythiaid", &m_pythiaid, "pythiaid/I");

  // particle level
  slimtree->Branch("nparticles", &nparticles, "nparticles/I");
  slimtree->Branch("particle_E", particle_E, "particle_E[nparticles]/F");
  slimtree->Branch("particle_Pt", particle_Pt, "particle_Pt[nparticles]/F");
  slimtree->Branch("particle_Eta", particle_Eta, "particle_Eta[nparticles]/F");
  slimtree->Branch("particle_Phi", particle_Phi, "particle_Phi[nparticles]/F");
  slimtree->Branch("particle_pid", particle_pid, "particle_pid[nparticles]/I");
  slimtree->Branch("particle_trkid", particle_trkid, "particle_trkid[nparticles]/I");
  slimtree->Branch("particle_photonclass", particle_photonclass, "particle_photonclass[nparticles]/I");
  slimtree->Branch("particle_photon_mother_pid", particle_photon_mother_pid, "particle_photon_mother_pid[nparticles]/I");
  slimtree->Branch("particle_truth_iso_02", particle_truth_iso_02, "particle_truth_iso_02[nparticles]/F");
  slimtree->Branch("particle_truth_iso_03", particle_truth_iso_03, "particle_truth_iso_03[nparticles]/F");
  slimtree->Branch("particle_truth_iso_04", particle_truth_iso_04, "particle_truth_iso_04[nparticles]/F");
  slimtree->Branch("particle_converted", particle_converted, "particle_converted[nparticles]/I");

  //daughter level
  slimtree->Branch("ndaughter", ndaughter, "ndaughter/I");
  slimtree->Branch("daughter_pid", daughter_pid, "daughter_pid[ndaughter]/I");
  slimtree->Branch("daughter_E", daughter_E, "daughter_E[ndaughter]/F");
  slimtree->Branch("daughter_Pt", daughter_Pt, "daughter_Pt[ndaughter]/F");
  slimtree->Branch("daughter_Eta", daughter_Eta, "daughter_Eta[ndaughter]/F");
  slimtree->Branch("daughter_Phi", daughter_Phi, "daughter_Phi[ndaughter]/F");

  for (int i = 0; i < nclustercontainer; i++)
  {
    slimtree->Branch(Form("ncluster_%s", clusternamelist[i].c_str()), &ncluster[i], Form("ncluster_%s/I", clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_E_%s", clusternamelist[i].c_str()), cluster_E[i], Form("cluster_E_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_Et_%s", clusternamelist[i].c_str()), cluster_Et[i], Form("cluster_Et_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_Eta_%s", clusternamelist[i].c_str()), cluster_Eta[i], Form("cluster_Eta_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_Phi_%s", clusternamelist[i].c_str()), cluster_Phi[i], Form("cluster_Phi_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_prob_%s", clusternamelist[i].c_str()), cluster_prob[i], Form("cluster_prob_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_CNN_prob_%s", clusternamelist[i].c_str()), cluster_CNN_prob[i], Form("cluster_CNN_prob_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_truthtrkID_%s", clusternamelist[i].c_str()), cluster_truthtrkID[i], Form("cluster_truthtrkID_%s[ncluster_%s]/I", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_pid_%s", clusternamelist[i].c_str()), cluster_pid[i], Form("cluster_pid_%s[ncluster_%s]/I", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_iso_02_%s", clusternamelist[i].c_str()), cluster_iso_02[i], Form("cluster_iso_02_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_iso_03_%s", clusternamelist[i].c_str()), cluster_iso_03[i], Form("cluster_iso_03_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_iso_04_%s", clusternamelist[i].c_str()), cluster_iso_04[i], Form("cluster_iso_04_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_iso_04_emcal_%s", clusternamelist[i].c_str()), cluster_iso_04_emcal[i], Form("cluster_iso_04_emcal_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_iso_04_hcalin_%s", clusternamelist[i].c_str()), cluster_iso_04_hcalin[i], Form("cluster_iso_04_hcalin_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_iso_04_hcalout_%s", clusternamelist[i].c_str()), cluster_iso_04_hcalout[i], Form("cluster_iso_04_hcalout_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e1_%s", clusternamelist[i].c_str()), cluster_e1[i], Form("cluster_e1_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e2_%s", clusternamelist[i].c_str()), cluster_e2[i], Form("cluster_e2_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e3_%s", clusternamelist[i].c_str()), cluster_e3[i], Form("cluster_e3_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e4_%s", clusternamelist[i].c_str()), cluster_e4[i], Form("cluster_e4_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_et1_%s", clusternamelist[i].c_str()), cluster_et1[i], Form("cluster_et1_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_et2_%s", clusternamelist[i].c_str()), cluster_et2[i], Form("cluster_et2_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_et3_%s", clusternamelist[i].c_str()), cluster_et3[i], Form("cluster_et3_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_et4_%s", clusternamelist[i].c_str()), cluster_et4[i], Form("cluster_et4_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_weta_%s", clusternamelist[i].c_str()), cluster_weta[i], Form("cluster_weta_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_wphi_%s", clusternamelist[i].c_str()), cluster_wphi[i], Form("cluster_wphi_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_ietacent_%s", clusternamelist[i].c_str()), cluster_ietacent[i], Form("cluster_ietacent_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_iphicent_%s", clusternamelist[i].c_str()), cluster_iphicent[i], Form("cluster_iphicent_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_detamax_%s", clusternamelist[i].c_str()), cluster_detamax[i], Form("cluster_detamax_%s[ncluster_%s]/I", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_dphimax_%s", clusternamelist[i].c_str()), cluster_dphimax[i], Form("cluster_dphimax_%s[ncluster_%s]/I", clusternamelist[i].c_str(), clusternamelist[i].c_str()));

    slimtree->Branch(Form("cluster_e11_%s", clusternamelist[i].c_str()), cluster_e11[i], Form("cluster_e11_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e22_%s", clusternamelist[i].c_str()), cluster_e22[i], Form("cluster_e22_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e13_%s", clusternamelist[i].c_str()), cluster_e13[i], Form("cluster_e13_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e15_%s", clusternamelist[i].c_str()), cluster_e15[i], Form("cluster_e15_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e17_%s", clusternamelist[i].c_str()), cluster_e17[i], Form("cluster_e17_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e31_%s", clusternamelist[i].c_str()), cluster_e31[i], Form("cluster_e31_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e51_%s", clusternamelist[i].c_str()), cluster_e51[i], Form("cluster_e51_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e71_%s", clusternamelist[i].c_str()), cluster_e71[i], Form("cluster_e71_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e33_%s", clusternamelist[i].c_str()), cluster_e33[i], Form("cluster_e33_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e35_%s", clusternamelist[i].c_str()), cluster_e35[i], Form("cluster_e35_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e37_%s", clusternamelist[i].c_str()), cluster_e37[i], Form("cluster_e37_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e53_%s", clusternamelist[i].c_str()), cluster_e53[i], Form("cluster_e53_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e73_%s", clusternamelist[i].c_str()), cluster_e73[i], Form("cluster_e73_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e55_%s", clusternamelist[i].c_str()), cluster_e55[i], Form("cluster_e55_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e57_%s", clusternamelist[i].c_str()), cluster_e57[i], Form("cluster_e57_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e75_%s", clusternamelist[i].c_str()), cluster_e75[i], Form("cluster_e75_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e77_%s", clusternamelist[i].c_str()), cluster_e77[i], Form("cluster_e77_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_w32_%s", clusternamelist[i].c_str()), cluster_w32[i], Form("cluster_w32_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e32_%s", clusternamelist[i].c_str()), cluster_e32[i], Form("cluster_e32_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_w72_%s", clusternamelist[i].c_str()), cluster_w72[i], Form("cluster_w72_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e72_%s", clusternamelist[i].c_str()), cluster_e72[i], Form("cluster_e72_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));

    slimtree->Branch(Form("cluster_ihcal_et_%s", clusternamelist[i].c_str()), cluster_ihcal_et[i], Form("cluster_ihcal_et_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_ohcal_et_%s", clusternamelist[i].c_str()), cluster_ohcal_et[i], Form("cluster_ohcal_et_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_ihcal_et22_%s", clusternamelist[i].c_str()), cluster_ihcal_et22[i], Form("cluster_ihcal_et22_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_ohcal_et22_%s", clusternamelist[i].c_str()), cluster_ohcal_et22[i], Form("cluster_ohcal_et22_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_ihcal_et33_%s", clusternamelist[i].c_str()), cluster_ihcal_et33[i], Form("cluster_ihcal_et33_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_ohcal_et33_%s", clusternamelist[i].c_str()), cluster_ohcal_et33[i], Form("cluster_ohcal_et33_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_ihcal_ieta_%s", clusternamelist[i].c_str()), cluster_ihcal_ieta[i], Form("cluster_ihcal_ieta_%s[ncluster_%s]/I", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_ihcal_iphi_%s", clusternamelist[i].c_str()), cluster_ihcal_iphi[i], Form("cluster_ihcal_iphi_%s[ncluster_%s]/I", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_ohcal_ieta_%s", clusternamelist[i].c_str()), cluster_ohcal_ieta[i], Form("cluster_ohcal_ieta_%s[ncluster_%s]/I", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_ohcal_iphi_%s", clusternamelist[i].c_str()), cluster_ohcal_iphi[i], Form("cluster_ohcal_iphi_%s[ncluster_%s]/I", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
  }

  h_tracking_radiograph = new TH3F("tracking_radiograph", "tracking_radiograph", 200, -100, 100, 200, -100, 100, 400, -200, 200);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloAna24::InitRun(PHCompositeNode *topNode)
{
  std::cout
      << "CaloAna24::InitRun(PHCompositeNode *topNode) Initializing for Run XXX"
      << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloAna24::process_event(PHCompositeNode *topNode)
{
  // trigger ana for data
  if (!isMC)
  {
    Gl1Packet *gl1PacketInfo =
        findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (!gl1PacketInfo)
    {
      std::cout << PHWHERE << "caloTreeGen::process_event: "
                << "Gl1Packet"
                << " node is missing. Output related to this node will be empty"
                << std::endl;
    }

    if (gl1PacketInfo)
    {
      uint64_t triggervec = gl1PacketInfo->getScaledVector();
      uint64_t triggervecraw = gl1PacketInfo->getLiveVector();
      for (int i = 0; i < 64; i++)
      {
        bool trig_decision = ((triggervec & 0x1U) == 0x1U);
        bool trig_decision_raw = ((triggervecraw & 0x1U) == 0x1U);

        if (i < 32)
        {
          // reset it just to be safe
          scaledtrigger[i] = false;
          scaledtrigger[i] = trig_decision;
          m_scaledtrigger[i] = trig_decision;
          // std::cout<<"Scaled trigger: "<<i<<" "<<trig_decision<<std::endl;
          livetrigger[i] = false;
          livetrigger[i] = trig_decision_raw;
          if (trig_decision)
            nscaledtrigger[i]++;
          if (trig_decision_raw)
            nlivetrigger[i]++;
          if (!initilized)
          {
            for (int j = 0; j < 3; j++)
            {
              initscaler[i][j] = gl1PacketInfo->lValue(i, j);
            }
            initilized = true;
          }
          for (int j = 0; j < 3; j++)
          {
            currentscaler[i][j] = gl1PacketInfo->lValue(i, j);
          }
        }

        triggervec = (triggervec >> 1U) & 0xffffffffU;
        triggervecraw = (triggervecraw >> 1U) & 0xffffffffU;
      }
    }
  }

  if (isMC && !isSingleParticle)
  {
    // hepmc record
    PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
    int ngenevents = 0;
    bool ispromptphoton = false;
    for (PHHepMCGenEventMap::Iter iter = genevtmap->begin(); iter != genevtmap->end(); ++iter)
    {
      PHHepMCGenEvent *genevt = iter->second;
      std::cout << "event embedded: " << genevt->get_embedding_id() << std::endl;
      HepMC::GenEvent *event = genevt->getEvent();
      // need to add the part to fish out the signal event
      singal_event = event;
      if (!event)
      {
        std::cout << PHWHERE << " no evt pointer under HEPMC Node found" << std::endl;
      }
      int process_id = event->signal_process_id();
      m_pythiaid = process_id;
      if (process_id > 200)
        ispromptphoton = true;
      // event->print();
      // std::cout<<"process_id: "<<process_id<<std::endl;
      ngenevents++;
    }
    // this is for photon jet sample, only look at the prompt photon process
    std::cout << "ispromptphoton: " << ispromptphoton << std::endl;
    // if(ispromptphoton) return Fun4AllReturnCodes::EVENT_OK;
    std::cout << "ngenevents: " << ngenevents << std::endl;

    // mbd trigger

    MbdPmtContainer *mbdtow = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");

    // bool mbdevent = false;
    if (mbdtow)
    {
      int northhit = 0;
      int southhit = 0;
      int sectormb = 128; // mbdtow->get_npmt();
      float mbenrgy[128] = {0};
      // if(_debug) cout << "Got " << sectormb << " mbd sectors in sim." << endl;
      for (int i = 0; i < sectormb; ++i)
      {
        MbdPmtHit *mbdhit = mbdtow->get_pmt(i);

        mbenrgy[i] = mbdhit->get_q();
        // std::cout<<mbenrgy[i]<<std::endl;
        if (mbenrgy[i] > 0.4 && i < 64)
          northhit += 1;
        if (mbenrgy[i] > 0.4 && i > 63)
          southhit += 1;
      }
      mbdnorthhit = northhit;
      mbdsouthhit = southhit;
    }
  }

  float m_vertex = -9999;
  std::vector<TLorentzVector> goodcluster;
  // GlobalVertexMap *vertexmap =
  //     findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!isSingleParticle)
  {
    MbdVertexMap *vertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");

    if (!vertexmap)
    {
      std::cout << "GlobalVertexMap node is missing" << std::endl;
    }
    if (vertexmap && !vertexmap->empty())
    {
      MbdVertex *vtx = vertexmap->begin()->second;
      if (vtx)
      {
        m_vertex = vtx->get_z();
        // std::cout<<m_vertex<<std::endl;
        //  if nan return
        /*
        if (m_vertex != m_vertex)
          return Fun4AllReturnCodes::EVENT_OK;
        if (abs(m_vertex) > 300)
          return Fun4AllReturnCodes::EVENT_OK;
        */
      }
      else
      {
        return Fun4AllReturnCodes::EVENT_OK;
      }
    }
  }

  vertexz = m_vertex;
  std::cout << "vertexz: " << vertexz << std::endl;
  // set of primary particles
  std::set<PHG4Particle *> primary_particles;
  std::set<PHG4Particle *> primary_photon_candidates;
  std::set<PHG4Particle *> photonsfrompi0;
  std::set<PHG4Particle *> photonsfrometa;
  std::set<PHG4Particle *> badphotons;
  std::set<PHG4Particle *> convertedphotons;
  std::map<PHG4Particle *, std::vector<float>> photontruthiso;
  if (isMC)
  {

    caloevalstack = new CaloEvalStack(topNode, "CEMC");
    clustereval = caloevalstack->get_rawcluster_eval();
    clustereval->set_usetowerinfo(true);
    clustereval->next_event(topNode);
    trutheval = caloevalstack->get_truth_eval();
    // CaloRawTowerEval *towereval = caloevalstack.get_rawtower_eval();
    //  clustereval->next_event(topNode);
    // CaloTruthEval *trutheval = m_caloevalstack->get_truth_eval();
    std::cout << "trutheval: " << trutheval << " clustereval: " << clustereval << std::endl;
    // trutheval->next_event(topNode);

    truthinfo =
        findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (!truthinfo)
    {
      std::cout << PHWHERE
                << "PHG4TruthInfoContainer node is missing, can't collect "
                   "some true information"
                << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    // std::cout<<"truthinfo: "<<truthinfo<<std::endl;

    int primaryvtxid = truthinfo->GetPrimaryVertexIndex();
    PHG4VtxPoint *primaryvtx = truthinfo->GetVtx(primaryvtxid);
    if (!primaryvtx)
    {
      std::cout << "primaryvtx is missing" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    vertexz_truth = primaryvtx->get_z();

    if (isSingleParticle)
    {
      vertexz = vertexz_truth;
      m_vertex = vertexz_truth;
    }

    //  loop over truth primary particles
    PHG4TruthInfoContainer::ConstRange range =
        truthinfo->GetPrimaryParticleRange();

    for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
         truth_itr != range.second; ++truth_itr)
    {
      PHG4Particle *truth = truth_itr->second;
      // std::cout<<trutheval->get_embed(truth)<<std::endl;
      if (trutheval->get_embed(truth) < 1)
        continue;
      primary_particles.insert(truth);
      int pid = truth->get_pid();
      // std::cout<<"truth pid: "<<pid<<std::endl;
      if (pid == 22)
      {
        primary_photon_candidates.insert(truth);
      }
    }

    std::cout << "primary_photon_candidates.size(): " << primary_photon_candidates.size() << std::endl;
    for (auto truth_photon : primary_photon_candidates)
    {
      // int trackid1 = truth_photon->get_track_id();
      //  std::cout<<"trackid1: "<<trackid1<<std::endl;
      //  check if the photon is from pi0 or eta
      if (photonsfrompi0.find(truth_photon) != photonsfrompi0.end())
        continue;
      if (photonsfrometa.find(truth_photon) != photonsfrometa.end())
        continue;
      // if we get the HEPMC tracking working we don't need the inv mass check

      TLorentzVector photon1;
      photon1.SetPxPyPzE(truth_photon->get_px(), truth_photon->get_py(), truth_photon->get_pz(), truth_photon->get_e());
      /*
      for (auto truth_photon_2 : primary_photon_candidates)
      {
        int trackid2 = truth_photon_2->get_track_id();
        // std::cout<<"trackid2: "<<trackid2<<std::endl;
        if (truth_photon == truth_photon_2)
          continue;
        if (trackid1 == trackid2)
          continue;
        TLorentzVector photon2;
        photon2.SetPxPyPzE(truth_photon_2->get_px(), truth_photon_2->get_py(), truth_photon_2->get_pz(), truth_photon_2->get_e());

        TLorentzVector diphoton = photon1 + photon2;
        float mass = diphoton.M();
        // std::cout<<"mass: "<<mass<<std::endl;
        // print the four momentum of the diphoton
        float pi0mass = 0.1349799931; // get from the print back
        float etamass = 0.5478500128;

        if (abs(mass - pi0mass) < 1E-8)
        {
          photonsfrompi0.insert(truth_photon);
          photonsfrompi0.insert(truth_photon_2);
          // h_pi0ET->Fill(diphoton.Et());
          // print with more precision
          // std::cout.precision(10);
          // std::cout<<"mass: "<<mass<<std::endl;
          break;
        }
        if (abs(mass - etamass) < 1E-8)
        {
          photonsfrometa.insert(truth_photon);
          photonsfrometa.insert(truth_photon_2);
          // std::cout.precision(10);
          // std::cout<<"mass: "<<mass<<std::endl;
          break;
        }
      }
      */
      // this part check for pair conversion
      int trackid = truth_photon->get_track_id();
      if (abs(photon1.Eta()) > 1.5)
        continue;
      if (photon1.Et() < 5)
        continue;
      // get PHG4Shower
      PHG4Shower *shower = truthinfo->GetShower(trackid);
      // std::cout<<"shower: "<<shower<<std::endl;
      if (!shower)
        continue;
      // loop over particles in the shower
      auto g4particle_ids = shower->g4particle_ids();
      for (auto g4particle_id : g4particle_ids)
      {
        PHG4Particle *g4particle = truthinfo->GetParticle(g4particle_id);
        if (!g4particle)
          continue;
        int vertexid = g4particle->get_vtx_id();
        PHG4VtxPoint *vtx = truthinfo->GetVtx(vertexid);
        if (!vtx)
          continue;
        float vertexr = sqrt(vtx->get_x() * vtx->get_x() + vtx->get_y() * vtx->get_y());
        if (vertexr < 93)
        {
          float momentum = sqrt(g4particle->get_px() * g4particle->get_px() + g4particle->get_py() * g4particle->get_py() + g4particle->get_pz() * g4particle->get_pz());
          if (momentum > 0.4 * photon1.E())
          {
            int g4particlepid = g4particle->get_pid();
            badphotons.insert(truth_photon);
            if (abs(g4particlepid) == 11)
              convertedphotons.insert(truth_photon);
          }
          h_tracking_radiograph->Fill(vtx->get_x(), vtx->get_y(), vtx->get_z(), momentum);
        }
      }
    }
    nparticles = 0;
    float merger = 0.001;
    float isor[3] = {0.2, 0.3, 0.4};
    // calculate truth iso
    for (auto truth : primary_particles)
    {
      // std::cout<<truth<<std::endl;
      float isoET[3] = {0, 0, 0};
      float clusterET = 0;
      int pid = truth->get_pid();
      int trackid = truth->get_track_id();

      TLorentzVector p1 = TLorentzVector(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
      // skip for soft stuff
      if (p1.E() < particlepTmin)
        continue;
      // cut on eta
      if (abs(p1.Eta()) > 1.5)
        continue;
      int barcode = truth->get_barcode();

      bool verbosephoton = false;
      // bool ispromptphoton = false;
      int photonclass = 0;
      int photonmotherpid = 0;
      if (pid == 22)
      {
        verbosephoton = true;
        // ispromptphoton = true;
        /*
        if (photonsfrompi0.find(truth) != photonsfrompi0.end())
        {
          verbosephoton = false;
          ispromptphoton = false;
        }
        if (photonsfrometa.find(truth) != photonsfrometa.end())
        {
          verbosephoton = false;
          ispromptphoton = false;
        }
        */
        if (!isSingleParticle)
        {
          photonclass = photon_type(barcode).first;
          photonmotherpid = photon_type(barcode).second;
        }
        if (photonclass <= 2)
          verbosephoton = true;
      }
      if (verbosephoton)
      {
        std::cout << "truth photon pid: " << pid << " eta: " << p1.Eta() << " phi: " << p1.Phi() << "ET: " << p1.Et() << std::endl;
      }
      for (auto truth2 : primary_particles)
      {
        // want to correlate with the same particle for clusterET
        // if(truth == truth2) continue;
        TLorentzVector p2 = TLorentzVector(truth2->get_px(), truth2->get_py(), truth2->get_pz(), truth2->get_e());
        /*
        int p2pid = truth2->get_pid();
        //if not eta, pi0 photon neuron anti-neutron, apply a 0.5 pT cut

        if (p2pid != 22 && p2pid != 111 && p2pid != -2112 && p2pid != 2112 ) {
          if(p2.Pt() < 1) continue;
        }
        */
        float dr = p1.DeltaR(p2);
        for (int i = 0; i < 3; i++)
        {
          if (dr < isor[i])
          {
            isoET[i] += p2.Et();
          }
        }

        if (dr < merger)
        {
          clusterET += p2.Et();
        }
      }
      for (int i = 0; i < 3; i++)
      {
        isoET[i] -= clusterET;
      }

      // bool isisophoton = false;
      // if photon, electron pi0, positron, eta
      // if (pid == 22 || pid == 11 || pid == -11 || pid == 111 || pid == 221) {
      if (pid != 0)
      {

        if (photonsfrompi0.find(truth) != photonsfrompi0.end())
        {
          pid = 111;
        }
        if (photonsfrometa.find(truth) != photonsfrometa.end())
        {
          pid = 221;
        }
        if (p1.Pt() > particlepTmin)
        {
          int converted = 0;
          if (badphotons.find(truth) != badphotons.end())
          {
            converted = 2;
            if (convertedphotons.find(truth) != convertedphotons.end())
              converted = 1;
          }
          // std::cout<<"nparticles: "<<nparticles<<std::endl;
          // std::cout<<"nparticles: "<<nparticles<<" E: "<<p1.E()<<" Pt: "<<p1.Pt()<<" Eta: "<<p1.Eta()<<" Phi: "<<p1.Phi()<<" pid: "<<pid<<" trackid: "<<trackid<<" converted: "<<converted<<" isoET: "<<isoET[0]<<" "<<isoET[1]<<" "<<isoET[2]<<std::endl;
          particle_E[nparticles] = p1.E();
          particle_Pt[nparticles] = p1.Pt();
          particle_Eta[nparticles] = p1.Eta();
          particle_Phi[nparticles] = p1.Phi();
          particle_pid[nparticles] = pid;
          particle_trkid[nparticles] = trackid;
          particle_photonclass[nparticles] = photonclass;
          particle_photon_mother_pid[nparticles] = photonmotherpid;
          particle_converted[nparticles] = converted;
          particle_truth_iso_02[nparticles] = isoET[0];
          particle_truth_iso_03[nparticles] = isoET[1];
          particle_truth_iso_04[nparticles] = isoET[2];
          nparticles++;
          if (nparticles > nparticlesmax)
          {
            std::cout << "nparticles exceed the max range: " << nparticles << std::endl;
            return Fun4AllReturnCodes::ABORTEVENT;
          }
        }
      }
    }
  }

  //only for single particles, get the first daughters from the primary 
  if (isSingleParticle)
  {
    ndaughter = 0;
    PHG4TruthInfoContainer::ConstRange range =
        truthinfo->GetSecondaryParticleRange();
    for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
         truth_itr != range.second; ++truth_itr)
    {
      PHG4Particle *truth = truth_itr->second;
      int pid = truth->get_pid();
      int parentid = truth->get_parent_id();
      if(parentid < 0) continue;
      TLorentzVector p1 = TLorentzVector(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());

      daughter_E[ndaughter] = p1.E();
      daughter_Pt[ndaughter] = p1.Pt();
      daughter_Eta[ndaughter] = p1.Eta();
      daughter_Phi[ndaughter] = p1.Phi();
      daughter_pid[ndaughter] = pid;
      daughter_parent_trackid[ndaughter] = parentid;
      ndaughter++;
      if (ndaughter > ndaughtermax)
      {
        std::cout << "ndaughter exceed the max range: " << ndaughter << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }

    }
  }

  // geom nodes:
  geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  int IHCALsize = 1;
  for (int i = 0; i < IHCALsize; i++)
  {

    const RawTowerDefs::keytype key = TowerInfoDefs::get_hcalin_geokey_at_channel(i);
    float tower_phi = geomIH->get_tower_geometry(key)->get_phi();
    std::cout << "tower_phi: " << tower_phi << std::endl;
  }

  std::string towerNodeName = "TOWERINFO_CALIB_CEMC";
  emcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, towerNodeName);
  if (!emcTowerContainer)
  {
    std::cout << "RawClusterCNNClassifier::process_event Could not locate tower node " << towerNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  std::string ihcalTowerNodeName = "TOWERINFO_CALIB_HCALIN";
  ihcalTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, ihcalTowerNodeName);
  if (!ihcalTowerContainer)
  {
    std::cout << "RawClusterCNNClassifier::process_event Could not locate tower node " << ihcalTowerNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  std::string ohcalTowerNodeName = "TOWERINFO_CALIB_HCALOUT";
  ohcalTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, ohcalTowerNodeName);
  if (!ohcalTowerContainer)
  {
    std::cout << "RawClusterCNNClassifier::process_event Could not locate tower node " << ohcalTowerNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (!geomEM || !geomIH || !geomOH)
  {
    std::cout << PHWHERE
              << "CaloAna24::process_event - missing tower geometry node"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  std::cout << "size of photonsfrompi0: " << photonsfrompi0.size() << std::endl;
  for (int i = 0; i < nclustercontainer; i++)
  {
    RawClusterContainer *clusterContainer =
        findNode::getClass<RawClusterContainer>(topNode,
                                                clusternamelist[i]);
    if (!clusterContainer)
    {
      std::cout << PHWHERE
                << "CaloAna24::process_event - missing cluster node: "
                << clusternamelist[i] << std::endl;

      return Fun4AllReturnCodes::ABORTEVENT;
    }
    ncluster[i] = 0;

    RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
    RawClusterContainer::ConstIterator clusterIter;
    float maxclusterpt = -1;

    for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; ++clusterIter)
    {
      RawCluster *recoCluster = clusterIter->second;
      // std::cout<<"recoCluster: "<<recoCluster<<std::endl;
      float prob = recoCluster->get_prob();

      CLHEP::Hep3Vector vertex(0, 0, m_vertex);

      CLHEP::Hep3Vector E_vec_cluster_Full = RawClusterUtility::GetEVec(*recoCluster, vertex);

      float ecalib = 1.00;
      float ET = E_vec_cluster_Full.perp() * ecalib;
      float E = E_vec_cluster_Full.mag();
      float phi = E_vec_cluster_Full.phi();
      float eta = E_vec_cluster_Full.eta();
      // std::cout<<E<<" "<<ET<<" "<<eta<<" "<<phi<<std::endl;
      // if (abs(eta) > 1.0)
      // continue;
      // if (ET > 1) h_ET->Fill(ET);
      if (ET < clusterpTmin)
        continue;

      // Array for storing the isolation energy for different radii

        float clusteriso[nRadii];

        // Loop to calculate the isolation energy for each radius
        // std::cout<<"nRadii: "<<nRadii<<std::endl;
        for (int i = 0; i < nRadii; ++i)
        {
          clusteriso[i] = recoCluster->get_et_iso(2 + i, false, true);
        }
      

      float emcalET_04 = calculateET(eta, phi, 0.4, 0);
      float ihcalET_04 = calculateET(eta, phi, 0.4, 1);
      float ohcalET_04 = calculateET(eta, phi, 0.4, 2);

      if (ET > maxclusterpt)
      {
        maxclusterpt = ET;
      }
      int trackid = -1;
      float clusterE = E_vec_cluster_Full.mag();
      int pid = 0;
      if (isMC)
      {
        // std::cout << clustereval << std::endl;
        // trutheval = caloevalstack->get_truth_eval();
        // PHG4Shower* max_shower = clustereval->max_truth_primary_shower_by_energy(recoCluster);
        // std::cout<<"max_shower: "<<max_shower<<std::endl;
        // int parentid = max_shower->get_parent_particle_id();
        // std::cout<<"parentid: "<<parentid<<std::endl;
        // PHG4Particle* showerparticle = truthinfo->GetParticle(parentid);
        // std::cout<<"showerparticle: "<<showerparticle<<std::endl;
        PHG4Particle *maxPrimary = clustereval->max_truth_primary_particle_by_energy(recoCluster);
        // if (ET > 1)
        std::cout << "maxPrimary: " << maxPrimary << std::endl;

        if (!maxPrimary)
          continue;
        if (trutheval->get_embed(maxPrimary) != 1)
          continue;
        pid = maxPrimary->get_pid();
        trackid = maxPrimary->get_track_id();
        float bestprimaryenergy = maxPrimary->get_e();
        float bestprimaryincluster = clustereval->get_energy_contribution(recoCluster, maxPrimary);

        // if from pi0 or eta
        if (photonsfrompi0.find(maxPrimary) != photonsfrompi0.end())
          pid = 111;
        if (photonsfrometa.find(maxPrimary) != photonsfrometa.end())
          pid = 221;

        if (ET > 5)
          std::cout << "pid: " << pid << " bestprimaryenergy: " << bestprimaryenergy
                    << " bestprimaryincluster: " << bestprimaryincluster << " clusterE: " << clusterE << std::endl;
      }

      std::vector<float> showershape = recoCluster->get_shower_shapes(0.070);
      std::pair<int, int> leadtowerindex = recoCluster->get_lead_tower();
      //-------------------------------------------------------------------------------------
      // filling the 7x7 matrix
      std::cout << "finding showershapes in 7x7" << std::endl;
      int maxieta = leadtowerindex.first;
      int maxiphi = leadtowerindex.second;

      int maxtowerieta = maxieta;
      int maxtoweriphi = maxiphi;
      float CNNprob = -1;
      std::vector<float> input;
      const int inputDimx = 5;
      const int inputDimy = 5;
      const int inputDimz = 1;
      const int outputDim = 1;
      // resize to inputDimx * inputDimy
      int vectorSize = inputDimx * inputDimy;
      input.resize(vectorSize, 0);
      // loop for classification
      if (ET > 0)
      {
        int xlength = int((inputDimx - 1) / 2);
        int ylength = int((inputDimy - 1) / 2);
        if (maxtowerieta - ylength < 0 || maxtowerieta + ylength >= 96)
        {
          continue;
        }
        for (int ieta = maxtowerieta - ylength; ieta <= maxtowerieta + ylength; ieta++)
        {
          for (int iphi = maxtoweriphi - xlength; iphi <= maxtoweriphi + xlength; iphi++)
          {
            int mappediphi = iphi;

            if (mappediphi < 0)
            {
              mappediphi += 256;
            }
            if (mappediphi > 255)
            {
              mappediphi -= 256;
            }
            unsigned int towerinfokey = TowerInfoDefs::encode_emcal(ieta, mappediphi);
            TowerInfo *towerinfo = emcTowerContainer->get_tower_at_key(towerinfokey);
            if (!towerinfo)
            {
              // should not happen
              std::cout << "No towerinfo for tower key " << towerinfokey << std::endl;
              std::cout << "ieta: " << ieta << " iphi: " << mappediphi << std::endl;
              continue;
            }
            int index = (ieta - maxtowerieta + ylength) * inputDimx + iphi - maxtoweriphi + xlength;
            input.at(index) = towerinfo->get_energy();
          }
        }
      }
      std::vector<float> probresult = onnxInference(onnxmodule, input, 1, inputDimx, inputDimy, inputDimz, outputDim);

      CNNprob = probresult[0];
      // to find the
      float avg_eta = showershape[4] + 0.5;
      float avg_phi = showershape[5] + 0.5;
      // don't use max tower use the center of the cluster
      maxieta = std::floor(avg_eta);
      maxiphi = std::floor(avg_phi);

      // here skip if we are too close to the edge +-3
      if (maxieta < 3 || maxieta > 92)
        continue;
      float E77[7][7] = {0};
      // loop over 7 by 7 around it
      for (int ieta = maxieta - 3; ieta < maxieta + 4; ieta++)
      {
        for (int iphi = maxiphi - 3; iphi < maxiphi + 4; iphi++)
        {
          int temp_ieta = ieta;
          int temp_iphi = iphi;
          shift_tower_index(ieta, iphi, 96, 256);

          unsigned int towerinfokey = TowerInfoDefs::encode_emcal(ieta, iphi);
          ieta = temp_ieta;
          iphi = temp_iphi;
          if (ieta < 0 || ieta > 95)
            continue;
          TowerInfo *towerinfo = emcTowerContainer->get_tower_at_key(towerinfokey);

          if (!towerinfo)
          {

            // should not happen
            std::cout << "No towerinfo for tower key " << towerinfokey << std::endl;
            std::cout << "ieta: " << ieta << " iphi: " << iphi << std::endl;
            continue;
          }

          E77[ieta - maxieta + 3][iphi - maxiphi + 3] = towerinfo->get_energy();
        }
      }
      //-------------------------------------------------------------------------------------
      float e11 = E77[3][3];
      float e22 = showershape[8] + showershape[9] + showershape[10] + showershape[11];
      float e13 = 0;
      float e15 = 0;
      float e17 = 0;
      float e31 = 0;
      float e51 = 0;
      float e71 = 0;

      float e33 = 0;
      float e35 = 0;
      float e37 = 0;
      float e53 = 0;
      float e73 = 0;

      float e55 = 0;
      float e57 = 0;
      float e75 = 0;
      float e77 = 0;
      // here also need to calculate second moment in eta for e32 and e72
      float w32 = 0;
      float e32 = 0;
      float w72 = 0;
      float e72 = 0;

      int signphi = (avg_phi - std::floor(avg_phi)) > 0.5 ? 1 : -1;

      for (int i = 0; i < 7; i++)
      {
        for (int j = 0; j < 7; j++)
        {
          int di = abs(i - 3);
          int dj = abs(j - 3);

          e77 += E77[i][j];
          if (di <= 1 && (dj == 0 || j == (3 + signphi)))
          {
            w32 += E77[i][j] * (i - 3) * (i - 3);
            e32 += E77[i][j];
          }
          if (di <= 3 && (dj == 0 || j == (3 + signphi)))
          {
            w72 += E77[i][j] * (i - 3) * (i - 3);
            e72 += E77[i][j];
          }

          if (di <= 0 && dj <= 1)
          {
            e13 += E77[i][j];
          }
          if (di <= 0 && dj <= 2)
          {
            e15 += E77[i][j];
          }
          if (di <= 0 && dj <= 3)
          {
            e17 += E77[i][j];
          }
          if (di <= 1 && dj <= 0)
          {
            e31 += E77[i][j];
          }
          if (di <= 2 && dj <= 0)
          {
            e51 += E77[i][j];
          }
          if (di <= 3 && dj <= 0)
          {
            e71 += E77[i][j];
          }

          if (di <= 1 && dj <= 1)
          {
            e33 += E77[i][j];
          }
          if (di <= 1 && dj <= 2)
          {
            e35 += E77[i][j];
          }
          if (di <= 1 && dj <= 3)
          {
            e37 += E77[i][j];
          }
          if (di <= 2 && dj <= 1)
          {
            e53 += E77[i][j];
          }
          if (di <= 3 && dj <= 1)
          {
            e73 += E77[i][j];
          }

          if (di <= 2 && dj <= 2)
          {
            e55 += E77[i][j];
          }
          if (di <= 2 && dj <= 3)
          {
            e57 += E77[i][j];
          }
          if (di <= 3 && dj <= 2)
          {
            e75 += E77[i][j];
          }
        }
      }
      w32 = e32 > 0 ? w32 / e32 : 0;
      w72 = e72 > 0 ? w72 / e72 : 0;

      w32 = sqrt(w32);
      w72 = sqrt(w72);

      //-------------------------------------------------------------------------------------
      // for detamax and dphimax
      std::cout << "for detamax and dphimax" << std::endl;
      int detamax = 0;
      int dphimax = 0;
      // loop over all towers in the cluster
      const RawCluster::TowerMap tower_map =
          recoCluster->get_towermap();
      for (auto tower_iter : tower_map)
      {
        RawTowerDefs::keytype tower_key = tower_iter.first;
        // get ieta iphi
        float eta = RawTowerDefs::decode_index1(tower_key);
        float phi = RawTowerDefs::decode_index2(tower_key);
        int totalphibins = 256;
        auto dphiwrap = [totalphibins](float towerphi, float maxiphi)
        {
          float idphi = towerphi - maxiphi;
          float idphiwrap = totalphibins - std::abs(idphi);
          if (std::abs(idphiwrap) < std::abs(idphi))
          {
            return (idphi > 0) ? -idphiwrap : idphiwrap;
          }
          return idphi;
        };

        float deta = eta - maxieta;
        float dphi = dphiwrap(phi, maxiphi);

        if (abs(deta) > detamax)
          detamax = abs(deta);
        if (abs(dphi) > dphimax)
          dphimax = abs(dphi);
      }
      //-------------------------------------------------------------------------------------

      //-------------------------------------------------------------------------------------
      // finding the ihcal ohcal energy behind the cluster
      std::cout << "finding the ihcal ohcal energy behind the cluster" << std::endl;
      std::vector<int> ihcal_tower = find_closest_hcal_tower(eta, phi, geomIH, ihcalTowerContainer, vertexz, true);
      std::vector<int> ohcal_tower = find_closest_hcal_tower(eta, phi, geomOH, ohcalTowerContainer, vertexz, false);

      float ihcal_et = 0;
      float ohcal_et = 0;
      float ihcal_et22 = 0;
      float ohcal_et22 = 0;
      float ihcal_et33 = 0;
      float ohcal_et33 = 0;

      int ihcal_ieta = ihcal_tower[0];
      int ihcal_iphi = ihcal_tower[1];
      float ihcalEt33[3][3] = {0};

      int ohcal_ieta = ohcal_tower[0];
      int ohcal_iphi = ohcal_tower[1];
      float ohcalEt33[3][3] = {0};

      // need to calculate ET from eta, sin(theta) =  sech(eta)
      std::cout << "ihcal_eta: " << ihcal_tower[0] << " ihcal_phi: " << ihcal_tower[1] << std::endl;

      for (int ieta = ihcal_ieta - 1; ieta <= ihcal_ieta + 1; ieta++)
      {
        for (int iphi = ihcal_iphi - 1; iphi <= ihcal_iphi + 1; iphi++)
        {
          int temp_ieta = ieta;
          int temp_iphi = iphi;
          shift_tower_index(ieta, iphi, 24, 64);
          // std::cout<<"ieta: "<<ieta<<" iphi: "<<iphi<<std::endl;
          if (ieta < 0)
          {
            ieta = temp_ieta;
            iphi = temp_iphi;
            continue;
          }
          unsigned int towerinfokey = TowerInfoDefs::encode_hcal(ieta, iphi);
          TowerInfo *towerinfo = ihcalTowerContainer->get_tower_at_key(towerinfokey);
          const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);

          ieta = temp_ieta;
          iphi = temp_iphi;
          if (!towerinfo)
          {
            std::cout << "No towerinfo for tower key " << towerinfokey << std::endl;
            std::cout << "ieta: " << ieta << " iphi: " << iphi << std::endl;
            continue;
          }
          if (!(towerinfo->get_isGood()))
            continue;
          RawTowerGeom *tower_geom = geomIH->get_tower_geometry(key);
          if (!tower_geom)
          {
            std::cout << "No tower geometry for tower key " << key << std::endl;
            continue;
          }
          float energy = towerinfo->get_energy();
          float eta = getTowerEta(tower_geom, 0, 0, m_vertex);
          float sintheta = 1 / cosh(eta);
          float Et = energy * sintheta;

          ihcalEt33[ieta - ihcal_ieta + 1][iphi - ihcal_iphi + 1] = Et;
        }
      }
      std::cout << "ohcal_eta: " << ohcal_tower[0] << " ohcal_phi: " << ohcal_tower[1] << std::endl;
      for (int ieta = ohcal_ieta - 1; ieta <= ohcal_ieta + 1; ieta++)
      {
        for (int iphi = ohcal_iphi - 1; iphi <= ohcal_iphi + 1; iphi++)
        {
          int temp_ieta = ieta;
          int temp_iphi = iphi;
          shift_tower_index(ieta, iphi, 24, 64);
          if (ieta < 0)
          {
            ieta = temp_ieta;
            iphi = temp_iphi;
            continue;
          }
          unsigned int towerinfokey = TowerInfoDefs::encode_hcal(ieta, iphi);
          TowerInfo *towerinfo = ohcalTowerContainer->get_tower_at_key(towerinfokey);
          const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
          ieta = temp_ieta;
          iphi = temp_iphi;
          if (!towerinfo)
          {
            std::cout << "No towerinfo for tower key " << towerinfokey << std::endl;
            std::cout << "ieta: " << ieta << " iphi: " << iphi << std::endl;
            continue;
          }

          RawTowerGeom *tower_geom = geomOH->get_tower_geometry(key);
          float energy = towerinfo->get_energy();
          float eta = getTowerEta(tower_geom, 0, 0, m_vertex);
          float sintheta = 1 / cosh(eta);
          float Et = energy * sintheta;

          ohcalEt33[ieta - ohcal_ieta + 1][iphi - ohcal_iphi + 1] = Et;
        }
      }
      // std::cout<<"ihcal_et33: "<<ihcalEt33[1][1]<<" ohcal_et33: "<<ohcalEt33[1][1]<<std::endl;
      ihcal_et = ihcalEt33[1][1];
      ohcal_et = ohcalEt33[1][1];

      for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3; j++)
        {
          ihcal_et33 += ihcalEt33[i][j];
          ohcal_et33 += ohcalEt33[i][j];
          if (i == 1 || j == 1 + ihcal_tower[2])
          {
            if (j == 1 || i == 1 + ihcal_tower[3])
            {
              ihcal_et22 += ihcalEt33[i][j];
            }
          }
          if (i == 1 || j == 1 + ohcal_tower[2])
          {
            if (j == 1 || i == 1 + ohcal_tower[3])
            {
              ohcal_et22 += ohcalEt33[i][j];
            }
          }
        }
      }

      cluster_E[i][ncluster[i]] = E;
      cluster_Et[i][ncluster[i]] = ET;
      cluster_Eta[i][ncluster[i]] = eta;
      cluster_Phi[i][ncluster[i]] = phi;
      cluster_prob[i][ncluster[i]] = prob;
      cluster_CNN_prob[i][ncluster[i]] = CNNprob;
      cluster_truthtrkID[i][ncluster[i]] = trackid;
      cluster_pid[i][ncluster[i]] = pid;
      cluster_iso_02[i][ncluster[i]] = clusteriso[0];
      cluster_iso_03[i][ncluster[i]] = clusteriso[1];
      cluster_iso_04[i][ncluster[i]] = clusteriso[2];
      cluster_iso_04_emcal[i][ncluster[i]] = emcalET_04 - ET;
      cluster_iso_04_hcalin[i][ncluster[i]] = ihcalET_04;
      cluster_iso_04_hcalout[i][ncluster[i]] = ohcalET_04;
      cluster_e1[i][ncluster[i]] = showershape[8];
      cluster_e2[i][ncluster[i]] = showershape[9];
      cluster_e3[i][ncluster[i]] = showershape[10];
      cluster_e4[i][ncluster[i]] = showershape[11];
      cluster_ietacent[i][ncluster[i]] = showershape[4];
      cluster_iphicent[i][ncluster[i]] = showershape[5];
      cluster_weta[i][ncluster[i]] = sqrt(showershape[6]);
      cluster_wphi[i][ncluster[i]] = sqrt(showershape[7]);
      cluster_detamax[i][ncluster[i]] = detamax;
      cluster_dphimax[i][ncluster[i]] = dphimax;
      cluster_et1[i][ncluster[i]] = showershape[0];
      cluster_et2[i][ncluster[i]] = showershape[1];
      cluster_et3[i][ncluster[i]] = showershape[2];
      cluster_et4[i][ncluster[i]] = showershape[3];
      cluster_e11[i][ncluster[i]] = e11;
      cluster_e22[i][ncluster[i]] = e22;
      cluster_e13[i][ncluster[i]] = e13;
      cluster_e15[i][ncluster[i]] = e15;
      cluster_e17[i][ncluster[i]] = e17;
      cluster_e31[i][ncluster[i]] = e31;
      cluster_e51[i][ncluster[i]] = e51;
      cluster_e71[i][ncluster[i]] = e71;
      cluster_e33[i][ncluster[i]] = e33;
      cluster_e35[i][ncluster[i]] = e35;
      cluster_e37[i][ncluster[i]] = e37;
      cluster_e53[i][ncluster[i]] = e53;
      cluster_e73[i][ncluster[i]] = e73;
      cluster_e55[i][ncluster[i]] = e55;
      cluster_e57[i][ncluster[i]] = e57;
      cluster_e75[i][ncluster[i]] = e75;
      cluster_e77[i][ncluster[i]] = e77;
      cluster_w32[i][ncluster[i]] = w32;
      cluster_e32[i][ncluster[i]] = e32;
      cluster_w72[i][ncluster[i]] = w72;
      cluster_e72[i][ncluster[i]] = e72;
      cluster_ihcal_et[i][ncluster[i]] = ihcal_et;
      cluster_ohcal_et[i][ncluster[i]] = ohcal_et;
      cluster_ihcal_et22[i][ncluster[i]] = ihcal_et22;
      cluster_ohcal_et22[i][ncluster[i]] = ohcal_et22;
      cluster_ihcal_et33[i][ncluster[i]] = ihcal_et33;
      cluster_ohcal_et33[i][ncluster[i]] = ohcal_et33;
      cluster_ihcal_ieta[i][ncluster[i]] = ihcal_ieta;
      cluster_ihcal_iphi[i][ncluster[i]] = ihcal_iphi;
      cluster_ohcal_ieta[i][ncluster[i]] = ohcal_ieta;
      cluster_ohcal_iphi[i][ncluster[i]] = ohcal_iphi;

      ncluster[i]++;
      if (ncluster[i] > nclustermax)
      {
        std::cout << "ncluster exceed the max range: " << ncluster[i] << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }

      // std::cout << "prob: " << prob << " clusteriso: " << clusteriso[0] << " ET: " << ET << std::endl;
    }
    std::cout << "done with cluster container: " << clusternamelist[i] << std::endl;
  }
  bool saveevent = false;
  for (int i = 0; i < nclustercontainer; i++)
  {
    if (ncluster[i] > 0)
      saveevent = true;
  }
  if (nparticles > 0)
    saveevent = true;
  if (saveevent)
  {
    slimtree->Fill();
  }

  std::cout << "done with cluster" << std::endl;

  // delete clustereval;
  // delete trutheval;

  return Fun4AllReturnCodes::EVENT_OK;
}

std::vector<int> CaloAna24::find_closest_hcal_tower(float eta, float phi, RawTowerGeomContainer *geom, TowerInfoContainer *towerContainer, float vertex_z, bool isihcal)
{
  int matchedieta = -1;
  int matchediphi = -1;
  double matchedeta = -999;
  double matchedphi = -999;
  unsigned int ntowers = towerContainer->size();

  float minR = 999;

  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    TowerInfo *tower = towerContainer->get_tower_at_channel(channel);
    if (!tower)
    {
      continue;
    }
    unsigned int towerkey = towerContainer->encode_key(channel);

    int ieta = towerContainer->getTowerEtaBin(towerkey);
    int iphi = towerContainer->getTowerPhiBin(towerkey);
    RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
    if (!isihcal)
    {
      key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
    }
    RawTowerGeom *tower_geom = geom->get_tower_geometry(key);
    double this_phi = tower_geom->get_phi();
    double this_eta = getTowerEta(tower_geom, 0, 0, vertex_z);
    double dR = deltaR(eta, this_eta, phi, this_phi);
    if (dR < minR)
    {
      minR = dR;
      matchedieta = ieta;
      matchediphi = iphi;
      matchedeta = this_eta;
      matchedphi = this_phi;
    }
  }
  float deta = eta - matchedeta;
  float dphi = phi - matchedphi;
  float dphiwrap = 2 * M_PI - std::abs(dphi);
  if (std::abs(dphiwrap) < std::abs(dphi))
  {
    dphi = (dphi > 0) ? -dphiwrap : dphiwrap;
  }
  int dphisign = (dphi > 0) ? 1 : -1;
  int detasign = (deta > 0) ? 1 : -1;

  std::vector<int> result = {matchedieta, matchediphi, detasign, dphisign};
  return result;
}

std::pair<int, int> CaloAna24::photon_type(int barcode)
{
  // check if the Genevent is null
  std::pair<int, int> photonclasspair = {-1, -1};
  if (!singal_event)
  {
    std::cout << "Genevent is null" << std::endl;
    return photonclasspair;
  }
  HepMC::GenParticle *particle = singal_event->barcode_to_particle(barcode);
  // check if pid is 22
  if (particle->pdg_id() != 22)
  {
    std::cout << "particle is not photon" << std::endl;
    return photonclasspair;
  }
  // find the production vertex
  HepMC::GenVertex *vertex = particle->production_vertex();
  if (!vertex)
  {
    std::cout << "vertex is null" << std::endl;
    return photonclasspair;
  }
  // find the incoming particles
  HepMC::GenVertex::particles_in_const_iterator inItr = vertex->particles_in_const_begin();
  std::vector<HepMC::GenParticle *> incoming_particles;
  for (; inItr != vertex->particles_in_const_end(); ++inItr)
  {
    incoming_particles.push_back(*inItr);
  }
  // check if there is only one incoming particle and pid is 22, noticed that in pythia there are vertex with photon in and photon out
  while (incoming_particles.size() == 1 && incoming_particles[0]->pdg_id() == 22)
  {
    // find the production vertex
    vertex = incoming_particles[0]->production_vertex();
    if (!vertex)
    {
      std::cout << "vertex is null" << std::endl;
      return photonclasspair;
    }
    // find the incoming particles
    inItr = vertex->particles_in_const_begin();
    incoming_particles.clear();
    for (; inItr != vertex->particles_in_const_end(); ++inItr)
    {
      incoming_particles.push_back(*inItr);
    }
  }
  std::vector<HepMC::GenParticle *> outgoing_particles;
  // find the outgoing particles
  HepMC::GenVertex::particles_out_const_iterator outItr = vertex->particles_out_const_begin();
  for (; outItr != vertex->particles_out_const_end(); ++outItr)
  {
    outgoing_particles.push_back(*outItr);
  }
  // direct photon:1, fragmentation photon:2 , decayed photon:3, can't identify: 0;
  int photonclass = 0;
  std::set<int> outgoing_pid;
  for (auto particle : outgoing_particles)
  {
    outgoing_pid.insert(particle->pdg_id());
  }
  // make sure there is photon in it
  if (outgoing_pid.find(22) == outgoing_pid.end())
  {
    std::cout << "no photon in the outgoing particles" << std::endl;
    return photonclasspair;
  }
  int incoming_pid = 0;
  if (incoming_particles.size() > 0)
    incoming_pid = incoming_particles.at(0)->pdg_id();

  // direct photon 2->2 both incoming are quark or gluons
  if (incoming_particles.size() == 2 && outgoing_particles.size() == 2)
  {
    if (abs(incoming_particles.at(0)->pdg_id()) <= 22 && abs(incoming_particles.at(1)->pdg_id()) <= 22)
    {
      if (abs(outgoing_particles.at(0)->pdg_id()) <= 22 && abs(outgoing_particles.at(1)->pdg_id()) <= 22)
      {
        photonclass = 1;
      }
    }
  }
  // fragmentation photon and  decayed photon should only have one incoming
  else if (incoming_particles.size() == 1)
  {
    // fragmentation photon 1->2
    if (abs(incoming_particles.at(0)->pdg_id()) <= 11 && outgoing_particles.size() == 2)
    {

      // both 22 and incoming pid should be in
      if (outgoing_pid.find((incoming_particles.at(0)->pdg_id())) != outgoing_pid.end())
      {
        photonclass = 2;
      }
    }
    if (abs(incoming_particles.at(0)->pdg_id()) > 37)
    {
      photonclass = 3;
    }
  }

  if (photonclass == 0)
  {
    std::cout << "can't identify the photon type" << std::endl;
    // debug print
    std::cout << "incoming particles: ";
    for (auto particle : incoming_particles)
    {
      std::cout << particle->pdg_id() << " ";
    }
    std::cout << std::endl;
    std::cout << "outgoing particles: ";
    for (auto particle : outgoing_particles)
    {
      std::cout << particle->pdg_id() << " ";
    }
    std::cout << std::endl;
  }
  photonclasspair = {photonclass, incoming_pid};
  return photonclasspair;
}

//____________________________________________________________________________..
int CaloAna24::End(PHCompositeNode *topNode)
{

  fout->cd();

  fout->Write();
  fout->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloAna24::Reset(PHCompositeNode *topNode)
{
  std::cout << "CaloAna24::Reset(PHCompositeNode *topNode) being Reset"
            << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void CaloAna24::Print(const std::string &what) const
{
  std::cout << "CaloAna24::Print(const std::string &what) const Printing "
               "info for "
            << what << std::endl;
}

double CaloAna24::getTowerEta(RawTowerGeom *tower_geom, double vx, double vy, double vz)
{
  float r;
  if (vx == 0 && vy == 0 && vz == 0)
  {
    r = tower_geom->get_eta();
  }
  else
  {
    double radius = sqrt((tower_geom->get_center_x() - vx) * (tower_geom->get_center_x() - vx) + (tower_geom->get_center_y() - vy) * (tower_geom->get_center_y() - vy));
    double theta = atan2(radius, tower_geom->get_center_z() - vz);
    r = -log(tan(theta / 2.));
  }
  return r;
}

float CaloAna24::calculateET(float eta, float phi, float dR, int layer) // layer: 0 EMCal, 1 IHCal, 2 OHCal
{
  float ET = 0;
  RawTowerGeomContainer *geomcontainer = nullptr;
  TowerInfoContainer *towercontainer = nullptr;
  RawTowerDefs::CalorimeterId caloid = RawTowerDefs::CalorimeterId::CEMC;

  if (layer == 0)
  {
    geomcontainer = geomEM;
    towercontainer = emcTowerContainer;
    caloid = RawTowerDefs::CalorimeterId::CEMC;
  }
  else if (layer == 1)
  {
    geomcontainer = geomIH;
    towercontainer = ihcalTowerContainer;
    caloid = RawTowerDefs::CalorimeterId::HCALIN;
  }
  else if (layer == 2)
  {
    geomcontainer = geomOH;
    towercontainer = ohcalTowerContainer;
    caloid = RawTowerDefs::CalorimeterId::HCALOUT;
  }
  else
  {
    std::cout << "Invalid layer" << std::endl;
    return ET;
  }
  float ntowers = towercontainer->size();
  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    TowerInfo *tower = towercontainer->get_tower_at_channel(channel);
    if (!tower)
    {
      continue;
    }
    if (tower->get_isGood() == false)
    {
      continue;
    }
    unsigned int towerkey = towercontainer->encode_key(channel);
    int ieta = towercontainer->getTowerEtaBin(towerkey);
    int iphi = towercontainer->getTowerPhiBin(towerkey);
    RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(caloid, ieta, iphi);
    RawTowerGeom *tower_geom = geomcontainer->get_tower_geometry(key);
    double this_phi = tower_geom->get_phi();
    double this_eta = getTowerEta(tower_geom, 0, 0, vertexz);
    if (deltaR(eta, this_eta, phi, this_phi) < dR)
    {
      ET += tower->get_energy()/cosh(this_eta);
    }
  }
  return ET;
}
