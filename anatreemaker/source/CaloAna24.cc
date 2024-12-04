

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
  slimtree->Branch("particle_truth_iso_02", particle_truth_iso_02, "particle_truth_iso_02[nparticles]/F");
  slimtree->Branch("particle_truth_iso_03", particle_truth_iso_03, "particle_truth_iso_03[nparticles]/F");
  slimtree->Branch("particle_truth_iso_04", particle_truth_iso_04, "particle_truth_iso_04[nparticles]/F");
  slimtree->Branch("particle_converted", particle_converted, "particle_converted[nparticles]/I");
  for (int i = 0; i < nclustercontainer; i++)
  {
    slimtree->Branch(Form("ncluster_%s", clusternamelist[i].c_str()), &ncluster[i], Form("ncluster_%s/I", clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_E_%s", clusternamelist[i].c_str()), cluster_E[i], Form("cluster_E_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_Et_%s", clusternamelist[i].c_str()), cluster_Et[i], Form("cluster_Et_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_Eta_%s", clusternamelist[i].c_str()), cluster_Eta[i], Form("cluster_Eta_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_Phi_%s", clusternamelist[i].c_str()), cluster_Phi[i], Form("cluster_Phi_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_prob_%s", clusternamelist[i].c_str()), cluster_prob[i], Form("cluster_prob_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_truthtrkID_%s", clusternamelist[i].c_str()), cluster_truthtrkID[i], Form("cluster_truthtrkID_%s[ncluster_%s]/I", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_pid_%s", clusternamelist[i].c_str()), cluster_pid[i], Form("cluster_pid_%s[ncluster_%s]/I", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_iso_02_%s", clusternamelist[i].c_str()), cluster_iso_02[i], Form("cluster_iso_02_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_iso_03_%s", clusternamelist[i].c_str()), cluster_iso_03[i], Form("cluster_iso_03_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_iso_04_%s", clusternamelist[i].c_str()), cluster_iso_04[i], Form("cluster_iso_04_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e1_%s", clusternamelist[i].c_str()), cluster_e1[i], Form("cluster_e1_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e2_%s", clusternamelist[i].c_str()), cluster_e2[i], Form("cluster_e2_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e3_%s", clusternamelist[i].c_str()), cluster_e3[i], Form("cluster_e3_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_e4_%s", clusternamelist[i].c_str()), cluster_e4[i], Form("cluster_e4_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_et1_%s", clusternamelist[i].c_str()), cluster_et1[i], Form("cluster_et1_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_et2_%s", clusternamelist[i].c_str()), cluster_et2[i], Form("cluster_et2_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_et3_%s", clusternamelist[i].c_str()), cluster_et3[i], Form("cluster_et3_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
    slimtree->Branch(Form("cluster_et4_%s", clusternamelist[i].c_str()), cluster_et4[i], Form("cluster_et4_%s[ncluster_%s]/F", clusternamelist[i].c_str(), clusternamelist[i].c_str()));
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

  if (isMC)
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
        if (mbenrgy[i] > 0.5 && i < 64)
          northhit += 1;
        if (mbenrgy[i] > 0.5 && i > 63)
          southhit += 1;
      }
      mbdnorthhit = northhit;
      mbdsouthhit = southhit;
    }
  }

  float m_vertex = -9999;
  std::vector<TLorentzVector> goodcluster;
  GlobalVertexMap *vertexmap =
      findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

  if (!vertexmap)
  {
    std::cout << "GlobalVertexMap node is missing" << std::endl;
  }
  if (vertexmap && !vertexmap->empty())
  {
    GlobalVertex *vtx = vertexmap->begin()->second;
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

  vertexz = m_vertex;
  // set of primary particles
  std::set<PHG4Particle *> primary_particles;
  std::set<PHG4Particle *> primary_photon_candidates;
  std::set<PHG4Particle *> photonsfrompi0;
  std::set<PHG4Particle *> photonsfrometa;
  std::set<PHG4Particle *> badphotons;
  std::map<PHG4Particle *, std::vector<float>> photontruthiso;
  if (isMC)
  {

    CaloEvalStack caloevalstack(topNode, "CEMC");
    clustereval = caloevalstack.get_rawcluster_eval();
    clustereval->set_usetowerinfo(true);
    clustereval->next_event(topNode);
    trutheval = caloevalstack.get_truth_eval();
    // CaloRawTowerEval *towereval = caloevalstack.get_rawtower_eval();
    //  clustereval->next_event(topNode);
    // CaloTruthEval *trutheval = m_caloevalstack->get_truth_eval();
    std::cout << "trutheval: " << trutheval << " clustereval: " << clustereval << std::endl;
    // trutheval->next_event(topNode);

    PHG4TruthInfoContainer *truthinfo =
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
      //int trackid1 = truth_photon->get_track_id();
      // std::cout<<"trackid1: "<<trackid1<<std::endl;
      // check if the photon is from pi0 or eta
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
      if (abs(photon1.Eta()) > 1.1)
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
            badphotons.insert(truth_photon);
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
      if (abs(p1.Eta()) > 1.1)
        continue;
      int barcode = truth->get_barcode();

      bool verbosephoton = false;
      //bool ispromptphoton = false;
      int photonclass = 0;
      if (pid == 22)
      {
        verbosephoton = true;
        //ispromptphoton = true;
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
        photonclass = photon_type(barcode);
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
        if (p1.E() > 5)
        {
          int converted = false;
          if (badphotons.find(truth) != badphotons.end())
          {
            converted = true;
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
      if (abs(eta) > 1.0)
        continue;
      // if (ET > 1) h_ET->Fill(ET);
      if (ET < clusterpTmin)
        continue;

      // Array for storing the isolation energy for different radii
      float clusteriso[nRadii];

      // Loop to calculate the isolation energy for each radius
      for (int i = 0; i < nRadii; ++i)
      {
        clusteriso[i] = recoCluster->get_et_iso(2 + i, false, true);
        // std::cout << "clusteriso: " << clusteriso[i] << std::endl;
      }

      if (ET > maxclusterpt)
      {
        maxclusterpt = ET;
      }
      int trackid = -1;
      float clusterE = E_vec_cluster_Full.mag();
      int pid = 0;
      if (isMC)
      {
        PHG4Particle *maxPrimary = clustereval->max_truth_primary_particle_by_energy(recoCluster);
        if (ET > 5)
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

      cluster_E[i][ncluster[i]] = E;
      cluster_Et[i][ncluster[i]] = ET;
      cluster_Eta[i][ncluster[i]] = eta;
      cluster_Phi[i][ncluster[i]] = phi;
      cluster_prob[i][ncluster[i]] = prob;
      cluster_truthtrkID[i][ncluster[i]] = trackid;
      cluster_pid[i][ncluster[i]] = pid;
      cluster_iso_02[i][ncluster[i]] = clusteriso[0];
      cluster_iso_03[i][ncluster[i]] = clusteriso[1];
      cluster_iso_04[i][ncluster[i]] = clusteriso[2];
      cluster_e1[i][ncluster[i]] = showershape[0];
      cluster_e2[i][ncluster[i]] = showershape[1];
      cluster_e3[i][ncluster[i]] = showershape[2];
      cluster_e4[i][ncluster[i]] = showershape[3];
      cluster_et1[i][ncluster[i]] = showershape[8];
      cluster_et2[i][ncluster[i]] = showershape[9];
      cluster_et3[i][ncluster[i]] = showershape[10];
      cluster_et4[i][ncluster[i]] = showershape[11];
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

int CaloAna24::photon_type(int barcode)
{
  // check if the Genevent is null
  if (!singal_event)
  {
    std::cout << "Genevent is null" << std::endl;
    return -1;
  }
  HepMC::GenParticle *particle = singal_event->barcode_to_particle(barcode);
  // check if pid is 22
  if (particle->pdg_id() != 22)
  {
    std::cout << "particle is not photon" << std::endl;
    return -1;
  }
  // find the production vertex
  HepMC::GenVertex *vertex = particle->production_vertex();
  if (!vertex)
  {
    std::cout << "vertex is null" << std::endl;
    return -1;
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
      return -1;
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
    outgoing_pid.insert(abs(particle->pdg_id()));
  }
  // make sure there is photon in it
  if (outgoing_pid.find(22) == outgoing_pid.end())
  {
    std::cout << "no photon in the outgoing particles" << std::endl;
    return -1;
  }

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
  return photonclass;
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
