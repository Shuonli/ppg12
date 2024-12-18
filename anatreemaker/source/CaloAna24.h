// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOANA24_H
#define CALOANA24_H

#include <fun4all/SubsysReco.h>

#include <phool/onnxlib.h>

#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile2D.h>
#include <TTree.h>
#include <string>
#include <TLorentzVector.h>

class PHCompositeNode;
class CaloEvalStack;
class CaloRawClusterEval;
class CaloTruthEval;
class RawTowerGeom;
class RawTowerGeomContainer;
class TowerInfoContainer;
class PHG4TruthInfoContainer;
class CaloEvalStack;

namespace HepMC
{
  class GenEvent;
}

class CaloAna24 : public SubsysReco
{
public:
  CaloAna24(const std::string &name = "CaloAna24");

  ~CaloAna24() override;

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  void set_isMC(bool isMC_) { isMC = isMC_; }

  void set_isSingleParticle(bool isSingleParticle_)
  {
    isSingleParticle = isSingleParticle_;
    isMC = true;
  }

private:
  int ievent = 0;
  Ort::Session *onnxmodule{nullptr};
  std::string m_modelPath{"/sphenix/u/shuhang98/core_patch/coresoftware/offline/packages/CaloReco/functional_model_single.onnx"};

  TFile *fout;

  TTree *slimtree;

  bool isMC{true};
  bool isSingleParticle{false};
  int m_scaledtrigger[32] = {0};
  bool initilized = false;
  long long initscaler[32][3] = {0};
  long long currentscaler[32][3] = {0};
  bool scaledtrigger[32] = {false};
  bool livetrigger[32] = {false};
  int nscaledtrigger[32] = {0};
  int nlivetrigger[32] = {0};

  float vertexz{-9999};
  int mbdnorthhit{0};
  int mbdsouthhit{0};
  float vertexz_truth{-9999};
  int m_pythiaid{-9999};
  float particlepTmin{1};
  static const int nparticlesmax = 10000;
  int nparticles{0};
  float particle_E[nparticlesmax] = {0};
  float particle_Pt[nparticlesmax] = {0};
  float particle_Eta[nparticlesmax] = {0};
  float particle_Phi[nparticlesmax] = {0};
  int particle_pid[nparticlesmax] = {0};
  int particle_trkid[nparticlesmax] = {0};
  int particle_photonclass[nparticlesmax] = {0};
  int particle_photon_mother_pid[nparticlesmax] = {0};
  float particle_truth_iso_02[nparticlesmax] = {0};
  float particle_truth_iso_03[nparticlesmax] = {0};
  float particle_truth_iso_04[nparticlesmax] = {0};

  int particle_converted[nparticlesmax] = {0};

  static const int ndaughtermax = 100;
  int ndaughter{0};
  int daughter_pid[ndaughtermax] = {0};
  int daughter_parent_trackid[ndaughtermax] = {0};
  float daughter_E[ndaughtermax] = {0};
  float daughter_Pt[ndaughtermax] = {0};
  float daughter_Eta[ndaughtermax] = {0};
  float daughter_Phi[ndaughtermax] = {0};


  std::vector<std::string> clusternamelist = {"CLUSTERINFO_CEMC"};
  static const int nclustercontainer = 1;
  // cluster wise stuff
  float clusterpTmin{1};
  static const int nclustermax = 10000;
  int ncluster[nclustercontainer] = {0};
  float cluster_E[nclustercontainer][nclustermax] = {0};
  float cluster_Et[nclustercontainer][nclustermax] = {0};
  float cluster_Eta[nclustercontainer][nclustermax] = {0};
  float cluster_Phi[nclustercontainer][nclustermax] = {0};
  float cluster_prob[nclustercontainer][nclustermax] = {0};
  float cluster_CNN_prob[nclustercontainer][nclustermax] = {0};
  int cluster_truthtrkID[nclustercontainer][nclustermax] = {0};
  int cluster_pid[nclustercontainer][nclustermax] = {0};
  float cluster_iso_02[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03[nclustercontainer][nclustermax] = {0};
  float cluster_iso_04[nclustercontainer][nclustermax] = {0};
  float cluster_iso_04_emcal[nclustercontainer][nclustermax] = {0};
  float cluster_iso_04_hcalin[nclustercontainer][nclustermax] = {0};
  float cluster_iso_04_hcalout[nclustercontainer][nclustermax] = {0};

  // shower shapes
  float cluster_e1[nclustercontainer][nclustermax] = {0};
  float cluster_e2[nclustercontainer][nclustermax] = {0};
  float cluster_e3[nclustercontainer][nclustermax] = {0};
  float cluster_e4[nclustercontainer][nclustermax] = {0};
  float cluster_et1[nclustercontainer][nclustermax] = {0};
  float cluster_et2[nclustercontainer][nclustermax] = {0};
  float cluster_et3[nclustercontainer][nclustermax] = {0};
  float cluster_et4[nclustercontainer][nclustermax] = {0};
  float cluster_ietacent[nclustercontainer][nclustermax] = {0};
  float cluster_iphicent[nclustercontainer][nclustermax] = {0};
  float cluster_weta[nclustercontainer][nclustermax] = {0};
  float cluster_wphi[nclustercontainer][nclustermax] = {0};
  int cluster_detamax[nclustercontainer][nclustermax] = {0};
  int cluster_dphimax[nclustercontainer][nclustermax] = {0};

  float cluster_e11[nclustercontainer][nclustermax] = {0};
  float cluster_e22[nclustercontainer][nclustermax] = {0};
  float cluster_e13[nclustercontainer][nclustermax] = {0};
  float cluster_e15[nclustercontainer][nclustermax] = {0};
  float cluster_e17[nclustercontainer][nclustermax] = {0};
  float cluster_e31[nclustercontainer][nclustermax] = {0};
  float cluster_e51[nclustercontainer][nclustermax] = {0};
  float cluster_e71[nclustercontainer][nclustermax] = {0};
  float cluster_e33[nclustercontainer][nclustermax] = {0};
  float cluster_e35[nclustercontainer][nclustermax] = {0};
  float cluster_e37[nclustercontainer][nclustermax] = {0};
  float cluster_e53[nclustercontainer][nclustermax] = {0};
  float cluster_e73[nclustercontainer][nclustermax] = {0};
  float cluster_e55[nclustercontainer][nclustermax] = {0};
  float cluster_e57[nclustercontainer][nclustermax] = {0};
  float cluster_e75[nclustercontainer][nclustermax] = {0};
  float cluster_e77[nclustercontainer][nclustermax] = {0};
  float cluster_w32[nclustercontainer][nclustermax] = {0};
  float cluster_e32[nclustercontainer][nclustermax] = {0};
  float cluster_w72[nclustercontainer][nclustermax] = {0};
  float cluster_e72[nclustercontainer][nclustermax] = {0};

  float cluster_ihcal_et[nclustercontainer][nclustermax] = {0};
  float cluster_ohcal_et[nclustercontainer][nclustermax] = {0};
  float cluster_ihcal_et22[nclustercontainer][nclustermax] = {0};
  float cluster_ohcal_et22[nclustercontainer][nclustermax] = {0};
  float cluster_ihcal_et33[nclustercontainer][nclustermax] = {0};
  float cluster_ohcal_et33[nclustercontainer][nclustermax] = {0};
  int cluster_ihcal_ieta[nclustercontainer][nclustermax] = {0};
  int cluster_ihcal_iphi[nclustercontainer][nclustermax] = {0};
  int cluster_ohcal_ieta[nclustercontainer][nclustermax] = {0};
  int cluster_ohcal_iphi[nclustercontainer][nclustermax] = {0};
  // Number of radii
  static const int nRadii = 3;

  TH3F *h_tracking_radiograph;

  const int truthisocut = 4;

  int process_cluster(std::vector<TLorentzVector> goodcluster);

  std::pair<int, int> photon_type(int barcode);

  void shift_tower_index(int &ieta, int &iphi, int maxeta, int maxphi)
  {
    if (ieta < 0)
      ieta = -1;
    if (ieta >= maxeta)
      ieta = -1;
    if (iphi < 0)
      iphi += maxphi;
    if (iphi >= maxphi)
      iphi -= maxphi;
  }

  double getTowerEta(RawTowerGeom *tower_geom, double vx, double vy, double vz);

  float calculateET(float eta, float phi, float dR, int layer); // layer: 0 EMCal, 1 IHCal, 2 OHCal

  std::vector<int> find_closest_hcal_tower(float eta, float phi, RawTowerGeomContainer *rawtowergeom, TowerInfoContainer *towercontainer, float vertex_z, bool isihcal);

  inline /*const*/ float deltaR(float eta1, float eta2, float phi1, float phi2)
  {
    float deta = eta1 - eta2;
    float dphi = phi1 - phi2;
    if (dphi > M_PI)
      dphi -= 2 * M_PI; // corrects to keep range -pi to pi
    if (dphi < -1 * M_PI)
      dphi += 2 * M_PI; // corrects to keep range -pi to pi
    return sqrt(deta * deta + dphi * dphi);
  }

  float DeltaR(TLorentzVector photon1, TLorentzVector photon2)
  {
    float deta = photon1.PseudoRapidity() - photon2.PseudoRapidity();
    float dphi = abs(photon1.Phi() - photon2.Phi());

    if (dphi > M_PI)
      dphi = 2 * M_PI - dphi;
    float dr = sqrt(deta * deta + dphi * dphi);
    return dr;
  }

  std::unique_ptr<CaloEvalStack> m_caloevalstack;
  CaloRawClusterEval *clustereval{nullptr};
  CaloTruthEval *trutheval{nullptr};
  HepMC::GenEvent *singal_event{nullptr};
  PHG4TruthInfoContainer *truthinfo{nullptr};
  CaloEvalStack *caloevalstack{nullptr};

  RawTowerGeomContainer *geomEM{nullptr};
  RawTowerGeomContainer *geomIH{nullptr};
  RawTowerGeomContainer *geomOH{nullptr};

  TowerInfoContainer *emcTowerContainer{nullptr};
  TowerInfoContainer *ihcalTowerContainer{nullptr};
  TowerInfoContainer *ohcalTowerContainer{nullptr};
};

#endif // CALOANA24_H
