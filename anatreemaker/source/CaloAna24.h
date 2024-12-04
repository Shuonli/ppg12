// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOANA24_H
#define CALOANA24_H

#include <fun4all/SubsysReco.h>

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

class CaloAna24 : public SubsysReco {
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

private:
  int ievent = 0;
  TFile *fout;

  TTree* slimtree;

  bool isMC{true};
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
  static const int nparticlesmax = 10000;
  int nparticles{0};
  float particle_E[nparticlesmax] = {0};
  float particle_Pt[nparticlesmax] = {0};
  float particle_Eta[nparticlesmax] = {0};
  float particle_Phi[nparticlesmax] = {0};
  int particle_pid[nparticlesmax] = {0};
  int particle_trkid[nparticlesmax] = {0};
  int particle_isprompt_photon[nparticlesmax] = {0};
  float particle_truth_iso_02[nparticlesmax] = {0};
  float particle_truth_iso_03[nparticlesmax] = {0};
  float particle_truth_iso_04[nparticlesmax] = {0};
  int particle_converted[nparticlesmax] = {0};


  std::vector<std::string> clusternamelist = {"CLUSTERINFO_CEMC"};
  static const int nclustercontainer = 1;
  //cluster wise stuff
  static const int nclustermax = 10000;
  int ncluster[nclustercontainer] = {0};
  float cluster_E[nclustercontainer][nclustermax] = {0};
  float cluster_Et[nclustercontainer][nclustermax] = {0};
  float cluster_Eta[nclustercontainer][nclustermax] = {0};
  float cluster_Phi[nclustercontainer][nclustermax] = {0};
  float cluster_prob[nclustercontainer][nclustermax] = {0};
  int cluster_truthtrkID[nclustercontainer][nclustermax] = {0};
  int cluster_pid[nclustercontainer][nclustermax] = {0};
  float cluster_iso_02[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03[nclustercontainer][nclustermax] = {0};
  float cluster_iso_04[nclustercontainer][nclustermax] = {0};
  //shower shapes
  float cluster_e1[nclustercontainer][nclustermax] = {0};
  float cluster_e2[nclustercontainer][nclustermax] = {0};
  float cluster_e3[nclustercontainer][nclustermax] = {0};
  float cluster_e4[nclustercontainer][nclustermax] = {0};
  float cluster_et1[nclustercontainer][nclustermax] = {0};
  float cluster_et2[nclustercontainer][nclustermax] = {0};
  float cluster_et3[nclustercontainer][nclustermax] = {0};
  float cluster_et4[nclustercontainer][nclustermax] = {0};
  // Number of radii
  static const int nRadii = 3;
  



  
  TH3F* h_tracking_radiograph;

  const int truthisocut = 4;

  int process_cluster(std::vector<TLorentzVector> goodcluster);

  float DeltaR(TLorentzVector photon1, TLorentzVector photon2) {
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
};

#endif // CALOANA24_H
