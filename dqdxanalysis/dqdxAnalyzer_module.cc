////////////////////////////////////////////////////////////////////////
// Class:       dqdxAnalyzer
// Plugin Type: analyzer (art v2_05_01)
// File:        dqdxAnalyzer_module.cc
//
// Generated at Wed Feb 14 10:17:20 2018 by Nicolo Foppiani using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"

#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibService.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "EnergyHelper.h"
#include "GeometryHelper.h"
#include "McPfpMatch.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "TTree.h"

class dqdxAnalyzer;

class dqdxAnalyzer : public art::EDAnalyzer
{
public:
  explicit dqdxAnalyzer(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  dqdxAnalyzer(dqdxAnalyzer const &) = delete;
  dqdxAnalyzer(dqdxAnalyzer &&) = delete;
  dqdxAnalyzer &operator=(dqdxAnalyzer const &) = delete;
  dqdxAnalyzer &operator=(dqdxAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;
  void reconfigure(fhicl::ParameterSet const &p) override;
  void clear();

  void recoTrueMatching(art::Event const &evt, art::Ptr<recob::PFParticle> const &pfparticle);
  void shiftTruePosition(double true_point[3], double true_point_shifted[3]);
  void shiftTruePosition(double true_point[3], double true_time, double true_point_shifted[3]);
  void shiftTruePosition_onlySCE(double true_point[3], double true_point_shifted[3]);
  void shiftTruePosition_onlyGeant(double true_point[3], double true_time, double true_point_shifted[3]);

  void trueNeutrinoInformation(art::Event const &evt);
  // double GetChargeCorrection(int plane, double x, double y, double z);

private:
  std::string _pfp_producer;

  // reco-true matching tools
  ubana::McPfpMatch _mcpfpMatcher;
  std::string _spacepointLabel;
  std::string _hitfinderLabel;
  std::string _geantModuleLabel;
  std::string _mcpHitAssLabel;
  bool _use_premade_ass;
  std::string _mctruthLabel;

  // dqdx tools
  lee::EnergyHelper energyHelper;
  lee::GeometryHelper geoHelper;

  TTree *fChargeTree;
  int fRun, fSubrun, fEvent;
  int fPdgCode, fPdgCodeParent;
  double fStartx, fStarty, fStartz;
  double fEndx, fEndy, fEndz;
  double fDirectionx, fDirectiony, fDirectionz;
  double fDistanceFromTrue;

  double fTrue_vx, fTrue_vy, fTrue_vz;
  double fTrue_v[3];
  int fTrue_nu_pdg, fTrue_ccnc;
  double fTrue_nu_energy, fTrue_theta, fTrue_vx_sce, fTrue_vy_sce, fTrue_vz_sce;
  double fTrue_v_sce[3];

  double fAngleZXplanePF;
  int fNhits, fNclusters;
  int fNhitsU, fNhitsV, fNhitsY;
  double fLength, fSpLength;

  // const lariov::TPCEnergyCalibProvider &energyCalibProvider = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();

  // reco-true matching information
  int fMatchedPdgCode;
  double fMatchedE;
  double fMatchedPx, fMatchedPy, fMatchedPz;
  double fMatchedVx, fMatchedVy, fMatchedVz;
  double fMatchedEndx, fMatchedEndy, fMatchedEndz;
  double fDistanceFromMatched, fMC_reco_costheta;
  double fFractionMatchedHits;
  unsigned int fMatchedTrackID;
  double fMatchedT;
  double fMatchedStartX_shifted, fMatchedStartY_shifted, fMatchedStartZ_shifted;
  double fMatchedStartX_shifted_onlySCE, fMatchedStartY_shifted_onlySCE, fMatchedStartZ_shifted_onlySCE;
  double fMatchedStartX_shifted_onlyGeant, fMatchedStartY_shifted_onlyGeant, fMatchedStartZ_shifted_onlyGeant;

  // dqdx information
  double _dQdxRectangleLength;
  double _dQdxRectangleWidth;
  double _dQdxRectangleLength_end;

  std::vector<double> fDQdx_hits_start;
  std::vector<int> fDQdx_wires_start;
  double fDQdx_start[3];
  double fDQdx_U_start, fDQdx_V_start, fDQdx_Y_start;
  int fn_hits_dQdx_U_start, fn_hits_dQdx_V_start, fn_hits_dQdx_Y_start;
  double fBox_start_z_start, fBox_start_x_start, fBox_direction_z_start, fBox_direction_x_start;
  double fReco_energy_U_start, fReco_energy_V_start, fReco_energy_Y_start;
  double fAngleZXplaneCluster_start, fDistance_starts_start;
  double fDistance_recotrue_starts_start, fDistance_recotrue_shifted_starts_start;

  double fMCShowerStartX, fMCShowerStartY, fMCShowerStartZ;
  double fMCShowerStartX_shifted, fMCShowerStartY_shifted, fMCShowerStartZ_shifted;
  double fMCShowerStartX_shifted_onlySCE, fMCShowerStartY_shifted_onlySCE, fMCShowerStartZ_shifted_onlySCE;
  double fMCShowerStartX_shifted_onlyGeant, fMCShowerStartY_shifted_onlyGeant, fMCShowerStartZ_shifted_onlyGeant;
  double fDistance_recotrue_mc_shower_starts_start, fDistance_recotrue_shifted_mc_shower_starts_start;

  std::vector<double> fDQdx_hits_end;
  std::vector<int> fDQdx_wires_end;
  double fDQdx_end[3];
  double fDQdx_U_end, fDQdx_V_end, fDQdx_Y_end;
  int fn_hits_dQdx_U_end, fn_hits_dQdx_V_end, fn_hits_dQdx_Y_end;
  double fBox_start_z_end, fBox_start_x_end, fBox_direction_z_end, fBox_direction_x_end;
  double fReco_energy_U_end, fReco_energy_V_end, fReco_energy_Y_end;
  double fAngleZXplaneCluster_end, fDistance_starts_end;

  double fReco_energy_U, fReco_energy_V, fReco_energy_Y;
  double fPitch_U, fPitch_V, fPitch_Y;
};

dqdxAnalyzer::dqdxAnalyzer(fhicl::ParameterSet const &p)
    : EDAnalyzer(p) // ,
                    // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;

  this->reconfigure(p);

  fChargeTree = tfs->make<TTree>("Charge", "Charge PF Particles Tree");
  fChargeTree->Branch("event", &fEvent, "event/i");
  fChargeTree->Branch("run", &fRun, "run/i");
  fChargeTree->Branch("subrun", &fSubrun, "subrun/i");
  fChargeTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
  fChargeTree->Branch("pdg_code_parent", &fPdgCodeParent, "pdg_code_parent/i");

  fChargeTree->Branch("true_nu_pdg", &fTrue_nu_pdg, "true_nu_pdg/d");
  fChargeTree->Branch("true_ccnc", &fTrue_ccnc, "true_ccnc/d");
  fChargeTree->Branch("true_nu_energy", &fTrue_nu_energy, "true_nu_energy/d");
  fChargeTree->Branch("true_theta", &fTrue_theta, "true_theta/d");

  fChargeTree->Branch("true_vx", &fTrue_vx, "true_vx/d");
  fChargeTree->Branch("true_vy", &fTrue_vy, "true_vy/d");
  fChargeTree->Branch("true_vz", &fTrue_vz, "true_vz/d");

  fChargeTree->Branch("true_vx_sce", &fTrue_vx_sce, "true_vx_sce/d");
  fChargeTree->Branch("true_vy_sce", &fTrue_vy_sce, "true_vy_sce/d");
  fChargeTree->Branch("true_vz_sce", &fTrue_vz_sce, "true_vz_sce/d");
  fChargeTree->Branch("distance_from_true", &fDistanceFromTrue, "distance_from_true/d");

  fChargeTree->Branch("start_x", &fStartx, "start_x/d");
  fChargeTree->Branch("start_y", &fStarty, "start_y/d");
  fChargeTree->Branch("start_z", &fStartz, "start_z/d");
  fChargeTree->Branch("end_x", &fEndx, "end_x/d");
  fChargeTree->Branch("end_y", &fEndy, "end_y/d");
  fChargeTree->Branch("end_z", &fEndz, "end_z/d");
  fChargeTree->Branch("dir_x", &fDirectionx, "dir_x/d");
  fChargeTree->Branch("dir_y", &fDirectiony, "dir_y/d");
  fChargeTree->Branch("dir_z", &fDirectionz, "dir_z/d");
  fChargeTree->Branch("angle_zxplane_pf", &fAngleZXplanePF, "angle_zxplane_pf/d");
  fChargeTree->Branch("n_hits", &fNhits, "n_hits/i");
  fChargeTree->Branch("n_hitsU", &fNhitsU, "n_hitsU/i");
  fChargeTree->Branch("n_hitsV", &fNhitsV, "n_hitsV/i");
  fChargeTree->Branch("n_hitsY", &fNhitsY, "n_hitsY/i");
  fChargeTree->Branch("n_clusters", &fNclusters, "n_clusters/i");
  fChargeTree->Branch("length", &fLength, "length/d");

  // reco-true matching information
  fChargeTree->Branch("matched_pdg_code", &fMatchedPdgCode, "matched_pdg_code/i");
  fChargeTree->Branch("matched_E", &fMatchedE, "matched_E/d");
  fChargeTree->Branch("matched_px", &fMatchedPx, "matched_px/d");
  fChargeTree->Branch("matched_py", &fMatchedPy, "matched_py/d");
  fChargeTree->Branch("matched_pz", &fMatchedPz, "matched_pz/d");
  fChargeTree->Branch("matched_vx", &fMatchedVx, "matched_vx/d");
  fChargeTree->Branch("matched_vy", &fMatchedVy, "matched_vy/d");
  fChargeTree->Branch("matched_vz", &fMatchedVz, "matched_vz/d");
  fChargeTree->Branch("matched_endx", &fMatchedEndx, "matched_endx/d");
  fChargeTree->Branch("matched_endy", &fMatchedEndy, "matched_endy/d");
  fChargeTree->Branch("matched_endz", &fMatchedEndz, "matched_endz/d");
  fChargeTree->Branch("distance_from_matched", &fDistanceFromMatched, "distance_from_matched/d");
  fChargeTree->Branch("mc_reco_costheta", &fMC_reco_costheta, "mc_reco_costheta/d");
  fChargeTree->Branch("fraction_matched_hits", &fFractionMatchedHits, "fraction_matched_hits/d");

  fChargeTree->Branch("matched_startX_shifted", &fMatchedStartX_shifted, "matched_startX_shifted/d");
  fChargeTree->Branch("matched_startY_shifted", &fMatchedStartY_shifted, "matched_startY_shifted/d");
  fChargeTree->Branch("matched_startZ_shifted", &fMatchedStartZ_shifted, "matched_startZ_shifted/d");
  fChargeTree->Branch("matched_startX_shifted_onlySCE", &fMatchedStartX_shifted_onlySCE, "matched_startX_shifted_onlySCE/d");
  fChargeTree->Branch("matched_startY_shifted_onlySCE", &fMatchedStartY_shifted_onlySCE, "matched_startY_shifted_onlySCE/d");
  fChargeTree->Branch("matched_startZ_shifted_onlySCE", &fMatchedStartZ_shifted_onlySCE, "matched_startZ_shifted_onlySCE/d");
  fChargeTree->Branch("matched_startX_shifted_onlyGeant", &fMatchedStartX_shifted_onlyGeant, "matched_startX_shifted_onlyGeant/d");
  fChargeTree->Branch("matched_startY_shifted_onlyGeant", &fMatchedStartY_shifted_onlyGeant, "matched_startY_shifted_onlyGeant/d");
  fChargeTree->Branch("matched_startZ_shifted_onlyGeant", &fMatchedStartZ_shifted_onlyGeant, "matched_startZ_shifted_onlyGeant/d");

  // dqdx information
  fChargeTree->Branch("reco_energy_U", &fReco_energy_U, "reco_energy_U/d");
  fChargeTree->Branch("reco_energy_V", &fReco_energy_V, "reco_energy_V/d");
  fChargeTree->Branch("reco_energy_Y", &fReco_energy_Y, "reco_energy_Y/d");

  fChargeTree->Branch("dQdx_hits_start", "std::vector<double>", &fDQdx_hits_start);
  fChargeTree->Branch("dQdx_wires_start", "std::vector<int>", &fDQdx_wires_start);
  // fChargeTree->Branch("dQdx_start", "std::vector<double>", &fDQdx_start);
  fChargeTree->Branch("dQdx_U_start", &fDQdx_U_start, "dQdx_U_start/d");
  fChargeTree->Branch("dQdx_V_start", &fDQdx_V_start, "dQdx_V_start/d");
  fChargeTree->Branch("dQdx_Y_start", &fDQdx_Y_start, "dQdx_Y_start/d");
  fChargeTree->Branch("n_hits_dQdx_U_start", &fn_hits_dQdx_U_start, "n_hits_dQdx_U_start/i");
  fChargeTree->Branch("n_hits_dQdx_V_start", &fn_hits_dQdx_V_start, "n_hits_dQdx_V_start/i");
  fChargeTree->Branch("n_hits_dQdx_Y_start", &fn_hits_dQdx_Y_start, "n_hits_dQdx_Y_start/i");
  fChargeTree->Branch("box_start_z_start", &fBox_start_z_start, "box_start_z_start/d");
  fChargeTree->Branch("box_start_x_start", &fBox_start_x_start, "box_start_x_start/d");
  fChargeTree->Branch("box_direction_z_start", &fBox_direction_z_start, "box_direction_z_start/d");
  fChargeTree->Branch("box_direction_x_start", &fBox_direction_x_start, "box_direction_x_start/d");
  fChargeTree->Branch("angle_ZXplane_cluster_start", &fAngleZXplaneCluster_start, "angle_ZXplane_cluster_start/d");
  fChargeTree->Branch("distance_starts_start", &fDistance_starts_start, "distance_starts_start/d");

  fChargeTree->Branch("distance_recotrue_starts_start", &fDistance_recotrue_starts_start, "distance_recotrue_starts_start/d");
  fChargeTree->Branch("distance_recotrue_shifted_starts_start", &fDistance_recotrue_shifted_starts_start, "distance_recotrue_shifted_starts_start/d");

  fChargeTree->Branch("mc_shower_start_x", &fMCShowerStartX, "mc_shower_start_x/d");
  fChargeTree->Branch("mc_shower_start_y", &fMCShowerStartY, "mc_shower_start_y/d");
  fChargeTree->Branch("mc_shower_start_z", &fMCShowerStartZ, "mc_shower_start_z/d");

  fChargeTree->Branch("mc_shower_startX_shifted", &fMCShowerStartX_shifted, "mc_shower_startX_shifted/d");
  fChargeTree->Branch("mc_shower_startY_shifted", &fMCShowerStartY_shifted, "mc_shower_startY_shifted/d");
  fChargeTree->Branch("mc_shower_startZ_shifted", &fMCShowerStartZ_shifted, "mc_shower_startZ_shifted/d");
  fChargeTree->Branch("mc_shower_startX_shifted_onlySCE", &fMCShowerStartX_shifted_onlySCE, "mc_shower_startX_shifted_onlySCE/d");
  fChargeTree->Branch("mc_shower_startY_shifted_onlySCE", &fMCShowerStartY_shifted_onlySCE, "mc_shower_startY_shifted_onlySCE/d");
  fChargeTree->Branch("mc_shower_startZ_shifted_onlySCE", &fMCShowerStartZ_shifted_onlySCE, "mc_shower_startZ_shifted_onlySCE/d");
  fChargeTree->Branch("mc_shower_startX_shifted_onlyGeant", &fMCShowerStartX_shifted_onlyGeant, "mc_shower_startX_shifted_onlyGeant/d");
  fChargeTree->Branch("mc_shower_startY_shifted_onlyGeant", &fMCShowerStartY_shifted_onlyGeant, "mc_shower_startY_shifted_onlyGeant/d");
  fChargeTree->Branch("mc_shower_startZ_shifted_onlyGeant", &fMCShowerStartZ_shifted_onlyGeant, "mc_shower_startZ_shifted_onlyGeant/d");

  fChargeTree->Branch("distance_recotrue_mc_shower_starts_start",
                      &fDistance_recotrue_mc_shower_starts_start,
                      "distance_recotrue_mc_shower_starts_start/d");
  fChargeTree->Branch("distance_recotrue_shifted_mc_shower_starts_start",
                      &fDistance_recotrue_shifted_mc_shower_starts_start,
                      "distance_recotrue_shifted_mc_shower_starts_start/d");

  fChargeTree->Branch("dQdx_hits_end", "std::vector<double>", &fDQdx_hits_end);
  fChargeTree->Branch("dQdx_wires_end", "std::vector<int>", &fDQdx_wires_end);
  // fChargeTree->Branch("dQdx_end", "std::vector<double>", &fDQdx_end);
  fChargeTree->Branch("dQdx_U_end", &fDQdx_U_end, "dQdx_U_end/d");
  fChargeTree->Branch("dQdx_V_end", &fDQdx_V_end, "dQdx_V_end/d");
  fChargeTree->Branch("dQdx_Y_end", &fDQdx_Y_end, "dQdx_Y_end/d");
  fChargeTree->Branch("n_hits_dQdx_U_end", &fn_hits_dQdx_U_end, "n_hits_dQdx_U_end/i");
  fChargeTree->Branch("n_hits_dQdx_V_end", &fn_hits_dQdx_V_end, "n_hits_dQdx_V_end/i");
  fChargeTree->Branch("n_hits_dQdx_Y_end", &fn_hits_dQdx_Y_end, "n_hits_dQdx_Y_end/i");
  fChargeTree->Branch("box_start_z_end", &fBox_start_z_end, "box_start_z_end/d");
  fChargeTree->Branch("box_start_x_end", &fBox_start_x_end, "box_start_x_end/d");
  fChargeTree->Branch("box_direction_z_end", &fBox_direction_z_end, "box_direction_z_end/d");
  fChargeTree->Branch("box_direction_x_end", &fBox_direction_x_end, "box_direction_x_end/d");
  fChargeTree->Branch("angle_ZXplane_cluster_end", &fAngleZXplaneCluster_end, "angle_ZXplane_cluster_end/d");
  fChargeTree->Branch("distance_starts_end", &fDistance_starts_end, "distance_starts_end/d");

  fChargeTree->Branch("pitch_U", &fPitch_U, "pitch_U/d");
  fChargeTree->Branch("pitch_V", &fPitch_V, "pitch_V/d");
  fChargeTree->Branch("pitch_Y", &fPitch_Y, "pitch_Y/d");
}

// double dqdxAnalyzer::GetChargeCorrection(int plane, double x, double y, double z)
// {
//   double x_correction, yz_correction, correction;
//   yz_correction = energyCalibProvider.YZdqdxCorrection(plane, y, z);
//   x_correction = energyCalibProvider.XdqdxCorrection(plane, x);
//   if (!yz_correction) yz_correction = 1.0;
//   if (!x_correction) x_correction = 1.0;
//   correction = yz_correction * x_correction;
//   return correction;
// }

void dqdxAnalyzer::recoTrueMatching(art::Event const &evt, art::Ptr<recob::PFParticle> const &pfparticle)
{
  bool _is_data = evt.isRealData();
  if (_is_data)
  {
    std::cout << "[RecoTrueMatcher] Running on a real data file. No MC-PFP matching will be attempted." << std::endl;
    fMatchedPdgCode = 1000000.;
    fMatchedE = 1000000.;
    fMatchedPx = 1000000.;
    fMatchedPy = 1000000.;
    fMatchedPz = 1000000.;
    fMatchedVx = 1000000.;
    fMatchedVy = 1000000.;
    fMatchedVz = 1000000.;
    fMatchedEndx = 1000000.;
    fMatchedEndy = 1000000.;
    fMatchedEndz = 1000000.;
    fFractionMatchedHits = -1;
    fMatchedTrackID = -1;
    fMatchedT = 1000000;
    fMatchedStartX_shifted = 1000000;
    fMatchedStartY_shifted = 1000000;
    fMatchedStartZ_shifted = 1000000;
    fMatchedStartX_shifted_onlySCE = 1000000;
    fMatchedStartY_shifted_onlySCE = 1000000;
    fMatchedStartZ_shifted_onlySCE = 1000000;
    fMatchedStartX_shifted_onlyGeant = 1000000;
    fMatchedStartY_shifted_onlyGeant = 1000000;
    fMatchedStartZ_shifted_onlyGeant = 1000000;
  }
  else
  {
    if (_use_premade_ass)
      _mcpfpMatcher.Configure(evt, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel, _mcpHitAssLabel, lar_pandora::LArPandoraHelper::kAddDaughters);
    else
      _mcpfpMatcher.Configure(evt, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel);

    // This is a map: PFParticle to matched MCParticle: std::map<art::Ptr<recob::PFParticle>, art::Ptr<simb::MCParticle> >
    lar_pandora::PFParticlesToMCParticles matched_pfp_to_mcp_map;
    lar_pandora::PFParticlesToFractionMCParticleHits matchedHitsFraction;
    _mcpfpMatcher.GetRecoToTrueMatches(matched_pfp_to_mcp_map, matchedHitsFraction);

    auto iter = matched_pfp_to_mcp_map.find(pfparticle);
    if (iter == matched_pfp_to_mcp_map.end())
    {
      fMatchedPdgCode = 1000000.;
      fMatchedE = 1000000.;
      fMatchedPx = 1000000.;
      fMatchedPy = 1000000.;
      fMatchedPz = 1000000.;
      fMatchedVx = 1000000.;
      fMatchedVy = 1000000.;
      fMatchedVz = 1000000.;
      fMatchedEndx = 1000000.;
      fMatchedEndy = 1000000.;
      fMatchedEndz = 1000000.;
      fFractionMatchedHits = -1;
      fMatchedTrackID = -1;
      fMatchedT = 1000000;
      fMatchedStartX_shifted = 1000000;
      fMatchedStartY_shifted = 1000000;
      fMatchedStartZ_shifted = 1000000;
      fMatchedStartX_shifted_onlySCE = 1000000;
      fMatchedStartY_shifted_onlySCE = 1000000;
      fMatchedStartZ_shifted_onlySCE = 1000000;
      fMatchedStartX_shifted_onlyGeant = 1000000;
      fMatchedStartY_shifted_onlyGeant = 1000000;
      fMatchedStartZ_shifted_onlyGeant = 1000000;
    }
    else
    {
      art::Ptr<simb::MCParticle> mc_part = iter->second;
      fMatchedPdgCode = mc_part->PdgCode();
      fMatchedE = mc_part->E();
      fMatchedPx = mc_part->Px();
      fMatchedPy = mc_part->Py();
      fMatchedPz = mc_part->Pz();
      fMatchedVx = mc_part->Vx();
      fMatchedVy = mc_part->Vy();
      fMatchedVz = mc_part->Vz();
      fMatchedEndx = mc_part->EndX();
      fMatchedEndy = mc_part->EndY();
      fMatchedEndz = mc_part->EndZ();
      double fStart_true[3] = {fMatchedVx, fMatchedVy, fMatchedVz};
      fDistanceFromMatched = geoHelper.distance(fStart_true, fTrue_v);

      auto iter2 = matchedHitsFraction.find(pfparticle);
      fFractionMatchedHits = iter2->second;

      fMatchedTrackID = mc_part->TrackId();
      fMatchedT = mc_part->T();

      double matched_start_shifted[3];
      shiftTruePosition(fStart_true, fMatchedT, matched_start_shifted);
      fMatchedStartX_shifted = matched_start_shifted[0];
      fMatchedStartY_shifted = matched_start_shifted[1];
      fMatchedStartZ_shifted = matched_start_shifted[2];

      double matched_start_shifted_onlySCE[3];
      shiftTruePosition_onlySCE(fStart_true, matched_start_shifted_onlySCE);
      fMatchedStartX_shifted_onlySCE = matched_start_shifted_onlySCE[0];
      fMatchedStartY_shifted_onlySCE = matched_start_shifted_onlySCE[1];
      fMatchedStartZ_shifted_onlySCE = matched_start_shifted[2];

      double matched_start_shifted_onlyGeant[3];
      shiftTruePosition_onlyGeant(fStart_true, fMatchedT, matched_start_shifted_onlyGeant);
      fMatchedStartX_shifted_onlyGeant = matched_start_shifted_onlyGeant[0];
      fMatchedStartY_shifted_onlyGeant = matched_start_shifted_onlyGeant[1];
      fMatchedStartZ_shifted_onlyGeant = matched_start_shifted_onlyGeant[2];
    }
  }
}

void dqdxAnalyzer::shiftTruePosition(double true_point[3], double true_point_shifted[3])
{
  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto offset = SCE->GetPosOffsets(true_point[0], true_point[1], true_point[2]);
  if (offset.size() == 3)
  {
    true_point_shifted[0] =
        true_point[0] - offset[0] + 0.7;
    true_point_shifted[1] =
        true_point[1] + offset[1];
    true_point_shifted[2] =
        true_point[2] + offset[2];
  }
}

void dqdxAnalyzer::shiftTruePosition(double true_point[3], double true_time, double true_point_shifted[3])
{
  true_point_shifted[0] = true_point[0];
  true_point_shifted[1] = true_point[1];
  true_point_shifted[2] = true_point[2];

  ::detinfo::DetectorProperties const* _detector_properties;
  ::detinfo::DetectorClocks const* _detector_clocks;

  _detector_properties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _detector_clocks = lar::providerFrom<detinfo::DetectorClocksService>();

  double g4Ticks = _detector_clocks->TPCG4Time2Tick(true_time)
                       + _detector_properties->GetXTicksOffset(0,0,0)
                       - _detector_properties->TriggerOffset();
  double xOffset = _detector_properties->ConvertTicksToX(g4Ticks, 0, 0, 0);

  true_point_shifted[0] += xOffset;

  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto offset = SCE->GetPosOffsets(true_point[0], true_point[1], true_point[2]);
  if (offset.size() == 3)
  {
    true_point_shifted[0] -= offset[0];
    true_point_shifted[1] += offset[1];
    true_point_shifted[2] += offset[2];
  }
}

void dqdxAnalyzer::shiftTruePosition_onlySCE(double true_point[3], double true_point_shifted[3])
{
  true_point_shifted[0] = true_point[0];
  true_point_shifted[1] = true_point[1];
  true_point_shifted[2] = true_point[2];

  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto offset = SCE->GetPosOffsets(true_point[0], true_point[1], true_point[2]);
  if (offset.size() == 3)
  {
    true_point_shifted[0] -= offset[0];
    true_point_shifted[1] += offset[1];
    true_point_shifted[2] += offset[2];
  }
}

void dqdxAnalyzer::shiftTruePosition_onlyGeant(double true_point[3], double true_time, double true_point_shifted[3])
{
  true_point_shifted[0] = true_point[0];
  true_point_shifted[1] = true_point[1];
  true_point_shifted[2] = true_point[2];

  ::detinfo::DetectorProperties const* _detector_properties;
  ::detinfo::DetectorClocks const* _detector_clocks;

  _detector_properties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _detector_clocks = lar::providerFrom<detinfo::DetectorClocksService>();

  double g4Ticks = _detector_clocks->TPCG4Time2Tick(true_time)
                       + _detector_properties->GetXTicksOffset(0,0,0)
                       - _detector_properties->TriggerOffset();
  double xOffset = _detector_properties->ConvertTicksToX(g4Ticks, 0, 0, 0);

  true_point_shifted[0] += xOffset;
}

void dqdxAnalyzer::clear()
{
  fDQdx_hits_start.clear();
  fDQdx_wires_start.clear();

  fDQdx_hits_end.clear();
  fDQdx_wires_end.clear();
  for (size_t i = 0; i < 3; i++)
  {
    fDQdx_start[i] = -1;
    fDQdx_end[i] = -1;
  }
}

void dqdxAnalyzer::trueNeutrinoInformation(art::Event const &evt)
{
  auto const &generator_handle =
      evt.getValidHandle<std::vector<simb::MCTruth>>(_mctruthLabel);
  auto const &generator(*generator_handle);

  // bool there_is_a_neutrino = false;
  //std::cout << "[PandoraLEEAnalyzer] Generator size " << generator.size() << std::endl;
  if (generator.size() > 0)
  {
    for (auto &gen : generator)
    {
      //std::cout << "[PandoraLEEAnalyzer] Generator origin " << gen.Origin() << std::endl;

      if (gen.Origin() == simb::kBeamNeutrino)
      {
        // there_is_a_neutrino = true;
        fTrue_nu_pdg = gen.GetNeutrino().Nu().PdgCode();
        fTrue_nu_energy = gen.GetNeutrino().Nu().E();
        fTrue_ccnc = gen.GetNeutrino().CCNC();
        // _qsqr = gen.GetNeutrino().QSqr();
        fTrue_theta = gen.GetNeutrino().Theta();

        // if (_ccnc == simb::kNC)
        // {
        //   fcategory = k_nc;
        // }

        fTrue_vx = gen.GetNeutrino().Nu().Vx();
        fTrue_vy = gen.GetNeutrino().Nu().Vy();
        fTrue_vz = gen.GetNeutrino().Nu().Vz();

        fTrue_v[0] = fTrue_vx;
        fTrue_v[1] = fTrue_vy;
        fTrue_v[2] = fTrue_vz;

        try
        {
          shiftTruePosition(fTrue_v, fTrue_v_sce);

          fTrue_vx_sce = fTrue_v_sce[0];
          fTrue_vy_sce = fTrue_v_sce[1];
          fTrue_vz_sce = fTrue_v_sce[2];
        }
        catch (...)
        {
          std::cout << "[PandoraLEEAnalyzer] "
                    << "Space Charge service offset size < 3" << std::endl;
          continue;
        }
      }
    }
  }
}

void dqdxAnalyzer::reconfigure(fhicl::ParameterSet const &p)
{
  _pfp_producer = p.get<std::string>("PFParticleProducer", "pandoraNu");
  _hitfinderLabel = p.get<std::string>("HitProducer", "pandoraCosmicHitRemoval");
  _geantModuleLabel = p.get<std::string>("GeantModule", "largeant");
  _spacepointLabel = p.get<std::string>("SpacePointProducer", "pandoraNu");
  _mcpHitAssLabel = p.get<std::string>("MCPHitAssProducer", "crHitRemovalTruthMatch");
  _use_premade_ass = p.get<bool>("UsePremadeMCPHitAss", true);
  _mctruthLabel = p.get<std::string>("MCTruthLabel", "generator");

  _dQdxRectangleWidth = p.get<double>("dQdxRectangleWidth", 1);
  _dQdxRectangleLength = p.get<double>("dQdxRectangleLength", 4);
  _dQdxRectangleLength_end = p.get<double>("dQdxRectangleLength_end", 8);
}

void dqdxAnalyzer::analyze(art::Event const &evt)
{
  clear();
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();
  // std::cout << "Event: " << fEvent << ", Run: " << fRun << ", subRun: " << fSubrun << std::endl;

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);

  std::vector<art::Ptr<recob::PFParticle>> pfp_v;
  art::fill_ptr_vector(pfp_v, pfparticle_handle);

  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  art::FindManyP<recob::Shower> showers_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Track> tracks_per_pfpart(pfparticle_handle, evt, _pfp_producer);

  art::FindManyP<recob::Vertex> vertices_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_cluster(cluster_handle, evt, _pfp_producer);
  auto const &spacepoint_handle =
      evt.getValidHandle<std::vector<recob::SpacePoint>>(_pfp_producer);

  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, _pfp_producer);

  // for (size_t i_pfp = 0; i_pfp < pfp_v.size(); i_pfp++)
  // {
  //   art::Ptr<recob::PFParticle> const pfparticle = pfp_v.at(i_pfp);
  //   if (pfparticle->Self() != i_pfp)
  //   {
  //     std::cout << "Self != pfparticle, " << pfparticle->Self() << " " <<  i_pfp << std::endl;
  //     std::cout << "Event: " << fEvent << ", Run: " << fRun << ", subRun: " << fSubrun << std::endl;
  //   }
  // }
  // std::cout << "event ok " << std::endl;
  trueNeutrinoInformation(evt);

  for (size_t i_pfp = 0; i_pfp < pfp_v.size(); i_pfp++)
  {
    art::Ptr<recob::PFParticle> const pfparticle = pfp_v.at(i_pfp);
    fPdgCode = pfparticle->PdgCode();

    if (fPdgCode != 11 && fPdgCode != 13)
    {
      continue;
    }

    try
    {
      double i_parent = pfparticle->Parent();
      fPdgCodeParent = pfp_v.at(i_parent)->PdgCode();
    }
    catch (...)
    {
      fPdgCodeParent = 1000000.;
    }

    // std::cout << "\n\nPfparticle:  " << i_pfp << ", pdgcode: " << fPdgCode << std::endl;
    // starting point and end point
    if (fPdgCode == 11)
    {
      std::vector<art::Ptr<recob::Shower>> pf_objs = showers_per_pfpart.at(i_pfp);
      // std::cout << "Found shower or tracks per pf part with lenght " << fNshowers_or_tracks << ", for pdg = " << fPdgCode << std::endl;
      if (pf_objs.size() != 0)
      {
        fStartx = pf_objs[0]->ShowerStart().X();
        fStarty = pf_objs[0]->ShowerStart().Y();
        fStartz = pf_objs[0]->ShowerStart().Z();
        fEndx = 1000000.;
        fEndy = 1000000.;
        fEndz = 1000000.;
        fLength = pf_objs[0]->Length();
        fDirectionx = pf_objs[0]->Direction().X();
        fDirectiony = pf_objs[0]->Direction().Y();
        fDirectionz = pf_objs[0]->Direction().Z();
        fAngleZXplanePF = atan2(fDirectionx, fDirectionz);
        // std::cout << "shower dir: " << pf_objs[0]->Direction().Z() << " , " << pf_objs[0]->Direction().X() << " , " << pf_objs[0]->Direction().X()/pf_objs[0]->Direction().Z() << std::endl;
      }
      else
      {
        continue;
      }
    }
    else if (fPdgCode == 13)
    {
      std::vector<art::Ptr<recob::Track>> pf_objs = tracks_per_pfpart.at(i_pfp);
      // std::cout << "Found shower or tracks per pf part with lenght " << fNshowers_or_tracks << ", for pdg = " << fPdgCode << std::endl;
      if (pf_objs.size() != 0)
      {
        fStartx = pf_objs[0]->Start().X();
        fStarty = pf_objs[0]->Start().Y();
        fStartz = pf_objs[0]->Start().Z();
        fEndx = pf_objs[0]->End().X();
        fEndy = pf_objs[0]->End().Y();
        fEndz = pf_objs[0]->End().Z();
        fLength = pf_objs[0]->Length();
        fDirectionx = pf_objs[0]->StartDirection().X();
        fDirectiony = pf_objs[0]->StartDirection().Y();
        fDirectionz = pf_objs[0]->StartDirection().Z();
        fAngleZXplanePF = atan2(fDirectionx, fDirectionz);
        // std::cout << "track dir: " << pf_objs[0]->StartDirection().Z() << " , " << pf_objs[0]->StartDirection().X() << " , " << pf_objs[0]->StartDirection().X()/pf_objs[0]->StartDirection().Z() << std::endl;
      }
      else
      {
        continue;
      }
    }
    else
    {
      fStartx = 1000000.;
      fStarty = 1000000.;
      fStartz = 1000000.;
      fEndx = 1000000.;
      fEndy = 1000000.;
      fEndz = 1000000.;
      fLength = 1000000.;
      fDirectionx = 1000000.;
      fDirectiony = 1000000.;
      fDirectionz = 1000000.;
    }

    // reco true matching
    recoTrueMatching(evt, pfparticle);
    double reco_dir[3] = {fDirectionx, fDirectiony, fDirectionz};
    double mc_dir[3] = {fMatchedPx, fMatchedPy, fMatchedPz};
    fMC_reco_costheta = geoHelper.costheta(reco_dir, mc_dir);

    // store pitch
    TVector3 t_reco_dir = TVector3(reco_dir);
    fPitch_U = geoHelper.getPitch(t_reco_dir, 0);
    fPitch_V = geoHelper.getPitch(t_reco_dir, 1);
    fPitch_Y = geoHelper.getPitch(t_reco_dir, 2);

    double fStart[3] = {fStartx, fStarty, fStartz};
    fDistanceFromTrue = geoHelper.distance(fStart, fTrue_v_sce);
    // std::cout << "Startz: " << fStartz <<  ", Startx:  " << fStartx << std::endl;
    // clusters
    fNhits = 0;
    fNhitsU = 0;
    fNhitsV = 0;
    fNhitsY = 0;
    fNclusters = 0;

    std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(i_pfp);
    fNclusters = clusters.size();
    for (art::Ptr<recob::Cluster> &cluster : clusters)
    {
      int aux_n_hits = cluster->NHits();
      fNhits += aux_n_hits;
      int fPlane = cluster->Plane().Plane;
      if (fPlane == 0)
      {
        fNhitsU += aux_n_hits;
      }
      else if (fPlane == 1)
      {
        fNhitsV += aux_n_hits;
      }
      else if (fPlane == 2)
      {
        fNhitsY += aux_n_hits;
      }
    }

    // reco energy
    std::vector<double> this_energy;
    std::vector<int> this_nhits;
    energyHelper.energyFromHits(*pfparticle, this_nhits, this_energy, evt, _pfp_producer);
    fReco_energy_U = this_energy[0];
    fReco_energy_V = this_energy[1];
    fReco_energy_Y = this_energy[2];

    // dqdx start
    // fDQdx_start = {-1., -1., -1.};
    double box_start_start[3][2], box_direction_start[3][2];
    for (size_t i = 0; i < 3; i++)
    {
      box_start_start[i][0] = 1000000;
      box_start_start[i][1] = 1000000;
      box_direction_start[i][0] = 1000000;
      box_direction_start[i][1] = 1000000;
    }
    std::vector<std::vector<double>> aux_dqdx_hits_start;
    aux_dqdx_hits_start.resize(3, std::vector<double>(0));
    std::vector<std::vector<int>> aux_dqdx_wires_start;
    aux_dqdx_wires_start.resize(3, std::vector<int>(0));
    std::string box_position_start = "start";
    energyHelper.dQdx(i_pfp,
                      evt,
                      fDQdx_start,
                      aux_dqdx_hits_start,
                      aux_dqdx_wires_start,
                      box_start_start,
                      box_direction_start,
                      box_position_start,
                      _dQdxRectangleLength,
                      _dQdxRectangleWidth,
                      _pfp_producer);
    fDQdx_U_start = fDQdx_start[0];
    fDQdx_V_start = fDQdx_start[1];
    fDQdx_Y_start = fDQdx_start[2];
    fn_hits_dQdx_U_start = aux_dqdx_hits_start[0].size();
    fn_hits_dQdx_V_start = aux_dqdx_hits_start[1].size();
    fn_hits_dQdx_Y_start = aux_dqdx_hits_start[2].size();
    fDQdx_hits_start = aux_dqdx_hits_start[2];
    fDQdx_wires_start = aux_dqdx_wires_start[2];
    fBox_start_z_start = box_start_start[2][0];
    fBox_start_x_start = box_start_start[2][1];
    fBox_direction_z_start = box_direction_start[2][0];
    fBox_direction_x_start = box_direction_start[2][1];
    fAngleZXplaneCluster_start = atan2(fBox_direction_x_start, fBox_direction_z_start);
    fDistance_starts_start = sqrt(pow((fBox_start_x_start - fStartx), 2) + pow((fBox_start_z_start - fStartz), 2));

    fDistance_recotrue_starts_start = sqrt(pow((fBox_start_x_start - fMatchedVx), 2) + pow((fBox_start_z_start - fMatchedVz), 2));
    double true_start[3] = {fMatchedVx, fMatchedVy, fMatchedVz};
    double true_start_shifted[3];
    shiftTruePosition(true_start, fMatchedT, true_start_shifted);
    fDistance_recotrue_shifted_starts_start = sqrt(pow((fBox_start_x_start - true_start_shifted[0]), 2) + pow((fBox_start_z_start - true_start_shifted[2]), 2));

    if (fMatchedPdgCode == 11 || fMatchedPdgCode == 22)
    {
      auto const &mc_shower_handle = evt.getValidHandle<std::vector<sim::MCShower>>("mcreco");

      for (size_t i_mc_shower = 0; i_mc_shower < mc_shower_handle->size(); i_mc_shower++)
      {
        sim::MCShower const mc_shower = mc_shower_handle->at(i_mc_shower);
        if (mc_shower.TrackID() != fMatchedTrackID)
        {
          continue;
        }
        else
        {
          sim::MCStep mc_step_shower_start = mc_shower.DetProfile();
          fMCShowerStartX = mc_step_shower_start.X();
          fMCShowerStartY = mc_step_shower_start.Y();
          fMCShowerStartZ = mc_step_shower_start.Z();
          double mc_shower_start[3] = {fMCShowerStartX, fMCShowerStartY, fMCShowerStartZ};

          double mc_shower_start_shifted[3];
          shiftTruePosition(mc_shower_start, fMatchedT, mc_shower_start_shifted);
          fMCShowerStartX_shifted = mc_shower_start_shifted[0];
          fMCShowerStartY_shifted = mc_shower_start_shifted[1];
          fMCShowerStartZ_shifted = mc_shower_start_shifted[2];

          double mc_shower_start_shifted_onlySCE[3];
          shiftTruePosition_onlySCE(mc_shower_start, mc_shower_start_shifted_onlySCE);
          fMCShowerStartX_shifted_onlySCE = mc_shower_start_shifted_onlySCE[0];
          fMCShowerStartY_shifted_onlySCE = mc_shower_start_shifted_onlySCE[1];
          fMCShowerStartZ_shifted_onlySCE = mc_shower_start_shifted[2];

          double mc_shower_start_shifted_onlyGeant[3];
          shiftTruePosition_onlyGeant(mc_shower_start, fMatchedT, mc_shower_start_shifted_onlyGeant);
          fMCShowerStartX_shifted_onlyGeant = mc_shower_start_shifted_onlyGeant[0];
          fMCShowerStartY_shifted_onlyGeant = mc_shower_start_shifted_onlyGeant[1];
          fMCShowerStartZ_shifted_onlyGeant = mc_shower_start_shifted_onlyGeant[2];
          
          fDistance_recotrue_mc_shower_starts_start =
              sqrt(pow((fBox_start_x_start - fMCShowerStartX), 2) + pow((fBox_start_z_start - fMCShowerStartZ), 2));
          fDistance_recotrue_shifted_mc_shower_starts_start =
              sqrt(pow((fBox_start_x_start - mc_shower_start_shifted[0]), 2) + pow((fBox_start_z_start - mc_shower_start_shifted[2]), 2));
          break;
        }
      }
    }
    else
    {
      fMCShowerStartX = 1000000;
      fMCShowerStartY = 1000000;
      fMCShowerStartZ = 1000000;
      fDistance_recotrue_mc_shower_starts_start = 1000000;
      fDistance_recotrue_shifted_mc_shower_starts_start = 1000000;
      fMCShowerStartX_shifted = 1000000;
      fMCShowerStartY_shifted = 1000000;
      fMCShowerStartZ_shifted = 1000000;
      fMCShowerStartX_shifted_onlySCE = 1000000;
      fMCShowerStartY_shifted_onlySCE = 1000000;
      fMCShowerStartZ_shifted_onlySCE = 1000000;
      fMCShowerStartX_shifted_onlyGeant = 1000000;
      fMCShowerStartY_shifted_onlyGeant = 1000000;
      fMCShowerStartZ_shifted_onlyGeant = 1000000;
    }

    // dqdx end
    // fDQdx_end = {-1., -1., -1.};
    double box_start_end[3][2], box_direction_end[3][2];
    for (size_t i = 0; i < 3; i++)
    {
      box_start_end[i][0] = 1000000;
      box_start_end[i][1] = 1000000;
      box_direction_end[i][0] = 1000000;
      box_direction_end[i][1] = 1000000;
    }
    std::vector<std::vector<double>> aux_dqdx_hits_end;
    aux_dqdx_hits_end.resize(3, std::vector<double>(0));
    std::vector<std::vector<int>> aux_dqdx_wires_end;
    aux_dqdx_wires_end.resize(3, std::vector<int>(0));
    std::string box_position_end = "end";
    energyHelper.dQdx(i_pfp,
                      evt,
                      fDQdx_end,
                      aux_dqdx_hits_end,
                      aux_dqdx_wires_end,
                      box_start_end,
                      box_direction_end,
                      box_position_end,
                      _dQdxRectangleLength_end,
                      _dQdxRectangleWidth,
                      _pfp_producer);
    fDQdx_U_end = fDQdx_end[0];
    fDQdx_V_end = fDQdx_end[1];
    fDQdx_Y_end = fDQdx_end[2];
    fn_hits_dQdx_U_end = aux_dqdx_hits_end[0].size();
    fn_hits_dQdx_V_end = aux_dqdx_hits_end[1].size();
    fn_hits_dQdx_Y_end = aux_dqdx_hits_end[2].size();
    fDQdx_hits_end = aux_dqdx_hits_end[2];
    fDQdx_wires_end = aux_dqdx_wires_end[2];
    fBox_start_z_end = box_start_end[2][0];
    fBox_start_x_end = box_start_end[2][1];
    fBox_direction_z_end = box_direction_end[2][0];
    fBox_direction_x_end = box_direction_end[2][1];
    fAngleZXplaneCluster_end = atan2(fBox_direction_x_end, fBox_direction_z_end);
    fDistance_starts_end = sqrt(pow((fBox_start_x_end - fStartx), 2) + pow((fBox_start_z_end - fStartz), 2));

    fChargeTree->Fill();
    clear();
  }
}

DEFINE_ART_MODULE(dqdxAnalyzer)
