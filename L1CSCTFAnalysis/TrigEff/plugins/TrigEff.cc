#include "TrigEff.h"
#include <map>
#include <vector>
#include <iostream>
#include <typeinfo>
#include "DataFormats/Common/interface/Ref.h"
#include <memory>

// Probably most of this includes are not useful or redundant
// clean them in another life
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/HcalIsolatedTrack/interface/IsolatedPixelTrackCandidate.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include <L1Trigger/CSCTrackFinder/src/CSCTFDTReceiverLUT.h>

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <cassert>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Lindsey's Objects
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonCandidateTrackFwd.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonCandidateTrack.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrackFwd.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrack.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonRegionalTracksFwd.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTrackSeed.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitiveFwd.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitive.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"


using namespace trigger;
using namespace reco;
using namespace std;
using namespace edm;

TrigEff::TrigEff(const edm::ParameterSet& pset):edm::EDAnalyzer(),
                                                _matchBox(4,4,4,0,true)
{
  L1extraTag = pset.getParameter<InputTag>("L1extraTag");
  muonsTag   = pset.getParameter<InputTag>("muonsTag");
  outputFile = pset.getParameter<string>("outputFile");

  csctfTag      = pset.getParameter<InputTag>("csctfTag");
  csctfLctsTag  = pset.getParameter<InputTag>("csctfLctsTag");
  dtSegTag      = pset.getParameter<InputTag>("dtSegTag");
  cscSegTag     = pset.getParameter<InputTag>("cscSegTag");


  // L1GTUtils
  //m_nameAlgTechTrig = pset.getParameter<std::string> ("AlgorithmName");
    
  // printLevel
  //  0 - No printout
  //  1 - Basic info printout
  //  2 - Verbose
  printLevel = pset.getUntrackedParameter<int>("printLevel",0);
  
  _matchBox . setPrintLevel( printLevel );

  file  = new TFile(outputFile.c_str(),"RECREATE");

  // These collections are added after Ivan's email about we need
  // to do trigger efficiency. They contain all the info about 
  // reco'd muons: global, standalone and tracker ones
  // reco Muons Collection
  recoMuons = new TTree("recoMuons", "recoMuons");
  
  //---------------------------------------------------------------------
  // general information Booking
  //--------------------------------------------------------------------- 
  recoMuons->Branch("Run",   &Run,   "Run/I"  );
  recoMuons->Branch("Event", &Event, "Event/I");
  recoMuons->Branch("Lumi",  &Lumi,  "Lumi/I" );
  recoMuons->Branch("Bx",    &Bx,    "Bx/I"   );
  recoMuons->Branch("Orbit", &Orbit, "Orbit/I");
  
  recoMuons->Branch("muonSize", &muonSize, "muonSize/I");
   
  recoMuons->Branch("isGlobalMuon"           , &isGlobalMuon       );
  recoMuons->Branch("isTrackerMuon"          , &isTrackerMuon      );
  recoMuons->Branch("isStandAloneMuon"       , &isStandAloneMuon   );
  recoMuons->Branch("isMuonAllArbitrated"    , &isMuonAllArbitrated);
  recoMuons->Branch("isTMLastStationAngTight", &isTMLastStationAngTight);
  recoMuons->Branch("isGlobalMuonPromptTight", &isGlobalMuonPromptTight);

  recoMuons->Branch("isEnergyValid"    , &isEnergyValid);
  recoMuons->Branch("caloCompatibility", &caloCompatibility);
  recoMuons->Branch("em",    &em    ); 
  recoMuons->Branch("emS9",  &emS9  ); 
  recoMuons->Branch("emS25", &emS25 );
  recoMuons->Branch("emMax", &emMax );
  recoMuons->Branch("had",   &had   );
  recoMuons->Branch("hadS9", &hadS9 );
  recoMuons->Branch("hadMax",&hadMax);
  recoMuons->Branch("ho",    &ho    );
  recoMuons->Branch("hoS9",  &hoS9  );
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // Global Muon Block Booking
  //---------------------------------------------------------------------
  //recoMuons->Branch("gmrEnergy"           , &gmrEnergy           );
  //recoMuons->Branch("gmrDEnergy"          , &gmrDEnergy          );
  recoMuons->Branch("gmrPt"               , &gmrPt               );
  recoMuons->Branch("gmrEta"              , &gmrEta              );
  recoMuons->Branch("gmrPhi"              , &gmrPhi              );
  recoMuons->Branch("gmrP"                , &gmrP                );
  recoMuons->Branch("gmrPx"               , &gmrPx               );
  recoMuons->Branch("gmrPy"               , &gmrPy               );
  recoMuons->Branch("gmrPz"               , &gmrPz               );
  recoMuons->Branch("gmrTheta"            , &gmrTheta            );
  recoMuons->Branch("gmrVx"               , &gmrVx               );
  recoMuons->Branch("gmrVy"               , &gmrVy               );
  recoMuons->Branch("gmrVz"               , &gmrVz               );
  recoMuons->Branch("gmrCharge"           , &gmrCharge           );
  recoMuons->Branch("gmrNDoF"             , &gmrNDoF             );
  recoMuons->Branch("gmrChi2"             , &gmrChi2             );
  recoMuons->Branch("gmrChi2Norm"         , &gmrChi2Norm         );
  recoMuons->Branch("gmrDXY"              , &gmrDXY              );
  recoMuons->Branch("gmrDTheta"           , &gmrDTheta           );
  recoMuons->Branch("gmrDPt"              , &gmrDPt              );
  recoMuons->Branch("gmrDEta"             , &gmrDEta             );
  recoMuons->Branch("gmrDPhi"             , &gmrDPhi             );
  recoMuons->Branch("gmrDDXY"             , &gmrDDXY             );
  recoMuons->Branch("gmrIso03nTracks"     , &gmrIso03nTracks     );
  recoMuons->Branch("gmrIso03sumPt"       , &gmrIso03sumPt       );
  recoMuons->Branch("gmrDz"               , &gmrDz               );
  recoMuons->Branch("gmrD0"               , &gmrD0               );
  recoMuons->Branch("gmrDsz"              , &gmrDsz              );
  recoMuons->Branch("gmrDDz"              , &gmrDDz              );
  recoMuons->Branch("gmrDD0"              , &gmrDD0              );
  recoMuons->Branch("gmrDDsz"             , &gmrDDsz             );
  recoMuons->Branch("gmrInnerX"           , &gmrInnerX           );
  recoMuons->Branch("gmrInnerY"           , &gmrInnerY           );
  recoMuons->Branch("gmrInnerZ"           , &gmrInnerZ           );
  recoMuons->Branch("gmrOuterX"           , &gmrOuterX           );
  recoMuons->Branch("gmrOuterY"           , &gmrOuterY           );
  recoMuons->Branch("gmrOuterZ"           , &gmrOuterZ           );
  recoMuons->Branch("gmrValHits"          , &gmrValHits          );

  //---------------------------------------------------------------------
  // Standalone Muon Block Booking
  //---------------------------------------------------------------------
  //recoMuons->Branch("stdEnergy"   , &stdEnergy   );
  //recoMuons->Branch("stdDEnergy"  , &stdDEnergy  );
  recoMuons->Branch("stdPt"       , &stdPt       );
  recoMuons->Branch("stdEta"      , &stdEta      );
  recoMuons->Branch("stdPhi"      , &stdPhi      );
  recoMuons->Branch("stdPx"       , &stdPx       );
  recoMuons->Branch("stdPy"       , &stdPy       );
  recoMuons->Branch("stdPz"       , &stdPz       );
  recoMuons->Branch("stdVx"       , &stdVx       );
  recoMuons->Branch("stdVy"       , &stdVy       );
  recoMuons->Branch("stdVz"       , &stdVz       );
  recoMuons->Branch("stdCharge"   , &stdCharge   );
  recoMuons->Branch("stdDPt"      , &stdDPt      );
  recoMuons->Branch("stdDEta"     , &stdDEta     );
  recoMuons->Branch("stdDPhi"     , &stdDPhi     );
  recoMuons->Branch("stdDz"       , &stdDz       );
  recoMuons->Branch("stdD0"       , &stdD0       );
  recoMuons->Branch("stdNDoF"     , &stdNDoF     ); 
  recoMuons->Branch("stdChi2"     , &stdChi2     );
  recoMuons->Branch("stdChi2Norm" , &stdChi2Norm );
  recoMuons->Branch("stdDXY"      , &stdDXY      );
  recoMuons->Branch("stdTheta"    , &stdTheta    );
  recoMuons->Branch("stdDTheta"   , &stdDTheta   );
  recoMuons->Branch("stdDDz"      , &stdDDz      );
  recoMuons->Branch("stdDD0"      , &stdDD0      );     
  recoMuons->Branch("stdValHits"  , &stdValHits  );
   
  //---------------------------------------------------------------------
  // Tracker Muon Block Booking
  //---------------------------------------------------------------------
  //recoMuons->Branch("trkEnergy"   , &trkEnergy   );
  //recoMuons->Branch("trkDEnergy"  , &trkDEnergy  );
  recoMuons->Branch("trkPt"       , &trkPt       );
  recoMuons->Branch("trkEta"      , &trkEta      );
  recoMuons->Branch("trkPhi"      , &trkPhi      );
  recoMuons->Branch("trkPx"       , &trkPx       );
  recoMuons->Branch("trkPy"       , &trkPy       );
  recoMuons->Branch("trkPz"       , &trkPz       );
  recoMuons->Branch("trkVx"       , &trkVx       );
  recoMuons->Branch("trkVy"       , &trkVy       );
  recoMuons->Branch("trkVz"       , &trkVz       );
  recoMuons->Branch("trkCharge"   , &trkCharge   );
  recoMuons->Branch("trkDPt"      , &trkDPt      );
  recoMuons->Branch("trkDEta"     , &trkDEta     );
  recoMuons->Branch("trkDPhi"     , &trkDPhi     );
  recoMuons->Branch("trkDz"       , &trkDz       );
  recoMuons->Branch("trkD0"       , &trkD0       );
  recoMuons->Branch("trkNDoF"     , &trkNDoF     ); 
  recoMuons->Branch("trkChi2"     , &trkChi2     );
  recoMuons->Branch("trkChi2Norm" , &trkChi2Norm );
  recoMuons->Branch("trkDXY"      , &trkDXY      );
  recoMuons->Branch("trkTheta"    , &trkTheta    );
  recoMuons->Branch("trkDTheta"   , &trkDTheta   );
  recoMuons->Branch("trkDDz"      , &trkDDz      );
  recoMuons->Branch("trkDD0"      , &trkDD0      );     
  recoMuons->Branch("trkValHits"  , &trkValHits  );

  // CSC segment for the tracker muon
  //chamber
  recoMuons->Branch("trkNchambers"    , &trkNchambers    );
  recoMuons->Branch("trkNofMatches"   , &trkNofMatches   );
  recoMuons->Branch("trkIsMatchValid" , &trkIsMatchValid );
  
  recoMuons->Branch("trkNSegs"           , &trkNSegs);

  recoMuons->Branch("trkSegChamberId"    , trkSegChamberId    ,"trkSegChamberId[muonSize][100]/I");
  recoMuons->Branch("trkSegRing"         , trkSegRing         ,"trkSegRing[muonSize][100]/I");
  recoMuons->Branch("trkSegStation"      , trkSegStation      ,"trkSegStation[muonSize][100]/I");
  recoMuons->Branch("trkSegEndcap"       , trkSegEndcap       ,"trkSegEndcap[muonSize][100]/I");
  recoMuons->Branch("trkSegTriggerSector", trkSegTriggerSector,"trkSegTriggerSector[muonSize][100]/I");
  recoMuons->Branch("trkSegTriggerCscId" , trkSegTriggerCscId ,"trkSegTriggerCscId[muonSize][100]/I");
  recoMuons->Branch("trkSegXfromMatch"   , trkSegXfromMatch   ,"trkSegXfromMatch[muonSize][100]/F");
  recoMuons->Branch("trkSegYfromMatch"   , trkSegYfromMatch   ,"trkSegYfromMatch[muonSize][100]/F");
  recoMuons->Branch("trkSegZfromMatch"   , trkSegZfromMatch   ,"trkSegZfromMatch[muonSize][100]/F");
  recoMuons->Branch("trkSegRfromMatch"   , trkSegRfromMatch   ,"trkSegRfromMatch[muonSize][100]/F");
  recoMuons->Branch("trkSegPhifromMatch" , trkSegPhifromMatch ,"trkSegPhifromMatch[muonSize][100]/F");
  recoMuons->Branch("trkSegEtafromMatch" , trkSegEtafromMatch ,"trkSegEtafromMatch[muonSize][100]/F");

  //segment
  recoMuons->Branch("trkSegIsArb", &trkSegIsArb, "trkSegIsArb[muonSize][100]/I");
  recoMuons->Branch("trkSegX"    , &trkSegX    , "trkSegX[muonSize][100]/F");
  recoMuons->Branch("trkSegY"    , &trkSegY    , "trkSegY[muonSize][100]/F");
  recoMuons->Branch("trkSegZ"    , &trkSegZ    , "trkSegZ[muonSize][100]/F");
  recoMuons->Branch("trkSegR"    , &trkSegR    , "trkSegR[muonSize][100]/F");
  recoMuons->Branch("trkSegPhi"  , &trkSegPhi  , "trkSegPhi[muonSize][100]/F");
  recoMuons->Branch("trkSegEta"  , &trkSegEta  , "trkSegEta[muonSize][100]/F");
  recoMuons->Branch("trkSegDxDz"    , &trkSegDxDz    , "trkSegDxDz[muonSize][100]/F");
  recoMuons->Branch("trkSegDyDz"    , &trkSegDyDz    , "trkSegDyDz[muonSize][100]/F");
  recoMuons->Branch("trkSegDxDzErr" , &trkSegDxDzErr , "trkSegDxDzErr[muonSize][100]/F");
  recoMuons->Branch("trkSegDyDzErr" , &trkSegDyDzErr , "trkSegDyDzErr[muonSize][100]/F");

  //---------------------------------------------------------------------
  // RECHIT information: only for standalone/global muons!
  //---------------------------------------------------------------------
  recoMuons->Branch("rchCSCtype" , &rchCSCtype ); 
  recoMuons->Branch("rchEtaLocal", &rchEtaLocal); 
  recoMuons->Branch("rchPhiLocal", &rchPhiLocal); 
  recoMuons->Branch("rchEta"     , &rchEta     ); 
  recoMuons->Branch("rchPhi"     , &rchPhi     ); 
  recoMuons->Branch("rchPhi_02PI", &rchPhi_02PI);

  recoMuons->Branch("rchStation", &rchStation);
  recoMuons->Branch("rchChamber", &rchChamber);
  recoMuons->Branch("rchRing"   , &rchRing   );
  recoMuons->Branch("rchLayer"  , &rchLayer  );

  recoMuons->Branch("rchMuonSize",      &rchMuonSize     );
 
  recoMuons->Branch("rchEtaMatrixLocal", rchEtaMatrixLocal, "rchEtaMatrixLocal[muonSize][35]/F");
  recoMuons->Branch("rchPhiMatrixLocal", rchPhiMatrixLocal, "rchPhiMatrixLocal[muonSize][35]/F");
  recoMuons->Branch("rchEtaMatrix"    , rchEtaMatrix    , "rchEtaMatrix[muonSize][35]/F");
  recoMuons->Branch("rchPhiMatrix"    , rchPhiMatrix    , "rchPhiMatrix[muonSize][35]/F");
  recoMuons->Branch("rchPhi02PIMatrix", rchPhi02PIMatrix, "rchPhi02PIMatrix[muonSize][35]/F");
  recoMuons->Branch("rchStationMatrix", rchStationMatrix, "rchStationMatrix[muonSize][35]/I");
  recoMuons->Branch("rchChamberMatrix", rchChamberMatrix, "rchChamberMatrix[muonSize][35]/I");
  recoMuons->Branch("rchRingMatrix"   , rchRingMatrix   , "rchRingMatrix[muonSize][35]/I");
  recoMuons->Branch("rchLayerMatrix"  , rchLayerMatrix  , "rchLayerMatrix[muonSize][35]/I");
  recoMuons->Branch("rchTypeMatrix"   , rchTypeMatrix   , "rchTypeMatrix[muonSize][35]/I");
   
  //--------------------------------------------------------------------- 
  // old format: keep it until you are sure the TMatrixF works
  recoMuons->Branch("nMu_nCscRchHits", &nMu_nCscRchHits, "nMu_nCscRchHits/I"); 
  recoMuons->Branch("rchEtaList",       rchEtaList,      "rchEtaList[nMu_nCscRchHits]/D"); 
  recoMuons->Branch("rchPhiList",       rchPhiList,      "rchPhiList[nMu_nCscRchHits]/D");
  recoMuons->Branch("rchPhiList_02PI",  rchPhiList_02PI, "rchPhiList_02PI[nMu_nCscRchHits]/D");
  resizeRchHits(1);  
  //--------------------------------------------------------------------- 

  //---------------------------------------------------------------------
  // Propagation block booking
  //---------------------------------------------------------------------
  // propagation to ME1/1
  recoMuons->Branch(  "muons_x_me11",  &muons_x_me11);
  recoMuons->Branch(  "muons_y_me11",  &muons_y_me11);
  recoMuons->Branch(  "muons_z_me11",  &muons_z_me11);
  recoMuons->Branch("muons_phi_me11",&muons_phi_me11);
  recoMuons->Branch("muons_eta_me11",&muons_eta_me11);

  // propagation to ME1
  recoMuons->Branch(  "muons_x_me1",  &muons_x_me1);
  recoMuons->Branch(  "muons_y_me1",  &muons_y_me1);
  recoMuons->Branch(  "muons_z_me1",  &muons_z_me1);
  recoMuons->Branch("muons_phi_me1",&muons_phi_me1);
  recoMuons->Branch("muons_eta_me1",&muons_eta_me1);

  // propagation to ME2
  recoMuons->Branch(  "muons_x_me2",  &muons_x_me2);
  recoMuons->Branch(  "muons_y_me2",  &muons_y_me2);
  recoMuons->Branch(  "muons_z_me2",  &muons_z_me2);
  recoMuons->Branch("muons_phi_me2",&muons_phi_me2);
  recoMuons->Branch("muons_eta_me2",&muons_eta_me2);

  // propagation to ME3
  recoMuons->Branch(  "muons_x_me3",  &muons_x_me3);
  recoMuons->Branch(  "muons_y_me3",  &muons_y_me3);
  recoMuons->Branch(  "muons_z_me3",  &muons_z_me3);
  recoMuons->Branch("muons_phi_me3",&muons_phi_me3);
  recoMuons->Branch("muons_eta_me3",&muons_eta_me3);

  //---------------------------------------------------------------------
  // all segment information
  //---------------------------------------------------------------------
  recoMuons->Branch("segsSize"     ,&segsSize     ,"segsSize/I");
  recoMuons->Branch("cscsegs_loc_x",&cscsegs_loc_x);
  recoMuons->Branch("cscsegs_loc_y",&cscsegs_loc_y);
  recoMuons->Branch("cscsegs_loc_z",&cscsegs_loc_z);

  recoMuons->Branch("cscsegs_loc_theta",&cscsegs_loc_theta);
  recoMuons->Branch("cscsegs_loc_eta",  &cscsegs_loc_eta);
  recoMuons->Branch("cscsegs_loc_phi",  &cscsegs_loc_phi);

  recoMuons->Branch("cscsegs_loc_dir_theta",&cscsegs_loc_dir_theta);
  recoMuons->Branch("cscsegs_loc_dir_eta",  &cscsegs_loc_dir_eta);
  recoMuons->Branch("cscsegs_loc_dir_phi",  &cscsegs_loc_dir_phi);

  recoMuons->Branch("cscsegs_gbl_x",&cscsegs_gbl_x);
  recoMuons->Branch("cscsegs_gbl_y",&cscsegs_gbl_y);
  recoMuons->Branch("cscsegs_gbl_z",&cscsegs_gbl_z);

  recoMuons->Branch("cscsegs_gbl_theta",&cscsegs_gbl_theta);
  recoMuons->Branch("cscsegs_gbl_eta",  &cscsegs_gbl_eta);
  recoMuons->Branch("cscsegs_gbl_phi",  &cscsegs_gbl_phi);

  recoMuons->Branch("cscsegs_gbl_dir_theta",&cscsegs_gbl_dir_theta);
  recoMuons->Branch("cscsegs_gbl_dir_eta",  &cscsegs_gbl_dir_eta);
  recoMuons->Branch("cscsegs_gbl_dir_phi",  &cscsegs_gbl_dir_phi);

  recoMuons->Branch("cscsegs_sector" ,&cscsegs_sector );
  recoMuons->Branch("cscsegs_endcap" ,&cscsegs_endcap );
  recoMuons->Branch("cscsegs_station",&cscsegs_station);
  recoMuons->Branch("cscsegs_ring"   ,&cscsegs_ring   );
  recoMuons->Branch("cscsegs_chamber",&cscsegs_chamber);

  recoMuons->Branch("isdtseg" , &isdtseg);
  recoMuons->Branch("cscsegs_wheel" ,&cscsegs_wheel);

  //---------------------------------------------------------------------
  // segments belonging to the muon
  //---------------------------------------------------------------------
  recoMuons->Branch("muonNsegs",&muonNsegs);

  recoMuons->Branch("muon_cscsegs_loc_x"      , muon_cscsegs_loc_x      ,"muon_cscsegs_loc_x[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_loc_y"      , muon_cscsegs_loc_y      ,"muon_cscsegs_loc_y[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_loc_eta"    , muon_cscsegs_loc_eta    ,"muon_cscsegs_loc_eta[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_loc_phi"    , muon_cscsegs_loc_phi    ,"muon_cscsegs_loc_phi[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_loc_dir_eta", muon_cscsegs_loc_dir_eta,"muon_cscsegs_loc_dir_eta[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_loc_dir_phi", muon_cscsegs_loc_dir_phi,"muon_cscsegs_loc_dir_phi[muonSize][16]/F");

  recoMuons->Branch("muon_cscsegs_gbl_x"      , muon_cscsegs_gbl_x      ,"muon_cscsegs_gbl_x[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_gbl_y"      , muon_cscsegs_gbl_y      ,"muon_cscsegs_gbl_y[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_gbl_eta"    , muon_cscsegs_gbl_eta    ,"muon_cscsegs_gbl_eta[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_gbl_phi"    , muon_cscsegs_gbl_phi    ,"muon_cscsegs_gbl_phi[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_gbl_dir_eta", muon_cscsegs_gbl_dir_eta,"muon_cscsegs_gbl_dir_eta[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_gbl_dir_phi", muon_cscsegs_gbl_dir_phi,"muon_cscsegs_gbl_dir_phi[muonSize][16]/F");

  recoMuons->Branch("muon_cscsegs_dxdz"   , muon_cscsegs_dxdz   ,"muon_cscsegs_dxdz[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_dydz"   , muon_cscsegs_dydz   ,"muon_cscsegs_dydz[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_dxdzErr", muon_cscsegs_dxdzErr,"muon_cscsegs_dxdzErr[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_dydzErr", muon_cscsegs_dydzErr,"muon_cscsegs_dydzErr[muonSize][16]/F");

  recoMuons->Branch("muon_cscsegs_endcap" , muon_cscsegs_endcap ,"muon_cscsegs_endcap[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_station", muon_cscsegs_station,"muon_cscsegs_station[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_ring"   , muon_cscsegs_ring   ,"muon_cscsegs_ring[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_chamber", muon_cscsegs_chamber,"muon_cscsegs_chamber[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_nhits"  , muon_cscsegs_nhits  ,"muon_cscsegs_nhits[muonSize][16]/I");
  
  recoMuons->Branch("muon_cscsegs_islctable", muon_cscsegs_islctable,"muon_cscsegs_islctable[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_ismatched"   , muon_cscsegs_ismatched,   "muon_cscsegs_ismatched[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_lctId"       , muon_cscsegs_lctId    ,   "muon_cscsegs_lctId[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_nmatched" , muon_cscsegs_nmatched ,"muon_cscsegs_nmatched[muonSize][16]/I");

  recoMuons->Branch("muon_isdtseg"          , muon_isdtseg,            "muon_isdtseg[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_wheel"    , muon_cscsegs_wheel,     "muon_cscsegs_wheel[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_sector"   , muon_cscsegs_sector,    "muon_cscsegs_sector[muonSize][16]/I");

  //---------------------------------------------------------------------
  // l1 extra muon collection booking
  //---------------------------------------------------------------------
  l1extraMuons = new TTree("l1extraMuons", "l1extraMuons");
  l1extraMuons->Branch("l1Size", &l1Size, "l1Size/I");
    
  l1extraMuons->Branch("l1Eta", &l1Eta);
  l1extraMuons->Branch("l1Pt" , &l1Pt );
  l1extraMuons->Branch("l1Phi", &l1Phi);
   
  l1extraMuons->Branch("isIsolated"  , &isIsolated  );
  l1extraMuons->Branch("isMip"       , &isMip       );
  l1extraMuons->Branch("isForward"   , &isForward   );
  l1extraMuons->Branch("isRPC"       , &isRPC       );
  l1extraMuons->Branch("detectorType", &detectorType);
  l1extraMuons->Branch("rank"        , &rank        );
  //---------------------------------------------------------------------

  csctfTTree = new TTree("csctfTTree","csctfTTree");
  // csctf
  csctfTTree->Branch("SizeTrk"       , &SizeTrk,      "SizeTrk/I");
  csctfTTree->Branch("isdtlct"       , &isdtlct      );
  csctfTTree->Branch("EndcapTrk"     , &EndcapTrk     );
  csctfTTree->Branch("SectorTrk"     , &SectorTrk     );
  csctfTTree->Branch("BxTrk"  	     , &BxTrk         );
  csctfTTree->Branch("me1ID"  	     , &me1ID         );
  csctfTTree->Branch("me2ID"         , &me2ID 	      );
  csctfTTree->Branch("me3ID"  	     , &me3ID         );
  csctfTTree->Branch("me4ID"  	     , &me4ID         );
  csctfTTree->Branch("mb1ID"  	     , &mb1ID         );
  csctfTTree->Branch("OutputLinkTrk" , &OutputLinkTrk );
  csctfTTree->Branch("ModeTrk"       , &ModeTrk       );
  csctfTTree->Branch("EtaTrk"  	     , &EtaTrk        );
  csctfTTree->Branch("PhiTrk"  	     , &PhiTrk        );
  csctfTTree->Branch("PhiTrk_02PI"   , &PhiTrk_02PI   );
  csctfTTree->Branch("PtTrk"  	     , &PtTrk         );
  csctfTTree->Branch("ChargeTrk"     , &ChargeTrk     );
  csctfTTree->Branch("ChargeValidTrk", &ChargeValidTrk);
  csctfTTree->Branch("QualityTrk"    , &QualityTrk    );
  csctfTTree->Branch("ForRTrk"       , &ForRTrk       );
  csctfTTree->Branch("Phi23Trk"      , &Phi23Trk      );
  csctfTTree->Branch("Phi12Trk"      , &Phi12Trk      );
  csctfTTree->Branch("PhiSignTrk"    , &PhiSignTrk    );
  csctfTTree->Branch("EtaBitTrk"     , &EtaBitTrk     );
  csctfTTree->Branch("PhiBitTrk"     , &PhiBitTrk     );
  csctfTTree->Branch("PtBitTrk"      , &PtBitTrk      );

  csctfTTree->Branch("NumLCTsTrk"       , &NumLCTsTrk       );
  csctfTTree->Branch("trLctEndcap"      , trLctEndcap      , "trLctEndcap[SizeTrk][4]/I");
  csctfTTree->Branch("trLctSector"      , trLctSector      , "trLctSector[SizeTrk][4]/I");
  csctfTTree->Branch("trLctSubSector"   , trLctSubSector   , "trLctSubSector[SizeTrk][4]/I");
  csctfTTree->Branch("trLctBx"          , trLctBx          , "trLctBx[SizeTrk][4]/I");
  csctfTTree->Branch("trLctBx0"         , trLctBx0         , "trLctBx0[SizeTrk][4]/I");
  csctfTTree->Branch("trLctStation"     , trLctStation     , "trLctStation[SizeTrk][4]/I");
  csctfTTree->Branch("trLctRing"        , trLctRing        , "trLctRing[SizeTrk][4]/I");
  csctfTTree->Branch("trLctChamber"     , trLctChamber     , "trLctChamber[SizeTrk][4]/I");
  csctfTTree->Branch("trLctTriggerCSCID", trLctTriggerCSCID, "trLctTriggerCSCID[SizeTrk][4]/I");
  csctfTTree->Branch("trLctFpga"        , trLctFpga        , "trLctFpga[SizeTrk][4]/I");                                                                                     
  csctfTTree->Branch("trLctlocalPhi"    , trLctlocalPhi    , "trLctlocalPhi[SizeTrk][4]/I");
  csctfTTree->Branch("trLctglobalPhi"   , trLctglobalPhi   , "trLctglobalPhi[SizeTrk][4]/I");
  csctfTTree->Branch("trLctglobalEta"   , trLctglobalEta   , "trLctglobalEta[SizeTrk][4]/I");                                                                                  
  csctfTTree->Branch("trLctstripNum"    , trLctstripNum    , "trLctstripNum[SizeTrk][4]/I");
  csctfTTree->Branch("trLctwireGroup"   , trLctwireGroup   , "trLctwireGroup[SizeTrk][4]/I");
  
  csctfTTree->Branch("isdttrlct"           , &isdttrlct          , "dt_trLctEndcap[SizeTrk][4]");
  csctfTTree->Branch("dt_NumLCTsTrk"       , &dt_NumLCTsTrk       );
  csctfTTree->Branch("dt_trLctEndcap"      , dt_trLctEndcap      , "dt_trLctEndcap[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctSector"      , dt_trLctSector      , "dt_trLctSector[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctSubSector"   , dt_trLctSubSector   , "dt_trLctSubSector[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctBx"          , dt_trLctBx          , "dt_trLctBx[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctBx0"         , dt_trLctBx0         , "dt_trLctBx0[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctStation"     , dt_trLctStation     , "dt_trLctStation[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctRing"        , dt_trLctRing        , "dt_trLctRing[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctChamber"     , dt_trLctChamber     , "dt_trLctChamber[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctTriggerCSCID", dt_trLctTriggerCSCID, "dt_trLctTriggerCSCID[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctFpga"        , dt_trLctFpga        , "dt_trLctFpga[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctlocalPhi"    , dt_trLctlocalPhi    , "dt_trLctlocalPhi[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctglobalPhi"   , dt_trLctglobalPhi   , "dt_trLctglobalPhi[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctglobalEta"   , dt_trLctglobalEta   , "dt_trLctglobalEta[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctstripNum"    , dt_trLctstripNum    , "dt_trLctstripNum[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctwireGroup"   , dt_trLctwireGroup   , "dt_trLctwireGroup[SizeTrk][4]");

  // all lcts
  csctfTTree->Branch("SizeLCTs"       , &SizeLCTs       ,"SizeLCTs/I");
  csctfTTree->Branch("lctEndcap"      , &lctEndcap      );
  csctfTTree->Branch("lctSector"      , &lctSector      );
  csctfTTree->Branch("lctSubSector"   , &lctSubSector   );
  csctfTTree->Branch("lctBx"          , &lctBx          );
  csctfTTree->Branch("lctBx0"         , &lctBx0         );
  csctfTTree->Branch("lctStation"     , &lctStation     );
  csctfTTree->Branch("lctRing"        , &lctRing        );
  csctfTTree->Branch("lctChamber"     , &lctChamber     );
  csctfTTree->Branch("lctTriggerCSCID", &lctTriggerCSCID);
  csctfTTree->Branch("lctFpga"        , &lctFpga        );
  csctfTTree->Branch("lctlocalPhi"    , &lctlocalPhi    );
  csctfTTree->Branch("lctglobalPhi"   , &lctglobalPhi   );
  csctfTTree->Branch("lctglobalEta"   , &lctglobalEta   );
  csctfTTree->Branch("lctstripNum"    , &lctstripNum    );
  csctfTTree->Branch("lctwireGroup"   , &lctwireGroup   );
  csctfTTree->Branch("lctWheel"       , &lctWheel       );
  csctfTTree->Branch("trig_prim_phi"   , &trig_prim_phi );
  csctfTTree->Branch("trig_prim_eta"   , &trig_prim_eta );
  
  csctfTTree->Branch("dt_SizeLCTs"       , &dt_SizeLCTs       ,"dt_SizeLCTs/I");
  csctfTTree->Branch("dt_lctEndcap"      , &dt_lctEndcap      );
  csctfTTree->Branch("dt_lctSector"      , &dt_lctSector      );
  csctfTTree->Branch("dt_lctSubSector"   , &dt_lctSubSector   );
  csctfTTree->Branch("dt_lctBx"          , &dt_lctBx          );
  csctfTTree->Branch("dt_lctBx0"         , &dt_lctBx0         );
  csctfTTree->Branch("dt_lctStation"     , &dt_lctStation     );
  csctfTTree->Branch("dt_lctRing"        , &dt_lctRing        );
  csctfTTree->Branch("dt_lctChamber"     , &dt_lctChamber     );
  csctfTTree->Branch("dt_lctTriggerCSCID", &dt_lctTriggerCSCID);
  csctfTTree->Branch("dt_lctFpga"        , &dt_lctFpga        );
  csctfTTree->Branch("dt_lctlocalPhi"    , &dt_lctlocalPhi    );
  csctfTTree->Branch("dt_lctglobalPhi"   , &dt_lctglobalPhi   );
  csctfTTree->Branch("dt_lctglobalEta"   , &dt_lctglobalEta   );
  csctfTTree->Branch("dt_lctstripNum"    , &dt_lctstripNum    );
  csctfTTree->Branch("dt_lctwireGroup"   , &dt_lctwireGroup   );
    
  // -----------------------------------------------------------------------------
  
  bzero(srLUTs_ , sizeof(srLUTs_));
  int sector=1;    // assume SR LUTs are all same for every sector
  bool TMB07=true; // specific TMB firmware
  // Create a pset for SR/PT LUTs: if you do not change the value in the 
  // configuration file, it will load the default minitLUTs
  edm::ParameterSet srLUTset;
  srLUTset.addUntrackedParameter<bool>("ReadLUTs", false);
  srLUTset.addUntrackedParameter<bool>("Binary",   false);
  srLUTset.addUntrackedParameter<std::string>("LUTPath", "./");
  
  // positive endcap
  int endcap = 1; 
  for(int station=1,fpga=0; station<=4 && fpga<5; station++)
    {
      if(station==1)
	for(int subSector=0; subSector<2 && fpga<5; subSector++)
	  srLUTs_[fpga++][1] = new CSCSectorReceiverLUT(endcap,sector,subSector+1,
							station, srLUTset, TMB07);
      else
	srLUTs_[fpga++][1]   = new CSCSectorReceiverLUT(endcap,  sector,   0, 
							station, srLUTset, TMB07);
    }
  
  // negative endcap
  endcap = 2; 
  for(int station=1,fpga=0; station<=4 && fpga<5; station++)
    {
      if(station==1)
	for(int subSector=0; subSector<2 && fpga<5; subSector++)
	  srLUTs_[fpga++][0] = new CSCSectorReceiverLUT(endcap,sector,subSector+1,
							station, srLUTset, TMB07);
      else
	srLUTs_[fpga++][0]   = new CSCSectorReceiverLUT(endcap,  sector,   0, 
							station, srLUTset, TMB07);
    }
  // -----------------------------------------------------------------------------
  
  edm::ParameterSet serviceParameters = pset.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);

}

void TrigEff::endJob() {
  //free the CSCTF array of pointers
  for(int j=0; j<2; j++) 
    for(int i=0; i<5; i++) 
      delete srLUTs_[i][j]; 
  
  // delete ts;
  // delete tpts;
}

// destructor
TrigEff::~TrigEff(void){ 
  file->Write(); 
  file->Close(); 
}

// analyze
void TrigEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //Get the Magnetic field from the setup
  iSetup.get<IdealMagneticFieldRecord>().get(theBField);
  // Get the GlobalTrackingGeometry from the setup
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
  theService->update(iSetup);

  // Get the propagators
  iSetup.get<TrackingComponentsRecord>().get("SmartPropagatorAnyRK", propagatorAlong   );
  iSetup.get<TrackingComponentsRecord>().get("SmartPropagatorAnyOpposite", propagatorOpposite);


  // General event info
  Run   = iEvent.id().run();
  Event = iEvent.id().event();
  Bx    = iEvent.bunchCrossing();   
  Lumi  = iEvent.luminosityBlock();
  Orbit = iEvent.orbitNumber();     
  
  // Get the CSC Geometry 
  iSetup.get<MuonGeometryRecord>().get(cscGeom);
  iSetup.get<MuonGeometryRecord>().get(dtGeom);
 
  // ============================================================================
  // C S C T F   R A W   B L O C K
  if (printLevel > 1 ) {
    cout << "\n\n================================== NEW EVENT =====================================\n";
    cout << " Event # " << Event << endl;
    cout << " Run   # " << Run << endl; 
    
    cout << "\n===================== FILLING CSCTF ========================\n"
	 <<   "============================================================\n";
      }



  edm::Handle< vector<L1TMuon::InternalTrack> > CSCTFtracks;
  if ( csctfTag.label() != "null" ) {
    iEvent.getByLabel(csctfTag, CSCTFtracks);
  }

  if( CSCTFtracks.isValid() ){
    csctfInit();

    if( iSetup.get< L1MuTriggerScalesRcd > ().cacheIdentifier() != m_scalesCacheID ||
        iSetup.get< L1MuTriggerPtScaleRcd >().cacheIdentifier() != m_ptScaleCacheID ){
      
      ESHandle< L1MuTriggerScales > scales;
      iSetup.get< L1MuTriggerScalesRcd >().get(scales);
      ts = scales.product();
      ESHandle< L1MuTriggerPtScale > ptscales;
      iSetup.get< L1MuTriggerPtScaleRcd >().get(ptscales);
      tpts = ptscales.product();
      m_scalesCacheID  = iSetup.get< L1MuTriggerScalesRcd  >().cacheIdentifier();
      m_ptScaleCacheID = iSetup.get< L1MuTriggerPtScaleRcd >().cacheIdentifier();
      
      //std::cout  << "Changing triggerscales and triggerptscales...\n";
    }    
    
    //fill the csctf information
      fillCSCTF(CSCTFtracks,
		//DTTFtracks,
		ts,
		tpts,
		srLUTs_
		);
    
  
    if (printLevel > 1 )
      cout<< "\n==================== FILLING ALL LCTS ======================\n"
	  <<   "============================================================\n"
	  << endl;
    
    Handle< vector<L1TMuon::TriggerPrimitive> > CSCTFlcts;
    
    if( csctfLctsTag.label() != "null" ) {
      iEvent.getByLabel(csctfLctsTag, CSCTFlcts);
    }

    if( CSCTFlcts.isValid() ) {
      
	fillAllLCTs( CSCTFlcts,
		   srLUTs_
		   );
    }    
    else cout << "Invalid CSCCorrelatedLCTDigiCollection... skipping it\n";
    
    
    // fill the ttree
    csctfTTree->Fill();

    //clean-up the pointers
    csctfDel();
  }
  else cout << "Invalid L1CSCTrackCollection... skipping it\n";
  
  
  // ============================================================================
  // R E C O   M U O N   B L O C K
  // Second: get global muons, stand alone muons, and track
  Handle<MuonCollection>  muons;
  if( muonsTag.label() != "null" ) iEvent.getByLabel(muonsTag, muons);
  
  // Check if the MuonCollection is valid  
  if( muons.isValid() ){
    
    // initialize the variables to dummy values
    muonsInit();
    
    //fill the muons
    if (printLevel >= 1) 
      cout<< "\n===================== FILLING MUONS ========================\n"
	  <<   "============================================================\n"
	  << endl;
    
    fillMuons(muons);
    
    // ==========================================================================
    //
    // look at SEGMENTs (from the CSC Validation, A. Kubik)
    //
    // ==========================================================================
    // get CSC segment collection
    edm::Handle<CSCSegmentCollection> cscSegments;
    edm::Handle<DTRecSegment4DCollection> dtSegments;
    //cscSegTag = cms.InputTag("cscSegments"),
    if( cscSegTag.label() != "null" ) iEvent.getByLabel(cscSegTag, cscSegments);
    if( dtSegTag.label()  != "null" ) iEvent.getByLabel(dtSegTag, dtSegments);
    if( cscSegments.isValid()) {
     
      //fill the segments variables
      if (printLevel >= 1) 
        cout<< "\n==================== FILLING SGMTS ========================\n"
	    <<   "===========================================================\n\n";
      
      fillSegments(cscSegments, dtSegments, cscGeom, dtGeom);
      
      edm::Handle< std::vector<L1TMuon::TriggerPrimitive> > CSCTFlcts;
      edm::Handle< std::vector<L1TMuon::InternalTrack> > CSCTFtracks;
      if( csctfLctsTag.label() != "null" ) iEvent.getByLabel(csctfLctsTag, CSCTFlcts);
      if( csctfTag.label() != "null" ) iEvent.getByLabel(csctfTag, CSCTFtracks);
      

      if( cscSegments.isValid() && CSCTFlcts.isValid() && CSCTFtracks.isValid() ) {
        if (printLevel >= 1) 
          cout<< "\n============== FILLING MUONS/SGMTS/CSCTF ==================\n"
	      <<   "===========================================================\n\n";
	
	fillSegmentsMuons(muons, cscSegments, dtSegments, cscGeom, CSCTFlcts, CSCTFtracks);
      }
    }
    
    recoMuons->Fill();
    muonsDel();
  }
  else {
    cout << "Invalid MuonCollection with label "
	 << muonsTag.label() 
	 << ": please check the root file with edmDumpEventContent "
	 << " to verify the collection with such label exists\n";
  }
  
  
   
  // ============================================================================
  // L 1 E X T R A 
  Handle<l1extra::L1MuonParticleCollection> l1muons;
  if( L1extraTag.label() != "null" ) iEvent.getByLabel(L1extraTag, l1muons);

  // Fetch all Level-1 muon candidates
  if( l1muons.isValid() ){
    
    if (printLevel > 2)
      std::cout<< "\n=================== FILLING L1EXTRA =======================\n"
               <<   "===========================================================\n\n";

    l1Size = l1muons->size();
    l1extraInit();
    
    //fill the l1 extra muons
    for(l1extra::L1MuonParticleCollection::const_iterator iter=l1muons->begin(); iter!=l1muons->end(); iter++) {
      
      // where the filling is happening
      fillExtra(iter);
    }
    
    l1extraMuons->Fill();
    l1extraDel();
  } 
  else{
    std::cout << "Invalid L1MuonParticleCollection with label "
              << L1extraTag.label() 
              << ": please check the root file with edmDumpEventContent "
              << " to verify the collection with such label exists\n";
  }
  
  
}


//-------------------------------------------------------------------------
// methods for the reco muons collection
void TrigEff::fillMuons(const edm::Handle<reco::MuonCollection> muons)
{
  
  if (printLevel > 2) {
    printf("%s: %d\n\n", "muons->size()", int(muons->size()) );
  }
  
  muonSize = muons->size();
  
  // Once you know how many muons we have in the event,
  // we resize the rchEtaList and rchPhiList
  nMu_nCscRchHits = muonSize * MAX_CSC_RECHIT;
  if (printLevel > 2) 
    std::cout << "nMu_nCscRchHits: " << nMu_nCscRchHits << std::endl;
  resizeRchHits(nMu_nCscRchHits);
  
    int whichMuon =-1;
  for (MuonCollection::const_iterator muon=muons->begin();
       muon!=muons->end(); muon++) {
    
    whichMuon++;

    if (printLevel > 2 ) {
      printf("************************************************\n");
      printf("M U O N  #%d\n", whichMuon);
      printf("************************************************\n\n");
      printf("%s\n"    , "--------------------------------");
      printf("%s: %d\n", "isGlobalMuon    ()", muon->isGlobalMuon    ());
      printf("%s: %d\n", "isTrackerMuon   ()", muon->isTrackerMuon   ());
      printf("%s: %d\n", "isStandAloneMuon()", muon->isStandAloneMuon());
      printf("%s: %d\n", "combinedMuon    ().isNonnull()", muon->combinedMuon  ().isNonnull());
      printf("%s: %d\n", "track           ().isNonnull()", muon->track         ().isNonnull());
      printf("%s: %d\n", "standAloneMuon  ().isNonnull()", muon->standAloneMuon().isNonnull());
      printf("%s\n\n"  , "--------------------------------");
    }
    
    // old left over from a previous ntuplizer
    if (muon->isIsolationValid()) {
      MuonIsolation isoR03 = muon->isolationR03();
      gmrIso03nTracks->push_back(isoR03.nTracks);
      gmrIso03sumPt  ->push_back(isoR03.sumPt  );
    }
    else {
      gmrIso03nTracks->push_back(-999);
      gmrIso03sumPt  ->push_back(-999);
    }
    //---------------------------------------------

    isGlobalMuon        -> push_back(muon->isGlobalMuon    ()); 
    isTrackerMuon       -> push_back(muon->isTrackerMuon   ()); 
    isStandAloneMuon    -> push_back(muon->isStandAloneMuon());
    int isAllArb = muon::isGoodMuon(*muon, muon::AllArbitrated);
    isMuonAllArbitrated -> push_back(isAllArb);
    int isTMLSAT = muon::isGoodMuon(*muon, muon::TMLastStationAngTight);
    isTMLastStationAngTight -> push_back(isTMLSAT);
    int isGBLPT = muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight);
    isGlobalMuonPromptTight -> push_back(isGBLPT);

    isEnergyValid       -> push_back(muon->isEnergyValid());
    caloCompatibility   -> push_back(muon->caloCompatibility());
    em                  -> push_back(muon->calEnergy().em    ); 
    emS9                -> push_back(muon->calEnergy().emS9  );
    emS25               -> push_back(muon->calEnergy().emS25 );
    emMax               -> push_back(muon->calEnergy().emMax );
    had                 -> push_back(muon->calEnergy().had   );
    hadS9               -> push_back(muon->calEnergy().hadS9 );
    hadMax              -> push_back(muon->calEnergy().hadMax);
    ho                  -> push_back(muon->calEnergy().ho    );
    hoS9                -> push_back(muon->calEnergy().hoS9  );

    // global muon
    if (muon->combinedMuon().isNonnull()) {

      TrackRef trackRef = muon->combinedMuon();

      if (printLevel > 1 ) 
        printf("(GBL) muon->pt(): %10.5f, muon->eta(): %10.5f, muon->phi(): %10.5f\n",
               trackRef->pt(), trackRef->eta(), trackRef->phi());
      
      gmrPt      -> push_back(trackRef->pt    ());
      gmrPhi     -> push_back(trackRef->phi   ());
      gmrEta     -> push_back(trackRef->eta   ());
      gmrP       -> push_back(trackRef->p     ());
      gmrPx      -> push_back(trackRef->px    ());
      gmrPy      -> push_back(trackRef->py    ());
      gmrPz      -> push_back(trackRef->pz    ());
      gmrTheta   -> push_back(trackRef->theta ());
      gmrCharge  -> push_back(trackRef->charge());
      gmrVx      -> push_back(trackRef->vx    ());
      gmrVy      -> push_back(trackRef->vy    ());
      gmrVz      -> push_back(trackRef->vz    ());
      gmrDXY     -> push_back(trackRef->dxy              ());
      gmrDDXY    -> push_back(trackRef->dxyError         ());
      gmrDPhi    -> push_back(trackRef->phiError         ());
      gmrDEta    -> push_back(trackRef->etaError         ());
      gmrDPt     -> push_back(trackRef->ptError          ());
      gmrChi2    -> push_back(trackRef->chi2             ());
      gmrNDoF    -> push_back(trackRef->ndof             ());
      gmrChi2Norm-> push_back(trackRef->normalizedChi2   ());
      gmrDTheta  -> push_back(trackRef->thetaError       ());
      gmrDz      -> push_back(trackRef->dz               ());
      gmrD0      -> push_back(trackRef->d0               ());
      gmrDsz     -> push_back(trackRef->dsz              ());
      gmrDDz     -> push_back(trackRef->dzError          ());
      gmrDD0     -> push_back(trackRef->d0Error          ());
      gmrDDsz    -> push_back(trackRef->dszError         ());
      gmrValHits -> push_back(trackRef->numberOfValidHits());

      gmrInnerX  -> push_back(trackRef->innerPosition().X());
      gmrInnerY  -> push_back(trackRef->innerPosition().Y());
      gmrInnerZ  -> push_back(trackRef->innerPosition().Z());   
  
      gmrOuterX  -> push_back(trackRef->outerPosition().X());
      gmrOuterY  -> push_back(trackRef->outerPosition().Y());
      gmrOuterZ  -> push_back(trackRef->outerPosition().Z());   

    }
    else {
      // fill with dummy values
      gmrPt      -> push_back(-999);
      gmrPhi     -> push_back(-999);
      gmrEta     -> push_back(-999);
      gmrP       -> push_back(-999);
      gmrPx      -> push_back(-999);
      gmrPy      -> push_back(-999);
      gmrPz      -> push_back(-999);
      gmrTheta   -> push_back(-999);
      gmrCharge  -> push_back(-999);
      gmrVx      -> push_back(-999);
      gmrVy      -> push_back(-999);
      gmrVz      -> push_back(-999);
      gmrDXY     -> push_back(-999);
      gmrDDXY    -> push_back(-999);
      gmrDPhi    -> push_back(-999);
      gmrDEta    -> push_back(-999);
      gmrDPt     -> push_back(-999);
      gmrChi2    -> push_back(-999);
      gmrNDoF    -> push_back(-999);
      gmrChi2Norm-> push_back(-999);
      gmrDTheta  -> push_back(-999);
      gmrDz      -> push_back(-999);
      gmrD0      -> push_back(-999);
      gmrDsz     -> push_back(-999);
      gmrDDz     -> push_back(-999);
      gmrDD0     -> push_back(-999);
      gmrDDsz    -> push_back(-999);
      gmrValHits -> push_back(-999);

      gmrInnerX  -> push_back(-999);
      gmrInnerY  -> push_back(-999);
      gmrInnerZ  -> push_back(-999);   
                                  
      gmrOuterX  -> push_back(-999);
      gmrOuterY  -> push_back(-999);
      gmrOuterZ  -> push_back(-999);   
    }

    
    // standalone muon
    if (muon->standAloneMuon().isNonnull()) {
      TrackRef trackRefStd = muon->standAloneMuon();

      if (printLevel > 1 ) 
        printf("(STA) muon->pt(): %10.5f, muon->eta(): %10.5f, muon->phi(): %10.5f\n",
               trackRefStd->pt(), trackRefStd->eta(), trackRefStd->phi());
    
      stdPt      -> push_back(trackRefStd->pt       ());
      stdEta     -> push_back(trackRefStd->eta      ());
      stdPhi     -> push_back(trackRefStd->phi      ());
      stdPx      -> push_back(trackRefStd->px       ());
      stdPy      -> push_back(trackRefStd->py       ());
      stdPz      -> push_back(trackRefStd->pz       ());
      stdVx      -> push_back(trackRefStd->vx       ());
      stdVy      -> push_back(trackRefStd->vy       ());
      stdVz      -> push_back(trackRefStd->vz       ());
      stdCharge  -> push_back(trackRefStd->charge   ());
      stdDPt     -> push_back(trackRefStd->ptError  ());
      stdDEta    -> push_back(trackRefStd->etaError ());
      stdDPhi    -> push_back(trackRefStd->phiError ());
      stdDz      -> push_back(trackRefStd->dz       ());
      stdD0      -> push_back(trackRefStd->d0       ());

      stdChi2    -> push_back(trackRefStd->chi2             ());
      stdNDoF    -> push_back(trackRefStd->ndof             ());
      stdChi2Norm-> push_back(trackRefStd->normalizedChi2   ());
      stdTheta   -> push_back(trackRefStd->theta            ());
      stdDTheta  -> push_back(trackRefStd->thetaError       ());
      stdDDz     -> push_back(trackRefStd->dzError          ());
      stdDD0     -> push_back(trackRefStd->d0Error          ());
      stdValHits -> push_back(trackRefStd->numberOfValidHits());

      stdInnerX  -> push_back(trackRefStd->innerPosition().X());
      stdInnerY  -> push_back(trackRefStd->innerPosition().Y());
      stdInnerZ  -> push_back(trackRefStd->innerPosition().Z());   

      stdOuterX  -> push_back(trackRefStd->outerPosition().X());
      stdOuterY  -> push_back(trackRefStd->outerPosition().Y());
      stdOuterZ  -> push_back(trackRefStd->outerPosition().Z());   
    }
    else {
      // fill with dummy values
      stdPt      -> push_back(-999);
      stdEta     -> push_back(-999);
      stdPhi     -> push_back(-999);
      stdPx      -> push_back(-999);
      stdPy      -> push_back(-999);
      stdPz      -> push_back(-999);
      stdVx      -> push_back(-999);
      stdVy      -> push_back(-999);
      stdVz      -> push_back(-999);
      stdCharge  -> push_back(-999);
      stdDPt     -> push_back(-999);
      stdDEta    -> push_back(-999);
      stdDPhi    -> push_back(-999);
      stdDz      -> push_back(-999);
      stdD0      -> push_back(-999);

      stdChi2    -> push_back(-999);
      stdNDoF    -> push_back(-999);
      stdChi2Norm-> push_back(-999);
      stdTheta   -> push_back(-999);
      stdDTheta  -> push_back(-999);
      stdDDz     -> push_back(-999);
      stdDD0     -> push_back(-999);
      stdValHits -> push_back(-999);

      stdInnerX  -> push_back(-999);
      stdInnerY  -> push_back(-999);
      stdInnerZ  -> push_back(-999);   
                              
      stdOuterX  -> push_back(-999);
      stdOuterY  -> push_back(-999);
      stdOuterZ  -> push_back(-999);   
      
    }
    
    //tracker muon    
    if (muon->track().isNonnull()) {
      
      TrackRef trackRefTrk = muon->track();

      if (printLevel > 1 ) 
        printf("(TRK) muon->pt(): %10.5f, muon->eta(): %10.5f, muon->phi(): %10.5f\n",
               trackRefTrk->pt(), trackRefTrk->eta(), trackRefTrk->phi());
      
      trkPt     ->push_back(trackRefTrk->pt       ());
      trkEta    ->push_back(trackRefTrk->eta      ());
      trkPhi    ->push_back(trackRefTrk->phi      ());
      trkPx     ->push_back(trackRefTrk->px       ());
      trkPy     ->push_back(trackRefTrk->py       ());
      trkPz     ->push_back(trackRefTrk->pz       ());
      trkVx     ->push_back(trackRefTrk->vx       ());
      trkVy     ->push_back(trackRefTrk->vy       ());
      trkVz     ->push_back(trackRefTrk->vz       ());
      trkCharge ->push_back(trackRefTrk->charge   ());
      trkDPt    ->push_back(trackRefTrk->ptError  ());
      trkDEta   ->push_back(trackRefTrk->etaError ());
      trkDPhi   ->push_back(trackRefTrk->phiError ());
      trkDz     ->push_back(trackRefTrk->dz       ());
      trkD0     ->push_back(trackRefTrk->d0       ());
      
      trkChi2    ->push_back(trackRefTrk->chi2             ());
      trkNDoF    ->push_back(trackRefTrk->ndof             ());
      trkChi2Norm->push_back(trackRefTrk->normalizedChi2   ());
      trkTheta   ->push_back(trackRefTrk->theta            ());
      trkDTheta  ->push_back(trackRefTrk->thetaError       ());
      trkDDz     ->push_back(trackRefTrk->dzError          ());
      trkDD0     ->push_back(trackRefTrk->d0Error          ());
      trkValHits ->push_back(trackRefTrk->numberOfValidHits());
    }
    else {
      // fill with dummy values
      trkPt     ->push_back(-999);
      trkEta    ->push_back(-999);
      trkPhi    ->push_back(-999);
      trkPx     ->push_back(-999);
      trkPy     ->push_back(-999);
      trkPz     ->push_back(-999);
      trkVx     ->push_back(-999);
      trkVy     ->push_back(-999);
      trkVz     ->push_back(-999);
      trkCharge ->push_back(-999);
      trkDPt    ->push_back(-999);
      trkDEta   ->push_back(-999);
      trkDPhi   ->push_back(-999);
      trkDz     ->push_back(-999);
      trkD0     ->push_back(-999);

      trkChi2    ->push_back(-999);
      trkNDoF    ->push_back(-999);
      trkChi2Norm->push_back(-999);
      trkTheta   ->push_back(-999);
      trkDTheta  ->push_back(-999);
      trkDDz     ->push_back(-999);
      trkDD0     ->push_back(-999);
      trkValHits ->push_back(-999);
    }


    // --------------------------------------------------------------------------
    // P R O P A G A T I O N   T O   T H E   E N D C A P
    // --------------------------------------------------------------------------
    if (printLevel > 2) 
      cout<< "\n\n"
	  <<"*******************************************\n"
	  <<"Projecting the track to the ME+X/Y surfaces\n"
	  <<"*******************************************\n";
    
    TrajectoryStateOnSurface tsos;

    // which TrackRef to use?
    TrackRef trackRef;
    // get the tracker muon
    if (muon->combinedMuon  ().isNonnull() ||
        muon->track         ().isNonnull()  )    trackRef = muon->track();
    else if (muon->standAloneMuon().isNonnull()) trackRef = muon->standAloneMuon();

    
    // z planes
    int endcapPlane=0;
    if (trackRef->eta() > 0) endcapPlane=1; 
    if (trackRef->eta() < 0) endcapPlane=-1; 

    float zzPlaneME11 = endcapPlane*585;  
    float zzPlaneME1  = endcapPlane*615;  
    float zzPlaneME2  = endcapPlane*830;  
    float zzPlaneME3  = endcapPlane*935;  

    // ... to ME1/1
    // track at ME1/1 surface, +/-5.85 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, zzPlaneME11);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_me11->push_back(xx);	
      muons_y_me11->push_back(yy);	
      muons_z_me11->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_me11->push_back(acos(cosphi));
      else       muons_phi_me11->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_me11->push_back(abspseta);
      else       muons_eta_me11->push_back(-abspseta);

      if (printLevel > 2 ) {
        cout<<"\nProjecting the track to the ME";
        if (endcapPlane>0) cout << "+";
        if (endcapPlane<0) cout << "-";
        cout<< "1/1 surface:" << endl;
        cout << " * muons_x_me11 = " << xx << " cm\n";
        cout << " * muons_y_me11 = " << yy << " cm\n";
        cout << " * muons_z_me11 = " << zz << " cm\n";      

        if (yy>=0) cout << " * muons_phi_me11 = " << acos(cosphi) << endl;
        else       cout << " * muons_phi_me11 = " << 2*PI-acos(cosphi) << endl; 
        if (zz>=0) cout << " * muons_eta_me11 = " <<  abspseta << endl;
        else       cout << " * muons_eta_me11 = " << -abspseta << endl;
      }
    }
    else {
      if (printLevel > 2) cout << "\nExtrapolation to ME1/1 NOT valid\n";

      muons_x_me11->push_back(-999);	
      muons_y_me11->push_back(-999);	
      muons_z_me11->push_back(-999);	
      muons_phi_me11->push_back(-999);
      muons_eta_me11->push_back(-999);
    }


    // ... to ME1
    // track at ME1 surface, +/-6.15 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, zzPlaneME1);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_me1->push_back(xx);	
      muons_y_me1->push_back(yy);	
      muons_z_me1->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_me1->push_back(acos(cosphi));
      else       muons_phi_me1->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_me1->push_back(abspseta);
      else       muons_eta_me1->push_back(-abspseta);

      if (printLevel > 2) {
          cout<<"\nProjecting the track to the ME";
        if (endcapPlane>0) cout << "+";
        if (endcapPlane<0) cout << "-";
        cout<< "1 surface:" << endl;
        cout << " * muons_x_me1 = " << xx << " cm\n";
        cout << " * muons_y_me1 = " << yy << " cm\n";
        cout << " * muons_z_me1 = " << zz << " cm\n";      
        if (yy>=0) cout << " * muons_phi_me1 = " << acos(cosphi)      << endl;
        else       cout << " * muons_phi_me1 = " << 2*PI-acos(cosphi) << endl;
        if (zz>=0) cout << " * muons_eta_me1 = " <<  abspseta << endl;
        else       cout << " * muons_eta_me1 = " << -abspseta << endl;
      }
    }
    else {
      if (printLevel > 2) cout << "\nExtrapolation to ME1 NOT valid\n";

      muons_x_me1->push_back(-999);	
      muons_y_me1->push_back(-999);	
      muons_z_me1->push_back(-999);	
      muons_phi_me1->push_back(-999);
      muons_eta_me1->push_back(-999);
    }


    // ... to ME2
    // track at ME2 surface, +/-8.30 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, zzPlaneME2);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_me2->push_back(xx);	
      muons_y_me2->push_back(yy);	
      muons_z_me2->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_me2->push_back(acos(cosphi));
      else       muons_phi_me2->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_me2->push_back(abspseta);
      else       muons_eta_me2->push_back(-abspseta);

      if (printLevel > 2 ) {
        std::cout<<"\nProjecting the track to the ME";
        if (endcapPlane>0) cout << "+";
        if (endcapPlane<0) cout << "-";
        std::cout<< "2 surface:" << endl;
        std::cout << " * muons_x_me2 = " << xx << " cm\n";
        std::cout << " * muons_y_me2 = " << yy << " cm\n";
        std::cout << " * muons_z_me2 = " << zz << " cm\n";      
        if (yy>=0) cout << " * muons_phi_me2 = " << acos(cosphi)      << endl;
        else       cout << " * muons_phi_me2 = " << 2*PI-acos(cosphi) << endl;
        if (zz>=0) cout << " * muons_eta_me2 = " << abspseta << endl;
        else       cout << " * muons_eta_me2 = " << -abspseta << endl;
      }
    }
    else {
      if (printLevel > 2) cout << "\nExtrapolation to ME2 NOT valid\n";
  
      muons_x_me2->push_back(-999);	
      muons_y_me2->push_back(-999);	
      muons_z_me2->push_back(-999);	
      muons_phi_me2->push_back(-999);
      muons_eta_me2->push_back(-999);
    }
    
    // ... to ME3
    // track at ME3 surface, +/-9.35 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, zzPlaneME3);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_me3->push_back(xx);	
      muons_y_me3->push_back(yy);	
      muons_z_me3->push_back(zz);	
      
      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_me3->push_back(acos(cosphi));
      else       muons_phi_me3->push_back(2*PI-acos(cosphi));
      
      
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_me3->push_back(abspseta);
      else       muons_eta_me3->push_back(-abspseta);
      
      if (printLevel > 2) {
        cout<<"\nProjecting the track to the ME";
        if (endcapPlane>0) cout << "+";
        if (endcapPlane<0) cout << "-";
        cout<< "3 surface:" << endl;
        cout << " * muons_x_me3 = " << xx << " cm\n";
        cout << " * muons_y_me3 = " << yy << " cm\n";
        cout << " * muons_z_me3 = " << zz << " cm\n";      
        if (yy>=0) cout << " * muons_phi_me3 = " << acos(cosphi) << endl;
        else       cout << " * muons_phi_me3 = " << 2*PI-acos(cosphi) << endl;
        if (zz>=0) cout << " * muons_eta_me3 = " << abspseta << endl;
        else       cout << " * muons_eta_me3 = " << -abspseta << endl;
      }
    }
    else {
      if (printLevel > 2) cout << "\nExtrapolation to ME3 NOT valid\n";

      muons_x_me3->push_back(-999);	
      muons_y_me3->push_back(-999);	
      muons_z_me3->push_back(-999);	
      muons_phi_me3->push_back(-999);
      muons_eta_me3->push_back(-999);
    }


    //---------------------------------------------------------------------
    // RECHIT information in CSC: only for standalone/global muons!
    //---------------------------------------------------------------------
    if (printLevel > 2) 
      cout<< "\n\n"
	  <<"*******************************************\n"
	  <<"R E C H I T S\n"
	  <<"*******************************************\n";
    
    // An artificial rank to sort RecHits 
    // the closer RecHit to key station/layer -> the smaller the rank
    // KK's original idea
    // !!! Verify the logic is still in place !!!
    const int lutCSC[4][6] = { {26,24,22,21,23,25} , //ME1
			       { 6, 4, 2, 1, 3, 5} , //ME2 
			       {16,14,12,11,13,15} , //ME3
			       {36,34,32,31,33,35} };//ME4
    
    if( muon->isGlobalMuon()     ||
        muon->isStandAloneMuon()  ){

      // var needed to register the better ranked muon
      int   globalTypeRCH = -999;
      float  localEtaRCH = -999; 
      float  localPhiRCH = -999;
      float globalEtaRCH = -999; 
      float globalPhiRCH = -999;
      int globalStationRCH = -999;
      int globalChamberRCH = -999;
      int globalRingRCH = -999;
      int globalLayerRCH = -999;
      
      int iRechit = 0;
      for(trackingRecHit_iterator segm  = muon->outerTrack()->recHitsBegin(); 
 	  segm != muon->outerTrack()->recHitsEnd(); 
 	  segm++){
        
        // some basic protection
	if ( !((*segm)->isValid()) ) continue;
        
	// Hardware ID of the RecHit (in terms of wire/strips/chambers)
	DetId detid = (*segm)->geographicalId();
        
	// Interested in muon systems only
        if( detid.det() != DetId::Muon ) continue;

        //Look only at CSC Hits (CSC id is 2)
        if (detid.subdetId() != MuonSubdetId::CSC) continue;
          
        CSCDetId id(detid.rawId());
          
        // another sanity check
        if  (id.station() < 1) continue;
          
        // get the CSCSegment
        const CSCSegment* cscSegment = dynamic_cast<const CSCSegment*>(&**segm);
        // check the segment is not NULL
        if (cscSegment == NULL) continue;

	
        if (printLevel > 2) {
          std::cout << "\n * Segment in ME";
          if (id.endcap() > 0) std::cout << "+";
          if (id.endcap() < 0) std::cout << "-";          
          std::cout << id.station() << "/" << id.ring() << "/" << id.chamber() << std::endl;           
        }
       
        // try to get the CSC recHits that contribute to this segment.
        std::vector<CSCRecHit2D> theseRecHits = (*cscSegment).specificRecHits();
        // loop over the rechits
        for ( std::vector<CSCRecHit2D>::const_iterator iRH = theseRecHits.begin(); 
              iRH != theseRecHits.end(); iRH++) {
          
          // get the rechit ID
          CSCDetId idRH = (CSCDetId)(*iRH).cscDetId();
   
          // CSC chamber
          const CSCChamber* cscchamber = cscGeom->chamber(idRH);
          if (!cscchamber) continue;

          // local position
          LocalPoint rhitlocal = iRH->localPosition();
          
          // translate into global position
          GlobalPoint gp = GlobalPoint(0.0, 0.0, 0.0);
          gp = cscchamber->toGlobal(rhitlocal);
        
          // identify the rechit position
          int pos = ( ((idRH.station()-1)*6) + (idRH.layer()-1) ) + (MAX_CSC_RECHIT*whichMuon);
        
          rchEtaList[pos] = gp.eta();
          rchPhiList[pos] = gp.phi();
          float phi02PI   = gp.phi();
          if (gp.phi() < 0) phi02PI += (2*PI); 
        
          rchPhiList_02PI[pos] = phi02PI;
        
          if (printLevel > 2) {
            cout << "iRechit: " << iRechit+1 
                 << " -> pos: " << pos;
            cout << "   - RCH Type: "         << lutCSC[idRH.station()-1][idRH.layer()-1]
                 <<      " RCH Eta: "         << gp.eta()	    
                 <<      " RCH Phi: "         << gp.phi()
                 <<      " RCH Phi [0,2PI]: " << phi02PI 
                 <<      " RCH X: "           << gp.x()  
                 <<      " RCH Y: "           << gp.y()  
                 <<      " RCH Z: "           << gp.z()  
                 << endl;
          }
        
          // -------------------------------------------------- 
          // new Matrix block
          if (whichMuon < MAX_MUONS && iRechit < MAX_CSC_RECHIT) {
            rchEtaMatrixLocal[whichMuon][iRechit] = rhitlocal.eta(); 
            rchPhiMatrixLocal[whichMuon][iRechit] = rhitlocal.phi();
            rchEtaMatrix     [whichMuon][iRechit] = gp.eta(); 
            rchPhiMatrix     [whichMuon][iRechit] = gp.phi();
            rchPhi02PIMatrix [whichMuon][iRechit] = phi02PI;
            rchStationMatrix [whichMuon][iRechit] = idRH.station();
            rchChamberMatrix [whichMuon][iRechit] = idRH.chamber();
            rchRingMatrix    [whichMuon][iRechit] = idRH.ring();
            rchLayerMatrix   [whichMuon][iRechit] = idRH.layer();
            rchTypeMatrix    [whichMuon][iRechit] = lutCSC[idRH.station()-1][idRH.layer()-1];
          }
          else
            cout << "ERROR: too many segment or muons "
                 << "MAX_MUONS is currently set to "      << MAX_MUONS      << endl
                 << "whichMuon loop is currently at "     << whichMuon      << endl
                 << "MAX_CSC_RECHIT is currently set to " << MAX_CSC_RECHIT << endl
                 << "iRechit loop is currently at "       << iRechit        << endl;
        
          iRechit++;
          // --------------------------------------------------
        
        
          // --------------------------------------------------
          // See if this hit is closer to the "key" position
          if( lutCSC[idRH.station()-1][idRH.layer()-1]<globalTypeRCH || globalTypeRCH<0 ){
            globalTypeRCH    = lutCSC[idRH.station()-1][idRH.layer()-1];
            localEtaRCH      = rhitlocal.eta();
            localPhiRCH      = rhitlocal.phi();
            globalEtaRCH     = gp.eta();
            globalPhiRCH     = gp.phi();
            globalStationRCH = idRH.station();
            globalChamberRCH = idRH.chamber();
            globalRingRCH    = idRH.ring();
            globalLayerRCH   = idRH.layer();
          }
          // -------------------------------------------------- 
        
        }//end loop rechit
      }

      rchMuonSize->push_back(iRechit);

      // at the end of the loop, write only the best rechit
      rchCSCtype  -> push_back(globalTypeRCH);       
      rchEtaLocal -> push_back(localEtaRCH);     
      rchPhiLocal -> push_back(localPhiRCH);     
      rchEta      -> push_back(globalEtaRCH);     
      rchPhi      -> push_back(globalPhiRCH);     
      rchStation  -> push_back(globalStationRCH);
      rchChamber  -> push_back(globalChamberRCH);
      rchRing     -> push_back(globalRingRCH);
      rchLayer    -> push_back(globalLayerRCH);
    
      if (printLevel > 2) {
        cout << "\n######### CLOSER #########";
        cout << "\n - RCH Type: " << globalTypeRCH
             <<     " (RCH Eta: " << globalEtaRCH ;	    
      }
      
      if ( globalPhiRCH < 0) {
        rchPhi_02PI -> push_back(globalPhiRCH + (2*PI));
 	if (printLevel > 2) 
 	  cout << ", RCH Phi [0,2PI]: " << globalPhiRCH + (2*PI) << ")" << endl;
      }
      else {
        rchPhi_02PI -> push_back(globalPhiRCH         );
 	if (printLevel > 2) 
 	  cout << ", RCH Phi [0,2PI]: " << globalPhiRCH << ")" << endl;
      }
      
      if (printLevel > 2) 
	cout << "##########################\n\n";
      
    }//isGlobal || isStandAlone
    else {
      // if only tracker muons 
      rchCSCtype  -> push_back(-999);       
      rchEta      -> push_back(-999);     
      rchPhi      -> push_back(-999);     
      rchPhi_02PI -> push_back(-999);
      rchStation  -> push_back(-999);
      rchChamber  -> push_back(-999);
      rchRing     -> push_back(-999);
      rchLayer    -> push_back(-999);
    }
    
    
    //---------------------------------------------------------------------
    // tracker muons: collect the segment associated to it
    //---------------------------------------------------------------------
    if( //(muon->isGlobalMuon()     == 0) ||
        (muon->isTrackerMuon()    == 1)  ){
      
      
      if (printLevel > 2) {
        
        std::cout<< "\n\n"
                 <<"*******************************************\n"
                 <<"T R A C K E R   M U O N   B L O C K\n"
                 <<"*******************************************\n";

        cout << "Number of Chambers:" << muon->numberOfChambers() << endl;
        cout << "Number of Matches:"  << muon->numberOfMatches()  << endl;
        cout << "StationMask: "       << muon->stationMask()      << endl;
        cout << "Is Matches Valid?:"  << muon->isMatchesValid()   << endl;
      }
      
      trkNchambers    -> push_back(muon->numberOfChambers());
      trkNofMatches   -> push_back(muon->numberOfMatches() );
      trkIsMatchValid -> push_back(muon->isMatchesValid()  );

      if (printLevel > 2)
        cout << "Is Muon All Arbitrated? " 
             << muon::isGoodMuon(*muon, muon::AllArbitrated) << endl;
      
      // fill the vectors only if the muon is arbitrated
      if (muon::isGoodMuon(*muon, muon::AllArbitrated)) {

        int iChamber=0;
        int iSegment=0;
        for( std::vector<MuonChamberMatch>::const_iterator chamber = muon->matches().begin();
             chamber != muon->matches().end(); ++chamber) {
          
          if (printLevel > 2) cout << "Chamber:" << iChamber+1 << endl;
          
          // we are interested only in CSC chamber
          if( chamber->detector() != MuonSubdetId::CSC ) continue; 

          //get the handle to the info, CSCDetId            
          CSCDetId cscId(chamber->id.rawId());
            
          // extract the information
          int chamberId     = cscId.chamber();
          int ring          = cscId.ring();
          int station       = cscId.station();
          int endcap        = cscId.endcap();
          int triggerSector = cscId.triggerSector();
          int triggerCscId  = cscId.triggerCscId();
          

          GlobalPoint gpChb = theService->trackingGeometry()->idToDet(cscId)->surface().toGlobal(LocalPoint(chamber->x,chamber->y,0));
          
          float xx = gpChb.x();
          float yy = gpChb.y();
          float zz = gpChb.z();
        
          float rr = sqrt(xx*xx + yy*yy);
          float phi = gpChb.phi();
          if (phi < 0   ) phi += 2*PI;
          if (phi > 2*PI) phi -= 2*PI;
          float eta = gpChb.eta();


          //look at the segments
          for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->segmentMatches.begin();
                segment != chamber->segmentMatches.end(); ++segment ){
            
            if (printLevel > 2) cout << "Segment:" << iSegment+1 << endl;
   
            int segIsArb=(segment->Arbitrated)>>8;

            if (!segIsArb) continue;
            
            // global coordinates
            GlobalPoint gpSeg = theService->trackingGeometry()->idToDet(cscId)->surface().toGlobal(LocalPoint(segment->x,segment->y,0));

            float segX=gpSeg.x();
            float segY=gpSeg.y();
            float segZ=gpSeg.z();
            
            float segR = sqrt(segX*segX + segY*segY);
            float segPhi = gpSeg.phi();
            if (segPhi < 0   ) segPhi += 2*PI;
            if (segPhi > 2*PI) segPhi -= 2*PI;
            float segEta = gpSeg.eta();

            float segDxDz = segment->dXdZ;
            float segDyDz = segment->dYdZ;
            float segDxDzErr = segment->dXdZErr;
            float segDyDzErr = segment->dYdZErr;


            if (printLevel > 2) {
              cout << "###### CSC SEGMENTS ########"          << endl;
              cout << "trkSegChamberId:"     << chamberId     << endl;
              cout << "trkSegRing:"          << ring          << endl;
              cout << "trkSegStation:"       << station       << endl;
              cout << "trkSegEndcap:"        << endcap        << endl;
              cout << "trkSegTriggerSector:" << triggerSector << endl;
              cout << "trkSegTriggerCscId :" << triggerCscId  << endl;
            
              cout<< "trkSegXfromMatch:"     << xx   << endl;
              cout<< "trkSegYfromMatch:"     << yy   << endl;
              cout<< "trkSegZfromMatch:"     << zz   << endl;
              cout<< "trkSegRfromMatch:"     << rr   << endl;
              
              cout << "trkSegPhifromMatch:" << phi << endl;
              cout << "trkSegEtafromMatch:" << eta << endl;

              cout << "SEGMENT:"                  << endl;
              cout << "trkSegIsArb: " << segIsArb << endl; 
              cout << "trkSegX: "     << segX     << endl; 
              cout << "trkSegY: "     << segY     << endl; 
              cout << "trkSegZ: "     << segZ     << endl; 
              cout << "trkSegR: "     << segR     << endl; 
              
              cout << "segPhi:"       << segPhi       << endl;
              cout << "segEta:"       << segEta       << endl;

              cout << "segDxDz:"      << segDxDz      << endl;
              cout << "segDyDz:"      << segDyDz      << endl;
              cout << "segDxDzErr:"   << segDxDzErr   << endl;
              cout << "segDyDzErr:"   << segDyDzErr   << endl;
            }
            
            
            // fill the matrices
            if (whichMuon < MAX_MUONS && iSegment<MAX_TRK_SEGS) {
              trkSegChamberId    [whichMuon][iSegment] = chamberId;
              trkSegRing         [whichMuon][iSegment] = ring;
              trkSegStation      [whichMuon][iSegment] = station;
              trkSegEndcap       [whichMuon][iSegment] = endcap;
              trkSegTriggerSector[whichMuon][iSegment] = triggerSector;
              trkSegTriggerCscId [whichMuon][iSegment] = triggerCscId;
            
              trkSegXfromMatch[whichMuon][iSegment] = xx;
              trkSegYfromMatch[whichMuon][iSegment] = yy;
              trkSegZfromMatch[whichMuon][iSegment] = zz;
              trkSegRfromMatch[whichMuon][iSegment] = rr;
              
              trkSegPhifromMatch[whichMuon][iSegment] = phi;
              trkSegEtafromMatch[whichMuon][iSegment] = eta;

              
              trkSegIsArb[whichMuon][iSegment]=segIsArb;
              trkSegX[whichMuon][iSegment]=segX;
              trkSegY[whichMuon][iSegment]=segY;
              trkSegZ[whichMuon][iSegment]=segZ;
              trkSegR[whichMuon][iSegment]=segR;
            
              trkSegPhi[whichMuon][iSegment]= segPhi;
              trkSegEta[whichMuon][iSegment]= segEta;

              trkSegDxDz[whichMuon][iSegment]   = segDxDz;
              trkSegDyDz[whichMuon][iSegment]   = segDyDz;
              trkSegDxDzErr[whichMuon][iSegment]= segDxDzErr;
              trkSegDyDzErr[whichMuon][iSegment]= segDyDzErr;


            }         
            else
              cout << "ERROR: too many segment or muons "
                   << "MAX_MUONS is currently set to "    << MAX_MUONS    << endl
                   << "whichMuon loop is currently at "   << whichMuon    << endl
                   << "MAX_TRK_SEGS is currently set to " << MAX_TRK_SEGS << endl
                   << "iSegment loop is currently at "    << iSegment     << endl;

            iSegment++;        
          }//end of segment loop 
          iChamber++;        
        }//end loop on the chambers
        if (printLevel > 2) cout << "trkNSegs="<<iSegment<<endl;
        trkNSegs->push_back(iSegment);
      }// is muon good?
      else trkNSegs->push_back(-999);
    }//isTrackerMuon
    else trkNSegs->push_back(-999);
  }

}


void TrigEff::muonsInit()
{
  
  isGlobalMuon        = new vector<int>;	  
  isTrackerMuon	      = new vector<int>;
  isStandAloneMuon    = new vector<int>; 
  isMuonAllArbitrated = new vector<int>; 
  isTMLastStationAngTight = new vector<int>; 
  isGlobalMuonPromptTight = new vector<int>; 

  isEnergyValid  = new vector<int>; 
  caloCompatibility  = new vector<float>; 
  em     = new vector<float>;
  emS9   = new vector<float>;
  emS25  = new vector<float>;
  emMax  = new vector<float>;
  had    = new vector<float>;
  hadS9  = new vector<float>;
  hadMax = new vector<float>;
  ho     = new vector<float>;
  hoS9   = new vector<float>;

  //gmrEnergy            = new vector<float>;
  //gmrDEnergy           = new vector<float>;
  gmrPt	               = new vector<float>;
  gmrEta	       = new vector<float>;
  gmrPhi	       = new vector<float>;
  gmrP	               = new vector<float>;
  gmrPx	               = new vector<float>;
  gmrPy	               = new vector<float>;
  gmrPz	               = new vector<float>;
  gmrTheta	       = new vector<float>;
  gmrVx	               = new vector<float>;
  gmrVy	               = new vector<float>;
  gmrVz	               = new vector<float>;
  gmrCharge            = new vector<float>;
  gmrNDoF	       = new vector<float>;
  gmrChi2	       = new vector<float>;
  gmrChi2Norm          = new vector<float>;
  gmrDXY	       = new vector<float>;
  gmrDTheta            = new vector<float>;
  gmrDPt	       = new vector<float>;
  gmrDEta	       = new vector<float>;
  gmrDPhi	       = new vector<float>;
  gmrDDXY	       = new vector<float>;
  gmrIso03nTracks      = new vector<float>;
  gmrIso03sumPt        = new vector<float>;
  gmrDz	               = new vector<float>;
  gmrD0	               = new vector<float>;
  gmrDsz	       = new vector<float>;
  gmrDDz	       = new vector<float>;
  gmrDD0	       = new vector<float>;
  gmrDDsz	       = new vector<float>;
  gmrInnerX            = new vector<float>;
  gmrInnerY            = new vector<float>;
  gmrInnerZ            = new vector<float>;
  gmrOuterX            = new vector<float>;
  gmrOuterY            = new vector<float>;
  gmrOuterZ            = new vector<float>;
  gmrValHits           = new vector<int>;

  //stdEnergy    = new vector<float>;
  //stdDEnergy   = new vector<float>;
  stdPt        = new vector<float>;
  stdEta       = new vector<float>;
  stdPhi       = new vector<float>;
  stdPx        = new vector<float>;
  stdPy        = new vector<float>;
  stdPz        = new vector<float>;
  stdVx        = new vector<float>;
  stdVy        = new vector<float>;
  stdVz        = new vector<float>;
  stdCharge    = new vector<float>;
  stdDPt       = new vector<float>;
  stdDEta      = new vector<float>;
  stdDPhi      = new vector<float>;
  stdDz        = new vector<float>;
  stdD0        = new vector<float>;
  stdNDoF      = new vector<float>;  
  stdChi2      = new vector<float>;
  stdChi2Norm  = new vector<float>;
  stdDXY       = new vector<float>;
  stdTheta     = new vector<float>;
  stdDTheta    = new vector<float>;
  stdDDz       = new vector<float>;
  stdDD0       = new vector<float>;     
  stdValHits   = new vector<int>;
  stdInnerX    = new vector<float>;
  stdInnerY    = new vector<float>;
  stdInnerZ    = new vector<float>;
  stdOuterX    = new vector<float>;
  stdOuterY    = new vector<float>;
  stdOuterZ    = new vector<float>;
   
  //trkEnergy            = new vector<float>;
  //trkDEnergy           = new vector<float>;
  trkPt                = new vector<float>;
  trkEta               = new vector<float>;
  trkPhi               = new vector<float>;
  trkPx                = new vector<float>;
  trkPy                = new vector<float>;
  trkPz                = new vector<float>;
  trkVx                = new vector<float>;
  trkVy                = new vector<float>;
  trkVz                = new vector<float>;
  trkCharge            = new vector<float>;
  trkDPt               = new vector<float>;
  trkDEta              = new vector<float>;
  trkDPhi              = new vector<float>;
  trkDz                = new vector<float>;
  trkD0                = new vector<float>;
  trkNDoF              = new vector<float>; 
  trkChi2              = new vector<float>;
  trkChi2Norm          = new vector<float>;
  trkDXY               = new vector<float>;
  trkTheta             = new vector<float>;
  trkDTheta            = new vector<float>;
  trkDDz               = new vector<float>;
  trkDD0               = new vector<float>;     
  trkValHits           = new vector<int>;

  // ------------------------------------------------------
  //segment
  trkNchambers        = new vector<int>;
  trkNofMatches       = new vector<int>;
  trkIsMatchValid     = new vector<int>;

  trkNSegs            = new vector<int>;
  for (int row=0; row < MAX_MUONS; row++) 
    for (int col=0; col < MAX_TRK_SEGS; col++) {
      trkSegChamberId[row][col] = -999;   
      trkSegRing[row][col] = -999;        
      trkSegStation[row][col] = -999;     
      trkSegEndcap[row][col] = -999;      
      trkSegTriggerSector[row][col] = -999;
      trkSegTriggerCscId[row][col] = -999;
         
      trkSegXfromMatch[row][col] = -999;  
      trkSegYfromMatch[row][col] = -999;  
      trkSegZfromMatch[row][col] = -999;  
      trkSegRfromMatch[row][col] = -999;  
      trkSegPhifromMatch[row][col] = -999;
      trkSegEtafromMatch[row][col] = -999;

      trkSegIsArb[row][col] = -999;
      trkSegX[row][col]     = -999;
      trkSegY[row][col]     = -999;
      trkSegZ[row][col]     = -999;
      trkSegR[row][col]     = -999;
      trkSegPhi[row][col]   = -999;
      trkSegEta[row][col]   = -999;

      trkSegDxDz[row][col]     = -999;
      trkSegDyDz[row][col]     = -999;
      trkSegDxDzErr[row][col]  = -999;
      trkSegDyDzErr[row][col]  = -999;
    }
  // ------------------------------------------------------  


  rchCSCtype  = new vector<int>  ; 
  rchEtaLocal = new vector<float>;     
  rchPhiLocal = new vector<float>;     
  rchEta      = new vector<float>;     
  rchPhi      = new vector<float>;     
  rchPhi_02PI = new vector<float>;
  rchStation  = new vector<int>;
  rchChamber  = new vector<int>;
  rchRing     = new vector<int>;
  rchLayer    = new vector<int>;

  rchMuonSize = new vector<int>;
  for (int row=0; row < MAX_MUONS; row++) 
    for (int col=0; col < MAX_CSC_RECHIT; col++) {
      rchEtaMatrixLocal[row][col] = -999;
      rchPhiMatrixLocal[row][col] = -999;
      rchEtaMatrix[row][col]     = -999;
      rchPhiMatrix[row][col]     = -999;
      rchPhi02PIMatrix[row][col] = -999;
      rchStationMatrix[row][col] = -999;
      rchChamberMatrix[row][col] = -999;
      rchRingMatrix[row][col]    = -999;
      rchLayerMatrix[row][col]   = -999;
      rchTypeMatrix[row][col]    = -999;
    }
  

  // propagation to ME1/1
    muons_x_me11 = new vector<float>;
    muons_y_me11 = new vector<float>;
    muons_z_me11 = new vector<float>;
  muons_phi_me11 = new vector<float>;
  muons_eta_me11 = new vector<float>;

  // propagation to ME1
    muons_x_me1 = new vector<float>;
    muons_y_me1 = new vector<float>;
    muons_z_me1 = new vector<float>;
  muons_phi_me1 = new vector<float>;
  muons_eta_me1 = new vector<float>;
                 
  // propagation to ME2
    muons_x_me2 = new vector<float>;
    muons_y_me2 = new vector<float>;
    muons_z_me2 = new vector<float>;
  muons_phi_me2 = new vector<float>;
  muons_eta_me2 = new vector<float>;

  // propagation to ME3                 
    muons_x_me3 = new vector<float>;
    muons_y_me3 = new vector<float>;
    muons_z_me3 = new vector<float>;
  muons_phi_me3 = new vector<float>;
  muons_eta_me3 = new vector<float>;
                 
  // csc segments
  cscsegs_loc_x = new vector<float>;
  cscsegs_loc_y = new vector<float>;
  cscsegs_loc_z = new vector<float>;

  cscsegs_loc_theta = new vector<float>;
  cscsegs_loc_eta   = new vector<float>;
  cscsegs_loc_phi   = new vector<float>;

  cscsegs_loc_dir_theta = new vector<float>;
  cscsegs_loc_dir_eta   = new vector<float>;
  cscsegs_loc_dir_phi   = new vector<float>;

  cscsegs_gbl_x = new vector<float>;
  cscsegs_gbl_y = new vector<float>;
  cscsegs_gbl_z = new vector<float>;

  cscsegs_gbl_theta = new vector<float>;
  cscsegs_gbl_eta   = new vector<float>;
  cscsegs_gbl_phi   = new vector<float>;

  cscsegs_gbl_dir_theta = new vector<float>;
  cscsegs_gbl_dir_eta   = new vector<float>;
  cscsegs_gbl_dir_phi   = new vector<float>;

  isdtseg         = new vector<int>;
  cscsegs_sector  = new vector<int>;
  cscsegs_endcap  = new vector<int>;
  cscsegs_station = new vector<int>;
  cscsegs_ring    = new vector<int>;
  cscsegs_chamber = new vector<int>;
  cscsegs_wheel   = new vector<int>;

  /*dtsegs_gbl_phi  = new vector<float>;
  dtsegs_gbl_eta  = new vector<float>;
  dtsegs_chamber  = new vector<int>;
  dtsegs_wheel    = new vector<int>;
  dtsegs_station  = new vector<int>;
  dtsegs_sector   = new vector<int>;
  */

  muonNsegs = new std::vector<int>; 

  for (int row=0; row < MAX_MUONS; row++) 
    for (int col=0; col < MAX_SEGS_STD; col++) {
 
      muon_cscsegs_loc_x[row][col] = -999;
      muon_cscsegs_loc_y[row][col] = -999;
      muon_cscsegs_loc_eta[row][col] = -999;
      muon_cscsegs_loc_phi[row][col] = -999;
      muon_cscsegs_loc_dir_eta[row][col] = -999;
      muon_cscsegs_loc_dir_phi[row][col] = -999;

      muon_cscsegs_gbl_x[row][col] = -999;
      muon_cscsegs_gbl_y[row][col] = -999;
      muon_cscsegs_gbl_eta[row][col] = -999;
      muon_cscsegs_gbl_phi[row][col] = -999;
      muon_cscsegs_gbl_dir_eta[row][col] = -999;
      muon_cscsegs_gbl_dir_phi[row][col] = -999;

      muon_cscsegs_dxdz[row][col] = -999;
      muon_cscsegs_dydz[row][col] = -999;
      muon_cscsegs_dxdzErr[row][col] = -999;
      muon_cscsegs_dydzErr[row][col] = -999;

      muon_cscsegs_endcap[row][col] = -999;
      muon_cscsegs_station[row][col] = -999;
      muon_cscsegs_ring[row][col] = -999;
      muon_cscsegs_chamber[row][col] = -999;
      muon_cscsegs_nhits[row][col] = -999;

      muon_cscsegs_islctable[row][col] = -999;
      muon_cscsegs_ismatched[row][col] = -999;
      muon_cscsegs_lctId[row][col] = -999;
      muon_cscsegs_nmatched[row][col] = -999;
      muon_isdtseg[row][col] = -999;

      muon_dtsegs_lctId[row][col] = -999;    
      muon_dtsegs_ismatched[row][col] = -999;
      muon_dtsegs_chamber[row][col] = -999;
      muon_cscsegs_wheel[row][col] = -999;
      muon_cscsegs_sector[row][col] = -999;
      muon_dtsegs_station[row][col] = -999;
      muon_dtsegs_gbl_eta[row][col] = -999;
      muon_dtsegs_gbl_phi[row][col] = -999;

    }
}


void TrigEff::muonsDel() {
  
  vector<int>().swap(*isGlobalMuon);	  
  vector<int>().swap(*isTrackerMuon);
  vector<int>().swap(*isStandAloneMuon); 
  vector<int>().swap(*isMuonAllArbitrated); 
  vector<int>().swap(*isTMLastStationAngTight); 
  vector<int>().swap(*isGlobalMuonPromptTight); 

  vector<int>().swap(*isEnergyValid); 
  vector<float>().swap(*caloCompatibility); 
  vector<float>().swap(*em); 
  vector<float>().swap(*emS9);
  vector<float>().swap(*emS25);
  vector<float>().swap(*emMax);
  vector<float>().swap(*had);
  vector<float>().swap(*hadS9);
  vector<float>().swap(*hadMax);
  vector<float>().swap(*ho);
  vector<float>().swap(*hoS9);

  //vector<float>().swap(*gmrEnergy);
  //vector<float>().swap(*gmrDEnergy);
  vector<float>().swap(*gmrPt);
  vector<float>().swap(*gmrEta);
  vector<float>().swap(*gmrPhi);
  vector<float>().swap(*gmrP);
  vector<float>().swap(*gmrPx);
  vector<float>().swap(*gmrPy);
  vector<float>().swap(*gmrPz);
  vector<float>().swap(*gmrTheta);
  vector<float>().swap(*gmrVx);
  vector<float>().swap(*gmrVy);
  vector<float>().swap(*gmrVz);
  vector<float>().swap(*gmrCharge);
  vector<float>().swap(*gmrNDoF);
  vector<float>().swap(*gmrChi2);
  vector<float>().swap(*gmrChi2Norm);
  vector<float>().swap(*gmrDXY);
  vector<float>().swap(*gmrDTheta);
  vector<float>().swap(*gmrDPt);
  vector<float>().swap(*gmrDEta);
  vector<float>().swap(*gmrDPhi);
  vector<float>().swap(*gmrDDXY);
  vector<float>().swap(*gmrIso03nTracks);
  vector<float>().swap(*gmrIso03sumPt);
  vector<float>().swap(*gmrChi2Norm);
  vector<float>().swap(*gmrDz);
  vector<float>().swap(*gmrD0);
  vector<float>().swap(*gmrDsz);
  vector<float>().swap(*gmrDDz);
  vector<float>().swap(*gmrDD0);
  vector<float>().swap(*gmrDDsz);
  vector<float>().swap(*gmrInnerX);
  vector<float>().swap(*gmrInnerY);
  vector<float>().swap(*gmrInnerZ);
  vector<float>().swap(*gmrOuterX);
  vector<float>().swap(*gmrOuterY);
  vector<float>().swap(*gmrOuterZ);
  vector<int>().swap(*gmrValHits);


  //vector<float>().swap(*stdEnergy);
  //vector<float>().swap(*stdDEnergy);
  vector<float>().swap(*stdPt);
  vector<float>().swap(*stdEta);
  vector<float>().swap(*stdPhi);
  vector<float>().swap(*stdPx);
  vector<float>().swap(*stdPy);
  vector<float>().swap(*stdPz);
  vector<float>().swap(*stdVx);
  vector<float>().swap(*stdVy);
  vector<float>().swap(*stdVz);
  vector<float>().swap(*stdCharge);
  vector<float>().swap(*stdDPt);
  vector<float>().swap(*stdDEta);
  vector<float>().swap(*stdDPhi);
  vector<float>().swap(*stdDz);
  vector<float>().swap(*stdD0);
  vector<float>().swap(*stdNDoF); 
  vector<float>().swap(*stdChi2);
  vector<float>().swap(*stdChi2Norm);
  vector<float>().swap(*stdDXY);
  vector<float>().swap(*stdTheta);
  vector<float>().swap(*stdDTheta);
  vector<float>().swap(*stdDDz);
  vector<float>().swap(*stdDD0);     
  vector<int>().swap(*stdValHits);
  vector<float>().swap(*stdInnerX);
  vector<float>().swap(*stdInnerY);
  vector<float>().swap(*stdInnerZ);
  vector<float>().swap(*stdOuterX);
  vector<float>().swap(*stdOuterY);
  vector<float>().swap(*stdOuterZ);
   

  //vector<float>().swap(*trkEnergy);
  //vector<float>().swap(*trkDEnergy);
  vector<float>().swap(*trkPt);
  vector<float>().swap(*trkEta);
  vector<float>().swap(*trkPhi);
  vector<float>().swap(*trkPx);
  vector<float>().swap(*trkPy);
  vector<float>().swap(*trkPz);
  vector<float>().swap(*trkVx);
  vector<float>().swap(*trkVy);
  vector<float>().swap(*trkVz);
  vector<float>().swap(*trkCharge);
  vector<float>().swap(*trkDPt);
  vector<float>().swap(*trkDEta);
  vector<float>().swap(*trkDPhi);
  vector<float>().swap(*trkDz);
  vector<float>().swap(*trkD0);
  vector<float>().swap(*trkNDoF); 
  vector<float>().swap(*trkChi2);
  vector<float>().swap(*trkChi2Norm );
  vector<float>().swap(*trkDXY);
  vector<float>().swap(*trkTheta);
  vector<float>().swap(*trkDTheta);
  vector<float>().swap(*trkDDz);
  vector<float>().swap(*trkDD0);     
  vector<int>().swap(*trkValHits);

  // CSC segment for the tracker muon
  //chamber
  vector<int>().swap(*trkNchambers); 
  vector<int>().swap(*trkNofMatches);
  vector<int>().swap(*trkIsMatchValid);

  vector<int>().swap(*trkNSegs);

  vector<int>  ().swap(*rchCSCtype ); 
  vector<float>().swap(*rchEtaLocal);     
  vector<float>().swap(*rchPhiLocal);     
  vector<float>().swap(*rchEta     );     
  vector<float>().swap(*rchPhi     );     
  vector<float>().swap(*rchPhi_02PI);

  vector<int>().swap(*rchStation);
  vector<int>().swap(*rchChamber);
  vector<int>().swap(*rchRing);
  vector<int>().swap(*rchLayer);
  
  vector<int>().swap(*rchMuonSize);

  delete [] rchEtaList;
  delete [] rchPhiList;
  delete [] rchPhiList_02PI;

    vector<float>().swap(*muons_x_me11);
    vector<float>().swap(*muons_y_me11);
    vector<float>().swap(*muons_z_me11);
  vector<float>().swap(*muons_phi_me11);
  vector<float>().swap(*muons_eta_me11);
                                    
    vector<float>().swap(*muons_x_me1);
    vector<float>().swap(*muons_y_me1);
    vector<float>().swap(*muons_z_me1);
  vector<float>().swap(*muons_phi_me1);
  vector<float>().swap(*muons_eta_me1);
                                    
    vector<float>().swap(*muons_x_me2);
    vector<float>().swap(*muons_y_me2);
    vector<float>().swap(*muons_z_me2);
  vector<float>().swap(*muons_phi_me2);
  vector<float>().swap(*muons_eta_me2);
                       
    vector<float>().swap(*muons_x_me3);
    vector<float>().swap(*muons_y_me3);
    vector<float>().swap(*muons_z_me3);
  vector<float>().swap(*muons_phi_me3);
  vector<float>().swap(*muons_eta_me3);

  vector<float>().swap(*cscsegs_loc_x);
  vector<float>().swap(*cscsegs_loc_y);
  vector<float>().swap(*cscsegs_loc_z);

  vector<float>().swap(*cscsegs_loc_theta);
  vector<float>().swap(*cscsegs_loc_eta);
  vector<float>().swap(*cscsegs_loc_phi);

  vector<float>().swap(*cscsegs_loc_dir_theta);
  vector<float>().swap(*cscsegs_loc_dir_eta);
  vector<float>().swap(*cscsegs_loc_dir_phi);

  vector<float>().swap(*cscsegs_gbl_x);
  vector<float>().swap(*cscsegs_gbl_y);
  vector<float>().swap(*cscsegs_gbl_z);

  vector<float>().swap(*cscsegs_gbl_theta);
  vector<float>().swap(*cscsegs_gbl_eta);
  vector<float>().swap(*cscsegs_gbl_phi);

  vector<float>().swap(*cscsegs_gbl_dir_theta);
  vector<float>().swap(*cscsegs_gbl_dir_eta);
  vector<float>().swap(*cscsegs_gbl_dir_phi);

  vector<int>().swap(*cscsegs_wheel);
  vector<int>().swap(*cscsegs_sector);
  vector<int>().swap(*cscsegs_endcap);
  vector<int>().swap(*cscsegs_station);
  vector<int>().swap(*cscsegs_ring);
  vector<int>().swap(*cscsegs_chamber);

  //vector<float>().swap(*dtsegs_gbl_phi);
  //vector<float>().swap(*dtsegs_gbl_eta);
  vector<int>().swap(*isdtseg);
  //vector<int>().swap(*dtsegs_chamber);
  vector<int>().swap(*cscsegs_wheel);
  //vector<int>().swap(*dtsegs_sector);
  //vector<int>().swap(*dtsegs_station);
  //vector<int>().swap(*muonNsegs); 
}


//-------------------------------------------------------------------------
// fill the l1 extra collection
//-------------------------------------------------------------------------
void TrigEff::fillExtra(l1extra::L1MuonParticleCollection::const_iterator l1muon) {
  
  l1Eta->push_back(l1muon->eta());
  l1Pt ->push_back(l1muon->pt());
  l1Phi->push_back(l1muon->phi());
          
  isIsolated->push_back(l1muon->isIsolated());
  isMip     ->push_back(l1muon->isMip()     );
  isForward ->push_back(l1muon->isForward() );
  isRPC     ->push_back(l1muon->isRPC()     );
  detectorType ->push_back(l1muon->gmtMuonCand().detector());
  rank ->push_back(l1muon->gmtMuonCand().rank());
}
  

void TrigEff::l1extraInit() {
  
  l1Eta = new vector<float>; 
  l1Pt  = new vector<float>; 
  l1Phi = new vector<float>; 
                               
  isIsolated = new vector<int>;   
  isMip      = new vector<int>;   
  isForward  = new vector<int>;   
  isRPC      = new vector<int>;   

  detectorType = new vector<int>;  
  rank         = new vector<int>;
}

void TrigEff::l1extraDel() {
  
  vector<float>().swap(*l1Eta); 
  vector<float>().swap(*l1Pt ); 
  vector<float>().swap(*l1Phi); 
                               
  vector<int>().swap(*isIsolated);   
  vector<int>().swap(*isMip     );   
  vector<int>().swap(*isForward );   
  vector<int>().swap(*isRPC     );   

  vector<int>().swap(*detectorType);  
  vector<int>().swap(*rank);
}


void TrigEff::resizeRchHits(int nMu_nMaxCscRchHits){
  
  if (!rchEtaList)      delete [] rchEtaList;
  if (!rchPhiList)      delete [] rchPhiList;
  if (!rchPhiList_02PI) delete [] rchPhiList_02PI;

  rchEtaList      = new Double_t[nMu_nMaxCscRchHits];
  rchPhiList      = new Double_t[nMu_nMaxCscRchHits];
  rchPhiList_02PI = new Double_t[nMu_nMaxCscRchHits];
  
  recoMuons->SetBranchAddress("rchEtaList"     , rchEtaList); 
  recoMuons->SetBranchAddress("rchPhiList"     , rchPhiList);
  recoMuons->SetBranchAddress("rchPhiList_02PI", rchPhiList_02PI);
  
  for (int i = 0; i < nMu_nMaxCscRchHits; i++) {
    rchEtaList[i]      = -999;
    rchPhiList[i]      = -999;
    rchPhiList_02PI[i] = -999;
  }
}

void TrigEff::csctfInit() {	
  
  EndcapTrk     = new vector<int>;
  SectorTrk     = new vector<int>;  
  BxTrk         = new vector<int>;  
     	    
  me1ID         = new vector<int>;  
  me2ID         = new vector<int>;
  me3ID         = new vector<int>;  
  me4ID         = new vector<int>;  
  mb1ID         = new vector<int>;  

  OutputLinkTrk = new vector<int>;      

  ModeTrk       = new vector<int>  ;
  EtaTrk        = new vector<float>;  
  PhiTrk        = new vector<float>;  
  PhiTrk_02PI   = new vector<float>;  
  PtTrk         = new vector<float>;  
  
  ChargeValidTrk = new vector<int>;
  ChargeTrk      = new vector<int>;
  QualityTrk     = new vector<int>;
  ForRTrk        = new vector<int>;
  Phi23Trk       = new vector<int>;
  Phi12Trk       = new vector<int>;
  PhiSignTrk     = new vector<int>;

  EtaBitTrk      = new vector<int>;
  PhiBitTrk      = new vector<int>;
  PtBitTrk       = new vector<int>;

  NumLCTsTrk    = new vector<int>; 
    
  /*
    trLctEndcap.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctSector.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctSubSector.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctBx.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctBx0.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
  
    trLctStation.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctRing.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctChamber.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctTriggerCSCID.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctFpga.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK);	  
    
    trLctlocalPhi.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctglobalPhi.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK);   
    trLctglobalEta.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
  
    trLctstripNum.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK);   
    trLctwireGroup.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK);     
  */

  // all LCTs
  lctEndcap       = new vector<int>; 
  lctSector       = new vector<int>; 
  lctSubSector    = new vector<int>; 
  lctBx           = new vector<int>; 
  lctBx0          = new vector<int>; 
  lctStation      = new vector<int>; 
  lctRing         = new vector<int>; 
  lctChamber      = new vector<int>; 
  lctTriggerCSCID = new vector<int>; 
  lctFpga         = new vector<int>; 
  lctlocalPhi     = new vector<int>; 
  lctglobalPhi    = new vector<int>; 
  lctglobalEta    = new vector<int>; 
  lctstripNum     = new vector<int>; 
  lctwireGroup    = new vector<int>;   
  lctWheel    = new vector<int>;
  isdtlct         = new vector<int>;
  trig_prim_phi   = new vector<double>;
  trig_prim_eta   = new vector<double>;
}

void TrigEff::csctfDel() {	

  vector<int>().swap(*EndcapTrk); 
  vector<int>().swap(*SectorTrk); 
  vector<int>().swap(*BxTrk    ); 

  vector<int>().swap(*me1ID); 
  vector<int>().swap(*me2ID); 
  vector<int>().swap(*me3ID); 
  vector<int>().swap(*me4ID); 
  vector<int>().swap(*mb1ID); 

  vector<int>().swap(*OutputLinkTrk);

  vector<int>  ().swap(*ModeTrk    ); 
  vector<float>().swap(*EtaTrk     ); 
  vector<float>().swap(*PhiTrk     ); 
  vector<float>().swap(*PhiTrk_02PI); 
  vector<float>().swap(*PtTrk      ); 

  vector<int>().swap(*ChargeTrk      ); 
  vector<int>().swap(*ChargeValidTrk ); 
  vector<int>().swap(*QualityTrk     ); 
  vector<int>().swap(*ForRTrk        ); 
  vector<int>().swap(*Phi23Trk       ); 
  vector<int>().swap(*Phi12Trk       ); 
  vector<int>().swap(*PhiSignTrk     ); 

  vector<int>().swap(*EtaBitTrk);
  vector<int>().swap(*PhiBitTrk);
  vector<int>().swap(*PtBitTrk );

  vector<int>().swap(*NumLCTsTrk );

  vector<int>().swap(*lctEndcap); 
  vector<int>().swap(*lctSector); 
  vector<int>().swap(*lctSubSector); 
  vector<int>().swap(*lctBx); 
  vector<int>().swap(*lctBx0); 
  vector<int>().swap(*lctStation); 
  vector<int>().swap(*lctRing); 
  vector<int>().swap(*lctChamber); 
  vector<int>().swap(*lctTriggerCSCID); 
  vector<int>().swap(*lctFpga);     
  vector<int>().swap(*lctlocalPhi); 
  vector<int>().swap(*lctglobalPhi);   
  vector<int>().swap(*lctglobalEta); 
  vector<int>().swap(*lctstripNum);   
  vector<int>().swap(*lctwireGroup);   
  vector<int>().swap(*lctWheel);
  vector<int>().swap(*isdtlct);
  vector<double>().swap(*trig_prim_phi);
  vector<double>().swap(*trig_prim_eta);
  
}

void TrigEff::fillCSCTF(const edm::Handle< vector<L1TMuon::InternalTrack> > tracks,
                        const L1MuTriggerScales  *ts, 
                        const L1MuTriggerPtScale *tpts, 
                        CSCSectorReceiverLUT* srLUTs_[5][2]
			) {
  
  // loop over CSCTF tracks
  int nTrk=0; 
  for( auto trk = tracks->cbegin(); trk < tracks->cend(); trk++){
    if (printLevel > 1) cout << "Looping over track # " << nTrk << endl;
    
    nTrk++;

    // Access track variables

    // Standard Pt LUTs	  
    edm::ParameterSet ptLUTset;
    ptLUTset.addParameter<bool>("ReadPtLUT", false);
    ptLUTset.addParameter<bool>("isBinary",  false);
    CSCTFPtLUT ptLUT(ptLUTset, ts, tpts);

    // Access the track variables in bit form
    int track_pt         = trk -> pt_packed();
    float trEta_bit      = trk -> eta_packed();
    unsigned long trMode = trk -> cscMode(); 

    // convert the variables in human readable values
    double trPt  = tpts->getPtScale()->getLowEdge(track_pt);
    double trEta = ts->getRegionalEtaScale(2)->getCenter(trEta_bit);
    
    if (printLevel > 1) {
      cout << "===========================" << endl;
      cout << " Track Pt   = " << trPt   << endl;
      cout << " Track mode = " << trMode << endl;
      cout << " Track Eta  = " << trEta  << endl;
    }
    
    // fill the track vectors
    PtTrk   -> push_back(trPt);
    EtaTrk  -> push_back(trEta);
    ModeTrk -> push_back(trMode);

    // Comment out rest of this info as it is not needed at this time

    //    int    trSector = trk -> getGlobalSector();
      /* int    trSector = 6*(trk->first.endcap()-1)+trk->first.sector();
    int    trBX     = trk->first.BX();
    int    trMe1ID = trk->first.me1ID();
    int    trMe2ID = trk->first.me2ID();
    int    trMe3ID = trk->first.me3ID();
    int    trMe4ID = trk->first.me4ID();
    int    trMb1ID = trk->first.mb1ID();
    
    int    trEtaBit = trk->first.eta_packed();
    int    trPhiBit = trk->first.localPhi();
    int    trCharge = trk->first.chargeValue();
    
    int    trRank   = trk->first.rank();    
      
    // PtAddress gives an handle on other parameters
    ptadd thePtAddress(trk->first.ptLUTAddress());
    
    int trPhiSign = thePtAddress.delta_phi_sign;
    int trPhi12   = thePtAddress.delta_phi_12;
    int trPhi23   = thePtAddress.delta_phi_23;
    int trMode    = thePtAddress.track_mode;
    int trForR    = thePtAddress.track_fr;
    
    
    int trOutputLink = trk->first.outputLink();
    
    //Pt needs some more workaround since it is not in the unpacked data
    //ptdat thePtData  = ptLUT->Pt(thePtAddress);
    ptdat thePtData  = ptLUT.Pt(thePtAddress);
    
    int trPtBit       = 0;
    int trQuality     = 0;
    int trChargeValid = 0;
    
    // front or rear bit? 
    if (thePtAddress.track_fr) {
      trPtBit = (thePtData.front_rank&0x1f);
      trQuality = ((thePtData.front_rank>>5)&0x3);
      trChargeValid = thePtData.charge_valid_front;
    } else {
      trPtBit = (thePtData.rear_rank&0x1f);
      trQuality = ((thePtData.rear_rank>>5)&0x3);
      trChargeValid = thePtData.charge_valid_rear;
    }
      
    if (printLevel >=1) {
      cout << " ===== CSCTF BIT VALUES ====\n"
           << " Track #"        << nTrk
           << "\n Endcap: "     << trEndcap
           << "\n Sector: "     << trSector 
           << "\n BX: "         << trBX 
           << "\n me1ID: "      << trMe1ID 			
	   << "\n me2ID: "      << trMe2ID
	   << "\n me3ID: "      << trMe3ID
	   << "\n me4ID: "      << trMe4ID 
           << "\n mb1ID: "      << trMb1ID 
           << "\n OutputLink: " << trOutputLink
        
           << "\n Charge: "     << trCharge
           << "\n DPhiSign: "   << trPhiSign
	   << "\n DPhi12: "     << trPhi12
	   << "\n DPhi23: "     << trPhi23
	   << "\n ForR: "       << trForR
                                
           << "\n Mode: "       << trMode 
           << "\n Quality: "    << trQuality 
           << "\n Rank : "      << trRank 
           << "\n ChargeValid: "<< trChargeValid
           << "\n Eta: "        << trEtaBit 
           << "\n Phi: "        << trPhiBit
           << "\n Pt: "         << trPtBit
           << endl;         
    }

    //... in radians
    // Type 2 is CSC
    double trEta = ts->getRegionalEtaScale(2)->getCenter(trk->first.eta_packed());
    double trPhi = ts->getPh;'.l/;"


    bnmiScale()->getLowEdge(trk->first.localPhi());
    //Phi in one sector varies from [0,62] degrees -> Rescale manually to global coords.
    double trPhi02PI = fmod(trPhi + 
                            ((trSector-1)*TMath::Pi()/3) + //sector 1 starts at 15 degrees 
                            (TMath::Pi()/12) , 2*TMath::Pi());
    
    // convert the Pt in human readable values (GeV/c)
    double trPt = tpts->getPtScale()->getLowEdge(trPtBit); 


    if (printLevel >=1) {
      cout << " ===== CSCTF TRK SCALES ====\n"
           << "\n Eta(scales): " << trEta 
           << "\n Phi(scales): " << trPhi02PI
           << "\n Pt(scales): "  << trPt
           << endl;
    }

    // fill the track vectors
    EndcapTrk      -> push_back(trEndcap);
    SectorTrk      -> push_back(trSector); 
    BxTrk          -> push_back(trBX    ); 
    me1ID          -> push_back(trMe1ID);  
    me2ID          -> push_back(trMe2ID);
    me3ID          -> push_back(trMe3ID); 
    me4ID          -> push_back(trMe4ID); 
    mb1ID          -> push_back(trMb1ID);  
    OutputLinkTrk  -> push_back(trOutputLink);  

    ModeTrk        -> push_back(trMode   );
    EtaTrk         -> push_back(trEta    ); 
    PhiTrk         -> push_back(trPhi    ); 
    PhiTrk_02PI    -> push_back(trPhi02PI); 
    PtTrk          -> push_back(trPt     );  

    ChargeValidTrk -> push_back(trChargeValid);
    ChargeTrk      -> push_back(trCharge     );
    QualityTrk     -> push_back(trQuality    );
    ForRTrk        -> push_back(trForR       );
    Phi23Trk       -> push_back(trPhi23      );
    Phi12Trk       -> push_back(trPhi12      );
    PhiSignTrk     -> push_back(trPhiSign    );
    
    EtaBitTrk      -> push_back(trEtaBit); 
    PhiBitTrk      -> push_back(trPhiBit); 
    PtBitTrk       -> push_back(trPtBit );  

    */
   
    


    
    // For each trk, get the list of its LCTs.  Use Lindsey's matched objects
    if (printLevel > 1) cout << "---------- Loop over LCTs matched to Track -----------\n" ;
    
    auto lct_map = trk -> getStubs();  // Get map for LCTs to stations
    //unsigned csc_mode = trk->cscMode();  // map needs to know what subsystem
    
    if (printLevel > 2) cout << "  Size of Lct Map is " << lct_map.size() << endl; // Tells us how many stations actually have lcts
    
    int LctTrkId_ = 0;   
    // Loop over stations and get all stubs in each one
    for( unsigned station = 1; station <= 4; ++station ) {
      if (printLevel > 1) cout << " \n Looping over station " << station << " to get LCTs" << endl;
      
      const unsigned id    = 4*L1TMuon::InternalTrack::kCSC+station-1; // unique identifier for each csc station
      const unsigned dt_id = 4*L1TMuon::InternalTrack::kDT+station-1; // unique identifier for each dt station
      
      if( !lct_map.count(id) ) {
        if (printLevel > 1) cout << "  CSC Station map empty " << endl;
      }
      if ( !lct_map.count(dt_id) ) {
	if (printLevel > 1) cout << "  DT Station map empty " << endl;
      }
      
      // ----- First loop over csc stubs, then dt later ----- 
      auto x_LCTs = lct_map[id];  // access the matched lcts for a given station id
      // Loop over lcts in each station
      for ( auto t_lcts = x_LCTs.cbegin(); t_lcts != x_LCTs.cend(); t_lcts++ ) {
	auto lcts = *t_lcts; // dereference the edm:Ref object	
	CSCDetId id = lcts->detId<CSCDetId>();
	
	auto trlct_station       = id.station();
        auto trlct_endcap        = id.endcap();
        double trlct_phi         = lcts->getCMSGlobalPhi();
        double trlct_eta         = lcts->getCMSGlobalEta();
        uint16_t trlct_bx        = lcts->getCSCData().bx;
	//unsigned trlct_sector    = lcts->getGlobalSector();
	//unsigned trlct_subsector = lcts->getSubSector();
	int trlct_sector         = CSCTriggerNumbering::triggerSectorFromLabels(id);
	int trlct_subsector      = CSCTriggerNumbering::triggerSubSectorFromLabels(id);
	uint16_t trlct_bx0       = lcts->getCSCData().bx0; 
	uint16_t trlct_cscID     = lcts->getCSCData().cscID;
	uint16_t trlct_strip     = lcts->getCSCData().strip;
	uint16_t trlct_pattern   = lcts->getCSCData().pattern;
	uint16_t trlct_bend      = lcts->getCSCData().bend;
	uint16_t trlct_quality   = lcts->getCSCData().quality;
	uint16_t trlct_keywire   = lcts->getCSCData().keywire;
	
        if (printLevel > 1) { 
	  cout << "  CSC LCT #     = "     << LctTrkId_ << endl;
	  cout << "  Phi       = "      << trlct_phi << endl;
	  cout << "  Eta       = "      << trlct_eta << endl;
	  cout << "  Bunch X   = "  << trlct_bx  << endl;
	  cout << "  Sector    = "   << trlct_sector  << endl;
	  cout << "  SubSector = " << trlct_subsector << endl;
	  cout << "  Endcap    = "   << trlct_endcap  << endl;
	  cout << "  Station   = "  << trlct_station  << endl;
	  cout << "  cscID     = "    << trlct_cscID  << endl << endl;
	}
	//------------- Fill ntuple with LCT values ------------
	
	if ( (nTrk-1) > (MAX_CSCTF_TRK-1) || LctTrkId_ > (MAX_LCTS_PER_TRK-1) ) {
          cout << "the track has " << nTrk-1 << " tracks, and " << LctTrkId_+1 << "Lcts, but the MAX allowed tracks is "
               << MAX_CSCTF_TRK << " and the max LCTs is " << MAX_LCTS_PER_TRK << " , -> Skipping the rest of tracks... " << endl;
          continue;
        }

	// Check if DetId is within range
        if(   trlct_sector < 1    || trlct_sector > 12 || 
	      trlct_station < 1   || trlct_station >  4 || 
	      trlct_cscID < 1     || trlct_cscID >  9 ||
	      trlct_endcap > 2    || trlct_endcap < 1  ) {
          if (printLevel > 1) cout << "  TRACK ERROR: CSC digi are out of range! " << endl;
          continue;
        }
	
	isdttrlct[nTrk-1][LctTrkId_] = 0;

	if (printLevel > 3) cout << "Filling endcap branch" << endl;
	trLctEndcap[nTrk-1][LctTrkId_] = trlct_endcap;
        
	// sector (end 1: 1->6, end 2: 7 -> 12)
        if (printLevel > 1) cout << "Filling sector branch" << endl;
	if ( trlct_endcap == 1)
          trLctSector[nTrk-1][LctTrkId_] = trlct_sector;
        else
          trLctSector[nTrk-1][LctTrkId_] = 6+trlct_sector;
	
	if (printLevel > 3) cout << "Filling subsector branch" << endl;
        trLctSubSector[nTrk-1][LctTrkId_] = trlct_subsector;
	
	if (printLevel > 3) cout << "Filling bx branch" << endl;
        trLctBx[nTrk-1][LctTrkId_] = trlct_bx;
        trLctBx0[nTrk-1][LctTrkId_] = trlct_bx0;
	
	if (printLevel > 3) cout << "Filling station branch" << endl;
        trLctStation[nTrk-1][LctTrkId_] = trlct_station;
	
	//trLctRing[nTrk-1][LctTrkId_] = (*lctOfTrks).first.ring();
	
        //trLctChamber[nTrk-1][LctTrkId_] = (*lctOfTrks).first.chamber();
	
	if (printLevel > 3) cout << "Filling cscID branch" << endl;
        trLctTriggerCSCID[nTrk-1][LctTrkId_] = trlct_cscID;
		
        // handles not to overload the method: mostly for readability	      
        if (trlct_endcap == 2) trlct_endcap = 0; 
        
	int FPGALctTrk = ( trlct_subsector ? trlct_subsector-1 : trlct_station );
	// local Phi
        lclphidat lclPhi;	
        try {
          
	  trLctstripNum[nTrk-1][LctTrkId_] = trlct_strip;
          lclPhi = srLUTs_[FPGALctTrk][trlct_endcap] -> localPhi(trlct_strip, 
								 trlct_pattern, 
								 trlct_quality, 
								 trlct_bend );
          
          trLctlocalPhi[nTrk-1][LctTrkId_] = lclPhi.phi_local;
        } 
        catch(...) { 
          bzero(&lclPhi,sizeof(lclPhi)); 
          trLctlocalPhi[nTrk-1][LctTrkId_] = -999;
        }


	  // Fill phi bit
	  if (printLevel > 3) cout << "Filling global phi in bit form " << endl;
	  
	  gblphidat gblPhi;
	  
	  try {
	    
	    lctwireGroup->push_back(trlct_keywire);
	    gblPhi = srLUTs_[FPGALctTrk][trlct_endcap] -> globalPhiME(lclPhi.phi_local,
								      trlct_keywire,
								      trlct_cscID);
 	    trLctglobalPhi[nTrk-1][LctTrkId_] = gblPhi.global_phi;
	    
	  } catch(...) {
	    bzero(&gblPhi,sizeof(gblPhi));
	
	    trLctglobalPhi[nTrk-1][LctTrkId_]  = -999;
	  }
	  
	  // Fill eta bit
	  if (printLevel > 3) cout << "Filling global eta in bit form " << endl;
	  gbletadat gblEta;
	  
	  try {
	    gblEta = srLUTs_[FPGALctTrk][trlct_endcap] -> globalEtaME(lclPhi.phi_bend_local,
								      lclPhi.phi_local,
								      trlct_keywire,
								      trlct_cscID);
	    trLctglobalEta[nTrk-1][LctTrkId_] = gblEta.global_eta;
	  }
	  catch(...) {
	    bzero(&gblEta,sizeof(gblEta));
	    trLctglobalEta[nTrk-1][LctTrkId_] = -999;
	  }
	  
	  if ( printLevel > 1 ) {
	    cout << "gblPhi = " << gblPhi.global_phi << endl;
	    cout << "gblEta = " << gblEta.global_eta << endl;
	  }
	
	LctTrkId_ += 1;

      } // end loop over CSC LCTs

      // ------- Now loop over any DT Lcts -------
      

      auto y_LCTs = lct_map[dt_id]; // access the matched lcts for a given station id
      // Loop over lcts in each station
      for ( auto t_lcts = y_LCTs.cbegin(); t_lcts != y_LCTs.cend(); t_lcts++ ) {
        auto lcts = *t_lcts; // dereference the edm:Ref object
        DTChamberId id = lcts->detId<DTChamberId>();
	
	int dt_trlct_station        = lcts->getDTData().station;
        double dt_trlct_phi         = lcts->getCMSGlobalPhi();
        double dt_trlct_eta         = lcts->getCMSGlobalEta();
        uint16_t dt_trlct_bx        = lcts->getDTData().bx;
        unsigned dt_trlct_sector    = lcts->getDTData().sector;
	auto dt_trlct_phi_bit       = lcts->getDTData().radialAngle;

	// convert dt phi bit into csc phi bit
	if (dt_trlct_phi_bit < 0) dt_trlct_phi_bit += 4096;

	if (dt_trlct_phi_bit > 4096) {
	  if (printLevel > 1) cout << "AAAAAAAAAAGH TOO BIG PHI:" << dt_trlct_phi_bit << endl;
	  continue;
	}
	if (dt_trlct_phi_bit < 0) {
	  if (printLevel > 1) cout << "AAAAAAAAH NEG PHI" << dt_trlct_phi_bit << endl;
	  continue;
	}

	int dt_to_csc_phi_bit = CSCTFDTReceiverLUT::lut[dt_trlct_phi_bit];

        if (printLevel > 1) {
          cout << "  DT LCT #  = "  << LctTrkId_         << endl;
	  cout << "  ChamberId = "  << id                << endl;
	  cout << "  Phi       = "  << dt_trlct_phi      << endl;
          cout << "  Phi Bit(before convert)       = "  << dt_trlct_phi_bit << endl;
	  cout << "  Phi Bit(after convert)       = "  << dt_to_csc_phi_bit << endl;
	  cout << "  Eta       = "  << dt_trlct_eta      << endl;
          cout << "  Bunch X   = "  << dt_trlct_bx       << endl;
          cout << "  Sector    = "  << dt_trlct_sector   << endl;
          cout << "  Station   = "  << dt_trlct_station  << endl;
	  cout << "  Subsystem  = " << lcts->subsystem() << endl;
	}
        //------------- Fill ntuple with LCT values ------------

	if ( (nTrk-1) > (MAX_CSCTF_TRK-1) || LctTrkId_ > (MAX_LCTS_PER_TRK-1) ) {
          cout << "the track has " << nTrk-1 << " tracks, and " << LctTrkId_+1 << "Lcts, but the MAX allowed tracks is "
               << MAX_CSCTF_TRK << " and the max LCTs is " << MAX_LCTS_PER_TRK << " , -> Skipping the rest of tracks... " << endl;
          continue;
        }

	isdttrlct[nTrk-1][LctTrkId_] = 1;

        // sector 
        if (printLevel > 3) cout << "Filling sector branch" << endl;        
	dt_trLctSector[nTrk-1][LctTrkId_] = dt_trlct_sector;
               
        if (printLevel > 3) cout << "Filling bx branch" << endl;
	dt_trLctBx[nTrk-1][LctTrkId_] = dt_trlct_bx;
       
        if (printLevel > 3) cout << "Filling station branch" << endl;
        dt_trLctStation[nTrk-1][LctTrkId_] = dt_trlct_station;

	dt_trLctChamber[nTrk-1][LctTrkId_] = id;

	trLctglobalPhi[nTrk-1][LctTrkId_] = dt_to_csc_phi_bit;

	LctTrkId_ += 1;
      } // end loop over DT LCTs

    } // end loop over stations
    NumLCTsTrk->push_back(LctTrkId_);
    if (printLevel > 1) cout << "# of Lcts in this track = " << LctTrkId_ << endl; 
    
  } // End for loop on tracks
  SizeTrk = nTrk;
  
} // End void FillCSCTF

void TrigEff::fillAllLCTs(const edm::Handle< vector<L1TMuon::TriggerPrimitive> > lcts,
                          CSCSectorReceiverLUT* srLUTs_[5][2] 
			  ) {
  
  int lctId = 0; // count number of lcts in event
  for ( auto Lct = lcts->cbegin(); Lct != lcts->cend(); Lct++) {
    // ------------- Fill all LCT variables ------------------
    //cout << "subsytem type is " << Lct->subsystem() << endl;
    if ( Lct->subsystem() == 1 ) {
      CSCDetId id = Lct->detId<CSCDetId>();
     
      isdtlct -> push_back(0);
     
      auto lct_station           = id.station();
      auto lct_endcap            = id.endcap();
      uint16_t lct_bx            = Lct->getCSCData().bx;
      int lct_ring               = id.ring();
      int lct_sector             = CSCTriggerNumbering::triggerSectorFromLabels(id);
      int lct_subsector          = CSCTriggerNumbering::triggerSubSectorFromLabels(id);
      uint16_t lct_bx0           = Lct->getCSCData().bx0;
      uint16_t lct_cscID         = Lct->getCSCData().cscID;
      uint16_t lct_strip         = Lct->getCSCData().strip;
      uint16_t lct_pattern       = Lct->getCSCData().pattern;
      uint16_t lct_bend          = Lct->getCSCData().bend;
      uint16_t lct_quality       = Lct->getCSCData().quality;
      uint16_t lct_keywire       = Lct->getCSCData().keywire;
      double lct_phi             = Lct->getCMSGlobalPhi();

      
      if ( printLevel > 1 ) {
	cout << "\nLCT CSC " << lctId
	     << "\n======\n";
	cout <<"lctEndcap       = " << lct_endcap << endl;
	cout <<"lctSector       = " << lct_sector<< endl;
	cout <<"lctSubSector    = " << lct_subsector << endl;
	cout <<"lctStation      = " << lct_station << endl;
	cout <<"lctRing         = " << lct_ring << endl;
	//cout <<"lctChamber      = " << (*corrLct).first.chamber() << std::endl;
	cout <<"lctTriggerCSCID = " << lct_cscID << endl;
	cout <<"lctBx           = " << lct_bx << endl;
	cout <<"lctKeyWire      = " << lct_keywire << endl;
	cout <<"lctStrip        = " << lct_strip << endl;      
	cout <<"lctBend         = " << lct_bend << endl;
	cout <<"lctQuality      = " << lct_quality << endl;
      }
      
      //    Check if DetId is within range
      if( lct_sector < 1   || lct_sector > 12 || 
	  lct_station < 1  || lct_station >  4 || 
	  lct_cscID < 1    || lct_cscID >  9 ||
	  lct_endcap > 2   || lct_endcap < 1 )
	{
	  if (printLevel > 1) cout << "  LCT ERROR: CSC digi are out of range! " << endl;;
	  continue;
	}

      // Fill with dummy values for csc Lcts. 
      lctChamber->push_back(-999);
      lctWheel->push_back(-999);
      trig_prim_phi->push_back(lct_phi);
      trig_prim_eta->push_back(-999);

      //if (printLevel > 1) cout << "Filling endcap branch" << endl;
      lctEndcap->push_back( lct_endcap );
      
      //if (printLevel > 1) cout << "Filling sector branch" << endl;
      if ( lct_endcap == 1)
	lctSector->push_back( lct_sector );
      else
	lctSector->push_back( 6+(lct_sector) );
      
      //if (printLevel > 1) cout << "Filling subsector branch" << endl;
      lctSubSector->push_back(lct_subsector);
      
      //if (printLevel > 1) cout << "Filling bx branch" << endl;
      lctBx->push_back(lct_bx);
      lctBx0->push_back(lct_bx0);
      
      //if (printLevel > 1) cout << "Filling station branch" << endl;
      lctStation->push_back( lct_station );
      
      // lctRing->push_back((*corrLct).first.ring());
      
      // lctChamber->push_back((*corrLct).first.chamber());
      
      //if (printLevel > 1) cout << "Filling cscID branch" << endl;
      lctTriggerCSCID->push_back( lct_cscID );
      
      int FPGALct = ( lct_subsector ? lct_subsector-1 : lct_station );
      if (lct_endcap == 2) lct_endcap = 0;
      
      // local Phi
      lclphidat lclPhi;
      
      try {
	lclPhi = srLUTs_[FPGALct][lct_endcap] -> localPhi(lct_strip, 
							  lct_pattern, 
							  lct_quality, 
							  lct_bend);
      
	lctlocalPhi->push_back(lclPhi.phi_local);
      } 
      catch(...) { 
	bzero(&lclPhi,sizeof(lclPhi)); 
	lctlocalPhi->push_back(-999);
      }	
      
      // Fill phi bit
      
      //if (printLevel > 1) cout << "Filling global phi in bit form " << endl;  
      gblphidat gblPhi;
      
      try {
	
	lctwireGroup->push_back(lct_keywire);
	
	gblPhi = srLUTs_[FPGALct][lct_endcap] -> globalPhiME(lclPhi.phi_local,
							     lct_keywire,
							     lct_cscID);
	
	lctglobalPhi->push_back(gblPhi.global_phi);
	
      } catch(...) {
	bzero(&gblPhi,sizeof(gblPhi));
	lctglobalPhi->push_back(-999);
      }
      
      // Fill eta bit
      //if (printLevel > 1) cout << "Filling global eta in bit form " << endl;
      gbletadat gblEta;
      
      try {
	gblEta = srLUTs_[FPGALct][lct_endcap] -> globalEtaME(lclPhi.phi_bend_local,
							     lclPhi.phi_local,
							     lct_keywire,
							     lct_cscID);
	lctglobalEta->push_back(gblEta.global_eta);
      }
      catch(...) {
	bzero(&gblEta,sizeof(gblEta));
	lctglobalEta->push_back(-999);
      }
      
      if ( printLevel > 1 ) {
	cout << "gblPhi = " << gblPhi.global_phi << endl;
	cout << "gblEta = " << gblEta.global_eta << endl << endl;
      }
 
      lctId += 1;
    } // end CSC LCT 
    
    if ( Lct->subsystem() == 0 ) {

      // ----- Fill all DT LCTs -------------------------------------------------------------------------------
      DTChamberId id = Lct -> detId<DTChamberId>();
      
      int dt_lct_station        = Lct->getDTData().station;
      double dt_lct_phi         = Lct->getCMSGlobalPhi();
      double dt_lct_eta         = Lct->getCMSGlobalEta();
      uint16_t dt_lct_bx        = Lct->getDTData().bx;
      unsigned dt_lct_sector    = Lct->getDTData().sector;
      auto dt_lct_wheel     = Lct->getDTData().wheel;
      auto dt_lct_phi_bit   = Lct->getDTData().radialAngle;
      
      // convert dt phi bit into csc phi bit
      if (dt_lct_phi_bit < 0) dt_lct_phi_bit += 4096;
      
      if (dt_lct_phi_bit > 4096) { 
	if (printLevel > 1) cout << "AAAAAAAAAAGH TOO BIG PHI:" << dt_lct_phi_bit << endl;
	continue; 
      }
      if (dt_lct_phi_bit < 0) {
	if (printLevel > 1) cout << "AAAAAAAAH NEG PHI" << dt_lct_phi_bit << endl;
	continue;
      }
      
      int dt_to_csc_phi_bit = CSCTFDTReceiverLUT::lut[dt_lct_phi_bit];
      

      if (printLevel > 1) {
	cout << "\nLCT DT " << lctId
	     << "\n======\n";
	cout << "  ChamberId = "  << id              << endl;
	cout << "  Phi       = "  << dt_lct_phi      << endl;
	cout << "  Phi Bit(before convert)   = "     << dt_lct_phi_bit  << endl;
	cout << "  Phi Bit(after convert)    = "     <<  dt_to_csc_phi_bit << endl;
	cout << "  Bunch X   = "  << dt_lct_bx       << endl;
	cout << "  Sector    = "  << dt_lct_sector   << endl;
	cout << "  Station   = "  << dt_lct_station  << endl;
	cout << "  Wheel     = "  << dt_lct_wheel    << endl;
      }
      
      //------------- Fill ntuple with LCT values ------------
      lctEndcap -> push_back(-999);
      
      isdtlct -> push_back(1); 
      
      trig_prim_phi ->push_back(dt_lct_phi);
      trig_prim_eta ->push_back(dt_lct_eta);
      
      //if (printLevel > 1) cout << "Filling sector branch" << endl;
      lctSector->push_back(dt_lct_sector);
       
      //if (printLevel > 1) cout << "Filling bx branch" << endl;
      lctBx->push_back(dt_lct_bx);
      
      //if (printLevel > 1) cout << "Filling station branch" << endl;
      lctStation->push_back(dt_lct_station);
      
      lctChamber->push_back(id);
      
      lctWheel->push_back(dt_lct_wheel);
      
      lctglobalPhi->push_back(dt_to_csc_phi_bit);
      lctglobalEta->push_back(-999);
      
      lctId += 1;
    } // end loop over DT LCTs
    
  } // End loop over lcts
  SizeLCTs = lctId;
  
}  //End FillLCTs

// this code snippet was taken from 
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/L1TriggerDPG/src/L1MuonRecoTreeProducer.cc?revision=1.6&view=markup
// to get track position at a particular (xy) plane given its z
TrajectoryStateOnSurface TrigEff::surfExtrapTrkSam(reco::TrackRef track, double z)
{
  Plane::PositionType pos(0, 0, z);
  Plane::RotationType rot;
  Plane::PlanePointer myPlane = Plane::build(pos, rot);

  FreeTrajectoryState recoStart = freeTrajStateMuon(track);
  TrajectoryStateOnSurface recoProp;
  recoProp = propagatorAlong->propagate(recoStart, *myPlane);
  if (!recoProp.isValid()) {
    recoProp = propagatorOpposite->propagate(recoStart, *myPlane);
  }
  return recoProp;
}

FreeTrajectoryState TrigEff::freeTrajStateMuon(reco::TrackRef track)
{
  GlobalPoint  innerPoint(track->innerPosition().x(),  track->innerPosition().y(),  track->innerPosition().z());
  GlobalVector innerVec  (track->innerMomentum().x(),  track->innerMomentum().y(),  track->innerMomentum().z());  
  
  FreeTrajectoryState recoStart(innerPoint, innerVec, track->charge(), &*theBField);
  
  return recoStart;
}


void TrigEff::fillSegments(edm::Handle<CSCSegmentCollection> cscSegments, 
			   edm::Handle<DTRecSegment4DCollection> dtSegments,
			   edm::ESHandle<CSCGeometry> cscGeom,
			   edm::ESHandle<DTGeometry> dtGeom){
  
  // get CSC segment collection
  int csc_segsSize = cscSegments->size();
  int dt_segsSize  = dtSegments ->size();
  
  if (printLevel>=1) {
    cout << "Found " << csc_segsSize << " CSC segments" << endl;
    cout << "Found " << dt_segsSize  << " Dt  segments" << endl;
  }
  
  int sumsegs = csc_segsSize + dt_segsSize;
  segsSize = sumsegs;  
  int iSegment = 0;  
  // -----------------------
  // loop over CSC segments
  // -----------------------
  for (auto dSiter = cscSegments->begin(); dSiter != cscSegments->end(); dSiter++, iSegment++) {
    
    // ----- Fill CSC Seg values -----
    CSCDetId id  = (CSCDetId)(*dSiter).cscDetId();
    int kEndcap  = id.endcap();
    int kRing    = id.ring();
    int kStation = id.station();
    int kChamber = id.chamber();
    int notdt = 0;

    isdtseg         ->push_back(notdt);
    cscsegs_endcap  ->push_back(kEndcap);
    cscsegs_station ->push_back(kRing);
    cscsegs_ring    ->push_back(kStation);
    cscsegs_chamber ->push_back(kChamber);
    cscsegs_wheel   ->push_back(-999);
    cscsegs_sector  -> push_back(-999);

    LocalPoint localPos = (*dSiter).localPosition();
    float segX = localPos.x();
    float segY = localPos.y();
    float segZ = localPos.z();
    float segPhi = localPos.phi();
    float segEta = localPos.eta();
    float segTheta = localPos.theta();
      
    cscsegs_loc_x->push_back(segX);
    cscsegs_loc_y->push_back(segY);
    cscsegs_loc_z->push_back(segZ);

    cscsegs_loc_theta->push_back(segTheta);
    cscsegs_loc_eta  ->push_back(segEta);
    cscsegs_loc_phi  ->push_back(segPhi);
      
    LocalVector segDir = (*dSiter).localDirection();
    double dirTheta = segDir.theta();
    double dirEta   = segDir.eta();
    double dirPhi   = segDir.phi();

    cscsegs_loc_dir_theta->push_back(dirTheta);
    cscsegs_loc_dir_eta  ->push_back(dirEta);
    cscsegs_loc_dir_phi  ->push_back(dirPhi);

    // global transformation
    float globX = 0.;
    float globY = 0.;
    float globZ = 0.;
    float globSegPhi   = 0.;
    float globSegTheta = 0.;
    float globSegEta   = 0.;
    float globDirPhi   = 0.;
    float globDirTheta = 0.;
    float globDirEta   = 0.;
    
    const CSCChamber* cscchamber = cscGeom->chamber(id);
    if (cscchamber) {
      GlobalPoint globalPosition = cscchamber->toGlobal(localPos);
      globX        = globalPosition.x();
      globY        = globalPosition.y();
      globZ        = globalPosition.z();
      globSegPhi   = globalPosition.phi();
      globSegEta   = globalPosition.eta();
      globSegTheta = globalPosition.theta();

      GlobalVector globalDirection = cscchamber->toGlobal(segDir);
      globDirTheta = globalDirection.theta();
      globDirEta   = globalDirection.eta();
      globDirPhi   = globalDirection.phi();

      cscsegs_gbl_x->push_back(globX);
      cscsegs_gbl_y->push_back(globY);
      cscsegs_gbl_z->push_back(globZ);

      cscsegs_gbl_theta->push_back(globSegTheta);
      cscsegs_gbl_eta  ->push_back(globSegEta  );
      cscsegs_gbl_phi  ->push_back(globSegPhi  );

      cscsegs_gbl_dir_theta->push_back(globDirTheta);
      cscsegs_gbl_dir_eta  ->push_back(globDirEta  );
      cscsegs_gbl_dir_phi  ->push_back(globDirPhi  );

      if (printLevel > 1) {
        cout << "\nCSC S E G M E N T #" << iSegment << endl;
        cout << "*******************"<< endl;
        
        cout << "Endcap : "<< kEndcap  << endl;
        cout << "Ring   : "<< kRing    << endl;
        cout << "Station: "<< kStation << endl;
        cout << "Chamber: "<< kChamber << endl;

        cout << "globX  : "    << globX      << endl;
        cout << "globY  : "    << globY      << endl;
        cout << "globZ  : "    << globZ      << endl;
        cout << "globSegPhi: " << globSegPhi << endl;
        cout << "globSegEta: " << globSegEta << endl;
        
        if (printLevel > 3) {
          cout << "segX: "<< segX << endl;
          cout << "segY: "<< segY << endl;
          cout << "segZ: "<< segZ << endl;
          cout << "segPhi: "<< segPhi << endl;
          cout << "segEta: "<< segEta << endl;
          cout << "segTheta: "<< segTheta << endl;
          cout << "dirTheta: "<< dirTheta << endl;
          cout << "dirEta: "  << dirEta   << endl;
          cout << "dirPhi: "  << dirPhi   << endl;
        }
      }
    } 
    else {
      cscsegs_gbl_x->push_back(-999);
      cscsegs_gbl_y->push_back(-999);
      cscsegs_gbl_z->push_back(-999);

      cscsegs_gbl_theta->push_back(-999);
      cscsegs_gbl_eta  ->push_back(-999);
      cscsegs_gbl_phi  ->push_back(-999);

      cscsegs_gbl_dir_theta->push_back(-999);
      cscsegs_gbl_dir_eta  ->push_back(-999);
      cscsegs_gbl_dir_phi  ->push_back(-999);
    }
  } // end loop over csc segs
  
  // -----------------------
  // loop over DT segments
  // -----------------------
  for (auto dSiter = dtSegments->begin(); dSiter != dtSegments->end(); dSiter++, iSegment++) {

    // ----- Fill DT Seg values -----
    // Note:  We will be filling cscsegs with DT values.  This is ok and best way to access data.
    DTChamberId id = (DTChamberId)(*dSiter).chamberId();
    int kStation   = id.station();
    int kSector    = id.sector();
    int kWheel     = id.wheel();
    int isdt       = 1;

    isdtseg        -> push_back(isdt);
    cscsegs_chamber -> push_back(id);
    cscsegs_station -> push_back(kStation);
    cscsegs_wheel   -> push_back(kWheel);
    cscsegs_sector  -> push_back(kSector);
    
    LocalPoint localPos = (*dSiter).localPosition();
    float segX = localPos.x();
    float segY = localPos.y();
    float segZ = localPos.z();
    float segPhi = localPos.phi();
    float segEta = localPos.eta();

    /*LocalVector segDir = (*dSiter).localDirection();
    double dirTheta = segDir.theta();
    double dirEta   = segDir.eta();
    double dirPhi   = segDir.phi();*/

    // global transformation
    float globX = 0.;
    float globY = 0.;
    float globZ = 0.;
    float globSegPhi   = 0.;
    float globSegTheta = 0.;
    float globSegEta   = 0.;
    //float globDirPhi   = 0.;
    //float globDirTheta = 0.;
    //float globDirEta   = 0.;

    const DTChamber* dtchamber = dtGeom->chamber(id);
    if (dtchamber) {
      GlobalPoint globalPosition = dtchamber->toGlobal(localPos);
      globX        = globalPosition.x();
      globY        = globalPosition.y();
      globZ        = globalPosition.z();
      globSegPhi   = globalPosition.phi();
      globSegEta   = globalPosition.eta();
      globSegTheta = globalPosition.theta();
      
      /*GlobalVector globalDirection = cscchamber->toGlobal(segDir);
	globDirTheta = globalDirection.theta();
	globDirEta   = globalDirection.eta();
	globDirPhi   = globalDirection.phi(); */
      
      cscsegs_gbl_x -> push_back(globX);
      cscsegs_gbl_y -> push_back(globY);
      cscsegs_gbl_z -> push_back(globZ);
      cscsegs_gbl_theta -> push_back(globSegTheta);
      cscsegs_gbl_eta   -> push_back(globSegEta);
      cscsegs_gbl_phi   -> push_back(globSegPhi);
      
      
      if (printLevel > 1) {
	cout << "\nDT S E G M E N T #" << iSegment << endl;
	cout << "*******************"<< endl;
	cout << "ChamberId: "  << id    << endl;
	cout << "Station: "    << kStation << endl;
	cout << "Sector: "     << kSector << endl;
	cout << "Wheel: "      << kWheel << endl;
	cout << "globSegPhi: " << globSegPhi << endl;
	cout << "globSegEta: " << globSegEta << endl;
	
	if (printLevel > 3) {
	  cout << "segX: "<< segX << endl;
	  cout << "segY: "<< segY << endl;
	  cout << "segZ: "<< segZ << endl;
	  cout << "segPhi: "<< segPhi << endl;
	  cout << "segEta: "<< segEta << endl;
	} 
      }
      
    } // end if(dtchamber)    
  } //end loop over DT segs  
} //end function fillSegs

   
void TrigEff::fillSegmentsMuons ( const edm::Handle<reco::MuonCollection> muons,
                                  edm::Handle<CSCSegmentCollection> cscSegments, 
				  edm::Handle<DTRecSegment4DCollection> dtSegments,
				  edm::ESHandle<CSCGeometry> cscGeom,
				  const edm::Handle< vector<L1TMuon::TriggerPrimitive> > CSCTFlcts,
				  const edm::Handle< vector<L1TMuon::InternalTrack> > CSCTFtracks) {
    
  if (printLevel > 1) {
    printf("Found %d muons in the event\n\n", int(muons->size()) );
  }
  muonSize = muons->size();
  
  // loop over the muons
  int whichMuon =-1;
  for (auto muon=muons->begin(); muon!=muons->end(); muon++) {
    
    if (whichMuon > (MAX_MUONS-1) ) {
      cout << "the muon has " << whichMuon << ", but the MAX allowed is "
           << MAX_MUONS << " -> Skipping the muon... " << endl;
      continue;
    }

    whichMuon++;
    if (printLevel > 1) {
      printf("************************************************\n");
      printf("M U O N  #%d\n", whichMuon);
      printf("************************************************\n\n");
      printf("%s\n"    , "--------------------------------");
      printf("%s: %d\n", "isGlobalMuon    ()", muon->isGlobalMuon    ());
      printf("%s: %d\n", "isTrackerMuon   ()", muon->isTrackerMuon   ());
      printf("%s: %d\n", "isStandAloneMuon()", muon->isStandAloneMuon());
      printf("%s: %d\n", "combinedMuon    ().isNonnull()", muon->combinedMuon  ().isNonnull());
      printf("%s: %d\n", "track           ().isNonnull()", muon->track         ().isNonnull());
      printf("%s: %d\n", "standAloneMuon  ().isNonnull()", muon->standAloneMuon().isNonnull());
      printf("%s\n\n"  , "--------------------------------");
    }
    
    // given a muon candidate, loop over the segments and 
    // find the segments belonging to the muon and see 
    // if they are "LCTAble" (i.e. could generate an LCT)
    if (!muon->combinedMuon())   continue;
    if (!muon->standAloneMuon()) continue;

    // get the segments which match the muon candidate.  Comes from class defined below
    SegmentVector    *segVect    = SegmentsInMuon( &(*muon), &(*cscSegments) );
    dt_SegmentVector *dt_segVect = dt_SegmentsInMuon( &(*muon), &(*dtSegments) );
    
    if (printLevel > 1){
      cout << "The muon has " << segVect -> size()    << " csc segments"     << endl;
      cout << "The muon has " << dt_segVect -> size() << " dt segments"     << endl;
    }
    
    muonNsegs->push_back( segVect->size() + dt_segVect->size() );
    if (printLevel > 1) cout << "The total segs is " <<  segVect->size() + dt_segVect->size() << endl;
    
    // sanity check
    if (segVect -> size() == 0 && dt_segVect == 0){
      delete segVect;
      delete dt_segVect;
      continue;
    }
    
    int iSegment=-1;
    // loop over the CSC segments
    if (printLevel > 1) cout << "Looping over the CSC segments\n\n";
    for(auto segmentCSC = segVect -> begin(); segmentCSC != segVect->end(); segmentCSC++){
      iSegment++;
      
      if (printLevel > 1) {
        printf("#############\n");
        printf("Segment  #%d\n", iSegment);
        printf("#############\n\n");
      
	
	if ( iSegment > (MAX_SEGS_STD-1) ) {
	  cout << "the muon has " << iSegment+1 << ", but the MAX allowed is "
	       << MAX_SEGS_STD << " -> Skipping the segment... " << endl;
	  continue;
	}
      }
      // basic info
      CSCDetId id  = (CSCDetId) (*segmentCSC)->cscsegcand.cscDetId();
      const CSCChamber* cscchamber = cscGeom->chamber(id);
      
      if (!cscchamber) { 
        cout << "cscchamber not valid" << endl;
        continue;
      }
      
      muon_isdtseg       [whichMuon][iSegment] = 0;
      muon_cscsegs_endcap [whichMuon][iSegment] = id.endcap(); 
      muon_cscsegs_station[whichMuon][iSegment] = id.station();
      muon_cscsegs_ring   [whichMuon][iSegment] = id.ring();   
      muon_cscsegs_chamber[whichMuon][iSegment] = id.chamber();
      muon_cscsegs_nhits  [whichMuon][iSegment] = (*segmentCSC)->cscsegcand.nRecHits();

      if (printLevel > 2) {
        std::cout << "Endcap  = " << muon_cscsegs_endcap [whichMuon][iSegment] << std::endl;
        std::cout << "Ring    = " << muon_cscsegs_station[whichMuon][iSegment] << std::endl;
        std::cout << "Station = " << muon_cscsegs_ring   [whichMuon][iSegment] << std::endl;
        std::cout << "Chamber = " << muon_cscsegs_chamber[whichMuon][iSegment] << std::endl;
        std::cout << "nRecHits= " << muon_cscsegs_nhits  [whichMuon][iSegment] << std::endl << std::endl;
      }
      
      // local segment position
      LocalPoint localPos = (*segmentCSC)->cscsegcand.localPosition();
      LocalVector  segDir = (*segmentCSC)->cscsegcand.localDirection();    
      
      muon_cscsegs_loc_x      [whichMuon][iSegment] = localPos.x();
      muon_cscsegs_loc_y      [whichMuon][iSegment] = localPos.y();
      muon_cscsegs_loc_eta    [whichMuon][iSegment] = localPos.eta();
      muon_cscsegs_loc_phi    [whichMuon][iSegment] = localPos.phi();
      muon_cscsegs_loc_dir_eta[whichMuon][iSegment] = segDir.eta();
      muon_cscsegs_loc_dir_phi[whichMuon][iSegment] = segDir.phi();

      if (printLevel > 2) {
        std::cout << "localPos.x()   = " << muon_cscsegs_loc_x      [whichMuon][iSegment] << " cm"<< std::endl; 
        std::cout << "localPos.y()   = " << muon_cscsegs_loc_y      [whichMuon][iSegment] << " cm"<< std::endl; 
        std::cout << "localPos.eta() = " << muon_cscsegs_loc_eta    [whichMuon][iSegment] << std::endl; 
        std::cout << "localPos.phi() = " << muon_cscsegs_loc_phi    [whichMuon][iSegment] << std::endl; 
        std::cout << "segDir.eta()   = " << muon_cscsegs_loc_dir_eta[whichMuon][iSegment] << std::endl; 
        std::cout << "segDir.phi()   = " << muon_cscsegs_loc_dir_phi[whichMuon][iSegment] << std::endl << std::endl; 
      }
      
     
      // global segment position
      GlobalPoint globalPosition   = cscchamber->toGlobal(localPos);
      GlobalVector globalDirection = cscchamber->toGlobal(segDir);      


      muon_cscsegs_gbl_x      [whichMuon][iSegment] = globalPosition.x();
      muon_cscsegs_gbl_y      [whichMuon][iSegment] = globalPosition.y();
      muon_cscsegs_gbl_eta    [whichMuon][iSegment] = globalPosition.eta();
      muon_cscsegs_gbl_phi    [whichMuon][iSegment] = globalPosition.phi();
      muon_cscsegs_gbl_dir_eta[whichMuon][iSegment] = globalDirection.eta();
      muon_cscsegs_gbl_dir_phi[whichMuon][iSegment] = globalDirection.phi();

      
      if (printLevel > 2) {
        std::cout << "globalPosition.x()    = " << muon_cscsegs_gbl_x      [whichMuon][iSegment] << " cm"<< std::endl; 
        std::cout << "globalPosition.y()    = " << muon_cscsegs_gbl_y      [whichMuon][iSegment] << " cm"<< std::endl; 
        std::cout << "globalPosition.eta()  = " << muon_cscsegs_gbl_eta    [whichMuon][iSegment] << std::endl; 
        std::cout << "globalPosition.phi()  = " << muon_cscsegs_gbl_phi    [whichMuon][iSegment] << std::endl; 
        std::cout << "globalDirection.eta() = " << muon_cscsegs_gbl_dir_eta[whichMuon][iSegment] << std::endl; 
        std::cout << "globalDirection.phi() = " << muon_cscsegs_gbl_dir_phi[whichMuon][iSegment] << std::endl << std::endl; 
      }

      
      // not sure we need them for anything: probably they should be deprecated...
      muon_cscsegs_dxdz   [whichMuon][iSegment] = -999;//(*segmentCSC)->cscsegcand.dXdZ;     
      muon_cscsegs_dydz   [whichMuon][iSegment] = -999;//(*segmentCSC)->cscsegcand.dYdZ;     
      muon_cscsegs_dxdzErr[whichMuon][iSegment] = -999;//(*segmentCSC)->cscsegcand.dXdZErr;  
      muon_cscsegs_dydzErr[whichMuon][iSegment] = -999;//(*segmentCSC)->cscsegcand.dYdZErr;  

      // is the segment triggerable?
      bool isTriggerAble = _matchBox.isLCTAble( (*segmentCSC)->cscsegcand, 0);
 
      // is the segment matched to an LCT?
      bool isLCTMatched  = _matchBox.isMatched( (*segmentCSC)->cscsegcand, CSCTFlcts , 0 );
      
      if (printLevel > 1) cout <<"isMatched?=" << isLCTMatched  << endl; 
      
      muon_cscsegs_islctable[whichMuon][iSegment] = isTriggerAble;
      muon_cscsegs_ismatched[whichMuon][iSegment] = isLCTMatched;
      
      vector<int> whichLCT;  // find the corresponding LCT in the list.  This seems to be done by matching CSCDetId for segment lct and list lct
      if (isLCTMatched) {
	int iLCT=-1;
	
        CSCDetId *segDetId = 0;
        const CSCDetId &origId = (*segmentCSC)->cscsegcand.cscDetId();
	
        // if we're in ME11a, we have to worry about triple-ganging of strips.
        if (origId.ring() == 4){
          segDetId = new CSCDetId ( origId.endcap(), origId.station(), 1,origId.chamber());
        } else {
          segDetId = new CSCDetId ( origId );
        }
	
	// loop over CSC Lcts 
	int match_count = 0; // keep track of how many lcts can be matched to the same segment
	for( auto corrLct = CSCTFlcts->cbegin(); corrLct != CSCTFlcts->cend(); corrLct++) {
          iLCT++;
	  if (printLevel > 2) cout << "Looping over Lct # " << iLCT << endl;
          
	  // get CSCDetId for Lcts
	  //cout << "setting LctId " << endl;
	  if ( corrLct->subsystem() != 1 ) continue;
	  CSCDetId LctId = corrLct->detId<CSCDetId>();
	  
	  // to show to Ivan: he was right, as always
          // segments have cscDetId == 4 for ME1/1a
          // while LCTs do not have it!!!
          //std::cout << "  lctRange1: " <<   (*segmentCSC)->cscsegcand.cscDetId() << endl;   	
          //std::cout << "  lctRange2: " <<   (*corrLct).first << endl; 
	  
	  // Compare type of segDetId and LctId
	  //	  cout << "LctId type is " << typeid(LctId).name() << endl;
	  //cout << "SegDetId type is " << typeid(*segDetId).name() << endl;
          
          // find the matching one
	  if ( (*segDetId) == LctId) {
	    match_count ++;
            whichLCT.push_back(iLCT);
	    if (printLevel > 1 ) cout << "Match is found. Corresponds to LCT number:" << iLCT << endl;
	  }

	}  // end loop CSC lcts
	
	// If match has been made to only one or no lcts do the following
	if (match_count == 1) {
	  muon_cscsegs_lctId[whichMuon][iSegment] = whichLCT.front();
	  if (printLevel > 1) cout << "Fill muon_cscsegs_lctId with " << whichLCT.front()  << endl;
	}
	
	if (match_count == 0) muon_cscsegs_lctId[whichMuon][iSegment] = -999;
	
	// --------------- If two lcts are matched to segment choose the one that is in CSCTF list, if at all. ---------------------
	if (match_count > 1) {
	  if (printLevel > 1) cout << "/n------- More than one lct was matched to segment -------" << endl;
	  // First loop over all Lcts in event
	  int iLCT = 0;
	  bool isMatch = 0;
	  for( auto Lct = CSCTFlcts->cbegin(); Lct != CSCTFlcts->cend(); Lct++, iLCT++) {
	    if (isMatch) continue;
	    if (printLevel > 2) cout << "Looping over event LCT # " << iLCT << endl;
	
	    // only look at Lcts that have been matched to a segment
	    for ( auto lct = whichLCT.cbegin(); lct < whichLCT.cend(); lct++) {
	      if (isMatch) continue;
	      if (printLevel > 2) cout << "LCT that was matched to seg is " << *lct << ". Comparing to event LCT # " << iLCT << endl;
	      if ( iLCT != *lct) continue;
	      if (printLevel > 2) cout << "Match was found" << endl;
	      
	      // access phi/eta values
	      if ( Lct->subsystem() != 1 ) continue;
	      CSCDetId id = Lct->detId<CSCDetId>();
	      auto lct_station           = id.station();
	      auto lct_endcap            = id.endcap();
	      //int lct_sector             = CSCTriggerNumbering::triggerSectorFromLabels(id);
	      int lct_subsector          = CSCTriggerNumbering::triggerSubSectorFromLabels(id);
	      uint16_t lct_strip         = Lct->getCSCData().strip;
	      uint16_t lct_pattern       = Lct->getCSCData().pattern;
	      uint16_t lct_bend          = Lct->getCSCData().bend;
	      uint16_t lct_quality       = Lct->getCSCData().quality;
	      uint16_t lct_keywire       = Lct->getCSCData().keywire;
	      uint16_t lct_cscID         = Lct->getCSCData().cscID;
	      int FPGALct = ( lct_subsector ? lct_subsector-1 : lct_station );
	      if (lct_endcap == 2) lct_endcap = 0;

	      // local Phi
	      lclphidat lclPhi;
	      lclPhi = srLUTs_[FPGALct][lct_endcap] -> localPhi(lct_strip,
								lct_pattern,
								lct_quality,
								lct_bend);
	      
	      // Fill phi bit
	      gblphidat gblPhi;
	      gblPhi = srLUTs_[FPGALct][lct_endcap] -> globalPhiME(lclPhi.phi_local,
								   lct_keywire,
								   lct_cscID);
	      
	      auto evt_lctPhi = gblPhi.global_phi;
	      
	      // Fill eta bit
	      gbletadat gblEta;
	      gblEta = srLUTs_[FPGALct][lct_endcap] -> globalEtaME(lclPhi.phi_bend_local,
								   lclPhi.phi_local,
								   lct_keywire,
								   lct_cscID);
	      auto evt_lctEta = gblEta.global_eta;
	     
	      // Loop over track lcts. Access phi/eta and compare to matched event lct phi/eta.
	      int nTrk=0;
	      int itrkLCT = 0;	
	      for( auto trk = CSCTFtracks->cbegin(); trk < CSCTFtracks->cend(); trk++){
		if (isMatch) continue;	
		if (printLevel > 2) cout << "Looping over track # " << nTrk << endl;
		auto lct_map = trk -> getStubs();

		for( unsigned station = 1; station <= 4; ++station ) {
		  const unsigned id = 4*L1TMuon::InternalTrack::kCSC+station-1; // unique identifier for each csc station
		  auto x_LCTs = lct_map[id];  // access the matched lcts for a given station id
		  
		  for ( auto t_lcts = x_LCTs.cbegin(); t_lcts != x_LCTs.cend(); t_lcts++, itrkLCT++ ) {
		    auto trkLct = *t_lcts; // dereference the edm:Ref object
		    if (isMatch) continue;
		    if (printLevel > 2) cout << "Looping over track lct # " << itrkLCT << endl;
		  
		  // access track lct phi/eta
		  if ( trkLct->subsystem() != 1 ) continue;
		  CSCDetId id = trkLct->detId<CSCDetId>();
		  auto trklct_station           = id.station();
		  auto trklct_endcap            = id.endcap();
		  //int trklct_sector             = CSCTriggerNumbering::triggerSectorFromLabels(id);
		  int trklct_subsector          = CSCTriggerNumbering::triggerSubSectorFromLabels(id);
		  uint16_t trklct_strip         = trkLct->getCSCData().strip;
		  uint16_t trklct_pattern       = trkLct->getCSCData().pattern;
		  uint16_t trklct_bend          = trkLct->getCSCData().bend;
		  uint16_t trklct_quality       = trkLct->getCSCData().quality;
		  uint16_t trklct_keywire       = trkLct->getCSCData().keywire;
		  uint16_t trklct_cscID         = trkLct->getCSCData().cscID;
		  int FPGALct = ( trklct_subsector ? trklct_subsector-1 : trklct_station );
		  if (trklct_endcap == 2) trklct_endcap = 0;

		  // local Phi
		  lclphidat lclPhi;
		  lclPhi = srLUTs_[FPGALct][trklct_endcap] -> localPhi(trklct_strip,
								    trklct_pattern,
								    trklct_quality,
								    trklct_bend);

		  // Fill phi bit
		  gblphidat gblPhi;
		  gblPhi = srLUTs_[FPGALct][trklct_endcap] -> globalPhiME(lclPhi.phi_local,
								       trklct_keywire,
								       trklct_cscID);

		  auto trk_lctPhi = gblPhi.global_phi;

		  // Fill eta bit
		  gbletadat gblEta;
		  gblEta = srLUTs_[FPGALct][trklct_endcap] -> globalEtaME(lclPhi.phi_bend_local,
								       lclPhi.phi_local,
								       trklct_keywire,
								       trklct_cscID);
		  auto trk_lctEta = gblEta.global_eta;

		  if (printLevel > 2) { 
		    cout << " Event LCT phi = " << evt_lctPhi << " , Trk LCT phi = " << trk_lctPhi << endl;
		    cout << " Event LCT phi = " << evt_lctEta << " , Trk LCT eta = " << trk_lctEta << endl; 
		  }

		  // Compare trk lct eta/phi to event lct eta/phi
		  if ( evt_lctPhi == trk_lctPhi && evt_lctEta == trk_lctEta) {
		    if (printLevel > 1) cout << "Match is found for event LCT. " << "Fill muon_cscsegs_lctId with " << iLCT << endl;		    
		    isMatch = 1;
		    muon_cscsegs_lctId[whichMuon][iSegment] = iLCT;
		    continue;
		  }		  		  		  
		  
		  } // end loop over station track lcts
		} // end loop over track stations
	      } // end loop over track 	      	      
	    } // end loop over matched lcts
	  } // end loop over event lcts
	  // If no match was found in Lct list then just choose the first whichLCT element
	  if (!isMatch) {
	    if (printLevel > 1) cout << "No match from event list.  Choose first lct then. muon_cscsegs_lctId = " << whichLCT.front() << endl;
	    muon_cscsegs_lctId[whichMuon][iSegment] = whichLCT.front();
	  }

	} // end if match_count > 1	
      } // end isMatched
    } // end loop over segments
    delete segVect;

    //---------- Begin Loop over DT Segments ------------
      iSegment++;
      if (printLevel > 1) cout << "Looping over the DT segments\n\n";
    for(auto segmentDT = dt_segVect->begin(); segmentDT != dt_segVect->end(); segmentDT++, iSegment++) {
      int matchcount = 0;
      vector<int> whichLCT;
      
      if (printLevel > 1) {
	printf("#############\n");
	printf("Segment  #%d\n", iSegment);
	printf("#############\n\n");
		
	if ( iSegment > (MAX_SEGS_STD-1) ) {
	  cout << "the muon has " << iSegment+1 << ", but the MAX allowed is "
	       << MAX_SEGS_STD << " -> Skipping the segment... " << endl;
	  continue;
	}
      }
      // basic info
      //      cout << " dt seg iterator is type " << typeid(*segmentDT->dtsegcand).name() << endl;
      
      DTChamberId seg_id = (*segmentDT)->dtsegcand.chamberId();
      const DTChamber* dtchamber = dtGeom->chamber(seg_id);
      
      if (!dtchamber) {
	cout << "dtchamber not valid" << endl;
	continue;
      }
  
      LocalPoint localPos = (*segmentDT)->dtsegcand.localPosition();
      
      GlobalPoint globalPosition = dtchamber->toGlobal(localPos);
      //globX        = globalPosition.x();
      //globY        = globalPosition.y();
      //globZ        = globalPosition.z();
      //globSegPhi   = globalPosition.phi();
      //globSegEta   = globalPosition.eta();
      //globSegTheta = globalPosition.theta();

      muon_isdtseg       [whichMuon][iSegment] = 1;
      muon_cscsegs_gbl_eta [whichMuon][iSegment] = globalPosition.eta();
      muon_cscsegs_gbl_phi [whichMuon][iSegment] = globalPosition.phi();
      muon_cscsegs_station [whichMuon][iSegment] = seg_id.station();
      muon_cscsegs_chamber [whichMuon][iSegment] = seg_id;
      muon_cscsegs_wheel   [whichMuon][iSegment] = seg_id.wheel();
      muon_cscsegs_sector  [whichMuon][iSegment] = seg_id.sector();
      //muon_dtsegs_nhits  [whichMuon][iSegment] = (*segmentCSC)->cscsegcand.nRecHits();

      if (printLevel > 1) {
	cout << "muon_dtsegs_gbl_eta = "       << muon_cscsegs_gbl_eta [whichMuon][iSegment] << endl;
	cout << "muon_dtsegs_gbl_phi = "       << muon_cscsegs_gbl_phi [whichMuon][iSegment] << endl;
	cout << "muon_dtsegs_gbl_chamberId = " << muon_cscsegs_chamber [whichMuon][iSegment] << endl;
	cout << "muon_dtsegs_gbl_station = "   << muon_cscsegs_station [whichMuon][iSegment] << endl;
	cout << "muon_dtsegs_gbl_wheel = "     << muon_cscsegs_wheel [whichMuon][iSegment]   << endl;
	cout << "muon_dtsegs_gbl_sector = "    << muon_cscsegs_sector [whichMuon][iSegment]  << endl;
      }

      // is DT seg matched to an lct in this event. Skip matchbox. Use, chamberid, wheel, station, sector, phi to make match.
      // loop over DT TrigPrimitives
      int iLCT=0;
      for ( auto Lct = CSCTFlcts->cbegin(); Lct != CSCTFlcts->cend(); Lct++, iLCT++) {
	
	// ----- Fill all DT LCTs ------
	if ( Lct->subsystem() != 0 ) continue;
	if (printLevel > 1) cout << "Looping over Lct # " << iLCT << endl;
	  
	DTChamberId id = Lct -> detId<DTChamberId>();
	int dt_lct_station        = Lct->getDTData().station;
	double dt_lct_phi         = Lct->getCMSGlobalPhi();
	double dt_lct_eta         = Lct->getCMSGlobalEta();
	int dt_lct_sector         = Lct->getDTData().sector;
	int dt_lct_wheel          = Lct->getDTData().wheel;
	
	if (printLevel > 1) {
	  cout << "dt_lct_eta = "     << dt_lct_eta << endl;
	  cout << "dt_lct_phi = "     << dt_lct_phi << endl;
	  cout << "chamberId = "      << id          << endl;
	  cout << "dt_lct_station = " << dt_lct_station << endl;
	  cout << "dt_lct_wheel = "   << dt_lct_wheel << endl;
	  cout << "dt_lct_sector = "  << dt_lct_sector << endl;
	}

	// Find delta eta/phi		       
	float deta = muon_cscsegs_gbl_eta [whichMuon][iSegment] - dt_lct_eta;
	float dphi = muon_cscsegs_gbl_phi [whichMuon][iSegment] - dt_lct_phi;

	// The dphi cutoffs come from plotting and choosing a good distribution width
	if ( abs(deta) < 0.1    &&    abs(dphi) < 0.03   &&
	     muon_cscsegs_station [whichMuon][iSegment] == dt_lct_station &&
	     muon_cscsegs_wheel  [whichMuon][iSegment] == dt_lct_wheel 
	     //muon_dtsegs_sector [whichMuon][iSegment] == dt_lct_sector 
	     //muon_dtsegs_chamber [whichMuon][iSegment] == id
	     ) {
	 
	  muon_cscsegs_ismatched[whichMuon][iSegment] = 1;
	  whichLCT.push_back(iLCT);
	  matchcount++;
	  if (printLevel > 1) cout << "-----> DT Seg-Lct Match IS found! " << endl;
	}
	else {  if (printLevel > 1) cout << "-----> DT Seg-Lct Match is NOT found! " << endl; }
	 	
      } // end loop over DT lcts
      
      // Check to see if segment was matched to more than one lct. If so take special care
      if ( matchcount == 1 ) {
	muon_cscsegs_lctId[whichMuon][iSegment] = whichLCT.front();
	if (printLevel > 1) cout << "Fill muon_cscsegs_lctId with " << whichLCT.front() << endl;
      }
      if ( matchcount == 0 ) muon_cscsegs_lctId[whichMuon][iSegment] = -999;
      
      // More than one seg/lct match
      if (matchcount > 1) {
	if (printLevel >1 ) cout << "------- More than one Lct matched to DT segment --------" << endl;
	
	// First loop over all Lcts in event
	int iLCT = 0;
	bool isMatch = 0;
	for( auto Lct = CSCTFlcts->cbegin(); Lct != CSCTFlcts->cend(); Lct++, iLCT++) {
	  if ( Lct->subsystem() != 0 ) continue;
	  if (isMatch) continue;
	  if (printLevel > 1) cout << "Looping over event LCT # " << iLCT << endl;

	  // only look at Lcts that have been matched to a segment
	  for ( auto lct = whichLCT.cbegin(); lct < whichLCT.cend(); lct++) {
	    if (isMatch) continue;
	    if (printLevel > 1) cout << "LCT that was matched to seg is " << *lct << ". Comparing to event LCT # " << iLCT << endl;
	    if ( iLCT != *lct) continue;
	    if (printLevel > 1) cout << "Lct from list found" << endl;
	    
	    DTChamberId lct_id        = Lct->detId<DTChamberId>();
	    int dt_lct_station        = Lct->getDTData().station;
	    double dt_lct_phi         = Lct->getCMSGlobalPhi();
	    double dt_lct_eta         = Lct->getCMSGlobalEta();
	    unsigned dt_lct_sector    = Lct->getDTData().sector;
	    int dt_lct_wheel          = Lct->getDTData().wheel;

	    // Loop over track lcts. Access phi/eta and compare to matched event lct phi/eta.
	    int nTrk=0;
	    for( auto trk = CSCTFtracks->cbegin(); trk < CSCTFtracks->cend(); trk++){
	      if (isMatch) continue;
	      if (printLevel > 1) cout << "Looping over track # " << nTrk << endl;
	      auto lct_map = trk -> getStubs();
	      
	     
	      for( unsigned station = 1; station <= 4; ++station ) {
		const unsigned id = 4*L1TMuon::InternalTrack::kDT+station-1; // unique identifier for each dt station
		auto x_LCTs = lct_map[id];  // access the matched lcts for a given station id

		int itrkLCT = 0;
		for ( auto t_lcts = x_LCTs.cbegin(); t_lcts != x_LCTs.cend(); t_lcts++, itrkLCT++ ) {
		  auto trkLct = *t_lcts; // dereference the edm:Ref object
		  if (isMatch) continue;
		  if (printLevel > 1) cout << "Looping over track lct # " << itrkLCT << endl;

 		  DTChamberId trid            = trkLct->detId<DTChamberId>();
		  int dt_trlct_station        = trkLct->getDTData().station;
		  double dt_trlct_phi         = trkLct->getCMSGlobalPhi();
		  double dt_trlct_eta         = trkLct->getCMSGlobalEta();
		  //uint16_t dt_trlct_bx        = trkLct->getDTData().bx;
		  unsigned dt_trlct_sector    = trkLct->getDTData().sector;
		  int dt_trlct_wheel          = trkLct->getDTData().wheel;
		  
		  // compare track lct to event lct
		  if ( lct_id == trid &&
		       dt_lct_station == dt_trlct_station &&
 		       dt_lct_phi     == dt_trlct_phi     && 
		       dt_lct_eta     == dt_trlct_eta     &&
		       dt_lct_sector  == dt_trlct_sector  &&
		       dt_lct_wheel   == dt_trlct_wheel ) {
		    
		    if (printLevel > 1) {
		      cout << "Match is found between event LCT " << iLCT << " and track LCT " << itrkLCT << endl;
		      cout << "Fill muon_cscsegs_lctId with " << iLCT << endl;
		    }
		      isMatch = 1;
		      muon_cscsegs_lctId[whichMuon][iSegment] = iLCT;
		      continue;		    
		  } 
		       
		} // endl loop over DT track lcts
	      } // end loop over stations
	    } // end loop over tracks
	  } // end loop over matched lcts
	} // end loop over lcts
	// If no match was found in Lct list then just choose the first whichLCT element
	if (!isMatch) {
	  if (printLevel > 1) cout << "No match from event list.  Choose first lct then" << endl;
	  muon_cscsegs_lctId[whichMuon][iSegment] = whichLCT.front();
	}
	
      } // end matchcount > 1    
    } // end Loop over DT segemnts
  } // end loop over muons
  
}  // --------- end FillSegemntsMuon ----------------

TrigEff::SegmentVector* TrigEff::SegmentsInMuon(const reco::Muon* muon, 
                                                const CSCSegmentCollection* segments) {
  
  TrigEff::SegmentVector *retVal    = new TrigEff::SegmentVector();
  
  if (printLevel > 1) cout << "\nSegmentsInMuon Method\n";

  bool isMuonStd=false; // has the muon a standalone component
  if (muon->combinedMuon().isNonnull() || muon->standAloneMuon().isNonnull()) 
    isMuonStd=true;
 
  // return empty vector if the muon is not interesting
  if (!isMuonStd) return retVal;
  int nMuonMatchedHits=0;
  int iSegment=0;
  // --------- Loop over the CSC segments -------------
  for (auto segIter = segments->begin(); segIter != segments->end(); segIter++){

    int nHits=segIter -> nRecHits();

    if (printLevel > 2) {
      cout << " ======================================== " << endl;
      cout << "Segment in CSC:" << iSegment++ << endl;
      cout << "# segs hits:"    << nHits      << endl;
    }

    const std::vector<CSCRecHit2D>& theHits = segIter -> specificRecHits();
    std::vector<CSCRecHit2D>::const_iterator hitIter;
    
    int iHit=0;
    // loop over the segments hits
    for (hitIter = theHits.begin(); hitIter!= theHits.end(); hitIter++){

      if (printLevel>2) std::cout << "iHit:" << iHit++ << ", ";
      
      // check if the hit will match the standalone muon component
      bool isHitMatched=false;
      
      LocalPoint seghitlocal = hitIter -> localPosition();

      double segHitX = seghitlocal.x();      
      double segHitY = seghitlocal.y();      

      if (printLevel > 2)
        std::cout << "segHitX="<<segHitX <<  ", "
                  << "segHitY="<<segHitY;
        
      
      // The muon now returns segments (2012), while in 2010 it was returning hits...
      for(trackingRecHit_iterator segm  = muon->outerTrack()->recHitsBegin(); 
          segm != muon->outerTrack()->recHitsEnd(); 
          segm++){


        // some basic protection
	if ( !((*segm)->isValid()) ) continue;
        
	// Hardware ID of the RecHit (in terms of wire/strips/chambers)
	DetId detid = (*segm)->geographicalId();
        
	// Interested in muon systems only
	if( detid.det() != DetId::Muon ) continue;

        //Look only at CSC Hits (CSC id is 2)
        if (detid.subdetId() != MuonSubdetId::CSC) continue;
          
        CSCDetId id(detid.rawId());
          
        // another sanity check
        if  (id.station() < 1) continue;
        
        // get the CSCSegment
        const CSCSegment* cscSegment = dynamic_cast<const CSCSegment*>(&**segm);
        // check the segment is not NULL
        if (cscSegment == NULL) continue;

        // try to get the CSC recHits that contribute to this segment.
        std::vector<CSCRecHit2D> theseRecHits = (*cscSegment).specificRecHits();

        // loop over the rechits
        for ( std::vector<CSCRecHit2D>::const_iterator iRH = theseRecHits.begin(); 
              iRH != theseRecHits.end(); iRH++) {
          
          // get the rechit ID
          CSCDetId idRH = (CSCDetId)(*iRH).cscDetId();
   
          // CSC chamber
          const CSCChamber* cscchamber = cscGeom->chamber(idRH);
          if (!cscchamber) continue;

          // local position
          LocalPoint rhitlocal = iRH->localPosition();

          if (segHitX==rhitlocal.x() &&
              segHitY==rhitlocal.y()  )
            isHitMatched=true;
        } // end loop over hits of a segment

      } // end loop trackingRecHit_iterator (segments of a muon)

      if (printLevel > 2 ){
        if (isHitMatched) cout<< " -> Matched" << endl;
        else              cout<< " -> NOT Matched" << endl;
      }

      if (isHitMatched) nMuonMatchedHits++;

    }
    
    if (printLevel > 2) cout<< "segment has "  << nMuonMatchedHits
                               << " hits out of " << nHits 
                               << " matched"      << endl;
  

    // fill the the vector with the matching segments
    if (nMuonMatchedHits!=0) {
      Segment* segment = new Segment(*segIter, nMuonMatchedHits);
      retVal->push_back(segment);
    }
    
  }  // end loop on the CSC segments
  return retVal;

} // end function


TrigEff::dt_SegmentVector* TrigEff::dt_SegmentsInMuon(const reco::Muon* muon,                                                
						      const DTRecSegment4DCollection* dt_segments) {

  TrigEff::dt_SegmentVector *dt_retVal    = new TrigEff::dt_SegmentVector();

  if (printLevel > 1) cout << "dt_SegmentsInMuon Method\n\n";

  bool isMuonStd=false; // has the muon a standalone component
  if (muon->combinedMuon().isNonnull() || muon->standAloneMuon().isNonnull())
    isMuonStd=true;

  // return empty vector if the muon is not interesting
  if (!isMuonStd) return dt_retVal;
  int nMuonMatchedHits=0;
  int iSegment=0;

  // -------Loop over DT segments -------------
  for (auto segIter = dt_segments->begin(); segIter != dt_segments->end(); segIter++){

    int nHits = segIter->recHits().size();

    if (printLevel > 2) {
      cout << " ======================================== " << endl;
      cout << "Segment in DT:" << iSegment++ << endl;
      cout << "# segs hits:"    << nHits      << endl;
    }
    
    auto theHits = segIter->recHits();
    int iHit=0;
    // loop over the segments hits
    for (auto hitIter = theHits.begin(); hitIter!= theHits.end(); hitIter++, iHit++){
      if (printLevel > 2) cout << "iHit:" << iHit << endl;
   
      // check if the hit will match the standalone muon component
      bool isHitMatched=false;
      auto hit = *hitIter;     
      LocalPoint seghitlocal = hit->localPosition();
      double segHitX = seghitlocal.x();
      double segHitY = seghitlocal.y();
      
      if (printLevel > 2)
	cout << "segHitX=" << segHitX <<  ", "
	     << "segHitY=" << segHitY;
      
      // loop over muon segments 
      for( auto segm = muon->outerTrack()->recHitsBegin(); segm != muon->outerTrack()->recHitsEnd(); segm++) {
	
	// some basic protection
        if ( !((*segm)->isValid()) ) continue;

        // Hardware ID of the RecHit (in terms of wire/strips/chambers)
        DetId detid = (*segm)->geographicalId();

        // Interested in muon systems only
        if( detid.det() != DetId::Muon ) continue;

        //Look only at DT Hits (DT id is 0)
        if (detid.subdetId() != MuonSubdetId::DT) continue;

	//auto id(detid.rawId());

        // get the DT Segment
	const DTRecSegment4D* dtSegment = dynamic_cast<const DTRecSegment4D*>(&**segm);
	// check the segment is not NULL
        if (dtSegment == NULL) continue;
	
	// try to get the DT recHit vector that contribute to this segment.
	auto theseRecHits = (*dtSegment).recHits();
	
	// loop over the rechits
        for ( auto iRH = theseRecHits.begin(); iRH != theseRecHits.end(); iRH++) {
	  
          // get the rechit ID
          auto muon_hit = *iRH;
	  //DTChamberId idRH = muon_hit->rawId();
	  
	  // local position
          LocalPoint rhitlocal = muon_hit->localPosition();
	  
          if (segHitX==rhitlocal.x() &&
              segHitY==rhitlocal.y()  )
            isHitMatched=true;
	  
	} // end loop over muon seg recHits
      } // end loop over muon segments
      
      if (printLevel > 2 ){
	if (isHitMatched) cout << " -> Matched" << endl;
	else              cout << " -> NOT Matched" << endl;
      }
      
      if (isHitMatched) nMuonMatchedHits++;
      
    } // end loop over segment recHits
    
    if (printLevel > 2) cout<< "segment has "  << nMuonMatchedHits
			    << " hits out of " << nHits
			    << " matched"      << std::endl;
    
    // fill the the vector with the matching segments
    if (nMuonMatchedHits!=0) {
      dt_Segment* segment = new dt_Segment(*segIter, nMuonMatchedHits);
      dt_retVal->push_back(segment);
    }
    
  } // end loop over DT segments
  
  return dt_retVal;
}

