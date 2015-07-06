#include "include/EffiEvent.h"

EffiEvent::EffiEvent():
  isGlobalMuon(0),isStandAloneMuon(0), isTrackerMuon(0), isTMLastStationAngTight(0), isGlobalMuonPromptTight(0),
  ptReco (0), etaReco(0), phiReco(0), gmrChi2Norm(0), gmrDz(0), gmrD0(0),
  stdEta(0), stdPt (0), trkEta(0), trkPhi(0), trkPt(0), trkValHits(0), trkChi2Norm(0), trkDz(0), trkD0(0),
  rchEta(0), rchPhi(0), rchCSCtype(0), trkNSegs(0),
  muons_x_me11(0),muons_y_me11(0),muons_z_me11(0),muons_phi_me11(0),muons_eta_me11(0),
  muons_x_me1(0),muons_y_me1(0),muons_z_me1(0),muons_phi_me1(0),muons_eta_me1(0),
  muonNsegs(0),
  cscsegs_gbl_phi(0),cscsegs_gbl_eta(0),cscsegs_chamber(0),cscsegs_station(0),cscsegs_wheel(0),cscsegs_sector(0),
  isdtseg(0),
  SizeTrk(0),EndcapTrk(0),SectorTrk(0),BxTrk(0), me1ID(0), me2ID(0), me3ID(0), me4ID(0), mb1ID(0),     
  OutputLinkTrk(0), ModeTrk(0), EtaTrk(0), PhiTrk(0), PtTrk(0),ChargeTrk(0), ChargeValidTrk(0), 
  QualityTrk(0), ForRTrk(0), Phi23Trk(0), Phi12Trk(0), PhiSignTrk(0),EtaBitTrk(0), PhiBitTrk(0), PtBitTrk(0),
  NumLCTsTrk(0), lctEndcap(0), lctSector(0), lctSubSector(0), lctBx(0), lctBx0(0), lctStation(0), lctRing(0), 
  lctChamber(0), lctTriggerCSCID(0), lctFpga(0), lctWheel(0), lctTheta_bti(0),
  lctlocalPhi(0), lctglobalPhi(0), lctglobalEta(0), lctstripNum(0), lctwireGroup(0), trig_prim_eta(0), 
  trig_prim_phi(0), isdtlct(0), dt_lctSector(0), dt_lctStation(0), dt_lctBx(0), dt_lctChamber(0), 
  dt_lctglobalPhi(0), dt_lctglobalEta(0)
{



}

void EffiEvent::Initialize(){}


void EffiEvent::AttachToFile(TFile *file){

  if (!file) return;

  // get the Nutple  

  recoMuons  = (TTree*) file -> Get("recoMuons");
  csctfTTree = (TTree*) file -> Get("csctfTTree");

  //--------------------------------------------------------------------------

  recoMuons->SetBranchAddress("Run"  , &Run  );
  recoMuons->SetBranchAddress("Event", &Event);
  recoMuons->SetBranchAddress("Bx"   , &Bx   );
  recoMuons->SetBranchAddress("Lumi" , &Lumi );

  recoMuons->SetBranchAddress("muonSize"           , &muonSize        );
  recoMuons->SetBranchAddress("isGlobalMuon"       , &isGlobalMuon    );
  recoMuons->SetBranchAddress("isStandAloneMuon"   , &isStandAloneMuon);
  recoMuons->SetBranchAddress("isTrackerMuon"      , &isTrackerMuon   );
  recoMuons->SetBranchAddress("isTMLastStationAngTight", &isTMLastStationAngTight   );
  recoMuons->SetBranchAddress("isGlobalMuonPromptTight", &isGlobalMuonPromptTight);
  
  recoMuons->SetBranchAddress("gmrPt" , &ptReco );
  recoMuons->SetBranchAddress("gmrEta", &etaReco);
  recoMuons->SetBranchAddress("gmrPhi", &phiReco);
  recoMuons->SetBranchAddress("gmrChi2Norm", &gmrChi2Norm);
  recoMuons->SetBranchAddress("gmrD0",  &gmrD0);
  recoMuons->SetBranchAddress("gmrDz",  &gmrDz);
  
  recoMuons->SetBranchAddress("stdEta", &stdEta);
  recoMuons->SetBranchAddress("stdPt" , &stdPt);
  recoMuons->SetBranchAddress("trkEta", &trkEta);
  recoMuons->SetBranchAddress("trkPhi", &trkPhi);
  recoMuons->SetBranchAddress("trkPt",  &trkPt);
  recoMuons->SetBranchAddress("trkValHits",  &trkValHits);
  recoMuons->SetBranchAddress("trkChi2Norm", &trkChi2Norm);
  recoMuons->SetBranchAddress("trkD0", &trkD0);
  recoMuons->SetBranchAddress("trkDz", &trkDz);
  
  recoMuons->SetBranchAddress("rchEta", &rchEta);
  recoMuons->SetBranchAddress("rchPhi", &rchPhi);
  recoMuons->SetBranchAddress("rchCSCtype", &rchCSCtype);
  
  recoMuons->SetBranchAddress("trkNSegs"           , &trkNSegs          );
  recoMuons->SetBranchAddress("trkSegChamberId"    , trkSegChamberId    );
  recoMuons->SetBranchAddress("trkSegRing"         , trkSegRing         );
  recoMuons->SetBranchAddress("trkSegStation"      , trkSegStation      );
  recoMuons->SetBranchAddress("trkSegEndcap"       , trkSegEndcap       );
  recoMuons->SetBranchAddress("trkSegTriggerSector", trkSegTriggerSector);
  recoMuons->SetBranchAddress("trkSegTriggerCscId" , trkSegTriggerCscId );
  recoMuons->SetBranchAddress("trkSegXfromMatch"   , trkSegXfromMatch   );
  recoMuons->SetBranchAddress("trkSegYfromMatch"   , trkSegYfromMatch   );
  recoMuons->SetBranchAddress("trkSegPhifromMatch" , trkSegPhifromMatch );

  //segment
  recoMuons->SetBranchAddress("trkSegIsArb", trkSegIsArb);
  recoMuons->SetBranchAddress("trkSegX"    , trkSegX    );
  recoMuons->SetBranchAddress("trkSegY"    , trkSegY    );
  recoMuons->SetBranchAddress("trkSegR"    , trkSegR    );
  recoMuons->SetBranchAddress("trkSegPhi"  , trkSegPhi  );
  recoMuons->SetBranchAddress("trkSegEta"  , trkSegEta  );
   
  // propagation to ME1/1
  recoMuons->SetBranchAddress(  "muons_x_me11",  &muons_x_me11);
  recoMuons->SetBranchAddress(  "muons_y_me11",  &muons_y_me11);
  recoMuons->SetBranchAddress(  "muons_z_me11",  &muons_z_me11);
  recoMuons->SetBranchAddress("muons_phi_me11",&muons_phi_me11);
  recoMuons->SetBranchAddress("muons_eta_me11",&muons_eta_me11);
   
  recoMuons->SetBranchAddress(  "muons_x_me1",  &muons_x_me1);
  recoMuons->SetBranchAddress(  "muons_y_me1",  &muons_y_me1);
  recoMuons->SetBranchAddress(  "muons_z_me1",  &muons_z_me1);
  recoMuons->SetBranchAddress("muons_phi_me1",&muons_phi_me1);
  recoMuons->SetBranchAddress("muons_eta_me1",&muons_eta_me1);

  //---------------------------------------------------------------------
  // segments belonging to the muon
  //---------------------------------------------------------------------
  recoMuons->SetBranchAddress("muonNsegs",&muonNsegs);
  recoMuons->SetBranchAddress("segsSize",&segsSize);
  recoMuons->SetBranchAddress("muon_cscsegs_loc_x"      , muon_cscsegs_loc_x      );
  recoMuons->SetBranchAddress("muon_cscsegs_loc_y"      , muon_cscsegs_loc_y      );
  recoMuons->SetBranchAddress("muon_cscsegs_loc_eta"    , muon_cscsegs_loc_eta    );
  recoMuons->SetBranchAddress("muon_cscsegs_loc_phi"    , muon_cscsegs_loc_phi    );
  recoMuons->SetBranchAddress("muon_cscsegs_loc_dir_eta", muon_cscsegs_loc_dir_eta);
  recoMuons->SetBranchAddress("muon_cscsegs_loc_dir_phi", muon_cscsegs_loc_dir_phi);

  recoMuons->SetBranchAddress("muon_cscsegs_gbl_x"      , muon_cscsegs_gbl_x      );
  recoMuons->SetBranchAddress("muon_cscsegs_gbl_y"      , muon_cscsegs_gbl_y      );
  recoMuons->SetBranchAddress("muon_cscsegs_gbl_eta"    , muon_cscsegs_gbl_eta    );
  recoMuons->SetBranchAddress("muon_cscsegs_gbl_phi"    , muon_cscsegs_gbl_phi    );
  recoMuons->SetBranchAddress("muon_cscsegs_gbl_dir_eta", muon_cscsegs_gbl_dir_eta);
  recoMuons->SetBranchAddress("muon_cscsegs_gbl_dir_phi", muon_cscsegs_gbl_dir_phi);

  recoMuons->SetBranchAddress("muon_cscsegs_dxdz"   , muon_cscsegs_dxdz   );
  recoMuons->SetBranchAddress("muon_cscsegs_dydz"   , muon_cscsegs_dydz   );
  recoMuons->SetBranchAddress("muon_cscsegs_dxdzErr", muon_cscsegs_dxdzErr);
  recoMuons->SetBranchAddress("muon_cscsegs_dydzErr", muon_cscsegs_dydzErr);

  recoMuons->SetBranchAddress("muon_cscsegs_endcap" , muon_cscsegs_endcap );
  recoMuons->SetBranchAddress("muon_cscsegs_station", muon_cscsegs_station);
  recoMuons->SetBranchAddress("muon_cscsegs_ring"   , muon_cscsegs_ring   );
  recoMuons->SetBranchAddress("muon_cscsegs_chamber", muon_cscsegs_chamber);
  recoMuons->SetBranchAddress("muon_cscsegs_nhits"  , muon_cscsegs_nhits  );
  
  recoMuons->SetBranchAddress("muon_cscsegs_islctable", muon_cscsegs_islctable);
  recoMuons->SetBranchAddress("muon_cscsegs_ismatched", muon_cscsegs_ismatched);
  recoMuons->SetBranchAddress("muon_cscsegs_lctId"    , muon_cscsegs_lctId    );
  recoMuons->SetBranchAddress("muon_cscsegs_nmatched" , muon_cscsegs_nmatched );
  recoMuons->SetBranchAddress("muon_isdtseg"         , muon_isdtseg );

  recoMuons->SetBranchAddress("cscsegs_sector", &cscsegs_sector );
  recoMuons->SetBranchAddress("cscsegs_gbl_eta", &cscsegs_gbl_eta );
  recoMuons->SetBranchAddress("cscsegs_gbl_phi", &cscsegs_gbl_phi );
  recoMuons->SetBranchAddress("cscsegs_chamber", &cscsegs_chamber );
  recoMuons->SetBranchAddress("cscsegs_station", &cscsegs_station );
  recoMuons->SetBranchAddress("cscsegs_wheel"  , &cscsegs_wheel );

  recoMuons->SetBranchAddress("isdtseg", &isdtseg);

  csctfTTree->SetBranchAddress("SizeTrk"       , &SizeTrk       );
  csctfTTree->SetBranchAddress("EndcapTrk"     , &EndcapTrk     );
  csctfTTree->SetBranchAddress("SectorTrk"     , &SectorTrk     );
  csctfTTree->SetBranchAddress("BxTrk"         , &BxTrk         );
  csctfTTree->SetBranchAddress("me1ID"         , &me1ID         );
  csctfTTree->SetBranchAddress("me2ID"         , &me2ID         );
  csctfTTree->SetBranchAddress("me3ID"         , &me3ID         );
  csctfTTree->SetBranchAddress("me4ID"         , &me4ID         );
  csctfTTree->SetBranchAddress("mb1ID"         , &mb1ID         );
  csctfTTree->SetBranchAddress("OutputLinkTrk" , &OutputLinkTrk );
  csctfTTree->SetBranchAddress("ModeTrk"       , &ModeTrk       );
  csctfTTree->SetBranchAddress("EtaTrk"        , &EtaTrk        );
  csctfTTree->SetBranchAddress("PhiTrk_02PI"   , &PhiTrk        );
  csctfTTree->SetBranchAddress("PtTrk"         , &PtTrk         );
  csctfTTree->SetBranchAddress("ChargeTrk"     , &ChargeTrk     );
  csctfTTree->SetBranchAddress("ChargeValidTrk", &ChargeValidTrk);
  csctfTTree->SetBranchAddress("QualityTrk"    , &QualityTrk    );
  csctfTTree->SetBranchAddress("ForRTrk"       , &ForRTrk       );
  csctfTTree->SetBranchAddress("Phi23Trk"      , &Phi23Trk      );
  csctfTTree->SetBranchAddress("Phi12Trk"      , &Phi12Trk      );
  csctfTTree->SetBranchAddress("PhiSignTrk"    , &PhiSignTrk    );
  csctfTTree->SetBranchAddress("EtaBitTrk"     , &EtaBitTrk     );
  csctfTTree->SetBranchAddress("PhiBitTrk"     , &PhiBitTrk     );
  csctfTTree->SetBranchAddress("PtBitTrk"      , &PtBitTrk      );

  csctfTTree->SetBranchAddress("isdttrlct"        , isdttrlct         );
  csctfTTree->SetBranchAddress("NumLCTsTrk"       , &NumLCTsTrk       );
  csctfTTree->SetBranchAddress("trLctEndcap"      ,  trLctEndcap      );
  csctfTTree->SetBranchAddress("trLctSector"      ,  trLctSector      );
  csctfTTree->SetBranchAddress("trLctSubSector"   ,  trLctSubSector   );
  csctfTTree->SetBranchAddress("trLctBx"          ,  trLctBx          );
  csctfTTree->SetBranchAddress("trLctBx0"         ,  trLctBx0         );
  csctfTTree->SetBranchAddress("trLctStation"     ,  trLctStation     );
  csctfTTree->SetBranchAddress("trLctRing"        ,  trLctRing        );
  csctfTTree->SetBranchAddress("trLctChamber"     ,  trLctChamber     );
  csctfTTree->SetBranchAddress("trLctTriggerCSCID",  trLctTriggerCSCID);
  csctfTTree->SetBranchAddress("trLctFpga"        ,  trLctFpga        );
  //csctfTTree->SetBranchAddress("trLctlocalEta"    ,  trLctlocalEta    );
  csctfTTree->SetBranchAddress("trLctlocalPhi"    ,  trLctlocalPhi    );
  csctfTTree->SetBranchAddress("trLctglobalPhi"   ,  trLctglobalPhi   );
  csctfTTree->SetBranchAddress("trLctglobalEta"   ,  trLctglobalEta   );
  csctfTTree->SetBranchAddress("trLctstripNum"    ,  trLctstripNum    );
  csctfTTree->SetBranchAddress("trLctwireGroup"   ,  trLctwireGroup   );

  //csctfTTree->Branch("dt_NumLCTsTrk"        , &dt_NumLCTsTrk       );
  //csctfTTree->Branch("dt_trLctEndcap"      , dt_trLctEndcap      , "dt_trLctEndcap[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctSector"      , dt_trLctSector      , "dt_trLctSector[SizeTrk][4]");
  //csctfTTree->Branch("dt_trLctSubSector"   , dt_trLctSubSector   , "dt_trLctSubSector[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctBx"          , dt_trLctBx          , "dt_trLctBx[SizeTrk][4]");
  //csctfTTree->Branch("dt_trLctBx0"         , dt_trLctBx0         , "dt_trLctBx0[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctStation"     , dt_trLctStation     , "dt_trLctStation[SizeTrk][4]");
  //csctfTTree->Branch("dt_trLctRing"        , dt_trLctRing        , "dt_trLctRing[SizeTrk][4]");
  csctfTTree->Branch("dt_trLctChamber"     , dt_trLctChamber     , "dt_trLctChamber[SizeTrk][4]");
  //csctfTTree->Branch("dt_trLctTriggerCSCID", dt_trLctTriggerCSCID, "dt_trLctTriggerCSCID[SizeTrk][4]");
  //csctfTTree->Branch("dt_trLctFpga"        , dt_trLctFpga        , "dt_trLctFpga[SizeTrk][4]");
  //csctfTTree->Branch("dt_trLctlocalPhi"    , dt_trLctlocalPhi    , "dt_trLctlocalPhi[SizeTrk][4]");
  //csctfTTree->Branch("dt_trLctglobalPhi"   , dt_trLctglobalPhi   , "dt_trLctglobalPhi[SizeTrk][4]");
  //csctfTTree->Branch("dt_trLctglobalEta"   , dt_trLctglobalEta   , "dt_trLctglobalEta[SizeTrk][4]");
  //csctfTTree->Branch("dt_trLctstripNum"    , dt_trLctstripNum    , "dt_trLctstripNum[SizeTrk][4]");
  //csctfTTree->Branch("dt_trLctwireGroup"   , dt_trLctwireGroup   , "dt_trLctwireGroup[SizeTrk][4]");

  csctfTTree->SetBranchAddress("SizeLCTs"       , &SizeLCTs       );
  csctfTTree->SetBranchAddress("lctEndcap"      , &lctEndcap      );
  csctfTTree->SetBranchAddress("lctSector"      , &lctSector      );
  csctfTTree->SetBranchAddress("lctSubSector"   , &lctSubSector   );
  csctfTTree->SetBranchAddress("lctBx"          , &lctBx          );
  csctfTTree->SetBranchAddress("lctBx0"         , &lctBx0         );
  csctfTTree->SetBranchAddress("lctStation"     , &lctStation     );
  csctfTTree->SetBranchAddress("lctRing"        , &lctRing        );
  csctfTTree->SetBranchAddress("lctChamber"     , &lctChamber     );
  csctfTTree->SetBranchAddress("lctWheel"       , &lctWheel     );  
  csctfTTree->SetBranchAddress("lctTheta_bti"   , &lctTheta_bti     );
  csctfTTree->SetBranchAddress("lctTriggerCSCID", &lctTriggerCSCID);
  csctfTTree->SetBranchAddress("lctFpga"        , &lctFpga        );
  csctfTTree->SetBranchAddress("lctlocalPhi"    , &lctlocalPhi    );
  csctfTTree->SetBranchAddress("lctglobalPhi"   , &lctglobalPhi   );
  csctfTTree->SetBranchAddress("lctglobalEta"   , &lctglobalEta   );
  csctfTTree->SetBranchAddress("lctstripNum"    , &lctstripNum    );
  csctfTTree->SetBranchAddress("lctwireGroup"   , &lctwireGroup   );
  csctfTTree->SetBranchAddress("trig_prim_phi"   , &trig_prim_phi );
  csctfTTree->SetBranchAddress("trig_prim_eta"   , &trig_prim_eta );

  csctfTTree->SetBranchAddress("isdtlct"           , &isdtlct           );
  csctfTTree->SetBranchAddress("dt_lctSector"      , &dt_lctSector      );
  csctfTTree->SetBranchAddress("dt_lctBx"          , &dt_lctBx          );
  csctfTTree->SetBranchAddress("dt_lctStation"     , &dt_lctStation     );
  csctfTTree->SetBranchAddress("dt_lctChamber"     , &dt_lctChamber     );
  csctfTTree->SetBranchAddress("dt_lctglobalPhi"   , &dt_lctglobalPhi   );
  csctfTTree->SetBranchAddress("dt_lctglobalEta"   , &dt_lctglobalEta   );



}

void EffiEvent::GetEntry( long int iEvt ){

  recoMuons ->GetEntry(iEvt);
  csctfTTree->GetEntry(iEvt);
  
}

bool EffiEvent::IsGlobalRecoMuon( int iReco ){

  if (!isGlobalMuon->at(iReco))            return false;
  if (!isGlobalMuonPromptTight->at(iReco)) return false; 
  if (!isStandAloneMuon->at(iReco))        return false; 
  if (!isTrackerMuon->at(iReco))           return false; 
  if (gmrChi2Norm->at(iReco) > 10)         return false; 
  if (fabs(gmrD0->at(iReco)) >  2)         return false; 
  if ((rchEta->at(iReco)==-999))           return false; 

  return true;

}

bool EffiEvent::IsStandaloneMuon( int iReco ){

  if (isGlobalMuon->at(iReco))      return false;
  if (isTrackerMuon->at(iReco))     return false;
  if (!isStandAloneMuon->at(iReco)) return false;
  if ((rchEta->at(iReco)==-999))    return false;

  return true;

}

bool EffiEvent::IsTrackerMuon ( int iReco ){

  if (isGlobalMuon->at(iReco))             return false;
  if (isStandAloneMuon->at(iReco))         return false;
  if (muons_eta_me11->at(iReco) == -999)   return false;
  if (fabs(trkEta->at(iReco)) < 0.9)       return false;
  if (trkNSegs->at(iReco) < 1)             return false;
  // quarkonia tracker muons selections   
  //if (trkValHits->at(iReco) < 12)        return false;
  if (!isTMLastStationAngTight->at(iReco)) return false;
  if (trkChi2Norm->at(iReco) > 5.0)        return false;
  if (fabs(trkD0->at(iReco)) > 5.0)        return false;
  if (fabs(trkDz->at(iReco)) > 20)         return false;
  
  return true;

}

bool EffiEvent::IsStdTrkNotGlb( int iReco ){

  if (isGlobalMuon->at(iReco))       return false;
  if (!isStandAloneMuon->at(iReco))  return false;
  if (!isTrackerMuon->at(iReco))     return false;
  if ((rchEta->at(iReco)==-999))     return false;
  
  return true;

}
