#ifndef __EFFIEVENT_H__
#define __EFFIEVENT_H__

#include <vector>
#include "include/effi_constants.h"
#include "TFile.h"
#include "TTree.h"



using namespace std;

class EffiEvent{

 public:
  
  EffiEvent();

  //----------------------------------------------------------------------------
  // Useful Methods
  //----------------------------------------------------------------------------

  void Initialize();
  
  void AttachToFile(TFile *file);

  long int GetEntries(){ return   recoMuons ->GetEntries(); }

  void GetEntry(long int iEvt);

  bool IsGlobalRecoMuon( int iReco );
  bool IsStandaloneMuon( int iReco );
  bool IsTrackerMuon   ( int iReco );
  bool IsStdTrkNotGlb  ( int iReco );

  TTree *csctfTTree;
  TTree *recoMuons;

  //----------------------------------------------------------------------------
  // Access the needed variables
  //----------------------------------------------------------------------------
  int    Run, Event, Bx, Lumi, muonSize;

  vector<int>*   isGlobalMuon;
  vector<int>*   isStandAloneMuon;
  vector<int>*   isTrackerMuon;
  vector<int>*   isTMLastStationAngTight;
  vector<int>*   isGlobalMuonPromptTight;

  vector<float>* ptReco ;
  vector<float>* etaReco;
  vector<float>* phiReco;
  vector<float>* gmrChi2Norm;
  vector<float>* gmrDz;
  vector<float>* gmrD0;
  
  vector<float>* stdEta;
  vector<float>* stdPt ;
  vector<float>* trkEta;
  vector<float>* trkPhi;
  vector<float>* trkPt;
  vector<int>*   trkValHits;
  vector<float>* trkChi2Norm;
  vector<float>* trkDz;
  vector<float>* trkD0;
  
  vector<float>* rchEta;
  vector<float>* rchPhi;
  vector<int>*   rchCSCtype;
  
  // tracker muon variables
  vector<int>*    trkNSegs;
  
  int   trkSegChamberId    [MAX_MUONS][MAX_TRK_SEGS]; 
  int   trkSegRing         [MAX_MUONS][MAX_TRK_SEGS];    
  int   trkSegStation      [MAX_MUONS][MAX_TRK_SEGS]; 
  int   trkSegEndcap       [MAX_MUONS][MAX_TRK_SEGS];  
  int   trkSegTriggerSector[MAX_MUONS][MAX_TRK_SEGS];
  int   trkSegTriggerCscId [MAX_MUONS][MAX_TRK_SEGS]; 
  float trkSegXfromMatch   [MAX_MUONS][MAX_TRK_SEGS];
  float trkSegYfromMatch   [MAX_MUONS][MAX_TRK_SEGS];
  float trkSegPhifromMatch [MAX_MUONS][MAX_TRK_SEGS];

  int trkSegIsArb          [MAX_MUONS][MAX_TRK_SEGS];
  float trkSegX            [MAX_MUONS][MAX_TRK_SEGS];
  float trkSegY            [MAX_MUONS][MAX_TRK_SEGS];
  float trkSegR            [MAX_MUONS][MAX_TRK_SEGS];
  float trkSegPhi          [MAX_MUONS][MAX_TRK_SEGS];
  float trkSegEta          [MAX_MUONS][MAX_TRK_SEGS];

  std::vector<float>*    muons_x_me11;
  std::vector<float>*    muons_y_me11;
  std::vector<float>*    muons_z_me11;
  std::vector<float>*  muons_phi_me11;
  std::vector<float>*  muons_eta_me11;

  std::vector<float>*    muons_x_me1;
  std::vector<float>*    muons_y_me1;
  std::vector<float>*    muons_z_me1;
  std::vector<float>*  muons_phi_me1;
  std::vector<float>*  muons_eta_me1;
    
  //--------------------------------------------------------------------------
  // Record information about segments belonging to the STD muon component
  //--------------------------------------------------------------------------
  // how many segments are associated to the muon candidate
  int segsSize;
  std::vector<int>* muonNsegs; 
  
  // segment position information, local
  float muon_cscsegs_loc_x       [MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_loc_y       [MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_loc_eta     [MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_loc_phi     [MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_loc_dir_eta [MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_loc_dir_phi [MAX_MUONS][MAX_SEGS_STD];

  // segment position information, global
  float muon_cscsegs_gbl_x       [MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_gbl_y       [MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_gbl_eta     [MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_gbl_phi     [MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_gbl_dir_eta [MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_gbl_dir_phi [MAX_MUONS][MAX_SEGS_STD];

  // more on segment direction
  float muon_cscsegs_dxdz        [MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_dydz        [MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_dxdzErr     [MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_dydzErr     [MAX_MUONS][MAX_SEGS_STD];

  // general segment information
  int muon_cscsegs_endcap        [MAX_MUONS][MAX_SEGS_STD];
  int muon_cscsegs_station       [MAX_MUONS][MAX_SEGS_STD];
  int muon_cscsegs_ring          [MAX_MUONS][MAX_SEGS_STD];
  int muon_cscsegs_chamber       [MAX_MUONS][MAX_SEGS_STD];
  int muon_cscsegs_nhits         [MAX_MUONS][MAX_SEGS_STD];

  // isLCTAble
  int muon_cscsegs_islctable     [MAX_MUONS][MAX_SEGS_STD];
  int muon_cscsegs_ismatched     [MAX_MUONS][MAX_SEGS_STD];

  // lctId is the position of the lct in the all LCT collection
  // look also in fillAllLCTs
  int muon_cscsegs_lctId         [MAX_MUONS][MAX_SEGS_STD];
  int muon_isdtseg               [MAX_MUONS][MAX_SEGS_STD];

  // number of hits belonging to the STD fit
  int muon_cscsegs_nmatched      [MAX_MUONS][MAX_SEGS_STD];

  vector<float>* cscsegs_gbl_phi;
  vector<float>* cscsegs_gbl_eta;
  vector<int>*   cscsegs_chamber;
  vector<int>*   cscsegs_station;
  vector<int>*   cscsegs_wheel;
  vector<int>*   cscsegs_sector;

  // DT info
  vector<int>* isdtseg;

  /*vector<int>* dtsegs_chamber;
  vector<int>* dtsegs_station;
  vector<int>* dtsegs_sector ;
  vector<float>* dtsegs_gbl_x;
  vector<float>* dtsegs_gbl_y;
  vector<float>* dtsegs_gbl_z;

  vector<float>* dtsegs_gbl_theta;
  vector<float>* dtsegs_gbl_eta  ;
  vector<float>* dtsegs_gbl_phi  ;
  */

  // --- start new insertion of variables -- may be duplicates? //

  int SizeTrk;
  std::vector<int>*    EndcapTrk;  
  std::vector<int>*    SectorTrk; 
  std::vector<int>*    BxTrk    ;  

  std::vector<int>*    me1ID; 
  std::vector<int>*    me2ID; 
  std::vector<int>*    me3ID; 
  std::vector<int>*    me4ID; 
  std::vector<int>*    mb1ID;     

  std::vector<int>*    OutputLinkTrk;  

  std::vector<int>*    ModeTrk;  
  std::vector<float>*   EtaTrk;
  std::vector<float>*   PhiTrk;
  std::vector<float>*    PtTrk;

  std::vector<int>*    ChargeTrk;  
  std::vector<int>*    ChargeValidTrk;  
  std::vector<int>*    QualityTrk;  
  std::vector<int>*    ForRTrk;  
  std::vector<int>*    Phi23Trk;  

  std::vector<int>*    Phi12Trk;    
  std::vector<int>*    PhiSignTrk;    

  std::vector<int>*    EtaBitTrk;    
  std::vector<int>*    PhiBitTrk;    
  std::vector<int>*    PtBitTrk;    
 
  std::vector<int>* NumLCTsTrk;  

  int isdttrlct[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctEndcap[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctSector[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctSubSector[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctBx[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctBx0[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctStation[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctRing[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctChamber[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctTriggerCSCID[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctFpga[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];         
  // note: the SPs return them in bits 
  int trLctlocalPhi[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  //int trLctlocalEta[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctglobalPhi[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];   
  int trLctglobalEta[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctstripNum[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];   
  int trLctwireGroup[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];   

  int dt_trLctSector[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int dt_trLctStation[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int dt_trLctBx[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int dt_trLctChamber[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];

 // --------------------------------------------------------------------------- 
  int SizeLCTs;    
  vector<int>* lctEndcap; 
  vector<int>* lctSector; 
  vector<int>* lctSubSector; 
  vector<int>* lctBx; 
  vector<int>* lctBx0; 
  vector<int>* lctStation; 
  vector<int>* lctRing; 
  vector<int>* lctChamber; 
  vector<int>* lctTriggerCSCID; 
  vector<int>* lctFpga;     
  vector<int>* lctWheel;
  vector<int>* lctTheta_bti;
  
 // note: the SPs return them in bits 
  vector<int>* lctlocalPhi; 
  vector<int>* lctglobalPhi;   
  vector<int>* lctglobalEta; 
  vector<int>* lctstripNum;   
  vector<int>* lctwireGroup;   
  vector<double>* trig_prim_eta;
  vector<double>* trig_prim_phi;
  
  vector<int>* isdtlct      ;
  vector<int>* dt_lctSector ;
  vector<int>* dt_lctStation;
  vector<int>* dt_lctBx     ;
  vector<int>* dt_lctChamber;
  vector<float>* dt_lctglobalPhi;
  vector<float>* dt_lctglobalEta;



};



#endif