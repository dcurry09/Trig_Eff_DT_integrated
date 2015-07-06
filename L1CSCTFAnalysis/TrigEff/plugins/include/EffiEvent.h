#ifndef __EFFIEVENT_H__
#define __EFFIEVENT_H__

#include <vector>
#include "TFile.h"
#include "TTree.h"

using namespace std;

class EffiEvent{

 public:
  
  EffiEvent();

  //----------------------------------------------------------------------------
  // Useful Methods
  //----------------------------------------------------------------------------

  
  void AttachToFile(TFile *file);
  
  void CloseFile(TFile *file);

  long int GetEntries() { return CSCtree ->GetEntries(); }

  void GetEntry(long int iEvt);

  TTree *CSCtree;
  TTree *RPCtree; 
  
  //----------------------------------------------------------------------------
  // Access the needed variables
  //----------------------------------------------------------------------------
  
#define MAX_CSCTF_TRK 36    // max # of CSCTF tracks per BX
#define MAX_LCTS_PER_TRK 4  // max # of LCTS which form a CSCTF track
#define MAX_RPC_LCTS_EVT 36

  // CSC objects ------------------------------------------------------------
  int SizeTrk;


  float PtTrk   [MAX_CSCTF_TRK];
  double EtaTrk [MAX_CSCTF_TRK];
  double PhiTrk [MAX_CSCTF_TRK];
  int EtaBitTrk [MAX_CSCTF_TRK];
  int PhiBitTrk [MAX_CSCTF_TRK];
  int PtBitTrk  [MAX_CSCTF_TRK];
  int ModeTrk   [MAX_CSCTF_TRK];

  // Info about LCTs forming the track
  // Stored in 2D array [track x] [LCT y] - each event has x tracks containing y Lcts
  int NumLctsTrk[MAX_CSCTF_TRK];  // How many Lcts made up track

  double trLctglobalPhi [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  double trLctglobalEta [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctPhiBit       [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctEtaBit       [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctEtaBitMatt   [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctSector       [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctSubSector    [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctStation      [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctChamber      [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctEndcap       [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctBx           [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctCSCId        [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];

  // RPC objects ------------------------------------------------------------

  // track variables
  int rpc_SizeTrk;

  float rpc_gblEtaTrk [MAX_CSCTF_TRK];
  
  // Trigger Primitve variables
  int Size_rpc_LCTs;
  
  double rpc_gblEta  [MAX_RPC_LCTS_EVT];
  double rpc_gblPhi  [MAX_RPC_LCTS_EVT];
  unsigned rpc_strip [MAX_RPC_LCTS_EVT];
  unsigned rpc_layer [MAX_RPC_LCTS_EVT];
  int rpc_bx         [MAX_RPC_LCTS_EVT];

};



#endif
