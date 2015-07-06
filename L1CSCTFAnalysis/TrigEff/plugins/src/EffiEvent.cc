#include "../include/EffiEvent.h"


EffiEvent::EffiEvent():
  SizeTrk(0), Size_rpc_LCTs(0)
{
}


void EffiEvent::CloseFile(TFile *file) { file -> Close(); }

void EffiEvent::AttachToFile(TFile *file){

  if (!file) return;

  // get the Nutple  

  CSCtree = (TTree*) file -> Get("CSCtree");
  RPCtree = (TTree*) file -> Get("RPCtree");
  
  // CSCTF track variables
  CSCtree -> SetBranchAddress("SizeTrk"   ,  &SizeTrk  );
  CSCtree -> SetBranchAddress("EtaTrk"    ,  EtaTrk    );
  CSCtree -> SetBranchAddress("PhiTrk"    ,  PhiTrk    );
  CSCtree -> SetBranchAddress("PtTrk"     ,  PtTrk     );
  CSCtree -> SetBranchAddress("PtBitTrk"  ,  PtBitTrk  );
  CSCtree -> SetBranchAddress("EtaBitTrk" ,  EtaBitTrk );
  CSCtree -> SetBranchAddress("PhiBitTrk" ,  PhiBitTrk );
  CSCtree -> SetBranchAddress("ModeTrk"   ,  ModeTrk   );
  CSCtree -> SetBranchAddress("NumLctsTrk",  NumLctsTrk);

  // These are the variables LCTs that belong to the CSCTF track
  CSCtree -> SetBranchAddress("trLctglobalPhi" , trLctglobalPhi );
  CSCtree -> SetBranchAddress("trLctglobalEta" , trLctglobalEta );
  CSCtree -> SetBranchAddress("trLctPhiBit"    , trLctPhiBit    );
  CSCtree -> SetBranchAddress("trLctEtaBit"    , trLctEtaBit    );
  CSCtree -> SetBranchAddress("trLctEtaBitMatt", trLctEtaBitMatt);  
  CSCtree -> SetBranchAddress("trLctSector"    , trLctSector    );
  CSCtree -> SetBranchAddress("trLctSubSector" , trLctSubSector );
  CSCtree -> SetBranchAddress("trLctStation"   , trLctStation   );
  CSCtree -> SetBranchAddress("trLctChamber"   , trLctChamber   );
  CSCtree -> SetBranchAddress("trLctEndcap"    , trLctEndcap    );
  CSCtree -> SetBranchAddress("trLctBx"        , trLctBx        );
  CSCtree -> SetBranchAddress("trLctCSCId"    , trLctCSCId   );

  // --------- RPCTF Branches ------------------------------------------------------

  // RPC Trigger Primitive Branches(all hits in event record)
  RPCtree -> SetBranchAddress("Size_rpc_LCTs" ,  &Size_rpc_LCTs );
  RPCtree -> SetBranchAddress("rpc_gblEta"    ,  rpc_gblEta    );
  RPCtree -> SetBranchAddress("rpc_gblPhi"    ,  rpc_gblPhi    );
  RPCtree -> SetBranchAddress("rpc_layer"     ,  rpc_layer     );
  RPCtree -> SetBranchAddress("rpc_strip"     ,  rpc_strip     );
  RPCtree -> SetBranchAddress("rpc_bx"        ,  rpc_bx        );

}

void EffiEvent::GetEntry( long int iEvt ){

  CSCtree -> GetEntry(iEvt);
  RPCtree -> GetEntry(iEvt);


}
