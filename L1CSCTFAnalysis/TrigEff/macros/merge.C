#include "TChain.h"
#include "TFile.h"
#include <iostream>
#include "TFileMerger.h"
#include "TSystem.h"


using namespace std;

void merge() {

  hadd -f x.root f1.root f2.root f3.root
  


  /*  TChain* chain1 = new TChain("recoMuons");
  TChain* chain2 = new TChain("csctfTTree");


  chain1 -> AddFile("/cms/data/store/user/dcurry/trigeff/2012D/MinimumBiasTrigEffNtuple_new_CD_allevt_320_1_Ynm.root");
  chain1 -> AddFile("/cms/data/store/user/dcurry/trigeff/2012D/MinimumBiasTrigEffNtuple_new_CD_allevt_321_1_Duf.root");

  chain2 -> AddFile("/cms/data/store/user/dcurry/trigeff/2012D/MinimumBiasTrigEffNtuple_new_CD_allevt_320_1_Ynm.root");
  chain2 -> AddFile("/cms/data/store/user/dcurry/trigeff/2012D/MinimumBiasTrigEffNtuple_new_CD_allevt_321_1_Duf.root");


  cout << "Merging into MinBias_2012D.root (it may take a while)\n" << endl;
  chain1 -> Merge("MinBias_2012D_reco.root"); 
  chain2 -> Merge("MinBias_2012D_csctf.root");
  */


}


