########################################################
## me11.py   A script to analyze me11 lcts information
##
## By David Curry
##
########################################################

import sys
from ROOT import *
import numpy as np


# Set the print level. Default = 0
if len(sys.argv) is 1: printLevel = 0
else: printLevel = sys.argv[1]

# ====== Import Root file
#file = TFile("GlobalRun_CSCTF_222608.root")
file = TFile('/afs/cern.ch/work/d/dcurry/public/trigEff/withPier/CMSSW_6_2_0_pre5/src/L1CSCTFAnalysis/TrigEff/test/GlobalRun_CSCTF.root')
#file = TFile('/exports/uftrig01a/dcurry/data/trig_eff/FC9D8210-E2A7-E211-8E8B-485B39800C34.root')

# ========================


# Define Hist file to be saved
newfile = new TFile("all_twoHit_RPC_pt.root","recreate")


# Set the branch address of TTree in Tfile
evt  = file.Get("csctfTTree")
reco = file.Get("recoMuons")

# bind methods to improve speed
getEntry  = evt.GetEntry
getREntry = reco.GetEntry

# ============== Initialize Histograms =================

hdphi = TH1F('hdphi', '', 50, -10, 10)
hdeta = TH1F('hdeta', '', 50, -10, 10)


# ======================================================


# Loop over over events in TFile
for iEvt in range(evt.GetEntries()):

    # for testing
    if iEvt > 10000: break
    
    getEntry(iEvt)
    getREntry(iEvt)
    
    if iEvt % 10000 is 0: print 'Event #', iEvt

    if printLevel > 0:
        print '\n============== New Event # ', reco.Event, ' =================\n'
        print '  Tracks in Event = ', evt.SizeTrk
        print '  Event Lcts      = ', evt.SizeLCTs
        
        
    # Only look at 1 track events, Eta in Me1 range(eta < 1.1)
    if evt.SizeTrk is not 1: continue
    
    if abs(evt.EtaTrk[0]) < 2.1: continue
    
    if evt.NumLCTsTrk[0] is 1: continue
    
    if printLevel > 0: print '---> One CSCTF Track in present in ME11 range'
    
    # Loop over tracks
    for iCSCTrk in range(evt.SizeTrk):
        
        if printLevel > 0: print '\n==== Looping over Track #', iCSCTrk, 'with', evt.NumLCTsTrk[iCSCTrk], 'Lcts' 
        
        # Loop over lcts in track
        for iLct in range(evt.NumLCTsTrk[iCSCTrk]):
            
            # We only want Lcts in ME11
            if evt.trLctStation[iCSCTrk*4 + iLct] is not 1: continue
            
            isLctMatched = False

            if printLevel > 0:
                print '\n   ===== Looping over Track Lct #', iLct
                print ' trLct Station = ', evt.trLctStation[iCSCTrk*4 + iLct]
                print ' trLct Sector  = ', evt.trLctSector[iCSCTrk*4 + iLct]
                print ' trLct SSector = ', evt.trLctSubSector[iCSCTrk*4 + iLct]
                print ' trLct Endcap  = ', evt.trLctEndcap[iCSCTrk*4 + iLct]
                print ' trLct Gbl Phi = ', evt.trLctglobalPhi[iCSCTrk*4 + iLct]
                print ' trLct Gbl Eta = ', evt.trLctglobalEta[iCSCTrk*4 + iLct]
                print ' trLct CSCID   = ', evt.trLctTriggerCSCID[iCSCTrk*4 + iLct]

            
            # Loop over Lcts in event to find match to current track lct
            for iCsc in range(evt.SizeLCTs):

                if isLctMatched: break
                
                if printLevel > 0:
                    print '\n     ==== Looping over Event Lct #', iCsc
                    print ' Lct Station = ', evt.lctStation[iCsc]
                    print ' Lct Sector  = ', evt.lctSector[iCsc]
                    print ' Lct SSector = ', evt.lctSubSector[iCsc]
                    print ' Lct Endcap  = ', evt.lctEndcap[iCsc]
                    print ' Lct Station = ', evt.lctStation[iCsc]
                    print ' Lct Gbl Phi = ', evt.lctglobalPhi[iCsc]
                    print ' Lct Gbl Eta = ', evt.lctglobalEta[iCsc]
                    print ' Lct CSCID   = ', evt.lctTriggerCSCID[iCsc]
                    
                    
                # compare event and track lct variables
                if evt.trLctStation[iCSCTrk*4 + iLct]      is not evt.lctStation[iCsc]   or \
                   evt.trLctSector[iCSCTrk*4 + iLct]       is not evt.lctSector[iCsc]    or \
                   evt.trLctSubSector[iCSCTrk*4 + iLct]    is not evt.lctSubSector[iCsc] or \
                   evt.trLctTriggerCSCID[iCSCTrk*4 + iLct] is not evt.lctTriggerCSCID[iCsc] or \
                   evt.trLctEndcap[iCSCTrk*4 + iLct]       is not evt.lctEndcap[iCsc] : continue
                
                
                isLctMatched = True;
                
                if printLevel > 0: print '------> Track Lct and Event Lct are matched'
                
                # If trLct and Lct match then we compare their eta/phi in a hist
                hdphi -> Fill(evt.lctglobalPhi[iCsc] - evt.trLctglobalPhi[iCSCTrk*4 + iLct])
                hdeta -> Fill(evt.lctglobalEta[iCsc] - evt.trLctglobalEta[iCSCTrk*4 + iLct]) 
                
                    
            # end Lct loop
                
        # end loop over track Lcts
                        
    # end loop over tracks
            
# end loop over events


#  ======== Write Hists ==========

hdphi -> Write()
hdeta -> Write()

delete newfile

# ================================



            