import FWCore.ParameterSet.Config as cms

process = cms.Process("TrigEff")

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('Configuration.StandardSequences.ReconstructionCosmics_cff')

process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")
process.load("RecoMuon.TrackingTools.MuonTrackLoader_cff")

# global tag
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_P_V41::All'

# message logger

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
##process.load('FWCore/MessageService/MessageLogger_cfi')
##process.MessageLogger.cout.placeholder = cms.untracked.bool(False)
##process.MessageLogger.cout.threshold = cms.untracked.string('WARNING')
##process.MessageLogger.cout.threshold = cms.untracked.string('INFO')
##process.MessageLogger.cerr.threshold = cms.untracked.string('ERROR')
##process.MessageLogger.debugModules = cms.untracked.vstring('*')

#Tracer: uncomment helps debugging
#process.Tracer = cms.Service("Tracer")

# how many events?
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("Configuration.StandardSequences.Generator_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")


# ------------ Trig Eff ----------------
process.te = cms.EDAnalyzer("TrigEff",
                            process.MuonServiceProxy,
                            genTag       = cms.InputTag("genParticles"),
                            L1extraTag   = cms.InputTag("l1extraParticles"),
                            muonsTag     = cms.InputTag("muons"),
                            csctfTag     = cms.InputTag("csctfDigis"),
                            csctfLctsTag = cms.InputTag("csctfDigis"),
                            outputFile   = cms.string("MinimumBiasTrigEffNtuple-OneFileLocal.root"),
                            cscSegTag    = cms.InputTag("cscSegments"),
                            printLevel   = cms.untracked.int32(2)
                            )
# --------- Trig Eff END ----------------


#===============================================================================
process.load( "HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi" )

# triggerSummaryAnalyzerAOD,
# singleMuNoIso & hltEnd,
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#from PhysicsTools.HepMCCandAlgos.genParticleCandidatesFast_cfi import *
#process.p = cms.Path(process.genParticleCandidates*process.te)

### L1 Extra needed for the analysis
##process.load('L1Trigger.Configuration.L1Extra_cff')
##
### path
##process.p = cms.Path(process.L1Extra*process.te)

process.p = cms.Path(process.csctfDigis * 
                     process.te)

# ------------ PoolSource -------------
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource",
                             fileNames = readFiles,
                             secondaryFileNames = secFiles
                             )
readFiles.extend([

    #'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RECO/PromptReco-v2/000/201/196/142C80D9-AAEB-E111-88DA-5404A638869C.root'
    'file:142C80D9-AAEB-E111-88DA-5404A638869C.root'
 ]

)

secFiles.extend( [

##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/084C3CD1-8AE9-E111-BD4F-001D09F25109.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/0C431673-64E9-E111-870C-BCAEC5364CED.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/14692F98-66E9-E111-AD69-0025B320384C.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/1C561F44-7BE9-E111-B99F-BCAEC5329709.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/1CB7061C-79E9-E111-B201-0025901D62A6.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/241D90DB-80E9-E111-8178-001D09F28E80.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/2EFC8481-70E9-E111-9D8B-BCAEC5329709.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/4094BC9A-61E9-E111-A0B3-BCAEC532971F.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/5AB03D68-73E9-E111-BD40-E0CB4E4408E7.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/76EF915E-6EE9-E111-AD64-0025901D5D9A.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/7C656956-7DE9-E111-BD56-E0CB4E553673.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/987A3FB7-83E9-E111-8593-0025B32445E0.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/A61C268E-75E9-E111-8FDC-BCAEC5329716.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/A6771F4D-5DE9-E111-834F-BCAEC53296F7.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/BE56B775-5FE9-E111-89D4-5404A63886B4.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/C4247142-87E9-E111-BFBB-001D09F23C73.root',
##    'rfio:/castor/cern.ch/cms/store/data/Run2012C/MinimumBias/RAW/v1/000/201/196/F22F5DE4-6AE9-E111-B4AB-BCAEC518FF52.root'

    'file:084C3CD1-8AE9-E111-BD4F-001D09F25109.root',
    'file:0C431673-64E9-E111-870C-BCAEC5364CED.root',
    'file:14692F98-66E9-E111-AD69-0025B320384C.root',
    'file:1C561F44-7BE9-E111-B99F-BCAEC5329709.root',
    'file:1CB7061C-79E9-E111-B201-0025901D62A6.root',
    'file:241D90DB-80E9-E111-8178-001D09F28E80.root',
    'file:2EFC8481-70E9-E111-9D8B-BCAEC5329709.root',
    'file:4094BC9A-61E9-E111-A0B3-BCAEC532971F.root',
    'file:5AB03D68-73E9-E111-BD40-E0CB4E4408E7.root',
    'file:76EF915E-6EE9-E111-AD64-0025901D5D9A.root',
    'file:7C656956-7DE9-E111-BD56-E0CB4E553673.root',
    'file:987A3FB7-83E9-E111-8593-0025B32445E0.root',
    'file:A61C268E-75E9-E111-8FDC-BCAEC5329716.root',
    'file:A6771F4D-5DE9-E111-834F-BCAEC53296F7.root',
    'file:BE56B775-5FE9-E111-89D4-5404A63886B4.root',
    'file:C4247142-87E9-E111-BFBB-001D09F23C73.root',
    'file:F22F5DE4-6AE9-E111-B4AB-BCAEC518FF52.root'
    
    ]
                 
)
# -------- PoolSource END -------------
