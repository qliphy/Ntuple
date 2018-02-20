import FWCore.ParameterSet.Config as cms

process = cms.Process( "TEST" )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
runOnMC = True
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
if runOnMC:
   process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
elif not(runOnMC):
   process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v4'


option = 'RECO'
process.load("VAJets.PKUCommon.goodElectrons_cff")
if option == 'RECO':
    process.goodElectrons.src = "slimmedElectrons"
    process.mediumElectrons.src = "slimmedElectrons"

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff']
for idmod in my_id_modules:
   setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.leptonSequence = cms.Sequence(
                                      process.egmGsfElectronIDSequence*process.eleSequence 
                                      )

process.treeDumper = cms.EDAnalyzer("EDBRTreeMaker",
                                    genSrc =  cms.InputTag("prunedGenParticles"),
                                    EleSrc0= cms.InputTag("slimmedElectrons"),  
                                    EleSrc1 = cms.InputTag("mediumElectrons"),
                                    EleSrc2 = cms.InputTag("goodElectrons"),
                                    )

process.analysis = cms.Path(process.leptonSequence +
                            process.treeDumper)

process.load("VAJets.PKUCommon.data.RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8")
process.source.fileNames = [
"/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/D2E83075-55BB-E611-9525-0025905B85F6.root",
]

process.maxEvents.input = 400000
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 2000
process.MessageLogger.cerr.FwkReport.limit = 99999999

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("check_EleId.root")
                                   )
