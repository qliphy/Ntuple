from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName   = 'ZA-AQGC-Ntuple'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
#config.JobType.generator = 'lhe'
config.JobType.inputFiles = ['Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt','Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt','Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt']
# Name of the CMSSW configuration file
config.JobType.psetName    = 'Zanalysis.py'


config.section_("Data")
# This string determines the primary dataset of the newly-produced outputs.
# For instance, this dataset will be named /CrabTestSingleMu/something/USER
config.Data.inputDataset = '/ZA-AQGC-GEN/melu-ZA-AQGC-MiniAOD-9c82be0544462ee57448ab37ee74ae1a/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 472
config.Data.publication = False 
config.Data.outputDatasetTag = 'ZA-AQGC-Ntuple'

config.section_("Site")
# Where the output files will be transmitted to
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T3_US_FNALLPC'
