from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = 'DMu16H-v33'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.inputFiles = ['Spring16_25nsV6_DATA_L1FastJet_AK4PFchs.txt','Spring16_25nsV6_DATA_L2Relative_AK4PFchs.txt','Spring16_25nsV6_DATA_L3Absolute_AK4PFchs.txt','Spring16_25nsV6_DATA_L2L3Residual_AK4PFchs.txt']
# Name of the CMSSW configuration file
config.JobType.psetName    = 'Zanalysis_data.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/DoubleMuon/Run2016H-PromptReco-v3/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 18
config.Data.lumiMask = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
#config.Data.runRange = '246908-258750'
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'DMu16H-v33'

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.storageSite = 'T2_CH_CERN'



