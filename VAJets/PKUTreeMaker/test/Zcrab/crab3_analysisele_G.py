from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = 'DEle16G-v1'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.inputFiles = ['Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt','Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFchs.txt','Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFchs.txt','Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFchs.txt','Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFpuppi.txt','Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFpuppi.txt','Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFpuppi.txt','Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFpuppi.txt']
# Name of the CMSSW configuration file
config.JobType.psetName    = 'Zanalysis_data_G.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/DoubleEG/Run2016G-03Feb2017-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 80
config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#config.Data.runRange = '246908-258750'
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'DEle16G-v1'

config.section_("Site")
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T3_US_FNALLPC'


