from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TagAndProbe_DY'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'TagAndProbe_Zmumu.py'
config.JobType.allowUndistributedCMSSW = False

config.section_("Data")
config.Data.inputDataset = '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.publication = False
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.publishDataName = 'TagAndProbe'
config.Data.outLFNDirBase = '/store/group/lpcsusyhad/PHYS14_720_Mar14_2014_v2/hua/'

config.Data.ignoreLocality = False

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
