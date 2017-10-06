from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import datetime

#dataset = ['SingleMuon', 'Run2017B-PromptReco-v1']
#dataset = ['SingleMuon', 'Run2017B-PromptReco-v2']
#dataset = ['SingleMuon', 'Run2017C-PromptReco-v1']
#dataset = ['SingleMuon', 'Run2017C-PromptReco-v2']
#dataset = ['SingleMuon', 'Run2017C-PromptReco-v3']
#dataset = ['SingleMuon', 'Run2017D-PromptReco-v1']
#dataset = ['SingleMuon', 'Run2017E-PromptReco-v1']

#dataset = ['JetHT', 'Run2017B-PromptReco-v1']
#dataset = ['JetHT', 'Run2017B-PromptReco-v2']
#dataset = ['JetHT', 'Run2017C-PromptReco-v1']
#dataset = ['JetHT', 'Run2017C-PromptReco-v2']
#dataset = ['JetHT', 'Run2017C-PromptReco-v3']
#dataset = ['JetHT', 'Run2017D-PromptReco-v1']
#dataset = ['JetHT', 'Run2017E-PromptReco-v1']

name = 'TriggerStudies-{}-{}'.format(*dataset)

config = config()
config.General.requestName = name
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'process_data.py'
config.JobType.outputFiles = ['fakeroot_csv.root']
#config.JobType.scriptExe = 'userscript.sh'
#config.Data.outputPrimaryDataset = 'VBSWZ_H7VBFNLO_OF_LO_shower'
config.Data.splitting = 'LumiBased'
config.Data.inputDataset=   '/{}/{}/MINIAOD'.format(*dataset)

config.Data.lumiMask= '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-302654_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.unitsPerJob = 30
#config.Data.totalUnits = config.Data.unitsPerJob * 100

config.Data.outLFNDirBase = '/store/user/{}/TriggerStudies/'.format(getUsernameFromSiteDB())

#config.Data.publication = False
#config.Data.outputDatasetTag = 'July17_100k_LO'
config.Site.storageSite = 'T2_DE_DESY'
config.Site.blacklist = ['T2_US_Caltech',
                        'T2_US_Florida',
                        'T2_US_MIT',
                        'T2_US_Nebraska',
                        'T2_US_Purdue',
                        'T2_US_UCSD',
                        'T2_US_Vanderbilt',
                        'T2_US_Wisconsin']