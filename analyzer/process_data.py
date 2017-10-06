import FWCore.ParameterSet.Config as cms

process = cms.Process("TriggerAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source ("PoolSource",
                            fileNames = cms.untracked.vstring(
                            '/store/data/Run2017B/JetHT/MINIAOD/23Jun2017-v1/00000/004DBDB2-C859-E711-8DD0-002590D0B042.root',
                            '/store/data/Run2017B/JetHT/MINIAOD/23Jun2017-v1/00000/00C46816-E359-E711-9063-002590D0B066.root'
                            '/store/data/Run2017B/JetHT/MINIAOD/23Jun2017-v1/00000/027DC76A-0759-E711-8D64-0025902008F4.root'),
                            )

process.analyzer = cms.EDAnalyzer('TriggerAnalyzer',
                                DoPreselection  = cms.untracked.bool(True),
                                DoFilter        = cms.untracked.bool(True),
                                DoTightJetID    = cms.untracked.bool(True),

                                jets            = cms.InputTag('slimmedJetsAK8'),
                                htjets            = cms.InputTag('slimmedJets'),
                                vertices        = cms.InputTag("offlineSlimmedPrimaryVertices"),

                                TriggerResults  = cms.InputTag('TriggerResults','','HLT'),
                                TriggerObjects  = cms.InputTag("slimmedPatTrigger"),

                                NoiseFilter     = cms.InputTag('TriggerResults','', 'RECO'),
                                noiseFilterSelection_HBHENoiseFilterLoose = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Loose"),
                                noiseFilterSelection_HBHENoiseFilterTight = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Tight"),
                                noiseFilterSelection_HBHENoiseIsoFilter = cms.InputTag("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult"),

                                noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
                                noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
                                noiseFilterSelection_globalTightHalo2016Filter = cms.string('Flag_globalTightHalo2016Filter'),
                                noiseFilterSelection_muonBadTrackFilter = cms.string('Flag_muonBadTrackFilter'),
                                noiseFilterSelection_chargedHadronTrackResolutionFilter = cms.string('Flag_chargedHadronTrackResolutionFilter'),

                                # parameter
                                EtaCut          = cms.untracked.double(2.5),
                                PtCut           = cms.untracked.double(200.0),
                                dEtaCut         = cms.untracked.double(1.3),
                                HTEtaCut        = cms.untracked.double(2.4),
                                HTPtCut         = cms.untracked.double(30.0),
                                target          = cms.untracked.string('fakeroot_csv.root'),
                                triggernames    = cms.untracked.string( 'HLT_PFHT1050_v,'+
                                                                        'HLT_PFHT890_v,'+
                                                                        'HLT_PFHT780_v,'+
                                                                        'HLT_PFHT430_v,'+
                                                                        'HLT_AK8PFHT750_TrimMass50_v,'+
                                                                        'HLT_AK8PFHT800_TrimMass50_v,'+
                                                                        'HLT_AK8PFHT850_TrimMass50_v,'+
                                                                        'HLT_AK8PFHT900_TrimMass50_v,'+
                                                                        'HLT_PFJet140_v,'+
                                                                        'HLT_PFJet200_v,'+
                                                                        'HLT_PFJet400_v,'+
                                                                        'HLT_PFJet450_v,'+
                                                                        'HLT_PFJet500_v,'+
                                                                        'HLT_PFJet550_v,'+
                                                                        'HLT_AK8PFJet140_v,'+
                                                                        'HLT_AK8PFJet200_v,'+
                                                                        'HLT_AK8PFJet400_v,'+
                                                                        'HLT_AK8PFJet450_v,'+
                                                                        'HLT_AK8PFJet500_v,'+
                                                                        'HLT_IsoMu20_v,'+
                                                                        'HLT_IsoMu27_v,'+
                                                                        'HLT_Mu20_v,'+
                                                                        'HLT_Mu50_v'),
                                )


# MET Filter
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
##___________________________BadChargedCandidate_Noise_Filter________________________________|| 
process.load('Configuration.StandardSequences.Services_cff')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.debug = cms.bool(False)
process.BadChargedCandidateSequence = cms.Sequence (process.BadChargedCandidateFilter)





process.p = cms.Path()
process.p += process.HBHENoiseFilterResultProducer
process.p += process.BadChargedCandidateSequence
process.p += process.analyzer