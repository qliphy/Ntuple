import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/data/Run2016B/DoubleMuon/MINIAOD/23Sep2016-v3/00000/9C5757BC-569B-E611-BE0F-0019B9CACF1A.root'),
    secondaryFileNames = cms.untracked.vstring()
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)

process.pfJetIDSelector = cms.PSet(
    quality = cms.string('LOOSE'),
    version = cms.string('RUNIISTARTUP')
)

process.Ztoee = cms.EDProducer("CandViewCombiner",
    cut = cms.string(''),
    decay = cms.string('goodElectrons@+ goodElectrons@-')
)


process.Ztomumu = cms.EDProducer("CandViewCombiner",
    cut = cms.string(''),
    decay = cms.string('goodMuons@+ goodMuons@-')
)


process.cleanAK4Jets = cms.EDProducer("PATJetCleaner",
    checkOverlaps = cms.PSet(
        electrons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.4),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(True),
            src = cms.InputTag("goodElectrons")
        ),
        muons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.4),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(True),
            src = cms.InputTag("goodMuons")
        ),
        photons = cms.PSet(

        ),
        taus = cms.PSet(

        ),
        tkIsoElectrons = cms.PSet(

        )
    ),
    finalCut = cms.string('pt > 20 & abs(eta) < 4.7'),
    preselection = cms.string(''),
    src = cms.InputTag("goodAK4Jets")
)


process.goodElectrons = cms.EDProducer("PATElectronIdSelector",
    effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
    idLabel = cms.string('tight'),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedElectrons"),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.goodMuons = cms.EDProducer("PATMuonIdSelector",
    idLabel = cms.string('tight'),
    src = cms.InputTag("slimmedMuons"),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.leptonicV = cms.EDProducer("CandViewMerger",
    cut = cms.string(''),
    src = cms.VInputTag("Ztoee", "Ztomumu")
)


process.looseMuons = cms.EDProducer("PATMuonIdSelector",
    idLabel = cms.string('loose'),
    src = cms.InputTag("slimmedMuons"),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.photonIDValueMapProducer = cms.EDProducer("PhotonIDValueMapProducer",
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
    esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedESRecHits"),
    particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons"),
    pfCandidates = cms.InputTag("particleFlow"),
    pfCandidatesMiniAOD = cms.InputTag("packedPFCandidates"),
    src = cms.InputTag("gedPhotons"),
    srcMiniAOD = cms.InputTag("slimmedPhotons","","@skipCurrentProcess"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    verticesMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.vetoElectrons = cms.EDProducer("PATElectronIdSelector",
    effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
    idLabel = cms.string('veto'),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedElectrons"),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.BadChargedCandidateFilter = cms.EDFilter("BadChargedCandidateFilter",
    PFCandidates = cms.InputTag("packedPFCandidates"),
    debug = cms.bool(False),
    innerTrackRelErr = cms.double(1.0),
    maxDR = cms.double(1e-05),
    minMuonPt = cms.double(100.0),
    minMuonTrackRelErr = cms.double(2.0),
    minPtDiffRel = cms.double(1e-05),
    muons = cms.InputTag("slimmedMuons"),
    segmentCompatibility = cms.double(0.3),
    taggingMode = cms.bool(False)
)


process.BadPFMuonFilter = cms.EDFilter("BadPFMuonFilter",
    PFCandidates = cms.InputTag("packedPFCandidates"),
    algo = cms.int32(14),
    debug = cms.bool(False),
    innerTrackRelErr = cms.double(1.0),
    minDZ = cms.double(0.1),
    minMuPt = cms.double(100),
    minPtError = cms.double(2.0),
    muons = cms.InputTag("slimmedMuons"),
    segmentCompatibility = cms.double(0.3),
    taggingMode = cms.bool(False)
)


process.goodAK4Jets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
    filterParams = cms.PSet(
        quality = cms.string('LOOSE'),
        version = cms.string('RUNIISTARTUP')
    ),
    src = cms.InputTag("slimmedJets")
)


process.goodPhotons = cms.EDFilter("PATPhotonSelector",
    cut = cms.string('pt > 15 && abs(eta) < 2.5'),
    src = cms.InputTag("slimmedPhotons")
)


process.leptonicVFilter = cms.EDFilter("CandViewCountFilter",
    filter = cms.bool(False),
    minNumber = cms.uint32(0),
    src = cms.InputTag("leptonicV")
)


process.leptonicVSelector = cms.EDFilter("CandViewSelector",
    cut = cms.string('pt > 0.0'),
    filter = cms.bool(False),
    src = cms.InputTag("leptonicV")
)


process.treeDumper = cms.EDAnalyzer("ZPKUTreeMaker",
    PKUChannel = cms.string('VW_CHANNEL'),
    RunOnMC = cms.bool(False),
    ak4jetsSrc = cms.InputTag("cleanAK4Jets"),
    beamSpot = cms.InputTag("offlineBeamSpot","","RECO"),
    conversions = cms.InputTag("reducedEgamma","reducedConversions","RECO"),
    crossSectionPb = cms.double(1),
    effAreaChHadFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring15/effAreaPhotons_cone03_pfChargedHadrons_25ns_NULLcorrection.txt'),
    effAreaNeuHadFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring15/effAreaPhotons_cone03_pfNeutralHadrons_25ns_90percentBased.txt'),
    effAreaPhoFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring15/effAreaPhotons_cone03_pfPhotons_25ns_90percentBased.txt'),
    elPaths1 = cms.vstring('HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*'),
    elPaths2 = cms.vstring('HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*'),
    electrons = cms.InputTag("slimmedElectrons"),
    full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
    genJet = cms.InputTag("slimmedGenJets"),
    genSrc = cms.InputTag("prunedGenParticles"),
    generator = cms.InputTag("generator"),
    hltToken = cms.InputTag("TriggerResults","","HLT"),
    isGen = cms.bool(False),
    jecAK4PayloadNames = cms.vstring('Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFchs.txt', 
        'Summer16_23Sep2016HV3_DATA_L2Relative_AK4PFchs.txt', 
        'Summer16_23Sep2016HV3_DATA_L3Absolute_AK4PFchs.txt', 
        'Summer16_23Sep2016HV3_DATA_L2L3Residual_AK4PFchs.txt'),
    jecAK4chsPayloadNames = cms.vstring('Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFchs.txt', 
        'Summer16_23Sep2016HV3_DATA_L2Relative_AK4PFchs.txt', 
        'Summer16_23Sep2016HV3_DATA_L3Absolute_AK4PFchs.txt', 
        'Summer16_23Sep2016HV3_DATA_L2L3Residual_AK4PFchs.txt'),
    leptonicVSrc = cms.InputTag("leptonicV"),
    looseelectronSrc = cms.InputTag("vetoElectrons"),
    loosemuonSrc = cms.InputTag("looseMuons"),
    metSrc = cms.InputTag("slimmedMETs"),
    muPaths1 = cms.vstring('HLT_Mu17_TrkIsoVVL_v*'),
    muPaths2 = cms.vstring('HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*', 
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*'),
    muPaths3 = cms.vstring('HLT_IsoMu24_v*'),
    muPaths4 = cms.vstring('HLT_Mu17_v*'),
    noiseFilter = cms.InputTag("TriggerResults","","RECO"),
    noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
    noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),
    noiseFilterSelection_HBHENoiseIsoFilter = cms.string('Flag_HBHENoiseIsoFilter'),
    noiseFilterSelection_badChargedHadron = cms.InputTag("BadChargedCandidateFilter"),
    noiseFilterSelection_badMuon = cms.InputTag("BadPFMuonFilter"),
    noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
    noiseFilterSelection_globalTightHaloFilter = cms.string('Flag_globalTightHalo2016Filter'),
    noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),
    originalNEvents = cms.int32(1),
    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
    phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer","phoNeutralHadronIsolation"),
    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
    photonSrc = cms.InputTag("slimmedPhotons"),
    pileup = cms.InputTag("slimmedAddPileupInfo"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    t1jetSrc = cms.InputTag("slimmedJets"),
    t1muSrc = cms.InputTag("slimmedMuons"),
    targetLumiInvPb = cms.double(1.0),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.muSequence = cms.Sequence(process.goodMuons+process.looseMuons)


process.eleSequence = cms.Sequence(process.goodElectrons+process.vetoElectrons)


process.leptonicVSequence = cms.Sequence(process.Ztomumu+process.Ztoee+process.leptonicV)


process.NJetsSequence = cms.Sequence(process.goodAK4Jets+process.cleanAK4Jets)


process.photonSequence = cms.Sequence(process.goodPhotons)


process.metfilterSequence = cms.Sequence(process.BadPFMuonFilter+process.BadChargedCandidateFilter)


process.jetSequence = cms.Sequence(process.NJetsSequence)


process.leptonSequence = cms.Sequence(process.muSequence+process.eleSequence+process.leptonicVSequence+process.leptonicVSelector+process.leptonicVFilter)


process.analysis = cms.Path(process.leptonSequence+process.jetSequence+process.metfilterSequence+process.photonIDValueMapProducer+process.treeDumper)


process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring('FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary'),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(99999999),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(200)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    destinations = cms.untracked.vstring('warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        optionalPSet = cms.untracked.bool(True),
        placeholder = cms.untracked.bool(True)
    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    )
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('ZtreePKU.root')
)


process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring('HCAL', 
        'ZDC', 
        'CASTOR', 
        'EcalBarrel', 
        'EcalEndcap', 
        'EcalPreshower', 
        'TOWER')
)


process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


process.CaloTowerGeometryFromDBEP = cms.ESProducer("CaloTowerGeometryFromDBEP",
    applyAlignment = cms.bool(False),
    hcalTopologyConstants = cms.PSet(
        maxDepthHB = cms.int32(2),
        maxDepthHE = cms.int32(3),
        mode = cms.string('HcalTopologyMode::LHC')
    )
)


process.CaloTowerTopologyEP = cms.ESProducer("CaloTowerTopologyEP")


process.CastorDbProducer = cms.ESProducer("CastorDbProducer")


process.CastorGeometryFromDBEP = cms.ESProducer("CastorGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.EcalBarrelGeometryFromDBEP = cms.ESProducer("EcalBarrelGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder")


process.EcalEndcapGeometryFromDBEP = cms.ESProducer("EcalEndcapGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")


process.EcalPreshowerGeometryFromDBEP = cms.ESProducer("EcalPreshowerGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


process.GlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


process.HcalAlignmentEP = cms.ESProducer("HcalAlignmentEP")


process.HcalGeometryFromDBEP = cms.ESProducer("HcalGeometryFromDBEP",
    applyAlignment = cms.bool(True),
    hcalTopologyConstants = cms.PSet(
        maxDepthHB = cms.int32(2),
        maxDepthHE = cms.int32(3),
        mode = cms.string('HcalTopologyMode::LHC')
    )
)


process.MuonDetLayerGeometryESProducer = cms.ESProducer("MuonDetLayerGeometryESProducer")


process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    compatibiltyWith11 = cms.untracked.bool(True),
    useDDD = cms.untracked.bool(False)
)


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer")


process.ZdcGeometryFromDBEP = cms.ESProducer("ZdcGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.fakeForIdealAlignment = cms.ESProducer("FakeAlignmentProducer",
    appendToDataLabel = cms.string('fakeForIdeal')
)


process.hcalDDDRecConstants = cms.ESProducer("HcalDDDRecConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalDDDSimConstants = cms.ESProducer("HcalDDDSimConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalTopologyIdeal = cms.ESProducer("HcalTopologyIdealEP",
    Exclude = cms.untracked.string(''),
    appendToDataLabel = cms.string('')
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)


process.idealForDigiCSCGeometry = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.idealForDigiDTGeometry = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.idealForDigiTrackerGeometry = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(cms.PSet(
        record = cms.string('SiPixelQualityFromDbRcd'),
        tag = cms.string('')
    ), 
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        ))
)


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(cms.PSet(
        Label = cms.untracked.string(''),
        NormalizationFactor = cms.untracked.double(1.0),
        Record = cms.string('SiStripApvGainRcd')
    ), 
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(cms.PSet(
        record = cms.string('SiStripDetVOffRcd'),
        tag = cms.string('')
    ), 
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        )),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.trackerGeometryDB = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.trackerNumberingGeometryDB = cms.ESProducer("TrackerGeometricDetESModule",
    appendToDataLabel = cms.string(''),
    fromDDD = cms.bool(False)
)


process.trackerTopology = cms.ESProducer("TrackerTopologyEP",
    appendToDataLabel = cms.string('')
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('80X_dataRun2_2016SeptRepro_v4'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string(''),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)


process.eegeom = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd')
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    GainWidthsForTrigPrims = cms.bool(False),
    HERecalibration = cms.bool(False),
    HEreCalibCutoff = cms.double(20.0),
    HFRecalibration = cms.bool(False),
    iLumi = cms.double(-1.0),
    testHFQIE10 = cms.bool(False),
    toGet = cms.untracked.vstring('GainWidths')
)


process.prefer("es_hardcode")

