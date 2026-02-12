import FWCore.ParameterSet.Config as cms

process = cms.Process("GUN")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1000))

process.source = cms.Source("EmptySource")

process.RandomNumberGeneratorService = cms.Service(
    "RandomNumberGeneratorService",
    generator=cms.PSet(
        initialSeed=cms.untracked.uint32(12345),
        engineName=cms.untracked.string("HepJamesRandom"),
    ),
)

process.generator = cms.EDProducer(
    "JPsiDimuonGunProducer",
    MinPt=cms.double(20.0),
    MaxPt=cms.double(100.0),
    MinEta=cms.double(-2.4),
    MaxEta=cms.double(2.4),
    MinPhi=cms.double(-3.14159265359),
    MaxPhi=cms.double(3.14159265359),
    Vx=cms.double(0.0),
    Vy=cms.double(0.0),
    Vz=cms.double(0.0),
    JPsiPdgId=cms.int32(443),
    MuMinusPdgId=cms.int32(13),
    MuPlusPdgId=cms.int32(-13),
    JPsiMass=cms.double(3.0969),
    MuonMass=cms.double(0.105658),
    Verbosity=cms.untracked.int32(0),
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName=cms.untracked.string("jpsi_dimuon_gun.root"),
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep HepMCProduct_*_*_*",
    ),
)

process.p = cms.Path(process.generator)
process.e = cms.EndPath(process.out)
