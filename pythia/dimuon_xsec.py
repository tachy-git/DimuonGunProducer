import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from FWCore.ParameterSet.VarParsing import VarParsing
import math

options = VarParsing('analysis')
options.register('process', 'Inclusive', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'DY, SoftQCD, SoftQCDall, QCDShowerGamma, Inclusive, or Brem')
options.register('beam', 'proton', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'proton, pi+, or pi-')
options.register('target', 'proton', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'proton, neutron, or Cu')
options.register('eBeam', 137.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, 'Beam energy in GeV')
options.register('filter', 'none', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'none, dimuon, gamma, meson')
options.register('muMinP', 10.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, 'Minimum momentum cut for each muon')
options.register('iSeed', 42, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'Random seed')
options.register('outFile', 'gen.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Output GEN EDM ROOT file')
options.parseArguments()

N_EVENTS = 10000

def fixed_target_boost(Ebeam, m_beam, m_target):
    p_beam = math.sqrt(max(Ebeam * Ebeam - m_beam * m_beam, 0.0))
    s = m_beam * m_beam + m_target * m_target + 2.0 * Ebeam * m_target
    sqrts = math.sqrt(max(s, 0.0))
    beta_cm = p_beam / (Ebeam + m_target)
    gamma_cm = (Ebeam + m_target) / sqrts
    return p_beam, sqrts, beta_cm, gamma_cm

def cm_to_lab(pstar, costh, mass, beta_cm, gamma_cm):
    estar = math.sqrt(pstar * pstar + mass * mass)
    sinth = math.sqrt(max(0.0, 1.0 - costh * costh))
    pt_lab = pstar * sinth
    pz_lab = gamma_cm * (pstar * costh + beta_cm * estar)
    p_lab = math.sqrt(pt_lab * pt_lab + pz_lab * pz_lab)
    e_lab = gamma_cm * (estar + beta_cm * pstar * costh)
    return pt_lab, pz_lab, p_lab, e_lab

beam_info = {
    'proton': (2212, 0.9382720813),
    'pi+':    (211,  0.13957039),
    'pi-':    (-211, 0.13957039),
}

target_info = {
    'proton':  (2212,       0.9382720813),
    'neutron': (2112,       0.9395654133),
    'Cu':      (1000290630, 63.0 * 0.9314941024),
}

if options.beam not in beam_info:
    raise ValueError("Unknown beam '{}'. Valid: {}".format(options.beam, list(beam_info.keys())))

if options.target not in target_info:
    raise ValueError("Unknown target '{}'. Valid: {}".format(options.target, list(target_info.keys())))

BEAM_PDG, m_beam = beam_info[options.beam]
TARGET_PDG, m_target = target_info[options.target]

p_beam, sqrts, beta_cm, gamma_cm = fixed_target_boost(options.eBeam, m_beam, m_target)

m_mu = 0.105658
example_pstar = 5.0
example_lines = []
for costh in (-1.0, 0.0, 1.0):
    pt_lab, pz_lab, p_lab, e_lab = cm_to_lab(example_pstar, costh, m_mu, beta_cm, gamma_cm)
    example_lines.append("costh={:+.1f} -> p_lab={:.6f} GeV".format(costh, p_lab))

print("[Config] process={} beam={}({}) target={}({}) eBeam={} GeV filter={} muMinP={} GeV nEvents={} outFile={}".format(
    options.process, options.beam, BEAM_PDG, options.target, TARGET_PDG, options.eBeam, options.filter, options.muMinP, N_EVENTS, options.outFile
))
print("[Kinematics] p_beam={:.6f} GeV sqrt(s)={:.6f} GeV beta_cm={:.8f} gamma_cm={:.8f}".format(
    p_beam, sqrts, beta_cm, gamma_cm
))
print("[Example CM->lab] for p* = {:.3f} GeV: {}".format(example_pstar, " ; ".join(example_lines)))

process_params = {
    'SoftQCD': [
        'SoftQCD:nonDiffractive = on',
    ],

    'SoftQCDall': [
        'SoftQCD:all = on',
    ],

    'DY': [
        'WeakSingleBoson:ffbar2gmZ = on',
        '23:onMode = off',
        '23:onIfMatch = 13 -13',
    ],

    'QCDShowerGamma': [
        'HardQCD:all = on',
        'PartonLevel:ISR = on',
        'PartonLevel:FSR = on',
        'SpaceShower:QEDshowerByQ = on',
        'TimeShower:QEDshowerByQ = on',
        'SpaceShower:QEDshowerByL = off',
        'TimeShower:QEDshowerByL = off',
        'PhaseSpace:pTHatMin = 1.0',
    ],

    'Inclusive': [
        'SoftQCD:all = on',
        'PartonLevel:ISR = on',
        'PartonLevel:FSR = on',
        'SpaceShower:QEDshowerByQ = on',
        'TimeShower:QEDshowerByQ = on',
        'SpaceShower:QEDshowerByL = off',
        'TimeShower:QEDshowerByL = off',
    ],

    'Brem': [
        'PhotonCollision:gmgm2mumu = on',
        'PDF:beamA2gamma = on',
        'PDF:beamB2gamma = on',
        'PhaseSpace:mHatMin = 0.2',
        'PhaseSpace:mHatMax = {:.6f}'.format(sqrts),
    ],
}

if options.process not in process_params:
    raise ValueError("Unknown process '{}'. Valid: {}".format(options.process, list(process_params.keys())))

pythia_process_params = list(process_params[options.process])

angantyr_params = []
if options.target == 'Cu':
    angantyr_params = ['HeavyIon:mode = 2']

process = cms.Process('GEN')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(N_EVENTS))
process.source = cms.Source("EmptySource")
process.source.firstRun = cms.untracked.uint32(1)

process.genstepfilter.triggerConditions = cms.vstring("generation_step")

process.generator = cms.EDFilter(
    "Pythia8GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    comEnergy = cms.double(max(sqrts, 1.0)),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        processParameters = cms.vstring(
            'Beams:idA = {}'.format(BEAM_PDG),
            'Beams:idB = {}'.format(TARGET_PDG),
            'Beams:frameType = 2',
            'Beams:eA = {:.8f}'.format(options.eBeam),
            'Beams:eB = 0.0',
            'PartonLevel:MPI = off',
            *angantyr_params,
            *pythia_process_params
        ),
        parameterSets = cms.vstring(
            'pythia8CommonSettings',
            'processParameters',
        )
    )
)

process.mesonFilter = cms.EDFilter(
    "MCSmartSingleParticleFilter",
    ParticleID = cms.untracked.vint32(113),
    MinPt = cms.untracked.vdouble(0.0),
    MinEta = cms.untracked.vdouble(-9.9),
    MaxEta = cms.untracked.vdouble(9.9),
    Status = cms.untracked.vint32(2),
)

process.gammaFilter = cms.EDFilter(
    "MCSmartSingleParticleFilter",
    ParticleID = cms.untracked.vint32(22),
    MinPt = cms.untracked.vdouble(0.0),
    MinEta = cms.untracked.vdouble(-9.9),
    MaxEta = cms.untracked.vdouble(9.9),
    Status = cms.untracked.vint32(1),
)

process.dimuonFilter = cms.EDFilter(
    "MCParticlePairFilter",
    ParticleID1 = cms.untracked.vint32(13),
    ParticleID2 = cms.untracked.vint32(13),
    ParticleCharge = cms.untracked.int32(-1),
    MinPt = cms.untracked.vdouble(0.0, 0.0),
    MinP = cms.untracked.vdouble(options.muMinP, options.muMinP),
    MinEta = cms.untracked.vdouble(-9.9, -9.9),
    MaxEta = cms.untracked.vdouble(9.9, 9.9),
    Status = cms.untracked.vint32(1, 1),
    MinInvMass = cms.untracked.double(0.0),
    MaxInvMass = cms.untracked.double(999999.0),
)

if options.filter == 'none':
    active_filter = None
elif options.filter == 'dimuon':
    active_filter = process.dimuonFilter
elif options.filter == 'meson':
    active_filter = process.mesonFilter
elif options.filter == 'gamma':
    active_filter = process.gammaFilter
else:
    raise ValueError("Unknown filter '{}'. Valid: none, dimuon, meson, gamma".format(options.filter))

process.xsecAnalyzer = cms.EDAnalyzer(
    "GenXSecAnalyzer",
    genFilterInfoTag = cms.InputTag("generator"),
)

if active_filter is None:
    process.generation_step = cms.Path(process.generator * process.pgen)
else:
    process.generation_step = cms.Path(process.generator * active_filter * process.pgen)

process.GENoutput = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string(options.outFile),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep edmHepMCProduct_generatorSmeared__GEN',
        'keep edmGenEventInfoProduct_generator__GEN',
        'keep GenRunInfoProduct_generator__GEN',
        'keep recoGenParticles_genParticles__GEN',
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN'),
        filterName = cms.untracked.string('')
    )
)

process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.xsec_step = cms.EndPath(process.xsecAnalyzer)
process.out_step = cms.EndPath(process.GENoutput)
process.endjob_step = cms.EndPath(process.endOfProcess)

process.schedule = cms.Schedule(
    process.generation_step,
    process.genfiltersummary_step,
    process.xsec_step,
    process.out_step,
    process.endjob_step,
)

process.RandomNumberGeneratorService.generator.initialSeed = options.iSeed
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
