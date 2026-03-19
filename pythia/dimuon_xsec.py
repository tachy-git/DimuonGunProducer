import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
from FWCore.ParameterSet.VarParsing import VarParsing
import math

options = VarParsing('analysis')
options.register('process', 'SoftQCD',  VarParsing.multiplicity.singleton, VarParsing.varType.string, 'DY, SoftQCD, SoftQCDall, QCDShowerGamma, or Brem')
options.register('beam',    'proton',   VarParsing.multiplicity.singleton, VarParsing.varType.string, 'proton, pi+, or pi-')
options.register('target',  'proton',   VarParsing.multiplicity.singleton, VarParsing.varType.string, 'proton, neutron, or Cu')
options.register('eBeam',   137.,       VarParsing.multiplicity.singleton, VarParsing.varType.float,  'Beam energy in GeV')
options.register('filter',  'none',     VarParsing.multiplicity.singleton, VarParsing.varType.string, 'none, dimuon, gamma, rho0, omega, phi, eta, etap, pi0')
options.register('iSeed',   42,         VarParsing.multiplicity.singleton, VarParsing.varType.int,    'Random seed')
options.parseArguments()

# ======================================================
# Beam info: (PDG ID, mass in GeV)
# ======================================================
beam_info = {
    'proton': ( 2212,  0.93827),
    'pi+':    (  211,  0.13957),
    'pi-':    ( -211,  0.13957),
}

# ======================================================
# Target info: (PDG ID, mass in GeV)
# ======================================================
target_info = {
    'proton':  ( 2212,       0.93827),
    'neutron': ( 2112,       0.93827),
    'Cu':      ( 1000290630, 0.93827),
}

if options.beam not in beam_info:
    raise ValueError(f"Unknown beam '{options.beam}'. Valid: {list(beam_info.keys())}")
if options.target not in target_info:
    raise ValueError(f"Unknown target '{options.target}'. Valid: {list(target_info.keys())}")

BEAM_PDG,   m_beam   = beam_info[options.beam]
TARGET_PDG, m_target = target_info[options.target]

# ======================================================
# Fixed-target kinematics
# ======================================================
sqrts = math.sqrt(m_beam**2 + m_target**2 + 2.0 * options.eBeam * m_target)

print(
    f"[Config] process={options.process}  beam={options.beam}({BEAM_PDG})  "
    f"target={options.target}({TARGET_PDG})  eBeam={options.eBeam} GeV  "
    f"sqrts={sqrts:.3f} GeV  filter={options.filter}"
)

# ======================================================
# Process parameter sets
# ======================================================
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

    'Brem': [
        'PhotonCollision:gmgm2mumu = on',
        'PDF:beamA2gamma = on',
        'PDF:beamB2gamma = on',
        f'PhaseSpace:mHatMin = 0.2',
        f'PhaseSpace:mHatMax = {sqrts:.4f}',
    ],
}

if options.process not in process_params:
    raise ValueError(f"Unknown process '{options.process}'. Valid: {list(process_params.keys())}")

pythia_process_params = list(process_params[options.process])

# ======================================================
# Angantyr for Cu nuclear target
# ======================================================
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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))
process.source    = cms.Source("EmptySource")
process.genstepfilter.triggerConditions = cms.vstring("generation_step")

process.generator = cms.EDFilter(
    "Pythia8GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency      = cms.untracked.double(1.0),
    pythiaHepMCVerbosity  = cms.untracked.bool(False),
    comEnergy             = cms.double(sqrts),
    maxEventsToPrint      = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
        processParameters = cms.vstring(
            f'Beams:idA = {BEAM_PDG}',
            f'Beams:idB = {TARGET_PDG}',
            #f'Beams:eA = 137',
            #f'Beams:eB = 1',
            'PartonLevel:MPI = off',
            *angantyr_params,
            *pythia_process_params,
        ),
        parameterSets = cms.vstring(
            'pythia8CommonSettings',
            'pythia8CUEP8M1Settings',
            'processParameters',
        )
    )
)

# ======================================================
# One filter for each meson
# ======================================================
process.rho0Filter = cms.EDFilter(
    "MCSmartSingleParticleFilter",
    ParticleID = cms.untracked.vint32(113),
    MinPt      = cms.untracked.vdouble(0.),
    MinEta     = cms.untracked.vdouble(-9.9),
    MaxEta     = cms.untracked.vdouble( 9.9),
    Status     = cms.untracked.vint32(2),
)

process.omegaFilter = cms.EDFilter(
    "MCSmartSingleParticleFilter",
    ParticleID = cms.untracked.vint32(223),
    MinPt      = cms.untracked.vdouble(0.),
    MinEta     = cms.untracked.vdouble(-9.9),
    MaxEta     = cms.untracked.vdouble( 9.9),
    Status     = cms.untracked.vint32(2),
)

process.phiFilter = cms.EDFilter(
    "MCSmartSingleParticleFilter",
    ParticleID = cms.untracked.vint32(333),
    MinPt      = cms.untracked.vdouble(0.),
    MinEta     = cms.untracked.vdouble(-9.9),
    MaxEta     = cms.untracked.vdouble( 9.9),
    Status     = cms.untracked.vint32(2),
)

process.etaFilter = cms.EDFilter(
    "MCSmartSingleParticleFilter",
    ParticleID = cms.untracked.vint32(221),
    MinPt      = cms.untracked.vdouble(0.),
    MinEta     = cms.untracked.vdouble(-9.9),
    MaxEta     = cms.untracked.vdouble( 9.9),
    Status     = cms.untracked.vint32(2),
)

process.etapFilter = cms.EDFilter(
    "MCSmartSingleParticleFilter",
    ParticleID = cms.untracked.vint32(331),
    MinPt      = cms.untracked.vdouble(0.),
    MinEta     = cms.untracked.vdouble(-9.9),
    MaxEta     = cms.untracked.vdouble( 9.9),
    Status     = cms.untracked.vint32(2),
)

process.pi0Filter = cms.EDFilter(
    "MCSmartSingleParticleFilter",
    ParticleID = cms.untracked.vint32(111),
    MinPt      = cms.untracked.vdouble(0.),
    MinEta     = cms.untracked.vdouble(-9.9),
    MaxEta     = cms.untracked.vdouble( 9.9),
    Status     = cms.untracked.vint32(2),
)

process.gammaFilter = cms.EDFilter(
    "MCSmartSingleParticleFilter",
    ParticleID = cms.untracked.vint32(22),
    MinPt      = cms.untracked.vdouble(0.),
    MinEta     = cms.untracked.vdouble(-9.9),
    MaxEta     = cms.untracked.vdouble( 9.9),
    Status     = cms.untracked.vint32(1),
)

# ======================================================
# Opposite-sign dimuon filter
# ======================================================
process.dimuonFilter = cms.EDFilter(
    "MCParticlePairFilter",
    ParticleID1   = cms.untracked.vint32(13),
    ParticleID2   = cms.untracked.vint32(13),
    ParticleCharge = cms.untracked.int32(-1),
    MinPt         = cms.untracked.vdouble(0., 0.),
    MinEta        = cms.untracked.vdouble(-9.9, -9.9),
    MaxEta        = cms.untracked.vdouble( 9.9,  9.9),
    Status        = cms.untracked.vint32(1, 1),
    MinInvMass    = cms.untracked.double(0.0),
    MaxInvMass    = cms.untracked.double(999.),
)

# ======================================================
# Select active filter
# ======================================================
if options.filter == 'none':
    active_filter = None
elif options.filter == 'dimuon':
    active_filter = process.dimuonFilter
elif options.filter == 'rho0':
    active_filter = process.rho0Filter
elif options.filter == 'omega':
    active_filter = process.omegaFilter
elif options.filter == 'phi':
    active_filter = process.phiFilter
elif options.filter == 'eta':
    active_filter = process.etaFilter
elif options.filter == 'etap':
    active_filter = process.etapFilter
elif options.filter == 'pi0':
    active_filter = process.pi0Filter
elif options.filter == 'gamma':
    active_filter = process.gammaFilter
else:
    raise ValueError("Unknown filter '{}'. Valid: none, dimuon, rho0, omega, phi, eta, etap, pi0".format(options.filter))

process.xsecAnalyzer = cms.EDAnalyzer(
    "GenXSecAnalyzer",
    genFilterInfoTag = cms.InputTag("generator"),
)

if active_filter is None:
    process.generation_step = cms.Path(process.generator * process.pgen)
else:
    process.generation_step = cms.Path(process.generator * active_filter * process.pgen)

process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.xsec_step             = cms.EndPath(process.xsecAnalyzer)
process.endjob_step           = cms.EndPath(process.endOfProcess)

process.schedule = cms.Schedule(
    process.generation_step,
    process.genfiltersummary_step,
    process.xsec_step,
    process.endjob_step,
)

process.source.firstRun = cms.untracked.uint32(1)
process.RandomNumberGeneratorService.generator.initialSeed = options.iSeed
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
