#include <cmath>
#include <memory>

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/AbstractServices/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/Units.h"

#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

class JPsiDimuonGunProducer : public edm::global::EDProducer<> {
public:
  explicit JPsiDimuonGunProducer(const edm::ParameterSet&);
  ~JPsiDimuonGunProducer() override = default;

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  double minPt_, maxPt_;
  double minEta_, maxEta_;
  double minPhi_, maxPhi_;
  double vx_, vy_, vz_;

  int jpsiId_;
  int muMinusId_;
  int muPlusId_;
  double mJpsi_;
  double mMu_;

  int verbosity_;
};

JPsiDimuonGunProducer::JPsiDimuonGunProducer(const edm::ParameterSet& pset)
    : minPt_(pset.getParameter<double>("MinPt")),
      maxPt_(pset.getParameter<double>("MaxPt")),
      minEta_(pset.getParameter<double>("MinEta")),
      maxEta_(pset.getParameter<double>("MaxEta")),
      minPhi_(pset.getParameter<double>("MinPhi")),
      maxPhi_(pset.getParameter<double>("MaxPhi")),
      vx_(pset.getParameter<double>("Vx")),
      vy_(pset.getParameter<double>("Vy")),
      vz_(pset.getParameter<double>("Vz")),
      jpsiId_(pset.getParameter<int>("JPsiPdgId")),
      muMinusId_(pset.getParameter<int>("MuMinusPdgId")),
      muPlusId_(pset.getParameter<int>("MuPlusPdgId")),
      mJpsi_(pset.getParameter<double>("JPsiMass")),
      mMu_(pset.getParameter<double>("MuonMass")),
      verbosity_(pset.getUntrackedParameter<int>("Verbosity", 0)) {
  produces<edm::HepMCProduct>("unsmeared");
  produces<GenEventInfoProduct>();
}

void JPsiDimuonGunProducer::produce(edm::StreamID sid, edm::Event& e, const edm::EventSetup&) const {
  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &rng->getEngine(sid);

  // 1) Sample J/psi kinematics
  const double pt = CLHEP::RandFlat::shoot(engine, minPt_, maxPt_);
  const double eta = CLHEP::RandFlat::shoot(engine, minEta_, maxEta_);
  const double phi = CLHEP::RandFlat::shoot(engine, minPhi_, maxPhi_);

  const double px = pt * std::cos(phi);
  const double py = pt * std::sin(phi);
  const double pz = pt * std::sinh(eta);
  const double eJ = std::sqrt(px * px + py * py + pz * pz + mJpsi_ * mJpsi_);

  TLorentzVector p4J(px, py, pz, eJ);

  // 2) Decay J/psi -> mu+ mu- using TGenPhaseSpace
  TGenPhaseSpace phsp;
  double masses[2] = {mMu_, mMu_};
  phsp.SetDecay(p4J, 2, masses);
  phsp.Generate();

  TLorentzVector* p4d0 = phsp.GetDecay(0);
  TLorentzVector* p4d1 = phsp.GetDecay(1);

  // 3) Build HepMC event
  auto genEvent = std::make_unique<HepMC::GenEvent>(HepMC::Units::GEV, HepMC::Units::MM);
  genEvent->set_event_number(e.id().event());
  genEvent->set_signal_process_id(20);

  const double time = std::sqrt(vx_ * vx_ + vy_ * vy_ + vz_ * vz_);

  // Production vertex
  auto* vProd = new HepMC::GenVertex(HepMC::FourVector(vx_, vy_, vz_, time));

  // J/psi (status 2)
  auto* jpsi = new HepMC::GenParticle(HepMC::FourVector(px, py, pz, eJ), jpsiId_, 2);
  vProd->add_particle_out(jpsi);
  genEvent->add_vertex(vProd);

  // Decay vertex (same position for now)
  auto* vDec = new HepMC::GenVertex(HepMC::FourVector(vx_, vy_, vz_, time));
  vDec->add_particle_in(jpsi);

  // Assign mu- / mu+ (ordering doesn't matter; PDG IDs define charge)
  auto* muMinus = new HepMC::GenParticle(
      HepMC::FourVector(p4d0->Px(), p4d0->Py(), p4d0->Pz(), p4d0->E()), muMinusId_, 1);
  auto* muPlus = new HepMC::GenParticle(
      HepMC::FourVector(p4d1->Px(), p4d1->Py(), p4d1->Pz(), p4d1->E()), muPlusId_, 1);

  vDec->add_particle_out(muMinus);
  vDec->add_particle_out(muPlus);
  genEvent->add_vertex(vDec);

  if (verbosity_ > 0) {
    genEvent->print();
  }

  // GenEventInfoProduct: make it BEFORE transferring ownership
  auto genInfo = std::make_unique<GenEventInfoProduct>(genEvent.get());
  e.put(std::move(genInfo));

  // HepMCProduct takes ownership of GenEvent*
  auto out = std::make_unique<edm::HepMCProduct>();
  out->addHepMCData(genEvent.release());
  e.put(std::move(out), "unsmeared");
}

void JPsiDimuonGunProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<double>("MinPt", 20.0);
  desc.add<double>("MaxPt", 100.0);
  desc.add<double>("MinEta", -2.4);
  desc.add<double>("MaxEta", 2.4);
  desc.add<double>("MinPhi", -M_PI);
  desc.add<double>("MaxPhi", M_PI);

  desc.add<double>("Vx", 0.0);
  desc.add<double>("Vy", 0.0);
  desc.add<double>("Vz", 0.0);

  desc.add<int>("JPsiPdgId", 443);
  desc.add<int>("MuMinusPdgId", 13);
  desc.add<int>("MuPlusPdgId", -13);

  desc.add<double>("JPsiMass", 3.0969);
  desc.add<double>("MuonMass", 0.105658);

  desc.addUntracked<int>("Verbosity", 0);

  descriptions.add("jpsiDimuonGunProducer", desc);
}

DEFINE_FWK_MODULE(JPsiDimuonGunProducer);
