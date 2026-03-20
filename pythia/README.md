# Dimuon production via Pythia
`dimuon_xsec.py` is the configuration file that mimics the dimuon production
where proton(or pion) collides with the detector system (proton, neutron, or heavy necleus like copper).
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd CMSSW_15_1_0/src
cmsenv
cmsRun dimuon_xsec.py
```

## beam
[Pythia beam parameter](https://pythia.org/latest-manual/BeamParameters.html)

- beam A: The beam energy is set according to a logarithmic scale
- beam B: Particle in the detector system
- [PDF](https://pythia.org/latest-manual/PDFSelection.html)
  - Not sure whether we can convert from lab frame to com frame, due to the PDF. Anyway...
  - proton: (default)NNPDF2.3 QCD+QED LO alpha_s(M_Z) = 0.130
  - pion: (default)GRS 99 L. ? / anyway, PDF exists also for pion beam
  - nucleus: provided as the number of proton and neutrons?
- `Beams:frameType = 2` : fixed-target configurations
- [heavy collision](https://pythia.org/latest-manual/HeavyIons.html)
  - The default heavy ion model(Angantyr) stacks parton level events, according to individual nucleon-nucleon
sub-collisions, on top of each other.

## process
- SoftQCD

- gammaShower

- DY

## filter


## 
