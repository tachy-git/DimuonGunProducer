#!/usr/bin/env bash

PT="${1:?PT missing}"
RADIUS="${2:?RADIUS missing}"
PHI="${3:?PHI missing}"
NEVENTS="${4:-1000}"
# Condor sets $PROCESS (0-based job index within the cluster)
PROCESS="${PROCESS:-0}"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cms/ldap_home/taehee/CMSSW_15_1_0/src
eval "$(scram runtime -sh)"

BASEDIR="/cms_scratch/taehee/DimuonGun"
CURRDIR="/cms/ldap_home/taehee/CMSSW_15_1_0/src/MyAnalyzer/JPsiDimuonGunProducer/test"
mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

CFGDIR="${BASEDIR}/cfg"
ROOTDIR="${BASEDIR}/root"
LOGDIR="${BASEDIR}/log"

mkdir -p "${CFGDIR}" "${ROOTDIR}" "${LOGDIR}"

PHI_LABEL="${PHI//./p}"
LABEL="pT-${PT}_R-${RADIUS}_PHI-${PHI_LABEL}_${PROCESS}"
OUTROOT="${ROOTDIR}/DimuonGun_${LABEL}.root"
CFG="${CFGDIR}/DimuonGun_${LABEL}_cfg.py"
LOG="${LOGDIR}/DimuonGun_${LABEL}.log"

echo "[INFO] PT=${PT}  RADIUS=${RADIUS} PHI=${PHI}  PROCESS=${PROCESS}"

cp "$CURRDIR/DimuonGun_cfg_template.py" "${CFG}"

sed -i \
  -e "s/__PT__/${PT}/g"           \
  -e "s/__RADIUS__/${RADIUS}/g"   \
  -e "s/__PHI__/${PHI}/g"   \
  -e "s/__NEVENTS__/${NEVENTS}/g" \
  -e "s|__OUTFILE__|${OUTROOT}|g" \
  "${CFG}"

echo "[INFO] About to run: cmsRun ${CFG} at $(date) on $(hostname)"
cmsRun "${CFG}" 2>&1 | tee ${LOG}
echo "[INFO] Finished ${OUTROOT}"
