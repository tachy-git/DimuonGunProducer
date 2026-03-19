#!/bin/bash
# ======================================================
# run_all.sh
# SoftQCDall: filters rho0, omega, phi, eta, etap, pi0
# QCDShowerGamma: filter dimuon
# ======================================================

ENERGIES="5 30 100 300 500 700 1000 1200 1500 2000"
BEAMS="proton pi+ pi-"
TARGETS="proton neutron Cu"
SOFTQCD_FILTERS="rho0 omega phi eta etap pi0"
ISEED=42
CSV="results.csv"

mkdir -p logs

echo "process,beam,target,eBeam_GeV,filter,xsec_pb" > ${CSV}

for EBEAM in $ENERGIES; do
    for BEAM in $BEAMS; do
        for TARGET in $TARGETS; do

            # beam label for filename (pi+ -> piplus, pi- -> piminus)
            BEAM_LABEL=$(echo ${BEAM} | sed 's/+/plus/g' | sed 's/-/minus/g')

            # ==========================================
            # SoftQCDall: one job per meson filter
            # ==========================================
            for FILTER in $SOFTQCD_FILTERS; do

                LOG="logs/SoftQCDall_${BEAM_LABEL}_${TARGET}_${EBEAM}GeV_${FILTER}.log"
                echo "[RUN] process=SoftQCDall  beam=${BEAM}  target=${TARGET}  eBeam=${EBEAM} GeV  filter=${FILTER}"

                cmsRun dimuon_xsec.py  \
                    process=SoftQCDall \
                    beam=${BEAM}       \
                    target=${TARGET}   \
                    eBeam=${EBEAM}     \
                    filter=${FILTER}   \
                    iSeed=${ISEED}     \
                    >& ${LOG}

                XSEC=$(grep "After filter: final cross section" ${LOG} \
                       | tail -1 | awk '{print $7}')
                [ -z "$XSEC" ] && XSEC="FAILED"

                echo "SoftQCDall,${BEAM},${TARGET},${EBEAM},${FILTER},${XSEC}" >> ${CSV}
                echo "       xsec = ${XSEC} pb"
                echo ""

            done

            # ==========================================
            # QCDShowerGamma: dimuon filter
            # ==========================================
            LOG="logs/QCDShowerGamma_${BEAM_LABEL}_${TARGET}_${EBEAM}GeV_dimuon.log"
            echo "[RUN] process=QCDShowerGamma  beam=${BEAM}  target=${TARGET}  eBeam=${EBEAM} GeV  filter=dimuon"

            cmsRun dimuon_xsec.py       \
                process=QCDShowerGamma  \
                beam=${BEAM}            \
                target=${TARGET}        \
                eBeam=${EBEAM}          \
                filter=dimuon           \
                iSeed=${ISEED}          \
                >& ${LOG}

            XSEC=$(grep "After filter: final cross section" ${LOG} \
                   | tail -1 | awk '{print $7}')
            [ -z "$XSEC" ] && XSEC="FAILED"

            echo "QCDShowerGamma,${BEAM},${TARGET},${EBEAM},dimuon,${XSEC}" >> ${CSV}
            echo "       xsec = ${XSEC} pb"
            echo ""

        done
    done
done

echo "========================================================"
echo " Results saved to ${CSV}"
echo "========================================================"
cat ${CSV}
echo "========================================================"
