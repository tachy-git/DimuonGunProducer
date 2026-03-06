#!/usr/bin/env bash
# Arrays
PTS=(40)
PHIS=(0 0.785375 1.57075 3.1415 3.926875 4.71225)
#RADII=(15 40 80 130 200 300 500 600 800 1000 1200 1500) #mm
RADII=(15 1500)
NEVENTS=1000
QUEUE=1

# Base submit description
cat > submit_dimuongun_base.sub <<'EOT'
universe              = vanilla
executable            = runDimuonGun.sh

output                = joblogs/job.$(Cluster).$(Process).out
error                 = joblogs/job.$(Cluster).$(Process).err
log                   = joblogs/job.$(Cluster).$(Process).log
accounting_group      = group_cms

stream_output         = True
stream_error          = True
should_transfer_files = YES
request_memory = 4GB
EOT

mkdir -p joblogs

# Submit
for pt in "${PTS[@]}"; do
  for r in "${RADII[@]}"; do
    for phi in "${PHIS[@]}"; do
    batch="DimuonGun_pT-${pt}_R-${r}_PHI-${phi//./p}"
    echo "Submitting: ${batch}"

    condor_submit submit_dimuongun_base.sub \
      -append "arguments = ${pt} ${r} ${phi} ${NEVENTS}" \
      -append "JobBatchName = ${batch}" \
      -append "environment = PROCESS=\$(Process)"  \
      -append "queue $QUEUE"
      done
  done
done
