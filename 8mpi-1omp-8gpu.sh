cat << EOF > select_gpu
#!/bin/bash
export TRACE_ALL_EVENTS=yes
export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
exec \$*
EOF

CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"

