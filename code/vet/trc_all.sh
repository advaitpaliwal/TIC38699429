#!/bin/zsh
# 10 runs, 3 at a time, single-threaded each
cd ~/Scratch/planet-hunt/v4/vet_38699429
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
seq 0 9 | xargs -P 3 -I{} sh -c 'taskpolicy -a ~/Scratch/planet-hunt/hunt4/.venv/bin/python -I glue/trc_run.py {} > trc/log_{}.txt 2>&1'
