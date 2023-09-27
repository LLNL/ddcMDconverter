#!/bin/bash

echo "Run cretsims for patch $1 struct $2, in local dir $3 using inpath $4"
echo "Init"
date

CODE_ROOT=/g/g90/helgi/mummi
echo "CODE_ROOT = $CODE_ROOT"
source $CODE_ROOT/setup/setup.env.sh

echo "Begin"
date

pwd

#jsrun -n 1 --bind=none -c 24 -a 1 python $CODE_ROOT/mummi/transformations/createsim/createsims_RAS_RAF_hack.py \
#  --gromacs gmx --loglevel 1 --fstype taridx \
#  --logpath ./ -c --mpi "gmx mdrun" \
#  --patch 'pfpatch_000000000300' --mdrunopt " -nt 96 -rdd 2.0 -ntomp 4 -dd 4 3 2" \
#  --inpath '/p/gpfs1/bhatia4/pilot2/campaign3/macro_patches_from_sr4' --outpath ./ --simnum test2 >> creatsims_out_test2.out 2>&1

# Used for RAS-RAF patches
#jsrun -n 1 --bind=none -c 24 -a 1 python $CODE_ROOT/mummi/transformations/createsim/createsims_RAS_RAF_hack.py \
#  --gromacs gmx --loglevel 1 --fstype taridx \
#  --logpath ./ -c --mpi "gmx mdrun" \
#  --patch $1 --mdrunopt " -nt 96 -rdd 2.0 -ntomp 4 -dd 4 3 2" \
#  --inpath '/p/gpfs1/bhatia4/pilot2/campaign1star/patches/' --outpath ./ --simnum $2 >> creatsims_${1}_${2}.out  2>&1

# Used for RAS only patches
#jsrun -n 1 --bind=none -c 24 -a 1 python $CODE_ROOT/mummi/transformations/createsim/createsims.py \
#  --gromacs gmx --loglevel 1 --fstype taridx \
#  --logpath ./ -c --mpi "gmx mdrun" \
#  --patch $1 --mdrunopt " -nt 96 -rdd 2.0 -ntomp 4 -dd 4 3 2" --outlocal "${3}" \
#  --inpath '/p/gpfs1/bhatia4/pilot2/campaign1star/patches/' --simnum "${2}" >> creatsims_${1}_${2}.out  2>&1

# Run creatsims

# TEMP intry point is not working on summit - change back when fixed
#jsrun -n 1 --bind=none -c 24 -a 1 python /gpfs/alpine/lrn005/scratch/helgi/mummi/mummi/transformations/createsim/createsims.py \
jsrun -n 1 --bind=none -c 24 -a 1 createsim \
  --gromacs gmx --loglevel 1 --fstype taridx \
  --logpath ./ -c --mpi "gmx mdrun" \
  --patch "create" --mdrunopt " -nt 96 -rdd 2.0 -ntomp 4 -dd 4 3 2" --outlocal ./ \
  --inpath '/p/gpfs1/bhatia4/pilot2/campaign1star/patches/' --simnum "RAS-RAF" >> creatsims.out  2>&1

echo "End"
date
