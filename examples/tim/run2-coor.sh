#!/bin/bash

source /usr/gapps/kras/spack6/share/spack/setup-env.sh

module load gcc/7.3.1

spack load -r ddcmdconverter@1.0.4

    cDir=`pwd`

    cd ddcmd
    pdbmartini2obj -p $cDir/gromacs/equil.pdb -t $cDir/gromacs/topol.top -m martini.data -f ConsAtom.data  -o atoms#000000

    cp $cDir/para/resItpList .
    cp $cDir/para/POPX_Martini_v2.0_lipid.itp .
    restraint -i atoms#000000 -o restraint.data -p resItpList

    mkdir -p snapshot.mem
    mv atoms#000000  restart snapshot.mem
    ln -s snapshot.mem/restart
    cp $cDir/para/object.data .
  
