#!/bin/bash

    cDir=`pwd`

    cd ddcmd
    pdbmartini2obj -p $cDir/gromacs/lipids-water-eq4.pdb -t $cDir/gromacs/system.top -m martini.data -f ConsAtom.data  -o atoms#000000

    cp $cDir/para/resItpList .
    cp $cDir/para/POPX_Martini_v2.0_lipid.itp .
    restraint -i atoms#000000 -o restraint.data -p resItpList

    mkdir -p snapshot.mem
    mv atoms#000000  restart snapshot.mem
    ln -s snapshot.mem/restart
    cp $cDir/para/object.data .
  
