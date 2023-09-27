#!/bin/bash

source /usr/gapps/kras/spack6/share/spack/setup-env.sh

module load gcc/7.3.1
spack compiler find
spack compiler add

spack load -r ddcmdconverter@1.0.4

cDir=`pwd`
cd $cDir/para

martini2obj  -i itpList -t proItpList -p  martini_v2.1-dna.itp -o martini.data -l molecule.data

# copy the martini.data molecule.data ConsAtom.data to ddcmd simulation directory
cp martini.data molecule.data ConsAtom.data $cDir/ddcmd
