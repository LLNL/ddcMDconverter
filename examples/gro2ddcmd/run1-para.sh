#!/bin/bash

cDir=`pwd`
cd $cDir/para

martini2obj  -i itpList -t proItpList -p  martini_v2.1-dna.itp -o martini.data -l molecule.data

# copy the martini.data molecule.data ConsAtom.data to ddcmd simulation directory
cp martini.data molecule.data ConsAtom.data $cDir/ddcmd
