# ddcMD Converter: a Python tool to convert GROMACS files to ddcMD inputs

Authors: Xiaohua Zhang, Francesco Di Natale, and Jim Glosli

The code uses python3 syntax after version v1.0.0.

## Convert GROMACS inputs to ddcMD inputs
Here is a sample script to Convert GROMACS inputs to ddcMD inputs.
ddcMD requires the residue name to be unique. So ION residue name is replaced with the actual name
"+" and "-" are reserved key words inn ddcMD that cannot be used in residue/atom names.

```
gmx editconf -f ../topol.tpr -o system.pdb > gmx-editconf.log

sed 's/ZN+ ION/ZN   ZN/' system.pdb | \
sed 's/MG+ ION/MG   MG/' | \
sed 's/NA+ ION/NA   NA/' | \
sed 's/CL- ION/CL   CL/'  \
 > system_fix.pdb

sed 's/NA+ /NA  /' ../topol.top | sed 's/CL- /CL  /' | sed 's/TERNARY /RAS_RAF /' > topol.top

ln -s $sDir/martini.data
ln -s $sDir/ConsAtom.data
ln -s $sDir/object.data
ln -s $sDir/resItpList
ln -s $sDir/POPX_Martini_v2.0_lipid.itp
ln -s $sDir/molecule.data

pdbmartini2obj -p system_fix.pdb -t topol.top -m martini.data -f ConsAtom.data  -o atoms#000000
restraint -i atoms#000000 -o restraint.data -p resItpList

mkdir snapshot.mem
mv atoms#000000 restart snapshot.mem
ln -s snapshot.mem/restart
```

