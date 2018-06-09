# Python tools for ddcMD 

The code uses python3 syntax after version v1.0.0.

## pdb2obj
Convert pdb file to ddcMD anatomy file.

```
[zhang30@soifon ddcMDconverter]$ python pdb2obj.py -h
usage: pdb2obj.py [-h] [-t TOPFILE] [-p PDBFILE] [-o OBJFILE] [-s SPEFILE]
                  [-l SPLFILE] [-x X] [-y Y] [-z Z]

optional arguments:
  -h, --help            show this help message and exit
  -t TOPFILE, --top TOPFILE
                        CHARMM topology file (default=top_all22_prot.inp).
  -p PDBFILE, --pdb PDBFILE
                        PDB input file (default=test.pdb).
  -o OBJFILE, --obj OBJFILE
                        ddcMD object output file (default=atom#.data).
  -s SPEFILE, --spe SPEFILE
                        ddcMD species all output file (default=species.data).
  -l SPLFILE, --spl SPLFILE
                        ddcMD species less output file (default=speless.data).
  -x X, --bx X          Boundary x.
  -y Y, --by Y          Boundary y.
  -z Z, --bz Z          Boundary z.

```

An usual way to run the script:
```
python ../../pdb2obj.py -t top_all36_prot.rtf -p prot-mini.pdb  -o atom#.data
```
In the PDB file header, the box size should be defined (starting with "CRYST1").
```
CRYST1   24.00   24.00   24.00  90.00  90.00  90.00  P 21 21 21   4
ATOM      1  N   GLY     1      -5.514  -0.055  -0.000  1.00  0.00      PROT
ATOM      2  HT1 GLY     1      -6.344   0.570  -0.000  1.00  0.00      PROT
ATOM      3  HT2 GLY     1      -5.536  -0.678   0.836  1.00  0.00      PROT
ATOM      4  HT3 GLY     1      -5.535  -0.678  -0.837  1.00  0.00      PROT
ATOM      5  CA  GLY     1      -4.216   0.717   0.000  1.00  0.00      PROT
ATOM      6  HA1 GLY     1      -4.166   1.303   0.908  1.00  0.00      PROT
ATOM      7  HA2 GLY     1      -4.166   1.304  -0.907  1.00  0.00      PROT
...
```
Otherwise the box size should be set in the options.

```
python ../../pdb2obj.py -t top_all36_prot.rtf -p prot-mini.pdb  -o atom#.data -x 24 -y 24 -z 24
```

## obj2pdb

```
[zhang30@soifon ddcMDconverter]$ python obj2pdb.py -h
usage: obj2pdb.py [-h] [-i OBJFILE] [-o PDBFILE] [-c CUTOFF] [-x X] [-y Y]
                  [-z Z]

optional arguments:
  -h, --help            show this help message and exit
  -i OBJFILE, --obj OBJFILE
                        ddcMD object input file (default=atom#.data).
  -o PDBFILE, --pdb PDBFILE
                        PDB output file (default=test.pdb).
  -c CUTOFF, --cut CUTOFF
                        Cutoff for bond.
  -x X, --bx X          Boundary x.
  -y Y, --by Y          Boundary y.
  -z Z, --bz Z          Boundary z.

```


## CHARMM Psf reader



## install the software

```
ssh://git@cz-bitbucket.llnl.gov:7999/xzr/ddcmdconvertor.git
cd ddcmdconvertor
pip3 install -e .
```
