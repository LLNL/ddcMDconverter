__author__ = 'zhang30'


import argparse
import ddcmdconverter.base.Pdb as Pdb
from ddcmdconverter.charmm.CharmmTop import CharmmTop
from ddcmdconverter.base.Obj import Obj
from ddcmdconverter.base.Specie import Specie


def getArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--top', action='store', dest='topfile', default='top_all22_prot.inp', help='CHARMM topology file (default=top_all22_prot.inp).')
    parser.add_argument('-p', '--pdb', action='store', dest='pdbfile', default='test.pdb', help='PDB input file (default=test.pdb).')
    parser.add_argument('-o', '--obj', action='store', dest='objfile', default='atom#.data', help='ddcMD object output file (default=atom#.data).')
    parser.add_argument('-s', '--spe', action='store', dest='spefile', default='species.data', help='ddcMD species all output file (default=species.data).')
    parser.add_argument('-l', '--spl', action='store', dest='splfile', default='speless.data', help='ddcMD species less output file (default=speless.data).')
    parser.add_argument('-x', '--bx', action='store', dest='x', type=float, default=64, help='Boundary x.')
    parser.add_argument('-y', '--by', action='store', dest='y', type=float, default=64, help='Boundary y.')
    parser.add_argument('-z', '--bz', action='store', dest='z', type=float, default=64, help='Boundary z.')
    args = parser.parse_args()

    return args


def main():
    args=getArgs()
    print ("Default inputs: ",args.topfile, args.pdbfile, args.objfile, args.spefile, args.splfile)

    print ("Reading in CHARMM topology file ", args.topfile)
    charmmTop=CharmmTop.CharmmTop()
    charmmTop.parse(args.topfile)

    print ("Reading in pdb file ", args.pdbfile)
    comPDB=Pdb.ComPDB()
    comPDB.parse(args)
    comPDB.assignGid(charmmTop)

    print ("Generating ddcMD object file ", args.objfile)
    obj=Obj.Obj()
    obj.toObj(args, comPDB)

    print ("Generating species file ", args.splfile)
    specie=Specie(charmmTop, comPDB)
    specie.toSpeData(args.splfile)


if __name__ == '__main__':
    main()
