__author__ = 'zhang30'


import argparse
import Charmm
import Pdb
import Obj


def getArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--top', action='store', dest='topfile', default='top_all22_prot.inp', help='CHARMM topology file (default=top_all22_prot.inp).')
    parser.add_argument('-p', '--pdb', action='store', dest='pdbfile', default='test.pdb', help='PDB input file (default=test.pdb).')
    parser.add_argument('-o', '--obj', action='store', dest='objfile', default='atom#.data', help='ddcMD object output file (default=atom#.data).')
    parser.add_argument('-s', '--spe', action='store', dest='spefile', default='species.data', help='ddcMD species all output file (default=species.data).')
    parser.add_argument('-l', '--spl', action='store', dest='splfile', default='speless.data', help='ddcMD species less output file (default=speless.data).')
    args = parser.parse_args()

    return args

if __name__ == '__main__':

    args=getArgs()
    print "Default inputs: ",args.topfile, args.pdbfile, args.objfile, args.spefile, args.splfile

    charmmTop=Charmm.CharmmTop()
    charmmTop.parse(args.topfile)

    comPDB=Pdb.ComPDB()
    comPDB.parse(args.pdbfile)
    comPDB.assignGid(charmmTop)

    obj=Obj.Obj()
    obj.toObj(args.objfile, comPDB)


