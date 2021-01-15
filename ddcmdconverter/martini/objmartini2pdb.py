__author__ = 'zhang30'

import argparse
from ddcmdconverter.base import Obj
from ddcmdconverter.martini import ITP

def getArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--obj', action='store', dest='objfile', default='atom#.data', help='ddcMD object input file (default=atom#.data).')
    parser.add_argument('-o', '--pdb', action='store', dest='pdbfile', default='test.pdb', help='PDB output file (default=test.pdb).')
    parser.add_argument('-t', '--pro', action='store', dest='profile', default=None,
                        help='Protein ITP file list (default=proItpList).')
    parser.add_argument('-c', '--cut', action='store', dest='cutoff', type=float, help='Cutoff for bond.')
    parser.add_argument('-x', '--bx', action='store', dest='x', type=float, help='Boundary x.')
    parser.add_argument('-y', '--by', action='store', dest='y', type=float, help='Boundary y.')
    parser.add_argument('-z', '--bz', action='store', dest='z', type=float, help='Boundary z.')
    #parser.add_argument('-hx', '--hex', action='store_true', dest='hex', default=False,
    #                    help='gid is in hex format')
    args = parser.parse_args()

    return args


def main():
    args=getArgs()
    print ("Default inputs: ",args.objfile, args.pdbfile)

    itpList=[]

    if args.profile != None:
        with open(args.profile, "r") as f:
            for line in f:
                tipFileName=line.rstrip("\n\r")
                itp = ITP.ITP(tipFileName)
                itpList.append(itp)


    atomNameMapCollection={}
    # fix the atom name in the
    for itp in itpList:
        count = 0
        atomNameMap={}
        newResName=itp.header.moleculetype.data['name']
        for atom in itp.header.moleculetype.atoms.data:
            newAtomname= 'P' + str(count)
            residueAtom=atom['resname'].decode("utf-8")+":"+str(atom['resnr'])+":"+atom['atomname'].decode("utf-8")
            atomNameMap[newAtomname]=residueAtom
            #print newAtomname, residueAtom
            count = count + 1
        atomNameMapCollection[newResName]=atomNameMap

    obj=Obj.Obj()
    obj.parseHeader(args)
    if obj.datatype=='FIXRECORDBINARY':
        obj.toBinaryMartiniPDB(args, atomNameMapCollection)
    else:
        obj.toMartiniPDB(args, atomNameMapCollection)

if __name__ == '__main__':
    main()
