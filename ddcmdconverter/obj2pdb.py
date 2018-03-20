__author__ = 'zhang30'

import argparse
import Obj

def getArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--obj', action='store', dest='objfile', default='atom#.data', help='ddcMD object input file (default=atom#.data).')
    parser.add_argument('-o', '--pdb', action='store', dest='pdbfile', default='test.pdb', help='PDB output file (default=test.pdb).')
    parser.add_argument('-c', '--cut', action='store', dest='cutoff', type=float, help='Cutoff for bond.')
    parser.add_argument('-x', '--bx', action='store', dest='x', type=float, help='Boundary x.')
    parser.add_argument('-y', '--by', action='store', dest='y', type=float, help='Boundary y.')
    parser.add_argument('-z', '--bz', action='store', dest='z', type=float, help='Boundary z.')
    args = parser.parse_args()

    return args


def main():
    args=getArgs()
    print "Default inputs: ",args.objfile, args.pdbfile

    obj=Obj.Obj()
    obj.toPDB(args)


if __name__ == '__main__':
    main()
