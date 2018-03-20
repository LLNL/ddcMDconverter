__author__ = 'zhang30'


import argparse
import CharmmTop
import Pdb
import Obj
import Specie


def getArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--top', action='store', dest='topfile', default='top_all22_prot.inp', help='CHARMM topology file (default=top_all22_prot.inp).')
    parser.add_argument('-o', '--obj', action='store', dest='objfile', default='atom#.data', help='ddcMD object output file (default=atom#.data).')

    args = parser.parse_args()

    return args


def main():
    args=getArgs()
    print "Default inputs: ", args.topfile,  args.objfile

    print "Reading in CHARMM topology file ", args.topfile
    charmmTop=CharmmTop.CharmmTop()
    charmmTop.parse(args.topfile)


if __name__ == '__main__':
    main()
