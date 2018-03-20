__author__ = 'zhang30'

import argparse
from ddcmdconverter.Obj import Obj
from ddcmdconverter.ITP import ITP

def getArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--obj', action='store', dest='objfile', default='atom#.data', help='ddcMD object input file (default=atom#.data).')
    parser.add_argument('-o', '--res', action='store', dest='resfile', default='restraint.data', help='Restraint file (default=restraint.data).')
    parser.add_argument('-p', '--itp', action='store', dest='itpfile', default=None,
                        help='ITP file list (default=itpList).')

    args = parser.parse_args()

    return args

if __name__ == '__main__':

    args=getArgs()
    print "Default inputs: ", args.objfile, args.resfile, args.itpfile

    itpList=[]

    if args.itpfile != None:
        with open(args.itpfile, "r") as f:
            for line in f:
                tipFileName=line.rstrip("\n\r")
                itp = ITP.ITP(tipFileName)
                itpList.append(itp)

    restraintMap={}
    for itp in itpList:
        resName = itp.header.moleculetype.data['name']
        sectionKeys=itp.header.moleculetype.sections.keys()

        if 'position_restraints' in sectionKeys:
            restraintList = []
            restraintSize = len(itp.header.moleculetype.position_restraints.data)
            for i in range(restraintSize):
                restraint=itp.header.moleculetype.position_restraints.data[i]
                restraintList.append(restraint)
            restraintMap[resName]=restraintList

    obj=Obj.Obj()
    obj.toRestraint(args, restraintMap)