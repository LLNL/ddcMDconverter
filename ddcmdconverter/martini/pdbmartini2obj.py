__author__ = 'zhang30'


import argparse
import sys

import ddcmdconverter.base.Pdb as Pdb
import ddcmdconverter.base.Gro as Gro
from ddcmdconverter.martini.ITP import ITP
#import Specie


def getArgs(in_options):

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pdb', action='store', dest='pdbfile', default='test.pdb', help='PDB input file (default=test.pdb).')
    parser.add_argument('-g', '--gro', action='store', dest='grofile', default=None,
                        help='GRO input file (default=None).')
    parser.add_argument('-t', '--pro', action='store', dest='profile', default=None,
                        help='Protein ITP file list (default=proItpList).')
    parser.add_argument('-f', '--fre', action='store', dest='frefile', default=None,
                        help='Constrant atom name for setting free in group.')
    parser.add_argument('-o', '--obj', action='store', dest='objfile', default='atom#.data',
                        help='ddcMD object output file (default=atom#.data).')
    parser.add_argument('-c', '--cut', action='store', dest='cutoff', type=float, help='Cutoff for bond.')
    parser.add_argument('-r', '--ori', action='store', dest='origin', default=True,
                        help='Move origin to -1/2L.')

    args = parser.parse_args(in_options)

    return args

def renameProt(comPDB):
    for molID, molPDB in enumerate(comPDB.molList):
        if molPDB.isPortein:
            count=0
            for resID, resPDB in enumerate(molPDB.resList):
                resPDB.resName=molPDB.proteinName
                for atmNum, atmPDB in enumerate(resPDB.atmList):
                    atmPDB.name="P"+str(count)
                    count=count+1


def assignGid(comPDB):
    for molID, molPDB in enumerate(comPDB.molList):
        #lastResID = len(molPDB.resList) - 1
        if molPDB.isPortein:
            count=0
            for resID, resPDB in enumerate(molPDB.resList):
                grpID = 0
                resID = 0
                for atmNum, atmPDB in enumerate(resPDB.atmList):
                    atmPDB.gid = (molID << 32) + (resID << 16) + (grpID << 8) + count #for RAS atmNum exceed 8 bits.
                    count = count + 1
        else:
            for resID, resPDB in enumerate(molPDB.resList):
                grpID = 0
                for atmNum, atmPDB in enumerate(resPDB.atmList):
                    atmPDB.gid = (molID << 32) + (resID << 16) + (grpID << 8) + atmNum #for RAS atmNum exceed 8 bits.

def getUniqueConsAtomList(frefile):

    consAtomDict={}

    with open(frefile, "r") as f:
        for line in f:
            strs=line.split('{')
            resName=strs[0].strip()
            subStrs=strs[1].split('}')
            nameLine=subStrs[0].strip()
            names=nameLine.split()
            consAtomDict[resName]=names

    return consAtomDict


def toObj(args, comPDB, vList):
    totAtmNum=comPDB.getTotAtmNum()

    if args.grofile:
        if totAtmNum !=len(vList):
            raise Exception("toObj - Total number of atoms doesn't match the number of velocity")

    if args.frefile:
        consAtomDict=getUniqueConsAtomList(args.frefile)

    filename=args.objfile
    outFh=open(filename, "w")
    header=getHeader(totAtmNum, args)
    outFh.write(header)

    xLHalf = args.x / 2
    yLHalf = args.y / 2
    zLHalf = args.z / 2

    zeroV=0.0
    count=0

    for molID, molPDB in enumerate(comPDB.molList):
        for resID, resPDB in enumerate(molPDB.resList):
            for atmNum, atmPDB in enumerate(resPDB.atmList):
                coor = atmPDB.coor
                speciename = resPDB.resName + "x" + atmPDB.name

                x = coor.x
                y = coor.y
                z = coor.z

                if args.origin:
                    x = x - xLHalf
                    y = y - yLHalf
                    z = z - zLHalf

                group="group"
                if args.frefile:
                    if resPDB.resName in consAtomDict:
                        atomNames=consAtomDict[resPDB.resName]
                        if atmPDB.name in atomNames:
                            group = "free"

                if args.grofile:
                    vel=vList[count]
                    count=count+1
                    outLine = "%14d ATOM %11s %s %21.13e %21.13e %21.13e %21.13e %21.13e %21.13e\n" \
                          % (atmPDB.gid, speciename, group, x, y, z, vel.vx, vel.vy, vel.vz)
                else:
                    outLine = "%14d ATOM %11s %s %21.13e %21.13e %21.13e %21.13e %21.13e %21.13e\n" \
                          % (atmPDB.gid, speciename, group, x, y, z, zeroV, zeroV, zeroV)

                outFh.write(outLine)


def getHeader(totAtmNum, args):
    header="particle FILEHEADER {type=MULTILINE; datatype=VARRECORDASCII; checksum=NONE;\n"
    header=header+"loop=0; time=0.000000;\n"
    header=header+"nfiles=1; nrecord="+str(totAtmNum)+"; nfields=10;\n"
    header=header+"field_names=id class type group rx ry rz vx vy vz;\n"
    header=header+"field_types=u s s s f f f f f f;\n"
    header=header+"h=     "+ str(args.x)+"      0.00000000000000      0.00000000000000\n"
    header=header+"       0.00000000000000     "+ str(args.y)+"       0.00000000000000\n"
    header=header+"       0.00000000000000      0.00000000000000     "+ str(args.z)+" ;\n"
    header=header+"groups = group ;\n"
    header=header+"types = ATOM ;\n"
    header=header+"} \n\n"

    return header

def getRestart(args, comPDB):
    totAtmNum = comPDB.getTotAtmNum()

    line="simulate SIMULATE { loop=0; time=0.000000 ;}\n"
    line=line+"box BOX {\n"
    line=line+"h=     "+ str(args.x)+"         0.00000000000000      0.00000000000000\n"
    line = line + "       0.00000000000000      "+ str(args.y)+"         0.00000000000000\n"
    line = line + "       0.00000000000000      0.00000000000000     "+ str(args.z)+"  ;\n"
    line=line+"}\n"
    line=line+"collection COLLECTION { mode=VARRECORDASCII; size="+str(totAtmNum)+"; files=snapshot.mem/atoms#;}\n"

    outFh = open("restart", "w")
    outFh.write(line)

def main():
    run_converter(sys.argv[1:])

def run_converter(in_options):

    args=getArgs(in_options)
    print ("Default inputs: ", args.pdbfile, args.grofile)

    itpList=[]

    if args.profile != None:
        with open(args.profile, "r") as f:
            for line in f:
                tipFileName=line.rstrip("\n\r")
                itp = ITP(tipFileName)
                itpList.append(itp)

    # fix the atom name in the
    for itp in itpList:
        count = 1
        for atom in itp.header.moleculetype.atoms.data:
            atom['atomname'] = 'P' + str(count)
            count = count + 1

    print( "Reading in pdb file ", args.pdbfile)
    comPDB=Pdb.ComPDB()
    comPDB.parse(args)
    renameProt(comPDB)
    assignGid(comPDB)
    gro = Gro.GRO()
    if args.grofile:
        gro.parse(args.grofile)
    toObj(args, comPDB, gro.vList)
    getRestart(args, comPDB)



if __name__ == '__main__':
    main()
