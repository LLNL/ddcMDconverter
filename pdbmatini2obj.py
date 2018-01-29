__author__ = 'zhang30'


import argparse

import Pdb
import Obj
import Specie


def getArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pdb', action='store', dest='pdbfile', default='test.pdb', help='PDB input file (default=test.pdb).')
    parser.add_argument('-o', '--obj', action='store', dest='objfile', default='atom#.data',
                        help='ddcMD object output file (default=atom#.data).')

    args = parser.parse_args()

    return args

def assignGid(comPDB):
    for molID, molPDB in enumerate(comPDB.molList):
        lastResID = len(molPDB.resList) - 1
        for resID, resPDB in enumerate(molPDB.resList):

            grpID = 0
            for atmNum, atmPDB in enumerate(resPDB.atmList):
                atmPDB.gid = (molID << 32) + (resID << 16) + (grpID << 8) + atmNum #for RAS atmNum exceed 8 bits.

def toObj(args, comPDB):
    totAtmNum=comPDB.getTotAtmNum()

    filename=args.objfile
    outFh=open(filename, "w")
    header=getHeader(totAtmNum, args)
    outFh.write(header)

    zeroV=0.0

    for molID, molPDB in enumerate(comPDB.molList):
        for resID, resPDB in enumerate(molPDB.resList):
            for atmNum, atmPDB in enumerate(resPDB.atmList):
                coor = atmPDB.coor
                speciename = resPDB.resName + "x" + atmPDB.name
                outLine = "%14d ATOM %11s group %21.13e %21.13e %21.13e %21.13e %21.13e %21.13e\n" \
                          % (atmPDB.gid, speciename, coor.x, coor.y, coor.z, zeroV, zeroV, zeroV)
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



if __name__ == '__main__':

    args=getArgs()
    print "Default inputs: ", args.pdbfile


    print "Reading in pdb file ", args.pdbfile
    comPDB=Pdb.ComPDB()
    comPDB.parse(args)
    assignGid(comPDB)
    toObj(args, comPDB)



