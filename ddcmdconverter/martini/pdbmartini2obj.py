__author__ = 'zhang30'


import argparse
import sys

import ddcmdconverter.base.Object as Object

def getArgs(in_options):

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pdb', action='store', dest='pdbfile', default='system_fix.pdb',
                        help='PDB input file (default=system_fix.pdb).')
    parser.add_argument('-t', '--top', action='store', dest='topfile', default='topol.top',
                        help='GROMACS top file (default=topol.top).')
    parser.add_argument('-m', '--mar', action='store', dest='marfile', default='martini.data',
                        help='ddcMD martini.data file (default=martini.data).')
    parser.add_argument('-f', '--fre', action='store', dest='frefile', default=None,
                        help='Constrant atom name for setting free in group.')
    parser.add_argument('-o', '--obj', action='store', dest='objfile', default='atom#000000.data',
                        help='ddcMD object output file (default=atom#000000.data).')

    args = parser.parse_args(in_options)

    return args

def getUniqueConsAtomList(args):

    consAtomDict={}
    contents=""
    with open(args.frefile, "r") as f:
        for line in f:
            contents=contents+line

        contentList=contents.split("}")
        for item in contentList[:-1]: # skip last one which is empty
            itemList=item.split('{')
            if len(itemList) != 2:
                raise Exception("Object curly brackets do not match: "+item+"}")

            itemList = item.replace('\n', ' ').split('{')
            name=itemList[0].strip()
            consList=itemList[1].split()
            consAtomDict[name]=consList

    return consAtomDict


def getHeader(totNumAtoms, boxSize):
    (x,y,z)=boxSize
    header="particle FILEHEADER {type=MULTILINE; datatype=VARRECORDASCII; checksum=NONE;\n"
    header=header+"loop=0; time=0.000000;\n"
    header=header+"nfiles=1; nrecord="+str(totNumAtoms)+"; nfields=10;\n"
    header=header+"field_names=id class type group rx ry rz vx vy vz;\n"
    header=header+"field_types=u s s s f f f f f f;\n"
    header=header+"h=     "+ str(x)+"      0.00000000000000      0.00000000000000\n"
    header=header+"       0.00000000000000     "+ str(y)+"       0.00000000000000\n"
    header=header+"       0.00000000000000      0.00000000000000     "+ str(z)+" ;\n"
    header=header+"groups = group ;\n"
    header=header+"types = ATOM ;\n"
    header=header+"} \n\n"

    return header

def getRestart(totNumAtoms, boxSize):
    (x, y, z) = boxSize

    line="simulate SIMULATE { loop=0; time=0.000000 ;}\n"
    line=line+"box BOX {\n"
    line=line+"h=     "+ str(x)+"         0.00000000000000      0.00000000000000\n"
    line = line + "       0.00000000000000      "+ str(y)+"         0.00000000000000\n"
    line = line + "       0.00000000000000      0.00000000000000     "+ str(z)+"  ;\n"
    line=line+"}\n"
    line=line+"collection COLLECTION { mode=VARRECORDASCII; size="+str(totNumAtoms)+"; files=snapshot.mem/atoms#;}\n"

    outFh = open("restart", "w")
    outFh.write(line)

def martiniDatafile(args):
    objectList = Object.ObjectList(args.marfile)
    martiniObject=objectList.getObject('martini')
    martiniAttrNames=martiniObject.getAttrNames()
    #print(martiniAttrNames)
    resiParms=martiniObject.getAttr('resiParms')
    print("Molecules in martini.data: "+resiParms)
    resList=resiParms.split()
    resDict={}
    for res in resList:
        resObject=objectList.getObject(res)
        isProteinStr = resObject.getAttr('isProtein')
        isProtein=int(isProteinStr)
        numAtomsStr = resObject.getAttr('numAtoms')
        numAtoms=int(numAtomsStr)
        groupListStr = resObject.getAttr('groupList')
        groupList=groupListStr.split()

        atomListStr=""
        for grp in groupList:
            grpObject = objectList.getObject(grp)
            atomListStr=atomListStr+grpObject.getAttr('atomList')+" "

        atomList=atomListStr.split()
        if len(atomList) != numAtoms:
            raise Exception("Molecule "+res+" size of atomList is not equal to numAtoms "+numAtomsStr)


        resDict[res]=(isProtein, numAtoms, atomList)

    return resDict

def groTopfile(args):
    contents=""
    with open(args.topfile, 'r') as f:
        for line in f:
            if line[0]==';' or line[0]=='#':
                continue
            contents=contents+line
    sections=contents.split('[')

    molTopList=[]
    for secRaw in sections[1:]:
        sec=secRaw.split(']')
        secName=sec[0].strip()
        if secName=='molecules':
            itemList=sec[1].strip().split('\n')
            for item in itemList:
                molinfo=item.split()
                if len(molinfo)>=2:
                    molTopList.append((molinfo[0],int(molinfo[1])))
                else:
                    print("Warning Gromacs Top file molecules wrong format - " + item)
    return molTopList

def groTopfileNumAtoms(molTopList, resDict):
    totNumAtoms=0
    for (name, numMol) in molTopList:
        (isProtein, numAtoms, atomList)=resDict[name]
        totNumAtoms=totNumAtoms+numMol*numAtoms
    return totNumAtoms

def pdbfile(args):
    pdbatomList=[]
    with open(args.pdbfile, 'r') as f:
        for line in f:
            if line[:4] == "ATOM" or line[:6]=="HETATM":
                pdbatomList.append(line)
    return pdbatomList

def pdbBoxSize(args):
    with open(args.pdbfile, 'r') as f:
        for line in f:
            if line[:6] == "CRYST1":
                strs=line.split()
                if len(strs)>4:
                    x = float(strs[1])
                    y = float(strs[2])
                    z = float(strs[3])
                    return (x,y,z)
    return (0,0,0)

def main():
    run_converter(sys.argv[1:])

def run_converter(in_options):

    args=getArgs(in_options)
    print ("Default inputs: ", args.pdbfile, args.topfile, args.marfile, args.objfile)

    resDict=martiniDatafile(args)
    molTopList=groTopfile(args)
    consAtomDict=getUniqueConsAtomList(args)
    pdbatomList=pdbfile(args)

    # number of atoms in GROMACS top file
    totNumAtoms=groTopfileNumAtoms(molTopList, resDict)
    if totNumAtoms != len(pdbatomList):
        raise Exception('Number of atoms in GROMACS top file does NOT match that in PDB file')

    boxSize=pdbBoxSize(args)
    if boxSize[1] == 0 :
        raise Exception('Box size in PDB file is NOT correct')

    header=getHeader(totNumAtoms, boxSize)

    with open(args.objfile, 'w') as f:
        f.write(header)
        count=-1
        molID=-1
        zeroV=0.0
        for (name, numMol) in molTopList:
            (isProtein, numAtoms, atomList)=resDict[name]
            consAtomList=None
            if name in consAtomDict:
                consAtomList=consAtomDict[name]
            for mol in range(numMol):
                molID=molID+1
                for atmID, atom in enumerate(atomList):
                    # figure out the resname and atomnem in atomList from martini.data
                    molnameSize=len(name)
                    atomNameMartini=atom[molnameSize+1:]
                    # found out information in pdb
                    count=count+1
                    pdbatomline=pdbatomList[count]
                    resName = pdbatomline[17:21].strip()
                    atmName = pdbatomline[12:17].strip()
                    resID =   pdbatomline[22:27].strip()
                    x = float(pdbatomline[30:38])
                    y = float(pdbatomline[38:46])
                    z = float(pdbatomline[46:54])

                    if isProtein == 0:
                        if atomNameMartini != atmName:
                            raise Exception("Atom Name in martini.data ["+atomNameMartini+
                                            "] and pdb file ["+atmName+"] are different\n for residue "
                                            +resName+" residue ID "+ resID)
                        if name != resName:
                            raise Exception("Residue Name in martini.data [" + name +
                                            "] and pdb file [" + resName + "] are different\n for residue ID " + resID)

                    speice=name+"x"+atomNameMartini
                    group = "group"
                    if consAtomList and atomNameMartini in consAtomList:
                        group = "free"
                    gid = (molID << 32) + atmID
                    outLine = "%16d ATOM %15s %6s %21.13e %21.13e %21.13e %21.13e %21.13e %21.13e\n" \
                              % (gid, speice, group, x, y, z, zeroV, zeroV, zeroV)

                    f.write(outLine)


    getRestart(totNumAtoms, boxSize)



if __name__ == '__main__':
    main()
