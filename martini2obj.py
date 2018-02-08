__author__ = 'zhang30'

# CHARMM energy terms are 1 based and GROMACS' are 1/2 based. k(r-r0)^2 v.s. 1/2 k(r-r0)^2

import argparse
import ITP
import math
import MartiniFF


def getArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--itp', action='store', dest='itpfile', default=None, help='ITP file list (default=itpList).')
    parser.add_argument('-t', '--pro', action='store', dest='profile', default=None,
                        help='Protein ITP file list (default=proItpList).')
    parser.add_argument('-p', '--par', action='store', dest='parfile', default='martini_v2.2.itp',
                        help='Overall ITP file (default=martini_v2.2.itp).')
    parser.add_argument('-o', '--obj', action='store', dest='objfile', default='martini.data', help='Martini object output file (default=martini.data).')
    parser.add_argument('-l', '--spl', action='store', dest='splfile', default='speless.data',
                        help='ddcMD species less output file (default=speless.data).')

    args = parser.parse_args()

    return args

if __name__ == '__main__':

    args=getArgs()
    print "Default inputs: ", args.itpfile,  args.objfile

    # Interaction Matrix from JPC B. Vol 111, No 27, 2007
    martiniFF=MartiniFF.MartiniFF()

    # Mass
    par_itp = ITP.ITP(args.parfile)

    massDict = {}
    for atomtype in par_itp.header.atomtypes.data:
        massDict[atomtype['name']]=atomtype['mass']

    # constant
    DEG2RAD=math.pi/180

    itpList=[]

    if args.profile != None:
        with open(args.profile, "r") as f:
            for line in f:
                tipFileName=line.rstrip("\n\r")
                itp = ITP.ITP(tipFileName)
                itpList.append(itp)

    # fix the atom name in the
    for itp in itpList:
        count = 0
        for atom in itp.header.moleculetype.atoms.data:
            atom['atomname'] = 'P' + str(count)
            count = count + 1

        #print itp.header.moleculetype.atoms

    if args.itpfile != None:
        with open(args.itpfile, "r") as f:
            for line in f:
                tipFileName=line.rstrip("\n\r")
                itp = ITP.ITP(tipFileName)
                itpList.append(itp)

    """
    martini
    MMFF
    {
        resiParms = CHOL
    POPC;
    }    
    """

    fh=open(args.objfile, "w")

    fh.write("martini MMFF\n{\n");

    fh.write("resiParms=");
    for itp in itpList:
        fh.write(itp.header.moleculetype.data['name'])
        fh.write(" ");

    fh.write(";\n");

    atmTypeDict={}

    for itp in itpList:
        for atom in itp.header.moleculetype.atoms.data:
            atomType=atom['atomtype']
            atmTypeDict[atomType]=1   # try to find all atom type

    # Print out atom type
    atmTypeList=atmTypeDict.keys()

    listLine="atomTypeList="
    for atomType in atmTypeList:
        listLine=listLine+atomType+" "

    listLine=listLine+";\n"

    nbLine=""
    listLine=listLine+"ljParms="
    for nonbond in par_itp.header.nonbond_params.data:
        typeI=nonbond['i']
        typeJ=nonbond['j']
        if typeI in atmTypeList:
            if typeJ in atmTypeList:
                ljname=typeI+"_"+ typeJ
                listLine=listLine+ljname+" "
                (eps, sigma)=martiniFF.getLJparm(typeI, typeJ)
                nbLine=nbLine+ljname+" LJPARMS{atomtypeI="+typeI+"; indexI="+str(atmTypeList.index(typeI))+"; atomtypeJ="+typeJ+"; indexJ="+str(atmTypeList.index(typeJ))+"; sigma="+str(sigma)+" nm; eps="+str(eps)+" kJ*mol^-1;}\n"

    listLine=listLine+";\n"+"}\n\n";
    fh.write(listLine);
    nbLine = nbLine +"\n"

    #MASSPARM
    listLine=""
    for index, atomType in enumerate(atmTypeList):
        listLine=listLine+atomType+" MASSPARMS { atomType="+atomType+"; atomTypeID="+str(index)+"; mass="+str(massDict[atomType])+"M_p ; }\n"

    listLine = listLine+"\n"
    fh.write(listLine)

    resID=1
    for itp in itpList:
        #['atoms', 'bonds', 'constraints', 'angles', 'dihedrals']
        sectionKeys=itp.header.moleculetype.sections.keys()
        # itp must have atoms, but not necessory bonds, angles, dihedrals or constraints
        hasBond=False
        hasAngle=False
        hasDihedral=False
        hasConstraints=False
        hasExclusions=False
        if 'bonds' in sectionKeys:
            hasBond=True
        if 'angles' in sectionKeys:
            hasAngle=True
        if 'dihedrals' in sectionKeys:
            hasDihedral = True
        if 'constraints' in sectionKeys:
            hasConstraints=True
        if 'exclusions' in sectionKeys:
            hasExclusions=True

        line = ""
        # RESIPARMS
        resName=itp.header.moleculetype.data['name']
        line = line + resName+" RESIPARMS\n{\n"
        line = line + "  resID="+str(resID)+";\n"
        resID=resID+1
        line = line + "  resType=0;\n"
        line = line + "  resName="+resName+";\n"
        resCharge=sum(itp.header.moleculetype.atoms.data.charge)
        line = line + "  charge=" + str(resCharge) + ";\n"
        line = line + "  groupList=" + resName + "_g0;\n"
        if hasBond:
            bondSize=len(itp.header.moleculetype.bonds.data)
            line = line + "  bondList="
            for i in range(bondSize):
                line = line + resName + "_b" + str(i)+" "
            line = line + ";\n"
        if hasConstraints:
            constraintSize = len(itp.header.moleculetype.constraints.data)
            line = line + "  constraintList="
            for i in range(constraintSize):
                line = line + resName + "_c" + str(i)+" "
            line = line + ";\n"
        if hasExclusions:
            exclusionSize = len(itp.header.moleculetype.exclusions.data)
            line = line + "  exclusionList="
            for i in range(exclusionSize):
                line = line + resName + "_e" + str(i)+" "
            line = line + ";\n"
        if hasAngle:
            angleSize = len(itp.header.moleculetype.angles.data)
            line = line + "  angleList="
            for i in range(angleSize):
                line = line + resName + "_a" + str(i)+" "
            line = line + ";\n"
        if hasDihedral:
            dihedralSize = len(itp.header.moleculetype.dihedrals.data)
            line = line + "  dihedralList="
            for i in range(dihedralSize):
                line = line + resName + "_d" + str(i)+" "
            line = line + ";\n"
        line = line + "}\n\n"

        # GROUPPARMS
        line = line + resName + "_g0 GROUPPARMS{\n"
        line = line + "  groupID=0;\n"
        line = line + "  atomList="
        for atom in itp.header.moleculetype.atoms.data:
            line = line + resName + "_" + atom['atomname'] +" "
        line = line + ";\n"
        line = line + "}\n\n"

        #ATOMPARMS
        atomTypeList=[]
        for atom in itp.header.moleculetype.atoms.data:
            atomName = atom['atomname']
            atomID = atom['atomnr'] - 1 # 0 based
            atomType=atom['atomtype']
            atomTypeList.append(atomType)
            atomCharge=atom['charge']
            line = line + resName + "_" + atomName +" ATOMPARMS{"
            line = line + "atomID="+str(atomID)+"; atomName="+atomName+"; atomType="+atomType+"; atomTypeID="+str(atmTypeList.index(atomType))+"; charge="+str(atomCharge)+"; mass="+str(massDict[atomType]) +" M_p ; }\n"
        line = line +"\n"

        #BONDPARMS
        if hasBond:
            for i in range(bondSize):
                bond=itp.header.moleculetype.bonds.data[i]
                atomI = bond['ai']-1 # 0 based
                atomJ = bond['aj']-1 # 0 based
                atomTypeI = atomTypeList[atomI]
                atomTypeJ = atomTypeList[atomJ]
                func=bond['func']
                b0=bond['b0']
                kb=bond['kb']
                line = line + resName + "_b" + str(i) + " BONDPARMS{"
                line = line + "atomI="+str(atomI)+"; atomTypeI="+atomTypeI+"; atomJ="+str(atomJ) + "; atomTypeJ="+atomTypeJ \
                       + "; func=" + str(func)+"; kb="+str(0.5*kb)+" kJ*mol^-1*nm^-2; b0="+str(b0)+" nm; }\n"   # 1/2 k -> K
            line = line + "\n"

        # CONSPARMS
        if hasConstraints:
            for i in range(constraintSize):
                constraint=itp.header.moleculetype.constraints.data[i]
                atomI = constraint['ai']-1 # 0 based
                atomJ = constraint['aj']-1 # 0 based
                r0=constraint['r0']
                line = line + resName + "_c" + str(i) + " CONSPARMS{"
                line = line + "atomI="+str(atomI)+"; atomJ="+str(atomJ)+"; r0="+str(r0)+" nm; }\n"
            line = line + "\n"

        # EXCLUDEPARMS
        if hasExclusions:
            for i in range(exclusionSize):
                exclusion=itp.header.moleculetype.exclusions.data[i]
                atomI = exclusion['ai']-1 # 0 based
                atomJ = exclusion['aj']-1 # 0 based
                atomTypeI = atomTypeList[atomI]
                atomTypeJ = atomTypeList[atomJ]
                line = line + resName + "_e" + str(i) + " EXCLUDEPARMS{"
                line = line + "atomI="+str(atomI)+"; atomTypeI="+atomTypeI+"; atomJ="+str(atomJ)+ "; atomTypeJ="+atomTypeJ+"; }\n"
            line = line + "\n"

        # ANGLEPARMS
        if hasAngle:
            for i in range(angleSize):
                angle=itp.header.moleculetype.angles.data[i]
                atomI = angle['ai']-1 # 0 based
                atomJ = angle['aj']-1 # 0 based
                atomK = angle['ak'] - 1  # 0 based
                theta0=angle['theta0']
                ktheta=angle['cth']
                func = angle['func']
                line = line + resName + "_a" + str(i) + " ANGLEPARMS{"
                if func==1:
                    line = line + "atomI="+str(atomI)+"; atomJ="+str(atomJ)+"; atomK="+str(atomK) + "; func=" + str(func)\
                           +"; theta0="+str(theta0*DEG2RAD)+"; ktheta="+str(0.5*ktheta)+" kJ*mol^-1; }\n"
                if func==2:
                    line = line + "atomI="+str(atomI)+"; atomJ="+str(atomJ)+"; atomK="+str(atomK) + "; func=" + str(func)\
                           +"; theta0="+str(math.cos(theta0*DEG2RAD))+"; ktheta="+str(0.5*ktheta)+" kJ*mol^-1; }\n" # cosine based
            line = line + "\n"
        # TORSPARMS
        if hasDihedral:
            for i in range(dihedralSize):
                dihedral=itp.header.moleculetype.dihedrals.data[i]
                atomI = dihedral['ai']-1 # 0 based
                atomJ = dihedral['aj']-1 # 0 based
                atomK = dihedral['ak'] - 1  # 0 based
                atomL = dihedral['al'] - 1  # 0 based
                delta=dihedral['delta']
                kchi=dihedral['kchi']
                func=dihedral['func']
                if func==1:
                    delta=-delta
                if func==2:
                    kchi=kchi/2.0
                line = line + resName + "_d" + str(i) + " TORSPARMS{"
                line = line + "atomI="+str(atomI)+"; atomJ="+str(atomJ)+"; atomK="+str(atomK)+"; atomL="+str(atomL) \
                       + "; func=" + str(func) +"; delta="+str(delta*DEG2RAD)+"; kchi="+str(kchi)+" kJ*mol^-1; n=1;}\n"
            line = line + "\n"

        fh.write(line)


    fh.write(nbLine)

    # Print out species

    sysLine="system SYSTEM\n{species =\n"
    speLine=""
    for itp in itpList:
        resName=itp.header.moleculetype.data['name']

        for atom in itp.header.moleculetype.atoms.data:
            specie=resName+"x"+atom['atomname']
            sysLine=sysLine+specie+"\n"
            atomType=atom['atomtype']
            speLine=speLine+specie+" SPECIES { type = ATOM ; charge ="+str(atom['charge'])+"; id="+str(atmTypeList.index(atomType))+ "; mass ="+str(massDict[atomType])+" M_p ; }\n"

    sysLine=sysLine+"   ;\n}\n\n"

    speFh=open(args.splfile, "w")
    speFh.write(sysLine)
    speFh.write(speLine)
