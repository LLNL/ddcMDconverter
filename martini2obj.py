__author__ = 'zhang30'


import argparse
import ITP


def getArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--top', action='store', dest='itpfile', default='itpList', help='ITP file list (default=itpList).')
    parser.add_argument('-o', '--obj', action='store', dest='objfile', default='martini.data', help='Martini object output file (default=martini.data).')

    args = parser.parse_args()

    return args

if __name__ == '__main__':

    args=getArgs()
    print "Default inputs: ", args.itpfile,  args.objfile

    itpList=[]

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

    fh.write(";\n}\n\n");

    resID=1
    for itp in itpList:
        #['atoms', 'bonds', 'constraints', 'angles', 'dihedrals']
        sectionKeys=itp.header.moleculetype.sections.keys()
        hasAngle=False
        hasDihedral=False
        hasConstraints=False
        if 'angles' in sectionKeys:
            hasAngle=True
        if 'dihedrals' in sectionKeys:
            hasDihedral = True
        if 'constraints' in sectionKeys:
            hasConstraints=True

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
        line = line + "  groupList=" + resName + "-g0;\n"
        bondSize=len(itp.header.moleculetype.bonds.data)
        line = line + "  bondList="
        for i in range(bondSize):
            line = line + resName + "-b" + str(i)+" "
        line = line + ";\n"
        if hasConstraints:
            constraintSize = len(itp.header.moleculetype.constraints.data)
            line = line + "  constraintList="
            for i in range(constraintSize):
                line = line + resName + "-c" + str(i)+" "
            line = line + ";\n"
        if hasAngle:
            angleSize = len(itp.header.moleculetype.angles.data)
            line = line + "  angleList="
            for i in range(constraintSize):
                line = line + resName + "-a" + str(i)+" "
            line = line + ";\n"
        if hasDihedral:
            dihedralSize = len(itp.header.moleculetype.dihedrals.data)
            line = line + "  dihedralList="
            for i in range(constraintSize):
                line = line + resName + "-d" + str(i)+" "
            line = line + ";\n"
        line = line + "}\n\n"

        # TGROUP_PARMS
        line = line + resName + "-g0 TGROUP_PARMS{\n"
        line = line + "  atomList="
        for atom in itp.header.moleculetype.atoms.data:
            line = line + resName + "-" + atom['atomname'] +" "
        line = line + ";\n"
        line = line + "}\n\n"

        #TATOM_PARMS
        for atom in itp.header.moleculetype.atoms.data:
            atomName = atom['atomname']
            atomID = atom['atomnr'] - 1 # 0 based
            atomType=atom['atomtype']
            atomCharge=atom['charge']
            line = line + resName + "-" + atomName +" TATOM_PARMS{"
            line = line + "atmID="+str(atomID)+"; atmName="+atomName+"; atmType="+atomType+"; charge="+str(atomCharge)+"; }\n"
        line = line +"\n"

        #BOND_CONN
        for i in range(bondSize):
            bond=itp.header.moleculetype.bonds.data[i]
            atomI = bond['ai']-1 # 0 based
            atomJ = bond['aj']-1 # 0 based
            b0=bond['b0']
            kb=bond['kb']
            line = line + resName + "-b" + str(i) + " BOND_CONN{"
            line = line + "atmI="+str(atomI)+"; atomJ="+str(atomJ)+"; kb="+str(kb)+"; b0="+str(b0)+"; }\n"
        line = line + "\n"

        # CONSTRAINT_CONN
        if hasConstraints:
            for i in range(constraintSize):
                constraint=itp.header.moleculetype.constraints.data[i]
                atomI = constraint['ai']-1 # 0 based
                atomJ = constraint['aj']-1 # 0 based
                r0=constraint['r0']
                line = line + resName + "-c" + str(i) + " CONSTRAINT_CONN{"
                line = line + "atmI="+str(atomI)+"; atomJ="+str(atomJ)+"; r0="+str(r0)+"; }\n"
            line = line + "\n"

        # Angle_CONN
        if hasAngle:
            for i in range(angleSize):
                angle=itp.header.moleculetype.angles.data[i]
                atomI = angle['ai']-1 # 0 based
                atomJ = angle['aj']-1 # 0 based
                atomK = angle['ak'] - 1  # 0 based
                theta0=angle['theta0']
                ktheta=angle['cth']
                line = line + resName + "-a" + str(i) + " ANGLE_CONN{"
                line = line + "atmI="+str(atomI)+"; atomJ="+str(atomJ)+"; atomK="+str(atomK)+"; theta0="+str(theta0)+"; ktheta="+str(ktheta)+"; }\n"
            line = line + "\n"
        # TORS_CONN
        if hasDihedral:
            for i in range(dihedralSize):
                dihedral=itp.header.moleculetype.dihedrals.data[i]
                atomI = dihedral['ai']-1 # 0 based
                atomJ = dihedral['aj']-1 # 0 based
                atomK = dihedral['ak'] - 1  # 0 based
                atomL = dihedral['al'] - 1  # 0 based
                delta=dihedral['delta']
                kchi=dihedral['kchi']
                line = line + resName + "-d" + str(i) + " TORS_CONN{"
                line = line + "atmI="+str(atomI)+"; atomJ="+str(atomJ)+"; atomK="+str(atomK)+"; atomL="+str(atomL)+"; delta="+str(delta)+"; kchi="+str(kchi)+"; }\n"
            line = line + "\n"

        fh.write(line)