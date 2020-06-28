__author__ = 'zhang30'

# CHARMM energy terms are 1 based and GROMACS' are 1/2 based. k(r-r0)^2 v.s. 1/2 k(r-r0)^2

import argparse
import math
import numpy as np

from ddcmdconverter.martini.ITP import ITP
from ddcmdconverter.martini.MartiniFF import MartiniFF
from ddcmdconverter.martini.MartiniInput import martini_input


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
    parser.add_argument('-c', '--c2b', action='store', dest='cons2bond', default=False,
                        help='Add constraints to the bondList(default=True).')

    args = parser.parse_args()

    return args

def uniqueConsAtom(itp):

    resName=itp.header.moleculetype.data['name']

    atomList=[]
    for atom in itp.header.moleculetype.atoms.data:
        atomName = atom['atomname'].decode('UTF-8')
        atomList.append(atomName)

    consData=itp.header.moleculetype.constraints.data

    atomDict={}
    for i in range(len(consData)):
        constraint = consData[i]
        atomI = constraint['ai'] - 1   # 0 based
        atomJ = constraint['aj'] - 1  # 0 based
        atomDict[atomI] = 1
        atomDict[atomJ] = 1

    atomIndexes=atomDict.keys()
    #print resName, atomIndexes

    uniqueNameList=resName+" { "

    for i in atomIndexes:
        name=atomList[i]
        uniqueNameList=uniqueNameList+name+" "

    uniqueNameList = uniqueNameList + "}\n"

    return uniqueNameList


def clusterConstraint(data):
    consCluster=[]
    consList=[]

    for i in range(len(data)):
        addcons = False
        constraint = data[i]
        atomI = constraint['ai']   # 0 based
        atomJ = constraint['aj']   # 0 based
        for consTemp in consList:
            atomItemp = consTemp['ai']   # 0 based
            atomJtemp = consTemp['aj']   # 0 based
            if atomItemp==atomI:
                consList.append(constraint)
                addcons = True
                break
            if atomItemp==atomJ:
                consList.append(constraint)
                addcons = True
                break
            if atomJtemp==atomI:
                consList.append(constraint)
                addcons = True
                break
            if atomJtemp==atomJ:
                consList.append(constraint)
                addcons = True
                break

        if not addcons:
            consList = []
            consList.append(constraint)
            consCluster.append(consList)


    clustering=True

    while clustering:
        clustering = False
        clusterSize=len(consCluster)
        firstSave=0
        secondSave=0
        for first in range(clusterSize):
            for second in range(first + 1, clusterSize):
                consListFirst  = consCluster[first]
                consListSecond = consCluster[second]
                for consFirst in consListFirst:
                    for consSecond in consListSecond:
                        consFirstI = consFirst['ai']  # 0 based
                        consFirstJ = consFirst['aj']  # 0 based
                        consSecondI = consSecond['ai']  # 0 based
                        consSecondJ = consSecond['aj']  # 0 based
                        if consFirstI==consSecondI or consFirstI==consSecondJ or consFirstJ==consSecondI or consFirstJ==consSecondJ :
                            firstSave=first
                            secondSave=second
                            clustering = True
                            break
                    if clustering:
                        break
                if clustering:
                    break
            if clustering:
                break

        if clustering:
            consListFirst = consCluster[firstSave]
            consListSecond = consCluster[secondSave]
            for cons in consListSecond:
                consListFirst.append(cons)
            del consCluster[secondSave]

    return consCluster

def rmDuplicateBonds(itp):
    bondDict = {}
    new_constraints=None
    new_bonds=None

    sectionKeys = itp.header.moleculetype.sections.keys()
    if 'constraints' in sectionKeys:
        constraintSize = len(itp.header.moleculetype.constraints.data)
        constraintRMlist = []
        for i in range(constraintSize):
            constraint = itp.header.moleculetype.constraints.data[i]
            atomI = constraint['ai']  # 0 based
            atomJ = constraint['aj']  # 0 based

            if (atomI, atomJ) in bondDict.keys() or (atomJ, atomI) in bondDict.keys():
                if (atomI, atomJ) in bondDict.keys():
                    bondDict[(atomI, atomJ)] = bondDict[(atomI, atomJ)] + 1
                if (atomJ, atomI) in bondDict.keys():
                    bondDict[(atomJ, atomI)] = bondDict[(atomJ, atomI)] + 1
                constraintRMlist.append(i)
                print("Warning: Remove constraint ", constraint)
            else:
                bondDict[(atomI, atomJ)] = 1

        print("Remove constraints: ", len(constraintRMlist))
        new_constraints = np.delete(itp.header.moleculetype.constraints.data, constraintRMlist)

    if 'bonds' in sectionKeys:

        bondSize = len(itp.header.moleculetype.bonds.data)
        bondRMlist = []
        for i in range(bondSize):
            bond = itp.header.moleculetype.bonds.data[i]
            atomI = bond['ai']  # 0 based
            atomJ = bond['aj']  # 0 based

            if (atomI, atomJ) in bondDict.keys() or (atomJ, atomI) in bondDict.keys():
                if (atomI, atomJ) in bondDict.keys():
                    bondDict[(atomI, atomJ)] = bondDict[(atomI, atomJ)] + 1
                if (atomJ, atomI) in bondDict.keys():
                    bondDict[(atomJ, atomI)] = bondDict[(atomJ, atomI)] + 1
                bondRMlist.append(i)
                print("Warning: Remove bond ", bond)
            else:
                bondDict[(atomI, atomJ)] = 1

        # print(bondRMlist)
        # bondRMlist.sort(reverse=True)
        # print(bondRMlist)
        print("Original total bonds: ", len(itp.header.moleculetype.bonds.data))
        print("Remove bonds: ", len(bondRMlist))

        new_bonds = np.delete(itp.header.moleculetype.bonds.data, bondRMlist)

        print("Final total number of bonds: ", len(new_bonds))

    return new_constraints, new_bonds

def main():
    args=getArgs()
    print( "Default inputs: ", args.itpfile,  args.objfile)

    # Interaction Matrix from JPC B. Vol 111, No 27, 2007
    martiniFF=MartiniFF()

    # Mass
    par_itp = ITP(args.parfile)

    massDict = {}
    for atomtype in par_itp.header.atomtypes.data:
        massDict[atomtype['name'].decode('UTF-8')]=atomtype['mass']

    # constant
    DEG2RAD=math.pi/180

    itpList=[]
    proteinList=[]

    if args.profile != None:
        with open(args.profile, "r") as f:
            for line in f:
                tipFileName=line.rstrip("\n\r")
                itp = ITP(tipFileName)
                itpList.append(itp)
                resName = itp.header.moleculetype.data['name']
                proteinList.append(resName)

    # fix the atom name in the
    for itp in itpList:
        count = 0
        for atom in itp.header.moleculetype.atoms.data:
            atom['atomname'] = 'P' + str(count)
            count = count + 1

        #print(itp.header.moleculetype.atoms)

    if args.itpfile != None:
        with open(args.itpfile, "r") as f:
            for line in f:
                tipFileName=line.rstrip("\n\r")
                itp = ITP(tipFileName)
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
            atomType=atom['atomtype'].decode('UTF-8')
            atmTypeDict[atomType]=1   # try to find all atom type

    # Print out atom type
    atmTypeList=[x for x in atmTypeDict.keys()]

    listLine="atomTypeList="
    for atomType in atmTypeList:
        listLine=listLine+atomType+" "

    listLine=listLine+";\n"

    nbLine=""
    listLine=listLine+"ljParms="
    for nonbond in par_itp.header.nonbond_params.data:
        typeI=nonbond['i'].decode('UTF-8')
        typeJ=nonbond['j'].decode('UTF-8')
        if typeI in atmTypeList:
            if typeJ in atmTypeList:
                ljname=typeI+"_"+ typeJ
                listLine=listLine+ljname+" "
                (eps, sigma)=martiniFF.getLJparm(typeI, typeJ)
                if eps == None:
                    print("Warning: Not LJ paramters for ", typeI, typeJ)
                nbLine=nbLine+ljname+" LJPARMS{atomtypeI="+typeI+"; indexI="+str(atmTypeList.index(typeI))+"; atomtypeJ="+typeJ+"; indexJ="+str(atmTypeList.index(typeJ))+"; sigma="+str(sigma)+" nm; eps="+str(eps)+" kJ*mol^-1;}\n"

    listLine=listLine+";\n"+"}\n\n";
    fh.write(listLine);
    nbLine = nbLine +"\n"

    #MASSPARM
    listLine=""
    for index, atomType in enumerate(atmTypeList):
        try:
            listLine=listLine+atomType+" MASSPARMS { atomType="+atomType+"; atomTypeID="+str(index)+"; mass="+str(massDict[atomType])+"M_p ; }\n"
        except:
            print("atomType="+atomType)
            print("atomTypeID="+str(index))
            print("mass="+str(massDict[atomType]))

    listLine = listLine+"\n"
    fh.write(listLine)

    consAtomFh = open("ConsAtom.data", "w")
    resID=1
    for itp in itpList:
        print("###########  ", itp.header.moleculetype.data["name"], "  ###########",)
        cons2bond=args.cons2bond   # TODO: also print out the constraint pair as bond - hacking the code - should be a better solution
        #['atoms', 'bonds', 'constraints', 'angles', 'dihedrals']
        sectionKeys=itp.header.moleculetype.sections.keys()
        # itp must have atoms, but not necessory bonds, angles, dihedrals or constraints
        hasBond=False
        hasAngle=False
        hasDihedral=False
        hasConstraints=False
        hasExclusions=False
        hasRestraints = False
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
        if 'position_restraints' in sectionKeys:
            hasRestraints=True

        #Remove duplicate constraints and bonds
        newConstraints, newBonds=rmDuplicateBonds(itp)

        if hasConstraints:
            consCluster=clusterConstraint(newConstraints)
            uniqueConsNameList=uniqueConsAtom(itp)
            consAtomFh.write(uniqueConsNameList)

        line = ""
        # RESIPARMS
        resName=itp.header.moleculetype.data['name']
        line = line + resName+" RESIPARMS\n{\n"
        line = line + "  resID="+str(resID)+";\n"
        resID=resID+1
        line = line + "  resType=0;\n"
        line = line + "  resName="+resName+";\n"
        if resName in proteinList:
            line = line + "  isProtein=1;\n"
        else:
            line = line + "  isProtein=0;\n"
        numAtoms=len(itp.header.moleculetype.atoms.data)
        line = line + "  numAtoms=" + str(numAtoms) + ";\n"
        resCharge=sum(itp.header.moleculetype.atoms.data.charge)
        line = line + "  charge=" + str(resCharge) + ";\n"
        line = line + "  groupList=" + resName + "_g0;\n"
        if hasBond:
            bondSize=len(newBonds)
            totalBondSize = bondSize
            if cons2bond and hasConstraints:  # TODO: also print out the constraint pair as bond - hacking the code - should be a better solution
                constraintSize = len(newConstraints)
                totalBondSize=bondSize+constraintSize
            line = line + "  bondList="
            for i in range(totalBondSize):
                line = line + resName + "_b" + str(i)+" "
            line = line + ";\n"
        if hasConstraints:
            constraintSize = len(consCluster)
            line = line + "  constraintList="
            for i in range(constraintSize):
                line = line + resName + "_cl" + str(i)+" "
            line = line + ";\n"
        if hasExclusions:
            exclusionSize = len(itp.header.moleculetype.exclusions.data)
            line = line + "  exclusionList="
            for i in range(exclusionSize):
                line = line + resName + "_e" + str(i)+" "
            line = line + ";\n"
        if hasRestraints:
            restraintSize = len(itp.header.moleculetype.position_restraints.data)
            line = line + "  restraintList="
            for i in range(restraintSize):
                line = line + resName + "_r" + str(i)+" "
            line = line + ";\n"
        if hasAngle:
            angleSize = len(itp.header.moleculetype.angles.data)
            line = line + "  angleList="
            for i in range(angleSize):
                line = line + resName + "_a" + str(i)+" "
            line = line + ";\n"
        if hasDihedral:
            try:
                dihedralSize = len(itp.header.moleculetype.dihedrals.data)
            except:
                print(itp.filename())
            line = line + "  dihedralList="
            for i in range(dihedralSize):
                line = line + resName + "_d" + str(i)+" "
            line = line + ";\n"
        if resName=="KRAS":
            line = line + "  centerAtom=160;\n"
        else:
            line = line + "  centerAtom=0;\n"

        line = line + "}\n\n"

        # GROUPPARMS
        line = line + resName + "_g0 GROUPPARMS{\n"
        line = line + "  groupID=0;\n"
        line = line + "  atomList="
        for atom in itp.header.moleculetype.atoms.data:
            line = line + resName + "_" + atom['atomname'].decode('UTF-8') +" "
        line = line + ";\n"
        line = line + "}\n\n"

        #ATOMPARMS
        atomTypeList=[]
        for atom in itp.header.moleculetype.atoms.data:
            atomName = atom['atomname'].decode('UTF-8')
            atomID = atom['atomnr'] - 1 # 0 based
            atomType=atom['atomtype'].decode('UTF-8')
            atomTypeList.append(atomType)
            atomCharge=atom['charge']
            line = line + resName + "_" + atomName +" ATOMPARMS{"
            line = line + "atomID="+str(atomID)+"; atomName="+atomName+"; atomType="+atomType+"; atomTypeID="+str(atmTypeList.index(atomType))+"; charge="+str(atomCharge)+"; mass="+str(massDict[atomType]) +" M_p ; }\n"
        line = line +"\n"

        #BONDPARMS
        if hasBond:
            for i in range(bondSize):
                bond=newBonds[i]
                atomI = bond['ai']-1 # 0 based
                atomJ = bond['aj']-1 # 0 based
                atomTypeI = atomTypeList[atomI]
                atomTypeJ = atomTypeList[atomJ]
                func=bond['func']
                b0=bond['b0']
                kb=bond['kb']
                if b0<0.63:
                    func=1
                line = line + resName + "_b" + str(i) + " BONDPARMS{"
                line = line + "atomI="+str(atomI)+"; atomTypeI="+atomTypeI+"; atomJ="+str(atomJ) + "; atomTypeJ="+atomTypeJ \
                       + "; func=" + str(func)+"; kb="+str(0.5*kb)+" kJ*mol^-1*nm^-2; b0="+str(b0)+" nm; }\n"   # 1/2 k -> K
            line = line + "\n"
            count=bondSize
            if cons2bond and hasConstraints:  # TODO: also print out the constraint pair as bond - hacking the code - should be a better solution
                constraintSize = len(consCluster)
                for i in range(constraintSize):
                    consList = consCluster[i]
                    for j in range(len(consList)):
                        constraint = consList[j]
                        atomI = constraint['ai'] - 1  # 0 based
                        atomJ = constraint['aj'] - 1  # 0 based
                        atomTypeI = atomTypeList[atomI]
                        atomTypeJ = atomTypeList[atomJ]
                        b0 = constraint['r0']
                        line = line + resName + "_b" + str(count) + " BONDPARMS{"
                        line = line + "atomI=" + str(atomI) + "; atomTypeI=" + atomTypeI + "; atomJ=" + str(atomJ) + "; atomTypeJ=" + atomTypeJ \
                               + "; func=1; kb=10000 kJ*mol^-1*nm^-2; b0=" + str(b0) + " nm; }\n"  # 1/2 k -> K
                        count=count+1
                    line = line + "\n"

        # CONSPARMS
        if hasConstraints:
            constraintSize = len(consCluster)
            for i in range(constraintSize):
                line = line + resName + "_cl" + str(i) + " CONSLISTPARMS{ constraintSubList="
                consList=consCluster[i]
                for j in range(len(consList)):
                    line = line + resName + "_cl" + str(i) +"_c"+str(j)+" "
                line = line + "; }\n"

                for j in range(len(consList)):
                    constraint=consList[j]
                    atomI = constraint['ai']-1 # 0 based
                    atomJ = constraint['aj']-1 # 0 based
                    atomTypeI = atomTypeList[atomI]
                    atomTypeJ = atomTypeList[atomJ]
                    r0=constraint['r0']
                    func = constraint['func']
                    line = line + resName + "_cl" + str(i) +"_c"+str(j) + " CONSPARMS{"
                    line = line + "atomI="+str(atomI)+"; atomTypeI="+atomTypeI+"; atomJ="+str(atomJ)+ "; atomTypeJ="+atomTypeJ + "; func=" + str(func) +"; r0="+str(r0)+" nm; }\n"
                line = line + "\n"


        #if hasConstraints:
        #    constraintSize = len(consCluster)
        #    for i in range(constraintSize):
        #        constraint=itp.header.moleculetype.constraints.data[i]
        #        atomI = constraint['ai']-1 # 0 based
        #        atomJ = constraint['aj']-1 # 0 based
        #        atomTypeI = atomTypeList[atomI]
        #        atomTypeJ = atomTypeList[atomJ]
        #        r0=constraint['r0']
        #        func = constraint['func']
        #        line = line + resName + "_c" + str(i) + " CONSPARMS{"
        #        line = line + "atomI="+str(atomI)+"; atomTypeI="+atomTypeI+"; atomJ="+str(atomJ)+ "; atomTypeJ="+atomTypeJ + "; func=" + str(func) +"; r0="+str(r0)+" nm; }\n"
        #    line = line + "\n"

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

        # RESTRAINT
        if hasRestraints:
            for i in range(restraintSize):
                restraint=itp.header.moleculetype.position_restraints.data[i]
                atomI = restraint['ai']-1 # 0 based
                atomTypeI = atomTypeList[atomI]
                func=restraint['func']
                fcx = restraint['fcx']
                fcy = restraint['fcy']
                fcz = restraint['fcz']
                line = line + resName + "_r" + str(i) + " RESTRAINTPARMS{"
                line = line + "atomI="+str(atomI)+"; atomTypeI="+atomTypeI + "; func=" + str(func) \
                       +"; fcx="+str(fcx)+"; fcy="+str(fcy)+"; fcz="+str(fcz)+"; kb= 1.0 kJ*mol^-1*nm^-2; }\n"   # 1/2 k -> K
            line = line + "\n"

        # ANGLEPARMS
        if hasAngle:
            for i in range(angleSize):
                if i==418 :
                    print(itp.header.moleculetype.angles.data[i])
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
                if func==2 or func==10:
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
                ntor = dihedral['ntor']
                if not ntor:
                    ntor=1
                if func==1:
                    delta=-delta
                if func==2:
                    kchi=kchi/2.0
                line = line + resName + "_d" + str(i) + " TORSPARMS{"
                line = line + "atomI="+str(atomI)+"; atomJ="+str(atomJ)+"; atomK="+str(atomK)+"; atomL="+str(atomL) \
                       + "; func=" + str(func) +"; delta="+str(delta*DEG2RAD)+"; kchi="+str(kchi)+" kJ*mol^-1; n="+str(ntor)+";}\n"
            line = line + "\n"

        fh.write(line)


    fh.write(nbLine)

    # Print out species

    # comment out old specie input
    #sysLine="system SYSTEM\n{species =\n"
    #speLine=""
    #for itp in itpList:
    #    resName=itp.header.moleculetype.data['name']

    #    for atom in itp.header.moleculetype.atoms.data:
    #        specie=resName+"x"+atom['atomname']
    #        sysLine=sysLine+specie+"\n"
    #        atomType=atom['atomtype']
    #        speLine=speLine+specie+" SPECIES { type = ATOM ; charge ="+str(atom['charge'])+"; id="+str(atmTypeList.index(atomType))+ "; mass ="+str(massDict[atomType])+" M_p ; }\n"

    #sysLine=sysLine+"   ;\n}\n\n"


    molecLine="moleculeClass MOLECULECLASS { molecules = "
    moleSpecie="\n"
    speLine = "\n"
    for itp in itpList:
        resName = itp.header.moleculetype.data['name']
        molName=resName+"x"
        molecLine=molecLine+" "+molName

        moleSpecie =moleSpecie+molName+" MOLECULE {ownershipSpecies = "

        count=0
        for atom in itp.header.moleculetype.atoms.data:
            specie = molName + atom['atomname'].decode('UTF-8')

            atomType=atom['atomtype'].decode('UTF-8')
            speLine=speLine+specie+" SPECIES { type = ATOM ; charge ="+str(atom['charge'])+"; id="+str(atmTypeList.index(atomType))+ "; mass ="+str(massDict[atomType])+" M_p ; }\n"

            if count==0:
                moleSpecie=moleSpecie+specie+"; species = "+specie
            else:
                moleSpecie=moleSpecie+" "+specie
            count=count+1

        moleSpecie = moleSpecie + ";}\n"


    molecLine=molecLine+"; }\n"

    speFh=open(args.splfile, "w")
    speFh.write(martini_input)
    speFh.write(molecLine)
    speFh.write(moleSpecie)
    speFh.write(speLine)


if __name__ == '__main__':
    main()
