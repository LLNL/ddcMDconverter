from ddcmdconverter.martini.ITP import ITP

#itp=ITP.ITP("martini_v2.0_POPC_02.itp")
#itp=ITP.ITP("martini_v2.1-dna.itp")
#itp=ITP("KRAS-GTP-04-HVR-best-guess-M22-CYFpos.itp")
itp=ITP("TERNARY.itp")
#print(itp.header.moleculetype.data['name'])
#itp=ITP.ITP("martini_v2.0_CHOL_01.itp")
#print itp.header
print(itp.header.moleculetype)

#print itp.header.moleculetype.atoms

count=1
for atom in itp.header.moleculetype.atoms.data:
    atom['atomname']='P'+str(count)
    count=count+1

#print(itp.header.moleculetype.atoms)

#print itp.header.moleculetype.data
#print(itp.header.moleculetype.data['name'])

#print itp.header.moleculetype.sections.keys()

#print itp.header.moleculetype.exclusions.data

#print itp.header.atomtypes.data

#print itp.header.nonbond_params.data

#print sum(itp.header.moleculetype.atoms.data.charge)

#print itp.header.moleculetype.bonds.data

bondDict={}
bondSize = len(itp.header.moleculetype.bonds.data)
for i in range(bondSize):
    bond = itp.header.moleculetype.bonds.data[i]
    atomI = bond['ai']   # 0 based
    atomJ = bond['aj']   # 0 based

    if (atomI, atomJ) in bondDict.keys() or (atomJ, atomI) in bondDict.keys():
        if (atomI, atomJ) in bondDict.keys():
            bondDict[(atomI, atomJ)] = bondDict[(atomI, atomJ)] + 1
        if (atomJ, atomI) in bondDict.keys():
            bondDict[(atomJ, atomI)] = bondDict[(atomJ, atomI)] + 1
    else:
        bondDict[(atomI, atomJ)] = 1

constraintSize = len(itp.header.moleculetype.constraints.data)
for i in range(constraintSize):
    constraint = itp.header.moleculetype.constraints.data[i]
    atomI = constraint['ai']   # 0 based
    atomJ = constraint['aj']  # 0 based

    if (atomI, atomJ) in bondDict.keys() or (atomJ, atomI) in bondDict.keys():
        if (atomI, atomJ) in bondDict.keys():
            bondDict[(atomI, atomJ)] = bondDict[(atomI, atomJ)] + 1
        if (atomJ, atomI) in bondDict.keys():
            bondDict[(atomJ, atomI)] = bondDict[(atomJ, atomI)] + 1
    else:
        bondDict[(atomI, atomJ)] = 1

exclusionSize = len(itp.header.moleculetype.exclusions.data)
for i in range(exclusionSize):
    exclusion = itp.header.moleculetype.exclusions.data[i]
    atomI = exclusion['ai']   # 0 based
    atomJ = exclusion['aj']  # 0 based


for key, value in bondDict.items():
    if value> 1:
        print(key, value)


angleDict={}
angleSize = len(itp.header.moleculetype.angles.data)
for i in range(angleSize):
    angle=itp.header.moleculetype.angles.data[i]
    atomI = angle['ai'] # 0 based
    atomJ = angle['aj'] # 0 based
    atomK = angle['ak'] # 0 based

    if (atomI, atomJ, atomK) in angleDict.keys() or (atomJ, atomI, atomK) in angleDict.keys():
        if (atomI, atomJ, atomK) in angleDict.keys():
            angleDict[(atomI, atomJ, atomK)] = angleDict[(atomI, atomJ, atomK)] + 1
        if (atomJ, atomI, atomK) in angleDict.keys():
            angleDict[(atomJ, atomI, atomK)] = angleDict[(atomJ, atomI, atomK)] + 1
    else:
        angleDict[(atomI, atomJ)] = 1

    if (atomI, atomK) in bondDict.keys() or (atomK, atomI) in bondDict.keys():
        if (atomI, atomK) in bondDict.keys():
            bondDict[(atomI, atomK)] = bondDict[(atomI, atomK)] + 1
        if (atomK, atomI) in bondDict.keys():
            bondDict[(atomK, atomI)] = bondDict[(atomK, atomI)] + 1
    else:
        bondDict[(atomI, atomJ)] = 1

for key, value in bondDict.items():
    if value> 1:
        print(key, value)


for key, value in angleDict.items():
    if value> 1:
        print(key, value)

with open('duplicate', 'w') as f:
    for key, value in bondDict.items():
        if value > 1:
            f.write("("+' '.join(str(s) for s in key) +") "+ str(value)+"\n")

    for key, value in angleDict.items():
        if value > 1:
            f.write("("+' '.join(str(s) for s in key) +") "+ str(value)+"\n")


#print itp.header.moleculetype.constraints.data

#print itp.header.moleculetype.angles.data

#print itp.header.moleculetype.dihedrals.data

#for atom in itp.header.moleculetype.atoms.data:
#    print atom

#for bond in itp.header.moleculetype.bonds.data:
#    print bond


#print itp.header.moleculetype.angles
#for angle in itp.header.moleculetype.angles.data:
#    print angle