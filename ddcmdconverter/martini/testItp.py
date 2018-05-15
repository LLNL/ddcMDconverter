from ddcmdconverter.martini.ITP import ITP

#itp=ITP.ITP("martini_v2.0_POPC_02.itp")
#itp=ITP.ITP("martini_v2.1-dna.itp")
itp=ITP("KRAS-GTP-04-HVR-best-guess-M22-CYFpos.itp")
print(itp.header.moleculetype.data['name'])
#itp=ITP.ITP("martini_v2.0_CHOL_01.itp")
#print itp.header
print(itp.header.moleculetype)

#print itp.header.moleculetype.atoms

count=1
for atom in itp.header.moleculetype.atoms.data:
    atom['atomname']='P'+str(count)
    count=count+1

print(itp.header.moleculetype.atoms)

#print itp.header.moleculetype.data
print(itp.header.moleculetype.data['name'])

#print itp.header.moleculetype.sections.keys()

#print itp.header.moleculetype.exclusions.data

#print itp.header.atomtypes.data

#print itp.header.nonbond_params.data

#print sum(itp.header.moleculetype.atoms.data.charge)

#print itp.header.moleculetype.bonds.data

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