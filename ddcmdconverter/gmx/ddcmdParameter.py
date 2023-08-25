import logging
from ddcmdconverter.gmx.gmxTopology import Toplogy
from ddcmdconverter.gmx.gmxParameter import Parameter
from ddcmdconverter.gmx.gmxTPRLog import TPRLog
from ddcmdconverter.gmx.gmxCoor import TprCoors, GroCoors
import math
import re
import os
import argparse

# set up logger
LOGLEVEL = 1
LOG_FMT = '%(asctime)s - %(name)s:%(funcName)s:%(lineno)s - %(levelname)s - %(message)s'
logger = logging.getLogger(__name__)
logger.setLevel(LOGLEVEL)

sh = logging.StreamHandler()
sh.setLevel(LOGLEVEL)
sh.setFormatter(logging.Formatter(LOG_FMT))
logger.addHandler(sh)

class ddcMDpara():
    def __init__(self, topology, parameter):
        self.top=topology
        self.par=parameter
        self.types= {}
        self.gatherTypes()
        self.fileHandle = None
        self.proteinList=["RAS", "RAF", "POY", "XZB", "WS4", "BVG", "JKI", "B2AR"]

    def gatherTypes(self):
        for moltype in self.par.moltypes:
            for atom in moltype.atomtypes.atoms:
                self.types[atom['typename'].strip('"')]=atom

    def add_proteinList(self, newProtList):
        if len(newProtList)>0:
            for prot in newProtList.split():
                self.proteinList.append(prot)
        logger.info("The current protein list is: ["+", ".join(self.proteinList)+"]")

    def _writeHeader(self):
        #print MMFF struct
        self.fileHandle.write("martini MMFF\n{\n")
        self.fileHandle.write("\tversion =1;\n")
        # print residue list
        self.fileHandle.write("\tresiParms=")

        for moltype in self.par.moltypes:
            # rename RAS and RAS/RAF
            if 'KRAS' in moltype.name:
                moltype.name = 'RAS'
            if 'RAF' in moltype.name:
                moltype.name = 'RAF'
            self.fileHandle.write(moltype.name+" ")
        self.fileHandle.write(";\n")
        # print atom type
        self.fileHandle.write("\tatomTypeList=")
        for type in self.types.keys():
            self.fileHandle.write(type + " ")
        self.fileHandle.write(";\n")
        # print lj list
        self.fileHandle.write("\tljParms=")
        typeList=list(self.types.keys())
        for i in range(len(typeList)):
            for j in range(i, len(typeList)):
                self.fileHandle.write(typeList[i] + "_" + typeList[j] + " ")
        self.fileHandle.write(";\n")
        self.fileHandle.write("}\n\n")

        for (type, typeifo)in self.types.items():
            self.fileHandle.write(type+"  MASSPARMS { atomType="+type+"; atomTypeID="+str(typeifo['type'])+"; mass="+str(typeifo['m'])+"M_p ; }\n")
        self.fileHandle.write("\n")

    def _writeResidues(self):
        # TODO move proteinList to arg list of program

        functypes=self.top.functype
        resID=0
        for moltype in self.par.moltypes:
            self.fileHandle.write(moltype.name+" RESIPARMS\n{\n")
            self.fileHandle.write("\tresID=" + str(resID) + ";\n")
            self.fileHandle.write("\tresType=" + str(0) + ";\n")
            self.fileHandle.write("\tresName=" + moltype.name + ";\n")
            isProtein=0
            if moltype.name in self.proteinList or moltype.name == 'CRS' or moltype.name =='RBS' or 'syn' in moltype.name:
                isProtein = 1
            self.fileHandle.write("\tisProtein=" + str(isProtein) + ";\n")
            numAtoms = len(moltype.atomtypes.atoms)
            self.fileHandle.write("\tnumAtoms=" + str(numAtoms) + ";\n")
            totCharge = 0
            for atom in moltype.atomtypes.atoms:
                totCharge += atom['q']
            self.fileHandle.write("\tcharge=" + str(totCharge) + ";\n")
            self.fileHandle.write("\tgroupList=" + moltype.name + "_g0" + ";\n")

            if len(moltype.bond.bonds) > 0:
                count=0
                self.fileHandle.write("\tbondList=")
                for bond in moltype.bond.bonds:
                    self.fileHandle.write(moltype.name + "_b" +str(count)+" ")
                    bond['id']=count
                    count=count+1
                self.fileHandle.write(";\n")
                #self.fileHandle.write(bondstr)

            if len(moltype.angle.angles) > 0:
                count=0
                self.fileHandle.write("\tangleList=")
                for angle in moltype.angle.angles:
                    self.fileHandle.write(moltype.name + "_a" + str(count) + " ")
                    angle['id']=count
                    count = count + 1
                self.fileHandle.write(";\n")
                #self.fileHandle.write(anglestr)

            if len(moltype.dihedral.dihedrals) >0:
                count=0
                self.fileHandle.write("\tdihedralList=")
                for dihedral in moltype.dihedral.dihedrals:
                    self.fileHandle.write(moltype.name+"_d" +str(count) +" ")
                    dihedral['id']=count
                    count = count + 1
                self.fileHandle.write(";\n")
                #self.fileHandle.write(dihedralstr)

            if moltype.name == "POPX":
                self.fileHandle.write("\trestraintList=POPX_r0 ;\n")

            #constraintStr="\tconstraints="
            if len(moltype.constraint.constraints)>0:
                """
                count=0
                constraintStr="\tconstraintList="
                for consList in moltype.constraint.consCluster:
                    constraintStr=constraintStr+moltype.name+"_cl"+str(count) +" "
                    count = count + 1
                constraintStr = constraintStr + ";\n"
                self.fileHandle.write(constraintStr)
                """
                count = 0

            if len(moltype.excls.exclsions)>0:
                count=0
                self.fileHandle.write("\texclusionList=")
                for exclusion in moltype.excls.exclsions:
                    self.fileHandle.write(moltype.name+"_e"+str(count) +" ")
                    count = count + 1
                self.fileHandle.write(";\n")
                #self.fileHandle.write(exclusionStr)

            if moltype.virtual and len(moltype.virtual.vSites)>0:
                """
                count=0
                vSiteStr="\tvsiteList="
                for vs in moltype.virtual.vSites:
                    vSiteStr=vSiteStr+moltype.name+"_vs"+str(count)+" "
                    vs['id']=count
                    count=count+1
                vSiteStr=vSiteStr+";\n"
                self.fileHandle.write(vSiteStr)
                """
                count = 0
            """
            constraintStr = constraintStr +";\n"
            if len(constraintStr) > 15:
                self.fileHandle.write(constraintStr)
            """

            self.fileHandle.write("\tcenterAtom=0;\n")
            self.fileHandle.write("}\n\n")
            resID += 1

            self.fileHandle.write(moltype.name+"_g0 GROUPPARMS{\n\tgroupID=0;\n\tatomList=")
            # if it is protein change atom name save original to oldname
            if isProtein == 1:
                count = 0
                for atom in moltype.atomtypes.atoms:
                    atom['oldname']=atom['name']
                    atom['name'] = 'P'+str(count)
                    count=count+1

            for atom in moltype.atomtypes.atoms:
                self.fileHandle.write(moltype.name+"_"+atom['name']+" ")
            self.fileHandle.write(";\n}\n\n")

            #atomID=0
            for atom in moltype.atomtypes.atoms:
                #atom['id']=atomID
                self.fileHandle.write(moltype.name + "_" + atom['name'] +
                        " ATOMPARMS{atomID="+str(atom['id'])+
                        "; atomName="+atom['name'] +
                        "; atomType="+atom['typename']+
                        "; atomTypeID="+str(atom['type'])+
                        "; charge="+str(atom['qB'])+
                        "; mass="+str(atom['m'])+
                        " M_p ; }\n")
                #atomID=atomID+1
            self.fileHandle.write("\n")

            #bond
            atoms=moltype.atomtypes.atoms
            for bond in moltype.bond.bonds:
                indexI=bond['atom1']
                indexJ=bond['atom2']
                atomI=atoms[indexI]
                atomJ=atoms[indexJ]
                func=functypes[bond['type']]
                kb=(func['cbA']+func['cbB'])/4  # average and ddcMD uses 1/2 k
                b0=(func['b0A']+func['b0B'])/2  # average

                self.fileHandle.write(moltype.name + "_b"+str(bond['id'])+
                        " BONDPARMS{atomI="+str(indexI)+
                        "; atomTypeI="+atomI['typename']+
                        "; atomJ="+str(indexJ)+
                        "; atomTypeJ="+atomJ['typename']+
                        "; func=1; kb="+str(kb)+
                        " kJ*mol^-1*nm^-2; b0="+str(b0)+
                        " nm; }\n")

            self.fileHandle.write("\n")

            #angle
            for angle in moltype.angle.angles:
                indexI = angle['atom1']
                indexJ = angle['atom2']
                indexK = angle['atom3']
                # atomI = atoms[indexI]
                # atomJ = atoms[indexJ]
                # atomK = atoms[indexK]
                func = functypes[angle['type']]
                try:
                    if func['type'] =='ANGLES' or func['type'] =='G96ANGLES':
                        ktheta=func['ctA']/2 # ddcMD uses 1/2 k
                        theta0=func['thA'] #
                    elif func['type'] =='RESTRANGLES':
                        #FIXED: `gmx dump -s topol.tpr > tpr.log` output wrong output wrong RESTRANGLES parameters.
                        #FIXED: After gromacs fix bug, the code show be changed back.
                        #ktheta=func['kthetaA']/2 # ddcMD uses 1/2 k
                        #theta0=func['costheta0A']
                        ktheta=func['costheta0A']/2 # ddcMD uses 1/2 k
                        theta0=func['kthetaA']
                except:
                    print(moltype.name, indexI, indexJ, indexK, angle['type'], func)
                    continue
                if func['type'] == 'ANGLES':
                    funcID = 1
                    #theta0 = theta0 * math.pi/180  # deg to rad
                elif func['type'] == 'G96ANGLES':
                    funcID = 2
                    theta0 = round(math.acos(theta0) * 180 / math.pi, 2) # convert to deg and round up to 2 decimal places
                elif func['type'] == 'RESTRANGLES':
                    funcID = 10
                    #theta0 = theta0 * math.pi / 180  # deg to rad
                else:
                    funcID = 0
                self.fileHandle.write(moltype.name + "_a"+str(angle['id'])+
                        " ANGLEPARMS{atomI="+str(indexI)+
                        "; atomJ="+str(indexJ)+
                        "; atomK=" + str(indexK) +
                        "; func="+ str(funcID) +
                        "; ktheta="+str(ktheta))
                if(funcID==2 or funcID==10):
                    self.fileHandle.write(" kJ*mol^-1;")
                else:
                    self.fileHandle.write(" kJ*mol^-1*rad^-2;")
                self.fileHandle.write(" theta0="+str(theta0)+
                        " deg; }\n")

            self.fileHandle.write("\n")

            #dihedral
            for dihedral in moltype.dihedral.dihedrals:
                indexI = dihedral['atom1']
                indexJ = dihedral['atom2']
                indexK = dihedral['atom3']
                indexL = dihedral['atom4']
                func=functypes[dihedral['type']]

                unit="kJ*mol^-1"
                if func['type'] == 'IDIHS': #IMPROPER
                    funcID = 2
                    kchi=(func['cxA']+func['cxB'])/4 # average and ddcMD IMPROPER uses 1/2 k
                    delta=(func['xiA']+func['xiB'])/2 # average
                    unit = "kJ*mol^-1*rad^-2"
                    #delta = math.pi * delta / 180
                elif func['type'] =='PDIHS': #TORSION
                    funcID = 1
                    kchi=(func['cpA']+func['cpB'])/2 # average and ddcMD TORSION uses k
                    delta=(func['phiA']+func['phiB'])/2 # average
                    #delta = math.pi * delta / 180  # deg to rad and ddcMD use opposite sign for dihedral
                else:
                    logger.error("Dihedral type is neither IDIHS nor PDIHS, skipping")
                    raise

                #RAS_d0 TORSPARMS{atomI=26; atomJ=29; atomK=31; atomL=32; func=1; delta=2.0943951023931953; kchi=400.0 kJ*mol^-1; n=1;}
                self.fileHandle.write(moltype.name + "_d"+str(dihedral['id'])+
                        " TORSPARMS{atomI="+str(indexI)+
                        "; atomJ="+str(indexJ)+
                        "; atomK=" + str(indexK) +
                        "; atomL=" + str(indexL) +
                        "; func=" + str(funcID) +
                        "; kchi="+str(kchi) +
                        " "+unit+"; delta="+str(delta)+
                        " deg;  n=1; }\n")

            self.fileHandle.write("\n")

            #restraint for POPX
            #if moltype.name=="POPX":
            #    for atom in moltype.atomtypes.atoms:
            #        if atom['name']=="PO4":
            #            self.fileHandle.write("POPX_r0 RESTRAINTPARMS{atomI="+str(atom['id'])+
            #                    "; atomTypeI="+atom['typename']+
            #                    "; func=1; fcx=0; fcy=0; fcz=2; kb= 1.0 kJ*mol^-1*nm^-2; }\n\n")
            #    self.fileHandle.write("\n")

            #constraint
            #CHOL_cl0 CONSLISTPARMS{ constraintSubList=CHOL_cl0_c0 CHOL_cl0_c1 CHOL_cl0_c2 CHOL_cl0_c3 CHOL_cl0_c4 CHOL_cl0_c5 ; }
            #CHOL_cl0_c0 CONSPARMS{atomI=0; atomTypeI=SP1; atomJ=2; atomTypeJ=SC3; func=1; r0=0.493 nm; }
            """
            count = 0
            for consCl in moltype.constraint.consCluster:
                clName = moltype.name + "_cl" + str(count)
                clStr = clName + " CONSLISTPARMS{ constraintSubList="
                count2 = 0
                cStr=""
                for cons in consCl:
                    cName=clName + "_c" + str(count2)
                    clStr=clStr+cName+" "
                    indexI = cons['atom1']
                    indexJ = cons['atom2']
                    atomI = atoms[indexI]
                    atomJ = atoms[indexJ]
                    func = functypes[cons['type']]
                    r0=(func['dA']+func['dB'])/2
                    cStr=cStr+cName+"  CONSPARMS{atomI="+str(indexI)+\
                                    "; atomTypeI="+atomI['typename']+\
                                    "; atomJ="+str(indexJ)+\
                                    "; atomTypeJ="+atomJ['typename']+\
                                    "; func=1; r0="+str(r0)+" nm; }\n"
                    count2 = count2 + 1
                clStr=clStr+";}\n"
                cStr = cStr + "\n"
                self.fileHandle.write(clStr)
                self.fileHandle.write(cStr)
                count = count + 1
            """

            # Exclusion list
            count = 0
            for exclusion in moltype.excls.exclsions:
                exName = moltype.name + "_e"+ str(count)
                indexI = exclusion['atom1']
                indexJ = exclusion['atom2']
                atomI = atoms[indexI]
                atomJ = atoms[indexJ]
                exStr = exName + " EXCLUDEPARMS{atomI="+str(indexI)+\
                                    "; atomTypeI="+atomI['typename']+\
                                    "; atomJ="+str(indexJ)+\
                                    "; atomTypeJ="+atomJ['typename']+\
                                    ";}\n"
                self.fileHandle.write(exStr)
                count = count + 1
            self.fileHandle.write("\n")
            
            # Virtual site
            if moltype.virtual:
                """
                for vs in moltype.virtual.vSites:
                    funct=functypes[vs['type']]
                    vsName=moltype.name+"_vs"+str(vs['id'])
                    vsType=funct['type']
                    vsStr=vsName+" VSITEPARMS{"
                    for k, v in funct.items():
                        if vsType == 'VSITE3OUT' and k=='c':
                            v = v*0.1
                        vsStr=vsStr + k + "=" +str(v) +"; "
                    count=0
                    for i in vs['idx']:
                        count=count+1
                        vsStr = vsStr + "atom"+str(count)+"="+str(i)+"; "
                        atomI = atoms[int(i)]
                        vsStr = vsStr + "atomType" + str(count) + "=" + atomI['typename'] + "; "
                    vsStr = vsStr + "}\n"
                    self.fileHandle.write(vsStr)
                self.fileHandle.write("\n")
                """
                count = 0

    def _writeLJPARMS(self):
        numTypes = self.top.atnr
        functypes = self.top.functype
        ljPairs=[]
        for funct in functypes:
            if funct['type'] == 'LJ_SR':
                ljPairs.append(funct)
        if len(ljPairs)!=numTypes*numTypes:
            logger.warning("Number of LJ pairs in functypes doesn't match number of Type")
            raise

        lj2dParms=[]
        count=0
        for i in range(numTypes):
            ljRowParm=[]
            for j in range(numTypes):
                lj={}
                ljPair=ljPairs[count]
                c6 = ljPair['c6']
                c12 = ljPair['c12']
                if c6 != 0 or c12 != 0 :
                    sigma = (c12 / c6) ** (1 / 6.0)
                    eps = c6 * c6 / (4 * c12)
                else:
                    sigma = 0
                    eps = 0
                lj['sigma']=round(sigma,3)
                lj['eps']  =round(eps, 3)
                ljRowParm.append(lj)
                count=count+1
            lj2dParms.append(ljRowParm)

        atomTypeIDmap={}
        #print(self.types)
        for key, value in self.types.items():
            atomTypeIDmap[value['type']]=key

        #Q0_Q0 LJPARMS{atomtypeI=Q0; indexI=19; atomtypeJ=Q0; indexJ=19; sigma=0.47 nm; eps=3.5 kJ*mol^-1;}
        count=0
        for i in range(numTypes):
            for j in range(i, numTypes):
                self.fileHandle.write(atomTypeIDmap[i] + "_" + atomTypeIDmap[j] +
                                      " LJPARMS{atomtypeI=" + atomTypeIDmap[i] +
                                      "; indexI=" + str(i) +
                                      "; atomtypeJ=" + atomTypeIDmap[j] +
                                      "; indexJ=" + str(j) +
                                      "; sigma=" +str(lj2dParms[i][j]['sigma'])+
                                      " nm; eps="+str(lj2dParms[i][j]['eps'])+
                                      " kJ*mol^-1;}\n")


    def toFile(self, filename):
        with open(filename, 'w') as self.fileHandle:
            logger.info('write martini.data header')
            self._writeHeader()
            logger.info('write martini.data residues')
            self._writeResidues()
            logger.info('write martini.data LJparms')
            self._writeLJPARMS()

class ddcMDobj():
    def __init__(self, topology, parameter, coordinates):
        self.top=topology
        self.par=parameter
        self.coordinates=coordinates

    def getHeader(self):
        box = self.coordinates.box
        totNumAtoms=self.coordinates.numAtoms
        header = "particle FILEHEADER {type=MULTILINE; datatype=VARRECORDASCII; checksum=NONE;\n"
        header = header + "loop=0; time=0.000000;\n"
        header = header + "nfiles=1; nrecord=" + str(totNumAtoms) + "; nfields=10;\n"
        header = header + "field_names=id class type group rx ry rz vx vy vz;\n"
        header = header + "field_types=u s s s f f f f f f;\n"
        header = header + "h=     " + str(box[0][0]) + "   " + str(box[0][1]) + "   " + str(box[0][2]) + "\n"
        header = header + "       " + str(box[1][0]) + "   " + str(box[1][1]) + "   " + str(box[1][2]) + "\n"
        header = header + "       " + str(box[2][0]) + "   " + str(box[2][1]) + "   " + str(box[2][2]) + ";\n"
        header = header + "groups = group ;\n"
        header = header + "types = ATOM ;\n"
        header = header + "} \n\n"

        return header

    def getRestart(self):
        box = self.coordinates.box
        totNumAtoms = self.coordinates.numAtoms
        line = "simulate SIMULATE { loop=0; time=0.000000 ;}\n"
        line = line + "box BOX {\n"
        line = line + "h=     " + str(box[0][0]) + "   " + str(box[0][1]) + "   " + str(box[0][2]) + "\n"
        line = line + "       " + str(box[1][0]) + "   " + str(box[1][1]) + "   " + str(box[1][2]) + "\n"
        line = line + "       " + str(box[2][0]) + "   " + str(box[2][1]) + "   " + str(box[2][2]) + ";\n"
        line = line + "}\n"
        line = line + "collection COLLECTION { mode=VARRECORDASCII; size=" + str(
            totNumAtoms) + "; files=snapshot.mem/atoms#;}\n"

        outFh = open("restart", "w")
        outFh.write(line)

    def toFile(self, filename):
        self.getRestart()
        count=-1
        moltypes=self.par.moltypes
        molID=-1
        with open(filename, 'w') as f:
            header = self.getHeader()
            f.write(header)
            for molBlk in self.top.molBlks:
                for i in range(molBlk.numMols):
                    molID=molID+1
                    moltype = moltypes[molBlk.moltypeID]
                    assert molBlk.moltypeName == moltype.name, molBlk.moltypeName+" and "+moltype.name+" molBlk name doesn't match moltype name"
                    for atom in moltype.atomtypes.atoms:
                        count=count+1
                        speice = moltype.name + 'x' + atom['name']
                        (idx, x, y, z) = self.coordinates.coors[count]
                        assert idx == count, "index in coordinates.coorlines doesn't match count"
                        (idx, vx, vy, vz) = self.coordinates.velocity[count]
                        assert idx == count, "index in coordinates.velocitylines doesn't match count"

                        gid = (molID << 32) + atom['id']
                        outLine = "%16d ATOM %15s %6s %21.13e %21.13e %21.13e %21.13e %21.13e %21.13e\n" \
                                  % (gid, speice, atom['group'], x, y, z, vx, vy, vz)
                        f.write(outLine)

    def toRestraint(self, filename):
        functypes = self.top.functype
        box = self.top.box # use box from tpr.log not from gro
        xbox = box[0][0]
        ybox = box[1][1]
        zbox = box[2][2]
        func = 1
        restraintCount=0
        moltypes=self.par.moltypes
        molID=-1
        with open(filename, 'w') as f:
            header = "restraint RESTRAINTLIST{\n  xbox = " + str(xbox) \
                     + ";\n  ybox = " + str(ybox) \
                     + ";\n  zbox = " + str(zbox) + ";\n"
            outLine = ""
            for molBlk in self.top.molBlks:
                count = 0
                for i in range(molBlk.numMols):
                    molID=molID+1
                    moltype = moltypes[molBlk.moltypeID]
                    assert molBlk.moltypeName == moltype.name, molBlk.moltypeName+" and "+moltype.name+" molBlk name doesn't match moltype name"
                    if len(moltype.restraint.restraints)>0:
                        for restraint in moltype.restraint.restraints:
                            funct=functypes[restraint['type']]
                            atomI = restraint['atom1']
                            gid = (molID << 32) + atomI
                            restraintLine = molBlk.posResxA[count+atomI]
                            restrStrs = re.split('\{|\}', restraintLine)
                            xyzStrs = restrStrs[1].split(',')
                            xyz =[10*float(x) for x in xyzStrs] # convert nm to angstrom
                            outLine = outLine + "r_" + str(restraintCount) \
                                      + " RESTRAINTPARMS{gid=" + str(gid) \
                                      + "; atomI=" + str(atomI) \
                                      + "; func=" + str(func) \
                                      + "; fcx=" + str(funct['fcx']) \
                                      + "; fcy=" + str(funct['fcy']) \
                                      + "; fcz=" + str(funct['fcz']) \
                                      + "; x0=" + str(xyz[0] / xbox) \
                                      + "; y0=" + str(xyz[1] / ybox) \
                                      + "; z0=" + str(xyz[2] / zbox) \
                                      + "; kb= 1.0 kJ*mol^-1*nm^-2; }\n"
                            restraintCount = restraintCount + 1
                            #print(outLine)

                        count = count + len(moltype.atomtypes.atoms)

            header = header+"restraintList="
            for i in range(restraintCount):
                header = header + "r_" + str(i) +" "
            header = header + ";\n}\n\n"
            f.write(header)
            f.write(outLine)

class ddcMDinput():
    def __init__(self, topology, parameter, ddcmdpara):
        self.top=topology
        self.par=parameter
        self.ddc=ddcmdpara

    def tofile(self, filename):
        molecLineFH=open('molecule_1.data', 'w')
        moleSpecieFH=open('molecule_2.data', 'w')
        speLineFH=open('molecule_3.data', 'w')
        functypes = self.top.functype
        types=self.ddc.types
        molecLineFH.write("moleculeClass MOLECULECLASS { molecules = ")
        moleSpecieFH.write("\n")
        speLineFH.write("\n")
        for moltype in self.par.moltypes:
            resName = moltype.name
            molName = resName + "x"
            molecLineFH.write(" " + molName)

            moleSpecieFH.write( molName + " MOLECULE {ownershipSpecies = ")

            count = 0
            for atom in moltype.atomtypes.atoms:
                specie = molName + atom['name']

                atomType = atom['typename']
                speLineFH.write( specie + " SPECIES { type = ATOM ; charge =" + str(atom['q']) + "; id=" + str(
                    types[atomType]['type']) + "; mass =" + str(atom['m']) + " M_p ; }\n")

                if count == 0:
                    moleSpecieFH.write(specie + "; species = " + specie)
                else:
                    moleSpecieFH.write(" " + specie)
                count = count + 1

            moleSpecieFH.write(";")

            hasCons=len(moltype.constraint.constraints) > 0
            hasVsite=moltype.virtual and len(moltype.virtual.vSites) > 0

            if hasCons or hasVsite:
                moleSpecieFH.write("\n\tconstraints=")
            if hasCons:
                count = 0
                for cons in moltype.constraint.constraints:
                    cName = moltype.name + "_C" + str(count)
                    moleSpecieFH.write(cName + " ")
                    count = count + 1

            if hasVsite:
                count = 0
                for vs in moltype.virtual.vSites:
                    vsName = moltype.name + "_VS" + str(count)
                    moleSpecieFH.write(vsName + " ")
                    count = count + 1

            if hasCons or hasVsite:
                moleSpecieFH.write(";")
            #if len(constraintStr) > 15:
            #    moleSpecie = moleSpecie + constraintStr

            moleSpecieFH.write("}\n")

            # constraint
            count = 0
            #cStr=""
            for cons in moltype.constraint.constraints:
                cName = moltype.name + "_C" + str(count)
                moleSpecieFH.write(cName + " ")
                indexI = cons['atom1']
                indexJ = cons['atom2']
                func = functypes[cons['type']]
                r0 = (func['dA'] + func['dB']) / 2
                moleSpecieFH.write("  CONSTRAINT { type=bondLength; offsets=" + \
                        str(indexI) +" " + str(indexJ) + ";" + \
                       " r0=" + str(r0) + " nm; }\n")
                count = count + 1
            #moleSpecie = moleSpecie+cStr

            #virtual site
            if moltype.virtual:
                count = 0
                for vs in moltype.virtual.vSites:
                    funct=functypes[vs['type']]
                    vsName=moltype.name+"_VS"+str(count)
                    vsType=funct['type']
                    if vsType =='VSITEN':
                        num = len(vs['idx']) -1
                        funct['type']='VSITE'+str(num)+'N'

                    vsStr=vsName+" CONSTRAINT{"
                    for k, v in funct.items():
                        if vsType == 'VSITE3OUT' and k=='c':
                            v = v*0.1
                        vsStr=vsStr + k + "=" +str(v) +"; "
                    vsStr=vsStr + " offsets="
                    for i in vs['idx']:
                        vsStr = vsStr +str(i) + " "
                    vsStr = vsStr + ";}\n"
                    count = count + 1
                    moleSpecieFH.write(vsStr)
                moleSpecieFH.write("\n")

        molecLineFH.write("; }\n")

        molecLineFH.close()
        moleSpecieFH.close()
        speLineFH.close()
        os.system('cat molecule_1.data molecule_2.data molecule_3.data >' + filename )
        #with open(filename, "w") as f:
        #    f.write(molecLine)
        #    f.write(moleSpecie)
        #    f.write(speLine)

def getArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--intpr', action='store', dest='inTprFile', default="tpr.log",
                        help='Input file uses Gromacs TPR log file (default=tpr.log).')
    parser.add_argument('-g', '--ingro', action='store', dest='inGroFile', default="prod.gro",
                        help='Input file uses Gromacs Gro file (default=prod.gro).')
    parser.add_argument('-m', '--outm', action='store', dest='outMartiniFile', default="martini.data",
                        help='Output martini.data file (default=martini.data).')
    parser.add_argument('-l', '--outl', action='store', dest='outMoleculeFile', default="molecule.data",
                        help='Output molecule.data file (default=molecule.data).')
    parser.add_argument('-r', '--outr', action='store', dest='outRestraintFile', default="restraint.data",
                        help='Output restraint.data file (default=restraint.data).')
    parser.add_argument('-p', '--prot', action='store', dest='newProtList', default="",
                        help='Output restraint.data file (default=RAS RAF POY XZB WS4 BVG JKI B2AR).')

    args = parser.parse_args()

    return args

def main():
    logger.info("Please run 'gmx dump -s topol.tpr > tpr.log' command if tpr.log input file is not available")
    args = getArgs()
    # 1. Parse tpr.log file
    logger.info("Parse tpr.log")
    tprLog=TPRLog(args.inTprFile)
    logger.info("Parse topology")
    topology=Toplogy(tprLog)
    logger.info("Parse parameter")
    parameter=Parameter(tprLog)
    # 2. Generate martini.data file
    logger.info("Generate martini.data")
    ddcmdPar=ddcMDpara(topology, parameter)
    ddcmdPar.add_proteinList(args.newProtList)
    ddcmdPar.toFile(args.outMartiniFile)
    # 3. Parse Gromacs coordinate file
    logger.info("Parse coordinate")
    #coordinates=TprCoors(tprLog, topology)
    #coordinates = GroCoors("start.gro")
    coordinates = GroCoors(args.inGroFile)
    # 4. Generate ddcMD coordinate file
    logger.info("Generate atoms#000000 and restart")
    ddcmdobj=ddcMDobj(topology, parameter, coordinates)
    ddcmdobj.toFile('atoms#000000')
    # 5. Generate ddcMD restraint file
    logger.info("Generate restraint.data")
    ddcmdobj.toRestraint(args.outRestraintFile)
    # 6. Generate ddcMD molecule file
    logger.info("Generate molecule.data")
    ddcmdinput=ddcMDinput(topology, parameter, ddcmdPar)
    ddcmdinput.tofile(args.outMoleculeFile)
    logger.info("End")

if __name__ == '__main__':
    main()
