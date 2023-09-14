import logging
import re

LOGLEVEL =2
LOG_FMT = '%(asctime)s - %(name)s:%(funcName)s:%(lineno)s - %(levelname)s - %(message)s'
logger = logging.getLogger(__name__)
logger.setLevel(LOGLEVEL)

sh = logging.StreamHandler()
sh.setLevel(LOGLEVEL)
sh.setFormatter(logging.Formatter(LOG_FMT))
logger.addHandler(sh)


def getSecName(line, level=2):
    idx=3*level
    if len(line)>0 and line[0] == ',':
        return None
    if line[idx:idx+1] == " " or line[0:idx] != " "*idx:
        return None
    else:
        strs=re.split(':|=', line[idx:])
        return strs[0]

class Parameter():
    def __init__(self, tprLog):
        self.sections = ['atomtypes', 'moltype', 'end']
        self.atomtypes=[]
        self.moltypes = []
        self.lines = tprLog.getCmap()

        self._parse()

    def _parse(self):
        secIdx=1
        secLen=len(self.sections[secIdx])
        count=0
        for line in self.lines:
            count+=1
            if line[3:3+secLen]==self.sections[secIdx]:
                secIdx+=1
                secLen = len(self.sections[secIdx])
            if secIdx==1:
                self.atomtypes.append(line)
            if secIdx>1:
                break

        hasMoltype=False
        for line in self.lines[count-1:]:
            #print(line)
            #print("count", count, 'section', self.sections[secIdx])
            if line[3:10]=='moltype':
                #print(line)
                hasMoltype=True
                moltype=MolType()
                self.moltypes.append(moltype)
            if hasMoltype:
                moltype.lines.append(line)

        logger.info("Parse moltype ...")
        for moltype in self.moltypes:
            moltype.parse()
            print(moltype.name, end =" ")
        print(" ")

class MolType():
    def __init__(self):
        '''
        self.sections={'name', 'atoms', 'cgs', 'excls', 'Bond','G96Bond',
                        'Harmonic', 'FENE','Restraint Pot.','Angle','G96Angle',
                        'Restricted', 'Proper', 'Ryckaert-Bell', 'Improper',
                        'Tab. Dih.', 'Constraint', 'Constr. No Conn.', 'Position Rest.',
                        'Virtual site 2', 'Virtual site 3', 'Virtual site 3fd',
                        'Virtual site 3fad', 'Virtual site 3out', 'Virtual site 4fd',
                        'Virtual site 4fdn', 'Virtual site N', 'COM Pull En', 'end'}
        '''
        self.lines=[]
        self.name=''
        self.atomtypes=None
        #self.cgs=[]
        self.excls=None
        self.bond=None
        self.angle=None
        self.dihedral=None
        self.constraint=None
        self.restraint=None
        self.virtual = None

    def parse(self):
        section=""
        for line in self.lines:
            secName=getSecName(line, level=2)
            if secName:
                section=secName

            if 'name' in section:
                self.name=line.split('=')[1].strip('"').strip('+').strip('-')

            if 'atoms'in section:
                if not self.atomtypes:
                    self.atomtypes=AtomTypes()
                self.atomtypes.lines.append(line)

            #if section == 'cgs':
            #    self.cgs.append(line)

            if 'excls' in section:
                if not self.excls:
                    self.excls=Exclusions()
                self.excls.lines.append(line)

            if 'Bond' in section:
                if not self.bond:
                    self.bond=Bonds()
                self.bond.lines.append(line)

            if 'Harmonic' in section:
                if not self.bond:
                    self.bond=Bonds()
                self.bond.lines.append(line)

            if 'Angle' in section:
                if not self.angle:
                    self.angle=Angles()
                self.angle.lines.append(line)
            # add G96Angle to angle list
            #if section == 'G96Angle':
            #    if not self.angle:
            #        self.angle=Angles()
            #    self.angle.lines.append(line)

            if 'Dih.' in section:
                if not self.dihedral:
                    self.dihedral=Dihedral()
                self.dihedral.lines.append(line)

            #if section == 'Improper':
            #    if not self.dihedral:
            #       self.dihedral = Dihedral()
            #    self.dihedral.lines.append(line)

            if 'Constraint' in section:
                if not self.constraint:
                    self.constraint = Constraint()
                self.constraint.lines.append(line)

            if 'Position Rest.' in section:
                if not self.restraint:
                    self.restraint = Restraint()
                self.restraint.lines.append(line)

            #'Virtual site 2', 'Virtual site 3', 'Virtual site 3fd',
            #'Virtual site 3fad', 'Virtual site 3out', 'Virtual site 4fd',
            #'Virtual site 4fdn', 'Virtual site N',
            if 'Virtual' in section:
                if not self.virtual:
                    self.virtual=VirtualSite()
                self.virtual.lines.append(line)

        logger.info("Parse atomtypes ...")
        self.atomtypes.parse()
        self.excls.parse()
        logger.info("Parse bond ...")
        self.bond.parse()
        logger.info("Parse angle ...")
        self.angle.parse()
        logger.info("Parse dihedral ...")
        self.dihedral.parse()
        logger.info("Parse constraint ...")
        self.constraint.parse()
        logger.info("Parse restraint ...")
        self.restraint.parse()
        logger.info("Parse virtual ...")
        if self.virtual:
            self.virtual.parse()
        logger.info("Set group ...")
        self.setGroup()

    def setGroup(self):
        logger.info("SetGroup set group...")
        for atom in self.atomtypes.atoms:
            atom['group']='group'

        logger.info("SetGroup constraint ...")
        for cons in self.constraint.constraints:
            atom1=self.atomtypes.atoms[cons['atom1']]
            atom1['group'] = 'free'
            atom2=self.atomtypes.atoms[cons['atom2']]
            atom2['group'] = 'free'

        logger.info("SetGroup Virtual ...")
        if self.virtual:
            for vs in self.virtual.vSites:
                vIdx=int(vs['idx'][0])
                atom = self.atomtypes.atoms[vIdx]
                atom['group'] = 'vsite'

        '''
        consAtoms=[]
        logger.info("SetGroup constraint ...")
        for cons in self.constraint.constraints:
            consAtoms.append(cons['atom1'])
            consAtoms.append(cons['atom2'])

        logger.info("Setgroup constraint unique ...")
        uniqConsAtoms=[]
        for x in consAtoms:
            # check if exists in unique_list or not
            if x not in uniqConsAtoms:
                uniqConsAtoms.append(x)
        logger.info("Setgroup Virtual ...")
        vSiteAtoms=[]
        if self.virtual:
            for vs in self.virtual.vSites:
                vIdxStr=vs['idx'][0]
                vSiteAtoms.append(int(vIdxStr))
        logger.info("Setgroup Virtual unique...")
        uniqVsAtoms=[]
        for x in vSiteAtoms:
            if x not in uniqVsAtoms:
                uniqVsAtoms.append(x)
        logger.info("Setgroup set group...")
        for atom in self.atomtypes.atoms:
            if atom['id'] in uniqConsAtoms:
                atom['group']='free'
            elif atom['id'] in uniqVsAtoms:
                atom['group'] = 'vsite'
            else:
                atom['group']='group'
        '''

class AtomTypes():
    def __init__(self):
        self.sections = ['atom', '']
        self.lines = []
        self.atom1 = []
        self.atom2 = []
        self.atom3 = []
        self.atoms = []

    def isfloat(self, value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    def toDictionary(self, str):
        strs=re.split('\{|\}', str)
        if len(strs) !=3:
            logger.error("Format of "+str+" doesn't have match curly brackets")

        atomstrs=re.split('\[|\]', strs[0])
        if len(atomstrs) != 3:
            logger.error("Format of " + strs[0] + " doesn't have match square brackets")
        idx = int(atomstrs[1])

        vDict={}
        valuestrs=strs[1].split(',')
        for vstr in valuestrs:
            vstrsplits=vstr.split('=')
            if len(vstrsplits) != 2:
                logger.error("Format of " + vstr + " is missing = sign")
            value=vstrsplits[1].strip().strip('"')
            if value.isdigit():
                value = int(value)
            elif self.isfloat(value):
                value = float(value)

            vDict[vstrsplits[0].strip()] =value

        return (idx, vDict)


    def parse(self):
        count=0
        logger.info("Parse atoms1 ... number line ="+str(len(self.lines)))
        for line in self.lines:
            count+=1
            if line[9:13] =='atom':
                numAtoms1=int(re.split('\(|\)', line)[1])
                self.atom1=self.lines[count:count+numAtoms1]
                break

        logger.info("numAtoms1 =" + str(numAtoms1))
        logger.info("Parse atoms2 ...")
        count=count+numAtoms1
        for line in self.lines[count:]:
            count += 1
            if line[9:13] =='atom':
                numAtoms2=int(re.split('\(|\)', line)[1])
                self.atom2=self.lines[count:count+numAtoms2]
                break
        logger.info("numAtoms2 =" + str(numAtoms2))
        logger.info("Parse atoms3 ...")
        count = count + numAtoms2
        for line in self.lines[count:]:
            count += 1
            if line[9:13] =='type':
                numAtoms3=int(re.split('\(|\)', line)[1])
                self.atom3=self.lines[count:count+numAtoms3]
                break
        logger.info("numAtoms3 =" + str(numAtoms3))
        if numAtoms1 != numAtoms2:
            logger.error("Size of atom list 1 doesn't equal to 2")

        if numAtoms1 != numAtoms3:
            logger.error("Size of atom list 1 doesn't equal to type list")

        logger.info("Mege dictionary ...")
        for i in range(numAtoms1):
            (idx1, vDict1) = self.toDictionary(self.atom1[i])
            (idx2, vDict2) = self.toDictionary(self.atom2[i])
            (idx3, vDict3) = self.toDictionary(self.atom3[i])
            if idx1 != idx2:
                logger.error("Indexs in"+self.atom1[i] + "and "+ self.atom2[i]+" not equal")
            if idx1 != idx3:
                logger.error("Indexs in"+self.atom1[i] + "and "+ self.atom3[i]+" not equal")

            mergeDict={**vDict1, **vDict2}
            mergeDict['typename'] = vDict3['name']
            self.atoms.append(mergeDict)

        logger.info("Set atom id ...")
        #set atom id
        atomID = 0
        for atom in self.atoms:
            atom['id'] = atomID
            atomID = atomID + 1

class Exclusions():
    def __init__(self):
        self.lines = []
        self.exclsions=[]
        self.numExcls=0

    def parse(self):
        pairList=[]
        line_iter = iter(self.lines)
        for line in line_iter:
            if 'numLists' in line:
                self.numExcls = int(line.split('=')[1])
            elif 'excls' in line:
                if 'excls:' in line:
                    continue
                strs=re.split('\{|\}', line)
                if len(strs)!=3:
                    if len(strs) == 2:
                        nextline = next(line_iter)
                        line = line + nextline
                        while ("}" not in nextline):
                            #print("2.. ", nextline)
                            nextline = next(line_iter)
                            line = line + nextline
                        strs = re.split('\{|\}', line)
                        #print("Multiline excls:   "+line)
                    else:
                        logger.warning("Incorrect excls format: " + line)
                        continue

                headers=re.split('\[|\]', strs[0])
                headIndex = int(headers[1])
                substrs = strs[1].split(",")
                subList = [int(i) for i in substrs]
                for i in subList:
                    if i != headIndex:
                        pairList.append((headIndex, i))

        uniqList0=set(pairList)
        uniqList1=set((a, b) if a <= b else (b, a) for a, b in uniqList0)

        for pair in uniqList1:
            exclusion = {}
            exclusion['atom1'] = pair[0]
            exclusion['atom2'] = pair[1]
            self.exclsions.append(exclusion)


class Bonds():
    def __init__(self):
        self.lines = []
        self.bonds =[]

    def parse(self):
        for line in self.lines:
            if 'BONDS' in line or 'HARMONIC' in line:
                strs=re.split('\(|\)', line)
                if len(strs)!=3:
                    logger.warning("Incorrect bond format: " + line)
                    continue
                bond={}
                bond['type'] = int(strs[0].split('=')[1])
                pair=strs[2].split()
                bond['atom1'] =int(pair[0])
                bond['atom2'] = int(pair[1])
                self.bonds.append(bond)

class Angles():
    def __init__(self):
        self.lines = []
        self.angles =[]

    def parse(self):
        for line in self.lines:
            if 'ANGLES' in line:
                strs=re.split('\(|\)', line)
                if len(strs)!=3:
                    logger.warning("Incorrect bond format: " + line)
                    continue
                angle={}
                angle['type'] = int(strs[0].split('=')[1])
                triple=strs[2].split()
                angle['atom1'] = int(triple[0])
                angle['atom2'] = int(triple[1])
                angle['atom3'] = int(triple[2])
                self.angles.append(angle)

class Dihedral():
    def __init__(self):
        self.lines=[]
        self.dihedrals=[]
    def parse(self):
        for line in self.lines:
            if 'PDIHS' in line or 'IDIHS' in line:
                strs=re.split('\(|\)', line)
                if len(strs)!=3:
                    logger.warning("Incorrect bond format: " + line)
                    continue
                dihedral={}
                dihedral['type'] = int(strs[0].split('=')[1])
                triple=strs[2].split()
                dihedral['atom1'] = int(triple[0])
                dihedral['atom2'] = int(triple[1])
                dihedral['atom3'] = int(triple[2])
                dihedral['atom4'] = int(triple[3])
                self.dihedrals.append(dihedral)

class Constraint():
    def __init__(self):
        self.lines=[]
        self.constraints=[]
        self.consCluster=[]

    def parse(self):
        logger.info("Constraint begin")
        for line in self.lines:
            if 'CONSTR' in line:
                strs=re.split('\(|\)', line)
                if len(strs)!=3:
                    logger.warning("Incorrect bond format: " + line)
                    continue
                cons={}
                cons['type'] = int(strs[0].split('=')[1])
                pair=strs[2].split()
                cons['atom1'] =int(pair[0])
                cons['atom2'] = int(pair[1])
                self.constraints.append(cons)

        #logger.info("Constraint cluster begin")
        #self.clusterConstraintOLD()
        #logger.info("Constraint cluster end")
    def clusterConstraint(self):
        setCluster=[]
        for cons in self.constraints:
            setCluster.append({cons['atom1'], cons['atom2']})

        clustering = True
        while clustering:
            clustering = False
            clusterSize=len(setCluster)
            clusterMask=[0] * clusterSize
            setClusterTmp=[]
            i=0
            for i in range(clusterSize):
                for j in range(i+1, clusterSize):
                    if clusterMask[i] == 1 or clusterMask[j] == 1 :
                        continue
                    if setCluster[i] & setCluster[j] != {}:
                        setClusterTmp.append(setCluster[i] | setCluster[j])
                        clusterMask[i] = 1
                        clusterMask[j] = 1
                        clustering = True
                        break
                    else:
                        setClusterTmp.append(setCluster[i])
                        setClusterTmp.append(setCluster[j])
            for idx, mask in enumerate(clusterMask):
                if mask == 0:
                    setClusterTmp.append(setCluster[idx])
            setCluster=setClusterTmp

        for set in setCluster:
            consList=[]
            for cons in self.constraints:
                if cons['atom1'] in set:
                    consList.append(cons)
                elif cons['atom2'] in set:
                    consList.append(cons)
            self.consCluster.append(consList)

    def clusterConstraintOLD(self):
        consList = []

        for cons in self.constraints:
            addcons = False
            atomI = cons['atom1']  # 0 based
            atomJ = cons['atom2']  # 0 based
            for consTemp in consList:
                atomItemp = consTemp['atom1']  # 0 based
                atomJtemp = consTemp['atom2']  # 0 based
                if atomItemp == atomI:
                    consList.append(cons)
                    addcons = True
                    break
                if atomItemp == atomJ:
                    consList.append(cons)
                    addcons = True
                    break
                if atomJtemp == atomI:
                    consList.append(cons)
                    addcons = True
                    break
                if atomJtemp == atomJ:
                    consList.append(cons)
                    addcons = True
                    break

            if not addcons:
                consList = []
                consList.append(cons)
                self.consCluster.append(consList)

        clustering = True

        while clustering:
            clustering = False
            clusterSize = len(self.consCluster)
            firstSave = 0
            secondSave = 0
            for first in range(clusterSize):
                for second in range(first + 1, clusterSize):
                    consListFirst = self.consCluster[first]
                    consListSecond = self.consCluster[second]
                    for consFirst in consListFirst:
                        consFirstI = consFirst['atom1']  # 0 based
                        consFirstJ = consFirst['atom2']  # 0 based
                        for consSecond in consListSecond:
                            consSecondI = consSecond['atom1']  # 0 based
                            consSecondJ = consSecond['atom2']  # 0 based
                            if consFirstI == consSecondI or consFirstI == consSecondJ or consFirstJ == consSecondI or consFirstJ == consSecondJ:
                                firstSave = first
                                secondSave = second
                                clustering = True
                                break
                        if clustering:
                            break
                    if clustering:
                        break
                if clustering:
                    break

            if clustering:
                consListFirst = self.consCluster[firstSave]
                consListSecond = self.consCluster[secondSave]
                for cons in consListSecond:
                    consListFirst.append(cons)
                del self.consCluster[secondSave]

class Restraint():
    def __init__(self):
        self.lines = []
        self.restraints = []

    def parse(self):
        for line in self.lines:
            if 'POSRES' in line:
                strs = re.split('\(|\)', line)
                if len(strs) != 3:
                    logger.warning("Incorrect restraint format: " + line)
                    continue
                restraint = {}
                restraint['type'] = int(strs[0].split('=')[1])
                restraint['atom1'] = int(strs[2].split()[0])
                self.restraints.append(restraint)

class VirtualSite():
    def __init__(self):
        self.lines = []
        self.vSites = []

    def parse(self):
        siteType=None
        vsiten_list=[]
        for line in self.lines:
            secName=getSecName(line, level=2)
            if secName:
                section=secName
                strs=line.split()
                siteType=strs[2].strip(':')
            # VSITEN has to treat differently
            if 'VSITEN' in line:
                vsiten_list.append(line)
            elif 'VSITE' in line:
                strs=re.split('\(|\)', line)
                if len(strs)!=3:
                    logger.warning("Incorrect Virtual Site format: " + line)
                    continue
                vS={}
                tStrs=strs[0].split('=')
                vS['type'] = int(tStrs[1])
                iStrs=strs[2].split()
                vS['idx'] = iStrs
                self.vSites.append(vS)

        if len(vsiten_list) > 0:
            vsiten_dcit = {}
            for line in vsiten_list:
                strs = re.split('\(|\)', line)
                if len(strs) != 3:
                    logger.warning("Incorrect Virtual Site format: " + line)
                    continue
                tStrs = strs[0].split('=')
                iStrs = strs[2].split()
                vAtomStr = iStrs[0]
                if vAtomStr not in vsiten_dcit:
                    vsiten_dcit[vAtomStr] = {'type': int(tStrs[1]), 'idx':[iStrs[1]]}
                else:
                    vsiten_dcit[vAtomStr]['idx'].append(iStrs[1])
            for key, value in vsiten_dcit.items():
                vS = {}
                vS['type'] = value['type']
                vS['idx'] = [key] + value['idx']
                self.vSites.append(vS)
