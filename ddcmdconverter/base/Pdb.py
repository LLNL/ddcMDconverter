__author__ = 'zhang30'

from ddcmdconverter.charmm.CharmmTop import ResTop

class Coor:
    def __init__(self, x, y, z):
        self.x=x
        self.y=y
        self.z=z

    def dist2(self, coor):
        dx=self.x-coor.x
        dy=self.y-coor.y
        dz=self.z-coor.z
        return (dx*dx+dy*dy+dz*dz)

class AtomPDB:
    def __init__(self):
        self.name=""
        self.gid=0
        self.coor=0

class GroupPDB:
    def __init__(self):
        self.grpID=0
        self.atmGrpList=[]

class ResPDB:
    def __init__(self):
        self.resID=0
        self.resName=""
        #self.nTer=0
        #self.cTer=0
        self.grpList=[]
        self.atmList=[]

    def isStdAA(self):
        stdAAList=["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        return self.resName in stdAAList

class MolPDB:
    def __init__(self):
        self.molID=0
        self.resList=[]
        self.isPortein=False
        self.proteinName=""

    def parse(self, molstrs):

        oldResName=""
        oldResID=-1
        for line in molstrs:
            if line[:7] == "PROTEIN":
                strs=line.split()
                if len(strs) > 2:
                    self.isPortein = True
                    self.proteinName = strs[2]
                continue

            if line[:4]=="ATOM" :
                atmID=int(line[4:11])
            elif line[:6]=="HETATM" :
                atmID=int(line[6:11])

            atmName=line[12:17].strip()
            resName=line[17:21].strip()

            resID=int(line[22:27], 16)
            x=float(line[30:38])
            y=float(line[38:46])
            z=float(line[46:54])

            if resID != oldResID or resName !=oldResName :
                resPDB=ResPDB()
                resPDB.resID=resID
                resPDB.resName=resName
                self.resList.append(resPDB)

            oldResID=resID
            oldResName=resName

            coor=Coor(x,y,z)
            atmPDB=AtomPDB()
            atmPDB.coor=coor
            atmPDB.name=atmName
            resPDB.atmList.append(atmPDB)

class ComPDB:
    def __init__(self):
        #self.totalAtm=0
        self.molList=[]

    def parse(self, args):
        molsList=[]
        aMol=[]

        oldResName = ""
        oldResID = -1
        isNewRes=False

        filename=args.pdbfile
        hasBox=False

        with open(filename, "r") as f:
            for line in f:
                if line[:4]=="ATOM":
                    resName = line[17:21].strip()
                    resID = int(line[22:27], 16)
                    if resID != oldResID or resName != oldResName:
                        isNewRes = True
                        oldResID = resID
                        oldResName = resName
                    else:
                        isNewRes = False

                if isNewRes:
                    isRes=ResTop.isResidue(resName)
                    if isRes==0:
                        if (len(aMol) > 0):
                            molsList.append(aMol)
                        aMol = []

                if line[:3]=="TER" or line[:3]=="END":
                    if(len(aMol)>0):
                        molsList.append(aMol)
                    aMol=[]
                elif line[:4]=="ATOM" or line[:6]=="HETATM" or line[:7]=="PROTEIN":
                    aMol.append(line)
                elif line[:6]=="CRYST1":
                    hasBox=self.getBoxSize(args, line)

        if hasBox:
            print("The box size is: ", args.x, args.y, args.z)
        else:
            print("Dosen't have box size info in pdb, use default box size.")

        for aMol in molsList:
            mol=MolPDB()
            mol.parse(aMol)
            self.molList.append(mol)

    def getTotAtmNum(self):
        totalAtm=0
        for molPDB in self.molList:
            for resPDB in molPDB.resList:
                totalAtm=totalAtm+len(resPDB.atmList)

        return totalAtm

    def getBoxSize(self, args, line):
        strs=line.split()
        if len(strs)>4:
            args.x = float(strs[1])
            args.y = float(strs[2])
            args.z = float(strs[3])
            return True
        return False

    def assignGid(self, charmmTop):
        #for aMol in self.molList:
        for molID, molPDB in enumerate(self.molList):
            lastResID=len(molPDB.resList)-1
            for resID, resPDB in enumerate(molPDB.resList):
                nTER=0
                if resID==0 and resPDB.isStdAA():
                    if resPDB.resName=="GLY":
                        nTER = charmmTop.findResiParm("GLYP")
                    elif resPDB.resName=="PRO":
                        nTER = charmmTop.findResiParm("PROP")
                    else:
                        nTER=charmmTop.findResiParm("NTER")
                cTER=0
                if resID==lastResID and resPDB.isStdAA():
                    cTER=charmmTop.findResiParm("CTER")

                curResiParm=charmmTop.findResiParm(resPDB.resName)

                oldGrpID=-10
                for atmNum, atmPDB in enumerate(resPDB.atmList):
                    (success, grpID, atmID)=ComPDB.findGrpAtmID(curResiParm, nTER, cTER, atmPDB)
                    if success<0:
                        print("Cannot find PDB atom name: ", atmPDB.name, " in residuce parm", curResiParm.resName)
                    atmPDB.gid=(molID<<32)+(resID<<16)+(grpID<<8)+atmID
                    if grpID !=oldGrpID:
                        grpPDB=GroupPDB()
                        grpPDB.grpID=grpID
                        resPDB.grpList.append(grpPDB)
                        oldGrpID=grpID

                    grpPDB.atmGrpList.append(atmPDB)


    @staticmethod
    def findGrpAtmID(resiParm, nTER, cTER, atmPDB):
        if nTER !=0:
            (success, grpID, atmID)=ComPDB.findAtmIDinGrp(nTER, atmPDB)
            if success>0:
                grpID=0
                return (success, grpID, atmID)

        if cTER !=0:
            (success, grpID, atmID)=ComPDB.findAtmIDinGrp(cTER, atmPDB)
            if success> 0:
                grpID=len(resiParm.groupList)-1
                return (success, grpID, atmID)

        (success, grpID, atmID)=ComPDB.findAtmIDinGrp(resiParm, atmPDB)
        return (success, grpID, atmID)


    @staticmethod
    def findAtmIDinGrp(resiParm, atmPDB):
        atmName=atmPDB.name
        #print("atmPDB.name=", atmName, " resiParm", resiParm.resName)
        for grpID, grpTop in enumerate(resiParm.groupList):
            for atmID, atmTop in enumerate(grpTop.grpAtoms):
                if atmTop.atmName==atmName:
                    success=1
                    return (success, grpID, atmID)

        return (-1, -1, -1)

    def findUniqRes(self):
        aaList={}
        nterList={}
        cterList={}

        for molPDB in self.molList:
            lastResID=len(molPDB.resList)-1
            for resID, resPDB in enumerate(molPDB.resList):
                if CharmmTop.ResTop.isStdAA(resPDB.resName):
                    if resID==0:
                        nterList[resPDB.resName] = 1
                    elif resID==lastResID:
                        cterList[resPDB.resName] = 1
                    else:
                        aaList[resPDB.resName] = 1
                else:
                    aaList[resPDB.resName] = 1

        aaKeys=aaList.keys()
        nterKeys=nterList.keys()
        cterKeys=cterList.keys()

        return (aaKeys, nterKeys, cterKeys)
