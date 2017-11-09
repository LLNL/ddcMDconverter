__author__ = 'zhang30'

import sys

class LineTop:
    @staticmethod
    def getLineType(line):
        if line[0]=='\0' or line[0]==' ' or line[0]=='!':
            return "SKIP"

        if line[:4]=="MASS":
            return "MASS"
        elif line[:4]=="RESI" or line[:4]=="PRES" :
            return "RESI"
        elif line[:5]=="GROUP":
            return "GROUP"
        elif line[:4]=="ATOM":
            return "ATOM"
        elif line[:4]=="BOND":
            return "BOND"
        elif line[:4]=="IMPR":
            return "IMPR"
        elif line[:4]=="CMAP":
            return "CMAP"
        elif line[:2]=="IC":
            return "IC"
        elif line[:5]=="DONOR":
            return "DONOR"
        elif line[:8]=="ACCEPTOR":
            return "ACCEPTOR"
        elif line[:3]=="END":
            return "END"
        else:
            return "SKIP"

    @staticmethod
    def removeComment(line):
        if '!' in line:
            i=line.index('!')-1
            return line[:i]
        else:
            return line

class AtomTypeTop:
    def __init__(self):
        self.atmType=""
        self.element=""
        self.mass=0.0
    def parse(self, line):
        strs=line.split()
        self.atmType=strs[2]
        self.mass=float(strs[3])
        self.element=strs[4]


class AtomTop:
    def __init__(self):
        self.atmID=0
        self.atmName=""
        self.charge=0.0
        self.atmTpyeTop=0

    @staticmethod
    def getAtmType(atmType, atmTypeTopList):
        for atmTypeTop in atmTypeTopList:
            if atmType== atmTypeTop.atmType :
                return atmTypeTop

        return 0
        print "Cannot find atom type for ", atmType

    def parse(self, line, count, atmTypeTopList):
        strs=line.split()
        self.atmID=count
        self.atmName=strs[1]
        self.charge=float(strs[3])
        atmType=strs[2]
        self.atmTpyeTop=AtomTop.getAtmType(atmType, atmTypeTopList)


class GroupTop:
    def __init__(self):
        self.grpID=0
        self.grpAtoms=[]

class BondTop:
    def __init__(self, atmI, atmJ, c=0.0):
        self.atomI=atmI
        self.atomJ=atmJ
        self.const=c

class DihTop:
    def __init__(self, atmI, atmJ, atmK, atmL):
        self.atomI=atmI
        self.atomJ=atmJ
        self.atomK=atmK
        self.atomL=atmL

class ICTop:
    def __init__(self):
        self.dihTop=0
        self.kconst1=0.0
        self.angle1=0.0
        self.torsion=0.0
        self.angle2=0.0
        self.kconst2=0.0
    def parse(self, line):
        strs=line.split()
        if len(strs) <10:
            print "IC field number is wrong ", line
        self.dihTop=DihTop(strs[1],strs[2], strs[3],strs[4])
        self.kconst1=float(strs[5])
        self.angle1=float(strs[6])
        self.torsion=float(strs[7])
        self.angle2=float(strs[8])
        self.kconst2=float(strs[9])

class ResTop:
    def __init__(self):
        self.resID=0
        self.resType=0
        self.resName=""
        self.charge=0.0

        self.atomList=[]
        self.groupList=[]
        self.bondList=[]
        self.imprList=[]
        self.cmapList=[]
        self.donorList=[]
        self.accList=[]
        self.icList=[]

    @staticmethod
    def isStdAA(resName):
        stdAA=["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HSD", "HSE", "HSP",
                "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        if resName in stdAA:
            return 1
        else:
            return 0

    def parseName(self, line, count):
        strs=line.split()
        self.resID=count
        self.resName=strs[1]
        self.resType=ResTop.isStdAA(self.resName)
        self.charge=float(strs[2])

    def parse(self, strList, count, atmTypeTopList):
        grpCount=0
        atmCount=0
        for line in strList:
            lt=LineTop.getLineType(line)
            if lt=="RESI":
                self.parseName(line, count)
            elif lt=="GROUP":
                grpTop=GroupTop()
                grpTop.grpID=grpCount
                grpCount=grpCount+1
                atmCount=0
                self.groupList.append(grpTop)
            elif lt=="ATOM":
                atmTop=AtomTop()
                atmTop.parse(line, atmCount, atmTypeTopList)
                atmCount=atmCount+1
                self.atomList.append(atmTop)
                #print atmTop.atmID, atmTop.atmName, self.resName
                if 'grpTop' not in locals():  # If there is no group before ATOM
                    grpTop = GroupTop()
                    grpTop.grpID = grpCount
                    grpCount=grpCount+1
                    atmCount=0
                grpTop.grpAtoms.append(atmTop)
            elif lt=="BOND":
                ncline=LineTop.removeComment(line)
                strs=ncline.split()
                numStr=len(strs)-1
                if numStr%2 != 0 :
                    print "Odd bond atom number", line
                nPair=numStr/2
                for i in range(nPair):
                    bndTop=BondTop(strs[2*i+1], strs[2*i+2])
                    self.bondList.append(bndTop)
            elif lt=="IMPR":
                ncline=LineTop.removeComment(line)
                strs=ncline.split()
                numStr=len(strs)-1
                if numStr%4 != 0 :
                    print "Odd IMPR atom number", line
                nPair=numStr/4
                for i in range(nPair):
                    dihTop=DihTop(strs[4*i+1], strs[4*i+2],strs[4*i+3], strs[4*i+4])
                    self.imprList.append(dihTop)
            elif lt=="CMAP":
                ncline=LineTop.removeComment(line)
                strs=ncline.split()
                numStr=len(strs)-1
                if numStr%4 != 0 :
                    print "Odd CMAP atom number", line
                nPair=numStr/4
                for i in range(nPair):
                    dihTop=DihTop(strs[4*i+1], strs[4*i+2],strs[4*i+3], strs[4*i+4])
                    self.cmapList.append(dihTop)
            elif lt=="DONOR":
                ncline=LineTop.removeComment(line)
                strs=ncline.split()
                numStr=len(strs)-1
                for i in range(numStr):
                    self.donorList.append(strs[i+1])
            elif lt=="ACCEPTOR":
                ncline=LineTop.removeComment(line)
                strs=ncline.split()
                numStr=len(strs)-1
                for i in range(numStr):
                    self.accList.append(strs[i+1])
            elif lt=="IC":
                icTop=ICTop()
                icTop.parse(line)
                self.icList.append(icTop)


class CharmmTop:
    def __init__(self):
        self.atmTypeTopList=[]
        self.resTopList=[]

    def parse(self, filename):
        firstRes=0

        massList=[]
        resiLists=[]
        aResiList=[]

        with open(filename, "r") as f:
            for line in f:
                lt=LineTop.getLineType(line)
                if lt=="SKIP":
                    continue
                elif lt== "MASS":
                    massList.append(line)
                else:
                    if lt=="RESI" or lt=="END":
                        if firstRes!=0:
                            resiLists.append(aResiList)
                            if lt=="END":
                                break
                        aResiList=[]
                        firstRes=1
                    aResiList.append(line)

        for line in massList:
            #sys.stdout.write(line)
            atomTypeTop=AtomTypeTop()
            atomTypeTop.parse(line)
            self.atmTypeTopList.append(atomTypeTop)

        for count, aResiList in enumerate(resiLists):
            resTop=ResTop()
            resTop.parse(aResiList, count, self.atmTypeTopList)
            #if resTop.resName == 'GLUP':
            #    print resTop.resName
            self.resTopList.append(resTop)

    def findResiParm(self, resName):
        for resTop in self.resTopList:
            if resName==resTop.resName:
                return resTop

        print "Cannot find residue name: ", resName, " in CharmmTop.resTopList"
        return 0