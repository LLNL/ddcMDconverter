import re


class Index:
    def __init__(self):
        self.ids=[]

class AtmPsf:
    def __init__(self):
        self.id=0
        self.atmName=""
        self.resName=""

class Psf:
    def __init__(self):
        self.atomList=[]
        self.angleList=[]
        self.diheList=[]

    def parse(self, filename):
        atomRe=re.compile('!NATOM')
        angleRe=re.compile('!NTHETA: angles')
        diheRe=re.compile('!NPHI: dihedrals')

        with open(filename, "r") as f:
            iterLine=iter(f)
            for line in iterLine:
                if atomRe.search(line):
                    strs=line.split()
                    numAtom=int(strs[0])
                    count=0
                    while count<numAtom:
                        linestr=next(iterLine)
                        strs=linestr.split()
                        atmPsf=AtmPsf()
                        atmPsf.id=int(strs[0])
                        atmPsf.atmName=strs[4]
                        atmPsf.resName = strs[4]
                        self.atomList.append(atmPsf)
                        count=count+1
                if angleRe.search(line):
                    strs=line.split()
                    numAngle=int(strs[0])
                    count=0
                    while count<numAngle:
                        linestr=next(iterLine)
                        strs=linestr.split()
                        inc=len(strs)/3
                        count = count + inc
                        for i in range(inc):
                            ind=Index()
                            self.angleList.append(ind)
                            for j in range(3):
                                ind.ids.append(int(strs[i*3+j]))
                if diheRe.search(line):
                    strs=line.split()
                    numDihe=int(strs[0])
                    count=0
                    while count < numDihe:
                        linestr=next(iterLine)
                        strs=linestr.split()
                        inc=len(strs)/4
                        count = count + inc
                        for i in range(inc):
                            ind=Index()
                            self.diheList.append(ind)
                            for j in range(4):
                                ind.ids.append(int(strs[i*4+j]))

    def printAngle(self):
         for i, angle in enumerate(self.angleList):
             print("Angle: ", i+1,  end=' ')
             for j, id in enumerate(angle.ids):
                 print(self.atomList[id-1].atmName, end='')
                 if j !=2:
                     print(" - ", end='')
             print(" ")

    def printDihe(self):
        for i, dihe in enumerate(self.diheList):
            print("Dihe: ", i + 1, end=' ')
            for j, id in enumerate(dihe.ids):
                print(self.atomList[id - 1].atmName, end='')
                if j != 3:
                    print(" - ", end='')
            print(" ")
