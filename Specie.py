__author__ = 'zhang30'

import Charmm

class Specie:
    def __init__(self, charmmtop, compdb):
        self.charmmTop=charmmtop
        self.comPDB=compdb

    def toSpeData(self, filename):
        aaList={}
        nterList={}
        cterList={}

        for molPDB in self.comPDB.molList:
            lastResID=len(molPDB.resList)-1
            for resID, resPDB in enumerate(molPDB.resList):
                if resID==0:
                    nterList[resPDB.resName] = 1
                elif resID==lastResID:
                    cterList[resPDB.resName] = 1
                else:
                    aaList[resPDB.resName] = 1

        aaKeys=aaList.keys()
        nterKeys=nterList.keys()
        cterKeys=cterList.keys()

        sysOutLine=""
        defOutLine=""

        for resTop in self.charmmTop.resTopList:
            resname=resTop.resName

            if resname in aaKeys:
                for atmTop in resTop.atomList:
                    speciename=resname+"x"+atmTop.atmName
                    sysStr="   %11s \n" % speciename
                    sysOutLine=sysOutLine+sysStr
                    defStr="%11s SPECIES { type = ATOM ; charge = %f ; mass = %f M_p ; }\n"\
                           % (speciename, atmTop.charge, atmTop.atmTpyeTop.mass)
                    defOutLine=defOutLine+defStr

            if resname in nterKeys:
                for atmTop in resTop.atomList:
                    speciename=resname+"n"+atmTop.atmName
                    sysStr="   %11s \n" % speciename
                    sysOutLine=sysOutLine+sysStr
                    defStr="%11s SPECIES { type = ATOM ; charge = %f ; mass = %f M_p ; }\n"\
                           % (speciename, atmTop.charge, atmTop.atmTpyeTop.mass)
                    defOutLine=defOutLine+defStr

            if resname in cterKeys:
                for atmTop in resTop.atomList:
                    speciename=resname+"c"+atmTop.atmName
                    sysStr="   %11s \n" % speciename
                    sysOutLine=sysOutLine+sysStr
                    defStr="%11s SPECIES { type = ATOM ; charge = %f ; mass = %f M_p ; }\n"\
                           % (speciename, atmTop.charge, atmTop.atmTpyeTop.mass)
                    defOutLine=defOutLine+defStr

        outFh=open(filename, "w")

        outFh.write("system SYSTEM\n{species =\n")
        outFh.write(sysOutLine)
        outFh.write("   ;\n}\n\n")
        outFh.write(defOutLine)


