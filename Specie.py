__author__ = 'zhang30'


class Specie:
    def __init__(self, charmmtop, compdb):
        self.charmmTop=charmmtop
        self.comPDB=compdb

    def toSpeData(self, filename):
        (aaKeys, nterKeys, cterKeys) = self.comPDB.findUniqRes()

        sysOutLine=""
        defOutLine=""

        nTerTop=None
        glypTop=None
        propTop=None
        cTerTop=None
        newnTop=None

        for resTop in self.charmmTop.resTopList:
            resname=resTop.resName
            if resname == 'NTER':
                nTerTop=resTop
            if resname == 'GLYP':
                glypTop=resTop
            if resname == 'PROP':
                propTop=resTop
            if resname == 'CTER':
                cTerTop=resTop

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
                # Use the NTER/GLYP/PROP for the first group
                if resname == 'GLY':
                    newnTop =glypTop
                elif resname == 'PRO':
                    newnTop =propTop
                else:
                    newnTop = nTerTop

                for atmTop in newnTop.atomList:
                    speciename=resname+"n"+atmTop.atmName
                    sysStr="   %11s \n" % speciename
                    sysOutLine=sysOutLine+sysStr
                    defStr="%11s SPECIES { type = ATOM ; charge = %f ; mass = %f M_p ; }\n"\
                           % (speciename, atmTop.charge, atmTop.atmTpyeTop.mass)
                    defOutLine=defOutLine+defStr
                # Skip the first group
                iterGrps=iter(resTop.groupList)
                next(iterGrps)
                for grpTop in iterGrps:
                    for atmTop in grpTop.grpAtoms:
                        speciename=resname+"n"+atmTop.atmName
                        sysStr="   %11s \n" % speciename
                        sysOutLine=sysOutLine+sysStr
                        defStr="%11s SPECIES { type = ATOM ; charge = %f ; mass = %f M_p ; }\n"\
                               % (speciename, atmTop.charge, atmTop.atmTpyeTop.mass)
                        defOutLine=defOutLine+defStr

            if resname in cterKeys:
                # Skip the last group
                for grpTop in resTop.groupList[:-1]:
                    for atmTop in grpTop.grpAtoms:
                        speciename=resname+"c"+atmTop.atmName
                        sysStr="   %11s \n" % speciename
                        sysOutLine=sysOutLine+sysStr
                        defStr="%11s SPECIES { type = ATOM ; charge = %f ; mass = %f M_p ; }\n"\
                               % (speciename, atmTop.charge, atmTop.atmTpyeTop.mass)
                        defOutLine=defOutLine+defStr
                # use the CTER for the last group
                for atmTop in cTerTop.atomList:
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


