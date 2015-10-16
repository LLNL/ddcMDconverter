__author__ = 'zhang30'


class Obj:
    def __init__(self):
        self.header=""

    def toObj(self, filename, comPDB):
        totAtmNum=comPDB.getTotAtmNum()

        outFh=open(filename, "w")
        self.getHeader(totAtmNum)
        outFh.write(self.header)

        zeroV=0.0

        for molPDB in comPDB.molList:
            lastResID=len(molPDB.resList)-1
            for resID, resPDB in enumerate(molPDB.resList):
                lastGrpID=len(resPDB.grpList)-1
                for grpID, grpPDB in enumerate(resPDB.grpList):
                    for atmPDB in grpPDB.atmGrpList:
                        coor=atmPDB.coor
                        if resID==0 and grpID==0:
                            speciename=resPDB.resName+"n"+atmPDB.name
                        elif resID==lastResID and grpID==lastGrpID:
                            speciename=resPDB.resName+"c"+atmPDB.name
                        else:
                            speciename=resPDB.resName+"x"+atmPDB.name

                        outLine="%14d ATOM %11s group %21.13e %21.13e %21.13e %21.13e %21.13e %21.13e\n" \
                                 %(atmPDB.gid, speciename, coor.x, coor.y, coor.z, zeroV, zeroV, zeroV)
                        outFh.write(outLine)

    def getHeader(self, totAtmNum):
        self.header=self.header+"particle FILEHEADER {type=MULTILINE; datatype=VARRECORDASCII; checksum=NONE;\n"
        self.header=self.header+"loop=0; time=0.000000;\n"
        self.header=self.header+"nfiles=1; nrecord="+str(totAtmNum)+"; nfields=10;\n"
        self.header=self.header+"field_names=id class type group rx ry rz vx vy vz;\n"
        self.header=self.header+"field_types=u s s s f f f f f f;\n"
        self.header=self.header+"h=    64.00000000000000      0.00000000000000      0.00000000000000\n"
        self.header=self.header+"       0.00000000000000     64.00000000000000      0.00000000000000\n"
        self.header=self.header+"       0.00000000000000      0.00000000000000     64.00000000000000;\n"
        self.header=self.header+"groups = group ;\n"
        self.header=self.header+"types = ATOM ;\n"
        self.header=self.header+"} \n\n"
