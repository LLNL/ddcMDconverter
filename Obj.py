__author__ = 'zhang30'

import Charmm
from Pdb import Coor

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
                isStdAA=Charmm.ResTop.isStdAA(resPDB.resName)
                for grpID, grpPDB in enumerate(resPDB.grpList):
                    for atmPDB in grpPDB.atmGrpList:
                        coor=atmPDB.coor
                        if isStdAA == 1:
                            if resID==0 and grpID==0:
                                speciename=resPDB.resName+"n"+atmPDB.name
                            elif resID==lastResID and grpID==lastGrpID:
                                speciename=resPDB.resName+"c"+atmPDB.name
                            else:
                                speciename=resPDB.resName+"x"+atmPDB.name
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


    def toPDB(self, args):
        resMask=(65535<<16)

        prtflg=0
        first=0
        fileID=0

        # index for fields
        fieldSize=0
        gidIndex=0
        typeIndex=0
        rxIndex=0
        ryIndex=0
        rzIndex=0

        outFh=open(args.pdbfile, "w")

        with open(args.objfile, "r") as f:
            for line in f:
                if prtflg==1:
                    strs=line.split()
                    if len(strs)>=fieldSize:
                        gid=int(strs[gidIndex])
                        resid=(gid&resMask)>>16
                        resid=resid+1

                        names=strs[typeIndex].replace('x', ' ').replace('n', ' ').replace('c', ' ').split()
                        if len(names) !=2:
                            print "Specie name ", strs[3], " dosen't split properly"
                        else:
                            resName=names[0]
                            atmName=names[1]

                        x=float(strs[rxIndex])
                        y=float(strs[ryIndex])
                        z=float(strs[rzIndex])

                        coor=Coor(x,y,z)
                        newcoor=coor
                        if first !=0:
                            newcoor=Obj.reduceImage(args, coor, oldcoor)
                            if newcoor==-1:
                                print "Fail to shift the atom coordinates"
                                print line
                                print "Use the original coordinates"
                                newcoor=coor

                        first=1
                        fileID=fileID+1
                        oldcoor=newcoor

                        if len(atmName)==4:
                            outLine="ATOM%7d %-4s%4s %5d    %8.3f%8.3f%8.3f\n" \
                                    % (fileID, atmName, resName, resid, newcoor.x, newcoor.y, newcoor.z)
                        else:
                            outLine="ATOM%7d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n" \
                                    % (fileID, atmName, resName, resid, newcoor.x, newcoor.y, newcoor.z)

                        outFh.write(outLine)

                if line[0]=='}':
                    prtflg=1

                if line[0:11]=='field_names':
                    fstrs=line.split('=')
                    fnameLine=fstrs[1]
                    strs=fnameLine.split()
                    # id class type group rx ry rz vx vy vz
                    gidIndex=strs.index("id")
                    typeIndex = strs.index("type")
                    rxIndex = strs.index("rx")
                    ryIndex = strs.index("ry")
                    rzIndex = strs.index("rz")
                    fieldSize=len(strs)



    @staticmethod
    def reduceImage(args, coor, oldcoor):
        cutoff2=args.cutoff*args.cutoff
        d2=coor.dist2(oldcoor)
        if d2<cutoff2:
            return coor

        dx=args.x
        dy=args.y
        dz=args.z

        for i in range(-2, 1, 2):
            for j in range(-2, 1, 2):
                for k in range(-2, 1, 2):
                    xnew=coor.x+i*dx
                    ynew=coor.y+i*dy
                    znew=coor.z+i*dz
                    newcoor=Coor(xnew, ynew, znew)

                    d2=newcoor.dist2(oldcoor)
                    if d2<cutoff2:
                        return coor

        return -1