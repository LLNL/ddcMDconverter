__author__ = 'zhang30'

from ddcmdconverter.CharmmTop import CharmmTop
from ddcmdconverter.Pdb import Coor
import re
import struct

class PDBItem:
    def __init__(self):
        self.atomName=""
        self.resName=""
        self.resid=0
        self.fileID=0
        self.x=0.0
        self.y=0.0
        self.z=0.0

class Obj:
    def __init__(self):
        self.header=""
        self.datatype=""
        # BOX
        self.lx = 0
        self.ly = 0
        self.lz = 0
        self.nfiles=0
        self.nrecord=0
        self.nfield=0
        self.fieldName=[]
        self.species=[]


    def toObj(self, args, comPDB):
        totAtmNum=comPDB.getTotAtmNum()

        filename=args.objfile
        outFh=open(filename, "w")
        self.getHeader(totAtmNum, args)
        outFh.write(self.header)

        zeroV=0.0

        for molPDB in comPDB.molList:
            lastResID=len(molPDB.resList)-1
            for resID, resPDB in enumerate(molPDB.resList):
                lastGrpID=len(resPDB.grpList)-1
                isStdAA=CharmmTop.ResTop.isStdAA(resPDB.resName)
                for grpID, grpPDB in enumerate(resPDB.grpList):
                    for atmPDB in grpPDB.atmGrpList:
                        coor=atmPDB.coor
                        if isStdAA == 1:
                            if resID==0:
                                speciename=resPDB.resName+"n"+atmPDB.name
                            elif resID==lastResID:
                                speciename=resPDB.resName+"c"+atmPDB.name
                            else:
                                speciename=resPDB.resName+"x"+atmPDB.name
                        else:
                            speciename=resPDB.resName+"x"+atmPDB.name
                        outLine="%14d ATOM %11s group %21.13e %21.13e %21.13e %21.13e %21.13e %21.13e\n" \
                                 %(atmPDB.gid, speciename, coor.x, coor.y, coor.z, zeroV, zeroV, zeroV)
                        outFh.write(outLine)

    def getHeader(self, totAtmNum, args):
        self.header=self.header+"particle FILEHEADER {type=MULTILINE; datatype=VARRECORDASCII; checksum=NONE;\n"
        self.header=self.header+"loop=0; time=0.000000;\n"
        self.header=self.header+"nfiles=1; nrecord="+str(totAtmNum)+"; nfields=10;\n"
        self.header=self.header+"field_names=id class type group rx ry rz vx vy vz;\n"
        self.header=self.header+"field_types=u s s s f f f f f f;\n"
        self.header=self.header+"h=     "+ str(args.x)+"      0.00000000000000      0.00000000000000\n"
        self.header=self.header+"       0.00000000000000     "+ str(args.y)+"       0.00000000000000\n"
        self.header=self.header+"       0.00000000000000      0.00000000000000     "+ str(args.z)+" ;\n"
        self.header=self.header+"groups = group ;\n"
        self.header=self.header+"types = ATOM ;\n"
        self.header=self.header+"} \n\n"


    def toPDB(self, args):
        resMask=(65535<<16)

        prtflg=0
        boxflg=0
        boxCount=0
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
                            print("Specie name ", strs[3], " dosen't split properly")
                        else:
                            resName=names[0]
                            atmName=names[1]

                        x=float(strs[rxIndex])
                        y=float(strs[ryIndex])
                        z=float(strs[rzIndex])

                        coor=Coor(x,y,z)
                        #coor=Obj.image(args, coor)
                        newcoor=coor
                        #if first !=0:
                        #    newcoor=Obj.reduceImage(args, coor, oldcoor)
                        #    if newcoor==-1:
                        #        print( "Fail to shift the atom coordinates")
                        #        print( line)
                        #        print( "Use the original coordinates")
                        #        newcoor=coor

                        #first=1
                        fileID=fileID+1
                        #oldcoor=newcoor

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

                if line[0:2]=='h=':
                    boxflg=1

                if boxflg==1:
                    strs = line.split()
                    if boxCount==0:
                        args.x=float(strs[1])
                    elif boxCount==1:
                        args.y=float(strs[boxCount])
                    elif boxCount==2:
                        endstr=strs[boxCount]
                        estrs=endstr.split(";")
                        args.z=float(estrs[0])
                        boxflg = 0

                    boxCount=boxCount+1


    def toMartiniPDB(self, args, atomNameMapCollection):
        resMask=(65535<<16)

        prtflg=0
        boxflg=0
        boxCount=0
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

        oldmolres=-1
        resid=0
        with open(args.objfile, "r") as f:
            for line in f:
                if prtflg==1:
                    strs=line.split()
                    if len(strs)>=fieldSize:
                        gid=int(strs[gidIndex])
                        molres=gid>>16
                        if molres!=oldmolres:
                            resid=resid+1
                            oldmolres=molres
                            first=1

                        names=strs[typeIndex].replace('x', ' ').replace('n', ' ').replace('c', ' ').split()
                        if len(names) !=2:
                            print ("Specie name ", strs[3], " dosen't split properly")
                        else:
                            resName=names[0]
                            atmName=names[1]
                            if resName in atomNameMapCollection:
                                atomNameMap=atomNameMapCollection[resName]
                                if atmName in atomNameMap:
                                    residueAtom=atomNameMap[atmName]
                                    residueAtomsplits=residueAtom.split(':')
                                    resName=residueAtomsplits[0]
                                    resid=int(residueAtomsplits[1])
                                    atmName=residueAtomsplits[2]

                        x=float(strs[rxIndex])
                        y=float(strs[ryIndex])
                        z=float(strs[rzIndex])

                        coor=Coor(x,y,z)
                        #coor=Obj.image(args, coor)
                        newcoor=coor
                        #if first ==0:
                        #    newcoor=Obj.reduceImage(args, coor, oldcoor)
                        #    if newcoor==-1:
                        #        print( "Fail to shift the atom coordinates")
                        #        print( line)
                        #        print( "Use the original coordinates")
                        #        newcoor=coor

                        first=0
                        fileID=fileID+1
                        #oldcoor=newcoor
                        if fileID==100000:
                            fileID=0

                        if resid==10000:
                            resid=0

                        #outLine = "ATOM%7d %-4s %-4s%5d    %8.3f%8.3f%8.3f\n"  \
                        #            % (fileID, atmName, resName, resid, newcoor.x, newcoor.y, newcoor.z)
                        if len(atmName)==4:
                            outLine="ATOM%7d %-4s %-4s%5d    %8.3f%8.3f%8.3f\n" \
                                    % (fileID, atmName, resName, resid, newcoor.x, newcoor.y, newcoor.z)
                        else:
                            outLine="ATOM%7d  %-3s %-4s%5d    %8.3f%8.3f%8.3f\n" \
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

                if line[0:2]=='h=':
                    boxflg=1

                if boxflg==1:
                    strs = line.split()
                    if boxCount==0:
                        args.x=float(strs[1])
                    elif boxCount==1:
                        args.y=float(strs[boxCount])
                    elif boxCount==2:
                        endstr=strs[boxCount]
                        estrs=endstr.split(";")
                        args.z=float(estrs[0])
                        boxflg = 0
                        outLine = "REMARK CONVERTED FROM %s BY ddcMDconvertor\n" % (args.objfile)
                        outFh.write(outLine)
                        outLine = "CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00 P 1           1\n" % (args.x, args.y, args.z)
                        outFh.write(outLine)

                    boxCount=boxCount+1

    def parseObj(self, object):
        replaceObj=object.replace("\n", " ")
        splitObjs=re.split('{|}', replaceObj)

        if len(splitObjs)<3:
            print("Object doesn't enclose with a pair of {}")

        contentList=splitObjs[1].split(";")

        contentDict={}

        for content in contentList:
            pair=content.split("=")
            if len(pair) == 2:
                contentDict[pair[0].strip()]=pair[1].strip()
            #else:
            #   print ("Wrong pair in object:", content)

        return contentDict

    def parseHeader(self, args):

        header=""

        with open(args.objfile, "r") as f:
            for line in f:
                header = header + line
                if line[0] == '}':
                    break

        headerDict=self.parseObj(header)

        requiredKeys=['datatype', 'h', 'field_names', 'nfields', 'nfiles', 'nrecord', 'species']
        for key in requiredKeys:
            if key not in headerDict:
                print ("Missing required key in header :", key)
                exit(-1)

        self.datatype = headerDict['datatype']

        boxSplit=headerDict['h'].split()
        if len(boxSplit)<9:
            print ("Item of box values should be 9 but it is ", len(boxSplit))
            exit(-1)
        self.lx = float(boxSplit[0])
        self.ly = float(boxSplit[4])
        self.lz = float(boxSplit[8])

        self.nfield = int(headerDict['nfields'])
        self.nfiles = int(headerDict['nfiles'])
        self.nrecord = int(headerDict['nrecord'])
        self.fieldName = headerDict['field_names'].split()
        self.species = headerDict['species'].split()


    def toBinaryMartiniPDB(self, args, atomNameMapCollection):

        # index for fields
        requiredFields=['id', 'species', 'rx', 'ry', 'rz']
        for field in requiredFields:
            if field not in self.fieldName:
                print ("Missing field name: ", field)
                exit(-1)

        # find out the offset for the first file
        contents=""
        skipCount=0
        isContent=False
        with open(args.objfile, "rb") as f:
            for line in f:
                contents=contents+line
                if isContent:
                    skipCount=skipCount+1
                    if skipCount==2:
                        # skip two empty line and start to read in binary data
                        break
                if line[0] == '}':
                    isContent = True

        offset=len(contents.encode('utf-8'))

        fileNames=[]
        fileBaseName=args.objfile.split("#")[0]
        for i in range(self.nfiles):
            numbering=str(i).zfill(6)
            newfileName=fileBaseName+"#"+numbering
            fileNames.append(newfileName)

        atomCount = 0
        oldmolres=-1
        resid=0
        fileID =0
        pdbItemList=[]

        for idx, filename in enumerate(fileNames):
            with open(filename, "rb") as f:
                if idx==0: # the first file has header
                    f.seek(offset)
                data = f.read(8)
                while data:
                    # get data from binary file
                    gid = struct.unpack('<Q', data)[0]
                    data = f.read(4)
                    specie = struct.unpack("<L", data)[0]
                    specieName=self.species[specie]
                    data = f.read(4)
                    rx = struct.unpack("<f", data)[0]
                    data = f.read(4)
                    ry = struct.unpack("<f", data)[0]
                    data = f.read(4)
                    rz = struct.unpack("<f", data)[0]

                    #create PDBItem object
                    pdbItem=PDBItem()
                    molres = gid >> 16
                    if molres != oldmolres:
                        resid = resid + 1
                        oldmolres = molres

                    names = specieName.replace('x', ' ').replace('n', ' ').replace('c', ' ').split()
                    if len(names) != 2:
                        print ("Specie name ", specieName, " dosen't split properly")
                    else:
                        resName = names[0]
                        atmName = names[1]
                        if resName in atomNameMapCollection:
                            atomNameMap = atomNameMapCollection[resName]
                            if atmName in atomNameMap:
                                residueAtom = atomNameMap[atmName]
                                residueAtomsplits = residueAtom.split(':')
                                resName = residueAtomsplits[0]
                                resid = int(residueAtomsplits[1])
                                atmName = residueAtomsplits[2]

                    fileID = fileID + 1

                    if fileID == 100000:
                        fileID = 1

                    if resid == 10000:
                        resid = 1

                    pdbItem.resid = resid
                    pdbItem.fileID = fileID
                    pdbItem.resName=resName
                    pdbItem.atomName=atmName
                    pdbItem.x = rx
                    pdbItem.y = ry
                    pdbItem.z = rz
                    pdbItemList.append(pdbItem)

                    #Goto next record
                    atomCount = atomCount + 1
                    data = f.read(8)

        self.toPDBbyItem(args, pdbItemList)

    def toPDBbyItem(self, args, pdbItemList):

        outFh = open(args.pdbfile, "w")
        outLine = "REMARK CONVERTED FROM %s BY ddcMDconvertor\n" % (args.objfile)
        outFh.write(outLine)
        outLine = "CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00 P 1           1\n" % (self.lx, self.ly, self.lz)
        outFh.write(outLine)

        for pdbItem in pdbItemList:
            if len(pdbItem.atomName) == 4:
                outLine = "ATOM%7d %-4s %-4s%5d    %8.3f%8.3f%8.3f\n" \
                          % (pdbItem.fileID, pdbItem.atomName, pdbItem.resName, pdbItem.resid, pdbItem.x, pdbItem.y, pdbItem.z)
            else:
                outLine = "ATOM%7d  %-3s %-4s%5d    %8.3f%8.3f%8.3f\n" \
                          % (pdbItem.fileID, pdbItem.atomName, pdbItem.resName, pdbItem.resid, pdbItem.x, pdbItem.y, pdbItem.z)
            outFh.write(outLine)

        outFh.write("END")


    def toRestraint(self, args, restraintMap):
        atmMask = 255

        prtflg=0
        boxflg=0
        boxCount=0
        first=0
        fileID=0

        # index for fields
        fieldSize=0
        gidIndex=0
        typeIndex=0
        rxIndex=0
        ryIndex=0
        rzIndex=0



        oldmolres=-1
        resid=0

        count=0
        entryLine = "restraint RESTRAINTLIST{\n"
        listLine = "  restraintList="
        outLine = ""
        with open(args.objfile, "r") as f:
            for line in f:
                if prtflg==1:
                    strs=line.split()
                    if len(strs)>=fieldSize:
                        gid=int(strs[gidIndex])
                        molres=gid>>16
                        if molres!=oldmolres:
                            resid=resid+1
                            oldmolres=molres

                        names=strs[typeIndex].replace('x', ' ').replace('n', ' ').replace('c', ' ').split()
                        if len(names) !=2:
                            print ("Specie name ", strs[3], " dosen't split properly")
                        else:
                            resName=names[0]
                            atmName=names[1]

                        if resName in restraintMap:

                            restraintList=restraintMap[resName]

                            atomI=gid&atmMask

                            for restraint in restraintList:
                                resAtomI = restraint['ai']-1 # 0 based
                                if resAtomI==atomI:
                                    func = restraint['func']
                                    fcx = restraint['fcx']
                                    if fcx>0:
                                        fcx=1
                                    fcy = restraint['fcy']
                                    if fcy>0:
                                        fcy=1
                                    fcz = restraint['fcz']
                                    if fcz>0:
                                        fcz=1
                                    x = float(strs[rxIndex])
                                    y = float(strs[ryIndex])
                                    z = float(strs[rzIndex])
                                    listLine = listLine+"r_" + str(count) + " "
                                    outLine = outLine+ "r_" + str(count) + " RESTRAINTPARMS{gid="+str(gid)+"; atomI=" + str(atomI) + "; func=" + str(func) \
                                              + "; fcx=" + str(fcx) + "; fcy=" + str(fcy) + "; fcz=" + str(fcz) \
                                              + "; x0=" + str(x/args.x) + "; y0=" + str(y/args.y) + "; z0=" + str(z/args.z) \
                                              + "; kb= 1.0 kJ*mol^-1*nm^-2; }\n"
                                    count = count + 1

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


                if line[0:2]=='h=':
                    boxflg=1

                if boxflg==1:
                    strs = line.split()
                    if boxCount==0:
                        args.x=float(strs[1])
                    elif boxCount==1:
                        args.y=float(strs[boxCount])
                    elif boxCount==2:
                        endstr=strs[boxCount]
                        estrs=endstr.split(";")
                        args.z=float(estrs[0])
                        boxflg = 0
                        entryLine = entryLine + "  xbox = " + str(args.x) + ";\n"
                        entryLine = entryLine + "  ybox = " + str(args.y) + ";\n"
                        entryLine = entryLine + "  zbox = " + str(args.z) + ";\n"

                    boxCount=boxCount+1

        entryLine=entryLine+listLine
        entryLine = entryLine + ";\n"
        entryLine = entryLine + "}\n\n"

        outLine = outLine +"\n"

        outFh=open(args.resfile, "w")
        outFh.write(entryLine)
        outFh.write(outLine)

    @staticmethod
    def reduceImage(args, coor, oldcoor):
        cutoff2=args.cutoff*args.cutoff
        d2=coor.dist2(oldcoor)
        if d2<cutoff2:
            return coor

        dx=args.x
        dy=args.y
        dz=args.z

        # search for -1, 0, 1
        for i in range(-1, 2, 1):
            for j in range(-1, 2, 1):
                for k in range(-1, 2, 1):
                    xnew=coor.x+i*dx
                    ynew=coor.y+j*dy
                    znew=coor.z+k*dz
                    newcoor=Coor(xnew, ynew, znew)

                    d2=newcoor.dist2(oldcoor)
                    if d2<cutoff2:
                        return newcoor

        return -1

    @staticmethod
    def image(args, coor):
        coor.x = coor.x % args.x
        coor.y = coor.y % args.y
        coor.z = coor.z % args.z

        return coor



