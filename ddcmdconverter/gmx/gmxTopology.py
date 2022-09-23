import logging
import re
from gmxParameter import getSecName

LOGLEVEL =2
LOG_FMT = '%(asctime)s - %(name)s:%(funcName)s:%(lineno)s - %(levelname)s - %(message)s'
logger = logging.getLogger(__name__)
logger.setLevel(LOGLEVEL)

sh = logging.StreamHandler()
sh.setLevel(LOGLEVEL)
sh.setFormatter(logging.Formatter(LOG_FMT))
logger.addHandler(sh)

class MolBlock():
    def __init__(self):
        self.moltypeID=None
        self.moltypeName=''
        self.numMols=0
        self.numPosResxA=0
        #self.posResxA =[]
        self.numPosResxB=0
        #self.posResxB =[]
        self.lines=[]

    def parse(self):
        section=""
        for line in self.lines:
            secName=getSecName(line, level=2)
            if secName:
                section=secName
            if 'moltype' in section:
                (id, name)=line.split('=')[1].split()
                self.moltypeID = int(id)
                self.moltypeName = name.strip('"').strip('+').strip('-')
                if 'KRAS' in self.moltypeName:
                    self.moltypeName='RAS'
                if 'RAF' in self.moltypeName:
                    self.moltypeName='RAS_RAF'
            if '#molecules' in section:
                self.numMols = int(line.split('=')[1])
            """
            self.numPosResxA = int(line.split('=')[1])
            nrecord=4
            if self.numPosResxA>0:
                self.posResxA = self.lines[5:5+self.numPosResxA]
                nrecord = 5+self.numPosResxA
            self.numPosResxB = int(self.lines[nrecord].split('=')[1])
            if self.numPosResxB > 0:
                nrecord+=2
                self.posResxB = self.lines[nrecord:nrecord+self.numPosResxA]
            """
class Toplogy():
    def __init__(self, tprLog):
        logger.info("Create Toplogy object")
        #'molblock', 'ffparams'
        self.sections = []
        self.lines = tprLog.getTopology().split('\n')
        self.boxLines = tprLog.getBox().split('\n')
        self.box=[]
        self.name = ''
        self.numAtoms = 0
        self.numMolBlk = 0
        self.molBlks =[]
        self.atnr=0
        self.ntypes=0
        self.functype=[]

        self._parse()

    def _parsebox(self):
        for b in self.boxLines:
            if '{' in b:
                strs = re.split('\{|\}', b)
                nums = strs[1].split(',')
                self.box.append([float(a)*10 for a in nums])
        assert len(self.box) == 3, 'Box dimesion does not equal to 3'

    def _parse(self):
        self._parsebox()
        molblockLines=[]
        ffparamsLines=[]
        section = ""
        for line in self.lines:
            secName=getSecName(line, level=1)
            if secName:
                section=secName
            if 'name' in section:
                self.name = self.lines[1].split('=')[1].strip('"').strip('+').strip('-')
            if '#atom' in section:
                self.numAtoms = int(self.lines[2].split('=')[1])
            if '#molblock' in section:
                self.numMolBlk = int(self.lines[3].split('=')[1])
            if 'molblock' in section:
                molblockLines.append(line)
            if 'ffparams' in section:
                ffparamsLines.append(line)

        for line in ffparamsLines:
            section=""
            secName=getSecName(line, level=2)
            if secName:
                section=secName
            if 'atnr' in section:
                self.atnr = int(line.split('=')[1])
            if 'types' in section:
                self.ntypes = int(line.split('=')[1])
            if 'functype' in line:
                funcType={}
                strs=line.split(',')
                try:
                    fStrs=strs[0].split('=')
                except:
                    logger.error("functype format is not correct")

                try:
                    funcType['type']=fStrs[1]
                    iStrs=re.split('\[|\]',fStrs[0])
                    funcType['index'] = int(iStrs[1])
                except:
                    logger.error("functype[...] = ... format is not correct")

                for item in strs[1:]:
                    dStrs=item.split('=')
                    try:
                        funcType[dStrs[0].strip()] = float(dStrs[1])
                    except:
                        logger.error("functype[] = , ... = ..., format is not correct: "+line)

                self.functype.append(funcType)

        if len(self.functype)!=self.ntypes:
            logger.error("Number of functype doesn't equal to ntypes")

        #self.molBlks = [None] * self.numMolBlk
        molblockComLines=""
        for line in molblockLines[1:]:
            molblockComLines=molblockComLines+line+"\n"
        molblockStrs = molblockComLines.split("molblock")

        for mbstr in molblockStrs[1:]:
            molblck = MolBlock()
            molblck.lines=mbstr.split('\n')
            self.molBlks.append(molblck)

        for molblck in self.molBlks:
            molblck.parse()



