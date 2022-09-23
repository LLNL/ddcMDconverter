import logging
from ddcmdconverter.gmx.gmxTopology import Toplogy
from ddcmdconverter.gmx.gmxParameter import Parameter
from gmxParameter import getSecName
from itertools import islice
import re

# set up logger
LOGLEVEL = 1
LOG_FMT = '%(asctime)s - %(name)s:%(funcName)s:%(lineno)s - %(levelname)s - %(message)s'
logger = logging.getLogger(__name__)
logger.setLevel(LOGLEVEL)

sh = logging.StreamHandler()
sh.setLevel(LOGLEVEL)
sh.setFormatter(logging.Formatter(LOG_FMT))
logger.addHandler(sh)

class TPRLog():
    def __init__(self, tprlog):
        logger.info("Create TPRLog object")
        #'header', 'topology', 'cmap','grp'
        #self.sections=['head', 'topo', 'cmap','grp[']
        #self.secflag=[False, False, False]
        self.header=''
        self.topology=''
        self.cmap=''
        self.boxlines=''
        self.coors=''
        self.velocity=''
        self.tprlog=tprlog
        self._parse()


    def _parse(self):
        logger.info("Start to parse tpr.log file ... ")
        with open(self.tprlog, 'r') as f:
            #secIdx=0
            section = ""
            for line in f:
                secName = getSecName(line, level=0)
                if secName:
                    section = secName
                if section == 'header':
                    self.header = self.header + line
                if section == 'topology':
                    self.topology = self.topology+line
                if section[0:4] == 'cmap':
                    self.cmap = self.cmap + line
                if section == 'box (3x3)':
                    self.boxlines = self.boxlines + line
                if section[0:3] == 'x (':
                    logger.info("Parse 3D coordinates ... ")
                    numberStr = re.split('\(|\)', section)[1]
                    numAtoms = int(re.split('x', numberStr)[0])
                    lines_gen = islice(f, numAtoms)
                    self.coors = ''.join(x for x in lines_gen)

                if section[0:3] == 'v (':
                    logger.info("Parse velocitys ... ")
                    numberStr = re.split('\(|\)', section)[1]
                    numAtoms = int(re.split('x', numberStr)[0])
                    lines_gen = islice(f, numAtoms)
                    self.velocity = ''.join(x for x in lines_gen)


        logger.info("End of parse tpr.log file")

    def getHeader(self):
        return  self.header

    def getTopology(self):
        return self.topology

    def getCmap(self):
        return self.cmap

    def getBox(self):
        return self.boxlines

    def getCoors(self):
        return self.coors

    def getVelocity(self):
        return self.velocity

def main():
    tprLog=TPRLog('tpr.log')
    topology=Toplogy(tprLog)
    parameter=Parameter(tprLog)
    pass

if __name__ == '__main__':
    main()




