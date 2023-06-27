import logging
import re
from itertools import islice
from ddcmdconverter.gmx.gmxParameter import getSecName

LOGLEVEL =2
LOG_FMT = '%(asctime)s - %(name)s:%(funcName)s:%(lineno)s - %(levelname)s - %(message)s'
logger = logging.getLogger(__name__)
logger.setLevel(LOGLEVEL)

sh = logging.StreamHandler()
sh.setLevel(LOGLEVEL)
sh.setFormatter(logging.Formatter(LOG_FMT))
logger.addHandler(sh)

class Coordinates():
    def __init__(self, box, numAtoms):
        logger.info("Create Toplogy object")
        #'molblock', 'ffparams'
        self.box = box
        self.numAtoms = numAtoms
        self.coors=[]
        self.velocity=[]

class TprCoors(Coordinates):
    def __init__(self, tprLog, topology):
        super().__init__(topology.box, topology.numAtoms)
        self.coorlines = tprLog.getCoors().split('\n')
        self.velocitylines = tprLog.getVelocity().split('\n')
        self.parse()

    def parseCoor(self, line):
        ''' x[    0]={ 2.57160e+01,  1.90230e+01,  1.48260e+01}'''
        strs = line.split('=')
        idx = int(re.split('\[|\]', strs[0])[1])
        tmpStr = re.split('\{|\}', strs[1])[1]
        numStrs = tmpStr.split(',')

        array = [float(x)*10 for x in numStrs]

        return (idx, array[0], array[1], array[2])

    def parse(self):
        for line in self.coorlines:
            if line=='':
                continue
            coor=self.parseCoor(line)
            self.coors.append(coor)

        for line in self.velocitylines:
            if line=='':
                continue
            vel = self.parseCoor(line)
            self.velocity.append(vel)

class GroCoors(Coordinates):
    def __init__(self, groFileName):
        numAtoms = 0
        box =[]
        coors=[]
        velocity=[]
        count=-1
        with open(groFileName, 'r') as f:
            line=f.readline()
            line = f.readline()
            numAtoms = int(line)
            lines_gen = islice(f, numAtoms)
            hasVelocity = len(line) > 67
            for line in lines_gen:
                #idx = int(line[15:20]) - 1
                count=count+1
                idx = count
                x = 10 * float(line[20:28])
                y = 10 * float(line[28:36])
                z = 10 * float(line[36:44])

                coors.append((idx, x, y, z))
                if hasVelocity:
                    vx = 10 * float(line[44:52])
                    vy = 10 * float(line[52:60])
                    vz = 10 * float(line[60:68])
                    velocity.append((idx, vx, vy, vz))
                else:
                    velocity.append((idx, 0.0, 0.0, 0.0))

            line = f.readline()
            strs = line.split()
            bx = [float(x) * 10 for x in strs]
            box.append([bx[0], 0.0, 0.0])
            box.append([0.0, bx[1], 0.0])
            box.append([0.0, 0.0, bx[2]])

        super().__init__(box, numAtoms)
        self.coors=coors
        self.velocity=velocity
