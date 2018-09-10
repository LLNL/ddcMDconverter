__author__ = 'zhang30'

class Velocity:
    def __init__(self, vx, vy, vz):
        self.vx = vx
        self.vy = vy
        self.vz = vz


class GRO:
    def __init__(self):
        self.vList=[]

    def parse(self, filename):
        scale=0.01
        totAtmNum=0
        count=0
        with open(filename, "r") as f:
            for line in f:
                if count==1:
                    totAtmNum=int(line.rstrip())
                if count>1 and count < totAtmNum+2:
                    #print(line[44:52], line[52:60], line[60:68])
                    vx = scale * float(line[44:52])
                    vy = scale * float(line[52:60])
                    vz = scale * float(line[60:68])
                    vel= Velocity(vx, vy, vz)
                    self.vList.append(vel)

                count=count+1

        if len(self.vList) != totAtmNum:
            raise Exception("GRO::parse - Number of total atoms doesn't match header record")

    def zeroVelocity(self, input, output):

        outFh=open(output, "w")
        totAtmNum = 0
        count = 0

        with open(input, "r") as f:
            for line in f:
                outLine = line
                if count == 1:
                    totAtmNum = int(line.rstrip())
                if count > 1 and count < totAtmNum + 2:
                    # print(line[44:52], line[52:60], line[60:68])
                    outLine=line[0:44]+"  0.0000  0.0000  0.0000\n"

                count = count + 1
                outFh.write(outLine)

        if (count-3) != totAtmNum:
            print("Number of atoms doesn't match the record")