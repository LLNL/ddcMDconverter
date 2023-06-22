import numpy as np
import logging
import re

class MartiniData:
    def __init__(self, mfile):
        self.martinifile = mfile
        self.keys = ['mass', 'atom', 'bond', 'cons', 'lj', 'exclusion']
        self.data = {}
        for k in self.keys:
            self.data[k]=''
        self.objs = {}
        for k in self.keys:
            self.objs[k]=[]

        self.parse()

    def parse(self):
        with open(self.martinifile) as f:
            for line in f:
                if 'MASSPARMS' in line:
                    self.data['mass'] = self.data['mass'] + line
                elif 'ATOMPARMS' in line:
                    self.data['atom'] = self.data['atom'] + line
                elif 'BONDPARMS' in line:
                    self.data['bond'] = self.data['bond'] + line
                elif 'CONSPARMS' in line:
                    self.data['cons'] = self.data['cons'] + line
                elif 'LJPARMS' in line:
                    self.data['lj'] = self.data['lj'] + line
                elif 'EXCLUDEPARMS' in line:
                    self.data['exclusion'] = self.data['exclusion'] + line

        for k in self.keys:
            lines = self.data[k].split('\n')
            for line in lines:
                success, object = self._parseLine(line)
                if success:
                    self.objs[k].append(object)

    def isfloat(self, value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    def _parseLine(self, line):
        obj={}
        #strs = line.split('{|}')
        strs = re.split('\{|\}', line)
        if len(strs) !=3 :
            logging.warning("Object line format error " + line)
            return False, obj
        name = strs[0].split()[0]
        obj['name'] = name
        kvs = strs[1].split(';')
        for kv in kvs:
            if '=' in kv:
                kvstrs=kv.split('=')
                key = kvstrs[0].strip()
                value = kvstrs[1].split()[0]
                if value.isdigit():
                    value = int(value)
                elif self.isfloat(value):
                    value = float(value)
                obj[key] = value
        return True, obj

    def save_npz(self, file):
        massAtomType=[]
        massAtomTypeID=[]
        mass =[]
        for val in self.objs['mass']:
            massAtomType.append(val['atomType'])
            massAtomTypeID.append(val['atomTypeID'])
            mass.append(val['mass'])
        massAtomType_np = np.array(massAtomType)
        massAtomTypeID_np = np.array(massAtomTypeID)
        mass_np = np.array(mass)
        types = np.rec.fromarrays((massAtomType_np, massAtomTypeID_np, mass_np), names=('atomType', 'atomTypeID', 'mass'))

        atomName=[]
        atomID = []
        atomName2=[]
        atomType=[]
        atomTypeID=[]
        atomCharge=[]
        for val in self.objs['atom']:
            atomName.append(val['name'])
            atomID.append(val['atomID'])
            atomName2.append(val['atomName'])
            atomType.append(val['atomType'])
            atomTypeID.append(val['atomTypeID'])
            atomCharge.append(val['charge'])

        atomName_np = np.array(atomName)
        atomID_np = np.array(atomID)
        atomName2_np = np.array(atomName2)
        atomType_np = np.array(atomType)
        atomTypeID_np = np.array(atomTypeID)
        atomCharge_np = np.array(atomCharge)
        atoms = np.rec.fromarrays((atomName_np, atomID_np, atomName2_np, atomType_np, atomTypeID_np, atomCharge_np),
                                  names=('name', 'atomID', 'atomName', 'atomType', 'atomTypeID', 'charge'))

        bondName=[]
        bondAtomI=[]
        bondAtomTypeI=[]
        bondAtomJ=[]
        bondAtomTypeJ=[]
        bondFuc=[]
        bondkb=[]
        bondb0=[]
        for val in self.objs['bond']:
            bondName.append(val['name'])
            bondAtomI.append(val['atomI'])
            bondAtomTypeI.append(val['atomTypeI'])
            bondAtomJ.append(val['atomJ'])
            bondAtomTypeJ.append(val['atomTypeJ'])
            bondFuc.append(val['func'])
            bondkb.append(val['kb'])
            bondb0.append(val['b0'])

        bondName_np = np.array(bondName)
        bondAtomI_np = np.array(bondAtomI)
        bondAtomTypeI_np = np.array(bondAtomTypeI)
        bondAtomJ_np = np.array(bondAtomJ)
        bondAtomTypeJ_np = np.array(bondAtomTypeJ)
        bondFuc_np = np.array(bondFuc)
        bondkb_np = np.array(bondkb)
        bondb0_np = np.array(bondb0)
        bonds = np.rec.fromarrays((bondName_np, bondAtomI_np, bondAtomTypeI_np, bondAtomJ_np, bondAtomTypeJ_np, bondFuc_np, bondkb_np, bondb0_np),
                                  names=('name', 'atomI', 'atomTypeI', 'atomJ', 'atomTypeJ', 'func', 'kb', 'b0'))

        consName=[]
        consAtomI=[]
        consAtomTypeI=[]
        consAtomJ=[]
        consAtomTypeJ=[]
        consFuc=[]
        consr0=[]
        for val in self.objs['cons']:
            consName.append(val['name'])
            consAtomI.append(val['atomI'])
            consAtomTypeI.append(val['atomTypeI'])
            consAtomJ.append(val['atomJ'])
            consAtomTypeJ.append(val['atomTypeJ'])
            consFuc.append(val['func'])
            consr0.append(val['r0'])

        consName_np = np.array(consName)
        consAtomI_np = np.array(consAtomI)
        consAtomTypeI_np = np.array(consAtomTypeI)
        consAtomJ_np = np.array(consAtomJ)
        consAtomTypeJ_np = np.array(consAtomTypeJ)
        consFuc_np = np.array(consFuc)
        consr0_np = np.array(consr0)
        cons = np.rec.fromarrays((consName_np, consAtomI_np, consAtomTypeI_np, consAtomJ_np, consAtomTypeJ_np, consFuc_np, consr0_np),
                                  names=('name', 'atomI', 'atomTypeI', 'atomJ', 'atomTypeJ', 'func', 'r0'))

        exclusionName=[]
        exclusionAtomI=[]
        exclusionAtomTypeI=[]
        exclusionAtomJ=[]
        exclusionAtomTypeJ=[]
        for val in self.objs['exclusion']:
            exclusionName.append(val['name'])
            exclusionAtomI.append(val['atomI'])
            exclusionAtomTypeI.append(val['atomTypeI'])
            exclusionAtomJ.append(val['atomJ'])
            exclusionAtomTypeJ.append(val['atomTypeJ'])

        exclusionName_np = np.array(exclusionName)
        exclusionAtomI_np = np.array(exclusionAtomI)
        exclusionAtomTypeI_np = np.array(exclusionAtomTypeI)
        exclusionAtomJ_np = np.array(exclusionAtomJ)
        exclusionAtomTypeJ_np = np.array(exclusionAtomTypeJ)
        exclusion = np.rec.fromarrays((exclusionName_np, exclusionAtomI_np, exclusionAtomTypeI_np, exclusionAtomJ_np, exclusionAtomTypeJ_np),
                                  names=('name', 'atomI', 'atomTypeI', 'atomJ', 'atomTypeJ'))

        ljName=[]
        ljAtomTypeI = []
        ljIndexI=[]
        ljAtomTypeJ = []
        ljIndexJ = []
        ljSigma = []
        ljEps =[]
        for val in self.objs['lj']:
            ljName.append(val['name'])
            ljAtomTypeI.append(val['atomtypeI'])
            ljIndexI.append(val['indexI'])
            ljAtomTypeJ.append(val['atomtypeJ'])
            ljIndexJ.append(val['indexJ'])
            ljSigma.append(val['sigma'])
            ljEps.append(val['eps'])

        ljName_np = np.array(ljName)
        ljAtomTypeI_np = np.array(ljAtomTypeI)
        ljIndexI_np = np.array(ljIndexI)
        ljAtomTypeJ_np = np.array(ljAtomTypeJ)
        ljIndexJ_np = np.array(ljIndexJ)
        ljSigma_np = np.array(ljSigma)
        ljEps_np = np.array(ljEps)
        ljs = np.rec.fromarrays(
            (ljName_np, ljAtomTypeI_np, ljIndexI_np, ljAtomTypeJ_np, ljIndexJ_np, ljSigma_np, ljEps_np),
            names=('name', 'atomtypeI', 'indexI', 'atomtypeJ', 'indexJ', 'sigma', 'eps'))

        np.savez(file, types=types, atoms=atoms, bonds=bonds, cons=cons, exclusion=exclusion, lj=ljs)


def main():
    martiniData = MartiniData('martini.data')
    martiniData.save_npz('martini.npz')
    pass


if __name__ == '__main__':
    main()
