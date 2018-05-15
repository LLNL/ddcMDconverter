__author__ = 'zhang30'


import argparse
from ddcmdconverter.charmm.Psf import Psf



def getArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--psf', action='store', dest='psffile',  help='CHARMM psf file.')
    args = parser.parse_args()

    return args


def main():
    args=getArgs()
    print("Default inputs: ", args.psffile)
    psf=Psf.Psf()
    psf.parse(args.psffile)
    psf.printAngle()
    psf.printDihe()


if __name__ == '__main__':
    main()
