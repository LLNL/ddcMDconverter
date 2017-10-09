__author__ = 'zhang30'


import argparse
import Psf



def getArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--psf', action='store', dest='psffile',  help='CHARMM psf file.')
    args = parser.parse_args()

    return args

if __name__ == '__main__':

    args=getArgs()
    print "Default inputs: ", args.psffile
    psf=Psf.Psf()
    psf.parse(args.psffile)
    psf.printAngle()
