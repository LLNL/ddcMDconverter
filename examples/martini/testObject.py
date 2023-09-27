import sys
sys.path.append('/Users/zhang30/workspace/PycharmProjects/ddcMDconverter')

from ddcmdconverter.base import Object


def main():
    objectList=Object.ObjectList("/Users/zhang30/workspace/PycharmProjects/ddcMDconverter/test/martini/test/martini.data")

    martiniObject=objectList.getObject('martini')
    martiniAttrNames=martiniObject.getAttrNames()
    print(martiniAttrNames)
    resiParms=martiniObject.getAttr('resiParms')
    print(resiParms)
    resList=resiParms.split()
    for res in resList:
        resObject=objectList.getObject(res)
        numAtoms = resObject.getAttr('numAtoms')
        print(res, numAtoms)
        break


if __name__ == '__main__':
    main()