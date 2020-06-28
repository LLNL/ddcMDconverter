from enum import Enum

class ObjectFormat(Enum):
    INT=1
    FLOAT=2
    STRING=3
    LIST=4

class ObjectItem:
    def __init__(self, name, type, attrs):
        self._name=name
        self._type=type
        self._attrs=self._parseAttrs(attrs)

    def _parseAttrs(self, attrs):
        attrDict={}
        attrList=attrs.split(";")
        for attr in attrList[:-1]:
            itemList=attr.split("=")
            if len(itemList) !=2:
                raise Exception("Obejct item is not paired with = :" +attr)
            key=itemList[0].strip()
            value=itemList[1].strip()
            attrDict[key]=value
        return attrDict

    def getAttrFormt(self, attrName, format):
        if not isinstance(format, ObjectFormat):
            raise TypeError('format must be an instance of ObjectFormat Enum')

    def getAttr(self, attrName):
        if attrName not in self._attrs:
            raise Exception("Object attribute name " + attrName + " is not in the list")
        return self._attrs[attrName]

    def getAttrNames(self):
        return self._attrs.keys()

class ObjectList:
    def __init__(self, fileName):
        self._objectDict={}
        self.readFile(fileName)

    def readFile(self, fileName):
        contents=''
        with open(fileName,'r') as f:
            for line in f:
                contents=contents+line

        contentList=contents.split("}")
        for item in contentList[:-1]: # skip last one which is empty
            itemList=item.split('{')
            if len(itemList) != 2:
                raise Exception("Object curly brackets do not match: "+item+"}")

            itemList = item.replace('\n', ' ').split('{')
            nameList=itemList[0].split()
            if len(nameList) !=2:
                raise Exception("Object missing name or type")
            name=nameList[0].strip()
            type=nameList[1].strip()
            attrs=itemList[1].strip()
            obj=ObjectItem(name, type, attrs)
            self._objectDict[name]=obj

    def getObject(self, name):
        if name not in self._objectDict:
            raise Exception("Object name " + name + " is not in the list")
        return self._objectDict[name]





