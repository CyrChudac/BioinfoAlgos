from os import path

class Task():
    def __init__(self):
        self.name = "undefined"
    def help(self, pre: str, mult : int):
        print(pre*mult + "help not defined in class " + self.name)
    def help(self):
        self.help("\t", 0)
    def run(self, params):
        print("run not defined for class " + self.name)

class Result():
    def __init__(self, value, displayF = lambda x: Result.printVal(x)):
        self.val = value
        self.show = displayF
    def display(self):
        self.show(self.val)

        
    def printList(list):
        for x in list:
            print(x)
    class ListPrinter():
        def __init__(self, func):
            self.func = func
        def show(self, list):
            for x in list:
                print(self.func(x))
    def printList_advanced(func):
        return Result.ListPrinter(func).show
    def printVal(v):
        print(v)

def isFile(x):
    return path.isfile(x)