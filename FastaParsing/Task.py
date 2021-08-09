class Task():
    def __init__(self):
        self.name = "undefined"
    def help(self, pre: str, mult : int):
        print(pre*mult + "help not defined in class " + self.name)
    def help(self):
        self.help("\t", 0)
    def run(self, params):
        print("run not defined for class " + self.name)

