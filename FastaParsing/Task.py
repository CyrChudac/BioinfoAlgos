class Task():
    def __init__(self):
        self.name = "undefined"
    def help(self):
        print("help not defined in class " + self.name)
    def run(self, ab):
        print("run not defined for class " + self.name)

