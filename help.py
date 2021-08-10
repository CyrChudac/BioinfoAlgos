import Task

class Help(Task.Task):
    def __init__(self, dic, longer):
        self.dic = dic
        self.longer = longer
        return super(Help, self).__init__()
    def help(self, pre: str, mult: int):
        print(pre*mult + "show this page")
    def run(self, params):
        print("TOOL HELP PAGE")
        for x in self.longer:
            print("\t" + self.longer[x] + " " + x)
            self.dic[self.longer[x]].help("\t",2)