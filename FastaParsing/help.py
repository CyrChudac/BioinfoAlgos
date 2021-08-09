import Task

class Help(Task.Task):
    def __init__(self, dic):
        self.dic = dic
        return super(Help, self).__init__()
    def help(self, pre: str, mult: int):
        print(pre*mult + "show this page")
    def run(self, params):
        print("TOOL HELP PAGE")
        for x in self.dic:
            print("\t" + x)
            self.dic[x].help("\t",2)