import Task

class Help(Task.Task):
    def __init__(self, dic, longer):
        self.dic = dic
        self.longer = longer
        return super(Help, self).__init__()
    def help(self, pre: str, mult: int):
        print(pre*mult + "show this help page")
    def run(self, params):
        print("TOOL HELP PAGE")
        pre = "\t"
        mult = 1
        print(pre*mult + "-h --help --crocodile --anything-unused")
        self.help(pre, mult+1)
        for x in self.longer:
            print(pre*mult + self.longer[x] + " " + x)
            self.dic[self.longer[x]].help(pre, mult+1)