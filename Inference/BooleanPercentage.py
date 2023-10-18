from glob import glob
from path import Path

from BooleanInference.MIBNILibrary.LoadNodes import LoadNodes


class RatioAnalysis:
    def __init__(self, path):
        self.nodes = LoadNodes.addNodes(path)
        self.nodeSize = len(self.nodes)

    def Percentage(self):
        high = 0
        for node in self.nodes:
            for i in node:
                if i == 2:
                    high = high + 1
                    break
        return high / self.nodeSize


folder = "C:/caocao/gnw-master/mydata/size200-noise-5percent/hybrid/values/*.txt"
x = [Path(f).abspath() for f in glob(folder)]
for f in x:
    net = RatioAnalysis(f)
    print(net.Percentage())

