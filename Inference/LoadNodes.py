import re

class LoadNodes:

    def __init__(self, path):
        self.nodes = LoadNodes.addNodes(path)
        self.nodeSize = self.getNodeSize
        self.nodeLength = self.getNodeLength()

    """
    This function returns the number of nodes of the network
    In the case of "NetworkTransition.txt" file, this value equals 10
    """

    def getNodeSize(self):
        if self.nodes is None:
            return -1
        else:
            return len(self.nodes)

    """
    This function returns the length of each node of the network
    In the case of "NetworkTransition.txt" file, this value equals 20 
    """

    def getNodeLength(self):
        if self.nodes is None:
            return -1
        else:
            return len(self.nodes[0])

    """
    This function returns the nodes itself
    """

    def getNodes(self):
        return self.nodes

    """
    This function reads the list of nodes which is saved in the file
    """

    @staticmethod
    def addNodes(path):
        f = open(path, "r")
        txt = []
        nodes = []
        for x in f:
            arr = re.split(' |\t', x)
            txt.append(arr)
        f.close()
        for i in range(len(txt[0])):
            rowsitem = []
            for j in range(len(txt)):
                try:
                    rowsitem.append(int(txt.__getitem__(j).__getitem__(i)))
                except:
                    continue
            if len(rowsitem) > 0:
                nodes.append(rowsitem)
        return nodes

    """
    This function just prints the nodes for checking correctness
    """

    def printNode(self):
        for node in self.nodes:
            print(node)


# n = LoadNodes("discretized_genes.txt")
# n.printNode()
