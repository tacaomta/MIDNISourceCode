from BooleanInference.MIBNISimple.LoadNodes import LoadNodes
import numpy as np


class MIBNIUpdateRules:

    def __init__(self, path):
        self.allNodes = LoadNodes(path)

    """
    This function returns inferred dynamic boolean matrix using Update Rules
    """

    def inferBooleanDynamics(self, targets, regulators):
        nodes = self.allNodes.getNodes()
        target = np.zeros(self.allNodes.getNodeLength() - 1, dtype=int)  # 1-D matrix
        sols = []  # 2-D matrix
        ret = []
        j = 1
        k = 0
        for i in range(len(targets)):
            for j in range(self.allNodes.getNodeLength()):
                target[j - 1] = nodes.__getitem__(targets[i] - 1)[j]
            sols = np.zeros((len(regulators[i]), self.allNodes.getNodeLength() - 1), dtype=int)
            for v in range(self.allNodes.getNodeLength() - 1):
                for k in range(len(regulators[i])):
                    sols[k][v] = nodes.__getitem__(regulators[i][k] - 1)[v]
            ret.insert(i, MIBNIUpdateRules.test(target, sols))
        return ret

    def calculateErrors(self, result):
        nodes = self.allNodes.getNodes()
        ret = np.zeros(len(nodes), dtype=int)
        j = 0
        for i in range(len(nodes)):
            originalNode = nodes.__getitem__(i)
            reconstructedNode = result[i]
            if reconstructedNode is None:
                continue
            for j in range(len(reconstructedNode)):
                if reconstructedNode[j] != originalNode[j + 1]:
                    ret[i] += 1
        return ret

    """
    This function returns a matrix of binary values based a given value - dim
    if dim = 3 then generated matrix is
    0 0 0
    0 0 1
    0 1 0
    0 1 1
    1 0 0
    1 0 1
    1 1 0
    1 1 1
    """

    @staticmethod
    def binaryPermutation(dim):
        perms = []
        permsl = []
        if dim == 0:
            return perms
        perms.append("0")
        perms.append("1")
        i = 1
        while i < dim:
            for s in perms:
                permsl.append(s + "0")
                permsl.append(s + "1")
            perms.clear()
            for x in permsl:
                perms.append(x)
            permsl.clear()
            i += 1
        return perms

    """
    This method for doing something
    where are:
     - target is 1-D matrix
     - solution is 2-D matrix
    """

    @staticmethod
    def test(target, solution):
        perms = MIBNIUpdateRules.binaryPermutation(len(solution))
        ands = np.zeros((len(perms), len(target)), dtype=int)
        ors = np.zeros((len(perms), len(target)), dtype=int)
        cAnds = np.zeros(len(perms), dtype=int)
        cOrs = np.zeros(len(perms), dtype=int)
        i = 0
        j = 0
        for i in range(len(target)):
            for j in range(len(perms)):
                str = perms.__getitem__(j)
                for k in range(len(str)):
                    if k == 0:
                        if str.__getitem__(k) == "0":
                            ands[j][i] = solution[k][i]
                            ors[j][i] = solution[k][i]
                        else:
                            ands[j][i] = solution[k][i] ^ 1
                            ors[j][i] = solution[k][i] ^ 1
                    else:
                        if str.__getitem__(k) == "0":
                            ands[j][i] &= solution[k][i]
                            ors[j][i] |= solution[k][i]
                        else:
                            ands[j][i] &= solution[k][i] ^ 1
                            ors[j][i] |= solution[k][i] ^ 1
                if target[i] == ands[j][i]:
                    cAnds[j] += 1
                if target[i] == ors[j][i]:
                    cOrs[j] += 1
        ret = []
        max = 0
        for i in range(len(cAnds)):
            if cAnds[i] > max:
                max = cAnds[i]
                ret = ands[i]
        for i in range(len(cOrs)):
            if cOrs[i] > max:
                max = cOrs[i]
                ret = ors[i]
        return ret
