from BooleanInference.MIDNI.LoadNodes import LoadNodes
from sklearn.metrics import mutual_info_score
from scipy.stats import entropy
import numpy as np
import math


class MIDNIMutualInformationCalculation:

    def __init__(self, path):
        self.nodes = LoadNodes.addNodes(path)
        self.nodeSize = len(self.nodes)
        self.nodeLength = len(self.nodes.__getitem__(0))
        self.setF = []
        self.setS = []
        self.setAS = []
        self.tabularMI = [[]]
        self.tabularMI2 = [[]]

    def prepare(self):
        # F Lưu tất cả các node của network
        self.setF = []
        self.setAS = []
        # S lưu các gene đã được chọn
        self.setS = []
        self.initializeSetF()
        # self.initializeTabularMI()

    def prepare_for_each_target_gene(self, target):
        # F Lưu tất cả các node của network
        self.setF = []
        self.setAS = []
        # S lưu các gene đã được chọn
        self.setS = []
        self.initialize_SetF_No_Self_Regulator(target)

    def initializeSetF(self):
        for i in range(self.nodeSize):
            self.setF.append(i)

    def initialize_SetF_No_Self_Regulator(self, target):
        for i in range(self.nodeSize):
            if target != i:
                self.setF.append(i)

    """
    This function can be written by another simpler way using numpy library
    """

    def initializeTabularMI(self, print_out=False):
        self.tabularMI = np.zeros((self.nodeSize, self.nodeSize), dtype=float)
        self.tabularMI2 = np.zeros((self.nodeSize, self.nodeSize), dtype=float)
        for i in range(self.nodeSize):
            for j in range(self.nodeSize):
                if i != j:
                    self.tabularMI[i][j] = self.calculatePairwiseMutualInformation(i, j)
                    self.tabularMI2[i][j] = self.calculatePairwiseMutualInformationl(i, j)
                else:
                    self.tabularMI[i][j] = -1
                    self.tabularMI2[i][j] = -1
        if print_out:
            MIDNIMutualInformationCalculation.display_matrix(self.tabularMI)
            print("====== APPROXIMATE DEPENDENCY TABULAR ======")
            MIDNIMutualInformationCalculation.display_matrix(self.tabularMI2)

    @staticmethod
    def display_matrix(matrix):
        index = 1
        for i in matrix:
            if index == 1:
                for v in i:
                    print('\t\t\t\tG%d' % index, end=' ')
                    index = index + 1
                index = 1
                print('')
            print('G%d' % index, end='<--'),
            for v in i:
                print(v, end=' ')
            index = index + 1
            print('')

    def getResult(self):
        r = []
        for i in self.setS:
            r.append(i + 1)
        return r

    def executeFeatureSelection(self, target_gene, k):
        if len(self.setF) == 0:
            raise Exception("The method prepare must be called first!")
        if k == 0:
            return
        # Check entropy of a target gene, H(v0) at first
        entropyValue = self.entropyCheckingOfTargetGene(target_gene)
        if entropyValue == 0:
            return
        # Third steps of MIFS - Choice the first feature find the feature v that
        # maximizes I(v0,w); set W <- W\{v}; set S <- S Union {v}
        # featureMax là index của gene có MI lớn nhất
        featureMax = self.getMax_MI_features(target_gene)
        print("index before ", featureMax)
        self.setF.remove(featureMax)
        self.setS.append(featureMax)
        # Step 4 - select next feature and it is assume Greedy Selection: repeat until |S| = k
        while len(self.setS) < k:
            featureMax = self.getMax_MI_features1()
            print("index after ", featureMax)
            self.setF.remove(featureMax)
            self.setS.append(featureMax)
        for i in range(1, k):
            featureMax = self._getNextMax(target_gene)
            self.setAS.append(featureMax)
        print("set s sau 1 gene=======")
        print(self.setS)

    """
    target_gene - chỉ số của target gene, bắt đầu từ 0
    k - maximum incoming link
    """

    def my_executeFeatureSelection_original_implementation(self, target_gene, k, print_out=False):
        self.prepare_for_each_target_gene(target_gene)
        if len(self.setF) == 0:
            raise Exception("The method prepare must be called first!")
        if k == 0:
            return
        # Check entropy of a target gene, H(v0) at first
        entropyValue = self.entropyCheckingOfTargetGene(target_gene)
        if entropyValue == 0:
            return
        # Third steps of MIFS - Choice the first feature find the feature v that
        # maximizes I(v0,w); set W <- W\{v}; set S <- S Union {v}
        # featureMax là index của gene có MI lớn nhất
        featureMax = self.getMax_MI_features(target_gene)
        self.setF.remove(featureMax)
        self.setS.append(featureMax)
        if print_out:
            print('============================== TARGET GENE %d ==============================' % (target_gene + 1))
            print("Gene %d has maximum MI with gene %d" % (featureMax + 1, target_gene + 1))
            print("========== SET F ==========")
            MIDNIMutualInformationCalculation.print_set(self.setF)
            print("========== SET S ==========")
            MIDNIMutualInformationCalculation.print_set(self.setS)
        dynamics_consistency = self.dynamics_consistency(target_gene)
        if dynamics_consistency == 1:
            if print_out:
                print("============================= PARTIALLY RESULT =============================")
                print("========== SET S ==========")
                MIDNIMutualInformationCalculation.print_set(self.setS)
                dynamics_consistency = self.dynamics_consistency(target_gene)
                print("dynamics_consistency = ", dynamics_consistency)
            return
        # Step 4 - select next feature and it is assume Greedy Selection: repeat until |S| = k
        while len(self.setS) < k:
            featureMax = self.select_next_node(target_gene, print_out)
            self.setF.remove(featureMax)
            self.setS.append(featureMax)
            # Tính dynamic accuracy
            dynamics_consistency = self.dynamics_consistency(target_gene)
            if print_out:
                print("Next gene %d has maximum approximate dependency with the gene selected" % (featureMax + 1))
                print("========== SET F ==========")
                MIDNIMutualInformationCalculation.print_set(self.setF)
                print("========== SET S ==========")
                MIDNIMutualInformationCalculation.print_set(self.setS)
                print("dynamics_consistency = ", dynamics_consistency)
            if dynamics_consistency == 1:
                break
            if self.SWAP(target_gene):
                break
        if print_out:
            print("============================= PARTIALLY RESULT =============================")
            print("========== SET S ==========")
            MIDNIMutualInformationCalculation.print_set(self.setS)
            dynamics_consistency = self.dynamics_consistency(target_gene)
            print("dynamics_consistency = ", dynamics_consistency)

    def my_executeFeatureSelection_gene_pair_wise(self, target_gene, k):
        if len(self.setF) == 0:
            raise Exception("The method prepare must be called first!")
        if k == 0:
            return
        # Check entropy of a target gene, H(v0) at first
        entropyValue = self.entropyCheckingOfTargetGene(target_gene)
        if entropyValue == 0:
            return
        # Third steps of MIFS - Choice the first feature find the feature v that
        # maximizes I(v0,w); set W <- W\{v}; set S <- S Union {v}
        # featureMax là index của gene có MI lớn nhất
        featureMax = self.getMax_MI_features(target_gene)
        print("MI max = ", featureMax)
        self.setF.remove(featureMax)
        self.setS.append(featureMax)
        # Step 4 - select next feature and it is assume Greedy Selection: repeat until |S| = k
        # while len(self.setS) < k:
        #     featureMax = self.getMax_MI_featuresModified()
        #     print("index after ", featureMax)
        #     self.setF.remove(featureMax)
        #     self.setS.append(featureMax)
        # for i in range(1, k):
        #     featureMax = self._getNextMax(target_gene - 1)
        #     self.setAS.append(featureMax)
        current_gene = {i: self.tabularMI[target_gene][i] for i in range(0, len(self.tabularMI[target_gene]))}
        sorted_gene = {}
        sorted_index = sorted(current_gene, key=current_gene.get)
        for w in sorted_index:
            sorted_gene[w] = current_gene[w]
        print(sorted_index)
        self.setS = []
        for y in range(k):
            self.setS.append(sorted_index[len(sorted_index) - 1 - y])
        print("set s sau 1 gene=======")
        print(self.setS)

    """
    Hàm lựa chọn next candidate
    target - index của target node
    """

    def select_next_node(self, target, print_out=False):
        score_set = []
        if print_out:
            print("=========================== SCORE_FOR_NEXT_CHOICE ===========================")
        for i in self.setF:
            mutual = self.tabularMI[target][i]
            sigma = self.sigma_approximate_dependency_next_candidate_set_S(i)
            score_set.append(mutual - sigma)
            if print_out:
                print(mutual - sigma, end=' ')
        if print_out:
            print('')
        max_score = score_set[0]
        selected_index = 0
        for i in range(len(score_set)):
            if max_score < score_set[i]:
                max_score = score_set[i]
                selected_index = i
        return self.setF[selected_index]

    """
    Hàm tính tổng approximate_dependency của node kế tiếp với tập S
    next - index của node bất kì trong tập chưa chọn
    """

    def sigma_approximate_dependency_next_candidate_set_S(self, next_node):
        sigma = 0
        for i in self.setS:
            sigma = sigma + self.tabularMI2[i][next_node]
        return sigma

    @staticmethod
    def print_set(s):
        for i in s:
            print('G%s' % (i + 1), end=' ')
        print('')

    def getMax_MI_features(self, index):
        imax = 0
        valormax = self.tabularMI[index][0]
        for i in range(1, self.nodeSize):
            if valormax < self.tabularMI[index][i]:
                valormax = self.tabularMI[index][i]
                imax = i
        return imax

    def getMax_MI_features1(self):
        imax = -1
        valormax = -1
        for index in self.setS:
            for i in self.setF:
                if valormax < self.tabularMI2[index][i]:
                    valormax = self.tabularMI2[index][i]
                    imax = i
        return imax

    def getMax_MI_featuresModified(self):
        imax = -1
        valormax = -1
        for index in self.setS:
            for i in self.setF:
                if valormax < self.tabularMI[index][i]:
                    valormax = self.tabularMI[index][i]
                    imax = i
        return imax

    def _getNextMax(self, index):
        imax = -1
        valormax = -1
        for i in range(0, self.nodeSize):
            found = False
            for s in self.setS:
                if s == i:
                    found = True
                    continue
            if found:
                continue
            for s in self.setAS:
                if s == i:
                    found = True
                    continue
            if found:
                continue
            if valormax < self.tabularMI[index][i]:
                valormax = self.tabularMI[index][i]
                imax = i
        return imax

    def entropyCheckingOfTargetGene(self, x):
        s = self.nodes[x]
        # firstVector = []
        # for i in range(self.nodeLength - 1):
        #     firstVector.append(s[i + 1])
        return entropy(s[1:])

    @staticmethod
    def getMutualInformation(value):
        return mutual_info_score(value[0], value[1])

    @staticmethod
    def getMutualInformation2(x, y):
        return mutual_info_score(x, y)

    def calculatePairwiseMutualInformationl(self, x, y):
        # s = self.nodes.__getitem__(x)
        # q = self.nodes.__getitem__(y)
        # firstVector = []
        # secondVector = []
        # tmp = []
        # for i in range(self.nodeLength):
        #     firstVector.append(s[i])
        #     secondVector.append(q[i])
        # tmp.append(s)
        # tmp.append(q)
        return MIDNIMutualInformationCalculation.getMutualInformation2(self.nodes[x], self.nodes[y])

    # x là index của target gene
    # y là index của regulator gene
    def calculatePairwiseMutualInformation(self, x, y):
        s = self.nodes.__getitem__(x)
        q = self.nodes.__getitem__(y)
        return MIDNIMutualInformationCalculation.getMutualInformation2(s[1:], q[0:self.nodeLength - 1])

    """:cvar
    Hàm tính dynamics consistency tại thời điểm sau khi chọn 1 node xong
    target - index của target gene
    """

    def dynamics_consistency(self, target):
        observed = []
        input = []
        for i in range(self.nodeLength - 1):
            step = ""
            for j in self.setS:
                step = step + str(self.nodes[j][i])
            input.append(step)
            step = step + str(self.nodes[target][i + 1])
            observed.append(step)
        key_observed = list(dict.fromkeys(observed))
        key_input = list(dict.fromkeys(input))
        observed_dict = {i: 0 for i in key_observed}
        input_dict = {i: 0 for i in key_input}
        for key in key_observed:
            for obr in observed:
                if key == obr:
                    observed_dict[key] = observed_dict[key] + 1
        maximums = 0
        for i in key_input:
            max = 1
            for key in key_observed:
                if i == key[:len(key) - 1]:
                    if max < observed_dict[key]:
                        max = observed_dict[key]
            maximums = maximums + max
            input_dict[i] = max
        return maximums / (self.nodeLength - 1)

    def SWAP_ONLY_LAST_ITEM(self, target):
        perfect_meet = False
        changing = self.setS[len(self.setS) - 1]
        for f in self.setF:
            self.setS[len(self.setS) - 1] = f
            if self.dynamics_consistency(target) == 1:
                perfect_meet = True
                break
        if not perfect_meet:
            self.setS[len(self.setS) - 1] = changing
        return perfect_meet

    def SWAP(self, target):
        perfect_meet = False
        for s in range(len(self.setS)):
            changing = self.setS[s]
            for f in range(len(self.setF)):
                self.setS[s] = self.setF[f]
                self.setF[f] = changing
                if self.dynamics_consistency(target) == 1:
                    perfect_meet = True
                    return perfect_meet
                else:
                    self.setF[f] = self.setS[s]
                    self.setS[s] = changing
        return perfect_meet

    def printSetF(self):
        print(self.setF)

    def printSetS(self):
        print(self.setS)

    def printSetAS(self):
        print(self.setAS)
