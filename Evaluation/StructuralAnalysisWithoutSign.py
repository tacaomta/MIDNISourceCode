import re


class StructuralAnalysis:
    def __init__(self, gold_standard, gold_standard_signed, prediction):
        self.gold_standard = gold_standard
        self.prediction = prediction
        self.gold_standard_signed = gold_standard_signed

    def analysis(self):
        if self.gold_standard != "":
            gold_standard_set = StructuralAnalysis.readfile(self.gold_standard)
            true_negative_set = []
            for i in gold_standard_set:
                if i[2] == "0":
                    true_negative_set.append(i)
        else:
            true_negative_set = []
        gold_standard_signed_set = StructuralAnalysis.readfile(self.gold_standard_signed)
        prediction_set = StructuralAnalysis.readfile(self.prediction)
        prediction_standard_set = []
        for i in prediction_set:
            if float(i[2]) > 0:
                i[2] = "+"
            else:
                i[2] = "-"
            prediction_standard_set.append(i)
        true_positive = 0
        true_negative = len(true_negative_set)
        false_negative = 0
        for i in prediction_standard_set:
            for j in gold_standard_signed_set:
                if i[0] == j[0] and i[1] == j[1] and i[2] == j[2]:
                    true_positive = true_positive + 1
            for n in true_negative_set:
                if i[0] == n[0] and i[1] == n[1]:
                    true_negative = true_negative - 1
                    false_negative = false_negative + 1
        precision = true_positive / len(prediction_standard_set)
        recall = true_positive / len(gold_standard_signed_set)
        if recall != 0:
            PR = precision / recall
        else:
            PR = 'NAN'
        false_positive = len(prediction_standard_set) - true_positive
        ROC = true_positive / false_positive
        print('Total prediction = {}'.format(len(prediction_standard_set)))
        print('Total actual relations = {}'.format(len(gold_standard_signed_set)))
        print('True positive = {}'.format(true_positive))
        print('False positive = {}'.format(len(prediction_standard_set) - true_positive))
        print('True negative = {}'.format(true_negative))
        print('False negative = {}'.format(false_negative))
        print("Precision = true positive/ (true_positive + false_positive) = {}/({}+{}) = {}".format(true_positive,
                                                                                                     true_positive,
                                                                                                     false_positive,
                                                                                                     precision))
        print("Recall = true positive/ (true_positive + false_negative) = {}/({}+{}) = {}".format(true_positive,
                                                                                                  true_positive,
                                                                                                  false_negative,
                                                                                                  recall))
        print("Structural Accuracy = {}%".format(
            (true_positive + true_negative) * 100 / (len(prediction_standard_set) + len(
                gold_standard_signed_set))))
        print("PR = ", PR)
        print("ROC = true positive/ false positive = {}/{} = {}".format(true_positive,
                                                                        len(prediction_standard_set) - true_positive,
                                                                        ROC))

    @staticmethod
    def readfile(path):
        f = open(path, "r")
        txt = []
        gs = []
        for x in f:
            arr = re.split('\t|\n', x)
            txt.append(arr)
        f.close()
        for i in range(len(txt)):
            rowsitem = []
            for j in range(len(txt[0]) - 1):
                a = txt[i][j]
                rowsitem.append(a)
            if len(rowsitem) > 0:
                gs.append(rowsitem)
        return gs


# struct = StructuralAnalysis("Ecoli_size10_noise_0%-1_goldstanddard.tsv",
#                             "Ecoli_size10_noise_0%-1_goldstandard_signed.tsv", "Ecoli_size10_01_result.txt")
# struct.analysis()
struct = StructuralAnalysis('',
                            "A_goldstandard_signed.tsv", "A_result_network.txt")
struct.analysis()
