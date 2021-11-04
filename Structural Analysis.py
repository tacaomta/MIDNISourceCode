import re

"Cái này dành cho Dream4, khi tên trong structure theo kiểu G1, G2..."
"và với dữ liệu tự sinh"


class StructuralAnalysis:
    def __init__(self, prediction, gold_standard, gold_standard_signed=""):
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

    def none_sign_analysis(self, detail=False):
        gold_standard_set = StructuralAnalysis.readfile(self.gold_standard)
        negative_gold_standard_set = []
        positive_gold_standard_set = []
        for i in gold_standard_set:
            if i[2] == "0":
                negative_gold_standard_set.append(i)
            else:
                positive_gold_standard_set.append(i)
        prediction_set = StructuralAnalysis.readfile(self.prediction)
        positive_prediction_set = []
        negative_prediction_set = []
        for i in prediction_set:
            if i != "":
                if i[2] == "0":
                    negative_prediction_set.append(i)
                else:
                    positive_prediction_set.append(i)
        true_positive = 0
        false_positive = 0
        true_negative = 0
        false_negative = 0
        # Calculate true positive
        true_positive_set = []
        for i in positive_prediction_set:
            for j in positive_gold_standard_set:
                if i[0] == j[0] and i[1] == j[1]:
                    true_positive = true_positive + 1
                    true_positive_set.append(i)
        for i in negative_prediction_set:
            for j in negative_gold_standard_set:
                if i[0] == j[0] and i[1] == j[1]:
                    true_negative = true_negative + 1
        false_positive = len(positive_prediction_set) - true_positive
        false_negative = len(negative_prediction_set) - true_negative
        precision = true_positive / (true_positive + false_positive)
        recall = true_positive / (true_positive + false_negative)
        print('========== PREDICTION INFO ==========')
        print('Positive predicted relations = {}'.format(len(positive_prediction_set)))
        print('Negative predicted relations = {}'.format(len(negative_prediction_set)))
        print('========== ACTUAL INFO ==========')
        print('Positive actual relations = {}'.format(len(positive_gold_standard_set)))
        print('Negative actual relations = {}'.format(len(negative_gold_standard_set)))
        print('========== REPORT ==========')
        print('True positive = {}'.format(true_positive))
        if detail:
            print("Detail: ", end='')
            for i in true_positive_set:
                print(i, end='')
            print('\n')
        print('False positive = {}'.format(false_positive))
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
            (true_positive + true_negative) * 100 / (true_positive + false_positive + false_negative + true_negative)))

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
# for index in range(11, 21, 1):
#     print("======================={}=======================".format(index))
#     struct = StructuralAnalysis(
#         "C:\caocao\Python Project\ABC\generated_data\size200\hybrid\\values\dbn_output\size200_{}_dbn_prediction.txt".format(
#             index),
#         'C:\caocao\Python Project\ABC\structure\size200\Ecoli_size200_noise_0%-{0}_goldstandard.txt'.format(index),
#     )
#     struct.none_sign_analysis(detail=True)

for index in range(1, 21, 1):
    print("======================={}=======================".format(index))
    index_st = ''
    if index < 10:
        index_st = '0{}'.format(index)
    else:
        index_st = str(index)
    struct = StructuralAnalysis(
        "C:\caocao\Python Project\ABC\generated_data\size300\hybrid\\values\dbn_output\size300_{}_dbn_prediction.txt".format(
            index_st),
        'C:\caocao\Python Project\ABC\structure\size300\Ecoli_size300-{}_goldstandard.txt'.format(index),
    )
    struct.none_sign_analysis(detail=True)
