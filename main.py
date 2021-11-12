from BooleanInference.MIDNI.LoadNodes import LoadNodes
from BooleanInference.MIDNI.MIDNIMutualInformationCalculation import MIDNIMutualInformationCalculation
import time
import random
import _thread
from path import Path
from glob import glob
import os


def readInput(input_path):
    regu = []
    f = open(input_path, "r")
    reg = f.readline().split(' ')
    for i in reg:
        if i != '':
            regu.append(int(i))
    f.close()
    return regu


def regulatorGenes(size):
    l = []
    for i in range(size):
        l.append(random.randrange(1, 5))
    return l


class MIDNI:
    def __init__(self, gene_file):
        self.gene_file = gene_file
        self.root_path = os.path.dirname(gene_file)
        self.filename_none_extension = Path(gene_file).stem

    @staticmethod
    def save_prediction(path, predictions):
        network_size = len(predictions)
        dict = {}
        for i in range(network_size):
            for j in range(network_size):
                if i == j: continue
                key = "G{}G{}".format(i + 1, j + 1)
                values = "G{}\tG{}\t0".format(i + 1, j + 1)
                dict[key] = values
        relations = []
        for i in range(len(predictions)):
            for j in predictions[i]:
                key = "G{}G{}".format(j, i + 1)
                relations.append(key)
                dict[key] = "G{}\tG{}\t1".format(j, i + 1)
        non_relations = [item for item in dict.keys() if item not in relations]
        with open(path, 'w') as f:
            for i in relations:
                f.write('%s\n' % dict[i])
            for i in non_relations:
                f.write('%s\n' % dict[i])

    def directory_making(self, folder):
        p = Path(self.root_path)
        access = 0o755
        _path = "%s/%s" % (p, folder)
        try:
            if not os.path.exists(_path):
                os.mkdir(_path, access)
                return _path
        except OSError:
            print("Creation of the directory %s failed" % _path)
            return ""
        else:
            print("Successfully created the directory %s " % _path)
        return _path

    def execute(self):
        start_time = time.time()
        mif = MIDNIMutualInformationCalculation(self.gene_file)
        print("%s is starting to be processed...." % self.gene_file)
        print("system params are initializing...")
        mif.initializeTabularMI(print_out=False)
        path_input = "%s/input/%s_input.txt" % (self.root_path, self.filename_none_extension)
        # demoInputRegulatorGenes = readInput(path_input)
        output_path = self.directory_making("midni_output/")
        print("Target genes      regulatory genes")
        bow = []
        for i in range(mif.nodeSize):
            mif.my_executeFeatureSelection_original_implementation(i, 4, False)
            # mif.my_executeFeatureSelection_original_implementation(i, 6, False)
            print(i + 1, end='')
            print("\t\t\t\t\t", end='')
            result = mif.getResult()
            bow.append(result)
            arr = []
            for j in range(len(result)):
                print(result.__getitem__(j), end='')
                print("\t\t", end='')
            print("\t\t")
        end_time = time.time()
        path_output = '%s/%s_midni_output.txt' % (output_path, self.filename_none_extension)
        with open(path_output, 'w') as filehandle:
            for listitem in bow:
                for item in listitem:
                    filehandle.write('%s\t' % item)
                filehandle.write('\n')
        path_prediction = '%s/%s_midni_prediction.txt' % (output_path, self.filename_none_extension)
        MIDNI.save_prediction(path_prediction, bow)
        print("running time = ", 1000 * (end_time - start_time))
        print("[%s] has been executed completely ..." % self.gene_file)


folder = r"C:\caocao\gnw-master\Extra study\size 600\hybrid\values\*.txt"
#folder = r"C:\caocao\gnw-master\Extra study\size 500\hybrid\values\*.txt"
#folder = r"C:\caocao\Python Project\ABC\generated_data\size300\hybrid\values\*.txt"
# folder = "C:\caocao\gnw-master\Generated\size10_noise0%\Discritized/G_size10_noise0%_01.txt"
# for i in range(10, 210, 10):
#     folder = "C:\caocao\gnw-master\Generated\size" + str(i) + "_noise0%\Discritized/*.txt"
#     x = [Path(f).abspath() for f in glob(folder)]
#     for f in x:
#         p = MIDNI(f)
#         p.execute()
# for i in range(1,9):
#     folder = "C:\caocao\gnw-master\Incomming_links50\\%s\Discritized/*.txt" %i
x = [Path(f).abspath() for f in glob(folder)]
for f in x:
    p = MIDNI(f)
    p.execute()
