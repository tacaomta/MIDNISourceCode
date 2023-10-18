import random
import os
from path import Path
from glob import glob

":cvar path_genes đây là nơi lưu trữ file .csv"


class GeneratedInput:
    def __init__(self, path_genes, size):
        self.root_path = os.path.dirname(path_genes)
        self.filename_none_extension = Path(path_genes).stem
        self.size = size

    def directory_making(self, folder):
        p = Path(self.root_path)
        access = 0o755
        _path = "%s/%s" % (p, folder)
        try:
            if not os.path.exists(_path):
                os.makedirs(_path, access)
                return _path
        except OSError:
            print("Creation of the directory %s failed" % _path)
            return ""
        else:
            print("Successfully created the directory %s " % _path)
        return _path

    def generate(self):
        input_path = self.directory_making("/hybrid/values/input/")
        l = []
        for i in range(self.size):
            l.append(random.randrange(1, 5))
        path = "%s/%s_input.txt" % (input_path, self.filename_none_extension)
        with open(path, 'a') as f:
            for item in l:
                f.write("%s " % item)
        f.close()
        print("successfully generated...")


#folder = "C:/caocao/gnw-master/mydata/size200-noise-10percent/*.csv"
folder = "C:\caocao\gnw-master\Generated\size10_noise0%\Discritized/*.txt"
x = [Path(f).abspath() for f in glob(folder)]
for f in x:
    p = GeneratedInput(f, 10)
    p.generate()
