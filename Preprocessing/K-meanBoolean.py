import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import os
from glob import glob
from path import Path

""":cvar
    The last version of the discretization process
"""


class Discretization:
    def __init__(self, path):
        self.path = path
        self.centroids = []
        self.data = pd.read_csv(path).values
        self.LengthOfGene = self.data.shape[0]
        self.NumberOfGenes = self.data.shape[1] - 1
        self.time = self.data[:, 0].reshape(-1, 1)
        self.colmap = {1: 'r', 2: 'g'}
        self.bit = {1: '0', 2: '1'}
        self.root_path = os.path.dirname(path)
        self.filename_none_extension = Path(path).stem

    def execute(self, with_header=True):
        discretization_graph_path = self.directory_making("k_mean2/graphs/%s/" % self.filename_none_extension)
        discretization_gene_path = self.directory_making("k_mean2/values/")
        result = pd.DataFrame({

        })
        header = ""
        for i in range(self.NumberOfGenes):
            x_value = self.data[:, i + 1].reshape(-1, 1)
            k = 2
            df = self.K_mean(i, x_value, k, discretization_graph_path)
            result['G%d' % (i + 1)] = df['closest'].map(lambda j: self.bit[j])
            if i == 0:
                header = 'G%d' % (i + 1)
            else:
                header = '%s\tG%d' % (header, i + 1)
        outpath = '%s/%s_k2.txt' % (discretization_gene_path, self.filename_none_extension)
        if with_header:
            np.savetxt(outpath, result, delimiter='\t', header=header, fmt='%s')
        else:
            np.savetxt(outpath, result, delimiter='\t', fmt='%s')
        text_display = '%s  -  export_2_text -  done' % self.path
        print(text_display)

    def K_mean(self, index, x, k, save_path):

        maxValueOfThisGene = Discretization.gene_max(x)
        minValueOfThisGene = Discretization.gene_min(x)
        wcss = []
        list_validity = []
        df = pd.DataFrame({
            'x': [self.time[v][0] for v in range(self.LengthOfGene)],
            'y': [x[j][0] for j in range(self.LengthOfGene)]
        })
        """
        centroids initialization step 
        """
        interval = (maxValueOfThisGene - minValueOfThisGene) / k
        centroids = {
            j + 1: [self.time_max() / 2, minValueOfThisGene + (j + 0.5) * interval]
            for j in range(k)
        }
        # print("Khỏi tạo ================")
        # print(centroids)

        """
        assignment of points to each centroids
        """
        df = self.assignment(centroids, df)
        """
        fit step
        """
        while True:
            closest_centroids = df['closest'].copy(deep=True)
            for j in centroids.keys():
                centroids[j][1] = np.mean(df[df['closest'] == j]['y'])
            df = self.assignment(centroids, df)
            if closest_centroids.equals(df['closest']):
                break
        plt.clf()
        plt.figure(figsize=(5, 5))
        plt.scatter(df['x'], df['y'], color=df['color'], alpha=0.5, edgecolor='k')
        for i in centroids.keys():
            plt.scatter(*centroids[i], color=self.colmap[i])
        plt.xlim(0, self.time_max())
        plt.ylim(0, self.gene_max(x))

        full_path = '%s/gene_%d' % (save_path, index + 1)
        plt.savefig(full_path)
        text_display = '%s  - save discretized figure  -  gene %d' % (self.path, index + 1)
        print(text_display)
        return df

    def export_2_text(self, path, df, with_header=True):
        result = pd.DataFrame({

        })
        header = ""
        for i in range(self.NumberOfGenes):
            result['G%d' % (i + 1)] = self.DFS[i]['closest'].map(lambda j: self.gene_bit[j])
            if i == 0:
                header = 'G%d' % (i + 1)
            else:
                header = '%s\tG%d' % (header, i + 1)
        if with_header:
            np.savetxt(path, result, delimiter='\t', header=header, fmt='%s')
        else:
            np.savetxt(path, result, delimiter='\t', fmt='%s')
        text_display = '%s  -  export_2_text -  done' % self.source_path
        print(text_display)

    def assignment(self, centroids, df):
        for j in centroids.keys():
            df['distance_from_{}'.format(j)] = (
                np.sqrt(
                    (df['y'] - centroids[j][1]) ** 2
                )
            )
            df['x_distance_from_{}'.format(j)] = (
                    (df['y'] - centroids[j][1]) ** 2
            )
        centroid_distance_cols = ['distance_from_{}'.format(k) for k in centroids.keys()]
        df['closest'] = df.loc[:, centroid_distance_cols].idxmin(axis=1)
        df['closest'] = df['closest'].map(lambda x: int(x.lstrip('distance_from_')))
        df['color'] = df['closest'].map(lambda x: self.colmap[x])
        return df

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

    def time_max(self):
        max = 0
        for i in range(self.LengthOfGene):
            if max < self.time[i][0]:
                max = self.time[i][0]
        return max

    @staticmethod
    def gene_max(x):
        max = 0
        for i in range(len(x)):
            if max < x[i][0]:
                max = x[i][0]
        return max

    @staticmethod
    def gene_min(x):
        min = 100000000
        for i in range(len(x)):
            if min > x[i][0]:
                min = x[i][0]
        return min

    @staticmethod
    def x_hat_calculation(x):
        x_mean = np.mean(x)
        x_hat = []
        for i in range(len(x)):
            x_hat.append(x[i] - x_mean)
        return x_hat

    @staticmethod
    def k_optimal(x):
        min = x[0]
        optimal = 10
        for i in range(len(x)):
            if x[i] > 0:
                if min >= x[i]:
                    min = x[i]
                    optimal = i + 1
        return optimal


folder = "C:\caocao\gnw-master\Extra study\size 600/*.csv"
x = [Path(f).abspath() for f in glob(folder)]
for f in x:
    p = Discretization(f)
    p.execute(with_header=False)


# p = Discretization("DREAM3/size100_5/Ecoli_size100-5.csv")
# p.execute(with_header=False)
