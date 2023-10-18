import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math


class Elbow:
    def __init__(self, path):
        self.path = path
        self.centroids = []
        self.data = pd.read_csv(path).values
        self.LengthOfGene = self.data.shape[0]
        self.NumberOfGenes = self.data.shape[1] - 1
        self.time = self.data[:, 0].reshape(-1, 1)

    @staticmethod
    def distance(m, n):
        d = 0
        for i in range(len(m)):
            d += (m[i] - n[i]) ** 2
        return math.sqrt(d)

    def gene_distance(self, path):
        """
        Initialization data frame
        """
        columns = []
        for i in range(self.NumberOfGenes):
            columns.append("G%d" % (i + 1))
        df = pd.DataFrame([[0 for x in range(self.NumberOfGenes)] for j in range(self.NumberOfGenes)],
                          columns=list(columns))
        df.insert(0, "", columns, True)
        """
        Calculate the distance between each pair of genes
        """
        for i in range(self.NumberOfGenes):
            for j in range(i + 1, self.NumberOfGenes):
                m = self.data[:, i + 1].reshape(-1, 1)
                n = self.data[:, j + 1].reshape(-1, 1)
                col_label = "G%d" % (j + 1)
                df.iloc[i, df.columns.get_loc(col_label)] = Elbow.distance(m, n)
        full_path = '%s/gene_distance.csv' % path
        df.to_csv(full_path, index=False, header=True)

    def max_min_difference(self, path):
        dif = []
        dx = [p + 1 for p in range(self.NumberOfGenes)]
        label = ['Gene %d' % (p + 1) for p in range(self.NumberOfGenes)]
        max_dif = 0
        min_dif = 1000000
        for i in range(self.NumberOfGenes):
            x = self.data[:, i + 1].reshape(-1, 1)
            maxValueOfThisGene = Elbow.gene_max(x)
            minValueOfThisGene = Elbow.gene_min(x)
            dif.append(maxValueOfThisGene - minValueOfThisGene)
            if max_dif < (maxValueOfThisGene - minValueOfThisGene):
                max_dif = maxValueOfThisGene - minValueOfThisGene
            if min_dif > (maxValueOfThisGene - minValueOfThisGene):
                min_dif = maxValueOfThisGene - minValueOfThisGene
        plt.plot(range(1, self.NumberOfGenes + 1), dif)
        for a, b in zip(range(1, self.NumberOfGenes + 1), dif):
            plt.text(a, b, str(round(b, 2)))
        plt.title('max_min_gap of genes')
        plt.xlabel('Gene')
        plt.ylabel('Gap values')
        plt.legend(['max = %f, min = %f' % (round(max_dif, 2), round(min_dif, 2))])

        full_path = '%s/max_min_gap' % path
        plt.savefig(full_path)

    @staticmethod
    def normalization(x, min, max):
        if max == min:
            return min
        return (x - min) / (max - min)

    def visualization_input_save(self, path):
        for i in range(self.NumberOfGenes):
            x = self.data[:, i + 1].reshape(-1, 1)
            maxValueOfThisGene = Elbow.gene_max(x)
            minValueOfThisGene = Elbow.gene_min(x)
            df = pd.DataFrame({
                'x': [self.time[v][0] for v in range(self.LengthOfGene)],
                'y': [Elbow.normalization(x[j][0], minValueOfThisGene, maxValueOfThisGene) for j in
                      range(self.LengthOfGene)]
            })
            plt.clf()
            plt.scatter(df['x'], df['y'], color='r')
            plt.xlim(0, 200)
            plt.ylim(0, 1)
            full_path = '%s/gene_%d' % (path, i + 1)
            plt.savefig(full_path)
            text_display = '%s  -  visualization_input_save  -  gene %d' % (self.path, i + 1)
            print(text_display)

    def execute(self, out_path):
        # self.NumberOfGenes = 1
        for i in range(self.NumberOfGenes):
            x = self.data[:, i + 1].reshape(-1, 1)
            maxValueOfThisGene = Elbow.gene_max(x)
            minValueOfThisGene = Elbow.gene_min(x)
            wcss = []
            for k in range(1, 10):
                df = pd.DataFrame({
                    'x': [self.time[v][0] for v in range(self.LengthOfGene)],
                    'y': [Elbow.normalization(x[j][0], minValueOfThisGene, maxValueOfThisGene) for j in
                          range(self.LengthOfGene)]
                })
                """
                centroids initialization step 
                """
                interval = 1 / k
                centroids = {
                    j + 1: [self.time_max() / 2, (j + 0.5) * interval]
                    for j in range(k)
                }
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

                # print("k ====================== %d" % k)
                # print(df)
                # print(centroids)
                wcss_item = 0
                for q in centroids.keys():
                    wcss_item += np.sum(df[df['closest'] == q]['x_distance_from_{}'.format(q)])
                wcss.append(wcss_item)

            plt.clf()
            plt.plot(range(1, 10), wcss)
            for a, b in zip(range(1, 10), wcss):
                plt.text(a, b, str(round(b, 5)))
            plt.title('Elbow Method')
            plt.xlabel('Number of clusters')
            plt.ylabel('WCSS')
            full_path = '%s/gene_%d' % (out_path, i + 1)
            plt.savefig(full_path)
            text_display = '%s  -  elbow graph building-  gene %d' % (self.path, i + 1)
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
        return df

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


p = Elbow("DREAM3/size100/trajectories04/Ecoli100-trajectories04.csv")
p.gene_distance("DREAM3/size100/trajectories04/gene_distance")
# p.max_min_difference("DREAM3/size100/trajectories04/max_min_gap")
# p.visualization_input_save("DREAM3/size100/trajectories04/original_normalization")
# p.execute("DREAM3/size100/trajectories04/elbow_normalization")
