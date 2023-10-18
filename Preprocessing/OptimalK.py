import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import os
import math


class OptimalK:
    """
    Automatically to determine the optimal number of cluster K
    """

    def __init__(self, source_path, k_max=4, original_or_mean=True):
        self.source_path = source_path
        self.data = pd.read_csv(source_path).values
        self.LengthOfGene = self.data.shape[0]
        self.NumberOfGenes = self.data.shape[1] - 1
        self.time = self.data[:, 0].reshape(-1, 1)
        self.k_max = k_max
        self.original_or_mean = original_or_mean

    def execute(self):
        """:creation directory for saving result"""
        mean_path = self.directory_making("mean_original")
        validity_path = self.directory_making("validity")
        elbow_path = self.directory_making("elbow_with_optimal_decision")
        """:mean subtract"""
        x_hat_list = []
        list_optimal_k = []
        for i in range(self.NumberOfGenes):
            x = self.data[:, i + 1].reshape(-1, 1)
            x_hat = OptimalK.x_hat_calculation(x)
            if self.original_or_mean:
                x_value = x
            else:
                x_value = x_hat
            maxValueOfThisGene = OptimalK.gene_max(x_value)
            minValueOfThisGene = OptimalK.gene_min(x_value)
            wcss = []
            list_validity = []
            for k in range(1, 10):
                df = pd.DataFrame({
                    'x': [self.time[v][0] for v in range(self.LengthOfGene)],
                    'y': [OptimalK.normalization(x_value[j][0], minValueOfThisGene, maxValueOfThisGene) for j in
                          range(self.LengthOfGene)]
                })
                """
                centroids initialization step 
                """
                interval = (maxValueOfThisGene - minValueOfThisGene) / k
                centroids = {
                    j + 1: [100, minValueOfThisGene + (j + 0.5) * interval]
                    for j in range(k)
                }
                """
                assignment of points to each centroids
                """
                df = OptimalK.assignment(centroids, df)
                """
                fit step
                """
                while True:
                    closest_centroids = df['closest'].copy(deep=True)
                    for j in centroids.keys():
                        centroids[j][1] = np.mean(df[df['closest'] == j]['y'])
                    df = OptimalK.assignment(centroids, df)
                    if closest_centroids.equals(df['closest']):
                        break
                # calculate intra distance, here
                wcss_item = 0
                intra_value = 0
                validity = -1
                for q in centroids.keys():
                    wcss_item += np.sum(df[df['closest'] == q]['x_distance_from_{}'.format(q)])
                    intra_value += np.sum(pow(df[df['closest'] == q]['x_distance_from_{}'.format(q)], 2))
                intra_value = intra_value / self.LengthOfGene
                wcss.append(wcss_item)
                # calculate inter distance
                inter_value = centroids[1][1] ** 2
                for w in centroids.keys():
                    for v in centroids.keys():
                        if w != v:
                            sub = (centroids[w][1] - centroids[v][1]) ** 2
                            if sub < inter_value:
                                inter_value = sub
                if inter_value != 0:
                    validity = intra_value / inter_value
                list_validity.append(validity)
            # find the optimal k based on the list of validity
            k_optimal = OptimalK.k_optimal(list_validity)
            full_path = '%s/gene_%d.txt' % (validity_path, i + 1)
            np.savetxt(full_path, list_validity)
            text_display = '%s  -  export validity values of gene %d' % (
                self.source_path, i + 1)
            print(text_display)
            # list_k_optimal.append(k_optimal)
            plt.clf()
            plt.plot(range(1, 10), wcss)
            for a, b in zip(range(1, 10), wcss):
                plt.text(a, b, str(round(b, 5)))
            plt.title('Elbow Method - K = %d is optimal' % k_optimal)
            plt.xlabel('Number of clusters')
            plt.ylabel('WCSS')
            full_path = '%s/elbow_%d' % (elbow_path, i + 1)
            plt.savefig(full_path)
            text_display = '%s  -  elbow graph building for the chosen gene %d' % (
                self.source_path, i + 1)
            print(text_display)

    def directory_making(self, folder):
        p = Path(self.source_path)
        path = "%s/%s" % (p.parent, folder)
        try:
            if not os.path.exists(path):
                os.mkdir(path)
                return path
        except OSError:
            print("Creation of the directory %s failed" % path)
            return ""
        else:
            print("Successfully created the directory %s " % path)
        return path

    """
    This procedure gives values by subtracting mean value from original values of the input data
    """

    @staticmethod
    def x_hat_calculation(x):
        x_mean = np.mean(x)
        x_hat = []
        for i in range(len(x)):
            x_hat.append(x[i] - x_mean)
        return x_hat

    """:The procedure calculates the distance between two matrices"""

    @staticmethod
    def distance(m, n):
        d = 0
        for i in range(len(m)):
            d += (m[i] - n[i]) ** 2
        return math.sqrt(d)

    @staticmethod
    def normalization(x, min, max):
        if max == min:
            return min
        return (x - min) / (max - min)

    @staticmethod
    def assignment(centroids, df):
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
    def k_optimal(x):
        min = x[0]
        optimal = 10
        for i in range(len(x)):
            if min >= x[i]:
                min = x[i]
                optimal = i + 1
        return optimal


c = OptimalK("DREAM3/size100/trajectories04/Ecoli100-trajectories04.csv", 4, False)
c.execute()
