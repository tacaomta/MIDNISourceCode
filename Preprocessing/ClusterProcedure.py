import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import os
import math


class MyCluster:
    """
    MyCluster that discretizes GE into k clusters
    source_path - data file path
    threshold_distance - threshold distance values between two genes
    threshold_gap_max_min - threshold for max min gap that larger this to consider for discretization
    """

    def __init__(self, source_path, threshold_distance=0.4, threshold_gap_max_min=0.3, min_top=1, k_max=4):
        self.source_path = source_path
        self.threshold_distance = threshold_distance
        self.threshold_gap_max_min = threshold_gap_max_min
        self.min_top = min_top
        self.data = pd.read_csv(source_path).values
        self.LengthOfGene = self.data.shape[0]
        self.NumberOfGenes = self.data.shape[1] - 1
        self.time = self.data[:, 0].reshape(-1, 1)
        self.k_max = k_max

    """
    Step by step
    Step 1 - normalization by using formula g^ = g - g*
    Step 2 - calculation distance between each pair of genes
    Step 3 - neighbor calculation - count how many neighbors each point has, then choose the point
    that has the largest number of neighbors as the center point
    Step 4 - combination with ma_min gap to consider that gene will be normalized or not
    Step 5 - K mean algorithm apply 
    """

    def execute(self):
        """:creation directory for saving result"""
        mean_path = self.directory_making("mean_original")
        distance_path = self.directory_making("gene_distance_new")
        validity_path = self.directory_making("validity")
        optimal_path = self.directory_making("optimal_path")
        """:mean subtract"""
        x_hat_list = []
        for i in range(self.NumberOfGenes):
            x = self.data[:, i + 1].reshape(-1, 1)
            x_hat = MyCluster.x_hat_calculation(x)
            x_hat_list.append(x_hat)
            plt.clf()
            plt.scatter(self.time, x_hat, color='r')
            full_path = '%s/gene_%d' % (mean_path, i + 1)
            plt.savefig(full_path)
            text_display = '%s  -  mean subtract and graph building  -  gene %d' % (self.source_path, i + 1)
            print(text_display)
        """:Calculation the distance between two genes"""
        """:=========================================="""
        """:Initialization data frame================="""
        columns = []
        for i in range(self.NumberOfGenes):
            columns.append("G%d" % (i + 1))
        column_original = columns.copy()
        for i in range(self.min_top):
            columns.append("min_top_%d" % (i + 1))
        df = pd.DataFrame([[0 for x in range(self.NumberOfGenes + self.min_top)] for j in range(self.NumberOfGenes)],
                          columns=list(columns))
        df.insert(0, "", column_original, True)
        """
        Calculate the distance between each pair of genes
        """
        histogram_top_x = self.directory_making("histogram_top_x")
        nearest_neighbor = []
        count = 0
        list_k_optimal = []
        for i in range(self.NumberOfGenes):
            count = 0
            min_list = []
            # for j in range(i + 1, self.NumberOfGenes):
            for j in range(self.NumberOfGenes):
                m = x_hat_list[i]
                n = x_hat_list[j]
                col_label = "G%d" % (j + 1)
                distance = MyCluster.distance(m, n)
                # if distance != 0:
                min_list.append(distance)
                df.iloc[i, df.columns.get_loc(col_label)] = distance
                if 0 < distance <= self.threshold_distance:
                    count += 1
            nearest_neighbor.append(count)
            min_list_sort = min_list.copy()
            min_list_sort.sort()
            """ 
            Build histogram distance of top x minimums
            """
            histogram_min_list_sort = min_list_sort.copy()
            plt.clf()
            fig, ax = plt.subplots()
            labels = ['top1', 'top2', 'top3', 'top4', 'top5', 'top6', 'top7', 'top8', 'top9', 'top10']
            x = np.arange(len(labels))  # the label locations
            width = 0.2
            rects1 = ax.bar(x, histogram_min_list_sort[1:11], width, label='top10')
            ax.set_ylabel('Distance')
            ax.set_title('Top 10 distance minimums')
            ax.set_xticks(x)
            ax.set_xticklabels(labels)
            ax.legend()
            full_path = '%s/top10_gene_%d' % (histogram_top_x, i + 1)
            plt.savefig(full_path)
            text_display = '%s  - export top distance 10 minimums of gene %d' % (self.source_path, i + 1)
            print(text_display)
            """ End of histogram"""
            min_list_sort.pop(0)
            for m in range(self.min_top):
                if len(min_list_sort) != 0:
                    dis = min_list_sort.pop(0)
                    ref = min_list.index(dis)
                    val = 'vs G%d - %f' % (ref + 1, dis)
                else:
                    val = '='
                df.iloc[i, df.columns.get_loc('min_top_%d' % (m + 1))] = val
        full_path = '%s/gene_distance.csv' % distance_path
        df.to_csv(full_path, index=False, header=True)
        text_display = '%s  - export distance to excel file  -  done' % self.source_path
        print(text_display)

        """Find the center point that has the largest number of neighbors which their distances smaller then 
        threshold """
        nearest_path = self.directory_making("nearest_genes")
        plt.clf()
        plt.plot(range(1, self.NumberOfGenes + 1), nearest_neighbor, color='g')
        plt.legend(
            ['threshold = %f, Gene - %d' % (
                self.threshold_distance, 1 + nearest_neighbor.index(np.max(nearest_neighbor)))])
        for a, b in zip(range(1, self.NumberOfGenes + 1), nearest_neighbor):
            plt.text(a, b, str(round(b, 2)))
        full_path = '%s/nearest_genes' % nearest_path
        plt.savefig(full_path)
        text_display = '%s  -  export the figure expressed the nearest neighbors to a gene  -  done' % self.source_path
        print(text_display)

        """From index_chosen_gene we start apply max_min gap"""
        nearest_neighbor_copy = nearest_neighbor.copy()
        nearest_neighbor_copy.sort()
        list_of_five_chosen_genes = [nearest_neighbor_copy[i] for i in
                                     range(len(nearest_neighbor_copy) - 5, len(nearest_neighbor_copy))]
        slice_index_start = -1
        index_chosen_gene = -1
        while slice_index_start >= -5:
            index_chosen_gene = nearest_neighbor.index(list_of_five_chosen_genes[slice_index_start])
            original_x = self.data[:, index_chosen_gene + 1].reshape(-1, 1)
            if np.max(original_x) - np.min(original_x) < self.threshold_gap_max_min:
                index_chosen_gene = -1
            else:
                break
            slice_index_start -= 1

        if index_chosen_gene == -1:
            print("Process exist without a chosen gene for discretization of this trajectory")
        else:
            print("Gene %s has been chosen for the next step" % (index_chosen_gene + 1))
            chosen_x = self.data[:, index_chosen_gene + 1].reshape(-1, 1)
            maxValueOfThisGene = MyCluster.gene_max(chosen_x)
            minValueOfThisGene = MyCluster.gene_min(chosen_x)
            wcss = []
            list_validity = []
            for k in range(1, 10):
                df = pd.DataFrame({
                    'x': [self.time[v][0] for v in range(self.LengthOfGene)],
                    'y': [MyCluster.normalization(chosen_x[j][0], minValueOfThisGene, maxValueOfThisGene) for j in
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
                df = MyCluster.assignment(centroids, df)
                """
                fit step
                """
                while True:
                    closest_centroids = df['closest'].copy(deep=True)
                    for j in centroids.keys():
                        centroids[j][1] = np.mean(df[df['closest'] == j]['y'])
                    df = MyCluster.assignment(centroids, df)
                    if closest_centroids.equals(df['closest']):
                        break
                # calculate intra distance, here
                wcss_item = 0
                intra_value = 0
                validity = -1
                for q in centroids.keys():
                    wcss_item += np.sum(df[df['closest'] == q]['x_distance_from_{}'.format(q)])
                    intra_value += np.sum(pow(df[df['closest'] == q]['x_distance_from_{}'.format(q)], 2))
                wcss.append(wcss_item)
                intra_value = intra_value / self.LengthOfGene
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
                list_validity.append(round(validity, 3))
            # find the optimal k based on the list of validity
            k_optimal = MyCluster.k_optimal(list_validity)
            print("The most optimal K for chosen gene %d  =  %d" % (index_chosen_gene + 1, k_optimal))
            full_path = '%s/gene_%d.txt' % (validity_path, index_chosen_gene + 1)
            np.savetxt(full_path, list_validity)
            text_display = '%s  -  export validity values of gene %d' % (
                self.source_path, index_chosen_gene + 1)
            print(text_display)
            # list_k_optimal.append(k_optimal)
            elbow_path = self.directory_making("elbow_best_gene")
            plt.clf()
            plt.plot(range(1, 10), wcss)
            for a, b in zip(range(1, 10), wcss):
                plt.text(a, b, str(round(b, 5)))
            plt.title('Elbow Method')
            plt.xlabel('Number of clusters')
            plt.ylabel('WCSS')
            full_path = '%s/elbow_%d' % (elbow_path, index_chosen_gene + 1)
            plt.savefig(full_path)
            text_display = '%s  -  elbow graph building for the chosen gene %d' % (
                self.source_path, index_chosen_gene + 1)
            print(text_display)
        # full_path = '%s/optimal.txt' % optimal_path
        # np.savetxt(full_path, list_k_optimal)
        # text_display = '%s  -  export the list of optimal k-s' % self.source_path
        # print(text_display)

    """:Creation a new directory for saving results"""

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
        for i in range(len(x)):
            if min > x[i]:
                min = x[i]
                optimal = i + 1
        return optimal


c = MyCluster("DREAM3/size100/trajectories04/Ecoli100-trajectories04.csv", 0.2, 0.13, 10, 4)
c.execute()
