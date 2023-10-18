import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


class K_Means_Cluster:
    def __init__(self, source_path, target_path, k, index):
        self.source_path = source_path
        self.target_path = target_path
        self.K = k
        self.centroids = []
        data = pd.read_csv(source_path).values
        self.LengthOfGene = data.shape[0]
        self.NumberOfGenes = data.shape[1]
        self.index = index
        x = data[:, self.index].reshape(-1, 1)
        max = 0
        for i in range(self.LengthOfGene):
            if max < x[i][0]:
                max = x[i][0]
        self.maxValue = max
        self.df = pd.DataFrame({
            'x': [i + 1 for i in range(self.LengthOfGene)],
            'y': [x[j][0] for j in range(self.LengthOfGene)]
        })
        self.colmap = {1: 'r', 2: 'g', 3: 'b', 4: 'y'}

    def visualization_input(self):
        fig = plt.figure(figsize=(5, 5))
        plt.scatter(self.df['x'], self.df['y'], color='r')
        plt.xlim(0, self.LengthOfGene)
        plt.ylim(0, self.maxValue)
        plt.show()
        print(self.NumberOfGenes)

    def elbow(self, path):
        wcss = []
        plt.clf()
        for i in range(1, 10):
            kmeans = KMeans(n_clusters=i, init='k-means++', max_iter=300, n_init=10, random_state=0)
            kmeans.fit(self.df)
            wcss.append(kmeans.inertia_)
        plt.plot(range(1, 10), wcss)
        plt.title('Elbow Method')
        plt.xlabel('Number of clusters')
        plt.ylabel('WCSS')
        plt.savefig(path)
        #plt.show()

    def centroids_init(self):
        interval = self.maxValue / self.K
        self.centroids = {
            # i + 1: [self.LengthOfGene / 2, np.random.uniform(i * interval, (i + 1) * interval)]
            i + 1: [self.LengthOfGene / 2, (i + 0.5) * interval]
            for i in range(self.K)
        }
       # print(self.centroids)

    def assignment(self):
        for i in self.centroids.keys():
            self.df['distance_from_{}'.format(i)] = (
                np.sqrt(
                    (self.df['y'] - self.centroids[i][1]) ** 2
                )
            )
        centroid_distance_cols = ['distance_from_{}'.format(i) for i in self.centroids.keys()]
        self.df['closest'] = self.df.loc[:, centroid_distance_cols].idxmin(axis=1)
        self.df['closest'] = self.df['closest'].map(lambda x: int(x.lstrip('distance_from_')))
        self.df['color'] = self.df['closest'].map(lambda x: self.colmap[x])
        return self.df

    def update(self):
        for i in self.centroids.keys():
            self.centroids[i][1] = np.mean(self.df[self.df['closest'] == i]['y'])

    def fit(self):
        while True:
            closest_centroids = self.df['closest'].copy(deep=True)
            self.update()
            self.df = self.assignment()
            if closest_centroids.equals(self.df['closest']):
                break

    def visualization_result(self, path):
        plt.clf()
        fig = plt.figure(figsize=(5, 5))
        plt.scatter(self.df['x'], self.df['y'], color=self.df['color'], alpha=0.5, edgecolor='k')
        for i in self.centroids.keys():
            plt.scatter(*self.centroids[i], color=self.colmap[i])
        plt.xlim(0, self.LengthOfGene)
        plt.ylim(0, self.maxValue)
        plt.savefig(path)
        # plt.show()

    def print(self):
        print(self.df)

    def path(self):
        return p.target_path


for i in range(2810):
    p = K_Means_Cluster("DREAM5/net2/net2_expression_data_avg.csv", "DREAM5/net2/elbow", 4, i)
    #p.visualization_input()
    # p.elbow()
   # p.centroids_init()
   # p.assignment()
    #p.fit()
   # name = '%s/fig_gene_%d' % (p.path(), i + 1)
    name = '%s/elbow_gene_%d' % (p.path(), i + 1)
    p.elbow(name)
    #p.visualization_result(name)
    print("=====================", i + 1)
    # p.print()
