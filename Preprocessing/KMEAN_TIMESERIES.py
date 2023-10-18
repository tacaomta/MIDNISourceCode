import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


class K_Means_Cluster:
    def __init__(self, source_path, k):
        self.source_path = source_path
        self.K = k
        self.centroids = []
        self.data = pd.read_csv(source_path).values
        self.LengthOfGene = self.data.shape[0]
        self.NumberOfGenes = self.data.shape[1] - 1
        self.time = self.data[:, 0].reshape(-1, 1)
        self.DFS = []
        self.CP = []
        if k == 2:
            self.colour_map = {1: 'r', 2: 'g'}
            self.gene_bit = {1: '0', 2: '1'}
        else:
            self.colour_map = {1: 'r', 2: 'g', 3: 'b'}
            self.gene_bit = {1: '-1', 2: ' 0', 3: ' 1'}

    def visualization_input_save(self, path):
        for i in range(self.NumberOfGenes):
            x = self.data[:, i + 1].reshape(-1, 1)
            df = pd.DataFrame({
                'x': [self.time[k][0] for k in range(self.LengthOfGene)],
                'y': [x[j][0] for j in range(self.LengthOfGene)]
            })
            plt.clf()
            # fig = plt.figure(figsize=(5, 5))
            plt.scatter(df['x'], df['y'], color='r')
            plt.xlim(0, self.time_max())
            plt.ylim(0, K_Means_Cluster.gene_max(x))
            full_path = '%s/gene_%d' % (path, i + 1)
            plt.savefig(full_path)
            text_display = '%s  -  visualization_input_save  -  gene %d' % (self.source_path, i + 1)
            print(text_display)

    def elbow(self, path):
        for i in range(self.NumberOfGenes):
            x = self.data[:, i + 1].reshape(-1, 1)
            df = pd.DataFrame({
                'x': [self.time[k][0] for k in range(self.LengthOfGene)],
                'y': [x[j][0] for j in range(self.LengthOfGene)]
            })
            wcss = []
            plt.clf()
            for j in range(1, 10):
                kmeans = KMeans(n_clusters=j, init='k-means++', max_iter=300, n_init=10, random_state=0)
                kmeans.fit(df)
                wcss.append(kmeans.inertia_)
            plt.plot(range(1, 10), wcss)
            for a, b in zip(range(1, 10), wcss):
                plt.text(a, b, str(round(b, 2)))
            plt.title('Elbow Method')
            plt.xlabel('Number of clusters')
            plt.ylabel('WCSS')
            full_path = '%s/gene_%d' % (path, i + 1)
            plt.savefig(full_path)
            text_display = '%s  -  elbow -  gene %d' % (self.source_path, i + 1)
            print(text_display)

    def execute(self):
        for i in range(self.NumberOfGenes):
            x = self.data[:, i + 1].reshape(-1, 1)
            df = pd.DataFrame({
                'x': [self.time[k][0] for k in range(self.LengthOfGene)],
                'y': [x[j][0] for j in range(self.LengthOfGene)]
            })
            """
            centroids initialization step 
            """
            maxValueOfThisGene = K_Means_Cluster.gene_max(x)
            minValueOfThisGene = K_Means_Cluster.gene_min(x)
            interval = (maxValueOfThisGene - minValueOfThisGene) / self.K
            centroids = {
                # i + 1: [self.LengthOfGene / 2, np.random.uniform(i * interval, (i + 1) * interval)]
                j + 1: [self.time_max() / 2, minValueOfThisGene + (j + 0.5) * interval]
                for j in range(self.K)
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
                    self.DFS.append(df)
                    self.CP.append(centroids)
                    break

    def visualization_result(self, path):
        for i in range(self.NumberOfGenes):
            x = self.data[:, i + 1].reshape(-1, 1)
            plt.clf()
            plt.scatter(self.DFS[i]['x'], self.DFS[i]['y'],
                        color=self.DFS[i]['color'], alpha=0.5, edgecolor='k')
            for j in self.CP[i].keys():
                plt.scatter(*self.CP[i][j], color=self.colour_map[j])
            plt.xlim(0, self.time_max())
            plt.ylim(0, K_Means_Cluster.gene_max(x))
            full_path = '%s/gene_%d' % (path, i + 1)
            plt.savefig(full_path)
            text_display = '%s  -  result visualization -  gene %d' % (self.source_path, i + 1)
            print(text_display)

    def export_2_text(self, path, with_header=True):
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
        centroid_distance_cols = ['distance_from_{}'.format(k) for k in centroids.keys()]
        df['closest'] = df.loc[:, centroid_distance_cols].idxmin(axis=1)
        df['closest'] = df['closest'].map(lambda x: int(x.lstrip('distance_from_')))
        df['color'] = df['closest'].map(lambda x: self.colour_map[x])
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


p = K_Means_Cluster("DREAM3/size100/trajectories01/Ecoli100-trajectories01.csv", 3)
p.visualization_input_save("DREAM3/size100/trajectories01/original")
p.elbow("DREAM3/size100/trajectories01/elbow")
p.execute()
p.visualization_result("DREAM3/size100/trajectories01/kmean3")
p.export_2_text("DREAM3/size100/trajectories01/result/trajectories01_k3.txt")
