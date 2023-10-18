import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class Elbow:
    def __init__(self, path):
        self.path = path
        self.centroids = []
        self.data = pd.read_csv(path).values
        self.LengthOfGene = self.data.shape[0]
        self.NumberOfGenes = self.data.shape[1] - 1
        self.time = self.data[:, 0].reshape(-1, 1)

    def execute(self, out_path):
       # self.NumberOfGenes = 1
        for i in range(self.NumberOfGenes):
            x_value = self.data[:,  i+ 1].reshape(-1, 1)
            #x_value = Elbow.x_hat_calculation(x)
            maxValueOfThisGene = Elbow.gene_max(x_value)
            minValueOfThisGene = Elbow.gene_min(x_value)
            wcss = []
            list_validity = []
            for k in range(2, 4):
                df = pd.DataFrame({
                    'x': [self.time[v][0] for v in range(self.LengthOfGene)],
                    'y': [x_value[j][0] for j in range(self.LengthOfGene)]
                })
                """
                centroids initialization step 
                """
                interval = (maxValueOfThisGene - minValueOfThisGene) / k
                centroids = {
                    j + 1: [self.time_max() / 2, minValueOfThisGene + (j + 0.5) * interval]
                    for j in range(k)
                }
                #print("Khỏi tạo ================")
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

                # print("k ====================== %d" % k)
                # print(df)
                # print(centroids)
                intra_value = 0
                validity = -1
                wcss_item = 0
                for q in centroids.keys():
                    wcss_item += np.sum(df[df['closest'] == q]['x_distance_from_{}'.format(q)])
                    intra_value += np.sum(pow(df[df['closest'] == q]['x_distance_from_{}'.format(q)], 2))
                intra_value = intra_value / self.LengthOfGene
                wcss.append(wcss_item)
                if k == 1:
                    inter_value = centroids[1][1] ** 2
                else:
                    inter_value = (centroids[1][1] - centroids[2][1]) ** 2
              #  print("cuoi cung================k = %d" % k)
               # print(centroids)
                for w in centroids.keys():
                    for v in centroids.keys():
                        if w != v:
                            sub = (centroids[w][1] - centroids[v][1]) ** 2
                            if sub < inter_value:
                                inter_value = sub
                if inter_value != 0:
                    validity = intra_value / inter_value / inter_value
               # print("inter_value = %f" % inter_value)
                #print("validity = %f" % validity)
                list_validity.append(validity)
            print("validity")
            print(list_validity)
            k_optimal = Elbow.k_optimal(list_validity)
            # plt.clf()
            # plt.plot(range(1, 10), wcss)
            # for a, b in zip(range(1, 10), wcss):
            #     plt.text(a, b, str(round(b, 5)))
            # plt.title('Elbow Method - K = %d is optimal' % k_optimal)
            # plt.xlabel('Number of clusters')
            # plt.ylabel('WCSS')
            # full_path = '%s/gene_%d' % (out_path, i + 1)
            # plt.savefig(full_path)
            # text_display = '%s  -  elbow graph building-  gene %d' % (self.path, i + 1)
            # print(text_display)

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


p = Elbow("DREAM3/size100/trajectories04/Ecoli100-trajectories04.csv")
p.execute("DREAM3/size100/trajectories04/elbow_original")
