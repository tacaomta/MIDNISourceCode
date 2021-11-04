from sklearn.metrics import mutual_info_score
from sklearn.metrics import adjusted_mutual_info_score
import math
import random
from scipy.stats import entropy

a = [0, 0, 0, 0, 1, 1, 1, 1]
b = [0, 0, 1, 1, 0, 0, 1, 1]
c = [0, 1, 0, 1, 0, 1, 0, 1]
aph = [0, 0, 1, 1, 0, 0, 1, 1]
bph = [0, 1, 0, 1, 1, 1, 1, 1]
cph = [0, 0, 0, 1, 0, 1, 1, 1]
def regulatorGenes(size):
    l = []
    for i in range(size):
        l.append(random.randrange(0, 3))
    return l
t = regulatorGenes(2000)
h = regulatorGenes(2000)
print(t)
print(t[1:])
print(h)
mul = mutual_info_score(t[1:], h[1:])
print(mul)
