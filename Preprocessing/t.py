import matplotlib
import matplotlib.pyplot as plt
import numpy as np


labels = ['G1', 'G2', 'G3', 'G4', 'G5']
top1 = [20, 34, 30, 35, 27]
top2 = [25, 32, 34, 20, 25]
top3 = [15, 27, 19, 22, 20]
y = top1[0:5]
print(y)
x = np.arange(len(labels))  # the label locations
width = 0.3  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x, y, width, label='top1')
#rects2 = ax.bar(2, top2, width, label='top2')
#ax.bar()
#rects3 = ax.bar(3, top3, width, label='top3')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Scores')
ax.set_title('Scores by group and gender')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


#autolabel(rects1)
#autolabel(rects2)
#autolabel(rects3)

fig.tight_layout()

plt.show()

fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
b = []
a = 0.01245124555555554
b.append(round(a, 3))
print(b)