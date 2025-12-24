import dill
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

ax = dill.load(open('/home/nap24/Documents/Pipeline/ccgd/STS0008b/160815_Cutaneous_Squam_B2_STS0008b_MutationsForCoMutationPlot_chart.pickle', 'rb'))
figure(num=1, figsize=(8, 6), dip=80, facecolor='w', edgecolor='k')
plt.savefig("chart.png")
