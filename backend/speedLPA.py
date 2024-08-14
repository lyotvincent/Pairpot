import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import make_interp_spline
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 14

df_sc = pd.read_csv("mistake-sc.tsv", sep='\t', index_col=0)
df_sp = pd.read_csv("mistake-sp.tsv", sep='\t', index_col=0)

fig, ax = plt.subplots(figsize=(3,3.5), constrained_layout=True)
df = pd.concat([df_sc['Time'], df_sp['Time']], axis=1)
df.columns = ['sc(n=3000)', 'sp(n=4598)']
bp = ax.boxplot(df, sym='.', patch_artist=True,widths=0.7)
plt.setp(bp['whiskers'], color='#666699')
plt.setp(bp['fliers'], markerfacecolor='#666699', marker='+', markersize=5)
for bplot in [bp]:
    for patch, color in zip(bplot['boxes'], ['pink','lightgreen']):
        patch.set_facecolor(color)
        patch.set_edgecolor('#666699')
ax.set_yticks(np.arange(0,11,2)*0.1)
ax.set_yticklabels(["{:.1f}".format(s) for s in np.arange(0,11,2)*0.1], fontsize=13)
ax.set_xticklabels(['single-cell', 'spatial'], fontsize=15)
ax.set_ylabel("Execution Time(s)")
ax.yaxis.grid(True, color='#e0e0e0')
ax.xaxis.grid(True, color='#e0e0e0')
ax.spines['top'].set_edgecolor('#e0e0e0')  # 设置上边坐标轴颜色
ax.spines['bottom'].set_edgecolor('#e0e0e0')  # 设置下边坐标轴颜色
ax.spines['left'].set_edgecolor('#e0e0e0')  # 设置左边坐标轴颜色
ax.spines['right'].set_edgecolor('#e0e0e0')
ax.tick_params(axis='both', color='#e0e0e0')
plt.savefig("speed.svg")