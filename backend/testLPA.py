import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import make_interp_spline
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 14


df = pd.read_csv("mistake-sp.tsv", sep='\t', index_col=0)
mistakes = pd.DataFrame()
mistakes_wo = pd.DataFrame()
groups = df.groupby('MR')
for group in groups:
    mistakes[group[0]] = group[1]['ARI_r'].reset_index(drop=True)
    mistakes_wo[group[0]] = group[1]['ARI_o'].reset_index(drop=True)
x = np.arange(len(mistakes.columns))+1
y = list(mistakes.mean(axis=0))
y_wo = list(mistakes_wo.mean(axis=0))
x_smooth = np.linspace(x.min(), x.max(), 200)
y_smooth = make_interp_spline(x, y)(x_smooth)
y_smooth_wo = make_interp_spline(x, y_wo)(x_smooth)
fig, ax = plt.subplots(figsize=(4,3.5), constrained_layout=True)
bp = ax.boxplot(mistakes, sym='.', patch_artist=True)
plt.setp(bp['whiskers'], color='#666699')
plt.setp(bp['fliers'], markerfacecolor='#666699', markeredgecolor='white', marker='.', markersize=5)
colors = ['pink' if s > 0.7 else 'lightgrey' for s in y]
for bplot in [bp]:
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('#666699')
l=ax.plot(x_smooth, y_smooth, label="LassoView")
l_wo=ax.plot(x_smooth, y_smooth_wo, label="LassoView w/o",c='#4ab311')
ax.set_yticks(np.arange(0,11,2)*0.1)
ax.set_yticklabels(["{:.1f}".format(s) for s in np.arange(0,11,2)*0.1], fontsize=13)
ax.set_xticklabels(list(mistakes.columns), fontsize=13)
ax.set_xticks(range(-1,20,2))
ax.set_xlim(0,21)
ax.set_ylim(0,1)
ax.set_xlabel("Error Rate")
ax.set_ylabel("ARI")
ax.legend(fontsize=12, handletextpad=0.5, loc="lower left",
              borderpad=0.3,
              columnspacing=1.3,
              handlelength=0.65,)
ax.yaxis.grid(True, color='#e0e0e0')
ax.xaxis.grid(True, color='#e0e0e0')
ax.spines['top'].set_edgecolor('#e0e0e0')  # 设置上边坐标轴颜色
ax.spines['bottom'].set_edgecolor('#e0e0e0')  # 设置下边坐标轴颜色
ax.spines['left'].set_edgecolor('#e0e0e0')  # 设置左边坐标轴颜色
ax.spines['right'].set_edgecolor('#e0e0e0')
ax.tick_params(axis='both', color='#e0e0e0')
plt.savefig("mistake-sp.svg")
