import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.signal import savgol_filter
import sys
sys.path.insert(0, '/Users/alantan/Library/Mobile Documents/com~apple~CloudDocs/Fourier Optics/1.PaperFig/SizeFormat.py')
import SizeFormat

plt.rcParams['font.family'] = 'Arial'


C_3 = 0.345 
C_5 = 92.8

E = 15010 
lamda = 6.6261e-34 / np.sqrt(2 * 1.6022e-19 * 9.1095e-31 * E)

delta_z_series_nac = []
R_G1_nac = []
R_FN_nac = []

delta_z_series_ac = []
R_G1_ac = []
R_FN_ac = []

with open('Run0_R(dz)_G1_nac.csv', 'r') as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            if row == []:
                continue
            else:
                delta_z_series_nac.append(float(row[0]))
                R_G1_nac.append(float(row[1]))    
    csvfile.close()

with open('Run0_R(dz)_FN_nac.csv', 'r') as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            if row == []:
                continue
            else:
                R_FN_nac.append(float(row[1]))    
    csvfile.close()
    
with open('Run0plot_R(dz)_G1_ac.csv', 'r') as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            if row == []:
                continue
            else:
                delta_z_series_ac.append(float(row[0]))
                R_G1_ac.append(float(row[1]))    
    csvfile.close()

with open('Run0plot_R(dz)_FN_ac.csv', 'r') as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            if row == []:
                continue
            else:
                R_FN_ac.append(float(row[1]))    
    csvfile.close()
    


delta_z_seriesGauss = np.array(delta_z_series_nac)
delta_z_seriesFN = np.array(delta_z_series_nac)

Scherzer3 = np.sqrt(C_3*lamda/2) 
Scherzer5 = (9/64*C_5*lamda**2)**(1/3)


smooth_R_G1_nac = savgol_filter(R_G1_nac, 61, 9)
smooth_R_FN_nac = savgol_filter(R_FN_nac, 61, 9)
smooth_R_G1_ac = savgol_filter(R_G1_ac, 25, 8)
smooth_R_FN_ac = savgol_filter(R_FN_ac, 25, 8)

########## Plotting the curves ############
size = SizeFormat.size+1
framesize = SizeFormat.frame_size
framewidth = SizeFormat.frame_width
tick_font = SizeFormat.tick_font
width = SizeFormat.width
minor_tick_size = SizeFormat.minor_tick_size
major_tick_size = SizeFormat.major_tick_size
tick_width = SizeFormat.tick_width
vline_width = SizeFormat.vline_width

fig, ax = plt.subplots(2,1)
fig.set_size_inches(5,5.5)
ax1,ax2=ax
fig.tight_layout(rect=[0.09,0.08,1.01,1.01])
fig.subplots_adjust(hspace=0.05)

for i in range(2):
    ax[i].set_xlim(-1.8,1.8)
    ax[i].minorticks_on()
    ax[i].tick_params(which='major', length=major_tick_size, width=tick_width)
    ax[i].tick_params(which='minor', length=minor_tick_size, width=tick_width)
    
    
y1min, y1max=1.2, 5.8
ax1.axhline(y=y1min+0.2, linestyle='-', c='k', linewidth=0.9)
ax1.axvline(x=Scherzer3/((C_3*lamda)**(1/2)), ymin=0.2/(y1max-y1min), ymax=1, linestyle = '--', c = 'b')

ax2.axvline(x=Scherzer5/((C_5*lamda**2)**(1/3)), ymin=0, ymax=1.09, linestyle = '--', c = 'r', clip_on=False)

ax2.set_ylabel("Resolution (nm)", fontsize=size+1)
ax2.yaxis.set_label_coords(-0.15,1)
ax1.spines["bottom"].set_visible(False)
ax2.spines["top"].set_visible(False)
ax1.set_xticks([])
ax1.plot(delta_z_seriesGauss/((C_3*lamda)**(1/2)), smooth_R_G1_nac, 'b--', linewidth=width, label = "NAC, Gaussian")
ax1.plot(delta_z_seriesFN/((C_3*lamda)**(1/2)), smooth_R_FN_nac, 'b-', linewidth=width, label = "NAC, FN")
ax2.plot(delta_z_series_ac/((C_5*lamda**2)**(1/3)), smooth_R_G1_ac, 'r--', linewidth=width, label = "AC, Gaussian")
ax2.plot(delta_z_series_ac/((C_5*lamda**2)**(1/3)), smooth_R_FN_ac, 'r-', linewidth=width, label = "AC, FN")

# ax1.scatter(delta_z_seriesGauss/((C_3*lamda)**(1/2)), R_G1_nac, color='b', linewidth=width, label = "NAC, Gaussian")
# ax1.scatter(delta_z_seriesFN/((C_3*lamda)**(1/2)), R_FN_nac, color='g', linewidth=width, label = "NAC, FN")
# ax2.scatter(delta_z_series_ac/((C_5*lamda**2)**(1/3)), R_G1_ac, color='r', linewidth=width, label = "AC, Gaussian")
# ax2.scatter(delta_z_series_ac/((C_5*lamda**2)**(1/3)), R_FN_ac, color='y', linewidth=width, label = "AC, FN")

ax1.set_ylim([y1min, y1max])
ax2.set_ylim([0.4,1.9])

yticks1=[2.0,3.0,4.0,5.0]
ax1.set_yticks(yticks1)
ax1.set_yticklabels(yticks1, fontsize=tick_font)
yticks2 = np.array([0.4, 0.6, 0.8, 1. , 1.2, 1.4, 1.6, 1.8])
# yticks2=[.5,1.0,1.5]
ax2.set_yticks(yticks2)
ax2.set_yticklabels(yticks2, fontsize=tick_font)
ax2.tick_params(axis='y', which='minor', left=False)

xticks = np.array([-1.5, -1. , -0.5,  0. ,  0.5,  1. ,  1.5])
ax2.set_xticks(xticks)
ax2.set_xticklabels(xticks, fontsize=tick_font)

ypos=-0.05
ax2.text(-1.1,ypos,'$\Delta z$', color='k', fontsize=size)
ax2.text(-0.77,ypos,'[', color='k', fontsize=size+3)
ax2.text(-0.69,ypos,'$(C_3 \lambda)^{1/2} $', color='b', fontsize=size)
ax2.text(0.09,ypos,', ', color='k', fontsize=size)
ax2.text(0.15,ypos,'$(C_5 \lambda^2)^{1/3} $', color='r', fontsize=size)
ax2.text(1.06,ypos,']', color='k', fontsize=size+3)

ax1.text(1.32, 1.55,'NAC', color='b', fontsize=size)
ax2.text(1.47,0.45,'AC', color='r', fontsize=size)

d = .01  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

# ax1.text(-2.28, 5.75, '(a)',fontsize=s)

for axis in ['top', 'bottom', 'left', 'right']:
    ax1.spines[axis].set_linewidth(framewidth)
    ax2.spines[axis].set_linewidth(framewidth)

plt.savefig("Figure_7a.png", dpi=1000)
