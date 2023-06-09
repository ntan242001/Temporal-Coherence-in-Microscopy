import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.signal import savgol_filter
import sys
sys.path.insert(0, '/Users/alantan/Library/Mobile Documents/com~apple~CloudDocs/Fourier Optics/1.PaperFig/SizeFormat.py')
import SizeFormat

plt.rcParams['font.family'] = 'Arial'

a_array_nacG = []
a_array_nacFN = []
R_G1_nac = []
R_FN_nac = []

a_array_acG = []
a_array_acFN = []
R_G1_ac = []
R_FN_ac = []


with open('Run1plot_R(a)_G1_nac.csv', 'r') as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            if row == []:
                continue
            else:
                a_array_nacG.append(float(row[0]))
                R_G1_nac.append(float(row[1]))   
    csvfile.close()

with open('Run1plot_R(a)_FN_nac.csv', 'r') as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            if row == []:
                continue
            else:
                a_array_nacFN.append(float(row[0]))
                R_FN_nac.append(float(row[1]))   
    csvfile.close()
a_array_nacG = np.array(a_array_nacG)*1e3
a_array_nacFN = np.array(a_array_nacFN)*1e3
    
with open('Run1plot_R(a)_G1_ac.csv', 'r') as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            if row == []:
                continue
            else:
                a_array_acG.append(float(row[0]))
                R_G1_ac.append(float(row[1]))   
    csvfile.close()
    

with open('Run1plot_R(a)_FN_ac.csv', 'r') as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            if row == []:
                continue
            else:
                a_array_acFN.append(float(row[0]))
                R_FN_ac.append(float(row[1]))   
    csvfile.close()
    
a_array_acG = np.array(a_array_acG)*1e3
a_array_acFN = np.array(a_array_acFN)*1e3
    
smooth_R_G1_nac = savgol_filter(R_G1_nac, 31, 9)
smooth_R_FN_nac = savgol_filter(R_FN_nac, 31, 9)
smooth_R_G1_ac = savgol_filter(R_G1_ac, 31, 9)
smooth_R_FN_ac = savgol_filter(R_FN_ac, 31, 9)
    
########## Plotting the curves ############
yminimum, ymaximum = 0.4, 10
size = SizeFormat.size
framesize = SizeFormat.frame_size
framewidth = SizeFormat.frame_width
tick_font = SizeFormat.tick_font
width = SizeFormat.width
minor_tick_size = SizeFormat.minor_tick_size
major_tick_size = SizeFormat.major_tick_size
tick_width = SizeFormat.tick_width

fig, ax = plt.subplots()
fig.tight_layout(rect=[0.09,0.09,1.02,1])
framesize=(5.38,4)
fig.set_size_inches(framesize)


ax.loglog(a_array_nacFN, smooth_R_FN_nac, 'b-', label = "NAC, FE", linewidth=width)
ax.loglog(a_array_nacG, smooth_R_G1_nac, 'b--', label = "NAC, Gaussian", linewidth=width)
ax.loglog(a_array_acFN, smooth_R_FN_ac, 'r-', label = "AC, FE", linewidth=width)
ax.loglog(a_array_acG, smooth_R_G1_ac, 'r--', label = "AC, Gaussian", linewidth=width)


import matplotlib.ticker as ticker
def myLogFormat(y,pos):
    # Find the number of decimal places required
    decimalplaces = int(np.maximum(-np.log10(y),0))   
    # Insert that number into a format string
    formatstring = '{{:.{:1d}f}}'.format(decimalplaces)
    # Return the formatted tick label
    return formatstring.format(y)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(myLogFormat))
ax.xaxis.set_major_formatter(ticker.FuncFormatter(myLogFormat))

yNAC = 2.9
yAC = 0.67
textx = 14.5
ax.text(textx-1.5, yNAC,'NAC', c='b', fontsize=size)
ax.text(textx, yAC,'AC', c='r', fontsize=size)

ax.set_xlim(0.5, 20)
ax.text(0.257,ymaximum-0.7, '(b)',fontsize=size)

xticks=[1,10]
ax.set_xticks(xticks)
ax.set_xticklabels(xticks,fontsize=tick_font)
ax.set_yticks(xticks)
ax.set_yticklabels(xticks,fontsize=tick_font)

for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(framewidth)

ax.tick_params(which='major', length=major_tick_size, width=tick_width)
ax.tick_params(which='minor', length=minor_tick_size, width=tick_width)

ax.set_xlabel("Aperture angle (mrad)", fontsize=size+1, labelpad=5)
ax.set_ylabel("Resolution (nm)", fontsize=size+1, labelpad=5)
ax.set_ylim(yminimum, ymaximum)
# ax.legend(frameon=False)
plt.savefig("Figure_6b.png", dpi=1000)