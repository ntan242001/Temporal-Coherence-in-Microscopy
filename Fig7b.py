import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.signal import savgol_filter

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


smooth_R_G1_nac = savgol_filter(R_G1_nac, 61, 9)
smooth_R_FN_nac = savgol_filter(R_FN_nac, 61, 9)
smooth_R_G1_ac = savgol_filter(R_G1_ac, 25, 8)
smooth_R_FN_ac = savgol_filter(R_FN_ac, 25, 8)

########## Plotting the curves ############
fig, ax = plt.subplots(2,1,sharex=True)
fig.set_size_inches(6.5,7)
ax1,ax2=ax
fig.tight_layout(rect=[0.03,0.08,1,1])
fig.subplots_adjust(hspace=0)
width = 2

for i in range(2):
    ax[i].set_ylabel("Resolution (nm)", fontsize=12)
    ax[i].set_xlim(-1.8,1.8)
    ax[i].minorticks_on()
    
ax1.plot(delta_z_seriesGauss/((C_3*lamda)**(1/2)), smooth_R_G1_nac, 'b--', linewidth=width, label = "NAC, Gaussian")
ax1.plot(delta_z_seriesFN/((C_3*lamda)**(1/2)), smooth_R_FN_nac, 'b-', linewidth=width, label = "NAC, FN")
ax2.plot(delta_z_series_ac/((C_5*lamda**2)**(1/3)), smooth_R_G1_ac, 'r--', linewidth=width, label = "AC, Gaussian")
ax2.plot(delta_z_series_ac/((C_5*lamda**2)**(1/3)), smooth_R_FN_ac, 'r-', linewidth=width, label = "AC, FN")

# ax1.scatter(delta_z_seriesGauss/((C_3*lamda)**(1/2)), R_G1_nac, color='b', linewidth=width, label = "NAC, Gaussian")
# ax1.scatter(delta_z_seriesFN/((C_3*lamda)**(1/2)), R_FN_nac, color='g', linewidth=width, label = "NAC, FN")
# ax2.scatter(delta_z_series_ac/((C_5*lamda**2)**(1/3)), R_G1_ac, color='r', linewidth=width, label = "AC, Gaussian")
# ax2.scatter(delta_z_series_ac/((C_5*lamda**2)**(1/3)), R_FN_ac, color='y', linewidth=width, label = "AC, FN")

# ax1.set_ylim([])

yticks1=[2.0,3.0,4.0,5.0, 6.0]
ax1.set_ylim([1.64,6])
tickLabels = map(str, yticks1)
ax1.set_yticks(yticks1)
ax1.set_yticklabels(tickLabels)
yticks2 = 0.2*np.arange(2,9)
ax2.set_yticks(yticks2)

ypos=0.1
ax2.text(-0.57,ypos,'$\Delta$z', color='k', fontsize=13)
ax2.text(-0.39,ypos,'[', color='k', fontsize=15)
ax2.text(-0.35,ypos,'$(C_3 \lambda)^{1/2} $', color='b', fontsize=12)
ax2.text(0.05,ypos,', ', color='k', fontsize=12)
ax2.text(0.12,ypos,'$(C_5 \lambda^2)^{1/3} $', color='r', fontsize=12)
ax2.text(0.59,ypos,']', color='k', fontsize=15)

ax1.text(1.5,1.8,'NAC', color='b', fontsize=12)
ax2.text(1.55,0.5,'AC', color='r', fontsize=12)

Scherzer3 = np.sqrt(C_3*lamda/2) 
Scherzer5 = (9/64*C_5*lamda**2)**(1/3)

ax1.axvline(x=Scherzer3/((C_3*lamda)**(1/2)), linestyle = '--', c = 'b')
ax2.axvline(x=Scherzer5/((C_5*lamda**2)**(1/3)), linestyle = '--', c = 'r')

# ax.plot(delta_z_seriesGauss*1e6, smooth_R_G1_nac, 'm-', linewidth=width, label = "NAC, Gaussian")
# ax.plot(delta_z_seriesFN*1e6, smooth_R_FN_nac, 'c-', linewidth=width, label = "NAC, FN")



# ax.set_ylim(0, 7)
# ax.set_xlim(-1.8,1.8)


# ax.set_xlabel(r"$\Delta $z $(\mathrm{\mu m})$")
# ax.set_xlabel('$\Delta z ((C_3 \lambda)^{1/2}, (C_5 \lambda^2)^{1/3})$')
# ax.set_xticks([-1,0,1])



# ax.axvline(x=6.90e-7/((C_3*lamda)**(1/2)), linestyle = '--', c = 'm', linewidth=2)
# ax.axvline(x=1.21e-6/((C_3*lamda)**(1/2)), linestyle = '--', c = 'c')
# ax.axvline(x=4.73e-8/((C_5*lamda**2)**(1/3)), linestyle = '--', c = 'r')
# ax.axvline(x=7.1e-8/((C_5*lamda**2)**(1/3)), linestyle = '--', c = 'b')

