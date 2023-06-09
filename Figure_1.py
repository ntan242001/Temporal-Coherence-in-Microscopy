# Comparing the single and double Gaussian distributions for different aperture angles
import numpy as np
import matplotlib.pyplot as plt
#from joblib import Parallel, delayed
import time
import sys
sys.path.insert(1, 'Users/alantan/Library/Mobile Documents/com~apple~CloudDocs/Fourier Optics/1.PaperFig/')
import SizeFormat

plt.rcParams['font.family'] = 'Arial'

# A function to choose different LEEM parameters
def choose_LEEM_type(LEEM_type_str, aberration_corrected_bool = False):
    global E, E_0, C_c, C_cc, C_3c, C_3, C_5, alpha_ap, alpha_ill, \
        delta_E, M_L, lamda, lamda_0, q_ill, q_ap, LEEM_type, aberration_corrected
    LEEM_type = LEEM_type_str
    aberration_corrected = aberration_corrected_bool
    
    if LEEM_type == "IBM":
        if aberration_corrected == False:
            E = 15010  # eV  Nominal Energy After Acceleration
            E_0 = 10  # eV  Energy at the sample
            
            C_c = -0.075  # m  Second Rank Chromatic Aberration Coefficient
            C_cc = 23.09 # m   Third Rank Chromatic Aberration Coefficient
            C_3c = -59.37  # m   Forth Rank Chromatic Aberration Coefficient
            
            C_3 = 0.345  # m  Third Order Spherical Aberration Coefficient
            C_5 = 39.4  # m  Fifth Order Spherical Aberration Coefficient
            
            alpha_ap = 2.34e-3  # rad Aperture angle
            alpha_ill = 0.055e-3  # rad Illumination Divergence Angle
            
            delta_E = 0.25  # eV  Energy Spread
            M_L = 0.653  # Lateral Magnification
            print("IBM nac chosen.")
            
        elif aberration_corrected == True:
            E = 15010  # eV  Nominal Energy After Acceleration
            E_0 = 10  # eV  Energy at the sample
            
            C_c = 0  # m   Second Rank Chromatic Aberration Coefficient
            C_cc = 27.9 # m   Third Rank Chromatic Aberration Coefficient
            C_3c = -67.4 # m   Forth Rank Chromatic Aberration Coefficient
            
            C_3 = 0  # m   Spherical Aberration Coefficient
            C_5 = 92.8
        
            alpha_ap = 7.37e-3  # rad Aperture angle
            alpha_ill = 0.055e-3  # rad Illumination Divergence Angle
        
            delta_E = 0.25  # eV  Energy Spread
            M_L = 0.653  # Lateral Magnification
            print("IBM ac chosen.")
            
        lamda = 6.6261e-34 / np.sqrt(2 * 1.6022e-19 * 9.1095e-31 * E)
        lamda_0 = 6.6261e-34 / np.sqrt(2 * 1.6022e-19 * 9.1095e-31 * E_0)
    
        q_ap = alpha_ap/lamda
        q_ill = alpha_ill/lamda
        
    elif LEEM_type == "Energy dependent":
        if aberration_corrected == False:
            E = 15010  # eV  Nominal Energy After Acceleration
            E_0 = 10 # eV  Energy at the sample ##########CUSTOMIZABLE INPUT##########
            kappa = np.sqrt(E/E_0)
            
            C_c = -0.0121 * kappa**(1/2) + 0.0029 # m  Second Rank Chromatic Aberration Coefficient
            C_cc = 0.5918 * kappa**(3/2) - 87.063 # m   Third Rank Chromatic Aberration Coefficient
            C_3c = -1.2141 * kappa**(3/2) + 169.41  # m   Forth Rank Chromatic Aberration Coefficient
            
            C_3 = 0.0297 * kappa**(1/2) + 0.1626  # m  Third Order Spherical Aberration Coefficient
            C_5 = 0.6223 * kappa**(3/2) - 79.305  # m  Fifth Order Spherical Aberration Coefficient
            
            delta_E = 0.25  # eV  Energy Spread
            alpha_ill = 0.055e-3  # rad Illumination divergence angle
            M_L = 0.653  # Lateral Magnification
            
            lamda = 6.6261e-34 / np.sqrt(2 * 1.6022e-19 * 9.1095e-31 * E) # in metre
            alpha_ap = (lamda/C_3)**(1/4) # rad Aperture angle for optimal resolution
            
            print("Custom nac LEEM at E_0 = " + str(E_0) + " eV chosen.")
            
        if aberration_corrected == True:
            E = 15010  # eV  Nominal Energy After Acceleration
            E_0 = 10 # eV  Energy at the sample
            kappa = np.sqrt(E/E_0)
            
            C_c = 0 # m  Second Rank Chromatic Aberration Coefficient
            C_cc = 0.5984 * kappa**(3/2) - 84.002 # m   Third Rank Chromatic Aberration Coefficient 
            C_3c = -1.1652 * kappa**(3/2) + 153.58  # m   Forth Rank Chromatic Aberration Coefficient  
            
            C_3 = 0  # m  Third Order Spherical Aberration Coefficient
            C_5 = 0.5624 * kappa**(3/2) - 16.541  # m  Fifth Order Spherical Aberration Coefficient
            
            delta_E = 0.25  # eV  Energy Spread
            alpha_ill = 0.055e-3  # rad Illumination divergence angle
            M_L = 0.653  # Lateral Magnification
            
            lamda = 6.6261e-34 / np.sqrt(2 * 1.6022e-19 * 9.1095e-31 * E) # in metre
            alpha_ap = (3/2*lamda/C_5)**(1/6) # rad Aperture angle for optimal resolution
            
            print("Custom ac LEEM at E_0 = " + str(E_0) + " eV chosen.")   
        
        lamda_0 = 6.6261e-34 / np.sqrt(2 * 1.6022e-19 * 9.1095e-31 * E_0) # in metre
        
        q_ap = alpha_ap/lamda
        q_ill = alpha_ill/lamda

object_size = 400            # simulating object size in nm
simulating_steps = 1 + 2**11 # total simulating steps
# An array of points in the x space
x_array = (np.linspace(-object_size/2, object_size/2, simulating_steps) + object_size/simulating_steps)*1e-9

choose_LEEM_type("IBM", aberration_corrected_bool = True)



t_0 = time.time()
print("Simulation start.")
# # The object image is reversed through the lens
# object_function_reversed = object_function[::-1] 

#epsilon_0_series = np.linspace(0.0008694, 0.0008694, 1)
epsilon_0_series = np.array([-0.2, 0, 0.2])

q = 1 / (simulating_steps* (x_array[1] - x_array[0])) * np.arange(0, simulating_steps, 1)

# Shifting the q array to centre at 0 
q = q - np.max(q) / 2
# Arrays for the calculation of the double integration 
Q, QQ = np.meshgrid(q, q)

# A function to calculate the image for single Gaussian distribution
def findEctot(LEEMType, epsilon_n):
    
    if LEEMType == 'NAC':
        C_c = -0.075  # m  Second Rank Chromatic Aberration Coefficient
        C_cc = 23.09 # m   Third Rank Chromatic Aberration Coefficient
        C_3c = -59.37  # m   Forth Rank Chromatic Aberration Coefficient
    if LEEMType == 'AC':
        C_c = 0  # m   Second Rank Chromatic Aberration Coefficient
        C_cc = 27.9 # m   Third Rank Chromatic Aberration Coefficient
        C_3c = -67.4 # m   Forth Rank Chromatic Aberration Coefficient
    
    E = 15010 # eV
    delta_E = 0.2424 # eV
    sigma_E = delta_E/(2*np.sqrt(2*np.log(2)))
    b_1 = 1/2*C_c*lamda*(Q**2 - QQ**2)/E + 1/4*C_3c*lamda**3*(Q**4 - QQ**4)/E
    b_2 = 1/2*C_cc*lamda*(Q**2 - QQ**2)/E**2
    
    d_n = b_1 - 1j*epsilon_n/(2*np.pi*sigma_E**2) 
    
    # The purely chromatic envelop functions
    E_cc = (1 - 1j*4*np.pi*b_2*sigma_E**2)**(-1/2)
    E_ct = E_cc * np.exp(-2*(np.pi*sigma_E*E_cc*d_n)**2 - epsilon_n**2/(2*sigma_E**2))
    
    E_ctot = E_ct[int(len(q)/2),:]
    
    return E_ctot

q_apNAC = 2.34e-3/lamda 
q_apAC = 7.37e-3/lamda 
xticks = [0, 0.5, 1, 1.5]
yticks = [-1,-0.5, 0, 0.5, 1]

size = SizeFormat.size + 3
framesize = SizeFormat.frame_size
framewidth = SizeFormat.frame_width
tick_font = SizeFormat.tick_font
width = SizeFormat.width
minor_tick_size = SizeFormat.minor_tick_size
major_tick_size = SizeFormat.major_tick_size
tick_width = SizeFormat.tick_width
vline_width = SizeFormat.vline_width

fig, ax = plt.subplots(nrows=1, ncols=3)
fig.set_size_inches(14, 3.7)
fig.tight_layout(rect=[0.04, 0.0, 1.01, 1.01])
fig.subplots_adjust(wspace=0.3)

ax1,ax2,ax3 = ax
# e=-0.2
tex_y = 1.03
tex_y_epsilon = 1.015
sub_label = ['(a)', '(b)', '(c)']

for i in range(3):
    ax[i].set_ylim(-1, 1.2)
    ax[i].set_xlim(0, 1.6)
    ax[i].axhline(y=0, color='k', linestyle='-', linewidth=vline_width)
    ax[i].set_xlabel(r'$q$ $\mathrm{(nm^{-1})}$', fontsize=size)
    ax[i].axvline(x=q_apNAC/1e9, color='r', linestyle=':',linewidth=2)
    ax[i].axvline(x=q_apAC/1e9, color='b', linestyle=':',linewidth=2)
    ax[i].minorticks_on()
    ax[i].set_ylabel('Amplitude', fontsize=size)
    ax[i].set_yticks(yticks)
    ax[i].set_xticks(xticks)
    ax[i].set_xticklabels(xticks, rotation=0, fontsize=tick_font+2)
    ax[i].set_yticklabels(yticks, rotation=0, fontsize=tick_font+2)
    E_ctotNAC = findEctot('NAC', epsilon_0_series[i])
    E_ctotAC = findEctot('AC', epsilon_0_series[i])
    ax[i].plot(q/(1e9), E_ctotAC.real, 'b-', linewidth=2, label = r'Re[$E_{C,n,\rm AC}(q,0)$]')
    ax[i].plot(q/(1e9), E_ctotAC.imag, 'b--', linewidth=2, label = r'Im[$E_{C,n,\rm AC}(q,0)$]')
    ax[i].plot(q/(1e9), E_ctotNAC.real, 'r-', linewidth=2, label = r'Re[$E_{C,n,\rm NAC}(q,0)$]')
    ax[i].plot(q/(1e9), E_ctotNAC.imag, 'r--', linewidth=2, label = r'Im[$E_{C,n,\rm NAC}(q,0)$]')
    ax[i].text(q_apNAC/1e9 + 0.02,tex_y,'NAC', c='r', fontsize=size-2)
    ax[i].text(q_apAC/1e9 + 0.02,tex_y,'AC', c='b', fontsize=size-2)
    ax[i].text(-0.34,1.1,sub_label[i],fontsize=size)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax[i].spines[axis].set_linewidth(framewidth)

    ax[i].tick_params(which='major', length=major_tick_size, width=tick_width)
    ax[i].tick_params(which='minor', length=minor_tick_size, width=tick_width)

ax1.text(0.99, tex_y_epsilon, r'$\epsilon_{\rm n} = $'+str(epsilon_0_series[0])+' eV',fontsize=size-3)
ax2.text(1.02, tex_y_epsilon, r'$\epsilon_{\rm n} = $'+str(epsilon_0_series[1])+' eV',fontsize=size-3)
ax3.text(1.02, tex_y_epsilon, r'$\epsilon_{\rm n} = $'+str(epsilon_0_series[2])+' eV',fontsize=size-3)
# ax2.text(q_apNAC/1e9 + 0.02,tex_y,'NAC', c='r', fontsize=s)
# ax2.text(q_apAC/1e9 + 0.02,tex_y,'AC', c='b', fontsize=s)

# ax1.text(x=1.35, y=1.08, s = r'$\epsilon_n = -0.2 eV$', color = 'k')

handles, labels = ax3.get_legend_handles_labels()
# fig.legend(handles, labels, loc=8, ncol = 4, frameon=False, fontsize=s)
fig.subplots_adjust(bottom=0.2)


plt.savefig('Figure_1.png', dpi=1000)



