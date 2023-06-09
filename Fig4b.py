import numpy as np
import matplotlib.pyplot as plt
import time
import sys
sys.path.insert(0, '/Users/alantan/Library/Mobile Documents/com~apple~CloudDocs/Fourier Optics/1.PaperFig/SizeFormat.py')
import SizeFormat

plt.rcParams['font.family'] = 'Arial'

t0=time.time()

E = 15010  # eV  Nominal Energy After Acceleration
E_0 = 10  # eV  Energy at the sample
alpha_ill = 0.055e-3  # rad Illumination Divergence Angle
lamda = 6.6261e-34 / np.sqrt(2 * 1.6022e-19 * 9.1095e-31 * E)
lamda_0 = 6.6261e-34 / np.sqrt(2 * 1.6022e-19 * 9.1095e-31 * E_0)
q_ill = alpha_ill/lamda
sigma_ill = q_ill/(2*np.sqrt(2*np.log(2)))
Scherzer3 = np.sqrt(0.345*lamda/2) 
Scherzer5 = (9/64*92.8*lamda**2)**(1/3)

object_size = 400            # simulating object size in nm
simulating_steps = 1 + 2**11 # total simulating steps
x_array = (np.linspace(-object_size/2, object_size/2, simulating_steps) + object_size/simulating_steps)*1e-9

# A function to choose different LEEM parameters
def choose_LEEM_type(LEEM_type):
    global C_c, C_cc, C_3c, C_3, C_5
    if LEEM_type == 'NAC':
        C_c = -0.075  # m  Second Rank Chromatic Aberration Coefficient
        C_cc = 23.09 # m   Third Rank Chromatic Aberration Coefficient
        C_3c = -59.37  # m   Forth Rank Chromatic Aberration Coefficient
        C_3 = 0.345  # m  Third Order Spherical Aberration Coefficient
        C_5 = 39.4  # m  Fifth Order Spherical Aberration Coefficient
        
    elif LEEM_type == 'AC':
        C_c = 0  # m   Second Rank Chromatic Aberration Coefficient
        C_cc = 27.9 # m   Third Rank Chromatic Aberration Coefficient
        C_3c = -67.4 # m   Forth Rank Chromatic Aberration Coefficient
        C_3 = 0  # m   Spherical Aberration Coefficient
        C_5 = 92.8


q = 1 / (simulating_steps* (x_array[1] - x_array[0])) * np.arange(0, simulating_steps, 1)

# Shifting the q array to centre at 0 
q = q - np.max(q) / 2
# Arrays for the calculation of the double integration 
Q, QQ = np.meshgrid(q, q)

def Ws(delta_z, LEEMType):
    choose_LEEM_type(LEEMType)
    Ws = np.exp(1j*2*np.pi*(C_3*lamda**3 * (Q**4 - QQ**4)/4 + C_5*lamda**5 *(
        Q**6 - QQ**6)/6 - delta_z*lamda*(Q**2 - QQ**2)/2))
    return Ws[int(len(q)/2),:]

def Es(delta_z):
    a_1 = C_3*lamda**3 *(Q**3 - QQ**3) + C_5*lamda**5 * (Q**5 - QQ**5) - delta_z*lamda*(Q - QQ)
    E_s = np.exp(-2*np.pi**2 *sigma_ill**2 *a_1**2)
    return E_s[int(len(q)/2),:]

# A function to calculate the image for single Gaussian distribution
def Ectot(epsilon_n, sigma_E):
    b_1 = 1/2*C_c*lamda*(Q**2 - QQ**2)/E + 1/4*C_3c*lamda**3*(Q**4 - QQ**4)/E
    b_2 = 1/2*C_cc*lamda*(Q**2 - QQ**2)/E**2
    
    d_n = b_1 - 1j*epsilon_n/(2*np.pi*sigma_E**2) 
    
    # The purely chromatic envelop functions
    E_cc = (1 - 1j*4*np.pi*b_2*sigma_E**2)**(-1/2)
    E_ct = E_cc * np.exp(-2*(np.pi*sigma_E*E_cc*d_n)**2 - epsilon_n**2/(2*sigma_E**2))
    
    E_ctot = E_ct[int(len(q)/2),:]
    return E_ctot

def EcFN():
    return 0.28682*Ectot(-0.03925, 0.0531)+ 0.33146*Ectot(-0.3382, 0.1991) + 0.38173*Ectot(-0.1438, 0.0962)

def R_G1(delta_z, LEEMType):
    choose_LEEM_type(LEEMType)
    delta_E = 0.2424 
    sigma_E = delta_E/(2*np.sqrt(2*np.log(2)))
    return Ws(delta_z, LEEMType)*Es(delta_z)*Ectot(0,sigma_E)

def R_FE(delta_z, LEEMType):
    choose_LEEM_type(LEEMType)
    return Ws(delta_z, LEEMType)*Es(delta_z)*EcFN()


q=q/1e9
q_apNAC_FE = 1.93e-3/lamda 
q_apNAC_G = 2.34e-3/lamda 
q_apAC_FE = 7.04e-3/lamda 
q_apAC_G = 7.37e-3/lamda 
# xticks = [0, 0.5, 1, 1.5]
yticks = 0.2*np.arange(-1,7)
alpha_ap_list = [2.51e-3, 2.34e-3]
# defocus_list = [0.75e-6,0]

CTF = R_FE(0, 'AC')
CTFG = R_G1(0, 'AC')
W_s = Ws(0, 'AC')

fig, ax = plt.subplots(2,1)
ax1,ax2 = ax
fig.set_size_inches(4.8, 5.2)
fig.tight_layout(rect=[0.09, 0.08, 1, 1])
fig.subplots_adjust(hspace=0)


###############################################################################
size = SizeFormat.size+1
framesize = SizeFormat.frame_size
framewidth = SizeFormat.frame_width
tick_font = SizeFormat.tick_font+1
width = SizeFormat.width
minor_tick_size = SizeFormat.minor_tick_size
major_tick_size = SizeFormat.major_tick_size
tick_width = SizeFormat.tick_width
vline_width = SizeFormat.vline_width

ax1.plot(q, W_s.real, 'k:', linewidth=width)
ax1.plot(q, CTFG.real, 'b--', linewidth=width)
ax1.plot(q, CTF.real, 'r-', linewidth=width)

ax2.plot(q, W_s.imag, 'k:', linewidth=width, label = r'$\mathrm{W_s}$')
ax2.plot(q, CTFG.imag, 'b--', linewidth=width, label = r'$\mathrm{R_G}$')
ax2.plot(q, CTF.imag, 'r-', linewidth=width, label = r'$\mathrm{R_{FE}}$')


ax2.axhline(y=0, color='k', linestyle='-', linewidth=vline_width)

ax2.set_xlabel(r'q $\mathrm{(nm^{-1})}$', fontsize=size)

yticks= np.array([-1. , -0.5,  0. ,  0.5,  1. ])

for i in range(2):
    ax[i].minorticks_on()
    ax[i].set_yticks(yticks)
    ax[i].set_yticklabels(yticks, fontsize=size-2)
    ax[i].axvline(x=q_apAC_FE/1e9, color='r', linestyle=':', linewidth=2)
    ax[i].axvline(x=q_apAC_G/1e9, color='b', linestyle='-.', linewidth=2)
    ax[i].set_ylabel('Amplitude', fontsize=size)
    ax[i].axhline(y=0, color='k', linestyle='-', linewidth=vline_width)
    ax[i].set_ylim(-1.2, 1.2)
    ax[i].set_xlim(0, 1.3)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax[i].spines[axis].set_linewidth(framewidth)

    ax[i].tick_params(which='major', length=major_tick_size, width=tick_width)
    ax[i].tick_params(which='minor', length=minor_tick_size, width=tick_width)

xticks=np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
# xticks=np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6])
ax1.set_xticks([])

ax2.set_xticks(xticks)
ax2.set_xticklabels(xticks, fontsize=size-2)

ax1.text(-.3, 1.1, '(b)',fontsize=size)

x0=0.056
ax1.text(x0, 0.4, 'AC',fontsize=size)
ax1.text(x0, -1.1, 'Real',fontsize=size)
ax2.text(x0, 0.9, 'Imaginary',fontsize=size)

ax2.legend(frameon=False, fontsize=size-2, bbox_to_anchor=(0.35,0.5), loc=0)

print(time.time() - t0)

plt.savefig("Figure_4b.png", dpi=1000)