# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                                                 #
# UnitCellPostProcessor.py                                                                        #
# Postprocessing tool                                                                             #
# Units: Kilograms-Meters-Seconds                                                                 #
# Initiated by Dennis Hagnas on Tue Jan 14 2025                                                   #
#                                                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import csv
from scipy.interpolate import interp1d
from numpy.linalg import eig

plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'Times',
    'font.size': 8,
    'axes.labelsize': 8,
    'axes.formatter.limits': (-2,2)
})
plt.rcParams['lines.linewidth'] = 0.8

#--------------------------------------------------------------------------------------------------
# Import parametric variables

with open('parameters.txt', 'r') as f:
    exec(f.read())
#--------------------------------------------------------------------------------------------------
# Define strain points to interpolate and extrapolate results to

strainPoints1 = np.array([
    -1.05e-02,-1.00e-02,-9.78e-03,-8.79e-03,-7.89e-03,-7.09e-03, 
    -6.36e-03,-5.72e-03,-5.14e-03,-4.61e-03,-4.14e-03,-3.72e-03,
    -3.34e-03,-3.00e-03,-2.69e-03,-2.41e-03,-2.17e-03,-1.95e-03,
    -1.75e-03,-1.57e-03,-1.41e-03,-1.26e-03,-1.13e-03,-1.02e-03,
    -9.12e-04,-8.18e-04,-7.33e-04,-6.57e-04,-5.89e-04,-5.28e-04,
    -4.72e-04,-4.23e-04,-3.78e-04,-3.38e-04,-3.02e-04,-2.69e-04,
    -2.41e-04,-2.14e-04,-1.91e-04,-1.70e-04,-1.51e-04,-1.34e-04,
    -1.04e-04,-9.25e-05,-8.14e-05,-7.14e-05,-6.25e-05,-5.44e-05,
    -4.72e-05,-4.07e-05,-3.49e-05,-2.96e-05,-2.49e-05,-2.07e-05,
    -1.69e-05,-1.35e-05,-1.04e-05,-7.70e-06,-5.22e-06,-3.00e-06,
    -1.00e-06,0.00e+00,1.00e-06,3.00e-06,5.22e-06,7.70e-06,
    1.04e-05,1.35e-05,1.69e-05,2.07e-05,2.49e-05,2.96e-05,
    3.49e-05,4.07e-05,4.72e-05,5.44e-05,6.25e-05,7.14e-05,
    8.14e-05,9.25e-05,1.04e-04,1.34e-04,1.51e-04,1.70e-04,
    1.91e-04,2.14e-04,2.41e-04,2.69e-04,3.02e-04,3.38e-04,
    3.78e-04,4.23e-04,4.72e-04,5.28e-04,5.89e-04,6.57e-04,
    7.33e-04,8.18e-04,9.12e-04,1.02e-03,1.13e-03,1.26e-03,
    1.41e-03,1.57e-03,1.75e-03,1.95e-03,2.17e-03,2.41e-03,
    2.69e-03,3.00e-03,3.34e-03,3.72e-03,4.14e-03,4.61e-03,
    5.14e-03,5.72e-03,6.36e-03,7.09e-03,7.89e-03,8.79e-03,
    1.00e-02
])
strainPoints2 = np.array([
    -5.10e-01,-5.00e-01,-4.97e-01,-4.41e-01,-3.90e-01,-3.46e-01,
    -3.06e-01,-2.71e-01,-2.40e-01,-2.13e-01,-1.88e-01,-1.67e-01,
    -1.47e-01,-1.31e-01,-1.15e-01,-1.03e-01,-9.07e-02,-8.03e-02,
    -7.11e-02,-6.30e-02,-5.57e-02,-4.93e-02,-4.36e-02,-3.86e-02,
    -3.41e-02,-3.02e-02,-2.67e-02,-2.36e-02,-2.09e-02,-1.84e-02,
    -1.63e-02,-1.44e-02,-1.12e-02,-9.92e-03,-8.75e-03,-7.71e-03,
    -6.79e-03,-5.97e-03,-5.25e-03,-4.61e-03,-4.05e-03,-3.55e-03,
    -3.10e-03,-2.71e-03,-2.36e-03,-2.05e-03,-1.78e-03,-1.54e-03,
    -1.32e-03,-1.13e-03,-9.67e-04,-8.18e-04,-6.86e-04,-5.69e-04,
    -4.65e-04,-3.37e-04,-2.92e-04,-2.20e-04,-1.56e-04,-1.00e-04,
    -5.00e-05,0.00e+00,5.00e-05,1.00e-04,1.56e-04,2.20e-04,
    2.92e-04,3.37e-04,4.65e-04,5.69e-04,6.86e-04,8.18e-04,
    9.67e-04,1.13e-03,1.32e-03,1.54e-03,1.78e-03,2.05e-03,
    2.36e-03,2.71e-03,3.10e-03,3.55e-03,4.05e-03,4.61e-03,
    5.25e-03,5.97e-03,6.79e-03,7.71e-03,8.75e-03,9.92e-03,
    1.12e-02,1.44e-02,1.63e-02,1.84e-02,2.09e-02,2.36e-02,
    2.67e-02,3.02e-02,3.41e-02,3.86e-02,4.36e-02,4.93e-02,
    5.57e-02,6.30e-02,7.11e-02,8.03e-02,9.07e-02,1.03e-01,
    1.15e-01,1.31e-01,1.47e-01,1.67e-01,1.88e-01,2.13e-01,
    2.40e-01,2.71e-01,3.06e-01,3.46e-01,3.90e-01,4.41e-01,
    5.00e-01
])
#--------------------------------------------------------------------------------------------------
# Create class for unit cell results

class Results:
    def __init__(self, csvPath, strainPoints):
        df = pd.read_csv(csvPath)
        threshold = 0
        self.frames = df.iloc[:, 0].to_numpy()
        self.N11 = np.where(np.abs(df.iloc[:,1].to_numpy())< threshold, 0,df.iloc[:, 1].to_numpy())
        self.N22 = np.where(np.abs(df.iloc[:,2].to_numpy())< threshold, 0,df.iloc[:, 2].to_numpy())
        self.N12 = np.where(np.abs(df.iloc[:,3].to_numpy())< threshold, 0,df.iloc[:, 3].to_numpy())
        self.M11 = np.where(np.abs(df.iloc[:,4].to_numpy())< threshold, 0,df.iloc[:, 4].to_numpy())
        self.M22 = np.where(np.abs(df.iloc[:,5].to_numpy())< threshold, 0,df.iloc[:, 5].to_numpy())
        self.M12 = np.where(np.abs(df.iloc[:,6].to_numpy())< threshold, 0,df.iloc[:, 6].to_numpy())
        self.E11 = df.iloc[:, 7].to_numpy()
        self.E22 = df.iloc[:, 8].to_numpy()
        self.E12 = df.iloc[:, 9].to_numpy()
        self.G11 = df.iloc[:, 10].to_numpy()
        self.G22 = df.iloc[:, 11].to_numpy()
        self.G12 = df.iloc[:, 12].to_numpy()
        
        self.interpN11 = None
        self.interpN22 = None
        self.interpN12 = None
        self.interpM11 = None
        self.interpM22 = None
        self.interpM12 = None
        self.interpolateExtrapolate(csvPath, strainPoints)

    def interpolateExtrapolate(self, csvPath, strainPoints):
        strainArrays = {
            'E11': self.E11,
            'E22': self.E22,
            'E12': self.E12,
            'G11': self.G11,
            'G22': self.G22,
            'G12': self.G12
        }
        stressArrays = {
            'N11': self.N11,
            'N22': self.N22,
            'N12': self.N12,
            'M11': self.M11,
            'M22': self.M22,
            'M12': self.M12
        }

        if csvPath == 'loadCase-1.csv':
            curStrain = strainArrays['E11']
        elif csvPath == 'loadCase-2.csv':
            curStrain = strainArrays['E22']
        elif csvPath == 'loadCase-3.csv':
            curStrain = strainArrays['E12']
        elif csvPath == 'loadCase-4.csv':
            curStrain = strainArrays['G11']
        elif csvPath == 'loadCase-5.csv':
            curStrain = strainArrays['G22']
        elif csvPath == 'loadCase-6.csv':
            curStrain = strainArrays['G12']
        else:
            print("No unit cell data found")
            
        interpResults = {}
        for key, stressArray in stressArrays.items():
            interpFunc = interp1d(curStrain, stressArray, fill_value="extrapolate")
            interpResults[key] = interpFunc(strainPoints)
        
        self.interpN11 = interpResults.get('N11')
        self.interpN22 = interpResults.get('N22')
        self.interpN12 = interpResults.get('N12')
        self.interpM11 = interpResults.get('M11')
        self.interpM22 = interpResults.get('M22')
        self.interpM12 = interpResults.get('M12')
#--------------------------------------------------------------------------------------------------
# Import and process unit cell results

loadCase1 = Results('loadCase-1.csv', strainPoints1)
loadCase2 = Results('loadCase-2.csv', strainPoints1)
loadCase3 = Results('loadCase-3.csv', strainPoints1)
loadCase4 = Results('loadCase-4.csv', strainPoints2)
loadCase5 = Results('loadCase-5.csv', strainPoints2)
loadCase6 = Results('loadCase-6.csv', strainPoints2)
#--------------------------------------------------------------------------------------------------
# Calculate stiffnesses

def stiffness(data, strainPoints):
    stif = [0.0]
    for n in range(1, len(strainPoints)):
        stif.append((data[n] - data[n-1]) / (strainPoints[n] - strainPoints[n-1]))
    return stif

A11 = stiffness(loadCase1.interpN11, strainPoints1) 
A21 = stiffness(loadCase1.interpN22, strainPoints1)
A12 = stiffness(loadCase2.interpN11, strainPoints1)
A22 = stiffness(loadCase2.interpN22, strainPoints1)
A33 = stiffness(loadCase3.interpN12, strainPoints1)

B11 = stiffness(loadCase4.interpN11, strainPoints2) 
B21 = stiffness(loadCase4.interpN22, strainPoints2)
B12 = stiffness(loadCase5.interpN11, strainPoints2)
B22 = stiffness(loadCase5.interpN22, strainPoints2)
B33 = stiffness(loadCase6.interpN12, strainPoints2)

C11 = stiffness(loadCase1.interpM11, strainPoints1) 
C21 = stiffness(loadCase1.interpM22, strainPoints1)
C12 = stiffness(loadCase2.interpM11, strainPoints1)
C22 = stiffness(loadCase2.interpM22, strainPoints1)
C33 = stiffness(loadCase3.interpM12, strainPoints1)

D11 = stiffness(loadCase4.interpM11, strainPoints2) 
D21 = stiffness(loadCase4.interpM22, strainPoints2)
D12 = stiffness(loadCase5.interpM11, strainPoints2)
D22 = stiffness(loadCase5.interpM22, strainPoints2)
D33 = stiffness(loadCase6.interpM12, strainPoints2)

open('stiffnessMatrix_'+str(panelName)+'.f', 'w').close()
def writeToFile(component, data):    
    with open('stiffnessMatrix_'+str(panelName)+'.f', 'a') as file: 
        values = data
        file.write(f'      {component}=(/')
        for i, value in enumerate(values):
            if i > 0 and i % 4 == 0: 
                file.write("\n     + ")
            file.write(f"{value: .6e}")
            if i < len(values) - 1:
                file.write(",")
        file.write("/)\n")

writeToFile('A11', A11)
writeToFile('A21', A21)
writeToFile('A12', A12)
writeToFile('A22', A22)
writeToFile('A33', A33)

writeToFile('B11', B11)
writeToFile('B21', B21)
writeToFile('B12', B12)
writeToFile('B22', B22)
writeToFile('B33', B33)

writeToFile('C11', C11)
writeToFile('C21', C21)
writeToFile('C12', C12)
writeToFile('C22', C22)
writeToFile('C33', C33)

writeToFile('D11', D11)
writeToFile('D21', D21)
writeToFile('D12', D12)
writeToFile('D22', D22)
writeToFile('D33', D33)

fig1, ([ax1, ax2], [ax3, ax4]) = plt.subplots(2, 2, figsize=(6.456, 3.937))
fig1.subplots_adjust(left=0.02, right=0.92, top=0.96, bottom=0.1, wspace=0.25, hspace=0.35)
ax1.plot(strainPoints1, loadCase1.interpN11, color='k', label='$A_{11}\\varepsilon^{(0)}_{xx}$')
ax1.plot(strainPoints1, loadCase2.interpN22, color='0.6', label='$A_{22}\\varepsilon^{(0)}_{yy}$')
ax1.plot(strainPoints1, loadCase2.interpN11, color='k', linestyle='-.', label='$A_{12}\\varepsilon^{(0)}_{yy}$')
ax1.plot(strainPoints1, loadCase1.interpN22, color='0.6', linestyle='-.', label='$A_{21}\\varepsilon^{(0)}_{xx}$')
ax1.plot(strainPoints1, loadCase3.interpN12, color='k', linestyle='--', label='$A_{33}\\gamma^{(0)}_{xy}$')
#ax1.grid()
ax1.set_xlabel('Extensional Strain')
ax1.set_ylabel('Stress Resultant [N/m]')
ax1.set_xlim([-0.01, 0.01])
ax1.set_ylim([-0.4e7, 0.4e7])
ax1.set_yticks([-0.4e7, -0.2e7, 0, 0.2e7, 0.4e7])
ax1.set_xticks([-0.01, -0.005, 0, 0.005, 0.01])

ax2.plot(strainPoints2, loadCase4.interpN11, color='k', label='$B_{11}\\varepsilon^{(1)}_{xx}$')
ax2.scatter([strainPoints2[19]], [loadCase4.interpN11[19]], marker='*', s=25, color='k')
ax2.plot(strainPoints2, loadCase5.interpN22, color='0.6', label='$B_{22}\\varepsilon^{(1)}_{yy}$')
ax2.plot(strainPoints2, loadCase5.interpN11, color='k', linestyle='-.', label='$B_{12}\\varepsilon^{(1)}_{yy}$')
ax2.plot(strainPoints2, loadCase4.interpN22, color='0.6', linestyle='-.', label='$B_{21}\\varepsilon^{(1)}_{xx}$')
ax2.plot(strainPoints2, loadCase6.interpN12, color='k', linestyle='--', label='$B_{33}\\gamma^{(1)}_{xy}$')
ax2.set_xlabel('Bending Strain')
ax2.set_ylabel('Stress Resultant [Nm/m]')
#ax2.grid()
ax2.set_xlim([-0.50, 0.50])
ax2.set_ylim([-1.2e6, 1.2e6])
ax2.set_yticks([-1.2e6, -0.6e6, 0, 0.6e6, 1.2e6])
ax2.set_xticks([-0.5, -0.25, 0, 0.25, 0.5])

ax3.plot(strainPoints1, loadCase1.interpM11, color='k', label='$C_{11}\\varepsilon^{(0)}_{xx}$')
ax3.plot(strainPoints1, loadCase2.interpM22, color='0.6', label='$C_{22}\\varepsilon^{(0)}_{yy}$')
ax3.plot(strainPoints1, loadCase2.interpM11, color='k', linestyle='-.', label='$C_{12}\\varepsilon^{(0)}_{yy}$')
ax3.plot(strainPoints1, loadCase1.interpM22, color='0.6', linestyle='-.', label='$C_{21}\\varepsilon^{(0)}_{xx}$')
ax3.plot(strainPoints1, loadCase3.interpM12, color='k', linestyle='--', label='$C_{33}\\gamma^{(0)}_{xy}$')
ax3.set_xlabel('Extensional Strain')
ax3.set_ylabel('Stress Resultant [N/m]')
#ax3.grid()
ax3.set_xlim([-0.01, 0.01])
ax3.set_ylim([-1.0e4, 1.0e4])
ax3.set_yticks([-1.0e4, -0.5e4, 0, 0.5e4, 1.0e4])
ax3.set_xticks([-0.01, -0.005, 0, 0.005, 0.01])

ax4.plot(strainPoints2, loadCase4.interpM11, color='k', label='$D_{11}\\varepsilon^{(1)}_{xx}$')
ax4.scatter([strainPoints2[19]], [loadCase4.interpM11[19]], marker='*', s=25, color='k')
ax4.plot(strainPoints2, loadCase5.interpM22, color='0.6', label='$D_{22}\\varepsilon^{(1)}_{yy}$')
ax4.plot(strainPoints2, loadCase5.interpM11, color='k', linestyle='-.', label='$D_{12}\\varepsilon^{(1)}_{yy}$')
ax4.plot(strainPoints2, loadCase4.interpM22, color='0.6', linestyle='-.', label='$D_{21}\\varepsilon^{(1)}_{xx}$')
ax4.plot(strainPoints2, loadCase6.interpM12, color='k', linestyle='--', label='$D_{33}\\gamma^{(1)}_{xy}$')
ax4.set_xlim([-0.5, 0.5])
ax4.set_ylim([-0.8e5, 0.8e5])
ax4.set_yticks([-0.8e5, -0.4e5, 0, 0.4e5, 0.8e5])
ax4.set_xticks([-0.5, -0.25, 0, 0.25, 0.5])
ax4.set_xlabel('Bending Strain')
ax4.set_ylabel('Stress Resultant [Nm/m]')
#ax4.grid()
labels = ['\\bf{a)}', '\\bf{b)}', '\\bf{c)}', '\\bf{d)}']
axes = [ax1, ax2, ax3, ax4]

fig1.align_labels()

for label, ax in zip(labels, axes):
    ax.text(-0.30, 1.08, label, transform=ax.transAxes,
            fontsize=11, va='top', ha='left')
for ax in axes:
    ax.set_box_aspect(.85)
    ax.legend(loc='lower right', bbox_to_anchor=(1.5, 0))

fig2, ([ax1, ax2], [ax3, ax4]) = plt.subplots(2, 2, figsize=(6.456, 3.937))
fig2.subplots_adjust(left=0.02, right=0.92, top=0.96, bottom=0.1, wspace=0.25, hspace=0.35)
ax1.plot(strainPoints1, A11, color='k', label='$A_{11}$')
ax1.plot(strainPoints1, A22, color='0.6', label='$A_{22}$')
ax1.plot(strainPoints1, A12, color='k', linestyle='-.', label='$A_{12}$')
ax1.plot(strainPoints1, A21, color='0.6', linestyle='-.', label='$A_{21}$')
ax1.plot(strainPoints1, A33, color='k', linestyle='--', label='$A_{33}$')
ax1.yaxis.set_major_locator(MaxNLocator(integer = True, nbins = 3))
ax1.scatter([strainPoints1[60]], [A11[60]], marker='o', facecolors='none', edgecolors='k', linewidths=0.8, s=20, color='k')
ax1.scatter([strainPoints1[60]], [A22[60]], marker='o', facecolors='none', edgecolors='0.6', linewidths=0.8, s=20, color='k')
#ax1.grid()
ax1.set_xlabel('Extensional Strain')
ax1.set_ylabel('$A_{ij}$ [N/m]')
ax1.set_xlim([-0.01, 0.01])
ax1.set_ylim([-0.5e9, 1.5e9])
ax1.set_yticks([-0.5e9, 0, 0.5e9, 1e9, 1.5e9])
ax1.set_xticks([-0.01, -0.005, 0, 0.005, 0.01])

ax2.plot(strainPoints2, B11, color='k', label='$B_{11}$')
ax2.scatter([strainPoints2[22]], [B11[22]], marker='*', s=25, color='k')
ax2.plot(strainPoints2, B22, color='0.6', label='$B_{22}$')
ax2.plot(strainPoints2, B12, color='k', linestyle='-.', label='$B_{12}$')
ax2.plot(strainPoints2, B21, color='0.6', linestyle='-.', label='$B_{21}$')
ax2.plot(strainPoints2, B33, color='k', linestyle='--', label='$B_{33}$')
ax2.yaxis.set_major_locator(MaxNLocator(integer = True, nbins = 3))
ax2.set_xlabel('Bending Strain')
ax2.set_ylabel('$B_{ij}$ [Nm]')
#ax2.grid()
ax2.set_xlim([-0.5, 0.5])
ax2.set_ylim([-0.6e7, 2e7])
ax2.set_yticks([-0.5e7, 0, 0.5e7, 1.0e7, 1.5e7, 2e7])
ax2.set_xticks([-0.5, -0.25, 0, 0.25, 0.5])

ax3.plot(strainPoints1, C11, color='k', label='$C_{11}$')
ax3.plot(strainPoints1, C22, color='0.6', label='$C_{22}$')
ax3.plot(strainPoints1, C12, color='k', linestyle='-.', label='$C_{12}$')
ax3.plot(strainPoints1, C21, color='0.6', linestyle='-.', label='$C_{21}$')
ax3.plot(strainPoints1, C33, color='k', linestyle='--', label='$C_{33}$')
ax3.yaxis.set_major_locator(MaxNLocator(integer = True, nbins = 3))
ax3.set_xlabel('Extensional Strain')
ax3.set_ylabel('$C_{ij}$ [N/m]')
#ax3.grid()
ax3.set_xlim([-0.01, 0.01])
ax3.set_ylim([-0.2e7, 2e7])
ax3.set_yticks([0, 0.5e7, 1e7, 1.5e7, 2e7])
ax3.set_xticks([-0.01, -0.005, 0, 0.005, 0.01])

ax4.plot(strainPoints2, D11, color='k', label='$D_{11}$')
ax4.scatter([strainPoints2[22]], [D11[22]], marker='*', s=25, color='k')
ax4.plot(strainPoints2, D22, color='0.6', label='$D_{22}$')
ax4.plot(strainPoints2, D12, color='k', linestyle='-.', label='$D_{12}$')
ax4.plot(strainPoints2, D21, color='0.6', linestyle='-.', label='$D_{21}$')
ax4.plot(strainPoints2, D33, color='k', linestyle='--', label='$D_{33}$')
ax4.yaxis.set_major_locator(MaxNLocator(integer = True, nbins = 3))
ax4.set_xlim([-0.5, 0.5])
ax4.set_ylim([-0.2e6, 1.6e6])
ax4.set_xlabel('Bending Strain')
ax4.set_ylabel('$D_{ij}$ [Nm]')
ax4.set_yticks([0, 0.5e6, 1.0e6, 1.5e6])
ax4.set_xticks([-0.5, -0.25, 0, 0.25, 0.5])
#ax4.grid()
labels = ['\\bf{a)}', '\\bf{b)}', '\\bf{c)}', '\\bf{d)}']
axes = [ax1, ax2, ax3, ax4]

fig2.align_labels()

for label, ax in zip(labels, axes):
    ax.text(-0.30, 1.08, label, transform=ax.transAxes,
            fontsize=11, va='top', ha='left')
for ax in axes:
    ax.set_box_aspect(.85)
    ax.legend(loc='lower right', bbox_to_anchor=(1.45, 0))

plt.show()
#--------------------------------------------------------------------------------------------------
# Calculate analytical transverse shear stiffnesses

h = tp + hw + tf

zNA = ((s * tp) * (0.5 * tp) + (tw * hw) * (tp + 0.5 * hw) + bf * tf * (tp + hw + 0.5 
    * tf)) / (s * tp + tw * hw + bf * tf)

Ip = (s * tp**3) / 12 + (s * tp) * (0.5 * tp - zNA)**2
Iw = (tw * hw**3) / 12 + (tw * hw) * (tp + 0.5 * hw - zNA)**2
If = (bf * tf**3) / 12 + (bf * tf) * (tp + hw + 0.5 * tf - zNA)**2

Iyy = Ip + Iw + If

kxz = (2 * Iyy * tw) / ((tw * tp + tw * tf + tw * hw) * (tw * zNA**2 + (2 * s * tp - 2 * tw * tp) 
    * zNA + (tw * tp**2 - s * tp * tp))) 
kyz = 5/6

G = E / (2 * (1 + nu))

DQx = kxz * G * (tp + (tw * hw) / s + (bf * tf) / s)
DQy = kyz * G * tp
#--------------------------------------------------------------------------------------------------
# Write to file

n = 62
with open('sectionDefinition_'+str(panelName)+'.txt', 'w') as file:   
    file.write(f"*SHELL GENERAL SECTION, ELSET = , USER, UNSYMM, VARIABLES = 52, "
        f"PROPERTIES = 36\n")
    file.write(f"{h:.4f},\n")
    file.write(f"{A11[n]:.4e}, {A12[n]:.4e}, {0.0000:.4e}, {B11[n]:.4e}, {B12[n]:.4e}, "
        f"{0.0000:.4e}, {A21[n]:.4e}, {A22[n]:.4e}\n")
    file.write(f"{0.0000:.4e}, {B21[n]:.4e}, {B22[n]:.4e}, {0.0000:.4e}, {0.0000:.4e}, "
        f"{0.0000:.4e}, {A33[n]:.4e}, {0.0000:.4e}\n")
    file.write(f"{0.0000:.4e}, {B33[n]:.4e}, {C11[n]:.4e}, {C12[n]:.4e}, {0.0000:.4e}, "
        f"{D11[n]:.4e}, {D12[n]:.4e}, {0.0000:.4e}\n")
    file.write(f"{C21[n]:.4e}, {C22[n]:.4e}, {0.0000:.4e}, {D21[n]:.4e}, {D22[n]:.4e}, "
        f"{0.0000:.4e}, {0.0000:.4e}, {0.0000:.4e}\n")
    file.write(f"{C33[n]:.4e}, {0.0000:.4e}, {0.0000:.4e}, {D33[n]:.4e}\n")
    file.write(f"*TRANSVERSE SHEAR STIFFNESS\n")
    file.write(f"{DQx:.4e}, {DQy:.4e}")
# ---------------------------------------- END OF SCRIPT ---------------------------------------- #

fig1.savefig('Figure_12.pdf', transparent=True)
fig2.savefig('Figure_13.pdf', transparent=True)