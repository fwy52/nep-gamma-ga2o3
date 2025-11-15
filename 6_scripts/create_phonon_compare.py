import os
import numpy as np
from pylab import *
from ase.io import read,write

cx,cy,cz = 6, 6, 2
npoints = 600       
special_points = {'G': [0, 0, 0], 'X': [0, 0.5, 0],'M': [0.5, 0.5, 0], 'G': [0, 0, 0], 'R': [0.5, 0.5, 0.5],'X': [0, 0.5, 0],'M': [0.5, 0.5, 0],'R': [0.5, 0.5, 0.5]}  
points_path = 'GXMGRXMR'   

uc = read('POSCAR-unitcell') 
struc = uc * (cx,cy,cz)
write("model.xyz", struc)

with open('basis.in', 'w') as f:
    f.write(f"{len(uc)}\n")
    for i, mass in enumerate(uc.get_masses()):
        f.write(f"{i} {mass}\n")
    for _ in range(cx*cy*cz): 
        for i in range(len(uc)):
            f.write(f"{i}\n")

path = uc.cell.bandpath(path=points_path, npoints = npoints, special_points=special_points) 
kpath, sym_points, labels = path.get_linear_kpoint_axis()
gpumd_kpts = np.matmul(path.kpts, uc.cell.reciprocal() * 2 * np.pi)
gpumd_kpts[np.abs(gpumd_kpts) < 1e-15] = 0.0
np.savetxt('kpoints.in', gpumd_kpts, header=str(npoints), comments='', fmt='%g')

data = np.loadtxt("omega2.out")

for i in range(len(data)):
    for j in range(len(data[0])):
        data[i, j] = np.sqrt(abs(data[i, j])) / (2 * np.pi) * np.sign(data[i, j])
nu = data

np.savetxt('y_omega.dat', nu, fmt='%g')
np.savetxt('x_omega.dat', kpath, fmt='%g')


import numpy as np
import matplotlib.pyplot as plt
import os

figure(figsize=(9, 8))
if os.path.exists('phonon.out'):
    data_vasp = np.loadtxt('phonon.out')
    vasp_path = data_vasp[:, 0] / max(data_vasp[:, 0]) * max(kpath)
    plot(vasp_path, data_vasp[:, 1], color='orange', linestyle='--', label='DFT')
else:
    None
plot(kpath, nu, color='red', lw=1, label='NEP')
xlim([0, max(kpath)])
for sym_point in sym_points[1:-1]:
    plt.axvline(sym_point, color='black', linestyle='--') 
gca().set_xticks(sym_points)
gca().set_xticklabels([r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'R', 'X', 'M', 'R'])


from matplotlib.ticker import MultipleLocator


plt.gca().yaxis.set_major_locator(MultipleLocator(3))

plt.gca().yaxis.set_minor_locator(MultipleLocator(1))

plt.ylim(0, 27) 


plt.tick_params(axis='y', which='major', length=6, width=1, colors='k') 
plt.tick_params(axis='y', which='minor', length=3, width=0.8, colors='k') 
# =======================================

ylabel(r'$\nu$ (THz)', fontsize=15)
'''
if os.path.exists('phonon.out'):
    legend(['DFT', 'NEP'])
else:
   legend(['NEP'])
'''
savefig('phonon.png', dpi=300, bbox_inches='tight')
