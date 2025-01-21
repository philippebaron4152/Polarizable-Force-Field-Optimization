import numpy as np
import matplotlib.pyplot as plt
import packmol_to_lammps as pl
import argparse

import os
import shutil
import glob
from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove

MM_Na = 22.989769
MM_Cl = 35.453
MM_Ar = 39.948
NA = 6.02214076 * 10 ** (23)

def replace(file_path, pattern, subst):
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    copymode(file_path, abs_path)
    remove(file_path)
    move(abs_path, file_path)

def create_NaCl_crystal(a, N_pairs, lattice, plot_flag=False):
    Cl_0 = []
    Na_0 = []
    N = 2 * N_pairs

    for i in range(N):
        y = i * a / 2 
        z = 0
        if i % 2 == 1:
            for j in range(int(N/2)):
                x = a / 2 + j * a
                coord = [x, y, z]
                Cl_0.append(coord)
        else:
            for j in range(int(N/2)):
                x = j * a
                coord = [x, y, z]
                Cl_0.append(coord)
    Cl_0 = np.array(Cl_0)

    for i in range(N):
        y = i * a / 2
        z = 0
        if i % 2 == 0:
            for j in range(int(N/2)):
                x = a / 2 + j * a
                coord = [x, y, z]
                Na_0.append(coord)
        else:
            for j in range(int(N/2)):
                x = j * a
                coord = [x, y, z]
                Na_0.append(coord)
    Na_0 = np.array(Na_0)

    if plot_flag:
        plt.figure(figsize=(8,8))
        plt.scatter(Na_0[:,0], Na_0[:,1], c = 'g', s = 50)
        plt.scatter(Cl_0[:,0], Cl_0[:,1], c = 'b', s = 300)

    Na = Na_0.copy()
    Cl = Cl_0.copy()
    for k in range(N):
        if k == 0:
            continue
        elif k % 2 == 1:
            arr1 = Cl_0.copy()
            arr1[:, 2] = arr1[:,2] + k * a / 2
            Na = np.vstack((Na, arr1))

            arr2 = Na_0.copy()
            arr2[:, 2] = arr2[:,2] + k * a / 2
            Cl = np.vstack((Cl, arr2))
        else:
            arr1 = Na_0.copy()
            arr1[:, 2] = arr1[:,2] + k * a / 2
            Na = np.vstack((Na, arr1))

            arr2 = Cl_0.copy()
            arr2[:, 2] = arr2[:,2] + k * a / 2
            Cl = np.vstack((Cl, arr2))
    
    min_coord = np.min([np.min(Na), np.min(Cl)])
    max_coord = np.max([np.max(Na), np.max(Cl)])
    lo = min_coord - a/4
    hi = max_coord + a/4
    L = hi - lo
    V = L ** 3
    
    N_Na = N_Cl = N ** 3 / 2
    
    if lattice == "FCC":
        rho = 1000 * ((N_Na * MM_Ar) / NA) / (V * 10 ** (-24))
    else:
        rho = 1000 * ((N_Na * MM_Na + N_Cl * MM_Cl) / NA) / (V * 10 ** (-24))
    return Na, Cl, lo, hi, rho

def lattice_setup(obj_rho, N_pairs, path, lattice, bonds, cutoff=False):
    N_atoms = (N_pairs * 2) ** 3
    a1 = 3; a2 = 7
    val = 0
    while abs(val - obj_rho) > 0.0001:
        _, _, _, _, rho1 = create_NaCl_crystal(a1, N_pairs, lattice)
        _, _, _, _, rho2 = create_NaCl_crystal(a2, N_pairs, lattice)

        val = (rho1 + rho2) / 2

        _, _, _, _, rho = create_NaCl_crystal((a1 + a2)/2, N_pairs, lattice)
        if rho > obj_rho:
            a1 = (a1 + a2)/2
        else:
            a2 = (a1 + a2)/2

    A = (a1 + a2) / 2
    print('a* = ' + str(A))

    Na, Cl, lo, hi, rho = create_NaCl_crystal(A, N_pairs, lattice)
    print('DENSITY = ' + str(rho) + ' kg/m^3')
    print('N_ATOMS = ' + str(N_atoms))
    print('DIMENSIONS = ', lo, hi)
    print('CUTOFF = ' + str((hi - lo) / 2))
    print('BOX LENGTH = ' + str(hi - lo))

    with open(path + '/crystal.xyz', 'w+') as f:
        f.write(str(len(Na)) + '\n')
        f.write('\n')
        for i, a in enumerate(Na):
                f.write('A   ' + str(a[0]) + '  ' + str(a[1]) + '  ' + str(a[2]) + '\n')
    f.close()
    print('CRYSTAL XYZ written to: ' + path)  

    pl.generate_data_lammps(path + '/crystal.xyz', path + '/data_crystal.lammps', '', bonds)

    replace(path + "/data_crystal.lammps", "replace1", str(lo))
    replace(path + "/data_crystal.lammps", "replace2", str(hi))

    if cutoff:
        dirs = os.listdir(path)
        for d in dirs:
            if d[:2] == "in":
                replace(path + "/" + d, "CUTOFF", str((hi - lo) / 2))

argParser = argparse.ArgumentParser()
argParser.add_argument("-p", "--path", type=str, help="save_path")
argParser.add_argument("-n", "--pairs", type=int, help="number of ion pairs per dimension")
argParser.add_argument("-r", "--rho", type=float, help="target density (kg/m^3)")
argParser.add_argument("-l", "--lattice", type=str, help="lattice type (FCC or NACL)")
argParser.add_argument("-b", "--bonds", type=str, help="bonds present")
argParser.add_argument('--cutoff', action=argparse.BooleanOptionalAction)

args = argParser.parse_args()

lattice_setup(args.rho,args.pairs,args.path, args.lattice, args.bonds, args.cutoff)
