import numpy as np
import matplotlib.pyplot as plt
import os
import lammps_logfile
import pymbar
from scipy.special import roots_legendre
from compute_chemical_potential import acf_analysis, extrapolation_error

import warnings
warnings.filterwarnings("ignore")

L_db = (9.9**3) * 10 ** (-30)                   # De Broglie Wavelength (consistency required between two phases)
kB_kJ_mol = 8.314462618 * 10 ** (-3)
kB_kcal_mol = 1.987204259 * 10 ** (-3)
NA = 6.0221408 * 10 ** 23
kB_J = 8.314462618 / NA
T = 50                                          # K
N_tot = 900                                     # Number of total particles

BETA = 1 / (kB_kcal_mol * T)
BETA2 = 1 / (kB_kJ_mol * T)

# Set parameters for Bennett Acceptance Ratio solver
solver_options = {"maximum_iterations":50000,"verbose":True}
solver_protocol = {"method":"adaptive","options":solver_options}

###### COMPUTE AQUEOUS CHEMICAL POTENTIAL ######

current_dir = os.getcwd()
lambdas = np.genfromtxt(current_dir + "/lambdas.txt")

path = current_dir + "/aqueous/solv/"
N_A = 1
xA = N_A / N_tot

nn = 1; corr_time = 100
while corr_time > 1.5:
    taus = np.zeros(len(lambdas))
    for i, l in enumerate(lambdas):
        p = path + "/lambda_" + str(i+1) + "/log_.lammps"
        U = np.array(lammps_logfile.File(p).get("PotEng"))[::nn]
        _, _, tau = acf_analysis(U, 250)
        taus[i] = tau
    corr_time = np.mean(taus)
    nn += 1

path_ideal = path + "/lambda_1/log_run.lammps"
V = np.mean(np.array(lammps_logfile.File(path_ideal).get("Volume")[1000:]) * (10 ** (-30)))
mu_ideal = np.log((N_A + 0.5) * L_db / V) 

mu = 0; mu_err = 0; t=100
for i, l in enumerate(lambdas[:-1]):
    p_00 = path + "/lambda_" + str(i+1) + "/log_.lammps"
    p_01 = path + "/lambda_" + str(i+1) + "/log_up.lammps"
    p_10 = path + "/lambda_" + str(i+2) + "/log_down.lammps"
    p_11 = path + "/lambda_" + str(i+2) + "/log_.lammps"

    U_00 = np.array(lammps_logfile.File(p_00).get("PotEng"))[t::nn] * BETA
    U_01 = np.array(lammps_logfile.File(p_01).get("PotEng"))[t::nn] * BETA
    U_10 = np.array(lammps_logfile.File(p_10).get("PotEng"))[t::nn] * BETA
    U_11 = np.array(lammps_logfile.File(p_11).get("PotEng"))[t::nn] * BETA
    L = np.min([len(U_00), len(U_01), len(U_10), len(U_11)])
    
    Ukn = np.array([np.hstack((U_00[:L], U_10[:L])), np.hstack((U_01[:L], U_11[:L]))])
    Nk = np.array([L, L])
    
    mbar = pymbar.MBAR(Ukn, Nk, solver_protocol = (solver_protocol,))
    results = mbar.compute_free_energy_differences()

    mu += results['Delta_f'][0,-1]
    mu_err += (results['dDelta_f'][0,-1]) ** 2

MU = mu_ideal + mu
MU_ERR = np.sqrt(mu_err)

Ns = np.array([4, 6])
mu_s = np.array(np.zeros_like(Ns), dtype=np.float64); mu_s_err= np.array(np.zeros_like(Ns), dtype=np.float64)
for jj, N in enumerate(Ns):
    N_crys = int((2 * N) ** 3 / 2); K = 5
    fep_path = current_dir + "/crystal/" + str(N_crys) + "/fep/"
    lattice_path = current_dir + "/crystal/" + str(N_crys) + "/lattice/"
    kB_J = 1.380649 * 10 ** (-23); lambda_s = 500
        
    beta = 1 / (kB_J * T) 
    lambda_s_a0 = lambda_s * kB_J * T * 10 ** (20)
    V =  np.mean(np.array(lammps_logfile.File(current_dir + "/crystal/" + str(N_crys) + "/density/log_run.lammps").get("Volume")[1000:])) * 10 ** (-30)
    A0_1 = np.log((N_crys * L_db) / V) 
    A0_2 = (3 / 2) * (N_crys - 1) * np.log((L_db ** (2/3) * beta * lambda_s_a0) / np.pi)
    A0 = (A0_1 + A0_2) / N_crys

    U_fep = np.array(lammps_logfile.File(fep_path + "log_rerun.lammps").get("PotEng"))[100:]
    U_lattice = np.array(lammps_logfile.File(lattice_path + "log_run.lammps").get("PotEng"))[1]
    dA1s = np.zeros(K); n = int(len(U_fep) / K)
    for ii in range(K):
        exp_avg = np.mean(np.exp(BETA * (U_lattice - U_fep)[ii*n:n*(ii+1)]))
        dA1s[ii] = BETA * U_lattice - np.log(exp_avg) 
    dA1 = np.mean(dA1s) / (N_crys)
    dA1_err = np.std(dA1s) / np.sqrt(K) / (N_crys)

    lambdas, weights = roots_legendre(18)
    lambdas = 0.5 * lambdas + 0.5
    dA2 = np.zeros(K)
    for i, l in enumerate(lambdas):
        path = current_dir + "/crystal/" + str(N_crys) + "/springs/k" + str(i+1) + "/msd.out"
        msd = np.genfromtxt(path)[5001:]
        MSD = np.reshape(msd[:,1], (K, int(len(msd)/K)))
        dA2 += -0.5 * weights[i] * lambda_s * np.mean(MSD, axis=1)
    dA2_err = np.std(dA2) / np.sqrt(K); dA2 = np.mean(dA2)

    mu_s[jj] = A0 + dA1 + dA2
    mu_s_err[jj] = np.sqrt(dA1_err ** 2 + dA2_err ** 2)

N_arr = (2 * Ns) ** 3 / 2
params, mu_solid_LJ, mu_solid_LJ_err = extrapolation_error(np.array(N_arr), mu_s, mu_s_err, 10)

m_s = 0.01 * np.abs(MU - mu_solid_LJ) ** 2
N_A_m_s = int(m_s * N_tot)
dN = int(N_A_m_s * 0.2)
N_A_start_val = N_A_m_s - 3*dN 
N_lim = 0.95 * N_tot

N_As = np.array([N_A_start_val + j*dN for j in range(6)])
N_As = N_As[(N_As > 0) & (N_As < N_lim)] 
while len(N_As) < 5:
    dN /= 1.5
    N_As = np.array([N_A_start_val + j*dN for j in range(6)])
    N_As = N_As[(N_As > 0) & (N_As < N_lim)]    

with open("conc.txt", "w+") as f:
    for N_A in N_As:
        f.write(str(int(N_A)) + "\n")
