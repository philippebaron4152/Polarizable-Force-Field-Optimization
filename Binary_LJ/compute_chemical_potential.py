import numpy as np
import matplotlib.pyplot as plt
import os
import lammps_logfile
import pymbar
from scipy.special import roots_legendre
from scipy.optimize import curve_fit, fsolve

import warnings
warnings.filterwarnings("ignore")

L_db = (9.9**3) * 10 ** (-30)                   # De Broglie Wavelength (consistency required between two phases)
kB_kJ_mol = 8.314462618 * 10 ** (-3)
kB_kcal_mol = 1.987204259 * 10 ** (-3)
NA = 6.0221408 * 10 ** 23
kB_J = 8.314462618 / NA
T = 50                                          # K
N_tot = 900                                       # Number of solvent molecules

BETA = 1 / (kB_kcal_mol * T)
BETA2 = 1 / (kB_kJ_mol * T)

def acf(data, N):
    
    """
    Compute the autocorrelation function of a quantity by Fourier Transform.
    """


    T = np.arange(0, N, 1)
    size = 2 ** np.ceil(np.log2(2*len(data) - 1)).astype('int')

    var = np.var(data)
    ndata = data - np.mean(data)
    fft = np.fft.fft(ndata, size)
    pwr = np.abs(fft) ** 2
    C = np.fft.ifft(pwr).real / var / len(data)
    
    C = C[:N]
    return T, C

def acf_analysis(data, N):

    """
    Average and integrate autocorrelation function to estimate correlation time. 

    data (double array) -> time-dependent trajectory of some quantity
    N (int) -> number of elements to consider for autocorrelation function
    """
     
    C = np.zeros(N); T = np.zeros(N)
    n = int(len(data) / N)
    for j in range(n):
        T, C_ = acf(data[j*N:], N)
        C += C_
    C /= n

    corr_time = 0
    for i in range(len(C)):
        corr_time += (1 - T[i]/len(C)) * C[i]
        
    return T, C, corr_time

def extrapolation_error(Xs, A, A_std, N_samples):

    """
    Estimate the extrapolated value and its uncertainty when linearly extrapolating to the thermodynamic limit (A vs. 1/N)

    Xs (double array) -> x-values of system size related parameter (L or N)
    A (double array) -> Function values being extrapolated
    A_std (double array) -> Uncertainties of function values being extrapolated
    N_samples (int) -> Number of samples from a normal distribution to be taken about each point
    """

    inf_dist = []
    slopes = []
    
    final_coeffs = np.polyfit(1/Xs, A, 1)

    for i in range(N_samples):
        points = []
        
        for j, a in enumerate(A):
            val = np.random.normal(a, A_std[j])
            points.append(val)
        
        coeffs = np.polyfit(1/Xs, points, 1)
        slopes.append(coeffs[0])
        inf_dist.append(coeffs[1])
    
    return final_coeffs, final_coeffs[1], np.std(inf_dist)

def fit_func(x, a, b, c, d):
    return np.log(a * x ** (1/2) + b * x ** 2 + c * x + d)

def mu_model(x, mu, V_f):
    bounds=([0, 0, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf])
    params, _ = curve_fit(fit_func, x, mu, bounds=bounds)
    MU_model = lambda y: fit_func(y,*params) + np.log((y * N_tot + 0.5) * L_db / V_f(y))
    return MU_model

if __name__ == "__main__":
    # Set parameters for Bennett Acceptance Ratio solver
    solver_options = {"maximum_iterations":50000,"verbose":True}
    solver_protocol = {"method":"adaptive","options":solver_options}

    ###### COMPUTE AQUEOUS CHEMICAL POTENTIAL ######

    current_dir = os.getcwd()
    lambdas = np.genfromtxt(current_dir + "/lambdas.txt")
    conc_dirs = [d for d in np.sort(os.listdir(current_dir + "/aqueous/")) if "A" in d]

    xA = np.zeros(len(conc_dirs)); Vs = np.zeros(len(conc_dirs))
    MU_ideal = np.zeros(len(conc_dirs))
    MU_ex = np.zeros(len(conc_dirs))
    MU_ERR = np.zeros(len(conc_dirs))

    for ii, dir in enumerate(conc_dirs):
        path = current_dir + "/aqueous/" + dir
        N_A = float(dir[:-1])
        xA[ii] = N_A / N_tot

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
        Vs[ii] = V
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

        MU_ideal[ii] = mu_ideal
        MU_ex[ii] = mu
        MU_ERR[ii] = np.sqrt(mu_err)

    MU = MU_ex + MU_ideal
    plt.figure(figsize=(8,5))
    plt.scatter(xA, MU, c='r', marker='s', s=70)
    plt.errorbar(xA, MU, MU_ERR, c='r', capsize=3, linestyle="")

    V_params = np.polyfit(xA, Vs, 1)
    V_f = lambda x: V_params[0] * x + V_params[1]
    f = mu_model(xA, MU_ex, V_f)
    X = np.linspace(np.min(xA)-0.005, np.max(xA)+0.05, 1000)
    plt.plot(X, f(X), "r--")


    ###### COMPUTE CRYSTAL CHEMICAL POTENTIAL ######

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
    plt.axhline(y=mu_solid_LJ - mu_solid_LJ_err, c='r')
    plt.axhline(y=mu_solid_LJ + mu_solid_LJ_err, c='r')

    ###### COMPUTE SOLUBILITY ######

    M = 20
    samples = np.random.normal(0.0, 1.0, (M, len(xA)))
    samples_crys = mu_solid_LJ_err*np.random.normal(0.0, 1.0, M) + mu_solid_LJ

    for i in range(len(xA)):
        samples[:,i] = samples[:,i] * MU_ERR[i] + MU_ex[i]

    MU_model = mu_model(xA, MU_ex, V_f)
    f = lambda x: MU_model(x) - mu_solid_LJ
    initial_guess = np.mean(xA); m_s = initial_guess
    while initial_guess == m_s:
        new_initial_guess = np.random.uniform(0.0, 1.0) * (xA[-1] - xA[0]) + xA[0]
        m_s = fsolve(f, new_initial_guess)
        initial_guess = new_initial_guess

    solubility = np.zeros(M)
    for i in range(M):
        MU_model = mu_model(xA, samples[i, :], V_f)
        mu_crys_i = samples_crys[i]
        f = lambda x: MU_model(x) - mu_crys_i
        solubility[i] = fsolve(f, new_initial_guess)

    m_s = np.mean(solubility)
    solubility_err = np.std(solubility) / np.sqrt(M)
    plt.axvline(x= m_s - solubility_err, c="k", linestyle=":")
    plt.axvline(x= m_s + solubility_err, c="k", linestyle=":")
    plt.text(m_s + 0.01, MU[0], r"$m_s = {}  \pm {}$".format(np.round(m_s,5), np.round(solubility_err,5)), fontsize=16)

    plt.xlabel(r"$x_A$", fontsize=16)
    plt.ylabel(r"$\mu/k_BT$", fontsize=16)
    plt.savefig("mu.png", dpi=700)