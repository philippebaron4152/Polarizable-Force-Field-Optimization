import numpy as np
import linecache
import matplotlib.pyplot as plt
import os
import shutil
import glob
from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove
from matplotlib.patches import Rectangle

def replace(file_path, pattern, subst):
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    copymode(file_path, abs_path)
    remove(file_path)
    move(abs_path, file_path)
    
def compute_rdf(box, coords, group1, group2, nbins, rlim):
    L = box[0][1] - box[0][0]
    Np = len(group2)
    rho = Np / (L ** 3)
    bins = np.linspace(0, rlim, nbins)
    dR = (bins[1] - bins[0])
    bin_centers = bins[:-1] + dR / 2

    rdf = np.zeros(nbins-1)
    for I in group1:
        Rs = []
        for jj, J in enumerate(group2): 
                
            rI = coords[I, 2:5]; rJ = coords[J, 2:5]
            dr = rI - rJ
            dr_pbc = dr - np.round(dr / L) * L
            if np.linalg.norm(dr_pbc) < 0.1:
                continue
            Rs.append(np.linalg.norm(dr_pbc))
        Rs = np.array(Rs)
        counts, _ = np.histogram(Rs[Rs < rlim], bins=bins)
        rdf += counts / (4 * np.pi * (bin_centers ** 2) * dR * rho)
        
    return bin_centers, rdf / len(group1)
 
def load_data(path):
    coords = []
    box = []
    Np = int(float(linecache.getline(path, 4)))
    buffer = 9
    with open(path, 'r') as f:
        lines = f.readlines()
        for ii in range(int(len(lines)/(Np + buffer))):
#             print(ii, int(len(lines)/(Np + buffer)))
            
            sub_box = []
            sub_box.append(np.array(lines[ii * (Np + buffer) + 5].split(), dtype=np.float64))
            sub_box.append(np.array(lines[ii * (Np + buffer) + 6].split(), dtype=np.float64))
            sub_box.append(np.array(lines[ii * (Np + buffer) + 7].split(), dtype=np.float64))
            box.append(sub_box)
            
            arr1 = lines[ii * (Np + buffer) + buffer : (ii + 1) * (Np + buffer)]
            arr2 = []
            for item in arr1:
                arr2.append(np.array(item.split(), dtype=np.float64))
            coords.append(arr2)
    print(path + " LOADED")
    return coords, box

def read_xyz(path):
    elements = []
    coords = []
    with open(path, 'r') as f:
        for i, line in enumerate(f):
            if i < 2:
                continue
            ele, x, y, z = line.strip().split()
            point = (float(x), float(y), float(z))
            elements.append(ele)
            coords.append(point)
    return np.array(elements), np.array(coords)

def extrapolation_error(Xs, A, A_std, N_samples):
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

def split_array(A, n):
    B = []
    for i in range(int(len(A)/n)):
        try:
            B.append(A[i*n : (i+1)*n])
        except:
            continue
    return B

def block_average(arr, plot_flag=True, N=300):
    delta_A_bar = np.var(arr,ddof=1)
    
    num_blocks = np.arange(N, 10, -1)
    tau_B = (len(arr) / num_blocks).astype(int)
    f_B = []

    for i, b in enumerate(tau_B):
        sub_means = []
        split_arr = split_array(arr, b)
        for sub_arr in split_arr:
            sub_means.append(np.mean(sub_arr))
        delta_Ab_bar = np.var(sub_means)
        f_B.append((b / 2) * (delta_Ab_bar / delta_A_bar))

    B = np.array(f_B)
    tau = np.max(B)

    if plot_flag == True:
        plt.figure()
        plt.scatter(tau_B, B)
        plt.axhline(y=tau)
        plt.ylabel(r'$\tau$')
        plt.xlabel(r'$\tau_B$')
        plt.show()

    variance = (delta_A_bar/len(arr)) * (1 + (2 * tau)) 
    error = np.sqrt(variance)
    
    return tau, error, tau_B, B

def mix_kong(A1, B1, A2, B2):
    k1 = -B1 / (B1 + B2)
    k2 = -B2 / (B1 + B2)
    Aij = 0.5 * (A1 * ((A1 * B1) / (A2 * B2)) ** k1 + A2 * ((A2 * B2) / (A1 * B1)) ** k2)
    Bij = 2 * B1 * B2 / (B1 + B2)
    return Aij, Bij
    
def mix_berthelot(C1, C2):
    Cij = np.sqrt(C1 * C2)
    return Cij

def n_salt(m, nw):
    return m * nw * (18.01528 / 1000)

def n_water(m, ns):
    return ns/(m * (18.01528 / 1000))

def molality(n_salt, nw):
    return (n_salt) / (nw * (18.01528 / 1000))

def acf(data, N):
    T = np.arange(0, N, 1)
    size = 2 ** np.ceil(np.log2(2*len(data) - 1)).astype('int')

    var = np.var(data)
    ndata = data - np.mean(data)
    fft = np.fft.fft(ndata, size)
    pwr = np.abs(fft) ** 2
    C = np.fft.ifft(pwr).real / var / len(data)
    
    C = C[:N]
    return T, C

def tau(T, C):
    tau = 0
    for i in range(len(C)):
        tau += (1 - T[i]/len(C)) * C[i]
    return tau

def acf_analysis(A, acf_len, plot_flag=True):
     
    C = np.zeros(acf_len); T = np.zeros(acf_len)
    n = int(len(A) / acf_len)
    for j in range(n):
        T, C_ = acf(A[j*acf_len:], acf_len)
        C += C_
    C /= n

    corr_time = tau(T, C)
    
    if plot_flag:
        print("CORRELATION TIME = " + str(corr_time) + " STEPS")
        plt.plot(T, C)
        
    return T, C, corr_time

def fermi_dirac(x, beta):
    return 1 / (1 + np.exp(beta * x))

def outfile_to_xyz(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        i0 = 0; i1 = 0; i2 = 0
        box_line = ''
        for k, line in enumerate(lines):
            if 'Atoms' in line:
                i0 = k
            if 'Bonds' in line:
                i1 = k
            if 'Velocities' in line:
                i2 = k
            if 'xlo xhi' in line:
                box_line = line
        
        box = np.array(box_line.split()[:2]).astype(np.float64)
        
        ind = int(np.min([i1, i2]))
        coords = lines[i0+2:ind-1]
        
        temp_arr = []
        for row in coords:
            arr = row.split()
            temp_arr.append([int(float(arr[0])), int(float(arr[2])), float(arr[4]), float(arr[5]), float(arr[6])]) 
        temp_arr = np.array(temp_arr)
        sorted_inds = np.argsort(temp_arr[:,0])
        A = temp_arr[sorted_inds, :]
        
        return A, box
    
def density_profile(C, B, START, STOP, sigma=3.405, plot_flag=False, species=[1., 2.]):
    coords = np.array(C[0]); box = np.array(B[0])
    solute_inds = np.where(coords[:,1] == 1.)[0]
    vals, edges = np.histogram(coords[solute_inds, 4], 250)
    X = edges[:-1] + (edges[1] - edges[0]) / 2
    
    midpoint = np.max(vals) / 2
    
    crystal_lim1 = np.min(X[vals > midpoint])
    crystal_lim2 = np.max(X[vals > midpoint])
    
    if np.abs(crystal_lim1 - crystal_lim2) > (np.max(X) / 2):
        vals2 = np.roll(vals, int(len(X)/2))
        crystal_lim1 = np.min(X[vals2 > midpoint])
        crystal_lim2 = np.max(X[vals2 > midpoint])
        midpoint_crystal = 0.5*(crystal_lim1 + crystal_lim2) - (np.max(X) / 2)
    else:
        midpoint_crystal = 0.5*(crystal_lim1 + crystal_lim2)
    
    ind = np.argmin(np.abs(coords[solute_inds, 4] - midpoint_crystal))
    tracked_particle = solute_inds[ind]
    
    nbins=150
    profiles = np.zeros((2, nbins))
    dzs = []

    for i in range(START, STOP):
        for jj, s in enumerate(species):
            coords = np.array(C[i]); box = np.array(B[i])
            inds = np.where((coords[:,1] == s))[0]

            reference_z = coords[tracked_particle, 4]
            Lz = box[2,1] - box[2,0]
            edges = np.linspace(0, Lz, nbins+1)
            dzs.append(edges[1]-edges[0])

            z_coords = coords[inds, 4] - reference_z
            z_coords[z_coords < 0] += Lz
            
            vals, _ = np.histogram(z_coords, edges)

            X = edges[:-1] + (edges[1] - edges[0]) / 2

            profiles[jj, :] += vals
    
    dz = np.mean(dzs)
    profiles /= (STOP-START)
    X = np.arange(0, nbins) * dz
    
    thresh = (np.max(profiles[0,:]) + np.min(profiles[0,:])) / 2

    peaks = np.where(profiles[0,:] > thresh)[0]
    left = np.where(peaks < (nbins)/2)[0]
    right = np.where(peaks > (nbins)/2)[0]
    
    ll_a = X[np.max(peaks[left])] + 3*sigma; rl_a = X[np.min(peaks[right])] - 3*sigma
    ll_b = X[np.max(peaks[left])] + 4*sigma; rl_b = X[np.min(peaks[right])] - 4*sigma
    
    n_A = np.mean(profiles[0, int(ll_a/dz):int(rl_a/dz)])
    n_B = np.mean(profiles[1, int(ll_b/dz):int(rl_b/dz)])

    if plot_flag:
        plt.plot(X, profiles[0,:], label='A')
        plt.plot(X, profiles[1,:], label='B')

        plt.axvline(x=rl_a, c='k', linestyle=':')
        plt.axvline(x=ll_a, c='k', linestyle=':')

        plt.axhline(y=n_A, c='k', linestyle='--')
        plt.axhline(y=n_B, c='k', linestyle='--')
        plt.legend()
    
    x_A = n_A / (n_B + n_A)
    return x_A

def conc_profile(Cs, Bs, label, species=[1., 2.], block_time=2, dt=1, save_rate = 0.06, avg_thresh=0.2, c=0):
    T0 = 0
    X = np.array([]); t= np.array([])
    for k, C in enumerate(Cs):
        frames_per_block = int(block_time / save_rate)
        tot_steps = np.floor(len(C) * save_rate / 10) * 10
        nblocks = int((tot_steps/save_rate) / frames_per_block)
        Xs = []
        for i in range(nblocks):
            x = density_profile(C, Bs[k], i*frames_per_block, (i+1)*frames_per_block, species=species)
            Xs.append(x)
        T = np.arange(T0, T0 + nblocks, 1)
        T0 = T0 + nblocks
        t = np.hstack((t, T*block_time/dt)); X = np.hstack((X, Xs))

    plt.plot(t*dt, X, '--', linewidth=2, label=label)
    plt.xlabel('time (ns)')
    plt.ylabel(r'$X_A$')

    I = int(len(X)*avg_thresh)
    xs = np.mean(X[I:])
    _, err, _, _ = block_average(X[I:], N=40, plot_flag=False)

    R = Rectangle([-10, xs-err], 660, 2*err, alpha=0.4, color='C' + str(c))
    ax = plt.gca()
    ax.add_patch(R)
    
    print("SOLUBILITY = " + str(xs) + u" \u00B1 " + str(err))
    return xs, err