import argparse
import numpy as np
import pandas as pd

esu_to_e = 1

MM = {
    'O': 0.01528,
    'H': 6.0,
    'Na': 22.989769,
    'Cl': 35.453,
    'NaI': 22.989769,
    'ClI': 35.453,
    'X': 6.0,
    'DP': 0.2,
    'DPs': 4.0,
    'Ar': 39.948,
}

q = {
    'O': 0.0,
    'H': 0.584 * esu_to_e,
    'Na': 0.0,
    'Cl': 0.0,
    'Nai': 0.0,
    'Cli': 0.0,
    'X': -1.168 * esu_to_e,
    'A' : 0,
    'B' : 0,
    'C' : 0,
}


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

def generate_data_lammps_polarizable(p, o, m, dp, s, b, a, mol, sdm=None, einstein=False):
    ri = []
    rf = []
    e, c = read_xyz(p)

    try:
        with open(m, 'r') as d:
            lines = d.readlines()
            for k, line in enumerate(lines):
                if line.find('inside box') != -1:
                    for i, z in enumerate(line.split()):
                        if 1 < i < 5:
                            ri.append(float(z))
                        elif i >= 5:
                            rf.append(float(z))
                    break
    except:
            print("WARNING: MIXTURE FILE DOESN'T EXIST")
        
    try:
        with open("box.xyz", 'r') as B:
            lines = B.readlines()
            for k, line in enumerate(lines):
                arr = np.array(line.split(), dtype=np.float64)
                ri.append(arr[0]); rf.append(arr[1])
    except:
        ri.append("replace1"); ri.append("replace1"); ri.append("replace1")
        rf.append("replace2"); rf.append("replace2"); rf.append("replace2")
    
    unique_ele = pd.unique(e)
    unique_ele = unique_ele.astype("object")
    E = unique_ele; count = 0
    drude_cores = np.array(dp.split())
    for i, ele in enumerate(unique_ele):
        if ele in drude_cores:
            E = np.insert(E, count+i+1, str('D' + ele))
            count+=1
    
    n_atoms = 0
    for e_ in e:
        if e_ in drude_cores:
            n_atoms += 2
        else:
            n_atoms += 1
  
    bonds = []
    b = b.split()
    for i in range(0, len(b)-1, 2):
        bonds.append((b[i], b[i+1]))
  
    for i, ele in enumerate(E):
        if ele in drude_cores:
            bonds.append((E[i], E[i+1]))
            
    angles = []
    a = a.split()
    for i in range(0, len(a)-2, 3):
        angles.append((a[i], a[i+1], a[i+2]))
        
    n_bonds = 0
    for bond in bonds:
        rep_ele = bond[0]
        instances = np.where(e == rep_ele)[0]
        if bond[0] == "O" and bond[1] == "H":
            n_bonds += 2*len(instances)
        else:
            n_bonds += len(instances)
        
    n_angles = 0
    for angle in angles:
        rep_ele = angle[0]
        instances = np.where(e == rep_ele)[0]
        n_angles += len(instances)

    with open(o, 'w') as f:
        f.write('LAMMPS datafile' + '\n')
        f.write('\n')

        f.write(str(n_atoms) + ' atoms\n')
        f.write(str(len(E)) + ' atom types\n')
        f.write(str(n_bonds) + ' bonds\n')
        f.write(str(len(bonds)) + ' bond types\n')
        f.write(str(n_angles) + ' angles\n')
        f.write(str(len(angles)) + ' angle types\n')
        f.write('\n')
        
        f.write(str(ri[0]) + ' ' + str(rf[0]) + '  xlo xhi\n')
        f.write(str(ri[1]) + ' ' + str(rf[1]) + '  ylo yhi\n')
        f.write(str(ri[2]) + ' ' + str(rf[2]) + '  zlo zhi\n')
        f.write('\n')

        f.write('Masses\n')
        f.write('\n')

        for i, ele in enumerate(E):
            
            if ele[0] == 'D':
                if ele[1:] in s.split():
                    if sdm == None:
                        molar_mass = MM['DPs']
                    else:
                        molar_mass = sdm
                else:
                    molar_mass = MM['DP']
            elif ele in drude_cores:
                try:
                    if ele in s.split():
                        if sdm == None:
                            molar_mass = MM[ele] - MM['DPs']
                        else:
                            molar_mass = MM[ele] - sdm
                    else:
                        molar_mass = MM[ele] - MM['DP']
                except:
                    if ele in s.split():
                        molar_mass = MM['Ar'] - MM['DP']
                    else:
                        molar_mass = MM['Ar'] - MM['DP']
            else:
                try:
                    molar_mass = MM[ele]
                except:
                    molar_mass = MM['Ar']
           
            f.write(str(i+1) + ' ' + str(molar_mass) + '\n')
        f.write('\n')

        f.write('Atoms\n')
        f.write('\n')
        
        mol_id = 0
        tot_bonds = []
        tot_angles = []
        atoms_added = 0
        e_ = []
        for j, coord in enumerate(c):
                    
            atom_id = np.where(E == e[j])[0][0] + 1
            if E[atom_id-1] in mol:
                mol_id += 1

            charge = 0.0
            
            dr = np.random.random(3)*0.05 
            if e[j] in drude_cores:
                drude_id = np.where(E == 'D' + e[j])[0][0] + 1
                f.write(str(j+1+atoms_added) + ' ' + str(mol_id) + ' ' + str(atom_id) + ' ' + str(charge) + ' ' + str(coord[0])
                    + ' ' + str(coord[1]) + ' ' + str(coord[2]) + '\n')
                f.write(str(j+2+atoms_added) + ' ' + str(mol_id) + ' ' + str(drude_id) + ' ' + str(charge) + ' ' 
                        + str(coord[0]+ dr[0]) + ' ' + str(coord[1] + dr[1]) + ' ' + str(coord[2] + dr[2]) + '\n')
                e_.append(e[j]); e_.append('D' + e[j])
                atoms_added += 1
            else:
                f.write(str(j+1+atoms_added) + ' ' + str(mol_id) + ' ' + str(atom_id) + ' ' + str(charge) + ' ' + str(coord[0])
                    + ' ' + str(coord[1]) + ' ' + str(coord[2]) + '\n')
                e_.append(e[j])
        f.write('\n')
        
        if n_bonds > 0:
            f.write('Bonds\n')
            f.write('\n')
            count=1
            for i, elei in enumerate(e_):
                for k, bond in enumerate(bonds):
                    if elei != bond[0]:
                        continue
                    bond_id = k+1
                    for j, elej in enumerate(e_[i:]):
                        if elej == bond[1]:
                            f.write(str(count) + " " + str(bond_id) + " " + str(i+1) + " " + str(i + j + 1) + '\n')
                            count += 1
                            if bond[0] == "O" and bond[1]=="H":
                                f.write(str(count) + " " + str(bond_id) + " " + str(i+1) + " " + str(i + j + 3) + '\n')
                                count += 1
                            break
            f.write('\n')
        
        if n_angles > 0:
            f.write('Angles\n')
            f.write('\n')
            count=1
            for i, elei in enumerate(e_):
                for k, angle in enumerate(angles):
                    if elei != angle[0]:
                        continue
                        
                    angle_id = k+1
                    id1 = 0; id2 = 0
                    for j, elej in enumerate(e_[i:]):
                        if elej == angle[1]:
                            id1 = j
                            break
                            
                    for j, coordj in enumerate(e_[(i+id1+1):]):
                        if e[j] == angle[2]:
                            id2 = j
                            f.write(str(count) + " " + str(angle_id) + " " + str(i+1) + " " + str(i + id1 + 1) + " " + str(i + id1 + id2 + 1) + '\n')
                            count += 1
                            break
                            
            f.write('\n')

if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument("-i", "--xyz", type=str, help="input xyz file")
    argParser.add_argument("-m", "--mixture", type=str, help="input packmol mixture file")
    argParser.add_argument("-o", "--out", type=str, help="output LAMMPS data file")
    argParser.add_argument("-d", "--drude", type=str, help="list of atoms that are drude cores")
    argParser.add_argument("-s", "--sym", type=str, help="list of atoms that have symmetrized interactions")
    argParser.add_argument("-b", "--bonds", type=str, help="list of bonds")
    argParser.add_argument("-a", "--angles", type=str, help="list of angles")
    argParser.add_argument("-sdm", "--sym_mass", type=float, help="symmetrized drude mass")
    argParser.add_argument("-mol", "--molecules", type=str, help="list of atoms characteristic of distinct molecules")
    args = argParser.parse_args()

    generate_data_lammps_polarizable(args.xyz, args.out, args.mixture, args.drude, args.sym, args.bonds, args.angles, args.molecules, sdm=args.sym_mass)
