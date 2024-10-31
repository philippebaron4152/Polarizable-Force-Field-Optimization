import argparse
import numpy as np
import pandas as pd

K = 1

MM = {
    'Ar': 39.948,
    'O': 15.9994,
    'H': 1.00784,
    'Na': 22.989769,
    'Cl': 35.453,
    'Nai': 22.989769,
    'Cli': 35.453,
}

q = {
    'O': -0.8476,
    'H': 0.4238,
    'Na': K*1.0000,
    'Cl': K*-1.0000,
    'Nai': K*1.0000,
    'Cli': K*-1.0000,
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


def generate_data_lammps(p, o, m, b="", a="", einstein=False):
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
                
    bonds = []
    b = b.split()
    for i in range(0, len(b)-1, 2):
        bonds.append((b[i], b[i+1]))
    
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

    unique_ele = pd.unique(e)
    with open(o, 'w') as f:
        f.write('LAMMPS datafile' + '\n')
        f.write('\n')

        f.write(str(len(e)) + ' atoms\n')
        f.write(str(len(unique_ele)) + ' atom types\n')
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
        for i, ele in enumerate(unique_ele):
            try:
                molar_mass = MM[ele]
            except:
                at_per_mol = 1
                for bond in bonds:
                    if ele in bond:
                        at_per_mol += 1
                molar_mass = MM['Ar'] / at_per_mol
            f.write(str(i+1) + ' ' + str(molar_mass) + '\n')
        f.write('\n')

        f.write('Atoms\n')
        f.write('\n')
        
        mol_id = 0
        tot_bonds = []
        tot_angles = []
        for j, coord in enumerate(c):
            
            if n_bonds > 0:
                for i, bond in enumerate(bonds):
                    if e[j] == bond[0]:
                        mol_id += 1
            else:
                mol_id = j + 1
                    
            atom_id = np.where(unique_ele == e[j])[0][0] + 1

            try:
                charge = q[e[j]]
            except:
                charge = 0.0
                
            f.write(str(j+1) + ' ' + str(mol_id) + ' ' + str(atom_id) + ' ' + str(charge) + ' ' + str(coord[0])
                    + ' ' + str(coord[1]) + ' ' + str(coord[2]) + '\n')
        f.write('\n')
        
        if n_bonds > 0:
            f.write('Bonds\n')
            f.write('\n')
            count=1
            for i, elei in enumerate(e):
                for k, bond in enumerate(bonds):
                    if elei != bond[0]:
                        continue
                    bond_id = k+1
                    for j, elej in enumerate(e[i:]):
                        if elej == bond[1]:
                            f.write(str(count) + " " + str(bond_id) + " " + str(i+1) + " " + str(i + j + 1) + '\n')
                            count += 1
                            if bond[0] == "O" and bond[1]=="H":
                                f.write(str(count) + " " + str(bond_id) + " " + str(i+1) + " " + str(i + j + 2) + '\n')
                                count += 1
                            break
            f.write('\n')
        
        if n_angles > 0:
            f.write('Angles\n')
            f.write('\n')
            count=1
            for i, elei in enumerate(e):
                for k, angle in enumerate(angles):
                    if elei != angle[0]:
                        continue
                        
                    angle_id = k+1
                    id1 = 0; id2 = 0
                    for j, elej in enumerate(e[i:]):
                        if elej == angle[1]:
                            id1 = j
                            break
                            
                    for j, coordj in enumerate(e[(i+id1+1):]):
                        if e[j] == angle[2]:
                            id2 = j
                            f.write(str(count) + " " + str(angle_id) + " " + str(i+id1+1) + " " + str(i + 1) + " " + str(i + id1 + id2 + 1) + '\n')
                            count += 1
                            break
                            
            f.write('\n')

if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument("-i", "--xyz", type=str, help="input xyz file")
    argParser.add_argument("-m", "--mixture", type=str, help="input packmol mixture file")
    argParser.add_argument("-o", "--out", type=str, help="output LAMMPS data file")
    argParser.add_argument("-b", "--bonds", type=str, help="pairs of bonded atoms")
    argParser.add_argument("-a", "--angles", type=str, help="atoms present in molecular angles")
    argParser.add_argument('--einstein', action=argparse.BooleanOptionalAction)
    args = argParser.parse_args()
    generate_data_lammps(args.xyz, args.out, a=args.angles, m=args.mixture, b=args.bonds, einstein=args.einstein)
