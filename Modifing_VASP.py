#!/bin/env python3
def Fix_Type_POSCAR(fix_type,atom_type,atom_num_tot,atom_num):
    l = 0
    atoms_relax1 = [[0 for _ in range(3)] for _ in range(atom_num_tot)] # creat 2-d 0 array   
    for i in range(len(atom_type)):
        for j in fix_type:
            if (j == atom_type[i]):
                for k in range(atom_num[i]):
                    atoms_relax1[k+l] = [-1,-1,-1]
        l = l + atom_num[i]
    return(atoms_relax1)

def Gen_Type_Show_POSCAR(atoms_relax,atom_type,atoms_coor,atom_num):
    atom_type = ['O','N']
    fixed = 0
    relaxed = 0
    atoms_coor_new = []
    for i in range(len(atoms_relax)):
        if -1 in atoms_relax[i]:            
            fixed = fixed + 1
            atoms_coor_new.append(atoms_coor[i])
    for i in range(len(atoms_relax)):
        if [0,0,0] == atoms_relax[i]:            
            relaxed = relaxed + 1
            atoms_coor_new.append(atoms_coor[i])
    atom_num = [fixed,relaxed]
    return(atom_type,atom_num,atoms_coor_new)

def Fix_Height_POSCAR(direct_ban,lattice_con,D_or_C,atom_num_tot,atoms_coor):
# direct_ban only for direct coordinate
    if D_or_C < 0: # C2D
        D_or_C_D2C,atoms_coor_D2C = D2C_POSCAR(lattice_con,atom_num_tot,D_or_C,atoms_coor)
        # use function D2C_POSCAR
        atoms_coor = atoms_coor_D2C        
    atoms_relax1 = [[0 for _ in range(3)] for _ in range(atom_num_tot)]
    for i in range(atom_num_tot):
        if atoms_coor[i][2] < direct_ban:
            atoms_relax1[i] = [-1,-1,-1]
    return(atoms_relax1)

def D2C_POSCAR(lattice_con,atom_num_tot,D_or_C,atoms_coor):
    import numpy as np
    atoms_coor1 = [[0 for _ in range(3)] for _ in range(atom_num_tot)]
    if D_or_C >= 0: # D2C
        D_or_C1 = -1
        for i in range(atom_num_tot):
            atoms_coor1[i][0] = atoms_coor[i][0]*lattice_con[0][0] + \
                atoms_coor[i][1]*lattice_con[1][0] + atoms_coor[i][2]*lattice_con[2][0]
            atoms_coor1[i][1] = atoms_coor[i][0]*lattice_con[0][1] + \
                atoms_coor[i][1]*lattice_con[1][1] + atoms_coor[i][2]*lattice_con[2][1]
            atoms_coor1[i][2] = atoms_coor[i][0]*lattice_con[0][2] + \
                atoms_coor[i][1]*lattice_con[1][2] + atoms_coor[i][2]*lattice_con[2][2]
    else:       # C2D
        D_or_C1 = 0
        for i in range(atom_num_tot):
            matrix_lc = np.array([[lattice_con[0][0],lattice_con[1][0],lattice_con[2][0]],\
                          [lattice_con[0][1],lattice_con[1][1],lattice_con[2][1]],\
                          [lattice_con[0][2],lattice_con[1][2],lattice_con[2][2]]])
            matrix_coor = np.array([[atoms_coor[i][0]],[atoms_coor[i][1]],[atoms_coor[i][2]]])
            x = np.linalg.solve(matrix_lc,matrix_coor)
            atoms_coor1[i][0] = x[0,0]
            atoms_coor1[i][1] = x[1,0]
            atoms_coor1[i][2] = x[2,0]
    return(D_or_C1, atoms_coor1)

def Surface_Perturbation(weight,which_fre,lattice_con,atom_num_tot,D_or_C,atoms_coor,modes_3d):
# num_fre tells which fre is used to apply perturbate the surface: -1 the last one, -2 the second to last...
    D_or_C1 = -1
    D_or_C2 = -1
    atoms_coor1 = [[0 for _ in range(3)] for _ in range(atom_num_tot)]
    atoms_coor2 = [[0 for _ in range(3)] for _ in range(atom_num_tot)]
    if D_or_C >= 0: # D2C
        D_or_C_D2C,atoms_coor_D2C = D2C_POSCAR(lattice_con,atom_num_tot,D_or_C,atoms_coor)
        # use function D2C_POSCAR
    else:
        atoms_coor_D2C = atoms_coor
    for i in range(atom_num_tot):
        atoms_coor1[i][0] = atoms_coor_D2C[i][0] + weight*modes_3d[which_fre][i][0]
        atoms_coor1[i][1] = atoms_coor_D2C[i][1] + weight*modes_3d[which_fre][i][1]
        atoms_coor1[i][2] = atoms_coor_D2C[i][2] + weight*modes_3d[which_fre][i][2]
    for i in range(atom_num_tot):
        atoms_coor2[i][0] = atoms_coor_D2C[i][0] - weight*modes_3d[which_fre][i][0]
        atoms_coor2[i][1] = atoms_coor_D2C[i][1] - weight*modes_3d[which_fre][i][1]
        atoms_coor2[i][2] = atoms_coor_D2C[i][2] - weight*modes_3d[which_fre][i][2]
    return(D_or_C1,D_or_C2,atoms_coor1,atoms_coor2)

def Atoms_Remove(rm_index,lattice_con,atom_type,atom_num,atom_num_tot,\
        D_or_C,atoms_coor,atoms_relax):
# add_coord only for direct coordinate
    import copy
    D_or_C1 = -1
    rm_index_tem = copy.deepcopy(rm_index)
    rm_index = []
    atoms_coor1 = []
    atom_type_tem = copy.deepcopy(atom_type)
    atom_type1 = []
    atom_num_tem = copy.deepcopy(atom_num)
    atom_num1 = []
    atoms_relax1 = []
    atom_num_tot1 = atom_num_tot
    if D_or_C >= 0: # D2C
        D_or_C_D2C,atoms_coor_D2C = D2C_POSCAR(lattice_con,atom_num_tot,D_or_C,atoms_coor)
    else:
        atoms_coor_D2C = atoms_coor
########################################### delete index (> tot_atom_num)
    for i in range(len(rm_index_tem)):
        if rm_index_tem[i] <= 0 or rm_index_tem[i] > atom_num_tot:
            with open('Error_message','a') as ff:
                ff.write('The ' + str(rm_index_tem[i]) + 'th assigned atoms are not found in POSCAR.\n')
            rm_index_tem[i] = 0
    for item in rm_index_tem:
        if item != 0:
            rm_index.append(item)
    rm_index_tem = copy.deepcopy(rm_index)
###########################################
    which_type = [0 for _ in range(len(rm_index))]
    for i in range(len(rm_index_tem)):        
        atom_num_tot1 = atom_num_tot1 - 1
        for j in range(len(atom_num)):
            rm_index_tem[i] = rm_index_tem[i] - atom_num_tem[j]
            if rm_index_tem[i] <= 0:
                which_type[i] = atom_type_tem[j]
                break
    for i in which_type:
        atom_num_tem[atom_type_tem.index(i)] = atom_num_tem[atom_type_tem.index(i)] - 1
    for i in range(len(atom_num_tem)):
        if atom_num_tem[i] == 0:
            atom_type_tem[i] = ''
    for item in atom_num_tem:
        if item != 0:
            atom_num1.append(item)   
    for item in atom_type_tem:
        if item != '':
            atom_type1.append(item)
    for i in range(len(atoms_coor_D2C)):
        if i+1 not in rm_index:
            atoms_coor1.append(atoms_coor_D2C[i])
    for i in range(len(atoms_relax)):
        if i+1 not in rm_index:
            atoms_relax1.append(atoms_relax[i])
    return(D_or_C1,atoms_coor1,atom_type1,atom_num1,atom_num_tot1,atoms_relax1)    
