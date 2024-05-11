#!/bin/env python3
import linecache

def Print_POSCAR(gen_POS_name,comment,sc_factor,lattice_con,atom_type,atom_num,\
               atom_num_tot,sel_or_not,D_or_C,atoms_coor,atoms_relax):
    linecache.clearcache()
    atoms_relax_TF = [[0 for _ in range(3)] for _ in range(atom_num_tot)]
    for i in range(atom_num_tot):
        for j in range(3):
            if atoms_relax[i][j] < 0:
                atoms_relax_TF[i][j] = 'F'
            else:
                atoms_relax_TF[i][j] = 'T'   
    with open(gen_POS_name,'w') as ff:
        ff.write(comment + '\n')
    with open(gen_POS_name,'a') as ff:
        ff.write(str(sc_factor) + '\n')
        ff.write(' '.join([str(x) for x in lattice_con[0]]) + '\n')
        ff.write(' '.join([str(x) for x in lattice_con[1]]) + '\n')
        ff.write(' '.join([str(x) for x in lattice_con[2]]) + '\n')
        ff.write(' '.join(atom_type) + ' \n')
        ff.write(' '.join([str(x) for x in atom_num]) + '\n')
        ff.write('Selective Dynamics\n')
        if D_or_C >=0:
            ff.write('Direct\n')
        else:
            ff.write('Cartesian\n')
    with open(gen_POS_name,'a') as ff:
        for i in range(atom_num_tot):
            ff.write(' '.join([str(x) for x in atoms_coor[i]]) + ' ' +\
                     ' '.join(atoms_relax_TF[i]) + '\n')
    linecache.clearcache()
    
def Print_MODECAR(which_fre,f_or_nf,modes_3d):
# generate MODECAR from *fre-calculation* OUTCAR, default choice is the *last* mode 
# num_fre tells which fre is used to apply perturbate the surface: -1 the last one, -2 the second to last...
    linecache.clearcache()
    if f_or_nf < 0:
        with open('Error_message','w') as ff:
            ff.write('Frequency test does not finished!')
    with open('MODECAR','w') as ff:
        ff.write(' '.join(str(x) for x in modes_3d[which_fre][0]) + '\n')
    with open('MODECAR','a') as ff:
        for i in range(len(modes_3d[-1])-1):
            ff.write(' '.join(str(x) for x in modes_3d[which_fre][i+1]) + '\n')
        modes_3d[-1]
    linecache.clearcache()
    
def Print_Figures(f_or_nf,Energy0,force,E_sigma_0,E_diff):
    import matplotlib.pyplot as plt
    linecache.clearcache()
    if f_or_nf < 0:
        with open('Error_message','w') as ff:
            ff.write('Calculation does not finished!')
    else:
        with open('Complete!','w') as ff:
            ff.write(str(E_sigma_0[-1]) + '\n')
            ff.write(str(force[-1]) + '\n')            
    with open('Energy_Data','w') as ff:
        for i in range(len(E_sigma_0)):
            ff.write(str(E_sigma_0[i]) + ' \n')
    with open('Force_Data','w') as ff:
        for i in range(len(force)):
            ff.write(str(force[i]) + ' \n')    
    with open('Energy_Diff_Data','w') as ff:
        for i in range(len(E_diff)):
            ff.write(str(E_diff[i]) + ' \n')            
    step = []
    for i in range(len(force)):
        step.append(i)
    plt.plot(step, force)
    plt.xlabel('Step')
    plt.ylabel('Force (eV/Ang)')
    plt.tight_layout()
    plt.savefig('./'+'Force')
    plt.cla()    
###############################################################################
    plt.plot(step, E_sigma_0) 
    plt.xlabel('Step')
    plt.ylabel('Energy (eV)')
    plt.tight_layout()
    plt.savefig('./'+'E_sigma_0')
    plt.cla()
###############################################################################
    del step[0]       
    plt.plot(step, E_diff)
    plt.xlabel('Step')
    plt.ylabel('Energydiff (eV)')
    plt.tight_layout()
    plt.savefig('./'+'Energydiff')
    linecache.clearcache()
