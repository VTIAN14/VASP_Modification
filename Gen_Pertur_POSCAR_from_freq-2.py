#!/bin/env python3
import linecache
import os
import Parsing_VASP
import Modifing_VASP
import Printing_VASP


linecache.clearcache()

###############################################################################
    
read_POS_name = 'POSCAR'
comment,sc_factor,lattice_con,atom_type,atom_num,atom_num_tot,sel_or_not,\
       D_or_C,atoms_coor,atoms_relax = Parsing_VASP.Parsing_POSCAR(read_POS_name)

###############################################################################

# Generate from POSCAR for freq calculation usually fix more atoms than we want
# This block aims to relax the top 2-3 layers atoms of the bulk
direct_ban = 0.20
atoms_relax1 = Modifing_VASP.Fix_Height_POSCAR(direct_ban,lattice_con,D_or_C,atom_num_tot,atoms_coor)

###############################################################################

read_OUT_name = 'OUTCAR'
f_or_nf,imagin,cm,mev,modes_3d,zpe,Energy0,force,E_sigma_0,E_diff = Parsing_VASP.Parsing_OUTCAR_Result(read_OUT_name)

###############################################################################

weight = 0.3
which_fre = -2
D_or_C1,D_or_C2,atoms_coor1,atoms_coor2 = Modifing_VASP.Surface_Perturbation(weight,which_fre,\
                                                             lattice_con,atom_num_tot,D_or_C,atoms_coor,modes_3d)  

###############################################################################

gen_POS_name = 'POSCAR1-2'    
Printing_VASP.Print_POSCAR(gen_POS_name,comment,sc_factor,lattice_con,atom_type,atom_num,\
               atom_num_tot,sel_or_not,D_or_C1,atoms_coor1,atoms_relax1)

gen_POS_name = 'POSCAR2-2'
Printing_VASP.Print_POSCAR(gen_POS_name,comment,sc_factor,lattice_con,atom_type,atom_num,\
               atom_num_tot,sel_or_not,D_or_C2,atoms_coor2,atoms_relax1)

###############################################################################
    
linecache.clearcache()
