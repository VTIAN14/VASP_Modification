#!/bin/env python3
import linecache
import os
import Parsing_VASP
import Modifing_VASP
import Printing_VASP


linecache.clearcache()

###############################################################################

os.rename('POSCAR','POSCARO')

###############################################################################
    
read_POS_name = 'POSCARO'
comment,sc_factor,lattice_con,atom_type,atom_num,atom_num_tot,sel_or_not,\
       D_or_C,atoms_coor,atoms_relax = Parsing_VASP.Parsing_POSCAR(read_POS_name)

###############################################################################

rm_index = []
D_or_C1,atoms_coor1,atom_type1,atom_num1,atom_num_tot1,atoms_relax1 \
    = Modifing_VASP.Atoms_Remove(rm_index,lattice_con,atom_type,atom_num,atom_num_tot,D_or_C,atoms_coor,atoms_relax)

###############################################################################

gen_POS_name = 'POSCAR'
Printing_VASP.Print_POSCAR(gen_POS_name,comment,sc_factor,lattice_con,atom_type1,atom_num1,\
               atom_num_tot1,sel_or_not,D_or_C1,atoms_coor1,atoms_relax1)

###############################################################################
    
linecache.clearcache()
