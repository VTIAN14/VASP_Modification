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

fix_type = ['Ag','Au','Cu','Ni','Pd','Pt','Rh']
atoms_relax1 = Modifing_VASP.Fix_Type_POSCAR(fix_type,atom_type,atom_num_tot,atom_num)

###############################################################################

gen_POS_name = 'POSCAR'
Printing_VASP.Print_POSCAR(gen_POS_name,comment,sc_factor,lattice_con,atom_type,atom_num,\
               atom_num_tot,sel_or_not,D_or_C,atoms_coor,atoms_relax1)

###############################################################################
    
linecache.clearcache()
