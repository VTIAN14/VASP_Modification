#!/bin/env python3
import linecache
import Parsing_VASP
import Modifing_VASP
import Printing_VASP


linecache.clearcache()

###############################################################################

read_OUT_name = 'OUTCAR'
f_or_nf,imagin,cm,mev,modes_3d,zpe,Energy0,force,E_sigma_0,E_diff = Parsing_VASP.Parsing_OUTCAR_Result(read_OUT_name)

read_OUT_name = 'OUTCAR'
nkps = Parsing_VASP.Parsing_OUTCAR_Kpoints(read_OUT_name)

###############################################################################

Printing_VASP.Print_Figures(f_or_nf,Energy0,force,E_sigma_0,E_diff)

read_POS_name = 'CONTCAR'
comment,sc_factor,lattice_con,atom_type,atom_num,atom_num_tot,sel_or_not,\
        D_or_C,atoms_coor,atoms_relax = Parsing_VASP.Parsing_POSCAR(read_POS_name)

atom_type,atom_num,atoms_coor = Modifing_VASP.Gen_Type_Show_POSCAR(atoms_relax,atom_type,atoms_coor,atom_num)

gen_POS_name = 'CONTCAR_TF'
Printing_VASP.Print_POSCAR(gen_POS_name,comment,sc_factor,lattice_con,atom_type,atom_num,atom_num_tot,sel_or_not,D_or_C,atoms_coor,atoms_relax)

Printing_VASP.Print_CONTCAR_Conf()

read_OUT_name = 'OUTCAR'
c_or_nc = Printing_VASP.Print_Check_Keywords(read_OUT_name,nkps)

Printing_VASP.Print_Summary_pdf(f_or_nf,c_or_nc,Energy0,force,E_sigma_0,E_diff)

###############################################################################
    
linecache.clearcache()
