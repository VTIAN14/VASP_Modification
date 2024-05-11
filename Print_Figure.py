#!/bin/env python3
import linecache
import os
import Parsing_VASP
import Modifing_VASP
import Printing_VASP


linecache.clearcache()

###############################################################################

read_OUT_name = 'OUTCAR'
f_or_nf,imagin,cm,mev,modes_3d,zpe,Energy0,force,E_sigma_0,E_diff = Parsing_VASP.Parsing_OUTCAR_Result(read_OUT_name)

###############################################################################

Printing_VASP.Print_Figures(f_or_nf,Energy0,force,E_sigma_0,E_diff)

###############################################################################
    
linecache.clearcache()