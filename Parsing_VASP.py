#!/bin/env python3
import linecache
import os

def Parsing_POSCAR(read_POS_name):
    
###############################################################################    
    linecache.clearcache()
    proc_message = 'proc_message'  
###############################################################################
    
###############################################################################
    with open(read_POS_name,'r') as ff:
        a = linecache.getline(read_POS_name,1)
    comment = a[:-1]# this comment involves "\n"
###############################################################################
    
###############################################################################
    with open(read_POS_name,'r') as ff:
        a = linecache.getline(read_POS_name,2)
    b = a.replace(' ','')
    a = b.replace('\n','')
    sc_factor = a
###############################################################################
    
###############################################################################
    lattice_con = [[],[],[]]
    for i in range(3):
        with open(read_POS_name,'r') as ff:
            a = linecache.getline(read_POS_name,i+3)
        b = a.replace(' ','\n')
        with open(proc_message,'w') as ff:
            ff.write(b)
        with open(proc_message,'r') as ff:
            a = ff.readlines()
        for j in range(len(a)):
            b = a[j].replace('\n','')
            a[j] = b
        while '' in a:
            a.remove('')
        for j in range(len(a)):
            b = eval(a[j])
            a[j] = b
        lattice_con[i] = a
###############################################################################
    
###############################################################################
    with open(read_POS_name,'r') as ff:
        a = linecache.getline(read_POS_name,6)
    b = a.replace(' ','\n')
    with open(proc_message,'w') as ff:
        ff.write(b)
    with open(proc_message,'r') as ff:
        a = ff.readlines()
    for i in range(len(a)):
        b = a[i].replace('\n','')
        a[i] = b
    while '' in a:
        a.remove('')
    atom_type = a
###############################################################################
    
###############################################################################
    with open(read_POS_name,'r') as ff:
        a = linecache.getline(read_POS_name,7)
    b = a.replace(' ','\n')
    with open(proc_message,'w') as ff:
        ff.write(b)
    with open(proc_message,'r') as ff:
        a = ff.readlines()
    for i in range(len(a)):
        b = a[i].replace('\n','')
        a[i] = b
    while '' in a:
        a.remove('')
    for i in range(len(a)):
        b = eval(a[i])
        a[i] = b
    atom_num = a # how many atoms for each type of elements
    atom_num_tot = sum(atom_num)
###############################################################################
    
###############################################################################
    with open(read_POS_name,'r') as ff:
        a = linecache.getline(read_POS_name,8)
        sel_or_not= a.find('S') # if we have Selective dynamics. 0:yes, -1:no
###############################################################################
    
###############################################################################
    with open(read_POS_name,'r') as ff:
        a = linecache.getline(read_POS_name,9+sel_or_not)
        D_or_C = a.find('D') # if it's Direct coordinate. 0.Direct, -1:Cartesian
###############################################################################
    
###############################################################################
    atoms_coor=[]
    a2 = [0,0,0]
    for i in range(atom_num_tot):
        with open(read_POS_name,'r') as ff:
            a = linecache.getline(read_POS_name,10+i+sel_or_not)
        b = a.replace(' ','\n')
        with open(proc_message,'w') as ff:
            ff.write(b)
        with open(proc_message,'r') as ff:
            a = ff.readlines()
        for j in range(len(a)):
            b = a[j].replace('\n','')
            a[j] = b
        while '' in a:
            a.remove('')
        for j in range(3):
            b = eval(a[j])
            a2[j] = b
        one_atom_coor = a2
        one_atom_coor2 = one_atom_coor.copy()
        atoms_coor.append(one_atom_coor2)
###############################################################################
    
###############################################################################
    atoms_relax=[]
    a2 = [0,0,0]
    for i in range(atom_num_tot):
        with open(read_POS_name,'r') as ff:
            a = linecache.getline(read_POS_name,10+i+sel_or_not)
        b = a.replace(' ','\n')
        with open(proc_message,'w') as ff:
            ff.write(b)
        with open(proc_message,'r') as ff:
            a = ff.readlines()
        for j in range(len(a)):
            b = a[j].replace('\n','')
            a[j] = b
        while '' in a:
            a.remove('')
        if len(a) == 3:
            a2 = [0,0,0]
        else:
            for j in range(3,6):
                b = a[j].find("T")
                a2[j-3] = b
        one_atom_relax = a2
        one_atom_relax2 = one_atom_relax.copy()
        atoms_relax.append(one_atom_relax2) # fixed:-1, relaxed:0
###############################################################################

###############################################################################
    os.remove(proc_message)
############################################################################### 
   
    return(comment,sc_factor,lattice_con,atom_type,atom_num,atom_num_tot,sel_or_not,\
           D_or_C,atoms_coor,atoms_relax)

def Parsing_OUTCAR_Result(read_OUT_name):

###############################################################################    
    linecache.clearcache()
    proc_message = 'proc_message'
    imagin = []
    cm = []
    mev = []
    mode_2d = [] # mode in one frequecy: 2D array
    modes_3d = [] # modes for all frequencies: 3D array
    force = []
    E_sigma_0 = []
    E_diff = []
    f_or_nf = -1
    all_lines = len(open(read_OUT_name).readlines())
###############################################################################
# stop criteria   
###############################################################################
    with open(read_OUT_name,'r') as ff:
        i = 1
        while True:
            a = linecache.getline(read_OUT_name,i)
            #print(i) # for debug purpose
            if a.find("Voluntary context switches") >=0:
                break
            if i > all_lines:
                break
###############################################################################
# reading frequency mode message
###############################################################################
            if a.find('f  =') >=0 or a.find('f/i=') >=0:
                b = a.replace(' ','\n')
                with open(proc_message,'w') as ff:
                    ff.write(b)
                with open(proc_message,'r') as ff:
                    a = ff.readlines()
                for j in range(len(a)):
                    b = a[j].replace('\n','')
                    a[j] = b
                while '' in a:
                    a.remove('')
                im_or_re = a[1].find('i')
                imagin.append(im_or_re) # if it is imaginry:2 real:-1
                if im_or_re >= 0:
                    cm.append(eval(a[6]))
                    mev.append(eval(a[8]))
                # interface here! if you want more parameters...
                if im_or_re <= 0:
                    cm.append(eval(a[7]))
                    mev.append(eval(a[9]))
                # interface here! if you want more parameters... 
                i = i + 2
                k = 0
                while True:
                    a = linecache.getline(read_OUT_name,i)
                    b = a.replace(' ','\n')
                    with open(proc_message,'w') as ff:
                        ff.write(b)
                    with open(proc_message,'r') as ff:
                        a = ff.readlines()
                    for j in range(len(a)):
                        b = a[j].replace('\n','')
                        a[j] = b
                    while '' in a:
                        a.remove('')
                    if (len(a) <= 1):
                        modes_3d.append(mode_2d)
                        mode_2d = []
                        break
                    for j in range(3,6):
                        b = eval(a[j])
                        a[j]=b
                    mode_2d.append(a[3:6])
                    i = i + 1
                    k = k + 1
            a = linecache.getline(read_OUT_name,i)
###############################################################################
# reading force/energy message
###############################################################################
            if a.find('FORCES:') >=0:
                b = a.replace(' ','\n')
                with open(proc_message,'w') as ff:
                    ff.write(b)
                with open(proc_message,'r') as ff:
                    a = ff.readlines()
                for j in range(len(a)):
                    b = a[j].replace('\n','')
                    a[j] = b
                while '' in a:
                    a.remove('')
                force.append(eval(a[4]))
            a = linecache.getline(read_OUT_name,i)
###############################################################################
            if a.find('energy  without entropy') >=0:
                b = a.replace(' ','\n')
                with open(proc_message,'w') as ff:
                    ff.write(b)
                with open(proc_message,'r') as ff:
                    a = ff.readlines()
                for j in range(len(a)):
                    b = a[j].replace('\n','')
                    a[j] = b
                while '' in a:
                    a.remove('')
                E_sigma_0.append(eval(a[6]))
            a = linecache.getline(read_OUT_name,i)
###############################################################################
            if a.find('d Energy =') >=0:
                b = a.replace(' ','\n')
                b = b.replace('=-','=\n-')  # divid '=-'
                with open(proc_message,'w') as ff:
                    ff.write(b)
                with open(proc_message,'r') as ff:
                    a = ff.readlines()
                for j in range(len(a)):
                    b = a[j].replace('\n','')
                    a[j] = b
                while '' in a:
                    a.remove('')
                b = a[a.index('Energy')+2]  # find the second element after element 'Energy'
                E_diff.append(eval(b[0:13]))# cut first part of '0.1496133E-03-0.151E-04'
            a = linecache.getline(read_OUT_name,i)
###############################################################################
            if a.find('reached required accuracy') >=0:
                f_or_nf = 0
            if a.find('Finite differences POTIM=') >=0:
                f_or_nf = 0 # finished:0 not finished: -1
            a = linecache.getline(read_OUT_name,i)
###############################################################################
            i = i + 1 # read the next line
###############################################################################
# out of while
###############################################################################
    # E_diff.insert(0,E_diff[0])
    Energy0 = E_sigma_0[-1]
    zpe = 0
    j = 0
    for i in mev:
        if imagin[j] <= 0:
            zpe = zpe + i
        j = j + 1
    zpe = zpe / 2000
    os.remove(proc_message)
############################################################################### 
       
    return(f_or_nf,imagin,cm,mev,modes_3d,zpe,Energy0,force,E_sigma_0,E_diff)




