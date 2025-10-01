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
            ff.write('Calculation finished!')       
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
    #plt.savefig('./Force {:.3f}_eVÅ-1.png'.format(force[-1]))
    plt.cla()    
###############################################################################
    plt.plot(step, E_sigma_0) 
    plt.xlabel('Step')
    plt.ylabel('Energy (eV)')
    plt.tight_layout()
    plt.savefig('./'+'E_sigma_0')
    #plt.savefig('./E_sigma_0 {:.3f}eV.png'.format(E_sigma_0[-1]))
    plt.cla()
###############################################################################
    del step[0]       
    plt.plot(step, E_diff)
    plt.xlabel('Step')
    plt.ylabel('Energydiff (eV)')
    plt.tight_layout()
    plt.savefig('./'+'Energy_diff')
    #plt.savefig('./Energy_diff {:.3f}eV.png'.format(E_diff[-1]))
    plt.cla()
    linecache.clearcache()

def Print_CONTCAR_Conf():
    import os
    import linecache
    from ase.io import read
    from ase.io import write  # Correct import for the write function
    
    linecache.clearcache() 
    
    # Define the main directory
    main_dir = os.getcwd()
    output_format = 'png'  # Change this to 'pdf' if you prefer PDF files
    
    # Loop through all subdirectories
    for root, dirs, files in os.walk(main_dir):
        if 'CONTCAR' in files:
            contcar_path = os.path.join(root, 'CONTCAR')
            atoms = read(contcar_path)
    
            #subfolder_name = os.path.basename(root)
            subfolder_name = 'CONTCAR'
    
            top_fn = f"{subfolder_name}_Top.{output_format}"
            top_path = os.path.join(root, top_fn)
            write(top_path, atoms,
                  format = 'png',
                  show_unit_cell=2,
                  rotation="0x,0y,0z",
                  scale=30.0)
            print(f"Saved {top_path}")
    
            front_fn = f"{subfolder_name}_Front.{output_format}"
            front_path = os.path.join(root, front_fn)
            write(front_path, atoms,
                  format = 'png',
                  show_unit_cell=2,
                  rotation="-90x,0y,0z",
                  scale=50.0)
            print(f"Saved {front_path}")
            
        if 'CONTCAR_TF' in files:
            contcar_path = os.path.join(root, 'CONTCAR_TF')
            atoms = read(contcar_path)
    
            #subfolder_name = os.path.basename(root)
            subfolder_name = 'CONTCAR_TF'
    
            top_fn = f"{subfolder_name}_Top.{output_format}"
            top_path = os.path.join(root, top_fn)
            write(top_path, atoms,
                  format = 'png',
                  show_unit_cell=2,
                  rotation="0x,0y,0z",
                  scale=50.0)
            print(f"Saved {top_path}")
    
            front_fn = f"{subfolder_name}_Front.{output_format}"
            front_path = os.path.join(root, front_fn)
            write(front_path, atoms,
                  format = 'png',
                  show_unit_cell=2,
                  rotation="-90x,0y,0z",
                  scale=50.0)
            print(f"Saved {front_path}")
    linecache.clearcache()

def Print_Check_Keywords(read_OUT_name, nkps):
    import os, re
    
    # RED = "\033[31m"
    # RESET = "\033[0m"

    CHECK_FILE = "check_keywords"
    if not os.path.exists(CHECK_FILE):   
        with open(CHECK_FILE, "w") as f:
            pass
    OUTCAR_FILE = read_OUT_name
    REPORT_FILE = "keyword_check_report.txt"
    RTOL, ATOL = 1e-6, 1e-6

    def is_number_token(tok: str) -> bool:
        t = tok.strip().strip(",").replace("d","e").replace("D","E")
        if not t: return False
        try: float(t); return True
        except: return False

    def to_float(tok: str) -> float:
        return float(tok.strip().strip(",").replace("d","e").replace("D","E"))

    def almost_equal(a, b, rtol=RTOL, atol=ATOL):
        return abs(a - b) <= max(atol, rtol * max(abs(a), abs(b)))

    def str_loose_equal(a: str, b: str) -> bool:
        """Case-insensitive, prefix-loose match."""
        a = a.strip().lower()
        b = b.strip().lower()
        # exact OR either is a prefix of the other
        return a == b or a.startswith(b) or b.startswith(a)

    def load_checklist(path: str):
    
        items = []
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                raw = ln.strip()
                if not raw or raw.startswith("#"):
                    continue

                s = raw.split("#", 1)[0].strip()
                parts = s.split()
                if not parts:
                    continue
                key = parts[0]
                rest = s[len(key):].strip()
    
                if key.upper() in {"KPOINTS", "KPOINT"}:

                    m = re.search(r'(\d+)[^\d]+(\d+)[^\d]+(\d+)', rest)
                    if m:
                        val = f"{m.group(1)} {m.group(2)} {m.group(3)}"
                    else:
 
                        val = rest
                else:

                    m = re.match(r'^\s*=?\s*(\S+)', rest)
                    val = m.group(1) if m else ""
    
                items.append((key, val))
        return items


    def extract_outcar_values(outcar_path: str, keys):
        regex_map = {key: re.compile(rf"(?<!\S){re.escape(key)}\s*=\s*([^\s;]+)")
                     for key in keys}
        found = {key: None for key in keys}
        with open(outcar_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                for key, rgx in regex_map.items():
                    for m in rgx.finditer(line):
                        found[key] = m.group(1)
        return found


    def parse_kmesh_list_any(x):
        if isinstance(x, (list, tuple)):
            if len(x) != 3: raise ValueError("kmesh_list must have 3 integers")
            return tuple(int(v) for v in x)
        if isinstance(x, str):
            nums = re.findall(r'\d+', x)
            if len(nums) != 3:
                raise ValueError(f"kmesh_list string not recognized: {x}")
            return tuple(int(n) for n in nums)
        raise TypeError("kmesh_list must be list/tuple of 3 ints or string like '4 5 1'")

    kmesh_list_tuple = parse_kmesh_list_any(nkps)


    if not os.path.isfile(CHECK_FILE):
        raise FileNotFoundError(f"Cannot find '{CHECK_FILE}' in current directory.")
    if not os.path.isfile(OUTCAR_FILE):
        raise FileNotFoundError(f"Cannot find '{OUTCAR_FILE}' in current directory.")

    checklist = load_checklist(CHECK_FILE)
    if not checklist:
        print(f"Warning '{CHECK_FILE}' is empty or has no valid lines.")


    normal_keys = [k for k, _ in checklist if k.upper() not in {"KPOINTS", "KPOINT"}]
    outcar_vals = extract_outcar_values(OUTCAR_FILE, normal_keys)

    lines_out = []
    c_or_nc = 1

    for key, check_val_str in checklist:
        U = key.upper()

        if U in {"KPOINTS", "KPOINT"}:

            try:
                kmesh_list_from_check = parse_kmesh_list_any(check_val_str)
            except Exception:
                c_or_nc = -1
                lines_out.append(f"KPOINTS {check_val_str} (OUTCAR: {kmesh_list_tuple[0]} {kmesh_list_tuple[1]} {kmesh_list_tuple[2]}) Not Consistent")
                continue

            exp_str = f"{check_val_str}"
            kmesh_list_str  = f"{kmesh_list_tuple[0]} {kmesh_list_tuple[1]} {kmesh_list_tuple[2]}"
            if kmesh_list_from_check == kmesh_list_tuple:
                lines_out.append(f"KPOINTS {exp_str} Consistent.")
            else:
                c_or_nc = -1
                lines_out.append(f"KPOINTS {exp_str} (OUTCAR: {kmesh_list_str}) Not Consistent")
            continue


        out_val_str = outcar_vals.get(key)
        if out_val_str is None:
            c_or_nc = -1
            lines_out.append(f"{key} {check_val_str} (OUTCAR: N/A) Not Consistent")
            continue

        if is_number_token(check_val_str) and is_number_token(out_val_str):
            a, b = to_float(check_val_str), to_float(out_val_str)
            if almost_equal(a, b):
                lines_out.append(f"{key} {check_val_str} Consistent.")
            else:
                c_or_nc = -1
                lines_out.append(f"{key} {check_val_str} (OUTCAR: {out_val_str}) Not Consistent")
        else:
            if str_loose_equal(check_val_str, out_val_str):
                lines_out.append(f"{key} {check_val_str} Consistent.")
            else:
                c_or_nc = -1
                lines_out.append(f"{key} {check_val_str} (OUTCAR: {out_val_str}) Not Consistent")

    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        for ln in lines_out:
            f.write(ln + "\n")

    print(f"Done. Report saved to: {os.path.abspath(REPORT_FILE)}")
    return c_or_nc

def Print_Summary_pdf(f_or_nf,c_or_nc,Energy0,force,E_sigma_0,E_diff):
    from PIL import Image
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.gridspec import GridSpec
    import os
    import numpy as np
    
    # ---- helpers ----
    KEYS = ["E_sigma_0", "Energydiff", "Force"]
    
    def short_label_from_path(path):
        base = os.path.splitext(os.path.basename(path))[0]
        for k in KEYS:
            if base.lower().startswith(k.lower()):
                return k
        return base
    
    def last_val(x):
        return float(np.ravel(x)[-1])  # robust for list/ndarray, 1D/ND
    

    series_map = {
        "E_sigma_0": (lambda: last_val(E_sigma_0), "eV"),
        "Energy_diff": (lambda: last_val(E_diff),   "eV"),
        "Force":      (lambda: last_val(force),    "eV/Å"),
    }
        
    # Define the main directory
    main_dir = os.getcwd()
    folder = os.getcwd()
    output_pdf = os.path.join(main_dir, 'Summary.pdf')  # Save PDF in the main directory
    # if f_or_nf < 0:
    #     output_pdf = os.path.join(main_dir, 'Summary (not Converged).pdf')  # Save PDF in the main directory
    # else:
    #     output_pdf = os.path.join(main_dir, 'Summary (Converged).pdf')  # Save PDF in the main directory
    

    desired_order = [
        "E_sigma_0.png",
        "Energy_diff.png",
        "Force.png",
        "CONTCAR_Top.png",
        "CONTCAR_Front.png",
        "CONTCAR_TF_Top.png",
        "CONTCAR_TF_Front.png",
    ]
    desired_lower = [n.lower() for n in desired_order]
    

    found = {name: None for name in desired_lower}
    for root, dirs, files in os.walk(main_dir):
        for f in files:
            fl = f.lower()
            if fl in found and found[fl] is None:
                found[fl] = os.path.join(root, f)
    
    top_row = [found[desired_lower[i]] for i in range(3)]
    bottom_row = [found[desired_lower[i]] for i in range(3, 7)]
    

    dpi_pdf = 600
    figsize = (8.27, 11.69)
    interp = 'nearest'
    
    with PdfPages(output_pdf) as pdf:
        fig = plt.figure(figsize=figsize, dpi=dpi_pdf)
        if f_or_nf < 0:
            if c_or_nc < 0:
                fig.text(0.5, 0.99, folder, ha='center', va='top', fontsize=9)
                fig.suptitle("Summary (Not Converged & Not Consistent)", fontsize=16, y=0.965)
            else:
                fig.text(0.5, 0.99, folder, ha='center', va='top', fontsize=9)
                fig.suptitle("Summary (Not Converged & Consistent)", fontsize=16, y=0.965)
        else:
            if c_or_nc < 0:
                fig.text(0.5, 0.99, folder, ha='center', va='top', fontsize=9)
                fig.suptitle("Summary (Converged & Not Consistent)", fontsize=16, y=0.965)
            else:
                fig.text(0.5, 0.99, folder, ha='center', va='top', fontsize=9)
                fig.suptitle("Summary (Converged & Consistent)", fontsize=16, y=0.965)
       

        report_path = os.path.join(main_dir, "keyword_check_report.txt")
        report_txt = None
        if os.path.isfile(report_path):
            with open(report_path, "r", encoding="utf-8", errors="ignore") as f:
                lines = [ln.strip() for ln in f if ln.strip()]
            if lines:

                report_txt = "\n".join(lines)

        
        if report_txt:

            fig.text(
                0.33, 0.9,
                report_txt,
                fontsize=10,
                family="monospace",
                va="top",
                ha="left"
            )
            

        gs = GridSpec(3, 12, figure=fig, wspace=0.08, hspace=0.15)

        for i in range(3):
            ax = fig.add_subplot(gs[1, 4*i:4*(i+1)])
            p = top_row[i]
            if p:
                img = Image.open(p)
                ax.imshow(img, interpolation=interp)
                label = short_label_from_path(p)
                if label in series_map:
                    getter, unit = series_map[label]
                    val = getter()
                    ax.set_title(f"{label}: {val:.3f} {unit}", fontsize=11)
                else:
                    ax.set_title(label, fontsize=11)
                img.close()
            ax.axis('off')
    

        for i in range(4):
            ax = fig.add_subplot(gs[2, 3*i:3*(i+1)])
            p = bottom_row[i]
            if p:
                img = Image.open(p)
                ax.imshow(img, interpolation=interp)
                ax.set_title(os.path.basename(p), fontsize=11)
                img.close()
            ax.axis('off')

        pdf.savefig(fig, dpi=dpi_pdf)
        plt.close(fig)
    
    print(f"Summary PDF created: {output_pdf}")


