#!/bin/env python3
import os
import shutil
import traceback
import linecache
from contextlib import contextmanager

# Your modules
import Parsing_VASP
import Modifing_VASP
import Printing_VASP

CHECK_FILE = "check_keywords"
if not os.path.exists(CHECK_FILE):   
    with open(CHECK_FILE, "w") as f:
        pass

@contextmanager
def pushd(path):
    prev = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev)

def run_pipeline_in_dir(workdir, src_check_keywords):
    """
    Execute the original workflow inside `workdir`:
      1) Copy check_keywords from root
      2) Parse OUTCAR (results and kpoints)
      3) Plot figures
      4) Parse/modify CONTCAR and write outputs
      5) Generate structure images
      6) Check keywords (using k-mesh list)
      7) Create per-directory summary.pdf
    """
    with pushd(workdir):
        linecache.clearcache()

        # 1) Copy check_keywords
        if src_check_keywords and os.path.isfile(src_check_keywords):
            try:
                shutil.copy2(src_check_keywords, os.path.join(workdir, "check_keywords"))
                print(f"[{workdir}] Copied check_keywords")
            except Exception as e:
                print(f"[{workdir}] [WARN] Failed to copy check_keywords: {e}")
        else:
            print(f"[{workdir}] [WARN] Root check_keywords not found; continuing")

        # 2) Parse OUTCAR
        outcar = "OUTCAR"
        if not os.path.isfile(outcar):
            print(f"[{workdir}] [SKIP] OUTCAR not found")
            return

        try:
            f_or_nf, imagin, cm, mev, modes_3d, zpe, Energy0, force, E_sigma_0, E_diff = \
                Parsing_VASP.Parsing_OUTCAR_Result(outcar)
            print(f"[{workdir}] Parsed OUTCAR result")
        except Exception:
            print(f"[{workdir}] [ERROR] Parsing_OUTCAR_Result failed:\n{traceback.format_exc()}")
            return

        # Parse KPOINTS (inferred from OUTCAR)
        try:
            nkps = Parsing_VASP.Parsing_OUTCAR_Kpoints(outcar)
        except Exception:
            print(f"[{workdir}] [WARN] Parsing_OUTCAR_Kpoints failed:\n{traceback.format_exc()}")
            nkps = None

        # 3) Plot line figures
        try:
            Printing_VASP.Print_Figures(f_or_nf, Energy0, force, E_sigma_0, E_diff)
            print(f"[{workdir}] Figures saved")
        except Exception:
            print(f"[{workdir}] [WARN] Print_Figures failed:\n{traceback.format_exc()}")

        # 4) Parse/modify CONTCAR and write
        poscar = "CONTCAR"
        if os.path.isfile(poscar):
            try:
                comment, sc_factor, lattice_con, atom_type, atom_num, atom_num_tot, sel_or_not, \
                    D_or_C, atoms_coor, atoms_relax = Parsing_VASP.Parsing_POSCAR(poscar)

                atom_type, atom_num, atoms_coor = Modifing_VASP.Gen_Type_Show_POSCAR(
                    atoms_relax, atom_type, atoms_coor, atom_num
                )

                gen_pos_name = "CONTCAR_TF"
                Printing_VASP.Print_POSCAR(
                    gen_pos_name, comment, sc_factor, lattice_con,
                    atom_type, atom_num, atom_num_tot, sel_or_not, D_or_C,
                    atoms_coor, atoms_relax
                )
                print(f"[{workdir}] Wrote {gen_pos_name}")
            except Exception:
                print(f"[{workdir}] [WARN] POSCAR processing failed:\n{traceback.format_exc()}")

            # 5) Export top/front structure images
            try:
                Printing_VASP.Print_CONTCAR_Conf()
                print(f"[{workdir}] Structure images saved")
            except Exception:
                print(f"[{workdir}] [WARN] Print_CONTCAR_Conf failed:\n{traceback.format_exc()}")
        else:
            print(f"[{workdir}] [WARN] CONTCAR not found; skipping structure steps")

        # 6) Keyword check (compare with kmesh_list where applicable)
        try:
            c_or_nc = Printing_VASP.Print_Check_Keywords(outcar, nkps)
            print(f"[{workdir}] keyword_check_report.txt generated (status={c_or_nc})")
        except Exception:
            print(f"[{workdir}] [WARN] Print_Check_Keywords failed:\n{traceback.format_exc()}")
            c_or_nc = -1

        # 7) Per-directory summary PDF
        try:
            Printing_VASP.Print_Summary_pdf(f_or_nf, c_or_nc, Energy0, force, E_sigma_0, E_diff)
            print(f"[{workdir}] summary.pdf saved")
        except Exception:
            print(f"[{workdir}] [WARN] Print_Summary_pdf failed:\n{traceback.format_exc()}")

        linecache.clearcache()

def merge_all_summaries(root, out_pdf_name="All_Summaries.pdf"):
    # Try pypdf first, fall back to PyPDF2
    try:
        from pypdf import PdfMerger, PdfReader
        PYPDF = True
    except Exception:
        from PyPDF2 import PdfMerger, PdfReader  # type: ignore
        PYPDF = False

    def is_bad_dir(p):
        name = os.path.basename(p)
        return name == "__pycache__" or name.startswith(".")

    # Collect all Summary.pdf / summary.pdf (case-insensitive)
    pdfs = []
    for dirpath, dirnames, filenames in os.walk(root):
        if dirpath == root or is_bad_dir(dirpath):
            continue
        for fn in filenames:
            if fn.lower() == "summary.pdf":
                pdfs.append(os.path.join(dirpath, fn))
                break

    if not pdfs:
        print("[MERGE] No Summary.pdf found")
        return

    pdfs = sorted(pdfs)
    merger = PdfMerger()
    page_offset = 0

    for p in pdfs:
        try:
            # how many pages this PDF has
            n_pages = len(PdfReader(p).pages)
            # append without deprecated bookmark arg
            merger.append(p)

            # add a top-level outline/bookmark at the starting page
            title = os.path.relpath(os.path.dirname(p), root)
            try:
                merger.add_outline_item(title, page_offset)
                # if PYPDF:
                    # merger.add_outline_item(title, page_offset)
                # else:
                    # merger.addBookmark(title, page_offset)  # PyPDF2
            except Exception as e:
                print(f"[MERGE] [WARN] cannot add outline for {p}: {e}")

            page_offset += n_pages
            print(f"[MERGE] appended {p} ({n_pages} pages)")
        except Exception as e:
            print(f"[MERGE] [WARN] skip {p}: {e}")

    out_path = os.path.join(root, out_pdf_name)
    try:
        with open(out_path, "wb") as f:
            merger.write(f)
        print(f"[MERGE] Done -> {out_path}")
    finally:
        try:
            merger.close()
        except Exception:
            pass



def collect_target_dirs(root):
    targets = []
    for dirpath, dirnames, filenames in os.walk(root):
        if dirpath == root:
            continue
        name = os.path.basename(dirpath)
        if name == "__pycache__" or name.startswith("."):
            continue
        if os.path.isfile(os.path.join(dirpath, "OUTCAR")):
            targets.append(dirpath)
    return sorted(targets)

def main():
    root = os.getcwd()
    src_check = os.path.join(root, "check_keywords")  # root-level check_keywords

    # Collect all subdirectories (recursive)
    subdirs = collect_target_dirs(root)
    # for dirpath, dirnames, filenames in os.walk(root):
        # if dirpath == root:
            # continue
        # subdirs.append(dirpath)

    if not subdirs:
        print("[INFO] No subdirectories found.")
        return

    print(f"[INFO] Found {len(subdirs)} subdirectories. Start processing...\n")

    for d in sorted(subdirs):
        print("=" * 80)
        print(f"[RUN] {d}")
        print("=" * 80)
        try:
            run_pipeline_in_dir(d, src_check_keywords=src_check)
        except Exception:
            print(f"[{d}] [FATAL] pipeline crashed:\n{traceback.format_exc()}")

    # Merge all per-directory summaries
    merge_all_summaries(root, out_pdf_name="All_Summaries.pdf")

    print("\n[ALL DONE]")

if __name__ == "__main__":
    main()