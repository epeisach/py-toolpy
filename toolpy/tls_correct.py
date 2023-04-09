#!/usr/bin/env python

import os, sys, shutil, math, re
from cifparse import cif2pdb


ERRLOG = []


def usage():
    content = """

###############################################################################
This script is to detect/correct various TLS problems generated from refmac/phenix/buster 
              (created 2010-06-20 and modified as necessary  )
-------------------------------------------------------------------------------
Usages: tls.py [option]  inputfile

1.  tls -pdb processed_file -dep deposited_file
    The most used option to do various corrections: input cif or pdb file!


2.  tls -pdb  pdbfile/ciffile/pdbid : (only for small correction)
    Correct various tls problems in pdb. If it is phenix, first convert to
    refmac, then get full B and do some correction.

3.  tls -all pdbfile/pdbid : (2010-10-10, system depended)
    Check all the rcsb0????? directory and find the deposited one and map
    all the TLS to the released one. All the mingled ligand/water chain 
    residue/id will be mapped to the processed one.

4.  tls -anis  pdbfile : (2010-11-08)
    get ANISOU records (if missing) from Biso_total.
    
5.  tls -bres  pdbfile : (2010-12-03)
    get B_residual of atoms in the tls group from Biso_total.
    
6.  tls -btls  pdbfile : (2010-12-03)
    get B_tls for atoms in tls group from Biso_total.
    
7.  tls -exttls  pdbfile : (2010-12-20)
    extract TLS from PDB. The tls file is the same format as PDB
    
8.  tls -newtls  pdbfile1 pdbfile2 : (2010-12-30)
    Put TLS groups from pdbfile2 to pdbfile1. 
    
9.  tls -tls2cif  pdbfile/tls.inp : (2011-01-08)
    Extract TLS from PDBfile or tls.inp (for refmac) and convert to CIF.

10. tls -ntls  deposited_pdbfile : (2011-02-18)
    get a new TLS cif file (taking care of HOH/LIG chainID/residue number change)

11. tls -xc  pdbfile : (2011-03-10)
    Calculated center of origin of TLS groups (moved from the auto)

12. tls -ph2ref  pdbfile : (2011-07-20)
    convert tls in phenix format to refmac
    
13. tls -ref2ph  pdbfile : (2011-07-20)
    convert tls in refmac to phenix format 
    
14. tls -bus2ref  pdbfile : (2011-09-02)
    convert tls in buster format to refmac
    
15. tls -ref2bus  pdbfile
    convert tls in buster format to refmac

16. tls -checktls  pdbfile : 
    check TLS problems if exist.


If add '-o output', user specify the output file name for the TLS group.
###############################################################################

"""
    print(content)
    sys.exit()


###############################################################################
def process(*files):
    """Add TLS contribution to the ATOM records to get the total B factors.
    (a correction for REFMAC that only export B_residual when TLS was used.

    """

    global ERRLOG

    arg = files[0]
    if len(arg) < 2:
        usage()

    pdbfile, pdbdep, outfile = "", "", ""
    opt, all_dep, dep, anis, ncs = {}, 0, 0, 0, 0

    for k in range(len(arg)):
        if arg[k].upper() == "-PDB" or arg[k].upper() == "-F":
            opt["pdb"] = 1
            n1 = is_cif(arg[k + 1])

            if n1 > 0:
                pdbfile = cif2pdb(arg[k + 1])
            else:
                pdbfile = arg[k + 1]

            if not os.path.exists(pdbfile):
                pdbfile = get_file_by_pdbid(pdbfile, 1)

        elif arg[k].upper() == "-DEP":
            opt["dep"] = 1
            n1 = is_cif(arg[k + 1])
            if n1 > 0:
                dep, pdbdep = 1, cif2pdb(arg[k + 1])
            else:
                dep, pdbdep = 1, arg[k + 1]

        elif arg[k].upper() == "-O":
            outfile = arg[k + 1]

        elif arg[k].upper() == "-NCS":
            ncs = 1

        elif arg[k].upper() == "-ALL_DEP" or arg[k].upper() == "-ALL":
            opt["all_dep"] = 1
            all_dep, pdbdep = 1, arg[k + 1]
            if not os.path.exists(pdbdep):
                pdbdep = get_file_by_pdbid(pdbdep, 1)

        elif arg[k].upper() == "-NTLS":
            cif = get_newtls_cif(arg[k + 1])
            out = "%s_tls2cif" % arg[k + 1]
            os.system("mv -f %s %s" % (cif, out))
            print("The output file = %s\n" % out)

            sys.exit()

        elif arg[k].upper() == "-ANIS":  # get ANISOU record from Biso_total.
            anis, pdbfile = 1, arg[k + 1]
            if is_cif(arg[k + 1]):
                pdbfile = cif2pdb(arg[k + 1])
            outfile = pdbfile + "_anis"

        elif arg[k].upper() == "-BRES":  # get B_residual from Biso_total
            anis, pdbfile = 2, arg[k + 1]
            if is_cif(arg[k + 1]):
                pdbfile = cif2pdb(arg[k + 1])
            outfile = pdbfile + "_bres"

        elif arg[k].upper() == "-BTLS":  # get B_tls from Biso_total
            anis, pdbfile = 3, arg[k + 1]
            if is_cif(arg[k + 1]):
                pdbfile = cif2pdb(arg[k + 1])
            outfile = pdbfile + "_btls"

        elif arg[k].upper() == "-EXTTLS":  # extract TLS from PDB (same format)
            pdbfile = arg[k + 1]
            if is_cif(arg[k + 1]):
                pdbfile = cif2pdb(arg[k + 1])
            cif = extract_tls_pdb(pdbfile)
            sys.exit()

        elif arg[k].upper() == "-TLS2CIF":  # convert TLS to cif(sequential)
            pdbfile = arg[k + 1]
            if is_cif(pdbfile):
                pdbfile = cif2pdb(arg[k + 1])
            cif = tls2cif(pdbfile)
            out = "%s_tls2cif" % arg[k + 1]
            os.system("mv -f %s %s" % (cif, out))
            print("\nThe output cif file of TLS groups = %s\n" % out)

            sys.exit()

        elif arg[k].upper() == "-NEWTLS":  # overwrite TLS group in PDB1 by pdb2
            wk, pdb_new = get_newtls(arg[k + 1], arg[k + 2])
            sys.exit()

        elif arg[k].upper() == "-XC":
            pdbfile = arg[k + 1]
            if is_cif(pdbfile):
                pdbfile = cif2pdb(arg[k + 1])

            tlsxc_origin(pdbfile, 1)
            sys.exit()

        elif arg[k].upper() == "-BUS2REF":
            convert_tls2refmac(arg[k + 1], "buster")
            sys.exit()

        elif arg[k].upper() == "-REF2BUS":
            convert_refmac2buster(arg[k + 1])
            sys.exit()

        elif arg[k].upper() == "-PH2REF":
            convert_tls2refmac(arg[k + 1], "phenix")
            sys.exit()

        elif arg[k].upper() == "-REF2PH":
            convert_refmac2phenix(arg[k + 1])
            out = arg[k + 1] + "_ref2ph"
            print("\nTLS is converted to phenix format. Output=%s\n" % out)
            sys.exit()

        elif arg[k].upper() == "-CHECKTLS":
            pdbfile = arg[k + 1]
            if is_cif(pdbfile):
                pdbfile = cif2pdb(arg[k + 1])
            check_tls(pdbfile, 0)
            sys.exit()

        elif arg[k].upper() == "-NUM":  # auto renumber the waters/ligands
            check_and_correct_tls(arg[k + 1])
            new_ligand_water_residue_number(arg[k + 1])
            sys.exit()

    # ---------------------------------------------------
    if len(outfile) == 0:
        if os.path.exists(pdbfile):
            outfile = pdbfile + "_new_tls"
        else:
            outfile = pdbdep + "_new_tls"

    prog = "REFMAC"
    tmp = os.popen('grep "REMARK   3   PROGRAM     :" -m 1  %s' % pdbfile).read()
    if tmp:
        prog = tmp.split(":")[1].strip()

    print("The deposited file format = %s" % prog)

    if not os.path.exists(pdbfile) and not os.path.exists(pdbdep):
        print("Error! no input pdb file.")
        sys.exit()

    if anis > 1:
        recover_anisou(anis, pdbfile, prog, outfile)

    elif all_dep == 1:  # search all file in the rcsb dir
        pdb_sort = reorder_pdb(pdbdep)
        search_rcsb_dir(pdb_sort, outfile, 1)
        delete_file(pdb_sort)

    elif dep == 1 and len(pdbfile) == 0:  # use only one file in the rcsb dir
        search_rcsb_dir(pdbdep, outfile, 0)

    elif len(pdbfile) > 0 and len(pdbdep) > 0 and "PHENIX" in prog.upper():
        map_tls_4_phenix(pdbfile, pdbdep, outfile)

    elif len(pdbfile) > 0 and len(pdbdep) > 0 and "BUSTER" in prog.upper():
        map_tls_4_buster(pdbfile, pdbdep, outfile)

    elif len(pdbfile) > 0 and len(pdbdep) > 0:
        print("Doing tls using %s %s " % (pdbfile, pdbdep))
        pdb_sort = reorder_pdb(pdbfile)
        wk, pdb_match = match_all_residues(pdb_sort, pdbdep)

        if wk == 1:
            tlspdb = outfile + ".tlspdb"
            proc_tlsanl(pdb_match, tlspdb)
            cif = tls2cif(tlspdb)

            shutil.move(cif, outfile)
            print("The new output PDB file = %s" % tlspdb)
            print("The new TLS groups in cif = %s\n" % outfile)

        delete_file(pdb_sort, pdb_match)

    elif len(pdbfile) > 0 and "PHENIX" in prog.upper():
        ph2ref = convert_tls2refmac(pdbfile, "phenix")
        proc_tlsanl(ph2ref, outfile)
        # cif = tls2cif(outfile)
        delete_file(ph2ref)

    elif len(pdbfile) > 0 and "BUSTER" in prog.upper():
        ph2ref = convert_tls2refmac(pdbfile, "buster")
        proc_tlsanl(ph2ref, outfile)
        # cif = tls2cif(outfile)
        delete_file(ph2ref)

    else:  # only one pdbfile
        pdb_sort = reorder_pdb(pdbfile)
        proc_tlsanl(pdb_sort, outfile)
        # cif = tls2cif(outfile)
        delete_file(pdb_sort)

    if os.path.exists("maxit.err"):
        delete_file("maxit.err")


##########################################################
def perror(info):
    """print error messages"""

    if info not in ERRLOG:
        ERRLOG.append(info)
        print(info.strip())


##########################################################
def check_tls(pdb, id):
    """
    id = 0, pdb is a file; id=1, pdb is a list
    """

    fp = pdb
    if id == 0 and os.path.exists(pdb):
        fp = open(pdb, "r").readlines()
    tls, pdbn = remove_pdb_item(fp, 1, "tls")
    if len(tls) < 15:
        return

    chain, ch_range = chain_res_list(pdbn, 1)
    ntls, prog = "?", "?"
    all_range = []  # [['ch', ni, nf, ntls] ...] ni/nf=first/last residue number
    prog = os.popen('grep "REMARK   3   PROGRAM     :" -m 1  %s' % pdb).read().split(":")[1].strip()

    for i, x in enumerate(tls):
        if "TLS GROUP :" in x[:25]:
            ntls = x.split(":")[1].strip()
            if not is_number(ntls):
                t1 = "Error: TLS group (=%s) is not a number \n" % (ntls)
                perror(t1)
                continue
            ntls = int(ntls)

        elif ("REMARK   3    RESIDUE RANGE :" in x[:29]) or ("REMARK   3    SELECTION:" in x[:29] and "buster" in prog.lower() and "{" not in x[20:]):
            if "SELECTION:" in x:  # buster changed set to  SELECTION
                xx = "REMARK   3    RESIDUE RANGE :" + x[24:48]
                t = parse_tls_range_refmac(xx, ntls)
            else:  # refmac
                t = parse_tls_range_refmac(x, ntls)
            t.append(ntls)
            all_range.append(t)
        #            print 'from refmac=',  t

        elif "REMARK   3    SELECTION:" in x[:29] and "phenix" in prog.lower():  # phenix,
            tls_str = tls_selection_one_string(tls, i, 25)
            t = parse_tls_range_phenix(tls_str, ntls, chain, ch_range, fp)
            for y in t:
                y.append(ntls)
                all_range.append(y)
            # print 'from phenix=', t

        elif "REMARK   3    SET :" in x[:20] or "SELECTION: {" in x[11:29]:  # buster in original format  # buster in modified format
            tls_str = tls_selection_one_string(tls, i, 20)
            # print tls_str
            t = parse_tls_range_buster(tls_str, ntls, chain, ch_range)
            for y in t:
                y.append(ntls)
                all_range.append(y)
            # print 'from buster=', t

    # print all_range
    if not all_range:
        return
    for x in all_range:
        if len(x) < 4:
            continue
        ch, n1, n2, nt = x[0], x[1], x[2], x[3]

        if ch not in chain.keys():
            t1 = "Error: chainID (%s) not in coordinate. TLS group=%d !" % (ch, nt)
            perror(t1)

        if ch in chain.keys() and n1 not in chain[ch]:
            t1 = "Warning: residue number (%d) not in coordinate. TLS group=%d." % (n1, nt)
            perror(t1)

        if ch in chain.keys() and n2 not in chain[ch]:
            t1 = "Warning: residue number (%d) not in coordinate. TLS group=%d." % (n2, nt)
            perror(t1)

    check_residue_range(all_range)
    if id == 0:
        tlsxc_origin(pdb, 0)  # check TLS origin, only pdbfile and refmac


##########################################################
def check_residue_range(all_range):
    """check residue range overlap :
    four type of overlaps for residue ranges (n1, n2) and (k1, k2)
    (n1, k1, k2, n2;   k1, n1, n2, k2;  n1, k1, n2, k2;   k1, n1, k2, n2;)
    """

    nr = len(all_range)
    for i in range(0, nr):
        if len(all_range[i]) < 4:
            continue
        n1, n2, nt1 = all_range[i][1], all_range[i][2], all_range[i][3]
        if n2 < n1:
            t = "Error! Residue range problem [%d %d] in TLS group %d." % (n1, n2, nt1)
            perror(t)
        for j in range(i + 1, nr):
            if all_range[i][0] != all_range[j][0]:
                continue
            k1, k2, nt2 = all_range[j][1], all_range[j][2], all_range[j][3]
            ch = all_range[j][0]
            if (n1 <= k1 and n2 >= k2) or (k1 <= n1 and k2 >= n2) or (k1 >= n1 and n2 >= k1 and k2 >= n2) or (n1 >= k1 and k2 >= n1 and n2 >= k2):  # 1  # 2  # 3  # 4
                t = "Error! TLS group overlaps: group=%d, residue (%d  %d) with group=%d, residue (%d  %d) ; chainID=%s " % (nt1, n1, n2, nt2, k1, k2, ch)
                perror(t)


##########################################################
def tlsxc_origin(pdbfile, info):
    """calculate the origin for TLS groups (occupancy weighted)
    refmac & buster: center of coordinate.
    phenix: center of mass
    info ==1, print all. info==0, only print error/warning

    """

    fp = open(pdbfile, "rU").readlines()

    pdbid, prog, tls, xyz = "", "?", 0, []
    for x in fp:
        if "HEADER" in x[:6]:
            pdbid = x[62:66].strip()

        elif "REMARK   3   PROGRAM     :" in x:
            prog = x.split(":")[1].strip()

        elif "S21:" in x and "S22:" in x:
            tls = tls + 1

        elif "ATOM  " in x[:6] or "HETATM" in x[:6]:
            xyz.append(x)

        elif "ENDMDL" in x[:6]:
            print("Warning. This is a multiple model entry (first model used!).\n")
            break

    if tls == 0:
        print("Note: No TLS groups are detected in file=(%s)." % pdbfile)
        return

    chain, ch_range = chain_res_list(fp, 1)

    nxyz, ntls, nph, nref, nbu = [], "?", 0, 0, 0
    for i, x in enumerate(fp):
        if "REMARK" in x[:6] and "BUSTER" in x[24:35]:
            nbu = 1

        elif "TLS GROUP :" in x[:25]:
            ntls = x.split(":")[1].strip()

        elif "REMARK   3    RESIDUE RANGE :" in x[:29]:  # refmac
            nref = nref + 1
            if nref == 1 and info == 1:
                print("\nNOTE: TLS format is REFMAC. Origin is OCC weighted center of coordinate.")

            t = parse_tls_range_refmac(x, ntls)
            if not t:
                continue
            tmp = [t[0], t[1], t[0], t[2]]

            txyz = get_xyz_tls(tmp, xyz)
            nxyz.extend(txyz)

        elif "REMARK   3    SELECTION:" in x[:29] and "phenix" in prog.lower():  # phenix,
            nph = nph + 1
            if nph == 1 and info == 1:
                print("\nNOTE: TLS format is PHENIX. Origin is OCC weighted center of mass.")

            tls_str = tls_selection_one_string(fp, i, 25)
            if not tls_str:
                continue
            t = parse_tls_range_phenix(tls_str, ntls, chain, ch_range, fp)
            for z in t:
                if len(z) < 3:
                    continue
                tmp = [z[0], z[1], z[0], z[2]]
                txyz = get_xyz_tls(tmp, xyz)
                nxyz.extend(txyz)

        elif ("REMARK   3    SET :" in x[:20]) or ("REMARK   3    SELECTION:" in x[:29] and "buster" in prog.lower()):  # buster
            nbu = nbu + 1
            if nbu == 1 and info == 1:
                print("\nNOTE: TLS format is BUSTER. Origin is OCC weighted center of coordinate.")
                print("        The origin column size is not the same as phenix & refmac. ")

            if "SELECTION:" in x and "{" not in x[20:]:  # buster (tls in phenix/refmac format)
                t = parse_tls_range_refmac(x, ntls)
                if not t:
                    continue
                tmp = [t[0], t[1], t[0], t[2]]
                txyz = get_xyz_tls(tmp, xyz)
                nxyz.extend(txyz)

            else:
                if "{" in x[20:]:
                    tls_str = tls_selection_one_string(fp, i, 25)
                else:
                    tls_str = tls_selection_one_string(fp, i, 20)

                t = parse_tls_range_buster(tls_str, ntls, chain, ch_range)
                for z in t:
                    if len(z) < 3:
                        continue
                    tmp = [z[0], z[1], z[0], z[2]]
                    txyz = get_xyz_tls(tmp, xyz)
                    nxyz.extend(txyz)

        elif "REMARK   3    ORIGIN FOR THE GROUP" in x:
            xc = [x[39:48], x[48:57], x[57:66]]
            if nbu > 0:
                xc = [x[39:49], x[49:59], x[59:69]]

            if not is_number(xc[0]) or not is_number(xc[1]) or not is_number(xc[2]):
                t = "Error: Wrong TLS origin (%s) for group=%s" % (xc, ntls)
                perror(t)
                continue

            nxc = []
            oxc = [float(y) for y in xc]
            if nref:
                nxc = calc_tlsxc(nxyz, "refmac")
            elif nph:
                nxc = calc_tlsxc(nxyz, "phenix")
            elif nbu:
                nxc = calc_tlsxc(nxyz, "buster")

            if len(nxc) < 1:
                t = "\nWarning: Origin not calculated. No xyz found in TLS group=%s" % ntls
                perror(t)
            else:
                dif = [math.fabs(nxc[0] - oxc[0]), math.fabs(nxc[1] - oxc[1]), math.fabs(nxc[2] - oxc[2])]
                if info == 1:
                    print("\nTLS origin (calculated): %8.3f %8.3f %8.3f" % (nxc[0], nxc[1], nxc[2]))
                    print("TLS origin (reported  ): %8.3f %8.3f %8.3f" % (oxc[0], oxc[1], oxc[2]))
                if max(dif) > 2:
                    t = "(%s:) Warning: TLS origin differs from calculated(max_dev=%.2f),group=%s." % (pdbid, max(dif), ntls)
                    perror(t)
                else:
                    if info == 1:
                        print("(%s:) TLS origin is similar to the calculated for group=%s." % (pdbid, ntls))
            #  print('tls_range= ', tmp )

            nxyz = []

        elif "MODEL" in x[:6]:
            print("Warning. This is a multiple model entry!\n")
        elif "ATOM" in x[:6]:
            break
    # for x in xyz : print(x)


##########################################################
def calc_tlsxc(nxyz, id):
    """get the center of coordinate (refmac) and center of mass (phenix)"""

    occup = 0.0
    for x in nxyz:
        if len(x.strip()) > 60:
            occup = occup + float(x[55:60])

    if occup <= 0:
        return []

    occmass, xt, yt, zt = 0.0, 0.0, 0.0, 0.0
    for ln in nxyz:
        x, y, z = float(ln[29:38]), float(ln[38:46]), float(ln[46:54])
        occ = float(ln[55:60])
        mass = 1
        if id.lower() == "phenix":
            atom_mass = 0.0
            if len(ln) > 78:
                atom_mass = mass_of_atom(ln[76:78].strip())
            if atom_mass > 0:
                mass = atom_mass

        xt = xt + x * occ * mass
        yt = yt + y * occ * mass
        zt = zt + z * occ * mass
        occmass = occmass + occ * mass

    return [xt / occmass, yt / occmass, zt / occmass]


##########################################################


def get_xyz_tls(res, xyz):
    """input res (ch, range1, ch, range2); xyz coordinates."""
    nxyz = []
    for x in xyz[:]:
        if not atom_record(x):
            continue
        if res[0] == x[20:22].strip() and int(res[1]) <= int(x[22:26]) <= int(res[3]):
            nxyz.append(x)
    #        xyz.remove(x)
    return nxyz


##########################################################
def tls2cif(pdbfile):
    """
    convert TLS groups into mmcif (group number, component number
    are re-generated.
    1. extract TLS from PDB file (strainght forward)
    2. extract TLS from tls.inp (for refmac)

    First determine if it is a pdb file or a tlsanl.inp file,
    then extract the tls group and convert it to mmcif format
    """

    pdbid = "unknown"
    if not os.path.exists(pdbfile) or os.path.getsize(pdbfile) < 100:
        print("Error: pdb file (%s) not exist!" % pdbfile)
        return ""

    fp = open(pdbfile, "rU")
    id = "?"
    for x in fp:
        if "HEADER" in x[:6]:
            pdbid = x[62:66]
        elif "REMARK   3   TLS GROUP" in x[:22]:
            id = "pdb"
            break
    fp.close()

    out = pdbfile + "2cif"
    fw = open(out, "w")
    fw.write("data_%s\n" % pdbid.lower())

    if "pdb" in id:
        ntls, res, par, detail = pdbtls_2cif(pdbfile)
    else:
        ntls, res, par, detail = inptls_2cif(pdbfile)

    if ntls < 1:
        print("Error! no TLS groups are extracted. ")
        fw.close()
        return out

    resh = """
# 
loop_
_pdbx_refine_tls_group.pdbx_refine_id 
_pdbx_refine_tls_group.id 
_pdbx_refine_tls_group.refine_tls_id 
_pdbx_refine_tls_group.beg_auth_asym_id 
_pdbx_refine_tls_group.beg_auth_seq_id 
_pdbx_refine_tls_group.end_auth_asym_id 
_pdbx_refine_tls_group.end_auth_seq_id 
_pdbx_refine_tls_group.selection_details 
"""
    fw.write(resh)
    n = 0
    for i in range(1, ntls + 1):
        for y in res[i]:
            n = n + 1
            if "'" in detail[i]:
                t = '"X-RAY DIFFRACTION" %2d %2d %2s %4s %2s %4s "%s"\n' % (n, i, y[0], y[1], y[2], y[3], detail[i])
            else:
                t = "'X-RAY DIFFRACTION' %2d %2d %2s %4s %2s %4s '%s'\n" % (n, i, y[0], y[1], y[2], y[3], detail[i])
            fw.write(t)

    parh = """
# 
loop_
_pdbx_refine_tls.pdbx_refine_id 
_pdbx_refine_tls.id 
_pdbx_refine_tls.details 
_pdbx_refine_tls.method 
_pdbx_refine_tls.origin_x 
_pdbx_refine_tls.origin_y 
_pdbx_refine_tls.origin_z 
_pdbx_refine_tls.T[1][1] 
_pdbx_refine_tls.T[2][2] 
_pdbx_refine_tls.T[3][3] 
_pdbx_refine_tls.T[1][2] 
_pdbx_refine_tls.T[1][3] 
_pdbx_refine_tls.T[2][3] 
_pdbx_refine_tls.L[1][1] 
_pdbx_refine_tls.L[2][2] 
_pdbx_refine_tls.L[3][3] 
_pdbx_refine_tls.L[1][2] 
_pdbx_refine_tls.L[1][3] 
_pdbx_refine_tls.L[2][3] 
_pdbx_refine_tls.S[1][1] 
_pdbx_refine_tls.S[2][2] 
_pdbx_refine_tls.S[3][3] 
_pdbx_refine_tls.S[1][2] 
_pdbx_refine_tls.S[1][3] 
_pdbx_refine_tls.S[2][3] 
_pdbx_refine_tls.S[2][1] 
_pdbx_refine_tls.S[3][1] 
_pdbx_refine_tls.S[3][2] 
"""

    fw.write(parh)
    for i in range(1, ntls + 1):
        x = par[i]

        t0 = "'X-RAY DIFFRACTION' %2d ? refined %8s %8s %8s\n" % (i, x[0], x[1], x[2])
        t1 = "%7s %7s %7s %7s %7s %7s\n" % (x[3], x[4], x[5], x[6], x[7], x[8])
        t2 = "%7s %7s %7s %7s %7s %7s\n" % (x[9], x[10], x[11], x[12], x[13], x[14])
        t3 = "%7s %7s %7s %7s %7s %7s %7s %7s %7s\n" % (x[15], x[16], x[17], x[18], x[19], x[20], x[21], x[22], x[23])

        t = t0 + t1 + t2 + t3
        fw.write(t)
    fw.close()

    return out


##########################################################
def tls_selection_one_string(fp, i, nc):
    """fp: a list;  i: the position for SELECTION; nc, start column of range.
    return a string for this selection.
    """

    tls = fp[i].split(":", 1)[1].strip() + " "  # the first one
    #    tlsn=tls.strip().replace('(', ' ( ').replace(')', ' ) ')
    tlsn = ""
    ss = ""
    for n in range(i + 1, len(fp)):
        if "REMARK   3    ORIGIN FOR " in fp[n][:25]:  # stop here
            tls = tls + ss
            tlsn = tls.strip().replace("(", " ( ").replace(")", " ) ")
            break
        ss = ss + (fp[n][nc:].strip() + " ")

    return tlsn


##########################################################
def pdbtls_2cif(pdbfile):
    """extract TLS from PDB file and convert it to mmcif"""

    pdbid, prog = "????", "?"
    fp = open(pdbfile, "rU").readlines()

    detail, par, res, ntls, rs, p = {}, {}, {}, 0, [], []
    s11, s12, s13, s21, s22, s23 = "?" * 6

    chain, ch_range = chain_res(pdbfile, 1)

    for i, x in enumerate(fp):
        x = x.strip()
        if "HEADER" in x[:6]:
            pdbid = x[62:66]
        elif "ATOM" in x[:4] or "HETA" in x[:4]:
            break

        elif "REMARK   3   PROGRAM     :" in x:
            prog = x.split(":")[1].strip()

        elif "REMARK   3   TLS GROUP :" in x:
            ntls = ntls + 1

        elif "REMARK   3    RESIDUE RANGE :" in x[:29]:  # refmac
            t = x[29:].split()
            if len(t) != 4 or len(t[0]) != 1 or len(t[2]) != 1:
                print("Error in residue range: %s" % (x[29:].strip()))
                sys.exit()
            t1 = [t[0], t[1], t[2], t[3]]
            rs.append(t1)

            detail[ntls] = "?"

        elif "REMARK   3    SELECTION:" in x[:26]:  # phenix
            tls_str = tls_selection_one_string(fp, i, 25)
            if not tls_str:
                continue

            tlsg = "%d" % ntls
            t = parse_tls_range_phenix(tls_str, tlsg, chain, ch_range, fp)

            detail[ntls] = tls_str

            for z in t:
                t1 = [z[0], z[1], z[0], z[2]]
                rs.append(t1)

        elif "REMARK   3    SET :" in x[:26]:  # buster
            tls_str = tls_selection_one_string(fp, i, 20)
            tlsg = "%d" % ntls
            t = parse_tls_range_buster(tls_str, tlsg, chain, ch_range)
            detail[ntls] = tls_str

            for z in t:
                t1 = [z[0], z[1], z[0], z[2]]
                rs.append(t1)

        elif "REMARK   3    ORIGIN FOR THE GROUP" in x:
            xx, yy, zz = origin(x)
            p = [xx, yy, zz]

        elif "REMARK   3 " in x[:12] and (
            ("T11:" in x[14:21] and "T22:" in x[25:])
            or ("T33:" in x[14:21] and "T12:" in x[25:])
            or ("T13:" in x[14:21] and "T23:" in x[25:])
            or ("L11:" in x[14:21] and "L22:" in x[25:])
            or ("L33:" in x[14:21] and "L12:" in x[25:])
            or ("L13:" in x[14:21] and "L23:" in x[25:])
        ):
            tt = x.split(":")
            p.append(tt[1].split()[0])
            p.append(tt[2].split()[0])

        elif "REMARK   3" in x[:12] and "S11:" in x[14:21]:
            tt = x.split(":")
            s11, s12, s13 = tt[1].split()[0], tt[2].split()[0], tt[3].split()[0]

        elif "REMARK   3" in x[:12] and "S21:" in x[14:21]:
            tt = x.split(":")
            s21, s22, s23 = tt[1].split()[0], tt[2].split()[0], tt[3].split()[0]

        elif "REMARK   3" in x[:12] and "S31:" in x[14:21]:
            tt = x.split(":")
            s31, s32, s33 = tt[1].split()[0], tt[2].split()[0], tt[3].split()[0]

            t = [s11, s22, s33, s12, s13, s23, s21, s31, s32]
            p.extend(t)
            par[ntls] = p
            res[ntls] = rs
            p, rs = [], []
    return ntls, res, par, detail


##########################################################
def inptls_2cif(pdbfile):
    """extract TLS from tls.inp (for refmac), input should have the follows:
    TLS    chain D   (example)
    RANGE  'D  14.' 'D 447.' ALL (example)
    ORIGIN   -22.199  -23.538  -30.994     (example)
    T11, T22, T33, T12, T13, T23
    L11, L22, L33, L12, L13, L23
    S22_S11, S11_S33, S12, S13, S23, S21, S31, S32
    Since there are 8 independent Sij, the following equations are used to get Sij
    1). S11 + S22 + S33 = 0
    2). S22 - S11 = a     (a=S22_S11)
    3). S11 - S33 = b     (b=S11_S33)

    solving the equation:
    S11 = (b-a)/3.0
    S22 = (b+2*a)/3.0
    S33 = -(a+2*b)/3.0
    """

    fp = open(pdbfile, "r").readlines()

    detail, par, res, ntls, rs, p = {}, {}, {}, 0, [], []
    s11, s12, s13, s21, s22, s23 = "?" * 6

    for x in fp:
        y = x.replace("'", " ").split()
        if len(y) == 0:
            continue

        if "TLS" == y[0].upper():
            ntls = ntls + 1

        elif "RANGE" == y[0].upper():
            if "'" not in x:
                continue

            s = x.split("'")
            s1, s2 = s[1].split(), s[3].split()
            y1, y2, y3, y4 = "?", "?", "?", "?"
            if len(s1) == 1:
                y1, y2 = s1[0][0], s1[0][1:]
            else:
                y1, y2 = s1[0], s1[1]
            if len(s2) == 1:
                y3, y4 = s2[0][0], s2[0][1:]
            else:
                y3, y4 = s2[0], s2[1]

            a, b = y2.split("."), y4.split(".")
            print(s, y1, a[0], y3, b[0])

            t1 = [y1, a[0], y3, b[0]]
            rs.append(t1)
            detail[ntls] = "?"

        elif "ORIGIN" == y[0].upper():
            if len(y) == 4:
                p = [y[1], y[2], y[3]]

        elif "T" == y[0].upper() or "L" == y[0].upper():
            if len(y) == 7:
                t = [y[1], y[2], y[3], y[4], y[5], y[6]]
                p.extend(t)

        elif "S" == y[0].upper():
            # S22_S11, S11_S33, S12, S13, S23, S21, S31, S32 in TLSANL.inp

            a, b = float(y[1]), float(y[2])
            s11 = (b - a) / 3.0
            s22 = (b + 2 * a) / 3.0
            s33 = -(a + 2 * b) / 3.0

            s11 = "%8.4f " % s11
            s22 = "%8.4f " % s22
            s33 = "%8.4f " % s33

            t = [s11, s22, s33, y[3], y[4], y[5], y[6], y[7], y[8]]
            p.extend(t)
            par[ntls] = p
            res[ntls] = rs
            p, rs = [], []

    return ntls, res, par, detail


##########################################################
def extract_tls_pdb(pdbfile):
    """extract tls groups from PDB and keep the PDB format (remark 3)
    also extract TLS for refmac refinement
    """

    arg = "tlsextract xyzin %s tlsout tls4refmac.inp >/dev/null" % pdbfile
    os.system(arg)

    str1 = "REMARK   3  TLS DETAILS"
    str2 = "REMARK   3  BULK SOLVENT MODELLING."

    tlsp, pdbp = cut_between_pattern(pdbfile, str1, str2)

    if len(tlsp) > 6:
        out = pdbfile + ".tls"
        fw = open(out, "w")
        for x in tlsp:
            fw.write(x)
        fw.close()
        print("The extracted tls group = %s : file for refmac = tls4refmac.inp" % out)
    else:
        print("Error: No tls groups were extracted from pdb=%s" % pdbfile)

    delete_file(pdbp)


##########################################################
def get_file_by_pdbid(pdbid_in, id):
    """id==0, both pdb & sf;  id==1, only pdb; id==2, only sf;"""

    pdb_path = "/data/remediation-alt/ftp-v4.0/pdb/data/structures/all/pdb/"
    sf_path = "/data/remediation-alt/ftp-v4.0/pdb/data/structures/all/structure_factors/"
    www_path = "http://www.rcsb.org/pdb/files"

    pdbid = pdbid_in.lower()

    pdb = pdb_path + "pdb" + pdbid + ".ent.gz"
    sf = sf_path + "r" + pdbid + "sf.ent.gz"

    pdbfile = "pdb" + pdbid + ".ent"
    sffile = "r" + pdbid + "sf.ent"

    if id == 0:
        os.system("zcat  %s > %s " % (pdb, pdbfile))
        os.system("zcat  %s > %s " % (sf, sffile))
    elif id == 1:
        os.system("zcat  %s > %s " % (pdb, pdbfile))
    elif id == 2:
        os.system("zcat  %s > %s " % (sf, sffile))

    if (id == 1 and not os.path.exists(pdbfile)) or (id == 2 and not os.path.exists(sffile)):  # try wget
        pdbfile = pdbid + ".pdb"
        sffile = pdbid + "-sf.cif"
        if id == 0:
            os.system("wget %s/%s" % (www_path, pdbfile))
            os.system("wget %s/%s" % (www_path, sffile))
        elif id == 1:
            os.system("wget %s/%s" % (www_path, pdbfile))
        elif id == 2:
            os.system("wget %s/%s" % (www_path, sffile))

    if id == 0:
        return (pdbfile, sffile)
    elif id == 1:
        return pdbfile
    elif id == 2:
        return sffile


##########################################################
def recover_anisou(anis, pdbfile, prog, outfile):
    """If the PDB file already has B_full, but the ANISOU records are missing,
    validation still show 1~3% larger than the deposited. This function is to
    recover the ANISOU (1. subtract TLS from B_full; 2. get B_full with ANISOU.)
    """

    tls_inp = pdbfile + "_tls.inp"

    if "buster" in prog.lower():
        convert_tls2refmac(pdbfile, "buster")
        pdb_tmp = pdbfile + "_buster2ref"
        pdb_new = correct_tls(pdb_tmp)  # do various corrections

    elif "phenix" in prog.lower():
        convert_tls2refmac(pdbfile, "phenix")
        pdb_tmp = pdbfile + "_phenix2ref"
        pdb_new = correct_tls(pdb_tmp)  # do various corrections

    else:
        pdb_new = correct_tls(pdbfile)  # do various corrections

    tls_inp = get_tlsinp(pdb_new)  # extract TLS groups into tlsanl
    if not check_file(50, tls_inp):
        return

    if anis == 1:  # get ANISOU record from Biso_total
        (pdb_tls, log1) = do_tlsanl(pdbfile, tls_inp, "f", "resi")
        (pdb_out, log2) = do_tlsanl(pdb_tls, tls_inp, "t", "full")
        delete_file(pdb_tls, log2)

    elif anis == 2:  # get B_residual from Biso_total
        (pdb_out, log1) = do_tlsanl(pdbfile, tls_inp, "f", "resi")
        tmp = pdb_out + "__tmp"
        os.system('grep -v "^ANISOU" %s  >%s' % (pdb_out, tmp))
        os.system("mv %s %s" % (tmp, pdb_out))

    elif anis == 3:  # get B_TLS from Biso_total
        (pdb_out, log1) = do_tlsanl(pdbfile, tls_inp, "f", "tlsc")

    delete_file(tls_inp, pdb_new, log1)

    shutil.move(pdb_out, outfile)
    finalize_pdb(pdbfile, outfile)
    print("\nThe final output =%s\n" % outfile)


##########################################################
def proc_tlsanl(pdbfile, outfile):
    """run tlsanl program: correct residue ranges if they are wrong."""

    delete_file(outfile)
    print("%s" % (60 * "-"))
    (pdb_tls1, log1, pdb_tls2, log2, tls_inp, tls_new_inp, pdb_0tls, pdb_new, pdb_newb) = ("TMP.LOG",) * 9

    tls = precheck_tls(pdbfile)
    if tls == 2:  # phenix
        print("Note:  TLS format is phenix (%s)." % pdbfile)
        return 0
    elif tls == 3:  # Buster
        print("Note:  TLS format is Buster (%s)." % pdbfile)
        return 0

    tls_inp = get_tlsinp(pdbfile)
    pdb_tls1, log1 = "", ""
    # if tls==1: (pdb_tls1, log1) = do_tlsanl(pdbfile, tls_inp, 't', 'full')

    if os.path.exists(pdb_tls1) and os.path.getsize(pdb_tls1) > 500:
        shutil.move(pdb_tls1, outfile)
        print("\n%s; TLSANL was successful (residue-range is ok).\n" % pdbfile)
    else:
        pdb_0tls = delete_0_origin(pdbfile)  # remove tls group with 0 origin
        pdb_new = correct_tls(pdb_0tls)  # do various corrections
        tls_new_inp = get_tlsinp(pdb_new)  # extract TLS groups into tlsanl

        (pdb_tls2, log2) = do_tlsanl(pdb_new, tls_new_inp, "t", "full")

        if os.path.exists(pdb_tls2) and os.path.getsize(pdb_tls2) > 500:
            shutil.move(pdb_tls2, outfile)
            print("\n%s; TLSANL was successful (residue-range corrected).\n" % pdbfile)
        else:
            t = "Error! TLS correction failed (%s).\n" % pdbfile
            perror(t)
            if os.path.exists(log2) and os.path.getsize(log2) > 100:
                fp = open(log2, "r").readlines()
                n = 0
                for ln in fp:
                    n = n + 1
                    if "ERROR:" in ln:
                        print(pdbfile + "; " + ln)
                        if len(fp[n].strip()) > 2:
                            print(pdbfile + "; " + fp[n])
                        break

    wk = 0
    if os.path.exists(outfile) and os.path.getsize(outfile) > 100:
        bigb, bval = check_bfactor(outfile)  # only show Biso status
        finalize_pdb(pdbfile, outfile)
        cif = tls2cif(outfile)  # get a tls cif file
        #        print('The corrected TLS in mmCIF format = %s\n' %cif)
        #        print('The pdbfile after tlsanl correction = %s\n' %outfile)
        wk = 1
        if bigb > 0:
            wk = 0
            # delete_file(outfile)

    delete_file(pdb_tls1, log1, pdb_tls2, log2, pdb_0tls, pdb_new, pdb_newb)
    delete_file(tls_inp, tls_new_inp)

    return wk


##########################################################
def parse_tls_range_refmac(x, ntls):
    """parse refmac tls range; return (ch,nres1,nres2) return none if bad!"""

    res = []
    tmp = x.split(":", 1)[1].strip().split()
    if len(tmp) != 4 or tmp[0] != tmp[2] or " NULL " in x or len(str(tmp[0])) > 1 or (is_number(tmp[1]) and is_number(tmp[3]) and float(tmp[3]) < float(tmp[1])):
        t = "Error: Wrong residue range for tls group=%s: (%s)" % (ntls, x[29:].strip())
        perror(t)

    elif not is_number(tmp[1]) or not is_number(tmp[3]):
        print("Note: maybe inserted code (not treated currently):", tmp)

    else:
        a, b = tmp[1].split("."), tmp[3].split(".")
        res = [tmp[0], int(a[0]), int(b[0])]
    #        res=[tmp[0], int(tmp[1]), int(tmp[3])]

    return res


##########################################################
def make_ph_str(resid):
    """write the dic into a string for the phenix format."""

    str1 = ""
    kk = len(resid)
    for k, x in enumerate(resid.keys()):
        str1 = str1 + "(CHAIN %s AND " % x
        nn = len(resid[x])
        if nn == 1:
            t = resid[x][0]
            str1 = str1 + "RESID %d:%d)" % (t[0], t[1])
        else:
            str1 = str1 + "("
            for i, y in enumerate(resid[x]):
                str1 = str1 + "RESID %d:%d" % (y[0], y[1])
                if i != nn - 1:
                    str1 = str1 + " OR "
            str1 = str1 + "))"

        if k != kk - 1:
            str1 = str1 + " OR "

    return str1


##########################################################
def make_tls_str(resid, id):
    """write the dic into a string for the ID format.
    for buster & phenix
    """

    str1 = ""
    kk = len(resid)
    if id == "buster":
        str1 = "{ "

        for k, v in resid.items():
            for y in v:
                str1 = str1 + "%s|%d - %s|%d " % (k, y[0], k, y[1])

        str1 = str1 + "}"
        return str1

    for k, x in enumerate(resid.keys()):
        str1 = str1 + "(CHAIN %s AND " % x
        nn = len(resid[x])
        if nn == 1:
            t = resid[x][0]
            str1 = str1 + "RESID %d:%d)" % (t[0], t[1])
        else:
            str1 = str1 + "("
            for i, y in enumerate(resid[x]):
                str1 = str1 + "RESID %d:%d" % (y[0], y[1])
                if i != nn - 1:
                    str1 = str1 + " OR "
            str1 = str1 + "))"

        if k != kk - 1:
            str1 = str1 + " OR "

    return str1


##########################################################
def make_ph_str_pdb1(str1):
    """write one line string into multiple lines of PDB format
    words will be broken.
    """

    s1 = "REMARK   3    SELECTION: "
    s2 = "REMARK   3               "

    str2 = s1
    n = 0
    for i, x in enumerate(str1):
        if i == 0:
            n = 1
            str2 = s1 + x

        elif i == 54 * n:
            str2 = str2 + "\n" + s2
            n = n + 1
        else:
            str2 = str2 + x

    return str2 + "\n"


##########################################################
def make_ph_str_pdb(str1):
    """write one line string into multiple lines of PDB format
    words are not broken.
    """

    s1 = "REMARK   3    SELECTION: "
    s2 = "REMARK   3               "

    str2 = s1
    ss = str1.split()
    n, tmps = 0, ""

    for i, x in enumerate(ss):
        tmps = tmps + " " + x

        if i == 0:
            n = 1
            str2 = s1 + x

        elif len(tmps) >= 50 * n:
            str2 = str2 + " " + x + "\n" + s2
            n = n + 1
        else:
            str2 = str2 + " " + x

    return str2 + "\n"


##########################################################
def make_tls_str_pdb(str1, id):
    """write one line string into multiple lines of PDB format
    words are not broken. (only for phenix & buster)
    """

    s1 = "REMARK   3    SELECTION: "  # phenix
    s2 = "REMARK   3               "

    if id == "buster":
        s1 = "REMARK   3    SET : "  # buster
        s2 = "REMARK   3          "

    str2 = s1
    ss = str1.split()
    n, tmps = 0, ""

    for i, x in enumerate(ss):
        tmps = tmps + " " + x

        if i == 0:
            n = 1
            str2 = s1 + x

        elif len(tmps) >= 54 * n:
            str2 = str2 + " " + x + "\n" + s2
            n = n + 1
        else:
            str2 = str2 + " " + x

    return str2 + "\n"


##########################################################
def map_tls_4_phenix(pdbfile, pdbdep, outfile):
    print("Mapping tls using %s %s " % (pdbfile, pdbdep))
    #    pdb_ref = convert_tls2refmac(pdbfile, 'phenix')
    dep_ref = convert_tls2refmac(pdbdep, "phenix")

    #    pdb_sort = reorder_pdb(pdb_ref)
    wk, pdb_match = match_all_residues(pdbfile, dep_ref)

    if wk == 1:
        outfile_tmp = outfile + "_tmp"
        proc_tlsanl(pdb_match, outfile_tmp)

        ref2ph = convert_refmac2phenix(outfile_tmp)
        cif = tls2cif(ref2ph)

        delete_file(outfile_tmp)

    if os.path.exists(ref2ph) and os.path.exists(cif):
        tlspdb = outfile + ".tlspdb"
        shutil.move(ref2ph, tlspdb)
        shutil.move(cif, outfile)

        print("The new output PDB file = %s" % tlspdb)
        print("The new TLS groups in cif = %s\n" % outfile)
    else:
        t = "Warning: TLS conversion to cif failed (%s)." % pdbfile
        perror(t)

    delete_file(dep_ref, pdb_match, cif)


##########################################################
def map_tls_4_buster(pdbfile, pdbdep, outfile):
    print("Mapping (buster) tls using %s %s " % (pdbfile, pdbdep))

    dep_ref = convert_tls2refmac(pdbdep, "buster")

    #    pdb_sort = reorder_pdb(dep_ref)
    #    wk, pdb_match = match_all_residues(pdbfile, pdb_sort)
    wk, pdb_match = match_all_residues(pdbfile, dep_ref)

    ref2bus = convert_refmac2buster(pdb_match)
    cif = tls2cif(ref2bus)

    if os.path.exists(ref2bus) and os.path.exists(cif):
        tlspdb = outfile + ".tlspdb"
        shutil.move(ref2bus, tlspdb)
        shutil.move(cif, outfile)

        print("The new output PDB file = %s" % tlspdb)
        print("The new TLS groups in cif = %s\n" % outfile)
    else:
        t = "Warning: TLS conversion to cif failed (%s)." % pdbfile
        perror(t)

    # delete_file(dep_ref, pdb_match, cif)


##########################################################
def convert_refmac2phenix(pdbfile):
    """Parse tls residue ranges of refmac output and convert them into phenix."""

    out = pdbfile + "_ref2ph"
    if not check_file(50, pdbfile):
        return out

    fw, fp = open(out, "w"), open(pdbfile, "r").readlines()

    resid, ntls = {}, "?"
    for i, x in enumerate(fp):
        if "TLS GROUP :" in x[11:25]:
            ntls = x.split(":")[1].strip()
            fw.write(x)

        elif "NUMBER OF COMPONENTS GROUP" in x[12:41] or "C SSSEQI" in x[31:41] or (len(x.strip()) < 12 and i < len(fp) - 1 and "TLS GROUP :" in fp[i + 1]):
            continue

        elif "RESIDUE RANGE :" in x[13:30]:
            res = parse_tls_range_refmac(x, ntls)
            if not res:
                print("Warning: Skipping residue range (%s) for TLS group (%s)" % (ntls, x[30:].strip()))
                continue
            ch, rrange = res[0], [res[1], res[2]]
            if ch in resid.keys():
                resid[ch].append(rrange)
            else:
                resid[ch] = [rrange]

        elif " ORIGIN FOR THE GROUP (A):" in x[12:40]:
            str1 = make_ph_str(resid)
            str2 = make_ph_str_pdb(str1)
            resid = {}
            fw.write(str2)
            fw.write(x)
        else:
            fw.write(x)

    fw.close()

    return out


##########################################################
def convert_refmac2buster(pdbfile):
    """Parse tls residue ranges of refmac output and convert them into buster."""

    out = pdbfile + "_ref2bus"
    fw, fp = open(out, "w"), open(pdbfile, "r").readlines()

    resid, ntls = {}, "?"
    for i, x in enumerate(fp):
        if "TLS GROUP :" in x[11:25]:
            ntls = x.split(":")[1].strip()
            fw.write(x)

        elif "NUMBER OF COMPONENTS GROUP" in x[12:41] or "C SSSEQI" in x[31:41] or (len(x.strip()) < 12 and i < len(fp) - 1 and "TLS GROUP :" in fp[i + 1]):
            continue

        elif "RESIDUE RANGE :" in x[13:30]:
            res = parse_tls_range_refmac(x, ntls)
            if not res:
                print("Warning: Skipping residue range (%s) for TLS group (%s)" % (ntls, x[30:].strip()))
                continue
            ch, rrange = res[0], [res[1], res[2]]
            if ch in resid.keys():
                resid[ch].append(rrange)
            else:
                resid[ch] = [rrange]

        elif " ORIGIN FOR THE GROUP (A):" in x[12:40]:
            str1 = make_tls_str(resid, "buster")
            str2 = make_tls_str_pdb(str1, "buster")
            resid = {}
            fw.write(str2)
            fw.write(x)
        else:
            fw.write(x)

    fw.close()
    print("\nTLS is converted to buster format. Output=%s\n" % out)

    return out


##########################################################
def convert_tls2refmac(pdbfile, id):
    """Parse tls residue ranges and convert them into refmac.
    id=phenix:  from phenix to refmac.
    id=buster:  from buster to refmac.
    """

    ref0 = "REMARK   3    NUMBER OF COMPONENTS GROUP : 1   \n"
    ref1 = "REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI    \n"
    ref2 = "REMARK   3    RESIDUE RANGE : "

    if id == "phenix":
        str1 = "REMARK   3    SELECTION: "
    elif id == "buster":
        str1 = "REMARK   3    SET : "
    str2 = "REMARK   3    ORIGIN FOR "

    chain, ch_range = chain_res(pdbfile, 1)
    tls, tlsg, ntls, nn = "", 0, 0, 0

    pdbnew = pdbfile + "_%s2ref" % id
    fw = open(pdbnew, "w")
    fp = open(pdbfile, "r")
    fpp = open(pdbfile, "r").readlines()
    for x in fp:
        if " NUMBER OF TLS GROUPS  :" in x[11:37]:
            tt = x.split(":", 1)[1].strip()
            if is_number(tt):
                ntls = int(tt)
                fw.write(x)
            else:
                print("Error: (%s) not integer after  TLS GROUP :" % tt)

        elif " TLS GROUP :" in x[11:24]:
            tt = x.split(":", 1)[1].strip()
            if is_number(tt):
                tlsg = int(tt)
                fw.write(x)
            else:
                print("Error: (%s) not integer after  TLS GROUP :" % tt)

        elif str1 in x or (id == "buster" and "SELECTION:" in x):
            fw.write(ref0)
            fw.write(ref1)

            if id == "buster" and "SELECTION:" in x and "{" not in x:  # refmac style
                t = parse_tls_range_refmac(x, ntls)
                tls = tls + t[0] + "|" + "%d" % t[1] + " - " + t[0] + "|" + "%d" % t[2]
            else:
                tls = tls + (x.split(":", 1)[1].strip() + " ")

            #            print 'x1=', x
            for y in fp:
                #                print 'x=', y
                if "CRYST1" in y[:6]:
                    break
                if str2 in y:  # write in refmac foramt
                    if id == "phenix":
                        tlsn = tls.strip().replace("(", " ( ").replace(")", " ) ")
                        tls_range = parse_tls_range_phenix(tlsn, tlsg, chain, ch_range, fpp)

                    elif id == "buster":
                        tlsn = tls.strip().replace("{", " ").replace("}", " ")
                        tls_range = parse_tls_range_buster(tlsn, tlsg, chain, ch_range)

                    for r in tls_range:
                        t = ref2 + "%3s  %4d %8s  %4d\n" % (r[0], r[1], r[0], r[2])
                        fw.write(t)

                    tls = ""
                    # fw.write(y)
                    nn = tlsg

                    t = (y.split(":", 1)[1]).split()
                    if len(t) != 3:
                        print("Warning: connected values in %s" % x)
                        fw.write(y)
                    else:
                        tt = "REMARK   3    ORIGIN FOR THE GROUP (A):%9s%9s%9s \n" % (t[0], t[1], t[2])
                        fw.write(tt)

                    break

                if id == "phenix":
                    tls = tls + (y[24:].strip() + " ")
                elif id == "buster" and "SET:" in y:
                    tls = tls + (y[19:].strip() + " ")
                    print("tls=%s" % tls)

                elif id == "buster" and "SELECTION:" in y:
                    t = parse_tls_range_refmac(y, ntls)
                    tls = tls + " " + t[0] + "|" + "%d" % t[1] + " - " + t[0] + "|" + "%d" % t[2]

        elif "REMARK   3    ORIGIN FOR THE GROUP (A):" in x:
            t = (x.split(":", 1)[1]).split()
            if len(t) != 3:
                print("Warning: connected values in %s" % x)
                fw.write(x)
            else:
                tt = "REMARK   3    ORIGIN FOR THE GROUP (A):%9s%9s%9s \n" % (t[0], t[1], t[2])
                fw.write(tt)

        else:
            fw.write(x)
            if nn > 0 and nn == ntls and " S31:" in x[13:20]:
                fw.write("REMARK   3                           \n")
                fw.write("REMARK   3  BULK SOLVENT MODELLING. \n")
                nn = 0

    fw.close(), fp.close()
    print("\nTLS is converted to refmac format. Output = %s\n" % pdbnew)
    return pdbnew


##########################################################
def parse_tls_range_buster(tls, ng, chain, chain_range):
    """tls: all the tls range as one string.  ng: the tls group number
    chain:  chainID with a list of residue number.
    chain_range: chainID with tls residual range.

    """

    tls_range = []

    tlsn = tls.replace(" |", "|").replace("| ", "|").replace("{", "").replace("}", "")

    if "|" not in tlsn or "null" in tlsn.lower():
        t = "Error: problem with TLS group=%s. \n" % (ng)
        perror(t)

    elif " - " not in tlsn:  # one residue or all in one chain
        st = tlsn.split()

        for x in st:
            if "|" in x:
                tmp = parse_tls_range_buster_1(x, chain_range)
                if x[0] == "*":
                    tls_range = tmp
                else:
                    tls_range.append(tmp)

    else:  # multiple ranges
        ss = tlsn.split()
        i = 0
        nn = len(ss)
        for k in range(i, nn):
            if "|" in ss[k] and "*" in ss[k] and nn - k > 1 and "-" not in ss[k + 1]:  # isolate entity
                tmp = parse_tls_range_buster_1(ss[k], chain_range)
                if tmp:
                    tls_range.append(tmp)
                else:
                    t = "Error: residue range (%s) not parsed (tls group=%s)" % (ss[k], ng)
                    perror(t)

            elif "|" in ss[k] and nn - k > 1 and "-" not in ss[k + 1]:  # isolate entity
                tmp = parse_tls_range_buster_1(ss[k], chain_range)

            elif "|" in ss[k] and nn - k > 2 and "-" in ss[k + 1]:
                ch1, res1 = parse_tls_range_buster_2(ss[k])

                if "|" not in ss[k + 2]:
                    if is_number(res1) and is_number(ss[k + 2]):
                        tls_range.append([ch1, int(res1), int(ss[k + 2])])
                    else:
                        t = "Error: residue range (%s - %s) not parsed (tls group=%s)" % (ss[k], ss[k + 2], ng)
                        perror(t)

                else:
                    ch2, res2 = parse_tls_range_buster_2(ss[k + 2])
                    if ch1 == ch2 and is_number(res1) and is_number(res2):
                        tls_range.append([ch1, int(res1), int(res2)])
                    else:
                        print("Warning: parse problem for TLS group=%s" % ng)

                i = i + 2

            i = i + 1

    check_parsed_tls_ranges(tls_range, chain, ng)

    return tls_range


##########################################################
def parse_tls_range_buster_2(x):
    s = x.split("|")
    ch, res = s[0], s[1]
    return ch, res


##########################################################
def parse_tls_range_buster_1(x, chain_range):
    """parse single items  (A|*,  A|45, A|12-45, A|12-A|45)"""

    tmp = []
    ch, res = parse_tls_range_buster_2(x)
    if "*" in ch and "*" in res:  # all chains
        for k, v in chain_range.items():
            tmp.append([k, v[0], v[1]])
        return tmp

    if ch not in chain_range.keys():
        t = "Error: chain ID (%s) not in coordinates " % (ch)
        perror(t)
        return tmp

    if "*" in x:
        tmp = [ch, chain_range[ch][0], chain_range[ch][1]]

    elif is_number(res):
        tmp = [ch, int(res), int(res)]

    else:  # for connected case (A|12-45, A|12-A|45)
        tt = x.split("|")
        if "-" in res and len(tt) == 2:  # A|12-45
            t = res.split("-")
            if len(t) != 2:
                return tmp
            t1, t2 = t[0], t[1]
            if is_number(t1) and is_number(t2):
                tmp = [ch, int(t1), int(t2)]

        elif "-" in res and len(tt) == 3:  # A|12-A|45
            t1 = x.split("-")
            if len(t1) != 2:
                return tmp

            if "|" in t1[0] and "|" in t1[1]:
                res1, res2 = t1[0].split("|"), t1[1].split("|")

                if is_number(res1[1]) and is_number(res2[1]) and res1[0] == res2[0]:
                    tmp = [ch, int(res1[1]), int(res2[1])]
                else:
                    t = "Error: The first chain ID (%s) differs from the second(%s)" % (res1[0], x)
                    perror(t)

    return tmp


##########################################################
def check_parsed_tls_ranges(tls_range, chain, ng):
    if not tls_range:
        return
    if not len(tls_range):
        t = "Error: failed to parse residue range (TLS group=%s)." % ng
        perror(t)
    #    elif ' not ' in tls_range:
    #        return
    else:  # check if residue exist
        for x in tls_range:
            if not x:
                continue
            k1, v1, v2 = x[0], x[1], x[2]

            if k1 not in chain.keys():
                t = "Error: chain ID (%s) not in coordinates (tls group=%s)" % (k1, ng)
                perror(t)
            else:
                if v1 not in chain[k1]:
                    t = "Warning: residue (%s_%d) not in coordinates (tls group=%s)" % (k1, v1, ng)
                    perror(t)
                if v1 != v2 and v2 not in chain[k1]:
                    t = "Warning: residue (%s_%d) not in coordinates (tls group=%s)" % (k1, v2, ng)
                    perror(t)


##########################################################
def keep_digit(x):
    """only keep the digit and '-' and ':' from x (string)"""

    t = "".join(filter(lambda y: y.isdigit() or y == ":" or y == "-", x))
    return t


##########################################################
def myreplace_case(before, after, ss):
    """before: string to be replaced by after. ss is the string"""
    beforein = "(?i)%s" % before
    return re.sub(beforein, after, ss)


##########################################################
def get_residue_range_from_resname(ch, res, fp):
    range = []
    for x in fp:
        if ("ATOM" in x[:4] or "HETA" in x[:4]) and (ch == x[20:22].strip() and res == x[17:20].lower()):
            range.append(x[22:26])
    if not range:
        perror("\nError: the chainid/resname (%s %s) is not in the coordinate." % (ch, res))
        return range

    nres1, nres2 = int(range[0]), int(range[-1])
    range1 = [ch, nres1, nres2]

    return range1


##########################################################
def parse_tls_range_phenix_special(tls, chain, chain_range, fp):
    """It involves in resname : for temp use. Modify later!"""

    range1 = []
    for i, x in enumerate(tls):
        if pattern(i, tls, "chain", "?", "and", "resname"):
            ch, resname = tls[2], tls[i + 4]
            range_t = get_residue_range_from_resname(ch, resname, fp)
            range1.append(range_t)
        #            print 'Got it  %s %s ' %(ch, resname), tls
        elif pattern(i, tls, "(", "chain", "?", "or"):
            chs = []
            tlstmp = tls[i:]
            for j, y in enumerate(tlstmp):
                if "chain" in y:
                    chs.append(tlstmp[j + 1])
                if "resname" in y and "(" not in tlstmp[j + 1]:
                    resname = tlstmp[j + 1]
                    for c in chs:
                        range_t = get_residue_range_from_resname(c, resname, fp)
                        range1.append(range_t)
                    break
    return range1


##########################################################
def parse_tls_range_phenix(tlsn, ng, chain, chain_range, fp):
    """tlsn: all the tls range as one string. ng: the tls group number
    chain:  chainID with a list of residue number.
    chain_range: chainID with tls residual range.
    fp: the list of coordinate
    """

    if len(tlsn.split("(")) != len(tlsn.split(")")):
        t = "\nError:  parentheses is not closed for TLS group = %s" % ng
        perror(t)

    tls, tls_range = [], []
    s0 = " ".join(tlsn.strip().replace("'", "").replace('"', "").split())
    s1 = myreplace_case("resseq", "resid", s0)
    s2 = myreplace_case(" through ", ":", s1)
    tmp = s2.replace(" : ", ":").replace(": ", ":").replace(" :", ":").split()
    #    print tmp
    for i, x in enumerate(tmp):
        t = x
        if x.endswith(":") and i < len(tmp) - 1 and is_number(tmp[i + 1]):
            t = x + tmp[i + 1]
        if i > 0 and is_number(x) and tmp[i - 1].endswith(":"):
            continue

        if i > 0 and is_number(x) and "resid" in tmp[i - 1]:
            tls.append("%s:%s" % (x, x))
            continue
        if t and len(t) != 1 and len(t) != 2 and not is_number(t[0]) and t[0] != ":":
            tls.append(t.lower())
        else:  # chid
            tls.append(t)
    #    print 'tmp=', s2, tmp, tls

    if "resname" in tls:
        spec_range = parse_tls_range_phenix_special(tls, chain, chain_range, fp)
        print("Special range = %s (tls group=%s)" % (spec_range, ng))
        return spec_range

    nn = len(tls)
    for i, x in enumerate(tls):  # error detect
        #        print 'helll=', i, x
        if "null" in x:
            t = "Error: NULL values for residue range with TLS group = %s" % ng
            perror(t)

        elif "chain" in x:
            if i > nn - 2:
                t = "Error: No chainID is given (tls group=%s)" % (ng)
                perror(t)

            else:
                chid = tls[i + 1]
                if len(tls[i + 1]) > 2:
                    t = "Error: chainID (%s) is wrong (tls group=%s)" % (tls[i + 1], ng)
                    perror(t)
                else:
                    if chid not in chain_range.keys():
                        t = "Error: chain ID (%s) not in coordinates (tls group=%s)" % (chid, ng)
                        perror(t)

        elif "segid" in x and i < nn - 1:
            chid = tls[i + 1]
            if len(tls[i + 1]) > 1:
                t = "Error: segidID (%s) is used (tls group=%s)" % (tls[i + 1], ng)
                perror(t)

    for i, x in enumerate(tls):
        if "all" in x:
            for y in chain_range.keys():
                tls_range.append([y, chain_range[y][0], chain_range[y][1]])

        elif pattern(i, tls, "chain", "?", "and", "not", "resid"):
            ch = tls[i + 1]
            if ch in chain.keys():
                tt = parse_tls_phenix_not(i, ch, tls, chain[ch], chain_range, 0)
                tls_range.extend(tt)

        elif pattern(i, tls, "(", "chain", "?", ")", "and", "not", "(", "resid"):
            ch = tls[i + 2]
            #            print 'I am here', i,ch, tls
            if ch in chain.keys():
                tt = parse_tls_phenix_not(i, ch, tls, chain[ch], chain_range, 1)
                tls_range.extend(tt)

        elif pattern(i, tls, "chain", "?", "and", "not", "(", "resid"):
            ch = tls[i + 1]
            if ch in chain.keys():
                tt = parse_tls_phenix_not(i, ch, tls, chain[ch], chain_range, 1)
                tls_range.extend(tt)

        elif pattern(i, tls, "chain", "?", "and", "(", "not", "("):
            ch = tls[i + 1]
            if ch in chain.keys():
                tt = parse_tls_phenix_not(i, ch, tls, chain[ch], chain_range, 1)
                tls_range.extend(tt)

        #        elif pattern(i, tls, 'chain', '?', 'and', '(', 'resname' , '?') :

        elif "chain" in x and i < nn - 1:
            chid = tls[i + 1]
            if chid not in chain_range.keys():
                continue

            #            elif (pattern(i, tls, 'chain', '?') and (nn==2 or nn==4)) : #one chain
            elif (
                pattern(i, tls, "chain", "?", "or")
                or pattern(i, tls, "chain", "?", ")", "or")
                or (pattern(i, tls, "chain", "?") and i > nn - 3 and nn > 5)
                or (pattern(i, tls, "chain", "?") and nn == 2)
                or pattern(i, tls, "chain", "?", ")")
            ):  # one chain
                n1, n2 = chain_range[chid][0], chain_range[chid][1]
                tls_range.append([chid, n1, n2])

            elif pattern(i, tls, "chain", "?", "or", "chain", "?"):  # two chain
                n1, n2 = chain_range[chid][0], chain_range[chid][1]
                tls_range.append([chid, n1, n2])

                ch2 = tls[i + 4]
                if ch2 not in chain_range.keys():
                    continue
                n1, n2 = chain_range[ch2][0], chain_range[ch2][1]
                tls_range.append([ch2, n1, n2])

            elif i < nn - 4 and (pattern(i, tls, "chain", "?", "and", "res")):
                n1, n2 = get_range_ph(chain_range[chid], tls[i + 4])
                tls_range.append([chid, n1, n2])

            elif pattern(i, tls, "chain", "?", "and", "(", "res"):
                get_range_phe(ng, chid, i + 4, 1, chain_range, tls, tls_range)

            elif pattern(i, tls, "chain", "?", "and", "(", "(", "res"):
                get_range_phe(ng, chid, i + 5, 2, chain_range, tls, tls_range)

        elif i == 0 and pattern(i, tls, "res", "?", "and", "chain", "?"):
            chid = tls[i + 4]
            if chid not in chain_range.keys():
                continue
            n1, n2 = get_range_ph(chain_range[chid], tls[i + 1])
            tls_range.append([chid, n1, n2])

        elif i == 0 and pattern(i, tls, "(", "res", "?", ")", "and", "chain", "?"):
            chid = tls[i + 6]
            if chid not in chain_range.keys():
                continue

            n1, n2 = get_range_ph(chain_range[chid], tls[i + 2])
            tls_range.append([chid, n1, n2])

        elif i == 0 and pattern(i, tls, "(", "res", "?", "or", "res", "?") and pattern(8, tls, "chain", "?"):
            chid = tls[tls.index("chain") + 1]
            if chid not in chain_range.keys():
                continue
            get_range_phe(ng, chid, i + 0, 1, chain_range, tls, tls_range)

        elif pattern(i, tls, "res") and "chain" not in tls:
            if len(chain_range.keys()) > 1:
                t = "Error: wrong residue selection, no chain identifer (tls group=%s)." % (ng)
                perror(t)
            else:
                chid = list(chain_range.keys())[0]
                n1, n2 = get_range_ph(chain_range[chid], tls[i + 1])
                tls_range.append([chid, n1, n2])

    check_parsed_tls_ranges(tls_range, chain, ng)

    return tls_range


##########################################################
def pattern(i, tls, *ss):
    """the pattern of TLS after position i.  ss is a list of input"""

    key, nn = 1, len(tls)

    if len(ss) > nn - i:
        return 0
    for k, x in enumerate(ss):
        n = i + k
        if n > nn - 1 or x == "?" or (x == "res" and x == tls[n][:3]):
            continue
        if x != tls[n]:
            key = 0
            break

    return key


##########################################################
def parse_tls_phenix_not(i, ch, tls, chain, chain_range, id):
    """parse tls range involve NOT;
    id=0 for style 'chain F and not resid 0'
    id=1 for style 'chain F and not  not (resseq 539:581 or ..)'
    """

    chain_range_not = []
    chain_range_yes = []
    chain_new = chain

    if id == 0:
        n1, n2 = get_range_ph(chain_range, tls[i + 5])
        chain_range_not.append([n1, n2])

    else:
        m, tmp_list = get_parenthetic_contents(i + 4, tls)
        for i, x in enumerate(tmp_list):
            if "resid" in x:
                n1, n2 = get_range_ph(chain_range, tmp_list[i + 1])
                chain_range_not.append([n1, n2])

    chain_range_not.sort()
    nn1 = len(chain_range_not)
    chain_new.sort()
    nn2 = len(chain_new)
    # print('here=', nn1, chain_range_not, nn2, chain_new)
    for i, x in enumerate(chain_range_not):
        if i == 0:  # first one
            if x[0] - chain_new[0] > 2:
                chain_range_yes.append([ch, chain_new[0], x[0] - 1])
        if i == nn1 - 1:  # last one
            if chain_new[nn2 - 1] - x[1] > 2:
                chain_range_yes.append([ch, x[1] + 1, chain_new[nn2 - 1]])

        if i < nn1 - 1 and chain_range_not[i + 1][0] - x[1] > 2:
            chain_range_yes.append([ch, x[1] + 1, chain_range_not[i + 1][0] - 1])

    #    print('chain_range_not=',chain_range_not)
    #    print('chain_range_yes=',chain_range_yes)
    return chain_range_yes


##########################################################
def get_parenthetic_contents(start, lists):
    """start:  the start position;  lists: input a list or string
    return the position of last ')' and the contents in ()
    """

    nn = len(lists)
    tmp_list = []
    m1, m2, n = 0, 0, 0
    for i in range(start, nn):
        if lists[i] == "(":
            m1 = m1 + 1
        if lists[i] == ")":
            m2 = m2 + 1
        if type(lists) == type(list()):  # a list
            tmp_list.append(lists[i])

        if m1 > 0 and m1 == m2:
            if type(lists) == type(list()):  # a list
                return n, tmp_list
            elif type(lists) == type(str()):  # a str
                return n, lists[start:n]
            break
        n = n + 1

    if m1 != m2:
        t = "Error! Parenthese not equal in line %s \n" " ".join(lists)
        perror(t)
    return 0, tmp_list


##########################################################
def get_range_phe(ng, chid, st, nb, chain_range, tls, tls_range):
    """st: start point;  nb: number of bracket;  chain_range: {ch:[n1,n2]}
    tls_range: return the value
    """
    nn = len(tls)
    kk = 0
    for m in range(st, nn):
        if (nb == 1 and tls[m] == ")") or (nb == 2 and tls[m] == ")" and m + 1 < nn and tls[m + 1] == ")"):
            kk = 1
            break
        if ":" in tls[m]:
            x = keep_digit(tls[m])  # strip the insertion code.
            n1, n2 = get_range_ph(chain_range[chid], x)
            tls_range.append([chid, n1, n2])

        elif tls[m][0] != "-" and "-" in tls[m]:
            t = tls[m].split("-")
            n1, n2 = int(t[0]), int(t[1])
            tls_range.append([chid, n1, n2])

    if kk == 0:
        ss = "Error: residue range is possibly truncated (TLS group=%s)!" % ng
        perror(ss)


##########################################################
def get_range_ph(chain_range, tls):
    """tls= resid n1:n2 (parse n1, and n2)"""
    n1, n2 = -9999, -9999

    if ":" in tls:
        t = tls.split(":")
        if len(t) == 2:
            #            n1, n2 = int(t[0]), int(t[1])
            if is_number(t[0]):
                n1 = int(t[0])
            if is_number(t[1]):
                n2 = int(t[1])
    else:
        if "-" in tls and len(tls.split("-")) == 2:
            t = tls.split("-")
            if is_number(t[0]):
                n1 = int(t[0])
            if is_number(t[1]):
                n2 = int(t[1])
        else:
            if is_number(tls):
                n1 = n2 = int(tls)

    if n1 == -9999:
        n1 = chain_range[0]
    if n2 == -9999:
        n2 = chain_range[1]

    return n1, n2


##########################################################
def parenthetic_contents(string):
    """Generate parenthesized contents in string as pairs (level, contents)."""
    stack = []
    for i, c in enumerate(string):
        if c == "(":
            stack.append(i)
        elif c == ")" and stack:
            start = stack.pop()
            yield (len(stack), string[start + 1 : i])


##########################################################
def precheck_tls(pdbfile):
    """pre-check the type of TLS and possible errors"""

    wk = 1  # default refmac

    fr = open(pdbfile, "r")
    for x in fr:
        if "REMARK   3    RESIDUE RANGE : " in x:
            tmp = x.split(":")[1].split()
            if len(tmp) != 4 or is_digit(tmp[1]) == 0 or is_digit(tmp[3]) == 0 or (int(tmp[1]) < 0 and (int(tmp[3]) == -1 or int(tmp[3]) >= 999)):
                wk = 0
                break

        elif "REMARK" in x[:6] and " SELECTION: " in x[13:25]:
            wk = 2
            break

        elif "REMARK" in x[:6] and " SET : " in x[12:20]:
            wk = 3
            break

        elif "ATOM" in x[:4] or "HETATM" in x[:6]:
            break

    fr.close()

    return wk


##########################################################
def finalize_pdb(pdbfile, outfile):
    """this step is only to replace REMARK 3 and coordinates in the pdbfile
    from outfile, so that nothing changes unless the modified parts.
    This function would be skiped if the TLSANL program does not change other
    things (such as  KEYWDS in, REMARK 500, HETNAM)  the original PDB file.
    """

    # print('Finalizing the output file')

    str1 = "REMARK   3 REFINEMENT."
    str2 = "REMARK   3  OTHER REFINEMENT REMARKS:"
    str3 = "SCALE3   "
    str4 = "MASTER"

    m = check_string_inpdb(pdbfile, str1, str2, str3, str4)
    if m == 0:
        return

    newpdb = replace_lines(pdbfile, outfile, str1, str2)
    newpdb1 = replace_lines(newpdb, outfile, str3, str4)

    delete_file(newpdb)

    if os.path.exists(newpdb1):
        shutil.move(newpdb1, outfile)


##########################################################
def check_string_inpdb(pdbfile, *items_inp):
    """check if items are in the pdbfile"""

    items = list(items_inp)
    if not os.path.exists(pdbfile):
        return 0
    fp = open(pdbfile, "r")
    npdb = 0
    for x in fp:
        if "ATOM" in x[:4] or "HETA" in x[:4] or "ANISOU" in x[:6]:
            continue
        for y in items[:]:
            if y in x[: len(y)]:
                items.remove(y)
                npdb = npdb + 1
                break

    fp.close()
    if len(items_inp) != npdb:
        #        print(items, ' not in the pdb file %s' %pdbfile)
        return 0

    return 1


##########################################################


def delete_file(*files):
    for x in files:
        os.system("rm -f " + x)


##########################################################
def reorder_pdb(pdbfile):
    """sort the ligand and waters by chain and residue number"""

    poly, ligand, water = separate_pdb(pdbfile)
    k1, k2 = 21, 22  # sort chain in order
    lignew1 = sort_column_pdb(ligand, k1, k2, 1)
    watnew1 = sort_column_pdb(water, k1, k2, 1)

    k1, k2 = 22, 26  # sort residue number from small to large
    lignew2 = sort_column_pdb(lignew1, k1, k2, 2)
    watnew2 = sort_column_pdb(watnew1, k1, k2, 2)

    pdbnew = pdbfile + "_sort"

    lig1 = pdbfile + "_lig1"
    lig2 = pdbfile + "_lig2"
    fw1, fw2 = open(lig1, "w"), open(lig2, "w")
    for x in open(lignew2, "r").readlines():
        if "ATOM" in x[:4] or "HETATM" in x[:6] or "ANISOU" in x[:6]:
            fw1.write(x)
        else:
            fw2.write(x)
    fw1.close(), fw2.close()

    os.system("cat %s %s %s  %s >%s" % (poly, lig1, watnew2, lig2, pdbnew))

    delete_file(poly, ligand, water, lignew1, watnew1, lignew2, watnew2, lig1, lig2)
    return pdbnew


##########################################################
def sort_column_pdb(pdb, k1, k2, id):
    """Sort the columns of PDB file in order.
    id=1, sort chain for all ligand/water;
    id=2, sort residue number in order.
    """

    fp = open(pdb, "r").readlines()
    pdb_new = pdb + "_sort"
    fw = open(pdb_new, "w")

    n, d = 0, []
    for ln in fp:
        n = n + 1
        if "ATOM" in ln[:4] or "HETATM" in ln[:6] or "ANISOU" in ln[:6]:
            ch = ln[20:22].strip()
            d.append(ln)
            if id == 1:
                cond = n == len(fp)
            if id == 2:
                cond = n == len(fp) or ch != fp[n][20:22].strip()
            if cond:
                d1 = sort_column(d, k1, k2)
                for x in d1:
                    fw.write("".join(x))
                d = []

        else:
            if len(d) > 0:
                d1 = sort_column(d, k1, k2)
                for x in d1:
                    fw.write("".join(x))
                d = []

            fw.write(ln)

    fw.close()
    return pdb_new


##########################################################
def sort_column(fp, k1, k2):
    """Sort the columns of a list (fp); k1,k2 column position)"""

    d1 = [[x[:k1], x[k1:k2], x[k2:]] for x in fp]  # separate 3 column
    d1.sort(key=lambda y: y[1])  # sort 2th column

    return d1


##########################################################
def sort_column_list(fp, k1, k2):
    """Sort the columns of PDB file in order.
    (input a list; k1,k2 column position)
    """

    fpnew, tmp = [], []
    for n, ln in enumerate(fp):
        if atom_record(ln):
            tmp.append(ln)
            if n == len(fp) - 1:  # last one
                d1 = sort_column(tmp, k1, k2)
                fpnew.extend(d1)
                tmp = []

        else:
            if len(tmp) > 0:
                d1 = sort_column(tmp, k1, k2)
                fpnew.extend(d1)
                tmp = []

    return fpnew


##########################################################
def get_tlsinp(pdbfile):
    """Get TLS input for TLSANL:
    This could be done by tlsextract, but data columns can be connected,
    then TLSANL will fail.
    """

    tlsinp = pdbfile + "_tls.inp"
    fr = open(pdbfile, "r")
    fo = open(tlsinp, "w")
    fo.write("REFMAC \n\n")

    n = 0
    for ln in fr:
        if ln[0:4] == "ATOM" or ln[0:6] == "HETATM":
            break
        if "REMARK   3" not in ln[:12]:
            continue

        if "TLS GROUP :" in ln:
            ntls = get_value_after_id(ln, ":")
            fo.write("TLS GROUP : " + ntls + "\n")

        elif "RESIDUE RANGE :" in ln:  # for refmac
            t1 = parse_tls_range_refmac(ln, ntls)
            if len(t1) < 2:
                continue
            fo.write("RANGE  '%s%4d.' '%s%4d.' ALL\n" % (t1[0], t1[1], t1[0], t1[2]))

        elif "  SELECTION: " in ln:  # for phenix
            if "NULL" in ln:
                print(pdbfile + "; Error! wrong residue range(%s)\n" % ln[25:].strip())
                break
            n = n + 1
            if n == 1:
                chain, chain_range, chain1, chain_range1, chain2, chain_range2 = chain_res_range(pdbfile)
            range = convert_range_phenix_to_refmac(ln, chain_range)
            for k, v in range.items():
                for x in v:
                    arg = "RANGE  '%s%4s.' '%s%4s.' ALL\n" % (k, x[0], k, x[1])
                    fo.write(arg)

        elif "REMARK   3    ORIGIN FOR THE GROUP (" in ln and "):" in ln:
            tmp = ln.split()

            if "NULL" in ln or "0.0000   0.0000   0.0000" in ln or len(tmp) > 10:
                print(pdbfile + "; Error! wrong TLS origin (%s).\n" % ln[40:].strip())
                break
            x, y, z = origin(ln)
            fo.write("ORIGIN  %8s %8s %8s \n" % (x, y, z))

        elif "T11:" in ln and "T22:" in ln:
            T11 = get_one_string_after_id(ln, "T11:")
            T22 = get_one_string_after_id(ln, "T22:")

        elif "T33:" in ln and "T12:" in ln:
            T33 = get_one_string_after_id(ln, "T33:")
            T12 = get_one_string_after_id(ln, "T12:")

        elif "T13:" in ln and "T23:" in ln:
            T13 = get_one_string_after_id(ln, "T13:")
            T23 = get_one_string_after_id(ln, "T23:")

            if "NULL" in T11 + T22 + T33 + T12 + T13 + T23:
                t = 'Error! "NULL" values in matrix Tij for tls group=%s\n' % ntls
                perror(t)
                break

            fo.write("T  %s %s %s %s %s %s\n" % (T11, T22, T33, T12, T13, T23))

        elif "L11:" in ln and "L22:" in ln:
            L11 = get_one_string_after_id(ln, "L11:")
            L22 = get_one_string_after_id(ln, "L22:")

        elif "L33:" in ln and "L12:" in ln:
            L33 = get_one_string_after_id(ln, "L33:")
            L12 = get_one_string_after_id(ln, "L12:")

        elif "L13:" in ln and "L23:" in ln:
            L13 = get_one_string_after_id(ln, "L13:")
            L23 = get_one_string_after_id(ln, "L23:")

            if "NULL" in L11 + L22 + L33 + L12 + L13 + L23:
                t = 'Error! "NULL" values in matrix Lij for tls group=%s\n' % ntls
                perror(t)

                break

            fo.write("L  %s %s %s %s %s %s\n" % (L11, L22, L33, L12, L13, L23))

        elif "S11:" in ln and "S12:" in ln and "S13:" in ln:
            S11 = get_one_string_after_id(ln, "S11:")
            S12 = get_one_string_after_id(ln, "S12:")
            S13 = get_one_string_after_id(ln, "S13:")

        elif "S21:" in ln and "S22:" in ln and "S23:" in ln:
            S21 = get_one_string_after_id(ln, "S21:")
            S22 = get_one_string_after_id(ln, "S22:")
            S23 = get_one_string_after_id(ln, "S23:")

        elif "S31:" in ln and "S32:" in ln and "S33:" in ln:
            S31 = get_one_string_after_id(ln, "S31:")
            S32 = get_one_string_after_id(ln, "S32:")
            S33 = get_one_string_after_id(ln, "S33:")

            if "NULL" in S11 + S22 + S33 + S12 + S13 + S23:
                t = 'Error! "NULL" values in matrix Sij for tls group=%s\n' % ntls
                perror(t)
                break

            S22_S11 = float(S22) - float(S11)
            S11_S33 = float(S11) - float(S33)

            fo.write("S  %8.4f %8.4f %s %s %s %s %s %s\n\n" % (S22_S11, S11_S33, S12, S13, S23, S21, S31, S32))

    fr.close(), fo.close()
    return tlsinp


##########################################################
def origin(ln):
    """only for tls origin"""

    t1 = ln.split(":", 1)[1].split()
    if len(t1) != 3:
        perror("Warning: origin connected in %s" % ln)
        return ln[39:48], ln[48:57], ln[57:66]
    return t1[0], t1[1], t1[2]


##########################################################
def convert_range_phenix_to_refmac(ln, chain_range):
    """Convert TLS residue range from phenix to refmac"""

    line = ln[24:].replace("(", "").replace(")", "").replace("'", "")
    tmp = line.upper().split()
    n = len(tmp)

    val = {}
    v = []

    if n == 2:
        if tmp[0] == "CHAIN" and len(tmp[1]) == 1:
            v.append(chain_range[tmp[1]])
            val[tmp[1]] = v

    elif n == 5:
        if "CHAIN" in tmp[0] and len(tmp[1]) == 1 and "AND" in tmp[2] and "RES" in tmp[3]:
            res = resi_range(tmp[4])
            v.append(res)
            val[tmp[1]] = v

        elif "RES" in tmp[0] and len(tmp[4]) == 1 and "AND" in tmp[2] and "CHAIN" in tmp[3]:
            val[tmp[4]] = v
            if ":" in tmp[1]:
                t1 = tmp[1].split(":")
                if not t1[0]:
                    t1[0] = "-10"
                v.append([t1[0], t1[1]])

    elif n == 8:
        if "CHAIN" in tmp[0] and "AND" in tmp[2] and "RES" in tmp[3] and "OR" in tmp[5] and "RES" in tmp[6]:
            chain = tmp[1]
            val[chain] = v
            res = resi_range(tmp[4])
            v.append(res)
            res = resi_range(tmp[7])
            v.append(res)

    #    print(tmp , n, val)

    return val


##########################################################
def resi_range(tmp):
    v = []
    if ":" in tmp:
        t1 = tmp.split(":")
        if not t1[0]:
            t1[0] = "-10"
        v = [t1[0], t1[1]]

    elif ":" not in tmp and "-" in tmp:
        t1 = tmp.split("-")
        v = [t1[0], t1[1]]

    else:
        v = [tmp, tmp]

    return v


##########################################################
def do_tlsanl(pdbfile, tls_new_inp, bresid, isoout):
    """execute csh script to get a new pdb file with modified B=Bres+Btls"""

    csh_script = """#!/bin/csh -f

############################################################### 
#   tlsanl: to correct B factors (B_total = B_residual + B_tls),
#     when the following conditions are satisfied
#     a). The PDB coordates were refined by REFMAC
#     b). When TLS was involved in refinement
#     c). When the residue ranges are correct (so that tlsanl works).
##################### Usage ################################## 
#  tlsanl_script.csh  pdbfile      (use tls.inp  by tlsextract)
#    or 
#  tlsanl_script.csh  pdbfile  tls.inp  (use the provided tls.inp)
##################### Do TLSAN ###############################

if( $#argv <1 ) then
    echo "Usage: tlsanl_script.csh  pdbfile"
    echo " OR "
    echo "Usage: tlsanl_script.csh  pdbfile  tls.inp"
    exit()
endif

set pdbfile=$1
set refmac=`grep "^REMARK   3   PROGRAM " -m 1 $pdbfile | awk '{print $5}' `  
if( $refmac !~ "REFMAC" ) then
echo "Entry ($1) not refined by REFMAC!"
#    exit()
endif

set aniso=`grep "^ANISOU" -m 2 $pdbfile | wc -l  | awk '{print $1}' ` 
if ( $aniso >1 ) then
echo "Note: ANISOU records exist in the file=$pdbfile."
#    cp $pdbfile  ${pdbfile}_tls
#    exit()
endif

if ($#argv <2) then    # no tls.inp, extract tls from PDB
    
    set tls_inp = "${pdbfile}_tls.inp"
    tlsextract XYZIN $pdbfile TLSOUT $tls_inp >/dev/null

    set size=`wc -l $tls_inp  | awk '{print $1}'`
    if( ! -e $tls_inp ||  $size <= 6 ) then
        /bin/rm -f $tls_inp
        echo "Error! TLS parameters were not extracted!"
        exit()
    endif
else  #tls.inp provided!
    set tls_inp = $2
endif


# do tlsanl to get a new pdb, memo below:
# BRESID [true | t | false | f]: if t, pdb contains B_partial; if f, B_full
# ISOOUT [FULL | RESI | TLSC]: RESI, B_partial; TLSC, TLS contribution

#echo "Doing TLSANL to get full B factors($1)"
set tlspdb="${pdbfile}_tls"
set tlslog="${pdbfile}_tls.log"
tlsanl XYZIN  $pdbfile  TLSIN  $tls_inp  XYZOUT  $tlspdb <<end >$tlslog
BRESID %s
#NUMERIC
ISOOUT %s 
end

""" % (
        bresid,
        isoout,
    )
    pdb_tls = pdbfile + "_tls"
    log = pdb_tls + ".log"
    if os.path.exists(pdb_tls):
        os.remove(pdb_tls)
    if os.path.exists(log):
        os.remove(log)

    script = "tlsanl_script.csh"
    fo = open(script, "w")
    fo.write(csh_script)
    fo.close()
    command = "chmod +x %s ; ./%s %s %s" % (script, script, pdbfile, tls_new_inp)
    os.system(command)
    #    sys.exit()

    return (pdb_tls, log)


##########################################################
def check_bfactor(pdbfile):
    """check B factors in the PDB. if B_full<0 add a constant"""

    bval, negb, bigb, wk = [], [], 0, 1
    if os.path.exists(pdbfile) and os.path.getsize(pdbfile) > 1000:
        fr = open(pdbfile, "r")
        for x in fr:
            if (x[0:5] == "ATOM " or x[0:6] == "HETATM") and "*" in x[60:66]:
                bigb = bigb + 1
                continue
            elif (x[0:5] == "ATOM " or x[0:6] == "HETATM") and float(x[60:66]) < 0:
                negb.append(float(x[60:66]))

            if x[0:4] == "ATOM" or x[0:4] == "HETA":
                bval.append(float(x[60:66]))

        fr.close()
    n = len(negb)
    if bigb > 0:
        print("Error (%s): %d atom with B_full>999." % (pdbfile, bigb))
    if n > 0:
        print("Warning (%s): (%d) atom with B_full<0." % (pdbfile, n))

    return bigb, bval


##########################################################
def delete_0_origin(pdbfile):
    """delete the tls group with zero origins.This group is empty"""

    fp = open(pdbfile, "r")
    pdb_new = pdbfile + "_000"
    fw = open(pdb_new, "w")
    i, id, ntls, nn = 0, 0, 0, []
    for x in fp:
        nn.append(x)
        if "REMARK   3   TLS GROUP :" in x:
            ntls = x.split(":")[1].strip()

        elif "ENDMDL" in x[:6]:
            print("This is a multi-model entry")
            # break

        elif "REMARK   3    ORIGIN FOR THE GROUP (" in x and "):" in x:
            val = x.split(":")[1].split()
            if "NULL" in val or not is_number(x[39:48]) or not is_number(x[48:57]):
                print("Error: wrong origin: " + x)
                continue
            if float(x[39:48]) == 0.000 and float(x[48:57]) == 0.000:
                n = len(nn)
                for m in range(n):  # move up (remove tls group)
                    k = n - m - 1
                    # print(i, k, nn[k])
                    if "REMARK   3   TLS GROUP :" in nn[k]:
                        nn.pop(k)
                        break
                    nn.pop(k)

                for x in fp:  # move down
                    id = id + 1
                    if "S31:" in x:
                        break
                    continue

        i = i + 1
    for x in nn:
        fw.write(x)
    if id > 3:
        t = "Warning: tls group (%s) is deleted, since all parameters are zero." % ntls
        perror(t)
    fp.close(), fw.close()
    return pdb_new


##########################################################
def correct_tls(pdbfile):
    """check PDB for possible corrections
    Correct TLS residue range if given as default (for TLSANL)
    """

    pdb_new = pdbfile + "_new"
    fw = open(pdb_new, "w")
    fr = open(pdbfile, "r")

    res_range = chain_res_range(pdbfile)  # dic with chainID and residue range
    print("Polymer= " % res_range[1])
    print("Ligands= " % res_range[3])
    print("Waters = " % res_range[5])

    range1 = {}
    prog, remark_b, tls_b, anisou, ntls = "?" * 5
    for line in fr:
        if "REMARK   3   PROGRAM     :" in line and "REFMAC" in line:
            prog = "REFMAC"

        elif "NUMBER OF TLS GROUPS  :" in line:
            ntls = get_value_after_id(line, ":")
        elif "ATOM RECORD CONTAINS RESIDUAL B FACTORS ONLY" in line:
            tls_b = "RESID"
        elif "ATOM RECORD CONTAINS SUM OF TLS AND RESIDUAL B " in line:
            tls_b = "FULL"
        elif "REMARK   3  U VALUES      : RESIDUAL ONLY" in line:
            remark_b = "RESID"
        elif "REMARK   3  U VALUES      : WITH TLS ADDED" in line:
            remark_b = "FULL"

        elif "ANISOU" in line[0:6]:
            anisou = "T"

        elif "REMARK   3    RESIDUE RANGE : " in line:
            tmp = line.split(":")[1].strip().split()
            if " NULL " not in line and len(tmp) == 4 and tmp[0] == tmp[2]:
                if is_digit(tmp[1]) <= 0 or is_digit(tmp[3]) <= 0:
                    print("Error! (%s; Residue-range: %s)" % (pdbfile, line[29:].strip()))
                    continue

                if tmp[0] not in range1.keys() and tmp[0] != " ":
                    range1[tmp[0]] = []
                range1[tmp[0]].append([int(eval(tmp[1])), int(eval(tmp[3]))])
                line = correct_tls_overlap(pdbfile, line, range1)
                if len(line) > 0:
                    line = correct_tls_overlap(pdbfile, line, range1)
                if len(line) < 5:
                    print("Warning! (%s)skiping this line. \n" % (pdbfile))
                    # continue
            line1 = correct_residue_range(ntls, pdbfile, line, res_range)

            if len(line1) > 10:
                fw.write(line1)
            continue

        elif line[:6] == "ENDMDL":
            print("Warning! The PDB file %s has multi-model. (1st model used)" % pdbfile)
            break

        fw.write(line)
    #    print(prog,ntls, tls_b,  remark_b, anisou)
    fr.close()
    fw.close()
    #    if remark_b == 'FULL' or tls_b == 'FULL' or anisou == 'T' :
    #        print(pdbfile + "; Warning! PDB file with B_full (or with ANISOU record).\n")
    #        sys.exit()
    return pdb_new


##########################################################
def correct_tls_overlap(pdbfile, line, range):
    tmp = line.split(":")[1].split()
    ch, n1, n2 = tmp[0], int(eval(tmp[1])), int(eval(tmp[3]))  # from current line
    n = 0
    ln = line
    # print(ch, range[ch])
    for x in reversed(range[ch]):
        n = n + 1
        if n == 1:
            continue
        k1, k2 = x[0], x[1]
        # print(n1, n2, '|', k1, k2)
        if n2 < n1 or k2 < k1:
            print("Error! Residue range error [%d %d] " % (n1, n2), x)
            return ln

        elif (n1 >= k1 and n2 <= k2) or (n1 < k1 and n2 > k2):
            print("Error! Total overlaps for [%d %d] " % (n1, n2), x)
            return ""

        elif n1 <= k2 and n1 >= k1 and n2 > k2:  # overlap current range with previous one
            t = "%4s%6d%9s%6d" % (tmp[0], k2 + 1, tmp[0], n2)
            t1 = "Warning: residue range corrected.(%s)->(%s)" % (line[30:60].strip(), t)
            perror(t1)
            ln = line[:29] + t + "   \n"
            break

        elif n2 >= k1 and n2 <= k2 and n1 < k1:  # overlap current range with previous one
            t = "%4s%6d%9s%6d" % (tmp[0], n1, tmp[0], k1 - 1)
            t1 = "Warning: residue range corrected.(%s)->(%s)" % (line[30:60].strip(), t)
            perror(t1)
            ln = line[:29] + t + "   \n"
            break

    return ln


##########################################################
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        # print('Warning! %s is not a number' %s)
        return False


#########################################################
def is_digit(s):
    n = 0
    if len(s) == 0:
        return 0
    if s[0:1] == "-":
        if s[1:].isdigit():
            n = 1
    else:
        if s.isdigit():
            n = 1

    return n


##########################################################
def correct_residue_range(ntls, pdbfile, line, res_range):
    """Correct various errors in TLS residue range"""

    line1 = line
    chain, chain_range, chain1, chain_range1, chain2, chain_range2 = res_range

    res = line[29:].split()
    nr = len(res)
    ch1, n1, ch2, n2 = (" ", 0, " ", 0)

    if "NULL" in line or nr < 2:
        print(pdbfile + "; Error! wrong residue range(%s)." % line[29:].strip())

    elif nr == 4:  # most common cases
        ch1 = ch2 = res[0]
        if res[0] != res[2]:
            print(pdbfile + "; Error! wrong residue range(%s).\n" % line[29:].strip())

        else:
            line1 = correct_range(pdbfile, res[1], res[3], line, ch1, res_range)

    elif nr == 3:
        ch1 = ch2 = res[0]

        line1 = correct_range(pdbfile, res[1], res[2], line, ch1, res_range)
        t = "Warning: residue range correction(%s)->(%s).\n" % (line[29:55].strip(), line1[29:55].strip())
        perror(t)

    elif nr == 2 and is_digit(ntls) > 0 and int(ntls) == 1 and is_digit(res[0]) > 0 or is_digit(res[1]) > 0:
        if "-" in res[0] and "-" in res[1]:
            line1 = ""
            if len(chain) > 0:
                line1 = correct_default(line1, line, chain_range, pdbfile)
            if len(chain1) > 0:
                line1 = correct_default(line1, line, chain_range1, pdbfile)
            if len(chain2) > 0:
                line1 = correct_default(line1, line, chain_range2, pdbfile)

        elif int(res[1]) > int(res[0]):
            ch1 = list(chain.keys())[0]
            line1 = correct_range(pdbfile, res[0], res[1], line, ch1, res_range)
            t = "Warning: residue range correction(%s)->(%s).\n" % (line[29:55].strip(), line1[29:55].strip())
            perror(t)

        else:
            ch1 = list(chain.keys())[0]
            line1 = correct_range(pdbfile, res[0], res[1], line, ch1, res_range)

    return line1


##########################################################
def correct_default(line1, line, chain_range, pdbfile):
    for x, y in chain_range.items():
        rrange = "%4s%6d%9s%6d" % (x, y[0], x, y[1])
        tmp = line[:29] + rrange + "   \n"
        line1 = tmp + line1
        t = "Warning: residue range correction (%s)->(%s).\n" % (line[34:55].strip(), tmp[30:55].strip())
        perror(t)
    return line1


##########################################################
def correct_range(pdbfile, res1, res2, ln, ch, res_range):
    """correct various residue range errors to pass the TLSANL"""

    chain, chain_range, chain1, chain_range1, chain2, chain_range2 = res_range

    line = ""

    if is_digit(res1) == 0 or is_digit(res2) == 0:
        print(pdbfile + "; Error! wrong residue range(%s)." % ln[29:].strip())
        return ln

    elif ch not in chain_range.keys() and ch not in chain_range1.keys() and ch not in chain_range2.keys():
        print(pdbfile + "; Error! Chain (%s) not in coordinate." % ch)
        return ln

    n1, n2 = int(res1), int(res2)  # tls residue number
    range1 = []
    if ch in chain_range.keys():
        range1.append([chain_range[ch][0], chain_range[ch][1]])
    if ch in chain_range1.keys():
        range1.append([chain_range1[ch][0], chain_range1[ch][1]])
    if ch in chain_range2.keys():
        range1.append([chain_range2[ch][0], chain_range2[ch][1]])

    np = 0
    for x in range1:  #
        if (n1 < x[0] and n2 < x[0]) or (n1 > x[1] and n2 > x[1]):
            np = np + 1

    if ch in chain_range.keys() and (n1 in chain[ch] or n2 in chain[ch] or (n1 < chain_range[ch][0] and n2 > chain_range[ch][1]) and n2 <= 999):  # polymer
        line = check_res_range(pdbfile, ln, chain, ch, n1, n2, 0)

    elif ch in chain_range1.keys() and (n1 in chain1[ch] or n2 in chain1[ch]):  # ligands
        line = check_res_range(pdbfile, ln, chain1, ch, n1, n2, 1)

    elif ch in chain_range2.keys() and (n1 in chain2[ch] or n2 in chain2[ch]):  # waters
        line = check_res_range(pdbfile, ln, chain2, ch, n1, n2, 2)

    elif ch in chain_range.keys() and ((n1 < 0 and n2 < 0) or (n1 < 1 and n2 > 999)):  # default (all in ch)
        poly1, poly2 = chain_range[ch][0], chain_range[ch][1]
        tmp1, tmp2, tmp3 = "", "", ""
        if ch in chain_range.keys():  # poly
            n1, n2, id = chain_range[ch][0], chain_range[ch][1], 1
            tmp1 = make_new_line(pdbfile, ln, ch, n1, n2, 1)

        if ch in chain_range1.keys():  # lig
            n1, n2, id = chain_range1[ch][0], chain_range1[ch][1], 1

            if n1 < poly1 and n2 > poly2:  # lig residue range span polymer
                tmp2 = make_two_newline(pdbfile, chain1[ch], ln, ch, poly1, poly2, 1)
            else:
                tmp2 = make_new_line(pdbfile, ln, ch, n1, n2, 1)

        if ch in chain_range2.keys():  # water
            n1, n2, id = chain_range2[ch][0], chain_range2[ch][1], 1

            if n1 < poly1 and n2 > poly2:  # lig residue range span polymer
                tmp3 = make_two_newline(pdbfile, chain2[ch], ln, ch, poly1, poly2, 1)
            else:
                tmp3 = make_new_line(pdbfile, ln, ch, n1, n2, 1)

        line = tmp1 + tmp2 + tmp3

    elif len(range1) == np and np > 0:  #
        line = " "
        t = "Warning: Residue range [%s %d-%d] not exist in coordinate. This range removed." % (ch, n1, n2)
        perror(t)

    if len(line) == 0:
        line = ln

    return line


##########################################################
def check_res_range(pdbfile, ln, chain, ch, n1, n2, id):
    line = ln
    if n1 <= n2:
        id, nd = 0, n2 - n1
        if n1 not in chain[ch]:
            n1 = new_number(chain[ch], n1, n2, 1)
        if n2 not in chain[ch]:
            n2 = new_number(chain[ch], n1, n2, -1)

        if nd != n2 - n1:
            id = 1
        line = make_new_line(pdbfile, ln, ch, n1, n2, id)

    elif n1 > n2:
        line = make_new_line(pdbfile, ln, ch, n2, n1, 1)

    return line


##########################################################
def make_two_newline(pdbfile, resn_in, ln, ch, n1, n2, id):
    resn = sorted(resn_in)
    # print(resn, ch, n1, n2, id)
    k1, k2, nn = 0, 0, len(resn) - 1
    for n in range(nn):
        if resn[n + 1] > n1:
            k2 = n
            break

    m1, m2, m3, m4 = resn[0], resn[k2], resn[k2 + 1], resn[nn]
    # print(m1, m2, m3, m4, ';', resn, ch, n1, n2, id)

    tmp1 = make_new_line(pdbfile, ln, ch, m1, m2, id)
    tmp2 = make_new_line(pdbfile, ln, ch, m3, m4, id)

    return tmp1 + tmp2


##########################################################
def make_new_line(pdbfile, ln, ch, n1, n2, id):
    """return a new line with correct residue ranges.
    id==1, corrected. id==0, original
    """

    tmp = "%4s%6d%9s%6d" % (ch, n1, ch, n2)
    if id > 0:
        t = "Warning: residue range correction(%s)->(%s).\n" % (ln[30:60].strip(), tmp.strip())
        perror(t)
    line1 = ln[:29] + tmp + "   \n"

    return line1


##########################################################
def new_number(chain, n1, n2, id):
    if id == 1:  # increase the number
        n = n1
        for i in range(n2 - n1 + 1):
            if n in chain or n >= n2:
                break
            n = n1 + i
        if (n - n1) > 20:
            t = "Warning: Residue range changed too much (%d).Watch out!" % (n - n1)
            perror(t)

    else:  # decrease the number
        n = n2
        for i in range(n2 - n1 + 1):
            if n in chain or n < n1:
                break
            n = n2 - i

        if (n2 - n) > 20:
            t = "Warning: Residue range changed too much (%d).Watch out!" % (n2 - n)
            perror(t)

    return n


##########################################################
def get_value_after_id(line, id):
    """get value after a given id"""

    if id not in line:
        return "?"
    li = line.split(id)
    value = li[1].strip()
    if len(value) <= 0 or value.upper() == "NULL":
        value = "?"
    return value


##########################################################
def get_one_string_after_id(line, id):
    """get one value after a given id"""

    if id not in line:
        return "?"
    value = line.split(id)[1].split()[0].strip()
    if not value:
        value = "?"
    # print(line, id, value)
    return value


##########################################################
def chain_res_range(pdbfile):
    """use dic to contain chain-ID and residue range
    separate poly-peptide and hetatoms:
    """

    pdb1, pdb2, pdb3 = separate_pdb(pdbfile)
    chain, chain_range = chain_res(pdb1, 0)  # poly-peptide
    chain1, chain_range1 = chain_res(pdb2, 0)  # ligands
    chain2, chain_range2 = chain_res(pdb3, 1)  # waters
    delete_file(pdb1, pdb2, pdb3)

    return chain, chain_range, chain1, chain_range1, chain2, chain_range2


##########################################################
def chain_res_range_list(pdb):
    """use dic to contain chain-ID and residue range
    separate poly-peptide and hetatoms. Input a list.
    """

    pdb1, pdb2, pdb3, pdb4 = separate_pdb_list(pdb)
    chain, chain_range = chain_res_list(pdb1, 0)  # poly-peptide
    chain1, chain_range1 = chain_res_list(pdb2, 0)  # ligands
    chain2, chain_range2 = chain_res_list(pdb3, 1)  # waters

    return chain, chain_range, chain1, chain_range1, chain2, chain_range2


##########################################################
def separate_pdb(pdbfile):
    """separate PDB file by polymer pdb1 and ligands pdb2, water pdb3"""
    pdb1 = pdbfile + "tmp1"
    pdb2 = pdbfile + "tmp2"
    pdb3 = pdbfile + "tmp3"

    fw1, fw2, fw3 = open(pdb1, "w"), open(pdb2, "w"), open(pdb3, "w")

    fr, fr1 = [], open(pdbfile, "r").readlines()
    for x in fr1:
        if atom_record(x) and ("HOH" in x[17:20] or "DOD" in x[17:20]):
            fw3.write(x)
        else:
            fr.append(x)

    length = len(fr)
    i, n = 0, 0
    for i in range(length):
        k = length - i - 1
        if atom_record(fr[k]) or "TER" in fr[k][:3]:  # go up until aa or na
            if residue(fr[k][17:20]) and "ATOM" in fr[k][:5]:
                n = k
                break

    if i > length - 3:
        n = length
    for i in range(length):
        if i <= n:
            fw1.write(fr[i])
        else:
            fw2.write(fr[i])

    fw1.close(), fw2.close(), fw3.close()
    return pdb1, pdb2, pdb3


##########################################################
def separate_pdb_list(pdb):
    """separate PDB file by polymer pdb1 and ligands pdb2, water pdb3,tail pdb4
    poly (pdb1) has all the remarks.
    """

    pdb1, pdb2, pdb3, pdb4, fr, tmp = [], [], [], [], [], []  # poly,lig,wat,tail

    for x in pdb:
        if atom_record(x) and ("HOH" in x[17:20] or "DOD" in x[17:20]):
            pdb3.append(x)  # wat
        else:
            fr.append(x)

    length = len(fr)
    i, n = 0, 0
    for i in range(length):
        k = length - i - 1
        if atom_record(fr[k]) or "TER" in fr[k][:3]:
            if residue(fr[k][17:20]) and "ATOM" in fr[k][:5]:  # go up until aa or na residue
                n = k
                break

    for i in range(length):
        if n == 0:
            tmp.append(fr[i])
        else:
            if i <= n:
                pdb1.append(fr[i])  # poly
            else:
                tmp.append(fr[i])

    for i, x in enumerate(tmp):
        if atom_record(x) or "TER" in x[:3]:
            pdb2.append(x)  # lig
        else:
            pdb4.append(x)

    return pdb1, pdb2, pdb3, pdb4


##########################################################
def separate_pdb_old(pdbfile):
    """separate PDB file by polymer pdb1 and ligands pdb2, water pdb3"""

    pdb1 = pdbfile + "tmp1"
    pdb2 = pdbfile + "tmp2"
    pdb3 = pdbfile + "tmp3"

    fw1, fw2, fw3 = open(pdb1, "w"), open(pdb2, "w"), open(pdb3, "w")

    n = 0
    fr, fr1 = [], open(pdbfile, "r").readlines()

    for x in fr1:
        if ("ATOM" in x[:4] or "HETA" in x[:4] or "ANISOU" in x[:6]) and ("HOH" in x[17:20] or "DOD" in x[17:20]):
            fw3.write(x)
        else:
            fr.append(x)

    length = len(fr)
    i = 0
    for i in range(length):
        k = length - i - 1
        if ("CONECT" in fr[k][:6] or "MASTER" in fr[k][:6] or "END" in fr[k][:3] or "HETATM" in fr[k][:6] or "TER" in fr[k][:3]) and (
            "ATOM" in fr[k - 1][:4] or "TER" in fr[k - 1][:3]
        ):
            n = k
            break
    if i > length - 3:
        n = length
    for i in range(length):
        if i < n:
            fw1.write(fr[i])
        else:
            fw2.write(fr[i])
            """
            if 'HOH' in fr[i]:
                fw3.write(fr[i])
            else:
                fw2.write(fr[i])
            """

    fw1.close(), fw2.close(), fw3.close()
    return pdb1, pdb2, pdb3


##########################################################
def str_after_id(line, id):
    """get string after a given id"""

    if id not in line:
        print("Warning! %s not in string (%s)" % (id, line))
        return line

    n = line.index(id)
    value = line[n + len(id) :].strip()
    if len(value) <= 0 or value.upper() == "NULL":
        value = "?"
    return value


##########################################################
def chain_res(pdbfile, id):
    """use dic to contain chain-ID and residue range
    id=0: not include waters; id=1: include waters.
    """

    fr = open(pdbfile, "r")

    chain = {}
    chain_range = {}
    n_old = -99999
    ch_old = "?"
    for ln in fr:
        if "ATOM" in ln[:4] or "HETATM" in ln[:6]:
            if id == 0 and "HOH" in ln[16:20]:
                continue
            ch = ln[20:22].strip()

            if ch not in chain.keys() and ch != " ":
                chain[ch] = []
            nres = ln[22:26]
            n = int(nres)
            if ch in chain.keys():
                if n == n_old and ch == ch_old:
                    continue
                chain[ch].append(n)
            n_old = n
            ch_old = ch

    for key in chain:
        if len(chain[key]) > 0:
            chain_range[key] = [min(chain[key]), max(chain[key])]

    fr.close()
    #    print  chain_range
    return chain, chain_range


##########################################################
def chain_res_list(fr, id):
    """use dic to contain chain-ID and residue range, pdb is a list.
    id=0: not include waters; id=1: include waters.
    """

    chain, chain_range = {}, {}
    ch_old, n_old = "?", -99999

    for ln in fr:
        if "ATOM" in ln[:4] or "HETATM" in ln[:6]:
            if id == 0 and "HOH" in ln[16:20]:
                continue
            ch = ln[20:22].strip()
            if ch not in chain.keys() and ch != " ":
                chain[ch] = []
            nres = ln[22:26]
            n = int(nres)
            if ch in chain.keys():
                if n == n_old and ch == ch_old:
                    continue
                chain[ch].append(n)
            n_old = n
            ch_old = ch

    for x in chain.keys():
        if len(chain[x]) > 0:
            chain_range[x] = [min(chain[x]), max(chain[x])]

    return chain, chain_range


##########################################################


def chain_resint_atom(pdbfile):
    """use dic to contain chain-ID and residue range, and atom
    like dic={chain:{resn:[atoms, ...]} , ...}
    """

    fr = open(pdbfile, "r")

    n_old, ch_old = -99999, "?"
    chain = {}
    for ln in fr:
        if ("ATOM" not in ln[:4] and "HETATM" not in ln[:6]) or len(ln) < 54:
            continue
        ch = ln[20:22].strip()
        if ch not in chain.keys() and ch != " ":  # 1st line of chain
            chain[ch] = [ln]
        elif ch in chain.keys():
            chain[ch].append(ln)
        ch_old = ch

    record = {}
    for ch in chain.keys():
        record[ch] = {}
        n_old = -99999
        res = {}
        for ln in chain[ch]:  # add residues to the chain
            n = int(ln[22:26])
            if n != n_old:
                res[n] = [ln]  # 1st residue

            if n == n_old:
                res[n].append(ln)  # add atoms to chain[n]
            n_old = n
        record[ch] = res

    fr.close()
    return record


##########################################################
def replace_lines(pdbfile, rcsbpdb, str1, str2):
    """replace the lines between str1 & str2
    return a new pdbfile with modified lines
    """

    newpdb = pdbfile + "_cut_add"

    tls1, pdb1 = cut_between_pattern(pdbfile, str1, str2)
    tls2, pdb2 = cut_between_pattern(rcsbpdb, str1, str2)
    if len(tls1) == 0 or len(tls2) == 0:
        print("Error! cut & paste between two strings failed.")
        return " "

    fw = open(newpdb, "w")
    fr = open(pdb1, "r")
    n = 0
    m = len(str1)
    for x in fr:
        fw.write(x)
        if str1 in x[: m + 1]:  # add tls2 to pdb1
            for y in tls2:
                fw.write(y)
                n = n + 1

    fr.close(), fw.close()
    delete_file(pdb1, pdb2)
    if n == 0:
        print("Error! Nothing added to the pdb=" + pdbfile)
    return newpdb


##########################################################
def cut_between_pattern(pdbfile, str1, str2):
    """Cut the lines between str1 and str2
    return a list for the cutted and the rest lines for a file
    """

    newpdb = pdbfile + "cut1"
    fr1 = open(pdbfile, "r")
    fw1 = open(newpdb, "w")

    tls1 = []
    for x in fr1:
        fw1.write(x)
        if str1 in x:
            for y in fr1:
                if str2 in y:
                    fw1.write(y)
                    break
                tls1.append(y)

    fr1.close()
    fw1.close()
    if len(tls1) == 0:
        print('Error! Nothing between "%s" & "%s" in file (%s)' % (str1, str2, pdbfile))
    return tls1, newpdb


##########################################################
def atom_list(pdbfile):
    """only return the ATOM/HETATM coordinates as list"""

    fr1 = open(pdbfile, "r")
    d1 = []
    for x in fr1:
        if (x[:4] == "ATOM" or x[:6] == "HETATM") and len(x.strip()) > 64:
            d1.append(x)

    fr1.close()
    return d1


##########################################################
def compare_two_pdb(pdbfile, pdbdep):
    """check if the two pdbfiles are the same (every 15th lines are checked)"""

    if not (os.path.exists(pdbfile) and os.path.getsize(pdbfile) > 100 and os.path.exists(pdbdep) and os.path.getsize(pdbdep) > 100):
        print("Error: Either (%s) or (%s) do not exist." % (pdbfile, pdbdep))
        return 0

    print("Comparing two PDBs (%s | %s)" % (pdbfile, pdbdep))
    pdbf1, pdbd1 = pdbfile + "tmp", pdbdep + "tmp"
    os.system('egrep "^ATOM" %s > %s ' % (pdbfile, pdbf1))
    os.system('egrep "^ATOM" %s > %s ' % (pdbdep, pdbd1))

    d1 = atom_list(pdbf1)  # only polymer
    d2 = atom_list(pdbd1)

    delete_file(pdbf1, pdbd1)

    n = len(d1)
    check_list = [int(x) for x in range(0, n, 10)]
    m = 0
    for i in check_list:
        for y in d2[:]:
            if d1[i][30:54] == y[30:54]:
                m = m + 1
                d2.remove(y)  # add this to speed up
                break
    val = 0
    n = len(check_list)
    rate = float(m) / float(n)

    if rate > 0.8 and m > 0:
        val = 1
    else:
        print("Error: files (%s | %s) are different (%d|%d)." % (pdbfile, pdbdep, m, n))
    return val


##########################################################
def match_pdb(tlsd, xx, pdb, dep, record, tmp):
    """use the deposited residue range to match the processed residue range.
    it is possible to have one-->many since ligands and waters are assigned
    to the closer chains.
    """

    ch, n1, n2, m = tmp[0], int(tmp[1]), int(tmp[3]), 0

    pdblist = []
    file = tmp[0] + "_" + tmp[1] + "_" + tmp[3] + ".txt"
    fw = open(file, "w")
    for n in range(n1, n2 + 1):
        if n not in record[ch].keys():
            continue
        for x in record[ch][n]:
            t1 = x[30:54]
            for y in pdb[:]:
                t2 = y[30:54]
                if t1 == t2:
                    fw.write(y)
                    pdblist.append(y)
                    pdb.remove(y)
                    m = m + 1

    fw.close()
    if m == 0:
        print("Error! No correspondence between dep/proc ligands.")
        return

    chain, ch_range = chain_res(file, 1)
    n = tlsd.index(xx)

    print("original=%s" % xx.strip())
    for k in ch_range.keys():
        t1 = "REMARK   3    RESIDUE RANGE : "
        s = t1 + "%3s  %4d %8s  %4d\n" % (k, ch_range[k][0], k, ch_range[k][1])
        tlsd.insert(n, s)
        print("new=%s" % s.strip())

    delete_file(file)


##########################################################
def match_ligand_water(pdbfile, pdbdep):
    """Ligands & water are re-order (closer to the polymer chain in the
    PDB file during deposition, However, tls also have to be moved around
    in order to do correct TLS correction
    """

    n = compare_two_pdb(pdbfile, pdbdep)
    if n == 0:
        return

    print("Matching tls groups between %s and %s ..." % (pdbfile, pdbdep))
    pdb1, pdb2, pdb3 = separate_pdb(pdbfile)  # polymer & ligand
    if os.path.getsize(pdb2) < 1:
        print("Warning: no ligands in pdb=" + pdbfile)
        delete_file(pdb1, pdb2)
        return

    pdb = open(pdb2, "r").readlines()  # for ligands
    dep = open(pdbdep, "r").readlines()

    pch, pch_range, pch1, pch1_range = chain_res_range(pdbfile)
    dch, dch_range, dch1, dch1_range = chain_res_range(pdbdep)

    str1 = "REMARK   3  TLS DETAILS"
    str2 = "REMARK   3  BULK SOLVENT MODELLING."

    tlsd, pdbd = cut_between_pattern(pdbdep, str1, str2)
    tlsp, pdbp = cut_between_pattern(pdbfile, str1, str2)

    if len(tlsd) < 16 or len(tlsp) < 16:
        print("Warning: No TLS records either in (%s) or (%s)" % (pdbfile, pdbdep))
        os.system("rm -f %s %s" % (pdbd, pdbp))
        return

    record = chain_resint_atom(pdbdep)

    tlsd_new = []
    for x in tlsd[:]:
        if "REMARK   3    RESIDUE RANGE :" in x:
            # print(x)
            tmp = str_after_id(x, ":").split()
            if len(tmp) < 4:
                continue
            ch = tmp[0]
            if ch not in pch and ch not in pch1:  # match
                match_pdb(tlsd, x, pdb, dep, record, tmp)
                tlsd.remove(x)

    m = len(tlsd)
    k = 0
    for i in range(m):  # get correct component number
        n = m - i - 1
        if " RESIDUE RANGE :" in tlsd[n]:
            k = k + 1
        if "NUMBER OF COMPONENTS GROUP :" in tlsd[n]:
            tlsd[n] = "REMARK   3    NUMBER OF COMPONENTS GROUP :  %2d\n" % k
            k = 0

    #    for s in tlsd:   print(s.strip())

    newpdb = pdbfile + "_newtls"

    fw = open(newpdb, "w")
    fr = open(pdbp, "r")
    delete_file(pdbd, pdbp)

    for x in fr:  # add tls2 to pdb1
        fw.write(x)
        if str1 in x[:24]:
            for y in tlsd:
                fw.write(y)

    fr.close(), fw.close()


##########################################################
def get_tls_group(fp):
    """input a list
    output: {group:[[residue ranges] [] ...]}
    """

    tls, res, tmp = {}, [], []
    for x in fp:
        if "REMARK   3   TLS GROUP :" in x:
            ntls = x.split(":")[1].strip()

        elif "REMARK   3    RESIDUE RANGE :" in x:
            tmp = parse_tls_range_refmac(x, ntls)
            if tmp:
                res.append(tmp)

        elif "REMARK   3    ORIGIN FOR THE GROUP (A):" in x:
            if is_number(ntls) and res:
                tls[int(ntls)] = res
                res, ntls = [], "?"
        elif "BULK SOLVENT" in x or atom_record(x):
            break

    return tls


##########################################################
def check_and_correct_tls(pdbfile):
    """check tls residue ranges; If tls residue ranges is mixed (polymer/lig/wat),
    it will be separated.

    """

    fp = open(pdbfile, "r").readlines()
    out = pdbfile + "__TLSNEW"
    fw = open(out, "w")

    tls = get_tls_group(fp)

    xyzt = []
    for x in fp:
        if atom_record(x):
            xyzt.append(x)

    # print(xyzt)
    for x in fp:
        if "REMARK   3   TLS GROUP :" in x:
            ntls = x.split(":")[1].strip()
            fw.write(x)
        elif "REMARK   3    RESIDUE RANGE :" in x:
            tmp = parse_tls_range_refmac(x, ntls)
            if not tmp:
                continue
            res = x.split(":")[1].split()
            xyz = get_xyz_tls(res, xyzt)  # get xyz in the tls range
            pdb1, pdb2, pdb3, pdb4 = separate_pdb_list(xyz)  # poly/lig/water
            if pdb1:
                write_tmp(pdbfile, fw, x, pdb1, res[0])
            if pdb2:
                write_tmp(pdbfile, fw, x, pdb2, res[0])
            if pdb3:
                write_tmp(pdbfile, fw, x, pdb3, res[0])

        else:
            fw.write(x)

    fw.close()
    return out


##########################################################
def write_tmp(pdbfile, fw, x, pdb1, ch):
    chain1, chain_range1 = chain_res_list(pdb1, 1)
    if ch not in chain_range1.keys():
        return
    n1, n2 = chain_range1[ch][0], chain_range1[ch][1]
    ln = make_new_line(pdbfile, x, ch, n1, n2, 1)
    fw.write(ln)


##########################################################
def match_residue(ntls, chain1_range, xx, pdb, record_dep):
    """use the deposited residue range to match the processed residue range.
    it is possible to have one-->many ranges since ligands and waters are
    assigned to the closer chains.
    ntls: tls group number; residue range for polymer; xx: current line
    pdb: the processed pdb;  record_dep: the deposited pdb
    """

    new_range, matched = [], []
    res = parse_tls_range_refmac(xx, ntls)
    if not res:
        return new_range

    # print('ntls=', ntls, res)
    ch, n1, n2 = res[0], res[1], res[2]
    if ch not in record_dep.keys():
        print("Error: chain (" + ch + ") not in deposited pdb file (TLS group=%s)!" % ntls, res)
        return []

    for n in range(n1, n2 + 1):
        if n not in record_dep[ch].keys():
            continue

        for x in record_dep[ch][n]:
            #            t1=x[32:56]
            t1 = x[29:54]
            for y in pdb[:]:
                #                t2=y[32:56]
                t2 = y[29:54]
                if t1 == t2:
                    matched.append(y)
                    pdb.remove(y)
                    break

    if not matched:
        print("Error! No correspondence between deposited & processed files (TLS group=%s)!" % ntls, res)
        return new_range

    ch1, ch1_range, ch2, ch2_range, ch3, ch3_range = chain_res_range_list(matched)

    new_range = write_resi_range(ch1, ch1_range, chain1_range, res, ntls)
    tmp = write_resi_range(ch2, ch2_range, chain1_range, res, ntls)
    if tmp:
        new_range.extend(tmp)
    tmp = write_resi_range(ch3, ch3_range, chain1_range, res, ntls)
    if tmp:
        new_range.extend(tmp)

    return new_range


##########################################################
def write_resi_range(chain, ch_range, chain1_range, res, ntls):
    ch, n1, n2 = res[0], res[1], res[2]
    new_range = []
    for c in sorted(ch_range.keys()):
        if c not in chain1_range:
            continue
        poly1, poly2 = chain1_range[c][0], chain1_range[c][1]

        k1, k2 = ch_range[c][0], ch_range[c][1]
        t1 = "REMARK   3    RESIDUE RANGE : "
        t2 = "%3s  %4d %8s  %4d\n" % (c, k1, c, k2)
        ln = t1 + t2

        if n1 != k1 or n2 != k2 or c != ch:
            print("RESIDUE RANGE (tls group=%2s): %s -> [%s]" % (ntls, res, t2.strip()))

        if k1 < poly1 and k2 > poly2:  # residue range span polymer
            ln = make_two_newline("!", chain[c], ln, c, poly1, poly2, 1)

        if n2 - n1 + 10 < k2 - k1:
            print("Warning: residue range change too much.", n2 - n1, k2 - k1)
        new_range.append(ln)

    return new_range


##########################################################
def is_cif(file):
    """if a pdb file, return 0; if a cif by pdb_extract, return 1;
    if a cif NOT by pdb_extract, return 2; if file not exist, return -1.
    """

    if not os.path.exists(file):
        t = "Error: file (%s) not exist. please check input file!" % file
        perror(t)
        return -1

    n, n1, n2 = 0, 0, 0
    for x in open(file, "r").readlines():
        t = x.strip()
        if len(t) < 4:
            continue
        if "data_" in t[:5] or "loop_" in t[:5] or ("_" == t.split()[0][0] and "." in t.split()[0]):
            n1 = n1 + 1
        elif "REMARK   3" in x[:10]:
            n2 = n2 + 1
        elif "HEADER " in x[:7] and n1 == 0:
            n = 0
            break
    if n1 > 1 and n2 > 10:  # by pdb_extract
        n = 2
    elif n1 > 1 and n2 < 10:
        n = 1

    return n


##########################################################
def get_newtls_cif(pdbfile):
    """match TLS of pdbfile with the processed one
    pdbfile is either the original (deposited) PDB file or the CIF file
    """

    n = is_cif(pdbfile)

    if n == 0:  # a pdb file
        cif = runmaxit(pdbfile, 1, " ")  # get a cif file
        ncif = runmaxit(cif, 35, " ")  # move (HOH/Lig) around the polymer

        # encif = runmaxit(pdbfile, 7, '') #convert pdb to encif
        # cif = runmaxit(encif, 8, ' -water_number ')  #convert encif to cif
        # ncif = runmaxit(cif, 33, ' -M %s '%pdbfile) #finalize cif(handle W/L)

    else:  # a mmcif file
        cif = pdbfile
        pdbfile = runmaxit(cif, 2, "-exchange_in ")
        ncif = runmaxit(cif, 35, " ")  # move (HOH/Lig) around the polymer

    npdb = runmaxit(ncif, 2, " ")  # new PDB file with HOH/Lig change

    wk, pdbf = match_all_residues(npdb, pdbfile)
    tlscif = tls2cif(pdbf)

    # delete_file(cif)
    if n == 0:
        delete_file(cif)
    # delete_file(npdb, pdbf)

    return tlscif


##########################################################
def match_all_residues(pdbfile, pdbdep):
    """Ligands & water are re-ordered (closer to the polymer chain in the
    PDB file) during deposition, However, residue range of tls groups also have
    to be changed to match the coordinates.

    This function will use the residue ranges (worked in pdbdep) and match all
    the residue ranges in the processed pdb (pdbfile).
    """

    chain1, chain1_range, chain2, chain2_range, chain3, chain3_range = chain_res_range(pdbfile)
    # n=compare_two_pdb(pdbfile, pdbdep)
    # if n==0 : return 0, pdbdep

    print("Matching tls groups between %s and %s" % (pdbfile, pdbdep))

    tlsd, pdbd = remove_pdb_item(pdbdep, 0, "tls")  # get a list & file
    if len(tlsd) < 16 or os.path.getsize(pdbd) < 20:
        print("Warning: No TLS records in (%s) or xyz not extracted.(Stop!)" % (pdbdep))
        delete_file(pdbd)
        sys.exit()
        # return 0, pdbdep

    pdb = atom_list(pdbfile)  # only ATOM/HETA in the list
    record_dep = chain_resint_atom(pdbd)
    ntls, tlsd_new = "?", []
    for x in tlsd:  # through deposited TLS groups
        tlsd_new.append(x)
        if "REMARK   3   TLS GROUP :" in x:
            ntls = x.split(":")[1].strip()

        elif "REMARK   3    RESIDUE RANGE :" in x:
            tmp = match_residue(ntls, chain1_range, x, pdb, record_dep)

            if len(tmp) > 0:
                tlsd_new.pop(-1)  # remove the last one and get new
                for x in tmp:
                    if len(x) > 0:
                        tlsd_new.append(x)

    m = len(tlsd_new)
    k = 0
    for i in range(m):  # get correct component number
        n = m - i - 1
        if " RESIDUE RANGE :" in tlsd_new[n]:
            k = k + 1
        if "NUMBER OF COMPONENTS GROUP :" in tlsd_new[n]:
            tlsd_new[n] = "REMARK   3    NUMBER OF COMPONENTS GROUP :  %2d\n" % k
            k = 0

    # for s in tlsd_new:   print(s.strip())

    pdbp, pdbf = remove_pdb_item(pdbfile, 0, "tls")  # get a list & file
    str1 = "REMARK   3  TLS DETAILS"
    newpdb = insert_lines2file(pdbf, str1, tlsd_new)
    delete_file(pdbd, pdbf)

    return 1, newpdb


##########################################################
def remove_pdb_item(fp_inp, id, item):
    """
    id=0, fp is a pdbfile & return a pdb file;
    id=1, fp is a list of pdb & return two lists
    item, is the items to be removed.
    """

    n1, n2, tls, pdb = 0, 0, [], []

    if id == 0:
        fp = open(fp_inp, "r").readlines()
    else:
        fp = fp_inp

    if item == "tls":
        for n, x in enumerate(fp):
            if "REMARK   3  TLS DETAILS" in x[:26]:
                n1 = n + 1
            elif "S31:" in x[13:20]:
                n2 = n
            elif n2 > 0 and "REMARK   3" not in x[:15]:
                break

    pdb.extend(fp[:n1])
    tls.extend(fp[n1 : n2 + 1])
    pdb.extend(fp[n2 + 1 :])

    # for x in tls: print (x.strip())

    if id == 0:
        ftmp = fp_inp + "tmptls"
        fw = open(ftmp, "w")
        for x in pdb:
            fw.write(x)
        fw.close()
        return tls, ftmp

    return tls, pdb


##########################################################
def insert_lines2file(file, str, lines):
    """lines (a list) are inserted in the file after string (str)"""

    newpdb = file + "_newtls"

    fw = open(newpdb, "w")
    fr = open(file, "r")
    for x in fr:
        fw.write(x)
        if str in x[:24]:
            for y in lines:
                fw.write(y)

    fr.close(), fw.close()

    return newpdb


##########################################################
def search_rcsb_dir(pdbfile, outfile, id):
    """Search all the files in deposited directory (/an../prot/rcsb0.)
    to find the correct TLS group
    """
    rcsbid = get_code_id(pdbfile).lower()
    cif = rcsbid + "-deposited.cif"
    if id == 1:
        files = get_deposited_file(rcsbid, 1)
        for x in files:
            print("file=", x)
        for x in files:
            val = 0
            pdbdep = get_local_file(x)
            print("Trying file=" + pdbdep)
            wk, pdb_match = match_all_residues(pdbfile, pdbdep)
            if wk == 1:
                val = proc_tlsanl(pdb_match, outfile)
            if val == 0 and wk == 1 and cif in pdbdep:
                wk1, pdb_match = get_newtls(pdbdep, cif)
                if wk1 == 1:
                    val = proc_tlsanl(pdb_match, outfile)

            # delete_file(pdb_match)
            if val == 1:
                break
    else:
        p = get_deposited_file(rcsbid, 0)  # path only
        x = p[0] + "/" + pdbfile
        print(x)
        pdbdep = get_local_file(x)
        wk, pdb_match = match_all_residues(pdbfile, pdbdep)
        if wk == 1:
            val = proc_tlsanl(pdb_match, outfile)
        # delete_file(pdb_match)


##########################################################
def get_newtls(pdbdep, cif):
    """Here, the pdbdep is converted from cif to pdb. But for some old
    entries, the TLS group is messed up. The function tries to put TLS
    from cif (remark) to pdbdep

    """

    print("Puting TLS from %s to %s " % (cif, pdbdep))

    str1 = "REMARK   3  TLS DETAILS"
    str2 = "REMARK   3  BULK SOLVENT MODELLING."
    tlsd, pdbd = cut_between_pattern(pdbdep, str1, str2)  # list & file
    tlsp, pdbp = cut_between_pattern(cif, str1, str2)

    if len(tlsd) < 16 or len(tlsp) < 16:
        print("Warning: No TLS records either in (%s) or (%s)" % (cif, pdbdep))
        delete_file(pdbd, pdbp)
        return 0, pdbdep

    newpdb = insert_lines2file(pdbd, str1, tlsp)

    delete_file(pdbd, pdbp)
    print("The new PDB file =%s" % newpdb)
    return 1, newpdb


##########################################################
def get_code_id(pdbfile):
    code = os.popen('grep " ID CODE IS " %s ' % pdbfile).read()
    if len(code) < 2:
        print("Error: No RCSBID was extracted from " + pdbfile)
    else:
        code = code.split()[-1].replace(".", "").lower()
    return code


##########################################################
def get_deposited_file(rcsbid, id):
    """get all the deposited pdb files or cif file by pdb_extract"""

    p1, p2 = "/annotation/prot/" + rcsbid, "/annotation/ndb/" + rcsbid
    files = []
    pdb = p1 + "/" + rcsbid + ".pdb"
    p = p1
    if not os.path.exists(pdb):
        p = p2  # ndb dir
    if id == 0:
        files.append(p)
        return files

    val = os.popen("find %s/*.*" % p).read().split()
    if len(val) == 0:
        print("Error! No files were obtained from " + p)

    for x in val:
        tmp = rcsbid + "-deposited.cif"
        if tmp in x or ("pdb" in x and "rcsb" in x):
            files.append(x)
    if len(files) == 0:
        print("Warning: no extra files in dir=%s" % p)
    return files


##########################################################
def runmaxit(file, option, other):
    """some options to run maxit
    use '-exchange_in' to get tls (from mmcif of pdb_extract)
    """
    #    print('Converting %s by Maxit with option %d ...' %(file, option))
    nfile = file + ".pdb"

    if option == 2:
        nfile = file + ".pdb"
    else:
        nfile = file + ".cif"

    if "pdb_extract" not in other.lower():  # regular
        arg = "maxit-v8.01-O  -i %s -o %d  %s " % (file, option, other)
        os.system(arg)
        return nfile
    else:
        arg = "maxit-v8.01-O  -i %s -o %d " % (file, option)
        os.system(arg)

        tls1, juck1 = remove_pdb_item(file, 0, "tls")  # original cif
        tls2, pdb2 = remove_pdb_item(nfile, 0, "tls")

        str1 = "REMARK   3  TLS DETAILS"
        newpdb = insert_lines2file(pdb2, str1, tls1)
        delete_file(juck1, pdb2)

        return newpdb


##########################################################
def get_local_file(x):
    pdbdep = x.split("/")[-1]  # remove path
    shutil.copy(x, pdbdep)

    if ".gz" in x[-3:]:
        pdbdep = pdbdep[0:-3]
        os.system("zcat %s > %s " % (x, pdbdep))

    elif ".Z" in x[-2:]:
        pdbdep = pdbdep[0:-2]
        os.system("zcat %s > %s " % (x, pdbdep))

    elif ".bz2" in x[-4:]:
        os.system("bunzip2 -f " + pdbdep)
        pdbdep = pdbdep[0:-4]

    if ".cif" in pdbdep:
        pdbdep = runmaxit(pdbdep, 2, "-exchange_in")
    return pdbdep


##########################################################


def change_nres(fw, chain_range1, chain_range2, fp2, shift):
    """change residue number"""

    for ch in sorted(chain_range2.keys()):
        if ch not in chain_range1:  # write all
            for x in fp2[:]:
                if ch == x[20:22].strip():
                    fw.write(x)
                    fp2.pop()
        else:
            nst = chain_range1[ch][1] + shift
            for i, x in enumerate(fp2):
                if ch == x[20:22].strip():
                    if i > 0 and fp2[i][22:26] != fp2[i - 1][22:26]:
                        nst = nst + 1
                    xx = x[:22] + "%4d" % nst + x[26:]
                    fw.write(xx)
                    # fp2.pop()


##########################################################
def new_ligand_water_residue_number(pdbfile):
    """re-assign the residue numbers for ligands/waters"""

    pdb = open(pdbfile, "r").readlines()
    pdb1, pdb2, pdb3, pdb4 = separate_pdb_list(pdb)  # poly/lig/water

    chain1, chain_range1 = chain_res_list(pdb1, 0)  # poly
    chain2, chain_range2 = chain_res_list(pdb2, 0)  # ligand
    chain3, chain_range3 = chain_res_list(pdb3, 1)  # wat

    out = pdbfile + "_new_num"
    fw = open(out, "w")

    print("Polymer=" % chain_range1)
    print("Ligands=" % chain_range2)
    print("Waters =" % chain_range3)

    nn = 0
    for x in chain2.values():
        if len(x) > nn:
            nn = len(x)  # max number of ligand

    for x in pdb1:
        fw.write(x)  # write polymer first

    change_nres(fw, chain_range1, chain_range2, pdb2, 2)
    change_nres(fw, chain_range1, chain_range3, pdb3, 3 + nn)

    for x in pdb4:
        fw.write(x)  # write tail

    fw.close()

    if len(open(out, "r").readlines()) != len(open(pdbfile, "r").readlines()):
        print("Warning: Size of the new/old pdb are different. Check the new pdb")

    print("\nNote: Residue number of Ligands/Waters have been re-assigned.")
    print("The new pdbfile=%s\n" % out)


##########################################################
def anisotropy(pdbfile):
    """Analyze anisoutropy of each atom."""

    fr = open(pdbfile, "r")
    output = pdbfile + ".anis"
    fo = open(output, "w")

    head = """
    Note: column 1 is atom id.
    column 2,3,4 are the principle axis (eigenvalues) of ellipsoid of
    each atom.
    column 5 is the anisotropy of each atom (defined as longest axis
    divided by shortest axis. The ratio is 1.0 for a perfectly isotropic
    (spherical) atom.

    use 'sort output -n -k 5' to sort anisotropy in order

"""

    fo.write(head)
    n = 0
    for ln in fr:
        if "ANISOU" in ln[:6]:
            n = n + 1
            atom = ln[12:26].lstrip().replace(" ", "_")

            u = [float(ln[27:35]), float(ln[35:42]), float(ln[42:49]), float(ln[49:56]), float(ln[56:63]), float(ln[63:70])]
            u = [x * 0.0001 for x in u]
            id, e1, e2, e3 = eigenvalue(u)
            if id == 0 or e2 == 0:
                print("Warning: anisotropy can not be determined for " + atom)
                continue

            anis = e2 / e1
            tmp = "%s %8.4f %8.4f %8.4f  %8.4f\n" % (atom, e1, e2, e3, anis)
            if e1 <= 0 or e2 <= 0 or e3 <= 0:
                print("Error: [%s] has negative eigenvalues" % ln[12:26])

            fo.write(tmp)

    #            all_ani.append(tmp)
    #            print(tmp)

    #    sorted=all_ani.sort(key= lambda ln: float(ln.split(' ')[3]))
    #    for ln in all_ani: fo.write('%s' %ln)
    #    fo.write('below is sorted\n')
    #    for ln in sorted: fo.write(ln)

    fo.close(), fr.close()
    if n == 0:
        print("Eigenvalues not calculated, since no ANISOU records in " + pdbfile)
        os.remove(output)
        return
    print("The output file = " + output)


##########################################################
def eigenvalue(u):
    """
    reference:http://en.wikipedia.org/wiki/Eigenvalue_algorithm


    get eigen values for a 3X3 matrix (u) that has 3 REAL eigen values

        a  b  c  u11 u12 u13 or u[0] u[3] u[4]
    u = d  e  f  u21 u22 u23    u[3] u[1] u[5]
        g  h  i  u31 u32 u33    u[4] u[5] u[2]

    det= -l**3 +L**2(a+e+i) + L(d*b+g*c+f*h-a*e-a*i-e*i) +
    (a*e*i -a*f*h -d*b*i +d*c*h + g*b*f -g*c*e)

    Eqn = a1L**3 + b1L**2 + c1L + d1 = 0
    """
    a, b, c = u[0], u[3], u[4]
    d, e, f = u[3], u[1], u[5]
    g, h, i = u[4], u[5], u[2]

    a1 = -1
    b1 = a + e + i
    c1 = d * b + g * c + f * h - a * e - a * i - e * i
    d1 = a * e * i - a * f * h - d * b * i + d * c * h + g * b * f - g * c * e

    a, b, c, d = a1, b1, c1, d1

    if a == 0:
        return 0, 0, 0, 0

    x = ((3 * c / a) - (b * b / a * a)) / 3.0
    y = ((2 * math.pow(b, 3) / math.pow(a, 3)) - 9 * b * c / (a * a) + (27 * d / a)) / 27.0
    z = y * y / 4.0 + math.pow(x, 3) / 27.0

    t1 = y * y / 4.0 - z
    if t1 <= 0:
        return 0, 0, 0, 0
    i = math.sqrt(t1)

    j = -math.pow(i, (1.0 / 3))
    t2 = -(y / (2 * i))
    if math.fabs(t2) > 1:
        return 0, 0, 0, 0
    k = math.acos(t2)

    m = math.cos(k / 3.0)
    n = math.sqrt(3) * math.sin(k / 3.0)
    p = -(b / 3 * a)

    Eig1 = -2 * j * m + p
    Eig2 = j * (m + n) + p
    Eig3 = j * (m - n) + p

    return 1, Eig1, Eig2, Eig3


##########################################################
def check_file(size, *files):
    n = 1
    for f in files:
        if not os.path.exists(f) or os.path.getsize(f) < size:
            print("Error: file (%s) not exist." % f)
            n = 0
    return n


##########################################################
def residue(resname):
    val = 0
    resid = [
        "GLY",
        "ALA",
        "VAL",
        "LEU",
        "ILE",
        "CYS",
        "MET",
        "PRO",
        "PHE",
        "TRP",
        "TYR",
        "HIS",
        "LYS",
        "ARG",
        "ASP",
        "GLU",
        "ASN",
        "GLN",
        "THR",
        "SER",
        "MSE",
        "UNK",
        "  A",
        "  G",
        "  C",
        "  T",
        " DA",
        " DG",
        " DC",
        " DT",
        "  U",
        "  I",
    ]

    if resname in resid:
        val = 1
    return val


##########################################################
def atom_record(x):
    val = 0
    if "ATOM" in x[:4] or "HETA" in x[:4] or "ANISOU" in x[:6]:
        val = 1

    return val


##########################################################
def mass_of_atom(atom):
    atom_mass = {
        "H": 1.0079,
        "HE": 4.0026,
        "LI": 6.941,
        "BE": 9.01218,
        "B": 10.81,
        "C": 12.011,
        "N": 14.0067,
        "O": 15.9994,
        "F": 18.9984,
        "NE": 20.179,
        "NA": 22.98977,
        "MG": 24.305,
        "AL": 26.98154,
        "SI": 28.0855,
        "P": 30.97376,
        "S": 32.06,
        "CL": 35.453,
        "AR": 39.948,
        "K": 39.0983,
        "CA": 40.08,
        "SC": 44.9559,
        "TI": 47.9,
        "V": 50.9415,
        "CR": 51.996,
        "MN": 54.938,
        "FE": 55.847,
        "CO": 58.9332,
        "NI": 58.7,
        "CU": 63.546,
        "ZN": 65.38,
        "GA": 69.72,
        "GE": 72.59,
        "AS": 74.9216,
        "SE": 78.96,
        "BR": 79.904,
        "KR": 83.8,
        "RB": 85.4678,
        "SR": 87.62,
        "Y": 88.9059,
        "ZR": 91.22,
        "NB": 92.9064,
        "MO": 95.94,
        "TC": 98.0,
        "RU": 101.07,
        "RH": 102.9055,
        "PD": 106.4,
        "AG": 107.868,
        "CD": 112.41,
        "IN": 114.82,
        "SN": 118.69,
        "SB": 121.75,
        "TE": 127.6,
        "I": 126.9045,
        "XE": 131.3,
        "CS": 132.9054,
        "BA": 137.33,
        "LA": 138.9055,
        "CE": 140.12,
        "PR": 140.9077,
        "ND": 144.24,
        "PM": 145.0,
        "SM": 150.4,
        "EU": 151.96,
        "GD": 157.25,
        "TB": 158.9254,
        "DY": 162.5,
        "HO": 164.9304,
        "ER": 167.26,
        "TM": 168.9342,
        "YB": 173.04,
        "LU": 174.967,
        "HF": 178.49,
        "TA": 180.9479,
        "W": 183.85,
        "RE": 186.207,
        "OS": 190.2,
        "IR": 192.22,
        "PT": 195.09,
        "AU": 196.9665,
        "HG": 200.59,
        "TL": 204.37,
        "PB": 207.2,
        "BI": 208.9804,
        "PO": 209.0,
        "AT": 210.0,
        "RN": 222.0,
        "FR": 223.0197,
        "RA": 226.0254,
        "AC": 227.0278,
        "TH": 232.0381,
        "PA": 231.0359,
        "U": 238.029,
        "NP": 237.0482,
        "PU": 244.0,
        "AM": 243.0,
        "CM": 247.0,
        "BK": 247.0,
        "CF": 251.0,
        "ES": 252.0,
        "FM": 257.0,
        "MD": 258.0,
        "NO": 259.0,
        "LR": 260.0,
        "RF": 261.0,
        "DB": 262.0,
        "SG": 263.0,
        "BH": 262.0,
        "HS": 265.0,
        "MT": 266.0,
        "CN": 277.0,
        "X": 12.011,
    }
    # for x in atom_mass.keys(): print(x, atom_mass[x])

    if atom.upper() in atom_mass.keys():
        return atom_mass[atom.upper()]
    else:
        print("Warning: The atom (%s) not exist in the list" % atom)
        return 0


##########################################################
if __name__ == "__main__":
    process(sys.argv)
##########################################################
