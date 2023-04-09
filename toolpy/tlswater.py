# ==================================================
#   This module is tls related
# ==================================================

import os, sys, math
import util, prog, parse


def compare_r_tls_water(pdbfile_inp, sffile):
    """compare the R factors with/without TLS correction for waters.
    Remove all the ANISOU and apply TLS again with/without water.
    """

    pdbfile = rid_of_anisou(pdbfile_inp)

    pdbfile_new = find_tls_water(pdbfile)
    out2 = prog.run_dcc(pdbfile_new, sffile, "")

    util.delete_file(pdbfile)

    return


##########################################################
def find_tls_water(pdbfile):
    """use dic to contain chain-ID and residue range"""
    fr = open(pdbfile, "r")

    ter, tmp, pdb, hoh, hoh_orig = [], [], [], [], []
    for ln in fr:
        if "ATOM" in ln[:4] or "HETATM" in ln[:6] or "TER" in ln[:3]:
            if "HOH" not in ln[17:20]:
                if "TER" in ln[:3]:  # use TER as one list
                    ter.append(tmp)
                    tmp = []

                tmp.append(ln)
                pdb.append(ln)

            elif "HOH" in ln[17:20]:
                hoh.append(ln)
                hoh_orig.append(ln)

    if len(ter) == 0:
        ter.append(pdb)
    fr.close()

    newpdb = newpdb_tls4water(pdbfile, pdb, hoh, hoh_orig)

    return newpdb


##########################################################
def newpdb_tls4water(pdbfile, pdb, hoh, hoh_orig):
    """
    pdb: only contains xyz of non-HOH atoms.
    hoh: only contains xyz of HOH atoms. (hoh_orig, a duplicate)

    get a new PDB file with the following modifications:
    1. sort waters around the residue range of each TLS group.The selection
    criteria is within 3.5A.
    2. Give a new chain W for waters and sort the waters in sequential order
    e.g. group 1 (W 1  W n1), group 2 (W n1+1 W n2) ....

    """

    pdbfile_new = pdbfile + "_water"

    fr = open(pdbfile, "r")
    fo = open(pdbfile_new, "w")

    wnew = []
    nwa = 1
    for ln in fr:
        if "REMARK   3    RESIDUE RANGE :" in ln:
            tmp = ln[29:].split()
            if len(tmp) != 4 or util.is_digit(tmp[1]) == 0 or util.is_digit(tmp[3]) == 0:
                print("Error! Residue range (%s) is wrong.\n" % ln[30:60].strip())
            else:
                if int(tmp[3]) < int(tmp[1]):
                    print("Error: Residue range (%s) is wrong." % ln[30:60].strip())
                    return pdbfile_new

                ch, rmin, rmax = tmp[0], int(tmp[1]), int(tmp[3])
                wres = waters_around_tls_group(ch, rmin, rmax, pdb, hoh)
                #                print (ln,wres )
                nw1 = nwa
                nw2 = nwa + len(wres) - 1
                for n1 in wres:
                    for n2 in hoh_orig:
                        if n1 == n2[22:26]:
                            s = n2[:21] + "W" + "%4d" % nwa + n2[26:]
                            wnew.append(s)
                            nwa = nwa + 1

                s1 = "REMARK   3    RESIDUE RANGE :" + "%4s%6d%9s%6d  \n" % ("W", nw1, "W", nw2)

                fo.write(s1)
            #                print(s1)
            fo.write(ln)
        else:
            if "HOH" in ln[17:20] or "END" in ln[:3]:
                continue
            fo.write(ln)

    #    fo.seek(0)
    for ln in wnew:
        fo.write(ln)
    for ln in hoh:
        fo.write(ln)

    fo.close(), fr.close()
    return pdbfile_new


##########################################################
def waters_around_tls_group(ch, rmin, rmax, pdb, hoh):
    """get all the waters around the TLS group (chain, rmin, rmax)
    selection criteria 3.5A around the TLS.
    wres contains the residue number.
    """

    wres = []
    for ln in pdb:
        chain, nres = ln[21:22], int(ln[22:26])
        if chain == ch and nres >= rmin and nres <= rmax:
            #            print (ch, rmin, rmax, chain, nres)
            if len(ln.strip()) < 54:
                continue
            check_water(wres, ln, hoh)

    return wres


##########################################################
def check_water(wres, line, hoh):
    """check water to see if they are arounding the heavy atoms"""

    dwater = 3.0

    xp, yp, zp = float(line[30:38]), float(line[38:46]), float(line[46:54])
    n = 0
    for ln in hoh:
        xh, yh, zh = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
        dx, dy, dz = xp - xh, yp - yh, zp - zh
        dist = math.sqrt(dx * dx + dy * dy + dz * dz)
        if dist <= dwater and dist > 1:  # water found
            wres.append(ln[22:26])
            hoh.pop(n)
            break
        n = n + 1


##########################################################
def rid_of_anisou(pdbfile):
    pdbnew = pdbfile + "_new"
    fr = open(pdbfile, "r")
    fo = open(pdbnew, "w")

    for ln in fr:
        if "ANISOU" not in ln[:6]:
            fo.write(ln)

    fo.close()
    fr.close()
    return pdbnew
