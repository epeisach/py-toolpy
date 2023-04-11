# ===========================================================
# this module is a interface to parse values from a log file
# ===========================================================

import util
import cifparse as cif


###############################################################################
def cell(flist):
    """return cell values as a list of float!"""

    mycell = [0, 0, 0, 0, 0, 0]
    items, values = cif.cifparse(flist, "_cell.")
    if not len(items):
        return cell

    a = cif.parse_values(items, values, "_cell.length_a")
    b = cif.parse_values(items, values, "_cell.length_b")
    c = cif.parse_values(items, values, "_cell.length_c")
    alpha = cif.parse_values(items, values, "_cell.angle_alpha")
    beta = cif.parse_values(items, values, "_cell.angle_beta")
    gamma = cif.parse_values(items, values, "_cell.angle_gamma")

    if not (a and b and c and alpha and beta and gamma):
        print("Warning: cells not extracted. Check ciftokens")

    for i, x in enumerate([a, b, c, alpha, beta, gamma]):
        if len(x) == 0 or not util.is_number(x[0]):
            print("Error: cell has wrong (%s) values" % x)
            continue
        mycell[i] = float(x[0].strip())

    return mycell


###############################################################################
def space_group_name(flist):
    """get space_group, return a string"""

    spg = ""
    items, values = cif.cifparse(flist, "_symmetry.")
    symm = cif.parse_values(items, values, "_symmetry.space_group_name_H-M")
    if symm:
        spg = symm[0].replace("'", "").replace('"', "").strip()
    else:
        print("Warning: space group not extracted. Check ciftokens")

    return spg


###############################################################################
def ncs_matrix(flist):
    """get NCS matrix, return a list of 3X4 matrix"""
    mtrix = []
    items, values = cif.cifparse(flist, "_struct_ncs_oper.")
    if not items:
        return []

    id = cif.parse_values(items, values, "_struct_ncs_oper.id")  # pylint: disable=redefined-builtin
    code = cif.parse_values(items, values, "_struct_ncs_oper.code")

    b11 = cif.parse_values(items, values, "_struct_ncs_oper.matrix[1][1]")
    b12 = cif.parse_values(items, values, "_struct_ncs_oper.matrix[1][2]")
    b13 = cif.parse_values(items, values, "_struct_ncs_oper.matrix[1][3]")

    b21 = cif.parse_values(items, values, "_struct_ncs_oper.matrix[2][1]")
    b22 = cif.parse_values(items, values, "_struct_ncs_oper.matrix[2][2]")
    b23 = cif.parse_values(items, values, "_struct_ncs_oper.matrix[2][3]")

    b31 = cif.parse_values(items, values, "_struct_ncs_oper.matrix[3][1]")
    b32 = cif.parse_values(items, values, "_struct_ncs_oper.matrix[3][2]")
    b33 = cif.parse_values(items, values, "_struct_ncs_oper.matrix[3][3]")

    t1 = cif.parse_values(items, values, "_struct_ncs_oper.vector[1]")
    t2 = cif.parse_values(items, values, "_struct_ncs_oper.vector[2]")
    t3 = cif.parse_values(items, values, "_struct_ncs_oper.vector[3]")

    if not (id and code and b11 and b12 and b13 and t1 and b21 and b22 and b23 and t2 and b31 and b32 and b33 and t3):
        return []

    for i in range(len(b11)):
        # idd = " "
        # if code[i] == "given":
        #     idd = "1"
        mt1 = [float(x) for x in (b11[i], b12[i], b13[i], t1[i]) if util.is_number(x)]
        mt2 = [float(x) for x in (b21[i], b22[i], b23[i], t2[i]) if util.is_number(x)]
        mt3 = [float(x) for x in (b31[i], b32[i], b33[i], t3[i]) if util.is_number(x)]
        if not (mt1 and mt2 and mt3):
            print("Error: NCS matrix has wrong values")
            continue
        mtrix.append([mt1, mt2, mt3])

    return mtrix


###############################################################################
def scale(flist):
    """get SCALE matrix, return a list of 3X4 matrix"""
    outscale = []
    items, values = cif.cifparse(flist, "_atom_sites.")
    if not items:
        return []

    b11 = cif.parse_values(items, values, "_atom_sites.fract_transf_matrix[1][1]")
    b12 = cif.parse_values(items, values, "_atom_sites.fract_transf_matrix[1][2]")
    b13 = cif.parse_values(items, values, "_atom_sites.fract_transf_matrix[1][3]")
    t1 = cif.parse_values(items, values, "_atom_sites.fract_transf_vector[1]")

    b21 = cif.parse_values(items, values, "_atom_sites.fract_transf_matrix[2][1]")
    b22 = cif.parse_values(items, values, "_atom_sites.fract_transf_matrix[2][2]")
    b23 = cif.parse_values(items, values, "_atom_sites.fract_transf_matrix[2][3]")
    t2 = cif.parse_values(items, values, "_atom_sites.fract_transf_vector[2]")

    b31 = cif.parse_values(items, values, "_atom_sites.fract_transf_matrix[3][1]")
    b32 = cif.parse_values(items, values, "_atom_sites.fract_transf_matrix[3][2]")
    b33 = cif.parse_values(items, values, "_atom_sites.fract_transf_matrix[3][3]")
    t3 = cif.parse_values(items, values, "_atom_sites.fract_transf_vector[3]")

    if not (b11 and b12 and b13 and t1 and b21 and b22 and b23 and t2 and b31 and b32 and b33 and t3):
        return []

    for i in range(len(b11)):
        mt1 = [float(x) for x in (b11[i], b12[i], b13[i], t1[i]) if util.is_number(x)]
        mt2 = [float(x) for x in (b11[i], b12[i], b13[i], t1[i]) if util.is_number(x)]
        mt3 = [float(x) for x in (b11[i], b12[i], b13[i], t1[i]) if util.is_number(x)]
        if not (mt1 and mt2 and mt3):
            print("Error: Scale matrix has wrong values")
            continue
        outscale.append([mt1, mt2, mt3])

    return outscale


###############################################################################
def atom_site(flist):
    """get all atom record. return a dictionary"""

    dic = {}
    chain, nres, comp, atom, symbol, alt, ins = [], [], [], [], [], [], []
    x, y, z, biso, occ = [], [], [], [], []

    items, values = cif.cifparse(flist, "_atom_site.")  # a loop

    # group1 = cif.parse_values(items, values, "_atom_site.group_PDB")
    natm = cif.parse_values(items, values, "_atom_site.id")
    symbol1 = cif.parse_values(items, values, "_atom_site.type_symbol")
    if not symbol1:
        symbol1 = cif.parse_values(items, values, "_atom_site.atom_type_symbol")
    atom1 = cif.parse_values(items, values, "_atom_site.label_atom_id")
    asym1 = cif.parse_values(items, values, "_atom_site.auth_asym_id")
    comp1 = cif.parse_values(items, values, "_atom_site.label_comp_id")
    nres1 = cif.parse_values(items, values, "_atom_site.auth_seq_id")
    x1 = cif.parse_values(items, values, "_atom_site.Cartn_x")
    y1 = cif.parse_values(items, values, "_atom_site.Cartn_y")
    z1 = cif.parse_values(items, values, "_atom_site.Cartn_z")
    biso1 = cif.parse_values(items, values, "_atom_site.B_iso_or_equiv")
    occ1 = cif.parse_values(items, values, "_atom_site.occupancy")
    ins1 = cif.parse_values(items, values, "_atom_site.pdbx_PDB_ins_code")
    if not ins1:
        ins1 = cif.parse_values(items, values, "_atom_site.ndb_ins_code")
    alt1 = cif.parse_values(items, values, "_atom_site.label_alt_id")

    if not (atom1 and comp1 and asym1 and nres1 and x1 and y1 and z1 and occ1 and biso1):
        print("Error: there is problem to parse atom_site.")
        return dic

    n = len(x1)

    for i in range(n):
        chain.append(asym1[i])
        nres.append(int(nres1[i]))
        comp.append(comp1[i])
        atom.append(atom1[i])
        if natm:
            natm.append(natm[i])  # newly added
        symbol.append(symbol1[i])
        if not alt1:
            alt.append(".")
        else:
            alt.append(alt1[i])
        if not ins1:
            ins.append(".")
        else:
            ins.append(ins1[i])

        x.append(float(x1[i]))
        y.append(float(y1[i]))
        z.append(float(z1[i]))
        biso.append(float(biso1[i]))
        occ.append(float(occ1[i]))

    dic = {"chain": chain, "nres": nres, "comp": comp, "atom": atom, "symbol": symbol, "alt": alt, "ins": ins, "x": x, "y": y, "z": z, "biso": biso, "occ": occ, "natm": natm}

    return dic


###############################################################################


def val_from_dcc(outf):
    """The format of outf can not be changed !!
    index (0,1,2,3,4,5,6,7...) corresponds to
    (resh,resl,Rwork,Rfree,Comp,FCC, Real_R, Dcc)
    res_rep: the reported; res_not: no TLS; res_tls: with TLS
    """

    if not util.check_file(500, outf):
        return "", ""

    flist = open(outf, "r").readlines()

    items, values = cif.cifparse(flist, "_pdbx_density.")
    # res = cif.parse_values(items, values, "_pdbx_density.ls_d_res_high")
    # rw = cif.parse_values(items, values, "_pdbx_density.R_value_R_work")
    # rf = cif.parse_values(items, values, "_pdbx_density.R_value_R_free")
    biso = cif.parse_values(items, values, "_pdbx_density.Biso_mean")
    bwil = cif.parse_values(items, values, "_pdbx_density.B_wilson")
    l2 = cif.parse_values(items, values, "_pdbx_density.Padilla-Yeates_L2_mean")
    z = cif.parse_values(items, values, "_pdbx_density.Z_score_L_test")
    fom = cif.parse_values(items, values, "_pdbx_density.fom")

    items, values = cif.cifparse(flist, "_pdbx_density_corr.")
    prog = cif.parse_values(items, values, "_pdbx_density_corr.program")
    resh = cif.parse_values(items, values, "_pdbx_density_corr.ls_d_res_high")
    rwork = cif.parse_values(items, values, "_pdbx_density_corr.ls_R_factor_R_work")
    rfree = cif.parse_values(items, values, "_pdbx_density_corr.ls_R_factor_R_free")
    fcc = cif.parse_values(items, values, "_pdbx_density_corr.correlation_coeff_Fo_to_Fc")
    rsr = cif.parse_values(items, values, "_pdbx_density_corr.real_space_R")
    dcc = cif.parse_values(items, values, "_pdbx_density_corr.correlation")
    detail = cif.parse_values(items, values, "_pdbx_density_corr.details")

    nr, nb = 0, 0
    for i, x in enumerate(detail):
        if "Best" in x:
            nb = i
            break

    rep = "%8s  %4s  %6s  %6s  %6s " % (prog[nr], resh[nr], rwork[nr], rfree[nr], biso[0])
    val = "%8s  %4s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s " % (prog[nb], resh[nb], rwork[nb], rfree[nb], bwil[0], l2[0], z[0], fom[0], fcc[nb], rsr[nb], dcc[nb])

    return rep, val


########################################################################
def values_from_dcc(outf):
    """The format of outf can not be changed !!
    index (0,1,2,3,4,5,6,7...) corresponds to
    (resh,resl,Rwork,Rfree,Comp,FCC, Real_R, Dcc)
    res_rep: the reported; res_not: no TLS; res_tls: with TLS
    """

    bstatus, prog, prog_rep, value = "?", "REFMAC", "?", "?"

    res_rep = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
    res_not = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
    res_tls = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]

    if not util.check_file(500, outf):
        return value, res_rep, res_not, res_tls, bstatus, prog_rep

    fr = open(outf, "r")
    for line in fr:
        if "'PDB reported'" in line:
            get_value(res_rep, line)
            prog_rep = line.split()[1].strip()

        elif prog in line and "Best " in line:
            get_value(res_not, line)

        elif "_pdbx_density.iso_B_value_type" in line:
            bstatus = "P"
            if " FULL" in line:
                bstatus = "T"

    fr.close()

    value = "%4.2f  %5.3f %5.3f %5.3f   %5.3f %5.3f %5.3f " % (res_rep[0], res_rep[2], res_rep[3], res_rep[5], res_not[2], res_not[3], res_not[5])

    return value, res_rep, res_not, res_tls, bstatus, prog_rep


########################################################################


def get_value(res_rep, line):
    tmp = line.split()
    if "." in tmp[2]:
        res_rep[0] = float(tmp[2])  # resh
    if "." in tmp[3]:
        res_rep[1] = float(tmp[3])  # resl
    if "." in tmp[5]:
        res_rep[2] = float(tmp[5])  # Rwork
    if "." in tmp[4] and "?" in tmp[5]:
        res_rep[2] = float(tmp[4])  # Rwork
    if "." in tmp[6]:
        res_rep[3] = float(tmp[6])  # Rfree
    if "." in tmp[8]:
        res_rep[4] = float(tmp[8])  # Comp
    if "." in tmp[10]:
        res_rep[5] = float(tmp[10])  # FCC
    if "." in tmp[11]:
        res_rep[6] = float(tmp[11])  # Real_R
    if "." in tmp[12]:
        res_rep[7] = float(tmp[12])  # Dcc

    return res_rep
