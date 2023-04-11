# ===========================================================
# this module contains the utility functions
# ===========================================================


import os
import math
import shutil
import sys

# from  config import *


##########################################################
def sort_column(data, col):
    data.sort(key=lambda y: y[col])


##########################################################
def gsort(data, col, idd):
    """A general function to sort a list (col=-1) or a list of list (col>=0)
    idd : 1 for sorting from small to large; -1 from large to small
    """

    if col < 0:  # a single list [..]
        if idd >= 0:
            data.sort()
        else:
            data.sort(reverse=True)

    else:  # a list of list [[..], [..]]
        if idd >= 0:
            data.sort(key=lambda y: y[col])
        else:
            data.sort(key=lambda y: y[col], reverse=True)


##########################################################
def chain_res_atom(pdb):
    """parse the pdb  into chain-residnum-atom list. a dictionary
    dd[x][y]: atom list, x:chainID, y:residue_number_in_string,
    dd={ch:{resn:[atom_lines]}}
    """

    if check_file(300, pdb) == 0:
        return

    fp = open(pdb, "r")

    id = 0  # pylint: disable=redefined-builtin
    dd = {}
    for x in fp:
        if not ("ATOM" in x[:4] or "HETATM" in x[:6] or (x.strip()) < 50):
            continue
        ch = x[20:22].strip()
        if ch == " ":
            id = 1

        if ch not in dd.keys() and ch != " ":
            dd[ch] = {}
        res = x[22:27]  # include inserted
        if ch != " " and ch in dd.keys() and res not in dd[ch].keys():
            dd[ch][res] = []
        if ch != " ":
            dd[ch][res].append(x)

    if id == 1:
        print("Warning: file (%s) has no ChainID." % pdb)
    return dd


##########################################################
def pdb2cif(xyz):
    """not sure why maxit can not get the refine. table
    add it from pdb_extract
    """

    os.system("maxit_local  -i %s  -o 1  " % xyz)
    cif = xyz + ".cif"
    if not check_file(200, cif):
        print("Error: pdb2cif conversion failed.")
        return " "

    cif_list = open(cif, "r").readlines()
    n = 0
    for x in cif_list:
        if "_refine." in x:
            n = n + 1
    if n < 4:  # add this table by pdb_extract
        out = xyz + "_TMP"
        os.system("pdb_extract -r refmac -ipdb %s -o %s >/dev/null" % (xyz, out))

        refine = []
        fp = open(out, "r").readlines()
        for x in fp:
            if "_refine.entry_id" in x or "_refine.details" in x:
                continue
            if "_refine." in x:
                refine.append(x)
        #        delete_file(out)

        fw = open(cif, "w")
        for y in cif_list:
            fw.write(y)
        for y in refine:
            fw.write(y)
        fw.close()

    return cif


##########################################################
def delete_file(*files):
    for x in files:
        os.system("rm -f " + x)


##########################################################
def check_file(size, *files):
    n = 1
    for f in files:
        if not os.path.exists(f) or os.path.getsize(f) < size:
            #            print('Error: file (%s) does not exist (or file size=0).' %f)
            n = 0
            break
    return n


##########################################################
def str_after_id(line, id):  # pylint: disable=redefined-builtin
    """get string after a given id"""

    if id not in line:
        print("Warning: %s not in string (%s)" % (id, line))
        return line

    n = line.index(id)
    value = line[n + len(id) :].strip()
    if len(value) <= 0 or value.upper() == "NULL":
        value = "?"
    return value


##########################################################


def float_after_id(line, id):  # pylint: disable=redefined-builtin
    """get float after a given id"""

    if id not in line:
        print("Warning: %s not in string." % (id))
        return 9999.0

    n = line.index(id)
    li = line[n + 1 :].strip().split()
    if len(li) <= 0:
        print("Error! No value after id (%s)." % (id))
        return 9999.0
    else:
        if isnum(li[0]):
            return float(li[0])
        else:
            print("Error! %s is not a numeric." % li[0])
            return 9999.0


##########################################################


def int_after_id(line, id):  # pylint: disable=redefined-builtin
    """get int after a given id"""

    if id not in line:
        print("Warning: %s not in string." % (id))
        return 9999

    n = line.index(id)
    li = line[n + 1 :].strip().split()
    if len(li) <= 0:
        print("Error! No value after id (%s)." % (id))
        return 9999
    else:
        if isnum(li[0]):
            return int(li[0])
        else:
            print("Error! %s is not a numeric." % li[0])
            return 9999


##########################################################
def str_between_id(line, id1, id2):
    """get string between the two ids"""

    if id1 not in line or id2 not in line:
        print("Warning: either %s or %s not in string" % (id1, id2))
        return "?"
    n1, n2 = line.index(id1) + len(id1), line.index(id2)
    return line[n1 + 1 : n2].strip()


##########################################################
def float_between_id(line, id1, id2):
    """get float value between the two ids"""

    if id1 not in line or id2 not in line:
        print("Warning: either %s or %s not in string" % (id1, id2))
        return 9999.0
    n1, n2 = line.index(id1), line.index(id2)
    return float(line[n1 + 1 : n2])


##########################################################
def int_between_id(line, id1, id2):
    """get int value between the two ids"""

    if id1 not in line or id2 not in line:
        print("Warning: either %s or %s not in string" % (id1, id2))
        return 9999
    n1, n2 = line.index(id1), line.index(id2)
    return int(line[n1 + 1 : n2])


##########################################################
def isnum(value):
    return str(value).replace(".", "").replace("-", "").isdigit()


##########################################################
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        # print('Error: %s is not a number.' %s)
        return False


##########################################################
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
def get_value_after_id(line, id):  # pylint: disable=redefined-builtin
    """get value after a given id"""

    if id not in line:
        return "?"
    li = line.split(id)
    value = li[1].strip()
    if len(value) <= 0 or value.upper() == "NULL":
        value = "?"
    return value


##########################################################
def is_cif(file):
    n, m = 0, 0
    if not os.path.exists(file):
        return n

    for x in open(file, "r").readlines():
        m = m + 1
        if m > 100:
            break
        t = x.strip()
        if len(t) < 4:
            continue

        if "data_" in t[:5] or "loop_" in t[:5] or ("_" == t.split()[0][0] and "." in t.split()[0]):
            n = 1
            break
        elif "HEADER " in x[:7]:
            break

    return n


##########################################################
def value_between_string(ln, s1, s2):
    """get the string between two strings (s1,s2) in a line(case sensitive)"""
    if s1 not in ln or s2 not in ln:
        print("Warning: string (%s or %s) not in %s" % (s1, s2, ln))
        return ""
    n = ln.index(s1) + len(s1)
    t1 = ln[n:]
    n = t1.index(s2)
    return t1[:n]


##########################################################
def residue():
    """non-ligand residues"""

    res = [
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
        "HOH",
        "DOD",
        "  A",
        "  G",
        "  C",
        "  T",
        " DA",
        " DG",
        " DC",
        " DT",
        " DI",
        "  U",
        "  I",
        "ADE",
        "THY",
        "GUA",
        "CYT",
        "URA",
    ]

    return res


# ##########################################################
# def perror(info):
#     """print error messages"""
#     if info not in ERRLOG:
#         ERRLOG.append(info)
#         print(info.strip())


##########################################################
def space_group_crystal_match(spg):
    """
    The dic spg_cell contains all the chiral space groups.

    http://www.ruppweb.org/Xray/comp/space_instr.htm
    """
    spg_cryst = {
        # TRICLINIC
        "A1": 1,
        "P1": 1,
        "P1-": 1,
        "P-1": 1,
        # MONOCLINIC
        "A121": 2,
        "A2": 2,
        "B112": 20,  # C is unique axis. alpha, beta, =90
        "B2": 20,  # C is unique axis. alpha, beta, =90
        "C2": 2,
        "C121": 2,
        "C21": 2,
        "C1211": 2,
        "C2(A112)": 2,
        "I121": 2,
        "I1211": 2,
        "I2": 2,
        "I21": 2,
        "P2": 2,
        "P121": 2,
        "P112": 2,
        "P21": 2,
        "P1211": 2,
        "P1121": 2,
        "P21(C)": 2,
        # ORTHORHOMBIC
        "P222": 3,
        "P2221": 3,
        "P2212": 3,
        "P21212": 3,
        "P212121": 3,
        "P22121": 3,
        "P21221": 3,
        "P21212A": 3,
        "B2212": 3,
        "C222": 3,
        "C2221": 3,
        "I222": 3,
        "I212121": 3,
        "F222": 3,
        # TETRAGONAL
        "P4": 4,
        "P41": 4,
        "P42": 4,
        "P43": 4,
        "P422": 4,
        "P4212": 4,
        "P4122": 4,
        "P41212": 4,
        "P4222": 4,
        "P42212": 4,
        "P4322": 4,
        "P43212": 4,
        "I4": 4,
        "I41": 4,
        "I422": 4,
        "I4122": 4,
        # TRIGONAL   #hexagonal axis a=b gamma=120
        "P3": 5,
        "P31": 5,
        "P32": 5,
        "P312": 5,
        "P321": 5,
        "P3112": 5,
        "P3121": 5,
        "P3212": 5,
        "P3221": 5,
        "R3": 50,  # (rhombohedral axis, a=b=c & alpha=beta=gamma)
        "R32": 50,  # (rhombohedral axis, a=b=c & alpha=beta=gamma)
        "H3": 5,
        "H32": 5,
        # HEXAGONAL
        "P6": 6,
        "P61": 6,
        "P62": 6,
        "P63": 6,
        "P64": 6,
        "P65": 6,
        "P622": 6,
        "P6122": 6,
        "P6222": 6,
        "P6322": 6,
        "P6422": 6,
        "P6522": 6,
        # CUBIC
        "C4212": 7,  # ?
        "F422": 7,  # ?
        "P23": 7,
        "F23": 7,
        "I23": 7,
        "P213": 7,
        "I213": 7,
        "P432": 7,
        "P4132": 7,
        "P4232": 7,
        "P4332": 7,
        "F432": 7,
        "F4132": 7,
        "I432": 7,
        "I4132": 7,
        # others
        "P121/c1": 2,
        "I41/a": 4,
    }
    tmp = spg.replace(" ", "")
    if tmp in spg_cryst.keys():
        return spg_cryst[tmp]
    else:
        print("Warning: The space group (%s) is not in the list" % spg)
        return 0


##########################################################


def get_file_by_pdbid(pdbid_in, idd):
    """pdbid_in : 4 char PDBID;   idd='cifid|pdbid|cif|pdb|sf'"""

    #    pth='/net/data/remediation-alt/ftp-v4.0/pdb/data/structures/divided'
    pth = "/net/ftp_tree_v5/ftp-v5.0/pdb/data/structures/divided/"
    # www_path = "http://www.rcsb.org/pdb/files"

    pdbid = pdbid_in.lower()
    if len(pdbid) != 4:
        print("Error: PDBID (%s) is not 4 characters. " % pdbid)
        sys.exit()

    hash = pdbid[1:3]  # pylint: disable=redefined-builtin

    cif = "%s/mmCIF/%s/%s.cif.gz" % (pth, hash, pdbid)
    pdb = "%s/pdb/%s/pdb%s.ent.gz" % (pth, hash, pdbid)
    sf = "%s/structure_factors/%s/r%ssf.ent.gz" % (pth, hash, pdbid)

    sffile = "r" + pdbid + "sf.ent"
    pdbfile = "pdb" + pdbid + ".ent"
    ciffile = pdbid + ".cif"

    if idd == "pdb":
        os.system("zcat  %s > %s " % (pdb, pdbfile))
        if not check_file(10, pdbfile):
            os.system("zcat  %s > %s " % (cif, pdbfile))  # for large entry
        return pdbfile, ""
    elif idd == "sf":
        os.system("zcat  %s > %s " % (sf, sffile))
        return "", sffile
    elif idd == "cif":
        os.system("zcat  %s > %s " % (cif, ciffile))
        return ciffile, ""
    elif idd == "cifid":
        os.system("zcat  %s > %s " % (cif, ciffile))
        os.system("zcat  %s > %s " % (sf, sffile))
        return ciffile, sffile
    elif idd == "pdbid":
        os.system("zcat  %s > %s " % (pdb, pdbfile))
        os.system("zcat  %s > %s " % (sf, sffile))
        return pdbfile, sffile
    else:
        print("Error: wrong file format. ")
        sys.exit()


##########################################################
def mean_dev(data_in, col):
    """get the min, max, average and deviation
    if col=-1, data_in is a list of data, else a list of list!
    """
    if len(data_in) == 0:
        return 0, 0, 0, 0
    data = []
    if col < 0:
        data = data_in
    else:
        for x in data_in:
            data.append(float(x[col]))
    # print(data)
    mini, maxi = min(data), max(data)
    n = len(data)
    if n == 0:
        return 0, 0, 0, 0
    avg = sum(data) / float(n)
    s1 = 0
    for x in data:
        a2 = (x - avg) * (x - avg)
        s1 = s1 + a2
    dev = math.sqrt(s1 / float(n))

    return avg, dev, mini, maxi


##########################################################
def disto(x1, x2):
    """x1,x2 must be in orthogonal"""

    d = math.sqrt((x1[0] - x2[0]) ** 2 + (x1[1] - x2[1]) ** 2 + (x1[2] - x2[2]) ** 2)
    return d


##########################################################
def distf(x1, x2, cell):
    """x1,x2 must be in fractional"""

    d2r = 3.141592654 / 180

    dx = x1[0] - x2[0]
    dy = x1[1] - x2[1]
    dz = x1[2] - x2[2]

    ca = math.cos(cell[3] * d2r)
    cb = math.cos(cell[4] * d2r)
    cg = math.cos(cell[5] * d2r)

    d2 = (
        (cell[0] * dx) ** 2
        + (cell[1] * dy) ** 2
        + (cell[2] * dz) ** 2
        + 2 * cell[1] * cell[2] * dy * dz * ca
        + 2 * cell[0] * cell[2] * dx * dz * cb
        + 2 * cell[0] * cell[1] * dx * dy * cg
    )

    d = math.sqrt(d2)
    return d


##########################################################
def frac_orth_matrix(cell):
    """get fractional and orthogonal matrix to apply xyz
    cell is a list containing of six floats
    """

    frac = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    orth = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

    if not cell:
        return frac, orth

    d2r = 3.141592654 / 180
    a = cell[0]
    b = cell[1]
    c = cell[2]

    sa = math.sin(cell[3] * d2r)
    sb = math.sin(cell[4] * d2r)
    sg = math.sin(cell[5] * d2r)

    ca = math.cos(cell[3] * d2r)
    cb = math.cos(cell[4] * d2r)
    cg = math.cos(cell[5] * d2r)

    vol = sg * sb * sa
    frac[0][0] = 1.0 / a
    frac[0][1] = (-cg) / (a * sg)
    frac[0][2] = (cg * ca - cb) / (a * vol * sg)
    frac[1][0] = 0.0
    frac[1][1] = 1.0 / (b * sg)
    frac[1][2] = (cg * cb - ca) / (b * vol * sg)
    frac[2][0] = 0.0
    frac[2][1] = 0.0
    frac[2][2] = (sg) / (c * vol)

    orth[0][0] = a
    orth[0][1] = b * cg
    orth[0][2] = c * cb
    orth[1][0] = 0.0
    orth[1][1] = b * sg
    orth[1][2] = c * (ca - cb * cg) / sg
    orth[2][0] = 0.0
    orth[2][1] = 0.0
    orth[2][2] = c * sb * sa

    return frac, orth


##########################################################
def matrix_prod_(A, B):
    """input matrix A=m1Xn ;  B=nXm2 ::  C = m1Xm2 is the output"""
    rows_A = len(A)
    cols_A = len(A[0])
    rows_B = len(B)
    cols_B = len(B[0])

    if cols_A != rows_B:
        print("Error: Incorrect dimensions.(must be A=m1Xn;B=nXm2)")
        return

    # Create the result matrix
    # Dimensions would be rows_A x cols_B
    C = [[0 for row in range(cols_B)] for col in range(rows_A)]
    #    print C

    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):
                C[i][j] += A[i][k] * B[k][j]
    return C


##########################################################
def matrix_prod(f, x):
    """f: a list of list [[],[],[]];  x: a list []"""

    xo = [0, 0, 0]

    xo[0] = x[0] * f[0][0] + x[1] * f[0][1] + x[2] * f[0][2]
    xo[1] = x[0] * f[1][0] + x[1] * f[1][1] + x[2] * f[1][2]
    xo[2] = x[0] * f[2][0] + x[1] * f[2][1] + x[2] * f[2][2]

    return [xo[0], xo[1], xo[2]]


##########################################################
def move(finp, fout):
    if os.path.exists(finp):
        shutil.move(finp, fout)


############################################################
def gnu_plot(datafile, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot):
    """some head infor for gnuplot (one data set plot by default)
        bar>0: add error bar on the histogram
        rot>0, xlabel rotation;
        key=0: key>0, on;  key=2,left;
        style=0, histogram ; =1 linepoints; =2, points; =3 dot

        ------


    lt=linetypes color	pt

    -1	black		-1	n/a
    0	black		0	dotted
    1	red		1	+
    2	green		2	x
    3	blue		3	*
    4	magenta		4	empty square
    5	cyan		5	filled square
    6	brown		6	empty circle
    7	light brown	7	filled circle
    8	light red	8	empty triangle
            -----------
    9	red		9	filled triangle
    10	green		10	empty nabla
    11	blue		11	filled nabla
    12	magenta		12	empty rhombus
    13	cyan		13	filled rhombus
                                    --------------
    14	brown		14	+
    15	light brown	15	x

    example:
    set palette
    plot sin(x) lt -1

    """

    #    psfix='png' #must be the same as the dic
    psfix = "svg"  # must be the same as the dic

    gnuscr = datafile + ".gnu"
    gnuout = datafile + ".%s" % psfix
    fw = open(gnuscr, "w")

    #    x1='set terminal %s large size 1000,800 # 840,640\n' %psfix
    x1 = "set terminal svg enhanced size 840,640  \n"
    x2 = "set output '%s'\n" % gnuout
    fw.write(x1 + x2)
    if style == 0:  # histogram
        fw.write("set boxwidth 0.5 relative\n")
        fw.write("set style histogram clustered gap 1 title offset 0, 0, 0\n")
        fw.write("set style data histograms\n")
        if bar:
            fw.write("set style histogram errorbars\n")
    elif style == 1:  # line
        fw.write("set style data linespoints\n")
    elif style == 2:  # point
        fw.write("set style data points\n")
        fw.write("set pointsize  1\n")
    elif style == 3:  # dot
        fw.write("set style data dots\n")
        fw.write("set pointsize  1\n")

    if len(xrange) > 1:
        fw.write("set xrange %s\n" % xrange)
    if len(yrange) > 1:
        fw.write("set yrange %s\n" % yrange)
    if len(title) > 0:
        fw.write('set title "%s" \n' % title)
    fw.write("set datafile missing '.'\n")
    fw.write("set style fill  solid 1.0 border -1\n")

    if key == 0:
        fw.write("set key off\n")
    else:
        fw.write("set key on\n")
        if key == 2:
            fw.write("set key reverse left Left\n")

    fw.write("set grid ytics\n")
    fw.write("set size 1.0,1.0\n")
    fw.write("set autoscale\n")

    if rot > 0:
        fw.write("set xtics rotate by 90 \n")
    else:
        fw.write("set xtics\n")

    fw.write("set ytics\n")

    if len(xlabel) > 0:
        if rot > 0:
            fw.write('set xlabel "%s" offset 0, -1, 0\n' % xlabel)  # offset xlabel down -1
        else:
            fw.write('set xlabel "%s" offset 0, 0, 0\n' % xlabel)
    else:
        fw.write('set xlabel " " \n')

    if len(ylabel) > 0:
        fw.write('set ylabel "%s" \n' % ylabel)
    else:
        fw.write('set ylabel " " \n')

    fw.write("set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 bdefault\n")

    fw.write(plot)

    fw.write("""\n# plot 'file' u 1:3:xtic(1) t "?" with linespoints,'' using 1:2 t "?" with boxes\n""")
    fw.write("""# plot 'filename' using 1:2:3:xtic(1) t "all" with  errorbars \n""")
    fw.write("# set term x11 \n")
    fw.write("# replot \n")
    fw.write("# pause -6 \n")

    #    fw.write("set term x11\nreplot \npause 4\n")

    fw.close()
    arg = "gnuplot %s " % gnuscr
    os.system(arg)  # if not given /bin/csh -c, use B shell

    return gnuscr, gnuout
