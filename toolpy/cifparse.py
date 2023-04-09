##########################################################
# A cif parser in python (2013-01-06 ), HY
##########################################################

import os, sys, shutil, string, math


##########################################################
def cifparse(flist, cate):
    """parse cif tokens. flist is a list of content; cate is the category"""

    items, values = [], []
    # remove the empty string and front end space
    line = [x.strip() for x in flist if not x.isspace()]

    m, loop, start = 0, 0, 0
    nline, clen = len(line), len(cate)
    for i in range(0, nline):  # look for the table
        if line[i][:clen] == cate:
            if i > 0 and line[i - 1][:5] == "loop_":
                loop = 1
            start = i
            m = m + 1
            break  # find start position and loop or nonloop

    if m == 0:
        # print("Error: cif table (%s) is not found." %cate)
        return items, values

    if loop == 0:
        for i in range(start, nline):
            val = line[i].split()
            m = len(val)
            if m > 1 and val[0][:clen] == cate:  # ciftoken & value in same line
                ss = ""
                if val[1][0] == "'":
                    ss = string_between_char(line[i], "'", "'")
                elif val[1][0] == '"':
                    ss = string_between_char(line[i], '"', '"')
                else:
                    ss = val[1]

                items.append(val[0])
                values.append(ss)

            elif m == 1 and i < nline - 1 and val[0][:clen] == cate:
                ss = ""
                if line[i + 1][0] == ";":  # go through lines until hit ; or
                    ssnew = []
                    n = 0
                    for j in range(i + 1, nline):
                        if n > 0 and line[j][0] == ";":
                            start = j
                            break
                        ssnew.append(line[j])
                        n = n + 1
                        if line[j] == "loop_" or line[j] == "#" or line[j][:clen] == cate:
                            print('Error: can not find the second ";" in (%s). Stop here!' % cate)
                            start = j
                            break

                    ss = "\n".join(ssnew)
                    if ss[0] == ";":
                        ss = ss[1:]

                elif line[i + 1][0] == "'" or line[i + 1][0] == '"':
                    cc = line[i + 1][0]
                    ss = line[i + 1][1:]
                    if ss[-1] != cc:
                        print("Error: in table (%s). End of string is not (%s)" % (cate, cc))
                    else:
                        ss = ss[:-1]
                elif cate not in line[i + 1][:clen]:
                    ss = line[i + 1]

                items.append(val[0])
                values.append(ss)

            elif line[i][0] == "#" or line[i][:5] == "loop_" or (line[i][0] == "_" and cate not in line[i][:clen]):
                break

    else:  # parse looped items
        for i in range(start, nline):
            val = line[i].split()

            m = len(val)
            if m == 1 and val[0][:clen] == cate:  # get items
                items.append(val[0])
            elif (len(items) > 0 and i < nline - 1 and cate not in line[i + 1]) or i == nline - 1:
                start = i
                break

        j = start
        while j < nline:  # through list from start position
            if line[j][0] == ";":
                m, ss = lines_between_char(j, line, ";")
                values.append(ss)
                j = m

                continue
            elif line[j][0] == "#" or line[j][:clen] == "loop_" or (line[j][0] == "_" and "." in line[j]):
                break
            else:
                tmp = clean_str(line[j])

                m1 = 0
                while m1 < 800:
                    m1 = m1 + 1

                    if tmp[0] == "'" or tmp[0] == '"':
                        cc = tmp[0]
                        m, ss = string_after_char(0, tmp, cc)
                        tmp = tmp[m:].lstrip()
                    else:
                        m, ss = string_after_char(0, tmp, " ")
                        tmp = tmp[m:].lstrip()

                    values.append(ss)
                    if m == 0 or not tmp:
                        break
            j = j + 1

    ###finshed parsing
    nitem, nval = len(items), len(values)
    n = nval % nitem
    nrow = nval / nitem

    if nval % nitem != 0:
        print("Parsing (%s) is NOT succesful!  Row=%d : nitem=%d : loop=%d " % (cate, nrow, nitem, loop))

    return items, values


##########################################################
def get_rows(items, values):
    """items: all the cif items for the table;  values: all the values in a big row.
    return the row of list [[...] ...[...]]
    """

    row, nval, nitem = [], len(values), len(items)
    for i in range(nval):
        if (i + 1) % nitem == 0:
            row.append(values[i + 1 - nitem : i + 1])

    return row


##########################################################
def parse_values(items, values, item):
    m = -1
    for i, x in enumerate(items):
        if x == item:
            m = i
            break
    if m == -1:
        # print("Error: the item (%s) is not in the cif file." %item)
        return []

    nitem = len(items)
    nrow = len(values) // nitem
    val = []
    for i in range(nrow):
        val.append(values[m + i * nitem])

    return val


##########################################################
#  Below are the applications
##########################################################
def cif_nonloop_format(rows, items):
    """get the maximum length in the items"""

    fmt = []  # the writing format
    for x in items:
        n = len(x)
        fmt.append(n)
    return max(fmt)


##########################################################
def cif_loop_format(rows):
    """the columns are formated (the max length is determined for each column)"""

    fmt = []  # the writing format
    nrow, ncol = len(rows), len(rows[0])
    for i in range(ncol):  # put column format in a list
        tmp = []
        for j in range(nrow):
            if " " in rows[j][i] or "'" in rows[j][i] or '"' in rows[j][i]:
                n = len(rows[j][i]) + 2  # add two " or ' in the item
            else:
                n = len(rows[j][i])
            tmp.append(n)
        fmt.append(max(tmp))

    return fmt


##########################################################
def string_between_char(ss, c1, c2):
    """return the string between c1 & c2; ss is the input."""

    return ss[ss.find(c1) + 1 : ss.rfind(c2)]


##########################################################
def string_after_char(k, ss, cc):
    """get one string starting from k and return index and the new string"""

    end = k
    ssnew = ss[k:]
    nss = len(ss)
    if cc.isspace():
        n = 0
        for i in range(k, nss):
            if n > 0 and (ss[i].isspace() or i == nss - 1):
                if ss[i].isspace():
                    ssnew = ss[k:i]
                elif i == nss - 1:  # the last one
                    ssnew = ss[k : i + 1]
                end = i + 1
                break
            elif not ss[i].isspace():
                n = n + 1

    else:
        n = 0
        for i in range(k + 1, nss):
            if n > 0 and ((ss[i] == cc and i < nss - 1 and ss[i + 1] == " ") or i == nss - 1):
                ssnew = ss[k + 1 : i]
                end = i + 1
                break
            elif not ss[i].isspace():
                n = n + 1
    return end - k, ssnew


##########################################################
def lines_between_char(k, lines, cc):
    """k: the starting position; return cc lines and end position"""

    ss, ssnew = "", []
    n = 0
    for i in range(k, len(lines)):
        if n > 0 and lines[i] == cc:
            end = i + 1
            break
        ssnew.append(lines[i])
        n = n + 1
    ss = "\n".join(ssnew)
    if ss[0] == cc:
        ss = ss[1:]

    return end, ss


##########################################################
def clean_str(ss):
    ssnew = ss.split()
    return " ".join(ssnew)


##########################################################
def cif2cif(file):
    """parse cif and rewrite in new cif"""

    if not (os.path.exists(file) and os.path.getsize(file) > 0):
        print("Error: file (%s) is empty or zero size")
        return ""

    blockid = "data_xxxx\n"
    nfile = file + "_new.cif"
    fw = open(nfile, "w")

    flist = open(file, "r").readlines()
    nblock = []
    for i, x in enumerate(flist):
        if "data_" in x.lstrip()[:5] and len(x.strip().split("data_")[1]) > 0:
            nblock.append(i)

    if nblock:  # has data block name
        nblock.append(len(flist))
        for i in range(1, len(nblock)):
            alist = flist[nblock[i - 1] : nblock[i]]
            blockid = flist[nblock[i - 1]]
            cif2cif_single_block(fw, blockid, alist)
    else:  # no data block name
        alist = flist[0 : len(flist)]
        cif2cif_single_block(fw, blockid, alist)

    return nfile


##########################################################
def cif2cif_single_block(fw, block, flist):
    """parse cif and rewrite in new cif. block: the line with block id.; flist: a list"""

    fw.write(block)
    cif = cif_table_items(flist)
    for x in cif:
        items, values = cifparse(flist, x[1])
        rows = get_rows(items, values)

        n = len(values) / len(items)
        if n == 1:  # non-looped
            write_cif_non_loop(fw, items, values, rows)
        elif n > 1:  # looped
            write_cif_loop(fw, items, rows)


##########################################################
def write_cif_non_loop(fw, items, values, rows):
    fmt = cif_nonloop_format(rows, items)
    fw.write("#\n")
    for i, p in enumerate(items):
        v = values[i].strip()
        if not v:
            v = "?"
        if "\n" in v:
            fw.write("%s  \n" % p)
            fw.write(";%s\n;\n" % v)
        elif " " in v or '"' in v or "'" in v:
            fw.write(p.ljust(fmt + 1))
            if "'" in v:
                t = '"%s"' % v
            else:
                t = "'%s'" % v
            fw.write("   %s\n" % t)
        else:
            fw.write(p.ljust(fmt + 1))
            fw.write("   %s\n" % v)


##########################################################
def write_cif_loop(fw, items, rows):
    fmt = cif_loop_format(rows)
    fw.write("\n#\nloop_\n")
    for p in items:
        fw.write("%s\n" % p)
    for row in rows:  # (left formated)
        for ii, zz in enumerate(row):
            yy = zz.strip()
            if not yy:
                yy = "?"
            if len(yy) > 60:
                fw.write("\n;%s\n;\n" % yy)
            else:
                if " " in yy or '"' in yy or "'" in yy:  # related to fmt
                    if "'" in yy:
                        t = '"%s"' % yy
                    else:
                        t = "'%s'" % yy
                    fw.write(t.ljust(fmt[ii] + 1))
                else:
                    fw.write(yy.ljust(fmt[ii] + 1))
        fw.write("\n")


##########################################################
def cif_table_items(flist_all):
    """get all the block_category_items from the file. Return a list.
    list=[[block, table, [items..], [block, table, [items..]...]
    """

    flist = []
    for x in flist_all:
        y = x.strip()
        if len(y) > 5 and "_" in y[0] or "loop_" in y[:5] or "data_" in y[:5]:
            flist.append(y)

    #    print('file= ', file, len(flist), flist)

    cif = []
    nlen = len(flist)
    st = 0
    block = "X"
    while st < nlen:
        x = flist[st].strip()
        if len(x) > 5 and "data_" in x[:5]:
            block = x[5:].split()[0]

        n, cate, item = check_item(x)
        if n > 0:  # potentional cif item
            cates, items = "", []
            items.append(item)
            if len(items) == 1:
                cates = cate
            if st == nlen - 1:
                cif.append([block, cates, items])

            for j in range(st + 1, nlen):
                y = flist[j].strip()
                if len(y) > 5 and "data_" in y[:5]:
                    block = y[5:].split()[0]
                m, cate, item = check_item(y)

                if len(items) > 0 and ("loop_" in y[:5] or (m > 0 and cate != cates)):
                    st = j - 1
                    cif.append([block, cates, items])
                    cates, items = "", []
                    break
                elif (m > 0 and cate == cates) and j == nlen - 1:
                    st = j
                    items.append(item)
                    cif.append([block, cates, items])
                    cates, items = "", []
                    break
                if m > 0:
                    items.append(item)
        st = st + 1
    return cif


##########################################################
def check_item(y):
    """y: a one line string: check if it is a possible cif item (n>0)"""

    n = 0
    ss = y.split()[0]
    nlen = len(ss)
    cate, item = "", ""
    if nlen > 4 and ss[0] == "_" and "." in ss[1:]:
        m = ss.find(".") + 1
        cate = ss[:m]
        item = ss[m:]
        if m > 1 and "." not in ss[m + 1 :] and len(cate) > 2 and len(item) > 0:
            n = 1

    return n, cate, item


##########################################################


def get_cell(flist):
    """flist is a list of the cif file
    return cell as a list of float!
    """

    cell = [0, 0, 0, 0, 0, 0]
    items, values = cifparse(flist, "_cell.")
    if not len(items):
        return cell
    a = parse_values(items, values, "_cell.length_a")
    b = parse_values(items, values, "_cell.length_b")
    c = parse_values(items, values, "_cell.length_c")
    alpha = parse_values(items, values, "_cell.angle_alpha")
    beta = parse_values(items, values, "_cell.angle_beta")
    gamma = parse_values(items, values, "_cell.angle_gamma")

    if not (a and b and c and alpha and beta and gamma):
        print("Warning: cells not extracted. Check ciftokens")

    #    print(a , b , c , alpha , beta , gamma)
    for i, x in enumerate([a, b, c, alpha, beta, gamma]):
        if len(x) == 0 or not is_a_number(x[0]):
            print("Error: cell has wrong (%s) values" % x)
            continue
        cell[i] = float(x[0].strip())

    return cell


##########################################################
def get_symm(flist):
    """get symmetry"""

    spg = ""
    items, values = cifparse(flist, "_symmetry.")
    symm = parse_values(items, values, "_symmetry.space_group_name_H-M")
    if symm:
        spg = symm[0].replace("'", "").replace('"', "").strip()
    #    else:
    #        print('Warning: space group not extracted. Check ciftokens')

    return spg


##########################################################
def asym2chain(asym):
    """If asym has length>2, assign it to 2 letters"""

    ch = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789"

    as2ch = {}
    #    print set(asym)
    uch = set(asym)
    id = 0
    for x in uch:
        if len(x) > 2:
            id = 1
            break
    nch = len(uch)
    if id == 0 or nch > 35 * 35:
        return as2ch

    n = 0
    ch2 = []
    for x in ch:
        for y in ch:
            n = n + 1
            xy = "%s%s" % (x, y)
            ch2.append(xy)
            if n > nch + 1:
                break
        if n > nch + 1:
            break

    for i, x in enumerate(uch):
        as2ch[x] = ch2[i]

    return as2ch


##########################################################
def cif2pdb(ciffile):
    """convert the cif to simple pdb file. (good enough for validation)"""

    pdbf = ciffile + ".PDB"
    fw = open(pdbf, "w")

    flist = open(ciffile, "r").readlines()
    method = add_header(fw, flist)
    add_remark3(fw, flist, method)
    add_twin(fw, flist)
    add_tls(fw, flist)
    add_remark200(fw, flist)
    add_remark285(fw, flist)
    add_ncs(fw, flist)

    c = get_cell(flist)
    s = get_symm(flist)

    ss = "CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f %-10s \n" % (c[0], c[1], c[2], c[3], c[4], c[5], s)
    fw.write(ss)
    add_scale(fw, flist)

    # for anisou records
    unatm, uatom, uasym, useq, ucomp, uins, ualt, u11, u22, u33, u12, u13, u23 = get_uij(flist)

    # for xyz
    group, natm, atom, asym, seq, comp, ins, alt, x, y, z, occ, biso, model, symbol = get_xyz(flist)
    if not (atom and comp and asym and seq and x and y and z and occ and biso and symbol):
        print("Error: Cif (%s) to pdb conversion failed. Missing the key cif tokens on atom_site." % ciffile)
        sys.exit()

    nmodel = 1
    if model and is_a_number(model[-1]) and int(model[-1]) > 1:
        nmodel = int(model[-1])

    asym2ch = asym2chain(asym)

    nline = len(x)
    if not group:
        group = ["ATOM" for i in x]
    if not alt:
        alt = ["." for i in x]
    if not ins:
        ins = ["." for i in x]
    if not natm:
        natm = ["%d" % i for i in range(nline)]

    n = 0
    for i in range(nline):
        nat = int(natm[i])
        nres = int(seq[i])
        if nat > 99999:
            nat = 99999
        if nres > 9999:
            nres = 9999

        if len(symbol[i].strip()) > 1:
            tmp1 = "%s    " % atom[i].strip()
        else:
            tmp1 = " %s    " % atom[i].strip()
        if len(atom[i].strip()) == 4:
            tmp1 = "%s    " % atom[i]

        atomname = tmp1[:4]
        if nmodel > 1:
            if int(model[i]) >= 180:
                break  # if too many atom, refmac breaks.
            if i == 0:
                fw.write("MODEL %8s \n" % model[i])
            elif 0 < i < nline - 1 and model[i - 1] != model[i]:
                fw.write("ENDMDL  \n")
                fw.write("MODEL %8s \n" % model[i])

        inst = ins[i]
        if ins[i] == "." or ins[i] == "?":
            inst = " "

        alter = alt[i]
        if alt[i] == "." or alt[i] == "?":
            alter = " "
        v = [float(xx) for xx in (x[i], y[i], z[i], occ[i], biso[i])]
        asym2 = asym[i]
        if asym2ch:
            asym2 = asym2ch[asym[i]]
        space = " "
        if v[4] < 1000:
            ss = "%-6s%5d %4s%1s%3s%2s%4d%1s   %8.3f%8.3f%8.3f%6.3f%6.2f%6s%-4s%2s  \n" % (
                group[i],
                nat,
                atomname,
                alter,
                comp[i],
                asym2,
                nres,
                inst,
                v[0],
                v[1],
                v[2],
                v[3],
                v[4],
                space,
                asym[i],
                symbol[i],
            )
        else:
            ss = "%-6s%5d %4s%1s%3s%2s%4d%1s   %8.3f%8.3f%8.3f%6.3f%6.1f%6s%-4s%2s  \n" % (
                group[i],
                nat,
                atomname,
                alter,
                comp[i],
                asym2,
                nres,
                inst,
                v[0],
                v[1],
                v[2],
                v[3],
                v[4],
                space,
                asym[i],
                symbol[i],
            )

        fw.write(ss)

        if n == len(uatom):
            continue
        if uatom and ucomp and uasym and useq:
            #            print(n, ucomp[n], comp[i], uatom[n], atom[i],uasym[n], asym[i])
            if n < len(ucomp) and ucomp[n] == comp[i] and uatom[n] == atom[i] and uasym[n] == asym[i] and useq[n] == seq[i]:
                u = [10000 * float(xx) for xx in (u11[n], u22[n], u33[n], u12[n], u13[n], u23[n])]
                # uu='%7.0f%7.0f%7.0f%7.0f%7.0f%7.0f' %(u[0], u[1],u[2],u[3],u[4],u[5])
                ss = "ANISOU%5d %4s%1s%3s%2s%4d%1s %7.0f%7.0f%7.0f%7.0f%7.0f%7.0f%2s%-4s%2s  \n" % (
                    nat,
                    atomname,
                    alter,
                    comp[i],
                    asym2,
                    nres,
                    inst,
                    u[0],
                    u[1],
                    u[2],
                    u[3],
                    u[4],
                    u[5],
                    space,
                    asym[i],
                    symbol[i],
                )
                fw.write(ss)
                n = n + 1

    if nmodel > 1:
        fw.write("ENDMDL  \n")
    fw.write("END\n")
    fw.close()

    return pdbf


##########################################################
def get_xyz(flist):
    items, values = cifparse(flist, "_atom_site.")  # a loop
    group = parse_values(items, values, "_atom_site.group_PDB")
    natm = parse_values(items, values, "_atom_site.id")
    symbol = parse_values(items, values, "_atom_site.type_symbol")
    if not symbol:
        symbol = parse_values(items, values, "_atom_site.atom_type_symbol")
    atom = parse_values(items, values, "_atom_site.label_atom_id")
    asym = parse_values(items, values, "_atom_site.auth_asym_id")
    comp = parse_values(items, values, "_atom_site.label_comp_id")
    seq = parse_values(items, values, "_atom_site.auth_seq_id")
    biso = parse_values(items, values, "_atom_site.B_iso_or_equiv")
    occ = parse_values(items, values, "_atom_site.occupancy")
    x = parse_values(items, values, "_atom_site.Cartn_x")
    y = parse_values(items, values, "_atom_site.Cartn_y")
    z = parse_values(items, values, "_atom_site.Cartn_z")
    ins = parse_values(items, values, "_atom_site.pdbx_PDB_ins_code")
    if not ins:
        ins = parse_values(items, values, "_atom_site.ndb_ins_code")
    alt = parse_values(items, values, "_atom_site.label_alt_id")
    model = parse_values(items, values, "_atom_site.pdbx_PDB_model_num")
    if not model:
        model = parse_values(items, values, "_atom_site.pdbx_model")
    if not model:
        model = parse_values(items, values, "_atom_site.ndb_model")

    return group, natm, atom, asym, seq, comp, ins, alt, x, y, z, occ, biso, model, symbol


##########################################################
def get_uij(flist):
    """pars the anisou records"""

    items, values = cifparse(flist, "_atom_site_anisotrop.")  # a loop

    unatm = parse_values(items, values, "_atom_site_anisotrop.id")

    uatom = parse_values(items, values, "_atom_site_anisotrop.pdbx_auth_atom_id")
    if not uatom:
        uatom = parse_values(items, values, "_atom_site_anisotrop.ndb_label_atom_id")
        if not uatom:
            uatom = parse_values(items, values, "_atom_site_anisotrop.ndb_PDB_atom_name")

    uasym = parse_values(items, values, "_atom_site_anisotrop.pdbx_auth_asym_id")
    if not uasym:
        uasym = parse_values(items, values, "_atom_site_anisotrop.ndb_PDB_strand_id")
        if not uasym:
            uasym = parse_values(items, values, "_atom_site_anisotrop.ndb_auth_asym_id")

    useq = parse_values(items, values, "_atom_site_anisotrop.pdbx_auth_seq_id")
    if not useq:
        useq = parse_values(items, values, "_atom_site_anisotrop.ndb_auth_seq_id")
        if not useq:
            useq = parse_values(items, values, "_atom_site_anisotrop.ndb_PDB_residue_no")

    ucomp = parse_values(items, values, "_atom_site_anisotrop.pdbx_auth_comp_id")
    if not ucomp:
        ucomp = parse_values(items, values, "_atom_site_anisotrop.ndb_label_comp_id")
        if not ucomp:
            ucomp = parse_values(items, values, "_atom_site_anisotrop.ndb_PDB_residue_name")

    uins = parse_values(items, values, "_atom_site_anisotrop.pdbx_PDB_ins_code")
    if not uins:
        uins = parse_values(items, values, "_atom_site_anisotrop.ndb_label_ins_code")

    ualt = parse_values(items, values, "_atom_site_anisotrop.pdbx_label_alt_id")
    if not ualt:
        ualt = parse_values(items, values, "_atom_site_anisotrop.ndb_label_alt_id")

    u11 = parse_values(items, values, "_atom_site_anisotrop.U[1][1]")
    u22 = parse_values(items, values, "_atom_site_anisotrop.U[2][2]")
    u33 = parse_values(items, values, "_atom_site_anisotrop.U[3][3]")
    u12 = parse_values(items, values, "_atom_site_anisotrop.U[1][2]")
    u13 = parse_values(items, values, "_atom_site_anisotrop.U[1][3]")
    u23 = parse_values(items, values, "_atom_site_anisotrop.U[2][3]")

    return unatm, uatom, uasym, useq, ucomp, uins, ualt, u11, u22, u33, u12, u13, u23


##########################################################
def add_twin(fw, flist):
    """Add twin records to PDB, if exist"""

    items, values = cifparse(flist, "_pdbx_reflns_twin.")  # a loop

    if not items:
        return

    domain_id = parse_values(items, values, "_pdbx_reflns_twin.domain_id")
    operator = parse_values(items, values, "_pdbx_reflns_twin.operator")
    fraction = parse_values(items, values, "_pdbx_reflns_twin.fraction")

    ndomain = len(domain_id)
    fw.write("REMARK   3  TWIN DETAILS \n")
    fw.write("REMARK   3   NUMBER OF TWIN DOMAINS  :   %s \n" % ndomain)

    for i in range(ndomain):
        if domain_id:
            fw.write("REMARK   3      TWIN DOMAIN   :  %s\n" % (domain_id[i]))
        if operator:
            fw.write("REMARK   3      TWIN OPERATOR :  %s\n" % (operator[i]))
        if fraction:
            fw.write("REMARK   3      TWIN FRACTION :  %s\n" % (fraction[i]))

    fw.write("REMARK   3    \n")


##########################################################
def get_prog(flist):
    """get program and version"""

    prog, vers = "?", "?"

    items, values = cifparse(flist, "_software.")
    p = parse_values(items, values, "_software.name")
    d = parse_values(items, values, "_software.classification")
    v = parse_values(items, values, "_software.version")
    if p and d:
        for i, x in enumerate(d):
            if "refinement" in x:
                prog = p[i].strip()
                vers = v[i].strip().replace(" ", "")
                if "PHENIX" in prog.upper() and "REFINE" not in prog.upper():
                    vers = "(phenix.refine: %s)" % vers
                break

    else:
        items, values = cifparse(flist, "_computing.")
        p = parse_values(items, values, "_computing.structure_refinement")
        if p:
            prog = p[0]

    return prog, vers


##########################################################
def add_tls(fw, flist):
    """Add TLS records to PDB, if exist"""

    prog, version = get_prog(flist)

    items, values = cifparse(flist, "_ccp4_refine_tls.")  # a loop
    if items:
        tab = "_ccp4_refine_tls"
    else:
        items, values = cifparse(flist, "_pdbx_refine_tls.")
        tab = "_pdbx_refine_tls"

    refid = parse_values(items, values, "%s.pdbx_refine_id" % tab)
    tlsid = parse_values(items, values, "%s.id" % tab)
    method = parse_values(items, values, "%s.method" % tab)
    xo = parse_values(items, values, "%s.origin_x" % tab)
    yo = parse_values(items, values, "%s.origin_y" % tab)
    zo = parse_values(items, values, "%s.origin_z" % tab)
    t11 = parse_values(items, values, "%s.T[1][1]" % tab)
    t22 = parse_values(items, values, "%s.T[2][2]" % tab)
    t33 = parse_values(items, values, "%s.T[3][3]" % tab)
    t12 = parse_values(items, values, "%s.T[1][2]" % tab)
    t13 = parse_values(items, values, "%s.T[1][3]" % tab)
    t23 = parse_values(items, values, "%s.T[2][3]" % tab)
    l11 = parse_values(items, values, "%s.L[1][1]" % tab)
    l22 = parse_values(items, values, "%s.L[2][2]" % tab)
    l33 = parse_values(items, values, "%s.L[3][3]" % tab)
    l12 = parse_values(items, values, "%s.L[1][2]" % tab)
    l13 = parse_values(items, values, "%s.L[1][3]" % tab)
    l23 = parse_values(items, values, "%s.L[2][3]" % tab)
    s11 = parse_values(items, values, "%s.S[1][1]" % tab)
    s12 = parse_values(items, values, "%s.S[1][2]" % tab)
    s13 = parse_values(items, values, "%s.S[1][3]" % tab)
    s21 = parse_values(items, values, "%s.S[2][1]" % tab)
    s22 = parse_values(items, values, "%s.S[2][2]" % tab)
    s23 = parse_values(items, values, "%s.S[2][3]" % tab)
    s31 = parse_values(items, values, "%s.S[3][1]" % tab)
    s32 = parse_values(items, values, "%s.S[3][2]" % tab)
    s33 = parse_values(items, values, "%s.S[3][3]" % tab)

    items, values = cifparse(flist, "%s_group." % tab)  # a loop
    refid1 = parse_values(items, values, "%s_group.pdbx_refine_id" % tab)
    groupid = parse_values(items, values, "%s_group.id" % tab)
    tlsid1 = parse_values(items, values, "%s_group.refine_tls_id" % tab)
    asymid1 = parse_values(items, values, "%s_group.beg_auth_asym_id" % tab)
    seqid1 = parse_values(items, values, "%s_group.beg_auth_seq_id" % tab)
    asymid2 = parse_values(items, values, "%s_group.end_auth_asym_id" % tab)
    seqid2 = parse_values(items, values, "%s_group.end_auth_seq_id" % tab)
    detail = parse_values(items, values, "%s_group.selection_details" % tab)

    if not (tlsid and tlsid1):
        return

    bpart1 = 0
    for x in flist:
        if "ATOM RECORD CONTAINS RESIDUAL B FACTORS ONLY" in x:
            bpart1 = 1
            break

    if not (is_a_number(tlsid[-1]) and is_a_number(tlsid1[-1]) and tlsid[-1] == tlsid1[-1]):
        print("Error: Number of TLS in %s & %s_group differ." % (tab, tab))
        return

    fw.write("REMARK   3  TLS DETAILS. \n")
    fw.write("REMARK   3   NUMBER OF TLS GROUPS: %s \n" % tlsid[-1])
    if bpart1:
        fw.write("REMARK   3   ATOM RECORD CONTAINS RESIDUAL B FACTORS ONLY\n")

    if "REFMAC" in prog:
        fw.write("REMARK   3  \n")
    else:
        fw.write("REMARK   3   ORIGIN: CENTER OF MASS \n")

    nlen = len(tlsid1)
    n = 0
    tmp = []
    for i in range(len(tlsid1)):
        if n >= len(tlsid):
            print("Error: Number of TLS group in _pdbx_refine_tls not consistent with _pdbx_refine_tls_group.")
            break

        if "REFMAC" not in prog and i < nlen - 1:
            if tlsid1[i] == tlsid1[i + 1]:
                continue
        elif "REFMAC" in prog and i < nlen - 1:
            if tlsid1[i] == tlsid1[i + 1]:
                tmp.append([asymid1[i], seqid1[i], asymid2[i], seqid2[i]])
                continue

        fw.write("REMARK   3   TLS GROUP : %2d \n" % (n + 1))
        ntmp = len(tmp) + 1

        if "REFMAC" in prog:
            fw.write("REMARK   3    NUMBER OF COMPONENTS GROUP : %d\n" % ntmp)
            fw.write("REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI \n")
            tmp.append([asymid1[i], seqid1[i], asymid2[i], seqid2[i]])
            for y in tmp:
                fw.write("REMARK   3    RESIDUE RANGE :   %s %6s    %s %6s \n" % (y[0], y[1], y[2], y[3]))
            tmp = []
        elif detail:
            details = detail[i].replace("\n", " ")
            write_select_phenix_buster(fw, details, prog)
        elif "PHENIX" in prog or "BUSTER" in prog:
            print("Error: problem in atom selection for the TLS group (%d)" % n)
        else:
            print("Warning: structure refined by program (%s), but TLS not supported by DCC." % prog)

        if is_a_number(xo[n]) and is_a_number(yo[n]) and is_a_number(zo[n]):
            if "BUSTER" in prog:
                orig = "%9.4f%10.4f%10.4f " % (float(xo[n]), float(yo[n]), float(zo[n]))

            else:
                orig = "%8.3f %8.3f %8.3f " % (float(xo[n]), float(yo[n]), float(zo[n]))
        else:
            orig = "   NULL     NULL     NULL "

        fw.write("REMARK   3    ORIGIN FOR THE GROUP (A): %s \n" % orig)
        fw.write("REMARK   3    T TENSOR \n")
        fw.write("REMARK   3      T11: %8s T22: %8s \n" % (t11[n], t22[n]))
        fw.write("REMARK   3      T33: %8s T12: %8s \n" % (t33[n], t12[n]))
        fw.write("REMARK   3      T13: %8s T23: %8s \n" % (t13[n], t23[n]))

        fw.write("REMARK   3    L TENSOR \n")
        fw.write("REMARK   3      L11: %8s L22: %8s \n" % (l11[n], l22[n]))
        fw.write("REMARK   3      L33: %8s L12: %8s \n" % (l33[n], l12[n]))
        fw.write("REMARK   3      L13: %8s L23: %8s \n" % (l13[n], l23[n]))

        fw.write("REMARK   3    S TENSOR \n")
        fw.write("REMARK   3      S11: %8s S12: %8s S13: %8s \n" % (s11[n], s12[n], s13[n]))
        fw.write("REMARK   3      S21: %8s S22: %8s S23: %8s \n" % (s21[n], s22[n], s23[n]))
        fw.write("REMARK   3      S31: %8s S32: %8s S33: %8s \n" % (s31[n], s32[n], s33[n]))
        n = n + 1

    fw.write("REMARK   3    \n")
    fw.write("REMARK   3  BULK SOLVENT MODELLING. \n")
    # print(items)


##########################################################
def write_select_phenix_buster(fw, detail, prog):
    select = "REMARK   3    SELECTION:"
    select1 = "REMARK   3             :"
    if "BUSTER" in prog:
        select = "REMARK   3    SET :"
        select1 = "REMARK   3        :"

    if len(detail) > 55:
        t = detail.split()
        s, ss = "", ""
        for j, x in enumerate(t):
            if j < len(t) - 1 and len(s) >= 46:
                ss = ss + s + "\n%s " % select1
                s = ""
            elif j == len(t) - 1:
                ss = ss + s + x
                break
            s = s + x + " "
            # print(s, ss)

        fw.write("%s %s\n" % (select, ss))
    else:
        fw.write("%s %s\n" % (select, detail))


##########################################################
def add_header(fw, flist):
    """Add simple stuff for the head"""

    items, values = cifparse(flist, "_database_2.")
    id = parse_values(items, values, "_database_2.database_id")
    code = parse_values(items, values, "_database_2.database_code")
    pdbid = ["XXXX"]
    if id and code:
        for i, x in enumerate(code):
            if id[i] == "PDB":
                pdbid = [x.strip()]
                break

    name = parse_values(items, values, "_struct_keywords.pdbx_keywords")
    if not name:
        name = [""]

    items, values = cifparse(flist, "_pdbx_database_related.")
    ids = parse_values(items, values, "_pdbx_database_related.db_id")
    type = parse_values(items, values, "_pdbx_database_related.content_type")
    pdbids = ""
    if ids and type:
        for i, x in enumerate(type):
            if x == "split":
                pdbids = pdbids + ids[i] + " "

    fw.write("HEADER    %-39s %16s    \n" % (name[0].upper(), pdbid[0].upper()))
    if len(pdbids) > 2:
        fw.write("SPLIT      %s  \n" % pdbids)  # legacy

    items, values = cifparse(flist, "_exptl.")
    method = parse_values(items, values, "_exptl.method")
    ss_method = ""
    if method:
        for x in method:
            ss_method = ss_method + x.strip() + "; "
    if ss_method.strip():
        ss = "EXPDTA    %-s\n" % ss_method
        fw.write(ss)
    return ss_method


##########################################################
def add_remark3(fw, flist, method):
    """Add simple stuff for REMARK3"""

    items, values = cifparse(flist, "_refine.")
    refid = parse_values(items, values, "_refine.pdbx_refine_id")

    nobs = parse_values(items, values, "_refine.ls_number_reflns_obs")
    if not nobs or "?" in nobs:
        nobs = ["NULL"]
    resh = parse_values(items, values, "_refine.ls_d_res_high")
    if not resh or "?" in resh:
        resh = ["NULL"]
    resl = parse_values(items, values, "_refine.ls_d_res_low")
    if not resl or "?" in resl:
        resl = ["NULL"]
    comp = parse_values(items, values, "_refine.ls_percent_reflns_obs")
    if not comp or "?" in comp:
        comp = ["NULL"]
    robs = parse_values(items, values, "_refine.ls_R_factor_obs")
    if not robs or "?" in robs:
        robs = ["NULL"]
    rwork = parse_values(items, values, "_refine.ls_R_factor_R_work")
    if not rwork or "?" in rwork:
        rwork = ["NULL"]
    if rwork[0] == "NULL" and robs[0] != "NULL":
        rwork = robs
    rfree = parse_values(items, values, "_refine.ls_R_factor_R_free")
    if not rfree or "?" in rfree:
        rfree = ["NULL"]
    compfr = parse_values(items, values, "_refine.ls_percent_reflns_R_free")
    if not compfr or "?" in compfr:
        compfr = ["NULL"]
    nfree = parse_values(items, values, "_refine.ls_number_reflns_R_free")
    if not nfree or "?" in nfree:
        nfree = ["NULL"]

    items, values = cifparse(flist, "_database_PDB_remark.")
    text_id = parse_values(items, values, "_database_PDB_remark.id")
    text = parse_values(items, values, "_database_PDB_remark.text")
    bres = ""
    for i, x in enumerate(text_id):
        if x == "3":
            if "U VALUES      : RESIDUAL ONLY" in text[i]:
                bres = "U VALUES      : RESIDUAL ONLY"

    prog, version = get_prog(flist)
    idd = ["REMARK   3 ", "REMARK   3 "]
    if refid and len(refid) == 2 and "X-RAY " in refid[0]:
        idd[0] = "REMARK   3  X-RAY DATA. "
        idd[1] = "REMARK   3  NEUTRON DATA. "
    elif refid and len(refid) == 2 and "NEUTRON " in refid[0]:
        idd[0] = "REMARK   3  NEUTRON DATA. "
        idd[1] = "REMARK   3  X-RAY DATA. "

    wds = """REMARK   2        
REMARK   3                   
REMARK   3 REFINEMENT.                                                          
REMARK   3   PROGRAM     : %s  %s 
REMARK   3                                                                      
""" % (
        prog,
        version,
    )
    fw.write(wds)

    wds = """REMARK   3      
%s
REMARK   3  
REMARK   3  REFINEMENT TARGET : NULL                                            
REMARK   3                                                                      
REMARK   3  DATA USED IN REFINEMENT.                                            
REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : %s    
REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : %s       
REMARK   3   COMPLETENESS (WORKING+TEST)   (%%) : %s        
REMARK   3   NUMBER OF REFLECTIONS             : %s                          
REMARK   3                                            
REMARK   3  FIT TO DATA USED IN REFINEMENT.          
REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT  
REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM   
REMARK   3   R VALUE            (WORKING SET) : %s    
REMARK   3   FREE R VALUE                     : %s    
REMARK   3   FREE R VALUE TEST SET SIZE   (%%) : %s  
REMARK   3   FREE R VALUE TEST SET COUNT      : %s    
REMARK   3  %s  
REMARK   3                                                                      
""" % (
        idd[0],
        resh[0],
        resl[0],
        comp[0],
        nobs[0],
        rwork[0],
        rfree[0],
        compfr[0],
        nfree[0],
        bres,
    )
    fw.write(wds)

    if "X-RAY" in method and "NEUTRON" in method and len(resh) > 1:
        wds = """REMARK   3
%s
REMARK   3  
REMARK   3  REFINEMENT TARGET : NULL                                            
REMARK   3                                                                      
REMARK   3  DATA USED IN REFINEMENT.                                            
REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : %s    
REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : %s       
REMARK   3   COMPLETENESS (WORKING+TEST)   (%%) : %s        
REMARK   3   NUMBER OF REFLECTIONS             : %s                          
REMARK   3                                            
REMARK   3  FIT TO DATA USED IN REFINEMENT.          
REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT  
REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM   
REMARK   3   R VALUE            (WORKING SET) : %s    
REMARK   3   FREE R VALUE                     : %s    
REMARK   3   FREE R VALUE TEST SET SIZE   (%%) : %s  
REMARK   3   FREE R VALUE TEST SET COUNT      : %s    
REMARK   3  %s  
REMARK   3                                                                      
""" % (
            idd[1],
            resh[1],
            resl[1],
            comp[1],
            nobs[1],
            rwork[1],
            rfree[1],
            compfr[1],
            nfree[1],
            bres,
        )

        fw.write(wds)


##########################################################
def add_remark285(fw, flist):
    """ """
    # for 285
    items, values = cifparse(flist, "_pdbx_database_remark.")
    text = parse_values(items, values, "_pdbx_database_remark.text")
    if text and "FOLLOWING TRANSFORMATION MATRIX" in text[0]:
        t = text[0].split("\n")
        for x in t:
            if "REMARK 285" not in x and len(x) < 72:
                fw.write("REMARK 285 %-s\n" % x)

        fw.write("REMARK 285 \n")

    # for 350
    items, values = cifparse(flist, "_pdbx_struct_oper_list.")
    if items:
        btype = parse_values(items, values, "_pdbx_struct_oper_list.type")
        b11 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[1][1]")
        b12 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[1][2]")
        b13 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[1][3]")
        t1 = parse_values(items, values, "_pdbx_struct_oper_list.vector[1]")

        b21 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[2][1]")
        b22 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[2][2]")
        b23 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[2][3]")
        t2 = parse_values(items, values, "_pdbx_struct_oper_list.vector[2]")

        b31 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[3][1]")
        b32 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[3][2]")
        b33 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[3][3]")
        t3 = parse_values(items, values, "_pdbx_struct_oper_list.vector[3]")

        if not (b11 and b12 and b13 and t1 and b21 and b22 and b23 and t2 and b31 and b32 and b33 and t3):
            return
        n = 0
        for i in range(len(b11)):
            if btype and "transform " in btype[i]:
                continue
            n = n + 1
            v = [float(x) for x in (b11[i], b12[i], b13[i], t1[i]) if is_a_number(x)]
            if len(v) == 4:
                fw.write("REMARK 350   BIOMT1%4d %9.6f %9.6f %9.6f %14.5f \n" % (n, v[0], v[1], v[2], v[3]))
            v = [float(x) for x in (b21[i], b22[i], b23[i], t2[i]) if is_a_number(x)]
            if len(v) == 4:
                fw.write("REMARK 350   BIOMT2%4d %9.6f %9.6f %9.6f %14.5f \n" % (n, v[0], v[1], v[2], v[3]))
            v = [float(x) for x in (b31[i], b32[i], b33[i], t3[i]) if is_a_number(x)]
            if len(v) == 4:
                fw.write("REMARK 350   BIOMT3%4d %9.6f %9.6f %9.6f %14.5f \n" % (n, v[0], v[1], v[2], v[3]))


##########################################################
def add_scale(fw, flist):
    """some programs need it."""

    items, values = cifparse(flist, "_atom_sites.")
    if items:
        b11 = parse_values(items, values, "_atom_sites.fract_transf_matrix[1][1]")
        b12 = parse_values(items, values, "_atom_sites.fract_transf_matrix[1][2]")
        b13 = parse_values(items, values, "_atom_sites.fract_transf_matrix[1][3]")
        t1 = parse_values(items, values, "_atom_sites.fract_transf_vector[1]")

        b21 = parse_values(items, values, "_atom_sites.fract_transf_matrix[2][1]")
        b22 = parse_values(items, values, "_atom_sites.fract_transf_matrix[2][2]")
        b23 = parse_values(items, values, "_atom_sites.fract_transf_matrix[2][3]")
        t2 = parse_values(items, values, "_atom_sites.fract_transf_vector[2]")

        b31 = parse_values(items, values, "_atom_sites.fract_transf_matrix[3][1]")
        b32 = parse_values(items, values, "_atom_sites.fract_transf_matrix[3][2]")
        b33 = parse_values(items, values, "_atom_sites.fract_transf_matrix[3][3]")
        t3 = parse_values(items, values, "_atom_sites.fract_transf_vector[3]")

        if not (b11 and b12 and b13 and t1 and b21 and b22 and b23 and t2 and b31 and b32 and b33 and t3):
            return

        for i in range(len(b11)):
            v = [float(x) for x in (b11[i], b12[i], b13[i], t1[i]) if is_a_number(x)]
            if not v:
                continue
            fw.write("SCALE1     %9.6f %9.6f %9.6f %14.5f \n" % (v[0], v[1], v[2], v[3]))
            v = [float(x) for x in (b21[i], b22[i], b23[i], t2[i])]
            fw.write("SCALE2     %9.6f %9.6f %9.6f %14.5f \n" % (v[0], v[1], v[2], v[3]))
            v = [float(x) for x in (b31[i], b32[i], b33[i], t3[i])]
            fw.write("SCALE3     %9.6f %9.6f %9.6f %14.5f \n" % (v[0], v[1], v[2], v[3]))


##########################################################
def add_ncs(fw, flist):
    """mainly for virus"""

    items, values = cifparse(flist, "_struct_ncs_oper.")
    if items:
        id = parse_values(items, values, "_struct_ncs_oper.id")
        code = parse_values(items, values, "_struct_ncs_oper.code")
        b11 = parse_values(items, values, "_struct_ncs_oper.matrix[1][1]")
        b12 = parse_values(items, values, "_struct_ncs_oper.matrix[1][2]")
        b13 = parse_values(items, values, "_struct_ncs_oper.matrix[1][3]")

        b21 = parse_values(items, values, "_struct_ncs_oper.matrix[2][1]")
        b22 = parse_values(items, values, "_struct_ncs_oper.matrix[2][2]")
        b23 = parse_values(items, values, "_struct_ncs_oper.matrix[2][3]")

        b31 = parse_values(items, values, "_struct_ncs_oper.matrix[3][1]")
        b32 = parse_values(items, values, "_struct_ncs_oper.matrix[3][2]")
        b33 = parse_values(items, values, "_struct_ncs_oper.matrix[3][3]")

        t1 = parse_values(items, values, "_struct_ncs_oper.vector[1]")
        t2 = parse_values(items, values, "_struct_ncs_oper.vector[2]")
        t3 = parse_values(items, values, "_struct_ncs_oper.vector[3]")

        if not (id and code and b11 and b12 and b13 and t1 and b21 and b22 and b23 and t2 and b31 and b32 and b33 and t3):
            return
        for i in range(len(b11)):
            idd = " "
            if code[i] == "given":
                idd = "1"
            v = [float(x) for x in (b11[i], b12[i], b13[i], t1[i]) if is_a_number(x)]
            if len(v) != 4:
                print("Error: matrix can not be parsed or not a numerical value.")
                continue
            fw.write("MTRIX1%4d %9.6f %9.6f %9.6f %14.5f %4s \n" % (i + 1, v[0], v[1], v[2], v[3], idd))
            v = [float(x) for x in (b21[i], b22[i], b23[i], t2[i])]
            fw.write("MTRIX2%4d %9.6f %9.6f %9.6f %14.5f %4s \n" % (i + 1, v[0], v[1], v[2], v[3], idd))
            v = [float(x) for x in (b31[i], b32[i], b33[i], t3[i])]
            fw.write("MTRIX3%4d %9.6f %9.6f %9.6f %14.5f %4s \n" % (i + 1, v[0], v[1], v[2], v[3], idd))


##########################################################
def add_remark200(fw, flist):
    """only for source"""

    items, values = cifparse(flist, "_diffrn_source.")
    syc = parse_values(items, values, "_diffrn_source.pdbx_synchrotron_y_n")
    if not syc:
        syc = parse_values(items, values, "_diffrn_source.ndb_synchrotron_y_n")
    if not syc:
        syc = ["NULL"]
    site = parse_values(items, values, "_diffrn_source.pdbx_synchrotron_site")
    if not site:
        site = parse_values(items, values, "_diffrn_source.ndb_synchrotron_site")
    if not site:
        site = ["NULL"]
    beam = parse_values(items, values, "_diffrn_source.pdbx_synchrotron_beamline")
    if not beam:
        beam = parse_values(items, values, "_diffrn_source.pdbx_synchrotron_beamline")
    if not beam:
        beam = ["NULL"]
    wave = parse_values(items, values, "_diffrn_source.pdbx_wavelength_list")
    if not wave:
        wave = parse_values(items, values, "_diffrn_source.rcsb_wavelength_list")
    if not wave:
        wave = ["NULL"]

    wds = """REMARK 200                                                                      
REMARK 200  SYNCHROTRON              (Y/N) : %s                                  
REMARK 200  RADIATION SOURCE               : %s                               
REMARK 200  BEAMLINE                       : %s                             
REMARK 200  X-RAY GENERATOR MODEL          : NULL                               
REMARK 200  MONOCHROMATIC OR LAUE    (M/L) : M                                  
REMARK 200  WAVELENGTH OR RANGE        (A) : %s                           
REMARK 200                                                                      
""" % (
        syc[0],
        site[0].strip(),
        beam[0].strip(),
        wave[0].strip(),
    )

    fw.write(wds)


##########################################################
def is_a_number(s):
    try:
        float(s)
        return True
    except ValueError:
        # print('Error: %s is not a number.' %s)
        return False


##########################################################
def usage():
    content = """
    
-------------------------------------------------------------------  
 A program for fast parsing mmcif files. (HY,2013-01-10)
-------------------------------------------------------------------
      How to use the cif parser (see the main below)
      

-------------------------------------------------------------------
              Defination of valid mmcif format:
1. for non-looped items, below are the valid syntaxes
_citation.title "The Structure and Evolution of the Major Capsid Protein"
or 
_citation.title 'The Structure and Evolution of the Major Capsid Protein'
or 
_citation.title
;The Structure and Evolution of
the Major Capsid Protein
;

2. for looped items, below are the valid syntaxes
loop_
_citation_author.citation_id 
_citation_author.name 
_citation_author.ordinal 
primary '???, N.' 1 
primary '???, A.'     2 

or
primary
"???, N." 1 
primary
'???, A.'     2 

or
primary
;???, N.
;
1 
primary
'???, A.'     2 

-------------------------------------------------------------------
"""
    print(content)
    sys.exit()

    """
##########################################################
if __name__ == '__main__':
    
    files=sys.argv[1]
    print('file=%s' %files)
    flist=open(files, 'r').readlines()
    items,values = cifparse(flist, '_cell.')
    v1=parse_values(items,values, '_cell.entry_id')
    v2=parse_values(items,values, '_cell.length_c')
    v3=parse_values(items,values, '_cell.angle_beta')
    print(v1,v2,v3)

    items,values = cifparse(flist, '_database_PDB_rev.')
    v1=parse_values(items,values, '_database_PDB_rev.num')
    v2=parse_values(items,values, '_database_PDB_rev.date')
    v3=parse_values(items,values, '_database_PDB_rev.status')

    print(v1,v2,v3)

##########################################################
    
    """
