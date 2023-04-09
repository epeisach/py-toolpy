import os, math
import cifparse as cif


def remove_alt_conformer(files):
    """keep the first conformer and make occ=1.0"""

    nlen = len(files)
    if nlen == 1:
        remove_allalt_conformer(files[0])
    else:
        out = files[0] + "_palt"
        file = files[0]  # give ch-res__conform.
        fp = open(file, "r")
        fw = open(out, "w")

        dic = {}
        for i in range(1, nlen):
            t = files[i].strip().split("_")
            dic[t[0]] = t[1]

        for x in fp:
            if x[16:17] != " ":
                ssk = "%s%s" % (x[20:22].strip(), x[22:26].strip())

                if ssk in dic.keys() and dic[ssk] == x[16:17]:
                    continue
                elif ssk in dic.keys() and dic[ssk] != x[16:17]:
                    if "ATOM" in x[:4] or "HETATM" in x[:6]:
                        t = "".join([x[:16], " ", x[17:56], "1.00", x[60:]])
                    else:
                        t = "".join([x[:16], " ", x[17:]])
                    fw.write(t)
                else:
                    fw.write(x)
            else:
                fw.write(x)

        fw.close()
        fp.close()

        print("The new pdb = %s" % out)


#######################################################
def remove_allalt_conformer(file):
    """remove all the alter conformers and make the OCC be one."""

    out = file + "_noalt"
    fw = open(out, "w")

    flist = open(file, "r").readlines()

    n = len(flist)
    tmp = []
    for i, x in enumerate(flist):
        if ("ATOM" in x[:4] or "HETATM" in x[:6] or "ANISOU" in x[:6]) and x[16:17] != " ":
            for j in range(i + 1, n):
                y = flist[j]
                if x[22:26] != y[22:26]:
                    break
                if x[11:16] == y[11:16] and x[17:26] == y[17:26] and x[16:17] != y[16:17] and y[16:17] != " ":
                    tmp.append(j)

    for i, x in enumerate(flist):
        if ("ATOM" in x[:4] or "HETATM" in x[:6] or "ANISOU" in x[:6]) and x[16:17] != " ":
            if i in tmp:
                continue
            t = x
            if "ATOM" in x[:4] or "HETATM" in x[:6]:
                t = x[:54] + "  1.00" + x[60:]
            t1 = t[:16] + " " + t[17:]
            fw.write(t1)

        else:
            fw.write(x)

    fw.close()
    print("The new pdb = %s" % out)


#######################################################
def remove_bad_water(pdbfile, sffile):
    npdb = pdbfile + "_rmwat"
    fw = open(npdb, "w")

    dccfile = pdbfile + "_dcc.cif"
    arg = "%s/bin/dcc  -pdb %s -sf %s -no_xtriage -o %s  " % (os.environ["DCCPY"], pdbfile, sffile, dccfile)
    os.system(arg)

    print("Removing bad waters from %s" % dccfile)
    flist = open(dccfile, "r").readlines()

    items, values = cif.cifparse(flist, "_pdbx_rscc_mapman_prob.")
    ch = cif.parse_values(items, values, "_pdbx_rscc_mapman_prob.auth_asym_id")
    comp = cif.parse_values(items, values, "_pdbx_rscc_mapman_prob.auth_comp_id")
    seq = cif.parse_values(items, values, "_pdbx_rscc_mapman_prob.auth_seq_id")
    rsr = cif.parse_values(items, values, "_pdbx_rscc_mapman_prob.real_space_R")
    wrsr = cif.parse_values(items, values, "_pdbx_rscc_mapman_prob.RsR_over_correlation")
    dcc = cif.parse_values(items, values, "_pdbx_rscc_mapman_prob.correlation")
    nn = len(ch)
    fp = open(pdbfile, "r")
    for x in fp:
        if ("ATOM" in x[:4] or "HETATM" in x[:6] or "ANISOU" in x[:6]) and x[17:20] == "HOH":
            ch1, seq1 = x[20:22].strip(), x[22:26].strip()

            id = 0
            for i, y in enumerate(comp):
                if "HOH" in y and ch1 == ch[i] and seq1 == seq[i]:
                    id = 1
                    break
            if id == 1:
                print("removing water (%s)" % x[:28])
                continue
            else:
                fw.write(x)
        else:
            fw.write(x)

    fw.close()
    fp.close()
    print("The new pdb file = %s" % npdb)


#######################################################
def add_water(mtz, pdbfile, id):
    map = mtz + "TMP.map"
    wpdb = ""
    if id == 1:  # mFo-dFc map
        mtz2map(mtz, pdbfile, map, "FO_FC")
        wpdb = water_peak(map, pdbfile, rms=3.0)
    else:
        mtz2map(mtz, pdbfile, map, "2FO_FC")
        wpdb = water_peak(map, pdbfile, rms=1.0)

    npdb = pdbfile + "_water"

    fw = open(npdb, "w")
    fp = open(wpdb, "r")
    fp1 = open(pdbfile, "r")

    nr, na = 0, 0
    for x in fp1:
        if "ATOM" in x[:4] or "HETATM" in x[:6]:
            na, nr = int(x[6:11]), int(x[22:26])
        fw.write(x)

    print("%s %s" % (nr, na))
    n = 1
    for x in fp:
        if "ATOM" in x[:4]:
            natom = n + na
            nres = n + nr + 1
            if natom > 99999:
                natom = 99999
            if nres > 9999:
                natom = 9999
            if float(x[54:60]) > 6:
                print("Warning: sigma value (%s ) is too large for water (%s)" % (x[54:60], x[22:26]))
            t = "".join(["HETATM", "%5d" % natom, x[11:20], " W", "%4d" % nres, x[26:54], "  1.00 30.00", x[66:]])
            #            t='HETATM' + '%5d'%natom  + x[11:20] + ' W' + '%4d'%nres + x[26:54]
            fw.write(t)
            n = n + 1
    fw.close()
    fp.close()
    fp1.close()

    print("The pdbfile with HOH = %s" % npdb)


#######################################################
def water_peak(map, pdbfile, rms):
    npdb = pdbfile + "_water_TMP"
    peak = map + ".PEAK"
    print("Adding waters ...")

    symm = os.popen('egrep "^CRYST1" %s -m 1 | cut -c 55-66' % pdbfile).read().strip()

    arg = """
#!/bin/sh

########################################################################
#  get the peaks of the density.
########################################################################

peakmax   MAPIN  %s  XYZOUT  %s << END-peakmax
NUMPEAKS 6000
THRESHOLD RMS %.2f
END-peakmax
#
########################################################################
# run atpeak to check distances within 3.6A
########################################################################
#

watpeak  PEAKS  %s XYZIN  %s  XYZOUT %s <<END-watpeak
TITLE  Toxd contacts
DISTANCE  5 2.4
#BFACTOR   30.0 1.0
SYMMETRY '%s'
!HETATOMONLY ! use if only distances from O and N wanted
END-watpeak
    """ % (
        map,
        peak,
        rms,
        peak,
        pdbfile,
        npdb,
        symm,
    )

    scr_name = "water_peak.sh"
    fw = open(scr_name, "w")
    fw.write(arg)
    fw.close()
    command = "chmod +x %s ; ./%s  >water_peak.log" % (scr_name, scr_name)
    os.system(command)

    return npdb


##########################################################
def mtz2map(mtz, pdbfile, mapout, maptype):
    """run fft to get the map: the mtz must be from refmac at the moment!"""

    if maptype == "2FO_FC":
        f1, phi = "FWT", "PHWT"
    elif maptype == "FO_FC":
        f1, phi = "DELFWT", "PHDELWT"
    elif maptype == "FC":
        f1, phi = "FC", "PHIC"
    elif maptype == "FEM":
        f1, phi = "FEM", "PHIFEM"

    scr = """#!/bin/tcsh  -f 

##################### mtz2map #####################
echo "Calculating %s map using %s %s ..."

fft  HKLIN  %s  MAPOUT  TMPMAP_BY_MTZ.map  <<end >/dev/null 
title mtz2map
 labin F1=%s  PHI=%s
end

# Extend the map to cover the volume of a model (about 4.1A)
mapmask MAPIN  TMPMAP_BY_MTZ.map  XYZIN  %s   MAPOUT %s <<end >/dev/null 
BORDER 5
end

rm -f TMPMAP_BY_MTZ.map

""" % (
        maptype,
        f1,
        phi,
        mtz,
        f1,
        phi,
        pdbfile,
        mapout,
    )

    scr_name = "get_mtz2map.csh"
    fw = open(scr_name, "w")
    fw.write(scr)
    fw.close()
    command = "chmod +x %s ; ./%s  " % (scr_name, scr_name)
    os.system(command)


###########################################
def scale_matrix(pdbfile):
    fp = open(pdbfile, "r").readlines()

    cell, scale = [], 0
    for x in fp:
        if "CRYST1" in x[:6]:
            v = x[6:].split()
            cell = [float(v[i]) for i in range(6)]
        elif "SCALE1" in x[:6]:
            scale = 1
        elif "ATOM" in x[:4] or "HETA" in x[:4]:
            break

    if len(cell) != 6:
        print("No Cell parameters in " + pdbfile)
        return [], [], []

    a, b, c = cell[0], cell[1], cell[2]
    sa = 3.141592654 / 180.0
    alpha, beta, gamma = sa * cell[3], sa * cell[4], sa * cell[5]

    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    cg, sg = math.cos(gamma), math.sin(gamma)

    vol = math.sqrt(1 - ca**2 - cb**2 - cg**2 + 2 * ca * cb * cg)
    sc1 = [sc11, sc12, sc13] = [1 / a, -cg / (a * sg), (ca * cg - cb) / (a * vol * sg)]
    sc2 = [sc21, sc22, sc23] = [0, 1.0 / (b * sg), (cb * cg - ca) / (b * vol * sg)]
    sc3 = [sc31, sc32, sc33] = [0, 0, sg / (c * vol)]

    return sc1, sc2, sc3


###########################################
def clean_file_4d3r(pdb):
    fp = open(pdb, "r").readlines()

    out = pdb + "_d3r"
    fw = open(out, "w")
    scale = 0
    occ, b = 0, 0
    for i, x in enumerate(fp):
        if "ATOM" in fp[i][:4] or "HETATM" in fp[i][:6]:
            occ = occ + float(x[54:60])
            b = b + float(x[60:66]) * float(x[54:60])
    bavg = b / occ
    #    bavg=0
    for i, x in enumerate(fp):
        if i > 0 and ("ATOM" in fp[i - 1][:4] or "HETATM" in fp[i - 1][:6]) and "ANISOU" in x[:6] and (len(x) < 77 or x[77:78] == " "):
            t = x[:70] + "      " + fp[i - 1][76:80] + "\n"
            fw.write(t)

        elif "CRYST1" in x and "SCALE1" not in fp[i + 1]:  # add scale
            m1, m2, m3 = scale_matrix(pdb)
            ss1 = "SCALE1    %10.6f%10.6f%10.6f        0.00000 \n" % (m1[0], m1[1], m1[2])
            ss2 = "SCALE2    %10.6f%10.6f%10.6f        0.00000 \n" % (m2[0], m2[1], m2[2])
            ss3 = "SCALE3    %10.6f%10.6f%10.6f        0.00000 \n" % (m3[0], m3[1], m3[2])

            fw.write(x)
            fw.write(ss1)
            fw.write(ss2)
            fw.write(ss3)
        elif x[17:20] == "UNK" and float(x[60:66]) == 0.0:
            t = "%s%6.2f%s" % (x[:60], bavg, x[66:])
            fw.write(t)
        else:
            fw.write(x)

    fw.close()
    print("The output = %s" % out)
