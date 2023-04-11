# ==================================================
#   This module is sf related
# ==================================================

import os
import sys
import shutil
import util
import prog
import parse
import cifparse as cif


def sf_quality(sffile, pdbfile):
    """run ctruncate and analysis the result"""

    sf_tmp = sffile + ".mtz"

    sf_tmp1 = prog.sf_convertor(pdbfile, sffile, "mmcif", "")
    sf_tmp = prog.sf_convertor(pdbfile, "SF_4_validate.cif", "mtz", "")
    if not util.check_file(500, sf_tmp):
        print("Error: Conversion failed! Maybe your sf file has no symmetry/cell.")
        sys.exit()

    if pdbfile == "-xtriage":
        outfile = prog.run_phenix(pdbfile, sf_tmp, "xtriage")

    elif pdbfile == "-ctruncate":  # do ctruncate
        outfile = prog.run_ccp4(pdbfile, sf_tmp, "ctruncate")

    else:  # do ctruncate
        outfile = prog.run_ccp4(pdbfile, sf_tmp, "pointless")

    if not util.check_file(500, outfile):
        print("Error: The log file (%s) is not generated" % outfile)
        sys.exit()

    fp = open(outfile, "r")
    flist = fp.readlines()
    for i, x in enumerate(flist):  # for xtriage/ctruncate
        t = x.strip()
        if "<|L|>:" in t[:6] or "<L^2>:" in t[:8] or "L-test suggests" in t[:30] or "twinning, twin fraction" in t:
            print(t)
        elif t == "Acentric reflections" or t == "Centric reflections" or t == "L test for acentric data":
            print(t)
        elif (
            ("<I^2>/<I>^2" in x or "<F>^2/<F^2>" in x or " <|E^2 - 1|> " in x or "Mean |L|   :" in x or "Mean  L^2  :" in x or "L statistic = " in x)
            and "untwinned" in x
            and " perfect twin" in x
        ):
            print(t)
        elif t == "$TABLE: Acentric moments of E:" or t == "$TABLE: Centric moments of E:":
            print(x.split(":")[1])
        elif "$GRAPHS: 4th moment of E (Expected " in x or ": 1st & 3rd moments of E (Expected " in x:
            print(x.split(":")[1])

        elif "ML estimate of overall B value of None:" in x:
            print("The ML estimate of overall B value = %s" % flist[i + 1])
        elif "[MaxAnisoB-MinAnisoB]/[MaxAnisoB" in x:
            print("%s" % x)
        elif "pseudo)merohedral twin laws" in x:
            print("%s" % x)

        elif "No significant pseudotranslation is detected" in x:
            print("No significant pseudotranslation is detected.\n")

    fp.seek(0)
    law = ""
    for x in fp:  # for plots
        if " Acentric_theory " in x and "_perfect_twin " in x:  # Xtriage
            fname, start, end = outfile + "_Ltest.data", "$$", "$$"
            extract_data4plot(fp, fname, start, end, "Ltest")
            plot_items(fname, "twin", "", "", "", "")

        elif "H test for possible twin law" in x:  # Xtriage
            law = x.split("law")[1].strip()
        elif "H " in x[:5] and "Observed_S(H)" in x:  # Xtriage
            fname, start, end = outfile + "_Htest.data", "$$", "$$"
            extract_data4plot(fp, fname, start, end, "Htest")
            title = "H test for Acentric data (twin law =%s)" % law
            plot = """  using 1:2 t "Observed data",  '' u 1:3 t "Untwinned by theory" ,  '' u 1:4 t "Fitted" """
            plot_items(fname, "twin", title, "|H|", "", plot)

        elif "1/resol^2   <I/sigma_I>;_Signal_to_noise   $$" in x:  # Xtriage
            fname, start, end = outfile + "_I_sigI.data", "$$", "$$"
            extract_data4plot(fp, fname, start, end, "Wilson")
            plot_items(fname, "I/sigI", "", "", "", "")

        if " Expected_untwinned " in x and "Expected_twinned" in x:  # Ctruncate
            fname, start, end = outfile + "_Ltest.data", "$$", "$$"
            extract_data4plot(fp, fname, start, end, "Ltest")
            plot_items(fname, "twin", "", "", "", "")

        elif "|L|       N(|L|)  Untwinned    Twinned" in x:  # Pointless
            fname, start, end = outfile + "_Ltest.data", "$$", "$$"
            extract_data4plot(fp, fname, start, end, "Ltest")
            plot_items(fname, "twin", "", "", "", "")

    fp.close()

    print("\nThe output file = %s\n" % outfile)

    util.delete_file("sf_information.txt", "sf_format_guess.text", sf_tmp, sf_tmp1)


##########################################################
def plot_items(plot_data, id, title, xlabel, ylabel, plot):  # pylint: disable=redefined-builtin
    """ """

    if id == "twin":
        if len(title) == 0:
            title = "The Padilla-Yeates cumulative probability (acentric data)"
        if len(xlabel) == 0:
            xlabel = "|L|"
        if len(ylabel) == 0:
            ylabel = "N(|L|)"
        if len(plot) == 0:
            plot = """  using 1:2 t "Observed data",  '' u 1:3 t "Untwinned by theory" ,  '' u 1:4 t "Twinned by theory" """
    elif id == "I/sigI":
        title = "The change of <I/sigI> with resolution (1/resol^2)"
        xlabel, ylabel = "1/resol^2", "<I/sigI>"
        plot = '  using 1:2 t "<I/sigI> vs resolution" '

    #    print(plot_data)
    prog.gnu_plot1(plot_data, xlabel, ylabel, title, plot)
    if not util.check_file(500, plot_data + ".png"):  # png not generate, try ps file
        print("Trying to plot in postscript")
        prog.gnu_plotps(plot_data, xlabel, ylabel, title, plot)


##########################################################
def extract_data4plot(fp, fname, start, end, idd):
    """data is between start and end"""

    fw = open(fname, "w")

    n = 0
    for y in fp:
        n = n + 1
        if start in y and n < 3:
            continue
        if end in y and n > 3:
            break
        fw.write(y)

    fw.close()

    # plug is ugly if using I/sigI vs resolution
    #
    if idd == "Htest":
        fp = open(fname, "r").readlines()
        fw = open(fname, "w")
        for x in fp:
            t = x.strip().split()
            ss = "  ".join([t[0], t[1], t[0], t[2]])
            ss = ss + "\n"
            fw.write(ss)

        fw.close()


##########################################################
def sf_symmetry(sffile, pdbfile):
    """get the best space group by pointless."""

    sfmtz = prog.sf_convertor(pdbfile, sffile, "mtz", "")

    if not util.check_file(500, sfmtz):
        print("Error: MTZ file not generated, check symmetry/cell in sf file.")
        sys.exit()

    print("Getting the best space group by pointless...")
    out = prog.run_ccp4(pdbfile, sfmtz, "pointless")

    print("For details, please see the output file =%s" % out)

    util.delete_file(sfmtz, "sf_format_guess.text", "sf_information.txt")


########################################################################
def test_sf_content(file_inp):
    """test if sf is mtz (id==1) or mmcif with pdbx_r_free_flag (id==2)
    or mmcif with maps (id==3)
    """

    id = 0  # pylint: disable=redefined-builtin
    n = 0
    fp = open(file_inp, "r")
    for x in fp:
        if n == 0 and x.split()[0] == "MTZ":
            id = 1
            break
        elif "_refln.pdbx_r_free_flag" in x:
            id = 2
        elif "_refln.pdbx_FWT" in x:
            id = 3

        n = n + 1
        if n > 1000:
            break

    fp.close()

    return id


########################################################################
def get_file_type(inp1, inp2):
    """determine xyz or sf from inp1&2"""

    xyz_inp, sf_inp = "", ""
    id = test_sf_content(inp1)  # pylint: disable=redefined-builtin
    if id:  # inp1 is SF file (cif or mtz)
        xyz_inp, sf_inp = inp2, inp1
    else:
        id = test_sf_content(inp2)  # inp2 is SF file (cif or mtz)
        if id:
            xyz_inp, sf_inp = inp1, inp2

    return id, xyz_inp, sf_inp


########################################################################
def display_map_coot(inp1, inp2):
    """display map by coordinate & SF file."""

    _mtz, xyz_inp, sf_inp = get_file_type(inp1, inp2)

    if not xyz_inp or not sf_inp:
        print("Error: File must be in MTZ or a mmcif file having _refln.pdbx_FWT\n")
        sys.exit()

    fp = open(sf_inp, "r")
    sf_new = sf_inp + "__new"
    fw = open(sf_new, "w")
    n = 0
    for x in fp:
        if "_refln.pdbx_r_free_flag" in x:  # change it to bypass a bug in CAD.
            fw.write("%s_x\n" % (x.rstrip()))
        else:
            if "_refln.pdbx_FWT" in x:
                n = n + 1
            fw.write(x)
    fw.close()

    add = ""
    if n > 0:  # a cif file with pdbx_FWT
        sf = prog.sf_convertor(xyz_inp, sf_new, "mtz", add)
        arg = "coot --pdb %s --mtz %s" % (xyz_inp, sf)
        os.system(arg)
    else:
        print("Error: SF file should be the mmcif having _refln.pdbx_FWT.")

    util.delete_file(sf_new, sf_new + ".mtz", "sf_information.cif", "sf_format_guess.text")


########################################################################
def get_freer_flag(inp1, inp2):
    """find the correct free set in the mtz or mmcif file by exaustive test."""

    idd, xyz_inp, sf_inp = get_file_type(inp1, inp2)

    if not xyz_inp or not sf_inp:
        print("Error: Stopped searching the free set.  ")
        print("       File must be in MTZ or a mmcif file having _refln.pdbx_r_free_flag\n")
        sys.exit()

    flist = open(sf_inp, "r").readlines()

    out = xyz_inp + "_test_Rf"
    fw = open(out, "w")
    tmp = "set  reso  R_rep  Rf_rep  CC_rep    R_cal  Rf_cal  CC_cal \n"
    print(tmp)
    fw.write(tmp)

    mark = " "
    for i in range(20):
        if idd == 1:  # sf is mtz
            add = " -freer %d" % i
            sf = prog.sf_convertor(xyz_inp, sf_inp, "mmcif", add)
        else:
            sf = cif2cif_sf(flist, i)

        dcc_out = prog.run_dcc(xyz_inp, sf, " -one -no_xtriage ")
        val, _rep, notls, _tls, _bstat, _prog_rep = parse.values_from_dcc(dcc_out)
        if (notls[3] - notls[2]) > 0.02:
            mark = "*"
        arg = " %2d  %s %s\n" % (i, val, mark)
        print(arg)
        fw.write(arg)
        if mark == "*":
            break
        util.delete_file(sf)
    fw.close()

    print("The statistics output file is =%s" % out)
    if mark == "*":
        print("The matched SF file is =%s\n" % sf)
    else:
        print("Warning: No proper free set is found! \n")


##########################################################
def cif2cif_sf(flist, num):
    sf = "SF-%d.cif" % num
    fw = open(sf, "w")

    n1 = 0
    n2 = len(flist)
    for i, x in enumerate(flist):
        if "loop_" in x.lstrip()[:5] and "_refln." in flist[i + 1].lstrip()[:7]:
            n1 = i
            break
    sflist = []
    for i in range(n1, n2):
        sflist.append(flist[i])
        t = flist[i].lstrip()
        if "#" in t[:1] or "data_" in t[:5] or ("_" in t[:1] and "_refln." not in t[:7]):
            n2 = i
            break

    for i in range(0, n1):
        fw.write(flist[i])  # write the front part
    items, values = cif.cifparse(flist[n1:n2], "_refln.")

    status_col = items.index("_refln.status")
    flag_col = items.index("_refln.pdbx_r_free_flag")

    rows = cif.get_rows(items, values)
    for i in range(len(rows)):  #
        if rows[i][status_col] == "f":
            rows[i][status_col] = "o"
    for i in range(len(rows)):
        if int(rows[i][flag_col]) == num and rows[i][status_col] == "o":
            rows[i][status_col] = "f"

    cif.write_cif_loop(fw, items, rows)

    for i in range(n2, len(flist)):
        fw.write(flist[i])  # write the end part

    fw.close()

    return sf


##########################################################
def test_swap_fi(pdb, sf):
    """Only for test_r (swap I->F or F->I)"""

    dcc_out = prog.run_dcc(pdb, sf, "")
    rep, calc1 = parse.val_from_dcc(dcc_out)
    (_swap, sf_new) = swap_sf(sf)
    shutil.move(dcc_out, dcc_out + "_orig")

    dcc_out = prog.run_dcc(pdb, sf_new, "")
    rep, calc2 = parse.val_from_dcc(dcc_out)
    shutil.move(dcc_out, dcc_out + "_swap")

    print("            prog    resh   Rwork  Rfree    Bwil     |L^2|     Z    FOM     FCC    RsR    Dcc ")
    print("Reported: %s" % rep)
    print("Original: %s" % calc1)
    print("Swapped : %s" % calc2)


##########################################################
def swap_sf(sffile):
    """exchange the cif items F-->I or I-->F"""

    sffile_new = sffile + "_swap"
    fr = open(sffile, "r")
    fw = open(sffile_new, "w")

    swap1, swap2 = "?->?", "?->?"
    n = 0
    for line in fr:
        if "_refln.F_meas_au" in line:
            swap1 = "F->I"
            line = "_refln.intensity_meas\n"
        elif "_refln.F_meas_sigma_au" in line:
            line = "_refln.intensity_sigma\n"
        elif "_refln.intensity_meas" in line:
            swap2 = "I->F"
            line = "_refln.F_meas_au\n"
        elif "_refln.intensity_sigma" in line:
            line = "_refln.F_meas_sigma_au\n"
        n = n + 1

        if n > 1000 and "data_r" in line:
            break

        fw.write(line)

    fr.close()
    fw.close()
    return swap1 + " " + swap2, sffile_new
