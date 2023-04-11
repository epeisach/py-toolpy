#!/usr/bin/env /apps/python-2.6.1/bin/python
# #!/usr/bin/env python


import MySQLdb    # pylint: disable=import-error # in apps-centos-5/python-2.6.6/bin/python (if not work, check .tcshrc)
import os
import sys
import math
import util

SQL = "mysql -u rcsbuser -prcsb0000 -h pdb-f-linux-2 cleanv1 -e "
YEAR = 2016  # fitting date upto 2013


######################################################
def proc_stat(args):
    print(args)

    if len(args) < 2:
        usage()

    ddir = initial()

    for _i, x in enumerate(args):
        if "-growth" in x:
            print("processing pdb-growth..")
            pdb_growth_fit(ddir, 1)  # fit growth data, 1 by polynom. 2 by exp.
        elif "-all" in x:
            get_stat()


######################################################
def get_stat():
    """all sort of stats (command: tool -stat >all.log )
    After running the command, mv pdb_each_item.html & DAT_STORE/ to
    /net/wwwdev/auto-check/html/dev/stat

    """

    ddir = initial()

    torsion_info(ddir)
    do_num_stat(ddir)  # stats for remark 3,200
    do_string_stat(ddir)  # all string items

    tls_growth(ddir)  # TLS related statistics
    twin_info(ddir)  # do twin related
    rfact_reso_growth(ddir)  # Rfactor changes with year
    do_reso_entity(ddir)  # resol. with protein,protein-NA, DNA, RNA,
    plot_entity(ddir)  # unique molecular weight/in asu
    rfactor_other(ddir)  # all the Rfactors

    pdb_growth_fit(ddir, 1)  # fit growth data, 1 by polynom. 2 by exp.
    sys.exit()


######################################################
def initial():
    # global SQL  # do not change it.
    # global YEAR
    ddir = "DAT_STORE/"  # create a dir to hold all the data/graphs
    os.system("mkdir %s" % ddir)  # make a fold to hold all the data
    get_html(ddir)  # A html file to link all the stat

    return ddir


######################################################
def usage():
    content = """
    ----------------------------------------------------------------------
    Do various statistics for the pdb items and pdb growth predictions

    1. pdb_stat -growth
    predict growth of deposited and released pdb entries.
    2. pdb_stat  -all
    do all the statistics
    ----------------------------------------------------------------------

"""
    print(content)
    sys.exit()


######################################################
def plot_data():
    ddir = "DATA_PLOT/"  # create a dir to hold all the data/graphs
    os.system("mkdir %s" % ddir)  # make a fold to hold all the data

    fp = open("rna-good.log", "r").readlines()
    data_in = []
    for x in fp:
        t = x.split()
        data_in.append(t)

    # item, xlabel, ylabel,data_type, min, max, steps
    y = ["RNA_residue_pop_Zscore", "RsR_Zscore", "Number of Residue", 1, -5, 5.0, 30]
    plot_pop_new(ddir, y, data_in, 7)

    y = ["RNA_residues_pop_RsR", "RsR", "Number of Residue", 1, 0, 0.7, 25]
    plot_pop_new(ddir, y, data_in, 6)

    y = ["RNA_residues_pop_DCC", "DCC", "Number of Residue", 1, 0.5, 1.2, 25]
    plot_pop_new(ddir, y, data_in, 5)

    #  correlations
    shell = data_bin(0.6, 5.0, 20)  # data change with resolution

    y = ["RNA_residue_cor_Zscor", "Resolution", "RsR_Zscore", 1, -5, 5.0, 30]
    plot_corr_new(ddir, y, data_in, 2, 7, shell)

    y = ["RNA_residues_cor_RsR", "Resolution", "RsR", 1, 0, 0.7, 25]
    plot_corr_new(ddir, y, data_in, 2, 6, shell)

    y = ["RNA_residues_cor_DCC", "Resolution", "DCC", 1, 0.5, 1.2, 25]
    plot_corr_new(ddir, y, data_in, 2, 5, shell)

    shell = data_bin(0.10, 0.5, 20)  # data change with Rwork
    #    ddir='DATA_PLOT_R/'  #create a dir to hold all the data/graphs
    #    os.system('mkdir %s' %ddir) #make a fold to hold all the data

    y = ["RNA_residue_cor_rwork_Zscor", "Rwork", "RsR_Zscore", 1, -5, 5.0, 30]
    plot_corr_new(ddir, y, data_in, 4, 7, shell)

    y = ["RNA_residue_cor_rwork_RsR", "Rwork", "RsR", 1, 0, 0.7, 25]
    plot_corr_new(ddir, y, data_in, 4, 6, shell)

    y = ["RNA_residue_cor_rwork_DCC", "Rwork", "DCC", 1, 0.1, 1.2, 25]
    plot_corr_new(ddir, y, data_in, 4, 5, shell)


######################################################
def plot_pop_new(ddir, x, data_all, ncol):
    """Plot the data populations (by manually selected boundary)
    data_all must be a list of list
    data_all: [[col1, col2.. ], [col1, col2 ..] ..]
    """

    item, xlabel, ylabel, _data_type, minv, maxv, nstep = x[0], x[1], x[2], x[3], x[4], x[5], x[6]  # noqa: F841
    shell = data_bin(minv, maxv, nstep)

    data1 = []
    for y in data_all:
        if ncol > len(y) or "?" in y[ncol] or not util.is_number(y[ncol]):
            print("Warning: problem with data row", y)
            continue
        data1.append([y[0], float(y[ncol])])

    data, outlier = filter_data(data1, shell, 1)  # only get the 1, 2th column
    fname = ddir + item + "_outlier.html"
    get_outlier_html(fname, outlier)

    fpop = ddir + item + "_pop.data"  # for data populations
    data_pop(fpop, data, shell, 1)  # 2th column
    avg, dev, _mini, _maxi = util.mean_dev(data, 1)

    title = "Population of %s (mean=%.2f; sigma=%.2f; entry=%d)" % (x[0], avg, dev, len(data))
    xrange, yrange, bar, rot, key, style = "", "", 0, 1, 0, 0
    plot = """plot '%s' using 3:xtic(1) lc rgb "blue" ,'' u 0:3:4 with labels offset 0, 0.5""" % (fpop)
    _gnuscr, _gnuout = gnu_plot(fpop, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)

    print("plot_pop=%s %s %s" % (x[1], len(data1), len(data_all)))


######################################################
def plot_corr_new(ddir, x, data_all, nc1, nc2, shell):
    """get the data correlations with resolution"""

    #    cate, token, data_type,  minv, maxv, nstep = x[0],x[1],x[2],x[3],x[4],x[5]
    item, xlabel, ylabel, minv, maxv = x[0], x[1], x[2], x[4], x[5]
    # data_type, nstep = x[3], x[6]

    data1 = []
    for y in data_all:
        if nc1 > len(y) or nc2 > len(y) or "?" in y[nc1] or "?" in y[nc2] or not util.is_number(y[nc1]) or not util.is_number(y[nc2]):
            print("Warning: problem with data row %s" % y)
            continue
        if float(y[nc2]) < minv or float(y[nc2]) > maxv:
            continue
        data1.append([y[0], float(y[nc1]), float(y[nc2])])
    #        print( nc1, nc2,y[nc1],y[nc2], y)

    fcorr = ddir + item + "_cor.data"
    data_corr(fcorr, data1, shell, 1)

    title = "Mean value of %s in each %s bin" % (x[2], x[1])
    xrange, yrange, bar, rot, key, style = "", "", 1, 1, 0, 0
    plot = """plot '%s' using 4:5:xtic(1) lc rgb "green" """ % (fcorr)
    _gnuscr, _gnuout = gnu_plot(fcorr, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)

    print("plot_corr=%s %s %s" % (x[1], len(data1), len(data_all)))


######################################################
def do_num_stat(ddir):
    """This is a large collection of numerically related plots"""
    data_all = []
    items = get_items()  # all the data items
    for x in items:  # stats for remark 3,200
        if "entity" in x[0]:
            continue  # special treatment
        if "torsion_outlier" in x[0] or "torsion_outlier" in x[1]:
            continue
        elif "start:" in x[0]:
            cate = x[0].split(":")[1]
            ofile = cate + ".all"
            get_file_from_sql(cate, ofile)
            data_all = from_file_to_list(ofile)

            print("%s %s" % (cate, len(data_all)))
            continue
        elif "end:" in x[0]:
            continue

        plot_pop(ddir, x, data_all)
        plot_corr(ddir, x, data_all)


######################################################
def plot_pop(ddir, x, data_all):
    """Plot the data populations (by manually selected boundary)
    data_all must be directly from SQL without change.
    """

    cate, token, minv, maxv, nstep = x[0], x[1], x[3], x[4], x[5]
    shell = data_bin(minv, maxv, nstep)

    if x[1] not in data_all[0]:
        print("Error, %s is not in the data list" % x[1])
        return

    i = data_all[0].index(x[1])
    data1 = []
    for y in data_all:
        if i >= len(y) or not util.is_number(y[i]):
            continue
        val = float(y[i])
        if val == 0.0:
            continue  # remove empty values
        data1.append([y[0], val])

    data, outlier = filter_data(data1, shell, 1)  # only get the 1, 2th column
    fname = ddir + cate + "_" + token + "_outlier.html"
    get_outlier_html(fname, outlier)

    fpop = ddir + cate + "_" + token + "_pop.data"  # for data populations
    data_pop(fpop, data, shell, 1)  # 2th column
    avg, dev, _mini, _maxi = util.mean_dev(data, 1)

    title = "population of _%s.%s (mean=%.2f; dev=%.2f; entry=%d)" % (x[0], x[1], avg, dev, len(data))
    xrange, yrange, xlabel, ylabel = "", "", "data range of %s" % x[1], "number of entry"
    bar, rot, key, style = 0, 1, 0, 0
    plot = """plot '%s' using 3:xtic(1) lc rgb "blue" ,'' u 0:3:4 with labels offset 0, 0.5""" % (fpop)
    _gnuscr, _gnuout = gnu_plot(fpop, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)

    print("plot_pop= %s %s %s %s" % (x[1], i, len(data1), len(data_all)))


######################################################
def plot_corr(ddir, x, data_all):
    """get the data correlations with resolution"""

    cate, token = x[0], x[1]
    # minv, maxv, nstep = x[3], x[4], x[5]

    res = "ls_d_res_high"
    if x[1] not in data_all[0] or res not in data_all[0]:
        print("Error, %s (or %s) is not in the data list" % (x[1], res))
        return

    i0 = data_all[0].index(res)
    i = data_all[0].index(x[1])
    data1 = []
    for y in data_all:
        if i >= len(y) or not util.is_number(y[i]) or not util.is_number(y[i0]):
            continue
        val = float(y[i])
        val0 = float(y[i0])
        if val == 0.0 or val0 == 0.0 or val < x[3] or val > x[4]:
            continue
        data1.append([y[0], val0, val])

    shell = data_bin(0.6, 5.0, 20)  # data change with resolution
    fcorr = ddir + cate + "_" + token + "_cor.data"
    data_corr(fcorr, data1, shell, 1)

    title = "The mean value of %s in each resolution bin" % (x[1])
    xrange, yrange, xlabel, ylabel = "", "", "resolution", "%s" % x[1]
    bar, rot, key, style = 1, 1, 0, 0
    plot = """plot '%s' using 4:5:xtic(1) lc rgb "green" """ % (fcorr)
    _gnuscr, _gnuout = gnu_plot(fcorr, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)

    print("plot_corr= %s %s %s %s" % (x[1], i, len(data1), len(data_all)))


######################################################
def from_file_to_list(outf):
    """input a file and output a list of list"""

    data_out = []
    fp = open(outf, "r")
    for x in fp:
        #    if 'structure_id' in x[:20] : continue
        t = x.strip().split()
        data_out.append(t)
    return data_out


######################################################
def get_file_from_sql(cate, outf):
    """get files from database (mainly for remark 3, 200)"""

    util.delete_file(outf)

    item_all = get_items()
    items = []
    for x in item_all:
        if cate != x[0]:
            continue
        items.append(x[1])

    print("items= %s" % items)

    y1 = '"select distinct r.structure_id, r.ls_d_res_high, '
    t = "".join(["x.%s, " % x for x in items])  # add more
    y2 = t.strip()[:-1]
    y3 = " from  refine as r , pdb_entry as p , %s as x  " % cate
    extra = " p.method like '%x-ray%' and p.status_code='REL' and p.method !='THEORETICAL MODEL' "
    y4 = " where  r.structure_id = x.structure_id and x.structure_id=p.pdb_id and %s " % extra
    y5 = 'order by  r.structure_id , x.%s" >%s' % (items[0], outf)

    s = "".join([SQL, y1, y2, y3, y4, y5])
    print("table=(%s); sql=%s" % (cate, s))
    os.system(s)


######################################################
def do_sql_gen(outf, tabs, type):  # pylint: disable=redefined-builtin
    """This is a general script to do SQL for any tables items
    outf: the output file from doing sql
    tabs: the tables (format : [[cate, item1,..], [cate, item1,..] ...]
    type: 'all' for all, 'x-ray' for XRAY, 'nmr' for NMR,
    'NEUTRON DIFFRACTION, X-RAY DIFFRACTION' for xray&neut ,
    'NEUTRON DIFFRACTION' for neut, 'POWDER' for POWDER DIFFRACTION
    'ELECTRON MICROSCOPY' for EM, 'FIBER' for FIBER DIFFRACTION
    """

    #    neut= " and  method not like '%neut%' " #exclude xry+neut
    pe = "pdb_entry"
    x1 = " %s.status_code='REL' and %s.method !='THEORETICAL MODEL'" % (pe, pe)
    if "all" in type:
        extra = x1
    else:
        extra = " pdb_entry.method like '%%%s%%' and " % type + x1

    cate = []
    sel = ""
    for x in tabs:
        for i, _y in enumerate(x):
            if i == 0:
                cate.append(x[0])
            else:
                sel = sel + "%s.%s, " % (x[0], x[i])
    sel = "select distinct " + sel
    # sel= 'select  ' + sel

    y0 = SQL
    y1 = sel.strip()[:-1]
    tmp = ""
    for x in cate:
        tmp = tmp + x + ","
    y2 = " from " + tmp.strip()[:-1] + " , pdb_entry "

    tmp = ""
    for i, x in enumerate(cate):
        if i == 0:
            continue
        tmp = tmp + " %s.structure_id= %s.structure_id and " % (cate[i - 1], cate[i])

    rev = " database_PDB_rev.num =1 and "
    if "database_PDB_rev" not in cate:
        rev = ""

    y3 = " where " + tmp + rev
    y4 = " %s.structure_id =pdb_entry.pdb_id and %s " % (cate[0], extra)
    arg = y0 + ' " ' + y1 + y2 + y3 + y4 + ' " ' + ">%s" % outf

    os.system(arg)
    print("arg=%s" % arg)


######################################################
def get_items():
    """all the annotated items
    numbers are float/int/str, low limit, upper limit, steps
    """

    items = [
        ["start:refine", "", 0, 0, 0, 0],
        ["refine", "ls_d_res_high", 1, 0.5, 6.0, 20],
        ["refine", "ls_d_res_low", 1, 5, 150, 20],
        ["refine", "ls_R_factor_R_work", 1, 0.08, 0.40, 20],
        ["refine", "ls_R_factor_R_free", 1, 0.08, 0.40, 20],
        ["refine", "overall_SU_R_Cruickshank_DPI", 1, 0.005, 0.95, 20],
        ["refine", "overall_SU_R_free", 1, 0.005, 0.95, 20],
        ["refine", "solvent_model_param_bsol", 1, 2, 250, 20],
        ["refine", "solvent_model_param_ksol", 1, 0.1, 0.95, 20],
        ["refine", "B_iso_mean", 1, 2, 200, 20],
        ["refine", "correlation_coeff_Fo_to_Fc", 1, 0.6, 1.1, 20],
        ["refine", "correlation_coeff_Fo_to_Fc_free", 1, 0.6, 1.1, 20],
        ["refine", "pdbx_overall_ESU_R", 1, 0.01, 0.95, 20],
        ["refine", "pdbx_overall_ESU_R_Free", 1, 0.01, 0.95, 20],
        ["refine", "overall_SU_B", 1, 0.1, 60, 20],
        ["refine", "overall_SU_ML", 1, 0.01, 0.95, 20],
        ["refine", "ls_percent_reflns_obs", 1, 50, 105, 20],
        ["refine", "ls_percent_reflns_R_free", 1, 0.5, 20, 20],
        ["refine", "ls_number_reflns_all", 1, 1000, 500000, 30],
        ["refine", "ls_number_reflns_obs", 1, 1000, 500000, 30],
        ["refine", "ls_number_reflns_R_free", 1, 20, 50000, 30],
        ["refine", "ls_number_reflns_R_work", 1, 1000, 150000, 30],
        ["refine", "ls_number_restraints", 1, 100, 50000, 30],
        ["refine", "occupancy_max", 1, 0.4, 1.1, 20],
        ["refine", "occupancy_min", 1, 0.05, 1.1, 20],
        ["refine", "pdbx_ls_sigma_I", 1, -6, 6, 20],
        ["refine", "pdbx_ls_sigma_F", 1, -6, 6, 20],
        ["refine", "pdbx_data_cutoff_high_absF", 1, 10000, 100000000, 30],
        ["refine", "pdbx_data_cutoff_low_absF", 1, 0, 4, 20],
        ["end:refine", "", 0, 0, 0, 0],
        ["start:refine_ls_shell", "", 1, 0.1, 0.95, 20],
        ["refine_ls_shell", "d_res_high", 1, 0.5, 5.5, 20],
        ["refine_ls_shell", "d_res_low", 1, 0.5, 6.5, 20],
        ["refine_ls_shell", "percent_reflns_obs", 1, 30, 106, 20],
        ["refine_ls_shell", "percent_reflns_R_free", 1, 0.5, 30, 20],
        ["refine_ls_shell", "R_factor_all", 1, 0.08, 0.70, 20],
        ["refine_ls_shell", "R_factor_R_free", 1, 0.08, 0.70, 20],
        ["refine_ls_shell", "R_factor_R_free_error", 1, 0.01, 0.2, 20],
        ["refine_ls_shell", "R_factor_R_work", 1, 0.08, 0.70, 20],
        ["refine_ls_shell", "pdbx_total_number_of_bins_used", 1, 2, 40, 20],
        ["refine_ls_shell", "number_reflns_all", 1, 100, 30000, 30],
        ["refine_ls_shell", "number_reflns_obs", 1, 100, 40000, 30],
        ["refine_ls_shell", "number_reflns_R_free", 1, 10, 2000, 30],
        ["refine_ls_shell", "number_reflns_R_work", 1, 100, 40000, 30],
        ["end:refine_ls_shell", "", 1, 0.1, 0.95, 20],
        ["start:reflns", "", 1, 0.1, 0.95, 20],
        ["reflns", "d_resolution_high", 1, 0.5, 6.5, 20],
        ["reflns", "d_resolution_low", 1, 5, 150, 20],
        ["reflns", "B_iso_Wilson_estimate", 1, 2, 200, 20],
        ["reflns", "observed_criterion_sigma_F", 1, -5, 5, 20],
        ["reflns", "observed_criterion_sigma_I", 1, -5, 5, 20],
        ["reflns", "percent_possible_obs", 1, 40, 105, 20],
        ["reflns", "pdbx_redundancy", 1, 0.5, 40, 20],
        ["reflns", "pdbx_netI_over_av_sigmaI", 1, 1.0, 50, 20],
        ["reflns", "pdbx_netI_over_sigmaI", 1, 1, 50, 20],
        ["reflns", "Rmerge_F_obs", 1, 0.01, 0.90, 20],
        ["reflns", "pdbx_Rmerge_I_obs", 1, 0.01, 1.0, 20],
        ["reflns", "pdbx_Rsym_value", 1, 0.01, 1.0, 20],
        ["reflns", "pdbx_Rrim_I_all", 1, 0.01, 1.0, 20],
        ["reflns", "pdbx_Rpim_I_all", 1, 0.01, 1.0, 20],
        ["reflns", "pdbx_number_measured_all", 1, 1000, 800000, 30],
        ["reflns", "number_all", 1, 1000, 800000, 30],
        ["reflns", "number_obs", 1, 1000, 800000, 30],
        ["end:reflns", "", 1, 0.1, 0.95, 20],
        ["start:reflns_shell", "", 1, 0.1, 0.95, 20],
        ["reflns_shell", "d_res_high", 1, 0.5, 5.0, 20],
        ["reflns_shell", "d_res_low", 1, 0.5, 6.0, 20],
        ["reflns_shell", "percent_possible_all", 1, 40, 105, 20],
        ["reflns_shell", "percent_possible_obs", 1, 40, 105, 20],
        ["reflns_shell", "meanI_over_sigI_obs", 1, 0.5, 50, 20],
        ["reflns_shell", "pdbx_redundancy", 1, 0.5, 50, 20],
        ["reflns_shell", "Rmerge_F_obs", 1, 0.02, 1.0, 20],
        ["reflns_shell", "Rmerge_I_all", 1, 0.02, 1.0, 20],
        ["reflns_shell", "Rmerge_I_obs", 1, 0.02, 1.0, 20],
        ["reflns_shell", "pdbx_Rsym_value", 1, 0.02, 0.95, 20],
        ["reflns_shell", "pdbx_Rrim_I_all", 1, 0.02, 0.95, 20],
        ["reflns_shell", "number_measured_all", 1, 100, 80000, 30],
        ["reflns_shell", "number_measured_obs", 1, 100, 80000, 30],
        ["reflns_shell", "number_unique_all", 1, 100, 40000, 30],
        ["reflns_shell", "number_unique_obs", 1, 100, 40000, 30],
        ["end:reflns_shell", "", 1, 0.1, 0.95, 20],
        ["start:exptl_crystal", "", 0, 0, 0, 0],
        ["exptl_crystal", "pdbx_mosaicity", 1, 0.1, 4.0, 20],
        ["exptl_crystal", "density_Matthews", 1, 1, 8.0, 20],
        ["exptl_crystal", "density_percent_sol", 1, 15, 92.0, 20],
        ["end:exptl_crystal", "", 0, 0, 0, 0],
        ["start:exptl_crystal_grow", "", 0, 0, 0, 0],
        ["exptl_crystal_grow", "pH", 1, 2, 11.0, 20],
        ["exptl_crystal_grow", "temp", 1, 273, 310.0, 30],
        ["end:exptl_crystal_grow", "", 1, 0.1, 4.0, 20],
        ["start:entity", "", 0, 0, 0, 0],
        ["entity", "water-weight-in-ASU", 1, 80, 60000, 20],
        ["entity", "polymer-weight", 1, 500, 100000, 30],
        ["entity", "non-polymer-weight", 1, 30, 3000, 30],
        ["entity", "total-weight-in-ASU", 1, 500, 800000, 30],
        ["end:entity", "", 1, 0.1, 4.0, 20],
        ["start:diffrn_radiation_wavelength", "", 0, 0, 0, 0],
        ["diffrn_radiation_wavelength", "wavelength", 1, 0.5, 2, 20],
        ["end:diffrn_radiation_wavelength", "", 1, 0.1, 4.0, 20],
        ["start:torsion_outlier", "", 0, 0, 0, 0],
        ["Amino_Acid", "torsion_outlier", 1, 0.0, 40, 20],  # percentage
        ["end: ", "", 0, 0, 0, 0],
    ]

    return items


######################################################
def get_outlier_html(fname, data):
    """write the outliers into a html file, data has two columns"""

    item = fname.replace("_outlier.html", "")

    url = "http://www.rcsb.org/pdb/explore/explore.do?structureId"
    fw = open(fname, "w")

    fw.write("<!DOCTYPE html> \n<html>\n<body>\n")
    fw.write('<style type="text/css">a {text-decoration: none}</style>\n')
    s1 = (
        "The outlier for the data item (%s) are excluded from the plot. They are not necessary wrong, \
    but the wrong data should be in the outlier values!<p>\n"
        % item
    )
    fw.write(s1)
    fw.write("<table>\n")

    data.sort(key=lambda y: y[1])
    fw.write("<tr> <td> PDBID </td> <td> OUTLIERS  </td></tr>\n")
    for x in data:
        t1 = '<tr> <td> <a href="%s=%s" target="dynamic"> %s </a> </td>' % (url, x[0], x[0])
        t2 = "<td> (%s) </td></tr>\n" % x[1]
        fw.write(t1 + t2)

    fw.write("</table> \n</body>\n</html>\n")
    fw.close()


######################################################
def html_head():
    """get the header of HTML"""

    head = """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
  <title>pdbitem statistics</title>
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
  <META NAME="description" CONTENT="pdbitem_stat is a tool to present the annotated pdb items in term of histogram. It shows how the data is populated in the database, and how the data is correlated with resolution, and how the data grows annually.">  # noqa: E501
  <META NAME="keywords" CONTENT=" statistics, population">
  <link href="/styles/general.css" rel="STYLESHEET" type="text/css">
</head>

<body>
<style type="text/css">a {text-decoration: none}</style> <!-- rid of underline for the link -->

<script language="javascript">
<!--
today = new Date();
/*document.write("<BR>The time now is: ", today.getHours(),":",today.getMinutes());*/
document.write("(Last update: ", today.getDate(),"/",today.getMonth()+1,"/",today.getFullYear(), ")");
//-->
</script>

<a name="top"></a>

<center><h1>A statistics tool for the annotated pdb items.</h1> </center> <!-- not underline -->


<li><a href="#info"> <b>A note for the distograms </b></a></li>
<li><a href="#refine"> <b>REFINE </b></a></li>
<li><a href="#refine_ls_shell"> <b> REFINE_LS_SHELL </b></a></li>
<li><a href="#reflns"> <b>REFLNS </b></a></li>
<li><a href="#reflns_shell"> <b>REFLNS_SHELL  </b></a></li>
<li><a href="#exptl_crystal"> <b>EXPTL_CRYSTAL  </b></a></li>
<li><a href="#exptl_crystal_grow"> <b>EXPTL_CRYSTAL_GROW </b></a></li>
<li><a href="#entity"> <b>ENTITY  </b></a></li>
<li><a href="#diffrn_radiation_wavelength"> <b>DIFFRN_RADIATION_WAVELENGTH </b></a></li>
<li><a href="#string"> <b>ITEMS WITH STRING VALUES </b></a></li>
<li><a href="#torsion_outlier"> <b>AMINO ACID TORSIONAL ANGLE OUTLIERS </b></a></li>
<li><a href="#other"> <b>OTHER </b></a></li>

<hr>

<h2>A note for the distograms  <a name="info"  href="#top" style="font-size:medium;">&nbsp;  (top)</a> </h2>


 <p>Tables below listed the statistics for some of the annotated cif
    items in terms of histograms. They show how the data is populated in the database,
    and how the data changes with  resolution.
    The outliers are excluded from the plot and are also given in the table.

    <p><b>For the graph of population:</b>
    The number on the top of each bar is the percentage of the population of the
    item. It was calculated by the number of the entries in the data range at the
    bottom of the graph divided by the total number of entries having proper values
    (outlier values excluded). The overall average is the summation of the values
    divided by the total number of entries (outliers excluded). The numbers on
    the left axis are the number of entries corresponding to each bar.

    <p><b>For the graph of correlation:</b>
    Each bar is the mean value (vertical axis) of the data item in the resolution
    range (horizontal axis). The error bars (the standard diviation) are given to each value.

    <p>

"""
    return head


######################################################
def get_html(ddir):
    """get a html table for the statistics tool"""

    fw = open("index.html", "w")

    tmp = """
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Frameset//EN"
   "http://www.w3.org/TR/html4/frameset.dtd">
<HTML>

<HEAD>
<TITLE>A simple frameset document</TITLE>
</HEAD>
<!--<frameset rows="200,*">-->
<frameset rows="30%, 70%">
  <frame name="top" scrolling="yes" noresize target="bottom" src="bottom.html">
  <frame name="bottom" scrolling="yes" noresize target="_self" src="temp.html">

</frameset>

</HTML>

    """
    fw.write(tmp)
    fw.close()

    fw = open("temp.html", "w")
    tmp = """
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Frameset//EN"
   "http://www.w3.org/TR/html4/frameset.dtd">
<HTML>

Please select an data item in the above panel to show the graphs.

</HTML>

    """
    fw.write(tmp)
    fw.close()

    # refer = {
    #     "refine": "refine",
    #     "refine_ls_shell": "refine_ls_shell",
    #     "reflns": "reflns",
    #     "refine_ls_shell": "refine_ls_shell",
    #     "exptl_crystal": "exptl_crystal",
    #     "exptl_crystal_grow": "exptl_crystal_grow",
    #     "entity": "entity",
    #     "diffrn_radiation_wavelength": "diffrn_radiation_wavelength",
    #     "string": "string",
    #     "other": "other",
    # }

    #    html='index.html'
    html = "bottom.html"
    fw = open(html, "w")

    items = get_items()
    head = html_head()
    fw.write(head)

    sp = "&nbsp;&nbsp; "
    nl = 0
    for x in items:
        if "start:" in x[0]:
            s1 = x[0].split(":")[1]
            s2 = (
                '<p><h3> Data population and correlation with resolution for the cif table of <i>%s</i> <a name="%s"  href="#top" style="font-size:medium;">&nbsp;  (top)</a>  </h3>\n'
                % (s1.upper(), s1)
            )
            fw.write(s2)
            fw.write("<table>\n")
            s3 = (
                "<tr><td><b>Data population %s</b></td> <td > data4plot  %s</td> <td> <b>Correlation with resolution %s</b></td> <td > data4plot %s</td> <td> outiliers </td></tr>\n"
                % (sp, sp, sp, sp)
            )
            fw.write(s3)

        elif "end:" in x[0]:
            fw.write("</table>\n")

        else:
            nl = nl + 1
            bgcolor = ' "" '
            if nl % 2 == 0:
                bgcolor = ' "#ccffff" '  # even number

            fw.write("<tr>\n")
            y1 = '<td bgcolor=%s ><a href="%s/%s_%s_pop.data.png" target="bottom">_%s.%s %s</a></td> \n' % (bgcolor, ddir, x[0], x[1], x[0], x[1], sp)
            y2 = '<td bgcolor=%s ><a href="%s/%s_%s_pop.data" target="bottom"> (data) %s </a></td> \n' % (bgcolor, ddir, x[0], x[1], sp)
            y3 = '<td bgcolor=%s ><a href="%s/%s_%s_cor.data.png" target="bottom">change with resolution  %s </a></td>\n' % (bgcolor, ddir, x[0], x[1], sp)
            y4 = '<td bgcolor=%s ><a href="%s/%s_%s_cor.data" target="bottom">(data) %s  </a></td> \n' % (bgcolor, ddir, x[0], x[1], sp)
            y5 = '<td bgcolor=%s ><a href="%s/%s_%s_outlier.html" target="bottom">outliers </a></td> \n' % (bgcolor, ddir, x[0], x[1])
            fw.write("".join([y1, y2, y3, y4, y5]))
            fw.write("</tr>\n")

    fw.write("\n<hr> \n")  # for the string items
    str_item = get_items_str()

    nl = 0
    for x in str_item:
        if "start:" in x[0]:
            s1 = x[0].split(":")[1]
            s2 = (
                '<p><h3> Data population and yearly growth and the growth rate for the <i>%s</i> items.  <a name="%s"  href="#top" style="font-size:medium;">&nbsp;  (top)</a> </h3>  <p> (Numbers on the top of the bar is the percentage calculated from the number of population and the total non-null entries. The yearly growth rate is calculated  by the entry number divided total number for the year.)  <p>\n'  # noqa: E501
                % (s1.upper(), s1)
            )
            fw.write(s2)
            fw.write("<table>\n")
            s3 = (
                "<tr><td><b>Data population %s</b></td>  <td >data4plot  %s</td>  <td><b>yearly growth %s</b></td>  <td> data4plot %s</td>   <td><b>yearly growth rate %s</b></td>   <td> data4plot </td></tr>\n"  # noqa: E501
                % (sp, sp, sp, sp, sp)
            )
            fw.write(s3)

        elif "end:" in x[0]:
            fw.write("</table>\n")

        else:
            nl = nl + 1
            bgcolor = ' "" '
            if nl % 2 == 0:
                bgcolor = ' "#ccffff" '  # even number

            fw.write("<tr>\n")
            y1 = '<td bgcolor=%s ><a href="%s/%s_%s.all_pop.data.png" target="bottom">_%s.%s %s</a></td> \n' % (bgcolor, ddir, x[0], x[1], x[0], x[1], sp)
            y2 = '<td bgcolor=%s ><a href="%s/%s_%s.all_pop.data" target="bottom"> (data) %s </a></td> \n' % (bgcolor, ddir, x[0], x[1], sp)
            y3 = '<td bgcolor=%s ><a href="%s/%s_%s.all.data.png" target="bottom"> yearly growth  %s </a></td>\n' % (bgcolor, ddir, x[0], x[1], sp)
            y4 = '<td bgcolor=%s ><a href="%s/%s_%s.all.data" target="bottom">(data) %s  </a></td> \n' % (bgcolor, ddir, x[0], x[1], sp)
            y5 = '<td bgcolor=%s ><a href="%s/%s_%s.all.data_rate.png" target="bottom"> yearly growth rate %s </a></td>\n' % (bgcolor, ddir, x[0], x[1], sp)
            y6 = '<td bgcolor=%s ><a href="%s/%s_%s.all.data" target="bottom">(data) %s  </a></td> \n' % (bgcolor, ddir, x[0], x[1], sp)
            fw.write("".join([y1, y2, y3, y4, y5, y6]))
            fw.write("</tr>\n")

    fw.write("\n<hr> \n")  # for other items

    fw.write('<p> <h3> Other statistics  <a name="other"  href="#top" style="font-size:medium;">&nbsp;  (top)</a> </h3> \n')

    s2 = "<p>PDB accumulate growth and prediction of the growth in next five years (DEPOSITED)."
    y1 = '<a href="%s/pdb_growth.txt_acum_DEP.data.png" target="bottom">%s %s</a> \n' % (ddir, s2, sp)
    y2 = '<a href="%s/pdb_growth.txt_acum_DEP.data" target="bottom"> (plot data) </a> \n' % (ddir)
    fw.write(y1 + y2)

    s2 = "<p>PDB accumulate growth and prediction of the growth in next five years (RELEASED)."
    y1 = '<a href="%s/pdb_growth.txt_acum_REL.data.png" target="bottom">%s %s</a> \n' % (ddir, s2, sp)
    y2 = '<a href="%s/pdb_growth.txt_acum_REL.data" target="bottom"> (plot data) </a> \n' % (ddir)
    fw.write(y1 + y2)

    s2 = "<p>PDB yearly growth and prediction of the growth in next five years (DEPOSITED)."
    y1 = '<a href="%s/pdb_growth.txt_year_DEP.data.png" target="bottom">%s %s</a> \n' % (ddir, s2, sp)
    y2 = '<a href="%s/pdb_growth.txt_year_DEP.data" target="bottom"> (plot data) </a> \n' % (ddir)
    fw.write(y1 + y2)

    s2 = "<p>PDB yearly growth and prediction of the growth in next five years (RELEASED)."
    y1 = '<a href="%s/pdb_growth.txt_year_REL.data.png" target="bottom">%s %s</a> \n' % (ddir, s2, sp)
    y2 = '<a href="%s/pdb_growth.txt_year_REL.data" target="bottom"> (plot data) </a> \n' % (ddir)
    fw.write(y1 + y2)

    s2 = "<p> Population of resolution with different entities"
    y1 = '<a href="%s/resolution_entity.png" target="bottom">%s %s</a> \n' % (ddir, s2, sp)
    y2 = '<a href="%s/resolution_entity.data" target="bottom"> (plot data) </a> \n' % (ddir)
    fw.write(y1 + y2)

    s2 = "<p> R_factor (different resolution groups) change in each year. "
    y1 = '<a href="%s/rfact_reso_growth.all.data.png" target="bottom">%s %s</a> \n' % (ddir, s2, sp)
    y2 = '<a href="%s/rfact_reso_growth.all.data" target="bottom"> (plot data) </a> \n' % (ddir)
    fw.write(y1 + y2)

    s2 = "<p> Resolution (different data ranges) change in each year. "
    y1 = '<a href="%s/rfact_reso_growth.all_res.data.png" target="bottom">%s %s</a> \n' % (ddir, s2, sp)
    y2 = '<a href="%s/rfact_reso_growth.all_res.data" target="bottom"> (plot data) </a> \n' % (ddir)
    fw.write(y1 + y2)

    s1 = "<p>yearly growth of twinned entry, "
    s2 = "space group dependance, "
    s3 = "comparison of R factors (twinned vs untwinned)"
    y1 = '<a href="%s/twin-info.all_growth.data.png" target="bottom">%s %s</a> \n' % (ddir, s1, sp)
    y2 = '<a href="%s/twin-info.all_sg_pop.data.png" target="bottom">%s %s</a> \n' % (ddir, s2, sp)
    y3 = '<a href="%s/twin-nontwin-rwork.data.png" target="bottom"> %s %s </a> \n' % (ddir, s3, sp)
    fw.write(y1 + y2 + y3)

    s1 = "<p>comparison of R factors produced by the popular programs "
    y1 = '<a href="%s/software_rwork.data.png" target="bottom">%s %s</a> \n' % (ddir, s1, sp)
    fw.write(y1)

    s2 = "<p>comparison of R_work and Rfree in various resolution range. "
    y2 = '<a href="%s/rwork-rfree.data.png" target="bottom">%s %s</a> \n' % (ddir, s2, sp)
    fw.write(y2)

    s2 = "<p>Growth of the entries refined with TLS "
    y2 = '<a href="%s/tls-growth.all.data.png" target="bottom">%s %s</a> \n' % (ddir, s2, sp)
    fw.write(y2)

    fw.write("\n<hr> </body>\n</html>\n")
    fw.close()


######################################################################
def data_pop(fpop, fp, shell, col):
    """get the data populations using the COL of the list FP"""

    nt = len(fp)
    if nt == 0:
        return

    fw = open(fpop, "w")
    fp.sort(key=lambda y: y[col])  # sort fp by the col
    s1 = "#data_range  avg_bin   entry_number  percentage\n"
    print(s1)
    fw.write(s1)

    nstep = len(shell)
    j = 0
    nnt = 0
    for i in range(nstep):
        if i == 0:
            continue
        bin = (shell[i - 1] + shell[i]) / 2.0  # pylint: disable=redefined-builtin
        n = 0
        for k in range(j, nt):
            if fp[k][col] > shell[i]:
                j = k
                break
            elif shell[i] > fp[k][col] >= shell[i - 1]:
                n = n + 1
                nnt = nnt + 1

        p = 100 * float(n) / nt
        print('pop:  "%8.2f %8.2f" %8.2f  %5d %6.1f ' % (shell[i - 1], shell[i], bin, n, p))
        fw.write('"%.2f %.2f" %8.2f  %5d %6.1f \n' % (shell[i - 1], shell[i], bin, n, p))

    fw.close()
    print("selected data (%s)= %d  %d " % (fpop, nt, nnt))


######################################################################
def data_corr(fcorr, fp, shell, col):
    """get data correlations between column of col and col+1
    fp: is  a list of list; fcorr is a out filename,
    shell is the value of column (col) break into bin.
    """

    fw = open(fcorr, "w")
    fp.sort(key=lambda y: y[col])

    nt = len(fp)
    nstep = len(shell)
    j = 0
    s1 = "#data_range  avg_bin  number  mean  deviation  minimum  maximum\n"
    print(s1)
    fw.write(s1)
    for i in range(nstep):
        if i == 0:
            continue

        bin = (shell[i - 1] + shell[i]) / 2.0  # pylint: disable=redefined-builtin
        n, subg = 0, []
        for k in range(j, nt):
            if fp[k][col] > shell[i]:
                j = k
                break
            elif shell[i] > fp[k][col] >= shell[i - 1]:
                subg.append(fp[k][col + 1])  # data after resolution
                n = n + 1

        if len(subg):
            avg, dev, mini, maxi = util.mean_dev(subg, -1)
            print('corr: "%8.2f %8.2f" %8.2f %5d %8.3f %8.3f %8.3f %8.3f' % (shell[i - 1], shell[i], bin, n, avg, dev, mini, maxi))
            fw.write('"%.2f %.2f" %8.2f %5d %8.3f %8.3f %8.3f %8.3f\n' % (shell[i - 1], shell[i], bin, n, avg, dev, mini, maxi))

    fw.close()


######################################################################
def clean_data(infile, idd):
    """infile has 3 columns (pdbid, resolution, items)
    idd=0: Do not clean data, idd==1: remove data with non-values
    it returns a list with proper values
    """

    print("Cleaning data (%s)." % infile)
    data = []
    fp = open(infile, "r")

    for x in fp:
        if "structure" in x:
            continue
        t = x.strip().split()
        n = len(t)
        if not util.is_number(t[1]) or (n > 2 and not util.is_number(t[2])):
            print("Error: value of (%s) is not a number." % x)
            continue
        if idd == 0:
            if n == 2:
                data.append([t[0], float(t[1])])
            elif n == 3:
                data.append([t[0], float(t[1]), float(t[2])])
        else:
            if n == 2:
                if float(t[1]) != 0.0:
                    data.append([t[0], float(t[1])])
            elif n == 3:
                if float(t[1]) != 0.0 and float(t[2]) != 0.0:
                    data.append([t[0], float(t[1]), float(t[2])])

    fp.close()

    return data


######################################################################
def filter_data(data_in, shell, col):
    """filter the data according to the shell & col values."""

    data_new, outlier = [], []
    for x in data_in:
        #        print(x[col], shell[0], shell[-1])
        if x[col] < shell[0] or x[col] > shell[-1]:
            outlier.append([x[0], x[col]])
        else:
            data_new.append(x)

    return data_new, outlier


######################################################################
def data_bin(lowv, upv, nstep):
    """Use the low and upper values and the nsteps to get a list"""
    shell = []
    if nstep == 0:
        return shell
    d = (upv - lowv) / float(nstep)
    a = lowv
    for x in range(1, nstep + 1):
        shell.append(a)
        a = lowv + d * x

    return shell


######################################################################
def plot_entity(ddir):
    """plot the entity"""

    out = "entity_data.all"
    refine = ["refine", "structure_id", "ls_d_res_high"]
    items = ["entity", "formula_weight", "id", "type", "pdbx_number_of_molecules "]
    tabs = [refine, items]
    do_sql_gen(out, tabs, "all")

    items = [
        ["entity", "water-weight-in-ASU", 1, 60, 50000, 20],
        ["entity", "polymer-weight", 1, 200, 100000, 30],
        ["entity", "non-polymer-weight", 1, 10, 1500, 25],
        ["entity", "total-weight-in-ASU", 1, 400, 600000, 30],
    ]

    data_all = from_file_to_list(out)
    data_water = [["structure_id", "ls_d_res_high", "water-weight-in-ASU"]]
    data_poly = [["structure_id", "ls_d_res_high", "polymer-weight"]]
    data_nonpoly = [["structure_id", "ls_d_res_high", "non-polymer-weight"]]
    data_allw = [["structure_id", "ls_d_res_high", "total-weight-in-ASU"]]
    nt = len(data_all)
    sm = 0
    for i, x in enumerate(data_all):
        if "structure_id" in x[0]:
            continue
        elif "water" == x[4]:
            wat = float(x[2]) * float(x[5])
            data_water.append([x[0], float(x[1]), wat])
        elif "non-polymer" == x[4]:
            data_nonpoly.append([x[0], float(x[1]), float(x[2])])
        elif "polymer" == x[4]:
            data_poly.append([x[0], float(x[1]), float(x[2])])

        sm = sm + float(x[2]) * float(x[5])
        if i < nt - 1 and x[0] not in data_all[i + 1][0]:
            data_allw.append([x[0], float(x[1]), sm])
            sm = 0

    for x in items:  #
        data_all = []
        if x[1] == "water-weight-in-ASU":
            data_all = data_water
        elif x[1] == "non-polymer-weight":
            data_all = data_nonpoly
        elif x[1] == "polymer-weight":
            data_all = data_poly
        elif x[1] == "total-weight-in-ASU":
            data_all = data_allw
        print(x)
        if not data_all:
            continue
        plot_pop(ddir, x, data_all)
        plot_corr(ddir, x, data_all)


######################################################################
def torsion_info(ddir):
    """The outliers of the psi and phi torsional angles
    return pdbid, resolution, percentage
    """
    connection = MySQLdb.connect(host="pdb-f-linux-2", user="rcsbuser", passwd="rcsb0000", db="cleanv1")
    cursor = connection.cursor()
    data = []

    aa = ["LEU", "ALA", "VAL", "GLU", "SER", "THR", "ASP", "ILE", "TYR", "PRO", "TRP", "ASN", "GLN", "ARG", "LYS", "PHE", "HIS", "CYS", "MET", "MSE"]

    ref = ["refine", "entry_id", "ls_d_res_high"]
    tabs = [ref]
    outf = ddir + "refine.all"
    do_sql_gen(outf, tabs, "x-ray")
    #    outf='testid.list'
    data.append(["pdbid", "ls_d_res_high", "torsion_outlier"])
    for i, x in enumerate(open(outf, "r").readlines()):
        naa = 0
        if i == 0:
            continue
        xt = x.split()
        pdbid, reso = xt[0], xt[1]

        query = "select  PDB_model_num, auth_comp_id from pdbx_validate_torsion where structure_id='%s'" % pdbid
        cursor.execute(query)
        rows = cursor.fetchall()
        num = 0
        for row in rows:
            if row[0] == 1 and row[1] in aa:
                num += 1
            elif "GLY" in row[1]:
                print("Note: %s has outlier for GLY!" % x.strip())
        if num == 0:
            continue

        query = "select asym_id, pdb_mon_id from pdbx_poly_seq_scheme where structure_id='%s'" % pdbid
        cursor.execute(query)
        rows = cursor.fetchall()
        nch = []
        for row in rows:
            if row[1] in aa:
                naa += 1
            if row[0] not in nch:
                nch.append(row[0])

        if naa == 0:
            print("Error: %s has no AA" % x.strip())
            continue
        per = 100 * float(num) / float(naa)
        #        print('%-15s  %3d  %5d %2d  %8.2f'%(x.strip(), num, naa, len(nch), per))
        data.append([pdbid, reso, "%.1f" % per])

    # #######ploting below
    y = ["Amino_Acid_torsion_outlier", "torsion_angle_outlier", "Number of Residue", 1, 0, 40.0, 20]
    plot_pop_new(ddir, y, data, 2)

    #  correlations
    shell = data_bin(0.6, 5.0, 20)  # data change with resolution
    y = ["Amino_Acid_torsion_outlier", "Resolution", "torsion_angle_outlier", 1, -5, 5.0, 30]
    plot_corr_new(ddir, y, data, 1, 2, shell)


#    return data


######################################################################
def tls_growth(ddir):
    tls = ["pdbx_refine_tls_group", "pdbx_refine_id"]
    pdbrev = ["database_PDB_rev", "structure_id", "date"]
    ref = ["refine", "ls_d_res_high"]
    tabs = [pdbrev, tls]
    outftls = ddir + "tls-growth.all"
    do_sql_gen(outftls, tabs, "x-ray")
    datatls = open(outftls, "r").readlines()

    tabs = [pdbrev, ref]
    outfref = ddir + "refine-growth.all"
    do_sql_gen(outfref, tabs, "x-ray")
    dataref = open(outfref, "r").readlines()

    plot_file = outftls + ".data"
    fw = open(plot_file, "w")
    s = "#year  number-of-tls-entry  percentage  number-of-xray-entry \n"
    fw.write(s)
    print(s.strip())

    for i in range(2000, YEAR):
        ntls, nref = 0, 0
        for x in datatls:
            if str(i) in x:
                ntls = ntls + 1
        for y in dataref:
            if str(i) in y:
                nref = nref + 1

        if nref == 0:
            continue
        p = 100 * float(ntls) / nref
        s = "%d  %5d %6.1f %5d\n" % (i, ntls, p, nref)
        fw.write(s)
        print(s.strip())
    fw.close()

    title = "Growth of entries refined with TLS(percentage on top bar)"
    xrange, yrange, xlabel, ylabel = "", "", "year", "number of entry with TLS"
    bar, rot, key, style = 0, 1, 0, 0

    plot = """plot '%s' using 2:xtic(1) , '' u 0:2:3 with labels offset 0, 0.5  \
      """ % (
        plot_file
    )
    _gnuscr, _gnuout = gnu_plot(plot_file, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)


######################################################################
def rfactor_other(ddir):
    ref = ["refine", "ls_d_res_high", "ls_R_factor_R_work"]
    comp = ["computing", "structure_refinement"]
    tabs = [ref, comp]

    outf = ddir + "software_rwork.all"
    do_sql_gen(outf, tabs, "x-ray")
    rawdata = from_file_to_list(outf)
    refmac, phenix, buster, cns = [], [], [], []
    for x in rawdata:
        if "d_res_high" in x[0] or float(x[1]) < 0.00001 or len(x) < 3:
            continue
        if "REFMAC" in x[2].upper():
            refmac.append([float(x[0]), float(x[1])])
        elif "PHENIX" in x[2].upper():
            phenix.append([float(x[0]), float(x[1])])
        elif "BUSTER" in x[2].upper():
            buster.append([float(x[0]), float(x[1])])
        elif "CNS" in x[2].upper():
            cns.append([float(x[0]), float(x[1])])

    shell = data_bin(0.9, 4.5, 20)  # data change with resolution
    prog = ["refmac", "phenix", "buster", "cns"]
    for i, x in enumerate([refmac, phenix, buster, cns]):
        fcorr = "rwork_%s.data" % (prog[i])
        data_corr(fcorr, x, shell, 0)

    title = "Comparison of R_work for different refinement programs"
    xrange, yrange, xlabel, ylabel = "", "", "resolution", "R_work"
    bar, rot, key, style = 0, 1, 2, 1
    plot = """plot 'rwork_%s.data' using 4:xtic(1) t "by REFMAC", \
    'rwork_%s.data' u 4 t "by PHENIX" , \
    'rwork_%s.data' u 4 t "by BUSTER" , \
    'rwork_%s.data' u 4 t "by CNS"  \
        """ % (
        prog[0],
        prog[1],
        prog[2],
        prog[3],
    )
    file = ddir + "software_rwork.data"
    _gnuscr, _gnuout = gnu_plot(file, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)

    # below is to plot the Rwork & Rfree together using previous data
    # If the file name changed, fwork/ffree must also change!!

    fwork = ddir + "refine_ls_R_factor_R_work_cor.data"
    ffree = ddir + "refine_ls_R_factor_R_free_cor.data"

    title = "Comparison of R_work & Rfree in different resolution ranges"
    xrange, yrange, xlabel, ylabel = "", "", "resolution", "R_work/Rfree"
    bar, rot, key, style = 0, 1, 2, 1
    plot = """plot '%s' using 4:xtic(1) t "R_work", \
    '%s' u 4 t "R_free"  """ % (
        fwork,
        ffree,
    )
    file = ddir + "rwork-rfree.data"
    _gnuscr, _gnuout = gnu_plot(file, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)


###################################################################
def twin_info(ddir):
    """population/stat of twin"""

    pdbrev = ["database_PDB_rev", "structure_id", "date_original"]
    twin = ["pdbx_reflns_twin", "fraction"]
    refine = ["refine", "ls_d_res_high", "ls_R_factor_R_work"]
    sym = ["symmetry", "space_group_name_H_M"]

    tabs = [pdbrev, twin, refine, sym]
    outf = ddir + "twin-info.all"
    do_sql_gen(outf, tabs, "x-ray")
    rawdata = from_file_to_list(outf)

    data1 = []
    for x in rawdata:
        # print(x)
        if "structur" in x[0] or x[3] == "0":
            continue
        s = "_".join(x[6:])
        data1.append([x[1][:4], x[4], x[5], s])

    data2 = uniq_list_of_list(data1)
    # print(len(data1), len(data2))
    data = []
    for x in data2:  # get year,sg, reso, Rfactor
        ss = [int(x[0]), x[3], float(x[1]), float(x[2])]
        data.append(ss)
        # print(ss)

    # print(data)
    shell = data_bin(1.0, 4.0, 15)  # data change with resolution
    fcorr = outf + ".data"
    data_corr(fcorr, data, shell, 2)

    fyear = outf + "_growth.data"
    fw = open(fyear, "w")
    s = "# year  number-of-twin-entry   percentage\n"
    fw.write(s)
    print(s)
    data3 = uniq_string_pop(data, 0)
    for x in data3:
        s = " ".join(["%5d" % x[0], "%5d" % x[1], "%6.1f" % x[2]]) + "\n"
        fw.write(s)
        print(s.strip())
    fw.close()

    fsg = outf + "_sg_pop.data"
    fw = open(fsg, "w")
    s = "# space-group  number-of-twin-entry   percentage\n"
    fw.write(s)
    print(s)
    data4 = uniq_string_pop(data, 1)
    for x in data4:
        s = " ".join(["%10s" % x[0], "%5d" % x[1], "%6.1f" % x[2]]) + "\n"
        fw.write(s)
        print(s.rstrip())

    fw.close()

    outf_reg = ddir + "regular-res-rfact.all"  # for general entry (twin-notwin)
    tabs = [["refine", "structure_id", "ls_d_res_high", "ls_R_factor_R_work"]]
    do_sql_gen(outf_reg, tabs, "x-ray")
    rawdata_reg = from_file_to_list(outf_reg)
    data_reg = []
    for x in rawdata_reg:
        if "structur" in x[0] or float(x[1]) == 0 or float(x[2]) == 0:
            continue
        data_reg.append([float(x[1]), float(x[2])])

    fcorr_reg = outf_reg + ".data"
    data_corr(fcorr_reg, data_reg, shell, 0)

    #  plot below
    title = "Yearly growth for twinned entry (percentage on each bar)"
    xrange, yrange, xlabel, ylabel = "", "", "Year", "number of twinned entry"
    bar, rot, key, style = 0, 1, 0, 0
    plot = """plot '%s' using 2:xtic(1)  lc rgb "green"  ,'' u 0:2:3 with labels offset 0, 0.5  \
        """ % (
        fyear
    )
    _gnuscr, _gnuout = gnu_plot(fyear, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)

    title = "Space group population with twinned data (percentage on each bar)"
    xrange, yrange, xlabel, ylabel = "", "", "space group", "number of twinned entry"
    bar, rot, key, style = 0, 1, 0, 0
    plot = """plot '%s' using 2:xtic(1) lc rgb "green"  ,'' u 0:2:3 with labels offset 0, 0.5  \
        """ % (
        fsg
    )
    _gnuscr, _gnuout = gnu_plot(fsg, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)

    title = "Comparison of R_work between twinned and non-twinned entry"
    xrange, yrange, xlabel, ylabel = "", "", "resolution", "R_work"
    bar, rot, key, style = 0, 1, 2, 1
    plot = """plot '%s' using 4:xtic(1) t "twinned", '%s' u 4 t "non-twinned"  \
        """ % (
        fcorr,
        fcorr_reg,
    )
    file = ddir + "twin-nontwin-rwork.data"
    _gnuscr, _gnuout = gnu_plot(file, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)


######################################################################
def uniq_list_of_list(data):
    """get uniq list"""
    data_new = [x for x in data if x not in locals()["_[1]"]]
    # data_new=[]
    # for x in data:
    #    if x not in data_new:  data_new.append(x)

    return data_new


######################################################################
def rfact_reso_growth(ddir):
    """test the mean resolution change with year"""

    out = "%srfact_reso_growth.all" % ddir

    plot_data = out + ".data"  # pylint: disable=redefined-outer-name
    plot_data1 = out + "_res.data"
    fw = open(plot_data, "w")
    fw1 = open(plot_data1, "w")

    arg = """%s " \
select distinct  d.structure_id, d.date_original , r.ls_d_res_high ,  \
r.ls_R_factor_R_work from database_PDB_rev as d , refine as r, pdb_entry as p \
where  d.structure_id = r.structure_id and r.ls_R_factor_R_work != '' \
and r.structure_id =p.pdb_id  and p.method like '%%x-ray%%' \
and p.status_code='REL' and p.method !='THEORETICAL MODEL'  \
order by d.date_original,  r.ls_d_res_high \
">%s
""" % (
        SQL,
        out,
    )
    os.system(arg)
    data1 = from_file_to_list(out)
    dic = {}

    reso = [0, 1.5, 2, 2.5, 3.0, 5.0]

    s = "#year  (0.0->1.5) dev (1.5->2.0) dev (2.0->2.5) dev (2.5->3.0) dev (3.0->4.0) dev "
    print(s)

    fw.write(s + "\n")

    for i in range(2000, YEAR):
        n = 0
        rf = []
        for x in data1[:]:
            if "structure_id" in x:
                continue
            if str(i) in x[1]:
                n = n + 1
                rf.append([float(x[3]), float(x[4])])
            elif str(i) not in x[1] and n > 0:
                break

        nrf = len(rf)

        tmp = [".", ".  . ", ".  .", ".  . ", ".  .", ".  ."]
        tmp1 = [0, 0, 0, 0, 0, 0]
        if rf:
            for k, _z in enumerate(reso):
                if k == 0:
                    continue
                t1 = []
                ns = 0
                for y in rf:
                    if reso[k - 1] <= y[0] < reso[k]:
                        t1.append([y[1]])
                        ns = ns + 1

                avg, dev = -1.0, -1.0
                if t1:
                    avg, dev, _mini, _maxi = util.mean_dev(t1, 0)
                    tmp[k] = "%.3f %.3f " % (avg, dev)

                    p1 = 100 * float(ns) / nrf
                    tmp1[k] = "%.2f " % (p1)

                    # s1 = "(%.2f->%.2f)  %.3f %.3f  %d" % (reso[k - 1], reso[k], avg, dev, len(t1))
                    # s2 = "%d  (%.2f->%.2f) %d  %d  %.2f" % (i, reso[k - 1], reso[k], ns, nrf, p1)
                    # print(s2)

            #  print(k, z, t1,avg, dev)
            ss = "%d %s  %s %s %s %s " % (i, tmp[1], tmp[2], tmp[3], tmp[4], tmp[5])
            print(ss)
            fw.write(ss + "\n")
            #            ss1='%d %.2f %.2f %.2f %.2f %.2f \n' %(i, tmp1[1], tmp1[2],tmp1[3],tmp1[4],tmp1[5])
            ss1 = "%d %s  %s %s %s %s \n" % (i, tmp1[1], tmp1[2], tmp1[3], tmp1[4], tmp1[5])
            fw1.write(ss1)
        dic[i] = rf

    fw.close()
    fw1.close()

    title = "The mean Rfactors change with year"
    xrange, yrange, xlabel, ylabel = "", "", "Year", "R_work"
    bar, rot, key, style = 0, 1, 2, 1
    plot = """plot '%s' using 2:xtic(1) t "reso: 0.0->1.5" , \
        '' u 4 t "reso: 1.5->2.0" , '' u 6 t "reso: 2.0->2.5", \
        '' u 8 t "reso: 2.5->3.0" , '' u 10 t "reso: 3.0->4.0"  \
        """ % (
        plot_data
    )

    _gnuscr, _gnuout = gnu_plot(plot_data, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)

    title = "The mean growth rate of resolution in each year"
    xrange, yrange, xlabel, ylabel = "", "", "Year", "percentage of resolution"
    bar, rot, key, style = 0, 1, 2, 1
    plot = """plot '%s' using 2:xtic(1) t "reso: 0.0->1.5" , \
        '' u 3 t "reso: 1.5->2.0" , '' u 4 t "reso: 2.0->2.5", \
        '' u 5 t "reso: 2.5->3.0" , '' u 6 t "reso: 3.0->4.0"  \
        """ % (
        plot_data1
    )

    _gnuscr, _gnuout = gnu_plot(plot_data1, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)


######################################################################
def get_items_str():
    """The format of the items must be fixed for the data base.
    col_1, col_2: the cif category, the cif item
    col_3: if 'x-ray', only get XRAY; if 'all', include all method
    col_4,5,6,7,8: if given, manual contral the first five strings
    """
    items = [
        ["start:string", "", "", ""],
        ["computing", "pdbx_data_reduction_ii", "x-ray", "DENZO", "HKL", "MOSFLM", "XDS/XSCALE", "D*TREK/DTREK"],
        ["computing", "pdbx_data_reduction_ds", "x-ray", "SCALEPACK", "HKL", "SCALA", "XDS/XSCALE", "D*TREK/DTREK"],
        ["computing", "structure_refinement", "x-ray", "REFMAC", "PHENIX", "BUSTER", "SHELX", "CNS/X-PLOR"],
        ["computing", "structure_solution", "x-ray", "PHASER", "CNS/X-PLOR", "MOLREP", "SOLVE", "AMORE"],
        ["computing", "pdbx_structure_refinement_method", "x-ray", "", "", "", "", ""],
        [
            "refine",
            "pdbx_method_to_determine_struct",
            "x-ray",
            "MOLECULAR_REPLACEMENT/MR",
            "SAD/SIRAS",
            "MAD/MIRAS",
            "MIR/MULTIPLE_ISOMORPHOUS_REPLACEMENT",
            "AB_INITIO/DIRECT_METHOD",
        ],
        ["diffrn_detector", "detector", "x-ray", "CCD", "IMAGE_PLATE", "AREA_DETECTOR", "PIXEL", "FILM"],
        ["diffrn_detector", "type", "x-ray", "MARRESEARCH/MAR", "ADSC", "RIGAKU", "SIEMENS", "PSI_PILATUS"],
        ["diffrn_source", "source", "x-ray", "", "", "", "", ""],
        ["diffrn_source", "type", "x-ray", "", "", "", "", ""],
        ["symmetry", "space_group_name_H_M", "x-ray", "", "", "", "", ""],
        ["exptl", "method", "all", "", "", "", "", ""],
        ["end:string", "", "", ""],
    ]
    return items


######################################################################
def do_string_stat(ddir):
    """The function is to do populations for string
    (date is from database_PDB_rev   fixed)
    """

    pdbrev = ["database_PDB_rev", "structure_id", "date_original"]
    items = get_items_str()

    for x in items:
        print(x)
        if "start:" in x[0] or "end:" in x[0]:
            continue
        tabs = [pdbrev, [x[0], x[1]]]
        outf = "%s%s_%s.all" % (ddir, x[0], x[1])
        do_sql_gen(outf, tabs, x[2])
        first5, file1, file2 = string_stat(outf, x)

        title = "The population of %s.%s" % (x[0], x[1])
        xrange, yrange, xlabel, ylabel = "", "", "Year", "PDB entry"
        bar, rot, key, style = 0, 1, 0, 0
        plot = """plot '%s' using 2:xtic(1) lc rgb "green" , \
        '' u 0:2:3 with labels offset 0, 0.5 """ % (
            file1
        )

        _gnuscr, _gnuout = gnu_plot(file1, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)

        title = "Growth of %s.%s" % (x[0], x[1])
        xrange, yrange, xlabel, ylabel = "", "", "Year", "PDB entry"
        bar, rot, key, style = 0, 1, 2, 0
        plot = """plot '%s' using 2:xtic(1) t "%s" ,'' u 4 t "%s" ,'' u 6 t "%s"  \
        ,'' u 8 t "%s" , '' u 10 t "%s" \
        """ % (
            file2,
            first5[0],
            first5[1],
            first5[2],
            first5[3],
            first5[4],
        )
        _gnuscr, _gnuout = gnu_plot(file2, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)

        title = "Growth rate of %s.%s" % (x[0], x[1])
        xrange, yrange, xlabel, ylabel = "", "", "Year", "percentage (%)"
        bar, rot, key, style = 0, 1, 2, 0
        plot = """plot '%s' using 3:xtic(1) t "%s" ,'' u 5 t "%s" ,'' u 7 t "%s"  \
        ,'' u 9 t "%s" ,'' u 11 t "%s" \
        """ % (
            file2,
            first5[0],
            first5[1],
            first5[2],
            first5[3],
            first5[4],
        )
        file3 = file2 + "_rate"
        _gnuscr, _gnuout = gnu_plot(file3, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)


######################################################################
def uniq_string_pop(data1, col):
    """get unique string populations
    data: a list,  col: the column that contains the string
    """

    data1.sort(key=lambda y: y[col])
    #    print(col, data1)
    nnt = len(data1)  # accumulate the unique strings
    n, data2 = 0, []
    for i, x in enumerate(data1):
        if i == nnt - 1 or (i < nnt - 1 and data1[i][col] != data1[i + 1][col]):
            pn = 100 * float(n) / nnt
            data2.append([x[col], n, pn])
            n = 0
        n = n + 1
    data2.sort(key=lambda y: y[2])

    return data2


######################################################################
def string_stat(outf, prog):
    """This is a general stat for string (plot population in order), growth and
    growth rate.
    outf contains more than 3 columns (id, date1, 0000 ..)
    prog is the list (column 4,5,6,7,8 contains programs)

    """

    fp = open(outf, "r")
    data1 = []  # remove null values
    for x in fp:
        t = x.split()
        if "structure_id" in x or len(t) <= 3:
            continue
        s = "_".join(t[3:])
        data1.append([t[1][:4], s.upper()])
    fp.close()
    data1.sort(key=lambda y: y[1])

    data2 = uniq_string_pop(data1, 1)  # accumulate the unique strings
    data2.reverse()

    plot_data1 = outf + "_pop.data"
    fw = open(plot_data1, "w")
    fw.write("# unique-intem  number-of-entry  percentage (n/n_total_per_year)\n")
    for i, x in enumerate(data2):
        s = "%20s  %d  %.1f \n" % (x[0], x[1], x[2])
        fw.write(s)
        print(s.strip())
        if i > 25:
            break
    fw.close()

    first5 = []  # plot the first top 5 string (programs)
    if len(prog[3]) and len(prog[4]):
        first5 = [prog[3], prog[4], prog[5], prog[6], prog[7]]
    else:
        for i in range(5):
            first5.append(data2[i][0])
    print("using %s for yearly growth" % first5)

    plot_data2 = outf + ".data"
    fw = open(plot_data2, "w")
    s = "#Year  " + "    ".join([" %s rate " % x for x in first5])
    print(s)
    fw.write(s + "\n")

    for i in range(1994, YEAR):
        num = []
        for y in first5:
            n, pn, ny = 0, 0, 0
            for x in data1:
                if i == int(x[0]):
                    ny = ny + 1
                    if len(y) > 2:
                        if "/" in y:
                            tt1 = y.split("/")
                            if tt1[0] in x[1] or tt1[1] in x[1]:
                                n = n + 1
                        else:
                            if y in x[1]:
                                n = n + 1

            #                    if y==x[1] :

            if ny > 0:
                pn = 100 * float(n) / ny
            num.append(n)
            num.append("%.1f" % pn)
            # print(i, y, n, ny)

        ss = "%d  " % i + "   ".join([str(x) for x in num])
        print(ss)
        fw.write(ss + "\n")

    fw.close()

    return first5, plot_data1, plot_data2


######################################################################
def do_reso_entity(ddir):
    xray_id = "select pdb_id from pdb_entry  where \
    method like '%x-ray%' and status_code='REL' and method !='THEORETICAL MODEL'"

    four = """%s "select  r.structure_id,  r.ls_d_res_high from refine as r where  \
r.structure_id  in ( select structure_id from entity_poly where lcase(type) like '%%polypeptide%%' ) \
and r.structure_id not in (select structure_id from entity_poly where lcase(type) like '%%polyribo%%' ) \
and r.structure_id not in (select structure_id from entity_poly where lcase(type) like '%%polydeoxy%%' ) \
and r.structure_id in (%s) order by r.ls_d_res_high \
">protein_reso.list


%s "select  r.structure_id,  r.ls_d_res_high from refine as r where  \
r.structure_id  in ( select structure_id from entity_poly where lcase(type) like '%%polyribo%%' ) \
and r.structure_id not in (select structure_id from entity_poly where lcase(type) like '%%polypeptide%%' ) \
and r.structure_id not in (select structure_id from entity_poly where lcase(type) like '%%polydeoxy%%' ) \
and r.structure_id in (%s) order by r.ls_d_res_high \
">rna_reso.list


%s "select  r.structure_id,  r.ls_d_res_high from refine as r where  \
r.structure_id  in ( select structure_id from entity_poly where lcase(type) like '%%polydeoxy%%' ) \
and r.structure_id not in (select structure_id from entity_poly where lcase(type) like '%%polypeptide%%' ) \
and r.structure_id not in (select structure_id from entity_poly where lcase(type) like '%%polyribo%%' ) \
and r.structure_id in (%s) order by r.ls_d_res_high \
">dna_reso.list


%s "select  r.structure_id,  r.ls_d_res_high from refine as r where  \
(r.structure_id  in ( select structure_id from entity_poly where lcase(type) like '%%polydeoxy%%' ) \
or r.structure_id in (select structure_id from entity_poly where lcase(type) like '%%polyribo%%' )) \
and r.structure_id in ( select structure_id from entity_poly where lcase(type) like '%%polypeptide%%' ) \
and r.structure_id in (%s) order by r.ls_d_res_high \
">protein_na_reso.list


""" % (
        SQL,
        xray_id,
        SQL,
        xray_id,
        SQL,
        xray_id,
        SQL,
        xray_id,
    )

    print(four)

    os.system(four)

    shell = data_bin(0.7, 4.8, 13)  # data change with resolution
    four = ["protein_reso.list", "protein_na_reso.list", "rna_reso.list", "dna_reso.list"]
    s1 = " ".join(four)
    os.system("cat %s > %s/resolution_entity.data" % (s1, ddir))

    res, dv = [], []
    col = 1  # data is in this column
    for x in four:
        data_clean = clean_data(x, col)
        data, _outlier = filter_data(data_clean, shell, col)
        fpop = x + "_pop.data"  # for data populations
        data_pop(fpop, data, shell, col)
        avg, dev, _mini, _maxi = util.mean_dev(data, col)
        #        print(avg, dev, mini, maxi)

        res.append(avg)
        dv.append(dev)

    title = "Population of resolution with different entity"
    xrange, yrange, xlabel, ylabel = "", "", "resolution", "percentage"
    bar, rot, key, style = 0, 1, 1, 1
    plot = """plot '%s_pop.data' using 4:xtic(1) t "protein (mean=%.2f, dev=%.2f)" , \
    '%s_pop.data' u 4:xtic(1) t "protein-NA (mean=%.2f, dev=%.2f)" ,  \
    '%s_pop.data' u 4:xtic(1) t "RNA (mean=%.2f, dev=%.2f)", \
    '%s_pop.data' u 4:xtic(1) t "DNA (mean=%.2f, dev=%.2f)" """ % (
        four[0],
        res[0],
        dv[0],
        four[1],
        res[1],
        dv[1],
        four[2],
        res[2],
        dv[2],
        four[3],
        res[3],
        dv[3],
    )

    fpop = ddir + "resolution_entity"

    _gnuscr, _gnuout = gnu_plot(fpop, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)


######################################################################
def fast_sort(file_in, k1, k2):
    """fast sort; k1 column for string, k2 column for number."""
    out = file_in + ".sorted"
    arg = "sort -k %d,%d -k %d,%dn %s > %s" % (k1, k1, k2, k2, file_in, out)
    os.system(arg)
    return out


######################################################################
def sort_multi_column(file_in):
    file_in = "tt"
    fp = open(file_in, "r").readlines()
    fp.sort(key=lambda y: y[0])
    print("puting it to dic")
    dic = {}
    for x in fp:
        t = x.split()
        if t[0] not in dic.keys():
            dic[t[0]] = []
        if t[0] in dic.keys():
            dic[t[0]].append([t[1], t[2]])

    print("printting")

    keys = list(dic.keys())
    keys.sort()
    for x in keys:
        dic[x].sort(key=lambda y: y[1])
        for m in dic[x]:
            print("%s %s" % (x, m))


######################################################################
def pdb_growth_fit(ddir, idd):
    """fit the data by equation, use the database."""

    yearly_dep, accumu_dep = pdb_growth(1976, YEAR, "DEP")
    fit_growth(ddir, idd, accumu_dep, "deposited")
    fit_growth(ddir, idd, yearly_dep, "deposited")

    yearly_rel, accumu_rel = pdb_growth(1976, YEAR, "REL")
    fit_growth(ddir, idd, accumu_rel, "released")
    fit_growth(ddir, idd, yearly_rel, "released")

    util.delete_file("pdb_growth.*__tmp* pdb_growth.*_DEP  pdb_growth.*_REL ")
    os.system("mv pdb_growth.txt_year_*.* %s" % (ddir))
    os.system("mv pdb_growth.txt_acum_*.* %s" % (ddir))


######################################################
def fit_growth(ddir, idd, file, type):  # pylint: disable=unused-argument,redefined-builtin
    """fit/plot data (from  2000 to 2011)
    ddir: the directory to hold graphs and data
    file: the data files
    type: is for released or for deposited.
    """

    x0_fit, xf_fit = 2000, 2013  # used to fit data. change them,if needed.

    fp = open(file, "r").readlines()
    d = [x for x in fp if len(x) > 1 and "year" not in x]
    # x0 = int(d[0].split()[0])  # first year
    xf = int(d[-1].split()[0])  # last year

    a, b, c, eq = fitted(file, x0_fit, xf_fit, idd)
    predict_year = 7
    pout = predict(file, a, b, c, x0_fit, xf, predict_year, idd)

    title = "Actual growth (%s) and prediction (by %s)." % (type, eq)
    xrange, yrange, xlabel, ylabel = "", "", "Year", "PDB entry"
    bar, rot, key, style = 0, 1, 1, 0

    plot = """plot '%s' using 1:3:xtic(1) t "predicted" with linespoints ,'' u 1:2 t "%s"  with boxes""" % (pout, type)
    _gnuscr, _gnuout = gnu_plot(pout, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)


######################################################
def fitted(datain, x0, xf, idd):
    """fit the equation by the data (two columns; 1: year, 2:data)
    idd=1: fit (a,b,c) for Polynomial  f(x) = a*(x-x0)**b + c
    idd=2: fit (a,b,c) for Exponential f(x) = a*exp(b*(x-x0)) + c
    x0, xf: the start & end data range used for fitting
    """

    data = datain + "__tmp"
    fw = open(data, "w")
    for x in open(datain, "r").readlines():
        t = x.split()
        if int(t[0]) < x0 or int(t[0]) > xf:
            continue
        fw.write(x)
    fw.close()

    fit = data + "_fit"
    fw = open(fit, "w")
    if idd == 1:
        fw.write("f(x)=a*((x-%d)**b) + c \n" % x0)
        fw.write("a=105; b=1.5;  c=1000; \n")  # initial values (adjust)
        fw.write('fit f(x) "%s" u 1:2 via a, b, c \n' % data)
    elif idd == 2:
        fw.write("f(x)=a*exp(b*(x-%d)) + c \n" % x0)
        fw.write("a=303.995; b= 0.51;  c=-100.165463;  \n")
        fw.write('fit f(x) "%s" u 1:2 via a, b, c \n' % data)

    elif idd == 3:
        fw.write("f(x)=a*exp(b*(x-%d)) \n" % x0)
        fw.write("a=303.995; b= 0.51;  \n")
        fw.write('fit f(x) "%s" u 1:2 via a, b \n' % data)

    elif idd == 4:
        fw.write("f(x)=a*exp(b*(x)) \n")
        fw.write("a=303.995; b= 0.51;  \n")
        fw.write('fit f(x) "%s" u 1:2 via a, b \n' % data)

    fw.close()

    log = data + "_fit.log"
    arg = "gnuplot  %s >& %s" % (fit, log)
    os.system(arg)

    util.delete_file("fit.log")
    #    print(x0,xf,d)

    a, b, c, da, db, dc = 0, 0, 0, 0, 0, 0
    fp = open(log, "r").readlines()  # get fitted values from log
    for x in fp:
        if "a" in x[:3] and "+/-" in x[30:45] and "(" in x[49:]:
            t = x.split()
            a, da = float(t[2]), float(t[4])
        elif "b" in x[:3] and "+/-" in x[30:45] and "(" in x[49:]:
            t = x.split()
            b, db = float(t[2]), float(t[4])
        elif "c" in x[:3] and "+/-" in x[30:45] and "(" in x[49:]:
            t = x.split()
            c, dc = float(t[2]), float(t[4])

    print("fitted values using data (%d->%d) a=%.4f  da=%.4f, b=%.4f db=%.4f, c=%.4f dc=%.4f" % (x0, xf, a, da, b, db, c, dc))
    eq = "f(x)=%.2f*((x-%d)**%.2f) + %.2f" % (a, x0, b, c)
    if idd == 2:
        eq = "f(x)=%.2f*exp(%.2f*(x-%d)) + %.2f" % (a, b, x0, c)
    return a, b, c, eq


############################################################
def predict(data, a, b, c, x0, xf, pred, idd):
    """predict the growth in next predict years
    data: the file containing two columns (year, number)
    a, b, c : the fitted parameters
    x0, xf : the start and end year for the data fitting
    pred : number of years ahead of the xf (predicting)
    idd : 1 for quadratic fitting, otherwise for exponential fitting
    """

    fr = open(data, "r").readlines()

    endy = xf + pred
    pout = data + ".data"
    fw = open(pout, "w")

    print("Year  deposited  Predicted  Accuracy(%)")
    fw.write("#Year  deposited  Predicted  Accuracy(%)\n")

    t1 = 0
    for x in fr:
        if "year" in x:
            continue
        t = x.strip().split()
        t1, t2 = int(t[0]), int(t[1])
        if t1 < x0:
            print("%4d %6d   ?      ? " % (t1, t2))
            fw.write("%4d %6d   ?      ? \n" % (t1, t2))
        else:
            y = fiteq(a, b, c, t1, x0, idd)
            percent = 100 * math.fabs(t2 - y) / t2
            print("%4d %6d %6d   %.1f  " % (t1, t2, y, percent))
            fw.write("%4d %6d %6d   %.1f  \n" % (t1, t2, y, percent))
    if endy >= t1 + 1:
        for x in range(t1 + 1, endy):
            y = fiteq(a, b, c, x, x0, idd)
            print("%4d     ?  %6d   ?  " % (x, y))
            fw.write("%4d     ?  %6d   ?  \n" % (x, y))

    fw.close()

    return pout


############################################################
def fiteq(a, b, c, x, x0, idd):
    """fited equation"""

    y = 0
    if idd == 1:
        y = a * ((x - x0) ** b) + c
    elif idd == 2:
        y = a * math.exp(b * (x - x0)) + c

    elif idd == 3:
        y = a * math.exp(b * (x - x0))

    elif idd == 4:
        y = a * math.exp(b * (x))

    return y


############################################################
def pdb_growth(start, end, type):  # pylint: disable=redefined-builtin
    """get growth of PDB using sql"""

    out = "pdb_growth.txt"

    if "REL" in type:  # only for release
        scr = """%s \
"select year(r.date), count(d.pdb_id) from pdb_entry_tmp d, database_PDB_rev r \
where r.structure_id = d.pdb_id and r.date >='%d' and r.date <= '%d' \
and r.num=1  and d.method!='theoretical model' and d.status_code ='REL'  \
group by year(r.date)" >%s \
""" % (
            SQL,
            start,
            end,
            out,
        )

    else:  # only for deposited
        scr = """%s \
"select year(initial_deposition_date),count(distinct structure_id) \
from pdb_entry_tmp where author_release_status_code not like '%%WDRN%%' and \
year(initial_deposition_date)<= %d group by year(initial_deposition_date); " >%s \
""" % (
            SQL,
            YEAR,
            out,
        )

    print("%s\n" % scr)
    os.system(scr)

    out_acum = out + "_acum_" + type
    out_year = out + "_year_" + type
    fp = open(out, "r")
    fw1 = open(out_acum, "w")
    fw2 = open(out_year, "w")
    n = 0
    for x in fp:
        if "year" in x.lower():
            continue
        t = x.split()
        year, val = int(t[0]), int(t[1])
        n = n + val
        fw1.write("%4d  %6d\n" % (year, n))
        fw2.write("%4d  %6d\n" % (year, val))

    fp.close()
    fw1.close()
    fw2.close()
    #    util.delete_file(out)
    return out_year, out_acum


############################################################
def gnu_plot(datafile, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot):
    """some head infor for gnuplot (one data set plot by default)
    bar>0: add error bar on the histogram
    rot>0, xlabel rotation;
    key=0: key>0, on;  key=2,left;
    style=0, histogram ; =1 linepoints; =1, points; =3 dot

    """

    gnuscr = datafile + ".gnu"
    gnuout = datafile + ".png"
    fw = open(gnuscr, "w")

    #    x1='set terminal  png large size 1000,740 # 840,640\n'
    x1 = "set terminal  gif large size 1000,740 # 840,640\n"
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
    fw.write("# plot 'filename' using 1:2:3:xtic(1) t \"all\" with  yerrorlines (error bar)\n")

    #    fw.write("set term x11\nreplot \npause 4\n")

    fw.close()
    #    os.system('/bin/csh -c "gnuplot %s" ' %gnuscr) #if not given /bin/csh -c, use B shell
    #    os.system('gnuplot %s ' %gnuscr) #if not given /bin/csh -c, use B shell
    os.system("gnuplot %s " % gnuscr)  # if not given /bin/csh -c, use B shell

    return gnuscr, gnuout


######################################################
def plot_data_fes():
    """only for tempary use"""

    ddir = "DATA_PLOT/"  # create a dir to hold all the data/graphs
    os.system("mkdir %s" % ddir)  # make a fold to hold all the data

    fp = open("Fe-Fe.log", "r").readlines()
    data_in = []
    for x in fp:
        data_in.append(x.split())

    # item, xlabel, ylabel,data_type, min, max, steps
    y = ["FE-FE", "Bond distance (Fe-Fe)", "percentage", 1, 1, 5.0, 20]
    plot_pop_new(ddir, y, data_in, 1)

    fp = open("Fe-noS-noFe.log", "r").readlines()
    data_in = []
    for x in fp:
        data_in.append(x.split())

    y = ["FE-other", "Bond distance (Fe-other)", "percentage", 1, 1, 5.0, 20]
    plot_pop_new(ddir, y, data_in, 1)

    fp = open("Fe-S.log", "r").readlines()
    data_in = []
    for x in fp:
        data_in.append(x.split())

    y = ["FE-S", "Bond distance (Fe-S)", "percentage", 1, 1, 5.0, 20]
    plot_pop_new(ddir, y, data_in, 1)


############################################################


if __name__ == "__main__":
    proc_stat(sys.argv)
#    get_stat()
#    plot_data_fes()
#    plot_data()
#    predict()
