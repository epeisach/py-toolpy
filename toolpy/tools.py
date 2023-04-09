#!/usr/bin/env python

import os, sys, math, shutil
from time import time
import prog, sf_util, util, matt, contact, ref_util, statis
import tlswater as tls, cifparse as cif

# import test_cython


def usage():
    content = """
    
###############################################################################
usage: tool  [option] file

This utility tool is to perform various jobs (created 2010-12-03)
    
*. tool -stat filename  xcol=?, ycol=?
      optional options:  xlabel=?, ylabel=?, xrange=? ?, yrange=? ?, title=?
      calculate statistics (min,max,mean,dev,skew, Kurtosis) and plot graphs.
      File must have columns of data.
         a). For data distribution: xcol=? [1,2,3, ..] must be given.
         b). For data correlations: xcol=?, ycol=? must be given.

*. tool -changeB pdbfile  const (or b11 b22 b33 b12 b13 b23)
       add a constant to Biso or (ANISOU) to pdb. (refmac V5.2 make B<0 with TLS) 
        
*. tool -contact pdbfile (id)  #calculate crystal contact in PDB (by ccp4)
       (id=0, use crystal frame (remove SCALE); id=1, use SCALE in PDB) 
    
*. tool -res  xyzfile  #reorder residue number sequentially
       if add A4_B5 after pdbfile, chain A 4 will be started from chain B 5

*. tool -pick file1 file2  #grep items from first column of file1 and 
      extract the lines in file2 containing the item.

*. tool -mr285  xyzfile  sffile  # (07/15/2011)
      coord. not in the crystal frame. Use MR to get new xyz and align the 
      original xyz to the new_xyz to get Matrix(X0). Then generate Remark 285.
    
*. tool -285  pdbfile1 pdbfile2  # (07/15/2011)
      coord. shift off crystal frame. pre_generate (rigid body + refine) pdbfile2
      Align pdbfile1 to pdbfile2 to get Matrix(X0). Then generate Remark 285.
    
*. tool -twin  sffile  ( -xtriage )  #(06/12/2012)
      analyze truncate results and make plots for twin/intesity/Wilson.
     (Add pdbfile after -twin , if no cell in the sffile)
    
*. tool -calc_sym  sffile    #(07/06/2012)
      use pointless to find best space group. (Add pdbfile, if no cell in sffile)

*. tool -tlswater  pdbfile  sffile   #(01/24/2013)
      Compare R/Rfree with/without waters(mainly for refmac5.5, waters added as
      default,but missing residue ranges).
      The pdbfile must be the deposited file with partial B and without water range.

*. tool -superpose  xyzfile -s "A" file1 -s "B" file2 ... -o output   #(2017-3)
       multipile structure alignment by superpose. If no output, default superpose.out

*. tool -mr prog=PROG, pdb=pdbfile, sf=sffile, nmol=1    #(2017-6-31)
       Do MR (program=epmr|phaser|molrep   , nmol=numb of molecule) 
       

*. tool -blob_map  pdbfile  map   #find blob of density from the map.
*. tool -anis pdbfile  #Analyze anisotropy of each atom. ANISOUs must exist.
*. tool -cont xyzfile #calculate atom contacts (by grid method) (2013-08).    
*. tool -ha  pdbfile       #remove H atoms,  (07/16/2012)
*. tool -mapdmp  mapfile   #look at the header of a map.(10/18/2012)
*. tool -rfree  pdb/pdbx  mtzfile/sfcif  #test free set.( 2016-10-10)
*. tool -swap_fi pdbfile sffile  #Test R/Rfree after swap F->I / I->F(2013-03)
*. tool -occ  pdb/pdbx     #(find/correct occupancy of atom on special pos.)
*. tool -occ_ph  pdb/pdbx  #(find/correct occupancy, use phenix)

*. tool -scale xyzfile # calculate scale matrix in pdb/pdbx file 
*. tool -matt pdb/pdbx  #calculate Matthew coefficient.(04/18/2013)
*. tool -pdf pdbid      #get wwPDB validation pdf full report (09/18/2014).
*. tool -xml pdbid      #get wwPDB validation xml report (09/18/2014).
*. tool -sf  pdbid      #get SF file from database (2017-6-23)
*. tool -pdb pdbid      #get PDB file from database (2017-6-23)
*. tool -cif pdbid      #get PDBX file from database (2017-6-23)
*. tool -pdbid pdbid    #get SF/PDB  files from database (2017-6-23)
*. tool -cifid pdbid    #get SF/PDBX files from database (2017-6-23)
*. tool -alt pdbfile    #Assign proper ALT id (by grid method) (09/10/2013) .
*. tool -rm_alt pdbfile  #remove the alt conformers (10/09/2014).
*. tool -rm_water pdbfile sffile  #remove water for bad RsR (10/09/2014).
*. tool -add_water mtzfile pdbfile  #add waters to pdbfile (10/09/2014).
*. tool -ext_pdbredo pdbid  #get fully optimised pdb/mtz file (03/02/2015).
*. tool -shift xyzfile   #shift xyz file to unit cell (03/31/2015).
*. tool -wwpdf xyzfile sffile outfile  #generate wwpdf report(API)(04/10/2017).
*. tool -dic file  #generate dic from a formated file (10/19/2015).
*. tool -clean file  #clean pdb file (add atom symbol to ANIS) (10/20/2015)
*. tool -map xyzfile sffile  #display map use FWT/PHWT,DELFWT,DELPHWT(2016-11)
*. tool -add_percentile file column_num  #add percentile to column(2017-5-22)
*. tool -geom pdbfile # use phenix to evaluate model based on molprobity,
###############################################################################

"""

    print(content)
    sys.exit()


###############################################################################


def process(*files):
    arg = files[0]
    if len(arg) < 2:
        usage()

    start_time = time()

    narg = len(arg)

    outfile = ""
    out_xyz = ""
    for k in range(narg):
        if arg[k].lower() == "-o":
            outfile = arg[k + 1]
        elif arg[k].lower() == "-oxyz":
            out_xyz = arg[k + 1]

    for k in range(narg):
        if arg[k].lower() == "-list":
            flist = arg[k + 1]

        elif arg[k].lower() == "-anis":
            ani = anisotropy(arg[k + 1])

        elif arg[k].lower() == "-stat":  # cal dev,mean,outliers, plot frequncy.
            statis.get_pop_corr(arg[k + 1], k + 2, narg, arg)
            sys.exit()

        elif arg[k].lower() == "-mr":  # Do MR by EPMR/Morep/Phaser
            prog.do_mr(k + 1, narg, arg)

        elif arg[k].lower() == "-add_percentile":  # add percentile rank for column #
            if narg - k != 3:
                print("Input error! Not enough argument (tool -add_perc file column).")
                return
            file, column = arg[k + 1], int(arg[k + 2])
            statis.add_percentile_rank(file, column)  #
            sys.exit()

        elif arg[k].lower() == "-scale":
            ani = calc_scale_matrix(arg[k + 1])

        elif arg[k].lower() == "-contact":
            id = 0  # using crystal frame
            if len(arg) == 4:
                id = int(arg[k + 2])
            ani = get_contact(arg[k + 1], id)

        elif arg[k].lower() == "-dev":
            min, max, mean, std, num = mean_std("", arg[k + 1], int(arg[k + 2]))

        elif arg[k].lower() == "-bin":
            t1, t2, t3, t4 = arg[k + 1], int(arg[k + 2]), int(arg[k + 3]), int(arg[k + 4])
            data = mean_std_bin("", t1, t2, t3, t4)

        elif arg[k].lower() == "-changeb":
            if len(arg) > 6:
                t = [float(arg[k + i + 2]) for i in range(6)]
            else:
                t = [float(arg[k + 2])]
            newpdb = change_bfactor(arg[k + 1], t)

        elif arg[k].lower() == "-res":  # rename residue number from start
            if len(arg) == 3:  # automatic
                newpdb = residue_num_seq(arg[k + 1])
            elif len(arg) > 3:
                newpdb = residue_num_seq(arg[k + 1], arg[k + 2])

        elif arg[k].lower() == "-pick":
            pickup_item_from_firstfile(arg[k + 1], arg[k + 2])

        elif arg[k].lower() == "-mr285":  # generate remark 285 & validate, by MR
            get_285_by_mr(arg[k + 1], arg[k + 2], 0)

        elif arg[k].lower() == "-285":  # generate remark 285 & validate
            get_285_by_mr(arg[k + 1], arg[k + 2], 1)

        elif arg[k].lower() == "-occ":  # find/correct occ
            contact.check_special_position(arg[k + 1], outfile, out_xyz)
        #            test_cython.check_special_position(arg[k+1], outfile)

        elif arg[k].lower() == "-occ_ph":  # find/correct occ, by phenix
            correct_occ_on_symm(arg[k + 1])

        elif arg[k].lower() == "-twin":  # analysis twin
            if k + 2 == narg:
                sf_util.sf_quality(arg[k + 1], "")
            elif k + 2 < narg:
                sf_util.sf_quality(arg[k + 1], arg[k + 2])

        elif arg[k].lower() == "-calc_sym":  # get space group by pointless
            if k + 2 == narg:
                sf_util.sf_symmetry(arg[k + 1], "")
            elif k + 2 < narg:
                sf_util.sf_symmetry(arg[k + 1], arg[k + 2])

        elif arg[k].lower() == "-ha":  # remove H atoms
            remove_h_atom(arg[k + 1])

        elif arg[k].lower() == "-mapdmp":  # header of a map
            map_header(arg[k + 1])

        elif arg[k].lower() == "-rfree" or arg[k].lower() == "-freer":  # test differ free sets
            sf_util.get_freer_flag(arg[k + 1], arg[k + 2])  # input pdb/cif & mtz/cif

        elif arg[k].lower() == "-swap_fi":  # test R/Rf after F->I
            sf_util.test_swap_fi(arg[k + 1], arg[k + 2])

        elif arg[k].lower() == "-tlswater":  # test tls with/without water
            tls.compare_r_tls_water(arg[k + 1], arg[k + 2])

        elif arg[k].lower() == "-matt":  # calculate Matthew coeffi.
            matt.matt_coeff(arg[k + 1], outfile)

        elif arg[k].lower() == "-cont":  # calculate contact by grid
            contact.contact(arg[k + 1])

        elif arg[k].lower() == "-alt":  # assign proper ALT id by grid
            contact.assign_alt_id(arg[k + 1], outfile)

        elif arg[k].lower() == "-pdf":  # follow pdbid, full report
            get_vtf_report(arg[k + 1], "fpdf")

        elif arg[k].lower() == "-xml":  # follow pdbid
            get_vtf_report(arg[k + 1], "xml")

        elif arg[k].lower() == "-wwpdf":  # follow xyz sf
            gen_vtf_report(arg[k + 1], arg[k + 2], arg[k + 3])

        elif arg[k].lower() == "-ext_pdbredo":  # follow pdbid
            get_pdbredo_report(arg[k + 1])

        elif arg[k].lower() == "-rm_alt":  # follow pdbid
            args = []
            for n in range(k + 1, narg):
                args.append(arg[n])

            ref_util.remove_alt_conformer(args)

        elif arg[k].lower() == "-rm_water":  # remove water if bad RsR
            ref_util.remove_bad_water(arg[k + 1], arg[k + 2])  # pdb, sf

        elif arg[k].lower() == "-add_water":  # remove water if bad RsR
            ref_util.add_water(arg[k + 1], arg[k + 2], 2)  # mtz, pdb

        elif arg[k].lower() == "-plot":  # plot histogram
            plot_file(arg[k + 1])

        elif arg[k].lower() == "-shift":  # plot histogram
            contact.shift_xyz_into_unitcell(arg[k + 1], "")

        elif arg[k].lower() == "-dic":  # make cif dictionary
            make_cif_dic(arg[k + 1])

        elif arg[k].lower() == "-clean":
            ref_util.clean_file_4d3r(arg[k + 1])

        elif arg[k].lower() == "-blob_sf":
            coord, sf = arg[k + 1], arg[k + 2]
            find_blob_density(coord, sf)

        elif arg[k].lower() == "-blob_map":
            coord, map = arg[k + 1], arg[k + 2]
            find_blob_density_map(coord, map)

        elif arg[k].lower() == "-map":  # display map by FWT/PHWT
            sf_util.display_map_coot(arg[k + 1], arg[k + 2])  # input pdb/cif & sfcif

        elif arg[k].lower() == "-superpose":
            run_superpose(k, arg, outfile)

        elif arg[k].lower() == "-pdb":
            pdbfile, sffile = util.get_file_by_pdbid(arg[k + 1], "pdb")
            print("pdbfile= %s " % pdbfile)

        elif arg[k].lower() == "-cif":
            ciffile, sffile = util.get_file_by_pdbid(arg[k + 1], "cif")
            print("ciffile= %s " % ciffile)

        elif arg[k].lower() == "-sf":
            pdbfile, sffile = util.get_file_by_pdbid(arg[k + 1], "sf")
            print("sffile= %s " % sffile)

        elif arg[k].lower() == "-pdbid":
            pdbfile, sffile = util.get_file_by_pdbid(arg[k + 1], "pdbid")
            print("pdbfile, sffile= %s  %s " % (pdbfile, sffile))

        elif arg[k].lower() == "-cifid":
            ciffile, sffile = util.get_file_by_pdbid(arg[k + 1], "cifid")
            print("ciffile, sffile=  %s  %s " % (ciffile, sffile))

        elif arg[k].lower() == "-geom":
            calc_geometry(arg[k + 1])

    print("Time elapsed : %.2f(sec)" % (time() - start_time))


##########################################################
def calc_geometry(ifile):
    print("Calculating geometry statistics by molprobtity..")
    #    ofile='molprobity.out'
    #    os.system('phenix.clashscore %s >%s'%(ifile, ofile))
    os.system("phenix.molprobity  %s > /dev/null" % (ifile))
    os.system("rm molprobity_coot.py molprobity_probe.txt")
    print("The output geometry statistics = molprobity.out")


##########################################################
def plot_file(file_in):
    file = file_in.replace("_", "-")

    title = "%s: before (red) and after (green) refinement" % file
    xrange, yrange, xlabel, ylabel = "", "", "entry", "%s" % file
    bar, rot, key, style = 0, 1, 0, 0
    #    style 0
    plot = "plot '%s' using 2:xtic(1) , '' using 3:xtic(1) " % (file_in)

    util.gnu_plot(file, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)


##########################################################
def get_285_by_mr(pdbfile, sffile, id):
    """coord. not in the crystal frame. Using MR to get new xyz and align the
    original xyz to the new_xyz to get Matrix(X0). Then generate Remark 285.
    """
    head = """REMARK 285                                                                      
REMARK 285 THE ENTRY COORDINATES                                                
REMARK 285 ARE NOT PRESENTED IN THE STANDARD CRYSTAL FRAME.                     
REMARK 285                                                                      
REMARK 285 IN ORDER TO GENERATE THE CRYSTAL AU, APPLY THE                       
REMARK 285 FOLLOWING TRANSFORMATION MATRIX OR MATRICES AND SELECTED             
REMARK 285 BIOMT RECORDS TO THE COORDINATES, AS SHOWN BELOW.                    
"""
    pdbnew = pdbfile + "_epmr"
    phelog = "phenix_superpose.log"
    delete_file(pdbnew, phelog)

    if id == 0:  # sf file
        arg = "auto -mr epmr -pdb %s -sf %s" % (pdbfile, sffile)
        os.system(arg)
    else:
        pdbnew = sffile

    arg = "auto -match %s %s " % (pdbnew, pdbfile)
    os.system(arg)

    if not os.path.exists(phelog):
        print("Error! problem of alignment  (%s & %s). Check phenix env." % (pdbnew, pdbfile))
        return

    fp = open(phelog, "r").readlines()

    chain, chain_range = chain_res(pdbfile, 1)
    print(chain_range.keys())

    r1, r2, r3, t1 = ([],) * 4
    matrix = ""
    for i, x in enumerate(fp):
        if "r={{" in x and "}," in x:
            r1 = x[3:].replace("{", "").replace("}", "").split(",")
            r2 = fp[i + 1][3:].replace("{", "").replace("}", "").split(",")
            r3 = fp[i + 2][3:].replace("{", "").replace("}", "").split(",")

        elif "t={{" in x and "}," in x:
            t1 = x[3:].replace("{", "").replace("}", "").split(",")
            tt1 = [float(x.strip()) for x in t1]
            print(head)

            rr1 = [float(x.strip()) for x in r1 if len(x.strip()) > 0]
            rr2 = [float(x.strip()) for x in r2 if len(x.strip()) > 0]
            rr3 = [float(x.strip()) for x in r3 if len(x.strip()) > 0]

            # print(rr1,rr2, rr3, tt1 )
            m1 = "REMARK 285 X0  1 %11.6f%11.6f%11.6f%11.6f \n" % (rr1[0], rr1[1], rr1[2], tt1[0])
            m2 = "REMARK 285 X0  2 %11.6f%11.6f%11.6f%11.6f \n" % (rr2[0], rr2[1], rr2[2], tt1[1])
            m3 = "REMARK 285 X0  3 %11.6f%11.6f%11.6f%11.6f \n" % (rr3[0], rr3[1], rr3[2], tt1[2])
            s1 = "REMARK 285 CRYSTAL AU =\n"
            s2 = "REMARK 285 (X0) * CHAINS "
            s3 = ""
            for x in chain.keys():
                s3 = s3 + x + ","
            matrix = m1 + m2 + m3 + s1 + s2 + s3[: len(s3) - 1] + "\n"

    print(matrix)

    pdb285 = pdbfile + "_r285"
    fw = open(pdb285, "w")
    pdb = open(pdbfile, "r").readlines()
    for i, x in enumerate(pdb):
        if "REMARK" in x[:6] and int(x[6:10]) < 290 and "REMARK 290" in pdb[i + 1]:
            fw.write(x)
            fw.write(head)
            fw.write(matrix)
        else:
            fw.write(x)
    fw.close()

    if id == 0:
        arg = "%s/bin/dcc -refmac -pdb %s -sf %s " % (os.environ["DCCPY"], pdb285, sffile)
        os.system(arg)

    print("The new pdb=%s" % pdb285)


##########################################################
def find_blob_density_map(coord, map):
    """id=0, map=SF; id=1, map=map"""

    # xyzlim=find_xyzlimt(0, coord)
    density = coord + "_peak.pdb"
    arg = """#!/bin/tcsh -f
#  Find position of significant peaks (>=1.8 rms) in the map
###################################################################
peakmax mapin %s  xyzout %s <<eof >peak.log
NUMPEAKS 8000000
THRESHOLD RMS 5.00
EXCLUDE EDGE
eof
""" % (
        map,
        density,
    )

    scr = "peakmax.csh"
    fw = open(scr, "w")
    fw.write(arg)
    fw.close()

    os.system("chmod +x %s; ./%s" % (scr, scr))
    nblob = get_blobs(coord, density)


# print xyzlim


##########################################################
def find_xyzlimt(extend, coord):
    """extend: extend the coordinate. return the xyzlimit in fract"""

    cell, xx, yy, zz = [], [], [], []
    xyzlim = ""

    if not util.check_file(100, coord):
        return xyzlim
    fp = open(coord, "r")
    for x1 in fp:
        if "CRYST1" in x1[:6]:
            cell = [float(p) for p in x1[8:54].split()]
        elif "ATOM" in x1[:4] or "HETATM" in x1[:6]:
            xx.append(float(x1[30:38]))
            yy.append(float(x1[38:46]))
            zz.append(float(x1[46:54]))
        elif "ENDMDL" in x1[:6]:
            break

    fp.close()

    if not xx or not yy or not zz:
        #  print('Error: %s can not be found in the coordinate. try a new id. ' %(compid))
        return xyzlim

    frac, orth = util.frac_orth_matrix(cell)  # get matrix
    xx_min, xx_max = min(xx) - extend, max(xx) + extend
    yy_min, yy_max = min(yy) - extend, max(yy) + extend
    zz_min, zz_max = min(zz) - extend, max(zz) + extend

    xf_min = util.matrix_prod(frac, [xx_min, yy_min, zz_min])
    xf_max = util.matrix_prod(frac, [xx_max, yy_max, zz_max])

    xyzlim = "%.3f %.3f  %.3f %.3f  %.3f %.3f" % (xf_min[0], xf_max[0], xf_min[1], xf_max[1], xf_min[2], xf_max[2])
    return xyzlim


##########################################################
def get_blobs(coord, density):
    #  print 'Finding blob of density (%s,%s) ...' %(coord, density)
    fp1 = open(density, "r")
    fp2 = open(coord, "r").readlines()

    blob = []
    for x in fp1:
        if "ATOM" not in x[:4] and "HETA" not in x[:4]:
            blob.append(x)
            continue
        # print x.strip()
        id = 0
        for y in fp2:
            if "ATOM" not in y[:4] and "HETA" not in y[:4]:
                continue
            dx = math.fabs(float(x[28:38]) - float(y[28:38]))
            if dx > 2:
                continue
            dy = math.fabs(float(x[38:46]) - float(y[38:46]))
            if dy > 2:
                continue
            dz = math.fabs(float(x[46:54]) - float(y[46:54]))
            if dz > 2:
                continue

            if dx < 2 and dy < 2 and dz < 2:  # possible overlap
                d = math.sqrt(dx * dx + dy * dy + dz * dz)
                if d < 2.2:  # overlap with coord
                    id = 1
                    break
                # print y.strip(), dx,dy,dz, d
        if id == 0:
            blob.append(x)

    for x in blob:
        print(x.strip())

    return blob


##########################################################
def get_pdbredo_report(pdbid):
    """get the pdb_redo generated file"""

    url = "http://www.cmbi.ru.nl/pdb_redo"

    id = pdbid.lower()
    pdb = "%s_final.pdb" % id
    mtz = "%s_final.mtz" % id
    fpdb = "%s/%s/%s/%s" % (url, id[1:3], id, pdb)
    fmtz = "%s/%s/%s/%s" % (url, id[1:3], id, mtz)

    os.system("wget  %s " % fpdb)
    os.system("wget  %s " % fmtz)

    print("\nThe pdb/mtz file = (%s  %s)" % (pdb, mtz))


##########################################################
def gen_vtf_report(xyzf, sff, outf):
    """generate wwpdf report from xyz and sf files. outf is the ID. (Use API)"""

    wwpdf_exec = "/net/techusers/hyang/prog-vari/bin/wwpdb-val.sh "

    cif = xyzf
    sf = sff
    if not util.is_cif(xyzf):
        cif = util.pdb2cif(xyzf)  # a pdbfile

    if not util.is_cif(sff):
        sf = prog.sf_convertor(xyzf, sff, "mmcif", "")

    os.system("%s  %s %s %s  >& all_wwpdf.log" % (wwpdf_exec, cif, sf, outf))

    print("The output files %s.pdf  %s.xml  %s.log  \n" % (outf, outf, outf))


##########################################################
def get_vtf_report(pdbid, id):
    """extract wwpdf report from database."""

    url = "http://ftp.wwpdb.org/pub/pdb/validation_reports"

    pdbid = pdbid.lower().strip()

    idd = "_validation.pdf"
    if id == "fpdf":
        idd = "_full_validation.pdf"
    elif id == "xml":
        idd = "_validation.xml.gz"

    mid = pdbid[1:3]
    zfile = "%s/%s/%s/%s%s" % (url, mid, pdbid, pdbid, idd)

    file = "%s%s" % (pdbid, idd)
    os.system("wget -q %s" % zfile)

    if id == "xml":
        os.system("gunzip -f %s" % file)
        file = "%s%s" % (pdbid, idd[:-3])

    print("\nThe VTF=%s " % file)


##########################################################
def count_atom(file):
    natm, natmh = 0, 0
    for x in open(file, "r").readlines():
        if "ATOM" in x[:4] or "HETA" in x[:4]:
            if " H" in x[76:78] or " D" in x[76:78]:
                natmh = natmh + 1
            else:
                natm = natm + 1
    return natm, natmh


##########################################################
def check_new_occ(fold, fnew, docc, y, nres2):
    tmp = ""
    for x in nres2:
        doccn = check_occ(x)
        if y in doccn.keys():
            t1 = ["%.2f" % t for t in docc[y]]
            t2 = ["%.2f" % t for t in doccn[y]]
            if len(t1) < 2:
                continue
            s1, s2 = sum(docc[y]), sum(doccn[y])
            tmp = "%s :: %s %s=%.2f %s  %s=%.2f\n" % (y, fold, str(t1), s1, fnew, str(t2), s2)
            print(tmp)
            break

    return tmp


##########################################################
def separate_pdb(pdbfile):
    """separate PDB file by polymer pdb1 and ligands pdb2, water pdb3"""

    pdb1 = pdbfile + "tmp1"
    pdb2 = pdbfile + "tmp2"
    pdb3 = pdbfile + "tmp3"

    fw1, fw2, fw3 = open(pdb1, "w"), open(pdb2, "w"), open(pdb3, "w")

    n = 0
    fr, fr1 = [], open(pdbfile, "r").readlines()

    for x in fr1:
        if ("ATOM" in x[:4] or "HETA" in x[:4]) and ("HOH" in x[17:20] or "DOD" in x[17:20]):
            fw3.write(x)
        else:
            fr.append(x)

    length = len(fr)
    i = 0
    for i in range(length):
        k = length - i - 1
        if ("CONECT" in fr[k][:6] or "MASTER" in fr[k][:6] or "END" in fr[k][:3] or "HETATM" in fr[k][:6]) and ("ATOM" in fr[k - 1][:4] or "TER " in fr[k - 1][:4]):
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
            ch = ln[21:22]
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
    return chain, chain_range


##########################################################
def compare_diff(fold, fnew, chain, atomo, reso, nres2):
    """
    reso: all the atoms in the residue([]); atomo: uniqe atoms;
    nres2: all for new PDB ([[res1], [res2] ...])
    """
    aa, ma, lig, hoh, hh = 0, 0, 0, 0, 0
    info = []
    for x in reso:
        tmp = x[11:16] + x[17:26]
        # print(tmp, atomo)
        m1, m2, m3 = x[17:20], x[29:54], x[54:60]
        if tmp == atomo:  # match uniq atom
            # print ('old:' , x)

            for y in nres2:  # new residue
                for z in y[:]:
                    n1, n2, n3 = z[17:20], z[29:54], z[54:60]
                    if m1 == n1 and m2 == n2 and m3 != n3:
                        if (m3.strip() == "1.00" and n3.strip() == "0.50") or n3.strip() == "1.00":
                            continue

                        t1 = "Error! files (%s : %s) differ in OCC\n" % (fold, fnew)
                        info.append(t1)
                        info.append("old: " + x)
                        info.append("new: " + z)
                        if "ATOM" in z[:4]:
                            aa = aa + 1
                            if " H " in z[76:79]:
                                hh = hh + 1

                        elif "HETA" in z[:4]:
                            if "HOH" in z[17:20]:
                                hoh = hoh + 1
                            else:  # either ligand or modified residue
                                ch, nres = z[21:22], int(z[22:26])
                                if ch in chain.keys():
                                    if nres >= chain[ch][0] and nres <= chain[ch][1]:
                                        ma = ma + 1
                                    else:
                                        lig = lig + 1

                        # print(t1,x,z)
                        y.remove(z)

    return info, aa, ma, lig, hoh, hh


##########################################################
def check_occ(pdb):
    """return a dic with key and a list of occ"""
    id, d1 = 0, {}
    for x in pdb:
        tmp = x[11:16] + x[17:27]
        if tmp not in d1.keys():
            d1[tmp] = []
        d1[tmp].append(float(x[54:60]))
    return d1


##########################################################
def chain_res_list(pdb):
    """put chain+residue in a list [[],[]...]"""

    n = len(pdb)
    if n == 0:
        print("Warning: empty PDB list")
        return []

    m, tmp, pdbn = 0, [], []
    for x in pdb:
        if "ATOM" in x[:4] or "HETA" in x[:4]:
            m = m + 1
            id0 = id = x[21:26]
            if m < n:
                id = pdb[m][21:26]
            tmp.append(x)

            if id0 != id or m == n:
                pdbn.append(tmp)
                tmp = []

    # print(pdbn)
    return pdbn


##########################################################
def sort_column(pdb, k1, k2):
    """Sort the columns (from k1 to k2) of PDB file in order."""

    id = 1
    n, d, pdbn = 0, [], []
    for ln in pdb:
        n = n + 1
        if "ATOM" in ln[:4] or "HETATM" in ln[:6] or "ANISOU" in ln[:6]:
            ch = ln[21:22]
            d.append(ln)
            if id == 1:
                cond = n == len(pdb)
            if id == 2:
                cond = n == len(pdb) or ch != pdb[n][21:22]
            if cond:
                d1 = [[x[:k1], x[k1:k2], x[k2:]] for x in d]  # separate 3 column
                d1.sort(key=lambda y: y[1])  # sort 2th column
                d = []
                for x in d1:
                    pdbn.append("".join(x))

        else:
            if len(d) > 0:
                d1 = [[x[:k1], x[k1:k2], x[k2:]] for x in d]  # separate 3 column
                d1.sort(key=lambda y: y[1])  # sort 2th column
                d = []
                for x in d1:
                    pdbn.append("".join(x))
            pdbn.append(ln)

    return pdbn


##########################################################
def pickup_item_from_firstfile(file1, file2):
    """grep items from first column of file1 and extract the
    lines in files containing the item.
    """

    script = """#!/bin/csh

set list = `cat $1 | awk '{print $1} '  `
set log = "${1}_picked"
rm -f $log 
foreach id ($list)
echo "========greping $id========="
grep -i $id $2 >>$log
end

"""
    scr = "grep_item.csh"
    fw = open(scr, "w")
    fw.write(script)
    fw.close()
    os.system("chmod +x %s; %s %s %s" % (scr, scr, file1, file2))


##########################################################
def residue_num_seq(*file):
    pdb = file[0]
    out = pdb + "_new"
    fw = open(out, "w")
    fr = open(pdb, "r")

    if len(file) > 1:
        pdb, ss = file[0], file[1]
        t = ss.split("_")
        s1, s2 = t[0], t[1]
        ch11, nres11, ch22, nres22 = s1[:1], s1[1:], s2[:1], s2[1:]
        print("input residue= %s %s %s %s" % (ch11, nres11, ch22, nres22))

        nr = int(nres22)

        res, ch, nres = "", "", ""
        for x in fr:
            if "ATOM" in x[:4] or "HETA" in x[:4] or "ANIS" in x[:4]:
                res, ch, nres = x[17:20], x[21:22], x[22:26]
                if ch.strip() == ch11 and nres.strip() == nres11:
                    s1, s2, s3, s4, s5 = x[:21], "%s" % ch22, "%4d" % nr, " ", x[27:]
                    s = s1 + s2 + s3 + s4 + s5
                    fw.write(s)
                    break
                else:
                    fw.write(x)
            else:
                fw.write(x)

        print("start residue=%s %s %s" % (res, ch, nres))

        for x in fr:
            res2, ch2, nres2 = x[17:20], x[21:22], x[22:26]
            if ("ATOM" in x[:4] or "HETA" in x[:4] or "ANIS" in x[:4]) and ch2 == ch11:
                if res2 != res:
                    nr = nr + 1
                    res = res2

                s1, s2, s3, s4, s5 = x[:21], "%s" % ch22, "%4d" % nr, " ", x[27:]
                s = s1 + s2 + s3 + s4 + s5
                fw.write(s)
            else:
                fw.write(x)
        fw.close()
        print("The new pdbfile=%s" % out)
        return

    ch_old, nres_old, nr = "?", -999, 0
    for x in fr:
        if "ATOM" in x[:4] or "HETA" in x[:4] or "ANIS" in x[:4]:
            # print(x[21:22], x[22:26])
            ch = x[21:22]
            s = ""
            if ch == ch_old:
                res = x[17:20]
                if res != res_old:
                    nr = nr + 1
                    res_old = res

                s1, s2, s3, s4 = x[:22], "%4d" % nr, " ", x[27:]
                s = s1 + s2 + s3 + s4
                # s=''.join(s1,s2,s3,s4)
            else:
                ch_old = ch
                res_old = x[17:20]
                y = x[22:26]
                if y[3:4].isalpha():
                    nres = 1
                else:
                    nres = int(y)
                nr = nres
                s1, s2, s3, s4 = x[:22], "%4d" % nr, " ", x[27:]
                s = s1 + s2 + s3 + s4

            fw.write(s)
        else:
            fw.write(x)

        fw.close()
        print("The new pdbfile=%s" % out)
        return


##########################################################
def delete_file(*files):
    for x in files:
        os.system("rm -f " + x)


##########################################################
def calc_scale_matrix(pdbfile_in):
    """compare  scale card in pdb with the calculated one
    refer to : http://en.wikipedia.org/wiki/Fractional_coordinates
    If the fractional coordinate system has the same origin as the cartesian
    coordinate system, the a-axis is collinear with the x-axis, and the b-axis
    lies in the xy-plane, fractional coordinates can be converted to cartesian
    coordinates through the following transformation matrix:

    """
    # acc=0.002  # used for detecting non-cyrstal frame
    acc = 0.000001
    be_cif = util.is_cif(pdbfile_in)
    pdbfile = pdbfile_in
    if be_cif:
        pdbfile = cif.cif2pdb(pdbfile_in)

    fp = open(pdbfile, "r")

    cell, s1, s2, s3 = [], [], [], []
    for x in fp:
        if "CRYST1" in x[:6]:
            v = x[6:].split()
            cell = [float(v[i]) for i in range(6)]
        elif "SCALE1" in x[:6]:
            v = x[6:].split()
            s1 = [float(v[i]) for i in range(3)]
        elif "SCALE2" in x[:6]:
            v = x[7:].split()
            s2 = [float(v[i]) for i in range(3)]
        elif "SCALE3" in x[:6]:
            v = x[6:].split()
            s3 = [float(v[i]) for i in range(3)]
        elif "ATOM" in x[:4] or "HETA" in x[:4]:
            break
    fp.close()

    if len(cell) != 6:
        print("No Cell parameters in " + pdbfile)
        return

    a, b, c = cell[0], cell[1], cell[2]
    sa = 3.141592654 / 180.0
    alpha, beta, gamma = sa * cell[3], sa * cell[4], sa * cell[5]

    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    cg, sg = math.cos(gamma), math.sin(gamma)

    vol = math.sqrt(1 - ca**2 - cb**2 - cg**2 + 2 * ca * cb * cg)
    #    vol = sg*sb*sa
    sc1 = [sc11, sc12, sc13] = [1 / a, -cg / (a * sg), (ca * cg - cb) / (a * vol * sg)]
    sc2 = [sc21, sc22, sc23] = [0, 1.0 / (b * sg), (cb * cg - ca) / (b * vol * sg)]
    sc3 = [sc31, sc32, sc33] = [0, 0, sg / (c * vol)]

    sct = []
    sct.extend(sc1)
    sct.extend(sc2)
    sct.extend(sc3)

    sst = []
    sst.extend(s1)
    sst.extend(s2)
    sst.extend(s3)

    df = 0
    if len(sct) == 9 and len(sst) == 9:
        for i in range(9):
            if math.fabs(sct[i] - sst[i]) > acc:
                df = 1
                break
    if not sst:
        print("Warning: Scale matrix is reported. The calculated is below:")
        print("%10.6f%10.6f%10.6f " % (sc11, sc12, sc13))
        print("%10.6f%10.6f%10.6f " % (sc21, sc22, sc23))
        print("%10.6f%10.6f%10.6f " % (sc31, sc32, sc33))

    elif df == 1:
        print("Error! Calculated and reported SCALE matrix is different(acc=%.6f)." % (acc))
        print("         calculated                  reported")
        print("%10.6f%10.6f%10.6f  %10.6f%10.6f%10.6f" % (sc11, sc12, sc13, s1[0], s1[1], s1[2]))
        print("%10.6f%10.6f%10.6f  %10.6f%10.6f%10.6f" % (sc21, sc22, sc23, s2[0], s2[1], s2[2]))
        print("%10.6f%10.6f%10.6f  %10.6f%10.6f%10.6f" % (sc31, sc32, sc33, s3[0], s3[1], s3[2]))
    else:
        print("Calculated and reported SCALE matrix is the same(acc=%.6f)." % (acc))

    if be_cif:
        util.delete_file(pdbfile)


##########################################################
def change_bfactor(pdbfile, t):
    """Add a constant to Biso or ANISOU which comes from overall scaling"""

    fp = open(pdbfile, "r").readlines()
    out = pdbfile + "_newB"
    fw = open(out, "w")
    n = len(t)
    v = t[0]
    if n == 6:
        v = (t[0] + t[1] + t[2]) / 3.0
    for x in fp:
        if "ATOM" in x[:4] or "HETA" in x[:4]:
            t1, t2, t3 = x[:60], float(x[60:66]) + v, x[66:].rstrip()
            tmp = "%s%6.2f%s\n" % (t1, t2, t3)
            fw.write(tmp)
        elif "ANISOU" in x[:6] and n == 6:
            u = [float(x1) for x1 in x[29:71].split()]
            v1 = [(x1 * 10**4) / (8 * 3.14156**2) for x1 in t]
            ut = [int(u[i] + v1[i]) for i in range(6)]
            # print(x)
            y1, y3 = x[:28], x[70:]
            y2 = "".join("%7d" % ut[i] for i in range(6))
            tmp = y1 + y2 + y3
            fw.write(tmp)
        else:
            fw.write(x)
    fw.close()
    print("The new pdb = " + out)
    return out


##########################################################
def mean_std_bin(list1, file, m, n, k):
    """get min,max,number,mean,standard deviation of nth column  for each
    resolution bins. (resolution & data in m & n columns;  k,number of bins)

    """

    all_data = []

    fp = list1
    if len(list1) == 0:
        fp = open(file, "r").readlines()
    res = []
    for x in fp:  # put resol in a list
        t = float(x.strip().split()[m - 1])
        res.append(t)
    min1, max1 = min(res), max(res)
    shell = reso_shell(min1, max1, k)

    for y in shell:
        s1, data = 0, []
        for x in fp:
            t = x.strip().split()
            y1, x1 = float(t[m - 1]), float(t[n - 1])
            if y1 >= y[0] and y1 < y[1]:
                data.append(x1)
                s1 = s1 + x1
        min1, max1, mean, std, num = mean_std(data, "", 1)
        avg = (y[0] + y[1]) / 2.0
        tmp = "%6.2f %6.2f %6.2f  %8.3f %8.3f %8.3f %8.3f %5d" % (y[0], y[1], avg, min1, max1, mean, std, num)
        all_data.append(tmp)

    print("\nres1, res2, avg, min1, max1, mean, std, num")
    for x in all_data:
        print(x)
    return all_data


##########################################################
def reso_shell(min1, max1, nstep):
    """get resolution shell from the range min-max and step"""

    step = (max1 - min1) / nstep
    shell = []
    y1 = min1

    for x in range(nstep):
        y2 = y1 + step
        shell.append([y1, y2])
        y1 = y2

    return shell


##########################################################
def mean_std(list1, file, n):
    """get min,max,mean,standard deviation of nth column in the file
    id=0, use file; id=1, use list1
    """

    id = len(list1)
    res = []
    if id == 0:
        fp = open(file, "r").readlines()
        for x in fp:
            v = float(x.strip().split()[n - 1])
            res.append(v)
    else:
        res = list1

    nlist = len(res)
    min1, max1, mean = min(res), max(res), sum(res) / float(nlist)

    s = 0
    for x in res:
        a = x - mean
        a2 = a * a
        s = s + a2
    std = math.sqrt(s / float(nlist))

    print("min, max, mean, std, num = %8.3f %8.3f %8.3f %8.3f %6d" % (min1, max1, mean, std, nlist))

    return min1, max1, mean, std, nlist


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
        print("Error: No ANISOU records in " + pdbfile)
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
def check_pdb_occ(pdbfile):
    tf1 = pdbfile + "_tmp"

    os.system('egrep "^ATOM|^HETATM"  %s  >%s' % (pdbfile, tf1))  # keep ATOM/HETA

    ddf1 = open(tf1, "r").readlines()
    nres1 = chain_res_list(ddf1)  # put chain+residue in a list [[],[]...]

    for x in nres1:  # for each residue in pdb
        docc = check_occ(x)  # return uniq atom(as key) and list of occ as val
        # print(docc)
        for y in docc.keys():
            s = sum(docc[y])

            if s > 1.2:  #
                print("Error: Sum of occupancy with atom(%s) >1.2" % y)

    delete_file(tf1)


##########################################################
def correct_occ_on_symm(file):
    """if file is cif, export corrected cif file; if pdb, then export pdb"""

    if not util.is_cif(file):  # pdbfile
        find_atom_on_symm(file, "pdb")
    else:  # cif
        pdb = prog.run_maxit(file, 2, "")
        site = find_atom_on_symm(pdb, "cif")
        newcif = parse_cif(file, site)
        delete_file(pdb)
        os.system("diff %s %s" % (file, newcif))


##########################################################
def find_atom_on_symm(pdbfile, mode):
    """find atoms on special position and return site = {id: fold, }
    mode = pdb; correct pdbfile. Otherwise,
    """

    check_pdb_occ(pdbfile)

    arg = "phenix.pdbtools " + pdbfile
    lines = os.popen(arg).read().split("\n")
    delete_file(pdbfile + "_modified.pdb", "maxit.err", "CifParser.log")

    if len(lines) < 50:
        print("\nError: phenix environment not set.")
        sys.exit()

    site = {}
    for i, x in enumerate(lines):
        if "Number of sites at special positions:" in x:
            print(x.strip())
        elif "pdb=" in x and "(" in x and ")" in x and "original" in x:
            n = 0
            tmp = x.split("pdb=")[1].split('"')[1]
            if "site sym" in lines[i + 1] and "(" in lines[i + 1] and "exact" in lines[i + 1]:
                n = int(lines[i + 1].split("sym")[1].split("(")[0])

            site[tmp] = n

    rmk375 = os.popen('egrep "^REMARK 375.*ON A SPECIAL POSITION" %s' % pdbfile).read().split("\n")
    for x in sorted(site.keys()):
        nn = 0
        for y in rmk375:
            if len(y[16:25]) > 0 and y[16:25].strip() in x:
                nn = 1
                break

        if nn == 0:
            print("Warning: SKIPPED! (%s : symmetry fold=%d), not match PDB criteria(min. dist=0.15)" % (x, site[x]))
            site.pop(x)
        else:
            print("%s :   symmetry fold=%d  " % (x, site[x]))

    if len(site) == 0 or mode == "cif":
        return site
    print("\nDoing occupancy correction if necessary\n")

    out = pdbfile + "_occ"
    fp, fw = open(pdbfile, "r"), open(out, "w")

    for x in fp:
        if "ATOM" in x[:4] or "HETA" in x[:4]:
            line = check_occ_site(x, site)
            if len(line) > 10:
                print("old:" + x.strip())
                print("new:" + line.strip())
                fw.write(line)
            else:
                fw.write(x)
        else:
            fw.write(x)

    fp.close(), fw.close()

    return site


##########################################################
def make_cif_dic(file):
    """make a cif dic from the file
        description after '==' is for category and '::' is for data items .  below is an example
    _pdbx_phasing_MR.d_res_high_fit ::The highest resolution limit used for rigid body refinement after molecular replacement (MR) solution.  :: float
    _pdbx_phasing_MR.d_res_low_fit :: The lowest resolution limit used for rigid body refinement after molecular replacement (MR) solution.:: float

    """

    fp = open(file, "r").readlines()
    outf = file + ".dic"
    fw = open(outf, "w")

    sp1 = "   "
    sp2 = "       "
    for x in fp:
        if "==" in x:  # category name
            t = x.split("==")
            cate = t[0].strip()
            cat = cate[1:]
            ss = t[1].strip()
            fw.write("\n######### Category=%s #######\n" % cat)
            fw.write("save%s\n" % cate)
            fw.write("%s_category.description\n" % sp1)
            fw.write(";%s%s\n;\n" % (sp2, ss))

            fw.write("%s_category.id              %s \n" % (sp1, cat))
            fw.write("%s_category.mandatory_code  no \n" % (sp1))
            fw.write("%ssave_\n\n" % (sp1))

        elif "::" in x:  # item names
            t = x.split("::")
            if len(t) < 3:
                print("problem with=%s" % x)
                continue
            item = t[0].strip()
            cat = item.split(".")[0][1:]
            ss = t[1].strip()
            code = t[2].strip()
            fw.write("save_%s\n" % item)
            fw.write("%s_item_description.description\n" % sp1)
            fw.write(";%s%s\n;\n" % (sp2, ss))
            fw.write("%s_item.name           '%s'\n" % (sp1, item))
            fw.write("%s_item.category_id     %s\n" % (sp1, cat))
            fw.write("%s_item.mandatory_code  no\n" % sp1)
            fw.write("%s_item_type.code       %s\n" % (sp1, code))
            fw.write("%ssave_\n\n" % sp1)

    print("\nThe output file = %s" % outf)


##########################################################


def check_occ_site(x, site):
    occ = {2: 0.5, 3: 0.33, 4: 0.25, 6: 0.16}
    line = ""
    for y in site.keys():
        if y in x:
            oc = float(x[54:60])
            if site[y] not in occ.keys():
                print("Warning: symmetry fold (%d) not (2,3,4,6). site=%s " % (site[y], site))
                continue

            if oc != occ[site[y]]:  #
                line = x[:54] + "%6.2f" % occ[site[y]] + x[60:]
                site.pop(y)
                break

    return line


##########################################################


def parse_cif(pdb, site):
    occ_assign = {2: 0.5, 3: 0.33, 4: 0.25, 6: 0.16}
    table = "_atom_site"
    fp = open(pdb, "r").readlines()
    fp1, fp2, fp3, fp4 = [], [], [], []
    n1, n2, n3, n4 = 0, 0, 0, len(fp)
    out = pdb + "new"
    fw = open(out, "w")

    for i, x in enumerate(fp):
        cate, item = "", ""
        if "." in x:
            t = x.split(".")
            cate, item = t[0].strip(), t[1].strip()

        if cate == table and len(item) > 1 and i > 0 and table not in fp[i - 1]:
            n1 = i
        elif cate == table and len(item) > 1 and i < n4 - 1 and table not in fp[i + 1]:
            n2 = i + 1
        elif n1 > 0 and "#" in x:
            n3 = i
            break

    fp1 = fp[:n1]
    fp2 = fp[n1:n2]  # table
    fp3 = fp[n2:n3]  # coord
    fp4 = fp[n3:n4]

    tab = {}
    nfield = n2 - n1
    for i, x in enumerate(fp2):  # parse table
        t = x.split(".")
        cate, item = t[0].strip(), t[1].strip()
        tab[item] = i

    atom = "auth_atom_id"
    resid = "auth_comp_id"
    chain = "auth_asym_id"
    nres = "auth_seq_id"
    occ = "occupancy"

    # site={' O   HOH B 104 ': 2, 'MG    MG A 100 ': 2}

    fp3_new = []
    for i, x in enumerate(fp3):  # coord
        t = x.split()
        if len(t) != nfield:
            print("Warning: nfield != parsed column")

        occ_old = float(t[tab[occ]])
        ln = ""
        for y in site.keys():  # check through site
            atom1, resid1 = y[:4].strip(), y[4:8].strip()
            chain1, nres1 = y[9:10].strip(), y[10:14].strip()
            if site[y] in occ_assign.keys():
                occ_new = float(occ_assign[site[y]])

            if t[tab[atom]] == atom1 and t[tab[resid]] == resid1 and t[tab[chain]] == chain1 and t[tab[nres]] == nres1:
                if occ_new != occ_old:
                    t[tab[occ]] = str(occ_new)
                    # print( atom1, resid1, chain1, nres1, site[y], occ_new, occ_old)
                    ln = " ".join(t) + "\n"
                    break

        if not ln:
            fp3_new.append(x)
        else:
            fp3_new.append(ln)

    for x in fp1:
        fw.write(x)
    for x in fp2:
        fw.write(x)
    for x in fp3_new:
        fw.write(x)
    for x in fp4:
        fw.write(x)
    fw.close()
    return out


##########################################################
def remove_h_atom(pdbfile):
    out = pdbfile + "_noH"
    fp = open(pdbfile, "r")
    fw = open(out, "w")

    for x in fp:
        if ("ATOM " in x[:6] or "HETATM" in x[:6]) and " H" in x[76:78]:
            continue
        fw.write(x)
    fp.close(), fw.close()
    print("The output file without H atoms = %s" % out)


##########################################################
def map_header(mapfile):
    log = mapfile + "_header"
    scr = mapfile + ".sh"

    arg = "mapdump mapin %s <<eof >%s \neof\n" % (mapfile, log)

    fw = open(scr, "w")
    fw.write(arg)
    fw.close()
    os.system("chmod +x %s ; ./%s" % (scr, scr))
    print("The header file of the map = %s" % log)


##########################################################
def run_superpose(k, arg, outfile):
    """ """

    if not outfile:
        outfile = "superpose.out"
    narg = len(arg)
    ss = ""
    for x in range(k + 1, narg):
        if arg[x] == "-o":
            break
        ss = ss + arg[x] + " "

    print(ss)
    args = "superpose %s  >%s " % (ss, outfile)
    os.system(args)

    if not os.path.exists(outfile) or os.path.getsize(outfile) < 100:
        print("Alignment failed. Check input.")
        return

    fp = open(outfile, "r").readlines()
    fw = open(outfile, "w")

    space = ""
    for i, x in enumerate(fp):
        if "Residue alignment:" in x:
            space = "    "
        if not space:
            fw.write(" %s  %s" % (space, x))
        #            print space , x.strip()
        else:
            brk = " "
            if ":" in x:
                ln = x.split()
                aa = []
                for y in ln:
                    if ":" in y:
                        t = y.split(":")
                        if not t[1] in aa:
                            aa.append(t[1])

                if len(ln) > 2 and i < len(fp) - 1:
                    lln = fp[i + 1].split()  # the next line
                    if util.is_number(ln[1]) and util.is_number(lln[1]) and int(ln[1]) + 1 != int(lln[1]):
                        brk = "?"

                if len(aa) > 1:
                    space = "--> "
            fw.write("%s %s %s" % (brk, space, x))
            #            print space , x.strip()
            space = "    "

    fw.close()
    print("\nThe output file = %s\n" % outfile)


##########################################################

if __name__ == "__main__":
    process(sys.argv)


##########################################################
# some corrections
# Only report the wrong OCC with atoms on special position (2015-04-23)
