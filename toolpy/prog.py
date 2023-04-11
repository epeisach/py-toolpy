# ===========================================================
# this module is a interface for all the external packages
# ===========================================================

import os
import util


##########################################################
def sf_convertor(pdbfile, sffile, type1, add):
    """convert to any type of formats,
    type1: the output format;  add: extra options
    """

    print("\nConverting SF file to %s format ... " % type1)
    out = sffile + "." + type1
    if not util.check_file(500, pdbfile):
        arg = "sf_convert -o %s -sf %s -out %s %s >/dev/null " % (type1, sffile, out, add)
    else:
        arg = "sf_convert -o %s -sf %s -pdb %s -out %s %s >/dev/null " % (type1, sffile, pdbfile, out, add)

    os.system(arg)
    return out


##########################################################
def run_phenix(pdbfile, sffile, type1):  # pylint: disable=unused-argument
    """run sub-programs of phenix"""
    outfile = "phenix__%s.log" % type1
    util.delete_file(outfile)
    print("\nDoing PHENIX calculation for (%s) ..." % type1)

    if type1 == "xtriage":
        arg = "phenix.xtriage %s log=%s  >/dev/null" % (sffile, outfile)
        os.system(arg)

    return outfile


##########################################################
def run_ccp4(pdbfile, sffile, type1):  # pylint: disable=unused-argument
    """run sub-programs of CCP4"""
    outfile = "ccp4__%s.log" % type1
    util.delete_file(outfile)
    print("\nDoing %s ..." % type1)

    if type1 == "ctruncate":
        arg = 'ctruncate -hklin %s -amplitudes -colin "/*/*/[FP,SIGFP]" -hklout ctruncate-SF.mtz>& %s' % (sffile, outfile)
        #        arg='ctruncate -hklin %s -amplitudes -colano "/*/*/[I(+),SIGI(+),I(-),SIGI(-)]" -hklout ctruncate-SF.mtz>& %s' %(sffile,outfile)
        os.system(arg)

    elif type1 == "pointless":
        arg = "pointless hklin %s > %s" % (sffile, outfile)
        os.system(arg)
        os.system('grep "Best Solution" %s ' % outfile)

    return outfile


##########################################################
def run_ccp4_contact(pdbfile, id):  # pylint: disable=redefined-builtin,unused-argument
    """Refer contact:  http://www.ccp4.ac.uk/dist/html/contact.html
    <mode> = ALL,IRES,ISUB,IMOL or AUTO (default: MODE IRES).
    """
    limit = 3.0
    pdb_new = pdbfile + "_new"
    fp = open(pdbfile, "r")
    fw = open(pdb_new, "w")
    util.delete_file("ccp4_contact.csh")
    fwc = open("ccp4_contact.csh", "w")
    for x in fp:
        #        if id==0 and 'SCALE' in x[:6] : continue
        if "ENDMDL" in x[:6]:
            break
        fw.write(x)
    fw.close()
    fp.close()

    log = pdbfile + "_contact"
    script = """#!/bin/csh
    #<mode> = ALL,IRES,ISUB,IMOL or AUTO (default: MODE IRES).
    #ALL: for all interatomic distances for chosen residues.
    #ISUB: intersubunit contacts (must have different chain name)
    #IMOL: intermolecular contacts
    #AUTO: as IMOL, but additional (primitive) lattice translations are generated
    contact xyzin  %s <<eof >%s
    mode   AUTO
    ATYPE ALL
    limits 0.0  %f
    """ % (
        pdb_new,
        log,
        limit,
    )

    os.system(script)
    fwc.write(script)
    fwc.close()

    arg = 'egrep "\[.*\]"  %s |wc -l' % log  # noqa: W605  # pylint: disable=anomalous-backslash-in-string
    ncont = int(os.popen(arg).read())
    print("Crystal contacts for %s (<%.1fA) = %d" % (pdbfile, limit, ncont))
    util.delete_file(pdb_new)
    # if ncont < 30: delete_file(log)
    return ncont


##########################################################


def run_dcc(pdbfile, sffile, add):
    """Add: add extra options"""

    out = pdbfile + "_dcc.cif"
    print("Calculating R/Rfree for (%s & %s)" % (pdbfile, sffile))
    arg = '%s/bin/dcc.sh -o %s -pdb %s -sf %s %s |egrep -i  "Error|Warn" |awk \'{print "%s: ", $0}\' ' % (os.environ["DCCPY"], out, pdbfile, sffile, add, sffile)
    # print 'arg=', arg
    os.system(arg)
    return out


##########################################################
def run_maxit(file, option, other):
    """some options to run maxit
    use '-exchange_in' to get tls (from mmcif of pdb_extract)
    """
    print("Converting %s by Maxit with option %d ..." % (file, option))
    nfile = file + ".pdb"

    if option == 2:
        nfile = file + ".pdb"
    else:
        nfile = file + ".cif"

    arg = "maxit-v8.01-O  -i %s -o %d  %s " % (file, option, other)
    os.system(arg)

    return nfile


##########################################################
def gnu_plot1(file, xlabel, ylabel, title, plot):
    """Plot the data using gnuplot, controled by key.
    if key==0, plot all of them
    """
    gnu_scr = file + ".gnu"
    fw = open(gnu_scr, "w")

    plot_gnu = """
#set terminal jpeg large size 840,640 #transparent nocrop enhanced font arial 8 size 420,320
set terminal png
set output '%s.png'
#set boxwidth 0.5 relative
#set style histogram clustered gap 1 title offset 0, 0, 0
#set style data histograms
set style data linespoints
set datafile missing '.'
set style fill  solid 1.0 border -1
#set key on   #default right top
set key  left top
set grid ytics
set size 1.0,1.0
set autoscale
set xtics rotate by 90
set ytics
set xlabel "%s"
set ylabel "%s"
set title "%s"

#set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 bdefault
#plot 'struct_growth.data' using 2:xtic(1) t "all",'' u 3 t "altered"
#plot 'filename'   using 1:2:3:xtic(1) t "hell" with  yerrorlines, '' u 1:4:5 t "altered"   with  yerrorlines
#plot 'datafile'   using 1:2:xtic(1) t "hell",  '' u 1:4:5 t "altered"   with  yerrorlines

plot '%s'  %s

set term x11
#replot
#pause 10

""" % (
        file,
        xlabel,
        ylabel,
        title,
        file,
        plot,
    )

    fw.write(plot_gnu)
    fw.close()
    os.system("gnuplot %s" % gnu_scr)
    print("output image = " + file + ".png")
    os.system("display " + file + ".png &")


##########################################################
def gnu_plotps(file, xlabel, ylabel, title, plot):
    """Plot the data using gnuplot, controled by key.
    if key==0, plot all of them
    Only plot ps file (if the png failed
    """
    gnu_scr = file + ".gnu"
    fw = open(gnu_scr, "w")

    plot_gnu = """
set terminal postscript eps enhanced color
set output '%s.ps'
#set boxwidth 0.5 relative
#set style histogram clustered gap 1 title offset 0, 0, 0
#set style data histograms
set style data linespoints
set datafile missing '.'
set style fill  solid 1.0 border -1
#set key on   #default right top
set key  left top
set grid ytics
set size 1.0,1.0
set autoscale
#set xtics rotate by 90
set ytics
set xlabel "%s"
set ylabel "%s"
set title "%s"

#set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 bdefault
#plot 'struct_growth.data' using 2:xtic(1) t "all",'' u 3 t "altered"
#plot 'filename'   using 1:2:3:xtic(1) t "hell" with  yerrorlines, '' u 1:4:5 t "altered"   with  yerrorlines
#plot 'datafile'   using 1:2:xtic(1) t "hell",  '' u 1:4:5 t "altered"   with  yerrorlines

plot '%s'  %s

set term x11
#replot
#pause 10

""" % (
        file,
        xlabel,
        ylabel,
        title,
        file,
        plot,
    )

    fw.write(plot_gnu)
    fw.close()
    os.system("gnuplot %s" % gnu_scr)
    print("output image = " + file + ".png")
    os.system("display " + file + ".png &")


##########################################################
def do_mr(nn, narg, arg):
    """#Do MR by EPMR/Morep/Phaser.
    http://www.csb.yale.edu/userguides/datamanip/phaser/phaser-1.3.html
    http://www.epmr.info/
    http://www.ccp4.ac.uk/html/molrep.html

    """

    pdb, sfin, nmol, prog = "", "", 1, "phaser"
    ss = ""
    for i in range(nn, narg):
        ss = ss + " %s " % arg[i]

    ss_split = ss.split(",")
    for x in ss_split:
        s = x.strip().split("=")
        key = s[0].strip().lower()
        if key == "pdb":
            pdb = s[1]
        elif key == "sf":
            sfin = s[1]
        elif key == "nmol":
            nmol = s[1]
        elif key == "prog":
            prog = s[1]

    sf = sf_convertor(pdb, sfin, "mtz", "")
    identity = "100"
    #    print ss, ss_split, prog, pdb, sf, nmol, identity

    phaser_arg = """#!/bin/csh  -f

  phaser << eof
  TITLe beta blip automatic
  MODE MR_AUTO
  HKLIn %s
  LABIn F=FP SIGF=SIGFP
  ENSEmble beta PDB %s IDENtity %s
#  ENSEmble blip PDB blip.pdb IDENtity 100
#  COMPosition PROTein SEQuence beta.seq NUM 1 #beta
#  COMPosition PROTein SEQuence blip.seq NUM 1 #blip
  SEARch ENSEmble beta NUM %s
#  SEARch ENSEmble blip NUM 1
  ROOT phaser_mr # not the default
  eof
    """ % (
        sf,
        pdb,
        identity,
        nmol,
    )

    if prog.lower() == "epmr":
        print("Doing MR by epmr ...")
        arg = "epmr %s %s -w 1 -o eprm_mr -n %s  > epmr.log" % (sf, pdb, nmol)
        os.system(arg)
        print("EPMR OUTPUT = eprm_mr.best.pdb ")

    elif prog.lower() == "phaser":
        print("Doing MR by phaser ...")

        exect = "phaser.csh"
        fw = open(exect, "w")
        fw.write(phaser_arg)
        fw.close()
        os.system("chmod +x %s ; %s >phaser.log" % (exect, exect))
        print("PHASER OUTPUT = phaser_mr.pdb")

    elif prog.lower() == "molrep":
        print("Doing MR by molrep ...")
        arg = "molrep -f %s -m %s  >molrep.log " % (sf, pdb)
        os.system(arg)
        print("MOLREP OUTPUT   ")
    else:
        print("Only the programs (EPMR, MOLREP, PHASER) are supported!!")
