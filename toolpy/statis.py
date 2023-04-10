#!/usr/bin/python

import os
import math
import util


######################################################################
def histogram_plot(infile, xcol, xlabel, ylabel, title, nstep, xxrange):
    """The infile has value at xcol column.  A histogram can be plotted with nstep.
    xcol, the column number in the file (1,2,3..)
    """

    scale = 1.5

    flist = open(infile, "r").readlines()
    data_o = []  # original data
    for x in flist:
        t = x.split()
        if xcol > len(t) or not util.is_number(t[xcol - 1]):
            continue
        data_o.append([float(t[xcol - 1])])

    col = 0
    bx1, bx2 = limits_by_iqr_2set(data_o, col, scale)
    #    bx1, bx2 ,q2= limits_by_iqr(data_o, col, scale)
    #    bx1, bx2, q2 = limits_by_iqr_mc(data_o, col, scale=1.5) # remove it for full data
    #    bx1, bx2 = limits_by_iqr_median(data_o, col, scale)

    if xxrange:
        bx1, bx2 = get_range(xxrange)
    data_n, outlier = filter_data(data_o, [bx1, bx2], col)

    print("\nStatistics summary for the input data (column=%d)" % xcol)
    print_stat(data_o, 0, "alldata")
    print_stat(data_n, 0, "filter ")

    if nstep <= 0:
        nstep = 20
    if not title:
        title = "Data distribution (entry=%d)" % (len(data_n))
    if not xlabel:
        xlabel = "data (column=%d)" % xcol
    if not ylabel:
        ylabel = "Frequency"

    fpop = "data_pop.data"  # for data populations
    shell = data_bin(bx1, bx2, nstep)
    peak = file_pop(fpop, data_n, shell, col)  # 2th column

    vtmp = peak + peak / 20.0
    xxrange = ""
    yyrange = "[0:%s]" % auto_format([vtmp])  # yrange for gnuplot (left)
    yrate = auto_format([100 * vtmp / len(data_n)])  # y2range for gnuplot (right)
    ss = '\nset y2tics \nset y2label "percentage (%%)" \nset y2range [0:%s] \n' % (yrate)
    plot = """%s  plot '%s' using 3:xtic(2) lc rgb "blue" """ % (ss, fpop)
    #    plot = """plot '%s' using 3:xtic(2) lc rgb "blue" ,'' u 0:3:4 with labels offset 0, 0.5""" %(fpop)
    bar, rot, key, style = 0, 1, 0, 0
    gnuscr, gnuout = gnu_plot(fpop, title, xxrange, yyrange, xlabel, ylabel, bar, rot, key, style, plot)

    os.system("display %s.png  &" % fpop)

    print("The histogram plot file = %s.png" % fpop)


######################################################################
def correlation_plot(infile, xcol, ycol, xlabel, ylabel, title, nstep, xxrange, yyrange):
    """plot correlations and statistics
    xcol & ycol,  the actual column number (1,2,..)
    nstep, the number of steps for box plot
    """

    print("Doing correlation and statistics for column (%d  %d)" % (xcol, ycol))
    flist = open(infile, "r").readlines()
    data_o = []  # original data
    for x in flist:
        t = x.split()
        if len(t) < xcol or len(t) < ycol or (not util.is_number(t[xcol - 1])) or (not util.is_number(t[ycol - 1])):
            continue
        data_o.append([float(t[xcol - 1]), float(t[ycol - 1])])

    bx1, bx2 = limits_by_iqr_2set(data_o, 0, scale=2.2)
    by1, by2 = limits_by_iqr_2set(data_o, 1, scale=6.0)

    if xxrange:
        bx1, bx2 = get_range(xxrange)
    if yyrange:
        by1, by2 = get_range(yyrange)

    data_n = []  # filtered data
    for x in data_o:
        if bx1 <= x[0] <= bx2 and by1 <= x[1] <= by2:
            data_n.append(x)

    print("\nStatistics summary for the input data (column=%d : %s)" % (xcol, xlabel))
    print_stat(data_o, 0, "alldata")
    print_stat(data_n, 0, "filter ")
    print("\nStatistics summary for the input data (column=%d : %s)" % (ycol, ylabel))
    print_stat(data_o, 1, "alldata")
    print_stat(data_n, 1, "filter ")

    cc_o = correlation(data_o, 0, 1)
    cc_n = correlation(data_n, 0, 1)
    print("\nCorrelations between columns (%d & %d) : cc(alldata)=%.2f , cc(filtered)=%.2f" % (xcol, ycol, cc_o, cc_n))

    if nstep <= 0:
        nstep = 20
    shell = data_bin(bx1, bx2, nstep)  # data range for x axis

    if not xlabel:
        xlabel = "data (column=%d)" % xcol
    if not ylabel:
        ylabel = "data (column=%d)" % ycol

    fbox = "box_plot.data"  # for boxplot
    data_corr(fbox, data_o, shell, 0, "box")
    xrang = "[0:%d]" % len(shell)
    yrang = ""
    title = "The adjusted boxplot"
    bar, rot, key, style = 0, 1, 0, 4  # boxplot
    plot = (
        """plot '%s' using 1:4:3:7:6:xticlabels(2) with candlesticks title 'Quartiles' whiskerbars, \
    ''         using 1:5:5:5:5 with candlesticks lt -1 notitle, \
    ''         using 1:5 with linespoints lt 3 pt 13 notitle
    """
        % fbox
    )
    gnuscr, gnuout = gnu_plot(fbox, title, xrang, yrang, xlabel, ylabel, bar, rot, key, style, plot)

    fcorr = "correlation_plot"  # for scattered plot
    #    data_corr(fcorr,data_o, shell, 0, '')
    xext = (bx2 - bx1) / 20.0
    yext = (by2 - by1) / 20.0
    xrang = "[%.3f:%.3f]" % (bx1 - xext, bx2 + xext)
    yrang = "[%.3f:%.3f]" % (by1 - yext, by2 + yext)
    title = "The scatter plot (data size=%s)" % len(data_n)
    bar, rot, key, style = 0, 1, 0, 2  # scatt

    plot = "\nf(x) = a*x + b\n fit f(x) '%s' u %d:%d via a, b\n title_f(a,b) = sprintf('f(x) = %%.2fx + %%.2f', a, b) \nplot '%s' using %d:%d , f(x) t title_f(a,b) " % (
        infile,
        xcol,
        ycol,
        infile,
        xcol,
        ycol,
    )
    gnuscr, gnuout = gnu_plot(fcorr, title, xrang, yrang, xlabel, ylabel, bar, rot, key, style, plot)

    if os.path.exists("fit.log"):
        print("\nThe linear fit f(x)=a*x + b ")
        fp = open("fit.log", "r")
        for x in fp:
            if "a " in x[:3] and " +/- " in x and "(" in x:
                a = x
            if "b " in x[:3] and " +/- " in x and "(" in x:
                b = x
        print("%s %s" % (a, b))
    print("\nThe output: %s.png ; %s.png \n" % (fbox, fcorr))


##########################################################
def correlation(data, col1, col2):
    """
    data = [[?,?], [?,?] ...]
    Pearson correlation: sum(xy-n*<x><y>)/((n-1)*sx*sy)
    sx and sy are the sample standard deviations.
    <x> and  <y> are the average of x and y.
    """

    xavg, xdev, xmin, xmax = util.mean_dev(data, col1)
    yavg, ydev, ymin, ymax = util.mean_dev(data, col2)

    #    print xavg, xdev, xmin, xmax, yavg, ydev, ymin, ymax

    n = len(data)
    if n < 4:
        print("Error: too few data ")
        return 0

    xy = 0
    for x in data:
        xy = xy + x[0] * x[1]

    cc = (xy - n * xavg * yavg) / ((n - 1) * xdev * ydev)

    return cc


######################################################################
def print_stat(data, coln, idd):
    """data=[[],[]]; coln=0, 1;  idd='alldata' or filter"""

    xavg, xdev, xmin, xmax = util.mean_dev(data, coln)
    xskew, xkurt = skewness_kurtosis(data, coln, xavg, xdev)

    print("N_%s=%d; mean=%.3f; dev=%.3f; min=%.3f; max=%.3f; skew=%.2f; kurtosis=%.2f" % (idd, len(data), xavg, xdev, xmin, xmax, xskew, xkurt))


######################################################################
def data_corr(fcorr, fp, shell, col, id):
    """get data correlations between column of col and col+1
    fp: is  a list of list; fcorr is a out filename,
    shell is the value of column (col) break into bin.
    """

    scale = 1.5
    fw = open(fcorr, "w")
    fp.sort(key=lambda y: y[col])

    nt = len(fp)
    nstep = len(shell)

    s1 = "#data_range  avg_bin  number  mean  deviation  minimum  maximum\n"
    if id == "box":
        s1 = "#n   xtics   low_BD    Q1     median   Q3   upp_BD  Number\n"
    fw.write(s1)

    j = 0
    for i in range(nstep):
        if i == 0:
            continue

        bin = (shell[i - 1] + shell[i]) / 2.0
        n, subg = 0, []
        for k in range(j, nt):
            if fp[k][col] >= shell[i]:
                j = k
                break
            elif shell[i] > fp[k][col] >= shell[i - 1]:
                subg.append(fp[k][col + 1])  # data after resolution
                n = n + 1
        nlen = len(subg)
        if nlen:
            if id == "box":
                #                lowb, q1, med, q3, uppb=pstat.get_boxplot(subg, -1, scale)
                #                lowb, q1, med, q3, uppb=pstat.get_adjusted_boxplot_mc(subg, -1, scale=1.5)
                lowb, q1, med, q3, uppb = get_adjusted_boxplot(subg, -1, scale)
                ss = auto_format([i, bin, lowb, q1, med, q3, uppb, nlen])
                fw.write("%s \n" % ss)
            else:
                avg, dev, mini, maxi = util.mean_dev(subg, -1)
                s = '"%.3f %.3f"' % (shell[i - 1], shell[i])
                ss = auto_format([s, bin, n, avg, dev, mini, maxi])
                fw.write("%s \n" % ss)

    fw.close()


######################################################################
def get_range(range_in):
    bx1 = float(range_in.split()[0])
    bx2 = float(range_in.split()[1])
    return bx1, bx2


######################################################################
def scatter_plot(file, col1, col2):
    """given a file (list of data, and the col1 & col2), a histgram/scatter can be plotted."""

    flist = open(file, "r").readlines()
    data1 = []
    for x in flist:
        t = x.split()
        if col1 >= len(t) or col2 >= len(t):
            # print 'Error: the given column > the column in file'

            continue
        if not util.is_number(t[col1]) or not util.is_number(t[col2]):
            continue

        data1.append([float(t[col1]), float(t[col2])])

    #    print 'all data = ', len(data1)
    #    for x in data1: print x
    # col = 0

    xavg0, xdev0, xmin0, xmax0 = util.mean_dev(data1, col=0)
    xavg, xdev, xmin, xmax = util.mean_dev(data1, col=1)

    fcc = "data_cc.data"
    scat = "data_cc.scat"
    fw = open(fcc, "w")
    for x in data1:
        fw.write("%s %s \n" % (x[0], x[1]))
    fw.close()

    title = "The scatter plot of column (%d and %d : entry=%d)" % (col1, col2, len(data1))  #
    xrange, yrange, xlabel, ylabel = "", "", "(column=%d)" % col1, "(column=%d)" % col2
    bar, rot, key, style = 0, 1, 0, 2
    plot = "plot '%s' using 1:2 " % (fcc)
    gnuscr, gnuout = gnu_plot(scat, title, xrange, yrange, xlabel, ylabel, bar, rot, key, style, plot)

    print("column1=%d : mean=%.3f; dev=%.3f; min=%.3f; max=%.3f" % (col1, xavg0, xdev0, xmin0, xmax0))
    print("column2=%d : mean=%.3f; dev=%.3f; min=%.3f; max=%.3f" % (col2, xavg, xdev, xmin, xmax))

    os.system("display %s.gif  &" % scat)


######################################################################
def filter_data(data_in, shell, col):
    """filter the data according to the shell & col values."""

    data_new, outlier = [], []
    for x in data_in:
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
    shell.append(a)
    #    print lowv,upv,a, shell, '<br>'
    return shell


######################################################################
def auto_format(data):
    """data is a single list [ ,,]; return string"""
    ss = []
    for x in data:
        if util.is_number(x):
            if math.fabs(x) < 1:
                ss.append(" %.3f " % x)
            elif 1.0 <= math.fabs(x) < 10:
                ss.append(" %.2f " % x)
            elif 10 <= math.fabs(x) < 100.0:
                ss.append(" %.1f " % x)
            elif math.fabs(x) >= 100.0:
                ss.append(" %.0f " % (x + 0.0001))
        else:
            ss.append(" %s " % x)

    return " ".join(ss)


######################################################################
def data_in_shell(data, shell, col):
    """
    data: a input list of list. shell, the shell values; col the data col.
    return new data list.
    """

    data.sort(key=lambda y: y[col])  # sort fp by the col

    nt = len(data)
    nstep = len(shell)
    j = 0
    nnt = 0
    data_inshell = []
    for i in range(nstep):
        if i == 0:
            continue
        bin = (shell[i - 1] + shell[i]) / 2.0
        n = 0
        for k in range(j, nt):
            if data[k][col] > shell[i]:
                j = k
                break
            elif shell[i] >= data[k][col] >= shell[i - 1]:
                n = n + 1
                nnt = nnt + 1

        p = 100 * float(n) / nt
        #        print ['"%.3f %.3f"' %(shell[i-1], shell[i]), bin, n, p]
        data_inshell.append(['"%.3f %.3f"' % (shell[i - 1], shell[i]), bin, n, p])

    return data_inshell


######################################################################
def file_pop(fpop, data, shell, col):
    """get the data populations using the col of the data [[],[]..]
    data must be sorted!
    """

    nt = len(data)
    if nt == 0:
        print("<br>Error: No data was found")
        return

    fw = open(fpop, "w")
    s1 = "#data_range  avg_bin   entry_number  percentage\n"
    fw.write(s1)

    data_inshell = data_in_shell(data, shell, col)
    mm = []
    for x in data_inshell:
        ss = auto_format([x[0], x[1], x[2], x[3]])
        mm.append(x[2])
        fw.write("%s\n" % (ss))

    fw.close()
    return max(mm)


######################################################################
def limits_by_iqr_2set(datain, col, scale):
    """get the boundaries [LB, UB] from the two data sets [-, 30] and
    [70, -]:  This could be a problem if too many wrong values.
    """

    datain.sort(key=lambda y: y[col])
    nt = len(datain)
    trim = 0.20  # use 20% of data on both sides
    n1, n3 = int(trim * nt), int((1 - trim) * nt)

    data1 = datain[:n1]
    data2 = datain[n3:]

    bb1, bb2, qb2 = limits_by_iqr(data1, col, scale)
    ba1, ba2, qa2 = limits_by_iqr(data2, col, scale)
    #   return bb1,ba2
    b1, b2 = adjust_outlier(bb1, ba2, datain, col)
    # print  b1,b2, bb1,ba2
    return b1, b2


######################################################################
def adjust_outlier(b1, b2, data, col):
    bb1, bb2 = data[0][col], data[-1][col]
    n = len(data)
    trim = 0.0001

    x1, x2 = b1, b2
    if b1 < bb1:  # over look outliers
        n1 = int(trim * n)
        x1 = data[n1][col]

    if b2 > bb2:
        n2 = int((1 - trim) * n)
        x2 = data[n2][col]

    return x1, x2


######################################################################
def limits_by_iqr_median(datain, col, scale):
    """get the boundaries [LB, UB] from the two data sets [-, median] and
    [median, -]
    """

    maxi = get_median(datain, col)
    data1, data2 = [], []
    for z in datain:  # use the raw data
        if z[col] <= maxi:
            data1.append(z)
        if z[col] >= maxi:
            data2.append(z)

    bb1, bb2, qb2 = limits_by_iqr(data1, col, scale)
    ba1, ba2, qa2 = limits_by_iqr(data2, col, scale)

    return bb1, ba2  # the [LB, UB]


######################################################################
def get_median(datain, col):
    """get the median value:
    data is a list [[],[], ] if col>=0; is single list [ ], if col<0
    """

    data = datain
    if col >= 0:  #
        data = []
        for x in datain:
            data.append(x[col])

    n = len(data)
    if n == 0:
        return 0
    if n == 1:
        return data[0]
    util.gsort(data, -1, 0)

    if n % 2 == 0:  # even number
        med = (data[n / 2] + data[n / 2 - 1]) / 2.0
    else:
        med = data[n / 2]

    return med


######################################################################
def limits_by_iqr(alist, col, scale):
    """Using the quartile method (IQR) to select the outliers
    col: the column of [[col,..], [...]]
    scale : a scale factor (1.5, 2.2, 3.0). scale=2.22 is equivalent to 3 standard deviations
    """

    if not alist:
        return 0, 0, 0

    iqr, q1, q2, q3 = get_quarts(alist, col)

    b1 = q1 - iqr * scale
    b2 = q3 + iqr * scale

    return b1, b2, q2


######################################################################
def get_quarts(alist, col):
    """
    alist : a single list (col<0) or list of list (col>=0).

    """

    nt = len(alist)
    n1 = int(nt / 4.0)
    n2 = int(nt / 2.0)  # the median
    n3 = int(nt * 3 / 4.0)

    if col < 0:
        util.gsort(alist, -1, 0)
        q1, q2, q3 = float(alist[n1]), float(alist[n2]), float(alist[n3])
    else:
        alist.sort(key=lambda y: y[col])  # sort fp by the col
        q1, q2, q3 = float(alist[n1][col]), float(alist[n2][col]), float(alist[n3][col])

    iqr = q3 - q1

    return iqr, q1, q2, q3


######################################################################
def get_adjusted_boxplot(alist, col, scale):
    """
    if col>0: the alist is [[col,..], [...]]; else: [....]
    """

    if not alist:
        return 0, 0, 0, 0, 0

    iqr, q1, med, q3 = get_quarts(alist, col)

    nt = len(alist)
    trim = 0.25  # use 10% of data on both sides
    n1, n3 = int(trim * nt), int((1 - trim) * nt)
    #    print iqr, q1,med,q3, n1,n3, nt

    data1 = alist[:n1]
    data2 = alist[n3:]
    #    print data1, '<br>',data2, '<br>'

    bb1, bb2, qb2 = limits_by_iqr(data1, col, scale)
    ba1, ba2, qa2 = limits_by_iqr(data2, col, scale)

    return bb1, q1, med, q3, ba2


############################################################
def limits_by_iqr_mc(alist, col, scale):
    """Using the quartile method (IQR) to select the outliers
    col: the column of [[col,..], [...]]
    scale : a scale factor (1.5, 2.2, 3.0)
    mc : the medcouple value (robust skewness)
    """

    mc = get_mc(alist, col)

    scale = 1.5
    if not alist:
        return 0, 0, 0

    if mc >= 0:
        left = math.exp(-3.5 * mc)
        right = math.exp(4 * mc)
    else:
        left = math.exp(-4.0 * mc)
        right = math.exp(3.5 * mc)

    iqr, q1, q2, q3 = get_quarts(alist, col)
    b1 = q1 - iqr * scale * left
    b2 = q3 + iqr * scale * right
    return b1, b2, q2


######################################################################
def get_mc(datain, col):
    """This is a robust alternative to classical skewness introduced by
    Brys et al. (2003), Another type of skewness is the medcouple (MC).

    MC = med (h( xi , xj )) , where the kernel function h is given by:

    h( xi , xj ) = ( xj + xi - 2*med_k )/(xj - xi)

    , where med_k is the median of data, and i and j have to satisfy
    xi <= med_k <= xj , and xi != xj . The value of the MC ranges between
    -1 and 1. If MC=0,the data is symmetric. If MC>0, the data has a right
    skewed distribution, whereas if MC<0, the data has a left skewed distribution.

    This method is proved to be useless for big data (>40000).
    The matrix takes too much memory!! natrux N*N

    datain : a aingle or a list of list [[..], [..]..];  col: the column of [..]
    """

    med = get_median(datain, col)
    data1, data2 = [], []  # separate the data
    if col < 0:  # single list
        for z in datain:
            if z <= med:
                data1.append(z)
            if z >= med:
                data2.append(z)
    else:  # list of list
        for z in datain:
            if z[col] <= med:
                data1.append(z[col])
            if z[col] >= med:
                data2.append(z[col])

    #    print 'med=', med  , len(data1),  len(data2)
    data = []
    for x in data1:  # make the matrix
        for y in data2:
            if x == y:
                continue
            mc1 = (x + y - 2 * med) / (y - x)
            #            print x, y, mc1
            data.append(mc1)
    util.gsort(data, -1, 0)
    mc = get_median(data, -1)

    return mc


############################################################


def skewness_kurtosis(data, col, xavg, xdev):
    """data = [], one column;  return skewness & kurtosis
    skewness_q = ((q3-md) - (md-q1))/(q3-q1) where md=median, q1=first quart
    or
    skewness_p = ((p90-md) - (md-p10))/(p90-p10) where md=median, p10=10%

    kurtosis = 0.5*(q3-q1)/(p90-p10)
    kurtosis>0.263, lepokurtic; ~ 0.263, mesokurtic;  <0.263 platy kurtic

    Other definations: http://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
    skewness = (sum(xi-xavg)^3)/((N-1)*dev^3)

    kurtosis = (sum(xi-xavg)^4)/((N-1)*dev^4)
    kurtosis = 3 for normal

    """

    if xdev == 0:
        return 0, 0

    n = len(data)
    # m = int(n * 0.5)
    # nq1, nq3 = int(n * 0.15), int(n * 0.75)
    # np10, np90 = int(n * 0.10), int(n * 0.90)

    # md = data[m][col]
    # q1, q3 = data[nq1][col], data[nq3][col]
    # p10, p90 = data[np10][col], data[np90][col]

    # skewness_q, skewness_p, kurtosis_p = 0, 0, 0

    # if (p90 - p10) != 0 and (q3 - q1) != 0:
    #     skewness_q = ((q3 - md) - (md - q1)) / (q3 - q1)
    #     skewness_p = ((p90 - md) - (md - p10)) / (p90 - p10)
    #     kurtosis_p = 0.5 * (q3 - q1) / (p90 - p10)

    #    print m,  nq1, nq3, np10, np90, q1, q3, p10, p90, md , '<br>'

    skew, kurt = 0, 0
    for x in data:
        d = x[col] - xavg
        skew = skew + d**3
        kurt = kurt + d**4

    skewness = skew / ((n - 1) * xdev**3)
    kurtosis = kurt / ((n - 1) * xdev**4)

    #    print 'skewness, kurtosis ', skewness, kurtosis, skewness_q, skewness_p, kurtosis_p

    return skewness, kurtosis


############################################################
def get_pop_corr(infile, nn, narg, arg):
    """nn, the nth argument; narg, the maximum argument;  arg, the argument"""

    xxrange, yyrange, title, xlabel, ylabel, xcol, ycol, nstep = "", "", "", "", "", 0, 0, 0

    ss = ""
    for i in range(nn, narg):
        ss = ss + "%s " % arg[i]

    ss_split = ss.split(",")
    for x in ss_split:
        s = x.strip().split("=")
        key = s[0].strip().lower()

        if ("xcol" == key or "ycol" == key) and not util.is_number(s[1]):
            print('Error: not an integer after xcol or ycol. both must be separated by ",".')
            return

        if "xrange" == key:
            xxrange = s[1]
        elif "yrange" == key:
            yyrange = s[1]
        elif "title" == key:
            title = s[1]
        elif "xlabel" == key:
            xlabel = s[1]
        elif "ylabel" == key:
            ylabel = s[1]
        elif "xcol" == key:
            xcol = int(s[1])
        elif "ycol" == key:
            ycol = int(s[1])
        elif "nstep" == key:
            nstep = int(s[1])

    print("Input file=%s" % infile)

    if xcol > 0 and ycol > 0:  # do correlations
        correlation_plot(infile, xcol, ycol, xlabel, ylabel, title, nstep, xxrange, yyrange)
    elif xcol > 0:
        histogram_plot(infile, xcol, xlabel, ylabel, title, nstep, xxrange)


##########################################################
def gnu_plot(datafile, title, xxrange, yrange, xlabel, ylabel, bar, rot, key, style, plot):
    """some head infor for gnuplot (one data set plot by default)
        bar>0: add error bar on the histogram
        rot>0, xlabel rotation;
        key=0: key>0, on;  key=2,left;
        style=0, histogram ; =1 linepoints; =2 point; =3 dot; =4 boxplot;

        ------

    lt=linetypes color      pt

    -1      black           -1      n/a
    0       black           0       dotted
    1       red             1       +
    2       green           2       x
    3       blue            3       *
    4       magenta         4       empty square
    5       cyan            5       filled square
    6       brown           6       empty circle
    7       light brown     7       filled circle
    8       light red       8       empty triangle
            -----------
    9       red             9       filled triangle
    10      green           10      empty nabla
    11      blue            11      filled nabla
    12      magenta         12      empty rhombus
    13      cyan            13      filled rhombus
                                    --------------
    14      brown           14      +
    15      light brown     15      x

    example:
    set palette
    plot sin(x) lt -1

    """

    psfix = "png"  # must be the same as the dic

    gnuscr = datafile + ".gnu"
    gnuout = datafile + ".%s" % psfix
    fw = open(gnuscr, "w")

    x1 = "set terminal %s large size 840,640 #1000,740 \n" % psfix
    x2 = "set output '%s'\n" % gnuout
    fw.write(x1 + x2)
    fw.write('set title "%s" \n' % title)
    fw.write("set datafile missing '.'\n")
    fw.write("set style fill  solid 1.0 border -1\n")

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
    elif style == 4:  # boxplot
        fw.write("set style fill empty\n")
        fw.write("set bars 2.0\nset boxwidth  0.5\n")
    elif style == 5:  # add linear regression
        fw.write("f(x) = a*x + b\n")
        fw.write("fit f(x) '%s' using '1:2' via a, b\n" % datafile)

    if len(xxrange) > 1:
        fw.write("set xrange %s\n" % xxrange)
    if len(yrange) > 1:
        fw.write("set yrange %s\n" % yrange)

    y2range, y2label = "", ""  # disabled at the moment
    if len(y2range) and len(y2label):
        fw.write("set y2tics\n")
        fw.write("set y2range %s\n" % yrange)
        fw.write('set y2label "%s" \n' % y2label)

    if key == 0:
        fw.write("set key off\n")
    else:
        fw.write("set key on\n")
        if key == 2:
            fw.write("set key reverse left Left\n")

    fw.write("set grid ytics\n")
    fw.write("set size 1.0,1.0\n")
    #    fw.write("set autoscale\n")

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

    #    fw.write("set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 bdefault\n")

    fw.write(plot)

    fw.write("""\n\n# plot 'file' u 1:3:xtic(1) t "?" with linespoints,'' using 1:2 t "?" with boxes\n""")
    fw.write("""# plot 'filename' using 1:2:3:xtic(1) t "all" with  errorbars \n""")
    fw.write("# set term x11 \n")
    fw.write("# replot \n")
    fw.write("# pause -6 \n")

    #    fw.write("set term x11\nreplot \npause 4\n")

    fw.close()
    arg = "gnuplot  %s >&/dev/null" % (gnuscr)
    os.system(arg)  # if not given /bin/csh -c, use B shell

    return gnuscr, gnuout


############################################################
def add_percentile_rank(infile, coln):
    """Add percentile after the column (col).
    This is the same as PERCENTRANK.INC(array,x,sig) in Microsoft Excel(differ PERCENTRANK.EXC)
    """

    print("Add percent_rank to the %dth column (same as PERCENTRANK.INC in Microsoft Excel)" % coln)

    col = coln - 1
    if col <= 0:
        print("Error: You selected wrong column. (The column should >0!)")
        return

    outf = "%s_new" % infile
    fw = open(outf, "w")

    flist = open(infile, "r").readlines()
    ulist = []  # must be a string!
    for x in flist:
        ss = x.strip().split()
        if not ss or len(ss) < coln or not ss[0].isdigit or "?" in ss[col]:
            continue
        ulist.append(ss[col])

    #    uniq=list(sorted(ulist, key=float))  # not work for string_number
    uniq = sorted_copy(ulist)
    #    uniq = sort_string_number(ulist) #not tested.
    #    sys.exit()

    nt = len(uniq)
    dic = {}
    for i, x in enumerate(uniq):  # same as excel
        p = (i) / float(nt - 1)
        if x not in dic:
            dic[x] = p

    for x in flist:  # add percentile rank to col
        ss = x.strip().split()
        if not ss:
            continue
        if "?" in ss[col]:
            ln = ln = "  ".join(ss[: col + 1]) + " ? " + "  ".join(ss[col + 1 :])
        else:
            prank = 100 * (1 - float(dic[ss[col]]))
            ln = "  ".join(ss[: col + 1]) + "  %6.1f " % (prank) + "  ".join(ss[col + 1 :])

        s1 = ln.split()
        tt = ""
        for x in s1:
            tt = tt + " %6s " % x.strip()
        fw.write("%s\n" % tt)
    #     print tt

    print("The new file = %s\n" % outf)
    fw.close()


############################################################
def sorted_copy(alist):
    """sort string in enumerical order. alist=['3.5', '3 str'...]"""
    indices = list(map(_generate_index, alist))
    decorated = list(zip(indices, alist))
    decorated.sort()
    return [item for index, item in decorated]


############################################################
def _generate_index(str):
    """
    Splits a string into alpha and numeric elements, which
    is used as an index for sorting
    """

    index = []  # the index is built progressively

    def _append(fragment, alist=index):
        if fragment.isdigit():
            fragment = int(fragment)
        alist.append(fragment)

    # initialize loop
    prev_isdigit = str[0].isdigit()
    current_fragment = ""
    # group a string into digit and non-digit parts
    for char in str:
        curr_isdigit = char.isdigit()
        if curr_isdigit == prev_isdigit:
            current_fragment += char
        else:
            _append(current_fragment)
            current_fragment = char
            prev_isdigit = curr_isdigit
    _append(current_fragment)
    return tuple(index)


############################################################
def sort_string_number(data):
    """sort int string in enumerical order. data = ['?', '2', '3','20'...]"""
    return sorted(data, key=lambda item: (int(item.partition(" ")[0]) if item[0].isdigit() else float("inf"), item))
