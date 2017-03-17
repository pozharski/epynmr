from windows import viewwindow, peakwindow, titrwindow, dualwindow
from matplotlib.pyplot import show
from scipy import array, arange, sqrt, log
from scipy.stats import scoreatpercentile
from nmrio import hsqc, peakset, bonferroni, iqr_sigma
import os, sys

def viewhsqc(args):
    viewwindow(args.input_file)
    show()

def dualhsqc(args):
    dualwindow(args.xargs[0], args.xargs[1], args.peakfile)
    show()

def peakhsqc(args):
    peakwindow(args.input_file, args.num_peaks, [args.bleft, args.bright, args.bbottom, args.btop])
    show()

def titrhsqc(args):
    titrwindow(args.aponum, eval(args.holonums), args.num_peaks, [args.bleft, args.bright, args.bbottom, args.btop])
    show()

def autocros(args):
    apo = hsqc(os.path.join(args.folder,str(args.aponum)+".nv"))
    lrbt = [args.bleft, args.bright, args.bbottom, args.btop]
    apo_peakdata = apo.peak_search(args.num_peaks, lrbt)
    apo_peaks = peakset(apo_peakdata)
    pmv = []
    for sample in eval(args.controls):
        ctrl_apo = hsqc(os.path.join(args.folder,str(sample)+".nv"))
        ctrl_peaks = peakset(ctrl_apo.peak_search(args.num_peaks, lrbt))
        pmv.append(apo_peaks.matchsignal(ctrl_apo,True))
    pmv = array(pmv)
    apo_good = []
    for (i, value) in enumerate(pmv.min(0)):
        if value>args.vcutoff:
            apo_good.append(apo_peakdata[i])
    print "Found %d good peaks in set of %d" % (len(apo_good),160)
    good_peaks = peakset(array(apo_good))
    for sample in eval(args.holonums):
        holo = hsqc(os.path.join(args.folder,str(sample)+".nv"))
        x = good_peaks.matchsignal(holo,True)
        print "Sample #" + str(sample).ljust(3) + ": " + " ".join(map(lambda x : "%5.2f" % x, map(lambda p : scoreatpercentile(x, p), range(0, 100+10, 10))))

def peakcros(args):
    with open(args.peakfile) as fin:
        peaks=map(lambda t : t.split(), fin)
    peaks = peakset(map(lambda t : [float(t[0]), float(t[1]), t[2]], peaks))
    lrbt = [args.bleft, args.bright, args.bbottom, args.btop]
    for sample in eval(args.holonums):
        holo = hsqc(os.path.join(args.folder,str(sample)+".nv"))
        holopeaks = peakset(holo.xyzconv(holo.peak_search(args.num_peaks, lrbt)))
        print "Sample #" + str(sample).ljust(3) + ": %d peaks matched" % holopeaks.matchcount(peaks, args.dcutoff, args.scutoff)

def peakmatch(args):
    with open(args.peakfile) as fin:
        peaks=map(lambda t : t.split(), fin)
    Np = len(peaks)
    peaks = peakset(map(lambda t : [float(t[0]), float(t[1]), t[2]], peaks))
    lrbt = [args.bleft, args.bright, args.bbottom, args.btop]
    allpms, hps = [], []
    for sample in eval(args.holonums):
        holo = hsqc(os.path.join(args.folder,str(sample)+".nv"))
        holopeaks = peakset(holo.xyzconv(holo.peak_search(args.num_peaks, lrbt)))
        samplenum = str(sample).ljust(3)
        print "Sample #" + samplenum + " processed"
        pms = holopeaks.peakmatch(peaks)
        allpms.append(list(pms[1]))
        hps.append(holopeaks)
    x = array(allpms).flatten()
    x.sort()
    N = len(x)
    n = arange(1,N).astype(float)
    sigmaD = sqrt(0.5*x[:N/2]**2/log(N/(N-n[:N/2]))).mean()
    print "Estimated peak position noise level %.3f" % ( sigmaD )
    zcut = bonferroni(args.zcutoff, Np)
    print "Requested Z-score cutoff (%.2f) adjusted to %.2f (Bonferroni correction)" % (args.zcutoff, zcut)
    dcut = zcut*sigmaD
    print "Peak mismatch cutoff set to %.3f ppm" % dcut
    for (i,sample) in enumerate(eval(args.holonums)):
        holopeaks = hps[i]
        pms = holopeaks.peakmatch(peaks)
        samplenum = str(sample).ljust(3)
        print "Sample #" + samplenum + ": %d peaks shifted" % (Np-holopeaks.matchcount(peaks, dcut, args.scutoff))
        if args.print_peakmatch:
            pms = sorted(zip(*(pms+[list(zip(*peaks.peaks)[2])])),key=lambda x : x[1],reverse=True)
            print '\n'.join(map(lambda x : samplenum+"%7.3f %6.2f %s" % x[1:], filter(lambda xx : xx[1]>=dcut,pms)))
        
def peakmxy(args):
    with open(args.peakfile) as fin:
        peaks=map(lambda t : t.split(), fin)
    Np = len(peaks)
    peaks = peakset(map(lambda t : [float(t[0]), float(t[1]), t[2]], peaks))
    lrbt = [args.bleft, args.bright, args.bbottom, args.btop]
    allxy, hps = [], []
    for sample in eval(args.holonums):
        holo = hsqc(os.path.join(args.folder,str(sample)+".nv"))
        holopeaks = peakset(holo.xyzconv(holo.peak_search(args.num_peaks, lrbt)))
        hps.append(holopeaks)
        allxy.append(holopeaks.peakmxy(peaks))
        samplenum = str(sample).ljust(3)
        sys.stderr.write("Sample #" + samplenum + " processed\n")
    xx, yy = array(allxy).T
    sx, sy = iqr_sigma(xx), iqr_sigma(yy)
    sys.stderr.write("Estimated peak position variations %.3f %.3f\n" % ( sx, sy ))
    zcut = bonferroni(args.zcutoff, Np)
    sys.stderr.write("Requested Z-score cutoff (%.2f) adjusted to %.2f (Bonferroni correction)\n" % (args.zcutoff, zcut))
    xcut, ycut = zcut*sx, zcut*sy
    sys.stderr.write("Peak mismatch cutoffs set to %.3f / %.3f ppm\n" % (xcut, ycut))
    for (i,sample) in enumerate(eval(args.holonums)):
        holopeaks = hps[i]
        pms = holopeaks.peakmatch(peaks)
        Nm, indm, xm, ym = holopeaks.matchcut_xy(peaks, xcut, ycut)
        samplenum = str(sample).ljust(3)
        print "Sample #" + samplenum + ": %d peaks shifted" % Nm
        if args.print_peakmatch:
            for (j, flag) in enumerate(indm):
                if flag or args.print_full:
                    print samplenum + " %7.3f %7.3f %s %7.3f %7.3f" % (holopeaks.peaks[pms[0][j]][0], holopeaks.peaks[pms[0][j]][1], peaks.peaks[j][2], xm[j], ym[j])
