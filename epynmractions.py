from windows import viewwindow, peakwindow, titrwindow, dualwindow
from matplotlib.pyplot import show
from scipy import array
from scipy.stats import scoreatpercentile
from nmrio import hsqc, peakset
import os

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
    peaks = peakset(map(lambda t : [float(t[0]), float(t[1]), t[2]], peaks))
    lrbt = [args.bleft, args.bright, args.bbottom, args.btop]
    for sample in eval(args.holonums):
        holo = hsqc(os.path.join(args.folder,str(sample)+".nv"))
        holopeaks = peakset(holo.xyzconv(holo.peak_search(args.num_peaks, lrbt)))
        print "Sample #" + str(sample).ljust(3) + ": %d peaks matched" % holopeaks.matchcount(peaks, args.dcutoff, args.scutoff)
        pms = sorted(zip(*(holopeaks.peakmatch(peaks)+[list(zip(*peaks.peaks)[2])])),key=lambda x : x[1],reverse=True)
        print '\n'.join(map(lambda x : "%7.3f %6.2f %s" % x[1:], pms))
        
