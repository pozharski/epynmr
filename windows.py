from nmrio import hsqc, peakset
from matplotlib.figure import Figure
from matplotlib.pyplot import figure, gca, show, grid, close
from scipy import logspace, log10, delete, array
import sys

def file_or_data(data):
    return data if data.__class__ is hsqc else hsqc(data)
    
def viewwindow(data, *args, **kwds):
    fig = figure(FigureClass=ViewWindow)
    fig.set_hsqc(file_or_data(data))
    fig.plot()
    fig.canvas.mpl_connect('key_press_event', fig.onkeypress)
    show()

def dualwindow(data1, data2, *args, **kwds):
    fig = figure(FigureClass=DualWindow)
    fig.set_hsqc(file_or_data(data1), file_or_data(data2))
    fig.plot()
    fig.canvas.mpl_connect('key_press_event', fig.onkeypress)
    show()

def peakwindow(data, num=50, box=None, *args, **kwds):
    fig = figure(FigureClass=PeakWindow)
    fig.set_hsqc(data if data.__class__ is hsqc else hsqc(data))
    fig.set_max_peaknum(num)
    fig.set_box(box)
    fig.peak_search()
    fig.plot()
    fig.canvas.mpl_connect('key_press_event', fig.onkeypress)
    show()

def titrwindow(aponum, holonums, peaknum, box=None, *args, **kwds):
    fig = figure(FigureClass=TitrWindow)
    fig.set_max_peaknum(peaknum)
    fig.set_box(box)
    fig.set_apo_hsqc(str(aponum)+".nv")
    fig.load_apo_hsqc()
    fig.set_holo_range(holonums)
    fig.peak_search()
    fig.plot()
    fig.canvas.mpl_connect('key_press_event', fig.onkeypress)
    show()

class NMRWindow(Figure):
    def __init__(self, *args, **kwds):
        Figure.__init__(self, *args, **kwds)
        self.zfloor = 1
    def log_contour_levels(self, z):
        return logspace(log10(z.mean()+self.zfloor*z.std()),log10(z.max()), 10)
    def onkeypress(self, event):
        if event.key == 'q':
            close(self)
            sys.exit()

class ViewWindow(NMRWindow):
    def set_hsqc(self, data):
        self.dataset = data
        self.set_contour_levels()
    def set_contour_levels(self):
        self.V = self.log_contour_levels(self.dataset.nv_data)
    def plot(self):
        z = self.dataset.nv_data
        self.add_subplot(111)
        self.axes[0].clear()
        self.contur = self.axes[0].contour(z, self.V, colors='blue')
        grid(True)
    def recontur(self, event):
        for coll in self.contur.collections:
            gca().collections.remove(coll)
        self.contur = self.axes[0].contour(self.dataset.nv_data, self.V)
        event.canvas.draw()
    def onkeypress(self, event):
        NMRWindow.onkeypress(self, event)
        if event.key == 'pageup':
            self.zfloor += 1
        if event.key == 'pagedown':
            self.zfloor = max(1, self.zfloor-1)
        self.set_contour_levels()
        self.recontur(event)

class DualWindow(NMRWindow):
    def __init__(self, *args, **kwds):
        Figure.__init__(self, *args, **kwds)
        self.zfloor = 1
        self.contur = []
    def set_hsqc(self, data1, data2):
        self.datasets = [ data1, data2 ]
        self.set_contour_levels()
    def set_contour_levels(self):
        self.V = map(lambda x : self.log_contour_levels(x.nv_data),  self.datasets)
    def plot(self):
        self.add_subplot(111)
        self.axes[0].clear()
        package = zip(*([self.datasets, self.V, ['red', 'blue']]))
        self.contur = map(lambda x : self.axes[0].contour(x[0].nv_data, x[1], colors=x[2]), package)
        grid(True)
    def recontur(self, event):
        for contur in self.contur:
            for coll in contur.collections:
                gca().collections.remove(coll)
        package = zip(*([self.datasets, self.V, ['red', 'blue']]))
        self.contur = map(lambda x : self.axes[0].contour(x[0].nv_data, x[1], colors=x[2]), package)
        event.canvas.draw()
    def onkeypress(self, event):
        NMRWindow.onkeypress(self, event)
        if event.key == 'pageup':
            self.zfloor += 1
        if event.key == 'pagedown':
            self.zfloor = max(1, self.zfloor-1)
        self.set_contour_levels()
        self.recontur(event)

class PeakWindow(Figure):
    def __init__(self, *args, **kwds):
        Figure.__init__(self, *args, **kwds)
        self.xzoom = 10
        self.yzoom = 10
        self.zfloor = 1
    def set_max_peaknum(self, peaknum):
        self.peaknum = peaknum
    def set_box(self, box=None):
        self.box = box
    def set_hsqc(self, data):
        self.dataset = data
        z = self.dataset.nv_data
        self.V = logspace(log10(z.mean()+self.zfloor*z.std()),log10(z.max()), 10)
    def peak_search(self):
        self.peaks = self.dataset.peak_search(self.peaknum, self.box)
        self.pcounter = 0
    def plot(self):
        z = self.dataset.nv_data
        self.add_subplot(111)
        self.axes[0].clear()
        self.contur = self.axes[0].contour(z, self.V, colors='blue')
        self.peakmarks = self.axes[0].plot(self.peaks.T[0],self.peaks.T[1],'r+')[0]
        grid(True)
    def recontur(self, event):
        for coll in self.contur.collections:
            gca().collections.remove(coll)
        self.contur = self.axes[0].contour(self.dataset.nv_data, self.V)
        event.canvas.draw()
    def limit_pcounter(self):
        self.pcounter = max(0, self.pcounter)
        self.pcounter = min(len(self.peaks)-1, self.pcounter)
    def limit_zoom(self):
        self.xzoom = max(10, self.xzoom)
        self.yzoom = max(10, self.yzoom)
    def onkeypress(self, event):
        if event.key == 'pageup':
            self.zfloor += 1
            z = self.dataset.nv_data
            self.V = logspace(log10(z.mean()+self.zfloor*z.std()),log10(z.max()), 10)
            self.recontur(event)
        if event.key == 'pagedown':
            self.zfloor = max(1, self.zfloor-1)
            z = self.dataset.nv_data
            self.V = logspace(log10(z.mean()+self.zfloor*z.std()),log10(z.max()), 10)
            self.recontur(event)
        if event.key == 'home':
            self.pcounter = 0
            self.onzoom(event)
        if event.key == '+':
            self.pcounter += 1
            self.limit_pcounter()
            self.onzoom(event)
        if event.key == '-':
            self.pcounter -= 1
            self.limit_pcounter()
            self.onzoom(event)
        if event.key == 'up':
            self.yzoom += 10
            self.onzoom(event)
        if event.key == 'down':
            self.yzoom -= 10
            self.limit_zoom()
            self.onzoom(event)
        if event.key == 'right':
            self.xzoom += 10
            self.onzoom(event)
        if event.key == 'left':
            self.xzoom -= 10
            self.limit_zoom()
            self.onzoom(event)
        if event.key == 'z':
            self.yzoom -= 10
            self.xzoom -= 10
            self.limit_zoom()
            self.onzoom(event)
        if event.key == 'Z':
            self.yzoom += 10
            self.xzoom += 10
            self.onzoom(event)
        if event.key == 'delete':
            self.peaks = delete(self.peaks, self.pcounter, 0)
            self.pcounter = min(len(self.peaks)-1, self.pcounter)
            self.peakmarks.set_data((self.peaks.T[0], self.peaks.T[1]))
            self.onzoom(event)
        if event.key == 'S':
            with open("epynmr.peaks", "w") as fout:
                for peak in self.peaks:
                    fout.write("%10f %10f %20f\n" % tuple(peak))
        if event.key == 'L':
            peaks = []
            with open("epynmr.peaks") as fin:
                for line in fin:
                    peaks.append(map(float, line.split()))
            self.peaks = array(peaks)
            self.plot()
    def onzoom(self, event):
        x,y,h = self.peaks[self.pcounter]
        event.canvas.figure.get_axes()[0].set_xlim((x-self.xzoom-1,x+self.xzoom))
        event.canvas.figure.get_axes()[0].set_ylim((y-self.yzoom-1,y+self.yzoom))
        event.canvas.draw()

class TitrWindow(PeakWindow):
    def __init__(self, *args, **kwds):
        PeakWindow.__init__(self, *args, **kwds)
    def set_apo_hsqc(self, fname):
        self.apo_hsqc = fname
    def load_apo_hsqc(self):
        self.set_hsqc(hsqc(self.apo_hsqc))
    def set_holo_range(self, holonums):
        self.holonums = holonums
        self.curholo = 0
    def load_current_holohsqc(self):
        self.set_hsqc(hsqc(str(self.holonums[self.curholo])+'.nv'))
    def set_apo_peaks(self, peaks):
        self.apo_peaks = peaks
        self.apo_dataset = self.dataset
        self.acounter = 0
    def limit_acounter(self):
        self.acounter = max(0, self.acounter)
        self.acounter = min(len(self.apo_peaks)-1, self.acounter)
    def onkeypress(self, event):
        PeakWindow.onkeypress(self, event)
        if event.key == 'A':
            self.set_apo_peaks(self.peaks)
            print str(len(self.peaks)) + " apo-peaks stored"
        if event.key == 'H':
            self.load_current_holohsqc()
            self.peak_search()
            self.plot()
            self.apoplot()
        if event.key == 'N':
            self.curholo = min(len(self.holonums)-1, self.curholo+1)
            print "Now loading sample #" + str(self.holonums[self.curholo])
            self.load_current_holohsqc()
            self.peak_search()
            self.plot()
            self.apoplot()
            print " ".join(map(lambda x : "%5.2f" % x, peakset(self.peaks).matchperc(peakset(self.apo_peaks))))
        if event.key == 'P':
            self.curholo = max(0, self.curholo-1)
            print "Now loading sample #" + str(self.holonums[self.curholo])
            self.load_current_holohsqc()
            self.peak_search()
            self.plot()
            self.apoplot()
            print " ".join(map(lambda x : "%5.2f" % x, peakset(self.peaks).matchperc(peakset(self.apo_peaks))))
        if event.key == ']':
            self.acounter += 1
            self.limit_acounter()
            self.onapozoom(event)
        if event.key == '[':
            self.acounter -= 1
            self.limit_acounter()
            self.onapozoom(event)
        if event.key == 'I':
            apox = peakset(self.apo_peaks)
            x = apox.matchsignal(self.dataset,True)
            print "Sample #" + str(self.holonums[self.curholo]).ljust(3) + ": " + " ".join(map(lambda x : "%5.2f" % x, map(lambda p : scoreatpercentile(x, p), range(0, 100+10, 10))))
        if event.key == 'D':
            apox = peakset(self.apo_peaks)
            f = apox.matchsignal(self.dataset,True)
            with open(str(self.holonums[self.curholo])+".out", 'w') as fout:
                fout.write("RPH: " + " ".join(map(lambda x : "%6.3f" % x, f)) + '\n')
        if event.key == 'Q':
            apox = peakset(self.apo_peaks)
            with open("holall.dat", 'w') as fout:
                for sample in self.holonums:
                    holo = hsqc(str(sample)+".nv")
                    f = apox.matchsignal(holo,True)
                    fout.write("RPHI " + str(sample).ljust(3) + " :" + " ".join(map(lambda x : "%6.3f" % x, f)) + '\n')
                    fout.write("RPHD " + str(sample).ljust(3) + " :" + " ".join(map(lambda x : "%5.2f" % x, map(lambda p : scoreatpercentile(f, p), range(0, 100+10, 10)))) + '\n')
                    print "RPHD "+ str(sample).ljust(3) + " :" + " ".join(map(lambda x : "%5.2f" % x, map(lambda p : scoreatpercentile(f, p), range(0, 100+10, 10))))
    def onapozoom(self, event):
        x,y,h = self.apo_peaks[self.acounter]
        event.canvas.figure.get_axes()[0].set_xlim((x-self.xzoom-1,x+self.xzoom))
        event.canvas.figure.get_axes()[0].set_ylim((y-self.yzoom-1,y+self.yzoom))
        event.canvas.draw()
    def apoplot(self):
        self.apocontur = self.axes[0].contour(self.apo_dataset.nv_data, self.V, colors='red')
        self.apopeakmarks = self.axes[0].plot(self.apo_peaks.T[0],self.apo_peaks.T[1],'g+')[0]
    def recontur(self, event):
        for coll in self.contur.collections:
            gca().collections.remove(coll)
        self.contur = self.axes[0].contour(self.dataset.nv_data, self.V, colors='blue')
        if 'apo_dataset' in dir(self):
            for coll in self.apocontur.collections:
                gca().collections.remove(coll)
            self.apocontur = self.axes[0].contour(self.apo_dataset.nv_data, self.V, colors='red')
        event.canvas.draw()
        self.canvas.draw()

