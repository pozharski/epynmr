from nmrio import hsqc, peakset
from matplotlib.figure import Figure
from matplotlib.pyplot import figure, gca, show, grid, close, annotate
from scipy import logspace, log10, delete, array
import sys

import Tkinter, tkFileDialog

def file_or_data(data):
    return data if data.__class__ is hsqc else hsqc(data)
    
def viewwindow(data, box=None, *args, **kwds):
    fig = figure(FigureClass=ViewWindow)
    fig.set_box(box)
    fig.set_hsqc(file_or_data(data))
    fig.plot()
    fig.canvas.mpl_connect('key_press_event', fig.onkeypress)
    show()

def dualwindow(data1, data2, peakfile, box=None, *args, **kwds):
    fig = figure(FigureClass=DualWindow)
    fig.set_box(box)
    fig.set_hsqc(file_or_data(data1), file_or_data(data2))
    fig.plot()
    if peakfile:
        fig.loadpeaks(peakfile)
    fig.canvas.mpl_connect('key_press_event', fig.onkeypress)
    return fig

def peakwindow(data, num=50, box=None, apfname=None, *args, **kwds):
    fig = figure(FigureClass=PeakWindow)
    fig.set_box(box)
    fig.set_hsqc(data if data.__class__ is hsqc else hsqc(data))
    fig.set_max_peaknum(num)
    fig.peak_search()
    fig.plot()
    if apfname is not None:
        fig.loadapeaks(apfname)
    fig.canvas.mpl_connect('key_press_event', fig.onkeypress)
    return fig

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

class TKWindow(object):
    def __init__(self, margs, xargs, *args, **kwds):
        self.tkroot = Tkinter.Tk()
        self.build(margs)
    def build(self, margs):
        pass
    def mainloop(self):
        show()
        self.tkroot.mainloop()
        close()

class TKPeaks(TKWindow):
    def build(self, margs):
        if margs.input_file is None:
            margs.input_file = tkFileDialog.askopenfilename(title = "Select NV file",filetypes = (("NMRview files","*.nv"),("all files","*.*")))
        self.gwindow = peakwindow(margs.input_file, margs.num_peaks, [margs.bleft, margs.bright, margs.bbottom, margs.btop], margs.auxpeakfile)
        self.label = Tkinter.Label(self.tkroot, text="HSQC peaks tool")
        self.label.grid()
        buttons =   [   ('pageup',    'Contour up', 1, 0),
                        ('pagedown',  'Contour down', 1, 1),
                        ('home','First peak', 2, 0),
                        ('+','Next peak', 2, 1),
                        ('-','Previous peak', 2, 2),
                        ('z','Zoom in', 3, 0),
                        ('Z','Zoom out', 3, 1),
                        ('delete','Delete current peak', 4, 0),
                    ]
        self.buttons = {}
        for key, label, row, column in buttons:
            self.add_button(key, label, row, column)
        Tkinter.Button(self.tkroot, text='Save peaks', command=self.savepeaks).grid(row=5, column=0)
        Tkinter.Button(self.tkroot, text='Load peaks', command=self.loadpeaks).grid(row=5, column=1)
        Tkinter.Button(self.tkroot, text='Load auxillary peaks', command=self.auxpeaks).grid(row=5, column=2)
        Tkinter.Button(self.tkroot, text='Quit', command=self.tkroot.quit).grid()
    def add_button(self, key, label, row, column):
        self.buttons[key] = Tkinter.Button(self.tkroot, text=label, command=lambda : self.gwindow.canvas.key_press_event(key))
        self.buttons[key].grid(row=row, column=column)
    def savepeaks(self):
        peakfile = tkFileDialog.asksaveasfilename(title = "Save as...",filetypes = (("Peak files","*.peaks"),("all files","*.*")))
        with open(peakfile, "w") as fout:
            for peak in self.gwindow.peaks:
                fout.write("%10f %10f %20f\n" % tuple(peak))
    def loadpeaks(self):
        peakfile = tkFileDialog.askopenfilename(title = "Load peaks",filetypes = (("Peak files","*.peaks"),("all files","*.*")))
        peaks = []
        with open(peakfile) as fin:
            for line in fin:
                peaks.append(map(float, line.split()))
        self.gwindow.peaks = array(peaks)
        self.gwindow.npks = self.gwindow.dataset.invconv(self.gwindow.peaks)
        self.gwindow.plot()
    def auxpeaks(self):
        self.gwindow.shiftsdat = tkFileDialog.askopenfilename(title = "Load auxillary peaks",filetypes = (("Peak files","*.dat"),("all files","*.*")))
        self.gwindow.canvas.key_press_event('A')

class NMRWindow(Figure):
    def __init__(self, *args, **kwds):
        Figure.__init__(self, *args, **kwds)
        self.zfloor = 1
    def log_contour_levels(self, z):
        return logspace(log10(z.mean()+self.zfloor*z.std()),log10(z.max()), 20)
    def onkeypress(self, event):
        if event.key == 'q':
            close(self)
            sys.exit()
    def set_box(self, box=None):
        self.box = box

class ViewWindow(NMRWindow):
    def set_hsqc(self, data):
        self.dataset = data
        self.set_contour_levels()
    def set_contour_levels(self):
        self.V = self.log_contour_levels(self.dataset.nv_data[self.box[3]:self.box[2],self.box[0]:self.box[1]])
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
        self.auxpeaks_loaded = False
    def set_hsqc(self, data1, data2):
        self.datasets = [ data1, data2 ]
        self.set_contour_levels()
    def set_contour_levels(self):
        self.V = map(lambda x : self.log_contour_levels(x.nv_data[self.box[3]:self.box[2],self.box[0]:self.box[1]]),  self.datasets)
    def plot(self):
        self.add_subplot(111)
        self.axes[0].clear()
        package = zip(*([self.datasets, self.V, ['red', 'blue']]))
        self.contur = map(lambda x : self.axes[0].contour(x[0].dim_gridvals(0), x[0].dim_gridvals(1), x[0].nv_data, x[1], colors=x[2]), package)
        grid(True)
        self.axes[0].set_xlim(sorted(self.axes[0].get_xlim(),reverse=True))
        self.axes[0].set_ylim(sorted(self.axes[0].get_ylim(),reverse=True))
        self.axlim = ['data', self.axes[0].get_xlim(), self.axes[0].get_ylim()]
    def recontur(self, event=None):
        for contur in self.contur:
            for coll in contur.collections:
                gca().collections.remove(coll)
        package = zip(*([self.datasets, self.V, ['red', 'blue']]))
        self.contur = map(lambda x : self.axes[0].contour(x[0].dim_gridvals(0), x[0].dim_gridvals(1), x[0].nv_data, x[1], colors=x[2]), package)
        if event:
            event.canvas.draw()
    def loadpeaks(self, fname):
        with open(fname) as fin:
            auxpeaks = zip(*map(lambda x : x.split(), fin))
        self.ax = array(auxpeaks[0]).astype(float)
        self.ay = array(auxpeaks[1]).astype(float)
        self.astepx = self.ax.ptp()*0.001
        self.astepy = self.ay.ptp()*0.001
        self.auxmarks = self.axes[0].plot(self.ax, self.ay, 'g+')[0]
        self.alabels = []
        for (i, name) in enumerate(auxpeaks[2]):
            self.alabels.append(annotate(name, (self.ax[i]-3*self.astepx, self.ay[i]-3*self.astepy)))
        self.marcounter=0
        self.labmark = None
        self.auxpeaks_loaded = True
    def limit_marcounter(self):
        self.marcounter = max(0, self.marcounter)
        self.marcounter = min(len(self.alabels)-1, self.marcounter)
    def onkeypress(self, event):
        NMRWindow.onkeypress(self, event)
        if event.key == 'pageup':
            self.zfloor += 1
        if event.key == 'pagedown':
            self.zfloor = max(1, self.zfloor-1)
        if event.key in ['pageup','pagedown']:
            self.set_contour_levels()
            self.recontur(event)
        if self.auxpeaks_loaded:
            if event.key == '+':
                self.marcounter += 1
            if event.key == '-':
                self.marcounter -= 1
            if event.key in '+-':
                self.limit_marcounter()
                if self.labmark:
                    self.labmark.remove()
                self.labmark = annotate(self.alabels[self.marcounter].get_text(), self.alabels[self.marcounter].xy)
                self.labmark.set_backgroundcolor('yellow')
                self.labmark.set_zorder(100)
                event.canvas.draw()
            if event.key == 'super':
                if self.labmark:
                    self.labmark.remove()
                    self.labmark = None
                else:
                    self.labmark = annotate(self.alabels[self.marcounter].get_text(), self.alabels[self.marcounter].xy)
                    self.labmark.set_backgroundcolor('yellow')
                    self.labmark.set_zorder(100)
                event.canvas.draw()
        if event.key == 'home':
            if self.axlim[0] == 'data':
                if self.auxpeaks_loaded:
                    self.axlim[0] = 'peak'
                    xpad, ypad = 0.05*self.ax.ptp(), 0.05*self.ay.ptp()
                    self.axes[0].set_xlim((max(self.ax)+xpad, min(self.ax)-xpad))
                    self.axes[0].set_ylim(max(self.ay)+ypad, min(self.ay)-ypad)
            elif self.axlim[0] == 'peak':
                self.axlim[0] = 'data'
                self.axes[0].set_xlim(self.axlim[1])
                self.axes[0].set_ylim(self.axlim[2])
            event.canvas.draw()
    def zfloor_down(self):
        self.zfloor = max(1, self.zfloor-1)
        self.set_contour_levels()
        self.recontur()
    def zfloor_up(self):
        self.zfloor += 1
        self.set_contour_levels()
        self.recontur()
    def zoom_peaks(self):
        if self.auxpeaks_loaded:
            self.axlim[0] = 'peak'
            xpad, ypad = 0.05*self.ax.ptp(), 0.05*self.ay.ptp()
            self.axes[0].set_xlim((max(self.ax)+xpad, min(self.ax)-xpad))
            self.axes[0].set_ylim(max(self.ay)+ypad, min(self.ay)-ypad)
    def aux_font_size(self, fs):
        for label in self.alabels:
            label.set_fontsize(fs)
        if self.labmark:
            self.labmark.set_fontsize(fs)

class PeakWindow(Figure):
    def __init__(self, *args, **kwds):
        Figure.__init__(self, *args, **kwds)
        self.xzoom = 10
        self.yzoom = 10
        self.zfloor = 1
        self.ax = []
        self.ay = []
        self.shiftsdat = 'shifts.dat'
        self.auxpeaks_loaded = False
    def set_max_peaknum(self, peaknum):
        self.peaknum = peaknum
    def set_box(self, box=None):
        self.box = box
    def set_hsqc(self, data):
        self.dataset = data
        self.set_vmatrix()
    def set_vmatrix(self):
        z = self.dataset.nv_data[self.box[3]:self.box[2],self.box[0]:self.box[1]]
        self.V = logspace(log10(z.mean()+self.zfloor*z.std()),log10(z.max()), 10)
    def peak_search(self):
        self.npks = self.dataset.peak_search(self.peaknum, self.box)
        self.peaks = self.dataset.xyzconv(self.npks)
        self.pcounter = 0
    def plot(self):
        x = self.dataset.dim_gridvals(0)
        y = self.dataset.dim_gridvals(1)
        z = self.dataset.nv_data
        self.add_subplot(111)
        self.axes[0].clear()
        self.contur = self.axes[0].contour(x, y, z, self.V, colors='blue')
        self.peakmarks = self.axes[0].plot(self.peaks.T[0],self.peaks.T[1],'r+')[0]
        self.auxmarks = self.axes[0].plot(self.ax, self.ay, 'go')[0]
        grid(True)
    def recontur(self, event):
        for coll in self.contur.collections:
            gca().collections.remove(coll)
        self.contur = self.axes[0].contour(self.dataset.dim_gridvals(0), self.dataset.dim_gridvals(1), self.dataset.nv_data, self.V)
        event.canvas.draw()
    def limit_pcounter(self):
        self.pcounter = max(0, self.pcounter)
        self.pcounter = min(len(self.peaks)-1, self.pcounter)
    def limit_auxcounter(self):
        self.auxcounter = max(0, self.auxcounter)
        self.auxcounter = min(len(self.ax)-1, self.auxcounter)
    def limit_zoom(self):
        self.xzoom = max(10, self.xzoom)
        self.yzoom = max(10, self.yzoom)
    def onkeypress(self, event):
        if event.key == 'pageup':
            self.zfloor += 1
            self.set_vmatrix()
            self.recontur(event)
        if event.key == 'pagedown':
            self.zfloor = max(1, self.zfloor-1)
            self.set_vmatrix()
            self.recontur(event)
        if event.key == 'home':
            self.pcounter = 0
            self.onzoom(event)
        if event.key == '+':
            self.pcounter += 1
        if event.key == '-':
            self.pcounter -= 1
        if event.key in '+-':
            self.limit_pcounter()
            self.onzoom(event)
        if event.key=='n':
            self.auxcounter += 1
        if event.key=='p':
            self.auxcounter -= 1
        if event.key in 'np':
            self.limit_auxcounter()
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
        if event.key == 'x':
            delta = array([self.dax,0])
            self.ax += self.dax
        if event.key == 'X':
            delta = array([-self.dax,0])
            self.ax -= self.dax
        if event.key == 'y':
            delta = array([0,self.day])
            self.ay += self.day
        if event.key == 'Y':
            delta = array([0,-self.day])
            self.ay -= self.day
        if event.key in 'xXyY':
            self.auxmarks.set_data((self.ax, self.ay))
            for alabel in self.alabels:
                alabel.set_position(delta+array(alabel.get_position()))
            self.onzoom(event)
        if event.key == 'delete':
            self.peaks = delete(self.peaks, self.pcounter, 0)
            self.npks = delete(self.npks, self.pcounter, 0)
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
        if event.key == 'A':
            self.loadapeaks(self.shiftsdat)
            self.onzoom(event)
        if event.key == 'M':
            if self.auxpeaks_loaded:
                self.automatch()
                self.onzoom(event)
        if event.key == 'q':
            close(self)
            sys.exit()
    def loadapeaks(self, fname):
        with open(fname) as fin:
            auxpeaks = zip(*map(lambda x : x.split()[:3], fin))
        self.ax = array(auxpeaks[0]).astype(float)
        self.ay = array(auxpeaks[1]).astype(float)
        self.dax = self.ax.ptp()*0.001
        self.day = self.ay.ptp()*0.001
        self.auxmarks.set_data((self.ax, self.ay))
        self.alabels = []
        for (i, name) in enumerate(auxpeaks[2]):
            self.alabels.append(annotate(name, (self.ax[i]-self.dax, self.ay[i]-self.day)))
        self.auxcounter = 0
        self.auxpeaks_loaded = True
    def automatch(self):
        axy = array([self.ax,self.ay]).T
        pxy = self.peaks[:,:2]
        sf=sqrt(axy.var(0)+pxy.var(0))
        nindex = [nonzero((abs((pxy-axy[i,:]))/sf<0.1).all(1))[0] for i in range(axy.shape[0])]
        # Continue here!!!
    def onzoom(self, event):
        if event.key in 'npAM':
            x,y = self.dataset.dim_convinv(0,self.ax[self.auxcounter]), self.dataset.dim_convinv(1,self.ay[self.auxcounter])
        else:
            x,y,h = self.npks[self.pcounter]
        xlims = map(lambda t : self.dataset.dim_conv(0,t), [x-self.xzoom-1,x+self.xzoom])
        ylims = map(lambda t : self.dataset.dim_conv(1,t), [y-self.yzoom-1,y+self.yzoom])
        event.canvas.figure.get_axes()[0].set_xlim(xlims)
        event.canvas.figure.get_axes()[0].set_ylim(ylims)
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

