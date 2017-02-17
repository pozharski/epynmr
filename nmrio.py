import array as binarray
from scipy import array, concatenate, sqrt, arange
from scipy.stats import scoreatpercentile
from scipy.interpolate import RectBivariateSpline

def readarray(fin, astype='I', num=1):
    '''
        Reads an element from a binary stream.
    '''
    var = binarray.array(astype)
    var.fromfile(fin, num)
    return var

class nvdata(object):
    def __init__(self, fname):
        with open(fname,'rb') as fin:
            self.nv_magic = readarray(fin, 'i')[0]
            self.nv_spare1 = readarray(fin, 'i')[0]
            self.nv_spare2 = readarray(fin, 'i')[0]
            self.nv_fheadersz = readarray(fin, 'i')[0]
            self.nv_bheadersz = readarray(fin, 'i')[0]
            self.nv_blkelems = readarray(fin, 'i')[0]
            self.nv_ndim = readarray(fin, 'i')[0]
            self.nv_temperature = readarray(fin, 'f')[0]
            self.nv_sequence = readarray(fin, 'c', 32).tostring().strip('\x00')
            self.nv_comment = readarray(fin, 'c', 160).tostring().strip('\x00')
            self.nv_month = readarray(fin, 'i')[0]
            self.nv_day = readarray(fin, 'i')[0]
            self.nv_year = readarray(fin, 'i')[0]
            self.nv_spare = readarray(fin, 'i',197).tolist()
            self.nv_dimheaders = []
            for i in range(8):
                self.nv_dimheaders.append(nvdim(fin))
            curpos = fin.tell()
            sizeofloat = readarray(fin, 'f').itemsize
            fin.seek(0,2)
            endpos = fin.tell()
            fin.seek(curpos,0)
            ndata = (endpos-curpos)/sizeofloat
            rawdata = array(readarray(fin, 'f',ndata).tolist())
            ms=map(lambda i : rawdata[self.nv_blkelems*i:self.nv_blkelems*(i+1)].reshape((self.nv_dimheaders[1].nv_blksize,self.nv_dimheaders[0].nv_blksize)), range(self.nv_dimheaders[1].nv_nblks*self.nv_dimheaders[0].nv_nblks))
            self.nv_data = concatenate(map(lambda i : concatenate(ms[i*self.nv_dimheaders[0].nv_nblks:(i+1)*self.nv_dimheaders[0].nv_nblks],1), range(self.nv_dimheaders[1].nv_nblks)))
    def __str__(self):
        return '\n'.join(["Magic number      : " + str(self.nv_magic),
                          "File header size  : " + str(self.nv_fheadersz),
                          "Block header size : " + str(self.nv_bheadersz),
                          "Block size        : " + str(self.nv_blkelems),
                          "Dimensions        : " + str(self.nv_ndim),
                          "Temperature       : " + str(self.nv_temperature),
                          "Sequence          : " + str(self.nv_sequence),
                          "Comment           : " + str(self.nv_comment),
                          "Month             : " + str(self.nv_month),
                          "Day               : " + str(self.nv_day),
                          "Year              : " + str(self.nv_year)
                         ])
    def dim_conv(self, dim, values):
        return self.nv_dimheaders[dim].dim_conv(values)
    def dim_gridvals(self, dim):
        return self.nv_dimheaders[dim].gridvals

class nvdim(object):
    def __init__(self, fin):
        self.nv_size, \
        self.nv_blksize, \
        self.nv_nblks, \
        self.nv_offblk, \
        self.nv_blkmask, \
        self.nv_offpt = readarray(fin, 'i', 6).tolist()
        self.nv_sf, \
        self.nv_sw, \
        self.nv_refpt, \
        self.nv_ref = readarray(fin, 'f', 4).tolist()
        self.nv_refunits = readarray(fin, 'i')[0]
        self.nv_foldup, \
        self.nv_folddown = readarray(fin, 'f', 2).tolist()
        self.nv_label = readarray(fin, 'c', 16).tostring().strip('\x00')
        self.nv_spares = readarray(fin, 'i',15).tolist()
        if self.nv_size:
            self.convfactor = self.nv_sw/self.nv_sf/self.nv_size
            self.gridvals = self.dim_conv(arange(self.nv_size))
    def dim_conv(self, values):
        return (self.nv_refpt-values)*self.convfactor+self.nv_ref
    def __str__(self):
        return '\n'.join(["Label           : " + self.nv_label,
                          "Size            : " + str(self.nv_size),
                          "Block size      : " + str(self.nv_blksize),
                          "Blocks          : " + str(self.nv_nblks),
                          "Offblk          : " + str(self.nv_offblk),
                          "Block mask      : " + str(self.nv_blkmask),
                          "Offpt           : " + str(self.nv_offpt),
                          "Frequency       : " + str(self.nv_sf),
                          "Spectral width  : " + str(self.nv_sw),
                          "Reference point : " + str(self.nv_refpt),
                          "Reference value : " + str(self.nv_ref),
                          "Reference units : " + str(self.nv_refunits),
                          "Fold up         : " + str(self.nv_foldup),
                          "Fold down       : " + str(self.nv_folddown)
                         ])
    def label(self):
        return self.nv_label

class hsqc(nvdata):
    def __init__(self, fname, *args, **kwds):
        nvdata.__init__(self, fname, *args, **kwds)
    def peakmatrix(self):
        z = self.nv_data
        zoo=z[1:-1,1:-1]
        zmp=z[0:-2,2:]
        zop=z[1:-1,2:]
        zpp=z[2:,2:]
        zmo=z[0:-2,1:-1]
        zpo=z[2:,1:-1]
        zmm=z[0:-2,0:-2]
        zom=z[1:-1,0:-2]
        zpm=z[2:,0:-2]
        a0  =  -8*zmp + 16*zop -  8*zpp + 16*zmo + 40*zoo + 16*zpo -  8*zmm + 16*zom -  8*zpm
        a1  = -12*zmp          + 12*zpp - 12*zmo          + 12*zpo - 12*zmm          + 12*zpm
        a2  =  12*zmp + 12*zop + 12*zpp                            - 12*zmm - 12*zom - 12*zpm
        a11 =  12*zmp - 24*zop + 12*zpp + 12*zmo - 24*zoo + 12*zpo + 12*zmm - 24*zom + 12*zpm
        a22 =  12*zmp + 12*zop + 12*zpp - 24*zmo - 24*zoo - 24*zpo + 12*zmm + 12*zom + 12*zpm
        a12 =  -9*zmp          +  9*zpp                            +  9*zmm          -  9*zpm
        D  = a11*a22-a12**2
        Dx = (a2*a12-a1*a22)/2
        Dy = (a1*a12-a2*a11)/2
        xe=Dx/D
        ye=Dy/D
        h = (a0 + a1*xe + a2*ye + a11*xe**2 + a22*ye**2 + 2*a12*xe*ye)*((xe**2+ye**2)<1).astype(int)*(D>0).astype(int)*(a11<0).astype(int)/72
        return h, xe, ye
    def xyzconv(self, xyz):
        x = self.dim_conv(0, xyz.T[0])
        y = self.dim_conv(1, xyz.T[1])
        z = xyz.T[2]
        return array([x,y,z]).T
    def peak_search(self, num, box=None, verbose=False):
        if verbose:
            sys.stdout.write("Looking for %d peaks\n" % num)
        h, xe, ye = self.peakmatrix()
        if box:
            if(box[0]):
                h[:,:box[0]] = 0
            if(box[1]):
                h[:,box[1]:] = 0
            if(box[2]):
                h[:box[2],:] = 0
            if(box[3]):
                h[box[3]:,:] = 0
        if verbose:
            sys.stdout.write("Got height matrix\n")
        peaks = []
        ymax, xmax = h.shape
        ymax -= 1
        xmax -= 1
        while len(peaks)< num:
            if h.max() > 0:
                x,y = h.max(0).argmax(), h.max(1).argmax()
                if verbose:
                    sys.stdout.write("Peak located at pixel %d,%d\n" % (x+1,y+1))
                    sys.stdout.write("Local shift %.3f %.3f\n" % (xe[y,x], ye[y,x]))
                    sys.stdout.write("Height: %g\n" % h[y][x])
                    sys.stdout.write(str(h.max())+'\n')
                peaks.append([x+ye[y,x]+1,y+xe[y,x]+1, h[y][x]])
                h[max(y-2,0):min(y+3,ymax),max(x-2,0):min(x+3,xmax)] = 0
            else:
                break
        if verbose:
            sys.stdout.write("Found %d peaks\n" % len(peaks))
        return array(peaks)

class peakset(object):
    def __init__(self, peaks, *args, **kwds):
        self.peaks = array(peaks)
    
    def distance2(self, peak):
        return ((self.peaks[:,:2].astype(float)-peak[:2].astype(float))**2).sum(1)

    def singlepeakmatch(self, peak):
        d2 = self.distance2(peak)
        ind = d2.argsort()
        return ind[0], sqrt(d2[ind[0]]), sqrt(d2[ind[1]]/d2[ind[0]])

    def peakmatch(self, other):
        return zip(*map(self.singlepeakmatch, other.peaks))

    def matchcount(self, other, d_cutoff=1.0, s_cutoff=2.0):
        pms = self.peakmatch(other)
        return sum((array(pms[2])>s_cutoff)*(array(pms[1])<d_cutoff))

    def matchperc(self, other, step=10):
        x = self.peakmatch(other)[1]
        return map(lambda p : scoreatpercentile(x, p), range(0, 100+step, step))

    def matchsignal(self, spec, normalize=False):
        a, b = spec.nv_data.shape
        f = RectBivariateSpline(range(a),range(b),spec.nv_data)
        if normalize:
            return f.ev(self.peaks.T[1], self.peaks.T[0])/self.peaks.T[2]
        else:
            return f.ev(self.peaks.T[1], self.peaks.T[0])
