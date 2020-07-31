import numpy
import os
import pickle as pkl
import geodesy
import matplotlib.pyplot as plt

# GLOBAL VARIABLES (treating this like a common block)
TLEN = dict()
IDXXX = dict()
IDXXY = dict()
IDXYX = dict()
IDXYY = dict()
DELTX = dict()
DELTY = dict()

def calcHorizontalAngle(nHorizon):
    """
    Given a number of horizons, calculate angle at which each horizon
    is defined
    """

    dAngle = 360.0 / nHorizon
    return numpy.deg2rad(numpy.arange(nHorizon)*dAngle)

def haversineStep(lat, dLon, dLat):
    """
    """

    Reff = geodesy.Reff(lat) # in [m]
    aLon = numpy.cos(numpy.deg2rad(lat))**2 * \
        numpy.sin(numpy.deg2rad(1./2.))**2
    dLat *= Reff*numpy.deg2rad(1.)
    dLon *= Reff*2. * \
        numpy.arctan2(numpy.sqrt(aLon), numpy.sqrt(1. - aLon))
    return dLat, dLon


class gridCell():
    """
    """

    def __init__(self, x1=None, x2=None, y1=None, y2=None):
        # bounding box
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        # length in km
        self.dy = None
        self.dx = None

    def getDxDy(self, lat, dLat=1., dLon=1.):
        self.dx, self.dy = haversineStep(lat / 2, dLon, dLat)

    def __str__(self):
        print("x1, x2 = ", self.x1, self.x2)
        print("y1, y2 = ", self.y1, self.y2)
        print("dx, dy = ", self.dx, self.dy)

# ======================================================================
#  define block of cells

class grid1D():
    """
    """

    def __init__(self, x1=None, x2=None, dX=None):
        print ("grid1D init", x1, x2, dX)
        self.x1 = x1
        self.x2 = x2
        self.dX = dX
        self.gridX = None
        self.NX = None
        self.idx = None

    def createGrid(self):
        """
        """

        self.gridX = numpy.linspace(
            self.x1, self.x2, (self.x2 - self.x1) // self.dX + 1)
        self.idx = 0

    def findStart(self, start, delta):
        """
        """

        if start < self.gridX[0]:
            self.NX = numpy.int32((self.gridX - start) / delta)
        else:
            self.NX = numpy.int32((start - self.gridX) / delta)

    def __iter__(self):
        return self

    def __next__(self):
        idx = self.idx
        self.idx += 1
        if self.idx >= len(self.gridX):
            raise StopIteration
        else:
            return idx, self.gridX[idx], self.gridX[idx+1], \
                self.NX[idx],  self.NX[idx+1]

    def __str__(self):
        """
        String that describes what is happening with object
        """

        out = '\n--grid1D-- \n x1, x2, dx:'
        out += ''.join( "%d "%(x) if x else ' None ' for x in \
            [self.x1, self.x2, self.dX])
        out += '\ngridX:      '
        out += 'None' if self.gridX is None else \
            ''.join('%d ' % (x) for x in self.gridX)

        out += '\nNX:         '
        out += 'None' if self.NX is None else ''.join(
            '%d ' % (x) for x in self.NX)
        return out

class gridBlock(gridCell):
    """
    """
    def __init__(self, x1=None, x2=None, y1=None, y2=None,
        delX=None, delY=None):

        print ("gridBlock init", delX, delY)
        super().__init__(x1=x1,  x2=x2, y1=y1, y2=y2)
        self.gridX = grid1D(x1, x2, delX)
        self.gridY = grid1D(y1, y2, delY)

    def readPickle(self, pklFile):
        with open(pklFile, 'rb') as fid:
            self.gridX.x1 = pkl.load(fid)
            self.gridX.x2 = pkl.load(fid)
            self.gridY.x1 = pkl.load(fid)
            self.gridY.x2 = pkl.load(fid)

    def writePickle(self, pklFile):
        with open(pklFile, 'wb') as fid:
            pkl.dump(self.gridX.x1, fid)
            pkl.dump(self.gridX.x2, fid)
            pkl.dump(self.gridY.x1, fid)
            pkl.dump(self.gridY.x2, fid)

    def createGrid(self):
        self.gridX.createGrid()
        self.gridY.createGrid()

    def findStartGrid(self, top):
        self.gridX.findStart(top.lons[0], top.dlon)
        self.gridY.findStart(top.lats[0], top.dlat)

class topographyBlock():
    """
    """

    def __init__(self, hh=None, lons=None, lats=None,
        dlon=None,  dlat=None):

        self.hh = hh
        self.lons = lons
        self.lats = lats
        self.dlon = dlon
        self.dlat = dlat

    def readPickle(self, pklFile):
        with open(pklFile, 'rb') as fid:
            self.hh = pkl.load(fid)
            self.dlat = pkl.load(fid)
            self.dlon = pkl.load(fid)
            self.lons = pkl.load(fid)
            self.lats = pkl.load(fid)

    def writePickle(self, pklFile):
        with open(pklFile, 'wb') as fid:
            pkl.dump(self.hh, fid)
            pkl.dump(self.dlat, fid)
            pkl.dump(self.dlon, fid)
            pkl.dump(self.lons, fid)
            pkl.dump(self.lats, fid)


    def show(self, gB=None, figFileName=None, isStatic=True):
        """
        """

        plt.pcolor(self.lons, self.lats, self.hh / 1e3)
        if gB:
            plt.axhline(gB.y1, color='k')
            plt.axhline(gB.y2, color='m')
            plt.axvline(gB.x1, color='k')
            plt.axvline(gB.x2, color='m')

        plt.xlabel('Longitude', fontsize=14)
        plt.ylabel('Latitude', fontsize=14)
        kk = plt.colorbar()
        kk.set_label('Altitude [km]')
        if isStatic:
            if figFileName:
                plt.savefig(figFileName)
            plt.show()

def preProcInit(dx, dy, dLon, dLat, maxDist=20000, nHz=16):
    """
    maxDist -- float, maximum distance in meters
    nHz -- int, number of horizons
    """

    horAngle = calcHorizontalAngle(nHz)

    for an, angle in enumerate(horAngle):
        maxY = - int(maxDist * numpy.sin(angle) / dy)
        maxX = + int(maxDist * numpy.cos(angle) / dx)
        ddY = numpy.sign(maxY)
        ddX = numpy.sign(maxX)
        TLEN[an] = numpy.array([])
        if ddX != 0:
            tx = 0.5 + numpy.arange(abs(maxX))* ddX
            ty = 0.5 + tx * numpy.tan(angle) * dLat / dLon*ddY
            TLEN[an] = numpy.append(
                TLEN[an], numpy.sqrt((tx * dx) ** 2 + (ty * dy) ** 2))
            iy = numpy.int32(numpy.floor(ty))
            DELTY[an] = ty - iy
            IDXXY[an] = iy
            IDXXX[an] = numpy.int32(numpy.floor(tx))

        if ddY != 0:
            ty = 0.5 + numpy.arange(abs(maxY)) * ddY
            tx = 0.5 + ty * numpy.cos(angle) / \
                numpy.sin(angle) * dLon / dLat
            TLEN[an] = numpy.append(
                TLEN[an], numpy.sqrt((tx * dx) ** 2 + (ty * dy) ** 2))
            ix = numpy.int32(numpy.floor(tx))
            DELTX[an] = tx - ix
            IDXYY[an] = numpy.int32(numpy.floor(ty))
            IDXYX[an] = ix

def Preprocess(dx, dy, N1, N2, M1, M2, top, cell, usepkl=False,
    nHz=16):
    """
    nHz -- int, number of horizons
    """

    pklFile = 'horizon_%05d_%05d.pkl'%(N1,M1)
    if usepkl and os.path.isfile(pklFile):
        with open(pklFile, 'rb') as fid:
            weight = pkl.load(fid)
            tanA = pkl.load(fid)
            S = pkl.load(fid)
            horz   = pkl.load(fid)
    else:
        weight = numpy.zeros((M2 - M1 + 1, N2 - N1 + 1))
        horz = numpy.zeros((M2 - M1 + 1, N2 - N1 + 1, nHz))
        grdY = ((top.hh[M1 + 1:M2 + 2, N1:N2 + 1] - \
                 top.hh[M1:M2 + 1, N1:N2 + 1]) + \
                (top.hh[M1 + 1:M2 + 2, N1 + 1:N2 + 2] - \
                 top.hh[M1:M2 + 1, N1 + 1:N2 + 2])) / 2. / dy
        grdX = ((top.hh[M1:M2 + 1, N1 + 1:N2 + 2] - \
                 top.hh[M1:M2 + 1, N1:N2 + 1]) + \
                (top.hh[M1 + 1:M2 + 2, N1 + 1:N2 + 2] - \
                 top.hh[M1 + 1:M2 + 2, N1:N2 + 1])) / 2. / dx
        cellAlt = ((top.hh[M1:M2 + 1, N1 + 1:N2 + 2] + \
                    top.hh[M1:M2 + 1, N1:N2 + 1]) + \
                   (top.hh[M1 + 1:M2 + 2, N1 + 1:N2 + 2] + \
                    top.hh[M1 + 1:M2 + 2, N1:N2 + 1])) / 4.
        for nk, k in enumerate(range(N1, N2+1)):
            xk=top.lons[k]
            if xk < cell.x1:
                wx = 1. - (cell.x1-xk)/top.dlon
            elif top.lons[k + 1] > cell.x2:
                wx = (cell.x2 - xk)/top.dlon
            else:
                wx = 1.0

            for nj, j in enumerate(range(M1, M2+1)):
                yj = top.lats[j]
                if yj > cell.y2:
                    wy = wx*(1. -  (yj - cell.y2) / top.dlat)
                elif top.lats[j+1] < cell.y1:
                    wy = wx*(yj -cell.y1) / top.dlat
                else:
                    wy = wx
                weight[nj, nk] = wy

                for an, angle in enumerate(horAngle):
                    tH = numpy.array([])
                    if an in DELTY:
                        ty = DELTY[an]
                        iy = IDXXY[an] + j
                        ix = IDXXX[an] + k + 1
                        tH = numpy.append(tH,
                             (top.hh[iy+1, ix]*ty + \
                              top.hh[iy, ix]*(1.-ty)) )
                    if an in DELTX:
                        tx = DELTX[an]
                        iy = IDXYY[an] + j + 1
                        ix = IDXYX[an]  + k
                        tH = numpy.append(tH,
                            (top.hh[iy, ix+1]*tx + \
                             top.hh[iy, ix]*(1.-tx)) )

                    horz[nj,nk,an] = numpy.pi/2. - \
                        max(0, numpy.arctan(numpy.max(
                            (tH - cellAlt[nj, nk])/TLEN[an])))
        tanA = numpy.sqrt(grdY ** 2 + grdX ** 2)
        S = numpy.arctan2(- grdY, grdX)
        if usepkl:
            with open(pklFile,'wb') as fid:
                pkl.dump(weight, fid)
                pkl.dump(tanA, fid)
                pkl.dump(S, fid)
                pkl.dump(horz, fid)

    return weight, tanA, S, horz

if __name__ == "__main__":
    print (os.path.basename(__file__), "cannot be run standalone")
