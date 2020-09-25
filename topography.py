import numpy
import os
import geodesy
import matplotlib.pyplot as plt

dbg = False
# GLOBAL VARIABLES (treating this like a common block)
TLENX = dict()
TLENY = dict()

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
    return azimuth angles in [rad]
    """

    dAngle = 360.0 / nHorizon
    return numpy.deg2rad(numpy.arange(nHorizon)*dAngle) # [rad]

def haversineStep(lat, dLon, dLat):
    """
        Given latitide returns dx, and dy corresponding to
        increment of lon and lat by dLon, dLat
        dLon and dLat in deg
        dx and dy in radians
    """
    Reff = geodesy.Reff(lat) # in [m]
    aLon = numpy.cos(numpy.deg2rad(lat))**2 * \
        numpy.sin(numpy.deg2rad(dLon/2.))**2

    dy = Reff*numpy.deg2rad(dLat)
    dx = Reff*2. * \
        numpy.arctan2(numpy.sqrt(aLon), numpy.sqrt(1. - aLon))
    return dx, dy

class gridCell():
    """
    class represent functionality for a computational Grid cell
    """
    def __init__(self, x1=None, x2=None, y1=None, y2=None):
        """
        initialization: bounding box
        """
        # bounding box in [degree]
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2

    def __str__(self):
        out  = "x1, x2 = %f, %f\n"%(self.x1, self.x2)
        out += "y1, y2 = %f, %f\n"%(self.y1, self.y2)
        return out

# ======================================================================
#  define block of cells

class grid1D():
    """
    class that represents 1D regular grid
    """
    def __init__(self, x1=None, x2=None, dX=None):
        if dbg: print ("grid1D init", x1, x2, dX)
        self.x1 = x1        # grid start
        self.x2 = x2        # grid end
        self.dX = dX        # grid step
        self.grid = None   # numpy array that contains all grid points
        self.idx = None     # selected grid index

    def createGrid(self):
        """
        creates uniform grid and points to the grid head
        """
        self.grid = numpy.linspace(
            self.x1, self.x2, num=int((self.x2 - self.x1) // self.dX + 1))
        self.idx = 0

    def floor(self, idx, start, delta):
        """
        for a regular sequence  z[k] = start + delta[k]
        return index that found index k0 that    z[k0] <= self.grid[idx] < z[k0+1]
        """
        if start < self.grid[0]:
            return numpy.int32(numpy.floor((self.grid[idx] - start) / delta))
        else:
            return numpy.int32(numpy.floor((start - self.grid[idx]) / delta))

    def ceil(self, idx, start, delta):
        """
        for a regular sequence  z[k] = start + delta[k]
        return index that found index k0 that    z[k0-1] <= self.grid[idx] < z[k0]
        """
        ind = self.floor(idx, start, delta) + 1
        if abs(start + delta*ind - self.grid[idx])< 1e-3 *delta:
            return ind
        else:
            return ind - 1

    def __iter__(self):
        """
        default iterator
        """
        return self

    def __next__(self):
        idx = self.idx
        self.idx += 1
        if self.idx >= len(self.grid):
            self.idx = 0
            raise StopIteration
        else:
            return idx, self.grid[idx], self.grid[idx+1]

    def __str__(self):
        """
        String that describes what is happening with object
        """
        out = '\n--grid1D-- \n x1, x2, dx:'
        out += ''.join( "%d "%(x) if x else ' None ' for x in \
            [self.x1, self.x2, self.dX])
        out += '\ngrid:      '
        out += 'None' if self.grid is None else \
            ''.join('%d ' % (x) for x in self.grid)

        return out


class gridBlock(gridCell):
    """
    class to represent computational 2D regular grid
    """
    def __init__(self, x1=None, x2=None, y1=None, y2=None,
        delX=None, delY=None):

        if dbg: print ("gridBlock init", delX, delY)
        super().__init__(x1=x1,  x2=x2, y1=y1, y2=y2)
        self.gridX = grid1D(x1, x2, delX)
        self.gridY = grid1D(y1, y2, delY)
        self.createGrid()

    def createGrid(self):
        if dbg: print ('creating gridX')
        self.gridX.createGrid()
        if dbg: print ('creating gridY')
        self.gridY.createGrid()


class topographyBlock():
    """
    define DEM subset with lats and lons representing its grid
    with resolution dlat, dlon
    """
    def __init__(self, dem=None, lons=None, lats=None,
        dlon=None,  dlat=None):

        self.dem = dem
        self.lons = lons
        self.lats = lats
        self.dlon = dlon
        self.dlat = dlat

    def show(self, gB=None, figFileName=None, isStatic=True):
        """
        display topography block on figure
        """
        plt.pcolor(self.lons, self.lats, self.dem / 1e3)

        # show a gridBlock is provided
        if gB:
            plt.plot([gB.x1, gB.x2, gB.x2, gB.x1, gB.x1],
                     [gB.y1, gB.y1, gB.y2, gB.y2, gB.y1], color='m')

        plt.xlabel('Longitude', fontsize=14)
        plt.ylabel('Latitude', fontsize=14)
        kk = plt.colorbar()
        kk.set_label('Altitude [km]')
        if isStatic:
            if figFileName:
                plt.savefig(figFileName)
            plt.show()

    def check(self, N1, N2, M1, M2, delta=1.):
        """
        check if there is a variability in  top[M1:M2, N1:N2]
        maxDist -- float, maximum distance used to compute horizon angles [m]
        return True if variability detected otherwise False
        """
        if dbg: print('evaluation of DEM ', N1, N2, M1, M2)
        maxAlt = numpy.max(self.dem[M1:M2, N1:N2])
        minAlt = numpy.min(self.dem[M1:M2, N1:N2])
        print ('   check altitude: max, min', maxAlt, minAlt)
        return  maxAlt-minAlt  > delta


def preProcInit(dx, dy, dLon, dLat, horAngle, maxDist=20000, debug=True):
    """
    maxDist -- float, maximum distance used to compute horizon angles [m]
    horAngle -- azimuthal angle at which the horizon angle are computed [rad]
                azimuth is computed from=n direction to the North clockwise

    (dx, dy)      - computational grid resolution in (lon, lat)
    (dLon, dLat)  - DEM grid resolution in (lon, lat)

    To compute the horizon angle for a give cell into a given azimuthal direction
    we find the intersections of view line with all vertical and horizontal lines
    connecting DEM cell centers distant less than maxDist to the given DEM cell center

    If computational area is not large this routine can be called just once.
    """
    for an, angle in enumerate(horAngle):
        # maxY, maxX is max number of computational cells along lon and lat
        if debug: print('\nprocessing: ', an, numpy.rad2deg(angle))
        maxX = + int(maxDist * numpy.sin(angle) / dx)
        ddX = numpy.sign(maxX)

        #  minus due to reverse DEM direction along latitudes
        maxY = - int(maxDist * numpy.cos(angle) / dy)
        ddY = numpy.sign(maxY)

        # if the view line has increment along longitudes
        if ddX != 0:
            tx = (1.0 + numpy.arange(abs(maxX)))* ddX
            #  corresponding increment along latitudinal direction
            ty = numpy.abs(tx * numpy.cos(angle) / numpy.sin(angle)) * dLat / dLon*ddY
            # represent distance between DEM center and intersections with
            # lines along latitudes
            TLENX[an] = numpy.sqrt((tx * dx) ** 2 + (ty * dy) ** 2)

            # index increment relative DEM cell latitude index
            iy = numpy.int32(numpy.floor(ty))
            IDXXY[an] = iy

            # the distance along latitudinal direction to the intersection point from
            # the nearest DEM cell (used in interpolation)
            DELTY[an] = ty - iy

            # index increment relative DEM cell longitude index
            IDXXX[an] = numpy.int32(numpy.floor(tx))
            if debug: print('IDXXX[an]: ', IDXXX[an])
            if debug: print('IDXXY[an]: ', IDXXY[an])
            if debug: print('TLENX[an]: ', TLENX[an])
            if debug: print('DELTY[an]: ', DELTY[an])
            if debug: print('')

        # if the view line has increment along latitude
        if ddY != 0:
            ty = (1.0 + numpy.arange(abs(maxY))) * ddY
            #  corresponding increment along longitudinal direction
            tx = numpy.abs(ty * numpy.sin(angle) / numpy.cos(angle)) *dLon / dLat * ddX

            # represent distance between DEM center and intersections with
            # lines along longitudes
            TLENY[an] =  numpy.sqrt((tx * dx) ** 2 + (ty * dy) ** 2)

            ix = numpy.int32(numpy.floor(tx))
            # index increment relative DEM cell longitude index
            IDXYX[an] = ix

            # the distance along longitudinal direction to the intersection point from
            # the nearest DEM cell (used in interpolation)
            DELTX[an] = tx - ix

            # index increment relative DEM cell latitude index
            IDXYY[an] = numpy.int32(numpy.floor(ty))

            if debug: print('IDXYY[an]: ', IDXYY[an])
            if debug: print('IDXYX[an]: ', IDXYX[an])
            if debug: print('TLENY[an]: ', TLENY[an])
            if debug: print('DELTX[an]: ', DELTX[an])
            if debug: print('')


def Preprocess(dx, dy, N1, N2, M1, M2, top, cell, horAzmAngle, debug=True):
    """
    dx, dy - resolution in [m] of DEM cell
    N1:N2 - DEM cell indices along longitude
    M1:M2 - DEM cell indices along latitude
    indices N2 and M2 indicates the upper boundary in pythonic sense

    horAzmAngle -- azimuthal angle at which the horizon angle computed [rad]
    output: for all DEM cell that are in computational cell
        weight - part of the area that DEM cell contributes to compuational cell
        tanSlopeAngle - tan of the DEM cell slope angle (with respect to slope normal)
        slopeAspect - DEM cell aspect angle (with respect to North direction [rad]
        horAngle    - horizon angles for horAzmAngle directions  [rad]
                    defined with respect to normal (flat corresponds to pi/2
    """

    # gradients along latitudinal direction
    grdY = ((top.dem[M1 + 1:M2 + 1, N1:N2] - top.dem[M1:M2, N1:N2]) + \
            (top.dem[M1 + 1:M2 + 1, N1 + 1:N2 + 1] - top.dem[M1:M2, N1 + 1:N2 + 1])) / 2. / dy

    # gradients along longitudinal direction
    grdX = ((top.dem[M1:M2, N1 + 1:N2 + 1] - top.dem[M1:M2, N1:N2]) + \
            (top.dem[M1 + 1:M2 + 1, N1 + 1:N2 + 1] -top.dem[M1 + 1:M2 + 1, N1:N2])) / 2. / dx

    # tan of slope angle
    #  more effective to store tan than angle
    tanSlopeAngle = numpy.sqrt(grdY ** 2 + grdX ** 2)

    # slope aspect
    #  minus due to DEM grid directions
    #  aspect angle defined with respect to the North direction
    slopeAspect = numpy.arctan2(grdX, -grdY) #  in [rad]

    if debug:
        print(cell)
        print(top.lons[N1], top.lons[N2])
        print(top.lats[M1], top.lats[M2])

    #  compute weight
    # initialization
    weight = numpy.ones((M2 - M1, N2 - N1))
    #  check the boundaries
    xW = top.lons[N1]
    if xW < cell.x1:
        weight[:,0] *= 1. + (xW - cell.x1) / top.dlon

    xE = top.lons[N2+1]
    if xE > cell.x2:
       weight[:,-1] *= 1. + (cell.x2 - xE) / top.dlon

    yN = top.lats[M1]
    if yN > cell.y2:
        weight[0,:] *= 1. + (cell.y2 - yN) / top.dlat

    yS = top.lats[M2 + 1]
    if yS < cell.y1:
        weight[-1,:] *= 1. + (yS - cell.y1) / top.dlat

    #  compute horizon angles
    horAngle = numpy.zeros((M2 - M1, N2 - N1, len(horAzmAngle)))
    for nk, k in enumerate(range(N1, N2)):
        for nj, j in enumerate(range(M1, M2)):
            for an, angle in enumerate(horAzmAngle):
                maxAngle = 0.
                # if the view line has increment along longitudes
                if an in IDXXX:
                    iy = IDXXY[an] + j
                    ix = IDXXX[an] + k
                    msk = numpy.bitwise_and(iy >= 0, ix >= 0)
                    tH = numpy.array((top.dem[iy[msk]+1, ix[msk]]*DELTY[an][msk] + top.dem[iy[msk], ix[msk]]*(1.-DELTY[an][msk])) - top.dem[j, k])
                    maxAngle = max(maxAngle, numpy.arctan(numpy.max(tH/TLENX[an])))

                # if the view line has increment along latitude
                if an in IDXYY:
                    iy = IDXYY[an] + j
                    ix = IDXYX[an] + k
                    msk = numpy.bitwise_and(iy >= 0, ix >= 0)
                    tH = numpy.array((top.dem[iy[msk]+1, ix[msk]]*DELTX[an][msk] + top.dem[iy[msk], ix[msk]]*(1.-DELTX[an][msk])) - top.dem[j, k])
                    maxAngle = max(maxAngle, numpy.arctan(numpy.max(tH/TLENY[an])))

                # in [rad]
                horAngle[nj,nk,an] = numpy.pi/2. - maxAngle

    return weight, tanSlopeAngle, slopeAspect, horAngle

if __name__ == "__main__":
    print (os.path.basename(__file__), "cannot be run standalone")
