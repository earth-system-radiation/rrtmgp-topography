import numpy
import matplotlib.pyplot as plt
import os, sys
import topography
numpy.set_printoptions(linewidth=200)
import matplotlib.colors as colors
from netCDF4 import Dataset

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to
        # make a simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return numpy.ma.masked_array(numpy.interp(value, x, y))

def readDEM(ncFile, numLat, numLon):
    '''
    Read in ETOPO1 Global Relief Model file and return its product
    (cell-registered, ice surface) as a NumPy array
    Call
        grid = readDEM(ncFile, numLat, numLon)

    Input
        ncFile -- string, ETOPO1 file from
            https://www.ngdc.noaa.gov/mgg/global/global.html
        numLat -- int, number of latitudes for output grid (x-axis)
        numLon -- int, number of longitudes for output grid (y-axis)

    Output
        grid -- NumPy array, surface heights
    '''

    print("Reading DEM database")
    with Dataset(ncFile) as fid:
        grid = numpy.reshape(fid.variables['z'][:], (numLat, numLon))

    print("Reading DEM database DONE")
    return grid

def topCalc(inGrid, numLat, numLon, saveFig=False,
    toPlot=False, nBlocks=100, add_aux=False, nZen=19, nHz=16):
    '''
    Topology parameter calculation and netcdf generation

    Call
        nxGridX, nxGridY = topCalc(inGrid, numLat, numLon, saveFig=False,
            toPlot=False, nBlocks=100, add_aux=False, nZen=19)

    Inputs
        inGrid -- NumPy array, grid from readDEM()
        numLat -- int, number of grid latitudes
        numLon -- int, number of grid longitudes

    Outputs
        nxGridX, nxGridY -- integer sets (list with unique elements),
            this is necessary for naming convention used and
            referenced throughout module; nxGridY is in descending order

    Keywords
        toPlot -- boolean, plot the topology (height, V) maps
        saveFig -- boolean, save height maps to PNG files
        nBlocks -- int,
        add_aux -- boolean, control information added to netCDF file
        nZen -- int, number of zenith angle that will be used in
            table preparations
        nHz -- int, number of horizons
    '''

    dLat = 180.0 / numLat
    dLon = 360.0 / numLon

    dLat = 180.0 / numLat
    dLon = 360.0 / numLon

    # lats and lons of the DEM pixels centers
    #  note that latitude changes from North to South
    mLat = 90.0 - numpy.arange(nLat) * dLat + dLat / 2.
    #  note that longitude changes from West to East
    mLon = -180.0 + numpy.arange(nLon) * dLon + dLon / 2.

    # lats and lons of the computational area
    # RP: Robert does not like magic numbers, these should be inputs
    # but i don't know what they are
    # (reference coordinates? what is 4?)
    lat1 = -18 - 4
    lat2 = -16 + 4

    lon1 = -68 - 4
    lon2 = -67 + 4

    # RP: what is "j"?
    jLon1 = int((lon1 + 180.0) / dLon)
    jLon2 = int((lon2 + 180.0) / dLon)

    jLat1 = int((90.0 - lat2) / dLat)
    jLat2 = int((90.0 - lat1) / dLat)

    # extra pixels to compute horizon angles
    print("Area boundaries to be extracted from DEM")
    print("Lat:", 90.0 - (jLat1 - 1) * dLat - dLat / 2., 90.0 -
        (jLat2) * dLat - dLat / 2.)
    print("Lon:", -180.0 + (jLon1) * dLon - dLon / 2., -180.0 +
        (jLon2 + 1) * dLon - dLon / 2.)

    hh = inGrid[jLat1 - 1 - nBlocks:jLat2 + nBlocks,
        jLon1 - nBlocks:jLon2 + 1 + nBlocks]

    #  ATTENTION
    # the sample datase use negative heights for ocean floor
    # so set them to 0 over ocean
    hh[hh < 0] = 0.

    top = topography.topographyBlock(
        hh=hh,
        lons=mLon[jLon1 - 1 - nBlocks:jLon2 + nBlocks],
        lats=mLat[jLat1 - 1 - nBlocks:jLat2 + nBlocks],
        dlon=dLon, dlat=dLat)

    # create computational area grid
    gB = topography.gridBlock(
        x1=lon1, x2=lon2, y1=lat1, y2=lat2, delX=1, delY=1)
    gB.createGrid()
    gB.findStartGrid(top)

    print("x=", gB.x1, gB.x2)
    print("y=", gB.y1, gB.y2)
    print("dLat, dLon=", top.dlat, top.dlon)

    cx = gB.x1
    cy = gB.y1
    print((gB.gridX.gridX - top.lons[0]) / top.dlon)
    print((top.lats[0] - gB.gridY.gridX ) / top.dlat)
    print('gridX', gB.gridX)
    print(gB.gridX.NX)
    print(top.lons[gB.gridX.NX])
    print(top.lons[gB.gridX.NX]+top.dlon)

    print('gridY', gB.gridY)
    print(gB.gridY.NX)
    print(top.lats[gB.gridY.NX])
    print(top.lats[gB.gridY.NX]-top.dlon)
    if toPlot:
        plt.figure(2, figsize=(12, 5))
        plt.figure(1)
        top.show(isStatic=False)

    nxGridY, nxGridX = [], []
    for ix, cx1 in enumerate(gB.gridX.gridX[:-1]):
        cx2 = gB.gridX.gridX[ix+1]
        N1, N2 = gB.gridX.NX[ix:ix + 2]
        nxGridX.append(N1)
        for iy, cy1 in enumerate(gB.gridY.gridX[:-1]):
            cy2 = gB.gridY.gridX[iy+1]
            M2, M1 = gB.gridY.NX[iy:iy + 2]
            nxGridY.append(M1)

            # RP: these variables have to be more descriptive
            ncFile = 'horizon_%05d_%05d.nc' % (N1, M1)
            if os.path.isfile(ncFile): continue
            print('Writing {}'.format(ncFile))

            dx, dy = topography.haversineStep(
                (cy1+cy2)/2, top.dlon, top.dlat)
            print(70*"-")
            print("Processing %d, %d : "%(ix, iy))
            print("x1, x2=", cx1, cx2,
                top.lons[N1:N1+2], top.lons[N2:N2+2], N1, N2)
            print("y1, y2=", cy1, cy2,
                top.lats[M1:M1+2], top.lats[M2:M2+2], M1, M2)
            print("dx,dy = ", dx, dy)

            topography.preProcInit(dx, dy, top.dlon, top.dlat)

            if toPlot:
                plt.figure(1)
                yy1=plt.axhline(cy1, color='k')
                yy2=plt.axhline(cy2, color='m')
                xx1=plt.axvline(cx1, color='k')
                xx2=plt.axvline(cx2, color='m')

                if saveFig:
                    plt.savefig('figure_spot.%d.%d.png' % (ix, iy))

                plt.draw()
                plt.pause(.5)
                yy1.remove()
                yy2.remove()
                xx1.remove()
                xx2.remove()

            weight, tanA, S, horz, = topography.Preprocess(
                dx, dy, N1, N2, M1, M2, top,
                topography.gridCell(x1=cx1, x2=cx2, y1=cy1, y2=cy2))

            oxAv = numpy.zeros((nZen, nHz))
            oxAvMask = numpy.zeros((nZen, nHz))
            solzen = numpy.deg2rad(numpy.arange(nZen) * 5.)
            for ja, sunA in enumerate(topography.horAngle):
                curHor = horz[:, :, ja]

                for kk, sunZ in enumerate(solzen):
                    mask = numpy.int32(curHor > sunZ)
                    oxAv[kk, ja] = numpy.mean(
                        weight * (1. + tanA * numpy.tan(sunZ) *
                        numpy.cos(sunA - S)))
                    oxAvMask[kk, ja] = numpy.mean(mask *
                        weight * (1. + tanA * numpy.tan(sunZ) *
                        numpy.cos(sunA - S)))

            # sky view factor
            V = numpy.zeros_like(weight)
            for an, angle in enumerate(topography.horAngle):
                tanTF = -numpy.arctan(1. / tanA / numpy.cos(angle - S))
                tanTF[tanA == 0] = numpy.pi / 2.
                tanTF[tanTF < 0] += numpy.pi
                V += numpy.sin(horz[:, :, an] +
                    numpy.pi / 2. - tanTF) ** 2

            V /= 16.

            print ("Writing topography netCDF")
            with Dataset(ncFile, 'w', format='NETCDF4_CLASSIC') as fid:
                nlat, nlon,nazm,  =  horz.shape
                # define axis size
                fid.createDimension('azm', nazm)

                # create azm axis
                tVar = fid.createVariable('azm', 'f4', ('azm',))
                tVar.long_name = 'horizon azimuth'
                tVar.units = 'rad'
                tVar[:] = topography.horAngle

                if add_aux:
                    fid.createDimension('lat', nlat)
                    fid.createDimension('lon', nlon)
                    # create latitude axis
                    tVar = fid.createVariable('lat', 'f8', ('lat'))
                    tVar.standard_name = 'latitude'
                    tVar.long_name = 'latitude'
                    tVar.units = 'degrees_north'
                    tVar[:] = top.lats[M1:M2 + 1] + top.dlat / 2.

                    # create longitude axis
                    tVar = fid.createVariable('lon', 'f8', ('lon'))
                    tVar.standard_name = 'longitude'
                    tVar.long_name = 'longitude'
                    tVar.units = 'degrees_east'
                    tVar[:] = top.lons[N1:N2 + 1] + top.dlon / 2.

                    # create weight
                    tVar = fid.createVariable(
                        'weight', 'f4', ('lat', 'lon'))
                    tVar.standard_name = 'weight'
                    tVar.long_name = 'weight'
                    tVar.units = 'none'
                    tVar[:] = weight

                    # create tanA
                    tVar = fid.createVariable(
                        'tanA', 'f4', ('lat', 'lon'))
                    tVar.standard_name = 'slope angle'
                    tVar.long_name = 'slope angle'
                    tVar.units = 'rad'
                    tVar[:] = numpy.rad2deg(tanA)

                    # create S
                    tVar = fid.createVariable(
                        'S', 'f4', ('lat', 'lon'))
                    tVar.standard_name = 'slope aspect'
                    tVar.long_name = 'slope aspect'
                    tVar.units = 'rad_north'
                    tVar[:] = numpy.rad2deg(S)

                    # create horz
                    tVar = fid.createVariable(
                        'horz', 'f4', ('lat', 'lon', 'azm'))
                    tVar.standard_name = 'horizon angle'
                    tVar.long_name = 'horizon angle'
                    tVar.units = 'rad'
                    tVar[:] = horz

                fid.createDimension('zen', nZen)
                # create azm axis
                tVar = fid.createVariable('zen', 'f4', ('zen',))
                tVar.long_name = 'solar zenith'
                tVar.units = 'deg'
                tVar[:] = numpy.rad2deg(solzen)

                # create SW albedo correction factor
                tVar = fid.createVariable(
                    'oxAvMask', 'f4', ('zen', 'azm'))
                tVar.standard_name = \
                    'SW surface albedo correction factor'
                tVar.long_name = \
                    'SW surface albedo correction factor'
                tVar.units = 'none'
                tVar[:] = oxAvMask

                # create LW sky view factor
                tVar = fid.createVariable('V', 'f4')
                tVar.standard_name = 'LW view factor'
                tVar.long_name = 'LW view factor'
                tVar.units = 'none'

                tVar[:] = numpy.mean(V)

                tVar = fid.createVariable('Vav', 'f4')
                tVar.standard_name = 'weighted LW view factor'
                tVar.long_name = 'weighted LW view factor'
                tVar.units = 'none'
                tVar[:] = numpy.mean(V * numpy.sqrt(1. + tanA ** 2))

                tVar = fid.createVariable('InvCosA', 'f4')
                tVar.standard_name = \
                    'averaged inverse cosine of slope angle'
                tVar.long_name = \
                    'averaged inverse cosine of slope angle'
                tVar.units = 'unit'
                tVar[:] = numpy.mean(numpy.sqrt(1. + tanA ** 2))

            if toPlot:
                plt.figure(2)
                plt.clf()
                plt.subplot(131)
                plt.pcolor(numpy.rad2deg(topography.horAngle),
                    numpy.rad2deg(solzen),oxAv , cmap='bwr',
                    norm=MidpointNormalize(midpoint=1.0))
                kk=plt.colorbar()
                plt.ylabel('solar zenith', fontsize=14)
                plt.xlabel('solar azimuth', fontsize=14)
                plt.tight_layout()
                plt.title('$f_{cor}^{no mask}$', fontsize=14)

                plt.subplot(132)
                plt.pcolor(numpy.rad2deg(topography.horAngle),
                    numpy.rad2deg(solzen),oxAvMask)
                kk=plt.colorbar()
                plt.xlabel('solar azimuth', fontsize=14)
                plt.tight_layout()
                plt.title('$f_{cor}^{masked}$', fontsize=14)

                plt.subplot(133)
                plt.pcolor(numpy.rad2deg(topography.horAngle),
                    numpy.rad2deg(solzen),oxAv - oxAvMask, cmap='bwr',
                    norm=MidpointNormalize(midpoint=0))
                kk=plt.colorbar()
                kk.set_label(r'$f_{cor}$')
                plt.xlabel('solar azimuth', fontsize=14)
                plt.tight_layout()
                plt.title('$f_{cor}^{no mask} - f_{cor}^{masked}$',
                    fontsize=14)
                if saveFig:
                    plt.savefig('./figure_f_cor.%d.%d.png'%(ix,iy))

                plt.draw()
                plt.pause(.1)

                plt.figure(2)
                plt.clf()
                plt.pcolor(top.lons[N1:N2 + 1], top.lats[M1:M2 + 1], V)
                kk = plt.colorbar()
                kk.set_label('V')
                plt.xlabel('Longitude', fontsize=14)
                plt.ylabel('Latitude', fontsize=14)
                plt.title('V=%.3f' % (numpy.mean(V)))
                if saveFig: plt.savefig('./figure_V.png')
                plt.draw()
                plt.pause(.1)

            print("writing to netcdf DONE")

    return sorted(set(nxGridX)), sorted(set(nxGridY))[::-1]

def prepareSW(nxGridX, nxGridY,
    outNC='rrtmgp-sw-topography-sample.nc', nCol=90,
    ncFormat='NETCDF4_CLASSIC'):
    '''
    Generate a netCDF file with minimal info required to include
    topography effect to shortwave radiative transfer

    Call
        prepareSW(outNC='rrtmgp-sw-topography-sample.nc', nCol=90, \
            ncFormat='NETCDF4_CLASSIC')

    Input
        nxGridX, nxGridY -- int, number of x points in two grids as
            determined in topCalc()

    Keywords
        outNC -- string, path of netCDF to which RT info is written
        nCol -- int, number of columns in RT calculation
        nCFormat -- string, format in which netCDF is written
    '''

    # topography
    fname = 'horizon_00100_00641.nc'
    with Dataset(fname) as fid:
        zenData = fid.variables['zen'][:]
        horAngle = fid.variables['azm'][:]

    nHorizon = len(horAngle)
    nZen =  len(zenData)
    oxAvMask = numpy.zeros(((nCol, nHorizon, nZen)))
    id = 0
    for nxx in nxGridX:
        for nxy in nxGridY:
            fname = 'horizon_{:05d}_{:05d}.nc'.format(nxx, nxy)
            oxAvMask[id] = numpy.transpose(
                Dataset(fname).variables['oxAvMask'][:])
            id += 1

    with Dataset(outNC, 'w', format=ncFormat) as oid:
        oid.description = 'Sample dataset created to illustrate ' + \
            'topography effect on SW radiation'
        oid.createDimension('zen', nZen)
        oid.createDimension('azm', nHorizon)
        oid.createDimension('col', nCol)

        tVar = oid.createVariable(
            'oxAvMask', 'f4', ('col', 'azm', 'zen'))
        tVar.standard_name = 'SW surface albedo correction factor'
        tVar.long_name = 'SW surface albedo correction factor'
        tVar.units = 'none'
        tVar[:] = oxAvMask

        tVar = oid.createVariable('solAzm', 'f4', ('azm',))
        tVar.long_name = 'solar azimuth'
        tVar.units = 'deg'
        tVar[:] = numpy.rad2deg(horAngle)

        tVar = oid.createVariable('solZen', 'f4', ('zen',))
        tVar.long_name = 'solar zenith'
        tVar.units = 'deg'
        tVar[:] = zenData

def prepareLW(nxGridX, nxGridY,
    outNC='rrtmgp-lw-topography-sample.nc', nCol=90,
    ncFormat='NETCDF4_CLASSIC'):
    '''
    Generate a netCDF file with minimal info required to include
    topography effect to longwave radiative transfer

    Call
        prepareSW(outNC='rrtmgp-lw-topography-sample.nc', nCol=90, \
            ncFormat='NETCDF4_CLASSIC')

    Input
        nxGridX, nxGridY -- int, number of x points in two grids as
            determined in topCalc()

    Keywords
        outNC -- string, path of netCDF to which RT info is written
        nCol -- int, number of columns in RT calculation
        nCFormat -- string, format in which netCDF is written
    '''

    # topography
    V = numpy.zeros(nCol)
    Vav = numpy.zeros(nCol)
    InvCosA = numpy.zeros(nCol)

    id = 0
    for nxx in nxGridX:
        for nxy in nxGridY:
            fname = 'horizon_{:05d}_{:05d}.nc'.format(nxx, nxy)
            with Dataset(fname) as tid:
                V[id] = tid.variables['V'][:]
                Vav[id] = tid.variables['Vav'][:]
                InvCosA[id] = tid.variables['InvCosA'][:]
                id += 1

    with Dataset(outNC, 'w', format=ncFormat) as oid:
        oid.description = 'Sample dataset created to illustrate ' + \
            'topography effect on LW radiation'
        oid.createDimension('col', nCol)

        # create LW sky view factor
        tVar = oid.createVariable('V', 'f4', ('col',))
        tVar.standard_name = 'LW view factor'
        tVar.long_name = 'LW view factor'
        tVar.units = 'none'
        tVar[:] = V

        tVar = oid.createVariable('Vav', 'f4', ('col',))
        tVar.standard_name = 'weighted LW view factor'
        tVar.long_name = 'weighted LW view factor'
        tVar.units = 'none'
        tVar[:] =Vav

        tVar = oid.createVariable('InvCosA', 'f4', ('col',))
        tVar.standard_name = 'averaged inverse cosine of slope angle'
        tVar.long_name = 'averaged inverse cosine of slope angle'
        tVar.units = 'unit'
        tVar[:] = InvCosA

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(\
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, \
        description='')
    parser.add_argument('--n_lat', '-n_lat', type=int, default=10800, \
        help='Number of latitude grid points for the ETOPO1 ' + \
        'grid reshape.')
    parser.add_argument('--n_lon', '-n_lon', type=int, default=21600, \
        help='Number of longitude grid points for the ETOPO1 ' + \
        'grid reshape.')
    parser.add_argument('--infile', '-i', type=str, \
        default='ETOPO1_Ice_c_gdal.grd', \
        help='ETOPO1 1 Arc-Minute Global Relief Model output netCDF')
    parser.add_argument('--savefig', '-s', action='store_true', \
        help='Save figures of V and altitude maps as PNG files.')
    parser.add_argument('--plot', '-p', action='store_true', \
        help='Generate V and altitude plots in addition to ' + \
        'the output netCDF files and either save (--savefig) or ' + \
        'display them.')
    parser.add_argument('--n_zenith', '-nz', default=19, type=int, \
        help='Number of zenith angles')
    parser.add_argument('--n_horizons', '-nh', default=16, type=int, \
        help='Number of horizons.')
    parser.add_argument('--n_blocks', '-nb', default=100, \
        help='Number of blocks in calculation into which the ' + \
        'topography grid is partitioned.')
    parser.add_argument('--auxiliary', '-a', action='store_true', \
        help='Add auxiliary information into netCDF.')
    args = parser.parse_args()

    nLat = args.n_lat
    nLon = args.n_lon
    inFile = args.infile

    demGrid = readDEM(inFile, nLat, nLon)
    nxx, nxy = topCalc(demGrid, nLat, nLon,
        saveFig=args.savefig, toPlot=args.plot,
        nBlocks=args.n_blocks, add_aux=args.auxiliary,
        nZen=args.n_zenith, nHz=args.n_horizons)
    prepareSW(nxx, nxy)
    prepareLW(nxx, nxy)
