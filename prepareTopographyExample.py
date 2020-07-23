import numpy
import matplotlib.pyplot as plt
import os
import topography
numpy.set_printoptions(linewidth=200)
import matplotlib.colors as colors
from netCDF4 import Dataset

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return numpy.ma.masked_array(numpy.interp(value, x, y))


# control whether the figures will be saved
saveFig = False

# control whether the figures will be generated
toPlot = True

# number of zenith angle that will be used in table preparations
NZEN = 19

# control information added to netCDF file
ADD_AUX = False

def prepareSW():
    '''
    netCDF file: minimal info required to include topography effect to SW RT
    '''
    myCol = 90
    # topography
    fname = 'horizon_00100_00641.nc'
    with Dataset(fname) as fid:
        zenData = fid.variables['zen'][:]
        horAngle = fid.variables['azm'][:]

    NHORIZON = len(horAngle)
    NZEN =  len(zenData)
    oxAvMask = numpy.zeros(((myCol, NHORIZON, NZEN)))
    id = 0
    for jj in range(9):
        for kk in range(10):
            fname = 'horizon_%05d_%05d.nc' % (100 + jj * 60, 641 - kk * 60)
            oxAvMask[id] = numpy.transpose(Dataset(fname).variables['oxAvMask'][:])
            id += 1

    with Dataset('rrtmgp-sw-topography-sample.nc', 'w', format='NETCDF4_CLASSIC') as oid:
        oid.createDimension('zen', NZEN)
        oid.createDimension('azm', NHORIZON)
        oid.createDimension('col', myCol)

        tVar = oid.createVariable('oxAvMask', 'f4', ('col', 'azm', 'zen'))
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

def prepareLW():
    '''
    netCDF file: minimal info required to include topography effect to LW RT
    '''
    myCol = 90
    # topography
    V = numpy.zeros(myCol)
    Vav = numpy.zeros(myCol)
    InvCosA = numpy.zeros(myCol)

    id = 0
    for jj in range(9):
        for kk in range(10):
            fname = 'horizon_%05d_%05d.nc' % (100 + jj * 60, 641 - kk * 60)
            with Dataset(fname) as tid:
                V[id] = tid.variables['V'][:]
                Vav[id] = tid.variables['Vav'][:]
                InvCosA[id] = tid.variables['InvCosA'][:]
                id += 1

    with Dataset('rrtmgp-lw-topography-sample.nc', 'w', format='NETCDF4_CLASSIC') as oid:
        oid.description = 'Sample dataset created to illustrate topography effect on LW radiation'
        oid.createDimension('col', myCol)

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
    # parameters for DEM dataset
    nLat = 10800
    nLon = 21600

    print("Reading DEM database")
    with Dataset("ETOPO1_Ice_c_gdal.grd") as fid:
        zz = numpy.reshape(fid.variables['z'][:], (nLat, nLon))

    print("Reading DEM database DONE")

    dLat = 180.0 / nLat
    dLon = 360.0 / nLon

    dLat = 180.0 / nLat
    dLon = 360.0 / nLon

    # lats and lons of the DEM pixels centers
    #  note that latitude changes from North to South
    mLat = 90.0 - numpy.arange(nLat) * dLat + dLat / 2.
    #  note that longitude changes from West to East
    mLon = -180.0 + numpy.arange(nLon) * dLon + dLon / 2.

    # lats and lons of the computational area
    lat1 = -18 - 4
    lat2 = -16 + 4

    lon1 = -68 - 4
    lon2 = -67 + 4

    jLon1 = int((lon1 + 180.0) / dLon)
    jLon2 = int((lon2 + 180.0) / dLon)

    jLat1 = int((90.0 - lat2) / dLat)
    jLat2 = int((90.0 - lat1) / dLat)


    # extra pixels to compute horizen angles
    Nblocks = 100
    print("Area boundaries to be extracted from DEM")

    print("Lat:", 90.0 - (jLat1 - 1) * dLat - dLat / 2., 90.0 - (jLat2) * dLat - dLat / 2.)
    print("Lon:", -180.0 + (jLon1) * dLon - dLon / 2., -180.0 + (jLon2 + 1) * dLon - dLon / 2.)


    hh = zz[jLat1 - 1 - Nblocks:jLat2 + Nblocks, jLon1 - Nblocks:jLon2 + 1 + Nblocks]

    #  ATTENTION
    # the sample datase use negative heights for ocean floor
    # so set them to 0 over ocean
    hh[hh < 0] = 0.

    top = topography.topographyBlock(hh=hh,
                                     lons=mLon[jLon1 - 1 - Nblocks:jLon2 + Nblocks],
                                     lats=mLat[jLat1 - 1 - Nblocks:jLat2 + Nblocks], dlon=dLon, dlat=dLat)


    # create computational area grid
    gB = topography.gridBlock(x1=lon1, x2=lon2, y1=lat1, y2=lat2, delX=1, delY=1)
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

    for ix, cx1 in enumerate(gB.gridX.gridX[:-1]):
        cx2 = gB.gridX.gridX[ix+1]
        N1, N2 = gB.gridX.NX[ix:ix + 2]
        for iy, cy1 in enumerate(gB.gridY.gridX[:-1]):
            cy2 = gB.gridY.gridX[iy+1]
            M2, M1 = gB.gridY.NX[iy:iy + 2]
            dx, dy = topography.haversineStep((cy1+cy2)/2, top.dlon, top.dlat)
            print("----------------------------------------------------------------")
            print("Processing %d, %d : "%(ix, iy))
            print("x1, x2=", cx1, cx2, top.lons[N1:N1+2], top.lons[N2:N2+2], N1, N2)
            print("y1, y2=", cy1, cy2, top.lats[M1:M1+2], top.lats[M2:M2+2], M1, M2)
            print("dx,dy = ", dx, dy)

            topography.preProcInit(dx, dy, top.dlon, top.dlat)

            if toPlot:
                plt.figure(1)
                yy1=plt.axhline(cy1, color='k')
                yy2=plt.axhline(cy2, color='m')
                xx1=plt.axvline(cx1, color='k')
                xx2=plt.axvline(cx2, color='m')

                if saveFig: plt.savefig('figure_spot.%d.%d.png' % (ix, iy))
                plt.draw()
                plt.pause(.5)
                yy1.remove()
                yy2.remove()
                xx1.remove()
                xx2.remove()


            ncFile = 'horizon_%05d_%05d.nc' % (N1, M1)
            if os.path.isfile(ncFile):
                continue

            weight, tanA, S, horz, = topography.Preprocess(dx, dy, N1, N2, M1, M2, top, topography.gridCell(x1=cx1, x2=cx2, y1=cy1, y2=cy2))

            oxAv = numpy.zeros((NZEN, topography.NHORIZON))
            oxAvMask = numpy.zeros((NZEN, topography.NHORIZON))
            solzen = numpy.deg2rad(numpy.arange(NZEN) * 5.)
            for ja, sunA in enumerate(topography.horAngle):
                curHor = horz[:, :, ja]

                for kk, sunZ in enumerate(solzen):
                    mask = numpy.int32(curHor > sunZ)
                    oxAv[kk, ja] = numpy.mean(weight * (1. + tanA * numpy.tan(sunZ) * numpy.cos(sunA - S)))
                    oxAvMask[kk, ja] = numpy.mean(mask * weight * (1. + tanA * numpy.tan(sunZ) * numpy.cos(sunA - S)))

            # sky view factor
            V = numpy.zeros_like(weight)
            for an, angle in enumerate(topography.horAngle):
                tanTF = -numpy.arctan(1. / tanA / numpy.cos(angle - S))
                tanTF[tanA == 0] = numpy.pi / 2.
                tanTF[tanTF < 0] += numpy.pi
                V += numpy.sin(horz[:, :, an] + numpy.pi / 2. - tanTF) ** 2

            V /= 16.

            print ("writing to netcdf")
            with Dataset(ncFile, 'w', format='NETCDF4_CLASSIC') as fid:
                nlat, nlon,nazm,  =  horz.shape
                # define axis size
                fid.createDimension('azm', nazm)

                # create azm axis
                tVar = fid.createVariable('azm', 'f4', ('azm',))
                tVar.long_name = 'horizon azimuth'
                tVar.units = 'rad'
                tVar[:] = topography.horAngle

                if ADD_AUX:
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
                    tVar = fid.createVariable('weight', 'f4', ('lat', 'lon'))
                    tVar.standard_name = 'weight'
                    tVar.long_name = 'weight'
                    tVar.units = 'none'
                    tVar[:] = weight

                    # create tanA
                    tVar = fid.createVariable('tanA', 'f4', ('lat', 'lon'))
                    tVar.standard_name = 'slope angle'
                    tVar.long_name = 'slope angle'
                    tVar.units = 'rad'
                    tVar[:] = numpy.rad2deg(tanA)

                    # create S
                    tVar = fid.createVariable('S', 'f4', ('lat', 'lon'))
                    tVar.standard_name = 'slope aspect'
                    tVar.long_name = 'slope aspect'
                    tVar.units = 'rad_north'
                    tVar[:] = numpy.rad2deg(S)

                    # create horz
                    tVar = fid.createVariable('horz', 'f4', ('lat', 'lon', 'azm'))
                    tVar.standard_name = 'horizon angle'
                    tVar.long_name = 'horizon angle'
                    tVar.units = 'rad'
                    tVar[:] = horz

                fid.createDimension('zen', NZEN)
                # create azm axis
                tVar = fid.createVariable('zen', 'f4', ('zen',))
                tVar.long_name = 'solar zenith'
                tVar.units = 'deg'
                tVar[:] = numpy.rad2deg(solzen)

                # create SW albedo correction factor
                tVar = fid.createVariable('oxAvMask', 'f4', ('zen', 'azm'))
                tVar.standard_name = 'SW surface albedo correction factor'
                tVar.long_name = 'SW surface albedo correction factor'
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
                tVar.standard_name = 'averaged inverse cosine of slope angle'
                tVar.long_name = 'averaged inverse cosine of slope angle'
                tVar.units = 'unit'
                tVar[:] = numpy.mean(numpy.sqrt(1. + tanA ** 2))

            plt.figure(2)
            plt.clf()
            plt.subplot(131)
            plt.pcolor(numpy.rad2deg(topography.horAngle),  numpy.rad2deg(solzen),oxAv , cmap='bwr', norm=MidpointNormalize(midpoint=1.0))
            kk=plt.colorbar()
            plt.ylabel('solar zenith', fontsize=14)
            plt.xlabel('solar azimuth', fontsize=14)
            plt.tight_layout()
            plt.title('$f_{cor}^{no mask}$', fontsize=14)

            plt.subplot(132)
            plt.pcolor(numpy.rad2deg(topography.horAngle),  numpy.rad2deg(solzen),oxAvMask)
            kk=plt.colorbar()
            plt.xlabel('solar azimuth', fontsize=14)
            plt.tight_layout()
            plt.title('$f_{cor}^{masked}$', fontsize=14)

            plt.subplot(133)
            plt.pcolor(numpy.rad2deg(topography.horAngle),  numpy.rad2deg(solzen),oxAv - oxAvMask, cmap='bwr', norm=MidpointNormalize(midpoint=0))
            kk=plt.colorbar()
            kk.set_label(r'$f_{cor}$')
            plt.xlabel('solar azimuth', fontsize=14)
            plt.tight_layout()
            plt.title('$f_{cor}^{no mask} - f_{cor}^{masked}$', fontsize=14)
            if saveFig: plt.savefig('/home/ipolonsk/_current/topography/figure_f_cor.%d.%d.png'%(ix,iy))
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
            if saveFig: plt.savefig('/home/ipolonsk/_current/topography/figure_V.png')
            plt.draw()
            plt.pause(.1)
            print("writing to netcdf DONE")

    prepareSW()
    prepareLW()