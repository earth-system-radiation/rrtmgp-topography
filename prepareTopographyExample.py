import numpy
import matplotlib.pyplot as plt
import os, sys
import topography
numpy.set_printoptions(linewidth=200)
import matplotlib.colors as colors
from netCDF4 import Dataset

import warnings
warnings.filterwarnings("ignore")
dbg = False
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

        numLat and numLon must be consistent with data provided with ncFile

    Output
        dem[numLat, numLon] -- NumPy array, surface heights
    '''

    print("Reading DEM database")
    with Dataset(ncFile) as fid:
        dem = numpy.reshape(fid.variables['z'][:], (numLat, numLon))

    print("Reading DEM database DONE")
    return dem


def wholeDEM(dem, nBlocks=100):
    '''
    return the whole DEM extended by nBlocks in all directions

    Input
        dem   -- NumPy array, surface heights

    Output
        topography.topographyBlock object
    '''
    numLat, numLon, = dem.shape
    dLat = 180.0 / numLat
    dLon = 360.0 / numLon

    # lats and lons of the DEM pixels centers
    #  note that latitude changes from North to South
    mLat = 90.0 - numpy.arange(numLat) * dLat + dLat / 2.
    #  note that longitude changes from West to East
    mLon = -180.0 + numpy.arange(numLon) * dLon + dLon / 2.

    # since DEM cells boundaries goes from 90, 90 - dLat, 90 -2*dLat, ...
    # latitude lat belong to the cell with index int((90.0 - lat) / dLat)
    # except lat = -90 which belongs to cell with index int((90.0 - lat) / dLat) - 1
    # and extend borders by nBlocks each sides

    #
    # # note that at the North and South pole vicinity we have to be creative
    # # and copy the rows similar if -180 longitude is close
    # #  flags that control this abnormal situation
    # addNorth = idxNorth < 0
    # addSouth = idxSouth >= numLat
    #
    # addWest = idxWest < 0
    # addEast = idxWest >= numLon
    #

    # allocate memory
    sizeLat = numLat + 2* nBlocks
    sizeLon = numLon + 2* nBlocks

    # row -> latitudes
    # col -> longitudes

    demSubset = numpy.zeros((sizeLat, sizeLon))

    demSubset[nBlocks:sizeLat-nBlocks, nBlocks:sizeLon-nBlocks] = dem
    # add to the North and South
    demSubset[:nBlocks,          nBlocks:numLon + nBlocks] = numpy.roll(dem[nBlocks - 1::-1, :], int(numLon / 2), axis=1)
    demSubset[numLat + nBlocks:, nBlocks:numLon + nBlocks] = numpy.roll(dem[:numLat - 1 - nBlocks:-1, :], int(numLon / 2), axis=1)

    # add to the East and West
    demSubset[:, :nBlocks] = demSubset[:, numLon:numLon + nBlocks]
    demSubset[:, numLon + nBlocks:] = demSubset[:, nBlocks:2 * nBlocks]

    mLat =   90.0 - numpy.arange(sizeLat) * dLat - dLat / 2. + nBlocks*dLon
    mLon = -180.0 + numpy.arange(sizeLon) * dLon + dLon / 2. - nBlocks*dLon

    return topography.topographyBlock(dem=demSubset, lons=mLon, lats=mLat, dlon=dLon, dlat=dLat)


def subsetDEM(dem, lonWest, lonEast, latSouth, latNorth, nBlocks=100):
    '''
    return a subset of DEM bounded by
    lower left  (lonWest, latSouth)
    upper right (lonEast, latNorth)

    Input
        dem[numLat, numLon] -- NumPy array, surface heights [m]
        lonWest, lonEast -- left and right longitudes       [deg]
        latSouth, latNorth -- lower and upper latitides     [deg]

    Output
        topography.topographyBlock object
    '''
    numLat, numLon, = dem.shape

    # compute DEM resolution
    dLat = 180.0 / numLat
    dLon = 360.0 / numLon

    numLon2 = int(numLon/2)
    # since DEM cells boundaries goes from noth to south 90, 90 - dLat, 90 -2*dLat, ...
    # latitude <lat> belong to the cell with index <int((90.0 - lat) / dLat)>
    # Special case <lat = -90> which belongs to cell with index <int((90.0 - lat) / dLat) - 1>
    # and extend borders by nBlocks each sides

    idxNorth = int((90.0 - latNorth) / dLat)

    # DEM cell should be above the lower boundary so
    # if latitude is coincide with DEM latitude grid
    #  the index is decreased by one
    idxSouth = int((90.0 - latSouth) / dLat)
    if abs(90.0 - latSouth - idxSouth*dLat) < 1e-8:
        idxSouth -= 1
    # python excludes the last index in expression [i1:i2]
    #  so add 1
    idxSouth += 1

    idxWest = int((lonWest + 180.0) / dLon)
    # DEM cell should be to the left at east boundary  so
    # if longitude  is coincide with DEM longitude grid
    #  the index is decreased by one
    idxEast = int((lonEast + 180.0) / dLon)
    if abs(lonEast + 180.0 - idxEast * dLon) < 1e-8:
        idxEast -= 1
    # python excludes the last index in expression [i1:i2]
    #  so add 1
    idxEast += 1


    # extend all boundaries by nBlocks
    idxNorth -= nBlocks
    idxSouth += nBlocks
    idxWest -= nBlocks
    idxEast += nBlocks

    # note that at the North and South pole vicinity we have to be creative
    # and copy the rows above 90 (below -90) from dem below 90 (above -90)
    # similar if -180 longitude is close
    # flags that control this possibility
    addNorth=idxNorth < 0
    addSouth=idxSouth > numLat

    addWest =idxWest < 0
    addEast =idxEast > numLon

    # allocate memory
    sizeLat = idxSouth - idxNorth
    sizeLon = idxEast - idxWest

    # we can have two type of cases
    # normal case lonEast > lonWest
    if lonWest > lonEast:
        # the array is filled in two steps
        # step 1 lonWest - -180,
        # step 2 -180 - lonEast
        sizeLon += numLon


    if dbg:
        print ("\nextended indices", sizeLat, sizeLon)
        print('sizeLat', sizeLat)
        print('sizeLon', sizeLon)

    demSubset = numpy.zeros( (sizeLat, sizeLon) )

    # mLat and mLon are for the DEM cell upper and left boundaries
    mLat =   90.0 - idxNorth*dLat -  numpy.arange(sizeLat) * dLat #- dLat / 2.
    mLon = -180.0 + idxWest *dLon +  numpy.arange(sizeLon) * dLon #+ dLon / 2.

    if dbg:
        print ("\nextended indices", sizeLat, sizeLon)
        print('idxNorth', idxNorth, mLat[ 0], '<>', latNorth + nBlocks*dLat)
        print('idxSouth', idxSouth, mLat[-1], '<>', latSouth - nBlocks*dLat)

        print('idxWest', idxWest, mLon[ 0],  '<>', lonWest- nBlocks*dLon)
        print('idxEast', idxEast, mLon[-1],  '<>', lonEast+ nBlocks*dLon)

    # <idx > are the indices for dem array
    #  however they can be unrealistic indicating that
    # data has to be extended through DEM array boundary
    # <indAdjusted.. >
    idxAdjNorth = 0       if addNorth else idxNorth
    idxAdjSouth = numLat  if addSouth else idxSouth

    idxAdjWest = 0      if addWest else idxWest
    idxAdjEast = numLon if addEast else idxEast

    # <dem... > are the indices for demSubset array
    # the index shows where dem data will be copied
    # in trivial case where all the data inside of dem
    # demNorth=0
    # demSouth=sizeLat
    # demEast =sizeLon
    # demWest =0
    # in the case it has to be extended to the South
    # width = idxSouth - numLat
    # demSouth = sizeLat - width

    demNorth = idxAdjNorth - idxNorth
    widthSouth = idxSouth - idxAdjSouth
    demSouth = sizeLat - widthSouth   # sizeLat - numLat if addSouth else sizeLat

    demWest = idxAdjWest - idxWest
    widthEast = idxEast- idxAdjEast
    demEast = sizeLon - widthEast

    if dbg:
        print ("\nadjustment", sizeLat, sizeLon)
        print("where %6s %6s %6s %6s"%('add', 'idxAdj', 'idx', "dem"))
        print('North %6r %6d %6d %6d'%(addNorth, idxAdjNorth, idxNorth, demNorth))
        print('South %6r %6d %6d %6d %6d'%(addSouth, idxAdjSouth, idxSouth, demSouth, widthSouth))

        print('West  %6r %6d %6d %6d'%(addWest, idxAdjWest, idxWest, demWest))
        print('East  %6r %6d %6d %6d %6d'%(addEast, idxAdjEast, idxEast, demEast, widthEast))

    # we can have two type of cases
    # normal case lonEast > lonWest
    idx = numpy.arange(idxAdjWest, idxAdjEast)
    idxMir = numpy.where(idx < numLon2, idx + numLon2, idx - numLon2)
    if addNorth:
        # copy rows if our extended area goes over the North
        # has to be mirrored  lon -> lon + 180
        demSubset[:demNorth, demWest:demEast] = dem[demNorth - 1::-1, idxMir]
        # demSubset[0, idx] = dem[demNorth-1, idxMir]
        # demSubset[1, idx] = dem[demNorth-2, idxMir]
        # ...
        # demSubset[demNorth-1, idx] = dem[0, idxMir]

    if addSouth:
        # copy rows if our extended area goes over the South
        # has to be mirrored  lon -> lon + 180
        demSubset[demSouth:, demWest:demEast] = dem[:-widthSouth-1:-1, idxMir]
        # demSubset[demSouth, idx] = dem[numLat-1, idxMir]
        # demSubset[demSouth+1, idx] = dem[numLat-2, idxMir]
        # ...
        # demSubset[sizeLat-1, idx] = dem[numLat-widthSouth, idxMir]

    if addWest:
        # copy columns if our extended area crosses -180
        # mirror west of -180 is the eastern part of dem
        # lon - > 360 - lon
        demSubset[demNorth:demSouth, :demWest] = dem[idxAdjNorth:idxAdjSouth, :numLon - 1 - demWest:-1]

        if addNorth:
            # copy rows if our extended area goes over the North
            # has to be mirrored  lon -> lon - 180
            # finally lon -> 180 - lon
            #  demWest < 1 degree
            demSubset[:demNorth, :demWest] = dem[demNorth - 1::-1, numLon2 - 1:numLon2 - 1 - demWest:-1]
            #  in lat index
            # demSubset[0, demWest-1] = dem[demNorth-1, numLon2]
            # demSubset[0, demWest-2] = dem[demNorth-1, numLon2-1]
            #  ...
            # demSubset[0, 0       ] = dem[demNorth-1, numLon2-demWest]
            # etc
            # demSubset[1, ...] = dem[demNorth-2, ...]
            # ...
            # demSubset[demNorth-1, ...] = dem[0, ...]

        if addSouth:
            # copy rows if our extended area goes over the South
            demSubset[demSouth:, :demWest] = dem[:-widthSouth-1:-1, numLon2 - 1:numLon2 - 1 - demWest:-1]

        # copy columns if our extended area crosses -180
    if addEast:
        demSubset[demNorth:demSouth, demEast:] = dem[idxAdjNorth:idxAdjSouth, :widthEast]
        if addNorth:
            # copy rows if our extended area goes over the North
            # has to be mirrored  lon -> lon - 180
            # finally lon -> 180 - lon
            #  demWest < 1 degree
            demSubset[:demNorth, demEast:] = dem[demNorth - 1::-1, numLon2:numLon2+widthEast]
            #  in lat index
            # demSubset[0, demWest-1] = dem[demNorth-1, numLon2]
            # demSubset[0, demWest-2] = dem[demNorth-1, numLon2-1]
            #  ...
            # demSubset[0, 0       ] = dem[demNorth-1, numLon2-demWest]
            # etc
            # demSubset[1, ...] = dem[demNorth-2, ...]
            # ...
            # demSubset[demNorth-1, ...] = dem[0, ...]

        if addSouth:
            # copy rows if our extended area goes over the South
            demSubset[demSouth:, demEast:] = dem[:-widthSouth-1:-1, numLon2:numLon2+widthEast]

    else:
        # the array is filled in two steps
        # -180 longitude index in demSubset
        dem180 = numLon - idxAdjWest

        # step 1 west of 180
        demSubset[demNorth:demSouth, demWest:dem180] = dem[idxAdjNorth:idxAdjSouth, idxAdjWest:]
        idx = numpy.arange(idxAdjWest, numLon)
        idxMir = numpy.where(idx < numLon2, idx + numLon2, idx - numLon2)
        if addNorth:
            # copy rows if our extended area goes over the North
            demSubset[:addNorth, demWest:dem180] = dem[demNorth - 1::-1, idxMir]
        if addSouth:
            demSubset[demSouth:, demWest:dem180] = dem[:-widthSouth - 1:-1, idxMir]

        # step 2 east of 180
        demSubset[demNorth:demSouth, dem180:] = dem[idxAdjNorth:idxAdjSouth, :idxAdjEast]
        idx = numpy.arange(idxAdjEast)
        idxMir = numpy.where(idx < numLon2, idx + numLon2, idx - numLon2)
        if addNorth:
            # copy rows if our extended area goes over the North
            demSubset[:addNorth, dem180:] = dem[demNorth - 1::-1, idxMir]

        # copy rows if our extended area goes over the South
        if addSouth:
            demSubset[demSouth:, dem180:] = dem[:-widthSouth - 1:-1, idxMir]

    #  ATTENTION
    # the sample datase use negative heights for ocean floor
    # so set them to 0 over ocean
    demSubset[demSubset < 0] = 0.

    return topography.topographyBlock(
                dem=demSubset, lons=mLon, lats=mLat, dlon=dLon, dlat=dLat)


def topCalc(dem, latSouth, latNorth, lonWest, lonEast, delLon=1, delLat=1,
    saveFig=False, toPlot=False, nBlocks=100, add_aux=False, nZen=19, nHz=16, maxDist=20000):
    '''
    Topology parameter calculation and netcdf generation

    Call
        nxGridX, nxGridY = topCalc(dem, numLat, numLon, saveFig=False,
            toPlot=False, nBlocks=100, add_aux=False, nZen=19)

    Inputs
        dem -- NumPy array, grid from readDEM()
        latSouth, latNorth, lonWest, lonEast latitudes and longitudes defining computational area

        maxDist --max distance [m] along which the horizon angle is estimated

    Outputs
        nxGridX, nxGridY -- integer sets (list with unique elements),
            this is necessary for naming convention used and
            referenced throughout module; nxGridY is in descending order

    Keywords
        delLon=1, delLat=1 the computational area resolution [deg]
        toPlot -- boolean, plot the topology (height, V) maps
        saveFig -- boolean, save height maps to PNG files
        nBlocks -- the number of DEM cell in latitudinal and longitudinal directions the
                -- computational area extended to read DEM database

        add_aux -- boolean, control information added to netCDF file
        nZen -- int, number of zenith angle that will be used in
            table preparations
        nHz -- int, number of azimuthal direction along which the horizon angles are estimated
    '''

    # check on input
    isError = False

    if abs(latSouth) > 90.0:
        print ("invalid latSouth: must be in range [-90, 90] : %f"%(latSouth))
        isError = True

    if abs(latNorth) > 90.0:
        print ("invalid latNorth: must be in range [-90, 90] : %f"%(latNorth))
        isError = True

    if latNorth <= latSouth:
        print ("invalid combination of latSouth and latNorth: latSouth < latNorth : %f"%(latSouth, latNorth))
        isError = True

    if abs(lonWest) > 180.0:
        print ("invalid lonWest: must be in range [-180, 180] : %f"%(lonWest))
        isError = True

    if abs(lonEast) > 180.0:
        print ("invalid lonEast: must be in range [-180, 180] : %f"%(lonEast))
        isError = True

    # if lonEast <= lonWest:
    #     print ("WARNING: -180 longitude is between lonWest and lonEast: %f"%(lonWest, lonEast))

    if isError:
        print ("ERROR: invalid input detected. please correct")
        exit(101)

    top = subsetDEM(dem, lonWest, lonEast, latSouth, latNorth, nBlocks=nBlocks)
    dx, dy = topography.haversineStep(min(abs(latSouth), abs(latNorth)), top.dlon, top.dlat)
    if maxDist / min(dx, dy) > nBlocks:
        print("")
        print("WARNING: increase nBlocks up to handle horizon angle computation at boundaries %d " % (
            numpy.ceil(maxDist / min(dx, dy))))
        print("")
        plt.pause(2)

    # initialize computational area grid
    # given user defined area with resolution delLon, delLat
    #  to handle situation lonWest > lonEast:
    if lonWest > lonEast:
        lonEast += 360.

    gB = topography.gridBlock(
        x1=lonWest, x2=lonEast, y1=latSouth, y2=latNorth, delX=delLon, delY=delLat)


    print("x=", gB.x1, gB.x2)
    print("y=", gB.y1, gB.y2)
    print("dLat, dLon=", top.dlat, top.dlon)

    if dbg: print('(gridX - top.lons[0]) / top.dlon', (gB.gridX.grid - top.lons[0]) / top.dlon)
    if dbg: print('(top.lats[0] - gridY) / top.dlat', (top.lats[0] - gB.gridY.grid ) / top.dlat)
    if dbg: print('gridX', gB.gridX)
    if dbg: print('gridY', gB.gridY)

    if toPlot:
        plt.figure(2, figsize=(12, 5))
        plt.figure(1)
        # draw topography subset extracted from DEM
        top.show(isStatic=False)
        blk=None
        # draw user defined area  over topography
        plt.plot([lonWest, lonEast, lonEast, lonWest, lonWest],
                 [latSouth, latSouth, latNorth, latNorth, latSouth], color='r')

    nxGridY, nxGridX = [], []

    # return azimuth angle to compute horizon angles [rad]
    horAzmAngle = topography.calcHorizontalAngle(nHz)

    print ("==============================================================================")
    print ("GRID cell")
    for ix, cx1, cx2 in gB.gridX:

        N1 = gB.gridX.floor(ix  , top.lons[0], top.dlon)
        N2 = gB.gridX.ceil (ix+1, top.lons[0], top.dlon)
        nxGridX.append(N1)

        for iy, cy1, cy2 in gB.gridY:
            # indices N2 and M2 indicates the upper boundary in the Pythonic sense
            #  DEM latitudinal grid and computational grids are in opposite directions
            M2 = gB.gridY.ceil (iy,     top.lats[0], top.dlat)
            M1 = gB.gridY.floor(iy + 1, top.lats[0], top.dlat)

            nxGridY.append(M1)
            ncFile = 'horizon_%05d_%05d.nc' % (N1, M1)
            if os.path.isfile(ncFile):
                print("SKIPPED: Processing %d, %d : " % (ix, iy))
                continue

            # compute the distance corresponding to DEM cell size in [deg]
            dx, dy = topography.haversineStep(
                (cy1+cy2)/2, top.dlon, top.dlat)

            print(70*"-")
            print("Processing %d, %d : "%(ix, iy))
            print("   x1, x2 = %.3f %.3f %.3f %.3f %d %d"%(cx1, cx2, top.lons[N1], top.lons[N2], N1, N2))
            print("   y1, y2 = %.3f %.3f %.3f %.3f %d %d"%(cy1, cy2, top.lats[M1], top.lats[M2], M1, M2))
            print("   dx, dy = %.1f %.1f"%(dx, dy))

            if toPlot:
                plt.figure(1)
                if blk:
                    blk.remove()
                # draw computational cell over topography
                lines=plt.plot([cx1, cx2, cx2, cx1, cx1],
                         [cy1, cy1, cy2, cy2, cy1], color='m')
                blk =lines.pop(0)
                if saveFig:
                    plt.savefig('figure_spot.%d.%d.png' % (ix, iy))
                plt.draw()
                plt.pause(.5)

            # To avoid computation for the flat terrain
            #  delta is max allowed altitude difference in [m]
            featureFound=top.check(N1, N2, M1, M2, delta=1)

            if featureFound:

                # If computational area is not large topography.preProcInit can be called just once
                topography.preProcInit(dx, dy, top.dlon, top.dlat, horAzmAngle, maxDist=maxDist, debug=False)

                weight, tanSlopeAngle, slopeAspect, horAngle, = topography.Preprocess(
                            dx, dy, N1, N2, M1, M2, top,
                            topography.gridCell(x1=cx1, x2=cx2, y1=cy1, y2=cy2), horAzmAngle, debug=False)
            else:
                if dbg: print ('computations skipped: no terrain features found')

                # default for flat surface
                weight = numpy.ones((M2 - M1, N2 - N1))
                horAngle = numpy.full((M2 - M1, N2 - N1, len(horAzmAngle)),fill_value=numpy.pi/2.)
                tanSlopeAngle = numpy.zeros_like(weight)
                slopeAspect = numpy.zeros_like(weight)

            # mean slope angle without shadow mask
            oxAv = numpy.ones((nZen, nHz))

            # mean slope angle with shadow mask
            oxAvMask = numpy.ones((nZen, nHz))

            # solar zenith grid for table
            solzen = numpy.deg2rad(numpy.linspace(0, 90, num=nZen, endpoint=True))
            solzen[-1]=numpy.deg2rad(89.9)

            # sky view factor
            V = numpy.ones_like(weight)

            # computation of mean slope angle
            if featureFound:
                for ja, sunA in enumerate(horAzmAngle):
                    curHor = horAngle[:, :, ja]

                    for kk, sunZ in enumerate(solzen):
                        # shadow mask
                        mask = numpy.int32(curHor > sunZ)
                        if toPlot:
                            oxAv[kk, ja] = numpy.mean(
                                weight * (1. + tanSlopeAngle * numpy.tan(sunZ) *
                                numpy.cos(sunA - slopeAspect)))/numpy.mean(weight)

                        oxAvMask[kk, ja] = numpy.mean(mask *
                            weight * (1. + tanSlopeAngle * numpy.tan(sunZ) *
                            numpy.cos(sunA - slopeAspect)))/numpy.mean(weight)

                # sky view factor
                for an, angle in enumerate(horAzmAngle):
                    tanTF = -numpy.arctan(1. / tanSlopeAngle / numpy.cos(angle - slopeAspect))
                    tanTF[tanSlopeAngle == 0] = numpy.pi / 2.
                    tanTF[tanTF < 0] += numpy.pi
                    V += numpy.sin(horAngle[:, :, an] +
                        numpy.pi / 2. - tanTF) ** 2

                V /= 16.

            print('   Writing {}'.format(ncFile), end='')

            with Dataset(ncFile, 'w', format='NETCDF4_CLASSIC') as fid:
                nlat, nlon, nazm,  =  horAngle.shape
                # define axis size
                fid.createDimension('azm', nazm)
                fid.createDimension('zen', nZen)

                # create azm axis
                tVar = fid.createVariable('azm', 'f4', ('azm',))
                tVar.long_name = 'horizon azimuth'
                tVar.units = 'rad'
                tVar[:] = horAzmAngle

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

                    # create tanSlopeAngle
                    tVar = fid.createVariable(
                        'tanA', 'f4', ('lat', 'lon'))
                    tVar.standard_name = 'slope angle'
                    tVar.long_name = 'slope angle'
                    tVar.units = 'rad'
                    tVar[:] = numpy.rad2deg(tanSlopeAngle)

                    # create slopeAspect
                    tVar = fid.createVariable(
                        'S', 'f4', ('lat', 'lon'))
                    tVar.standard_name = 'slope aspect'
                    tVar.long_name = 'slope aspect'
                    tVar.units = 'rad_north'
                    tVar[:] = numpy.rad2deg(slopeAspect)

                    # create horAngle
                    tVar = fid.createVariable(
                        'horz', 'f4', ('lat', 'lon', 'azm'))
                    tVar.standard_name = 'horizon angle'
                    tVar.long_name = 'horizon angle'
                    tVar.units = 'rad'
                    tVar[:] = horAngle

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
                tVar[:] = numpy.mean(weight*V)/numpy.mean(weight)

                tVar = fid.createVariable('Vav', 'f4')
                tVar.standard_name = 'weighted LW view factor'
                tVar.long_name = 'weighted LW view factor'
                tVar.units = 'none'
                tVar[:] = numpy.mean(weight*V * numpy.sqrt(1. + tanSlopeAngle ** 2))/numpy.mean(weight)

                tVar = fid.createVariable('InvCosA', 'f4')
                tVar.standard_name = \
                    'averaged inverse cosine of slope angle'
                tVar.long_name = \
                    'averaged inverse cosine of slope angle'
                tVar.units = 'unit'
                tVar[:] = numpy.mean(weight*numpy.sqrt(1. + tanSlopeAngle ** 2))/numpy.mean(weight)

            print(": DONE")

            if toPlot:
                # plot SW correction factor and
                #  LW view factor
                plt.figure(2)
                plt.clf()
                plt.subplot(121)
                plt.pcolor(numpy.rad2deg(horAzmAngle),
                    numpy.rad2deg(solzen),oxAvMask)
                kk=plt.colorbar()
                plt.xlabel('solar azimuth', fontsize=14)
                plt.tight_layout()
                plt.title('$f_{cor}^{masked}$', fontsize=14)

                plt.subplot(122)
                plt.pcolor(top.lons[N1:N2 + 1], top.lats[M1:M2 + 1], V)
                kk = plt.colorbar()
                kk.set_label('V')
                plt.xlabel('Longitude', fontsize=14)
                plt.ylabel('Latitude', fontsize=14)
                plt.title('V=%.3f' % (numpy.mean(weight*V)/numpy.mean(weight)))
                if saveFig:
                    plt.savefig('./figure_f_cor.%d.%d.png'%(ix,iy))
                plt.draw()
                plt.pause(.1)

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
    fname = 'horizon_{:05d}_{:05d}.nc'.format(nxGridX[0], nxGridY[0])
    with Dataset(fname) as fid:
        zenData = fid.variables['zen'][:]
        horAzmAngle = fid.variables['azm'][:]

    nHorizon = len(horAzmAngle)
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
        tVar[:] = numpy.rad2deg(horAzmAngle)

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
    parser.add_argument('--gridBox', '-grid', type=float, nargs=4, \
        help='Grid box bounday:\n The lower left corner geographical longitude and latitude, ' + \
             'The upper right corner geographical longitude and latitide\n' + \
              'Longitude: -180:+180; Latitude: -90: +90', default =None)

    parser.add_argument('--n_lat', '-n_lat', type=int, default=10800, \
        help='Number of latitude grid points for the ETOPO1 ' + \
        'grid reshape.')
    parser.add_argument('--n_lon', '-n_lon', type=int, default=21600, \
        help='Number of longitude grid points for the ETOPO1 ' + \
        'grid reshape.')
    parser.add_argument('--infile', '-i', type=str, \
        default='ETOPO1_Ice_c_gdal.grd', \
        help='ETOPO1 1 Arc-Minute Global Relief Model output netCDF')
    parser.add_argument('--savefig', '-s', action='store_true', default=False, \
        help='Save figures of V and altitude maps as PNG files.')
    parser.add_argument('--plot', '-p', action='store_true', default=True, \
        help='Generate V and altitude plots in addition to ' + \
        'the output netCDF files and either save (--savefig) or ' + \
        'display them.')
    parser.add_argument('--n_zenith', '-nz', default=19, type=int, \
        help='Number of zenith angles')
    parser.add_argument('--n_horizons', '-nh', default=16, type=int, \
        help='Number of horizons.')
    parser.add_argument('--n_blocks', '-nb', default=20, \
        help='The DEM projection on the user defined area is extended by n_blocks of DEM cells in all directions')
    parser.add_argument('--auxiliary', '-a', action='store_true', default=False, \
        help='Add auxiliary information into netCDF.')
    args = parser.parse_args()

    if args.gridBox is None:
        print ("User must enter computational grid box: the lower left corner geographical longitude and latitude,\n"\
             "the upper right corner geographical longitude and latitide.\n Example: prepareTopographyExample.py -grid -72 -22 -63 -12")
        exit(101)

    elif len(args.gridBox) != 4:
        print ("User must enter all 4 values for computational grid box: \n"
               "the lower left corner geographical longitude and latitude, "
               "the upper right corner geographical longitude and latitide")
        print ("Current input: ", args.gridBox)

    demGrid = readDEM(args.infile, args.n_lat, args.n_lon)
    print ("PrepareTopographyExample args:")
    print( '    gridBox', 	args.gridBox)
    print( '    saveFig', 	args.savefig)
    print( '    toPlot', 	args.plot)
    print( '    nBlocks', 	args.n_blocks)
    print( '    add_aux', 	args.auxiliary)
    print( '    nZen', 	args.n_zenith)
    print( '    nHz', 	args.n_horizons)

    nxx, nxy = topCalc(demGrid, args.gridBox[1], args.gridBox[3],
        args.gridBox[0], args.gridBox[2],
        delLon=1, delLat=1,  
        saveFig=args.savefig, toPlot=args.plot,
        nBlocks=args.n_blocks, add_aux=args.auxiliary,
        nZen=args.n_zenith, nHz=args.n_horizons, maxDist=20000)

    prepareSW(nxx, nxy)
    prepareLW(nxx, nxy)
