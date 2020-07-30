import numpy

# Earth equatorial radius
Re = 6378137.0

# flattening parameter
f  = 1.0 / 298.2572235630

# polar radius
Rp = Re * (1.0 - f)

# eccentricity
e  = (Re**2 - Rp**2)/(Re**2)

def geod2cart(rlat, rlon, height):
    """
    Geodetic to Cartesian coordinate conversion

    Call
        cart = geod2cart(rlat, rlon, height)

    Input
        rlat -- NumPy float array of Geodetic latitudes
        rlon -- NumPy float array of Geodetic longitudes
        height -- NumPy float array of heights (m)

    Output
        cart -- tuple, x, y, and z coordinates in meters
    """

    flatfn = (2.0 - f) * f
    funsq = (1.0 - f)**2
    # rlat = numpy.deg2rad(lat)
    # rlon = numpy.deg2rad(lon)

    gd = Re/numpy.sqrt(1.0 - flatfn*numpy.sin(rlat)**2)

    cart =  numpy.array([numpy.cos(rlat)*numpy.cos(rlon)*(gd + height),
            numpy.cos(rlat)*numpy.sin(rlon)*(gd + height),
            numpy.sin(rlat)*(gd*funsq + height)])
    return cart


def Reff(lat):
    """
    Calculate effective radius given a latitude
    """

    a_earth = 6378137.0 #     ! semi - major    axis(m)
    flatt = 0.003352811 #     ! flattening
    gm_ratio = 0.003449787 #  ! gravitational    ratio

    sinlat2 = numpy.sin(numpy.deg2rad(lat))**2
    return a_earth / (1.0 + flatt + gm_ratio - 2.0 * flatt * sinlat2)


def curvature(lat, lon, theta):
    '''
    :param lat:    Surface point latitude
    :param lon:    Surface point longitude
    :param theta: Azimuth direction
    :return:
    r_coc centre of curvature
    roc Radius of curvature
    '''
    rlat = numpy.deg2rad(lat)
    rlon = numpy.deg2rad(lon)
    x = geod2cart(rlat, rlon, 0.0)
    atoc2 = (Re / Rp) ** 2
    xy2 = x[0] ** 2 + x[1] ** 2
    z2 = (atoc2 * x[2]) ** 2
    cphi2 = xy2 / (xy2 + z2)
    sphi2 = z2 / (xy2 + z2)
    r_EW = 1.0 / numpy.sqrt( cphi2 + (sphi2/atoc2) )  # normalised wrt Re
    r_NS = r_EW**3 / atoc2                         # normalised wrt Re
    roc = Re / ( (numpy.cos(theta)**2/r_NS) + (numpy.sin(theta)**2/r_EW) )

    n = numpy.array([numpy.cos(rlat) * numpy.cos(rlon),
         numpy.cos(rlat) * numpy.sin(rlon),
         numpy.sin(rlat)])

    r_coc = x - roc * n
    return roc, r_coc
