#! /usr/bin/env python
"""WFC3 wavelength calibration
"""
from __future__ import print_function, division
import numpy as np


def WFC3Dispersion(xc, yc, subarray=256):
    """convert pixel coordinate to wavelength. Method and coefficient
    adopted from Kuntschner et al. (2009), Wilkins et al. (2014)

    =====
    xc -- X coordinate of direct image centroid
    yc -- Y coordinate of direct image centroid
    return wavelength mapping of x coordinate in angstrom

    Note that direct image and spectral image should be taken with the
    same aperture. If not, please adjust the centroid measurement
    according to table in
    http://www.stsci.edu/hst/observatory/apertures/wfc3.html

    e.g., if direct image is taken with aperture IRSUB256 and spectral
    image is taken wit aperture GRISM256, xc0 and yc0 are measurements
    in direct image. xc, and yc should be offset with
    xc = xc0 - (522 - 410)
    yc = yc0 - (522 - 532)
    """
    DLDP0 = [8949.40742544, 0.08044032819916265]
    DLDP1 = [44.97227893276267,
             0.0004927891511929662,
             0.0035782416625653765,
             -9.175233345083485e-7,
             2.2355060371418054e-7,  -9.258690000316504e-7]
    # calculate field dependent dispersion coefficient
    p0 = DLDP0[0] + DLDP0[1] * xc
    p1 = DLDP1[0] + DLDP1[1] * xc + DLDP1[2] * yc + \
         DLDP1[3] * xc**2 + DLDP1[4] * xc * yc + DLDP1[5] * yc**2
    dx = np.arange(1014) - xc
    wavelength = (p0 + dx * p1)
    if subarray < 1014:
        i0 = (1014 - subarray) // 2
        wavelength = wavelength[i0: i0 + subarray]
    return wavelength
