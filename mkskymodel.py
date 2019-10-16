#! /usr/bin/env python
from __future__ import division, print_function

import csv
import sys

from astropy.coordinates import SkyCoord
import astropy.units as units
import numpy as np

import astrobits.sourcecounts as sourcecounts
import astrobits.skymodel as skymodel

#center = SkyCoord('23h03m17.8869s -24d56m09.097s', unit=(units.hourangle, units.degree))
center = SkyCoord('0 -27', unit=(units.hourangle, units.degree))
ra0 = center.ra.radian
dec0 = center.dec.radian
lowerthreshold = 50e-6

bins = np.linspace(np.log10(lowerthreshold), np.log10(75), 10000000)
logS = (bins[:-1] + bins[1:]) / 2
ds = 10**bins[1:] - 10**bins[:-1]
logdN = sourcecounts.logdNdS(logS) + np.log10(ds)
dN = 10**logdN

radius = 15  # degrees
solidangle = 2 * np.pi * (1 - np.cos(np.radians(radius)))
print("Solid angle:", solidangle)
N = int(np.sum(dN) * solidangle)
print("Catalogue sources:", N)

fluxes = 10**np.random.choice(logS, size=N, p=dN/np.sum(dN))
alphas = np.random.normal(-0.8, 0.1, size=N)
# alphas = np.zeros(N)

minz = np.cos(np.radians(radius))
phis = 2 * np.pi * np.random.uniform(size=N)
thetas = np.arccos(np.random.uniform(low=minz, high=1, size=N))

# Rotate North pole into desired Ra, Dec

# First, convert to cartesian coordinates on unit sphere
xyz = np.array([np.sin(thetas) * np.cos(phis), np.sin(thetas) * np.sin(phis), np.cos(thetas)])

# Rotate North pole down to Dec, rotated around y axis
omega = np.pi / 2 - dec0
xyz = np.matmul([
    [np.cos(omega), 0, np.sin(omega)],
    [0, 1, 0],
    [-np.sin(omega), 0, np.cos(omega)],
    ], xyz)

# Rotate around z axis by Ra0
omega = ra0
xyz = np.matmul([
    [np.cos(omega), -np.sin(omega), 0],
    [np.sin(omega), np.cos(omega), 0],
    [0, 0, 1],
], xyz)

# Convert cartesian coordinates back to spherical
thetas = np.arccos(xyz[2])
phis = np.arctan2(xyz[1], xyz[0])

# Convert phis, thetas into Ra, Dec
ras = phis
decs = np.pi / 2 - thetas

coords = SkyCoord(ras, decs, unit=(units.radian, units.radian))

catalogue = np.empty((len(ras), 5)) # [[ra, dec, refreq, stokesI, alpha]]
print(catalogue.shape)
catalogue[:, [True, True, False, True, True]] = np.array([ras, decs, fluxes, alphas]).T
catalogue[:, 2] = 154.5E6
np.save('catalogue.npy', catalogue)

# with open('catalogue.txt', 'w') as f:
#     print("skymodel fileformat 1.1", file=f)
#     for i in range(len(coords)):
#         if i % 1000 == 0:
#             print("%d/%d" % (i+1, len(coords)), end="\r")
#             sys.stdout.flush()
#         f.write("source {\n")
#         f.write("  name \"source-%d\"\n" % i)
#         f.write("  component {\n")
#         f.write("    type point\n")
#         f.write("    position %s\n" % coords[i].to_string("hmsdms"))
#         f.write("    sed {\n")
#         f.write("      frequency 154.5 MHz\n")
#         f.write("      fluxdensity Jy %f 0 0 0\n" % fluxes[i])
#         f.write("      spectral-index { %f 0 }\n" % alphas[i])
#         f.write("    }\n")
#         f.write("  }\n")
#         f.write("}\n")

# # Create models
# models = [
#     skymodel.Model('source-%d' % i, [
#         skymodel.SEDComponent(coords[i], [154.5E6, fluxes[i], alphas[i], 0]),
#     ])
#     for i in range(len(coords))
# ]

# print("Writing...")

# with open('catalogue.txt', 'w') as f:
#     skymodel.writeto(f, models)
