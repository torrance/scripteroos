#! /usr/bin/env python
from __future__ import print_function, division

import argparse
import os.path

from astropy.io import fits
from astropy.wcs import WCS
from casacore.tables import table
from numba import njit, float64
import numpy as np


@njit([(float64[:, :], float64[:], float64[:], float64[:])])
def paint(data, xs, ys, fluxes):
    weights = np.zeros((2, 2), dtype=float64)

    for x, y, flux in zip(xs, ys, fluxes):
        x_int, y_int = int(x), int(y)

        # Calculate proportionate weights for the 4 nearby pixels
        for i in [0, 1]:
            for j in [0, 1]:
                xi, yj = x_int + i, y_int + j
                weights[i, j] = 1 / np.sqrt((x - xi)**2 + (y - yj)**2)

        weights /= weights.sum()

        # Set flux weighted by nearby pixels
        for i in [0, 1]:
            for j in [0, 1]:
                xi, yj = x_int + i, y_int + j
                data[yj, xi] += weights[i, j] * flux


def tukey(n, N, alpha):
    n = n + N/2  # Tukey defined on [0, N], shift for [-N/2, N/2]
    window = np.zeros_like(n)

    # Left
    idx = np.all([n >= 0, n < (alpha * N) / 2], axis=0)
    window[idx] = 0.5 * (1 + np.cos(np.pi * ((2 * n[idx]) / (alpha * N)) - 1))

    # Middle
    idx = np.all([n >= (alpha * N) / 2, n <= N * (1 - alpha / 2)], axis=0)
    window[idx] = 1

    # Right
    idx = np.all([n > N * (1 - alpha / 2), n <= N], axis=0)
    window[idx] = 0.5 * (1 + np.cos(np.pi * ((2 * n[idx]) / (alpha * N) - 2 / alpha + 1)))

    return window



parser = argparse.ArgumentParser()
parser.add_argument('--models', nargs='+', required=True)
parser.add_argument('--prefix', default='')
parser.add_argument('--templateimage', required=True)
parser.add_argument('--templatemset', required=True)
parser.add_argument('--channels-out', type=int, required=True)
args = parser.parse_args()

modelname = args.prefix + '-' if args.prefix else ''

models = []
for model in args.models:
    models.append(np.load(model))
#models[1][:, 3] = models[1] [:, 3] * 1000
#print("AMPLIFIED!!!!!!", models[1][:, 3].sum())
model = np.concatenate(models)
print("Total flux:", model[:, 3].sum())
print("Processing model (%d componenets)" % len(model))

ras, decs = model[:, 0], model[:, 1]

template = fits.open(args.templateimage)[0]
naxis1, naxis2 = template.header['NAXIS1'], template.header['NAXIS2']
ra0, dec0 = np.radians(template.header['CRVAL1']), np.radians(template.header['CRVAL2'])
# cdelt1, cdelt2 = np.radians(template.header['CDELT1']), np.radians(template.header['CDELT2'])
# crpix1, crpix2 = template.header['CRPIX1'], template.header['CRPIX2']

# Get channel requencies and reduce to channels-out
freqs = table(args.templatemset + '::SPECTRAL_WINDOW').getcell('CHAN_FREQ', 0)
assert(len(freqs) % args.channels_out == 0)
freqs = np.mean(np.reshape(freqs, (args.channels_out, -1)), axis=1)

# Calculate pixel coordinates of sources
coords = np.zeros((len(ras), 4))
coords[:, 0] = np.degrees(ras)
coords[:, 1] = np.degrees(decs)
wcs = WCS(template.header)
xs, ys, _, _ = wcs.wcs_world2pix(coords, 0).T

# Apply Tukey window
r = np.arccos(np.sin(model[:, 1]) * np.sin(dec0) + np.cos(model[:, 1]) * np.cos(dec0) * np.cos(model[:, 0] - ra0))
#model[:, 3] = model[:, 3] * tukey(r, N=np.radians(10), alpha=0.71)
#! /usr/bin/env python
from __future__ import print_function, division

from aplpy import FITSFigure



fig = FITSFigure('#! /usr/bin/env python
        from __future__ import print_function, division

        from aplpy import FITSFigure



        fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure('#! /usr/bin/env python
            from __future__ import print_function, division

            from aplpy import FITSFigure



            fig = FITSFigure(m.model[:, 3] = model[:, 3] * tukey(r, N=np.radians(1.5), alpha=0.71)

# Remove sources with 0 flux
zeroed = model[:, 3] == 0
print("Removing %d zero flux sources" % sum(zeroed))
model = model[~zeroed]
xs, ys = xs[~zeroed], ys[~zeroed]

# Remove out of bound sources
oob = np.any([xs < 0, xs >= naxis1, ys < 0, ys >= naxis2], axis=0)
print("Removing %d OOB sources" % sum(oob))

model = model[~oob]
xs, ys = xs[~oob], ys[~oob]


for i, freq in enumerate(freqs):
    fluxes = model[:, 3] * (freq / model[:, 2])**model[:, 4]
    print(fluxes)

    # Paint sources onto model image
    data = np.zeros(template.data.shape)
    paint(data[0, 0], xs, ys, fluxes)

    # Save image
    modelimage = template.copy()
    modelimage.header['CRVAL3'] = freq
    modelimage.data = data

    if len(freqs) > 1:
        modelimage.writeto("%s%04d-model.fits" % (modelname, i), overwrite=True)
    else:
        modelimage.writeto("%smodel.fits" % modelname, overwrite=True)


