#! /usr/bin/env python
from __future__ import division, print_function

import argparse

from astropy.io import fits
import matplotlib.pyplot as plt
from numba import njit, float64, prange
import numpy as np
from scipy.signal import correlate

from astrobits.coordinates import fits_to_radec


@njit([float64[:, :](float64[:, :], float64[:, :], float64[:, :], float64[:, :])], parallel=True)
def dister(ras1, decs1, ras2, decs2):
    N, M = ras1.shape
    corr = np.zeros((N, M))

    sindecs1 = np.sin(decs1)
    sindecs2 = np.sin(decs2)
    cosdecs1 = np.cos(decs1)
    cosdecs2 = np.cos(decs2)

    for k in prange(N):
        delta_i = k - N // 2

        for l in range(M):
            delta_j = l - M // 2

            i = N // 2 - 1
            j = M // 2 - 1
            corr[k, l] = np.arccos(
                sindecs1[i, j] * sindecs2[i - delta_i, j - delta_j]
                + cosdecs1[i, j] * cosdecs2[i - delta_i, j - delta_j] * np.cos(ras1[i, j] - ras2[i - delta_i, j - delta_j])
            )

    return corr


parser = argparse.ArgumentParser()
parser.add_argument('--images', nargs='+', required=True)
parser.add_argument('--galaxydensity', required=True)
args = parser.parse_args()

corrs = []
for image in args.images:
    image1 = fits.open(image)[0]
    image2 = fits.open(args.galaxydensity)[0]


    # Sanity checks
    assert(image1.data.shape == image2.data.shape)
    ras1, decs1 = fits_to_radec(image1)
    ras2, decs2 = fits_to_radec(image2)
    assert(np.allclose(ras1, ras2))
    assert(np.allclose(decs1, decs2))

    dists = np.degrees(dister(ras1[0, 0], decs1[0, 0], ras2[0, 0], decs2[0, 0]))

    tmp = np.empty(image1.data[0, 0].shape)
    tmp[:] = image1.data[0, 0]
    image1 = tmp
    tmp = np.empty(image2.data[0, 0].shape)
    tmp[:] = image2.data[0, 0]
    image2 = tmp

    Ns = correlate(np.ones_like(image1), np.ones_like(image2), mode='same')

    corr = correlate(image1 - np.mean(image1), image2 - np.mean(image2), mode='same')
    corr = corr / Ns
    corrs.append(corr)

# plt.subplot(1, 4, 1)
# plt.imshow(image1)
# plt.subplot(1, 4, 2)
# plt.imshow(image2)
# plt.subplot(1, 4, 3)
# plt.imshow(corr)

    bins = np.linspace(0, 0.4, 100)
    indices = np.digitize(dists, bins)

    mids, values = [], []
    for i, (left, right) in enumerate(zip(bins[:-1], bins[1:])):
        mids.append((left + right) / 2)
        values.append(np.mean(corr[indices == i]))

# np.save('radialCCF.npy', np.array([mids, values]))

    plt.plot(mids, values, label=image)
    plt.xlim([0, 0.4])
    plt.xlabel("Degrees")
    plt.legend()
    # plt.yscale('log', nonposy='clip')

plt.show()







