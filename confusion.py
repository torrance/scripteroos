#! /usr/bin/env python
from __future__ import print_function, division

from astropy.io import fits
from numba import njit, float64, prange
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from astrobits.sourcecounts import dNdS


@njit([float64[:](float64[:], float64[:, :], float64, float64)], parallel=True)
def R(xs, beam, pixelsize, minS):
    """
    pixelsize: area of cell in steradians
    """
    Rxs = np.empty_like(xs)
    beam = beam.flatten()

    for i in prange(0, len(xs)):
        x = xs[i]
        S = np.absolute(x / beam)
        _dNdS = dNdS(S) * pixelsize  # differential source count per pixel
        _dNdS[S < minS] = 0
        _dNdS[S > 75] = 0
        Rxs[i] = np.sum(_dNdS / beam)

    return Rxs


def p(Rxs):
    rws = np.fft.fft(Rxs)
    return np.exp(rws - rws[0])


def P(pws):
    return np.fft.ifft(pws)


# Set up Gaussian beam
xs, ys = np.meshgrid(range(-200, 201), range(-200, 201))
scale = 1  #arcseconds
xs, ys = xs * scale, ys * scale
FWHM = 70  # arcseconds
sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
beam = np.exp(-(xs**2 + ys**2) / (2 * sigma**2))

pixelsize = scale / 3600  # degrees
pixelarea = pixelsize**2 * (np.pi**2 / 180**2)  # steradians


for minS in [0.15e-3, 0.1e-3, 0.05e-3, 0.025e-3, 0.01e-3, 0.005e-3, 0.001e-3]: # [1e-3, 0.75e-3, 0.5e-3, 0.25e-3, 0.1e-3]:
    # Calculate R(x) for a logarithmic spacing of xs
    print("Calculating Rlogxs")
    logxs = np.logspace(-7, 1.7, num=10000)
    dlogxs = logxs[1:] - logxs[:-1]
    logxs = (logxs[1:] + logxs[:-1]) / 2
    Rlogxs = R(logxs, beam, pixelarea, minS=minS)

    # Interpolate Rxs onto regular xs grid
    print("Interpolating onto linear array")
    xs = np.linspace(logxs[0], logxs[-1], 1e8)
    dx = (xs[-1] - xs[0]) / (len(xs) - 1)
    Rxs = interp1d(logxs, Rlogxs, 'linear')(xs) * dx
    print("Done")

    pws = p(Rxs)
    PD = P(pws)

    # Normalise PD = 1 across distribution
    PD = PD / np.sum(PD * dx)

    # Calculate the 'width' of the distribution
    PD_integral = np.cumsum(PD * dx)
    low_idx = np.argsort(abs(PD_integral - 0.25))[0]
    high_idx = np.argsort(abs(PD_integral - 0.75))[0]
    width = (xs[high_idx] - xs[low_idx]) / 1.349
    print("Width:", width)

    print("Done. Plotting...")
    idx = np.all([xs < 5e-3], axis=0)
    plt.plot(xs[idx], PD[idx], label="S > %g Jy; $\sigma$ %g Jy" % (minS, width))

plt.legend()
plt.xlabel("D [Jy/beam]")
plt.ylabel("P(D) [(Jy/beam)$^{-1}$]")
plt.savefig("confusion-%d-arcseconds.pdf" % FWHM, transparent=True)
plt.show()
exit()

