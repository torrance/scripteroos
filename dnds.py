from __future__ import division, print_function

import glob

import matplotlib.pyplot as plt
import numpy as np


plt.figure(figsize=(7, 6))

fluxes = np.array([])

for f, label in [
    ('fullrange_galaxies_sf2_150.npy', 'SF2'),
    ('fullrange_galaxies_agn1_150.npy', 'AGN1'),
    ('fullrange_galaxies_agn23_150.npy', 'AGN2'),
    ('fullrange_galaxies_sf1_150.npy', 'SF1'),
]:
    fs = np.load(f)
    fluxes = np.concatenate([fluxes, fs])
    print("Sources:",  len(fs))

# SKADS is at 151 MHz - use spectral index of -0.8 to shift this to 154 MHz
fluxes = fluxes + 0.8 * np.log10(151e6) - 0.8 * np.log10(154e6)

fluxes = fluxes + np.log10(1.15) # 1.2 flux multiplier from Frazen
# Filter out sources beyond 1e-8 Jy < S < 1 Jy
idx = np.all([fluxes > -7.5, fluxes < 2], axis=0)
fluxes = fluxes[idx].copy()

bins = np.linspace(-7.5, 2, 1000)
dxs = 10**bins[1:] - 10**bins[:-1]
mids = (bins[:-1] + bins[1:]) / 2
hist, _ = np.histogram(fluxes, bins)
hist = hist / (400 * (np.pi**2 / 180**2))  # per 400 degrees squared  -> per steradian
hist = hist / dxs

# Remove bins with zeros
mids = mids[hist != 0]
dxs = dxs[hist != 0]
hist = hist[hist != 0]
loghist = np.log10(hist) + 2.5 * mids

# Plot raw SKADS data
plt.scatter(mids, loghist, s=2, alpha=0.2, label="SKADS catalogue @ 154 MHz x 1.15", c='black')

# Plot Franzen (2019) model
franzen = [0.006, 0.0351, -0.0404, -0.388, 0.307, 3.52]
idx = mids > -3
ys = np.zeros_like(mids[idx])
for i, const in enumerate(franzen):
    power = len(franzen) - i - 1
    ys += const * mids[idx]**power

plt.plot(mids[idx], ys, label="Franzen (2019)", color='red')

# Plot SKADS fit beneath 10^-3 Jy
idx = mids < -3
fit = np.polyfit(mids[idx], loghist[idx], 5)
ys = np.zeros_like(mids[idx])
for i, const in enumerate(fit):
    power = len(fit) - i - 1
    ys += const * (mids[idx])**power

print(np.flip(fit))

plt.plot(mids[idx], ys, label="SKADS fit @154 MHz x 1.17 ($10^{-7.5} < S < 10^{-3}$)", color='dodgerblue')

# # Franzen's SKADS fit
# fit = np.flip([-12.4489, -14.2704, -5.59825, -1.08487, -0.108262, -0.00433429])
# idx = mids < -3
# fit = np.polyfit(mids[idx], loghist[idx], 5)
# ys = np.zeros_like(mids[idx])
# for i, const in enumerate(fit):
#     power = len(fit) - i - 1
#     ys += const * (mids[idx])**power

# plt.plot(mids[idx], ys, label="Franzen's SKADS fit @154 MHz", color='green')

plt.xlim([-7.5, 2])
plt.legend()
plt.xlabel("154 MHz flux density (Jy)")
plt.ylabel("dN/dS S$^2.5$ (Jy$^{1.5}$ sr$^{-1}$)")
plt.savefig('dnds.pdf', transparent=True)
plt.show()



