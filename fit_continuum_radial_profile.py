from analysis_utils import (
    Gaussian1d,
    GaussianRing1d,
    PowerLaw,
    FourthPowerGaussianRing1d,
)
from scipy.optimize import curve_fit
import analysis_utils as au
import numpy as np
import matplotlib.pyplot as plt
import emcee
from qdisk.utils import is_within


def ring_model(r, I_c, sigma_c, r0_r, I_r, sigma_r, I_b, sigma_b):
    return (
        Gaussian1d(r, I_c, sigma_c)
        + GaussianRing1d(r, r0_r, I_r, sigma_r)
        + Gaussian1d(r, I_b, sigma_b)
    )


# load radial profile
profile = (
    au.VADPpath
    + "L1489IRS_SBLB_continuum_robust_1.0.image.tt0_radialProfileWedge90deg.txt"
)
r, I, dI = np.loadtxt(profile, unpack=True)

I *= 1e3
dI *= 1e3

fit_radii = (0, 2.8)
I = I[is_within(r, fit_radii)]
dI = dI[is_within(r, fit_radii)]
r = r[is_within(r, fit_radii)]


# ring model fit
# p0 = [-1.4, 0.05, 0.3, -1, 0.23, 0, 2]
# bounds = ([-3, 0.01, 0.05, -3, 0.05, -3, 0.5], [2, 0.1, 0.7, 2, 0.5, 2, 3])
# model = ring_model
# popt, pcov = curve_fit(model, r, I, p0=p0, sigma=dI, absolute_sigma=True, bounds=bounds)
# print(popt)

# fig, ax = plt.subplots()

# # data
# ax.plot(r, I, color="black")
# ax.fill_between(r, I - dI, I + dI, color="black", alpha=0.25, edgecolor=None)

# # best-fit model
# ax.plot(r, model(r, *popt), color="blue")
# # ax.plot(r, model(r, *p0))

# ax.set(xlim=(0, 3))

# plt.show()


# gap fit

# get ring peak as a baseline
# radii_region = (0.2, 0.8)
# baseline = np.max(I[is_within(r, radii_region)])

# print(baseline)


# def gap_model(r, r0, I, sigma):
#     return baseline - GaussianRing1d(r, r0, I, sigma)

def gap_model(r, I_c, sigma_c, I_b, p, r0, I_r, sigma_r):
    return Gaussian1d(r, I_c, sigma_c) + PowerLaw(r, I_b, p) - GaussianRing1d(r, r0, I_r, sigma_r)


# rmax = r[I == baseline]
# print(rmax)
# fit_radii = (0.1, rmax)
# p0 = [0.2, -1, 0.2]
# bounds = ([0.1, -3, 0.01], [rmax, 0, 0.3])


p0 = [-1, 0.05, -1, -1., 0.2, -1, 0.2]
bounds = ([-3, 0.01, -3, -2., 0.05, -3, 0.05], [1, 0.2, 1, -0.1, 0.5, 1, 0.5])
model = gap_model
popt, pcov = curve_fit(
    model,
    r,
    I,
    p0=p0,
    sigma=dI,
    absolute_sigma=True,
    bounds=bounds,
)
perr = np.sqrt(np.diag(pcov))
print(popt, perr)

fig, ax = plt.subplots()

# data
ax.plot(r, I, color="black")
ax.fill_between(r, I - dI, I + dI, color="black", alpha=0.25, edgecolor=None)

# best-fit model
ax.plot(r, model(r, *popt), color="tab:blue", lw=3)
# ax.plot(r, model(r, *popt) - Gaussian1d(r, *popt[:2]) - PowerLaw(r, *popt[2:4]))
# ax.plot(r, 0.3-GaussianRing1d(r, *popt[4:]))

# ax.set(yscale="log")

plt.show()
