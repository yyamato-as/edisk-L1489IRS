from qdisk.model import Keplerian_velocity
from qdisk.classes import FitsImage
import eDiskplot as eplot
from qdisk.utils import is_within
from qdisk.plot import ChannelMap, Map
import analysis_utils as au
from eDisk_source_dict import source_dict
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian2DKernel, convolve_fft
from analysis_utils import FWHM_to_sigma
from mcmc_tools import log_prior, emcee_run_wrapper, plot_corner, plot_walker
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.patheffects as pe
from astropy.visualization import ImageNormalize

source = "L1489IRS"
config = "SBLB"
line = "C18O"
robust = 1.0
center_coord = source_dict[source]["radec"]
PA = source_dict[source]["PA"]
incl = source_dict[source]["incl"]
distance = source_dict[source]["distance"]
vsys = source_dict[source]["v_sys"]

fit_vrange = (1.1, 13)
fit_vrange_blue = (1.1, 4.6)
fit_vrange_red = (9.8, 13)
# fit_region_x = (-0.25, 0.25)
# fit_region_y = (-0.2, 0.2)
fit_region = (-1.4, 1.4)

imagename = au.customimagepath + au.get_image_basename(source, config, line, robust=robust)

image = FitsImage(imagename, xlim=fit_region, ylim=fit_region, vlim=None, downsample=False)
image.shift_phasecenter_toward(center_coord)
image.estimate_rms(edgenchan=1)

image.data = image.data[is_within(image.v, fit_vrange_blue) | is_within(image.v, fit_vrange_red)]
image.v = image.v[is_within(image.v, fit_vrange_blue) | is_within(image.v, fit_vrange_red)]

norm = ImageNormalize(image.data*1e3, vmin=0.0)

kernel = Gaussian2DKernel(x_stddev=FWHM_to_sigma(image.bmaj)/image.dpix, y_stddev=FWHM_to_sigma(image.bmin)/image.dpix, theta=np.radians(90+image.bpa))
M8 = np.nanmax(image.data, axis=0)

image.spectrally_collapse(mode="mom0")
M0 = image.collapsed


def Keplerian_model(PA, sigma, Mstar, incl, vsys=vsys):
    r, t = image.get_disk_coord(incl=incl, PA=PA, frame="polar")
    vlos = - Keplerian_velocity(r, t, Mstar=Mstar, distance=distance, incl=incl) + vsys
    # vlos = convolve_fft(vlos, kernel)

    # modelcube = M0 / (np.sqrt(2 * np.pi) * sigma) * np.exp(- (image.v[:, None, None] - vlos[None, :, :]) ** 2 / (2 * sigma ** 2))

    modelcube = np.exp(- (image.v[:, None, None] - vlos[None, :, :]) ** 2 / (2 * sigma ** 2))
    modelcube = np.array([convolve_fft(im, kernel) for im in modelcube]) 
    modelcube *= M8 / np.nanmax(modelcube, axis=0)
    return modelcube

# set up parameters
param = [69, 0.5, 1.6, 70]
bound = [(0, 90), (0.01, 10), (0.1, 10), (0, 90)]

def log_likelihood(param):
    model = Keplerian_model(*param)
    # compute log likelihood
    ll = -0.5 * np.sum((image.data - model) ** 2 / image.rms ** 2)

    return ll

def log_probability(param):
    lp = log_prior(param, bound)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(param)
    return lp + ll

import multiprocessing

print("Going to sample MCMC...")
with multiprocessing.Pool(16) as pool:
    sampler, sample = emcee_run_wrapper(log_probability=log_probability, initial_state=param, nstep=500, pool=pool)


import pickle

with open("./channelmap_fit_C18O.pkl", "wb") as f:
    pickle.dump(sampler, f, protocol=pickle.HIGHEST_PROTOCOL)