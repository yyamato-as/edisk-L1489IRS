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
from mpl_toolkits.axes_grid1 import ImageGrid
from astropy.visualization import AsinhStretch
import astropy.io.fits as fits

source = "L1489IRS"
config = "SBLB"
lines = ["12CO", "13CO", "C18O", "SO"]
robust = {"12CO": 0.5, "13CO": 0.5, "C18O": 1.0, "SO": 1.0}
center_coord = source_dict[source]["radec"]
PA = source_dict[source]["PA"]
incl = source_dict[source]["incl"]
distance = source_dict[source]["distance"]
vsys = source_dict[source]["v_sys"]


dv = 2

for line in lines:
    mapname = au.VADPpath + au.get_image_basename(source, config, line, robust=robust[line]).replace(".fits", "_v0.fits")
    cubename = au.customimagepath + au.get_image_basename(source, config, line, robust=robust[line])

    print("Loading data...")
    cube = FitsImage(cubename)
    cube.shift_phasecenter_toward(center_coord)

    print("Generating threshold mask...")
    cube.estimate_rms(edgenchan=3)
    cube.get_threshold_mask(threshold=4)
    t_mask = np.sum(cube.threshold_mask, axis=0) > 0
    print("Done.")

    v0map = FitsImage(mapname)
    v0map.shift_phasecenter_toward(center_coord)
    v_mask = (v0map.data*1e-3 <= vsys + dv) & (v0map.data*1e-3 >= vsys - dv)
    mask = t_mask | v_mask

    v0map.data[~mask] = np.nan

    savename = mapname.replace("_v0.fits", "_v0_hybridmask.fits")

    fits.writeto(
        savename,
        v0map.data,
        v0map.header,
        overwrite=True,
        output_verify="silentfix",
    )


    


