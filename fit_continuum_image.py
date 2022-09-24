from fileinput import filename
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.visualization import AsinhStretch, ImageNormalize
import pickle
from qdisk.fit import imfit_wrapper
import analysis_utils as au
from analysis_utils import ContinuumNormalize
from qdisk.classes import FitsImage
from casatasks import importfits
import os

# image path etc.
source = "L1489IRS"
config = "SBLB"
imtype = "continuum"
robust = 1.0
path = "/works/yamato/eDisk/L1489IRS/try1/trial_image_backup/"
# imagename = au.imageproductpath + au.get_image_basename(
#     source, config, imtype, robust=robust
# )
imagename = path + au.get_image_basename(
    source, config, imtype, robust=robust
)

imagename = "./image_archival/L1489IRS_SBLB_continuum_robust_1.0.pbcor.tt0.fits"


image = FitsImage(imagename)
print(image.beam)


####################### MANUAL SETTINGS ###############################
# nominal source position from data reduction script
dir = "04h04m43.070001s +26d18m56.20011s"
maj = 8
min = 8
pa = 0.0
mask = au.get_casa_ellipse_region(dir, maj, min, pa)

# set the initial estimates file
# peak intensity, peak xpixel, peak ypixel, maj, min, pa
# values from Sai et al. 2020
modelname = "doubleGaussian"
est_str_list = [
    "0.003, 3000, 3000, 0.097arcsec, 0.037arcsec, 49deg\n",
    "0.001, 3000, 3000, 4.1arcsec, 1.2arcsec, 69deg\n",
]  # need \n
# est_str_list = [
#     "0.003, 800, 800, 0.097arcsec, 0.037arcsec, 49deg\n",
#     "0.001, 800, 800, 4.1arcsec, 1.2arcsec, 69deg\n",
# ]  # need \n
########################################################################

# filename_prefix = au.get_image_basename(source, config, imtype, robust=robust).replace(
#     ".fits", ".imfit_{:s}".format(modelname)
# )
filename_prefix = os.path.basename(imagename).replace(
    ".fits", ".imfit_{:s}".format(modelname)
)

estimates_filename = "./fit_archival/" + filename_prefix + ".estimates"
with open(estimates_filename, "w") as f:
    f.writelines(est_str_list)

# set model and residual file
model_filename = "./fit_archival/" + filename_prefix + ".model"
residual_filename = "./fit_archival/" + filename_prefix + ".residual"
logfile = "./fit_archival/" + filename_prefix + ".log"

# measure rms
image.estimate_rms(rmin=8, rmax=10)
print(image.rms)

importfits(fitsimage=imagename, imagename=imagename.replace(".fits", ""), overwrite=True)

# 2 component Gaussian fit
result, fig = imfit_wrapper(
    imagename.replace(".fits", ""),
    region=mask,
    model=model_filename,
    residual=residual_filename,
    estimates=estimates_filename,
    logfile=logfile,
    rms=image.rms,
    plot_result=True,
    plot_kw=dict(
        xlim=(-4, 4),
        ylim=(-4, 4),
        method="pcolorfast",
        # stretch=AsinhStretch(a=0.02)
    ),
)

# fig.savefig(
#     au.figurepath + filename_prefix + ".png",
#     dpi=500,
#     bbox_inches="tight",
#     pad_inches=0.01,
# )

plt.show()

# save the result
# savefile = au.analysisdatapath + filename_prefix + ".pkl"
# with open(savefile, "wb") as f:
#     pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
