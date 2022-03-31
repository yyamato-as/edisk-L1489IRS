"""
Script to plot the continuum image of L1489IRS for EA ALMA workshop 2022 presentation.
Images sent to Ohashi-san.
"""

import astropy.io.fits as fits
import pickle
from astropy.visualization import ImageNormalize, AsinhStretch
import numpy as np
import matplotlib.pyplot as plt
from analysis_utils import plot_2D_map
from astropy.wcs import WCS

prefix = (
    "/raid/work/yamato/eDisk_data/L1489IRS/v0_images/continuum/L1489IRS_SBLB_continuum_"
)
ext = ".pbcor.tt0.fits"

statfile = "/raid/work/yamato/eDisk_data/L1489IRS/analysis_data/L1489IRS_cont_statistics.pkl"
with open(statfile, "rb") as f:
    stat = pickle.load(f)

imaging_params = [
    "robust_2.0",
    "robust_1.0_taper_3000klambda",
]  # images that shows the clear visible ring


for ip in imaging_params:
    # get path and header
    imagepath = prefix + ip + ext
    header = fits.getheader(imagepath)

    # beam and scale setiing
    beam = (
        header["BMAJ"] / np.abs(header["CDELT1"]),
        header["BMIN"] / np.abs(header["CDELT1"]),
        90 + header["BPA"],
    )
    scale = (50 / 140 / 3600 / header["CDELT1"], "50 au")

    # get data
    data = fits.getdata(imagepath)[2100:3900, 2100:3900] * 1e3  # to mJy

    # normalization, contour levels, and color/contour setting
    norm = ImageNormalize(data, vmin=0.0, stretch=AsinhStretch(a=0.03))
    levels = (
        np.array([3, 5, 7, 10, 15, 20, 30, 50, 100, 150, 200])
        * stat.at[ip, "rms [mJy / beam]"]
    )
    imshow_kw = {"norm": norm, "cmap": "inferno"}
    contour_kw = {
        "levels": levels,
        "colors": "white",
        "linewidths": 0.2,
        "linestyles": "dashed",
    }

    # with cbar, label, beam and scalebar
    fig, ax = plt.subplots(subplot_kw={"projection": WCS(header)})
    plot_2D_map(
        data,
        ax=ax,
        beam=beam,
        scale=scale,
        imshow_kw=imshow_kw,
        contour_kw=contour_kw,
        cbar_kw={"label": "mJy / beam"},
    )
    ax.set(xlabel="Right Ascension (ICRS)", ylabel="Declination (ICRS)")
    plt.show()
    fig.savefig(
        "./figure/L1489IRS_continuum_{:s}.pdf".format(ip),
        bbox_inches="tight",
        pad_inches=0.01,
    )

    # no contour
    fig, ax = plt.subplots(subplot_kw={"projection": WCS(header)})
    plot_2D_map(
        data,
        ax=ax,
        beam=beam,
        scale=scale,
        contour=False,
        imshow_kw=imshow_kw,
        cbar_kw={"label": "mJy / beam"},
    )
    ax.set(xlabel="Right Ascension (ICRS)", ylabel="Declination (ICRS)")
    plt.show()
    fig.savefig(
        "./figure/L1489IRS_continuum_{:s}_noContour.pdf".format(ip),
        bbox_inches="tight",
        pad_inches=0.01,
    )

    # without cbar, label, beam and scalebar
    fig, ax = plt.subplots(subplot_kw={"projection": WCS(header)})
    plot_2D_map(data, ax=ax, imshow_kw=imshow_kw, contour_kw=contour_kw, colorbar=False)
    ax.axis("off")
    plt.show()
    fig.savefig(
        "./figure/L1489IRS_continuum_{:s}_noScale.pdf".format(ip),
        bbox_inches="tight",
        pad_inches=0.01,
    )  # with cbar, label, beam and scalebar

    # no contour
    fig, ax = plt.subplots(subplot_kw={"projection": WCS(header)})
    plot_2D_map(data, ax=ax, imshow_kw=imshow_kw, colorbar=False, contour=False)
    ax.axis("off")
    plt.show()
    fig.savefig(
        "./figure/L1489IRS_continuum_{:s}_noContour_noScale.pdf".format(ip),
        bbox_inches="tight",
        pad_inches=0.01,
    )  # with cbar, label, beam and scalebar
