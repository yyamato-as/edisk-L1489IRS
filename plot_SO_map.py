import eDiskplot as eplot
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import AsinhStretch
from eDisk_source_dict import source_dict
import analysis_utils as au
from qdisk.classes import FitsImage
from qdisk.plot import Map
from matplotlib.backends.backend_pdf import PdfPages

source = "L1489IRS"
config = "SBLB"
line = "SO"
robust = 1.0
vsys = source_dict[source]["v_sys"]
center_coord = source_dict[source]["radec"]
distance = source_dict[source]["distance"]
vrange = 5
linelabel = "SO ($J_N=$ $6_5$--$5_4$)"
imagename = au.customimagepath + au.get_image_basename(source, config, line, robust=robust)
contname = au.customimagepath + au.get_image_basename(source, config, "continuum", robust=1.0)

cont = FitsImage(contname, xlim=(-8, 8), ylim=(-8, 8))
cont.estimate_rms(rmin=8)

# integrated intensity map
fitsname = au.VADPpath + au.get_image_basename(source, config, line=line, robust=robust).replace(".fits", "_2sigma_M0.fits")
rmax = [20, 8, 2.5]

with PdfPages(au.figurepath + "integmap_{:s}_multiscale.pdf".format(line)) as pdf:
    # integrated intensity map
    for r in rmax:
        fig, ax = plt.subplots(1, 1, figsize=(6, 5), dpi=300)
        mom0map = Map(fitsname=fitsname, ax=ax, center_coord=center_coord, xlim=(-r, r), ylim=(-r, r))
        mom0map.plot_colormap(cmap=eplot.cmap["M0"], vmin=0.0, stretch=AsinhStretch(a=0.05))
        mom0map.add_beam(color="white")
        mom0map.add_scalebar(scale=200/distance, text="200 au", color="white")
        mom0map.add_colorbar(label="mJy beam$^{-1}$ km s$^{-1}$", labelpad=20)
        mom0map.overlay_contour(contname, levels=np.array([6, 12, 50, 200])*cont.rms, linewidth=1.0)
        mom0map.set_labels(xlabel=eplot.RAlabel, ylabel=eplot.Declabel)
        mom0map.set_ticker(minor=True, majornticks=7, minornticks=5)

        ax.set_title(linelabel)

        pdf.savefig(fig)


# rotation map

### plot
fitsname = au.VADPpath + au.get_image_basename(source, config, line=line, robust=robust).replace(".fits", "_v0.fits")
cubename = au.customimagepath + au.get_image_basename(source, config, line, robust=robust)
rmax = [20, 8, 2.5]

with PdfPages(au.figurepath + "rotmap_{:s}_mutltiscale.pdf".format(line)) as pdf:
    for r in rmax:

        print("Loading data...")
        cube = FitsImage(cubename, xlim=(-r, r), ylim=(-r, r))
        cube.shift_phasecenter_toward(center_coord)

        print("Generating threshold mask...")
        cube.estimate_rms(edgenchan=3)
        cube.get_threshold_mask(threshold=4)
        mask = np.sum(cube.threshold_mask, axis=0) > 0
        print("Done.")

        fig, ax = plt.subplots(1, 1, figsize=(6, 5), dpi=300)
        rotmap = Map(fitsname=fitsname, ax=ax, data_scaling_factor=1e-3, center_coord=center_coord, xlim=(-r, r), ylim=(-r, r))
        rotmap.mask(vmin=vsys-2, vmax=vsys+2, user_mask=mask, combine="or")
        rotmap.plot_colormap(cmap=eplot.cmap["M1"], vmin=vsys-vrange, vmax=vsys+vrange)
        rotmap.add_beam(color="black")
        rotmap.add_scalebar(scale=200/distance, text="200 au", color="black")
        rotmap.add_colorbar(label="km s$^{-1}$", labelpad=20)
        rotmap.overlay_contour(contname, levels=np.array([6, 12, 50, 200])*cont.rms, linewidth=1.0)
        rotmap.set_labels(xlabel=eplot.RAlabel, ylabel=eplot.Declabel)
        rotmap.set_ticker(minor=True, majornticks=7, minornticks=5)

        ax.set_title(linelabel)

        pdf.savefig(fig)

### peak intensity map
fitsname = au.VADPpath + au.get_image_basename(source, config, line=line, robust=robust).replace(".fits", "_Fnu.fits")
rmax = [20, 8, 2.5]

with PdfPages(au.figurepath + "peakmap_{:s}_mutltiscale.pdf".format(line)) as pdf:
    for r in rmax:

        fig, ax = plt.subplots(1, 1, figsize=(6, 5), dpi=300)
        peakmap = Map(fitsname=fitsname, ax=ax, data_scaling_factor=1e3, center_coord=center_coord, xlim=(-r, r), ylim=(-r, r))
        peakmap.convert_unit()
        peakmap.plot_colormap(cmap=eplot.cmap["M8"], vmin=0.0)
        peakmap.add_beam(color="white")
        peakmap.add_scalebar(scale=200/distance, text="200 au", color="white")
        peakmap.add_colorbar(label="K", labelpad=20)
        peakmap.overlay_contour(contname, levels=np.array([6, 12, 50, 200])*cont.rms, linewidth=1.0)
        peakmap.set_labels(xlabel=eplot.RAlabel, ylabel=eplot.Declabel)
        peakmap.set_ticker(minor=True, majornticks=7, minornticks=5)

        ax.set_title(linelabel)

        pdf.savefig(fig)