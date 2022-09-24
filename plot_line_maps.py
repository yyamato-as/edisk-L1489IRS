from qdisk.plot import ChannelMap, Map
from eDisk_source_dict import source_dict
import matplotlib.pyplot as plt
from matplotlib import ticker
from astropy.visualization import AsinhStretch, SinhStretch
import matplotlib.patheffects as pe
from matplotlib import patches
import eDiskplot as eplot
from qdisk.classes import FitsImage
import numpy as np
import analysis_utils as au

source = "L1489IRS"
config = "SBLB"
# lines = ["13CO", "C18O", "SO"]
line = "13CO"
robust = 1.0
dv = 0.2
center_coord = source_dict[source]["radec"]
PA = source_dict[source]["PA"]
incl = source_dict[source]["incl"]
distance = source_dict[source]["distance"]
vsys = source_dict[source]["vsys"]
vrange = 6
nu0 = source_dict[source]["cont_nu0"]

def arcsec2au(r):
    return r * distance

def au2arcsec(r):
    return r / distance


fig, axes = plt.subplots(
    1, 3, figsize=(9, 3.5), sharex=True, sharey=True
)  # , constrained_layout=True)
rmax = 3.5
scale = 100  # in au

# moment0
ax = axes[0]
mom0name = au.VADPpath + au.get_image_basename(
    source=source, baseline=config, line=line, robust=robust, dv=dv
).replace(".fits", "_M0.fits")
mom0map = Map(
    mom0name, ax=ax, center_coord=center_coord, xlim=(-rmax, rmax), ylim=(-rmax, rmax)
)
mom0map.estimate_rms(edgenchan=3)
mom0map.plot_colormap(cmap=eplot.cmap["M0"], vmin=0.0, stretch=AsinhStretch(a=0.4))
# mom0map.overlay_contour(contimagename, levels=np.array([8])*contrms, linewidth=0.6)
mom0map.add_colorbar(
    position="top", label="$I$ [mJy beam$^{-1}$ km s$^{-1}$]", rotation=0, labelpad=5
)
mom0map.add_beam()
mom0map.add_scalebar(scale=scale / distance, text=str(scale) + " au")
mom0map.colorbar.ax.minorticks_on()
mom0map.set_ticker(majornticks=5, minor=True, minornticks=4)
mom0map.set_labels(
    xlabel="$\Delta$R.A. [$^{\prime\prime}$]", ylabel="$\Delta$Dec. [$^{\prime\prime}$]"
)
# mom0map.ax.scatter(0, 0, marker="+", s=50, color="white", linewidth=1)
x, y = mom0map.get_disk_coord(PA=PA, incl=incl, frame="cartesian")
mom0map.overlay_contour(y, x=mom0map.x, y=mom0map.y, levels=[0.0], linewidth=0.6, linestyle="dashed", color="grey", alpha=0.5)
mom0map.overlay_contour(x, x=mom0map.x, y=mom0map.y, levels=[0.0], linewidth=0.6, linestyle="dashed", color="grey", alpha=0.5)
mom0map.ax.annotate(
    text="(a)",
    xy=(0.1, 0.9),
    xycoords="axes fraction",
    color="black",
    ha="center",
    va="center",
    path_effects=[pe.Stroke(linewidth=3, foreground="white"), pe.Normal()],
)

# v0map
ax = axes[1]
vrange = 4
cube = FitsImage(imagename, xlim=(-rmax, rmax), ylim=(-rmax, rmax))
cube.shift_phasecenter_toward(center_coord)
print("Generating threshold mask...")
cube.estimate_rms(edgenchan=3)
cube.get_threshold_mask(threshold=4)
mask = np.sum(cube.threshold_mask, axis=0) > 0
print("Done.")
v0name = au.VADPpath + au.get_image_basename(
    source=source, baseline=config, line=line, robust=robust, dv=dv
).replace(".fits", "_v0.fits")
v0map = Map(
    v0name,
    ax=ax,
    center_coord=center_coord,
    xlim=(-rmax, rmax),
    ylim=(-rmax, rmax),
    data_scaling_factor=1e-3,
    invert_xaxis=False,
)
v0map.mask(vmin=vsys - 2, vmax=vsys + 2, user_mask=mask, combine="or")
v0map.plot_colormap(cmap=eplot.cmap["M1"], vmin=vsys - vrange, vmax=vsys + vrange)
v0map.add_colorbar(position="top", label="$v_0$ [km s$^{-1}$]", rotation=0, labelpad=5)
v0map.add_beam(color="black")
v0map.add_scalebar(scale=scale / distance, text=str(scale) + " au", color="black")
# v0map.ax.scatter(0, 0, marker="+", s=50, color="black", linewidth=0.5)
x, y = v0map.get_disk_coord(PA=PA, incl=incl, frame="cartesian")
v0map.overlay_contour(y, x=mom0map.x, y=mom0map.y, levels=[0.0], linewidth=0.6, linestyle="dashed", color="0.2", alpha=0.5)
v0map.overlay_contour(x, x=mom0map.x, y=mom0map.y, levels=[0.0], linewidth=0.6, linestyle="dashed", color="0.2", alpha=0.5)
v0map.colorbar.ax.minorticks_on()
v0map.ax.annotate(
    text="(b)",
    xy=(0.1, 0.9),
    xycoords="axes fraction",
    color="black",
    ha="center",
    va="center",
    path_effects=[pe.Stroke(linewidth=3, foreground="white"), pe.Normal()],
)

# M8map
ax = axes[2]
mom8name = au.VADPpath + au.get_image_basename(
    source=source, baseline=config, line=line, robust=robust, dv=dv
).replace(".fits", "_Fnu.fits")
peakmap = Map(
    mom8name,
    ax=ax,
    center_coord=center_coord,
    xlim=(-rmax, rmax),
    ylim=(-rmax, rmax),
    invert_xaxis=False,
)
peakmap.convert_unit()
peakmap.plot_colormap(cmap=eplot.cmap["M8"], vmin=10.0, stretch=SinhStretch())
peakmap.add_colorbar(position="top", label="$T_\mathrm{B}$ [K]", rotation=0, labelpad=5)
peakmap.add_beam()
peakmap.add_scalebar(scale=scale / distance, text=str(scale) + " au")
# peakmap.ax.scatter(0, 0, marker="+", s=50, color="white", linewidth=0.5)
x, y = peakmap.get_disk_coord(PA=PA, incl=incl, frame="cartesian")
peakmap.overlay_contour(y, x=mom0map.x, y=mom0map.y, levels=[0.0], linewidth=0.6, linestyle="dashed", color="grey", alpha=0.5)
peakmap.overlay_contour(x, x=mom0map.x, y=mom0map.y, levels=[0.0], linewidth=0.6, linestyle="dashed", color="grey", alpha=0.5)
peakmap.colorbar.ax.minorticks_on()
peakmap.ax.annotate(
    text="(c)",
    xy=(0.1, 0.9),
    xycoords="axes fraction",
    color="black",
    ha="center",
    va="center",
    path_effects=[pe.Stroke(linewidth=3, foreground="white"), pe.Normal()],
)
