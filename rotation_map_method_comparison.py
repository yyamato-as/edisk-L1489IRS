import analysis_utils as au
import eDiskplot as eplot
from eDisk_source_dict import source_dict
from qdisk.classes import FitsImage
from qdisk.utils import plot_2D_map
import matplotlib.pyplot as plt
from qdisk.product import calculate_moment

source = "L1489IRS"
config = "SBLB"
line = "C18O"
robust = 1.0
center_coord = source_dict[source]["radec"]
vsys = source_dict[source]["v_sys"]
vrange = 5

plt.rcParams.update({
    "xtick.top": True,
    "ytick.right": True,
    "xtick.direction": "out",
    "ytick.direction": "out"})

imagename = au.customimagepath + au.get_image_basename(
    source, config, line, robust=robust
)

# traditional moment 1 with different sigma threshold
# fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(9, 3))

# sigma = [3, 4, 5]

# for i, s in enumerate(sigma):
#     calculate_moment(
#         imagename=imagename,
#         moments=[1],
#         threshold=s,
#         save=True,
#         savefilename=au.VADPpath
#         + au.get_image_basename(source, config, line, robust=robust),
#         nchunks=4,
#     )

#     ax = axes[i]
#     v0image = au.VADPpath + au.get_image_basename(
#         source, config, line, robust=robust
#     ).replace(".fits", "_M1.fits")
#     image = FitsImage(v0image)
#     image.get_directional_coord(center_coord=center_coord)
#     image.data *= 1e-3  # in km/s

#     # colorbar = True if i == 2 else False

#     plot_2D_map(
#         image.data,
#         X=image.x,
#         Y=image.y,
#         ax=ax,
#         cmap_method="pcolorfast",
#         contour=False,
#         cmap_kw=dict(cmap=eplot.cmap["M1"], vmin=vsys - vrange, vmax=vsys + vrange),
#         xlim=(-8, 8),
#         ylim=(-8, 8),
#     )

#     ax.set(
#         aspect=1.0 / ax.get_data_ratio(),
#         xlim=(8, -8),
#         ylim=(-8, 8),
#         title="threshold = {}$\sigma$".format(s),
#     )

#     if i == 0:
#         ax.set(xlabel="$\Delta$R.A. [arcsec]", ylabel="$\Delta$Dec. [arcsec]")

# fig.suptitle("moment 1")
# fig.savefig(
#     "./figure/"
#     + au.get_image_basename(source, config, line, robust=robust).replace(
#         ".fits", "_moment1_threshold_comparison.png"
#     ),
#     dpi=800,
#     bbox_inches="tight",
#     pad_inches=0.01,
# )
# plt.show()

# method comparison; traditional moment 1, moment 9, quadratic method
fig, axes = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(9, 5))
moment = [1, 9, "q"]
prefix = {1: "M1", 9: "M9", "q": "v0"}
cmap_kw = {
    "val": dict(cmap=eplot.cmap["M1"], vmin=vsys - vrange, vmax=vsys + vrange),
    "err": dict(cmap="viridis", vmin=0.0, vmax=0.5),
}
method = {1: "moment 1 (4$\sigma$)", 9: "moment 9", "q": "quadratic"}


for i, mom in enumerate(moment):
    calculate_moment(
        imagename=imagename,
        moments=[mom],
        threshold=4 if mom == 1 else None,
        # threshold=4,
        save=True,
        savefilename=au.VADPpath
        + au.get_image_basename(source, config, line, robust=robust),
        nchunks=4 if mom == 1 else None,
    )

    for j, pre in enumerate([prefix[mom], "d"+prefix[mom]]):
        ax = axes[j, i]
        v0image = au.VADPpath + au.get_image_basename(
            source, config, line, robust=robust
        ).replace(".fits", "_{:s}.fits".format(pre))
        image = FitsImage(v0image)
        image.get_directional_coord(center_coord=center_coord)
        image.data *= 1e-3  # in km/s

        # colorbar = True if i == 2 else False

        plot_2D_map(
            image.data,
            X=image.x,
            Y=image.y,
            ax=ax,
            cmap_method="pcolorfast",
            contour=False,
            cmap_kw=cmap_kw["err"] if "d" in pre else cmap_kw["val"],
            xlim=(-8, 8),
            ylim=(-8, 8),
        )

        ax.set(
            aspect=1.0 / ax.get_data_ratio(),
            xlim=(8, -8),
            ylim=(-8, 8),
        )

        if j == 0:
            ax.set_title(method[mom])

axes[1, 0].set(xlabel="$\Delta$R.A. [arcsec]", ylabel="$\Delta$Dec. [arcsec]")

fig.savefig(
    "./figure/"
    + au.get_image_basename(source, config, line, robust=robust).replace(
        ".fits", "_rotationmap_method_comparison_all4sigmaThresh.png"
    ),
    dpi=800,
    bbox_inches="tight",
    pad_inches=0.01,
)
