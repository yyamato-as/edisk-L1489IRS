import os
import analysis_utils as au
from analysis_utils import ContinuumNormalize
from qdisk.classes import FitsImage
from qdisk.plot import plot_2D_map
from eDisk_source_dict import source_dict
import matplotlib.pyplot as plt
import matplotlib.patches as patches

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Times"],
    "xtick.top": True,
    "ytick.right": True,
    "xtick.direction": "out",
    "ytick.direction": "out"})

############ MANUALLY SET #################
source = "L1489IRS"
config = "SBLB"
imtype = "continuum"
robust = 1.0
center_coord = source_dict[source]["radec"]
distance = source_dict[source]["distance"]
zoomed = False
imagename = au.imageproductpath + au.get_image_basename(
    source, config, imtype, robust=robust, type="pbcor"
)
###########################################

image = FitsImage(imagename)
image.get_directional_coord(center_coord=center_coord)
data = image.data * au.data_scaling_factor["cont"]

norm = ContinuumNormalize(data)
rmax = 2 if zoomed else 5

fig, ax = plt.subplots(figsize=(3, 3))
plot_2D_map(
    data,
    X=image.x,
    Y=image.y,
    ax=ax,
    cmap_method="pcolorfast",
    contour=False,
    beam=image.beam,
    scale=(50 / distance, "50 au"),
    xlim=(-rmax, rmax),
    ylim=(-rmax, rmax),
    cmap_kw=dict(cmap=au.cmap["cont"], norm=norm), 
    cbar_kw=au.cbar_kw["cont"],
    beam_kw=au.beam_kw["cont"]
)

# add pacthes to indicate ring position
# arc = 
# ax.add_patch


ax.set(xlabel=r"$\Delta$R.A. [$^{\prime\prime}$]", ylabel=r"$\Delta$Dec. [$^{\prime\prime}$]", aspect=1./ax.get_data_ratio())

figname = au.figurepath + os.path.basename(imagename).replace(".fits", "_zoomed.png" if zoomed else ".png")
fig.savefig(figname, dpi=500, bbox_inches="tight", pad_inches=0.01)
# plt.show()
