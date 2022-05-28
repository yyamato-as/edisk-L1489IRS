#%%
from fits_xarray import FitsDataArray
import matplotlib.pyplot as plt
from astropy.visualization import AsinhStretch, ImageNormalize
import numpy as np
import astropy.units as u
import astropy.io.fits as fits

#%%
robust = [-2.0, -1.0, -0.5, 0.5, 1.0, 2.0]
taper = [None, '1000klambda', '2000klambda', '3000klambda']
fitsimagepath = "/raid/work/yamato/eDisk_data/L1489IRS/v0_images/"

#%%
# different robust parameter
fig, axes = plt.subplots(1, len(robust), squeeze=False, figsize=(20,4))
for i, r in enumerate(robust):
    fitsimage = fitsimagepath + "L1489IRS_SBLB_continuum_robust_{:s}.pbcor.tt0.fits".format(str(r))

    # import fda
    fda = FitsDataArray(hdu=fits.open(fitsimage)[0])
    data = fda.da

    # unit conversion
    x_unit = data.coords["dRA cos(Dec)"].attrs["unit"]
    y_unit = data.coords["dDec"].attrs["unit"]
    data.coords["dRA cos(Dec)"].values = (data.coords["dRA cos(Dec)"].values * x_unit).to(u.arcsec).value
    data.coords["dDec"].values = (data.coords["dDec"].values * y_unit).to(u.arcsec).value

    # image normalization
    norm = ImageNormalize(data.values, vmin=0.0, stretch=AsinhStretch(a=0.03))

    # plot
    ax = axes[0,i]
    data.plot.pcolormesh(x="dRA cos(Dec)", y="dDec",  ax = ax, rasterized=True, xlim=(1, -1), ylim=(-1, 1), norm=norm, cmap='inferno')
    ax.set(aspect=1./ax.get_data_ratio())
plt.show()

# %%
import plotly.graph_objects as go

#fig, axes = plt.subplots(1, len(robust), squeeze=False, figsize=(20,4))
for i, r in enumerate(robust):
    fitsimage = fitsimagepath + "L1489IRS_SBLB_continuum_robust_{:s}.pbcor.tt0.fits".format(str(r))

    # import fda
    fda = FitsDataArray(hdu=fits.open(fitsimage)[0])
    data = fda.da

    # unit conversion
    x_unit = data.coords["dRA cos(Dec)"].attrs["unit"]
    y_unit = data.coords["dDec"].attrs["unit"]
    data.coords["dRA cos(Dec)"].values = (data.coords["dRA cos(Dec)"].values * x_unit).to(u.arcsec).value
    data.coords["dDec"].values = (data.coords["dDec"].values * y_unit).to(u.arcsec).value
    fig = go.Figure(data=go.Heatmap(z=data.values, x=data.coords['dRA cos(Dec)'].values, y=data.coords['dDec'].values))
    fig.show()
# %%
