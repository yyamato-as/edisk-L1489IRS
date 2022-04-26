from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.io.fits import header
from astropy.visualization.stretch import AsinhStretch
import numpy as np
import astropy.units as u
from astropy.time import Time
from radio_beam import Beam
from astropy.wcs import WCS
import itertools
import astropy.constants as ac
import xarray as xr

CTYPE_CELESTIAL = ["RA", "DEC"]
CTYPE_SPECTRAL = ["FREQ", "VELO-LSR"]
CTYPE_STOKES = ["STOKES"]


class FitsDataArray:
    def __init__(self, hdu, relative=True, ref_coord=None, scale_data=False):

        # get header
        self.header = hdu.header

        # get data
        self.data = hdu.data

        self.relative = relative
        self.ref_coord = ref_coord

        # data scaling based on 'BZERO' and 'BSCALE'
        if scale_data:
            self.data = self.header["BZERO"] + self.header["BSCALE"] * self.data

        # start reading the header
        # check standard (mandatory) keyword
        self._check_standard_keywords()

        # check data shape
        self._check_data_shape_consistency()

        # obs date
        self.obsdate = Time(
            self.header["DATE-OBS"], format="isot", scale="utc"
        )  # assume that the time format is isot and the observation date is no earlier than 1972-01-01

        # object name (todo: may not be needed?)
        self.object = self.header["OBJECT"]

        # beam information
        self.beam = Beam.from_fits_header(
            self.header
        )  # assume the unit is deg fro bmaj, bmin, and pa

        # WCS coordinate
        # WCS class do the transformation between the pixel coordinate and world coordinate. see https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf (https://staff.fukuoka-edu.ac.jp/kanamitu/fits/jdoc/fits_t60.pdf for Japanese version)
        self.wcs = WCS(self.header)
        # self._drop_degenerate_axes()
        # self.pixel_coords = np.meshgrid([np.arange(n) for n in self.data.shape[::-1]]) # inverse order as data.shape
        # self.world_coords = self.wcs.pixel_to_world(*self.pixel_coords) # transform from pixel coordinate to world coordinate

        # coordinates
        self._get_intermediate_world_coord()
        self._reorder_axes_data()
        self._calc_direction_axes()
        if self.has_spectral:
            self._calc_spectral_axes()
        if self.has_stokes:
            self._calc_stokes_axes()

    @property
    def da(self):
        self._fetch_data_unit()
        self._fetch_axes_unit()
        dims = tuple(self.axes_order)
        coords = {
            "Right Ascension": (
                tuple(self.celestial_ctype),
                self.alpha,
                {"unit": self.axesunit[0]},
            ),
            "Declination": (
                tuple(self.celestial_ctype),
                self.delta,
                {"unit": self.axesunit[1]},
            ),
        }
        if self.relative:
            coords["dRA cos(Dec)"] = (
                tuple(self.celestial_ctype),
                self.dalpha,
                {"unit": self.axesunit[0]},
            )
            coords["dDec"] = (
                tuple(self.celestial_ctype),
                self.ddelta,
                {"unit": self.axesunit[1]},
            )
        if self.has_spectral:
            coords["Spectral Axis"] = (
                tuple(self.spectral_ctype),
                self.specax,
                {"unit": self.axesunit[2]},
            )
        if self.has_stokes:
            raise NotImplementedError
        for i, d in enumerate(dims):
            coords[d] = (d, self.pixel_coord[i])

        attrs = {
            "unit": self.dataunit,
            "obs-date": self.obsdate,
            "object": self.object,
            "beam": self.beam,
            "radec-frame": (self.header["RADESYS"], self.header["EQUINOX"])
            if "EQUINOX" in self.header
            else self.header["RADESYS"],
            "center_coord": self.ref_coord,
            "wcs": self.wcs,
        }
        if self.has_spectral:
            self._fetch_restfreq()
            attrs["restfreq"] = self.nu0
            attrs["spec-frame"] = self.header["SPECSYS"]
            attrs["velocity-convention"] = self.header["VELREF"]

        return xr.DataArray(self.data, coords=coords, dims=dims, attrs=attrs)

    @property
    def has_spectral(self):
        return self.header["NAXIS"] > 2 and self.intermediate_world_coord[2].size != 1

    @property
    def has_stokes(self):
        return self.header["NAXIS"] > 3 and self.intermediate_world_coord[3].size != 1

    def _fetch_data_unit(self):
        """A function to get the appropriate unit of the data from the header. This is needed due to the discouraged usage of "BUNIT" keyword in the header. Even though the use of multiple "/" (slashes) is prohibited in the fits standard (see Section 4.3 of https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf), sometimes radio image fits files use (e.g., "Jy/beam km/s" for velocity-integrated intensity maps)"""
        if self.header["BUNIT"].count("/") > 1:
            unit_list = [u.Unit(s) for s in self.header["BUNIT"].split(" ")]
            self.dataunit = np.prod(unit_list)
        else:
            self.dataunit = u.Unit(self.header["BUNIT"], format="fits")

    def _fetch_axes_unit(self):
        self.axesunit = []
        for i in self.order_map:
            self.axesunit.append(u.Unit(self.header["CUNIT{:d}".format(i + 1)]))

    def _fetch_restfreq(self):
        if "RESTFRQ" in self.header:
            self.nu0 = self.header["RESTFRQ"]
        elif "RESTFREQ" in self.header:
            self.nu0 = self.header["RESTFREQ"]
        else:
            raise KeyError("Any rest frequency keywords not found.")

    def _check_standard_keywords(self):
        """A function to check the standard keyword of a fits header. If any standard keywords are lacked or invalid, warings will be appeared."""
        # primary HDU or extension
        if "SIMPLE" not in self.header.keys() and "XTENSION" not in self.header.keys():
            print(
                'Warning: "SIMPLE" nor "XTENSION" keyword found. Continuing anyway though this HDU may not satisfy the FITS Standard.'
            )
        elif "SIMPLE" in self.header.keys():
            if self.header["SIMPLE"] != "T":
                print(
                    'Warning: "SIMPLE" keyword is invalid. Continuing anyway though this HDU may not satisfy the FITS Standard.'
                )
            standard_keyword_list = ["BITPIX", "END"]
            self._keyword_absence_warning(keyword_list=standard_keyword_list)
        elif "XTENSION" in self.header.keys():
            standard_keyword_list = ["BITPIX", "PCOUNT", "GCOUNT", "END"]
            self._keyword_absence_warning(keyword_list=standard_keyword_list)

    def _keyword_absence_warning(self, keyword_list):
        """A function to generate the warning about the absence of the given keyword in the fits header.

        Args:
            keyword_list (list): The list of keywords one wants to check.
        """
        for keyword in keyword_list:
            if keyword in self.header.keys():
                print(
                    'Warning: "{:s}" keyword not found. Continuing anyway though this HDU may not satisfy the FITS Standard.'.format(
                        keyword
                    )
                )

    def _check_data_shape_consistency(self):
        # 'NAXIS' to specify the number of axis
        if "NAXIS" not in self.header.keys():
            raise KeyError(
                '"NAXIS" keyword not found. Unable to calculate coordinates.'
            )

        # check consistency between data and header
        if self.header["NAXIS"] != self.data.ndim:
            raise ValueError("Mismatch in the data dimension and the header.")
        if (
            tuple(
                [
                    self.header["NAXIS{:d}".format(i)]
                    for i in range(1, self.header["NAXIS"] + 1)
                ]
            )
            != self.data.shape
        ):
            raise ValueError("Mismatch in the data shape and the header.")

    # def _drop_degenerate_axes(self):
    #     self.data = self.data.squeeze()
    #     for i in range(1, self.header["NAXIS"]+1):
    #         if self.header['NAXIS{:d}'.format(i)] == 1:
    #             self.wcs.dropaxis(i)
    #     if self.data.shape != self.wcs.array_shape:
    #         raise ValueError("Mismatch in the data shape and the WCS shape after the dropping the degenerate axes.")

    # def _get_ctype(self):
    #     self.ctype = [self.header['CTYPE{:d}'.format(i)] for i in range(1, self.data.ndim+1)]
    #     self.celestial_ctype = []
    #     self.spectral_ctype = []
    #     self.stokes_ctype = []
    #     for c in self.ctype:
    #         c_name = c.split('-')[0]
    #         if c_name in CTYPE_CELESTIAL:
    #             self.celestial_ctype.append(c)
    #         elif c_name in CTYPE_SPECTRAL:
    #             self.spectral_ctype.append(c)
    #         elif c_name in CTYPE_STOKES:
    #             self.stokes_ctype.append(c)
    #         else:
    #             raise ValueError("Coordinate type {:s} not recognized.".format(c))

    @property
    def ctype(self):
        return [
            self.header["CTYPE{:d}".format(i)] for i in range(1, self.data.ndim + 1)
        ]

    @property
    def celestial_ctype(self):
        ctype = []
        for c in self.ctype:
            c_name = c.split("-")[0]
            if c_name.upper() in CTYPE_CELESTIAL:
                ctype.append(c)
        if len(ctype) > 2:
            raise ValueError("More than 2 axes for celestial coord.")
        return ctype

    @property
    def spectral_ctype(self):
        ctype = []
        for c in self.ctype:
            c_name = c.split("-")[0]
            if c_name.upper() in CTYPE_SPECTRAL:
                ctype.append(c)
        if len(ctype) > 1:
            raise ValueError("More than 1 axes for spectral coord.")
        return ctype

    @property
    def stokes_ctype(self):
        ctype = []
        for c in self.ctype:
            c_name = c.split("-")[0]
            if c_name.upper() in CTYPE_STOKES:
                ctype.append(c)
            if len(ctype) > 1:
                raise ValueError("More than 1 axes for stokes coord.")
        return ctype

    @property
    def celestial_axes_physical_types(self):
        type = [c.split("-")[0].upper() for c in self.celestial_ctype]
        if type != ["RA", "DEC"]:
            raise NotImplementedError("Coordinate type {:s} not supported.")
        return type

    @property
    def celestial_axes_projection(self):
        projection = np.unique([c.split("-")[-1] for c in self.celestial_ctype])
        if len(projection) != 1:
            raise ValueError("Different projection methods for longitude and latitude.")
        return projection[0]

    @property
    def pixel_coord(self):
        return [np.arange(n) for n in self.data.shape]

    def _get_pixcoord_conv_mat(self):
        self.conv_mat = np.diag([1.0] * self.data.ndim)
        try:
            # `CDi_j` keyword?
            for i, j in itertools.product(range(self.data.ndim), range(self.data.ndim)):
                self.conv_mat[i, j] = self.header["DC{:d}_{:d}".format(i + 1, j + 1)]
        except KeyError:
            try:
                # `PCi_j` keyword?
                for i, j in itertools.product(
                    range(self.data.ndim), range(self.data.ndim)
                ):
                    self.conv_mat[i, j] = (
                        self.header["CDELT{:d}".format(i + 1)]
                        * self.header["PC{:d}_{:d}".format(i + 1, j + 1)]
                    )
            except KeyError:
                # fall back to default: no transformation, just multiply `CDELTi`
                for i, j in itertools.product(
                    range(self.data.ndim), range(self.data.ndim)
                ):
                    self.conv_mat[i, j] *= self.header["CDELT{:d}".format(i + 1)]

    def _get_intermediate_world_coord(self):
        self._get_pixcoord_conv_mat()
        if np.trace(self.conv_mat) != np.sum(
            self.conv_mat
        ):  # trick to check if a mat is diagonal
            raise NotImplementedError(
                "Non-diagonal conversion matrix from pix coord to intermediate world coord not supported."
            )
        self.intermediate_world_coord = []
        for i in range(self.data.ndim):
            q = self.conv_mat[i, i] * (
                self.pixel_coord[i] + 1 - self.header["CRPIX{:d}".format(i + 1)]
            )  # important to +1 for FITS convention
            self.intermediate_world_coord.append(q)

    def _reorder_axes_data(self):
        self.axes_order = self.celestial_ctype + self.spectral_ctype + self.stokes_ctype
        self.order_map = [self.ctype.index(c) for c in self.axes_order]
        self.intermediate_world_coord = [
            self.intermediate_world_coord[i] for i in self.order_map
        ]
        self.data = self.data.transpose(
            self.order_map
        )  # data has intrinsically (nu, y, x) order

    def _calc_direction_axes(self):
        if self.celestial_axes_projection != "SIN":
            raise NotImplementedError(
                "Projection `{:s}` not supported.".format(
                    self.celestial_axes_projection
                )
            )

        xx, yy = np.meshgrid(*self.intermediate_world_coord[:2])
        # assume no non-linear correction (PVi_j keywords are 0.0)

        # set the native lon, lat of pole
        phi_P = self.header["LONPOLE"] if "LONPOLE" in self.header else 180.0
        theta_P = self.header["LATPOLE"] if "LATPOLE" in self.header else None

        # set the native celestial lon, lat of the pole
        if isinstance(self.ref_coord, SkyCoord):
            alpha_P = self.ref_coord.ra.deg
            delta_P = self.ref_coord.dec.deg
        elif isinstance(self.ref_coord, (tuple, list)):
            c = SkyCoord(*self.reference_coord[:-1], frame=self.ref_coord[-1])
            alpha_P = c.ra.deg
            delta_P = c.dec.deg
        elif self.ref_coord is None:
            alpha_P = self.header["CRVAL1"]
            delta_P = self.header["CRVAL2"]

        # step1 - intermediate coordinate to (phi, theta)
        phi = np.degrees(np.arctan2(xx, -yy))
        R_th = np.sqrt(xx**2 + yy**2)
        theta = np.degrees(np.arccos(np.pi / 180.0 * R_th))

        # step2 - (phi, theta) -> (alpha, delta) (i.e., RA, Dec)
        cos_th = np.cos(np.radians(theta))
        sin_th = np.sin(np.radians(theta))
        cos_delta_P = np.cos(np.radians(delta_P))
        sin_delta_P = np.sin(np.radians(delta_P))

        argy = -cos_th * np.sin(np.radians(phi - phi_P))
        argx = cos_delta_P * sin_th - sin_delta_P * cos_th * np.cos(
            np.radians(phi - phi_P)
        )
        alpha = np.degrees(np.arctan2(argy, argx)) + alpha_P

        delta = np.degrees(
            np.arcsin(
                sin_delta_P * sin_th
                + cos_delta_P * cos_th * np.cos(np.radians(phi - phi_P))
            )
        )

        # wrap into 0 ~ 360 deg
        self.alpha = alpha % 360.0
        self.delta = delta % 360.0

        if self.relative:
            alpha_0 = alpha_P
            delta_0 = delta_P
            self.dalpha = (self.alpha - alpha_0) * np.cos(np.radians(self.delta))
            self.ddelta = self.delta - delta_0

    def _calc_spectral_axes(self):
        # get the axis index for spectral axis
        ind = [
            k.replace("CTYPE", "")
            for k, v in self.header.items()
            if v == self.spectral_ctype[0]
        ][0]
        self.specax = self.intermediate_world_coord[2] + self.header["CRVAL" + ind]

    def _calc_stokes_axes(self):
        raise NotImplementedError

    # def freq_to_velo_radio(self):
    #     self.v = ac.c.to(u.km/u.s).value * (1. - self.nu / self.nu0)

    # def velo_to_freq_radio(self):
    #     self.nu = (self.nu0 * (1. - self.v / ac.c)).decompose()


if __name__ == "__main__":
    import numpy as np
    import astropy.units as u
    import astropy.io.fits as fits
    import matplotlib.pyplot as plt
    from astropy.visualization import MinMaxInterval, SqrtStretch, ImageNormalize
    import matplotlib.colors as col

    fitsimage = "/raid/work/yamato/eDisk_data/L1489IRS/v0_images/L1489IRS_SBLB_continuum_robust_1.0_taper_3000klambda.image.tt0.fits"

    fda = FitsDataArray(hdu=fits.open(fitsimage)[0])

    data = fda.da

    x_unit = data.coords["dRA cos(Dec)"].attrs["unit"]
    y_unit = data.coords["dDec"].attrs["unit"]

    data.coords["dRA cos(Dec)"].values = (
        (data.coords["dRA cos(Dec)"].values * x_unit).to(u.arcsec).value
    )
    data.coords["dDec"].values = (
        (data.coords["dDec"].values * y_unit).to(u.arcsec).value
    )

    norm = ImageNormalize(data.values, vmin=0.0, stretch=AsinhStretch(a=0.03))

    data.plot.pcolormesh(
        x="dRA cos(Dec)",
        y="dDec",
        rasterized=True,
        xlim=(1, -1),
        ylim=(-1, 1),
        norm=norm,
        cmap="inferno",
    )
    plt.show()
