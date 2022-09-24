import numpy as np
import astropy.constants as ac
import astropy.units as u
import emcee
from galario.double import get_image_size, check_image_size, sampleProfile
from uvplot import UVTable, export_uvtable, COLUMNS_V0
import casatools
import casatasks
import subprocess
import multiprocessing
import matplotlib.pyplot as plt
from scipy import stats

c = ac.c.to(u.m / u.s).value
pi = np.pi
deg_to_rad = pi / 180.0  # in rad
arcsec_to_rad = pi / 180.0 / 3600.0  # in rad

geometrical_param = ["PA", "incl", "dRA", "dDec"]
gp_default = {
    "PA": {"p0": 90.0, "bound": (0.0, 180.0), "fixed": True},
    "incl": {"p0": 45.0, "bound": (0.0, 90.0), "fixed": True},
    "dRA": {"p0": 0.0, "bound": (-5.0, 5.0), "fixed": True},
    "dDec": {"p0": 0.0, "bound": (-5.0, 5.0), "fixed": True},
}


tb = casatools.table()
ms = casatools.ms()


################## BASIC FUNCTIONS #####################
def stderr(x):
    return np.nanstd(x) / x.size

def mean_err(x):
    return np.sqrt(np.nansum(1/x)) / x.size

def weighted_mean_err(x):
    return 1. / np.sqrt(np.nansum(x))

def binned_statistic_1D_base(x, val, statistic="mean", bins=1):
    if np.iscomplexobj(val):
        rstat = stats.binned_statistic(x, val.real, statistic=statistic, bins=bins).statistic
        istat = stats.binned_statistic(x, val.imag, statistic=statistic, bins=bins).statistic
        return rstat + 1j*istat
    else:
        return stats.binned_statistic(x, val, statistic=statistic, bins=bins).statistic

def binned_statistic_1D(x, val, weight=None, statistic="mean", bins=1):
    if weight is not None and statistic == "weighted_mean":
        vw = binned_statistic_1D_base(x, val*weight, statistic=np.nansum, bins=bins)
        w = binned_statistic_1D_base(x, weight, statistic=np.nansum, bins=bins)
        return vw / w
    elif statistic == "weighted_mean" and weight is None:
        raise ValueError("Provide *weight* for 'weighted_mean' statistic.")
    else:
        return binned_statistic_1D_base(x, val, statistic=statistic, bins=bins)


############ function to load the data from msfiles ###############

def load_ms(filenames, datacolumn="data"):
    if isinstance(filenames, str):
        filenames = [filenames]
    
    # initialize the Visibility class
    visibility = Visibility()

    for msfile in filenames:
        print("Loading {}...".format(msfile))

        # open the msfile
        tb.open(msfile)
        
        # read quantities
        spwid = tb.getcol("DATA_DESC_ID") # shape (nvis)
        if len(np.unique(spwid)) != 1:
            print("Warning: Multiple SPWs detected.")
        ant1 = tb.getcol("ANTENNA1")  # shape (nvis,)
        ant2 = tb.getcol("ANTENNA2")  # shape (nvis,)
        uvw = tb.getcol("UVW")  # shape (3, nvis)
        weight = tb.getcol("WEIGHT")  # shape (npol, nvis)
        flag = tb.getcol("FLAG") # shape (npol, nchan, nvis)
        if datacolumn.upper() == "CORRECTED_DATA" and datacolumn.upper() not in tb.colnames():
            print("'CORRECTED_DATA' column not found. Will use 'DATA' column instead.")
            data = tb.getcol("DATA") # shape (npol, nchan, nvis)
        else:
            
            data = tb.getcol(datacolumn.upper())

        # close table instance
        tb.close()

        # channel info
        tb.open(msfile + "/SPECTRAL_WINDOW")
        try:
            nchan = np.unique(tb.getcol("NUM_CHAN")).item()
        except ValueError:
            raise AssertionError("Warning: SPWs have different number of channels (currently not supported). ")
        chan_freqs = tb.getcol("CHAN_FREQ") # shape (nchan, nspw)
        tb.close()

        # get frequency array
        repeats = [len(spwid[spwid == i]) for i in np.unique(spwid)]
        frequency = np.repeat(chan_freqs, repeats, axis=1) # shape (nchan, nvis)
        wavelength = c / frequency

        # broadcast and convert  to lambda
        u, v, w = uvw
        u = u * np.ones((nchan, 1)) / wavelength # shape (nchan, nvis)
        v = v * np.ones((nchan, 1)) / wavelength # shape (nchan, nvis)
        weight = weight[:, np.newaxis, :] # shape (npol, nchan, nvis)

        # average the polarization
        if data.shape[0] == 2:
            data = np.sum(data * weight, axis=0) / np.sum(weight, axis=0) # shape (nchan, nvis)
            weight = np.sum(weight, axis=0) # shape (nchan, nvis)
            flag = np.any(flag, axis=0) # shape (nchan, nvis)
        else:
            data = np.atleast_1d(data.squeeze()) # shape (nchan, nvis)
            weight = np.atleast_1d(weight.squeeze()) # shape (nchan, nvis)
            flag = np.atleast_1d(flag.squeeze()) # shape (nchan, nvis)

        # remove auto-correlation
        xc = ant1 != ant2
        u = u[:, xc]
        v = v[:, xc]
        data = data[:, xc]
        weight = weight[:, xc] 
        frequency = frequency[:, xc]
        flag = flag[:, xc]

        # add to Visibility class; all the quantities are in the shape of (nchan, nvis)
        visibility.add(u, v, data, weight, frequency, flag=flag)

    print("Loading visibility done.")

    return visibility


class Visibility:

    def __init__(self):

        """Store the 1D flattened visibiliy data.

        Parameters
        ----------
        u : 1D or 2D numpy array
            east-west spatial freuencies in lambda. Must be in a 2D shape of (nchan, nvis) or any 1D array.
        v : 1D or 2D numpy array
            north-south spatial freuencies in lambda. Must be in a 2D shape of (nchan, nvis) or any 1D array.
        data : 1D or 2D numpy array
            visibility data for each (u,v) point. Must be in the same shape as *u* and *v*. 
        weight : 1D or 2D numpy array
            visibility weight for each (u,v) point. Must be in the same shape as *u* and *v*. 
        frequency : 1D or 2D numpy array, optional
            frequency for each visibility channel. Must be in the same shape as *u* and *v*. 
        flag : 1D or 2D numpy array, optional
            flag for each visibility data. If None, all data are unflagged.
        """

        self.initialize()

    def initialize(self):
        # initialize
        self.size = 0
        self.u = np.empty(0, dtype=np.float64)
        self.v = np.empty(0, dtype=np.float64)
        self.data = np.empty(0, dtype=np.float64)
        self.weight = np.empty(0, dtype=np.float64)
        self.frequency = np.empty(0, dtype=np.float64)
        
    def add(self, u, v, data, weight, frequency, flag=None):
        
        # number of pol, chan, vis
        datashape = self._get_data_shape(u, v, data, weight, frequency, flag)

        # mask
        mask = ~np.atleast_1d(flag).flatten() if flag is not None else np.atleast_1d(np.ones(datashape).astype(bool)).flatten()

        # spatial frequencies
        u = np.atleast_1d(u).flatten()[mask]
        v = np.atleast_1d(v).flatten()[mask]

        # data 
        data = np.atleast_1d(data).flatten()[mask]

        # weight
        weight = np.atleast_1d(weight).flatten()[mask]

        # frequency
        frequency = np.atleast_1d(frequency).flatten()[mask] 

        # concat
        self.size += data.size
        self.u = np.concatenate((self.u, np.atleast_1d(u)))
        self.v = np.concatenate((self.v, np.atleast_1d(v)))
        self.data = np.concatenate((self.data, np.atleast_1d(data)))
        self.weight = np.concatenate((self.weight, np.atleast_1d(weight)))
        self.frequency = np.concatenate((self.frequency, np.atleast_1d(frequency)))

    def to_npz(self, filename):
        np.savez(filename, u=self.u, v=self.v, data=self.data, weight=self.weight, frequency=self.frequency)

    @staticmethod
    def _get_data_shape(*args):

        # get the reference data shape
        shape = args[0].shape

        # iterate over args to check data shape consistency
        for a in args:
            if a is None:
                continue
            if a.shape != shape:
                raise AssertionError("Mismatch in the shape of input data.")

        # data dimension
        ndim = len(shape)

        # get the data shapes
        nvis = shape[-1]
        nchan = shape[-2] if ndim >= 2 else 1
        npol = shape[-3] if ndim == 3 else 1

        return npol, nchan, nvis

    @property
    def uvdist(self):
        return np.hypot(self.u, self.v)

    def bin_1D(self, binsize=10e3, uvrange=None, stat_type="weighted_mean", uncertainty_type="errprop"):
        if uvrange is None:
            uvbins = np.arange(0, np.nanmax(self.uvdist), binsize)
        else:
            uvbins = np.arange(*uvrange, binsize)

        uvvals = np.average([uvbins[1:], uvbins[:-1]], axis=0)

        if stat_type == "mean":
            vals = binned_statistic_1D(self.uvdist, self.data, statistic=np.nanmean, bins=uvbins)
            if uncertainty_type == "stderr":
                errs = binned_statistic_1D(self.uvdist, self.data, statistic=stderr, bins=uvbins)
            elif uncertainty_type == "errprop":
                errs = binned_statistic_1D(self.uvdist, self.weight, statistic=mean_err, bins=uvbins)

        elif stat_type == "weighted_mean":
            vals = binned_statistic_1D(self.uvdist, self.data, weight=self.weight, statistic="weighted_mean", bins=uvbins)
            if uncertainty_type == "stderr":
                errs = binned_statistic_1D(self.uvdist, self.data, statistic=stderr, bins=uvbins)
            elif uncertainty_type == "errprop":
                errs = binned_statistic_1D(self.uvdist, self.weight, statistic=weighted_mean_err, bins=uvbins)

        return uvvals, vals, errs

    def plot_uvprofile(self, axes=None, binsize=10e3, uvrange=None, stat_type="weighted_mean", uncertainty_type="errprop", errorbar=True, **kwargs):

        if axes is None:
            fig, axes = plt.subplots(2, 1, figsize=(4, 8), sharex=True)

        uv, vals, errs = self.bin_1D(binsize=binsize, uvrange=uvrange, stat_type=stat_type, uncertainty_type=uncertainty_type)

        for ax, v, e in zip(axes, [vals.real, vals.imag], [errs.real, errs.imag]):
            if errorbar:
                ax.errorbar(uv, v, yerr=e, **kwargs)
            else:
                ax.plot(uv, v)

        ax[0].set(ylabel="Real [Jy]", xscale="log")
        ax[1].set(xlabel="Baseline [$\lambda$]", ylabel="Imaginary [Jy]", xscale="log")

        return axes

        

class VisFit1D:
    def __init__(self, msfilenames, verbose=True):

        self.obsvis = load_ms(msfilenames)
        self.verbose=verbose

    
    def _log_likelihood(self, param):

        param_dict = {name: p for name, p in zip(self.param_name, param)}

        # update fixed param
        param_dict.update({name: self.param_dict[name]["p0"] for name in self.fixed_param_name})

        # sample model visibility
        model_vis = self.sample_vis(param_dict)

        # compute log likelihood
        rms = np.sqrt(1.0 / self.obsvis.weight)
        ll = -0.5 * np.sum(
            ((model_vis.real - self.obsvis.data.real) ** 2 + (model_vis.imag - self.obsvis.data.imag) ** 2) / rms ** 2
            + np.log(2 * pi * rms ** 2)
        )

        return ll

    def _log_probability(self, param):
        lp = self._log_prior(param, self.bound)
        if not np.isfinite(lp):
            return -np.inf
        ll = self._log_likelihood(param)
        return lp + ll

    @staticmethod
    def _initialize_walker(p0, nwalker=200, scatter=1e-4):
        p0 = p0 + scatter * np.random.randn(nwalker, len(p0))
        return p0

    def _run_mcmc(self, p0, nwalker, nburnin, nstep, pool=None, **kwargs):

        p0 = self._initialize_walker(p0, nwalker=nwalker, scatter=kwargs.get("scatter", 1e-4))

        sampler = emcee.EnsembleSampler(nwalker, p0.shape[1], self._log_probability, pool=pool)
        
        progress = kwargs.pop("progress", True)
        sampler.run_mcmc(p0, nburnin + nstep, progress=progress, **kwargs)

        return sampler

    def fit_vis(
        self, model_func, param_dict, nwalker=32, nstep=500, nburnin=500, pool=None, progress=True
    ):
        # setup radial grid
        self.set_grid()

        self.param_dict = param_dict
        self.model_func = model_func

        # get lists of free/fixed parameter names and bounds
        self.param_name = [key for key in self.param_dict.keys() if not self.param_dict[key]["fixed"]]
        self.fixed_param_name = [key for key in self.param_dict.keys() if self.param_dict[key]["fixed"]]
        self.bound = [
            self.param_dict[key]["bound"] for key in self.param_dict.keys() if not self.param_dict[key]["fixed"]
        ]

        # initial state
        p0 = [self.param_dict[key]["p0"] for key in self.param_dict.keys() if not self.param_dict[key]["fixed"]]

        # run
        print(
            "starting to run the MCMC sampling with: \n \t initial state:",
            p0,
            "\n \t number of walkers:",
            nwalker,
            "\n \t number of steps:",
            nstep,
        )
        sampler = self._run_mcmc(p0, nwalker=nwalker, nburnin=nburnin, nstep=nstep, pool=pool, progress=progress)

        return sampler

    def sample_vis(self, param_dict):

        # retrieve geometrical params
        PA = param_dict.pop("PA", gp_default["PA"]["p0"])
        incl = param_dict.pop("incl", gp_default["incl"]["p0"])
        dRA = param_dict.pop("dRA", gp_default["dRA"]["p0"])
        dDec = param_dict.pop("dDec", gp_default["dDec"]["p0"])

        # get model array
        model = self.model_func(self.r, **param_dict)

        # sampling by GALARIO
        vis = sampleProfile(
            intensity=model,
            Rmin=self.rmin,
            dR=self.dr,
            nxy=self.nxy,
            dxy=self.dxy,
            u=self.obsvis.u,
            v=self.obsvis.v,
            dRA=dRA * arcsec_to_rad,
            dDec=dDec * arcsec_to_rad,
            PA=PA * deg_to_rad,
            inc=incl * deg_to_rad,
            check=False,
        )

        return vis

    def set_grid(self):
        # gridding parameter
        self.nxy, self.dxy = get_image_size(
            self.obsvis.u, self.obsvis.v, verbose=self.verbose, f_min=1.0, f_max=2.0
        )  # in rad

        # condition for GALARIO interpolation: dxy/2 - Rmin > 5 * dR
        # not to make dR too small, adopt dxy/2/1000 as Rmin
        self.rmin = self.dxy / 2.0 * 1e-3  # in rad

        # dR = (dxy/2 - Rmin) / 5.1
        self.dr = (self.dxy / 2.0 - self.rmin) / 5.1  # in rad

        r_pad = 2 * self.dxy  # padding parameter
        self.rmax = self.dxy * self.nxy / np.sqrt(2) + r_pad  # in rad

        # radial grid in image plane
        self.r = np.arange(self.rmin, self.rmax, self.dr)  # in rad

    @staticmethod
    def _condition(p, b):
        if (p >= b[0]) and (p <= b[1]):
            return True
        return False

    def _log_prior(self, param, bound):
        for p, b in zip(param, bound):
            if not self._condition(p, b):
                return -np.inf
        return 0.0


if __name__ == '__main__':

    def Gaussian1d(r, F, sigma):
        return 10 ** F / (np.sqrt(2 * np.pi) * sigma) * np.exp(-r**2 / (2 * sigma**2))

    def GaussianRing1d(r, r0, F, sigma):
        return 10 ** F / (np.sqrt(2 * np.pi) * sigma) * np.exp(-(r - r0)**2 / (2 * sigma**2))

    def model_func_1d(r, F_c, sigma_c, r0_r, F_r, sigma_r, F_b, sigma_b):
        return (
            Gaussian1d(r, F_c, sigma_c)
            + GaussianRing1d(r, r0_r, F_r, sigma_r)
            + Gaussian1d(r, F_b, sigma_b)
        )

    # dictionary of parameters including initial guess, bound, and fixed/free
    param_dict = {
        # "I_g": {"p0": 8.5, "bound": (-2.0, 20.0), "fixed": False},
        # "sigma_g": {"p0": 2.0, "bound": (0.1, 10), "fixed": False},
        "F_c": {"p0": 10.59, "bound": (5, 15), "fixed": False},
        "sigma_c": {"p0": 0.01, "bound": (1e-5, 0.2), "fixed": False},
        "r0_r": {"p0": 0.47, "bound": (0.01, 1.5), "fixed": False},
        "F_r": {"p0": 8.31, "bound": (3, 13), "fixed": False},
        "sigma_r": {"p0": 0.20, "bound": (1e-2, 1.5), "fixed": False},
        "F_b": {"p0": 9.14, "bound": (4, 14), "fixed": False},
        "sigma_b": {"p0": 1.59, "bound": (0.3, 5), "fixed": False},
        "PA": {"p0": 68.9, "bound": (0, 180), "fixed": False},
        "incl": {"p0": 72.86, "bound": (0, 90), "fixed": False},
        "dRA": {"p0": 0.0, "bound": (-2, 2), "fixed": False},
        "dDec": {"p0": 0.0, "bound": (-2, 2), "fixed": False},
    }

    path = "/works/yamato/eDisk/L1489IRS/ALMA_pipeline_calibrated_data/"
    msfilenames = [path + f"L1489IRS_{i}_continuum_shift.bin_30s.split.ms" for i in ["SB1", "SB2", "SB3", "LB1", "LB2"]]
    # vis = load_ms(msfilenames)
    # vis.to_npz("./testvis.npz")
    vf = VisFit1D(msfilenames)

    nwalker = 32
    nstep = 2500
    nburnin = 2500
    progress = True

    with multiprocessing.Pool(processes=8) as pool:
    # pool = None
        vf.fit_vis(model_func=model_func_1d, param_dict=param_dict, nwalker=nwalker, nstep=nstep, nburnin=nburnin, pool=pool)
