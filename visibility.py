import numpy as np
import casatools 
from scipy import stats
import astropy.constants as ac
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

c = ac.c.to(u.m/u.s).value
arcsec_to_deg = 1. / 3600.

ms = casatools.ms()

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

def binned_statistic_2D_base(x, y, val, statistic="mean", bins=1):
    if np.iscomplexobj(val):
        rstat = stats.binned_statistic_2d(x, y, val.real, statistic=statistic, bins=bins).statistic
        istat = stats.binned_statistic_2d(x, y, val.imag, statistic=statistic, bins=bins).statistic
        return rstat + 1j*istat
    else:
        return stats.binned_statistic_2d(x, y, val, statistic=statistic, bins=bins).statistic

def binned_statistic_2D(x, y, val, weight=None, statistic="mean", bins=1):
    if weight is not None and statistic == "weighted_mean":
        vw = binned_statistic_2D_base(x, y, val*weight, statistic=np.nansum, bins=bins)
        w = binned_statistic_2D_base(x, y, weight, statistic=np.nansum, bins=bins)
        return vw / w
    elif statistic == "weighted_mean" and weight is None:
        raise ValueError("Provide *weight* for 'weighted_mean' statistic.")
    else:
        return binned_statistic_2D_base(x, y, val, statistic=statistic, bins=bins)
        

class Visibility():

    def __init__(self, filenames):
        
        if isinstance(filenames, str) and filenames.endswith(".npz"):
            data = np.load(filenames)
            self.u = data["u"]
            self.v = data["v"]
            self.uvdist = np.hypot(self.u, self.v)
            self.V = data["V"]
            self.real = self.V.real
            self.imag = self.V.imag
            self.weight = data["weight"]

        else:
            if isinstance(filenames, str):
                self.filenames = [filenames]
            else:
                self.filenames = filenames
            
            self._load_ms()


    def _load_ms(self):

        
        self.u = []
        self.v = []
        self.V = []
        self.weight = []
        # self.freqs = []

        for filename in self.filenames:
            print("Loading {}...".format(filename))
            ms.open(filename)

            for i in self._get_spwids():
                print("Loading data in spw {}...".format(i))
                u, v = self._get_uv_axis(spwid=i)
                V, weight = self._get_visibility(spwid=i)
                # freqs = self._get_frequency_axis(spwid=i)
                uvdist = self._get_uvdistance(spwid=i)
                flag = self._get_flag(spwid=i)

                # remove flag
                u = u[~flag]
                v = v[~flag]
                V = V[~flag]
                weight = weight[~flag]
                uvdist = uvdist[~flag]
                # freqs = freqs[~flag]

                # remove autocorrelation
                cc = uvdist != 0.0
                u = u[cc]
                v = v[cc]
                V = V[cc]
                weight = weight[cc]
                # freqs = freqs[cc]

                # flatten and append
                self.u.append(u.flatten())
                self.v.append(v.flatten())
                self.V.append(V.flatten())
                self.weight.append(weight.flatten())
                # self.freqs.append(freqs.flatten())
            
            ms.close()

        # concatenate 
        self.u = np.ascontiguousarray(np.concatenate(self.u))
        self.v = np.ascontiguousarray(np.concatenate(self.v))
        self.uvdist = np.ascontiguousarray(np.hypot(self.u, self.v))
        self.V = np.ascontiguousarray(np.concatenate(self.V))
        self.real = np.ascontiguousarray(self.V.real)
        self.imag = np.ascontiguousarray(self.V.imag)
        self.weight = np.ascontiguousarray(np.concatenate(self.weight))
        # self.freqs = np.ascontiguousarray(np.concatenate(self.freqs))

        print("Done.")

    def deproject(self, incl=0.0, PA=0.0):
        incl = np.radians(incl)
        PA = np.radians(PA)

        u_dep = (self.u * np.cos(PA) - self.v * np.sin(PA)) * np.cos(incl)
        v_dep = self.u * np.sin(PA) + self.v * np.cos(PA)

        self.u, self.v = u_dep, v_dep
        self.uvdist = np.hypot(self.u, self.v)

    def shift_phase(self, dRA=0.0, dDec=0.0):
        if dRA == 0.0 and dDec == 0.0:
            return self.real, self.imag

        dRA = 2 * np.pi * np.radians(dRA*arcsec_to_deg)
        dDec = 2 * np.pi * np.radians(dDec*arcsec_to_deg)

        phi = self.u * dRA + self.v * dDec
        self.V = (self.real + 1j * self.imag) * (np.cos(phi) + 1j * np.sin(phi))

    # def shift_phase_toward(self, coord, frame="icrs"):

    #     pc_from = self._get_phasecenter()
    #     pc_to = SkyCoord(coord, frame=frame)
    #     dRA, dDec = pc_from.spherical_offsets_to(pc_to)

    #     return

    def export_vis(self, filename):
        np.savez(filename, u=self.u, v=self.v, V=self.V, weight=self.weight)

    def bin_1D(self, binsize=10e3, uvrange=None, stat_type="weighted_mean", uncertainty_type="errprop"):
        if uvrange is None:
            uvbins = np.arange(0, np.nanmax(self.uvdist), binsize)
        else:
            uvbins = np.arange(*uvrange, binsize)

        uvvals = np.average([uvbins[1:], uvbins[:-1]], axis=0)

        if stat_type == "mean":
            vals = binned_statistic_1D(self.uvdist, self.V, statistic=np.nanmean, bins=uvbins)
            if uncertainty_type == "stderr":
                errs = binned_statistic_1D(self.uvdist, self.V, statistic=stderr, bins=uvbins)
            elif uncertainty_type == "errprop":
                errs = binned_statistic_1D(self.uvdist, self.weight, statistic=mean_err, bins=uvbins)

        elif stat_type == "weighted_mean":
            vals = binned_statistic_1D(self.uvdist, self.V, weight=self.weight, statistic="weighted_mean", bins=uvbins)
            if uncertainty_type == "stderr":
                errs = binned_statistic_1D(self.uvdist, self.V, statistic=stderr, bins=uvbins)
            elif uncertainty_type == "errprop":
                errs = binned_statistic_1D(self.uvdist, self.weight, statistic=weighted_mean_err, bins=uvbins)

        return uvvals, vals, errs
    
    def bin_2D(self, binsize=10e3, uvrange=None, stat_type="weighted_mean", uncertainty_type="errprop", replace=True):
        if uvrange is None:
            uvmax = np.nanmax(np.abs([self.u, self.v]))
            uvbins = np.arange(-uvmax, uvmax, binsize)
        else:
            uvbins = np.arange(*uvrange, binsize)

        uvvals = np.average([uvbins[1:], uvbins[:-1]], axis=0)

        if stat_type == "mean":
            vals = binned_statistic_2D(self.v, self.u, self.V, statistic=np.nanmean, bins=uvbins)
            if uncertainty_type == "stderr":
                errs = binned_statistic_2D(self.v, self.u, self.V, statistic=stderr, bins=uvbins)
            elif uncertainty_type == "errprop":
                errs = binned_statistic_2D(self.v, self.u, self.weight, statistic=mean_err, bins=uvbins)

        elif stat_type == "weighted_mean":
            vals = binned_statistic_2D(self.v, self.u, self.V, weight=self.weight, statistic="weighted_mean", bins=uvbins)
            if uncertainty_type == "stderr":
                errs = binned_statistic_2D(self.v, self.u, self.V, statistic=stderr, bins=uvbins)
            elif uncertainty_type == "errprop":
                errs = binned_statistic_2D(self.v, self.u, self.weight, statistic=weighted_mean_err, bins=uvbins)

        # make arrays flatten
        u, v = np.meshgrid(uvvals, uvvals)
        uvals = u.flatten()
        vvals = v.flatten()
        vals = vals.flatten()
        errs = errs.flatten()

        if replace:
            self.u = np.ascontiguousarray(uvals)
            self.v = np.ascontiguousarray(vvals)
            self.uvdist = np.ascontiguousarray(np.hypot(self.u, self.v))
            self.V = np.ascontiguousarray(vals)
            self.real = np.ascontiguousarray(self.V.real)
            self.imag = np.ascontiguousarray(self.V.imag)
            self.weight = np.ascontiguousarray(1. / errs**2)
        else:
            return uvals, vvals, vals, errs

    def plot_uvprofile(self, axes=None, binsize=10e3, uvrange=None, stat_type="weighted_mean", uncertainty_type="errprop", errorbar=True, **kwargs):

        if axes is None:
            fig, axes = plt.subplots(2, 1, figsize=(4, 8), sharex=True)

        uv, vals, errs = self.bin_1D(binsize=binsize, uvrange=uvrange, stat_type=stat_type, uncertainty_type=uncertainty_type)

        for ax, v, e in zip(axes, [vals.real, vals.imag], [errs.real, errs.imag]):
            if errorbar:
                ax.errorbar(uv, v, yerr=e, **kwargs)
            else:
                ax.plot(uv, v)

        axes[0].set(ylabel="Real [Jy]", xscale="log")
        axes[1].set(xlabel="Baseline [$\lambda$]", ylabel="Imaginary [Jy]", xscale="log")

        return axes


    def _get_spwids(self):
        spwinfo = ms.getspectralwindowinfo()
        return [spwinfo[i]["SpectralWindowId"] for i in spwinfo.keys()]

    def _get_nchan(self, spwid):
        spwinfo = ms.getspectralwindowinfo()
        return spwinfo[str(spwid)]["NumChan"]

    def _get_frequency_axis(self, spwid):

        ms.selectinit(datadescid=int(spwid))
        data = ms.getdata(["axis_info"])

        ms.reset()

        return np.atleast_1d(data["axis_info"]["freq_axis"]["chan_freq"].squeeze())

    def _get_uv_axis(self, spwid, in_lambda=True):

        ms.selectinit(datadescid=int(spwid))

        data = ms.getdata(["u", "v"])

        nchan = self._get_nchan(spwid)
        u = np.broadcast_to(data["u"], (nchan, data["u"].size)).copy()
        v = np.broadcast_to(data["v"], (nchan, data["v"].size)).copy()

        if in_lambda:
            nu = self._get_frequency_axis(spwid)
            if nu.size != nchan:
                raise ValueError("Mismatch in the number of channels and frequency axis size.")
            u /= np.broadcast_to(c / nu[:, np.newaxis], u.shape).copy()
            v /= np.broadcast_to(c / nu[:, np.newaxis], v.shape).copy()

        ms.reset()

        return u, v

    def _get_uvdistance(self, spwid, in_lambda=True):

        ms.selectinit(datadescid=int(spwid))

        data = ms.getdata(["uvdist"])

        nchan = self._get_nchan(spwid)
        uvdist = np.broadcast_to(data["uvdist"], (nchan, data["uvdist"].size)).copy()

        if in_lambda:
            nu = self._get_frequency_axis(spwid)
            if nu.size != nchan:
                raise ValueError("Mismatch in the number of channels and frequency axis size.")
            uvdist /= np.broadcast_to(c / nu[:, np.newaxis], uvdist.shape).copy()

        ms.reset()

        return uvdist

    def _get_visibility(self, spwid, average_pol=True):

        ms.selectinit(datadescid=int(spwid))

        data = ms.getdata(["data", "weight"])

        V = data["data"]
        weight = np.broadcast_to(data["weight"][:, np.newaxis, :], V.shape).copy() # to make the shape of weight same as V

        if average_pol:
            if V.shape[0] == 2:
                V = np.sum(V * weight, axis=0) / np.sum(weight, axis=0)
                weight = np.sum(weight, axis=0)
            else:
                V = np.atleast_1d(V.squeeze())
                weight = np.atleast_1d(weight.squeeze())

        ms.reset()
        
        return V, weight

    def _get_flag(self, spwid, average_pol=True):

        ms.selectinit(datadescid=int(spwid))

        data = ms.getdata(["flag"])

        flag = data["flag"]

        if average_pol:
            if flag.shape[0] == 2:
                flag = np.any(flag, axis=0)
            else:
                flag = np.atleast_1d(flag.squeeze())

        ms.reset()

        return flag

    def _get_metadata(self):
        return ms.metadata()

    def _get_phasecenter(self):
        metadata = self._get_metadata()
        return self._read_phasecenter(metadata.phasecenter())

    def _read_phasecenter(self, record, repre=None):

        ra = record["m0"]["value"] * u.Unit(record["m0"]["unit"])
        dec = record["m1"]["value"] * u.Unit(record["m1"]["unit"])
        frame = record["refer"].lower()

        if repre is None:
            return SkyCoord(ra=ra, dec=dec, frame=frame)
        else:
            return SkyCoord(ra=ra, dec=dec, frame=frame).to_string(repre)



    








