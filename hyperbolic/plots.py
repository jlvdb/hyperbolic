import os

import numpy as np
import pandas as pd
import scipy.stats
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from . import Keys
from .logger import logger
from .magnitudes import compute_classic_magnitude
import hyperbolic


class SparseSampler:

    def __init__(self, total_length, samples=50000):
        sparse = np.arange(total_length)
        np.random.shuffle(sparse)
        self.sparse = sparse[:min(total_length, samples)]

    def apply(self, data, mask=None):
        try:
            data = data.to_numpy()
        except AttributeError:
            pass
        if mask is not None:
            data = data[mask]
        return data[self.sparse]


def add_left_edge(centers, edges):
    try:
        centers = centers.to_numpy()
    except AttributeError:
        pass
    try:
        edges = edges.to_numpy()
    except AttributeError:
        pass
    return np.append(edges[0], centers)


def compute_prediction(sn, b):
    # calculate normalised flux from relation of b and S/N
    norm_flux = sn * b / np.sqrt(hyperbolic.pogson)
    return hyperbolic.compute_magnitude(norm_flux, b)


def fig_add_xlabel(axes, label, offset=0):
    ncols = axes.shape[1]
    for i in range(axes.size+offset-ncols, axes.size+offset):
        axes.flatten()[i].set_xlabel(label)


def fig_add_ylabel(axes, label):
    for ax in axes[:, 0]:
        ax.set_ylabel(label)


class PlotLims:

    def __init__(self):
        self.lower = np.inf
        self.upper = -np.inf

    def update(self, lower=None, upper=None):
        if lower is not None:
            self.lower = min(lower, self.lower)
        if upper is not None:
            self.upper = max(upper, self.upper)

    def get(self, reverse=False):
        if reverse:
            return (self.upper, self.lower)
        else:
            return (self.lower, self.upper)


class Plotter:

    _magnitudes = None

    scatter_style = {
        "edgecolor": "none", "s": 3, "marker": ".",
        "alpha": 0.1, "rasterized": True}

    def __init__(self, config, data, stats, b):
        self.config = config
        self.data = data
        self.stats = stats
        self.b = b

    def make_figure(self, size=2.5, sharex=True, sharey=True, n_plot_offset=0):
        n_plots = len(self.config.filters) + n_plot_offset
        ncols = min(3, n_plots)
        nrows = n_plots // ncols
        if nrows * ncols < n_plots:
            nrows += 1
        fig, axes = plt.subplots(
            nrows, ncols, figsize=(0.5 + size*ncols, 0.5 + size*ncols),
            sharex=sharex, sharey=sharey)
        for i, ax in enumerate(axes.flatten()):
            if i >= n_plots:
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                plt.axis("off")
            else:
                for pos in ["top", "right"]:
                    ax.spines[pos].set_visible(False)
                ax.grid(alpha=0.33)
        return fig, axes

    def __enter__(self, *args, **kwargs):
        fname = os.path.splitext(self.config.outfile)[0] + ".pdf"
        logger.info(f"plotting to {fname}")
        self._backend = PdfPages(fname)
        return self

    def __exit__(self, *args, **kwargs):
        self._backend.close()

    def add_fig(self, fig):
        self._backend.savefig(fig)

    def plot_b(self, smoothing):
        # collect the b values
        b_median = {}
        b_all = {}
        for filt, data in smoothing.groupby(Keys.filter):
            b_all[filt] = data[Keys.b_abs].to_numpy()
            b_median[filt] = np.nanmedian(b_all[filt])
        # make figure
        fig, axes = self.make_figure(sharex=False)
        for i, filt in enumerate(self.config.filters):
            ax = axes.flatten()[i]
            _, bins, patches = ax.hist(b_all[filt], histtype="step")
            hand2 = ax.axvline(x=b_median[filt], color="k")
            # decorate
            ax.annotate(
                filt, (0.95, 0.95), xycoords="axes fraction",
                ha="right", va="top")
        fig_add_xlabel(axes, "Smoothing $b^\prime$")
        # add legend and fix layout
        fig.legend(
            handles=[patches[0], hand2],
            labels=["fields", "median (fiducial)"],
            ncol=2, frameon=False)
        fig.tight_layout()
        fig.subplots_adjust(top=0.94)
        self.add_fig(fig)

    def get_magnitudes(self):
        if self._magnitudes is not None:
            return self._magnitudes
        magnitudes = self.config.get_magnitudes(self.data)
        # compute classic magnitudes as fallback from fluxes and zeropoints
        if magnitudes is None:
            fields = self.config.get_fields(self.data)
            fluxes = self.config.get_fluxes(self.data)
            magnitudes = {}
            for filt in self.config.filters:
                zeropoint = self.stats[Keys.zp].loc[fields]
                zeropoint.index = fluxes[filt].index
                magnitudes = compute_classic_magnitude(fluxes[filt], zeropoint)
        self._magnitudes = magnitudes
        return magnitudes

    def plot_magnitudes(self):
        # get the required data
        fluxes = self.config.get_fluxes(self.data)
        errors = self.config.get_errors(self.data)
        magnitudes = self.get_magnitudes()
        # make figure
        xlims = [-3, 33]
        ylims = PlotLims()
        fig, axes = self.make_figure()
        for i, filt in enumerate(self.config.filters):
            # select all observed objects and create a sparse sampling
            ax = axes.flatten()[i]
            is_good = errors[filt] > 0.0
            sparse = SparseSampler(np.count_nonzero(is_good))
            SN = (sparse.apply(fluxes[filt], mask=is_good) /
                  sparse.apply(errors[filt], mask=is_good))
            # add prediction
            SN_theory = np.linspace(*xlims, 50)
            mag_theory = compute_prediction(SN_theory, self.b[filt])
            hand_theory = ax.plot(
                SN_theory, mag_theory, color="k", lw=0.7, ls="--", zorder=2)[0]
            ylims.update(  # update limits from y=mag(-b) to y=max(min)
                lower=mag_theory.min(),
                upper=hyperbolic.compute_magnitude(
                    -self.b[filt], self.b[filt]))
            # add classical magnitudes
            mags = sparse.apply(magnitudes[filt], mask=is_good)
            hand1 = ax.scatter(
                SN, mags, c="C0", zorder=-1, **self.scatter_style)
            # add hyperbolic magnitudes
            key_mag = self.config.outname[filt]
            mags = sparse.apply(self.data[key_mag], mask=is_good)
            hand2 = ax.scatter(
                SN, mags, c="C3", **self.scatter_style)
            # decorate
            ax.axhline(
                y=hyperbolic.compute_magnitude(0.0, self.b[filt]),
                color="k", lw=0.4)
            ax.axvline(x=0.0, color="k", lw=0.25)
            ax.annotate(
                filt, (0.1, 0.95), xycoords="axes fraction",
                ha="left", va="top")
            ax.set_xlim(*xlims)
            ax.set_ylim(*ylims.get(reverse=True))
        fig_add_xlabel(axes, "Signal-to-noise")
        fig_add_ylabel(axes, "Magnitude")
        # add legend and fix layout
        fig.legend(
            handles=[hand_theory, hand1, hand2],
            labels=["prediction", "classical", "hyperbolic"],
            markerscale=8, ncol=3, frameon=False)
        fig.tight_layout()
        fig.subplots_adjust(top=0.94)
        self.add_fig(fig)

    def plot_magnitude_distribution(self):
        # get the required data
        errors = self.config.get_errors(self.data)
        magnitudes = self.get_magnitudes()
        # make figure
        bins = np.arange(5, 40, 0.2)
        xlims = PlotLims()
        ylims = PlotLims()
        fig, axes = self.make_figure()
        for i, filt in enumerate(self.config.filters):
            # select all observed objects and create a sparse sampling
            ax = axes.flatten()[i]
            is_good = errors[filt] > 0.0
            has_mag = magnitudes[filt] < 90.0
            # add classical magnitudes
            c, _, patches1 = ax.hist(
                magnitudes[filt][is_good & has_mag], bins,
                log=True, color="C0", histtype="step")
            # add hyperbolic magnitudes
            key = self.config.outname[filt]
            c, _, patches2 = ax.hist(
                self.data[key][is_good], bins,
                log=True, color="C3", histtype="step")
            # decorate
            zeroflux = hyperbolic.compute_magnitude(0.0, self.b[filt])
            hand0 = ax.axvline(x=zeroflux, color="k", lw=0.7, ls="--")
            xlims.update(self.data[key][is_good].min(), zeroflux + 5)
            ax.annotate(
                filt, (0.05, 0.95), xycoords="axes fraction",
                ha="left", va="top")
            ax.set_xlim(*xlims.get())
            ylims.update(lower=1.0, upper=ax.get_ylim()[1])
            ax.set_ylim(*ylim.get())
        fig_add_xlabel(axes, "Magnitude")
        fig_add_ylabel(axes, "Frequency")
        # add legend and fix layout
        fig.legend(
            handles=[hand0, patches1[0], patches2[0]],
            labels=["zero-flux magnitude", "classical", "hyperbolic"],
            ncol=3, frameon=False)
        fig.tight_layout()
        fig.subplots_adjust(top=0.94)
        self.add_fig(fig)


    def plot_colour_distribution(self):
        # get the required data
        errors = self.config.get_errors(self.data)
        magnitudes = self.get_magnitudes()
        # make figure
        bins = np.linspace(-1.5, 2.5, 50)
        fig, axes = self.make_figure(n_plot_offset=-1)
        for i, (filt1, filt2) in enumerate(
                zip(self.config.filters[:-1], self.config.filters[1:])):
            ax = axes.flatten()[i]
            # select all observed objects and create a sparse sampling
            mag1 = magnitudes[filt1]
            mag2 = magnitudes[filt2]
            is_good = (errors[filt1] > 0.0) & (errors[filt2] > 0.0)
            has_color = np.isfinite(mag1) & (mag1 < 90.0)
            has_color &= np.isfinite(mag2) & (mag2 < 90.0)
            mask = is_good & has_color
            # add classical colours
            mag1 = magnitudes[filt1][mask].to_numpy()
            mag2 = magnitudes[filt2][mask].to_numpy()
            idx_sort = np.argsort(mag1)
            colours = mag1[idx_sort] - mag2[idx_sort]
            hand1a = ax.hist(
                colours[:len(colours)//4], bins,
                color="C0", alpha=0.3, zorder=-1)[-1][0]
            hand1b = ax.hist(
                colours[-len(colours)//4:], bins,
                color="C0", histtype="step")[-1][0]
            # add hyperbolic colours
            key_mag1 = self.config.outname[filt1]
            key_mag2 = self.config.outname[filt2]
            mag1 = self.data[key_mag1][mask].to_numpy()
            mag2 = self.data[key_mag2][mask].to_numpy()
            idx_sort = np.argsort(mag1)
            colours = (mag1[idx_sort] - mag2[idx_sort])
            hand2a = ax.hist(
                colours[:len(colours)//4], bins,
                color="C3", alpha=0.3, zorder=-1)[-1][0]
            hand2b = ax.hist(
                colours[-len(colours)//4:], bins,
                color="C3", histtype="step")[-1][0]
            # decorate
            ax.axvline(x=0.0, color="k", lw=0.25)
            ax.annotate(
                f"{filt1}$-${filt2}", (0.1, 0.95), xycoords="axes fraction",
                ha="left", va="top")
            ax.set_xlim(bins[0], bins[-1])
        fig_add_xlabel(axes, "Colour", offset=-1)
        fig_add_ylabel(axes, "Frequency")
        # add legend and fix layout
        axes.flatten()[-1].legend(
            loc="center",
            handles=[hand1a, hand1b, hand2a, hand2b],
            labels=[
                "classical\n(bright quartile)",
                "classical\n(faint quartile)",
                "hyperbolic\n(bright quartile)",
                "hyperbolic\n(faint quartile)"],
            ncol=1, frameon=False)
        fig.tight_layout()
        fig.subplots_adjust()
        self.add_fig(fig)

    def plot_magdiff(self):
        # get the required data
        errors = self.config.get_errors(self.data)
        magnitudes = self.get_magnitudes()
        # make figure
        fig, axes = self.make_figure(sharey=False)
        xlims = PlotLims()
        for i, filt in enumerate(self.config.filters):
            # select all observed objects and create a sparse sampling
            ax = axes.flatten()[i]
            is_good = (errors[filt] > 0.0) & (magnitudes[filt] < 90.0)
            # collect the data
            key_mag = self.config.outname[filt]
            df = pd.DataFrame({
                "mag": magnitudes[filt][is_good],
                "hyp": self.data[key_mag][is_good]})
            df["diff"] = df["mag"] - df["hyp"]
            sparse = SparseSampler(np.count_nonzero(is_good))
            # add the magnitude difference
            ax.scatter(
                sparse.apply(df["mag"]), sparse.apply(df["diff"]),
                **self.scatter_style)
            # plot statistics
            bins = np.linspace(*np.percentile(df["mag"], q=[0.5, 99.5]), 20)
            centers = (bins[1:] + bins[:-1]) / 2.0
            stats = df.groupby(pd.cut(df["mag"], bins)).agg([
                np.median, scipy.stats.median_abs_deviation])
            y = stats["diff"]["median"]
            ylow = y - stats["diff"]["median_abs_deviation"]
            yhigh = y + stats["diff"]["median_abs_deviation"]
            ax.plot(
                add_left_edge(centers, bins), add_left_edge(y, y),
                color="k", lw=0.7)
            ax.plot(
                add_left_edge(centers, bins), add_left_edge(ylow, ylow),
                color="k", ls="--", lw=0.7)
            ax.plot(
                add_left_edge(centers, bins), add_left_edge(yhigh, yhigh),
                color="k", ls="--", lw=0.7)
            # decorate
            zeroflux = hyperbolic.compute_magnitude(0.0, self.b[filt])
            xlims.update(bins[0], zeroflux)
            ax.set_xlim(*xlims.get())
            ax.set_ylim(-0.2, 0.6)
            ax.annotate(
                filt, (0.1, 0.95), xycoords="axes fraction",
                ha="left", va="top")
        fig_add_xlabel(axes, "Classical")
        fig_add_ylabel(axes, "Classical $-$ hyperbolic")
        # add legend and fix layout
        fig.tight_layout()
        fig.subplots_adjust()
        self.add_fig(fig)
