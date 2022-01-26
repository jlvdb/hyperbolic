import os

import numpy as np
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


class Plotter:

    def __init__(self, config):
        self.config = config

    def make_figure(self, size=2, sharex=True, sharey=True, n_plot_offset=0):
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

    def plot_b(self, stats):
        b_median = {}
        b_all = {}
        for filt, data in stats.groupby(Keys.filter):
            b_all[filt] = data[Keys.b_abs].to_numpy()
            b_median[filt] = np.nanmedian(b_all[filt])
        # make figure
        fig, axes = self.make_figure(sharex=False)
        for i, filt in enumerate(self.config.filters):
            ax = axes.flatten()[i]
            hand1 = ax.hist(b_all[filt], histtype="step")[-1][0]
            hand2 = ax.axvline(x=b_median[filt], color="k")
            # decorate
            ax.annotate(
                filt, (0.95, 0.95), xycoords="axes fraction",
                ha="right", va="top")
        fig_add_xlabel(axes, "Smoothing $b$")
        # add legend and fix layout
        fig.legend(
            handles=[hand1, hand2], labels=["fields", "median (fiducial)"],
            ncol=2, frameon=False)
        fig.tight_layout()
        fig.subplots_adjust(top=0.94)
        self.add_fig(fig)

    def plot_magnitudes(self, data, stats, b):
        # get the required data
        fluxes = self.config.get_fluxes(data)
        errors = self.config.get_errors(data)
        fields = self.config.get_fields(data)
        magnitudes = self.config.get_magnitudes(data)
        # compute classic magnitudes as fallback from fluxes and zeropoints
        if magnitudes is None:
            magnitudes = {}
            for filt in self.config.filters:
                zeropoint = stats[Keys.zp].loc[fields]
                zeropoint.index = fluxes[filt].index
                magnitudes = compute_classic_magnitude(fluxes[filt], zeropoint)
        # make figure
        style = {
            "edgecolor": "none", "s": 3, "marker": ".",
            "alpha": 0.3, "rasterized": True}
        xlims = [-3, 33]
        ylims = [0.0, 90.0]
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
            mag_theory = compute_prediction(SN_theory, b[filt])
            hand_theory = ax.plot(
                SN_theory, mag_theory, color="k", lw=0.7, ls="--", zorder=2)[0]
            ylims = [  # update limits from y=mag(-b) to y=max(min)
                max(ylims[0], hyperbolic.compute_magnitude(-b[filt], b[filt])),
                min(ylims[1], mag_theory.min())]
            # add classical magnitudes
            mags = sparse.apply(magnitudes[filt], mask=is_good)
            hand1 = ax.scatter(SN, mags, c="C0", zorder=-1, **style)
            # add hyperbolic magnitudes
            key_mag = self.config.outname[filt]
            mags = sparse.apply(data[key_mag], mask=is_good)
            hand2 = ax.scatter(SN, mags, c="C3", **style)
            # decorate
            ax.axhline(
                y=hyperbolic.compute_magnitude(0.0, b[filt]),
                color="k", lw=0.4)
            ax.axvline(x=0.0, color="k", lw=0.25)
            ax.annotate(
                filt, (0.1, 0.95), xycoords="axes fraction",
                ha="left", va="top")
            ax.set_xlim(*xlims)
            ax.set_ylim(*ylims)
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

    def plot_colours(self, data, stats):
        # get the required data
        fluxes = self.config.get_fluxes(data)
        errors = self.config.get_errors(data)
        fields = self.config.get_fields(data)
        magnitudes = self.config.get_magnitudes(data)
        # compute classic magnitudes as fallback from fluxes and zeropoints
        if magnitudes is None:
            magnitudes = {}
            for filt in self.config.filters:
                zeropoint = stats[Keys.zp].loc[fields]
                zeropoint.index = fluxes[filt].index
                magnitudes = compute_classic_magnitude(fluxes[filt], zeropoint)
        # make figure
        bins = np.linspace(-1.5, 2.5, 50)
        fig, axes = self.make_figure(n_plot_offset=-1)
        for i, (filt1, filt2) in enumerate(
                zip(self.config.filters[:-1], self.config.filters[1:])):
            # select all observed objects and create a sparse sampling
            ax = axes.flatten()[i]
            is_good = (errors[filt1] > 0.0) & (errors[filt2] > 0.0)
            # add classical colours
            idx_sort = np.argsort(magnitudes[filt1][is_good])
            mag1 = magnitudes[filt1][is_good].to_numpy()[idx_sort]
            mag2 = magnitudes[filt2][is_good].to_numpy()[idx_sort]
            colours = (mag1 - mag2)
            hand1a = ax.hist(
                colours[:len(colours)//4], bins,
                color="C0", alpha=0.3, zorder=-1)[-1][0]
            hand1b = ax.hist(
                colours[-len(colours)//4:], bins,
                color="C0", histtype="step")[-1][0]
            # add hyperbolic colours
            key_mag1 = self.config.outname[filt1]
            key_mag2 = self.config.outname[filt2]
            idx_sort = np.argsort(data[key_mag1][is_good])
            mag1 = data[key_mag1][is_good].to_numpy()[idx_sort]
            mag2 = data[key_mag2][is_good].to_numpy()[idx_sort]
            colours = (mag1 - mag2)
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
            loc="center", #prop={"size": 10}
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
