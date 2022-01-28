#!/usr/bin/env python3
import argparse

import numpy as np
import pandas as pd

import hyperbolic
import hyperbolic.config
import hyperbolic.plots


parser = argparse.ArgumentParser(
    description="Compute the hyperbolic magnitudes and add them into the "
                "input data catalogue.",
    add_help=False)
parser.add_argument(
    "infile", metavar="infile",
    help="input FITS file")
parser.add_argument(
    "--hdu", type=int, default=1,
    help="FITS HDU index to read (default: %(default)s)")
parser.add_argument(
    "outfile", metavar="outfile",
    help="output path for FITS file with added hyperbolic magnitudes")

parser.add_argument(
    "-c", "--config", required=True,
    help="JSON configuration file that specifies in- and output data "
         "(see --dump)")
parser.add_argument(
    "-s", "--stats", required=True,
    help="statistics file genereated with 'hyp_smoothing.py'")
parser.add_argument(
    "--smoothing",
    help="external statistics file genereated with 'hyp_smoothing.py' that is "
         "used to define the smoothing parameter b")
parser.add_argument(
    "-f", "--fields",
    help="column name that can uniquely indentify pointings")
parser.add_argument(
    "--b-global", action='store_true',
    help="compute the smoothing paramter globally for all filters")
parser.add_argument(
    "--plot", action='store_true',
    help="add a PDF file with summary plots alongside the output file")

group = parser.add_argument_group("help")
group.add_argument(
    "-h", "--help", action="help",
    help="show this help message and exit")
group.add_argument(
    "-d", "--dump", nargs=0, action=hyperbolic.config.DumpAction,
    help="dump an empty configuration file, 'filter_name' can be repeated")
parser.add_argument(
    "-v", "--verbose", action="store_true",
    help="show statistic summary per filter")


if __name__ == "__main__":
    config = hyperbolic.config.LoadConfigMagnitudes(parser.parse_args())
    # load the input data
    data = config.load_input()
    all_stats = config.load_stats()  # for the internal zeropoint corrections
    smooth_stats = config.load_smoothing()
    # get data columns
    fields = config.get_fields(data)
    fluxes = config.get_fluxes(data)
    errors = config.get_errors(data)

    # compute smoothing factor (from external data if --smoothing is provided)
    b = smooth_stats[hyperbolic.Keys.b].groupby(
        hyperbolic.Keys.filter).agg(np.nanmedian)
    if config.b_global:
        b[:] = np.median(b)
    if config.verbose:  # print in nicely looking dataframe
        _df = pd.DataFrame({hyperbolic.Keys.b: b})
        print(_df.loc[config.filters])  # maintain usual order
    b = b.to_dict()

    for filt in config.filters:
        hyperbolic.logger.logger.info(f"processing filter {filt}")
        stats = hyperbolic.fill_missing_stats(all_stats.loc[filt])
        # mask to valid observations
        is_good = errors[filt] > 0.0

        # compute normalised flux
        ref_flux = hyperbolic.fields_to_source(
            stats[hyperbolic.Keys.ref_flux], fields, index=fluxes[filt].index)
        norm_flux = fluxes[filt] / ref_flux
        norm_flux_err = errors[filt] / ref_flux

        # compute the hyperbolic magnitudes
        hyp_mag = hyperbolic.compute_magnitude(
            norm_flux, b[filt])
        hyp_mag_err = hyperbolic.compute_magnitude_error(
            norm_flux, b[filt], norm_flux_err)

        # add data to catalogue
        key_mag = config.outname[filt]
        key_mag_err = config.KiDS_aware_error_colname(key_mag)
        data[key_mag] = np.where(
            is_good, hyp_mag, -99.0).astype(np.float32)
        data[key_mag_err] = np.where(
            is_good, hyp_mag_err, -99.0).astype(np.float32)

    if config.plot:
        with hyperbolic.plots.Plotter(config, data, all_stats, b) as plotter:
            if config.fields is not None:
                plotter.plot_b(smooth_stats)
            plotter.plot_magnitudes()
            plotter.plot_magnitude_distribution()
            plotter.plot_colour_distribution()
            plotter.plot_magdiff()
    config.write_output(data)
