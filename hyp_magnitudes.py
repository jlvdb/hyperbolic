#!/usr/bin/env python3
import argparse

import astropandas as apd
import numpy as np
import pandas as pd

import hyperbolic
import hyperbolic.config


def KiDS_aware_error_colname(mag_col_name):
    if mag_col_name.startswith("HMAG_"):
        error_colname = f"HMAGERR_{mag_col_name[5:]}"
    else:
        error_colname = mag_col_name + hyperbolic.config.error_suffix
    return error_colname


parser = argparse.ArgumentParser(
    description="...",
    add_help=False)
parser.add_argument(
    "infile", metavar="path",
    help="input FITS file")
parser.add_argument(
    "outfile", metavar="path",
    help="output path for FITS file with added hyperbolic magnitudes")

parser.add_argument(
    "-c", "--config", required=True,
    help="JSON configuration file that specifies in- and output data "
         "(see --dump)")
parser.add_argument(
    "-s", "--stats", required=True,
    help="statistics file genereated with 'hyp_smoothing.py'")
parser.add_argument(
    "-f", "--fields",
    help="column name that can uniquely indentify pointings")
parser.add_argument(
    "--b-global", action='store_true',
    help="compute the smoothing paramter globally for all filters")

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
    all_stats = config.load_stats()
    # get data columns
    fields = config.get_fields(data)
    fluxes = config.get_fluxes(data)
    errors = config.get_errors(data)
    magnitudes = config.get_magnitudes(data)

    # compute smoothing factor
    global_stats = all_stats.groupby(
        hyperbolic.Keys.filter).agg(np.nanmedian)
    global_stats = global_stats.loc[config.filters]  # maintain order (print)
    b = global_stats[[hyperbolic.Keys.b]].copy()
    if config.b_global:
        b[hyperbolic.Keys.b] = np.median(b[hyperbolic.Keys.b])
    if config.verbose:
        print(b)
    b = b.to_dict()[hyperbolic.Keys.b]

    for filt in config.filters:
        hyperbolic.logger.logger.info(f"processing filter {filt}")
        stats = hyperbolic.fill_missing_stats(all_stats.loc[filt])
        # mask to valid observations
        is_good = errors[filt] > 0.0

        # compute normalised flux
        ref_flux = stats[hyperbolic.Keys.ref_flux].loc[fields]
        ref_flux.index = fluxes[filt].index
        norm_flux = fluxes[filt] / ref_flux
        norm_flux_err = errors[filt] / ref_flux

        # compute the hyperbolic magnitudes
        hyp_mag = hyperbolic.compute_magnitude(
            norm_flux, b[filt])
        hyp_mag_err = hyperbolic.compute_magnitude_error(
            norm_flux, b[filt], norm_flux_err)

        # add data to catalogue
        key_mag = config.outname[filt]
        key_mag_err = KiDS_aware_error_colname(key_mag)
        data[key_mag] = np.where(hyp_mag, is_good, -99.0)
        data[key_mag_err] = np.where(hyp_mag_err, is_good, -99.0)

    config.write_output(data)
