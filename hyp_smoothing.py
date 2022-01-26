#!/usr/bin/env python3
import argparse

import numpy as np
import pandas as pd

import hyperbolic
import hyperbolic.config


parser = argparse.ArgumentParser(
    description="...",
    add_help=False)
parser.add_argument(
    "infile", metavar="path",
    help="input FITS file")
parser.add_argument(
    "outfile", metavar="path",
    help="output path for data statistics")

parser.add_argument(
    "-c", "--config", required=True,
    help="JSON configuration file that specifies in- and output data "
         "(see --dump)")
parser.add_argument(
    "-f", "--fields",
    help="column name that can uniquely indentify pointings")
parser.add_argument(
    "--zeropoint", type=float,
    help="use this fixed magnitude zeropoint")

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
    config = hyperbolic.config.LoadConfigSmooting(parser.parse_args())
    # load the input data
    data = config.load_input()
    # get data columns
    fields = config.get_fields(data)
    fluxes = config.get_fluxes(data)
    errors = config.get_errors(data)
    magnitudes = config.get_magnitudes(data)

    all_stats = []
    for filt in config.filters:
        hyperbolic.logger.logger.info(f"processing filter {filt}")
        # mask to valid observations
        is_good = errors[filt] > 0.0

        # compute zeropoint and median flux error
        df = pd.DataFrame({
            hyperbolic.Keys.field: fields,
            "flux error": np.where(is_good, errors[filt], np.nan)})
        if config.zeropoint is None:
            df[hyperbolic.Keys.zp] = hyperbolic.estimate_zp(
                np.where(is_good, magnitudes[filt], np.nan),
                np.where(is_good, fluxes[filt], np.nan))
        else:
            df[hyperbolic.Keys.zp] = config.zeropoint
        stats = df.groupby(hyperbolic.Keys.field).agg(np.nanmedian)
        stats.index.name = hyperbolic.Keys.field
        stats[hyperbolic.Keys.ref_flux] = \
            hyperbolic.ref_flux_from_zp(stats[hyperbolic.Keys.zp])

        # compute b
        stats[hyperbolic.Keys.b] = hyperbolic.estimate_b(
            stats[hyperbolic.Keys.zp], stats["flux error"])
        stats[hyperbolic.Keys.b_abs] = \
            stats["ref. flux"] * stats[hyperbolic.Keys.b]

        # collect statistics
        stats[hyperbolic.Keys.filter] = filt
        stats = stats.reset_index().set_index([
            hyperbolic.Keys.filter, hyperbolic.Keys.field])
        if config.verbose:
            print(stats)
        all_stats.append(stats)

    config.write_output(pd.concat(all_stats))
