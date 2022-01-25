#!/usr/bin/env python3
import argparse

import astropandas as apd
import numpy as np
import pandas as pd
import hyperbolic

from hyperbolic.mags import DumpAction


parser = argparse.ArgumentParser(
    description="...")
parser.add_argument(
    "infile", metavar="path",
    help="input FITS file")
parser.add_argument(
    "outfile", metavar="path",
    help="output path for data statistics")
parser.add_argument(
    "-f", "--fields",
    help="column name that can uniquely indentify pointings")
parser.add_argument(
    "--zeropoint", type=float,
    help="use this fixed magnitude zeropoint")
parser.add_argument(
    "-c", "--config", required=True,
    help="JSON configuration file that specifies in- and output data "
         "(see --dump)")
parser.add_argument(
    "-d", "--dump", nargs=0, action=DumpAction,
    help="dump an empty configuration file, 'filter_name' can be repeated")


if __name__ == "__main__":
    config = hyperbolic.LoadConfig(parser.parse_args())
    # load the data file
    print(f"reading: {config.infile}")
    data = apd.read_fits(config.infile)
    # get data columns
    fields = config.get_fields(data)
    fluxes = config.get_fluxes(data)
    errors = config.get_errors(data)
    magnitudes = config.get_magnitudes(data)

    all_stats = []
    key_filt, key_field = "filter", "field ID"
    for filt in config.filters:
        # mask to valid observations
        is_good = errors[filt] > 0.0
        # compute zeropoint and median flux error
        df = pd.DataFrame({
            key_field: fields,
            "flux error": np.where(is_good, errors[filt], np.nan)})
        if config.zeropoint is None:
            df["zeropoint"] = hyperbolic.estimate_zp(
                np.where(is_good, magnitudes[filt], np.nan),
                np.where(is_good, fluxes[filt], np.nan))
        else:
            df["zeropoint"] = config.zeropoint
        stats = df.groupby(key_field).agg(np.nanmedian)
        stats.index.name = key_field
        stats["ref. flux"] = hyperbolic.ref_flux_from_zp(stats["zeropoint"])
        # compute b
        stats["b relative"] = hyperbolic.estimate_b(
            stats["zeropoint"], stats["flux error"])
        stats["b absolute"] = stats["ref. flux"] * stats["b relative"]
        # collect statistics
        stats[key_filt] = filt
        stats = stats.reset_index().set_index([key_filt, key_field])
        print(stats)
        all_stats.append(stats)

    # write statistics file
    print(f"writing: {config.outfile}")
    all_stats = pd.concat(all_stats)
    all_stats.to_csv(config.outfile)
