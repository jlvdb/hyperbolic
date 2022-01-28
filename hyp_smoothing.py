#!/usr/bin/env python3
import argparse
from types import new_class

import numpy as np
import pandas as pd

import hyperbolic
import hyperbolic.config


parser = argparse.ArgumentParser(
    description="Compute the statistics of a photometric data catalogue that "
                "are needed to compute the smoothing factor for the "
                "hyperbolic magnitudes.",
    add_help=False)
parser.add_argument(
    "infile", metavar="infile",
    help="input FITS file")
parser.add_argument(
    "--hdu", type=int, default=1,
    help="FITS HDU index to read (default: %(default)s)")
parser.add_argument(
    "outfile", metavar="statfile",
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

adpat = parser.add_argument_group(
    "flux adaptation",
    description="Adjust the flux noise level in the input data to match that "
                "of an external data catalogue before measuring statistics.")
adpat.add_argument(
    "-a", "--adapt",
    help="statistics file genereated with this script for the external data "
         "catalogue")
adpat.add_argument(
    "--suffix", default="_adapt",
    help="suffix that is added to the names of data columns with the adapted "
         "fluxes and flux errors (default: %(default)s)")
adpat.add_argument(
    "--adapt-file",
    help="file path where the output FITS file with added adapted fluxes is "
         "written (default derived from input: [infile]_adapted.*)")

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
    adapt_stats = config.load_adapt()  # statistics of the external data
    # get data columns
    fields = config.get_fields(data)
    fluxes = config.get_fluxes(data)
    errors = config.get_errors(data)
    magnitudes = config.get_magnitudes(data)

    # compute median flux error in the external catalogue for fixed ZP=0
    if config.adapt is not None:
        adapt_errors = hyperbolic.compute_flux_error_target(adapt_stats, 0.0)
        if config.verbose:  # print in nicely looking dataframe
            _df = pd.DataFrame({hyperbolic.Keys.flux_err: adapt_errors})
            print(_df.loc[config.filters])  # maintain usual order
        adapt_errors = adapt_errors.to_dict()
    else:
        adapt_errors = {filt: None for filt in config.filters}

    # compute the data statistics required for the smoothing
    all_stats = []
    for filt in config.filters:
        hyperbolic.logger.logger.info(f"processing filter {filt}")
        # mask to valid observations
        is_good = errors[filt] > 0.0

        # compute the additional noise required to match the external data
        if config.adapt is not None:
            # compute zeropoint and median flux error
            stats = hyperbolic.compute_flux_stats(
                fluxes[filt], errors[filt], fields, magnitudes[filt],
                config.zeropoint, is_good)
            stats = hyperbolic.fill_missing_stats(stats)
            zeropoint = hyperbolic.fields_to_source(
                stats[hyperbolic.Keys.zp], fields, index=fluxes[filt].index)
            # updated the flux and flux error and add to input data table
            fluxes[filt], errors[filt] = hyperbolic.adapt_flux(
                fluxes[filt], errors[filt], stats, adapt_errors[filt],
                fields, is_good)
            config.add_column_and_update(data, fluxes[filt], filt, "flux")
            config.add_column_and_update(data, errors[filt], filt, "error")
            # compute the classical adapted magnitudes
            if magnitudes[filt] is None:
                fill_mag, fill_err = None, None
            else:
                if np.isnan(magnitudes[filt]).any():
                    fill_mag, fill_err = None, None
                else:
                    fill_mag = magnitudes[filt].max()
                    fill_err = hyperbolic.compute_classic_magnitude(
                        errors[filt], zeropoint)
            mag_err = hyperbolic.compute_classic_magnitude_error(
                fluxes[filt], errors[filt], fill=fill_err)
            # update the internal values
            magnitudes[filt] = hyperbolic.compute_classic_magnitude(
                fluxes[filt], zeropoint, fill=fill_mag)
            # include them in the input data
            data[f"mag_classical_adapt_{filt}"] = \
                magnitudes[filt].astype(np.float32)
            data[f"magerr_classical_adapt_{filt}"] = mag_err.astype(np.float32)

        # compute zeropoint and median flux error
        stats = hyperbolic.compute_flux_stats(
            fluxes[filt], errors[filt], fields, magnitudes[filt],
            config.zeropoint, is_good)
        # compute b
        stats[hyperbolic.Keys.b] = hyperbolic.estimate_b(
            stats[hyperbolic.Keys.zp], stats[hyperbolic.Keys.flux_err])
        stats[hyperbolic.Keys.b_abs] = \
            stats[hyperbolic.Keys.ref_flux] * stats[hyperbolic.Keys.b]

        # collect statistics
        stats[hyperbolic.Keys.filter] = filt
        stats = stats.reset_index().set_index([
            hyperbolic.Keys.filter, hyperbolic.Keys.field])
        if config.verbose:
            print(stats)
        all_stats.append(stats)

    if config.adapt is not None:
        config.write_output(data)
    config.write_stats(pd.concat(all_stats))
