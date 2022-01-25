import argparse
import json
import warnings

import numpy as np


pogson = 2.5 * np.log10(np.e)


class DumpAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        default = json.dumps({
            "filter_name (e.g. 'r')": {
                "outname": "# output column name (suffix for errors: _err)",
                "flux": "# name of colum with flux",
                "error": "# name of column with flux error",
                "magnitude": "# name of column with magnitudes (ignored if "
                             "--zeropoint is given)",
            }
        }, indent=4)
        print(default)
        parser.exit()


class LoadConfig:

    def __init__(self, args):
        # check the commandline arguments
        self.infile = args.infile
        self.outfile = args.outfile
        self.fields = args.fields
        self.zeropoint = args.zeropoint
        # read the configuration file
        print(f"reading: {args.config}")
        with open(args.config) as f:
            data = json.load(f)
            self.filters = list(data.keys())
            self.outname = {
                fname: config["outname"] for fname, config in data.items()}
            self.flux = {
                fname: config["flux"] for fname, config in data.items()}
            self.error = {
                fname: config["error"] for fname, config in data.items()}
            self.magnitude = {
                fname: config["magnitude"] for fname, config in data.items()}
        if self.zeropoint is not None:
            self.magnitude = None

    def get_fields(self, df):
        try:
            if self.fields is None:
                return np.zeros(len(df), dtype=np.float)
            else:
                return df[self.fields].to_numpy()
        except KeyError as e:
            raise KeyError(f"fields column {e} not found")

    def get_fluxes(self, df):
        try:
            return {key: df[val] for key, val in self.flux.items()}
        except KeyError as e:
            raise KeyError(f"flux column {e} not found")

    def get_errors(self, df):
        try:
            return {key: df[val] for key, val in self.error.items()}
        except KeyError as e:
            raise KeyError(f"flux error column {e} not found")

    def get_magnitudes(self, df):
        try:
            if self.zeropoint is None:
                return {key: df[val] for key, val in self.magnitude.items()}
            else:
                return None
        except KeyError as e:
            raise KeyError(f"magnitude column {e} not found")


def ref_flux_from_zp(m_0):
    return np.exp(m_0 / pogson)


def zp_from_ref_flux(flux):
    return pogson * np.log(flux)


def estimate_zp(magnitude, flux):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return magnitude + pogson * np.log(flux)


def estimate_b(zp, flux_error):
    return np.sqrt(pogson) * np.exp(-zp / pogson) * flux_error
