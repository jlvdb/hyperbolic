import argparse
import json

import astropandas as apd
import numpy as np
import pandas as pd

from . import Keys
from .logger import logger


error_suffix = "_err"


class DumpAction(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        default = json.dumps({
            "filter_name (e.g. 'r')": {
                "outname": "# output column name (suffix appended for "
                          f"errors: {error_suffix})",
                "flux": "# name of colum with flux",
                "error": "# name of column with flux error",
                "magnitude": "# name of column with magnitudes (value ignored "
                             "if --zeropoint is given)",
            }
        }, indent=4)
        print(default)
        parser.exit()


class LoadConfig:

    def __init__(self, args):
        self.verbose = args.verbose
        # get the commandline arguments
        self.infile = args.infile
        self.outfile = args.outfile
        self.fields = args.fields
        # read the configuration file
        logger.info(f"reading configuration from {args.config}")
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
        self.zeropoint = None  # default, overwritten by subclasses

    def load_input(self):
        logger.info(f"reading data from {self.infile}")
        return apd.read_auto(self.infile)

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


class LoadConfigSmooting(LoadConfig):

    def __init__(self, args):
        super().__init__(args)
        # check the additional commandline arguments
        self.zeropoint = args.zeropoint
        if self.zeropoint is not None:
            self.magnitude = None

    def write_output(self, data):
        logger.info(f"writing statistics to {self.outfile}")
        data.to_csv(self.outfile)


class LoadConfigMagnitudes(LoadConfig):

    def __init__(self, args):
        super().__init__(args)
        # check the additional commandline arguments
        self.stats = args.stats
        self.smoothing = args.smoothing
        self.b_global = args.b_global
        self.plot = args.plot

    def load_stats(self):
        logger.info(f"reading statistics from {self.stats}")
        return pd.read_csv(self.stats, index_col=[Keys.filter, Keys.field])

    def load_smoothing(self):
        fpath = self.smoothing if self.smoothing is not None else self.stats
        logger.info(f"reading statistics from {fpath}")
        return pd.read_csv(fpath, index_col=[Keys.filter, Keys.field])

    @staticmethod
    def KiDS_aware_error_colname(mag_col_name):
        if mag_col_name.startswith("HMAG_"):
            error_colname = f"HMAGERR_{mag_col_name[5:]}"
        else:
            error_colname = mag_col_name + error_suffix
        return error_colname

    def write_output(self, data):
        logger.info(f"writing table data to {self.outfile}")
        apd.to_fits(data, self.outfile)
