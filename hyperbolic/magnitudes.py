import warnings

import numpy as np
import pandas as pd

from . import Keys
from .logger import logger


pogson = 2.5 * np.log10(np.e)


def fields_to_source(per_field_data, fields, index=None):
    per_source_data = per_field_data.loc[fields]
    if index is not None:
        per_source_data.index = index
    return per_source_data


def ref_flux_from_zp(m_0):
    return np.exp(m_0 / pogson)


def zp_from_ref_flux(flux):
    return pogson * np.log(flux)


def estimate_zp(magnitude, flux):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return magnitude + pogson * np.log(flux)


def convert_flux(flux, zeropoint, target_zp):
    return flux * np.exp((target_zp - zeropoint) / pogson)


def estimate_b(zp, flux_error):
    return np.sqrt(pogson) * np.exp(-zp / pogson) * flux_error


def compute_classic_magnitude(flux, zeropoint, fill=None):
    if fill is None:
        fill = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mag = np.where(flux > 0.0, -pogson * np.log(flux) + zeropoint, fill)
    return mag


def compute_classic_magnitude_error(flux, flux_error, fill=None):
    if fill is None:
        fill = np.nan
    mag_err = np.where(flux > 0.0, pogson * flux_error / flux, fill)
    return mag_err


def compute_magnitude(norm_flux, b):
    return -pogson * (np.arcsinh(0.5 * norm_flux / b) + np.log(b))


def compute_magnitude_error(norm_flux, b, norm_flux_err):
    num = pogson * norm_flux_err
    denom = np.sqrt(norm_flux**2 + 4 * b**2)
    return num / denom


def compute_flux_stats(
        flux, error, fields, magnitudes=None, zeropoint=None, is_good=None):
    if magnitudes is None and zeropoint is None:
        raise ValueError("either magnitudes or a zeropoint must be provided")
    if is_good is None:
        is_good = np.ones(len(error), dtype="bool")
    df = pd.DataFrame({
        Keys.field: fields,
        Keys.flux_err: np.where(is_good, error, np.nan)})
    if zeropoint is None:
        df[Keys.zp] = estimate_zp(
            np.where(is_good, magnitudes, np.nan),
            np.where(is_good, flux, np.nan))
    else:
        df[Keys.zp] = zeropoint
    stats = df.groupby(Keys.field).agg(np.nanmedian)
    stats.index.name = Keys.field
    stats[Keys.ref_flux] = ref_flux_from_zp(stats[Keys.zp])
    return stats


def fill_missing_stats(stats):
    fixed_stats = stats.copy()
    mean = stats.mean()
    row_good = np.isfinite(stats).all(axis=1).to_numpy()
    if not np.all(row_good):
        fields = ", ".join(str(r) for r in stats.index[~row_good])
        logger.warn(
            f"replaced statistics of fields with global mean: {fields}")
        for col in stats.columns:
            fixed_stats[col] = np.where(row_good, stats[col], mean[col])
    return fixed_stats


def compute_flux_error_target(adapt_stats, target_zp=0.0):
    adapt_errors = convert_flux(
        adapt_stats[Keys.flux_err],
        adapt_stats[Keys.zp], target_zp)
    adapt_errors = adapt_errors.groupby(
        Keys.filter).agg(np.nanmedian)
    return adapt_errors


def adapt_flux(flux, error, stats, adapt_error, fields, is_good=None):
    if is_good is None:
        is_good = np.ones(len(error), dtype="bool")
    # convert the targeted flux error to the correct zeropoint (originally ZP=0.0)
    adapt_error = adapt_error * np.exp(stats[Keys.zp] / pogson)
    # compute the additional noise that must be added to the flux
    add_variance = adapt_error**2 - stats[Keys.flux_err]**2
    is_modified = add_variance > 0.0  # find fields that cannot be adapted
    add_variance[:] = np.where(is_modified, add_variance, 0.0)
    if not is_modified.all():
        n_unchanged = np.count_nonzero(~is_modified)
        logger.warn(f"flux in {n_unchanged} fields too large to be adapted")
    # compute standard deviation for each object
    add_sigma = fields_to_source(
        np.sqrt(add_variance), fields, index=flux.index)
    # updated the flux and flux error
    new_flux = np.where(
        is_good, np.random.normal(flux, add_sigma), flux)
    new_error = np.where(
        is_good, np.sqrt(error**2 + add_sigma**2), error)
    return new_flux, new_error
