import warnings

import numpy as np

from .logger import logger


pogson = 2.5 * np.log10(np.e)


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


def compute_magnitude(norm_flux, b):
    return -pogson * (np.arcsinh(0.5 * norm_flux / b) + np.log(b))


def compute_magnitude_error(norm_flux, b, norm_flux_err):
    num = pogson * norm_flux_err
    denom = np.sqrt(norm_flux**2 + 4 * b**2)
    return num / denom
