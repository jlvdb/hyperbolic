import warnings
from typing import Optional

import numpy as np
import numpy.typing as npt
import pandas as pd

from . import Keys
from .logger import logger


pogson = 2.5 * np.log10(np.e)


def fields_to_source(
    per_field_data: pd.Series,
    fields: npt.NDArray,
    index: Optional[pd.Index] = None
) -> pd.Series:
    """
    Map field-wide statistics back to each individual source based on the field membership.

    Parameters:
    -----------
    per_field_data : pd.Series
        Series with field identifier as index and statistics as values.
    fields : array-like
        List of field identifiers that defines the field membership of each source.
    index : pd.Index
        Optional index to assign to the output.
    
    Returns:
    --------
    per_source_data : pd.Series
        Field-wide statistics mapped to each source.
    """
    per_source_data = per_field_data.loc[fields]
    if index is not None:
        per_source_data.index = index
    return per_source_data


def ref_flux_from_zp(m_0: npt.ArrayLike) -> npt.ArrayLike:
    """
    Convert the photometric zeropoint to the reference flux.
    """
    return np.exp(m_0 / pogson)


def zp_from_ref_flux(flux: npt.ArrayLike) -> npt.ArrayLike:
    """
    Convert the reference flux to the photometric zeropoint
    """
    return pogson * np.log(flux)


def estimate_zp(
    magnitude: npt.ArrayLike,
    flux: npt.ArrayLike
) -> npt.ArrayLike:
    """
    Estimate the photometric zeropoint from a set of magnitudes and corresponding fluxes.

    Parameters:
    -----------
    magnitude : array-like
        Magnitudes of objects.
    flux : array-like
        Fluxes of objects.
    
    Returns:
    --------
    zp : array-like
        Photometric zeropoint per object.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return magnitude + pogson * np.log(flux)


def convert_flux(
    flux: npt.ArrayLike,
    zeropoint: npt.ArrayLike,
    target_zp: npt.ArrayLike
) -> npt.ArrayLike:
    """
    Convert a flux measurement obtained with a given photometric zeropoint to a flux system with a
    different photometric zeropoint.

    Parameters:
    -----------
    flux : array-like
        Input fluxes of objects.
    zeropoint : array-like
        Photometric zeropoint (per source) corresponding to the input fluxes.
    target_zp : array-like
        Photometric zeropoint (per source) for the output fluxes.
    
    Returns:
    --------
    flux : array-like
        Output fluxes for the new photometric zeropoint.
    """
    return flux * np.exp((target_zp - zeropoint) / pogson)


def estimate_b(
    zp: npt.ArrayLike,
    flux_error: npt.ArrayLike
) -> npt.ArrayLike:
    """
    Estimate the smoothing factor for hyperbolic magnitudes in flux units for a given photometric
    zeropoint and flux errors.

    Parameters:
    -----------
    zp : array-like
        Photometric zeropoint (per source) for the output fluxes.
    flux_error : array-like
        Flux error (per object).

    Returns:
    --------
    b : array-like
        Estimate of the smoothing factor (per object).
    """
    return np.sqrt(pogson) * np.exp(-zp / pogson) * flux_error


def compute_classic_magnitude(
    flux: npt.ArrayLike,
    zeropoint: npt.ArrayLike,
    fill: Optional[float] = None
) -> npt.ArrayLike:
    """
    Compute the classical magnitude for a given flux measurement and photometric zeropoint.

    Parameters:
    -----------
    flux : array-like
        Input fluxes of objects.
    zeropoint : array-like
        Photometric zeropoint (per object).
    fill : array-like
        Fill value to use for objects with flux <= 0, defaults to NaN.
    
    Returns:
    --------
    mag : array-like
        Classical magnitude for objects.
    """
    if fill is None:
        fill = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mag = np.where(flux > 0.0, -pogson * np.log(flux) + zeropoint, fill)
    return mag


def compute_classic_magnitude_error(
    flux: npt.ArrayLike,
    flux_error: npt.ArrayLike,
    fill: Optional[float] = None
) -> npt.ArrayLike:
    """
    Compute the error of the classical magnitude for a given flux measurement and error.

    Parameters:
    -----------
    flux : array-like
        Input fluxes of objects.
    flux_error : array-like
        Input flux errors of objects.
    fill : array-like
        Fill value to use for objects with flux <= 0, defaults to NaN.
    
    Returns:
    --------
    mag : array-like
        Classical magnitude error for objects.
    """
    if fill is None:
        fill = np.nan
    mag_err = np.where(flux > 0.0, pogson * flux_error / flux, fill)
    return mag_err


def compute_magnitude(
    norm_flux: npt.ArrayLike,
    b: npt.ArrayLike
) -> npt.ArrayLike:
    """
    Compute the hyperbolic magnitude for a given flux measurement and smoothing parameter b.

    Parameters:
    -----------
    norm_flux : array-like
        Normalised fluxes (flux divided by reference flux) of objects.
    b : array-like
        Smoothing parameter (per object).
    
    Returns:
    --------
    mag : array-like
        Hyperbolic magnitude for objects.
    """
    return -pogson * (np.arcsinh(0.5 * norm_flux / b) + np.log(b))


def compute_magnitude_error(
    norm_flux: npt.ArrayLike,
    b: npt.ArrayLike,
    norm_flux_err: npt.ArrayLike
) -> npt.ArrayLike:
    """
    Compute the error of the hyperbolic magnitude for a given flux measurement and error, and
    smoothing parameter.

    Parameters:
    -----------
    norm_flux : array-like
        Normalised fluxes (flux divided by reference flux) of objects.
    b : array-like
        Smoothing parameter (per object).
    norm_flux_err : array-like
        Error of the normalised fluxes of objects.
    
    Returns:
    --------
    error : array-like
        Hyperbolic magnitude error for objects.
    """
    num = pogson * norm_flux_err
    denom = np.sqrt(norm_flux**2 + 4 * b**2)
    return num / denom


def compute_flux_stats(
    flux: npt.NDArray,
    error: npt.NDArray,
    fields: npt.NDArray,
    magnitudes: Optional[npt.NDArray] = None,
    zeropoint: Optional[npt.ArrayLike] = None,
    is_good: Optional[npt.NDArray] = None
) -> pd.DataFrame:
    """
    Compute the median flux error, photometric zeropoint, and reference flux per field (e.g. spatial
    regions, telescope pointings, etc.). The field identifiers can be set to a single value to
    compute the statistics globally. 

    Parameters:
    -----------
    flux : array-like
        Input fluxes of objects.
    error : array-like
        Input flux errors of objects.
    fields : array-like
        List of field identifiers that defines the field membership of each object.
    magnitudes : array-like
        Magnitudes of objects, used optionally to compute photometric zeropoints per object.
    zeropoint : array-like
        If no magnitudes are provided, use this photometric zeropoint (per object).
    is_good : array-like of bool
        Mask to exclude objects from computation by setting their mask value to False.

    Returns:
    --------
    stats : pd.DataFrame
        Dataframe with field identifier as index. The columns list the median flux error, the
        photometric zeropoint, and reference flux per field.
    """
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


def fill_missing_stats(stats: pd.DataFrame) -> pd.DataFrame:
    """
    Substitute missing values of per-field statistics (see `compute_flux_stats`) with the mean over
    all fields. Returns a new copy of the input statistics.
    """
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
