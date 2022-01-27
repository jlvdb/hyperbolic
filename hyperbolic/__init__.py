class Keys:
    filter = "filter"
    field = "field ID"
    b = "b relative"
    b_abs = "b absolute"
    ref_flux = "ref. flux"
    zp = "zeropoint"
    flux_err = "flux error"

from .magnitudes import (
    pogson, estimate_b,
    ref_flux_from_zp, zp_from_ref_flux, estimate_zp, convert_flux,
    compute_classic_magnitude, compute_magnitude, compute_magnitude_error,
    compute_flux_stats, fill_missing_stats,
    compute_flux_error_target, adapt_flux)
