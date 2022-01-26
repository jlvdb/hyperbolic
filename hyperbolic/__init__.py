class Keys:
    filter = "filter"
    field = "field ID"
    b = "b relative"
    b_abs = "b absolute"
    ref_flux = "ref. flux"
    zp = "zeropoint"

from .magnitudes import (
    pogson, fill_missing_stats, estimate_b,
    ref_flux_from_zp, zp_from_ref_flux, estimate_zp,
    compute_classic_magnitude, compute_magnitude, compute_magnitude_error)
