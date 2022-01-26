class Keys:
    filter = "filter"
    field = "field ID"
    b = "b relative"
    ref_flux = "ref. flux"

from .magnitudes import (
    pogson, fill_missing_stats, estimate_b,
    ref_flux_from_zp, zp_from_ref_flux, estimate_zp,
    compute_magnitude, compute_magnitude_error)
