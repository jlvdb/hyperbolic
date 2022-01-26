# hyperbolic

Implementation of hyperbolic magnitudes (Lupton et al. 1999).

## Installation

This code is written for `python3.6+`. All dependencies are listed in the
`requirements.txt`.

## Usage

The code is split in two parts. The first script, `hyp_smoothing.py`, reads an
input FITS table with photometric data (flux and flux errors and optional
classical magnitudes) and calculates statistics that determine the smoothing
paramter for the hyperbolical magnitudes, such as the photometric zeropoint
and the median flux uncertainty (see implementation section below). The input
data is specified in the configuration file. An empty configration file can be
generated with

   hyp_smoothing.py -d

An example usage might look like this:

```shell
hyp_smoothing.py -v \
   --field pointing_identifier_column_name \
   -c configration.json \
   input_catalogue.fits \
   field_statistics.csv
```

The second part of the code, `hyp_magnitudes.py`, computes the actual value for
the smoothing parameter and computes hyperbolic magnitudes consistently. The
code corrects automatically variations in the photometric zeropoint (if
magnitudes are provided in the configuration file) based on the
`field_statistics.csv` computed in the first step. An example usage might look
like this:

```shell
hyp_magnitudes.py -v \
   --field pointing_identifier_column_name \
   -c configration.json \
   -s field_statistics.csv \
   --plot \
   input_catalogue.fits \
   output_catalogue_with_hyperb_mags.fits
```

It is also possible to apply the smoothing parameter from an external data
cataloge. In this case, the field statistics of the external data should be
computed with `hyp_smoothing.py`. The output .csv file can be provided to
`hyp_magnitudes.py` by adding the `--smoothing` parameter to the call above.

## Implementation

Hyperbolical magnitudes approximate the classical magnitudes at high signal-to-
noise, but have linear behaviour around flux values of f=0 and are therefore
well defined if f<0. Lupton et al. (1999) suggest to define hyperbolical
magnitudes in terms of the normalised flux x=f/f0 as

    mu = a * [arcsinh(0.5*x/b) + log(b)] ,

with a=2.5log10(e). In this parameterisation f0 is the flux of an object with
magnitude zero and m0=a*log(f0) is the corresponding photometric zeropoint.

The free parameter b is a smoothing factor that determines the transition
between linear and logarithmic behaviour of mu. Lupton et al. (1999) show that
the optimal value for the smoothing parameter is

    b = sqrt(a) Dx ,

i.e. is determined by the variance of the normalised flux.

When applied to observational data, the hyperbolic magnitudes can be written as

    mu = a * [arcsinh(0.5*f/b') + log(b')] + m0 ,

where f is the measured flux and

    b' = f0 * b = sqrt(a) * Df

is determined from the variance of the measured fluxes. In this formulation, b
depends on the photometric zeropoint.

### Estimating the smoothing parameter

For a given set of flux measurements the hyperbolic magnitudes can be computed
once an appropriate smoothing parameter is chosen. This determines b globally
from the complete input data and for each photometric filter specified in the
configuration file:

1. It computes the photometric zeropoint m0 in each telescope pointing by
   comparing the individual flux measurements fi and magnitudes mi, since the
   latter may include a number of corrections, such as extinction or stellar
   locus regressions:

       m0 = median(mi + a*log(fi))

2. It computes the smoothing parameter for each pointing from the zeropoint
   according to
   
       b = b' / f0 = sqrt(a) * Df / f0 = sqrt(a) * e^(-m0/a) * Df ,

   where Df=median(fi) is the median of the measured flux errors.
3. It computes the global value for b in each filter by taking the median of
   all pointings.

### Computation

The code computes the hyperbolical magnitudes from their normalised flux
x=f/f0 with uncertainty Dx=Df/f0 based on these global values for b. In each
pointing it calcualtes the flux f0=e^(m0/a) from the zeropoint m0 (see above)
to compensate variations in the observing conditions.

## Maintainers

[Jan Luca van den Busch](jlvdb@astro.rub.de),
Ruhr-University Bochum, Astronomical Institute.
