# hyperbolic

Implementation of hyperbolic magnitudes (Lupton et al. 1999).

## Implementation

Hyperbolical magnitudes approximate the classical magnitudes at high signal-to-noise, but have linear behaviour around flux values of $f=0$ and are therefore well defined if $f<0$. Lupton et al. (1999) suggests to define hyperbolical magnitudes in terms of the normalised flux $x = f / f_0$ as
$$ \mu = a \left(\mathrm{arcsinh}\left(\frac{x}{2b}\right) + \log{b}\right) , $$
which has an uncertainty of
$$ \Delta\mu = \frac{a \, \Delta x}{\sqrt{x^2 + 4 b^2}} $$
with $a = 2.5 \log_{10}(e) \approx 1.086$. In this parameterisation $f_0$ is the flux of an object with magnitude zero and $m_0 = a \log{f_0}$ is the corresponding photometric zeropoint.

The free parameter $b$ is a smoothing factor that determines the transition between linear and logarithmic behaviour of $\mu$. Lupton et al. (1999) show that the optimal value for the smoothing parameter is
$$ b = \sqrt{a} \Delta x , $$
i.e. is determined by the variance of the normalised flux.

When applied to observational data, the hyperbolic magnitudes can be written as
$$ \mu = a \left(\mathrm{arcsinh}\left(\frac{f}{2b^\prime}\right) + \log{b^\prime}\right) + m_0 , $$
where $f$ is the measured flux and
$$ b^\prime = f_0 b = \sqrt{a} \Delta f $$
is determined from the variance of the measured fluxes. In this formulation, $b^\prime$ depends on the photometric zeropoint and can be converted according to
$$ b^\prime_1 = b^\prime_2 \exp{\left(\frac{m_{0,1} - m_{0,2}}{a}\right)} , $$
where $b^\prime_1$ and $b^\prime_2$ are the smoothing parameters for observations with zeropoints $m_{0,1}$ and $m_{0,2}$.

### Estimating the smoothing parameter

For a given set of flux measurements the hyperbolic magnitudes can be computed once an appropriate smoothing parameter is chosen. We determine $b$ globally from the whole survey data for each photometric band:

1. We compute the photometric zeropoint $m_0$ in each telescope pointing by comparing the individual flux measurements $f_i$ and magnitudes $m_i$, since the latter already include a number of corrections, such as extinction and stellar locus regressions:
   $$ m_0 = \mathrm{median}(m_i + a \log{f_i}) $$
2. We compute the smoothing parameter for each pointing from the zeropoint according to
   $$ b = \frac{b^\prime}{f_0} = \frac{\sqrt{a} \, \Delta f}{f_0} = \sqrt{a} \, e^{-m_0 / a} \Delta f , $$
   where $\Delta f = \mathrm{median}(f_i)$ is the median of the measured flux errors.
3. We compute the global value for $b$ in each band by taking the median of all pointings.

### Computation

We compute the hyperbolical magnitudes from their normalised flux $x=f/f_0$ with uncertainty $\Delta x = \Delta f / f_0$ based on these global values for $b$. In each pointing we calcualte the flux $f_0 = e^{m_0 / a}$ from the zero-point $m_0$ (see above) to compensate variations in the observing conditions.
