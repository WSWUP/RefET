Relationship to FAO-56
======================

A common question is how the ASCE standardized reference ET computed by this
package relates to the FAO-56 Penman-Monteith equation (Allen et al., 1998).
Short answer: for daily time steps, the ASCE short reference (ETos) is the
FAO-56 equation, so there is no need for a separate FAO-56 implementation.

Daily time step
---------------

The ASCE-EWRI (2005) standardized equation with the short reference surface
coefficients (Cn = 900, Cd = 0.34) has the same form and the same constants as
the FAO-56 Penman-Monteith ETo equation. Given the same inputs,
``refet.Daily(...).eto()`` is the FAO-56 daily reference ET. Small numerical
differences relative to another FAO-56 implementation can only come from
choices in the supporting calculations (for example the clear sky radiation
formulation used for cloudiness), not from the reference equation itself.

Hourly time step
----------------

The hourly equations differ. FAO-56 uses a constant denominator coefficient
(Cd = 0.34) for all hours, while the ASCE standardization refined the hourly
short reference to Cd = 0.24 for daytime and 0.96 for nighttime, with soil
heat flux ratios G/Rn of 0.1 (day) and 0.5 (night), where nighttime is defined
by Rn < 0. Hourly ETos from this package therefore follows ASCE-EWRI 2005 and
is not numerically identical to FAO-56 hourly ETo.

Surfaces at a glance
--------------------

=========  ============================  =====================  ==================
Surface    Description                   Daily Cn / Cd          Hourly Cn / Cd
=========  ============================  =====================  ==================
ETos       Short (clipped grass)         900 / 0.34             37 / 0.24 (day), 0.96 (night)
ETrs       Tall (full-cover alfalfa)     1600 / 0.38            66 / 0.25 (day), 1.7 (night)
=========  ============================  =====================  ==================

FAO-56 has no tall reference; ETr is specific to the ASCE standardization.

References
----------

Allen, R.G., Pereira, L.S., Raes, D., and Smith, M. (1998). Crop
evapotranspiration: Guidelines for computing crop water requirements. FAO
Irrigation and Drainage Paper 56. https://www.fao.org/4/x0490e/x0490e00.htm

ASCE-EWRI (2005). The ASCE standardized reference evapotranspiration equation.
https://ascelibrary.org/doi/book/10.1061/9780784408056
