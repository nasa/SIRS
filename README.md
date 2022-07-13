# Simplified Improved Reference Subtraction (SIRS)

Bernard J. Rauscher<br>
NASA Goddard Space Flight Center

Teledyne’s H4RG, H2RG, and H1RG near-infrared array detectors provide reference pixels embedded in their data streams. Although they do not respond to light, the reference pixels electronically mimic normal pixels and track correlated noise. SIRS uses linear algebra and least squares optimization to optimally reduce correlated noise. SIRS works with common detector clocking patterns and, when applicable, relies only on post-processing existing data. The resulting reference correction is optimal, in a least squares sense, when the embedded reference pixels are the only references available and they are treated as two timeseries of reference information (groups of reference columns on the left and right).

SIRS is written in Julia. Because python is more widely used in the astronomical community, SIRS proves a python-3 "back end" that can be used to reference correct data given a pre-computed SIRS calibration file.

## Modification History

### sirs_3.5.1
* Initial release

### sirs_3.5.2
This revision incorporates suggestions from Chris Willott of Hertzberg Astrophysics
* Improved rejection of alternating column noise (ACN)
* Support for operable pixel masks