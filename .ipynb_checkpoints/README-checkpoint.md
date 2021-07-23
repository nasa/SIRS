# Simplified Improved Reference Subtraction (SIRS)

Bernard J. Rauscher<br>
NASA Goddard Space Flight Center

Teledyne’s H4RG, H2RG, and H1RG near-infrared array detectors provide reference pixels embedded in their data streams. Although they do not respond to light, the reference pixels electronically mimic normal pixels and track correlated noise. In this paper, we describe how the reference pixels can be used with linear algebra and machine learning to optimally reduce correlated noise. Simple Improved Reference Subtraction (SIRS) works with common detector clocking patterns and, when applicable, relies only on post-processing existing data. The resulting reference correction is optimal, in a least squares sense, when the embedded reference pixels are the only references available. We demonstrate SIRS using H4RG ground test data from the Nancy Grace Roman Space Telescope Project. The SIRS software is freely available for download.

SIRS is written in Julia. Because python is more widely used in the astronomical community, SIRS proves a python-3 "back end" that can be used to reference correct data given a pre-computed SIRS calibration file.