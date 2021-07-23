# Simplified Improved Reference Subtraction (SIRS)

## Version 2 Release Note

The main difference between Rev. and Rev. 2 is how bad pixels and transients are handled. In Rev. 1, I used interpolation to fill in the gaps in the time series. This turned out to be computationally very expensive. The execution time was about 30 minutes per 60 frame up-the-ramp sampled exposure. Processing a full 105 exposure noise experiment required over 48 hours.

This version replaces interpolation with: (1) new row overhead columns are filled with the mirror of the last pixels read and (2) isolated bad pixels are replaced by the mean value of the same row. With these changes, execution time for one exposure dropped from about 30 minutes to about 1 minute and 20 seconds.