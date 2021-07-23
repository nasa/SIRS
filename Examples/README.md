# Simple Improved Reference Subtraction Examples

Bernard.J.Rauscher@nasa.gov<br>
NASA Goddard Space Flight Center

Simple Improved Reference Subtraction (SIRS) uses training data to compute frequency dependent weights that can be used to make reference corretions using an HxRG's reference columns. The training data are required to be a large set of up-the-ramp sampled darks.

Early on, this Julia language SIRS distribution contained both the "front-end" that solves for frequency dependent weights and the "back-end" that applies them to data. In practice, I have found that astronomers have a strong preference for python. For this reason, on 25 May 2021, I stopped working on the Julia back-end to focus all of my attention on the python-language sirspy back-end.

The examples below therefore show only how to use the SIRS front-end. See the python-3 sirspy distribution for examples of how to apply SIRS reference correction using the results of these notebooks.

Table 1. Example Notebooks

| File | Description |
| ------ | ------ |
| 20210415_find_sirs_weights_one_roman_sca.ipynb | Solve for SIRS frequency dependent weights for one Roman flight candidate  |
| 20210525_find_sirs_weights_all_roman_scas.ipynb | Solve for SIRS frequency dependent weights for all Roman flight candidates so far. |
| 20210428_find_sirs_weights_jpl_ppl.ipynb | Find SIRS weights for JPL Precision Projector Lab |
| 20210525_find_sirs_weights_jwst_nircam.ipynb | Find SIRS weights for JWST NIRCam |
