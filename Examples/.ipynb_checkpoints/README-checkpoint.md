# Simple Improved Reference Subtraction Examples

Bernard.J.Rauscher@nasa.gov<br>
NASA Goddard Space Flight Center

Simple Improved Reference Subtraction (SIRS) uses training data to compute frequency dependent weights that can be used to make reference corretions using an HxRG's reference columns. The training data are required to be a large set of up-the-ramp sampled darks.

SIRS consists of both a "front end" that computes the weights and a "back end" that applies them. Although the front end computation is currently available only in Julia, if the weights are known, there is a python-3 back end.

Table 1. Example Notebooks

| File | Description |
| ------ | ------ |
| 20210415_find_sirs_weights_one_roman_sca.ipynb | Solve for SIRS frequency dependent weights for one Roman flight candidate  |
| 20210525_find_sirs_weights_all_roman_scas.ipynb | Solve for SIRS frequency dependent weights for all Roman flight candidates so far. |
| 20210428_find_sirs_weights_jpl_ppl.ipynb | Find SIRS weights for JPL Precision Projector Lab |
| 20211203_find_sirs_weights_jwst_nircam.ipynb | Find SIRS weights for JWST NIRCam |
| 20211203_sirs_correct_jwst-nircam.ipynb | Apply SIRS correction top JWST NIRCam |

