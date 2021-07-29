# README for SIRSPY

Bernard J. Rauscher<br>
NASA Goddard Space Flight Center

This package provides the SIRS "back-end" implemented in python-3. It can be used to reference correct up-the-ramp sampled data given a SIRS calibration file. Sirspy also provides a few simple tools for fitting and modeling SIRS corrected data cubes.

## 1. Installation

It is easiest to install sirspy as a local (editable) package using pip. To do this, type

your_unix$ `pip install -e <path>`

where <path> is the path to sirspy on your machine.

## 2. Hints

### 2.1 Speeding it Up

Sirspy is computationally intensive. The best way to speed it up is to turn on the multithreaded linear algebra that is already built into numpy (by default, python and numpy start single-threaded). To do this, use `numpy.__config__.show()` to determine whether you are using MKL or OpenBLAS. On Intel machines, MKL is generally faster. If you are on an Intel machine, but numpy is using OpenBLAS, you may want to recompile it to use MKL if you do a lot of linear algebra.

Once you know what library numpy is using, figure out how many cores you want to allocate to it. Do this using the unix `lscpu` command. My machine has 8 cores, and I allow MKL to use all of them. You can do this in your python script by setting, `os.environ['MKL_NUM_THREADS'] = '8'`, or `os.environ['OPENBLAS_NUM_THREADS'] = '8'` (substitute however many threads you want to allow for '8').

Alternatively, you can set these envars in your `BASH` initialization script.

## 3. Revision History

26 July 2021, B.J.Rauscher, NASA/GSFC
* Now sirs_3.5.0
* Initial NASA GitHub release

29 July 2021, B.J.Rauscher, NASA/GSFC
* Now sirs_3.5.1
* Refined JWST NIRISS sirspy Example
* Refined JWST NIRCam sirspy Example
