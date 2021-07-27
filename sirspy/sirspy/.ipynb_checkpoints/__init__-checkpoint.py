# HxRG parameters
NAXIS1 = 4096          # Number of columns
NAXIS2 = 4096          # Number of rows
RLIM   = (4093,4094)   # Use these reference rows
NOUT   = 32            # Number of outputs
RB     = 4             # Reference pixel border width
XSIZE  = NAXIS1//NOUT  # Width of each output in columns
YSIZE  = NAXIS2        # Height of each output in columns

# Module imports
from .SIRS import SIRS
from .LegFit import LegFit
from .make_fourier_matrices import make_fourier_matrices
from .rowsonly import rowsonly
from .t2c import t2c
