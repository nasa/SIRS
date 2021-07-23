import numpy as np
import numba

@numba.jit(nopython=True)
def make_fourier_matrices(j, nvec, nstep):
    """
    Make the Incomplete Fourier basis matrix and its Moore-Penrose inverse. The
    Inverse is what is used for fitting.
    
    Parameters: j, Array
                  Array of reference column sample numbers. This will
                  be something like, [0,1,2,3,136,137,138,139,...]
                nvec, int
                  Number of Fourier vectors. This should be set to
                  (len(freq)-1)//2 + 1, where the final +1
                  allows for DC.
                nstep, int
                  Total number of time steps in frame including the new row
                  overhead.
    Notes:
      * 4/12/2021, B.J. Rauscher, NASA/GSFC
        - Making this a function. I am not expert, but my quick reading suggests
          that numba works on functions, not methods. This was unfortunately a bit
          slow in standard python.
    """
    B = np.empty((len(j),nvec), dtype=np.complex64)
    for row in np.arange(len(j)):
        for col in np.arange(nvec):
            # Note normalization just after = sign
            # B[row,col] = (1./len(j)) * np.exp(-2*np.pi*1J*j[row]*col/nstep)
            B[row,col] = (1/nstep) * np.exp(2*np.pi*1J*j[row]*col/nstep)
    pinv_B = np.linalg.pinv(B)
    return(B, pinv_B)