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
                  Number of low frequency Fourier vectors. This should be
                  set to (len(freq)-1)//2 + 1, where the final +1
                  allows for DC.
                nstep, int
                  Total number of time steps in frame including the new row
                  overhead.
    Notes:
      * 4/12/2021, B.J. Rauscher, NASA/GSFC
        - Making this a function. I am not expert, but my quick reading suggests
          that numba works on functions, not methods. This was unfortunately a bit
          slow in standard python.
      * 7/26/2021, B.J. Rauscher, NASA/GSFC
        - The initial implementation was for Roman Space Telescope H4RGs.
          These do not seem to have much alternating column noise. For this reason,
          we initially did not include a correction at high frequency. Unfortunately,
          JWST's H2RGs do have significant alternating column noise. We therefore have
          added the high frequencies back in.
    """
    nrows = len(j)       # Number of rows in Fourier basis matrix
    ncols = 2*(nvec-1)+1 # Number of columns in Fourier basis matrix. This includes
                         # f = 0 Hz, len(freq)//2 low frequencies, len(freq)//2 high frequencies
    B = np.empty((nrows,ncols), dtype=np.complex64) # Create the basis matrix
    for row in np.arange(len(j)):
        for col in np.arange(nvec):
            
            # Low frequency Fourier vectors
            B[row,col] = (1/nstep) * np.exp(2*np.pi*1J*j[row]*col/nstep)
            
            # High frequency Fourier vectors
            if col > 0:
                k_max = nstep-1  # Maximum basis vector index. See 
                                 # https://numpy.org/doc/stable/reference/routines.fft.html
                B[row,ncols-col] = (1/nstep) * np.exp(2*np.pi*1J*j[row]*(k_max-col+1)/nstep)
                
    # The Moore-Penrose inverse is needed for modeling            
    pinv_B = np.linalg.pinv(B)
    
    # Done!
    return(B, pinv_B)