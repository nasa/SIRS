import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

from .t2c import t2c
from .make_fourier_matrices import make_fourier_matrices


class SIRS():
    
    RB = 4 # Reference pixel border width
    
    def __init__(self, sirs_file):
        """
        __init__(sirs_file)
            
        Instantiate a SIRS object
        
        Parameters: sirs_file:string
                      Name of a SIRS weights file. The suffix will be .jld.
        """
        
        # Open the cal. file
        f = h5py.File(sirs_file, "r")
        
        # Recover just the parameters needed to
        # apply SIRS reference correction
        self.naxis1 = np.int64(f['naxis1'][...])       # Number of columns
        self.naxis2 = np.int64(f['naxis2'][...])       # Number of rows
        self.nout = np.int64(f['nout'][...])           # Number of outputs
        self.nroh = np.int64(f['nroh'][...])           # New row overhead in pixels
        self.xsize = np.int64(f['xsize'][...])         # x-size one output
        self.ysize = np.int64(f['ysize'][...])         # y-size one output)
        self.nstep = (self.xsize+self.nroh)*self.ysize # Total number of time steps
        self.freq = np.array(f['𝒇'][...],
                             dtype=np.float64)         # Incomplet Fourier transform frequencies
        self.α = f['α'][...]                           # SIRS weights
        self.β = f['β'][...]                           # SIRS weights
        
        # Parameters used for DC correction
        self.rowslim = (self.naxis2-3,self.naxis2-2) # 1st and last reference rows to use
        self.discard = np.int(.005 * (self.rowslim[1]-self.rowslim[0]+1) * self.xsize) # Discard this many elements 
                                                                                       # on top and bottom of distribution
                                                                                       # when robustly computing reference
                                                                                       # rows mean.
        
        # The SIRS calibration files contain the Fourier basis vectors, SFB_B,
        # that were used to compute the incomplete Fourier transform and their
        # Moore-Penrose inverse, SFT_Binv. Unfortunately, they unpack to non-standard
        # numpy arrays in python and I was not able to convert them quickly. I therefore 
        # re-compute the basis vectors using make_fourier_matrices().
        # self.B = f['SFT_B'][...]            # Fourier basis vectors
        # self.pinv_B = f['SFT_Binv'][...]    # Moore-Penrose inverse of above
        self.j = np.array(f['SFT_j'][...],
                          dtype=np.int64)     # A vector of j values, j in {1,2,... n}.
                                              # These specify which time steps are used
                                              # to compute the incomplete Fourier transform
        self.j -= 1                           # Python uses zero-offset arrays
        # Convert from arrays of tuples to complex arrays
        self.α = t2c(self.α)
        self.β = t2c(self.β)
        
        # Make the Fourier basis vectors and compute the Moore-Penrose inverse.
        # This is available in the JLD file, but it turns out to be faster in python to just
        # compute it anew. The basis vectors can be seen in the "Implementation Details" here:
        # 
        #     https://numpy.org/doc/stable/reference/routines.fft.html.
        #
        # Figure out how many basis vectors are needed to do just the low frequencies
        # shown in the JATIS article, Fig 3a
        #
        # ***** FOR JWST WE NEED TO KEEP THE HIGH FREQUENCIES TOO. THIS IS DONE IN THE
        #       JULIA CODE, BUT NOT SIRSPY AS OF NOW.
        self.nvec = (len(self.freq)-1)//2 + 1 # +1 is for zero frequency
        self.B, self.pinv_B = make_fourier_matrices(self.j,self.nvec,self.nstep)
        
    def plot(self, op, title="", mag=1.0):
        """
        Plot alpha and beta
        
        Parameters: op, int
                      Output number ∊ {0,1,... n_out-1}, where
                      n_out is the number of available outputs
                    title, string
                      Plot title
                    mag, float
                      Plot magnification
        """
        
        # Define ranges of frequencies to plot
        nfreq = len(self.freq)
        f_low_min  = self.freq[0]                       # Minimum frequency to plot
        f_low_max  = self.freq[(nfreq-1)//2] # Maximum frequency to plot
        f_high_min = self.freq[(nfreq-1)//2+1]
        f_high_max = self.freq[-1]
        
        # Make an over-under plot
        fig, ax = plt.subplots(2, 2, figsize=(mag*6.4, mag*4.8)) # Over-under plots
        
        # Plot low frequency amplitude
        ax[0,0].plot(self.freq[self.freq <= f_low_max],
                     np.abs(self.α[op, self.freq <= f_low_max]), '-', alpha=.7, label='$\\alpha$')
        ax[0,0].plot(self.freq[self.freq <= f_low_max],
                     np.abs(self.β[op, self.freq <= f_low_max]), '-', alpha=.7, label='$\\beta$')
        ax[0,0].set_ylim(-.05,1)
        ax[0,0].set_ylabel('Amplitude')
        ax[0,0].legend(loc='upper right')
        if title != "":
            ax1.set_title(title)
            
        # Plot high frequency amplitude
        ax[0,1].plot(self.freq[self.freq >= f_high_min] - self.freq[-1], np.abs(self.α[op, self.freq >= f_high_min]), '-', alpha=.7)
        ax[0,1].plot(self.freq[self.freq >= f_high_min] - self.freq[-1], np.abs(self.β[op, self.freq >= f_high_min]), '-', alpha=.7)
        ax[0,1].set_ylim(-.05,1)
        ax[0,1].yaxis.set_ticklabels([])
        
        # Plot low frequency phase
        ax[1,0].plot(self.freq[self.freq <= f_low_max], np.angle(self.α[op, self.freq <= f_low_max]), '.', alpha=.2)
        ax[1,0].plot(self.freq[self.freq <= f_low_max], np.angle(self.β[op, self.freq <= f_low_max]), '.', alpha=.2)
        ax[1,0].set_ylim(-np.pi,+np.pi)
        ax[1,0].set_xlabel('Frequency (Hz)')
        ax[1,0].set_ylabel('Phase')
        
        # Plot high frequency phase
        ax[1,1].plot(self.freq[self.freq >= f_high_min] - self.freq[-1], np.angle(self.α[op, self.freq >= f_high_min]), '.', alpha=.2)
        ax[1,1].plot(self.freq[self.freq >= f_high_min] - self.freq[-1], np.angle(self.β[op, self.freq >= f_high_min]), '.', alpha=.2)
        ax[1,1].set_ylim(-np.pi,+np.pi)
        ax[1,1].yaxis.set_ticklabels([])
        ax[1,1].set_xlabel('Frequency ($-f_{\\rm Ny}$ Hz)')
        
        # Done
        return(plt)
    
    def incomplete_ft(self, d):
        """
        Incomplete Fourier Transform.
        See Documentation/derivation_of_incomplete_ft.ipynb.
        
        Parameters: d, Data vector
                      The input data vector
        """
        return(np.matmul(self.pinv_B, d))
        
        
    def refcor(self, D, pplfix1=False, rowsonly=False):
        """
        SIRS reference correction
        
        Parameters: D, Datacube
                      The input datacube
                    pplfix1, Bool
                      The PPL detector has some bad reference pixels. Interpolate over them.
                    rowsonly, False
                      Optionally do a rowsonly correction. This uses the most stable reference
                      rows only. This is useful for studying the effect of SIRS.
        Notes:
          * This method overwrites the input data
        """
        nframes = D.shape[0] # Number of frames to correct
        
        # Work frame-by-frame...
        for z in np.arange(nframes):
            
            # Pick off reference columns
            l = D[z,:,:self.RB]
            r = D[z,:,-self.RB:]
            
            # Deal with bad reference pixels in the JPL PPL detector
            if pplfix1==True:
                x = np.arange(4096) # All possible row indices
                px = (x<1870) | (x>2224) # Bool
                _x = x[px]
                for col in np.arange(4):                
                    _y = r[px,col]
                    spl = interpolate.interp1d(_x, _y, kind='linear')
                    r[np.logical_not(px),col] = spl(x[np.logical_not(px)])
                
            # Go to Fourier space
            𝓵 = self.incomplete_ft(l.flatten()) # Project left refcols into Fourier space
            𝓻 = self.incomplete_ft(r.flatten()) # Project right refcols...
                        
            # Work output by output...
            for op in np.arange(self.nout):
                
                # We need the range of columns for this output
                x0 = op*self.xsize
                x1 = x0 + self.xsize

                # Do SIRS unless disabled using rowsonly
                if rowsonly == False:
        
                    # Work out reference correction for this output
                    ref = np.zeros(self.nstep//2+1, dtype=np.complex64) # Build full rfft here

                    # Low frequencies
                    ref[:self.ysize//2+1] = self.α[op,:self.ysize//2+1]*𝓵[:self.ysize//2+1] +\
                                            self.β[op,:self.ysize//2+1]*𝓻[:self.ysize//2+1]

                    # High frequencies
                    ref[-self.ysize//2:] = self.α[op,-self.ysize//2:]*𝓵[-self.ysize//2:] +\
                                           self.β[op,-self.ysize//2:]*𝓻[-self.ysize//2:]

                    # Invert the rfft
                    ref = np.fft.irfft(ref, n=self.nstep)

                    # Reformat as 2D image
                    ref = ref.reshape(self.ysize,self.xsize+self.nroh)

                    # Keep just real samples
                    ref = ref[:,:self.xsize]

                    # Flip odd numbered outputs
                    if np.mod(op,2)==1:
                        ref = np.fliplr(ref)


                    # A stub used for debugging. Leaving it in because
                    # it may be needed again. To use this, limit the number
                    # of frames and outputs to 1.
                    # _d = D[z,:,:4]
                    # _d -= np.mean(_d)
                    # _r = ref[:,:4]
                    # _r -= np.mean(_r)
                    # plt.plot(np.mean(_d, axis=1), '.', alpha=.2)
                    # plt.plot(np.mean(_r, axis=1), '-')
                    # return(plt)

                    # SIRS reference correct data
                    D[z,:,x0:x1] -= ref
                
                # Correct DC using reference rows
                ref = np.mean(np.sort(D[z,self.rowslim[0]:self.rowslim[1]+1,
                                        x0:x1].flatten())[self.discard:-self.discard])
                D[z,:,x0:x1] -= ref
                
