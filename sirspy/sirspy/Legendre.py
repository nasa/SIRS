import numpy as np

class Legendre():
    """
    Base Class for up-the-ramp Legendre fitting and modeling
    
    Parameters: nsamp, int
                  Number of equally spaced samples or groups up-the-ramp
                degree, int
                  Degree of Legendre polynomial fit
    """
    def __init__(self, nsamp, degree):
        
        # Pick off arguments
        self.nsamp = nsamp   # Number of up-the-ramp samples
        self.degree = degree # Fit degree
        
        # Build the Legendre basis matrix. In doing so, we allow
        # for the very first sample to be a "virtual" sample taken
        # immediately after reset. Although the Legendre polynomials
        # are defined over the interval x ∈ [-1,+1], our basis vectors
        # always have x > -1.
        self.x = (2*np.arange(nsamp+1)/nsamp-1)[1:]
        self.B = np.empty((nsamp,degree+1), dtype=np.float) # Empty basis matrix
        p = np.zeros((degree+1), dtype=np.int) # Used to pick out Legendre polynomials
        p[0] = 1 # Initialize it
        for col in np.arange(degree+1):
            self.B[:,col] = np.polynomial.legendre.legval(self.x,np.roll(p, col))
        
        # Compute its Moore-Penrose inverse, which does the fitting
        self.pinvB = np.linalg.pinv(self.B)
        
        # Compute the matrix that does the modeling
        self.B_x_pinvB = np.matmul(self.B, self.pinvB)
        
    def legfit(self, D):
        """
        Legendre fit a datacube
        
        Parameters: D, numpy.ndarray
                      The input datacube. This has been tested
                      with floating point inputs.
        Returns:
          * The Legendre fit. It can be converted to other things
            as follows.
              - Integrated DN = 2*λ1
              - "Slope" = 2*λ1/(NAXIS3-1)
        """
        return(np.einsum('ij,jkl',self.pinvB, D))