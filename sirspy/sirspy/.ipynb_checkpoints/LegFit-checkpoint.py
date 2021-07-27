import numpy as np

class LegFit():
    """
    Base Class for up-the-ramp Legendre Polynomial fitting
    
    Parameters: nsamp, int
                  Number of samples up-the-ramp
                degree, int
                  Legendre fit degree
    """
    def __init__(self, nsamp, degree):
        
        # Pick off arguments
        self.nsamp = nsamp   # Number of up-the-ramp samples
        self.degree = degree # Fit degree
        
        # Build the Legendre basis matrix
        self.x = 2*np.arange(nsamp)/(nsamp-1)-1 # x-values for computing Legendre polynomials
        self.B = np.empty((nsamp,degree+1), dtype=np.float) # The empty matrix
        p = np.zeros((degree+1), dtype=np.int) # Used to pick out Legendre polynomials
        p[0] = 1 # Initialize it
        for col in np.arange(degree+1):
            self.B[:,col] = np.polynomial.legendre.legval(self.x,np.roll(p, col))
        
        # Compute its Moore-Penrose inverse, which does the fitting
        self.pinvB = np.linalg.pinv(self.B)
        
    def fit(self, D):
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