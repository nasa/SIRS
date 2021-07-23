import numpy as np

from sirspy import XSIZE, RLIM, NOUT

def rowsonly(D, prop=.025):
    """
    Legacy rowsonly reference correction
    
    Parameters: D, ndarray
                  The input datacube
                prop, float
                  Discard prop % of low and high reference
                  row outliers. The total percent discarded
                  is = 2 * prop.
    """
    # sirspy globals
    global XSIZE, RLIM, NOUT
    
    # Get dimensions
    naxis3 = D.shape[0] # Number of samples up-the-ramp
    
    # We want very robust, but fast means. Discard the lowest
    # and highest pixels
    discard = np.int(np.round(prop*XSIZE*(RLIM[1]-RLIM[0]+1))) # Ignore this many pixel low and high
    
    # Do it
    for z in np.arange(naxis3):
        for op in np.arange(NOUT):
            x0 = op * XSIZE
            x1 = x0 + XSIZE
            ref = np.mean(np.sort(D[z,4093:4095,x0:x1].flatten())\
                                      [discard:-discard])
            # Reference correct
            D[z,:,x0:x1] -= ref
            
    # Done!
    return(D)
            
            