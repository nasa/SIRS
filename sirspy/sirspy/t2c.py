import numpy as np

def t2c(d):
    """
    Convert arrays of tuples containing complex numbers to 
    numpy complex arrays.
            
    Parameters: d, Complex array formatted as array of tuples
                  Input data
    Returns: A numpy complex array
            
    Notes:
      * This can be time consuming. Help speeding it up
        would be appreciated.
    """
    result = np.zeros(d.shape, dtype=np.complex64)
    for r in np.arange(d.shape[0]):
        for c in np.arange(d.shape[1]):
            result[r,c] = np.complex(d[r][c][0], d[r][c][1])
    return(result)