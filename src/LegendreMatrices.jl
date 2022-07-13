"""
    LegendreMatrices(size)

Build Legendre basis matrix and it's Moore-Penrose inverse. These
are used for fitting and modeling ramps.

    Parameters: size::Tuple{Int64,Int64})
                  Size of the basis matrix. Each columns will contain a
                  Legendre basis vector. The number of rows should be
                  set equal to the number of up-the-ramp samples. For
                  example, if we want to fit straight lines to a 60
                  frame sampled up-the-ramp exposure, size = (60,2). The
                  number of columns should be set to the degree of the fit
                  plus 1.
"""
struct LegendreMatrices

    x       # "x" values for computing basis vectors
    B       # The Legendre basis matrix. Basis vectors are columns
    pinv_B  # Its Moore-Penrose inverse
    
    function LegendreMatrices(size::Tuple{Int64,Int64})
        
        # Pick out parameters for readability
        nrows, ncols = size
        
        #= Make x. The first element is a hypothetical "virtual"
        sample taken immediately after reset. We discard this 
        and don't try to fit it. =#
        x = collect(-1:2/(nrows):+1)[2:end]
        
        # Make the basis vectors
        B = Array{Float64,2}(undef, size)
        for j in 1:ncols
           B[:,j] = legendre.(x, j-1) 
        end
        
        # Invert
        pinv_B = pinv(B)
        
        # Done!
        new(x, B, pinv_B)
        
    end
    
end