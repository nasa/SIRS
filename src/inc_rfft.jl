"""
    inc_rfft_pars(j, k, n [; restore=false])

This struct contains the information that is need to compute incomplete
rffts.

    Parameters: j, Array{Int64,1}
                  A vector of j values, j in {1,2,... n}. These specify
                  which time steps are used to compute the incomplete Fourier
                  transform.
                k, Array{Int64,1}
                  A vector of k values, k in {1,2,... n/2+1}; where
                  n/2 is rounded down. These are the columns in the basis
                  matrix.
                n, Int
                  Number of data elements. For SIRS, this is the
                  number of time steps per output per frame.
                restore, Bool
                  Set this if restoring from a save set
"""
struct inc_rfft_pars
   
    j    # Vector of j values
    k    # Vector of k values
    n    # Number of data elements
    B    # Incomplete Fourier basis matrix
    Binv # Moore-Penrose inverse of B
    
    function inc_rfft_pars(j::Array{Int64,1},
                    k::Array{Int64,1}, n::Int64; restore::Bool=false)
        
        if restore == false
            # Build the basis matrix. The basis vectors are the ones normally
            # used to compute irfft(). See Eqns. 4-5 in docs/fftw_rfft_basis_vectors.ipynb.
            B = Array{Complex{Float64}}(undef, (length(j), length(k)))
            Threads.@threads for pk in 1:length(k)
                B[:,pk] = n^(-1) .* exp.(2*π*sqrt(Complex(-1)) .* (j .- 1) .* ((k[pk]-1)/n))
            end
            # Compute it's Moore-Penrose inverse
            Binv = pinv(B)
        else
            # Just make place holders if restoring
            B = Array{Complex{Float64}}(undef, (length(j), length(k)))
            Binv = Array{Complex{Float64}}(undef, (length(k), length(j)))
        end
        
        # Instantiation
        new(j, k, n, B, Binv)
        
    end
        
end


"""
    inc_rfft(FT, d)

Compute the incomplete rfft of the input data.

    Parameters: FT::inc_rfft_pars
                  A inc_rfft_pars struct
                d::Array{Float64,1}
                  A data vector
"""
function inc_rfft(FT::inc_rfft_pars, d)
   return(FT.Binv * d) 
end
