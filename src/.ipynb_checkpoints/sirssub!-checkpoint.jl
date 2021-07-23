"""
    sirssub(SIRSCore, D; f_max=nothing, rowsonly=false)

SIRS reference subtraction.

    Parameters: sc::SIRSCore
                  A loaded SIRSCore structure
                D::Array{Float64,3}
                  A sampled up-the-ramp datacube
                f_max::Float64
                  Maximum frequency to use. Weights for higher
                  frequencies are zeroed out. If not set, all available
                  frequencies are used
                rowsonly::Bool
                  This is mostly for debugging. Use it to turn of SIRS
                  correction and apply only a legacy rowsonly correction
"""
function sirssub!(sc::SIRSCore, D::Array{Float64,3}; f_max=nothing, rowsonly=false)
    
    # Get necessary information
    naxis3 = size(D)[3]                    # Number of frames this ramp
    nsteps = (sc.xsize+sc.nroh) * sc.ysize # Number of time steps per frame
    
    # Loop over frames
    for z in 1:naxis3
        
        # Get the incomplete Fourier transforms of reference pixels on
        # left and right. Previous work has shown that there is no gain
        # from doing these in parallel.
        𝓵 = inc_rfft(sc.SFT, reshape(D[1:HXRG_RB,:,z], :))
        𝓻 = inc_rfft(sc.SFT, reshape(D[end-HXRG_RB+1:end,:,z], :))
        
        # Zero out high frequencies if requested
        if f_max != nothing
            𝓵[sc.𝒇 .> f_max] .= 0
            𝓻[sc.𝒇 .> f_max] .= 0
        end
        
        # Loop over outputs
        Threads.@threads for op in 1:sc.nout
            
            # Column range for this output
            x0 = (op-1)*sc.xsize + 1 # First column this output
            x1 = x0 + sc.xsize - 1 # Last column this output
            
            # Do SIRS correction unless turned off
            if rowsonly == false
            
                # Compute the Fourier transform of the reference correction
                # for this output. Recall, we are working with incomplete
                # Fourier transforms here.
                r = (sc.α[:,op] .* 𝓵) .+ (sc.β[:,op] .* 𝓻)

                # Invert the Fourier transform. We must invert the *complete*
                # Fourier transform including new row overheads.
                r2 = zeros(Complex{Float64}, nsteps÷2+1)    # Buffer for complete rfft
                r2[1:length(r)÷2+1] .= r[1:length(r)÷2+1]   # Low frequencies
                r2[end-(length(r)÷2-1):end] .= r[length(r)÷2+2:end] # High frequencies
                r = irfft(r2, nsteps) # Go back to time domain

                # Pick out steps that correspond  to normal pixels
                r = reshape(r, (sc.xsize+sc.nroh,sc.ysize))[1:sc.xsize,:]

                # flip even numbered outputs
                if mod(op,2) == 0
                    r = r[end:-1:1,:]
                end

                # Reference correct in place
                D[x0:x1,:,z] .-= r
                
            end
            
            # Use top 4 rows for DC correction. Just use the center 2 for
            # now. We will revisit this with LAML later. Use a more gentle clip
            # for reference pixels than for regular pixels.
            D[x0:x1,:,z] .-= mean(StatsBase.trim(reshape(D[x0:x1,sc.naxis2-2:sc.naxis2-1,z], :),
                    prop=TRIM_MED))
                        
        end
        
    end
    
end