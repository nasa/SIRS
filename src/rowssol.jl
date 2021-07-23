"""
    rowssol(input_filenames, sirscore_filename; retval)

Solve for LAML reference pixels in rows weights

    Parameters: input_filenames::Array{String,1}
                  List of input files
                sirscore_filename::String
                  The name of a SIRSCore save file
                retval::Bool
                  If true, returns the reference row weights. This
                  may be useful for development.
    Returns: Nothing
               Upon completion, this creates a JLD file in the same directory
               as sirscore_filename having the suffix _SIRS2.jld.
"""
function rowssol(LegendreMatrices::LegendreMatrices, input_filenames::Array{String,1},
                sirscore_filename::String; retval::Bool=false)
    
    # Load the sirscore as a dict and get needed information
    d = JLD.load(sirscore_filename)
    nout = d["nout"]            # Number of outputs
    xsize = d["xsize"]          # Width of each output
    rb = d["rb"]                # Reference border width
    naxis3 = d["naxis3"]        # Number of up-the-ramp samples
    
    # The results will go here
    R = Array{Float64,3}(undef, (length(input_filenames)*naxis3,2*xsize*rb,nout)) # Reference rows
    μ = Array{Float64,2}(undef, (length(input_filenames)*naxis3,nout)) # Regular pixel means
    
    idx = 1 # A counter
    for file in input_filenames
        
        # Display status
        println("rowssol() processing file: ", file)
        flush(stdout)
    
        # Read data
        f = FITS(file)
        D = -Float64.(read(f[2])[:,:,:,1]) # Multiply by -1x so that DCL
                                           # charge integrates up
        close(f)
    
        # Subtract out the line fit
        D .-= legval(LegendreMatrices, legfit(LegendreMatrices, D))
    
        # Compute needed information frame-by-frame
        for z in 1:naxis3
            Threads.@threads for op in 1:nout
            
                # Get columns of interest
                x1 = (op-1)*xsize + 1 # First col. this output
                x2 = x1 + xsize -1    # Last col. this output
            
                # Get robust mean of regular pixels.
                μ[idx,op] = mean(trim(reshape(D[x1:x2,rb+1:end-rb,z], :), prop=TRIM_HARD))
            
                # Get reference rows
                R[idx,:,op] = reshape(cat(D[x1:x2,1:rb,z], D[x1:x2,end-rb+1:end,z], dims=2), :)         
            end
            idx += 1 # Increment counter
        end
    
    end # Done looping over files
    
    # Compute reference pixel weights
    w = ein"ijk,jk->ik"(mapslices(pinv, R, dims=(1,2)),μ)
    
    # Save result
    output_filename = chsuf(sirscore_filename, "_SIRS.jld", "_SIRS2.jld")
    JLD.save(output_filename, merge(d, Dict("w"=>w)))
    
    # Return the weights if requested
    if retval == true
        return(w)
    end
    
end