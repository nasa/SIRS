"""
    coadd(sc, D; pplfix1)

Co-add a sampled up-the-ramp datacube into a SIRSCore.

    Parameters: sc::SIRSCore
                  A SIRSCore struct
                D::Array{Float64,3}
                  An HxRG datacube
                pplfix1::bool (optional)
                  JPL PPL detector has some bad rows. In Julia's
                  index convention, these are rows [1872,2226].
                  Interpolate over them.
"""
function coadd!(sc::SIRSCore, D::Array{Float64,3}; pplfix1::Bool=false)
    
    # JWST NIRCam seems to require a DC reference correction on account
    # of the SIDECAR resets every frame. I don't think it can hurt to do this for
    # everything, so leave it in.
    for z in 1:sc.naxis3
        Threads.@threads for op in 1:sc.nout
            x0 = (op-1)*sc.xsize + 1 # First col this output
            x1 = x0 + sc.xsize - 1   # Last col this output
            R = D[x0:x1,sc.ysize-2:sc.ysize-1,z] # These ones tend to be stable
            μ = mean(trim(R[:], prop=TRIM_MED))  # Robust mean
            D[x0:x1,:,z] .-= μ # Reference correct
        end
    end
    
    # Compute residuals by fitting and subtracting a straight line.
    # A bit of testing shows that using MKL to do the matrix multiplications is
    # faster than using mapslices. In other contexts, I have found
    # that using MKL is faster than tensor operations. Therefore, 
    # overwrite D so that it will work with MKL.
    D = permutedims(D, (3,1,2)) # This is expensive, on Racy ~14 seconds but done only once per exposure
    Δ = D - reshape(sc.L_x_Linv * reshape(D, (sc.naxis3,sc.naxis1*sc.naxis2)),
            (sc.naxis3,sc.naxis1,sc.naxis2))

    # Work in frames
    for z in 1:sc.naxis3
        
        # Show status
        # println("Processing frame ", z)
        # flush(stdout)
        
        # Pick out one frame
        frm = Δ[z,:,:]
        
        #=
        The JPL PPL detector has a bunch of bad rows on the right. These are rows [1872,2226]
        (unity offset). Interpolate over them.
        =#
        if pplfix1 == true
            x = collect(1:4096) # All x values
            badcols = 4093:4096    # The bad columns
            for col in badcols
                _x = x[(x .< 1871) .| (x .> 2225)]
                _y = frm[col, (x .< 1871) .| (x .> 2225)]
                spl = Spline1D(_x, _y, k=1) # Linear interpolation
                frm[col, (x .>= 1871) .& (x .<= 2225)] .= spl.(x[(x .>= 1871) .& (x .<= 2225)])
            end     
        end        
        
        # Get a copy of the good pixel mask. Make a copy so
        # as not to mess it up.
        gdpx = copy(sc.gdpx)
        
        # Compute incomplete Fourier transforms of reference pixels on left and right
        # Note added 2/21/22. I tried parallelizing these two lines on Racy.
        # There was no speed improvement. The overheads associated with starting
        # the additional threads were larger than the improvement.
        # From timing the code, these are actually very fast.
        #
        # On 4/20, it became apparent while working with JPL PPL data that there
        # needs to be quality control in the reference columns. We add that now.
        lft = frm[1:HXRG_RB,:][:]            # get left reference columns as vector
        rgt = frm[end:-1:end-HXRG_RB+1,:][:] # get right reference columns as vector
                                           # Note that this output gets flipped
        
        #=
        Try subtracting a smoothed version of the data from the data
        and filling with that. We are after transients here, not
        bad clusters of reference pixels. For this detector, there
        appears to be such a cluster on the right.
        =#
        sigrej = 4.0 # Sigma clipping threshold
        sg_ftr = SG_Filter(Float64, halfWidth=5, degree=3)
        # Left
        sm_lft= apply_SG_filter(lft, sg_ftr, derivativeOrder=0)
        dif = lft .- sm_lft
        trimmed = trim(dif, prop=TRIM_SOFT)
        μ = mean(trimmed)
        σ = std(trimmed)
        lft[(dif .< μ-sigrej*σ) .| (dif .> μ+sigrej*σ)] .=
            sm_lft[(dif .< μ-sigrej*σ) .| (dif .> μ+sigrej*σ)]
        # sum_lft = sum((dif .< μ-sigrej*σ) .| (dif .> μ+sigrej*σ)) # Stub
        # Right
        sm_rgt= apply_SG_filter(rgt, sg_ftr, derivativeOrder=0)
        dif = rgt .- sm_rgt
        trimmed = trim(dif, prop=TRIM_SOFT)
        μ = mean(trimmed)
        σ = std(trimmed)
        rgt[(dif .< μ-sigrej*σ) .| (dif .> μ+sigrej*σ)] .=
            sm_rgt[(dif .< μ-sigrej*σ) .| (dif .> μ+sigrej*σ)]
        # sum_rgt = sum((dif .< μ-sigrej*σ) .| (dif .> μ+sigrej*σ)) # Stub
        # println("L: ", sum_lft, "  R: ", sum_rgt) # Stub
        # flush(stdout) # Stub


        
        # Incomplete Fourier transforms after rejecting outliers
        𝓵 = inc_rfft(sc.SFT, lft)
        𝓻 = inc_rfft(sc.SFT, rgt)
        
        # Loop over outputs
        Threads.@threads for op in 1:sc.nout
        
            # Get column range for this output
            c1 = (op-1) * sc.xsize + 1
            c2 = c1 + sc.xsize - 1
            crng = c1:c2
        
            # Get data for this output.
            #   1) We overwrite all elements, so it is faster to leave this
            #      uninitialized.
            #   2) Include NROH additional columns for the new row overhead
            #      to eventually hold new row overhead
            #   3) Flip axis 1 of even numbered columns
            d = Array{Float64, 2}(undef, (sc.xsize+sc.nroh,sc.ysize))
            if mod(op,2) != 0
                d[1:sc.xsize,:] = frm[crng,:]             # Don't flip
            else 
                d[1:sc.xsize,:] = frm[crng,:][end:-1:1,:] # Flip
            end
            
            # Get the good pixel map for this output flipping
            # axis 1 as necessary
            _gdpx = BitArray{2}(Array{Bool,2}(undef, (sc.xsize,sc.ysize)))
            if mod(op,2) != 0
                _gdpx = gdpx[crng,:]
            else
                _gdpx = gdpx[crng,:][end:-1:1,:]
            end
            
            # Find and flag transients. This only makes sense
            # in the regular pixels. Known bad pixels are already
            # flagged. Trim thresholds are defined in SIRS.jl.
            if (op==1) || (op==32)
                regpix_mask = sc.regpix_edge
            else
                regpix_mask = sc.regpix_middle
            end
            regular_pixels = d[1:sc.xsize,:][regpix_mask]
            good_pixels = _gdpx[regpix_mask]
            sorted = sort(regular_pixels[good_pixels])
            nrej = Int64(round(TRIM_HARD * length(sorted))) # Figure out how many to trim
                                                            # on either side
            min_good = sorted[nrej]
            max_good = sorted[end-nrej+1]
            good_pixels[regular_pixels .<= min_good] .= 0
            good_pixels[regular_pixels .>= max_good] .= 0
            _gdpx[regpix_mask] = good_pixels
            
            # Output level quality control.
            # Compute the mean and standard deviation of remaining good pixels.
            # Use this to find any anomalous behavior in this output and frame.
            # Also compute the mean since it will be useful later.
            here = d[1:sc.xsize,HXRG_RB+1:end-HXRG_RB]
            μ = mean(here[_gdpx[:,HXRG_RB+1:end-HXRG_RB]])
            σ = std(here[_gdpx[:,HXRG_RB+1:end-HXRG_RB]])
            if ! ((GD_OP_SIG_MIN .< σ) * (σ .< GD_OP_SIG_MAX))
                println("Skipping frame ", z, " output ", op, " σ = ", σ)
                flush(stdout)
                continue
            end
            
            # Get the (robust) mean reference row value.
            # From looking at a lot of data, the rows selected here tend to be 
            # most indicative of the regular pixels in Roman H4RGs.
            μ_refrows = mean(trim(reshape((d[1:sc.xsize,sc.naxis2-2:sc.naxis2-1]), :),
                                prop=TRIM_HARD))
            
            # Replace known bad pixels and transients with the mean
            # of all good pixels in the same row. On 4/19/21, this part of the
            # code cause problems with JPL Precision Projector Lab data.
            # The data contained lines with few (or no) good pixels. The code
            # now replaces such lines with the output's mean value.
            for y in 1:sc.ysize
                # Only touch rows that contain bad pixels
                if !(0 ∈ _gdpx[:,y])
                    continue
                end
                here = d[1:sc.xsize,y]
                # Ensure that there is at least one good pixel
                if sum(_gdpx[:,y]) > 0
                    here[_gdpx[:,y] .== 0] .= mean(here[_gdpx[:,y] .== 1])
                    d[1:sc.xsize,y] = here
                else
                    # This handles the case of all pixels in the line being bad
                    d[1:sc.xsize,y] .= μ
                end
            end
            
            # Fill overhead columns by mirroring
            roi = d[end-2sc.nroh+1:end-sc.nroh,:]
            d[end-sc.nroh+1:end,:] = roi[end:-1:1,:]
        
            # Go to vectors. View the reference pixel stream as an xy-plot.
            # The x-axis is pixel index and the y-axis is signal in DN.
            d_vec = reshape(d,:)
             
            # Compute the FFT
            𝓷 = rfft(d_vec)
        
            # Keep just frequencies of interest. These are less than Nyquist
            # on the row rate and within the same frequency interval of Nyquist.
            𝓷 = cat(𝓷[1:sc.naxis2÷2+1],𝓷[length(𝓷)-(sc.naxis2÷2-1):end], dims=1)
        
            # Coadd sums for frequencies > 0 Hz
            sc.ℕ[:,op] .+= real.(𝓷 .* conj(𝓷))
            sc.𝕃[:,op] .+= real.(𝓵 .* conj(𝓵))
            sc.ℝ[:,op] .+= real.(𝓻 .* conj(𝓻))
            sc.𝕏[:,op] .+= 𝓷 .* conj(𝓻)
            sc.𝕐[:,op] .+= 𝓷 .* conj(𝓵)
            sc.ℤ[:,op] .+= 𝓻 .* conj(𝓵)
            
            # Coadd sums for f = 0 Hz
            sc.R[op]    += μ - μ_refrows
            sc.N[op]    += 1          # Serves as frame counter
            
        end
    end
end