"""
    SIRSCore()

Simple Improved Reference Subtraction (SIRS) core structure. This
contains the information that is needed to compute the SIRS alpha
and beta vectors from a set of up-the-ramp sampled darks.

    Parameters: naxis3::Int64
                  Number of up-the-ramp samples. If set =0, this tells
                  SIRSCore to 
                gdpx::BitArray{2}
                  The detector's good pixel map. Good pixels are =1. Bad pixels
                  and reference pixels are =0.
                nroh::Int64
                  New row overhead in pixels.
                restore::Bool
                  Set this true when restoreing a saved SIRSCore.

                The remaining kwargs are mostly FITS header keywords

    Notes:
      * 3/2/21: Computing SFT is expensive. Consider adding ability to save and reuse.
"""
struct SIRSCore

    # Constants that define the data
    hxrg_kind # Selected from {"h4rg", "h2rg"}
    naxis1    # Number of HxRG columns
    naxis2    # Number of HxRG rows
    naxis3    # Z-size of full array (number of UTR frames)
    xsize     # Number of columns per output
    ysize     # Number of rows per output (redundant with naxis2)
    nout      # Number of HxRG outputs
    nroh      # New row overhead in pixels
    n         # Number of time steps per output
    τ         # Pixel dwell time in seconds
    j_ref     # Reference pixel indices
    j_reg     # Regular pixel indices
    
    # An incomplete Fourier transform is used to compute the spectral properties
    # of the reference columns. See the inc_rfft() and inc_rfft_pars()
    # methods that are defined in this package.
    SFT         # Incomplete Fourier transform parameters
    𝒇           # The rfftfreqs for the values returned by the incomplete FT
    
    # These are used in the "key SIRS equations". See documenation folder.
    ℕ      # = Σ 𝓷 𝓷* (Elementwise sum over all available data)
    𝕃      # = Σ 𝓵 𝓵* ... 
    ℝ      # = Σ 𝓻 𝓻* ...
    𝕏      # = Σ 𝓷 𝓻* ...
    𝕐      # = Σ 𝓷 𝓵* ...
    ℤ      # = Σ 𝓻 𝓵* ...
    
    # These are used for computing reference row weights at f = 0 Hz
    R      # = Σ mean(refrows)
    N      # Number of frames included in sum. Some fail quality control.
    
    # Add the two resulting frequency dependent weight vectors.
    # alpha operates on the left-side reference columns. Beta
    # operates on the right-side reference columns.
    α
    β
    
    # For fitting ramps, we require the Legendre basis matrix
    # and its inverse
    L         # Each column is a naxis3 long Legendre polynomial
    Linv      # The Moore-Penrose inverse of B
    L_x_Linv  # = L * pinv(L)
    
    # This HxRG's good pixel map.
    gdpx
    
    # Masks for picking out regular pixels
    regpix_edge    # Mask for outputs containing reference columns
    regpix_middle  # Mask for all other outputs
        
    function SIRSCore(hxrg_kind::String, nout::Int64, nroh::Int64, τ::Float64,
                naxis3::Int64;gdpx=nothing, restore=false)
        
        # Definitions that depend on HxRG kind
        if hxrg_kind == "h4rg"
            naxis1 = 4096
            naxis2 = 4096
            xsize = naxis1 ÷ nout
            ysize = naxis2
        elseif hxrg_kind == "h2rg"
            naxis1 = 2048
            naxis2 = 2048
            xsize = naxis1 ÷ nout
            ysize = naxis2            
        else
           println("ERROR: hxrg_kind not supported!") 
        end
        
        # Definitions
        nleg   = 2             # Number of Legendre polynomials to fit to ramps
        n = (xsize+nroh)*ysize # Number of time steps per output
        
        # Declarations. These include the key SIRS sums
        ℕ = zeros(Float64, (naxis2+1,nout)) # These 3 are Fourier amplitudes
        𝕃 = zeros(Float64, (naxis2+1,nout)) # ...
        ℝ = zeros(Float64, (naxis2+1,nout)) # ...
        𝕏 = zeros(Complex{Float64}, (naxis2+1,nout)) # These 3 have phase info
        𝕐 = zeros(Complex{Float64}, (naxis2+1,nout)) # ...
        ℤ = zeros(Complex{Float64}, (naxis2+1,nout)) # ...
        
        # Declarations for reference rows
        R = zeros(Float64, (nout))
        N = zeros(Int64, (nout))
        
        # And the two resulting frequency dependent weight vectors.
        # alpha operates on the left-side reference columns. Beta
        # operates on the right-side reference columns.
        α = zeros(Complex{Float64}, (naxis2+1,nout))
        β = zeros(Complex{Float64}, (naxis2+1,nout))
        
        # Calculate pixel indices
        j = reshape(collect(1:n), (xsize+nroh,ysize))
        j_ref = reshape(j[1:HXRG_RB,:], :)      # Reference pixel indices
        j_reg = reshape(j[HXRG_RB+1:end,:], :)  # Regular pixel indices
        
        # Instantiate the incomplete Fourier transform for
        # the reference pixel streams. Keep
        # frequencies within about 1 / 2*τ_row
        # of DC. From looking at a few preliminary results,
        # there does not appear to be much going on near Nyquist.
        # Keep just a few to see if there is anything up there.
        𝒇 = Float64.(rfftfreq(n, τ^-1))  # All Fourier frequencies in Hz
        
        # We know that we are going to keep naxis/2+1 low frequencies and naxis2/2 high
        # frequencies. Just pick those out.
        lo_freq = collect(1:length(𝒇))[1:naxis2÷2+1]
        hi_freq = collect(1:length(𝒇))[end-naxis2÷2+1:end]
        k_ref = cat(lo_freq, hi_freq, dims=1)
        SFT = inc_rfft_pars(j_ref, k_ref, n, restore=restore)
        𝒇 = 𝒇[k_ref]
        
        # Build the Legendre basis matrix and its inverse.
        L = Array{Float64,2}(undef, (naxis3, nleg))
        for col in 1:nleg
            L[:,col] = legendre.(collect(-1:2/(naxis3-1):+1), col-1)
        end
        Linv = pinv(L)
        L_x_Linv = L * Linv
        
        # For SIRS, consider reference pixels to be operable. Do it so that we
        # don't overwrite the operability map
        if gdpx != nothing
            gdpx = copy(gdpx)
            gdpx[1:HXRG_RB,:]         .= 1 # Left
            gdpx[end-HXRG_RB+1:end,:] .= 1 # Right
            gdpx[:,1:HXRG_RB]         .= 1 # Bottom
            gdpx[:,end-HXRG_RB+1:end] .= 1 # Top
        else
            gdpx = BitArray{2}(ones(Bool, (naxis1,naxis2)))
        end
        
        #=
        Make masks that can be used for picking out just the regular pixels.
        These both operate on the data from one output
        with new row overheads trimmed. The software flips outputs to follow
        the HxRG clocking pattern, so we only need one mask for edge outputs.
        =#
        regpix_edge = BitArray{2}(zeros(Bool, (xsize,ysize))) # First & last outputs
        regpix_edge[HXRG_RB+1:end,HXRG_RB+1:end-HXRG_RB] .= true
        regpix_middle = BitArray{2}(zeros(Bool, (xsize,ysize))) # Middle outputs
        regpix_middle[:,HXRG_RB+1:end-HXRG_RB] .= true
        
        # Construct it
        new(hxrg_kind, naxis1, naxis2, naxis3, xsize, ysize, nout, nroh, n,
               τ, j_ref, j_reg, SFT, 𝒇, ℕ, 𝕃, ℝ, 𝕏, 𝕐, ℤ, R, N, α, β, L,
                    Linv, L_x_Linv, gdpx, regpix_edge, regpix_middle)
    end
    
end
