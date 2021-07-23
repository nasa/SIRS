"""
    sirssol(input_filenames, vdetbias; [,gdpx_filename=nothing], 
                [,output_directory="./"] [,savesol=false] 
                    [,nroh=HXRG_NROH])

Solve for SIRS alpha and beta coefficients

    Parameters: input_filenames::Array{String,1}
                  List of input files
                vdetbias::Float64
                  Photodiode bias. This information is not
                  directly contained in the ADAPT FITS headers but
                  is required to calibrate science data.
                gdpx_filename::String
                  Name of good pixel map file
                output_directory::String (optional)
                  Name of the output directory
                savesol::Bool (optional)
                  Save the solution
                nroh::Int (optional)
                  New row overhead in pixels. The default is
                  given in SIRS.jl
"""
function sirssol(input_filenames::Array{String,1},
                    vdetbias::Float64; gdpx_filename=nothing, 
                        output_directory::String=nothing,
                            savesol::Bool=false, nroh::Int64=HXRG_NROH)
    
    # Declarations
    D = []    # Data
    sc = []   # SIRSCore struct
    gdpx = [] # Good pixel mask. 1=good and 0=bad.
    
    FIRST = true # A flag
    for file in input_filenames

        # Print status
        println("SIRSSOL Processing: ", file)
        flush(stdout)
        
        # Get data. The only supported file format is 
        # currently ADAPT. Do this because we need
        # the number of frames to instantiate SIRSCore.
        # All files must have the same number of frames.
        f = FITS(file, "r")
        D = Float64.(dropdims(read(f[2]), dims=4))
        D .*= -1 # Rectify DCL data so that charge integrates up
        H = read_header(f[1])
        date_beg = H["DATE-BEG"]
        detector = H["DETECTOR"]
        dewar    = H["DEWAR"]
        fpapos   = H["FPAPOS"]
        lodfile  = H["LODFILE"]
        origin   = H["ORIGIN"]
        scatemp  = H["SCATEMP"] 
        close(f)
        
        # Do initializations on first one
        if FIRST == true
            FIRST = false # Clear flag
            
            # Get good pixel map if given
            if gdpx_filename != nothing
                f = FITS(gdpx_filename, "r")
                gdpx = BitArray{2}(Bool.(read(f[2])))
                close(f)
            else
                # Otherwise assume all good
                gdpx = BitArray{2}(zeros(Bool, (NAXIS1,NAXIS2)))
                gdpx[HXRG_RB+1:end-HXRG_RB,HXRG_RB+1:end-HXRG_RB] .= 1
            end
            
            # Instantiate SIRSCore
            coadd = size(D)[3] * length(input_filenames) # Number of coadded frames
            sc = SIRSCore(size(D)[3], gdpx, nroh=nroh, coadd=coadd, date_beg=date_beg,
                    detector=detector, dewar=dewar, files=input_filenames, 
                    fpapos=fpapos, lodfile=lodfile, origin=origin, scatemp=scatemp,
                    vdetbias=vdetbias)
        end
        
        # Coadd data
        coadd!(sc, D)
        
    end
    
    # Solve for alpha and beta
    solve!(sc)
    
    # Optionally save output file
    if savesol == true
        save(sc, output_directory=output_directory)
    end
    
    # Done
    # return(sc)
    
end