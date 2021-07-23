"""
    adapt_sirssub(fits_filenames, sc_filename, odir; [verbose=true])

Simplified Improved Reference Subtraction for Roman data in the ADAPT data format.

    Parameters: fits_filenames::Array{String,1}
                  Input FITS filenames
                sc_filename::String
                  SIRSCore filename
                odir::String
                  Output directory name. Output filenames are created 
                  automatically.
                verbose::Bool
                  Print diagnostic information

Modification History:
    3/10/21 BJR
      * In DCL charge integrates down. Multiply by -1x to integrate up.
"""
function adapt_sirssub(fits_filenames::Array{String,1},
                sc_filename::String, odir::String; verbose::Bool=true)
    
    # Ensure that odir terminates with a '/' character
    if odir[end] != '/'; odir = odir * '/'; end
    
    # Read in the SIRSCore
    sc = restore(sc_filename)
    
    # Reference correct files
    for file in fits_filenames
        
        # Print diagnostic information
        if verbose == true
            println("adapt_sirsub() processing file: ", file)
            flush(stdout)
        end
        
        # Read a file
        f = FITS(file, "r")
        H = read_header(f[1])
        D = Float64.(read(f[2])[:,:,:,1]) # Use for HyC data on ADAPT
        D .*= -1.0 # Make charge integrate up
        close(f)
        
        # Augment the header
        H["REFCOR"] = true
        set_comment!(H, "REFCOR", "Two-stream reference corrected by SIRS")
        H["SIRSBIAS"] = BZERO
        set_comment!(H, "SIRSBIAS", "DN; Constant added by SIRS to make values positive")
        
        # Reference correct in place
        sirssub!(sc, D)
        
        # Make output filename
        base_filename = basename(file)
        fout = base_filename[1:findlast(".fits", base_filename)[1]] * "sirs.fits"
        
        # Write result
        f = FITS(odir * fout, "w")
        D = Int64.(round.(D))    # Convert to integer
        D .+= BZERO              # Offset to make all the good stuff positive
        D[D .< 0] .= 0           # Anything less than zero is bad
        D[D .> 2^16-1] .= 2^16-1 # Anything over this is likewise bad
        write(f, UInt16.(D), header=H)
        close(f)
        
    end
end