"""
    export_to_sirspy(sc, file)

Export a SIRSCore for use by sirspy, the python-3 backend to SIRS.

    Parameters: sc::SIRSCore
                  A SIRSCore object
                file::String
                  The output filename. The file will be written in hdf5
                  format.

    Notes:
      * We do not export the entire SIRSCore. We only export the parameters
        that are used by sirspy.d
"""
function export_to_sirspy(sc::SIRSCore, file::String)
    h5open(file, "w") do file
        g = create_group(file, "SIRSCore") # create a group
        g["naxis1"] = sc.naxis1
        g["naxis2"] = sc.naxis2
        g["nout"] = sc.nout
        g["nroh"] = sc.nroh
        g["xsize"] = sc.xsize
        g["ysize"] = sc.ysize
        g["freq"] = sc.𝒇
        g["incft"] = copy(transpose(sc.SFT.Binv))  # copy(transpose()) converts to row major
        g["α"] = copy(transpose(sc.α))
        g["β"] = copy(transpose(sc.β))
    end
end
