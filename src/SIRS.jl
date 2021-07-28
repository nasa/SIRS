module SIRS

    using Dates          # Time functions
    using Dierckx        # Interpolation
    using FITSIO         # FITS files
    using FFTW           # FFTs
    using Glob           # Searching for files
    using HDF5           # Used when exporting for sirspy
    using Jacobi         # Legendre polynomials
    using JLD            # Save and load variables in Julia Data format (JLD)
                         #     JLD is a "dialect" of HDF5
    using LaTeXStrings   # Make nice labels in plots
    using LinearAlgebra  # Matrix maths
    using OMEinsum       # Tensor operations
    using Plots          # Plotting
    using Statistics     # Statistics
    using StatsBase      # Robust statistics
    using DirectConvolution # TBC, used in ref. col. quality control

    # Set defaults. These should be the same for any HxRG. Declare more specific
    # things in SIRSCore.jl
    const BZERO  = 4096                  # Offset when saving to UInt16 FITS files
    const HXRG_RB = 4                    # H4RG ref. pix. border width

    # Thresholds for trimming off outliers on top and on
    # bottom. These should be set to 1/2 the total amount to trim.
    const TRIM_HARD = .05/2              # Leaves 95%
    const TRIM_MED  = .01/2              # Leaves 99%
    const TRIM_SOFT = 3.16712e-5         # "4-sigma" clip

    # These values worked for flight candidates, but not for the JPL PPL data
    # const GD_OP_SIG_MIN =  2.0       # ...
    # const GD_OP_SIG_MAX =  14.0      # ...
    # 
    # These values worked for the JPL PPL and flight. The values are ad-hoc
    const GD_OP_SIG_MIN =  2.0   # Must be > 0
    const GD_OP_SIG_MAX =  28.0  # Must be > GD_OP_SIG_MIN. OK if big.

    export LegendreMatrices    # Legendre polynomial bases, used to fit and model
    export SIRSCore            # Core SIRS Structuture
    export inc_rfft_pars       # Incomplete Fourier transform parameters
    export adapt_sirssub       # SIRS subtraction for FITS files on ADAPT
    export chsuf               # Change a string's suffix
    export clear!              # Clear coadds from SIRSCore
    export coadd!              # Coadd a file into SIRSCore
    export export_to_sirspy    # Export SIRSCore parameters needed by sirspy
    export get_file_list       # Convenience function for data on ADAPT
    export legfit              # Fit Legendre polynomials
    export legval              # Model datacubes from Legendre polynomials
    export restore             # Load a SIRSCore from saved data
    export rowssol             # Solve for reference row weights (f = 0 Hz)
    export save                # Save a SIRSCore
    export sirssol             # SIRS solve for alpha and beta from files
    export sirssub!            # SIRS correct a data cube
    export solve!              # SIRS solve for alpha and beta from SIRSCore 
    export inc_rfft            # Incomplete Fourier transform
    export wplot               # Plot SIRS weights

    include("LegendreMatrices.jl")
    include("SIRSCore.jl")
    # include("inc_rfft_pars.jl")
    include("adapt_sirssub.jl")
    include("chsuf.jl")
    include("clear!.jl")
    include("coadd!.jl")
    include("export_to_sirspy.jl")
    include("get_file_list.jl")
    include("legfit.jl")
    include("legval.jl")
    include("restore.jl")
    include("rowssol.jl")
    include("save.jl")
    include("sirssol.jl")
    include("sirssub!.jl")
    include("solve!.jl")
    include("inc_rfft.jl")
    include("wplot.jl")

end
