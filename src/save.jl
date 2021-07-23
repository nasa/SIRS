"""
    save(SIRSCore, output_filename)

Save SIRSCore results to a JLD file. To ensure uniformity,
The filename is generated from information in the SIRSCore struct.

    Parameters: SIRSCore::SIRSCore
                  A coadded SIRSCore struct
                output_filename::String
                  The saveset goes here. It is required
                  to have the suffix .jld.
"""
function save(sc::SIRSCore, output_filename::String)
    
    # Check to be sure that the output filename's suffix is .jld
    if output_filename[end-3:end] != ".jld"
        println("Error: output_filename must have the suffix .jld")
        flush(stdout)
        return(-1)
    end
    
    # Get data from the SFT
    SFT_j = sc.SFT.j
    SFT_k = sc.SFT.k
    SFT_n = sc.SFT.n
    SFT_B = sc.SFT.B
    SFT_Binv = sc.SFT.Binv
    
    # Create a dictionary containing the information to save
    dict = Dict(
        "hxrg_kind" => sc.hxrg_kind,
        "naxis1" => sc.naxis1,
        "naxis2" => sc.naxis2,
        "naxis3" => sc.naxis3,
        "xsize" => sc.xsize,
        "ysize" => sc.ysize,
        "nout" => sc.nout,
        "nroh" => sc.nroh,
        "n" => sc.n,
        "τ" => sc.τ,
        "j_ref" => sc.j_ref,
        "j_reg" => sc.j_reg,
        "SFT_j" => SFT_j,
        "SFT_k" => SFT_k,
        "SFT_n" => SFT_n,
        "SFT_B" => SFT_B,
        "SFT_Binv" => SFT_Binv,
        "𝒇" => sc.𝒇,
        "ℕ" => sc.ℕ,
        "𝕃" => sc.𝕃,
        "ℝ" => sc.ℝ,
        "𝕏" => sc.𝕏,
        "𝕐" => sc.𝕐,
        "ℤ" => sc.ℤ,
        "R" => sc.R,
        "N" => sc.N,
        "α" => sc.α,
        "β" => sc.β,
        "L" => sc.L,
        "Linv" => sc.Linv,
        "L_x_Linv" => sc.L_x_Linv,
        "gdpx" => sc.gdpx,
        "regpix_edge" => sc.regpix_edge,
        "regpix_middle" => sc.regpix_middle
    )
    
    # Save results
    JLD.save(output_filename, dict)
    
end