"""
    restore(savefile)

Restore a previously saved SIRS calibration file.

    Parameters: savefile:String
                  The name of a SIRS savefile
"""
function restore(savefile)
   
    # Load the save file
    d = load(savefile)
    
    # Pick out all parameters. We could also do this directly from
    # the d Dict, but it may be more readable to mirror what is in
    # SIRSCore
    hxrg_kind = d["hxrg_kind"]
    naxis1 = d["naxis1"]
    naxis2 = d["naxis2"]
    naxis3 = d["naxis3"]
    xsize = d["xsize"]
    ysize = d["ysize"]
    nout = d["nout"]
    nroh = d["nroh"]
    n = d["n"]
    τ = d["τ"]
    d_jref = d["j_ref"]
    d_jreg = d["j_reg"]
    SFT_j = d["SFT_j"]
    SFT_k = d["SFT_k"]
    SFT_n = d["SFT_n"]
    SFT_B = d["SFT_B"]
    SFT_Binv = d["SFT_Binv"]
    𝒇 = d["𝒇"]
    ℕ = d["ℕ"]
    𝕃 = d["𝕃"]
    ℝ = d["ℝ"]
    𝕏 = d["𝕏"]
    𝕐 = d["𝕐"]
    ℤ = d["ℤ"]
    R = d["R"]
    N = d["N"]
    α = d["α"]
    β = d["β"]
    L = d["L"]
    Linv = d["Linv"]
    L_x_Linv = d["L_x_Linv"]
    gdpx = d["gdpx"]
    regpix_edge = d["regpix_edge"]
    regpix_middle = d["regpix_middle"]
    
    # Create the SIRSCore
    sc = SIRSCore(hxrg_kind, nout, nroh, τ, naxis3, gdpx=gdpx, restore=true)
    
    # Copy in computed arrays. This is possible for arrays in an unmutable
    # struct, but not for constants and variables.
    sc.SFT.B .= SFT_B
    sc.SFT.Binv .= SFT_Binv
    sc.ℕ .= ℕ
    sc.𝕃 .= 𝕃
    sc.ℝ .= ℝ
    sc.𝕏 .= 𝕏
    sc.𝕐 .= 𝕐
    sc.ℤ .= ℤ
    sc.R .= R
    sc.N .= N
    sc.α .= α
    sc.β .= β
    sc.𝒇 .= 𝒇
    sc.L .= L
    sc.Linv .= Linv
    sc.L_x_Linv .= L_x_Linv
    sc.gdpx .= gdpx
    sc.regpix_edge .= regpix_edge
    sc.regpix_middle .= regpix_middle
    
    return(sc)
    
end
