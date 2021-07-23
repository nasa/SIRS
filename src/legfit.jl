"""
    legfit(lm::LegendreMatrices, D::Array{Float64,3})

Legendre fit datacubes

    Parameters: lm::LegendreMatrices
                  A LegendreMatrices object
                D::Array{Float64,3}
                  The datacube to fit. The dimensions are
                    1) x-image dimension
                    2) y-image dimension
                    3) Up-the-ramp frame index dimension

    Returns: fit
               The legendre fit
"""
function legfit(lm::LegendreMatrices, D::Array{Float64,3})
    return(ein"ij,klj->kli"(lm.pinv_B, D))    
end