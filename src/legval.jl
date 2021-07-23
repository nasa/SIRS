"""
    legval(lm::LegendreMatrices, Λ::Array{Float64,3})

Given a Legendre fit, build a model of the up-the-ramp sampled data

    Parameters: lm::LegendreMatrices
                  A LegendreMatrices object
                Λ::Array{Float64,3}
                  A Legendre fit
                    1) x-image dimension
                    2) y-image dimension
                    3) Legendre coefficient dimension

    Returns: D
               A model of the up-the-ramp sampled data
"""
function legval(lm::LegendreMatrices, Λ::Array{Float64,3})
    return(ein"ij,klj->kli"(lm.B, Λ))    
end