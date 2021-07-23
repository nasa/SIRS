"""
    clear(sc::SIRSCore)

Reset SIRS sums to zero.

    parameters: None
"""
function clear!(sc::SIRSCore)
    sc.ℕ .*= 0
    sc.𝕃 .*= 0
    sc.ℝ .*= 0
    sc.𝕏 .*= 0
    sc.𝕐 .*= 0
    sc.ℤ .*= 0
    sc.R .*= 0
    sc.N .*= 0
end