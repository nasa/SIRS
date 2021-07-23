"""
    solve(sc::SIRSCore)

Solve a SIRCore for the SIRS alpha and beta arrays

    Parameters: sc::SIRSCore
                  A SIRCore containing coadded SIRS sums.
"""
function solve!(sc::SIRSCore)
    sc.α .= (sc.𝕐 .- sc.𝕏 .* sc.ℤ ./ sc.ℝ) ./
                (sc.𝕃 .-  sc.ℤ .* conj(sc.ℤ) ./ sc.ℝ)
    sc.β .= (sc.𝕏 .* sc.𝕃 ./ sc.ℝ .- sc.𝕐 .* conj(sc.ℤ) ./ sc.ℝ) ./
                (sc.𝕃 .-  sc.ℤ .* conj(sc.ℤ) ./ sc.ℝ)
end