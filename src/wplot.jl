"""
    wplot(sc; [,frange_low] [,frange_high] [,mag=1.0])

Plot frequency dependent weights

    Parameters: sc::SIRSCore
                  A SIRSCore containing frequency dependent weights
                op::Int64
                  Output to plot
                frange_low::Tuple
                  Range of low frequencies to plot
                frange_high::Tuple
                  Range of high frequencies to plot
                mag::Float64
                  Plot magnification
                nice::Bool
                  Use "nice" ticks on x-axis
"""
function wplot(sc::SIRSCore, op::Int64; frange_low=nothing, frange_high=nothing,
                    mag=1.0, nice=false)
    
    # Definitions
    τ_row = (sc.xsize+sc.nroh) * sc.τ # Row readout time
    f_ny_row = 1/2/τ_row # Nyquist frequeny for row sampling
    f_ny_pix = 1/sc.τ/2 # Nyquist frequency for the pixel sampling
    ylim_amp = (0.0 - .05, 1.0) # Y-limits for amplitude plots
    ylim_phs = (-pi,+pi)     # Y-limits for phase plots
    
    # Set defaults
    if frange_low == nothing; frange_low = (0., f_ny_row); end
    if frange_high == nothing; frange_high = (-f_ny_row, 0.); end
        
    # Make the plots
    p1 = plot(sc.𝒇, [abs.(sc.α)[:,op], abs.(sc.β)[:,op]], alpha=1,
        xlims=frange_low,
        ylims=ylim_amp,
        legend=:true,
        label = [L"\alpha" L"\beta"],
        titlefontsize=10, title = "Low Frequency",
        annotate=[(20,.95, text("Op: "*string(op), :left, 7, "Computer Modern"))])
    if nice==true
        p1 = plot!(xticks=60:120:frange_low[2], xminorticks=2,
        xminorgrid=:true)
    end


    p2 = plot(sc.𝒇, [angle.(sc.α)[:,op], angle.(sc.β)[:,op]], alpha=.7,
        seriestype=:scatter, ms=1.5, msw=0, ma=.5,
        xlims=frange_low,
        ylim = ylim_phs,
        legend=:false)
    if nice==true
        p2 = plot!(xticks=60:120:frange_low[2], xminorticks=2,
        xminorgrid=:true)
    end

    
    p3 = plot(sc.𝒇 .- f_ny_pix, [abs.(sc.α)[:,op], abs.(sc.β)[:,op]], alpha=1,
        xlims=frange_high,
        ylims=ylim_amp,
        legend=:false,
        titlefontsize=10, title = "High Frequency")
    
    p4 = plot(sc.𝒇 .- f_ny_pix, [angle.(sc.α)[:,op], angle.(sc.β)[:,op]],
        seriestype=:scatter, ms=1.5, msw=0, ma=.5,
        xlims=frange_high,
        ylim = ylim_phs,
        legend=:false)

    left_col  = plot(p1, p2, layout=(2,1), linked=:x,
                    xlabel=["" L"\textrm{Frequency (Hz)}"],
                    ylabel=[L"\textrm{Amplitude}" L"\textrm{Phase}"])
    right_col = plot(p3, p4, layout=(2,1), linked=:x,
                    xlabel=["" L"\textrm{Frequency} ~(+f_\textrm{Ny}~ \textrm{Hz})"])
    
    figure = plot(left_col, right_col, layout=(1,2), fontfamily="Computer Modern",
                size=(mag*640,mag*480))
    
    return(figure)
    
end
