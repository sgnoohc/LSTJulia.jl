using LSTJulia
using Arrow
using PlotlyJS
using Query
using DataFrames
using PlotlyJSWrapper

function getangledata(jobidx)
    # table = Arrow.Table(string("muon_gun_segmentpairangles",jobidx,".arrow"))
    table = Arrow.Table(string("muonhighe_gun_segmentpairangles",jobidx,".arrow"))
    df = DataFrame(table)
    as = []
    for item in 1:4
        angledata = df |> @filter(_.moddifftype == item) |> DataFrame
        push!(as, angledata)
    end
    return as
end

function makefigure(jobidx)
    table = Arrow.Table(string("muon_gun_segmentpairangles",jobidx,".arrow"))
    # table = Arrow.Table(string("muonhighe_gun_segmentpairangles",jobidx,".arrow"))
    df = DataFrame(table)
    for item in 1:4
        angledata = df |> @filter(_.moddifftype == item) |> DataFrame
        trace = histogram2d(x=abs.(angledata.Δϕ),
                            y=abs.(angledata.Δθ),
                            opacity=0.75,
                            xbins=attr(start=0, size=0.001),
                            xbins_end=0.02,
                            ybins=attr(start=0, size=0.001),
                            ybins_end=0.02,
                           )
        p = plot(trace)
        savefig(p, string("plots/moddifftype",item,"_",jobidx,".pdf"))
        savefig(p, string("plots/moddifftype",item,"_",jobidx,".png"))
        trace = histogram2d(x=abs.(angledata.Δϕ),
                            y=abs.(angledata.Δθ),
                            opacity=0.75,
                            xbins=attr(start=0, size=0.001),
                            xbins_end=0.02,
                            ybins=attr(start=0, size=0.001),
                            ybins_end=0.50,
                           )
        p = plot(trace)
        savefig(p, string("plots/moddifftypewide",item,"_",jobidx,".pdf"))
        savefig(p, string("plots/moddifftypewide",item,"_",jobidx,".png"))

        h = PlotlyJSWrapper.make_fhist1d(abs.(angledata.Δϕ), 0:0.0005:0.02)
        trace = PlotlyJSWrapper.build_hist1dtrace(h; witherror=true)
        # trace.fields[:marker][:color] = "blue"
        # trace.fields[:opacity] = 0.75
        p = plot(trace, Layout(title="Δϕ distribution diffmode=$item", xaxis_title="Δϕ [rad]", yaxis_title="Events", bargap=0, template="simple_white"))
        savefig(p, string("plots/moddifftypewide",item,"_",jobidx,"_dphi.pdf"))
        savefig(p, string("plots/moddifftypewide",item,"_",jobidx,"_dphi.png"))

        trace = PlotlyJSWrapper.build_hist1dtrace(PlotlyJSWrapper.make_fhist1d(abs.(angledata.Δθ), 0:0.0005:0.02), witherror=true)
        # trace.fields[:marker][:color] = "blue"
        # trace.fields[:opacity] = 0.75
        p = plot(trace, Layout(title="Δrz distribution diffmode=$item", xaxis_title="Δrz [rad]", yaxis_title="Events", bargap=0, template="simple_white"))
        savefig(p, string("plots/moddifftypewide",item,"_",jobidx,"_drz.pdf"))
        savefig(p, string("plots/moddifftypewide",item,"_",jobidx,"_drz.png"))
    end
end

function make_fishbone_result_figure(idx)

    # table = Arrow.Table("fishbone_result.arrow")
    # table = Arrow.Table("fishbonehighe_result$idx.arrow")
    table = Arrow.Table("fishbone_result$idx.arrow")
    df = DataFrame(table)

    trace1 = PlotlyJSWrapper.make_fhist1d(df.nsims, 0:1:20) |> PlotlyJSWrapper.build_hist1dtrace
    trace1.fields[:marker][:color] = "orange"
    trace1.fields[:opacity] = 0.5
    trace2 = PlotlyJSWrapper.make_fhist1d(df.nsims_fishboned, 0:1:20) |> PlotlyJSWrapper.build_hist1dtrace
    trace2.fields[:marker][:color] = "blue"
    trace2.fields[:opacity] = 0.5
    p = plot([trace1, trace2], Layout(barmode="overlay", bargap=0, template="simple_white"))
    savefig(p, "plots/nsims.pdf")
    savefig(p, "plots/nsims.png")

    trace1 = PlotlyJSWrapper.make_fhist1d(df.ndups, 0:1:20) |> PlotlyJSWrapper.build_hist1dtrace
    trace1.fields[:marker][:color] = "orange"
    trace1.fields[:opacity] = 0.5
    trace2 = PlotlyJSWrapper.make_fhist1d(df.ndups_fishboned, 0:1:20) |> PlotlyJSWrapper.build_hist1dtrace
    trace2.fields[:marker][:color] = "blue"
    trace2.fields[:opacity] = 0.5
    p = plot([trace1, trace2], Layout(barmode="overlay", bargap=0, template="simple_white"))
    savefig(p, "plots/ndups.pdf")
    savefig(p, "plots/ndups.png")

    trace1 = PlotlyJSWrapper.make_fhist1d(df.ndups./df.nsims, 0:0.05:1) |> PlotlyJSWrapper.build_hist1dtrace
    trace1.fields[:marker][:color] = "orange"
    trace1.fields[:opacity] = 0.5
    trace2 = PlotlyJSWrapper.make_fhist1d(df.ndups_fishboned./df.nsims_fishboned, 0:0.05:1) |> PlotlyJSWrapper.build_hist1dtrace
    trace2.fields[:marker][:color] = "blue"
    trace2.fields[:opacity] = 0.5
    p = plot([trace1, trace2], Layout(barmode="overlay", bargap=0, template="simple_white"))
    savefig(p, "plots/ndups_fraction.pdf")
    savefig(p, "plots/ndups_fraction.png")

end

function main()
    @time makefigure("12")
    # @time makefigure("23")
    # @time makefigure("34")
    # @time makefigure("45")
    # @time makefigure("56")
    @time make_fishbone_result_figure("12")
    # @time make_fishbone_result_figure("23")
    # @time make_fishbone_result_figure("34")
    # @time make_fishbone_result_figure("45")
    # @time make_fishbone_result_figure("56")
end

main()
